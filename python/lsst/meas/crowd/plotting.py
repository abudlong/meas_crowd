# This file is part of meas_crowd.
#
# Developed for the LSST Data Management System.
# See the COPYRIGHT file at the top-level directory of this distribution.

"""Diagnostic plots of a finished run: data | model | residual, residual/sigma, star zooms,
cross-sections, bright-star stacks and comparison tables.
"""

import numpy as np
from scipy.spatial import cKDTree
from . import fitLib as L
from .pipeline import (
    BAD_PIXEL_COLOR,
    EXCLUDED_COLOR,
    EXCLUDED_COLOR_GRAY,
    NEARSAT_COLOR,
    _params,
    _residual_stack,
    _satlike,
    sat_zones,
    saturation_level,
    select_psf_stars,
    unmodeled_near_mask,
)


def show_bright_star_stack(run, min_snr=200.0, half=20, plot=True):
    """Median of the bright, clean, unsaturated stars' residuals, each normalised by the star's
    flux and centred on its sub-pixel position -- once with the PSF-only joint model, once
    with the final model (incl. per-star corrections).  A pattern that survives the median is
    COMMON to all bright stars -> a PSF feature (e.g. diffraction spikes) missing from the PSF
    model; noise and star-specific errors average out.  Also shows the PSF kernel (log).
    Returns dict(psf_only=..., final=..., n=...)."""
    final = run.outer_history[-1]
    exp, res = final["exposure"], final["result"]
    S, fit = res["S"], res["snaps"][-1]["fit"]
    bb = exp.getBBox()
    c = bb.getCenter()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    data, var = exp.image.array, exp.variance.array
    psf = exp.getPsf()
    kern = psf.computeKernelImage(c).array
    fwhm = L._fwhm(psf, c)
    noise = float(np.sqrt(np.nanmedian(var)))
    P = dict(_params(run.config))
    P["psf_stack_min_snr"] = min_snr
    P["psf_stack_n"] = 60
    idx = select_psf_stars(exp, res, fwhm, kern.max(), noise, P, verbose=False)
    out = {}
    for name, model in (("psf_only", fit["model_full"]), ("final", res["model_full"])):
        out[name] = _residual_stack(data - model, fit["good"], S["xy"], fit["amps"], idx,
                                    X0, Y0, half)
    out["n"] = len(idx)
    if len(idx) >= 6:
        # brighter-fatter test: is the shape error the same for the brighter and fainter half?
        u = np.arange(-half, half + 1)
        R = np.hypot(*np.meshgrid(u, u))
        order = idx[np.argsort(-fit["amps"][idx])]
        halves = (("brighter half", order[:len(order) // 2]), ("fainter half", order[len(order) // 2:]))
        print("  core shape error by brightness (PSF-only model; median residual flux as % of the "
              "star's flux):")
        print(f"    {'':14s} {'median flux':>12s} {'r < 1 FWHM':>11s} {'1-2.5 FWHM':>11s}")
        rows = []
        for lab, sel in halves:
            st = _residual_stack(data - fit["model_full"], fit["good"], S["xy"], fit["amps"], sel,
                                 X0, Y0, half)
            cin = 100 * np.nansum(st[R < fwhm])
            cout = 100 * np.nansum(st[(R >= fwhm) & (R < 2.5 * fwhm)])
            rows.append((cin, cout))
            print(f"    {lab:14s} {np.median(fit['amps'][sel]):12.3g} {cin:+10.2f}% {cout:+10.2f}%")
        out["by_brightness"] = rows
        if np.sign(rows[0][0]) != np.sign(rows[1][0]) and min(abs(rows[0][0]), abs(rows[1][0])) > 0.2:
            print("    -> opposite signs: the PSF shape DEPENDS ON BRIGHTNESS (brighter-fatter);"
                  " one PSF cannot fit both")
    print(f"stacked {len(idx)} bright clean star(s) (peak S/N >= {min_snr:g})")
    _print_run_version(final)
    if plot and len(idx):
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        h = kern.shape[0] // 2
        ks = kern[max(h - half, 0):h + half + 1, max(h - half, 0):h + half + 1]
        lim = np.nanpercentile(np.abs(out["psf_only"]), 99) or 1e-6
        fig, ax = plt.subplots(1, 3, figsize=(17, 5.5))
        im = ax[0].imshow(np.clip(ks, ks.max() * 1e-5, None), origin="lower", cmap="magma",
                          norm=LogNorm(), extent=(-half - .5, half + .5, -half - .5, half + .5))
        plt.colorbar(im, ax=ax[0], fraction=0.046)
        ax[0].set_title("PSF model kernel (log)")
        for a, key, ttl in ((ax[1], "psf_only", "median residual / flux: PSF-only model"),
                            (ax[2], "final", "median residual / flux: final model")):
            im = a.imshow(out[key], origin="lower", cmap="RdBu_r", vmin=-lim, vmax=lim,
                          extent=(-half - .5, half + .5, -half - .5, half + .5))
            plt.colorbar(im, ax=a, fraction=0.046)
            a.set_title(ttl)
        for a in ax:
            a.set_xlabel("dx [px]")
        ax[0].set_ylabel("dy [px]")
        fig.suptitle(f"{len(idx)} bright clean stars: a pattern here is shared by all of them "
                     f"(a PSF feature); noise and per-star errors average out", fontsize=11)
        plt.tight_layout()
        plt.show()
    return out


def show_star_zooms(run, stars=None, n=8, half=12, sigma="model", grayscale=False,
                    min_snr=30.0, cuts=False):
    """Zoom on individual stars: data | model | residual | residual/sigma_eff, with every source
    in the window marked (o = Gaia-seeded, x = found in the residual) and, for Gaia sources,
    a line from the Gaia catalogue position (+) to the fitted position.
    stars: list of (x, y) detector positions, or None = the n worst stars of
    star_residual_table (highest chi2).  Prints each star's neighbours (offset, flux ratio,
    Gaia or not, distance from its Gaia position)."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    final = run.outer_history[-1]
    exp, res = final["exposure"], final["result"]
    S, fit = res["S"], res["snaps"][-1]["fit"]
    bb = exp.getBBox()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    data, model = exp.image.array, res["model_full"]
    m = np.clip(model, 0, None)
    var = exp.variance.array.copy()
    if sigma in ("effective", "model"):
        var = var + float(final.get("poisson_k", 0.0) or 0.0) * m
    if sigma == "model":
        var = var + (float(final.get("model_err", 0.0) or 0.0) * m) ** 2
    bad = ((exp.mask.array & run.bad_bitmask) != 0) | ~np.isfinite(data - model) | ~(var > 0)
    try:
        nearsat = (exp.mask.array & exp.mask.getPlaneBitMask("NEARSAT")) != 0
    except Exception:
        nearsat = np.zeros(bad.shape, bool)
    resid = data - model
    z = resid / np.sqrt(np.where(var > 0, var, np.nan))
    amps = fit["amps"]
    xy = S["xy"]
    targets = _pick_stars(run, stars, n, min_snr)
    if not targets:
        print("no stars to show")
        return None
    tree = cKDTree(xy)
    cm_res = plt.get_cmap("gray" if grayscale else "RdBu_r").copy()
    fig, axes = plt.subplots(len(targets), 4, figsize=(15, 3.6 * len(targets)), squeeze=False)
    for row, (j, chi) in zip(axes, targets):
        x, y = xy[j]
        xa, xb = max(int(round(x)) - half - X0, 0), min(int(round(x)) + half + 1 - X0, data.shape[1])
        ya, yb = max(int(round(y)) - half - Y0, 0), min(int(round(y)) + half + 1 - Y0, data.shape[0])
        ext = (X0 + xa - .5, X0 + xb - .5, Y0 + ya - .5, Y0 + yb - .5)
        sl = (slice(ya, yb), slice(xa, xb))
        d = data[sl]
        lo, hi = np.nanpercentile(d, 1), np.nanmax(d)
        from matplotlib.colors import SymLogNorm
        nrm = SymLogNorm(linthresh=max(3 * np.sqrt(np.nanmedian(exp.variance.array)), 1e-6),
                         vmin=lo, vmax=hi)
        rs = np.nanpercentile(np.abs(resid[sl][~bad[sl]]), 99) if (~bad[sl]).any() else 1.0
        imgs = [(d, dict(cmap=L.CMAP, norm=nrm), "data"),
                (model[sl], dict(cmap=L.CMAP, norm=nrm), "model"),
                (np.ma.masked_where(bad[sl], resid[sl]), dict(cmap=cm_res, vmin=-rs, vmax=rs),
                 "data − model"),
                (np.ma.masked_where(bad[sl], z[sl]), dict(cmap=cm_res, vmin=-5, vmax=5),
                 "residual / σ_eff")]
        near = [k for k in tree.query_ball_point((x, y), half * 1.42)]
        for c, (ax, (img, kw, ttl)) in enumerate(zip(row, imgs)):
            im = ax.imshow(img, origin="lower", extent=ext, interpolation="nearest", **kw)
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
            for msk, col in ((bad[sl] & ~nearsat[sl], BAD_PIXEL_COLOR), (nearsat[sl], NEARSAT_COLOR)):
                if msk.any():
                    ax.imshow(np.ma.masked_where(~msk, np.ones(msk.shape)), origin="lower",
                              extent=ext, cmap=ListedColormap([col]), interpolation="nearest")
            mc = "w" if c < 2 else ("#00b4d8" if grayscale else "k")
            for k in near:
                px, py = xy[k]
                if S["seed"][k]:
                    ax.plot(px, py, "o", mfc="none", mec=mc, ms=9 if k == j else 6, mew=1.2)
                    rx, ry = S["ref"][k]
                    if np.isfinite(rx):
                        ax.plot([rx, px], [ry, py], "-", color=mc, lw=0.8)
                        ax.plot(rx, ry, "+", color=mc, ms=6)
                else:
                    ax.plot(px, py, "x", color="#ff3b3b" if c < 2 or grayscale else "#d62728",
                            ms=7, mew=1.5)
            ax.set_xlim(ext[0], ext[1])
            ax.set_ylim(ext[2], ext[3])
            ax.set_title(ttl, fontsize=9)
        row[0].set_ylabel(f"({x:.0f}, {y:.0f})\nflux {amps[j]:.3g}" +
                          (f"\nchi2 {chi:.2f}" if np.isfinite(chi) else ""), fontsize=9)
        print(f"star ({x:.1f}, {y:.1f})  flux {amps[j]:.3g}  {'Gaia' if S['seed'][j] else 'non-Gaia'}"
              + (f"  chi2 {chi:.2f}" if np.isfinite(chi) else "")
              + (f"  [{np.hypot(*(xy[j] - S['ref'][j])):.2f}px from its Gaia position]"
                 if S["seed"][j] and np.isfinite(S["ref"][j, 0]) else ""))
        for k in sorted(near, key=lambda k: np.hypot(*(xy[k] - xy[j]))):
            if k == j:
                continue
            dx, dy = xy[k] - xy[j]
            off = (f", {np.hypot(*(xy[k] - S['ref'][k])):.2f}px from its Gaia position"
                   if S["seed"][k] and np.isfinite(S["ref"][k, 0]) else "")
            print(f"    neighbour dx {dx:+6.1f} dy {dy:+6.1f}  flux ratio {amps[k] / max(amps[j], 1e-30):6.3f}"
                  f"  {'Gaia' if S['seed'][k] else 'found in residual'}{off}")
    fig.suptitle("o = Gaia-seeded source (+ = its Gaia position, line = fitted shift), "
                 "x = source found in the residual; magenta = bad, orange = near-saturated core",
                 fontsize=10)
    plt.tight_layout()
    plt.show()
    if cuts:
        show_cross_sections(run, stars=[tuple(xy[j]) for j, _ in targets], half=half, sigma=sigma)
    return [j for j, _ in targets]


def _pick_stars(run, stars, n, min_snr):
    """[(index, chi2)] -- the n worst stars of star_residual_table, or the sources nearest to the
    given (x, y) positions."""
    import contextlib
    import io
    S = run.outer_history[-1]["result"]["S"]
    if stars is None:
        with contextlib.redirect_stdout(io.StringIO()):
            t = star_residual_table(run, min_snr=min_snr, plot=False)
        if t is None:
            return []
        o = np.argsort(-t["chi2_final"])[:n]
        return [(int(t["idx"][i]), float(t["chi2_final"][i])) for i in o]
    tree = cKDTree(S["xy"])
    return [(int(tree.query(p)[1]), np.nan) for p in stars]


def show_cross_sections(run, stars=None, n=4, half=15, width=1, sigma="model", min_snr=30.0,
                        log=True):
    """Cuts through stars: for each star, a cut along x (through its row), along y (through its
    column) and the azimuthally averaged RADIAL profile.  Top: data (points) and model (line);
    bottom: residual = data - model with its +-1 sigma_eff band (grey) -- a ring shows up as a
    bump in the radial residual at its radius.  Other sources lying in the cut are marked
    (dotted: Gaia, dashed: found in the residual).
    stars: list of (x, y), or None = the n worst stars of star_residual_table.
    width: number of rows / columns averaged for the x / y cuts (odd; 1 = single pixel).
    log: data/model on a symmetric-log axis (wings visible) instead of linear."""
    import matplotlib.pyplot as plt
    final = run.outer_history[-1]
    exp, res = final["exposure"], final["result"]
    S, fit = res["S"], res["snaps"][-1]["fit"]
    bb = exp.getBBox()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    data, model = exp.image.array, res["model_full"]
    m = np.clip(model, 0, None)
    var = exp.variance.array.copy()
    if sigma in ("effective", "model"):
        var = var + float(final.get("poisson_k", 0.0) or 0.0) * m
    if sigma == "model":
        var = var + (float(final.get("model_err", 0.0) or 0.0) * m) ** 2
    bad = ((exp.mask.array & run.bad_bitmask) != 0) | ~np.isfinite(data - model) | ~(var > 0)
    noise = float(np.sqrt(np.nanmedian(exp.variance.array)))
    targets = _pick_stars(run, stars, n, min_snr)
    if not targets:
        print("no stars to show")
        return None
    xy, amps = S["xy"], fit["amps"]
    w2 = max(int(width) // 2, 0)
    fig = plt.figure(figsize=(17, 4.6 * len(targets)))
    outer = fig.add_gridspec(len(targets), 1, hspace=0.32)
    out = []
    for row, (j, chi) in enumerate(targets):
        gs = outer[row].subgridspec(2, 3, height_ratios=[2, 1], hspace=0.06, wspace=0.25)
        x, y = xy[j]
        ix, iy = int(round(x)) - X0, int(round(y)) - Y0
        cuts = {}
        for name in ("x", "y"):
            if name == "x":
                lo, hi = max(ix - half, 0), min(ix + half + 1, data.shape[1])
                rows = slice(max(iy - w2, 0), iy + w2 + 1)
                D, Mo, V, B = (a[rows, lo:hi] for a in (data, model, var, bad))
                coord = X0 + np.arange(lo, hi) - x
                other = [(k, xy[k, 0] - x) for k in range(len(xy))
                         if k != j and abs(xy[k, 1] - y) <= w2 + 1.0 and abs(xy[k, 0] - x) <= half]
            else:
                lo, hi = max(iy - half, 0), min(iy + half + 1, data.shape[0])
                cols = slice(max(ix - w2, 0), ix + w2 + 1)
                D, Mo, V, B = (a[lo:hi, cols].T for a in (data, model, var, bad))
                coord = Y0 + np.arange(lo, hi) - y
                other = [(k, xy[k, 1] - y) for k in range(len(xy))
                         if k != j and abs(xy[k, 0] - x) <= w2 + 1.0 and abs(xy[k, 1] - y) <= half]
            ok = ~B
            nn = ok.sum(axis=0)
            with np.errstate(invalid="ignore", divide="ignore"):
                d = np.where(nn > 0, np.sum(np.where(ok, D, 0), axis=0) / nn, np.nan)
                mo = np.where(nn > 0, np.sum(np.where(ok, Mo, 0), axis=0) / nn, np.nan)
                sg = np.where(nn > 0, np.sqrt(np.sum(np.where(ok, V, 0), axis=0)) / nn, np.nan)
            cuts[name] = (coord, d, mo, sg, other)
        # radial profile
        R = half
        xa, xb = max(ix - R, 0), min(ix + R + 1, data.shape[1])
        ya, yb = max(iy - R, 0), min(iy + R + 1, data.shape[0])
        yy, xx = np.mgrid[ya + Y0:yb + Y0, xa + X0:xb + X0]
        rr = np.hypot(xx - x, yy - y)
        ok = ~bad[ya:yb, xa:xb]
        edges = np.arange(0, R + 1.0, 1.0)
        rc, rd, rm, rs = [], [], [], []
        for a_, b_ in zip(edges[:-1], edges[1:]):
            k_ = ok & (rr >= a_) & (rr < b_)
            if k_.sum() < 3:
                continue
            rc.append(np.mean(rr[k_]))
            rd.append(np.mean(data[ya:yb, xa:xb][k_]))
            rm.append(np.mean(model[ya:yb, xa:xb][k_]))
            rs.append(np.sqrt(np.sum(var[ya:yb, xa:xb][k_])) / k_.sum())
        cuts["r"] = (np.array(rc), np.array(rd), np.array(rm), np.array(rs), [])
        for col, (name, ttl) in enumerate((("x", "cut along x (row)"), ("y", "cut along y (column)"),
                                           ("r", "radial profile (azimuthal mean)"))):
            coord, d, mo, sg, other = cuts[name]
            a1 = fig.add_subplot(gs[0, col])
            a2 = fig.add_subplot(gs[1, col], sharex=a1)
            a1.plot(coord, d, "o", ms=3.5, color="k", label="data")
            a1.plot(coord, mo, "-", lw=1.5, color="C1", label="model")
            if log:
                a1.set_yscale("symlog", linthresh=max(3 * noise, 1e-6))
            r_ = d - mo
            a2.fill_between(coord, -sg, sg, color="0.85", step="mid", label="±1σ_eff")
            a2.axhline(0, color="k", lw=0.7)
            a2.plot(coord, r_, "o-", ms=3, lw=0.8, color="C3", label="data − model")
            good_ = np.isfinite(r_) & (sg > 0)
            chi_c = float(np.mean((r_[good_] / sg[good_]) ** 2)) if good_.any() else np.nan
            a2.text(0.01, 0.95, f"χ²/pt = {chi_c:.2f}", transform=a2.transAxes, va="top", fontsize=8)
            for k, off in other:
                ls = ":" if S["seed"][k] else "--"
                for a in (a1, a2):
                    a.axvline(off, color="C0" if S["seed"][k] else "C3", ls=ls, lw=0.9, alpha=0.8)
            a1.set_title(ttl + (f"  (width {2 * w2 + 1} px)" if name != "r" and w2 else ""), fontsize=9)
            a1.tick_params(labelbottom=False)
            a2.set_xlabel("offset from star centre [px]" if name != "r" else "radius [px]", fontsize=8)
            if col == 0:
                a1.set_ylabel(f"({x:.0f}, {y:.0f})  flux {amps[j]:.3g}" +
                              (f"\nchi2 {chi:.2f}" if np.isfinite(chi) else "") + "\ncounts", fontsize=8)
                a2.set_ylabel("residual", fontsize=8)
            if row == 0 and col == 2:
                a1.legend(fontsize=7, loc="upper right")
                a2.legend(fontsize=7, loc="upper right")
        out.append(dict(index=j, cuts=cuts))
    fig.suptitle("Cross-sections: data (points), model (line), residual with ±1σ_eff band; vertical "
                 "lines = other sources in the cut (dotted Gaia, dashed found in residual)",
                 fontsize=10, y=1.0)
    plt.show()
    return out


def compare_runs(variants, min_snr=200.0):
    """Side-by-side numbers for several runs of the same cutouts, e.g.
        compare_runs({"default": runs, "wider wing join": runs_b})
    For each cutout and variant: residual/sigma_eff (robust), sources, and the stacked
    PSF-only residual of the bright stars (median residual/flux x 1e5) in radial bins --
    a ring shows as a large value in one bin; closer to 0 everywhere = better PSF."""
    import contextlib
    import io
    from . import session as cn
    names = sorted({k for runs in variants.values() for k in runs})
    edges = (0, 4, 8, 12, 16, 20)
    for nm in names:
        print(f"===== {nm}")
        print(f"  {'variant':22s} {'res/σ_eff':>9s} {'N src':>6s}   stacked PSF-only residual/flux x1e5 "
              f"in r [px] = " + "  ".join(f"{a}-{b}" for a, b in zip(edges[:-1], edges[1:])))
        for vname, runs in variants.items():
            run = runs.get(nm)
            if run is None:
                continue
            final = run.outer_history[-1]
            exp, res = final["exposure"], final["result"]
            data, model = exp.image.array, res["model_full"]
            m = np.clip(model, 0, None)
            var = (exp.variance.array + float(final.get("poisson_k", 0) or 0) * m
                   + (float(final.get("model_err", 0) or 0) * m) ** 2)
            bad = ((exp.mask.array & run.bad_bitmask) != 0) | ~np.isfinite(data - model) | ~(var > 0)
            z = (data - model) / np.sqrt(np.where(var > 0, var, np.nan))
            st = cn.robust_stats(z, ~bad)
            with contextlib.redirect_stdout(io.StringIO()):
                stk = show_bright_star_stack(run, min_snr=min_snr, plot=False)
            prof = []
            if stk.get("psf_only") is not None:
                h = stk["psf_only"].shape[0] // 2
                u = np.arange(-h, h + 1)
                R = np.hypot(*np.meshgrid(u, u))
                for a, b in zip(edges[:-1], edges[1:]):
                    sel = (R >= a) & (R < b) & np.isfinite(stk["psf_only"])
                    prof.append(1e5 * np.nanmean(stk["psf_only"][sel]) if sel.any() else np.nan)
            print(f"  {vname:22s} {st['sigma_iqr']:9.3f} {len(res['peaks']):6d}   "
                  + "  ".join(f"{p:+6.2f}" for p in prof) + f"   ({stk.get('n', 0)} stars)")


def core_excess_diagnostic(run, min_snr=100.0, plot=True):
    """Why are bright stars' cores sharper than the model?  For each unsaturated, reasonably
    isolated star with peak S/N >= min_snr: the CENTRAL EXCESS = (data - model) summed over the
    central 3x3 px / model summed there (final model and PSF-only model), against
        - flux                      -> grows with brightness: sensor effects (non-linearity,
                                       brighter-fatter correction) or flux-dependent PSF
        - peak pixel / saturation   -> grows near saturation: non-linearity
        - sub-pixel phase           -> larger for stars centred near a pixel CORNER: the PSF
                                       stamp is blurred when shifted to sub-pixel positions
        - nothing                   -> the same fraction for all: the PSF core itself is too
                                       broad (e.g. blurred when the PSF was built)
    Returns dict of per-star arrays + correlations."""
    from scipy.stats import spearmanr
    final = run.outer_history[-1]
    exp, res = final["exposure"], final["result"]
    S, fit = res["S"], res["snaps"][-1]["fit"]
    bb = exp.getBBox()
    c = bb.getCenter()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    data = exp.image.array
    psf = exp.getPsf()
    kern = psf.computeKernelImage(c).array
    fwhm = L._fwhm(psf, c)
    noise = float(np.sqrt(np.nanmedian(exp.variance.array)))
    amps, xy = fit["amps"], S["xy"]
    snr = amps * kern.max() / noise
    P = dict(_params(run.config))
    P.update(psf_stack_min_snr=min_snr, psf_stack_n=1000)
    idx = select_psf_stars(exp, res, fwhm, kern.max(), noise, P, verbose=False)
    lvl, how = saturation_level(exp, P)
    rows = []
    for j in idx:
        x, y = xy[j]
        ix, iy = int(round(x)) - X0, int(round(y)) - Y0
        sl = (slice(iy - 1, iy + 2), slice(ix - 1, ix + 2))
        if iy < 1 or ix < 1 or iy + 2 > data.shape[0] or ix + 2 > data.shape[1]:
            continue
        if not fit["good"][sl].all():
            continue
        d = data[sl]
        mf, mp = res["model_full"][sl], fit["model_full"][sl]
        ef = float(np.sum(d - mf) / np.sum(mf))
        ep = float(np.sum(d - mp) / np.sum(mp))
        cen = float((d - mf)[1, 1] / mf[1, 1])
        phase = float(np.hypot(x - round(x), y - round(y)))       # 0 = pixel centre, 0.71 = corner
        err = float(np.sqrt(np.sum(exp.variance.array[sl])) / np.sum(mf))
        rows.append((j, amps[j], snr[j], d.max(), ef, ep, cen, phase, err))
    if len(rows) < 6:
        print(f"only {len(rows)} usable star(s) -- lower min_snr")
        return None
    T = np.array(rows)
    out = dict(idx=T[:, 0].astype(int), flux=T[:, 1], snr=T[:, 2], peak=T[:, 3], excess=T[:, 4],
               excess_psf_only=T[:, 5], centre_px=T[:, 6], phase=T[:, 7], err=T[:, 8])
    ex = out["excess"]
    se = 1.2533 * 1.4826 * np.median(np.abs(ex - np.median(ex))) / np.sqrt(len(ex))   # error of the median
    print(f"{len(rows)} isolated unsaturated stars with peak S/N >= {min_snr:g}")
    print(f"  central 3x3 excess (data - model)/model:  median {100 * np.median(out['excess']):+.2f}% "
          f"(final model),  {100 * np.median(out['excess_psf_only']):+.2f}% (PSF-only)")
    corr = {}
    tests = [("flux", np.log10(out["flux"])), ("sub-pixel phase", out["phase"])]
    if lvl:
        out["peak_frac_sat"] = out["peak"] / lvl
        tests.append(("peak / saturation", out["peak_frac_sat"]))
        print(f"  saturation level {lvl:.4g} ({how}); brightest peak = "
              f"{100 * out['peak_frac_sat'].max():.0f}% of it")
    for lab, v in tests:
        rho, pv = spearmanr(v, out["excess"])
        corr[lab] = (float(rho), float(pv))
        print(f"  excess vs {lab:18s} rho = {rho:+.2f}  (p = {pv:.2g})"
              + ("   <-- significant" if pv < 0.01 and abs(rho) > 0.3 else ""))
    o = np.argsort(out["flux"])
    k = max(len(o) // 3, 2)
    fa, br = o[:k], o[-k:]
    print(f"  faintest third: excess {100 * np.median(out['excess'][fa]):+.2f}%   "
          f"brightest third: {100 * np.median(out['excess'][br]):+.2f}%  "
          f"(noise per star ~{100 * np.median(out['err']):.2f}%; error of the median "
          f"{100 * se:.2f}%)")
    print("  VERDICT:")
    v = []
    # a trend only counts if it is large compared with the excess itself: a uniformly too-broad
    # PSF also shows a mild flux trend, because the fit weights bright cores less (Poisson noise)
    e_fa, e_br = np.median(out["excess"][fa]), np.median(out["excess"][br])
    big_trend = abs(e_br - e_fa) > 0.3 * max(abs(np.median(ex)), 3 * se)
    if corr["flux"][1] < 0.01 and corr["flux"][0] > 0.3 and big_trend:
        v.append("excess GROWS WITH FLUX -> a brightness-dependent effect (sensor non-linearity / "
                 "brighter-fatter correction in the pipeline); one PSF cannot fit it")
    if ("peak / saturation" in corr and corr["peak / saturation"][1] < 0.01
            and corr["peak / saturation"][0] > 0.3 and big_trend):
        v.append("excess grows toward SATURATION -> non-linearity of the brightest pixels")
    if corr["sub-pixel phase"][1] < 0.01 and corr["sub-pixel phase"][0] > 0.3:
        v.append("excess larger for stars near a pixel CORNER -> the PSF stamp is blurred when "
                 "shifted to sub-pixel positions (rendering)")
    if not v:
        if np.median(ex) > 3 * se:
            v.append("same fractional excess for all stars -> the PSF CORE is slightly too broad "
                     "(blurred when it was built); a core-sharpening of the PSF should fix it")
        else:
            v.append("no significant central excess")
    for t in v:
        print("   - " + t)
    out.update(correlations=corr, verdict=v)
    if plot:
        import matplotlib.pyplot as plt
        n = len(tests)
        fig, ax = plt.subplots(1, n, figsize=(5.5 * n, 4.3))
        for a, (lab, vv) in zip(ax, tests):
            a.scatter(vv, 100 * out["excess"], s=20, c=np.log10(out["flux"]), cmap="viridis")
            a.axhline(0, color="k", lw=0.7, ls=":")
            a.set_xlabel(lab + (" (log10)" if lab == "flux" else ""))
            a.set_ylabel("central 3x3 excess [% of model]")
            a.set_title(f"rho = {corr[lab][0]:+.2f}", fontsize=10)
        plt.suptitle("Central excess of bright stars (colour = log flux)", fontsize=11)
        plt.tight_layout()
        plt.show()
    return out


def plot_residual_gray(run, **kw):
    """The two residual panels of plot_final_comparison, larger and in GREYSCALE
    (black = negative, white = positive; 'not modelable' areas in blue)."""
    kw.setdefault("panels", "residuals")
    kw.setdefault("grayscale", True)
    return plot_final_comparison(run, **kw)


def _print_run_version(h):
    """Which code version produced this run, and was the brightness term used in its fit?"""
    v = h.get("version")
    if v is None:
        print("  (run made by an OLDER cutout_gaia without the brightness term -- restart the kernel "
              "and re-run)")
        return
    bf = h.get("bf")
    if bf is None:
        print(f"  (run made by cutout_gaia {v}; brightness term NOT active in the final fit)")
    else:
        print(f"  (run made by cutout_gaia {v}; brightness term active in the final fit, "
              f"L0 = 10^{bf['L0']:.2f})")


def star_residual_table(run, min_snr=30.0, r_fwhm=2.5, plot=True, n_worst=10):
    """Which stars carry the leftover residual, and what do they have in common?

    For every unsaturated point source with peak S/N >= min_snr: chi2/px of (data - final
    model)/sigma within r_fwhm*FWHM (1.0 = pure noise), the same for the PSF-only joint model,
    and the PSF-only CORE residual (flux within 1 FWHM / star flux; < 0 = star broader than the
    PSF).  Then rank correlations of the excess with
        flux            -> shape depends on brightness (brighter-fatter)
        distance / x / y-> PSF varies across the cutout
        crowding        -> blends / neighbours
    and how concentrated the excess is (a few stars = individual problems).
    Returns a dict of per-star arrays + correlations + verdict."""
    from scipy.stats import spearmanr
    final = run.outer_history[-1]
    exp, res = final["exposure"], final["result"]
    S, fit = res["S"], res["snaps"][-1]["fit"]
    bb = exp.getBBox()
    c = bb.getCenter()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    data = exp.image.array
    var = exp.variance.array + float(final.get("poisson_k", 0.0) or 0.0) * np.clip(res["model_full"], 0, None)
    good = fit["good"] & np.isfinite(var) & (var > 0)
    psf = exp.getPsf()
    kern = psf.computeKernelImage(c).array
    fwhm = L._fwhm(psf, c)
    noise = float(np.sqrt(np.nanmedian(exp.variance.array)))
    amps, xy = fit["amps"], S["xy"]
    snr = amps * kern.max() / noise
    cand = np.flatnonzero((snr >= min_snr) & (S["ext"] == 0) & (amps > 0))
    sat = _satlike(exp)
    if sat.any():
        from scipy.ndimage import binary_dilation
        near = binary_dilation(sat, iterations=int(np.ceil(3 * fwhm)))
        ix = np.clip(np.round(xy[cand, 0]).astype(int) - X0, 0, sat.shape[1] - 1)
        iy = np.clip(np.round(xy[cand, 1]).astype(int) - Y0, 0, sat.shape[0] - 1)
        cand = cand[~near[iy, ix]]
    tree = cKDTree(xy)
    r = r_fwhm * fwhm
    rows = []
    for j in cand:
        x, y = xy[j]
        xa, xb = int(np.floor(x - r)) - X0, int(np.ceil(x + r)) - X0 + 1
        ya, yb = int(np.floor(y - r)) - Y0, int(np.ceil(y + r)) - Y0 + 1
        if xa < 0 or ya < 0 or xb > data.shape[1] or yb > data.shape[0]:
            continue
        yy, xx = np.mgrid[ya + Y0:yb + Y0, xa + X0:xb + X0]
        R = np.hypot(xx - x, yy - y)
        g = good[ya:yb, xa:xb] & (R < r)
        if g.sum() < 0.8 * np.pi * r * r:              # too much of it masked
            continue
        d, v = data[ya:yb, xa:xb], var[ya:yb, xa:xb]
        rf = d - res["model_full"][ya:yb, xa:xb]
        rp = d - fit["model_full"][ya:yb, xa:xb]
        chi_f = float(np.mean(rf[g] ** 2 / v[g]))
        chi_p = float(np.mean(rp[g] ** 2 / v[g]))
        core = float(np.sum(rp[g & (R < fwhm)]) / amps[j])
        # flux-independent shape mismatch: rms residual flux beyond the noise, / star flux
        mis = float(np.sqrt(max(np.sum(rf[g] ** 2) - np.sum(v[g]), 0.0)) / amps[j])
        nb = [k for k in tree.query_ball_point(xy[j], 2.5 * fwhm) if k != j and amps[k] > 0.1 * amps[j]]
        rows.append((j, amps[j], snr[j], x, y, np.hypot(x - c.getX(), y - c.getY()), len(nb),
                     chi_f, chi_p, core, int(g.sum()), mis))
    if len(rows) < 8:
        print(f"only {len(rows)} usable star(s) with S/N >= {min_snr:g} -- lower min_snr")
        return None
    T = np.array(rows)
    names = ["idx", "flux", "snr", "x", "y", "dist", "n_nb", "chi2_final", "chi2_psf", "core_psf", "npix",
             "mismatch"]
    out = {n: T[:, i] for i, n in enumerate(names)}
    out["idx"] = out["idx"].astype(int)
    ex = out["chi2_final"] - 1.0
    mis = out["mismatch"]
    lf = np.log10(out["flux"])
    corr = {}
    for lab, v in (("flux", lf), ("distance from centre", out["dist"]), ("x", out["x"]),
                   ("y", out["y"]), ("crowding (n neighbours)", out["n_nb"])):
        rho, pv = spearmanr(v, mis) if np.ptp(v) > 0 else (0.0, 1.0)
        corr[lab] = (float(rho), float(pv))
    rho_core, p_core = spearmanr(lf, out["core_psf"])
    o = np.argsort(-ex)
    mmed = np.median(mis[mis > 0]) if np.any(mis > 0) else 0.0
    outl = (mis > 3 * mmed) & (out["chi2_final"] > 1.5)
    nb_ = len(lf) // 3
    fo = np.argsort(lf)
    _print_run_version(final)
    print(f"{len(rows)} unsaturated stars with peak S/N >= {min_snr:g}; chi2/px within "
          f"{r_fwhm:g} FWHM (1.0 = noise only):")
    print(f"  {'':16s} {'median flux':>11s} {'chi2 final':>11s} {'chi2 PSF-only':>14s} "
          f"{'core (PSF-only)':>16s} {'mismatch':>9s}")
    for lab, sel in (("faintest third", fo[:nb_]), ("middle third", fo[nb_:2 * nb_]),
                     ("brightest third", fo[2 * nb_:])):
        print(f"  {lab:16s} {np.median(out['flux'][sel]):11.3g} {np.median(out['chi2_final'][sel]):11.2f} "
              f"{np.median(out['chi2_psf'][sel]):14.2f} {100 * np.median(out['core_psf'][sel]):+15.2f}% "
              f"{100 * np.median(mis[sel]):8.2f}%")
    print("  (mismatch = rms residual flux beyond the noise / star flux: a flux-independent shape error;")
    print("   chi2 alone always grows with flux for the same fractional error)")
    print("  rank correlation of the mismatch with:")
    for lab, (rho, pv) in corr.items():
        flag = "  <-- significant" if pv < 0.01 and abs(rho) > 0.3 else ""
        print(f"    {lab:24s} rho = {rho:+.2f}  (p = {pv:.1g}){flag}")
    print(f"    core residual vs flux    rho = {rho_core:+.2f}  (p = {p_core:.1g})  "
          f"[negative = brighter stars broader than the PSF]")
    print(f"  {int(outl.sum())} star(s) with mismatch > 3x the median ({100 * mmed:.2f}%) and chi2 > 1.5")
    verdict = []
    bright_chi = np.median(out["chi2_final"][fo[2 * nb_:]])
    if bright_chi < 1.1:
        verdict.append("brightest stars are fitted to within ~noise (chi2 < 1.1): the splotches in the "
                       "data-model panel are mostly their PHOTON NOISE (sigma grows as sqrt(flux)); "
                       "judge the residual/sigma panel")
    if p_core < 0.01 and abs(rho_core) > 0.3:
        verdict.append(f"PSF-only core residual {'falls' if rho_core < 0 else 'rises'} with flux -> the PSF "
                       f"core width depends on brightness")
    if any(corr[k][1] < 0.01 and abs(corr[k][0]) > 0.3 for k in ("distance from centre", "x", "y")):
        verdict.append("shape mismatch depends on POSITION -> the PSF varies across the cutout")
    if corr["crowding (n neighbours)"][1] < 0.01 and corr["crowding (n neighbours)"][0] > 0.3:
        verdict.append("mismatch grows with CROWDING -> blends / missing or misplaced neighbours")
    if outl.sum() >= max(2, 0.05 * len(mis)):
        verdict.append(f"{int(outl.sum())} individual OUTLIER star(s) (mismatch > 3x typical) -> "
                       f"blends / wrong positions: see the list")
    if not verdict:
        verdict.append("no single cause stands out: excess is spread over many stars without a trend")
    print("  VERDICT:")
    for v in verdict:
        print("   - " + v)
    print(f"  worst {min(n_worst, len(o))} stars:")
    print(f"    {'x':>8s} {'y':>8s} {'flux':>9s} {'S/N':>6s} {'chi2':>6s} {'mismatch':>9s} {'nbrs':>5s} {'Gaia':>5s}")
    for i in o[:n_worst]:
        j = out["idx"][i]
        print(f"    {out['x'][i]:8.1f} {out['y'][i]:8.1f} {out['flux'][i]:9.3g} {out['snr'][i]:6.0f} "
              f"{out['chi2_final'][i]:6.2f} {100 * mis[i]:8.2f}% {int(out['n_nb'][i]):5d} {'yes' if S['seed'][j] else 'no':>5s}")
    out.update(correlations=corr, core_vs_flux=(float(rho_core), float(p_core)),
               outliers=out["idx"][outl], verdict=verdict)
    if plot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(18, 5.2))
        sc = ax[0].scatter(out["flux"], out["chi2_final"], c=out["n_nb"], cmap="viridis", s=18,
                           vmin=0, vmax=max(3, out["n_nb"].max()))
        ax[0].axhline(1, color="k", lw=0.8, ls=":")
        ax[0].set_xscale("log")
        ax[0].set_yscale("log")
        ax[0].set_xlabel("star flux")
        ax[0].set_ylabel(f"chi2/px within {r_fwhm:g} FWHM (final model)")
        ax[0].set_title("excess vs brightness (colour = neighbours)")
        plt.colorbar(sc, ax=ax[0], fraction=0.046, label="neighbours > 10% flux")
        ax[1].scatter(out["flux"], 100 * out["core_psf"], s=18, color="C1")
        ax[1].axhline(0, color="k", lw=0.8, ls=":")
        ax[1].set_xscale("log")
        ax[1].set_xlabel("star flux")
        ax[1].set_ylabel("PSF-only residual within 1 FWHM [% of flux]")
        ax[1].set_title("core residual vs brightness (brighter-fatter test)")
        z = (data - res["model_full"]) / np.sqrt(var)
        z = np.where(good, z, np.nan)
        ax[2].imshow(z, origin="lower", cmap="gray_r", vmin=-4, vmax=4,
                     extent=(X0 - .5, X0 + data.shape[1] - .5, Y0 - .5, Y0 + data.shape[0] - .5))
        sc = ax[2].scatter(out["x"], out["y"], c=out["chi2_final"], cmap="magma_r", s=40,
                           vmin=1, vmax=max(2, np.percentile(out["chi2_final"], 95)),
                           edgecolors="k", linewidths=0.4)
        plt.colorbar(sc, ax=ax[2], fraction=0.046, label="chi2/px")
        ax[2].set_title("where the excess is (on residual/sigma)")
        ax[2].set_xlabel("x [px]")
        ax[2].set_ylabel("y [px]")
        plt.tight_layout()
        plt.show()
    return out


def plot_final_comparison(run, view="interior", bad_color=BAD_PIXEL_COLOR, circles="images",
                          show_zones=True, exclude_near_mask=True, excluded_color=EXCLUDED_COLOR,
                          mask_pad_fwhm=1.5, z_clip=5.0, resid_clip_pct=99.0, figsize=(26, 7),
                          sigma="model", nearsat_color=NEARSAT_COLOR, grayscale=False,
                          panels="all", verbose=True):
    """Original | final model | residual | residual/sigma, bad pixels drawn OPAQUE in one
    colour (default magenta) in all four panels and excluded from the statistics.
    circles: where to draw the source circles -- 'images' (data + model only, default),
    'all', or 'none'.  Saturated-star zones: dashed outline.  Extended: dashed circles.
    exclude_near_mask: grey out (and leave out of the statistics) pixels within
    mask_pad_fwhm*FWHM of bad pixels plus residual blobs touching them -- flux of sources
    whose centre is under the mask, which the model cannot fit.
    view: 'interior' | None (whole fitted array) | (xmin, xmax, ymin, ymax).
    sigma: 'model' (default) = sqrt(variance + (eps*model)^2 + k*model) with the measured PSF
    model error eps (option A); 'effective' = without eps; 'pixel' = variance plane only.
    Near-saturated cores excluded by option B are drawn in nearsat_color (orange); they are
    modelled from their wings, so only significant residual next to them is greyed out.
    grayscale: residual panels in greyscale (black = negative, white = positive); the
    'not modelable' areas are then drawn in blue instead of grey.
    panels: 'all' (4 panels) or 'residuals' (only the two residual panels, larger).
    Returns dict(residual=..., residual_over_sigma=...) of robust statistics."""
    from . import session as cn
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.colors import ListedColormap
    from astropy.visualization import simple_norm
    final = run.outer_history[-1]
    exp, res = final["exposure"], final["result"]
    bb = exp.getBBox()
    data, model, var = exp.image.array, res["model_full"], exp.variance.array
    raw = run.raw_crop[bb].image.array
    bad = ((exp.mask.array & run.bad_bitmask) != 0) | ~np.isfinite(data - model)
    ny, nx = data.shape
    x0, y0 = bb.getMinX(), bb.getMinY()
    ext = (x0 - .5, x0 + nx - .5, y0 - .5, y0 + ny - .5)
    xs, ys = x0 + np.arange(nx), y0 + np.arange(ny)
    if isinstance(view, str) and view == "interior":
        b = run.interior_bbox
        view = (b.getMinX() - .5, b.getMaxX() + .5, b.getMinY() - .5, b.getMaxY() + .5)
    elif isinstance(view, str):
        view = None
    region = np.ones(data.shape, bool)
    if view is not None:
        region = ((xs[None, :] >= view[0]) & (xs[None, :] <= view[1]) &
                  (ys[:, None] >= view[2]) & (ys[:, None] <= view[3]))
    resid = data - model
    k_used, eps_used = 0.0, 0.0
    if sigma in ("effective", "model"):   # include the photon noise the variance plane misses
        k_used = float(final.get("poisson_k", 0.0) or 0.0)
        try:
            k_used = max(k_used, float(run.config.get("poisson_k", 0.0) or 0.0))
        except Exception:
            pass
        var = var + k_used * np.clip(model, 0, None)
    if sigma == "model":                  # + the measured PSF model error (option A)
        eps_used = float(final.get("model_err", 0.0) or 0.0)
        var = var + (eps_used * np.clip(model, 0, None)) ** 2
    varok = np.isfinite(var) & (var > 0)
    z = resid / np.sqrt(np.where(varok, var, np.nan))
    fw = L._fwhm(exp.getPsf(), bb.getCenter())
    try:
        nearsat = (exp.mask.array & exp.mask.getPlaneBitMask("NEARSAT")) != 0
    except Exception:
        nearsat = np.zeros(bad.shape, bool)
    excl = (unmodeled_near_mask(resid, var, bad, fw, mask_pad_fwhm, modeled=nearsat)
            if exclude_near_mask else np.zeros(bad.shape, bool))
    use = ~bad & ~excl & region
    st_r = cn.robust_stats(resid, use)
    st_z = cn.robust_stats(z, use & varok)
    sel = resid[use]
    rs = np.nanpercentile(np.abs(sel), resid_clip_pct) if sel.size else 1.0
    norm = simple_norm(raw, vmin=np.nanmin(raw), vmax=np.nanpercentile(raw, 99),
                       stretch="log", log_a=1000)
    cm_res = plt.get_cmap("gray" if grayscale else "RdBu_r").copy()
    cm_res.set_bad(bad_color)
    if grayscale and excluded_color == EXCLUDED_COLOR:
        excluded_color = EXCLUDED_COLOR_GRAY            # grey would vanish in a grey map
    zone_col = "#00b4d8" if grayscale else "k"
    badimg = np.ma.masked_where(~bad, np.ones(bad.shape))
    S = res.get("S")
    ext_src = (S["ext"] > 0) if S is not None and "ext" in S else np.zeros(len(res["peaks"]), bool)
    peaks = np.asarray(res["peaks"], float).reshape(-1, 2)
    zmask = None
    if show_zones:
        try:
            P = _params(run.config)
            import contextlib
            import io
            with contextlib.redirect_stdout(io.StringIO()):
                zmask = (sat_zones(exp, L._fwhm(exp.getPsf(), bb.getCenter()), P)[1]
                         if P["sat_zone_fwhm"] else None)
        except Exception:
            zmask = None
    radius = 0.75 * L._fwhm(exp.getPsf(), bb.getCenter())

    def sub(st):
        return f"median {st['median']:+.2f},  robust σ {st['sigma_iqr']:.2f}  (N={st['n']})"

    plist = [("Original visit image (cutout)", raw, dict(cmap=L.CMAP, norm=norm)),
              (f"Final model ({len(peaks)} sources)", model, dict(cmap=L.CMAP, norm=norm)),
              ("Residual = data − model\n" + sub(st_r), np.ma.masked_where(bad | excl, resid),
               dict(cmap=cm_res, vmin=-rs, vmax=rs)),
              (("Residual / σ_eff  (σ² + (%.2f%%·model)²%s)" % (100 * eps_used, " + %.3g·model" % k_used if k_used else "")
                if sigma == "model" else
                "Residual / σ_eff (σ² + %.3g·model)" % k_used if sigma == "effective" else "Residual / σ")
               + "\n" + sub(st_z), np.ma.masked_where(bad | excl | ~varok, z),
               dict(cmap=cm_res, vmin=-z_clip, vmax=z_clip))]
    show = [0, 1, 2, 3] if panels == "all" else [2, 3]
    if panels != "all" and figsize == (26, 7):
        figsize = (16, 7.5)
    fig, axes = plt.subplots(1, len(show), figsize=figsize, sharex=True, sharey=True)
    for ax, i in zip(axes, show):
        title, img, kw = plist[i]
        im = ax.imshow(img, extent=ext, origin="lower", **kw)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        if bad.any():
            ax.imshow(badimg, extent=ext, origin="lower", cmap=ListedColormap([bad_color]),
                      interpolation="nearest", zorder=3)
        if nearsat.any():
            ax.imshow(np.ma.masked_where(~nearsat, np.ones(bad.shape)), extent=ext, origin="lower",
                      cmap=ListedColormap([nearsat_color]), interpolation="nearest", zorder=3.5)
        if excl.any():
            if i >= 2:
                ax.imshow(np.ma.masked_where(~excl, np.ones(excl.shape)), extent=ext,
                          origin="lower", cmap=ListedColormap([excluded_color]),
                          interpolation="nearest", zorder=3)
            else:
                ax.contour(xs, ys, excl.astype(float), levels=[0.5], colors=excluded_color,
                           linewidths=1.2, zorder=4)
        if zmask is not None and zmask.any():
            ax.contour(xs, ys, zmask.astype(float), levels=[0.5], linestyles="--",
                       colors="w" if i < 2 else zone_col, linewidths=1.0, zorder=4)
        if len(peaks) and (circles == "all" or (circles == "images" and i < 2)):
            col = "w" if i < 2 else zone_col
            for (px, py), e in zip(peaks, ext_src):
                ax.add_patch(mpatches.Circle((px, py), radius * (1.6 if e else 1.0), fill=False,
                                             ec=col, lw=0.6 if not e else 0.9,
                                             ls="--" if e else "-",
                                             alpha=0.6 if i < 2 else 0.35, zorder=5))
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("x [px]")
        ax.grid(False)
        if i >= 2:
            st = st_r if i == 2 else st_z
            ax.text(0.02, 0.98, f"median = {st['median']:+.3f}\nσ_IQR = {st['sigma_iqr']:.3f}",
                    transform=ax.transAxes, va="top", ha="left", fontsize=10,
                    bbox=dict(fc="white", ec="k", alpha=0.85, pad=3), zorder=10)
    axes[0].set_ylabel("y [px]")
    axes[0].set_xlim(*(view[:2] if view else ext[:2]))
    axes[0].set_ylim(*(view[2:] if view else ext[2:]))
    handles = [mpatches.Patch(fc=bad_color, ec="k", label="bad pixel (excluded)")]
    if nearsat.any():
        handles.append(mpatches.Patch(fc=nearsat_color, ec="k",
                                      label="near-saturated core, fitted from wings (excluded)"))
    if excl.any():
        handles.append(mpatches.Patch(fc=excluded_color, ec="k",
                                      label="next to bad pixels, not modelable (excluded)"))
    if zmask is not None and zmask.any():
        handles.append(mpatches.Patch(fc="none", ec="k", ls="--",
                                      label="saturated-star zone (non-Gaia S/N ≥ 10)"))
    if ext_src.any() and circles != "none":
        handles.append(mpatches.Patch(fc="none", ec="k", ls="--",
                                      label=f"extended source ({int(ext_src.sum())})"))
    axes[show.index(2)].legend(handles=handles, loc="lower left", fontsize=8, framealpha=0.9)
    fig.text(0.5, -0.01, "Statistics exclude bad pixels" +
             (" and the grey unmodelable areas next to them" if excl.any() else "") +
             "; median and robust σ = (Q75−Q25)/1.349"
             + (" over the plotted view" if view is not None else " over the full array"),
             ha="center", fontsize=9)
    plt.tight_layout()
    stats = dict(residual=st_r, residual_over_sigma=st_z, excluded_px=int((excl & region).sum()))
    if verbose:
        cn._print_stats({k: v for k, v in stats.items() if isinstance(v, dict)})
        if excl.any():
            print(f"excluded next to bad pixels: {int((excl & region).sum())} px "
                  f"({100 * (excl & region).sum() / max(region.sum(), 1):.1f}% of the view)")
    return stats
