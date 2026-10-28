# This file is part of meas_crowd.
#
# Developed for the LSST Data Management System.
# See the COPYRIGHT file at the top-level directory of this distribution.

"""Per-source cutouts: data | model | residual | residual/sigma, x and y cross-sections and the
radial profile, labelled with the source's PSF and aperture magnitudes.

`source_cutout_data` returns the numbers behind one panel row (unit tested); `plot_source_cutouts`
draws them for a selection of sources (worst fitted, brightest, random, or given positions) and
can write every page to a PDF.
"""

__all__ = ["source_cutout_data", "select_sources", "plot_source_cutouts"]

import numpy as np
from scipy.spatial import cKDTree


def _effective_variance(final, model, var, sigma):
    m = np.clip(model, 0, None)
    v = var.copy()
    if sigma in ("effective", "model"):
        v = v + float(final.get("poisson_k", 0.0) or 0.0) * m
    if sigma == "model":
        v = v + (float(final.get("model_err", 0.0) or 0.0) * m) ** 2
    return v


def _cut(D, M, V, B):
    """Mean over the rows of a (k, n) strip of good pixels: data, model, sigma of the mean."""
    ok = ~B
    nn = ok.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        d = np.where(nn > 0, np.sum(np.where(ok, D, 0), axis=0) / nn, np.nan)
        m = np.where(nn > 0, np.sum(np.where(ok, M, 0), axis=0) / nn, np.nan)
        s = np.where(nn > 0, np.sqrt(np.sum(np.where(ok, V, 0), axis=0)) / nn, np.nan)
    return d, m, s


def source_cutout_data(data, model, var, bad, x, y, x0=0, y0=0, half=15, width=1):
    """Everything one cutout row shows, for the source at (x, y) (parent coordinates; arrays start
    at (x0, y0)).

    Returns a dict with
        ``window``  (xmin, xmax, ymin, ymax) parent pixel range shown,
        ``data``, ``model``, ``resid``, ``z`` (residual / sigma), ``bad`` image arrays,
        ``xcut`` / ``ycut``: dict(offset, data, model, resid, sigma) through the source's row /
        column (``width`` rows or columns averaged),
        ``radial``: dict(r, data, model, resid, sigma) azimuthal means in 1-px annuli,
        ``chi2`` dict(x=, y=, r=): chi2 per point of each cut.
    """
    ny, nx = data.shape
    ix, iy = int(round(x)) - x0, int(round(y)) - y0
    xa, xb = max(ix - half, 0), min(ix + half + 1, nx)
    ya, yb = max(iy - half, 0), min(iy + half + 1, ny)
    sl = (slice(ya, yb), slice(xa, xb))
    sig = np.sqrt(np.where(var > 0, var, np.nan))
    out = dict(window=(x0 + xa, x0 + xb - 1, y0 + ya, y0 + yb - 1), data=data[sl], model=model[sl],
               resid=(data - model)[sl], z=((data - model) / sig)[sl], bad=bad[sl])
    w2 = max(int(width) // 2, 0)
    rows = slice(max(iy - w2, 0), min(iy + w2 + 1, ny))
    cols = slice(max(ix - w2, 0), min(ix + w2 + 1, nx))
    d, m, s = _cut(data[rows, xa:xb], model[rows, xa:xb], var[rows, xa:xb], bad[rows, xa:xb])
    out["xcut"] = dict(offset=x0 + np.arange(xa, xb) - x, data=d, model=m, resid=d - m, sigma=s)
    d, m, s = _cut(data[ya:yb, cols].T, model[ya:yb, cols].T, var[ya:yb, cols].T, bad[ya:yb, cols].T)
    out["ycut"] = dict(offset=y0 + np.arange(ya, yb) - y, data=d, model=m, resid=d - m, sigma=s)
    yy, xx = np.mgrid[ya:yb, xa:xb]
    r = np.hypot(xx + x0 - x, yy + y0 - y)
    ok = ~bad[sl]
    rc, rd, rm, rs = [], [], [], []
    for a in range(half):
        k = ok & (r >= a) & (r < a + 1)
        if k.sum() < 3:
            continue
        rc.append(r[k].mean())
        rd.append(data[sl][k].mean())
        rm.append(model[sl][k].mean())
        rs.append(np.sqrt(var[sl][k].sum()) / k.sum())
    rd, rm = np.array(rd), np.array(rm)
    out["radial"] = dict(r=np.array(rc), data=rd, model=rm, resid=rd - rm, sigma=np.array(rs))
    chi = {}
    for key, cut in (("x", out["xcut"]), ("y", out["ycut"]), ("r", out["radial"])):
        g = np.isfinite(cut["resid"]) & (cut["sigma"] > 0)
        chi[key] = float(np.mean((cut["resid"][g] / cut["sigma"][g]) ** 2)) if g.any() else np.nan
    out["chi2"] = chi
    return out


def select_sources(phot, how="worst", n=6, positions=None, xy=None, chi2=None, min_snr=20.0,
                   seed=0):
    """Indices of the sources to show.

    how : 'worst' (highest local chi2 among S/N >= min_snr), 'brightest', 'faintest' (lowest
    S/N above min_snr), 'random', or 'positions' (nearest source to each (x, y) in
    ``positions``)."""
    if how == "positions" or positions is not None:
        return [int(cKDTree(xy).query(p)[1]) for p in positions]
    snr = np.asarray(phot["psf_snr"], float)
    ok = np.flatnonzero(np.isfinite(snr) & (snr >= min_snr))
    if not len(ok):
        return []
    if how == "worst" and chi2 is not None:
        return list(ok[np.argsort(-np.asarray(chi2)[ok])][:n])
    if how in ("brightest", "worst"):
        return list(ok[np.argsort(-snr[ok])][:n])
    if how == "faintest":
        return list(ok[np.argsort(snr[ok])][:n])
    if how == "random":
        rng = np.random.default_rng(seed)
        return list(rng.choice(ok, size=min(n, len(ok)), replace=False))
    raise ValueError(f"unknown selection {how!r}")


def _local_chi2(data, model, var, good, xy, x0, y0, r):
    """chi2/px within r px of each source."""
    out = np.full(len(xy), np.nan)
    ny, nx = data.shape
    for j, (x, y) in enumerate(xy):
        ix, iy = int(round(x)) - x0, int(round(y)) - y0
        xa, xb = max(ix - int(r) - 1, 0), min(ix + int(r) + 2, nx)
        ya, yb = max(iy - int(r) - 1, 0), min(iy + int(r) + 2, ny)
        if xb <= xa or yb <= ya:
            continue
        yy, xx = np.mgrid[ya:yb, xa:xb]
        m = good[ya:yb, xa:xb] & (np.hypot(xx + x0 - x, yy + y0 - y) < r) & (var[ya:yb, xa:xb] > 0)
        if m.sum() >= 5:
            out[j] = np.mean((data - model)[ya:yb, xa:xb][m] ** 2 / var[ya:yb, xa:xb][m])
    return out


def plot_source_cutouts(run, how="worst", n=6, positions=None, half=15, width=1, sigma="model",
                        phot=None, min_snr=20.0, log=True, per_page=4, pdf_path=None, show=True,
                        grayscale=False):
    """Cutouts and cross-sections of individual sources, labelled with their magnitudes.

    One row per source: data | model | residual | residual/sigma | x-cut | y-cut | radial
    profile.  The cuts show data (points) and model (line) on top, and the residual with its
    +-1 sigma band below; other sources in the window are marked (o Gaia, x found in the
    residual).

    Parameters
    ----------
    run : `CutoutRun`
        A finished run.
    how : `str`
        'worst' (default), 'brightest', 'faintest', 'random' -- or pass ``positions``.
    positions : list of (x, y), optional
        Show the sources nearest to these detector positions.
    phot : `astropy.table.Table`, optional
        Output of `measure_photometry` (computed if not given).
    pdf_path : `str`, optional
        Also write all pages to this PDF.

    Returns the photometry table rows of the sources shown.
    """
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import ListedColormap, SymLogNorm
    from .photometry import measure_photometry
    from .plotting import BAD_PIXEL_COLOR
    final = run.outer_history[-1]
    exp, res = final["exposure"], final["result"]
    S = res["S"]
    bb = exp.getBBox()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    data, model = exp.image.array, res["model_full"]
    var = _effective_variance(final, model, exp.variance.array, sigma)
    bad = ((exp.mask.array & run.bad_bitmask) != 0) | ~np.isfinite(data - model) | ~(var > 0)
    if phot is None:
        phot = measure_photometry(run, verbose=False)
    xy = S["xy"]
    fwhm = phot.meta.get("fwhm_px", 5.0)
    chi2 = _local_chi2(data, model, var, ~bad, xy, X0, Y0, 2.5 * fwhm)
    sel = select_sources(phot, how, n, positions, xy, chi2, min_snr)
    if not sel:
        print("no sources selected")
        return phot[:0]
    noise = float(np.sqrt(np.nanmedian(exp.variance.array)))
    tree = cKDTree(xy)
    cm_res = plt.get_cmap("gray" if grayscale else "RdBu_r").copy()
    pages = [sel[i:i + per_page] for i in range(0, len(sel), per_page)]
    pdf = PdfPages(pdf_path) if pdf_path else None
    try:
        for page in pages:
            fig = plt.figure(figsize=(26, 4.4 * len(page)))
            outer = fig.add_gridspec(len(page), 1, hspace=0.45)
            for row, j in enumerate(page):
                x, y = xy[j]
                cd = source_cutout_data(data, model, var, bad, x, y, X0, Y0, half, width)
                gs = outer[row].subgridspec(2, 7, height_ratios=[2, 1], hspace=0.05, wspace=0.35)
                w = cd["window"]
                ext = (w[0] - .5, w[1] + .5, w[2] - .5, w[3] + .5)
                lo, hi = np.nanpercentile(cd["data"], 1), np.nanmax(cd["data"])
                nrm = SymLogNorm(linthresh=max(3 * noise, 1e-6), vmin=lo, vmax=max(hi, lo + 1e-6))
                good_r = cd["resid"][~cd["bad"]]
                rs = np.nanpercentile(np.abs(good_r), 99) if good_r.size else 1.0
                imgs = [(cd["data"], dict(cmap="viridis", norm=nrm), "data"),
                        (cd["model"], dict(cmap="viridis", norm=nrm), "model"),
                        (np.ma.masked_where(cd["bad"], cd["resid"]),
                         dict(cmap=cm_res, vmin=-rs, vmax=rs), "data − model"),
                        (np.ma.masked_where(cd["bad"], cd["z"]),
                         dict(cmap=cm_res, vmin=-5, vmax=5), "residual / σ")]
                near = tree.query_ball_point((x, y), half * 1.42)
                for c, (img, kw, ttl) in enumerate(imgs):
                    ax = fig.add_subplot(gs[:, c])
                    im = ax.imshow(img, origin="lower", extent=ext, interpolation="nearest", **kw)
                    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
                    if cd["bad"].any():
                        ax.imshow(np.ma.masked_where(~cd["bad"], np.ones(cd["bad"].shape)),
                                  origin="lower", extent=ext, interpolation="nearest",
                                  cmap=ListedColormap([BAD_PIXEL_COLOR]))
                    mc = "w" if c < 2 else ("#00b4d8" if grayscale else "k")
                    for k in near:
                        px, py = xy[k]
                        if S["seed"][k]:
                            ax.plot(px, py, "o", mfc="none", mec=mc, ms=10 if k == j else 6, mew=1.2)
                        else:
                            ax.plot(px, py, "x", color="#ff3b3b", ms=7, mew=1.5)
                    ax.set_xlim(ext[:2])
                    ax.set_ylim(ext[2:])
                    ax.set_title(ttl, fontsize=9)
                    if c == 0:
                        p = phot[j]
                        lab = (f"#{j} ({x:.1f}, {y:.1f})  {'Gaia' if S['seed'][j] else 'non-Gaia'}"
                               f"{', extended' if S['ext'][j] > 0 else ''}\n"
                               f"PSF  {p['psf_mag']:.3f} ± {p['psf_mag_err']:.3f} mag  "
                               f"(S/N {p['psf_snr']:.0f})\n"
                               f"ap   {p['ap_mag']:.3f} ± {p['ap_mag_err']:.3f} mag  "
                               f"(ap−PSF {p['ap_minus_psf']:+.3f})\n"
                               + (f"Gaia G {p['gaia_G']:.2f}   " if np.isfinite(p["gaia_G"]) else "")
                               + f"χ²/px {chi2[j]:.2f}")
                        ax.set_ylabel(lab, fontsize=8, rotation=0, ha="right", va="center",
                                      labelpad=10)
                cuts = (("xcut", "cut along x (row)", "Δx [px]"),
                        ("ycut", "cut along y (column)", "Δy [px]"),
                        ("radial", "radial profile", "r [px]"))
                for c, (key, ttl, xl) in enumerate(cuts):
                    cut = cd[key]
                    xs_ = cut["r"] if key == "radial" else cut["offset"]
                    a1 = fig.add_subplot(gs[0, 4 + c])
                    a2 = fig.add_subplot(gs[1, 4 + c], sharex=a1)
                    a1.plot(xs_, cut["data"], "o", ms=3.5, color="k", label="data")
                    a1.plot(xs_, cut["model"], "-", lw=1.5, color="C1", label="model")
                    if log:
                        a1.set_yscale("symlog", linthresh=max(3 * noise, 1e-6))
                    a2.fill_between(xs_, -cut["sigma"], cut["sigma"], color="0.85", step="mid")
                    a2.axhline(0, color="k", lw=0.7)
                    a2.plot(xs_, cut["resid"], "o-", ms=3, lw=0.8, color="C3")
                    a2.text(0.01, 0.92, f"χ²/pt {cd['chi2'][key[0]]:.2f}", transform=a2.transAxes,
                            va="top", fontsize=8)
                    if key != "radial":
                        for k in near:
                            if k == j:
                                continue
                            off = (xy[k, 0] - x) if key == "xcut" else (xy[k, 1] - y)
                            perp = (xy[k, 1] - y) if key == "xcut" else (xy[k, 0] - x)
                            if abs(perp) <= max(width // 2, 0) + 1.0 and abs(off) <= half:
                                for a in (a1, a2):
                                    a.axvline(off, color="C0" if S["seed"][k] else "C3",
                                              ls=":" if S["seed"][k] else "--", lw=0.9)
                    a1.set_title(ttl, fontsize=9)
                    a1.tick_params(labelbottom=False)
                    a2.set_xlabel(xl, fontsize=8)
                    if row == 0 and c == 0:
                        a1.legend(fontsize=7, loc="upper right")
            fig.suptitle(f"{run.name if hasattr(run, 'name') else ''}: source cutouts "
                         f"(zero point {phot.meta.get('zero_point', np.nan):.3f}, "
                         f"{phot.meta.get('zero_point_source', '')}); grey band = ±1σ; "
                         f"o Gaia, x found in residual", fontsize=10)
            if pdf is not None:
                pdf.savefig(fig, bbox_inches="tight")
            if show:
                plt.show()
            else:
                plt.close(fig)
    finally:
        if pdf is not None:
            pdf.close()
            print(f"wrote {len(pages)} page(s) to {pdf_path}")
    return phot[sel]
