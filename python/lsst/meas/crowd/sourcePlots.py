# This file is part of meas_crowd.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Per-source cutouts, cross-sections and magnitudes.

For each source: data | model | residual | residual/sigma, x and y cross-sections and the
radial profile, labelled with the source's PSF and aperture magnitudes.

`sourceCutoutData` returns the numbers behind one panel row (unit tested); `plotSourceCutouts`
draws them for a selection of sources (worst fitted, brightest, random, or given positions) and
can write every page to a PDF.
"""

__all__ = ["sourceCutoutData", "selectSources", "plotSourceCutouts"]

import numpy as np
from scipy.spatial import cKDTree


def _effectiveVariance(final, model, var, sigma):
    """Return the variance for residual/sigma ('pixel', 'effective' or 'model')."""
    m = np.clip(model, 0, None)
    v = var.copy()
    if sigma in ("effective", "model"):
        v = v + float(final.get("poissonK", 0.0) or 0.0) * m
    if sigma == "model":
        v = v + (float(final.get("modelErr", 0.0) or 0.0) * m) ** 2
    return v


def _cut(D, M, V, B):
    """Return the mean data, model and sigma of the mean of a strip of good pixels."""
    ok = ~B
    nn = ok.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        d = np.where(nn > 0, np.sum(np.where(ok, D, 0), axis=0) / nn, np.nan)
        m = np.where(nn > 0, np.sum(np.where(ok, M, 0), axis=0) / nn, np.nan)
        s = np.where(nn > 0, np.sqrt(np.sum(np.where(ok, V, 0), axis=0)) / nn, np.nan)
    return d, m, s


def sourceCutoutData(data, model, var, bad, x, y, x0=0, y0=0, half=15, width=1):
    """Compute the arrays one cutout row shows for a source.

    Parameters
    ----------
    data : `numpy.ndarray`
        Image.
    model : `numpy.ndarray`
        Model image.
    var : `numpy.ndarray`
        Variance image.
    bad : `numpy.ndarray` [`bool`]
        Bad pixels.
    x, y : `float`
        Source position in parent pixel coordinates.
    x0, y0 : `int`, optional
        Parent pixel coordinates of ``data[0, 0]``.
    half : `int`, optional
        Half-size of the window [px].
    width : `int`, optional
        Rows / columns averaged in the x / y cuts.

    Returns
    -------
    cutout : `dict`
        ``window`` (xmin, xmax, ymin, ymax), the ``data``, ``model``, ``resid``, ``z``
        (residual / sigma) and ``bad`` arrays, ``xcut`` / ``ycut`` / ``radial`` (each a dict
        with the coordinate, ``data``, ``model``, ``resid`` and ``sigma``) and ``chi2`` (per
        point of each cut).
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


def selectSources(phot, how="worst", n=6, positions=None, xy=None, chi2=None, minSnr=20.0,
                  seed=0):
    """Choose the sources to show.

    Parameters
    ----------
    phot : `astropy.table.Table` or `dict`
        Photometry with ``psfSnr``.
    how : `str`, optional
        'worst' (highest ``chi2``), 'brightest', 'faintest' or 'random', among sources with S/N
        >= ``minSnr``.
    n : `int`, optional
        Number of sources.
    positions : `list` [`tuple` [`float`, `float`]], optional
        Show the nearest source to each position instead.
    xy : `numpy.ndarray`, optional
        Source positions (needed with ``positions``).
    chi2 : `numpy.ndarray`, optional
        Local chi2 per source (needed for 'worst').
    minSnr : `float`, optional
        Minimum S/N.
    seed : `int`, optional
        Random seed for 'random'.

    Returns
    -------
    indices : `list` [`int`]
        Selected source indices.

    Raises
    ------
    ValueError
        Raised if ``how`` is not a known selection.
    """
    if how == "positions" or positions is not None:
        return [int(cKDTree(xy).query(p)[1]) for p in positions]
    snr = np.asarray(phot["psfSnr"], float)
    ok = np.flatnonzero(np.isfinite(snr) & (snr >= minSnr))
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


def _localChi2(data, model, var, good, xy, x0, y0, r):
    """Return the chi2 per pixel within ``r`` px of each source."""
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


def plotSourceCutouts(run, how="worst", n=6, positions=None, half=15, width=1, sigma="model",
                      phot=None, minSnr=20.0, log=True, perPage=4, pdfPath=None, show=True,
                      grayscale=False):
    """Plot cutouts and cross-sections of individual sources, labelled with their magnitudes.

    One row per source: data | model | residual | residual/sigma | x-cut | y-cut | radial
    profile.  The cuts show data (points) and model (line), and the residual with its +-1 sigma
    band; other sources in the window are marked (o Gaia, x found in the residual).

    Parameters
    ----------
    run : `lsst.meas.crowd.session.CutoutRun`
        A finished run.
    how : `str`, optional
        'worst', 'brightest', 'faintest' or 'random' (see `selectSources`).
    n : `int`, optional
        Number of sources.
    positions : `list` [`tuple` [`float`, `float`]], optional
        Show the sources nearest these detector positions.
    half : `int`, optional
        Half-size of each window [px].
    width : `int`, optional
        Rows / columns averaged in the cuts.
    sigma : `str`, optional
        'model' (variance plus measured model error), 'effective' or 'pixel'.
    phot : `astropy.table.Table`, optional
        Output of `lsst.meas.crowd.photometry.measurePhotometry` (computed if not given).
    minSnr : `float`, optional
        Minimum S/N of the sources selected.
    log : `bool`, optional
        Symmetric-log axis for data and model.
    perPage : `int`, optional
        Sources per figure.
    pdfPath : `str`, optional
        Also write all pages to this PDF.
    show : `bool`, optional
        Display the figures.
    grayscale : `bool`, optional
        Residual panels in greyscale.

    Returns
    -------
    rows : `astropy.table.Table`
        Photometry rows of the sources shown.
    """
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import ListedColormap, SymLogNorm
    from .photometry import measurePhotometry
    from .plotting import BAD_PIXEL_COLOR
    final = run.outerHistory[-1]
    exp, res = final["exposure"], final["result"]
    S = res["S"]
    bb = exp.getBBox()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    data, model = exp.image.array, res["modelFull"]
    var = _effectiveVariance(final, model, exp.variance.array, sigma)
    bad = ((exp.mask.array & run.badBitmask) != 0) | ~np.isfinite(data - model) | ~(var > 0)
    if phot is None:
        phot = measurePhotometry(run, verbose=False)
    xy = S["xy"]
    fwhm = phot.meta.get("fwhmPx", 5.0)
    chi2 = _localChi2(data, model, var, ~bad, xy, X0, Y0, 2.5 * fwhm)
    sel = selectSources(phot, how, n, positions, xy, chi2, minSnr)
    if not sel:
        print("no sources selected")
        return phot[:0]
    noise = float(np.sqrt(np.nanmedian(exp.variance.array)))
    tree = cKDTree(xy)
    cmRes = plt.get_cmap("gray" if grayscale else "RdBu_r").copy()
    pages = [sel[i:i + perPage] for i in range(0, len(sel), perPage)]
    pdf = PdfPages(pdfPath) if pdfPath else None
    try:
        for page in pages:
            fig = plt.figure(figsize=(26, 4.4 * len(page)))
            outer = fig.add_gridspec(len(page), 1, hspace=0.45)
            for row, j in enumerate(page):
                x, y = xy[j]
                cd = sourceCutoutData(data, model, var, bad, x, y, X0, Y0, half, width)
                gs = outer[row].subgridspec(2, 7, height_ratios=[2, 1], hspace=0.05, wspace=0.35)
                w = cd["window"]
                ext = (w[0] - .5, w[1] + .5, w[2] - .5, w[3] + .5)
                lo, hi = np.nanpercentile(cd["data"], 1), np.nanmax(cd["data"])
                nrm = SymLogNorm(linthresh=max(3 * noise, 1e-6), vmin=lo, vmax=max(hi, lo + 1e-6))
                goodR = cd["resid"][~cd["bad"]]
                rs = np.nanpercentile(np.abs(goodR), 99) if goodR.size else 1.0
                imgs = [(cd["data"], dict(cmap="viridis", norm=nrm), "data"),
                        (cd["model"], dict(cmap="viridis", norm=nrm), "model"),
                        (np.ma.masked_where(cd["bad"], cd["resid"]),
                         dict(cmap=cmRes, vmin=-rs, vmax=rs), "data − model"),
                        (np.ma.masked_where(cd["bad"], cd["z"]),
                         dict(cmap=cmRes, vmin=-5, vmax=5), "residual / σ")]
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
                        if S["fromGaia"][k]:
                            ax.plot(px, py, "o", mfc="none", mec=mc, ms=10 if k == j else 6, mew=1.2)
                        else:
                            ax.plot(px, py, "x", color="#ff3b3b", ms=7, mew=1.5)
                    ax.set_xlim(ext[:2])
                    ax.set_ylim(ext[2:])
                    ax.set_title(ttl, fontsize=9)
                    if c == 0:
                        p = phot[j]
                        lab = (f"#{j} ({x:.1f}, {y:.1f})  {'Gaia' if S['fromGaia'][j] else 'non-Gaia'}"
                               f"{', extended' if S['ext'][j] > 0 else ''}\n"
                               f"PSF  {p['psfMag']:.3f} ± {p['psfMagErr']:.3f} mag  "
                               f"(S/N {p['psfSnr']:.0f})\n"
                               f"ap   {p['apMag']:.3f} ± {p['apMagErr']:.3f} mag  "
                               f"(ap−PSF {p['apMinusPsf']:+.3f})\n"
                               + (f"Gaia G {p['gaiaG']:.2f}   " if np.isfinite(p["gaiaG"]) else "")
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
                                    a.axvline(off, color="C0" if S["fromGaia"][k] else "C3",
                                              ls=":" if S["fromGaia"][k] else "--", lw=0.9)
                    a1.set_title(ttl, fontsize=9)
                    a1.tick_params(labelbottom=False)
                    a2.set_xlabel(xl, fontsize=8)
                    if row == 0 and c == 0:
                        a1.legend(fontsize=7, loc="upper right")
            fig.suptitle(f"{run.name if hasattr(run, 'name') else ''}: source cutouts "
                         f"(zero point {phot.meta.get('zeroPoint', np.nan):.3f}, "
                         f"{phot.meta.get('zeroPointSource', '')}); grey band = ±1σ; "
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
            print(f"wrote {len(pages)} page(s) to {pdfPath}")
    return phot[sel]


# old snake_case names still work (with a FutureWarning)
from ._compat import moduleGetattr  # noqa: E402
__getattr__ = moduleGetattr(globals())
