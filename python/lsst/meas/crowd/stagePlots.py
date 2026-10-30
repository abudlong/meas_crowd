# This file is part of meas_crowd.
#
# Developed for the LSST Data Management System.
# See the COPYRIGHT file at the top-level directory of this distribution.

"""Plots of the pipeline's intermediate stages.

The PSF stars and the PSF model, the background entering each pass, data | model | residual
round by round, and per-source pages whose coloured frames match boxes drawn on the data, model
and residual images.
"""

__all__ = ["psf_star_positions", "show_psf_stars", "show_backgrounds", "show_rounds",
           "show_linked_source_pages", "compare_passes"]

import numpy as np

from . import fitLib as L
from .pipeline import _params, select_psf_stars


def _xy_from_catalog(cat, flag=None):
    """Return (N, 2) positions of a source catalogue's rows, optionally only flagged ones."""
    if hasattr(cat, "getX"):                        # a single afw record
        cat = [cat]
    try:                                            # afw SourceCatalog
        rows = [r for r in cat if flag is None or r.get(flag)]
        return np.array([(r.getX(), r.getY()) for r in rows], float).reshape(-1, 2)
    except AttributeError:
        pass
    cols = list(getattr(cat, "colnames", None) or getattr(cat, "columns", []))
    for cx, cy in (("x", "y"), ("slot_Centroid_x", "slot_Centroid_y"),
                   ("centroid_x", "centroid_y"), ("base_SdssCentroid_x", "base_SdssCentroid_y")):
        if cx in cols and cy in cols:
            x, y = np.asarray(cat[cx], float), np.asarray(cat[cy], float)
            break
    else:
        raise ValueError(f"no x/y centroid columns in the catalogue (have {cols[:10]}...)")
    keep = np.ones(len(x), bool) if flag is None else np.asarray(cat[flag], bool)
    return np.column_stack([x[keep], y[keep]])


def psf_star_positions(session=None, catalog=None, run=None, outer=-1):
    """Return the candidate and used PSF stars of one of the pipeline's PSF determinations.

    Exactly one source is used, in this order: a catalogue with ``calib_psf_candidate`` /
    ``calib_psf_used`` flags (the visit PSF's own stars), a session (the detector-level
    ``psf_source='refit'`` PSF), or a run (the clean stars the residual-stack PSF update after
    pass ``outer`` used).

    Parameters
    ----------
    session : `lsst.meas.crowd.session.Session`, optional
        A session prepared with ``psf_source='refit'``.
    catalog : `lsst.afw.table.SourceCatalog` or `astropy.table.Table`, optional
        A source catalogue with PSF flags, e.g. the visit's ``single_visit_star_unstandardized``.
    run : `lsst.meas.crowd.session.CutoutRun`, optional
        A finished run.
    outer : `int`, optional
        Outer pass of ``run`` whose result the PSF update used.

    Returns
    -------
    stars : `dict`
        ``candidates`` and ``used`` (positions in detector px, shape (N, 2)) and ``source``
        (a description).

    Raises
    ------
    ValueError
        Raised if no source of PSF stars is given.
    """
    if catalog is not None:
        try:                                        # afw SourceCatalog
            cols = set(catalog.schema.getNames())
        except AttributeError:                      # astropy Table / pandas DataFrame
            cols = set(getattr(catalog, "colnames", None) or list(getattr(catalog, "columns", [])))
        cand = _xy_from_catalog(catalog, "calib_psf_candidate" if "calib_psf_candidate" in cols
                                else None)
        used = (_xy_from_catalog(catalog, "calib_psf_used") if "calib_psf_used" in cols
                else np.zeros((0, 2)))
        return dict(candidates=cand, used=used, source="catalogue (calib_psf_* flags)")
    if session is not None:
        if session.starCat is None:
            raise ValueError("this session used the visit PSF (psf_source='visit'): pass the "
                             "visit's source catalogue as `catalog` instead")
        used = session.psf_used_xy
        return dict(candidates=_xy_from_catalog(session.starCat),
                    used=np.zeros((0, 2)) if used is None else np.asarray(used, float),
                    source="detector-level refit (ObjectSizeStarSelector + PCA PSF)")
    if run is not None:
        h = run.outer_history[outer]
        exp, res = h["exposure"], h["result"]
        P = _params(run.config)
        c = exp.getBBox().getCenter()
        kernel = h["psf"].computeKernelImage(c).array
        noise = float(np.sqrt(np.nanmedian(exp.variance.array)))
        fwhm = L._fwhm(h["psf"], c)
        xy = np.asarray(res["S"]["xy"], float).reshape(-1, 2)
        amps = np.asarray(res["snaps"][-1]["fit"]["amps"], float)
        cand = xy[amps * kernel.max() / noise >= P["psf_stack_min_snr"]]
        idx = select_psf_stars(exp, res, fwhm, kernel.max(), noise, P, verbose=False)
        n = len(run.outer_history)
        return dict(candidates=cand, used=xy[np.asarray(idx, int)],
                    source=f"residual-stack PSF update after pass {outer % n + 1}")
    raise ValueError("give a catalogue, a session or a run")


def _kernel_profile(kernel):
    """Return the azimuthally averaged radial profile of a kernel image (r [px], value)."""
    ny, nx = kernel.shape
    yy, xx = np.mgrid[:ny, :nx]
    r = np.hypot(xx - nx // 2, yy - ny // 2)
    rb = np.round(r).astype(int)
    prof = np.bincount(rb.ravel(), kernel.ravel()) / np.maximum(np.bincount(rb.ravel()), 1)
    return np.arange(len(prof)), prof


def show_psf_stars(image, stars=None, psf=None, at=None, view=None, figsize=(19, 6), out=None):
    """Plot the PSF stars over an image, the PSF kernel and its radial profile.

    Parameters
    ----------
    image : `lsst.afw.image.Exposure`
        The visit image (or a cutout) to draw the stars on.
    stars : `dict`, optional
        Output of `psf_star_positions`; None = draw no stars.
    psf : `lsst.afw.detection.Psf`, optional
        PSF to show (default: the image's own PSF).
    at : `tuple` [`float`, `float`], optional
        Detector position at which to render the kernel (default: the image centre).
    view : `tuple` [`float`, `float`, `float`, `float`], optional
        (xmin, xmax, ymin, ymax) zoom in detector pixels.
    figsize : `tuple` [`float`, `float`], optional
        Figure size [inches].
    out : `str`, optional
        Also save the figure to this file.

    Returns
    -------
    info : `dict`
        ``kernel`` (the rendered kernel image), ``fwhm`` [px], ``n_candidates`` and ``n_used``.
    """
    import matplotlib.pyplot as plt
    from astropy.visualization import simple_norm
    psf = psf or image.getPsf()
    bb = image.getBBox()
    x0, y0 = bb.getMinX(), bb.getMinY()
    c = bb.getCenter() if at is None else L.geom.Point2D(*at)
    kernel = np.asarray(psf.computeKernelImage(c).array, float)
    fwhm = L._fwhm(psf, c)
    data = np.asarray(image.image.array, float)
    ny, nx = data.shape
    ext = (x0 - 0.5, x0 + nx - 0.5, y0 - 0.5, y0 + ny - 0.5)
    stars = stars or dict(candidates=np.zeros((0, 2)), used=np.zeros((0, 2)), source="")
    cand, used = np.asarray(stars["candidates"]).reshape(-1, 2), np.asarray(stars["used"]).reshape(-1, 2)

    fig, ax = plt.subplots(1, 3, figsize=figsize, gridspec_kw=dict(width_ratios=(1.6, 1, 1)))
    fin = data[np.isfinite(data)]
    norm = simple_norm(fin, stretch="asinh", vmin=np.percentile(fin, 1),
                       vmax=np.percentile(fin, 99.7)) if fin.size else None
    ax[0].imshow(data, origin="lower", extent=ext, cmap="gray_r", norm=norm, interpolation="nearest")
    ax[0].scatter(cand[:, 0], cand[:, 1], s=60, facecolors="none", edgecolors="orange", lw=1,
                  label=f"candidates ({len(cand)})")
    ax[0].scatter(used[:, 0], used[:, 1], s=90, marker="s", facecolors="none", edgecolors="lime",
                  lw=1.4, label=f"used ({len(used)})")
    if view is not None:
        ax[0].set_xlim(view[0], view[1])
        ax[0].set_ylim(view[2], view[3])
    ax[0].legend(loc="upper right", fontsize=9)
    ax[0].set_title(f"PSF stars: {stars['source']}", fontsize=10)
    ax[0].set_xlabel("x [px]")
    ax[0].set_ylabel("y [px]")

    h = kernel.shape[0] // 2
    pos = kernel[kernel > 0]
    knorm = simple_norm(kernel, stretch="log", vmin=pos.max() * 1e-5 if pos.size else None)
    im = ax[1].imshow(kernel, origin="lower", extent=(-h - .5, h + .5, -h - .5, h + .5),
                      norm=knorm, cmap="magma")
    plt.colorbar(im, ax=ax[1], fraction=0.046)
    ax[1].set_title(f"PSF kernel at ({c[0]:.0f}, {c[1]:.0f}), {kernel.shape[1]}x{kernel.shape[0]} px, "
                    f"FWHM {fwhm:.2f} px", fontsize=10)
    ax[1].set_xlabel("dx [px]")

    r, prof = _kernel_profile(kernel)
    ok = prof > 0
    ax[2].semilogy(r[ok], prof[ok] / prof.max(), "k.-")
    ax[2].axvline(fwhm / 2, color="C1", ls="--", label=f"HWHM {fwhm / 2:.2f} px")
    ax[2].set_xlabel("radius [px]")
    ax[2].set_ylabel("kernel / peak")
    ax[2].set_title("radial profile (log)", fontsize=10)
    ax[2].legend()
    plt.tight_layout()
    if out:
        fig.savefig(out, dpi=110)
    return dict(kernel=kernel, fwhm=fwhm, n_candidates=len(cand), n_used=len(used))


def _image_array(img, bbox):
    """Return an afw image's pixels over ``bbox`` (it may cover a larger area)."""
    arr = getattr(img, "array", None)
    if arr is not None and arr.shape == (bbox.getHeight(), bbox.getWidth()):
        return np.asarray(arr, float)
    return np.asarray(img[bbox].array, float)


def show_backgrounds(run, initial=None, figsize_per_panel=4.2, out=None):
    """Plot the background entering every outer pass, and the change from pass to pass.

    Pass 1 uses ``initial`` (the detector-level background of ``psf_source='refit'``; a visit
    image is already background-subtracted, so there is none).  Before each later pass the
    background is refitted on the raw cutout with the previous pass's sources masked.

    Parameters
    ----------
    run : `lsst.meas.crowd.session.CutoutRun`
        A finished run.
    initial : `lsst.afw.image.Image`, optional
        The background subtracted before pass 1 (``session.bgImage``).
    figsize_per_panel : `float`, optional
        Size of each panel [inches].
    out : `str`, optional
        Also save the figure to this file.

    Returns
    -------
    backgrounds : `list` [`numpy.ndarray` or `None`]
        Background entering each pass (None = nothing was subtracted).
    """
    import matplotlib.pyplot as plt
    hist = run.outer_history
    bb = hist[0]["exposure"].getBBox()
    ext = (bb.getMinX() - .5, bb.getMaxX() + .5, bb.getMinY() - .5, bb.getMaxY() + .5)
    bgs = [None if initial is None else _image_array(initial, bb)]
    bgs += [None if h.get("bg_image") is None else _image_array(h["bg_image"], bb) for h in hist[1:]]
    n = len(bgs)
    fig, ax = plt.subplots(2, n, figsize=(figsize_per_panel * n, 2 * figsize_per_panel), squeeze=False)
    have = [b for b in bgs if b is not None]
    lo, hi = ((np.nanpercentile(np.concatenate([b.ravel() for b in have]), [0.5, 99.5]))
              if have else (0, 1))
    for k, b in enumerate(bgs):
        a = ax[0, k]
        if b is None:
            a.text(.5, .5, "nothing subtracted\n(visit image is already\nbackground-subtracted)"
                   if k == 0 else "no refit", ha="center", va="center", transform=a.transAxes)
            a.set_title(f"background entering pass {k + 1}", fontsize=10)
        else:
            im = a.imshow(b, origin="lower", extent=ext, vmin=lo, vmax=hi, cmap="viridis")
            plt.colorbar(im, ax=a, fraction=0.046)
            a.set_title(f"background entering pass {k + 1}\nmedian {np.nanmedian(b):.4g}, "
                        f"range {np.nanmin(b):.4g}..{np.nanmax(b):.4g}", fontsize=10)
        d = ax[1, k]
        prev = bgs[k - 1] if k else None
        if b is not None and prev is not None:
            diff = b - prev
            v = np.nanpercentile(np.abs(diff), 99) or 1.0
            im = d.imshow(diff, origin="lower", extent=ext, vmin=-v, vmax=v, cmap="RdBu_r")
            plt.colorbar(im, ax=d, fraction=0.046)
            d.set_title(f"change from pass {k}: rms {np.nanstd(diff):.3g}", fontsize=10)
        else:
            d.axis("off")
    for a in ax[0]:
        a.set_xticks([])
        a.set_yticks([])
    plt.tight_layout()
    if out:
        fig.savefig(out, dpi=110)
    for k, b in enumerate(bgs):
        print(f"pass {k + 1}: " + ("no background subtracted" if b is None else
                                   f"background median {np.nanmedian(b):.4g}, rms {np.nanstd(b):.3g}"))
    return bgs


def _mad_std(v):
    """Return the robust (MAD-based) standard deviation."""
    return 1.4826 * float(np.nanmedian(np.abs(v - np.nanmedian(v)))) if v.size else np.nan


def _pick(n, max_rows):
    """Return up to ``max_rows`` indices out of ``n``, evenly spaced, first and last included."""
    if n <= max_rows:
        return list(range(n))
    return sorted({int(round(v)) for v in np.linspace(0, n - 1, max_rows)})


def show_rounds(run, level="outer", outer=-1, rounds=None, max_rows=8, view=None, z_clip=5.0,
                figsize_per_panel=3.6, out=None):
    """Plot data | model | residual | residual/sigma for each round of the fit.

    ``level='outer'``: one row per outer pass (each with its own background and PSF).
    ``level='inner'``: one row per inner iteration of outer pass ``outer`` (the source list
    growing as the detection threshold steps down).  Data and model share one colour scale, and
    the residual one scale for all rows, so the rows can be compared.

    Parameters
    ----------
    run : `lsst.meas.crowd.session.CutoutRun`
        A finished run.
    level : `str`, optional
        'outer' or 'inner'.
    outer : `int`, optional
        Outer pass whose inner iterations are shown (``level='inner'``).
    rounds : `list` [`int`], optional
        Rounds to show (default: up to ``max_rows``, evenly spaced, the last included).
    max_rows : `int`, optional
        Maximum number of rows when ``rounds`` is not given.
    view : `tuple` [`float`, `float`, `float`, `float`], optional
        (xmin, xmax, ymin, ymax) zoom in detector pixels.
    z_clip : `float`, optional
        Colour range of residual/sigma.
    figsize_per_panel : `float`, optional
        Size of each panel [inches].
    out : `str`, optional
        Also save the figure to this file.

    Returns
    -------
    rows : `list` [`dict`]
        Per row ``label``, ``n_sources``, ``resid_std`` (robust std of the residual) and
        ``chi_std`` (robust std of residual/sigma) over the good pixels in view.

    Raises
    ------
    ValueError
        Raised if ``level`` is not 'outer' or 'inner'.
    """
    import matplotlib.pyplot as plt
    hist = run.outer_history
    items = []
    if level == "outer":
        for k, h in enumerate(hist):
            res = h["result"]
            data = np.asarray(h["exposure"].image.array, float)
            items.append((f"outer pass {k + 1}", h["exposure"], data, np.asarray(res["model_full"]),
                          len(res["S"]["xy"])))
    elif level == "inner":
        h = hist[outer]
        exp = h["exposure"]
        data = np.asarray(exp.image.array, float)
        for s in h["result"]["snaps"]:
            model = (data - s["residual"]) if "residual" in s else np.asarray(s["fit"]["model_full"])
            thr = s.get("threshold")
            lab = (f"pass {outer % len(hist) + 1}, it {s.get('it', 0) + 1}"
                   + (f" ({thr:.2f}σ)" if thr is not None else ""))
            items.append((lab, exp, data, model, len(s.get("merged_arr", s["fit"]["amps"]))))
    else:
        raise ValueError("level must be 'outer' or 'inner'")
    sel = rounds if rounds is not None else _pick(len(items), max_rows)
    items = [items[i] for i in sel]

    bb = items[0][1].getBBox()
    x0, y0 = bb.getMinX(), bb.getMinY()
    ny, nx = items[0][2].shape
    xa, xb, ya, yb = 0, nx, 0, ny
    if view is not None:
        xa, xb = max(int(view[0]) - x0, 0), min(int(view[1]) - x0 + 1, nx)
        ya, yb = max(int(view[2]) - y0, 0), min(int(view[3]) - y0 + 1, ny)
    sl = (slice(ya, yb), slice(xa, xb))
    ext = (x0 + xa - .5, x0 + xb - .5, y0 + ya - .5, y0 + yb - .5)
    bad = ((items[0][1].mask.array & run.bad_bitmask) != 0)[sl]
    var = np.asarray(items[0][1].variance.array, float)[sl]
    norm = L._data_norm(np.where(bad, np.nan, items[-1][2][sl]))
    last_r = (items[-1][2] - items[-1][3])[sl]
    rlim = np.nanpercentile(np.abs(last_r[~bad]), 99) if (~bad).any() else 1.0

    fig, ax = plt.subplots(len(items), 4, figsize=(4 * figsize_per_panel, figsize_per_panel * len(items)),
                           squeeze=False)
    rows = []
    for r, (lab, exp, data, model, n) in enumerate(items):
        d, m = data[sl], model[sl]
        resid = d - m
        chi = resid / np.sqrt(np.maximum(var, 1e-12))
        good = ~bad & np.isfinite(resid)
        st, cs = _mad_std(resid[good]), _mad_std(chi[good])
        rows.append(dict(label=lab, n_sources=n, resid_std=float(st), chi_std=float(cs)))
        for c, (img, kw, t) in enumerate((
                (d, dict(norm=norm, cmap=L.CMAP), "data"),
                (m, dict(norm=norm, cmap=L.CMAP), f"model ({n} sources)"),
                (resid, dict(vmin=-rlim, vmax=rlim, cmap="RdBu_r"), f"residual (std {st:.3g})"),
                (chi, dict(vmin=-z_clip, vmax=z_clip, cmap="RdBu_r"), f"residual/σ (std {cs:.2f})"))):
            a = ax[r, c]
            im = a.imshow(np.where(bad, np.nan, img), origin="lower", extent=ext,
                          interpolation="nearest", **kw)
            plt.colorbar(im, ax=a, fraction=0.046)
            a.set_title(f"{lab}: {t}" if c == 0 else t, fontsize=9)
            a.set_xticks([])
            a.set_yticks([])
    plt.tight_layout()
    if out:
        fig.savefig(out, dpi=100)
    for row in rows:
        print(f"{row['label']:>26s}: {row['n_sources']:5d} sources, residual std "
              f"{row['resid_std']:.4g}, residual/σ std {row['chi_std']:.3f}")
    return rows


def _nearest(xy, positions):
    """Return the index of the source nearest each position."""
    xy = np.asarray(xy, float).reshape(-1, 2)
    return [int(np.argmin(np.hypot(xy[:, 0] - x, xy[:, 1] - y))) for x, y in positions]


def _pass_sources(h):
    """Return positions, fluxes, peak S/N and Gaia flags of one outer pass's sources."""
    res, exp = h["result"], h["exposure"]
    xy = np.asarray(res["S"]["xy"], float).reshape(-1, 2)
    amps = np.asarray(res["snaps"][-1]["fit"]["amps"], float)
    peak = float(h["psf"].computeKernelImage(exp.getBBox().getCenter()).array.max())
    noise = float(np.sqrt(np.nanmedian(exp.variance.array)))
    seed = np.asarray(res["S"].get("seed", np.zeros(len(xy), bool)), bool)
    return xy, amps, amps * peak / noise, seed


def _matched_flux(resid, var, good, kernel, x, y, x0, y0):
    """Return the PSF-weighted flux and its S/N of the residual at a position."""
    k = kernel / kernel.sum()
    hy, hx = k.shape[0] // 2, k.shape[1] // 2
    ix, iy = int(round(x)) - x0, int(round(y)) - y0
    ny, nx = resid.shape
    ya, yb, xa, xb = max(iy - hy, 0), min(iy + hy + 1, ny), max(ix - hx, 0), min(ix + hx + 1, nx)
    if yb <= ya or xb <= xa:
        return np.nan, np.nan
    kk = k[ya - iy + hy:yb - iy + hy, xa - ix + hx:xb - ix + hx]
    w = np.where(good[ya:yb, xa:xb], kk / np.maximum(var[ya:yb, xa:xb], 1e-12), 0.0)
    den = float(np.sum(w * kk))
    if den <= 0:
        return np.nan, np.nan
    f = float(np.sum(w * np.nan_to_num(resid[ya:yb, xa:xb]))) / den
    return f, f * np.sqrt(den)


def compare_passes(run, a=-2, b=-1, match_px=1.0, n_show=8, half=12, plot=True, out=None):
    """Find the sources that one outer pass lost or gained relative to another, and why.

    Sources are matched by position (within ``match_px``).  For every source of pass ``a``
    without a match in pass ``b`` the residual of pass ``b`` is measured at its position with a
    PSF-weighted (matched-filter) flux: if most of its flux is still there, a real source was
    lost; if a neighbour within 1.5 FWHM gained that flux, it was merged; if a halo or ring check
    removed it, it is also blocked from re-detection in later passes.

    Parameters
    ----------
    run : `lsst.meas.crowd.session.CutoutRun`
        A finished run with at least two outer passes.
    a, b : `int`, optional
        The two outer passes to compare (default: the last two).
    match_px : `float`, optional
        Matching radius [px].
    n_show : `int`, optional
        Number of lost sources drawn (brightest first).
    half : `int`, optional
        Half-size of each zoom [px].
    plot : `bool`, optional
        Draw the overview and the zooms.
    out : `str`, optional
        Also save the figure to this file.

    Returns
    -------
    table : `astropy.table.Table`
        One row per lost or gained source: ``kind`` ('lost' / 'gained'), ``x``, ``y``, ``flux``,
        ``snr`` (peak S/N in the pass it exists in), ``gaia``, ``flux_left`` and ``left_snr``
        (matched-filter flux and S/N in the residual of the pass without it), ``left_frac``,
        ``nearest_px`` (distance to the nearest source of the other pass), ``nearest_dflux``
        (that source's flux change) and ``verdict``.

    Raises
    ------
    ValueError
        Raised if the run has fewer than two outer passes.
    """
    from astropy.table import Table
    from scipy.spatial import cKDTree
    hist = run.outer_history
    if len(hist) < 2:
        raise ValueError("the run has only one outer pass")
    ia, ib = a % len(hist), b % len(hist)
    ha, hb = hist[ia], hist[ib]
    xa_, fa, sa, ga = _pass_sources(ha)
    xb_, fb, sb, gb = _pass_sources(hb)
    c = hb["exposure"].getBBox().getCenter()
    fwhm = L._fwhm(hb["psf"], c)
    da, _ = cKDTree(xb_).query(xa_) if len(xb_) else (np.full(len(xa_), np.inf), None)
    db, _ = cKDTree(xa_).query(xb_) if len(xa_) else (np.full(len(xb_), np.inf), None)
    rej_a = np.asarray(ha.get("rejects", np.zeros((0, 2))), float).reshape(-1, 2)
    rej_b = np.asarray(hb.get("rejects", np.zeros((0, 2))), float).reshape(-1, 2)
    new_rej = rej_b[len(rej_a):] if len(rej_b) >= len(rej_a) else rej_b

    rows = []
    for kind, src, other, flux, snr, gaia, dist, h_other in (
            ("lost", xa_, (xb_, fb, xa_, fa), fa, sa, ga, da, hb),
            ("gained", xb_, (xa_, fa, xb_, fb), fb, sb, gb, db, ha)):
        exp = h_other["exposure"]
        bb = exp.getBBox()
        resid = np.asarray(exp.image.array, float) - np.asarray(h_other["result"]["model_full"], float)
        var = np.asarray(exp.variance.array, float)
        good = ((exp.mask.array & run.bad_bitmask) == 0) & np.isfinite(resid)
        kern = np.asarray(h_other["psf"].computeKernelImage(bb.getCenter()).array, float)
        oxy, oflux, sxy, sflux = other
        tree = cKDTree(oxy) if len(oxy) else None
        stree = cKDTree(sxy) if len(sxy) else None
        for j in np.flatnonzero(dist > match_px):
            x, y = src[j]
            fl, fs = _matched_flux(resid, var, good, kern, x, y, bb.getMinX(), bb.getMinY())
            near, dfl = np.inf, np.nan
            if tree is not None:
                near, k = tree.query((x, y))
                # flux change of that neighbour between the passes
                if stree is not None:
                    d2, k2 = stree.query(oxy[k])
                    if d2 <= match_px:            # always pass b minus pass a
                        dfl = (oflux[k] - sflux[k2]) * (1 if kind == "lost" else -1)
            frac = fl / flux[j] if flux[j] > 0 else np.nan
            if kind == "lost":
                blocked = len(new_rej) and np.min(np.hypot(*(new_rej - (x, y)).T)) <= match_px
                if blocked:
                    verdict = "removed by halo/ring check; blocked from re-detection"
                elif np.isfinite(fs) and fs >= 5 and frac > 0.5:
                    verdict = "LOST: its flux is still in the residual"
                elif near <= 1.5 * fwhm and np.isfinite(dfl) and dfl >= 0.5 * flux[j]:
                    verdict = "merged into a neighbour"
                elif near <= match_px * 3:
                    verdict = "moved"
                else:
                    verdict = "flux not needed (likely spurious)"
            else:
                verdict = ("new: was missing from the residual" if np.isfinite(fs) and fs >= 5
                           else "new: faint / split off")
            rows.append((kind, x, y, flux[j], snr[j], bool(gaia[j]), fl, fs, frac, near, dfl,
                         verdict))
    t = Table(rows=rows, names=("kind", "x", "y", "flux", "snr", "gaia", "flux_left", "left_snr",
                                "left_frac", "nearest_px", "nearest_dflux", "verdict"),
              dtype=(str, float, float, float, float, bool, float, float, float, float, float, str))
    t.sort(["kind", "flux"], reverse=True)
    lost, gained = t[t["kind"] == "lost"], t[t["kind"] == "gained"]
    print(f"pass {ia + 1} -> pass {ib + 1}: {len(xa_)} -> {len(xb_)} sources; {len(lost)} lost, "
          f"{len(gained)} gained (match radius {match_px} px)")
    for v in sorted(set(lost["verdict"])):
        print(f"  lost, {v}: {int(np.sum(lost['verdict'] == v))}")
    for r in lost:
        print(f"    ({r['x']:7.1f}, {r['y']:7.1f})  flux {r['flux']:9.1f}  S/N {r['snr']:6.1f}  "
              f"{'Gaia' if r['gaia'] else '    '}  left {100 * r['left_frac']:5.0f}% "
              f"(S/N {r['left_snr']:5.1f})  nearest {r['nearest_px']:4.1f} px  -> {r['verdict']}")

    if plot and len(t):
        import matplotlib.pyplot as plt
        expb = hb["exposure"]
        bb = expb.getBBox()
        ext = (bb.getMinX() - .5, bb.getMaxX() + .5, bb.getMinY() - .5, bb.getMaxY() + .5)
        bad = (expb.mask.array & run.bad_bitmask) != 0
        rb = np.asarray(expb.image.array, float) - np.asarray(hb["result"]["model_full"], float)
        ra = (np.asarray(ha["exposure"].image.array, float)
              - np.asarray(ha["result"]["model_full"], float))
        rl = np.nanpercentile(np.abs(rb[~bad]), 99) if (~bad).any() else 1.0
        show = lost[:n_show]
        nrow = 1 + len(show)
        fig = plt.figure(figsize=(15, 4.6 * nrow))
        gs = fig.add_gridspec(nrow, 3)
        ax = fig.add_subplot(gs[0, :])
        ax.imshow(np.where(bad, np.nan, rb), origin="lower", extent=ext, cmap="RdBu_r",
                  vmin=-rl, vmax=rl, interpolation="nearest")
        ax.scatter(lost["x"], lost["y"], marker="x", c="k", s=60, lw=1.5, label=f"lost ({len(lost)})")
        ax.scatter(gained["x"], gained["y"], marker="+", c="g", s=70, lw=1.5,
                   label=f"gained ({len(gained)})")
        ax.legend(loc="upper right")
        ax.set_title(f"residual of pass {ib + 1}; sources lost / gained since pass {ia + 1}")
        ax.set_aspect("equal")
        norm = L._data_norm(np.where(bad, np.nan, np.asarray(ha["exposure"].image.array, float)))
        for r, row in enumerate(show, 1):
            x, y = row["x"], row["y"]
            win = (x - half, x + half, y - half, y + half)
            for cc, (img, kw, ttl) in enumerate((
                    (np.asarray(ha["exposure"].image.array, float), dict(norm=norm, cmap=L.CMAP),
                     f"pass {ia + 1} data: flux {row['flux']:.0f}, S/N {row['snr']:.0f}"),
                    (ra, dict(cmap="RdBu_r", vmin=-rl, vmax=rl), f"pass {ia + 1} residual (with it)"),
                    (rb, dict(cmap="RdBu_r", vmin=-rl, vmax=rl),
                     f"pass {ib + 1} residual: {100 * row['left_frac']:.0f}% left"))):
                a_ = fig.add_subplot(gs[r, cc])
                a_.imshow(np.where(bad, np.nan, img), origin="lower", extent=ext,
                          interpolation="nearest", **kw)
                a_.scatter([x], [y], marker="x", c="k", s=60, lw=1.5)
                a_.scatter(xb_[:, 0], xb_[:, 1], marker="+", c="w" if cc == 0 else "g", s=40, lw=1)
                a_.set_xlim(win[0], win[1])
                a_.set_ylim(win[2], win[3])
                a_.set_title(ttl if cc else f"({x:.1f}, {y:.1f}) {ttl}", fontsize=9)
            fig.axes[-1].set_xlabel(row["verdict"], fontsize=9)
        plt.tight_layout()
        if out:
            fig.savefig(out, dpi=100)
    return t


class _PngPages:
    """Stand-in for `PdfPages` that writes every page to its own PNG file at once."""

    def __init__(self, out_dir, prefix):
        import os
        os.makedirs(out_dir, exist_ok=True)
        self.out_dir, self.prefix, self.n, self.paths = out_dir, prefix, 0, []

    def savefig(self, fig, dpi=80, **kw):
        """Save one page and remember its file name."""
        import os
        self.n += 1
        path = os.path.join(self.out_dir, f"{self.prefix}_{self.n:03d}.png")
        fig.savefig(path, dpi=dpi, **kw)
        self.paths.append(path)


def show_linked_source_pages(run, outer=-1, snap=-1, sources=None, source_view="all",
                             sources_per_page=6, max_pages=1, view_half=20, show_pages=True,
                             out_dir=None, pdf=None, dpi=80, figsize=(21, 7.5), out=None):
    """Plot per-source pages whose coloured frames match boxes on the data, model and residual.

    For each page: an overview (data | model | residual of the whole cutout) with one coloured
    box per source on that page -- the window its row shows -- then the page itself: one row
    per source with data | model | residual | x cross-section | y cross-section, framed in the
    same colour.

    For many pages use ``out_dir``: every figure is written to its own PNG file and closed at
    once, so memory stays flat.  A `PdfPages` (``pdf``) keeps every image in memory until it is
    closed (matplotlib), which can exhaust a notebook's memory for hundreds of sources.

    Parameters
    ----------
    run : `lsst.meas.crowd.session.CutoutRun`
        A finished run.
    outer : `int`, optional
        Outer pass.
    snap : `int`, optional
        Inner-iteration snapshot of that pass; only the last one keeps the per-source stamps
        unless the run used ``keep_all_fits=True``.
    sources : `list`, optional
        Source indices, or (x, y) detector positions (the nearest source is used); None =
        choose with ``source_view``.
    source_view : `str`, optional
        'all' (brightest first), 'new' or 'changed' (see `fitLib.select_source_pages`).
    sources_per_page : `int`, optional
        Rows per page.
    max_pages : `int`, optional
        Maximum number of pages (None = all).
    view_half : `int`, optional
        Half-size of each source window [px] (the size of the boxes).
    show_pages : `bool`, optional
        Draw the pages (False = only the overviews).
    out_dir : `str`, optional
        Write ``overview_NNN.png`` and ``page_NNN.png`` here instead of showing them.
    pdf : `matplotlib.backends.backend_pdf.PdfPages`, optional
        Open PDF; the overviews and pages are saved into it instead of being shown.
    dpi : `int`, optional
        Resolution of the saved pages.
    figsize : `tuple` [`float`, `float`], optional
        Size of the overview [inches].
    out : `str`, optional
        Also save the first overview to this file.

    Returns
    -------
    linked : `dict`
        ``pages`` (source indices per page), ``colors`` (colour per source index) and ``files``
        (the PNG files written with ``out_dir``).

    Raises
    ------
    ValueError
        Raised if the snapshot no longer has the per-source stamps.
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    hist = run.outer_history
    h = hist[outer]
    exp = h["exposure"]
    sn = h["result"]["snaps"][snap]
    fit = sn["fit"]
    if "comps" not in fit:
        raise ValueError("this snapshot has no per-source stamps any more: use snap=-1 or rerun "
                         "with keep_all_fits=True")
    xy = np.asarray(sn.get("merged_arr", h["result"]["S"]["xy"]), float)
    sn = dict(sn, merged_arr=xy)
    sn.setdefault("prev_amp", np.full(len(xy), np.nan))
    sn.setdefault("is_new", np.zeros(len(xy), bool))
    sn.setdefault("bad_mask", (exp.mask.array & run.bad_bitmask) != 0)
    sn.setdefault("it", len(h["result"]["snaps"]) - 1 if snap == -1 else snap)
    if sources is not None:
        idx = [int(s) if np.ndim(s) == 0 else _nearest(xy, [s])[0] for s in sources]
        pages = [idx[i:i + sources_per_page] for i in range(0, len(idx), sources_per_page)]
        pages = pages[:max_pages] if max_pages is not None else pages
    else:
        pages = L.select_source_pages(sn, source_view, sources_per_page, max_pages)
    colors = L.source_colors(pages)

    bb = exp.getBBox()
    data = np.asarray(exp.image.array, float)
    model = np.asarray(fit["model_full"], float)
    bad = sn["bad_mask"]
    ext = (bb.getMinX() - .5, bb.getMaxX() + .5, bb.getMinY() - .5, bb.getMaxY() + .5)
    norm = L._data_norm(np.where(bad, np.nan, data))
    resid = data - model
    rlim = np.nanpercentile(np.abs(resid[~bad]), 99) if (~bad).any() else 1.0
    radius = L._circle_radius(exp)
    panels = ((data, dict(norm=norm, cmap=L.CMAP), "data"),
              (model, dict(norm=norm, cmap=L.CMAP), f"model ({len(xy)} sources)"),
              (resid, dict(vmin=-rlim, vmax=rlim, cmap="RdBu_r"), "residual"))
    png = _PngPages(out_dir, "page") if out_dir else None
    ovw = _PngPages(out_dir, "overview") if out_dir else None
    sink = png or pdf
    n_pass = outer % len(hist) + 1
    for p, page in enumerate(pages):
        fig, ax = plt.subplots(1, 3, figsize=figsize)
        for a, (img, kw, t) in zip(ax, panels):
            im = a.imshow(np.where(bad, np.nan, img), origin="lower", extent=ext,
                          interpolation="nearest", **kw)
            plt.colorbar(im, ax=a, fraction=0.046)
            for j in page:
                px, py = xy[j]
                ix, iy, xa, xb, ya, yb = L._view_window(bb, px, py, view_half)
                a.add_patch(patches.Rectangle((bb.getMinX() + xa - .5, bb.getMinY() + ya - .5),
                                              xb - xa, yb - ya, fill=False, ec=colors[j], lw=2,
                                              path_effects=L._OUTLINE))
                L._mark_source(a, px, py, colors[j], radius, label=j)
            a.set_title(f"pass {n_pass}, page {p + 1}/{len(pages)}: {t}; boxes = rows on the page",
                        fontsize=10)
        plt.tight_layout()
        if out and p == 0:
            fig.savefig(out, dpi=110)
        if sink is not None:
            (ovw or pdf).savefig(fig, dpi=dpi)
            plt.close(fig)
        else:
            plt.show()
        if show_pages:
            L.plot_source_pages(exp, sn, n_pass - 1, pages=[page], colors=colors,
                                view_half=view_half, pdf=sink, show=sink is None, dpi=dpi)
    files = sorted((png.paths if png else []) + (ovw.paths if ovw else []))
    if out_dir:
        print(f"{len(files)} files in {out_dir}/ ({len(pages)} page(s))")
    return dict(pages=pages, colors=colors, files=files)
