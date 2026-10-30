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

"""Residual check that decides whether sources are still missing.

A matched-filter S/N map of the residual is searched for local maxima and minima in a clean
region (away from the edge, bright stars, bad pixels and fitted sources).  Real missing sources
only add positive peaks, while noise gives positive and negative ones equally, so the excess
``N_pos - N_neg`` above an S/N threshold (error ``sqrt(N_pos + N_neg)``) estimates how many
sources are missing; the negative count is measured in the same field, so the test adapts to
crowding and correlated noise.  `residual_excess`, `is_missing`, `counts_at` and
`pick_threshold` are used by the pipeline's inner loop.

`install` adds the same check, with a ladder of relaxed thresholds, to
`lsst.meas.crowd.fitLib.run_cutout_pipeline`.
"""
import copy
import inspect

import numpy as np
from scipy.ndimage import distance_transform_edt, maximum_filter
from scipy.signal import fftconvolve
from scipy.spatial import cKDTree

ADAPT = dict(
    enabled=True,
    max_rounds=3,          # max extra (warm-started) rounds per outer pass
    extra_inner=3,         # inner iterations per extra round
    target_snr=None,       # S/N down to which the residual should be clean (None: run's min_snr)
    t_floor=3.0,           # never relax thresholds below this S/N (set = target_snr to forbid going deeper)
    fdr=0.25,              # tolerated  N_neg / N_pos  when choosing the lowest threshold
    z_sig=2.0,             # excess must exceed z_sig * sqrt(N_pos + N_neg) ...
    min_excess=3,          # ... and be at least this many sources
    edge_fwhm=2.0,         # clean region: farther than this from the cutout edge [FWHM]
    halo_fwhm=2.5,         # clean region: farther than this from bright sources [FWHM] (grows with flux)
    bright_sigma=100.0,    # 'bright' = amp*psf_peak/noise above this
    min_sep_fwhm=0.75,     # peaks closer than this to a fitted source are its own mismatch residual
    max_outer=6,           # hard cap on outer passes when auto-extending the outer loop
    min_outer_growth=0.02, # extend the outer loop only if the last pass added >= this fraction of sources
    viz_rounds=False,      # draw the per-iteration figures during the extra rounds
    verbose=True,
)
_STATE = dict(npeaks=[])

_DEFAULT_SHAPE = dict(min_snr_local=5.0, chi2_max=2.0, width_range=(0.75, 1.35),
                      dchi2_width=9.0, trust_snr=30.0, min_contrast=20.0,
                      weak_snr=10.0, dchi2_width_weak=4.0)


# =====================================================================
# Diagnostic (numpy / scipy only)
# =====================================================================
def _odd(x):
    """Return the nearest odd integer >= 3."""
    return max(3, int(round(x)) | 1)


def _matched_snr(resid, var, weight_ok, fwhm, calib):
    """Return the Gaussian matched-filter S/N of the residual.

    Pixels with ``weight_ok`` False carry zero weight.  The S/N is calibrated empirically (median
    offset, MAD scale >= 1) on the ``calib`` pixels so correlated noise does not inflate it.
    Returns (snr, den, scale, offset, ok).
    """
    sig = fwhm / 2.355
    h = int(np.ceil(3 * sig))
    u = np.arange(-h, h + 1)
    k = np.exp(-0.5 * (u[None, :] ** 2 + u[:, None] ** 2) / sig ** 2)
    good = weight_ok & np.isfinite(resid) & np.isfinite(var) & (var > 0)
    inv = np.where(good, 1.0 / np.where(good, var, 1.0), 0.0)
    num = fftconvolve(np.where(good, resid, 0.0) * inv, k, mode="same")
    den = np.sqrt(np.clip(fftconvolve(inv, k ** 2, mode="same"), 1e-30, None))
    snr = num / den
    ok = good & (den > 0.3 * np.median(den[good]))
    ref = calib & ok
    zz = snr[ref] if ref.sum() > 200 else snr[ok]
    med = float(np.median(zz))
    scale = max(1.0, float(1.4826 * np.median(np.abs(zz - med))))
    return np.where(ok, (snr - med) / scale, 0.0), den, scale, med, ok


def residual_excess(resid, var, good, peaks_xy, amps, x0, y0, fwhm, psf_peak, noise,
                    thresholds=(3.0, 3.5, 4.0, 4.5, 5.0, 6.0), edge_fwhm=2.0,
                    bright_sigma=100.0, halo_fwhm=2.5, min_sep_fwhm=0.75):
    """Count positive and negative peaks of the residual's matched-filter S/N.

    Only local maxima in the clean region count: away from the cutout edge, from bright sources
    (halo radius growing with flux) and from fitted sources themselves.  Real missing sources are
    positive only, so the positive excess over the negative count estimates how many are
    missing.

    Parameters
    ----------
    resid : `numpy.ndarray`
        Residual image.
    var : `numpy.ndarray`
        Variance image.
    good : `numpy.ndarray` [`bool`]
        Usable pixels.
    peaks_xy : array-like
        Fitted source positions (parent pixel coordinates), shape (N, 2).
    amps : array-like
        Fitted source amplitudes.
    x0, y0 : `int`
        Parent pixel coordinates of ``resid[0, 0]``.
    fwhm : `float`
        PSF FWHM [px].
    psf_peak : `float`
        Peak value of the unit-flux PSF.
    noise : `float`
        Typical pixel noise.
    thresholds : `tuple` [`float`], optional
        S/N thresholds at which to count.
    edge_fwhm : `float`, optional
        Ignore peaks within this many FWHM of the edge.
    bright_sigma : `float`, optional
        A source is bright if ``amp * psf_peak / noise`` exceeds this.
    halo_fwhm : `float`, optional
        Ignore peaks within this many FWHM of a bright source (grows with flux).
    min_sep_fwhm : `float`, optional
        Peaks closer than this to a fitted source are its own mismatch residual.

    Returns
    -------
    result : `dict`
        ``t`` (thresholds, ascending), ``pos`` / ``neg`` (number of peaks with |S/N| >= t),
        ``scale``, ``offset`` (S/N calibration), ``area_frac`` (fraction of the image in the
        clean region) and ``candidates`` ((x, y, snr, distance to the nearest fitted source) of
        the positive peaks).
    """
    ny, nx = resid.shape
    iy, ix = np.mgrid[0:ny, 0:nx]
    P = np.asarray(peaks_xy, float).reshape(-1, 2) - np.array([x0, y0])
    amps = np.asarray(amps, float)
    ts = np.array(sorted({round(float(t), 3) for t in thresholds}))

    region = good & (np.minimum.reduce([ix, nx - 1 - ix, iy, ny - 1 - iy]) >= edge_fwhm * fwhm)
    if len(P):
        sig_br = amps * psf_peak / noise
        for j in np.flatnonzero(sig_br > bright_sigma):
            r = halo_fwhm * fwhm * (1.0 + 0.5 * max(0.0, np.log10(sig_br[j] / bright_sigma)))
            region &= (ix - P[j, 0]) ** 2 + (iy - P[j, 1]) ** 2 > r * r
    if not good.all():
        region &= distance_transform_edt(good) >= 0.5 * fwhm
    snr, den, scale, offset, ok = _matched_snr(resid, var, good, fwhm, region)
    region &= ok

    size = _odd(fwhm)
    tree = cKDTree(P) if len(P) else None

    def extrema(sign):
        s = sign * snr
        m = np.where(region, s, -np.inf)
        pk = (m == maximum_filter(m, size=size, mode="constant", cval=-np.inf)) & (m >= ts[0])
        yy, xx = np.nonzero(pk)
        v = s[yy, xx]
        nn = np.full(len(v), np.inf)
        if tree is not None and len(v):
            nn = tree.query(np.column_stack([xx, yy]))[0]
            k = nn >= min_sep_fwhm * fwhm
            xx, yy, v, nn = xx[k], yy[k], v[k], nn[k]
        return xx, yy, v, nn

    xp, yp, vp, dp = extrema(+1)
    _, _, vn, _ = extrema(-1)
    order = np.argsort(-vp)
    return dict(t=ts,
                pos=np.array([int((vp >= t).sum()) for t in ts]),
                neg=np.array([int((vn >= t).sum()) for t in ts]),
                scale=scale, offset=offset, area_frac=float(region.mean()),
                candidates=[(float(xp[i] + x0), float(yp[i] + y0), float(vp[i]), float(dp[i]))
                            for i in order])


def _idx(d, t):
    """Return the index of the threshold in ``d['t']`` closest to ``t``."""
    return int(np.argmin(np.abs(d["t"] - t)))


def counts_at(d, t):
    """Return the positive and negative peak counts at an S/N threshold.

    Parameters
    ----------
    d : `dict`
        Output of `residual_excess`.
    t : `float`
        S/N threshold (the nearest one in ``d['t']`` is used).

    Returns
    -------
    n_pos : `int`
        Number of positive peaks.
    n_neg : `int`
        Number of negative peaks.
    """
    i = _idx(d, t)
    return int(d["pos"][i]), int(d["neg"][i])


def is_missing(d, t, z_sig=2.0, min_excess=3):
    """Decide whether the residual still contains missing sources at an S/N threshold.

    Parameters
    ----------
    d : `dict`
        Output of `residual_excess`.
    t : `float`
        S/N threshold.
    z_sig : `float`, optional
        The excess must exceed ``z_sig * sqrt(N_pos + N_neg)`` ...
    min_excess : `int`, optional
        ... and be at least this many sources.

    Returns
    -------
    missing : `bool`
        Whether sources are missing.
    excess : `int`
        Positive minus negative peak count.
    error : `float`
        Its uncertainty, sqrt(N_pos + N_neg).
    """
    p, n = counts_at(d, t)
    e, err = p - n, float(np.sqrt(p + n))
    return bool(e >= max(min_excess, z_sig * err)), e, err


def pick_threshold(d, fdr, t_floor, t_cap):
    """Return the lowest S/N threshold the field can support.

    The lowest S/N in [``t_floor``, ``t_cap``] whose estimated false fraction
    (N_neg + 1) / N_pos is at most ``fdr`` (with N_pos >= 3).

    Parameters
    ----------
    d : `dict`
        Output of `residual_excess`.
    fdr : `float`
        Tolerated false fraction.
    t_floor : `float`
        Lowest threshold allowed.
    t_cap : `float`
        Highest threshold; returned if no lower one qualifies.

    Returns
    -------
    threshold : `float`
        The chosen S/N threshold.
    """
    for t, p, n in zip(d["t"], d["pos"], d["neg"]):
        if t < t_floor - 1e-9 or t > t_cap + 1e-9:
            continue
        if p >= 3 and (n + 1.0) / p <= fdr:
            return float(t)
    return float(t_cap)


# =====================================================================
# Relaxation ladder.  Each rung: (cur_kwargs, t_new, target, A) ->
#   (new_kwargs, detection_threshold, eval_snr, description)  or  (None, ...) if a no-op
# =====================================================================
def _copy(cur):
    """Copy the inner-loop keyword arguments (with their shape_kwargs dict)."""
    new = dict(cur)
    if isinstance(cur.get("shape_kwargs"), dict):
        new["shape_kwargs"] = dict(cur["shape_kwargs"])
    return new


def _shape(cur):
    """Return the shape-cut keyword arguments with the defaults filled in."""
    sk = dict(_DEFAULT_SHAPE)
    sk.update(cur.get("shape_kwargs") or {})
    return sk


def _rung_thresholds(cur, t_new, target, A):
    """Relaxation rung: lower the S/N cuts and the detection threshold."""
    old, thr0 = cur["min_snr"], cur["threshold_min"]
    if t_new >= old - 1e-9 and t_new >= thr0 - 1e-9:
        return None, None, None, None
    new, ch = _copy(cur), []
    if t_new < old - 1e-9:
        new["min_snr"] = t_new
        new["near_bright_snr"] = max(t_new, cur["near_bright_snr"] * t_new / old)
        sk = _shape(cur)
        sk["min_snr_local"] = min(sk["min_snr_local"], t_new)
        new["shape_kwargs"] = sk
        ch.append(f"min_snr {old:g}->{t_new:g}")
    thr = min(thr0, t_new)
    ch.append(f"detect to {thr:g} sigma")
    return new, thr, min(target, t_new), "lower thresholds (" + ", ".join(ch) + ")"


def _rung_gates(cur, t_new, target, A):
    """Relaxation rung: loosen the shape-cut gates and the bright-star rejections."""
    new, sk = _copy(cur), _shape(cur)
    sk["chi2_max"] *= 1.5
    sk["dchi2_width"] *= 1.5
    sk["dchi2_width_weak"] *= 1.5
    if sk.get("min_contrast"):
        sk["min_contrast"] = min(sk["min_contrast"], 1.0)
    new["shape_kwargs"] = sk
    new["reject_fwhm"] = cur["reject_fwhm"] * 0.7
    new["near_bright_snr"] = max(cur["min_snr"], 0.75 * cur["near_bright_snr"])
    return (new, cur["threshold_min"], target,
            f"relax gates (chi2_max {sk['chi2_max']:.2f}, dchi2_width {sk['dchi2_width']:.1f}, "
            f"reject_fwhm {new['reject_fwhm']:.2f}, near_bright_snr {new['near_bright_snr']:.1f})")


def _rung_separations(cur, t_new, target, A):
    """Relaxation rung: allow closer pairs (merge / prune / reject separations)."""
    new = _copy(cur)
    new["merge_separation_sigma"] = max(0.8, cur["merge_separation_sigma"] * 0.8)
    new["prune_sep_sigma"] = max(0.8, cur["prune_sep_sigma"] * 0.8)
    new["reject_fwhm"] = cur["reject_fwhm"] * 0.7
    return (new, cur["threshold_min"], target,
            f"relax separations (merge {new['merge_separation_sigma']:.2f} sigma, "
            f"prune {new['prune_sep_sigma']:.2f} sigma, reject_fwhm {new['reject_fwhm']:.2f})")


_RUNGS = [_rung_thresholds, _rung_gates, _rung_separations]


# =====================================================================
# Controller
# =====================================================================
def _diagnose(exposure, psf, res, A, extra=()):
    """Run `residual_excess` on the last snapshot of an inner-loop result."""
    snap = res["snaps"][-1]
    fit = snap["fit"]
    bbox = exposure.getBBox()
    c = bbox.getCenter()
    shp = psf.computeShape(c)
    fwhm = 2.355 * np.sqrt(0.5 * (shp.getIxx() + shp.getIyy()))
    psf_peak = psf.computeKernelImage(c).array.max()
    var, good = fit["variance"], fit["good"]
    noise = float(np.sqrt(np.nanmedian(var[good])))
    resid = exposure.image.array - fit["model_full"]
    ts = [3.0, 3.5, 4.0, 4.5, 5.0, 6.0] + [float(e) for e in extra]
    return residual_excess(resid, var, good, snap["merged_arr"], fit["amps"],
                           bbox.getMinX(), bbox.getMinY(), fwhm, psf_peak, noise,
                           thresholds=ts, edge_fwhm=A["edge_fwhm"],
                           bright_sigma=A["bright_sigma"], halo_fwhm=A["halo_fwhm"],
                           min_sep_fwhm=A["min_sep_fwhm"])


def _print_diag(tag, d):
    """Print the positive / negative peak counts of a diagnostic."""
    print(f"{tag} residual check: clean region {100 * d['area_frac']:.0f}% of px, noise scale "
          f"{d['scale']:.2f}, S/N offset {d['offset']:+.2f}")
    print(f"{tag}   peaks >= S/N (positive/negative):  " +
          "  ".join(f"{t:g}: {p}/{n}" for t, p, n in zip(d["t"], d["pos"], d["neg"])))


def _accept(d_old, d_new, t_eval):
    """Decide whether an extra round improved the residual (fewer excess, not more negative peaks)."""
    po, no = counts_at(d_old, t_eval)
    pn, nn = counts_at(d_new, t_eval)
    eo, en = po - no, pn - nn
    improved = (eo - en) >= max(2.0, 0.25 * max(eo, 0))
    neg_ok = nn <= no + max(2.0, 0.5 * no)
    return bool(improved and neg_ok), dict(t_eval=t_eval, excess_before=eo, excess_after=en,
                                           neg_before=no, neg_after=nn)


def _merge_results(a, b):
    """Join a result with its warm-started continuation (all snapshots, joined histories)."""
    out = dict(b)
    out["snaps"] = list(a.get("snaps", [])) + list(b.get("snaps", []))
    for k in ("resid_std_history", "n_peaks_history", "n_zero_amp_history"):
        out[k] = list(a[k]) + list(b[k])[1:]
    fs = dict(b.get("first_seen", {}))
    fs.update(a.get("first_seen", {}))
    out["first_seen"] = fs
    return out


def _adapt(orig, exposure, psf, bad_bitmask, outer_idx, cur, res, A):
    """Run warm-started extra inner rounds, relaxing settings while sources are still missing."""
    tag = f"[adaptive o{outer_idx + 1}]"
    target = float(A["target_snr"] if A["target_snr"] is not None else cur["min_snr"])
    d = _diagnose(exposure, psf, res, A, extra=[target, A["t_floor"]])
    if A["verbose"]:
        _print_diag(tag, d)
    log, rung, rounds = [], 0, 0
    missing, e, err = is_missing(d, target, A["z_sig"], A["min_excess"])
    while missing and rounds < A["max_rounds"] and rung < len(_RUNGS):
        t_new = pick_threshold(d, A["fdr"], A["t_floor"], target)
        new = None
        while new is None and rung < len(_RUNGS):
            new, thr, t_eval, desc = _RUNGS[rung](cur, t_new, target, A)
            rung += 1
        if new is None:
            print(f"{tag} no further relaxations available")
            break
        rounds += 1
        print(f"{tag} excess {e:+d} +/- {err:.1f} at >={target:g} sigma -> round {rounds}: {desc}")
        cont = dict(new)
        cont.update(n_iterations=A["extra_inner"], initial_peaks=res["peaks"],
                    threshold_start=thr, threshold_step=0.0, threshold_min=thr,
                    convergence_frac=0.0, viz=bool(A["viz_rounds"] and cur.get("viz")))
        res2 = orig(exposure, psf, bad_bitmask, **cont)
        if not res2.get("snaps"):
            print(f"{tag} round {rounds} produced nothing -- stopping")
            break
        d2 = _diagnose(exposure, psf, res2, A, extra=[target, t_eval, A["t_floor"]])
        ok, info = _accept(d, d2, t_eval)
        log.append(dict(round=rounds, change=desc, n_before=len(res["peaks"]),
                        n_after=len(res2["peaks"]), accepted=ok, **info))
        if A["verbose"]:
            _print_diag(tag, d2)
        print(f"{tag} round {rounds}: {len(res['peaks'])} -> {len(res2['peaks'])} sources; excess at "
              f">={t_eval:g} sigma {info['excess_before']:+d} -> {info['excess_after']:+d}, negatives "
              f"{info['neg_before']} -> {info['neg_after']}  => {'KEPT' if ok else 'REVERTED'}")
        if not ok:
            break
        res, d, cur, target = _merge_results(res, res2), d2, new, t_eval
        missing, e, err = is_missing(d, target, A["z_sig"], A["min_excess"])
    if not missing:
        print(f"{tag} residual consistent with noise at >={target:g} sigma "
              f"(excess {e:+d} +/- {err:.1f}) -- done after {rounds} extra round(s)")
    else:
        print(f"{tag} residual still shows excess {e:+d} +/- {err:.1f} at >={target:g} sigma "
              f"after {rounds} extra round(s); top remaining positive peaks (x, y, S/N, dist):")
        for x, y, s, dd in [c for c in d["candidates"] if c[2] >= target][:8]:
            print(f"{tag}     ({x:7.1f}, {y:7.1f})  S/N {s:5.1f}  {dd:5.1f}px from nearest fitted")
    res["adaptive"] = dict(rounds=log, final=d, target=target, missing=bool(missing),
                           excess=int(e), params={k: cur[k] for k in
                           ("min_snr", "near_bright_snr", "threshold_min", "reject_fwhm",
                            "merge_separation_sigma", "prune_sep_sigma", "shape_kwargs")
                           if k in cur})
    return res


def _make_wrapper(L):
    """Build the adaptive replacement of fitLib's inner loop."""
    orig = L._orig_run_inner_loop_verbose
    sig = inspect.signature(orig)
    defaults = {k: v.default for k, v in sig.parameters.items()
                if v.default is not inspect.Parameter.empty}

    def run_inner_loop_adaptive(exposure, psf, bad_bitmask, outer_idx=0, **kw):
        A = ADAPT
        cur = dict(defaults)
        cur.update(kw)
        cur["outer_idx"] = outer_idx
        if outer_idx == 0:
            _STATE["npeaks"] = []
        res = orig(exposure, psf, bad_bitmask, **cur)
        if A["enabled"] and res.get("snaps"):
            res = _adapt(orig, exposure, psf, bad_bitmask, outer_idx, cur, res, A)
        _STATE["npeaks"].append(len(res["peaks"]))
        return res

    return run_inner_loop_adaptive


def _wants_more_outer(result, outer_idx):
    """Decide whether the outer loop should run another pass.

    Called by `fitLib.run_cutout_pipeline` once the configured outer passes are used up: True if
    the residual still shows missing sources, the hard cap is not reached and the last pass
    added enough sources.
    """
    A = ADAPT
    ad = result.get("adaptive") if isinstance(result, dict) else None
    if not (A["enabled"] and ad and ad["missing"]):
        return False
    if outer_idx + 1 >= A["max_outer"]:
        return False
    h = _STATE["npeaks"]
    if len(h) >= 2 and h[-1] < h[-2] * (1.0 + A["min_outer_growth"]):
        print(f"[adaptive o{outer_idx + 1}] excess remains but the last outer pass added "
              f"<{100 * A['min_outer_growth']:g}% sources -- not extending the outer loop")
        return False
    print(f"[adaptive o{outer_idx + 1}] excess remains -> running outer pass {outer_idx + 2}")
    return True


# =====================================================================
# Public API
# =====================================================================
def install(**overrides):
    """Route fitLib's inner loop through the adaptive controller.

    Also enables the adaptive outer loop.  Requires the hooks added by
    ``patch_cutout_fit_lib.py``.

    Parameters
    ----------
    **overrides
        Settings that replace entries of ``ADAPT``.

    Raises
    ------
    ValueError
        Raised if an override is not an ``ADAPT`` option.
    RuntimeError
        Raised if fitLib has not been patched with the hooks.
    """
    from . import fitLib as L
    unknown = set(overrides) - set(ADAPT)
    if unknown:
        raise ValueError(f"unknown option(s) {sorted(unknown)}; valid: {sorted(ADAPT)}")
    ADAPT.update(overrides)
    if not hasattr(L, "_orig_run_inner_loop_verbose"):
        L._orig_run_inner_loop_verbose = L.run_inner_loop_verbose
    if "initial_peaks" not in inspect.signature(L._orig_run_inner_loop_verbose).parameters:
        raise RuntimeError("cutout_fit_lib.py has not been patched: run\n"
                           "    python patch_cutout_fit_lib.py cutout_fit_lib.py\n"
                           "then restart the kernel (or importlib.reload) and install again.")
    L.run_inner_loop_verbose = _make_wrapper(L)
    if hasattr(L, "_wants_more_outer"):
        L._wants_more_outer = _wants_more_outer
    else:
        print("note: cutout_fit_lib has no _wants_more_outer hook -- outer loop will not auto-extend")
    print("cutout_adaptive installed:", {k: ADAPT[k] for k in sorted(ADAPT)})


def uninstall():
    """Restore fitLib's original inner and outer loops."""
    from . import fitLib as L
    if hasattr(L, "_orig_run_inner_loop_verbose"):
        L.run_inner_loop_verbose = L._orig_run_inner_loop_verbose
    if hasattr(L, "_wants_more_outer"):
        L._wants_more_outer = lambda result, outer_idx: False
    print("cutout_adaptive uninstalled")


def report(run_or_history):
    """Print what the adaptive controller did in every outer pass of a finished run.

    Parameters
    ----------
    run_or_history : `CutoutRun` or `list` [`dict`]
        A finished run or its outer-loop history.
    """
    hist = getattr(run_or_history, "outer_history", run_or_history)
    for i, h in enumerate(hist, 1):
        ad = h["result"].get("adaptive")
        if not ad:
            print(f"outer {i}: no adaptive record")
            continue
        p, n = counts_at(ad["final"], ad["target"])
        print(f"outer {i}: {h['n_peaks']} sources, {len(ad['rounds'])} extra round(s); final peaks "
              f">={ad['target']:g} sigma: +{p}/-{n} ({'EXCESS REMAINS' if ad['missing'] else 'clean'})")
        for r in ad["rounds"]:
            print(f"    round {r['round']} [{'kept' if r['accepted'] else 'reverted'}] {r['change']}: "
                  f"{r['n_before']} -> {r['n_after']} sources, excess {r['excess_before']:+d} -> "
                  f"{r['excess_after']:+d}")

