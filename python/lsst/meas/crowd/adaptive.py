# This file is part of meas_crowd.
#
# Developed for the LSST Data Management System.
# See the COPYRIGHT file at the top-level directory of this distribution.
#
# Vendored from the original notebook module cutout_adaptive.py (behaviour unchanged; only the
# imports of the sibling modules were made package-relative).

"""
cutout_adaptive.py
==================
Residual-driven adaptive control for the cutout pipeline.

IDEA
----
After every inner loop the pipeline looks at its own residual.  A matched-filter S/N map
of the residual is searched for local maxima (candidate missing sources) and local MINIMA
(pure noise, since real sources are positive).  In the *clean region* (away from the
cutout edge, bright-star halos, bad pixels, and fitted sources) the two counts should agree
if nothing is missing:

        excess(t) = N_pos(>= t) - N_neg(>= t)        error ~ sqrt(N_pos + N_neg)

The negative count is measured IN THE SAME FIELD, so the test adapts by itself to dense /
sparse regions and to correlated noise: in a crowded field the noise scale and the negative
count rise, and a significant excess needs more positive peaks.

If the excess is significant the inner loop is CONTINUED (warm start from the current
sources, a few extra iterations) with a ladder of relaxations, one rung per round:

    rung 0  lower thresholds  -> min_snr / near_bright_snr / shape min_snr_local / detection
                                 threshold go down to the lowest S/N at which the estimated
                                 false fraction  N_neg/N_pos  is still <= `fdr` (>= t_floor)
    rung 1  relax gates       -> shape-cut chi2 / width limits x1.5, min_contrast <= 1,
                                 reject_fwhm x0.7, near-bright S/N x0.75
    rung 2  relax separations -> merge / prune separations x0.8 (blends), reject_fwhm x0.7

A round is KEPT only if the excess at the evaluation S/N drops clearly and the negative
count does not grow (no noise fitting); otherwise it is reverted and the controller stops.
Optionally the OUTER loop is extended too (PSF / background refit + new inner loop) while
excess remains and the last outer pass still added sources.

SETUP (once)
------------
    python patch_cutout_fit_lib.py cutout_fit_lib.py      # 3 tiny hooks; makes a .bak
    # restart the kernel, then:
    from . import adaptive as ca
    ca.install(target_snr=None, max_rounds=3, max_outer=6)   # BEFORE cn.run_session(...)
    runs = cn.run_session(sess, cfg, ...)                    # unchanged
    ca.report(runs[name])                                    # what the controller did

`ca.install(enabled=False)` or `ca.uninstall()` gives back the original behaviour.

Everything in the diagnostic is plain numpy/scipy (importable and testable without LSST).
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
    return max(3, int(round(x)) | 1)


def _matched_snr(resid, var, weight_ok, fwhm, calib):
    """Gaussian matched-filter S/N of the residual.  Pixels with weight_ok False carry zero
    weight.  Calibrated empirically (median offset / MAD scale >= 1) on `calib` pixels so
    correlated noise does not inflate it.  Returns (snr, den, scale, offset, ok)."""
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
    """Positive / negative local-maximum counts of the matched-filter S/N of the residual,
    in the clean region.  Returns dict(t, pos, neg, scale, offset, area_frac, candidates)
    with t ascending and pos[i] / neg[i] = number of peaks with |S/N| >= t[i];
    candidates = [(x, y, snr, nearest_fitted_px), ...] for positive peaks >= min(t)."""
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
    return int(np.argmin(np.abs(d["t"] - t)))


def counts_at(d, t):
    i = _idx(d, t)
    return int(d["pos"][i]), int(d["neg"][i])


def is_missing(d, t, z_sig=2.0, min_excess=3):
    """(missing?, excess, error) at S/N threshold t."""
    p, n = counts_at(d, t)
    e, err = p - n, float(np.sqrt(p + n))
    return bool(e >= max(min_excess, z_sig * err)), e, err


def pick_threshold(d, fdr, t_floor, t_cap):
    """Lowest S/N in [t_floor, t_cap] whose estimated false fraction (N_neg + 1)/N_pos <= fdr
    (and N_pos >= 3); t_cap if there is none (the field cannot support going deeper)."""
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
    new = dict(cur)
    if isinstance(cur.get("shape_kwargs"), dict):
        new["shape_kwargs"] = dict(cur["shape_kwargs"])
    return new


def _shape(cur):
    sk = dict(_DEFAULT_SHAPE)
    sk.update(cur.get("shape_kwargs") or {})
    return sk


def _rung_thresholds(cur, t_new, target, A):
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
    print(f"{tag} residual check: clean region {100 * d['area_frac']:.0f}% of px, noise scale "
          f"{d['scale']:.2f}, S/N offset {d['offset']:+.2f}")
    print(f"{tag}   peaks >= S/N (positive/negative):  " +
          "  ".join(f"{t:g}: {p}/{n}" for t, p, n in zip(d["t"], d["pos"], d["neg"])))


def _accept(d_old, d_new, t_eval):
    po, no = counts_at(d_old, t_eval)
    pn, nn = counts_at(d_new, t_eval)
    eo, en = po - no, pn - nn
    improved = (eo - en) >= max(2.0, 0.25 * max(eo, 0))
    neg_ok = nn <= no + max(2.0, 0.5 * no)
    return bool(improved and neg_ok), dict(t_eval=t_eval, excess_before=eo, excess_after=en,
                                           neg_before=no, neg_after=nn)


def _merge_results(a, b):
    """b = warm-started continuation of a: keep all snapshots and the joined histories."""
    out = dict(b)
    out["snaps"] = list(a.get("snaps", [])) + list(b.get("snaps", []))
    for k in ("resid_std_history", "n_peaks_history", "n_zero_amp_history"):
        out[k] = list(a[k]) + list(b[k])[1:]
    fs = dict(b.get("first_seen", {}))
    fs.update(a.get("first_seen", {}))
    out["first_seen"] = fs
    return out


def _adapt(orig, exposure, psf, bad_bitmask, outer_idx, cur, res, A):
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
    """Hook called by run_cutout_pipeline once the configured outer passes are used up."""
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
    """Route cutout_fit_lib's inner loop through the adaptive controller (and enable the
    adaptive outer loop).  Requires the hooks from patch_cutout_fit_lib.py."""
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
    from . import fitLib as L
    if hasattr(L, "_orig_run_inner_loop_verbose"):
        L.run_inner_loop_verbose = L._orig_run_inner_loop_verbose
    if hasattr(L, "_wants_more_outer"):
        L._wants_more_outer = lambda result, outer_idx: False
    print("cutout_adaptive uninstalled")


def report(run_or_history):
    """Print what the controller did in every outer pass of a finished run."""
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

