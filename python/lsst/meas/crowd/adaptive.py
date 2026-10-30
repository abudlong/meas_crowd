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
crowding and correlated noise.  `residualExcess`, `isMissing`, `countsAt` and
`pickThreshold` are used by the pipeline's inner loop.

`install` adds the same check, with a ladder of relaxed thresholds, to
`lsst.meas.crowd.fitLib.runCutoutPipeline`.
"""
import copy
import inspect

import numpy as np
from scipy.ndimage import distance_transform_edt, maximum_filter
from scipy.signal import fftconvolve
from scipy.spatial import cKDTree

ADAPT = dict(
    enabled=True,
    maxRounds=3,          # max extra (warm-started) rounds per outer pass
    extraInner=3,         # inner iterations per extra round
    targetSnr=None,       # S/N down to which the residual should be clean (None: run's minSnr)
    tFloor=3.0,           # never relax thresholds below this S/N (set = targetSnr to forbid going deeper)
    fdr=0.25,              # tolerated  N_neg / N_pos  when choosing the lowest threshold
    zSig=2.0,             # excess must exceed zSig * sqrt(N_pos + N_neg) ...
    minExcess=3,          # ... and be at least this many sources
    edgeFwhm=2.0,         # clean region: farther than this from the cutout edge [FWHM]
    haloFwhm=2.5,         # clean region: farther than this from bright sources [FWHM] (grows with flux)
    brightSigma=100.0,    # 'bright' = amp*psfPeak/noise above this
    minSepFwhm=0.75,     # peaks closer than this to a fitted source are its own mismatch residual
    maxOuter=6,           # hard cap on outer passes when auto-extending the outer loop
    minOuterGrowth=0.02, # extend the outer loop only if the last pass added >= this fraction of sources
    vizRounds=False,      # draw the per-iteration figures during the extra rounds
    verbose=True,
)
_STATE = dict(npeaks=[])

_DEFAULT_SHAPE = dict(minSnrLocal=5.0, chi2Max=2.0, widthRange=(0.75, 1.35),
                      dchi2Width=9.0, trustSnr=30.0, minContrast=20.0,
                      weakSnr=10.0, dchi2WidthWeak=4.0)


# =====================================================================
# Diagnostic (numpy / scipy only)
# =====================================================================
def _odd(x):
    """Return the nearest odd integer >= 3."""
    return max(3, int(round(x)) | 1)


def _matchedSnr(resid, var, weightOk, fwhm, calib):
    """Return the Gaussian matched-filter S/N of the residual.

    Pixels with ``weightOk`` False carry zero weight.  The S/N is calibrated empirically (median
    offset, MAD scale >= 1) on the ``calib`` pixels so correlated noise does not inflate it.
    Returns (snr, den, scale, offset, ok).
    """
    sig = fwhm / 2.355
    h = int(np.ceil(3 * sig))
    u = np.arange(-h, h + 1)
    k = np.exp(-0.5 * (u[None, :] ** 2 + u[:, None] ** 2) / sig ** 2)
    good = weightOk & np.isfinite(resid) & np.isfinite(var) & (var > 0)
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


def residualExcess(resid, var, good, peaksXy, amps, x0, y0, fwhm, psfPeak, noise,
                    thresholds=(3.0, 3.5, 4.0, 4.5, 5.0, 6.0), edgeFwhm=2.0,
                    brightSigma=100.0, haloFwhm=2.5, minSepFwhm=0.75):
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
    peaksXy : array-like
        Fitted source positions (parent pixel coordinates), shape (N, 2).
    amps : array-like
        Fitted source amplitudes.
    x0, y0 : `int`
        Parent pixel coordinates of ``resid[0, 0]``.
    fwhm : `float`
        PSF FWHM [px].
    psfPeak : `float`
        Peak value of the unit-flux PSF.
    noise : `float`
        Typical pixel noise.
    thresholds : `tuple` [`float`], optional
        S/N thresholds at which to count.
    edgeFwhm : `float`, optional
        Ignore peaks within this many FWHM of the edge.
    brightSigma : `float`, optional
        A source is bright if ``amp * psfPeak / noise`` exceeds this.
    haloFwhm : `float`, optional
        Ignore peaks within this many FWHM of a bright source (grows with flux).
    minSepFwhm : `float`, optional
        Peaks closer than this to a fitted source are its own mismatch residual.

    Returns
    -------
    result : `dict`
        ``t`` (thresholds, ascending), ``pos`` / ``neg`` (number of peaks with |S/N| >= t),
        ``scale``, ``offset`` (S/N calibration), ``areaFrac`` (fraction of the image in the
        clean region) and ``candidates`` ((x, y, snr, distance to the nearest fitted source) of
        the positive peaks).
    """
    ny, nx = resid.shape
    iy, ix = np.mgrid[0:ny, 0:nx]
    P = np.asarray(peaksXy, float).reshape(-1, 2) - np.array([x0, y0])
    amps = np.asarray(amps, float)
    ts = np.array(sorted({round(float(t), 3) for t in thresholds}))

    region = good & (np.minimum.reduce([ix, nx - 1 - ix, iy, ny - 1 - iy]) >= edgeFwhm * fwhm)
    if len(P):
        sigBr = amps * psfPeak / noise
        for j in np.flatnonzero(sigBr > brightSigma):
            r = haloFwhm * fwhm * (1.0 + 0.5 * max(0.0, np.log10(sigBr[j] / brightSigma)))
            region &= (ix - P[j, 0]) ** 2 + (iy - P[j, 1]) ** 2 > r * r
    if not good.all():
        region &= distance_transform_edt(good) >= 0.5 * fwhm
    snr, den, scale, offset, ok = _matchedSnr(resid, var, good, fwhm, region)
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
            k = nn >= minSepFwhm * fwhm
            xx, yy, v, nn = xx[k], yy[k], v[k], nn[k]
        return xx, yy, v, nn

    xp, yp, vp, dp = extrema(+1)
    _, _, vn, _ = extrema(-1)
    order = np.argsort(-vp)
    return dict(t=ts,
                pos=np.array([int((vp >= t).sum()) for t in ts]),
                neg=np.array([int((vn >= t).sum()) for t in ts]),
                scale=scale, offset=offset, areaFrac=float(region.mean()),
                candidates=[(float(xp[i] + x0), float(yp[i] + y0), float(vp[i]), float(dp[i]))
                            for i in order])


def _idx(d, t):
    """Return the index of the threshold in ``d['t']`` closest to ``t``."""
    return int(np.argmin(np.abs(d["t"] - t)))


def countsAt(d, t):
    """Return the positive and negative peak counts at an S/N threshold.

    Parameters
    ----------
    d : `dict`
        Output of `residualExcess`.
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


def isMissing(d, t, zSig=2.0, minExcess=3):
    """Decide whether the residual still contains missing sources at an S/N threshold.

    Parameters
    ----------
    d : `dict`
        Output of `residualExcess`.
    t : `float`
        S/N threshold.
    zSig : `float`, optional
        The excess must exceed ``zSig * sqrt(N_pos + N_neg)`` ...
    minExcess : `int`, optional
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
    p, n = countsAt(d, t)
    e, err = p - n, float(np.sqrt(p + n))
    return bool(e >= max(minExcess, zSig * err)), e, err


def pickThreshold(d, fdr, tFloor, tCap):
    """Return the lowest S/N threshold the field can support.

    The lowest S/N in [``tFloor``, ``tCap``] whose estimated false fraction
    (N_neg + 1) / N_pos is at most ``fdr`` (with N_pos >= 3).

    Parameters
    ----------
    d : `dict`
        Output of `residualExcess`.
    fdr : `float`
        Tolerated false fraction.
    tFloor : `float`
        Lowest threshold allowed.
    tCap : `float`
        Highest threshold; returned if no lower one qualifies.

    Returns
    -------
    threshold : `float`
        The chosen S/N threshold.
    """
    for t, p, n in zip(d["t"], d["pos"], d["neg"]):
        if t < tFloor - 1e-9 or t > tCap + 1e-9:
            continue
        if p >= 3 and (n + 1.0) / p <= fdr:
            return float(t)
    return float(tCap)


# =====================================================================
# Relaxation ladder.  Each rung: (cur_kwargs, tNew, target, A) ->
#   (new_kwargs, detection_threshold, eval_snr, description)  or  (None, ...) if a no-op
# =====================================================================
def _copy(cur):
    """Copy the inner-loop keyword arguments (with their shapeKwargs dict)."""
    new = dict(cur)
    if isinstance(cur.get("shapeKwargs"), dict):
        new["shapeKwargs"] = dict(cur["shapeKwargs"])
    return new


def _shape(cur):
    """Return the shape-cut keyword arguments with the defaults filled in."""
    sk = dict(_DEFAULT_SHAPE)
    sk.update(cur.get("shapeKwargs") or {})
    return sk


def _rungThresholds(cur, tNew, target, A):
    """Relaxation rung: lower the S/N cuts and the detection threshold."""
    old, thr0 = cur["minSnr"], cur["thresholdMin"]
    if tNew >= old - 1e-9 and tNew >= thr0 - 1e-9:
        return None, None, None, None
    new, ch = _copy(cur), []
    if tNew < old - 1e-9:
        new["minSnr"] = tNew
        new["nearBrightSnr"] = max(tNew, cur["nearBrightSnr"] * tNew / old)
        sk = _shape(cur)
        sk["minSnrLocal"] = min(sk["minSnrLocal"], tNew)
        new["shapeKwargs"] = sk
        ch.append(f"minSnr {old:g}->{tNew:g}")
    thr = min(thr0, tNew)
    ch.append(f"detect to {thr:g} sigma")
    return new, thr, min(target, tNew), "lower thresholds (" + ", ".join(ch) + ")"


def _rungGates(cur, tNew, target, A):
    """Relaxation rung: loosen the shape-cut gates and the bright-star rejections."""
    new, sk = _copy(cur), _shape(cur)
    sk["chi2Max"] *= 1.5
    sk["dchi2Width"] *= 1.5
    sk["dchi2WidthWeak"] *= 1.5
    if sk.get("minContrast"):
        sk["minContrast"] = min(sk["minContrast"], 1.0)
    new["shapeKwargs"] = sk
    new["rejectFwhm"] = cur["rejectFwhm"] * 0.7
    new["nearBrightSnr"] = max(cur["minSnr"], 0.75 * cur["nearBrightSnr"])
    return (new, cur["thresholdMin"], target,
            f"relax gates (chi2Max {sk['chi2Max']:.2f}, dchi2Width {sk['dchi2Width']:.1f}, "
            f"rejectFwhm {new['rejectFwhm']:.2f}, nearBrightSnr {new['nearBrightSnr']:.1f})")


def _rungSeparations(cur, tNew, target, A):
    """Relaxation rung: allow closer pairs (merge / prune / reject separations)."""
    new = _copy(cur)
    new["mergeSeparationSigma"] = max(0.8, cur["mergeSeparationSigma"] * 0.8)
    new["pruneSepSigma"] = max(0.8, cur["pruneSepSigma"] * 0.8)
    new["rejectFwhm"] = cur["rejectFwhm"] * 0.7
    return (new, cur["thresholdMin"], target,
            f"relax separations (merge {new['mergeSeparationSigma']:.2f} sigma, "
            f"prune {new['pruneSepSigma']:.2f} sigma, rejectFwhm {new['rejectFwhm']:.2f})")


_RUNGS = [_rungThresholds, _rungGates, _rungSeparations]


# =====================================================================
# Controller
# =====================================================================
def _diagnose(exposure, psf, res, A, extra=()):
    """Run `residualExcess` on the last snapshot of an inner-loop result."""
    snap = res["snaps"][-1]
    fit = snap["fit"]
    bbox = exposure.getBBox()
    c = bbox.getCenter()
    shp = psf.computeShape(c)
    fwhm = 2.355 * np.sqrt(0.5 * (shp.getIxx() + shp.getIyy()))
    psfPeak = psf.computeKernelImage(c).array.max()
    var, good = fit["variance"], fit["good"]
    noise = float(np.sqrt(np.nanmedian(var[good])))
    resid = exposure.image.array - fit["modelFull"]
    ts = [3.0, 3.5, 4.0, 4.5, 5.0, 6.0] + [float(e) for e in extra]
    return residualExcess(resid, var, good, snap["mergedArr"], fit["amps"],
                           bbox.getMinX(), bbox.getMinY(), fwhm, psfPeak, noise,
                           thresholds=ts, edgeFwhm=A["edgeFwhm"],
                           brightSigma=A["brightSigma"], haloFwhm=A["haloFwhm"],
                           minSepFwhm=A["minSepFwhm"])


def _printDiag(tag, d):
    """Print the positive / negative peak counts of a diagnostic."""
    print(f"{tag} residual check: clean region {100 * d['areaFrac']:.0f}% of px, noise scale "
          f"{d['scale']:.2f}, S/N offset {d['offset']:+.2f}")
    print(f"{tag}   peaks >= S/N (positive/negative):  " +
          "  ".join(f"{t:g}: {p}/{n}" for t, p, n in zip(d["t"], d["pos"], d["neg"])))


def _accept(dOld, dNew, tEval):
    """Decide whether an extra round improved the residual (fewer excess, not more negative peaks)."""
    po, no = countsAt(dOld, tEval)
    pn, nn = countsAt(dNew, tEval)
    eo, en = po - no, pn - nn
    improved = (eo - en) >= max(2.0, 0.25 * max(eo, 0))
    negOk = nn <= no + max(2.0, 0.5 * no)
    return bool(improved and negOk), dict(tEval=tEval, excessBefore=eo, excessAfter=en,
                                           negBefore=no, negAfter=nn)


def _mergeResults(a, b):
    """Join a result with its warm-started continuation (all snapshots, joined histories)."""
    out = dict(b)
    out["snaps"] = list(a.get("snaps", [])) + list(b.get("snaps", []))
    for k in ("residStdHistory", "nPeaksHistory", "nZeroAmpHistory"):
        out[k] = list(a[k]) + list(b[k])[1:]
    fs = dict(b.get("firstSeen", {}))
    fs.update(a.get("firstSeen", {}))
    out["firstSeen"] = fs
    return out


def _adapt(orig, exposure, psf, badBitmask, outerIdx, cur, res, A):
    """Run warm-started extra inner rounds, relaxing settings while sources are still missing."""
    tag = f"[adaptive o{outerIdx + 1}]"
    target = float(A["targetSnr"] if A["targetSnr"] is not None else cur["minSnr"])
    d = _diagnose(exposure, psf, res, A, extra=[target, A["tFloor"]])
    if A["verbose"]:
        _printDiag(tag, d)
    log, rung, rounds = [], 0, 0
    missing, e, err = isMissing(d, target, A["zSig"], A["minExcess"])
    while missing and rounds < A["maxRounds"] and rung < len(_RUNGS):
        tNew = pickThreshold(d, A["fdr"], A["tFloor"], target)
        new = None
        while new is None and rung < len(_RUNGS):
            new, thr, tEval, desc = _RUNGS[rung](cur, tNew, target, A)
            rung += 1
        if new is None:
            print(f"{tag} no further relaxations available")
            break
        rounds += 1
        print(f"{tag} excess {e:+d} +/- {err:.1f} at >={target:g} sigma -> round {rounds}: {desc}")
        cont = dict(new)
        cont.update(nIterations=A["extraInner"], initialPeaks=res["peaks"],
                    thresholdStart=thr, thresholdStep=0.0, thresholdMin=thr,
                    convergenceFrac=0.0, viz=bool(A["vizRounds"] and cur.get("viz")))
        res2 = orig(exposure, psf, badBitmask, **cont)
        if not res2.get("snaps"):
            print(f"{tag} round {rounds} produced nothing -- stopping")
            break
        d2 = _diagnose(exposure, psf, res2, A, extra=[target, tEval, A["tFloor"]])
        ok, info = _accept(d, d2, tEval)
        log.append(dict(round=rounds, change=desc, nBefore=len(res["peaks"]),
                        nAfter=len(res2["peaks"]), accepted=ok, **info))
        if A["verbose"]:
            _printDiag(tag, d2)
        print(f"{tag} round {rounds}: {len(res['peaks'])} -> {len(res2['peaks'])} sources; excess at "
              f">={tEval:g} sigma {info['excessBefore']:+d} -> {info['excessAfter']:+d}, negatives "
              f"{info['negBefore']} -> {info['negAfter']}  => {'KEPT' if ok else 'REVERTED'}")
        if not ok:
            break
        res, d, cur, target = _mergeResults(res, res2), d2, new, tEval
        missing, e, err = isMissing(d, target, A["zSig"], A["minExcess"])
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
                           ("minSnr", "nearBrightSnr", "thresholdMin", "rejectFwhm",
                            "mergeSeparationSigma", "pruneSepSigma", "shapeKwargs")
                           if k in cur})
    return res


def _makeWrapper(L):
    """Build the adaptive replacement of fitLib's inner loop."""
    orig = L._origRunInnerLoopVerbose
    sig = inspect.signature(orig)
    defaults = {k: v.default for k, v in sig.parameters.items()
                if v.default is not inspect.Parameter.empty}

    def runInnerLoopAdaptive(exposure, psf, badBitmask, outerIdx=0, **kw):
        A = ADAPT
        cur = dict(defaults)
        cur.update(kw)
        cur["outerIdx"] = outerIdx
        if outerIdx == 0:
            _STATE["npeaks"] = []
        res = orig(exposure, psf, badBitmask, **cur)
        if A["enabled"] and res.get("snaps"):
            res = _adapt(orig, exposure, psf, badBitmask, outerIdx, cur, res, A)
        _STATE["npeaks"].append(len(res["peaks"]))
        return res

    return runInnerLoopAdaptive


def _wantsMoreOuter(result, outerIdx):
    """Decide whether the outer loop should run another pass.

    Called by `fitLib.runCutoutPipeline` once the configured outer passes are used up: True if
    the residual still shows missing sources, the hard cap is not reached and the last pass
    added enough sources.
    """
    A = ADAPT
    ad = result.get("adaptive") if isinstance(result, dict) else None
    if not (A["enabled"] and ad and ad["missing"]):
        return False
    if outerIdx + 1 >= A["maxOuter"]:
        return False
    h = _STATE["npeaks"]
    if len(h) >= 2 and h[-1] < h[-2] * (1.0 + A["minOuterGrowth"]):
        print(f"[adaptive o{outerIdx + 1}] excess remains but the last outer pass added "
              f"<{100 * A['minOuterGrowth']:g}% sources -- not extending the outer loop")
        return False
    print(f"[adaptive o{outerIdx + 1}] excess remains -> running outer pass {outerIdx + 2}")
    return True


# =====================================================================
# Public API
# =====================================================================
def install(**overrides):
    """Route fitLib's inner loop through the adaptive controller.

    Also enables the adaptive outer loop (through fitLib's ``_wantsMoreOuter`` hook).

    Parameters
    ----------
    **overrides
        Settings that replace entries of ``ADAPT``.

    Raises
    ------
    ValueError
        Raised if an override is not an ``ADAPT`` option.
    RuntimeError
        Raised if fitLib's inner loop has no ``initialPeaks`` argument.
    """
    from . import fitLib as L
    unknown = set(overrides) - set(ADAPT)
    if unknown:
        raise ValueError(f"unknown option(s) {sorted(unknown)}; valid: {sorted(ADAPT)}")
    ADAPT.update(overrides)
    if not hasattr(L, "_origRunInnerLoopVerbose"):
        L._origRunInnerLoopVerbose = L.runInnerLoopVerbose
    if "initialPeaks" not in inspect.signature(L._origRunInnerLoopVerbose).parameters:
        raise RuntimeError("fitLib.runInnerLoopVerbose has no initialPeaks argument; "
                           "the adaptive controller needs it")
    L.runInnerLoopVerbose = _makeWrapper(L)
    if hasattr(L, "_wantsMoreOuter"):
        L._wantsMoreOuter = _wantsMoreOuter
    else:
        print("note: fitLib has no _wantsMoreOuter hook -- outer loop will not auto-extend")
    print("adaptive controller installed:", {k: ADAPT[k] for k in sorted(ADAPT)})


def uninstall():
    """Restore fitLib's original inner and outer loops."""
    from . import fitLib as L
    if hasattr(L, "_origRunInnerLoopVerbose"):
        L.runInnerLoopVerbose = L._origRunInnerLoopVerbose
    if hasattr(L, "_wantsMoreOuter"):
        L._wantsMoreOuter = lambda result, outerIdx: False
    print("adaptive controller uninstalled")


def report(runOrHistory):
    """Print what the adaptive controller did in every outer pass of a finished run.

    Parameters
    ----------
    runOrHistory : `CutoutRun` or `list` [`dict`]
        A finished run or its outer-loop history.
    """
    hist = getattr(runOrHistory, "outerHistory", runOrHistory)
    for i, h in enumerate(hist, 1):
        ad = h["result"].get("adaptive")
        if not ad:
            print(f"outer {i}: no adaptive record")
            continue
        p, n = countsAt(ad["final"], ad["target"])
        print(f"outer {i}: {h['nPeaks']} sources, {len(ad['rounds'])} extra round(s); final peaks "
              f">={ad['target']:g} sigma: +{p}/-{n} ({'EXCESS REMAINS' if ad['missing'] else 'clean'})")
        for r in ad["rounds"]:
            print(f"    round {r['round']} [{'kept' if r['accepted'] else 'reverted'}] {r['change']}: "
                  f"{r['nBefore']} -> {r['nAfter']} sources, excess {r['excessBefore']:+d} -> "
                  f"{r['excessAfter']:+d}")



# old snake_case names still work (with a FutureWarning)
from ._compat import moduleGetattr  # noqa: E402
__getattr__ = moduleGetattr(globals())
