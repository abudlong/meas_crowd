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

"""Adaptive, iterative joint PSF-fitting pipeline starting from Gaia positions.

Inner and outer loops, saturated and bright-star handling, PSF refinement, extended sources
and residual diagnostics.
"""

import math
import os
import time
import traceback
import types
from collections import OrderedDict, defaultdict
import numpy as np
from scipy.linalg import cho_factor, solve_triangular
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from astropy.table import Table
from ._compat import camelKeys
from . import fitLib as L
from . import adaptive as ca
from .gaiaCatalog import (
    attachGaia,
    selectGaiaSources,
)


# Options of this driver on top of the runCutoutPipeline config keys.
PIPELINE_DEFAULTS = dict(
    # ---- adaptive inner loop
    maxInner=30,              # safety cap on inner iterations per outer pass
    gaiaThresholdStart=5.0,  # first detection threshold when the fit starts from Gaia
    thresholdStart=None,      # without Gaia (None: innerKwargs['thresholdStart'])
    thresholdStep=None,       # None: innerKwargs['thresholdStep']
    thresholdMin=None,        # None: innerKwargs['thresholdMin']
    growFrac=0.02,            # "still finding sources": added >= growFrac * N ...
    growMin=3,                # ... or >= growMin sources -> hold the threshold
    patience=2,                # stop after this many iterations without progress
    detectGain=1.0,           # detection runs on data - gain*model (1 = full joint model)
    targetSnr=None,           # residual must be clean down to this S/N (None: minSnr)
    tFloor=3.0,               # never detect below this S/N
    fdr=0.25,                  # tolerated N_neg/N_pos when going deeper
    zSig=2.0,                 # significant excess: >= zSig*sqrt(N_pos+N_neg) ...
    minExcess=3,              # ... and >= this many sources
    edgeFwhm=2.0, haloFwhm=2.5, minSepFwhm=0.75,   # residual-check clean region
    centroidRoundsAfterFirst=1,
    hopelessFrac=0.5,         # eliminateWeak: drop S/N < hopelessFrac*minSnr at once
    elimDropFrac=0.25,
    # ---- adaptive outer loop
    minOuter=2,               # at least this many passes (pass 2 refines PSF + shape cut)
    maxOuter=6,               # safety cap
    outerGrowthTol=0.01,     # source count changed by < 1% ...
    psfTol=0.01,              # ... and PSF kernel changed by < 1% (L2)  -> converged
    # ---- Gaia starting positions
    protectGaia=True,
    gaiaPadFwhm=1.5,         # use Gaia stars up to this far OUTSIDE the cutout
    gaiaDropOffDetector=True,  # ... except beyond the DETECTOR edge (the nearest cutout pixel is
    #   EDGE / NO_DATA): no data constrain their flux, which the fit would push to absurd values
    gaiaMaxDriftPx=1.5,     # max distance a Gaia source may move from its Gaia position
    gaiaMinSepPx=0.5,       # Gaia duplicates closer than this are merged
    fitGaiaOffset=True,      # remove a global WCS-vs-data offset first
    subtractInitialBg=True,  # subtract session.bgImage in pass 1
    keepAllFits=False,       # keep full fits (all PSF stamps) in every snapshot
    # ---- saturated / bright stars
    satGrowPx=2,             # also exclude pixels within this many px of SAT (non-linear,
                               # charge spill); 0 = off
    brightCorr=True,          # per-bright-star shift + width + flux correction terms
    brightCorrSigma=100.0,   # ... for stars with amp*psfPeak/noise above this (or SAT)
    brightCorrRadiusFwhm=4.0,
    brightCorrRadial=True,   # ... plus radial-profile terms (core/wing shape mismatch)
    brightCoreRingsPx=(0.0, 1.0, 2.0, 3.0, 4.0),  # bright (unsaturated) stars also get NARROW
    brightCoreRingSigma=0.6,  #   radial terms (a central spike + thin rings at these radii [px],
    #   Gaussian width this many px): the core of bright stars is sharper
    #   than the PSF on a 1-2 px scale (spike / dip / ring residual)
    brightCoreModes=2,       # ... and each narrow term may be lopsided / elongated (angular modes
    #   m <= this): the sharp core can sit a fraction of a pixel off the
    #   rest of the star (core dipoles)
    brightAngularModes=4,    # unsaturated bright stars: rings (<= 4 FWHM) may vary with
    #   angle (cos/sin m*theta, m = 1..this) -- removes the lopsided
    #   'dipole' residuals a shift term cannot.  0 = off
    satHaloModes=4,          # saturated stars: the halo rings may also vary with angle
    #   (cos/sin m*theta, m = 1..this) -- spikes, bleed charge and
    #   lopsided halos are not circular.  0 = radial only
    # ---- saturated-star zones
    satZoneFwhm=2.0,         # zone = pixels within satZoneFwhm*FWHM + satZoneScale*(half-
    satZoneScale=1.0,        #   width of the saturated region) of the SAT mask -- a band that
    #   follows the mask's shape (bleed trails stay narrow)
    satZoneSnr=10.0,         # non-Gaia sources inside a zone need S/N >= this (instead of
    #   minSnr): real stars survive, halo/spike speckle does not.
    #   The residual check ignores the zones. satZoneFwhm=0 -> off
    satZoneMinPx=4,         # ignore saturated blobs smaller than this (hot pixels, CRs)
    psfSpatial=False,         # let the PSF vary LINEARLY across the cutout (kernel gradient in x
    #   and y, measured from the clean stars); kept only if it lowers
    #   chi2 on those stars by >= psfSpatialMinGain
    psfSpatialMinGain=0.02,
    psfSpatialMinStars=10,
    psfCvMinGain=0.005,     # PSF updates are CROSS-VALIDATED: built from half the clean stars,
    #   judged on the other half (both ways); kept only if the held-out
    #   chi2 drops by >= this fraction (stops the PSF fitting noise)
    psfSmoothPx=0.6,         # Gaussian smoothing [px] of every PSF correction (a real PSF
    #   change is smooth on the PSF scale; pixel noise is not)
    psfSoftSigma=2.0,        # corrections below this many sigma per pixel are ignored
    psfStackN=25,            # PSF update between outer passes: up to this many CLEAN stars
    psfStackMinSnr=50.0,    #   (peak S/N >= this, no saturated pixel within 3 FWHM, no
    psfStackIsoFwhm=2.5,    #   brighter-than-10% neighbour within psfStackIsoFwhm*FWHM,
    #   not extended).  Saturated stars used to dominate the stack.
    maskPadFwhm=1.5,         # residual check ignores pixels within this many FWHM of bad
    #   pixels (flux of sources whose centre is masked)
    satRefineWings=True,     # saturated stars without a Gaia position: refine the position
    #   from the unsaturated wings (the core cannot be used)
    wingKeepAngular=False,   # wing extensions keep the kernel's angular pattern (spikes) beyond
    #   the blend radius.  OFF: the visit kernel's outer part is mostly
    #   noise and the pattern imprinted noisy rings on bright stars
    psfCoreMedian=False,     # (experimental, off: did not help on real data) after each outer pass: correct the PSF CORE with the MEDIAN residual
    psfCoreMinSnr=100.0,    #   of the bright clean stars (peak S/N >= this; what
    psfCoreRadiusFwhm=2.5,  #   showBrightStarStack shows), within this many FWHM, tapered,
    psfCoreCycles=2,         #   cross-validated on two halves of the stars
    bfTerm=False,             # (experimental, off: did not help on real data) BRIGHTNESS-DEPENDENT PSF (brighter-fatter): every star's PSF is
    bfMinSnr=30.0,           #   K + (log10 flux - L0) * B(r), with the radial profile B measured by
    bfRadiusFwhm=3.0,        #   regressing isolated stars' residual profiles on log flux (S/N >=
    #   bfMinSnr, within bfRadiusFwhm*FWHM); flux-neutral, smooth,
    #   cross-validated on two halves of the stars
    modelErrMeasure=True,    # (A) after each outer pass, MEASURE the PSF model error eps: the
    modelErrMinSnr=50.0,    #   residual variance around stars (S/N >= this) beyond the pixel
    #   variance, as (eps * model)^2.  plotFinalComparison then shows
    #   residual / sqrt(var + (eps*model)^2) (sigma="model", default)
    fitSysFromData=False,   # ... and also use eps as sysFrac in the FIT weights / residual check
    #   (default False: they keep sysFrac=0.03, which is more tolerant)
    nearSatFrac=0.0,         # (B, off by default; e.g. 0.7 to use) pixels above this fraction of the saturation level (not flagged
    nearSatLevel=None,       #   SAT) are excluded like saturated ones: those stars are fitted
    nearSatGrowPx=1,        #   from their wings.  nearSatLevel: saturation level in image
    #   units (None = estimated from the SAT pixels).  0/None frac = off
    psfRegularize=True,       # after EVERY PSF update: beyond psfRegCoreFwhm*FWHM the kernel is
    psfRegCoreFwhm=1.5,     #   rebuilt as a smooth, monotone radial profile (ring medians) times
    psfRegSigma=1.0,         #   a smoothed angular pattern -- no holes or isolated hot pixels
    psfExtBlend=(12.0, 16.0),  # initial 41->81 px wing extension: the kernel is replaced by a
                               #   fitted Moffat between these radii [px]; e.g. (17, 19.5) keeps
                               #   the PSF model's own profile almost to the kernel edge
    wingBlendFwhm=(2.0, 3.0),  # measured wing profile replaces the kernel between these radii
                               #   [FWHM]; try (3.5, 5.0) to keep more of the PSF model's own wing
    wingsFromStars=True,     # measure the PSF's outer wings from the bright + saturated stars'
    wingMinSnr=300.0,        #   radial profiles (unsaturated stars with peak S/N >= this set
                               #   the shape; saturated stars extend it outwards) instead of
                               #   extrapolating the kernel's noisy 6-14 px profile
    bgMode="plane",           # background refit between passes: "plane" (tilted plane, sigma-clipped),
                               #   "constant", "spline" (SubtractBackgroundTask surface on bgBinSize
                               #   bins -- can swing by 100s of counts in the corners of a small,
                               #   heavily masked cutout) or "none".  The visit image is already
                               #   background-subtracted, so only a low-order residual is left
    bgSatMaskFwhm=6.0,      # background refit ignores pixels within this many FWHM of
                               #   saturated pixels (their halos were absorbed into the sky)
    dropHaloMimics=True,     # end of each outer pass: remove non-Gaia sources next to bright
                               #   stars that checkBrightNeighbours() calls 'doubtful', and
                               #   do not re-detect them in later passes
    haloCheckSnr=300.0,      # ... bright = peak S/N >= this
    haloCheckFwhm=4.0,       # ... neighbours within this many FWHM
    haloSimilarMin=3,        # ... 'doubtful' if >= this many control angles show >= half
                               #   its flux (a repeating halo pattern), or significance < 3
    autoPoisson=True,         # after each outer pass: if the residual on bright stars is
                               #   PIXEL-SCALE noise (diagnoseResidual), measure the photon noise
                               #   the variance plane is missing (k per unit model flux) and use
                               #   var + k*model from then on (fit weights, S/N, residual check)
    oneSourcePerSat=True,   # keep only the brightest source on each saturated core (its
                               #   masked centre cannot separate 2-3 sources; they split the
                               #   flux and leave blotches) -- also applies to Gaia entries
    # ---- extended sources / shape cut
    keepExtended=True,        # a source the shape cut finds too BROAD is kept with a
                               # Gaussian-broadened PSF (small galaxy) instead of deleted,
                               # so its flux does not stay in the residual
    shapeSnrTest=False,      # the shape cut's own local-S/N test (eliminateWeak already
                               # enforces minSnr; the two disagreeing caused churn)
    laterShiftCap=0.5,       # centroid step cap after the first iteration of a pass
    dropRingMimics=True,     # end of every inner loop: a faint NON-Gaia source near a bright star
    ringCheckSnr=150.0,      #   (peak S/N >= ringCheckSnr, within ringCheckFwhm*FWHM, flux <
    ringCheckFwhm=4.0,       #   ringMaxRatio of the bright star) must stand out from the
    ringMaxRatio=0.05,       #   AZIMUTHAL MEDIAN of the residual at its distance from the star
    ringSnr=5.0,              #   with S/N >= ringSnr -- otherwise it is a fake absorbing the
                               #   star's wing excess: removed and not re-detected
    ringBgModes=4,           #   the wing level at each radius includes angular modes m <= this,
                               #   not just the azimuthal median
    ringSym=True,             #   ... and the source must also stand out from the same spot rotated
                               #   by 90/180/270 deg about the star (4-fold diffraction pattern):
                               #   a real star is in one place, a spike/lobe is in all four
    splitTest=True,           # end of every inner loop: a source (S/N splitMinSnr .. brightCorr
    splitMinSnr=20.0,        #   _sigma) with a positive residual peak (matched-filter S/N >=
    splitPeakSnr=5.0,        #   splitPeakSnr) 0.5-1.5 FWHM away gets a trial second source
    splitDchi2=25.0,          #   there; kept if chi2 within 2 FWHM drops by >= splitDchi2 and its
                               #   own S/N >= minSnr (unresolved pairs -> dipoles)
    ellipTest=True,           # end of every inner loop: a source (Gaia too) with S/N ellipMinSnr ..
    ellipMinSnr=15.0,        #   brightCorrSigma and local chi2/px >= ellipChi2Min is refitted,
    ellipChi2Min=1.3,        #   together with the non-Gaia points within ellipGroupFwhm*FWHM, as
    ellipGroupFwhm=1.5,      #   ONE elliptical extended source (PSF x elongated Gaussian); kept if
    ellipDchi2=10.0,          #   the local chi2 drops by >= ellipDchi2 (+5 per extra fitted number:
                               #   a LONE point source -> ellipse adds 3, so it needs >= 25)
    finalPolish=True,         # end of every inner loop (after all removals): re-optimise ALL
    finalPolishRounds=4,     #   positions jointly (Gauss-Newton steps capped at finalPolishCap
    finalPolishCap=0.3,      #   px, Gaia drift limit kept), kept only while the total chi2
                               #   drops -- removes dipoles from positions left stale when
                               #   neighbours were added / removed
    rememberRejects=True,     # do not re-add candidates already rejected in this pass
    # ---- speed
    fastStamps=True,          # render PSF stamps / derivatives with the C++ kernels (NumPy if not
                               #   built) for spatially CONSTANT PSFs; each PSF is first checked
                               #   against psf.computeImage and the fast path is only used if they
                               #   agree to 2e-3 of the peak
    stampFloorSigma=0.05,    # crop each source's PSF stamp where its wing drops below this
                               # many sigma per pixel (faint sources get small stamps);
                               # None = always the full psfWingSize stamp
)


# runCutoutPipeline defaults for the keys this driver uses
_LIB_DEFAULTS = dict(
    innerKwargs=None, supportFrac=0.01, refineCentroids=True, centroidMinAmp=50.0,
    nCentroidRounds=4, centroidShiftTol=0.03, centroidShiftSchedule=None,
    coarseCentroid=True, coarseCentroidHwSigma=2.5, mergeSeparationSigma=1.5,
    sysFrac=0.03, minSnr=5.0, pruneSepSigma=1.5, rejectFwhm=0.75,
    brightPeakSigma=100.0, freezeFwhm=1.0, poissonK=0.0, nearBrightFwhm=2.0,
    nearBrightSnr=10.0, edgeMarginFwhm=0.5, useEliminateWeak=True,
    useShapeCut=True, shapeKwargs=None, shapeInFirstPass=False,
    refitPsf=True, stackKwargs=None, psfPolishCycles=3, psfKernelSize=41,
    psfKernelSizes=None, bgBinSize=32, bgUseApprox=True, psfWingSize=81,
)


_IGNORED = {"nOuterIterations", "nInnerIterations", "outerConvergenceFrac",
            "plotEvery", "sourceView", "maxPagesPerIter", "sourcesPerPage",
            "viewHalf", "xRange", "yRange", "viz"}


# option names used before the "Gaia starting positions" naming -> current names
# (snake_case names are turned into camelCase first)
_OLD_KEYS = dict(seedThresholdStart="gaiaThresholdStart", seedPadFwhm="gaiaPadFwhm",
                 seedMaxDriftPx="gaiaMaxDriftPx", seedMinSepPx="gaiaMinSepPx",
                 protectSeeds="protectGaia")


def _params(cfg):
    """Merge a configuration with the fitLib and pipeline defaults (unknown keys raise)."""
    cfg = {_OLD_KEYS.get(k, k): v for k, v in camelKeys(cfg).items()}   # old names still work
    P = dict(_LIB_DEFAULTS)
    P.update(PIPELINE_DEFAULTS)
    unknown = set(cfg) - set(P) - _IGNORED
    if unknown:
        raise ValueError(f"unknown option(s) {sorted(unknown)}")
    P.update({k: v for k, v in cfg.items() if k in P})
    ik = dict(nIterations=10, thresholdStart=8.0, thresholdStep=-1.5, thresholdMin=4.5,
              minPixels=1)
    ik.update(P["innerKwargs"] or {})
    for k in ("thresholdStart", "thresholdStep", "thresholdMin"):
        if P[k] is None:
            P[k] = float(ik[k])
    P["minPixels"] = int(ik.get("minPixels", 1))
    if P["thresholdStep"] >= 0:
        P["thresholdStep"] = -abs(P["thresholdStep"]) or -1.0
    P["tFloor"] = min(P["tFloor"], P["thresholdMin"])
    return P


# =====================================================================
# Speedups (patched into fitLib)
# =====================================================================
STAMP_CACHE_MAX = 6000            # stamps kept per PSF (81x81 float64 = 52 kB each)


_STAMP_CACHE = {}                 # id(psf) -> (psf, OrderedDict)  (strong ref: no id reuse)


_HUGE = None


_ORIG = {}


_SIZE_RULE = dict(tree=None, half=None, tol=1.5)


_EXT_RULE = dict(tree=None, sig=None, q=None, th=None)   # positions -> Gaussian broadening [px]


__version__ = "2026-10-29a (meas_crowd package: fast C++ stamps, photometry, source cutouts)"


_BF = dict(r=None, prof=None, L0=None, tree=None, lf=None, tol=0.75, clip=1.5)


#   brightness-dependent PSF term: stamp of a star with flux f = K + (log10 f - L0) * prof(r)


def _setBfFluxes(xy, amps):
    """Tell the stamp builder each source's current flux (for the brightness term)."""
    if _BF["prof"] is None or len(xy) == 0:
        return
    a = np.asarray(amps, float)
    ok = a > 0
    _BF.update(tree=cKDTree(np.asarray(xy)[ok]) if ok.any() else None, lf=np.log10(a[ok]))


def _bfReset():
    """Clear the brightness term and the elliptical-source shapes, and the stamp cache."""
    _BF.update(r=None, prof=None, L0=None, tree=None, lf=None)
    _ELLIP.clear()
    _STAMP_CACHE.clear()


_ELLIP = {}       # source uid -> (axis ratio q <= 1, position angle [rad]) of an elliptical extended source


def _ellipKernel(sa, q, th, prof="gauss"):
    """Return a normalised elliptical 'gauss' or 'exp' profile.

    Major-axis scale ``sa`` [px], axis ratio ``q``, major axis at angle ``th``.  The
    exponential is integrated over 5x5 sub-pixels because its centre is a cusp.
    """
    sa = max(sa, 0.15)
    sb = max(q * sa, 0.15)
    h = int(np.ceil((4 if prof == "gauss" else 7) * sa)) + 1
    ns = 1 if prof == "gauss" else 5
    u = (np.arange(-h * ns, (h + 1) * ns) + 0.5) / ns - 0.5
    X, Y = np.meshgrid(u, u)
    c, s_ = np.cos(th), np.sin(th)
    a_ = (c * X + s_ * Y) / sa
    b_ = (-s_ * X + c * Y) / sb
    r2 = a_ ** 2 + b_ ** 2
    k = np.exp(-0.5 * r2) if prof == "gauss" else np.exp(-np.sqrt(r2))
    if ns > 1:
        k = k.reshape(2 * h + 1, ns, 2 * h + 1, ns).sum(axis=(1, 3))
    return k / k.sum()


def _setExtRule(S):
    """Register the extended sources' broadening (size, shape, profile) with the stamp builder."""
    m = S["ext"] > 0
    ell = [tuple(_ELLIP.get(int(u), (1.0, 0.0, "gauss"))) + ("gauss",) for u in S["uid"][m]]
    qs = np.array([e[0] for e in ell])
    ths = np.array([e[1] for e in ell])
    prs = [e[2] for e in ell]
    _EXT_RULE.update(tree=cKDTree(S["xy"][m]) if m.any() else None, sig=S["ext"][m], q=qs, th=ths,
                     prof=prs)


def clearStampCache():
    """Empty the PSF stamp cache and the stamp-size rule.

    Call after changing the PSF in place; a new PSF object gets its own cache automatically.
    """
    _STAMP_CACHE.clear()
    _SIZE_RULE["tree"] = None


def _cropStamp(arr, fx0, fy0, x, y, half):
    """Zero a stamp beyond radius ``half`` around (x, y), with a cosine taper over the outer quarter.

    The array keeps its shape (callers rely on it); the zeros never enter the fit.
    """
    yy, xx = np.mgrid[fy0:fy0 + arr.shape[0], fx0:fx0 + arr.shape[1]]
    w = np.clip((half - np.hypot(xx - x, yy - y)) / (0.25 * half), 0.0, 1.0)
    return arr * (0.5 - 0.5 * np.cos(np.pi * w))


def setStampSizes(xy, amps, psf, point, noise, fwhm, floorSigma=0.05, tol=1.5):
    """Choose a PSF stamp half-width per source.

    The smallest of a few tiers beyond which ``amp * PSF`` stays below ``floorSigma * noise``,
    so faint sources get small stamps.  Sources not in the list (new detections) get the full
    stamp.

    Parameters
    ----------
    xy : `numpy.ndarray`
        Source positions, shape (N, 2).
    amps : `numpy.ndarray`
        Source amplitudes.
    psf : `lsst.afw.detection.Psf`
        PSF model.
    point : `lsst.geom.Point2D`
        Position at which the PSF profile is evaluated.
    noise : `float`
        Typical pixel noise.
    fwhm : `float`
        PSF FWHM [px].
    floorSigma : `float`, optional
        Crop where the source's wing drops below this many noise sigma.
    tol : `float`, optional
        Positions within this many px of a listed source use its stamp size.

    Returns
    -------
    half : `numpy.ndarray` [`int`]
        Half-width per source (-1 = full stamp).
    """
    xy = np.asarray(xy, float).reshape(-1, 2)
    if floorSigma is None or not len(xy):
        _SIZE_RULE["tree"] = None
        return np.full(len(xy), -1)
    k = psf.computeKernelImage(point).array
    h = k.shape[0] // 2
    u = np.arange(-h, h + 1)
    rb = np.round(np.hypot(*np.meshgrid(u, u))).astype(int).ravel()
    prof = np.zeros(rb.max() + 2)
    np.maximum.at(prof, rb, np.abs(k).ravel())
    env = np.maximum.accumulate(prof[::-1])[::-1] / k.max()      # max |PSF| beyond radius r
    s = np.clip(np.asarray(amps, float), 0, None) * k.max() / noise
    minh = max(8, int(np.ceil(2.0 * fwhm)))
    half = np.full(len(s), -1)
    for t in sorted({max(minh, t) for t in (8, 12, 16, 24, 32)}, reverse=True):
        if t < h:
            half[s * env[t] < floorSigma] = t
    _SIZE_RULE.update(tree=cKDTree(xy), half=half, tol=tol)
    return half


def _psfStampCached(psf, x, y, clipBbox, taper=None):
    """Return the same result as `fitLib.psfStamp`, cached on the exact (x, y)."""
    ent = _STAMP_CACHE.get(id(psf))
    if ent is None or ent[0] is not psf:
        ent = (psf, OrderedDict())
        _STAMP_CACHE[id(psf)] = ent
    store = ent[1]
    half = -1
    if _SIZE_RULE["tree"] is not None:
        d, i = _SIZE_RULE["tree"].query((x, y))
        if d <= _SIZE_RULE["tol"]:
            half = int(_SIZE_RULE["half"][i])
    sig, eq, eth, epr = 0.0, 1.0, 0.0, "gauss"
    if _EXT_RULE["tree"] is not None:
        d, i = _EXT_RULE["tree"].query((x, y))
        if d <= 1e-6:
            sig, half = float(_EXT_RULE["sig"][i]), -1
            if _EXT_RULE.get("q") is not None and len(_EXT_RULE["q"]) > i:
                eq, eth = float(_EXT_RULE["q"][i]), float(_EXT_RULE["th"][i])
                epr = _EXT_RULE["prof"][i] if _EXT_RULE.get("prof") else "gauss"
    eps = 0.0
    if _BF["prof"] is not None and _BF["tree"] is not None and sig == 0:
        d, i = _BF["tree"].query((x, y))
        if d <= _BF["tol"]:
            eps = round(float(np.clip(_BF["lf"][i] - _BF["L0"], -_BF["clip"], _BF["clip"])) / 0.02) * 0.02
    key = (float(x), float(y), taper, half, sig, eps, eq, eth, epr)
    hit = store.get(key)
    if hit is None:
        arr, bb = _rawPsfStamp(psf, x, y, taper)
        if eps:
            yy, xx = np.mgrid[bb.getMinY():bb.getMinY() + arr.shape[0],
                              bb.getMinX():bb.getMinX() + arr.shape[1]]
            add = eps * np.interp(np.hypot(xx - x, yy - y), _BF["r"], _BF["prof"], right=0.0)
            arr = np.clip(arr + add, 0.0, None)
        if half > 0:
            arr = _cropStamp(arr, bb.getMinX(), bb.getMinY(), x, y, half)
        if sig > 0:
            tot = arr.sum()
            if eq < 0.999 or epr != "gauss":
                from scipy.signal import fftconvolve
                arr = fftconvolve(arr, _ellipKernel(sig, eq, eth, epr), mode="same")
            else:
                from scipy.ndimage import gaussian_filter
                arr = gaussian_filter(arr, sig, mode="constant")
            arr *= tot / max(arr.sum(), 1e-30)
        hit = (arr, bb.getMinX(), bb.getMinY(), bb.getMaxX(), bb.getMaxY())
        store[key] = hit
        if len(store) > STAMP_CACHE_MAX:
            store.popitem(last=False)
    else:
        store.move_to_end(key)
    arr, fx0, fy0, fx1, fy1 = hit
    bx0, by0 = max(fx0, clipBbox.getMinX()), max(fy0, clipBbox.getMinY())
    bx1, by1 = min(fx1, clipBbox.getMaxX()), min(fy1, clipBbox.getMaxY())
    if bx0 > bx1 or by0 > by1:
        return None, L.geom.Box2I()
    bb = L.geom.Box2I(L.geom.Point2I(bx0, by0), L.geom.Point2I(bx1, by1))
    return arr[by0 - fy0:by1 - fy0 + 1, bx0 - fx0:bx1 - fx0 + 1], bb


def _invDiag(M):
    """Return diag(M^-1) of a symmetric positive (semi)definite matrix."""
    try:
        c, low = cho_factor(M, lower=True, check_finite=False)
        Linv = solve_triangular(np.tril(c), np.eye(len(M)), lower=True, check_finite=False)
        return np.sum(Linv ** 2, axis=0)
    except (np.linalg.LinAlgError, ValueError):
        return np.diag(np.linalg.pinv(M))


def ampSnrFullFast(fit, noiseScale=1.0):
    """Return honest amplitude S/N values, computed block by block.

    Drop-in for `fitLib.ampSnrFull`: ``S/N = amp / sqrt(diag((A^T W A)^-1))``, inverted per
    independent group of overlapping sources with Cholesky instead of one dense pseudo-inverse.

    Parameters
    ----------
    fit : `dict`
        Joint fit (from `fitLib.jointPsfFitDiag`).
    noiseScale : `float`, optional
        Empirical noise scale dividing the S/N.

    Returns
    -------
    snr : `numpy.ndarray`
        S/N per source.
    """
    comps, var, good = fit["comps"], fit["variance"], fit["good"]
    n = len(comps)
    pix = -np.ones(good.shape, dtype=np.int64)
    idx = np.flatnonzero(good.ravel())
    pix.ravel()[idx] = np.arange(len(idx))
    R, C, V = [], [], []
    for j, c in enumerate(comps):
        if c is None:
            continue
        p = pix[c["ys"], c["xs"]]
        ok = (p >= 0) & (c["stamp"] != 0)
        w = 1.0 / np.sqrt(var[c["ys"], c["xs"]][ok])
        R.append(p[ok]); C.append(np.full(int(ok.sum()), j)); V.append(c["stamp"][ok] * w)
    if not R:
        return np.zeros(n)
    A = coo_matrix((np.concatenate(V), (np.concatenate(R), np.concatenate(C))),
                   shape=(len(idx), n)).tocsr()
    N = (A.T @ A).tocsr()
    _, lab = connected_components(N, directed=False)
    counts = np.bincount(lab)
    order = np.argsort(lab, kind="stable")
    starts = np.cumsum(counts) - counts
    d = N.diagonal()
    varAmp = np.full(n, np.inf)
    single = (counts[lab] == 1) & (d > 0)
    varAmp[single] = 1.0 / d[single]
    for c in np.flatnonzero(counts > 1):
        ii = order[starts[c]:starts[c] + counts[c]]
        varAmp[ii] = _invDiag(N[ii][:, ii].toarray())
    sig = np.sqrt(np.clip(varAmp, 1e-30, None)) * noiseScale
    return np.where(np.isfinite(sig), fit["amps"] / sig, 0.0)


NNLS_MIN_SOURCES = 300            # below this lsq_linear is already fast (often 0 iterations)


def _nnlsNormal(A, b):
    """Solve bounded (x >= 0) least squares exactly through the normal equations.

    With ``N = A^T A = R^T R``, ``||Ax - b||^2 = ||Rx - c||^2 + const`` (``c = R^-T A^T b``), so
    NNLS runs on an n x n problem instead of the (pixels x n) one.
    """
    from scipy.optimize import nnls
    N = (A.T @ A).toarray()
    r = A.T @ b
    s = np.sqrt(np.clip(np.diag(N), 1e-300, None))
    Ns = N / np.outer(s, s)
    Ns[np.diag_indices_from(Ns)] += 1e-10
    rs = r / s
    c, _ = cho_factor(Ns, lower=False, check_finite=False)
    Rm = np.triu(c)
    cc = solve_triangular(Rm.T, rs, lower=True, check_finite=False)
    x = solve_triangular(Rm, cc, lower=False, check_finite=False)
    if (x < 0).any():
        try:
            x, _ = nnls(Rm, cc, maxiter=50 * len(x))
        except TypeError:
            x, _ = nnls(Rm, cc)
    return np.where(np.diag(N) > 0, x / s, 0.0)


def _lsqLinearFast(A, b, bounds=(-np.inf, np.inf), **kw):
    """Drop-in for scipy's lsq_linear that solves large non-negative joint fits with `_nnlsNormal`."""
    from scipy.optimize import OptimizeResult
    from scipy.sparse import issparse
    try:
        isNonneg = tuple(bounds) == (0, np.inf)
    except Exception:
        isNonneg = False
    if A.shape[1] < NNLS_MIN_SOURCES or not issparse(A) or not isNonneg:
        return _ORIG["lsq_linear"](A, b, bounds=bounds, **kw)
    try:
        x = _nnlsNormal(A, b)
    except Exception:
        return _ORIG["lsq_linear"](A, b, bounds=bounds, **kw)
    res = A @ x - b
    return OptimizeResult(x=x, status=1, nit=0, cost=0.5 * float(res @ res), success=True)


def profile(fn, *args, top=30, **kw):
    """Run a function under cProfile and print where the time goes.

    For example ``runs = crowd.profile(crowd.runSessionGaia, sess, cfg, maxOuter=2)``.

    Parameters
    ----------
    fn : callable
        Function to profile.
    *args
        Positional arguments of ``fn``.
    top : `int`, optional
        Number of entries printed.
    **kw
        Keyword arguments of ``fn``.

    Returns
    -------
    result : `object`
        What ``fn`` returned.
    """
    import cProfile
    import pstats
    pr = cProfile.Profile()
    pr.enable()
    try:
        return fn(*args, **kw)
    finally:
        pr.disable()
        st = pstats.Stats(pr)
        st.strip_dirs()
        print("\n===== most time spent INSIDE these functions (tottime) =====")
        st.sort_stats("tottime").print_stats(top)
        print("\n===== pipeline steps incl. everything they call (cumulative) =====")
        st.sort_stats("cumulative").print_stats(r"cutout_|lsq|detect|computeImage|nnls", top)


_FAST = {"on": True, "cache": {}, "verbose": True}   # fast C++/NumPy stamps for constant PSFs


def _fastFor(psf):
    """Return a validated FastPsf for a spatially constant PSF, else None (cached per PSF)."""
    if not _FAST["on"]:
        return None
    ent = _FAST["cache"].get(id(psf))
    if ent is not None and ent[0] is psf:
        return ent[1]
    from .fastKernels import FastPsf, HAVE_CXX
    fast = None
    try:
        c = psf.getAveragePosition()
        cx, cy = c.getX(), c.getY()
        fast = FastPsf(psf.computeKernelImage(c).array)

        def ref(x, y):
            arr, bb = _ORIG["psfStamp"](psf, x, y, _HUGE, taper=None)
            return arr, (bb.getMinX(), bb.getMinY())

        pos = [(cx + 0.31, cy - 0.17), (cx - 57.45, cy + 33.62), (cx + 81.08, cy + 90.49),
               (cx - 120.5, cy - 64.27)]
        ok, worst = fast.validate(ref, pos)
        if _FAST["verbose"]:
            print(f"  fast PSF stamps ({'C++' if HAVE_CXX else 'NumPy'}): "
                  + (f"ON (max difference to psf.computeImage {worst:.1e} of the peak)" if ok else
                     f"off for this PSF (differs from psf.computeImage by {worst:.1e} of the peak "
                     f"-- spatially varying PSF?)"))
        if not ok:
            fast = None
    except Exception as e:                   # PSF types without these methods etc.
        if _FAST["verbose"]:
            print(f"  fast PSF stamps: off ({type(e).__name__}: {e})")
        fast = None
    _FAST["cache"][id(psf)] = (psf, fast)
    return fast


def _box(b):
    """Convert (xmin, ymin, xmax, ymax) into an `lsst.geom.Box2I`."""
    return L.geom.Box2I(L.geom.Point2I(int(b[0]), int(b[1])), L.geom.Point2I(int(b[2]), int(b[3])))


def _rawPsfStamp(psf, x, y, taper):
    """Return the full PSF stamp at (x, y): fast path for constant PSFs, else psf.computeImage."""
    fast = _fastFor(psf)
    if fast is not None:
        arr, b = fast.stamp(x, y, (-10 ** 8, -10 ** 8, 10 ** 8, 10 ** 8), taper)
        return arr, _box(b)
    return _ORIG["psfStamp"](psf, x, y, _HUGE, taper=taper)


def _psfDerivativesFast(psf, x, y, bboxClip, eps=0.25):
    """Drop-in for `fitLib.computePsfDerivatives` that uses the fast path when possible."""
    fast = _fastFor(psf)
    if fast is None:
        return _ORIG["computePsfDerivatives"](psf, x, y, bboxClip, eps)
    clip = (bboxClip.getMinX(), bboxClip.getMinY(), bboxClip.getMaxX(), bboxClip.getMaxY())
    dx, dy, b = fast.derivatives(x, y, clip, eps)
    if dx is None:
        return None, None, L.geom.Box2I()
    return dx, dy, _box(b)


def enableSpeedups():
    """Patch fitLib with the cached / fast PSF stamps, fast S/N and fast NNLS (idempotent).

    Called by the pipeline itself; also speeds up fitLib's own `fitLib.runCutoutPipeline`.
    """
    global _HUGE
    if "psfStamp" not in _ORIG:
        _ORIG["psfStamp"] = L.psfStamp
        _ORIG["ampSnrFull"] = L.ampSnrFull
        _ORIG["lsq_linear"] = L.lsq_linear
        _ORIG["computePsfDerivatives"] = L.computePsfDerivatives
        _HUGE = L.geom.Box2I(L.geom.Point2I(-10 ** 8, -10 ** 8), L.geom.Extent2I(2 * 10 ** 8, 2 * 10 ** 8))
    L.psfStamp = _psfStampCached
    L.ampSnrFull = ampSnrFullFast
    L.lsq_linear = _lsqLinearFast
    L.computePsfDerivatives = _psfDerivativesFast


def disableSpeedups():
    """Undo `enableSpeedups` and empty the caches."""
    if "psfStamp" in _ORIG:
        L.psfStamp = _ORIG["psfStamp"]
        L.ampSnrFull = _ORIG["ampSnrFull"]
        L.lsq_linear = _ORIG["lsq_linear"]
        L.computePsfDerivatives = _ORIG["computePsfDerivatives"]
    _FAST["cache"].clear()
    clearStampCache()


# =====================================================================
# Source state: positions + Gaia bookkeeping, kept aligned through every filter
# =====================================================================
def _state(xy=None, fromGaia=False, gaia=None, gmag=None):
    """Return a new source state.

    Keys: ``xy`` (N, 2) positions, ``fromGaia`` (started from a Gaia position), ``ext``
    (extended-source size, 0 = point source), ``g`` (Gaia G), ``ref`` (Gaia reference
    position), ``gaia`` (catalogue row), ``uid`` (unique id) and ``nextUid``; all arrays stay
    aligned through every filter.
    """
    xy = np.asarray(xy if xy is not None else np.zeros((0, 2)), float).reshape(-1, 2)
    n = len(xy)
    g = np.full(n, np.nan) if gmag is None else np.asarray(gmag, float).reshape(n)
    return dict(xy=xy.copy(), fromGaia=np.full(n, bool(fromGaia)), ext=np.zeros(n), g=g,
                ref=xy.copy() if fromGaia else np.full((n, 2), np.nan),
                gaia=(np.asarray(gaia, np.int64) if gaia is not None else np.full(n, -1, np.int64)),
                uid=np.arange(n, dtype=np.int64), nextUid=n)


def _take(S, keep):
    """Return the source state restricted to ``keep`` (index or boolean array)."""
    keep = np.asarray(keep)
    out = {k: S[k][keep] for k in ("xy", "fromGaia", "ext", "g", "ref", "gaia", "uid")}
    out["nextUid"] = S["nextUid"]
    return out


def _add(S, xyNew):
    """Return the source state with non-Gaia sources appended at ``xyNew``."""
    xyNew = np.asarray(xyNew, float).reshape(-1, 2)
    m = len(xyNew)
    out = dict(xy=np.vstack([S["xy"], xyNew]), fromGaia=np.concatenate([S["fromGaia"], np.zeros(m, bool)]),
               ext=np.concatenate([S["ext"], np.zeros(m)]),
               g=np.concatenate([S["g"], np.full(m, np.nan)]),
               ref=np.vstack([S["ref"], np.full((m, 2), np.nan)]),
               gaia=np.concatenate([S["gaia"], np.full(m, -1, np.int64)]),
               uid=np.concatenate([S["uid"], S["nextUid"] + np.arange(m)]),
               nextUid=S["nextUid"] + m)
    return out


def _xy(S):
    """Return the source positions as a list of (x, y) tuples."""
    return [tuple(p) for p in S["xy"]]


def _clampToGaia(S, maxDrift):
    """Pull sources started from Gaia back to within ``maxDrift`` px of their Gaia position (in place)."""
    if not maxDrift or not S["fromGaia"].any():
        return
    s = S["fromGaia"]
    d = S["xy"][s] - S["ref"][s]
    r = np.hypot(d[:, 0], d[:, 1])
    over = r > maxDrift
    d[over] *= (maxDrift / r[over])[:, None]
    S["xy"][s] = S["ref"][s] + d


def estimateGaiaOffset(exposure, S, fwhm, badBitmask, noise, minPeakSnr=20.0,
                       minStars=5, verbose=True):
    """Measure the global offset between the image and the Gaia positions.

    The median (data centroid - Gaia position) of bright, unsaturated, fairly isolated Gaia sources,
    from background-subtracted windowed centroids.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The cutout.
    S : `dict`
        Source state with the Gaia starting positions (see ``_state``).
    fwhm : `float`
        PSF FWHM [px].
    badBitmask : `int`
        Mask bits of pixels to ignore.
    noise : `float`
        Typical pixel noise.
    minPeakSnr : `float`, optional
        Minimum peak S/N of a star used.
    minStars : `int`, optional
        Minimum number of stars; otherwise no offset is applied.
    verbose : `bool`, optional
        Print the result.

    Returns
    -------
    dx, dy : `float`
        Offset (data - Gaia) [px]; 0 if too few stars.
    nUsed : `int`
        Number of stars used.
    """
    data = exposure.image.array
    bad = (exposure.mask.array & badBitmask) != 0
    bb = exposure.getBBox()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    ny, nx = data.shape
    hw = max(3, int(round(1.25 * fwhm)))
    d = []
    others = cKDTree(S["xy"]) if len(S["xy"]) > 1 else None
    for x, y in S["xy"][S["fromGaia"]]:
        ix, iy = int(round(x)) - X0, int(round(y)) - Y0
        if not (hw <= ix < nx - hw and hw <= iy < ny - hw):
            continue
        if others is not None and len(others.query_ball_point([x, y], 1.5 * fwhm)) > 1:
            continue
        st = data[iy - hw:iy + hw + 1, ix - hw:ix + hw + 1]
        if bad[iy - hw:iy + hw + 1, ix - hw:ix + hw + 1].any():
            continue
        border = np.concatenate([st[0], st[-1], st[:, 0], st[:, -1]])
        w = np.clip(st - np.median(border), 0, None)
        if w.max() < minPeakSnr * noise or w.sum() <= 0:
            continue
        yy, xx = np.mgrid[-hw:hw + 1, -hw:hw + 1]
        cx = (xx * w).sum() / w.sum() + ix + X0
        cy = (yy * w).sum() / w.sum() + iy + Y0
        d.append((cx - x, cy - y))
    if len(d) < minStars:
        if verbose:
            print(f"  Gaia offset: only {len(d)} usable stars -- not applied")
        return 0.0, 0.0, len(d)
    d = np.array(d)
    med = np.median(d, axis=0)
    mad = 1.4826 * np.median(np.abs(d - med), axis=0) / np.sqrt(len(d))
    if verbose:
        print(f"  Gaia offset (data - Gaia) from {len(d)} stars: dx={med[0]:+.3f}+/-{mad[0]:.3f}, "
              f"dy={med[1]:+.3f}+/-{mad[1]:.3f} px")
    return float(med[0]), float(med[1]), len(d)


# =====================================================================
# Index-preserving versions of the cleaning steps (Gaia sources protected)
# =====================================================================
def _fit(exposure, psf, S, badBitmask, varModel, fitkw):
    """Run the joint fit of the source state's positions (extended sources broadened)."""
    _setExtRule(S)
    return L.jointPsfFitDiag(exposure, psf, _xy(S), badBitmask, varModel=varModel, **fitkw)


def _refineSubset(exposure, psf, S, fit, idx, maxShift, minAmp, freezeSep, psfSigma):
    """Make one Gauss-Newton position step for the sources ``idx`` against the current residual.

    Sources fainter than ``minAmp``, extended ones and ones with a brighter neighbour within
    ``freezeSep`` stay put.  Returns the new positions of all sources.
    """
    eb = exposure.getBBox()
    data, var, good, amps = exposure.image.array, fit["variance"], fit["good"], fit["amps"]
    resid = data - fit["modelFull"]
    xy = S["xy"]
    new = xy.copy()
    tree = cKDTree(xy) if len(xy) > 1 else None
    R2 = (2.5 * psfSigma) ** 2
    for j in idx:
        a = amps[j]
        if a < minAmp or S["ext"][j] > 0:
            continue
        x, y = xy[j]
        if tree is not None and freezeSep:
            if any(k != j and amps[k] > a for k in tree.query_ball_point(xy[j], freezeSep)):
                continue
        dpdx, dpdy, bb = L.computePsfDerivatives(psf, x, y, eb)
        if dpdx is None:
            continue
        ys, xs = L.bboxToSlices(bb, eb)
        yy, xx = np.mgrid[bb.getMinY():bb.getMaxY() + 1, bb.getMinX():bb.getMaxX() + 1]
        g = good[ys, xs] & (((xx - x) ** 2 + (yy - y) ** 2) < R2)
        if g.sum() < 10:
            continue
        Pm, _ = L.psfStamp(psf, x, y, bb, taper=L.WING_TAPER)
        if Pm is None or Pm.shape != dpdx.shape:
            continue
        w = 1.0 / np.sqrt(var[ys, xs][g])
        J = np.column_stack([Pm[g], a * dpdx[g], a * dpdy[g]]) * w[:, None]
        try:
            delta, *_ = np.linalg.lstsq(J, resid[ys, xs][g] * w, rcond=None)
        except np.linalg.LinAlgError:
            continue
        new[j] = (x + np.clip(delta[1], -maxShift, maxShift),
                  y + np.clip(delta[2], -maxShift, maxShift))
    return new


def _pruneMerge(fit, S, minSnr, minSep, bright, nearSep, nearSnr, protect,
                inZone=None, zoneSnr=None):
    """Prune low-S/N sources and merge close pairs, keeping the bookkeeping; return (alive, positions).

    Protected (Gaia) sources are never pruned; a non-protected source within ``minSep`` of a
    protected one is absorbed by it.
    """
    var, amps, good = fit["variance"], fit["amps"], fit["good"]
    xy = S["xy"].copy()
    n = len(xy)
    snr = np.zeros(n)
    for j, c in enumerate(fit["comps"]):
        if c is None:
            continue
        g = good[c["ys"], c["xs"]]
        snr[j] = amps[j] * np.sqrt(np.sum(c["stamp"][g] ** 2 / var[c["ys"], c["xs"]][g]))
    thresh = np.full(n, float(minSnr))
    if nearSep and bright is not None and np.any(bright):
        dd, _ = cKDTree(xy[bright]).query(xy)
        near = (dd > 0.5) & (dd < nearSep)
        thresh[near] = np.maximum(thresh[near], nearSnr)
    nZoneLow = 0
    if inZone is not None and zoneSnr and inZone.any():
        thresh[inZone] = np.maximum(thresh[inZone], zoneSnr)
        nZoneLow = int((inZone & (snr < thresh) & ~S["fromGaia"]).sum())
    prot = S["fromGaia"] if protect else np.zeros(n, bool)
    alive = (snr >= thresh) | prot
    nLow = int((~alive).sum())
    tree = cKDTree(S["xy"]) if n else None
    for a in np.lexsort((-amps, ~prot)):          # protected first, then brightest first
        if not alive[a]:
            continue
        for b in tree.query_ball_point(S["xy"][a], minSep):
            if b == a or not alive[b] or prot[b]:
                continue
            if not prot[a] and amps[a] + amps[b] > 0:
                xy[a] = (amps[a] * xy[a] + amps[b] * xy[b]) / (amps[a] + amps[b])
            alive[b] = False
    print(f"  prune/merge: {n} -> {int(alive.sum())} ({nLow} low S/N"
          + (f", of which {nZoneLow} in saturated-star zones" if nZoneLow else "")
          + f", {n - nLow - int(alive.sum())} merged)")
    return alive, xy


def _eliminateWeak(exposure, psf, S, fit, badBitmask, P, fitkw, fwhm=None):
    """Backward elimination on the honest S/N for non-Gaia sources (Gaia sources are kept)."""
    n0 = len(S["xy"])
    rounds = 0
    scale = L.empiricalNoiseScale(fit, exposure)
    minSnr = P["minSnr"]
    for _ in range(40):
        n = len(S["xy"])
        if n <= 1:
            break
        snr = ampSnrFullFast(fit, scale)
        cand = ~S["fromGaia"] if P["protectGaia"] else np.ones(n, bool)
        weak = np.flatnonzero(cand & (snr < minSnr))
        if len(weak) == 0:
            break
        rounds += 1
        hopeless = snr[weak] < P["hopelessFrac"] * minSnr
        isolated = np.ones(len(weak), bool)          # no other weak source nearby: removing
        if fwhm and len(weak) > 1:                   # it barely changes anyone else's S/N
            cnt = np.array([len(l) for l in cKDTree(S["xy"][weak]).query_ball_point(
                S["xy"][weak], 2.0 * fwhm)])
            isolated = cnt <= 1
        rest = weak[~(hopeless | isolated)]
        rest = rest[np.argsort(snr[rest])][:int(np.ceil(P["elimDropFrac"] * len(rest)))]
        drop = np.concatenate([weak[hopeless | isolated], rest])[:n - 1]
        keep = np.ones(n, bool)
        keep[drop] = False
        S = _take(S, keep)
        fit = _fit(exposure, psf, S, badBitmask, fit["modelFull"], fitkw)
    print(f"  eliminateWeak: noiseScale={scale:.2f}, {n0} -> {len(S['xy'])} sources in "
          f"{rounds} round(s) (honest S/N>={minSnr:g}{', Gaia protected' if P['protectGaia'] else ''})")
    return S, fit


def _shapeCut(exposure, psf, S, fit, badBitmask, P, fitkw, dropFrac=0.25, maxRounds=15,
              judgeUids=None, psfSigma=None):
    """Remove non-PSF-like non-Gaia sources; too-broad ones become extended sources if allowed.

    ``judgeUids``: only these sources are judged (None = all non-Gaia point sources).
    """
    kw = dict(P["shapeKwargs"] or {})
    mkw = {k: kw.pop(k) for k in ("bgSigma", "widthScales", "minPix") if k in kw}
    if not P["shapeSnrTest"]:
        kw["minSnrLocal"] = 1e-3
    n0, rejected, nExt = len(S["xy"]), [], 0
    why = dict(snr=0, chi2=0, width=0, contrast=0)
    scale = L.empiricalNoiseScale(fit, exposure)
    pscale = L.empiricalNoiseScale(fit, exposure, matched=False)
    for _ in range(maxRounds):
        n = len(S["xy"])
        ok = (S["ext"] == 0) & (~S["fromGaia"] if P["protectGaia"] else True)
        if judgeUids is not None:
            ok &= np.isin(S["uid"], list(judgeUids))
        cand = np.flatnonzero(ok)
        if n <= 1 or len(cand) == 0:
            break
        sub = dict(fit, comps=[fit["comps"][j] for j in cand])
        m = L.shapeMetrics(exposure, psf, [tuple(S["xy"][j]) for j in cand], sub,
                           noiseScale=scale, chi2Scale=pscale, **mkw)
        fail, score, flags = L.judgeShape(m, **kw)
        bad = np.flatnonzero(fail)
        if len(bad) == 0:
            break
        broad = (flags["width"] & ~flags["chi2"] & ~flags["contrast"] &
                 (np.nan_to_num(m["bestScale"]) > 1.0))
        if P["keepExtended"] and psfSigma and broad[bad].any():
            conv = bad[broad[bad]]
            S["ext"][cand[conv]] = psfSigma * np.sqrt(np.clip(m["bestScale"][conv] ** 2 - 1, 0.05, None))
            nExt += len(conv)
            bad = bad[~broad[bad]]
            if len(bad) == 0:
                fit = _fit(exposure, psf, S, badBitmask, fit["modelFull"], fitkw)
                continue
        order = bad[np.argsort(-score[bad])]
        sel = order[:min(max(1, int(np.ceil(dropFrac * len(bad)))), n - 1)]
        for i in sel:
            for k in why:
                why[k] += int(flags[k][i])
        drop = cand[sel]
        rejected += [(S["xy"][j][0], S["xy"][j][1], float(score[i])) for j, i in zip(drop, sel)]
        keep = np.ones(n, bool)
        keep[drop] = False
        S = _take(S, keep)
        fit = _fit(exposure, psf, S, badBitmask, fit["modelFull"], fitkw)
    print(f"  shape cut: {n0} -> {len(S['xy'])} ({len(rejected)} removed, {nExt} kept as extended; "
          f"removed for chi2 {why['chi2']}, width {why['width']}, contrast {why['contrast']}"
          + (f", S/N {why['snr']}" if P["shapeSnrTest"] else "") + ")")
    return S, fit, rejected


def _satlike(exposure):
    """Return the SAT (and NEARSAT, if present) pixels of an exposure."""
    m = exposure.mask
    out = (m.array & m.getPlaneBitMask("SAT")) != 0
    try:
        out |= (m.array & m.getPlaneBitMask("NEARSAT")) != 0
    except Exception:
        pass
    return out


def saturationLevel(exposure, P):
    """Estimate the saturation level of an image.

    ``nearSatLevel`` if given; else the median of the SAT-flagged pixels, unless those look
    interpolated (below the brightest unflagged pixel), in which case the brightest unflagged
    pixel, a lower bound.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The image.
    P : `dict`
        Pipeline options (uses ``nearSatLevel``).

    Returns
    -------
    level : `float` or `None`
        Saturation level in image units (None if it cannot be estimated).
    how : `str`
        How it was obtained.
    """
    if P.get("nearSatLevel"):
        return float(P["nearSatLevel"]), "given (nearSatLevel)"
    data = exposure.image.array
    sat = (exposure.mask.array & exposure.mask.getPlaneBitMask("SAT")) != 0
    fin = np.isfinite(data)
    if (sat & fin).sum() < 5:
        return None, "no saturated pixels in this cutout -- set nearSatLevel to use option B"
    lvl = float(np.nanmedian(data[sat & fin]))
    top = float(np.nanmax(data[~sat & fin]))
    if lvl >= 0.9 * top:
        return lvl, "median of the SAT-flagged pixels"
    return top, "brightest unflagged pixel (SAT pixel values look interpolated; a lower bound)"


def measureModelError(exposure, res, P, verbose=True):
    """Measure the PSF model error as a fraction of the model per pixel.

    Around unsaturated point sources with peak S/N >= ``modelErrMinSnr``, the robust (MAD)
    variance of residual/sigma in bins of model value is fitted as ``1 + eps^2 * model^2 / var``.
    Values not significant at 2 sigma give eps = 0.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    res : `dict`
        Result of `runInnerAdaptive`.
    P : `dict`
        Pipeline options.
    verbose : `bool`, optional
        Print the result.

    Returns
    -------
    eps : `float`
        Model error (0 if not significant).
    table : `list` [`tuple`]
        Per model bin: (model, number of pixels, robust variance of residual/sigma, model^2/var).
    """
    fit = res["snaps"][-1]["fit"]
    S = res["S"]
    model = np.clip(res["modelFull"], 0, None)
    data = exposure.image.array
    var = exposure.variance.array + P["poissonK"] * model
    bb = exposure.getBBox()
    c = bb.getCenter()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    psf = exposure.getPsf()
    kern = psf.computeKernelImage(c).array
    fwhm = L._fwhm(psf, c)
    noise = float(np.sqrt(np.nanmedian(exposure.variance.array)))
    amps = fit["amps"]
    snr = amps * kern.max() / noise
    sel = np.flatnonzero((snr >= P["modelErrMinSnr"]) & (S["ext"] == 0))
    satl = _satlike(exposure)
    if satl.any() and len(sel):
        from scipy.ndimage import binary_dilation
        near = binary_dilation(satl, iterations=int(np.ceil(3 * fwhm)))
        ix = np.clip(np.round(S["xy"][sel, 0]).astype(int) - X0, 0, satl.shape[1] - 1)
        iy = np.clip(np.round(S["xy"][sel, 1]).astype(int) - Y0, 0, satl.shape[0] - 1)
        sel = sel[~near[iy, ix]]
    if len(sel) < 3:
        return 0.0, []
    use = _starMask(data.shape, X0, Y0, S["xy"][sel], 2.5 * fwhm) & fit["good"] & ~satl
    use &= np.isfinite(data) & (var > 0)
    m, r, v = model[use], (data - res["modelFull"])[use], var[use]
    lo, hi = 5 * noise, m.max() if m.size else 0
    if hi <= lo * 2:
        return 0.0, []
    edges = np.geomspace(lo, hi, 13)
    rows, num, den, vv = [], 0.0, 0.0, 0.0
    z = r / np.sqrt(v)
    for a, b in zip(edges[:-1], edges[1:]):
        k = (m >= a) & (m < b)
        if k.sum() < 30:
            continue
        # robust variance of residual/sigma: 1 + eps^2 * model^2 / var
        zv = (1.4826 * np.median(np.abs(z[k] - np.median(z[k])))) ** 2
        q = float(np.median(m[k] ** 2 / v[k]))
        rows.append((float(np.median(m[k])), int(k.sum()), zv, q))
        num += k.sum() * (zv - 1.0) * q
        den += k.sum() * q * q
        vv += k.sum() * q * q * 2.7          # var(MAD variance) ~ 2.7 / n per bin
    if den <= 0:
        return 0.0, rows
    e2, se2 = num / den, np.sqrt(vv) / den
    signif = e2 > 2 * se2
    eps = float(np.sqrt(e2)) if signif else 0.0
    if verbose and rows:
        if signif:
            msg = f"eps = {100 * eps:.2f}% of the model per pixel"
        else:
            msg = (f"not significant (eps < {100 * np.sqrt(max(e2, 0) + 2 * se2):.2f}%, "
                   f"consistent with 0) -> eps = 0")
        print(f"  PSF model error: {msg} ({len(sel)} stars, S/N >= {P['modelErrMinSnr']:g}); "
              f"fit weights use sysFrac = {P['sysFrac']:.3g}")
    return eps, rows


def satZones(exposure, fwhm, P):
    """Find the zones around saturated stars.

    The zone band follows the SAT mask's shape: pixels within
    ``w = satZoneFwhm * FWHM + satZoneScale * inscribed radius`` of each blob.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The cutout.
    fwhm : `float`
        PSF FWHM [px].
    P : `dict`
        Pipeline options (``sat_zone_*``).

    Returns
    -------
    centres : `numpy.ndarray`
        (x, y, halo radius) per saturated blob, shape (N, 3).
    mask : `numpy.ndarray` [`bool`]
        The zone band.
    """
    from scipy.ndimage import label, center_of_mass, distance_transform_edt
    shape = exposure.image.array.shape
    empty = (np.zeros((0, 3)), np.zeros(shape, bool))
    if not P["satZoneFwhm"]:
        return empty
    sat = (exposure.mask.array & exposure.mask.getPlaneBitMask("SAT")) != 0
    if not sat.any():
        return empty
    lab, n = label(sat, structure=np.ones((3, 3)))
    area = np.bincount(lab.ravel())
    inner = distance_transform_edt(sat)             # depth inside the saturated region
    bb = exposure.getBBox()
    mask = np.zeros(shape, bool)
    cen = []
    for l in range(1, n + 1):
        if area[l] < P["satZoneMinPx"]:
            continue
        blob = lab == l
        rIn = float(inner[blob].max())
        w = P["satZoneFwhm"] * fwhm + P["satZoneScale"] * rIn
        mask |= distance_transform_edt(~blob) <= w
        cy, cx = center_of_mass(blob)
        cen.append((cx + bb.getMinX(), cy + bb.getMinY(), 2.0 * rIn + w))
        print(f"  saturated region at ({cen[-1][0]:7.1f}, {cen[-1][1]:7.1f}): {int(area[l])} px, "
              f"half-width {rIn:4.1f}px -> zone band {w:4.1f}px around it")
    if cen:
        print(f"  saturated-star zones cover {100 * mask.mean():.1f}% of the cutout; non-Gaia "
              f"sources there need S/N >= {P['satZoneSnr']:g}")
    return np.array(cen).reshape(-1, 3), mask


def satCoreLabels(exposure, fwhm, minPx=4):
    """Label the saturated cores, each grown by max(3 px, FWHM/2).

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The cutout.
    fwhm : `float`
        PSF FWHM [px].
    minPx : `int`, optional
        Ignore saturated blobs smaller than this (hot pixels, cosmic rays).

    Returns
    -------
    labels : `numpy.ndarray` [`int`]
        Core label per pixel (0 = none).
    """
    from scipy.ndimage import label, binary_dilation
    sat = (exposure.mask.array & exposure.mask.getPlaneBitMask("SAT")) != 0
    lab = np.zeros(sat.shape, np.int32)
    if not sat.any():
        return lab
    l0, n = label(sat, structure=np.ones((3, 3)))
    area = np.bincount(l0.ravel())
    grow = max(3, int(round(0.5 * fwhm)))
    k = 0
    for l in range(1, n + 1):
        if area[l] < minPx:
            continue
        k += 1
        lab[binary_dilation(l0 == l, iterations=grow) & (lab == 0)] = k
    return lab


def _labelAt(xy, lab, X0, Y0):
    """Return the label image's value at each source position (0 outside)."""
    xy = np.asarray(xy, float).reshape(-1, 2)
    if not len(xy) or not lab.any():
        return np.zeros(len(xy), np.int32)
    ix = np.round(xy[:, 0]).astype(int) - X0
    iy = np.round(xy[:, 1]).astype(int) - Y0
    ok = (ix >= 0) & (ix < lab.shape[1]) & (iy >= 0) & (iy < lab.shape[0])
    out = np.zeros(len(xy), np.int32)
    out[ok] = lab[iy[ok], ix[ok]]
    return out


def _onePerCore(S, amps, lab, X0, Y0):
    """Return a keep-mask with one source per saturated core.

    The brightest Gaia star if Gaia has one there (its position is good to milliarcseconds; the
    faint Gaia entries around bright stars are usually spurious), otherwise the source closest to
    the core's centre -- never by fitted flux, which is degenerate when the centre is masked.
    """
    keep = np.ones(len(S["xy"]), bool)
    core = _labelAt(S["xy"], lab, X0, Y0)
    for l in np.unique(core[core > 0]):
        idx = np.flatnonzero(core == l)
        if len(idx) <= 1:
            continue
        keep[idx] = False
        g = S["g"][idx]
        if np.isfinite(g).any():
            keep[idx[np.nanargmin(g)]] = True
        else:
            cy, cx = np.argwhere(lab == l).mean(axis=0)
            d = np.hypot(S["xy"][idx, 0] - (cx + X0), S["xy"][idx, 1] - (cy + Y0))
            keep[idx[np.argmin(d)]] = True
    return keep


def _refineSaturated(exposure, psf, S, fit, coreLab, X0, Y0, fwhm, nIter=6, maxStep=2.0):
    """Fit the positions of saturated stars without Gaia positions on their unsaturated wings.

    Gauss-Newton on an annulus out to ~3 FWHM beyond the core, other sources held fixed.  Returns
    the new positions of all sources and the number moved.
    """
    eb = exposure.getBBox()
    data, var, good = exposure.image.array, fit["variance"], fit["good"]
    amps = fit["amps"]
    core = _labelAt(S["xy"], coreLab, X0, Y0)
    xy = S["xy"].copy()
    moved = 0
    for j in np.flatnonzero(core > 0):
        if S["fromGaia"][j] and np.isfinite(S["g"][j]):
            continue
        a = amps[j]
        if a <= 0 or fit["comps"][j] is None:
            continue
        c = fit["comps"][j]
        own = np.zeros_like(data)
        own[c["ys"], c["xs"]] = a * c["stamp"]
        others = fit["modelFull"] - own
        rcore = np.sqrt((coreLab == core[j]).sum() / np.pi)
        x, y = xy[j]
        for _ in range(nIter):
            dpdx, dpdy, bb = L.computePsfDerivatives(psf, x, y, eb)
            if dpdx is None:
                break
            ys, xs = L.bboxToSlices(bb, eb)
            Pm, _ = L.psfStamp(psf, x, y, bb, taper=L.WING_TAPER)
            if Pm is None or Pm.shape != dpdx.shape:
                break
            yy, xx = np.mgrid[bb.getMinY():bb.getMaxY() + 1, bb.getMinX():bb.getMaxX() + 1]
            g = good[ys, xs] & (np.hypot(xx - x, yy - y) < rcore + 3 * fwhm)
            if g.sum() < 30:
                break
            w = 1.0 / np.sqrt(var[ys, xs][g])
            r = (data[ys, xs] - others[ys, xs] - a * Pm)[g] * w
            J = np.column_stack([Pm[g], a * dpdx[g], a * dpdy[g]]) * w[:, None]
            try:
                d, *_ = np.linalg.lstsq(J, r, rcond=None)
            except np.linalg.LinAlgError:
                break
            a = a + d[0]
            step = np.clip(d[1:], -maxStep, maxStep)
            x, y = x + step[0], y + step[1]
            if np.hypot(*step) < 0.02:
                break
        if np.hypot(x - xy[j, 0], y - xy[j, 1]) > 0.05:
            moved += 1
        xy[j] = (x, y)
    return xy, moved


def _zoneFlags(xy, zoneMask, X0, Y0):
    """Return whether each source lies in the saturated-star zone band."""
    xy = np.asarray(xy, float).reshape(-1, 2)
    if not len(xy) or not zoneMask.any():
        return np.zeros(len(xy), bool)
    ix = np.clip(np.round(xy[:, 0]).astype(int) - X0, 0, zoneMask.shape[1] - 1)
    iy = np.clip(np.round(xy[:, 1]).astype(int) - Y0, 0, zoneMask.shape[0] - 1)
    return zoneMask[iy, ix]


def brightStarCorrections(exposure, S, fit, model, psfPeak, noise, fwhm, badBitmask,
                          sigmaMin=100.0, radiusFwhm=4.0, sweeps=2, radial=True,
                          zones=None, haloModes=0, brightModes=0, coreRings=(),
                          coreSigma=0.6, coreModes=0):
    """Fit local correction terms for every bright or saturated star.

    On top of the joint model, within ``radiusFwhm`` FWHM, each bright star's residual is fitted
    with ``amp * (a*P + bx*dP/dx + by*dP/dy + cxx*d2P/dx2 + cyy*d2P/dy2 + cxy*d2P/dxdy)`` (flux,
    position, size and ellipticity), optional radial rings with angular modes (halo and wing
    structure) and optional narrow core terms (core sharper than the PSF, possibly off-centre).
    Weighted by the pixel variance (not the sysFrac-inflated one) so the core, where the
    mismatch is, counts.  Two sweeps, brightest first, for overlapping bright stars.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    S : `dict`
        Source state.
    fit : `dict`
        Joint fit.
    model : `numpy.ndarray`
        Joint model.
    psfPeak : `float`
        Peak value of the unit-flux PSF.
    noise : `float`
        Typical pixel noise.
    fwhm : `float`
        PSF FWHM [px].
    badBitmask : `int`
        Mask bits of pixels to ignore.
    sigmaMin : `float`, optional
        Stars with peak S/N above this (or saturated) are corrected.
    radiusFwhm : `float`, optional
        Fit radius [FWHM].
    sweeps : `int`, optional
        Number of sweeps.
    radial : `bool`, optional
        Add radial ring terms.
    zones : `numpy.ndarray`, optional
        Saturated-star zones (see `satZones`); stars in them get halo terms out to the zone.
    haloModes : `int`, optional
        Angular modes of the saturated stars' rings.
    brightModes : `int`, optional
        Angular modes of the unsaturated stars' rings.
    coreRings : `tuple` [`float`], optional
        Radii [px] of the narrow core terms (unsaturated stars).
    coreSigma : `float`, optional
        Width [px] of the narrow core terms.
    coreModes : `int`, optional
        Angular modes of the narrow core terms.

    Returns
    -------
    corr : `numpy.ndarray`
        Correction image (add to the joint model).
    info : `dict` [`int`, `dict`]
        Per corrected star: ``x``, ``y``, ``snr``, ``nSat``, ``flux``, ``dfluxPct``, ``dx``,
        ``dy``, ``chi2Before``, ``chi2After``.
    """
    data, good = exposure.image.array, fit["good"]
    var = exposure.variance.array
    bb = exposure.getBBox()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    amps = fit["amps"]
    sat = (exposure.mask.array & exposure.mask.getPlaneBitMask("SAT")) != 0
    corr = np.zeros_like(model)
    snr = amps * psfPeak / noise
    picks = set(np.flatnonzero(snr > sigmaMin).tolist())
    if sat.any():
        for j, (x, y) in enumerate(S["xy"]):
            ix, iy = int(round(x)) - X0, int(round(y)) - Y0
            if 0 <= iy < sat.shape[0] and 0 <= ix < sat.shape[1] and sat[iy, ix] and amps[j] > 0:
                picks.add(j)
    order = sorted(picks, key=lambda j: -amps[j])
    R = radiusFwhm * fwhm
    imgs, info = {}, {}
    for sweep in range(sweeps):
        for j in order:
            c = fit["comps"][j]
            if c is None or amps[j] <= 0:
                continue
            ys, xs = c["ys"], c["xs"]
            Pm = c["stamp"]
            gy, gx = np.gradient(Pm)
            gyy, gyx = np.gradient(gy)
            _, gxx = np.gradient(gx)
            yy, xx = np.mgrid[ys.start:ys.stop, xs.start:xs.stop]
            px, py = S["xy"][j, 0] - X0, S["xy"][j, 1] - Y0
            rr = np.hypot(xx - px, yy - py)
            Rj, rks, isSat = R, (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5), False
            if zones is not None and len(zones):           # saturated star: model its halo
                dz = np.hypot(zones[:, 0] - S["xy"][j, 0], zones[:, 1] - S["xy"][j, 1])
                k = int(np.argmin(dz))
                if dz[k] < 0.5 * zones[k, 2]:
                    Rj = max(R, zones[k, 2])
                    rks = tuple(np.arange(0.5, Rj / fwhm + 0.01, 0.75))
                    isSat = True
            okpix = good[ys, xs] & (var[ys, xs] > 0)
            terms = [Pm, gx, gy, gxx, gyy, gyx]
            if radial:                  # smooth rings: core/halo profile differs from the PSF
                for rk in rks:
                    ring = np.exp(-0.5 * ((rr - rk * fwhm) / (0.5 * fwhm)) ** 2)
                    if (okpix & (ring > 0.5)).sum() >= 10:     # skip rings inside the mask
                        terms.append(Pm.max() * ring)
                        nm = (haloModes if isSat else (brightModes if rk <= 4.0 else 0))
                        if nm:                                 # non-circular structure
                            th = np.arctan2(yy - py, xx - px)
                            for mm in range(1, int(nm) + 1):
                                terms += [Pm.max() * ring * np.cos(mm * th),
                                          Pm.max() * ring * np.sin(mm * th)]
            if coreRings and not isSat:   # narrow core structure (1-2 px scale)
                thc = np.arctan2(yy - py, xx - px)
                for rk in coreRings:
                    ring = np.exp(-0.5 * ((rr - rk) / coreSigma) ** 2)
                    if (okpix & (ring > 0.5)).sum() >= 1:
                        terms.append(Pm.max() * ring)
                        base = ring * (rr / coreSigma if rk == 0 else 1.0)  # spike: its own dipole
                        for mm in range(1, int(coreModes) + 1):
                            terms += [Pm.max() * base * np.cos(mm * thc),
                                      Pm.max() * base * np.sin(mm * thc)]
            T = np.array(terms) * amps[j]
            if j in imgs:
                corr[ys, xs] -= imgs[j]
            sel = okpix & (rr < Rj)
            if sel.sum() < 30:
                continue
            w = 1.0 / np.sqrt(var[ys, xs][sel])
            r = (data[ys, xs] - model[ys, xs] - corr[ys, xs])[sel] * w
            coef, *_ = np.linalg.lstsq(T[:, sel].T * w[:, None], r, rcond=1e-8)
            imgs[j] = np.tensordot(coef, T, axes=1)
            corr[ys, xs] += imgs[j]
            if sweep == sweeps - 1:
                r2 = r - (imgs[j][sel] * w)
                ring = np.hypot(xx - px, yy - py) < R
                info[j] = dict(x=float(S["xy"][j, 0]), y=float(S["xy"][j, 1]), snr=float(snr[j]),
                               nSat=int((sat[ys, xs] & ring).sum()), flux=float(coef[0] * amps[j]),
                               dfluxPct=float(100 * coef[0]),
                               dx=float(-coef[1]), dy=float(-coef[2]),
                               chi2Before=float(np.mean(r ** 2)), chi2After=float(np.mean(r2 ** 2)))
    return corr, info


def printBrightCorrections(info, maxRows=15):
    """Print the table of bright-star corrections.

    Parameters
    ----------
    info : `dict` [`int`, `dict`]
        Second output of `brightStarCorrections`.
    maxRows : `int`, optional
        Number of (brightest) stars printed.
    """
    if not info:
        print("  no bright stars corrected")
        return
    print(f"  bright-star corrections ({len(info)} stars; chi2/px inside the fit radius with "
          f"the pixel variance, ~1 = residual at the noise level):")
    print(f"    {'x':>8s} {'y':>8s} {'peak S/N':>9s} {'SAT px':>6s} {'dx':>6s} {'dy':>6s} "
          f"{'dflux%':>7s} {'chi2 before':>11s} {'after':>7s}")
    for d in sorted(info.values(), key=lambda d: -d["snr"])[:maxRows]:
        print(f"    {d['x']:8.1f} {d['y']:8.1f} {d['snr']:9.0f} {d['nSat']:6d} {d['dx']:+6.2f} "
              f"{d['dy']:+6.2f} {d['dfluxPct']:+7.2f} {d['chi2Before']:11.2f} {d['chi2After']:7.2f}")


def _detect(workingResidual, threshold, minPixels):
    """Detect positive footprints on the working residual (None if detection fails)."""
    cfg = L.SourceDetectionTask.ConfigClass()
    cfg.thresholdValue = threshold
    cfg.thresholdType = "pixel_stdev"
    cfg.minPixels = minPixels
    cfg.doTempLocalBackground = True
    cfg.tempLocalBackground.binSize = 16
    cfg.reEstimateBackground = False
    schema = L.afwTable.SourceTable.makeMinimalSchema()
    task = L.SourceDetectionTask(schema=schema, config=cfg)
    try:
        return task.run(L.afwTable.SourceTable.make(schema), workingResidual).positive.getFootprints()
    except L.TooManyMaskedPixelsError:
        return None


def _slim(fit):
    """Return a fit without its large per-source arrays (for older snapshots)."""
    return {k: v for k, v in fit.items() if k not in ("comps", "chi2map", "groupMap")}


# =====================================================================
# Adaptive inner loop
# =====================================================================
def ringMimics(exposure, S, fit, modelFull, P, fwhm, psfPeak, noise):
    """Find faint non-Gaia sources that are fakes absorbing a bright star's wing excess.

    A faint source near a bright star must stand out from the residual's wing level at its
    distance from the star (including angular modes) and from the same spot rotated by 90, 180
    and 270 degrees (four-fold diffraction structure).

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    S : `dict`
        Source state.
    fit : `dict`
        Joint fit.
    modelFull : `numpy.ndarray`
        Model including the bright-star corrections.
    P : `dict`
        Pipeline options (``ring_*``).
    fwhm : `float`
        PSF FWHM [px].
    psfPeak : `float`
        Peak value of the unit-flux PSF.
    noise : `float`
        Typical pixel noise.

    Returns
    -------
    drop : `list` [`int`]
        Indices of the sources to remove.
    nChecked : `int`
        Number of sources checked.
    """
    data = exposure.image.array
    var, good, amps = fit["variance"], fit["good"], fit["amps"]
    xy = S["xy"]
    bb = exposure.getBBox()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    snr = amps * psfPeak / noise
    bright = np.flatnonzero((snr >= P["ringCheckSnr"]) & (S["ext"] == 0))
    if not len(bright) or len(xy) < 2:
        return [], 0
    tree = cKDTree(xy)
    rmax = P["ringCheckFwhm"] * fwhm
    drop, nChecked = set(), 0
    for b in bright:
        cand = [k for k in tree.query_ball_point(xy[b], rmax)
                if k != b and not S["fromGaia"][k] and S["ext"][k] == 0 and amps[k] > 0
                and amps[k] < P["ringMaxRatio"] * amps[b]
                and np.hypot(*(xy[k] - xy[b])) >= fwhm and fit["comps"][k] is not None]
        if not cand:
            continue
        x, y = xy[b]
        R = rmax + 2 * fwhm
        xa, xb = max(int(x - R) - X0, 0), min(int(x + R) + 2 - X0, data.shape[1])
        ya, yb = max(int(y - R) - Y0, 0), min(int(y + R) + 2 - Y0, data.shape[0])
        sl = (slice(ya, yb), slice(xa, xb))
        base = (data - modelFull)[sl]
        owns = {}
        for k in cand:
            ck = fit["comps"][k]
            own = np.zeros(data.shape)
            own[ck["ys"], ck["xs"]] = amps[k] * ck["stamp"]
            owns[k] = own[sl]
        res = base + sum(owns.values())              # ring level: ALL candidates out of the model
        yy, xx = np.mgrid[ya + Y0:yb + Y0, xa + X0:xb + X0]
        r = np.hypot(xx - x, yy - y)
        g = good[sl] & (var[sl] > 0)
        away = np.ones(r.shape, bool)
        for k in cand:
            away &= np.hypot(xx - xy[k, 0], yy - xy[k, 1]) > 1.0 * fwhm
        rb = np.floor(r).astype(int)
        nb = rb.max() + 1
        th = np.arctan2(yy - y, xx - x)
        M = int(P.get("ringBgModes", 0) or 0)
        nc = 1 + 2 * M
        coef = np.full((nb, nc), np.nan)
        for i in range(nb):
            m = g & (rb == i)
            use = m & away if (m & away).sum() >= 8 else m
            if use.sum() < 8:
                continue
            Mi = M if use.sum() >= 4 * nc else 0          # too few pixels: azimuthal level only
            X = np.column_stack([np.ones(use.sum())] +
                                [f(mm * th[use]) for mm in range(1, Mi + 1) for f in (np.cos, np.sin)])
            yv = res[use]
            w = np.ones(len(yv))
            for _ in range(3):                             # robust: down-weight outliers
                c_, *_ = np.linalg.lstsq(X * w[:, None], yv * w, rcond=None)
                e = yv - X @ c_
                s_ = 1.4826 * np.median(np.abs(e)) + 1e-12
                w = 1.0 / np.maximum(1.0, np.abs(e) / (2.5 * s_))
            coef[i, :] = 0.0
            coef[i, :X.shape[1]] = c_
        okb = np.isfinite(coef[:, 0])
        if okb.sum() < 3:
            continue
        rc = np.flatnonzero(okb) + 0.5
        bg = np.interp(r, rc, coef[okb, 0])
        for q in range(1, nc):
            mm = (q + 1) // 2
            f = np.cos if q % 2 else np.sin
            bg = bg + np.interp(r, rc, coef[okb, q]) * f(mm * th)
        for k in cand:
            nChecked += 1
            ck = fit["comps"][k]
            w = np.zeros(data.shape)
            w[ck["ys"], ck["xs"]] = ck["stamp"]
            w = w[sl]
            ap = g & (np.hypot(xx - xy[k, 0], yy - xy[k, 1]) < 1.5 * fwhm) & (w > 0)
            if ap.sum() < 5:
                continue
            v = var[sl][ap]
            rk = base + owns[k]                      # only this source out; neighbours stay modelled
            den = np.sum(w[ap] ** 2 / v)
            f = np.sum(w[ap] * (rk[ap] - bg[ap]) / v) / den
            snrK = f * np.sqrt(den)
            if P.get("ringSym", False):
                from scipy.ndimage import map_coordinates
                dx, dy = xx[ap] - x, yy[ap] - y
                rot = []
                for c_, s_ in ((0, 1), (-1, 0), (0, -1)):        # 90, 180, 270 deg
                    xr, yr = x + dx * c_ - dy * s_, y + dx * s_ + dy * c_
                    rot.append(map_coordinates(res, [yr - (ya + Y0), xr - (xa + X0)], order=1,
                                               mode="nearest"))
                bsym = np.median(np.array(rot), axis=0)
                fs = np.sum(w[ap] * (rk[ap] - bsym) / v) / den
                snrK = min(snrK, fs * np.sqrt(den))
            if snrK < P["ringSnr"]:
                drop.add(k)
    return sorted(drop), nChecked


def splitTest(exposure, psf, S, fit, modelFull, P, fwhm, psfPeak, noise, badBitmask, fitkw):
    """Split unresolved pairs.

    For each source with ``splitMinSnr <= peak S/N < brightCorrSigma``, a positive residual
    peak (matched-filter S/N >= ``splitPeakSnr``) 0.5-1.5 FWHM away and away from other sources
    gets a trial source.  All trials are fitted together; a trial is kept if the chi2 within
    2 FWHM drops by at least ``splitDchi2`` and its S/N is at least ``minSnr``.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    psf : `lsst.afw.detection.Psf`
        PSF model.
    S : `dict`
        Source state.
    fit : `dict`
        Joint fit.
    modelFull : `numpy.ndarray`
        Model including the bright-star corrections.
    P : `dict`
        Pipeline options (``split_*``).
    fwhm : `float`
        PSF FWHM [px].
    psfPeak : `float`
        Peak value of the unit-flux PSF.
    noise : `float`
        Typical pixel noise.
    badBitmask : `int`
        Mask bits of pixels to ignore.
    fitkw : `dict`
        Keyword arguments of the joint fit.

    Returns
    -------
    S : `dict`
        Source state with the kept second sources.
    fit : `dict`
        Joint fit (the input object if nothing was kept).
    n_tried : `int`
        Number of trial sources.
    nKept : `int`
        Number kept.
    """
    data = exposure.image.array
    var, good, amps = fit["variance"], fit["good"], fit["amps"]
    bb = exposure.getBBox()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    xy = S["xy"]
    resid = data - modelFull
    usable = good & np.isfinite(resid) & (var > 0)
    zmap = ca._matchedSnr(np.where(usable, resid, 0.0), np.where(usable, var, 1.0), usable, fwhm,
                          usable)[0]
    snr = amps * psfPeak / noise
    hosts = np.flatnonzero((snr >= P["splitMinSnr"]) & (snr < P["brightCorrSigma"])
                           & (S["ext"] == 0))
    tree = cKDTree(xy)
    trials, hostOf = [], []
    r1 = 1.5 * fwhm
    for j in hosts:
        x, y = xy[j]
        xa, xb = max(int(x - r1) - X0, 1), min(int(x + r1) + 2 - X0, data.shape[1] - 1)
        ya, yb = max(int(y - r1) - Y0, 1), min(int(y + r1) + 2 - Y0, data.shape[0] - 1)
        if xb <= xa or yb <= ya:
            continue
        yy, xx = np.mgrid[ya + Y0:yb + Y0, xa + X0:xb + X0]
        d = np.hypot(xx - x, yy - y)
        z = np.where((d >= 0.5 * fwhm) & (d <= r1) & usable[ya:yb, xa:xb], zmap[ya:yb, xa:xb], -np.inf)
        k = np.unravel_index(np.argmax(z), z.shape)
        if not np.isfinite(z[k]) or z[k] < P["splitPeakSnr"]:
            continue
        iy, ix = k[0] + ya, k[1] + xa
        w = np.clip(resid[iy - 1:iy + 2, ix - 1:ix + 2], 0, None)
        if w.sum() <= 0:
            continue
        oy, ox = np.mgrid[-1:2, -1:2]
        px = X0 + ix + float(np.sum(w * ox) / w.sum())
        py = Y0 + iy + float(np.sum(w * oy) / w.sum())
        if tree.query((px, py))[0] < 0.7 * fwhm:              # another source already there
            continue
        if trials and np.min(np.hypot(*(np.array(trials) - (px, py)).T)) < 0.7 * fwhm:
            continue
        trials.append((px, py))
        hostOf.append(j)
    if not trials:
        return S, fit, 0, 0
    S2 = _add(S, trials)
    f2 = _fit(exposure, psf, S2, badBitmask, fit["modelFull"], fitkw)
    n0 = len(xy)
    keep = np.ones(len(S2["xy"]), bool)
    yy, xx = np.mgrid[Y0:Y0 + data.shape[0], X0:X0 + data.shape[1]]
    for t, (px, py) in enumerate(trials):
        jn = n0 + t
        x0_, x1_ = int(px - 2 * fwhm) - X0, int(px + 2 * fwhm) + 2 - X0
        y0_, y1_ = int(py - 2 * fwhm) - Y0, int(py + 2 * fwhm) + 2 - Y0
        sl = (slice(max(y0_, 0), max(y1_, 0)), slice(max(x0_, 0), max(x1_, 0)))
        m = good[sl] & (np.hypot(xx[sl] - px, yy[sl] - py) < 2 * fwhm) & (var[sl] > 0)
        c0 = np.sum((data[sl] - fit["modelFull"][sl])[m] ** 2 / var[sl][m])
        c1 = np.sum((data[sl] - f2["modelFull"][sl])[m] ** 2 / f2["variance"][sl][m])
        cj = f2["comps"][jn]
        if cj is None or f2["amps"][jn] <= 0:
            keep[jn] = False
            continue
        g = f2["good"][cj["ys"], cj["xs"]]
        sNew = f2["amps"][jn] * np.sqrt(np.sum(cj["stamp"][g] ** 2 / f2["variance"][cj["ys"], cj["xs"]][g]))
        if (c0 - c1) < P["splitDchi2"] or sNew < P["minSnr"]:
            keep[jn] = False
    nKept = int(keep[n0:].sum())
    if not nKept:
        return S, fit, len(trials), 0
    S3 = _take(S2, keep)
    f3 = _fit(exposure, psf, S3, badBitmask, fit["modelFull"], fitkw)
    return S3, f3, len(trials), nKept


def _ellipModel(psf, x, y, sa, q, th, shape, X0, Y0, prof="gauss"):
    """Return the image of a unit-flux PSF convolved with an elliptical 'gauss' or 'exp' profile."""
    from scipy.signal import fftconvolve
    arr, bb = _rawPsfStamp(psf, x, y, None)
    arr = fftconvolve(arr, _ellipKernel(sa, q, th, prof), mode="same")
    arr = arr / max(arr.sum(), 1e-30)
    out = np.zeros(shape)
    ya, xa = bb.getMinY() - Y0, bb.getMinX() - X0
    y0_, x0_ = max(ya, 0), max(xa, 0)
    y1_, x1_ = min(ya + arr.shape[0], shape[0]), min(xa + arr.shape[1], shape[1])
    if y1_ > y0_ and x1_ > x0_:
        out[y0_:y1_, x0_:x1_] = arr[y0_ - ya:y1_ - ya, x0_ - xa:x1_ - xa]
    return out


def ellipTest(exposure, psf, S, fit, modelFull, P, fwhm, psfPeak, noise, badBitmask, fitkw,
              psfSigma, verboseEllip=True, only=None, explain=False):
    """Model badly fitted sources as small galaxies.

    A badly fitted source (local chi2/px >= ``ellipChi2Min``, S/N between ``ellipMinSnr`` and
    ``brightCorrSigma``) and the points within ``ellipGroupFwhm`` FWHM (at most one Gaia
    source) are refitted with alternative models: 'one' (one extended source replaces all of
    them) and 'star+gal' (the Gaia source stays a point, one extended source replaces the
    others), each with a Gaussian or exponential elliptical profile convolved with the PSF.  The
    best model, penalised by 5 per extra fitted number, must beat the point sources by at least
    ``ellipDchi2``.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    psf : `lsst.afw.detection.Psf`
        PSF model.
    S : `dict`
        Source state.
    fit : `dict`
        Joint fit.
    modelFull : `numpy.ndarray`
        Model including the bright-star corrections.
    P : `dict`
        Pipeline options (``ellip_*``).
    fwhm : `float`
        PSF FWHM [px].
    psfPeak : `float`
        Peak value of the unit-flux PSF.
    noise : `float`
        Typical pixel noise.
    badBitmask : `int`
        Mask bits of pixels to ignore.
    fitkw : `dict`
        Keyword arguments of the joint fit.
    psfSigma : `float`
        PSF sigma [px].
    verboseEllip : `bool`, optional
        Print every model tried.
    only : `list` [`int`], optional
        Only try these source indices.
    explain : `bool`, optional
        Print why sources are not tried, and tracebacks of failed fits.

    Returns
    -------
    S : `dict`
        Source state with the extended sources.
    fit : `dict`
        Joint fit.
    n_tried : `int`
        Number of sources tried.
    nKept : `int`
        Number of extended sources kept.
    """
    from scipy.optimize import least_squares
    data = exposure.image.array
    var, good, amps = fit["variance"], fit["good"], fit["amps"]
    bb = exposure.getBBox()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    xy = S["xy"]
    resid = data - modelFull
    snr = amps * psfPeak / noise
    tree = cKDTree(xy)
    yy, xx = np.mgrid[Y0:Y0 + data.shape[0], X0:X0 + data.shape[1]]
    hosts = np.flatnonzero((snr >= P["ellipMinSnr"]) & (snr < P["brightCorrSigma"])
                           & (S["ext"] == 0) & (amps > 0))
    if only is not None:
        for k in only:
            if explain and k not in hosts:
                print(f"    source {k} at ({xy[k, 0]:.1f}, {xy[k, 1]:.1f}) is not a candidate: peak S/N "
                      f"{snr[k]:.1f} (needs {P['ellipMinSnr']:g}..{P['brightCorrSigma']:g}), "
                      f"extended={S['ext'][k] > 0}, flux {amps[k]:.3g}")
        hosts = np.array([k for k in only if k in hosts], int)
    used, accepted, nTry = set(), [], 0
    for k0 in hosts[np.argsort(-amps[hosts])]:
        j = int(k0)
        if j in used:
            continue
        x, y = xy[j]
        m2 = good & (np.hypot(xx - x, yy - y) < 2 * fwhm) & (var > 0)
        c2Loc = np.mean(resid[m2] ** 2 / var[m2]) if m2.sum() else np.nan
        if m2.sum() < 20 or c2Loc < P["ellipChi2Min"]:
            if explain:
                print(f"    ({x:.1f}, {y:.1f}): local chi2/px {c2Loc:.2f} < {P['ellipChi2Min']:g} "
                      f"(fit variance incl. sysFrac) -> not tried")
            continue
        near = [k for k in tree.query_ball_point(xy[j], P["ellipGroupFwhm"] * fwhm)
                if k != j and k not in used and S["ext"][k] == 0 and amps[k] > 0
                and fit["comps"][k] is not None]
        gaiaNear = [k for k in near if S["fromGaia"][k]]
        if (S["fromGaia"][j] and gaiaNear) or len(gaiaNear) > 1:
            if explain:
                print(f"    ({x:.1f}, {y:.1f}): another Gaia source within {P['ellipGroupFwhm']:g} "
                      f"FWHM -> not merged (two Gaia stars)")
            continue
        gsrc = j if S["fromGaia"][j] else (gaiaNear[0] if gaiaNear else None)
        grp = list(dict.fromkeys(([gsrc] if gsrc is not None else []) +
                                 [k for k in [j] + near if not S["fromGaia"][k]]))
        if any(fit["comps"][k] is None for k in grp):
            continue
        stamps = {}
        own = np.zeros(data.shape)
        for k in grp:
            ck = fit["comps"][k]
            st = np.zeros(data.shape)
            st[ck["ys"], ck["xs"]] = ck["stamp"]
            stamps[k] = st
            own += amps[k] * st
        T = resid + own                                # the group's light alone
        f_ = amps[grp]
        cx, cy = np.sum(f_ * xy[grp, 0]) / f_.sum(), np.sum(f_ * xy[grp, 1]) / f_.sum()
        reg = good & (np.hypot(xx - cx, yy - cy) < 3 * fwhm) & (var > 0)
        if reg.sum() < 30:
            continue
        w = 1.0 / var[reg]
        Tr = T[reg]
        chiPt = float(np.sum((Tr - own[reg]) ** 2 * w))
        pPt = 3 * len(grp)
        nTry += 1
        others = [k for k in grp if k != gsrc]
        variants = [("one", None, grp)]
        if gsrc is not None and others:
            variants.append(("star+gal", gsrc, others))
        best = None
        for vname, keepPt, repl in variants:
            fr = amps[repl]
            gx = np.sum(fr * xy[repl, 0]) / fr.sum()
            gy = np.sum(fr * xy[repl, 1]) / fr.sum()
            fixed = stamps[keepPt][reg] if keepPt is not None else None
            for prof in ("gauss", "exp"):
                def comps_(p):
                    m = _ellipModel(psf, p[0], p[1], np.exp(p[2]), 0.1 + 0.9 / (1 + np.exp(-p[3])),
                                    p[4], data.shape, X0, Y0, prof)[reg]
                    A = [m] if fixed is None else [m, fixed]
                    A = np.array(A)
                    M = (A * w) @ A.T
                    try:
                        fl = np.linalg.solve(M + 1e-12 * np.eye(len(A)), (A * w) @ Tr)
                    except np.linalg.LinAlgError:
                        fl = np.zeros(len(A))
                    return A, fl

                def resFn(p):
                    A, fl = comps_(p)
                    return (Tr - fl @ A) * np.sqrt(w)

                s0 = np.log(0.7 * psfSigma if prof == "gauss" else 0.5 * psfSigma)
                try:
                    sol = least_squares(resFn, [gx, gy, s0, 1.0, 0.0],
                                        x_scale=[0.5, 0.5, 0.3, 0.5, 0.3],
                                        bounds=([gx - fwhm, gy - fwhm, np.log(0.15), -6, -10],
                                                [gx + fwhm, gy + fwhm, np.log(3 * fwhm), 6, 10]),
                                        max_nfev=200)
                except Exception as e:
                    print(f"    elliptical fit ({vname}, {prof}) at ({gx:.1f}, {gy:.1f}) FAILED: "
                          f"{type(e).__name__}: {e}")
                    if explain:
                        traceback.print_exc()
                    continue
                A, fl = comps_(sol.x)
                chi = float(np.sum((Tr - fl @ A) ** 2 * w))
                sa, q, th = float(np.exp(sol.x[2])), float(0.1 + 0.9 / (1 + np.exp(-sol.x[3]))), float(sol.x[4])
                pM = 6 + (3 if keepPt is not None else 0)
                need = P["ellipDchi2"] + 5.0 * max(0, pM - pPt)
                ok_ = (np.all(fl > 0) and chiPt - chi >= need and sa >= 0.3 * psfSigma
                       and q >= 0.15)
                if verboseEllip:
                    print(f"    try {vname:8s} {prof:5s} at ({sol.x[0]:.1f}, {sol.x[1]:.1f}): "
                          f"{len(grp)} point(s) chi2 {chiPt:.0f} -> {chi:.0f} (needs -{need:.0f}); "
                          f"scale {sa:.2f} px, q {q:.2f}, PA {np.degrees(th) % 180:.0f} deg"
                          f"{'' if np.all(fl > 0) else ', negative flux'} -> {'ok' if ok_ else 'no'}")
                score = chi - chiPt + need
                if ok_ and (best is None or score < best[0]):
                    best = (score, vname, keepPt, repl, sol.x[0], sol.x[1], sa, q, th, prof, chi)
        if best is not None:
            _, vname, keepPt, repl, ex, ey, sa, q, th, prof, chi = best
            accepted.append((vname, keepPt, repl, ex, ey, sa, q, th, prof, chiPt, chi))
            used.update(grp)
    if not accepted:
        return S, fit, nTry, 0
    S = dict(S, xy=S["xy"].copy(), ext=S["ext"].copy(), ref=S["ref"].copy())
    keep = np.ones(len(S["xy"]), bool)
    for vname, keepPt, repl, ex, ey, sa, q, th, prof, c0, c1 in accepted:
        host = repl[0] if keepPt is not None else (next((k for k in repl if S["fromGaia"][k]), repl[0]))
        S["xy"][host] = (ex, ey)
        if S["fromGaia"][host]:            # a galaxy's Gaia centroid is unreliable: no drift limit
            S["ref"][host] = (ex, ey)
        S["ext"][host] = sa
        _ELLIP[int(S["uid"][host])] = (float(q), float(th), prof)
        for k in repl:
            if k != host:
                keep[k] = False
        print(f"    extended source ({vname}, {prof}) at ({ex:.1f}, {ey:.1f}) replaces {len(repl)} point "
              f"source(s){' next to the Gaia star' if keepPt is not None else ''}: scale {sa:.2f} px, "
              f"axis ratio {q:.2f}, PA {np.degrees(th) % 180:.0f} deg; local chi2 {c0:.0f} -> {c1:.0f}")
    S = _take(S, keep)
    _clampToGaia(S, P["gaiaMaxDriftPx"])
    f2 = _fit(exposure, psf, S, badBitmask, fit["modelFull"], fitkw)
    return S, f2, nTry, len(accepted)


def tryElliptical(run, xy, chi2Min=None, dchi2=None):
    """Try the elliptical-source test on one source of a finished run, printing every step.

    Shows why the source is or is not tried, the fitted models, the chi2 change and any error.
    The run is not changed.

    Parameters
    ----------
    run : `lsst.meas.crowd.session.CutoutRun`
        A finished run.
    xy : `tuple` [`float`, `float`]
        Position; the nearest source is tested.
    chi2Min : `float`, optional
        Override the trigger (``ellipChi2Min``).
    dchi2 : `float`, optional
        Override the acceptance bar (``ellipDchi2``).

    Returns
    -------
    result : `dict` or `None`
        ``kept``, ``S`` and ``fit`` (with the extended source if it was kept); None if the run's
        last fit has no per-source stamps.
    """
    if not _ORIG or _HUGE is None:
        enableSpeedups()
    final = run.outerHistory[-1]
    exp, res = final["exposure"], final["result"]
    S, fit = res["S"], res["snaps"][-1]["fit"]
    if fit.get("comps") is None:
        print("this run's last fit has no per-source stamps -- cannot test")
        return None
    P = dict(_params(run.config))
    P["poissonK"] = float(final.get("poissonK", 0.0) or 0.0)
    if chi2Min is not None:
        P["ellipChi2Min"] = chi2Min
    if dchi2 is not None:
        P["ellipDchi2"] = dchi2
    psf = exp.getPsf()
    c = exp.getBBox().getCenter()
    fwhm = L._fwhm(psf, c)
    psfSigma = fwhm / 2.355
    psfPeak = psf.computeKernelImage(c).array.max()
    noise = float(np.sqrt(np.nanmedian(exp.variance.array)))
    j = int(cKDTree(S["xy"]).query(xy)[1])
    print(f"nearest source: index {j} at ({S['xy'][j, 0]:.1f}, {S['xy'][j, 1]:.1f}), "
          f"{'Gaia' if S['fromGaia'][j] else 'non-Gaia'}, flux {fit['amps'][j]:.3g}, "
          f"peak S/N {fit['amps'][j] * psfPeak / noise:.1f}; FWHM {fwhm:.2f} px")
    S = {k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in S.items()}
    fitkw = dict(supportFrac=P["supportFrac"], sysFrac=P["sysFrac"], poissonK=P["poissonK"])
    nBefore = len(_ELLIP)
    S2, f2, nTry, nKeep = ellipTest(exp, psf, S, fit, res["modelFull"], P, fwhm, psfPeak, noise,
                                    run.badBitmask, fitkw, psfSigma, only=[j], explain=True)
    if nKeep:                       # do not leave the what-if in the global table
        for u in list(_ELLIP)[nBefore:]:
            _ELLIP.pop(u, None)
    print(f"tried {nTry}, kept {nKeep}")
    return dict(kept=bool(nKeep), S=S2, fit=f2)


def finalPositionPolish(exposure, psf, S, fit, badBitmask, P, fitkw, fwhm, psfSigma):
    """Re-optimise all source positions jointly.

    Gauss-Newton steps (capped at ``finalPolishCap`` px, Gaia drift limit kept) and refit,
    repeated while the total chi2 drops.  Fixes dipoles from positions that went stale when
    neighbours were added or removed.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    psf : `lsst.afw.detection.Psf`
        PSF model.
    S : `dict`
        Source state.
    fit : `dict`
        Joint fit of ``S``.
    badBitmask : `int`
        Mask bits of pixels to ignore.
    P : `dict`
        Pipeline options (``final_polish_*``).
    fitkw : `dict`
        Keyword arguments of the joint fit.
    fwhm : `float`
        PSF FWHM [px].
    psfSigma : `float`
        PSF sigma [px].

    Returns
    -------
    S : `dict`
        Source state with the polished positions.
    fit : `dict`
        Joint fit (the input object if nothing improved).
    info : `dict`
        ``rounds``, ``moved``, ``chi0`` and ``chi1`` (total chi2 before / after).
    """
    data = exposure.image.array

    def chiAll(f):
        g = f["good"]
        return float(np.sum((data - f["modelFull"])[g] ** 2 / f["variance"][g]))

    S = dict(S)
    c0 = cPrev = chiAll(fit)
    best = fit
    allidx = np.arange(len(S["xy"]))
    moved, nR = 0, 0
    for r in range(1, P["finalPolishRounds"] + 1):
        before = S["xy"].copy()
        S["xy"] = _refineSubset(exposure, psf, S, best, allidx, P["finalPolishCap"],
                                P["centroidMinAmp"], P["freezeFwhm"] * fwhm, psfSigma)
        _clampToGaia(S, P["gaiaMaxDriftPx"])
        fNew = _fit(exposure, psf, S, badBitmask, best["modelFull"], fitkw)
        cNew = chiAll(fNew)
        if cNew >= cPrev:
            S["xy"] = before
            break
        d = np.hypot(*(S["xy"] - before).T)
        moved = max(moved, int((d > 0.01).sum()))
        best, cPrev, nR = fNew, cNew, r
        if d.max() < P["centroidShiftTol"]:
            break
    return S, best, dict(rounds=nR, moved=moved, chi0=c0, chi1=cPrev)


def runInnerAdaptive(exposure, psf, badBitmask, S, P, outerIdx=0, thresholdStart=None,
                     shapeNow=False, timer=None):
    """Run one outer pass of the adaptive pipeline.

    A warm fit of the current sources (Gaia starting positions or the previous pass), then detect -> add ->
    joint fit -> refine -> clean, repeated until the residual check says stop; then the halo,
    ring, position-polish, pair-split and galaxy tests.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Background-subtracted cutout.
    psf : `lsst.afw.detection.Psf`
        PSF model.
    badBitmask : `int`
        Mask bits of pixels to ignore.
    S : `dict`
        Starting source state.
    P : `dict`
        Pipeline options.
    outerIdx : `int`, optional
        Outer-pass index.
    thresholdStart : `float`, optional
        First detection threshold [sigma].
    shapeNow : `bool`, optional
        Apply the shape cut in this pass.
    timer : `lsst.meas.crowd.fitLib.StepTimer`, optional
        Timer collecting the step durations.

    Returns
    -------
    result : `dict`
        Compatible with `fitLib.runInnerLoopVerbose`'s result (``modelFull``,
        ``residualFull``, ``peaks``, ``amps``, histories, ``snaps``), plus ``S``,
        ``brightCorr``, ``adaptive`` (the loop decisions), ``isGaia`` and ``gaiaIndex``.
    """
    timer = timer or L.TIMER
    bbox = exposure.getBBox()
    ny, nx = bbox.getHeight(), bbox.getWidth()
    X0, Y0 = bbox.getMinX(), bbox.getMinY()
    data = exposure.image.array
    badMask = (exposure.mask.array & badBitmask) != 0
    sh = psf.computeShape(bbox.getCenter())
    psfSigma = math.sqrt(0.5 * (sh.getIxx() + sh.getIyy()))
    fwhm = 2.355 * psfSigma
    psfPeak = psf.computeKernelImage(bbox.getCenter()).array.max()
    noise = float(np.sqrt(np.nanmedian(exposure.variance.array)))
    fitkw = dict(supportFrac=P["supportFrac"], sysFrac=P["sysFrac"], poissonK=P["poissonK"])
    target = float(P["targetSnr"] if P["targetSnr"] is not None else P["minSnr"])
    tFloor, thrMin, step = P["tFloor"], P["thresholdMin"], P["thresholdStep"]
    thr = float(thresholdStart if thresholdStart is not None else P["thresholdStart"])
    edgePx = (P["edgeMarginFwhm"] or 0.0) * fwhm
    padPx = P["gaiaPadFwhm"] * fwhm
    mergeSep = (P["mergeSeparationSigma"] or 2.0) * psfSigma
    sched = list(P["centroidShiftSchedule"] or [1.5, 0.8, 0.4, 0.2, 0.1])
    sched += [sched[-1]] * max(0, P["nCentroidRounds"] - len(sched))
    ts = sorted({3.0, 3.5, 4.0, 4.5, 5.0, 6.0, target, tFloor})

    coreLab = (satCoreLabels(exposure, fwhm, P["satZoneMinPx"]) if P["oneSourcePerSat"]
               else np.zeros((ny, nx), np.int32))
    okpix = ~badMask & np.isfinite(data)

    def rstd(img):                  # residual std over usable pixels only
        return float(np.nanstd(img[okpix]))

    zones = P.get("_satZones", np.zeros((0, 3)))
    zoneMask = P.get("_satZoneMask")
    if zoneMask is None or zoneMask.shape != (ny, nx):
        zoneMask = np.zeros((ny, nx), bool)
    if P["maskPadFwhm"] and badMask.any():           # flux of sources centred in the mask
        from scipy.ndimage import distance_transform_edt
        zoneMask = zoneMask | (distance_transform_edt(~badMask) <= P["maskPadFwhm"] * fwhm)

    def check(fit, model):
        resid = data - model
        return ca.residualExcess(resid, fit["variance"], fit["good"] & ~zoneMask, S["xy"],
                                 fit["amps"],
                                 X0, Y0, fwhm, psfPeak, noise, thresholds=ts,
                                 edgeFwhm=P["edgeFwhm"], brightSigma=P["brightPeakSigma"],
                                 haloFwhm=P["haloFwhm"], minSepFwhm=P["minSepFwhm"])

    fit, modelFull = None, np.zeros((ny, nx))
    hist = dict(residStd=[rstd(data)], nPeaks=[0], nZero=[0])
    snaps, firstSeen, log = [], {}, []
    bcorr = {}
    stampPt = bbox.getCenter()

    def corrected(fit):
        """Return the model with the bright-star correction terms (updates the stamp sizes)."""
        nonlocal bcorr
        setStampSizes(S["xy"], fit["amps"], psf, stampPt, noise, fwhm, P["stampFloorSigma"])
        _setBfFluxes(S["xy"], fit["amps"])
        if not P["brightCorr"]:
            bcorr = {}
            return fit["modelFull"]
        corr, bcorr = brightStarCorrections(exposure, S, fit, fit["modelFull"], psfPeak,
                                            noise, fwhm, badBitmask, P["brightCorrSigma"],
                                            P["brightCorrRadiusFwhm"],
                                            radial=P["brightCorrRadial"], zones=zones,
                                            haloModes=P["satHaloModes"],
                                            brightModes=P["brightAngularModes"],
                                            coreRings=P["brightCoreRingsPx"],
                                            coreSigma=P["brightCoreRingSigma"],
                                            coreModes=P["brightCoreModes"])
        return fit["modelFull"] + corr
    prevAmpUid = {}
    diag, stop, noProg = None, "maxInner", 0
    prevN, prevE = None, None
    rejected = np.asarray(P.get("_haloRejects", np.zeros((0, 2))), float).reshape(-1, 2).copy()
    shapeJudgedOnce = False

    if len(S["xy"]):
        with timer.step("warm-start fit", outer=outerIdx):
            fit = _fit(exposure, psf, S, badBitmask, None, fitkw)
            fit = _fit(exposure, psf, S, badBitmask, fit["modelFull"], fitkw)
        modelFull = corrected(fit)
        prevAmpUid = dict(zip(S["uid"].tolist(), fit["amps"]))
        hist = dict(residStd=[rstd(data - modelFull)], nPeaks=[len(S["xy"])],
                    nZero=[fit["nAtZero"]])
        print(f"  warm start: {len(S['xy'])} sources ({int(S['fromGaia'].sum())} Gaia), residual "
              f"std {hist['residStd'][0]:.3f}")

    for it in range(P["maxInner"]):
        times = {}

        def T(label):
            return timer.step(label, outer=outerIdx, inner=it)

        print(f"\n--- outer {outerIdx + 1} | inner {it + 1} (cap {P['maxInner']}) | "
              f"threshold={thr:.2f}σ | sources={len(S['xy'])} ---")
        detIn = data - P["detectGain"] * modelFull
        wr = exposure.clone()
        wr.image.array[:] = detIn
        if fit is not None:
            wr.variance.array[:] = fit["variance"]
        wr.setPsf(psf)
        with T("detect on residual") as rec:
            fps = _detect(wr, thr, P["minPixels"])
        times["detect"] = rec["seconds"]
        if fps is None:
            stop = "detection failed (background estimation)"
            break
        rawPeaks = [(pk.getFx(), pk.getFy()) for fp in fps for pk in fp.getPeaks()]

        # ---- filter the new candidates
        with T("filter new peaks") as rec:
            new, _ = L.collapseSaturatedPeaks(rawPeaks, exposure)
            if len(S["xy"]) and len(new) and fit is not None:
                isBright = fit["amps"] * psfPeak > P["brightPeakSigma"] * noise
                new, nRej = L.rejectNearBright(new, _xy(S), isBright, P["rejectFwhm"] * fwhm)
            if P["coarseCentroid"] and new:
                hw = max(int(round(P["coarseCentroidHwSigma"] * psfSigma)), 3)
                new = [L.windowedCentroid(detIn, (X0, Y0), x, y, hw) for x, y in new]
            new, _ = L.dropMaskedPeaks(new, exposure)
            new = [(x, y) for x, y in new
                   if X0 + edgePx <= x <= bbox.getMaxX() - edgePx
                   and Y0 + edgePx <= y <= bbox.getMaxY() - edgePx]
            if len(new) > 1:
                new = L.mergeClosePeaks(new, exposure, minSeparation=mergeSep)
            if new and len(S["xy"]):
                dd, _ = cKDTree(S["xy"]).query(np.asarray(new))
                new = [p for p, di in zip(new, dd) if di >= mergeSep]
            nRem = 0
            if new and coreLab.any() and len(S["xy"]):     # core already has its source
                taken = set(_labelAt(S["xy"], coreLab, X0, Y0).tolist()) - {0}
                lc = _labelAt(new, coreLab, X0, Y0)
                new = [p for p, l in zip(new, lc) if l not in taken]
            if new and P["rememberRejects"] and len(rejected):
                dd, _ = cKDTree(rejected).query(np.asarray(new))
                nRem = int((dd < mergeSep).sum())
                new = [p for p, di in zip(new, dd) if di >= mergeSep]
        times["filter"] = rec["seconds"]
        print(f"  {len(rawPeaks)} raw peaks -> {len(new)} new candidate(s)"
              + (f" ({nRem} skipped: already rejected this pass)" if nRem else ""))

        nAdded = 0
        if (new or snaps == []) and (len(S["xy"]) + len(new)) > 0:
            firstNewUid = S["nextUid"]
            S = _add(S, new)
            with T("joint fit") as rec:
                fit = _fit(exposure, psf, S, badBitmask, modelFull if fit is not None else None,
                           fitkw)
            times["fit"] = rec["seconds"]

            if P["refineCentroids"]:
                first = not snaps
                nRounds = P["nCentroidRounds"] if first else P["centroidRoundsAfterFirst"]
                if first:
                    idx = np.arange(len(S["xy"]))
                else:                          # new sources + their neighbours only
                    isnew = S["uid"] >= firstNewUid
                    idx = np.arange(len(S["xy"]))
                    if isnew.any():
                        dd, _ = cKDTree(S["xy"][isnew]).query(S["xy"])
                        idx = np.flatnonzero(dd < 2.0 * fwhm)
                    else:
                        idx = idx[:0]
                r, shift = -1, 0.0
                with T("refine centroids") as rec:
                    for r in range(nRounds if len(idx) else 0):
                        before = S["xy"].copy()
                        cap = sched[r] if first else P["laterShiftCap"]
                        S["xy"] = _refineSubset(exposure, psf, S, fit, idx, cap,
                                                P["centroidMinAmp"], P["freezeFwhm"] * fwhm,
                                                psfSigma)
                        _clampToGaia(S, P["gaiaMaxDriftPx"])
                        fit = _fit(exposure, psf, S, badBitmask, fit["modelFull"], fitkw)
                        shift = float(np.max(np.hypot(*(S["xy"] - before).T))) if len(before) else 0.0
                        if shift < P["centroidShiftTol"]:
                            break
                if P["satRefineWings"] and coreLab.any():
                    S["xy"], nMv = _refineSaturated(exposure, psf, S, fit, coreLab, X0, Y0, fwhm)
                    if nMv:
                        fit = _fit(exposure, psf, S, badBitmask, fit["modelFull"], fitkw)
                        print(f"    saturated stars: {nMv} position(s) refined from the wings")
                times["centroids"] = rec["seconds"]
                print(f"    centroids: {len(idx)} source(s), {r + 1} round(s), last max shift "
                      f"{shift:.3f}px")

            with T("prune + merge + bounds") as rec:
                bright = fit["amps"] * psfPeak > P["brightPeakSigma"] * noise
                alive, S["xy"] = _pruneMerge(fit, S, P["minSnr"], P["pruneSepSigma"] * psfSigma,
                                             bright, P["nearBrightFwhm"] * fwhm,
                                             P["nearBrightSnr"], P["protectGaia"],
                                             inZone=_zoneFlags(S["xy"], zoneMask, X0, Y0),
                                             zoneSnr=P["satZoneSnr"])
                x, y = S["xy"][:, 0], S["xy"][:, 1]
                m = np.where(S["fromGaia"], padPx, -edgePx)        # Gaia sources may sit just outside
                inside = ((x >= X0 - m) & (x <= bbox.getMaxX() + m) &
                          (y >= Y0 - m) & (y <= bbox.getMaxY() + m))
                nEdge = int((alive & ~inside).sum())
                S = _take(S, alive & inside)
                if coreLab.any() and len(S["xy"]):
                    ampsNow = fit["amps"][np.flatnonzero(alive & inside)]
                    one = _onePerCore(S, ampsNow, coreLab, X0, Y0)
                    if not one.all():
                        print(f"  saturated cores: removed {int((~one).sum())} extra source(s) "
                              f"({int((~one & S['fromGaia']).sum())} of them Gaia entries)")
                        S = _take(S, one)
                if nEdge:
                    print(f"  dropped {nEdge} source(s) outside / within the edge margin")
                fit = _fit(exposure, psf, S, badBitmask, fit["modelFull"], fitkw) if len(S["xy"]) else None
            times["prune"] = rec["seconds"]
            if fit is None:
                stop = "everything pruned"
                break

            if P["useEliminateWeak"]:
                with T("eliminate weak") as rec:
                    S, fit = _eliminateWeak(exposure, psf, S, fit, badBitmask, P, fitkw, fwhm)
                times["eliminate weak"] = rec["seconds"]
            shapeRejected = []
            if shapeNow:
                with T("shape cut") as rec:
                    ju = None
                    if shapeJudgedOnce:          # only new sources and their neighbours
                        isnew = S["uid"] >= firstNewUid
                        ju = set()
                        if isnew.any():
                            dd, _ = cKDTree(S["xy"][isnew]).query(S["xy"])
                            ju = set(S["uid"][dd < 2.0 * fwhm].tolist())
                    S, fit, shapeRejected = _shapeCut(exposure, psf, S, fit, badBitmask, P, fitkw,
                                                      judgeUids=ju, psfSigma=psfSigma)
                    shapeJudgedOnce = True
                times["shape cut"] = rec["seconds"]

            modelFull = corrected(fit)
            resid = data - modelFull
            nAdded = int((S["uid"] >= firstNewUid).sum())
            if new:
                kept = set(S["uid"].tolist())
                lost = [p for k, p in enumerate(new) if firstNewUid + k not in kept]
                if lost:
                    rejected = np.vstack([rejected, np.asarray(lost, float).reshape(-1, 2)])
            isNew = S["uid"] >= firstNewUid
            prevAmp = np.array([prevAmpUid.get(u, np.nan) for u in S["uid"].tolist()])
            prevAmpUid = dict(zip(S["uid"].tolist(), fit["amps"]))
            for x, y in S["xy"]:
                firstSeen.setdefault((round(x, 1), round(y, 1)), it)
            oldStd, newStd = hist["residStd"][-1], rstd(resid)
            hist["residStd"].append(newStd)
            hist["nPeaks"].append(len(S["xy"]))
            hist["nZero"].append(fit["nAtZero"])
            if snaps and not P["keepAllFits"]:
                snaps[-1]["fit"] = _slim(snaps[-1]["fit"])
            snaps.append(dict(it=it, threshold=thr, detIn=detIn, rawPeaks=rawPeaks,
                              collapsedPeaks=new, mergedPeaks=_xy(S), mergedArr=S["xy"].copy(),
                              isNew=isNew, prevAmp=prevAmp,
                              footprintMask=L._footprintsToMask(fps, exposure), fit=fit,
                              residual=resid.copy(), stdBefore=oldStd, stdAfter=newStd,
                              times=times, badMask=badMask, shapeRejected=shapeRejected,
                              isGaia=S["fromGaia"].copy()))
            with T("residual check"):
                diag = check(fit, modelFull)
            print(f"  {len(S['xy'])} sources ({int(S['fromGaia'].sum())} Gaia, "
                  f"{int((S['ext'] > 0).sum())} extended, {nAdded} new kept), "
                  f"residual std {oldStd:.3f} -> {newStd:.3f}; time "
                  + ", ".join(f"{k}={v:.1f}s" for k, v in times.items()))
        elif not len(S["xy"]):
            print("  no sources yet")

        # ---- adaptive decision ---------------------------------------------------
        N = len(S["xy"])
        if diag is not None:
            missing, e, err = ca.isMissing(diag, target, P["zSig"], P["minExcess"])
            p_, n_ = ca.countsAt(diag, target)
        else:
            missing, e, err, p_, n_ = True, 0, 0.0, 0, 0
        growing = nAdded >= max(P["growMin"], P["growFrac"] * N)
        progress = nAdded > 0 and (prevN is None or N > prevN or (prevE is not None and e < prevE))
        atFloor = thr <= thrMin + 1e-9
        thrOld = thr
        print(f"  residual check at >={target:g}σ: +{p_}/-{n_} peaks, excess {e:+d} ± {err:.1f} "
              f"-> {'MISSING sources' if missing else 'consistent with noise'}")
        if growing and not atFloor:
            action = f"still finding sources ({nAdded} added) -> hold threshold {thr:.2f}σ"
        elif not atFloor:
            thr = max(thr + step, thrMin)
            action = f"lower threshold {thrOld:.2f} -> {thr:.2f}σ"
        elif growing and noProg < P["patience"]:
            action = f"at threshold floor, still adding ({nAdded}) -> continue"
        elif missing and thr > tFloor + 1e-9 and diag is not None:
            tNew = ca.pickThreshold(diag, P["fdr"], tFloor, thr)
            if tNew < thr - 1e-9:
                thr = tNew
                action = f"excess remains -> go deeper: threshold {thrOld:.2f} -> {thr:.2f}σ"
            else:
                stop, action = "excess remains, cannot go deeper (false-positive limit)", "STOP"
        else:
            stop = "residual clean" if not missing else "stalled (no progress)"
            action = "STOP"
        noProg = 0 if (progress or thr != thrOld) else noProg + 1
        if action != "STOP" and atFloor and noProg >= P["patience"] and thr == thrOld:
            stop, action = f"no progress for {noProg} iterations", "STOP"
        log.append(dict(it=it, threshold=thrOld, n=N, added=nAdded, pos=p_, neg=n_, excess=e,
                        missing=missing, action=action))
        print(f"  => {action}" + (f" ({stop})" if action == "STOP" else ""))
        prevN, prevE = N, e
        if action == "STOP":
            break
    else:
        print(f"  WARNING: hit maxInner={P['maxInner']} while still working -- raise maxInner")

    if P["dropHaloMimics"] and fit is not None and snaps and len(S["xy"]):
        import contextlib
        import io
        import types
        fake = types.SimpleNamespace(outerHistory=[dict(
            exposure=exposure, result=dict(S=S, snaps=[dict(fit=fit)], modelFull=modelFull))])
        with contextlib.redirect_stdout(io.StringIO()):
            tb = checkBrightNeighbours(fake, brightSnr=P["haloCheckSnr"],
                                       withinFwhm=P["haloCheckFwhm"], plot=False,
                                       verbose=False, similarMin=P["haloSimilarMin"])
        badIds = ([int(r["id"]) for r in tb if r["verdict"] == "doubtful"] if len(tb) else [])
        if badIds:
            keep = np.ones(len(S["xy"]), bool)
            keep[badIds] = False
            P["_haloRejects"] = np.vstack([np.asarray(P.get("_haloRejects", np.zeros((0, 2)))
                                                      ).reshape(-1, 2), S["xy"][~keep]])
            nTot = len(tb)
            S = _take(S, keep)
            fit = _fit(exposure, psf, S, badBitmask, fit["modelFull"], fitkw)
            modelFull = corrected(fit)
            diag = check(fit, modelFull)
            last = snaps[-1]
            last.update(fit=fit, mergedArr=S["xy"].copy(), mergedPeaks=_xy(S),
                        residual=(data - modelFull).copy(), isGaia=S["fromGaia"].copy(),
                        isNew=np.asarray(last["isNew"])[keep] if len(last["isNew"]) == len(keep)
                        else np.zeros(len(S["xy"]), bool),
                        prevAmp=np.asarray(last["prevAmp"])[keep] if len(last["prevAmp"]) == len(keep)
                        else np.full(len(S["xy"]), np.nan))
            print(f"  halo check: removed {len(badIds)} of {nTot} non-Gaia source(s) next to bright "
                  f"stars as halo mimics (repeat at >= {P['haloSimilarMin']} control angles or "
                  f"significance < 3); they will not be re-detected")
    if P["dropRingMimics"] and fit is not None and snaps and len(S["xy"]):
        badIds, nChk = ringMimics(exposure, S, fit, modelFull, P, fwhm, psfPeak, noise)
        if badIds:
            keep = np.ones(len(S["xy"]), bool)
            keep[badIds] = False
            P["_haloRejects"] = np.vstack([np.asarray(P.get("_haloRejects", np.zeros((0, 2)))
                                                      ).reshape(-1, 2), S["xy"][~keep]])
            S = _take(S, keep)
            fit = _fit(exposure, psf, S, badBitmask, fit["modelFull"], fitkw)
            modelFull = corrected(fit)
            diag = check(fit, modelFull)
            last = snaps[-1]
            last.update(fit=fit, mergedArr=S["xy"].copy(), mergedPeaks=_xy(S),
                        residual=(data - modelFull).copy(), isGaia=S["fromGaia"].copy(),
                        isNew=np.asarray(last["isNew"])[keep] if len(last["isNew"]) == len(keep)
                        else np.zeros(len(S["xy"]), bool),
                        prevAmp=np.asarray(last["prevAmp"])[keep] if len(last["prevAmp"]) == len(keep)
                        else np.full(len(S["xy"]), np.nan))
        print(f"  ring check: removed {len(badIds)} of {nChk} faint non-Gaia source(s) around bright "
              f"stars that do not stand out from the star's wing at their radius (S/N < "
              f"{P['ringSnr']:g}); they will not be re-detected")
    def _polish(label):
        nonlocal S, fit, modelFull, diag
        S, fNew, info = finalPositionPolish(exposure, psf, S, fit, badBitmask, P, fitkw, fwhm,
                                            psfSigma)
        if fNew is not fit:
            fit = fNew
            modelFull = corrected(fit)
            diag = check(fit, modelFull)
            snaps[-1].update(fit=fit, mergedArr=S["xy"].copy(), mergedPeaks=_xy(S),
                             residual=(data - modelFull).copy(), isGaia=S["fromGaia"].copy())
        print(f"  {label}: {info['rounds']} round(s), up to {info['moved']} source(s) "
              f"moved, total chi2 {info['chi0']:.4g} -> {info['chi1']:.4g} "
              f"({100 * (info['chi0'] - info['chi1']) / max(info['chi0'], 1e-12):+.2f}%)")

    if P["finalPolish"] and fit is not None and snaps and len(S["xy"]):
        _polish("final position polish")
    if P["splitTest"] and fit is not None and snaps and len(S["xy"]):
        S_new, fNew, nTry, nKeep = splitTest(exposure, psf, S, fit, modelFull, P, fwhm,
                                             psfPeak, noise, badBitmask, fitkw)
        print(f"  split test: {nTry} source(s) with a lopsided residual tried as pairs, "
              f"{nKeep} second source(s) kept")
        if nKeep:
            S, fit = S_new, fNew
            modelFull = corrected(fit)
            diag = check(fit, modelFull)
            snaps[-1].update(fit=fit, mergedArr=S["xy"].copy(), mergedPeaks=_xy(S),
                             residual=(data - modelFull).copy(), isGaia=S["fromGaia"].copy(),
                             isNew=np.zeros(len(S["xy"]), bool),
                             prevAmp=np.full(len(S["xy"]), np.nan))
            if P["finalPolish"]:
                _polish("position polish after split")
    if P["ellipTest"] and fit is not None and snaps and len(S["xy"]):
        S_new, fNew, nTry, nKeep = ellipTest(exposure, psf, S, fit, modelFull, P, fwhm,
                                             psfPeak, noise, badBitmask, fitkw, psfSigma)
        print(f"  elliptical test: {nTry} badly fitted source(s) tried as elliptical extended "
              f"sources, {nKeep} kept")
        if nKeep:
            S, fit = S_new, fNew
            modelFull = corrected(fit)
            diag = check(fit, modelFull)
            snaps[-1].update(fit=fit, mergedArr=S["xy"].copy(), mergedPeaks=_xy(S),
                             residual=(data - modelFull).copy(), isGaia=S["fromGaia"].copy(),
                             isNew=np.zeros(len(S["xy"]), bool),
                             prevAmp=np.full(len(S["xy"]), np.nan))
    finalMissing = bool(ca.isMissing(diag, target, P["zSig"], P["minExcess"])[0]) if diag else True
    ampsOut = np.array([]) if fit is None else fit["amps"].copy()
    for j, d in bcorr.items():
        if j < len(ampsOut):
            ampsOut[j] += d["flux"]
    if P["brightCorr"]:
        printBrightCorrections(bcorr)
    return dict(modelFull=modelFull, residualFull=data - modelFull, peaks=_xy(S),
                amps=ampsOut, brightCorr=bcorr,
                residStdHistory=hist["residStd"], nPeaksHistory=hist["nPeaks"],
                nZeroAmpHistory=hist["nZero"], firstSeen=firstSeen, snaps=snaps,
                S=S, isGaia=S["fromGaia"].copy(), gaiaIndex=S["gaia"].copy(),
                adaptive=dict(final=diag, target=target, missing=finalMissing, stop=stop,
                              log=log, nInner=len(log)))


# =====================================================================
# Adaptive outer loop
# =====================================================================
def selectPsfStars(exposure, res, fwhm, psfPeak, noise, P, verbose=True):
    """Select clean, bright, isolated, unsaturated point sources for the PSF stack.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    res : `dict`
        Result of `runInnerAdaptive`.
    fwhm : `float`
        PSF FWHM [px].
    psfPeak : `float`
        Peak value of the unit-flux PSF.
    noise : `float`
        Typical pixel noise.
    P : `dict`
        Pipeline options (``psf_stack_*``).
    verbose : `bool`, optional
        Print the selection.

    Returns
    -------
    idx : `numpy.ndarray` [`int`]
        Indices of the selected stars, brightest first.
    """
    S = res["S"]
    fit = res["snaps"][-1]["fit"]
    amps = fit["amps"]
    xy = S["xy"]
    bb = exposure.getBBox()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    snr = amps * psfPeak / noise
    ok = (snr >= P["psfStackMinSnr"]) & (S["ext"] == 0)
    sat = _satlike(exposure)
    if sat.any():
        from scipy.ndimage import binary_dilation
        nearSat = binary_dilation(sat, iterations=int(np.ceil(3 * fwhm)))
        ix = np.clip(np.round(xy[:, 0]).astype(int) - X0, 0, sat.shape[1] - 1)
        iy = np.clip(np.round(xy[:, 1]).astype(int) - Y0, 0, sat.shape[0] - 1)
        ok &= ~nearSat[iy, ix]
    h = 0.5 * fwhm * 4                                   # keep the stamp core inside the image
    ok &= ((xy[:, 0] > X0 + h) & (xy[:, 0] < bb.getMaxX() - h) &
           (xy[:, 1] > Y0 + h) & (xy[:, 1] < bb.getMaxY() - h))
    if len(xy) > 1:
        tree = cKDTree(xy)
        for j in np.flatnonzero(ok):
            nb = [k for k in tree.query_ball_point(xy[j], P["psfStackIsoFwhm"] * fwhm) if k != j]
            if any(amps[k] > 0.1 * amps[j] for k in nb):
                ok[j] = False
    rej = np.asarray(P.get("_haloRejects", np.zeros((0, 2))), float).reshape(-1, 2)
    if len(rej):                    # their light is back in the residual -> would bias the stack
        dd, _ = cKDTree(rej).query(xy)
        ok &= dd > 3.0 * fwhm
    idx = np.flatnonzero(ok)
    idx = idx[np.argsort(-amps[idx])][:P["psfStackN"]]
    if verbose:
        print(f"  PSF stack: {len(idx)} clean star(s) (of {int((snr >= P['psfStackMinSnr']).sum())} "
              f"with peak S/N >= {P['psfStackMinSnr']:g}; saturated / crowded / extended excluded)")
    return idx


def _starMask(shape, X0, Y0, pos, r):
    """Return the pixels within ``r`` px of any of the positions."""
    yy, xx = np.mgrid[Y0:Y0 + shape[0], X0:X0 + shape[1]]
    m = np.zeros(shape, bool)
    for x, y in pos:
        m |= (xx - x) ** 2 + (yy - y) ** 2 < r * r
    return m


def _chi2On(exposure, f, m):
    """Return the mean chi2 per pixel of a fit over a pixel mask."""
    m = m & f["good"]
    return float(np.mean((exposure.image.array - f["modelFull"])[m] ** 2
                         / exposure.variance.array[m])) if m.any() else np.nan


def refinePsfClean(exposure, psf, res, badBitmask, P, fitkw, nCycles=3, verbose=True):
    """Update the PSF from the residual stack of clean stars, cross-validated.

    The clean stars are split into two halves A and B (alternating in brightness).  Each cycle
    builds a correction from A and judges it on B's pixels, and vice versa; only if the held-out
    chi2 improves by at least ``psfCvMinGain`` is the correction from all stars applied.  This
    stops the PSF absorbing the pixel noise of its own stars, which shows up as speckle around
    every bright star.  Corrections are smoothed (``psfSmoothPx``) and soft-thresholded
    (``psfSoftSigma``).

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    psf : `lsst.afw.detection.Psf`
        Current PSF.
    res : `dict`
        Result of `runInnerAdaptive`.
    badBitmask : `int`
        Mask bits of pixels to ignore.
    P : `dict`
        Pipeline options.
    fitkw : `dict`
        Keyword arguments of the joint fit.
    nCycles : `int`, optional
        Maximum number of correction cycles.
    verbose : `bool`, optional
        Print each cycle.

    Returns
    -------
    psf : `lsst.afw.detection.Psf`
        The updated PSF (the input PSF if no cycle was kept).
    """
    from scipy.ndimage import gaussian_filter
    bb = exposure.getBBox()
    c = bb.getCenter()
    fwhm = L._fwhm(psf, c)
    psfPeak = psf.computeKernelImage(c).array.max()
    noise = float(np.sqrt(np.nanmedian(exposure.variance.array)))
    idx = selectPsfStars(exposure, res, fwhm, psfPeak, noise, P, verbose)
    if len(idx) < 6:
        print(f"  PSF stack: {len(idx)} clean star(s) -- too few to cross-validate, PSF not updated")
        return psf
    S = res["S"]
    fit = res["snaps"][-1]["fit"]
    data = exposure.image.array
    X0, Y0 = bb.getMinX(), bb.getMinY()
    A, B = idx[0::2], idx[1::2]
    mA = _starMask(data.shape, X0, Y0, S["xy"][A], 2.5 * fwhm)
    mB = _starMask(data.shape, X0, Y0, S["xy"][B], 2.5 * fwhm)
    skw = dict(P["stackKwargs"] or {})
    skw.setdefault("softSigma", P["psfSoftSigma"])
    sc = 0.5 * max(bb.getWidth(), bb.getHeight())

    def correctedKernel(kernel, sel, f):
        new, _ = L.estimatePsfCorrection(data, f["modelFull"], f["variance"], f["good"],
                                         (X0, Y0), S["xy"][sel], f["amps"][sel], kernel,
                                         **dict(skw, nUse=len(sel)))
        d = gaussian_filter(new - kernel, P["psfSmoothPx"]) if P["psfSmoothPx"] else new - kernel
        out = kernel + d
        return out / out.sum()

    def asPsf(k, like):
        prev = _SPATIAL.get(id(like))
        if prev is not None and prev[0] is like and prev[1].shape == k.shape:
            return spatialKernelPsf(k, prev[1], prev[2], c.getX(), c.getY(), sc)
        return L.kernelToLsstPsf(k)

    cur = psf
    for cyc in range(max(1, nCycles)):
        kernel = cur.computeKernelImage(c).array.copy()
        try:
            fA = _fit(exposure, asPsf(correctedKernel(kernel, A, fit), cur), S, badBitmask,
                      fit["modelFull"], fitkw)
            fB = _fit(exposure, asPsf(correctedKernel(kernel, B, fit), cur), S, badBitmask,
                      fit["modelFull"], fitkw)
        except ValueError as e:
            print(f"    PSF cycle {cyc + 1}: {e}")
            break
        before = _chi2On(exposure, fit, mB) + _chi2On(exposure, fit, mA)
        after = _chi2On(exposure, fA, mB) + _chi2On(exposure, fB, mA)
        gain = (before - after) / max(before, 1e-12)
        keep = gain >= P["psfCvMinGain"]
        if verbose:
            print(f"    PSF cycle {cyc + 1}/{nCycles}: held-out chi2/px on {len(idx)} PSF stars "
                  f"{before / 2:.2f} -> {after / 2:.2f} ({100 * gain:+.1f}%) "
                  f"({'kept' if keep else 'rejected, stopping'})")
        if not keep:
            break
        cand = asPsf(correctedKernel(kernel, idx, fit), cur)
        fit = _fit(exposure, cand, S, badBitmask, fit["modelFull"], fitkw)
        cur = cand
    return cur


def _residualStack(resid, good, xy, amps, idx, X0, Y0, half, withErr=False):
    """Return the median residual/flux of stars on a grid centred on each (and optionally its error).

    Sampled at each star's sub-pixel position on a (2*half+1)^2 grid; NaN where fewer than three
    stars contribute.
    """
    from scipy.ndimage import map_coordinates
    u = np.arange(-half, half + 1)
    V, U = np.meshgrid(u, u, indexing="ij")
    g = good.astype(float)
    st = []
    for j in idx:
        x, y = xy[j]
        co = [(y - Y0 + V).ravel(), (x - X0 + U).ravel()]
        r = map_coordinates(resid, co, order=3, mode="nearest").reshape(V.shape)
        ok = map_coordinates(g, co, order=1, mode="nearest").reshape(V.shape) > 0.99
        st.append(np.where(ok, r / amps[j], np.nan))
    if not st:
        return None
    st = np.array(st)
    med = np.nanmedian(st, axis=0)
    nn = np.sum(np.isfinite(st), axis=0)
    med = np.where(nn >= 3, med, np.nan)
    if not withErr:
        return med
    # error of the median: 1.25 * robust sigma / sqrt(n)
    err = 1.2533 * 1.4826 * np.nanmedian(np.abs(st - med), axis=0) / np.sqrt(np.maximum(nn, 1))
    return med, err


def refinePsfMedian(exposure, psf, res, badBitmask, P, fitkw, fit0=None, verbose=True):
    """Correct the PSF core with the median residual/flux of the bright clean stars.

    The median is robust to neighbours, cosmic rays and per-star noise, so it keeps the shape
    error shared by the bright stars.  Applied within ``psfCoreRadiusFwhm`` FWHM (cosine
    taper), soft-thresholded, smoothed and cross-validated on two halves of the stars.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    psf : `lsst.afw.detection.Psf`
        Current PSF.
    res : `dict`
        Result of `runInnerAdaptive`.
    badBitmask : `int`
        Mask bits of pixels to ignore.
    P : `dict`
        Pipeline options (``psf_core_*``).
    fitkw : `dict`
        Keyword arguments of the joint fit.
    fit0 : `dict`, optional
        Joint fit with ``psf`` (computed if not given).
    verbose : `bool`, optional
        Print each cycle.

    Returns
    -------
    psf : `lsst.afw.detection.Psf`
        The corrected PSF (the input PSF if no correction was kept).
    """
    from scipy.ndimage import gaussian_filter
    bb = exposure.getBBox()
    c = bb.getCenter()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    kernel = psf.computeKernelImage(c).array.astype(float).copy()
    fwhm = L._fwhm(psf, c)
    noise = float(np.sqrt(np.nanmedian(exposure.variance.array)))
    S = res["S"]
    fit = fit0 if fit0 is not None else _fit(exposure, psf, S, badBitmask,
                                             res["snaps"][-1]["fit"]["modelFull"], fitkw)
    Pq = dict(P, psfStackMinSnr=P["psfCoreMinSnr"], psfStackN=60)
    idx = selectPsfStars(exposure, dict(res, snaps=[dict(fit=fit)]), fwhm, kernel.max(),
                         noise, Pq, verbose=False)
    if len(idx) < 6:
        if verbose:
            print(f"  PSF core: {len(idx)} bright clean star(s) (S/N >= {P['psfCoreMinSnr']:g}) "
                  f"-- too few, not corrected")
        return psf
    rmax = P["psfCoreRadiusFwhm"] * fwhm
    half = int(np.ceil(rmax)) + 1
    H = kernel.shape[0] // 2
    half = min(half, H)
    u = np.arange(-half, half + 1)
    R = np.hypot(*np.meshgrid(u, u))
    taper = np.where(R <= 0.6 * rmax, 1.0,
                     np.where(R >= rmax, 0.0, 0.5 * (1 + np.cos(np.pi * (R - 0.6 * rmax) / (0.4 * rmax)))))
    data = exposure.image.array
    A_, B_ = idx[0::2], idx[1::2]
    mA = _starMask(data.shape, X0, Y0, S["xy"][A_], rmax)
    mB = _starMask(data.shape, X0, Y0, S["xy"][B_], rmax)

    def corrected(f, sel, kbase):
        d, e = _residualStack(data - f["modelFull"], f["good"], S["xy"], f["amps"], sel, X0, Y0,
                              half, withErr=True)
        d, e = np.nan_to_num(d), np.nan_to_num(e, nan=np.inf)
        # keep only significant pixels (soft threshold): the outer, noise-dominated part of the
        # median punched holes into the kernel
        d = np.sign(d) * np.clip(np.abs(d) - P["psfSoftSigma"] * e, 0.0, None) * taper
        if P["psfSmoothPx"]:
            d = gaussian_filter(d, P["psfSmoothPx"])
        k = kbase.copy()
        sub = k[H - half:H + half + 1, H - half:H + half + 1]
        # never remove more than 70% of the PSF at a pixel (no holes)
        k[H - half:H + half + 1, H - half:H + half + 1] = np.maximum(sub + d, 0.3 * sub)
        if P.get("psfRegularize"):
            k = smoothPsfKernel(k, fwhm, P["psfRegCoreFwhm"], P["psfRegSigma"])
        return k / k.sum()

    cur = psf
    for cyc in range(max(1, P["psfCoreCycles"])):
        try:
            fA = _fit(exposure, L.kernelToLsstPsf(corrected(fit, A_, kernel)), S, badBitmask,
                      fit["modelFull"], fitkw)
            fB = _fit(exposure, L.kernelToLsstPsf(corrected(fit, B_, kernel)), S, badBitmask,
                      fit["modelFull"], fitkw)
        except ValueError as e:
            print(f"    PSF core cycle {cyc + 1}: {e}")
            break
        before = _chi2On(exposure, fit, mB) + _chi2On(exposure, fit, mA)
        after = _chi2On(exposure, fA, mB) + _chi2On(exposure, fB, mA)
        gain = (before - after) / max(before, 1e-12)
        keep = gain >= P["psfCvMinGain"]
        if verbose:
            print(f"    PSF core cycle {cyc + 1}: median of {len(idx)} bright stars, held-out "
                  f"chi2/px {before / 2:.2f} -> {after / 2:.2f} ({100 * gain:+.1f}%) "
                  f"({'kept' if keep else 'rejected'})")
        if not keep:
            break
        kernel = corrected(fit, idx, kernel)
        cur = L.kernelToLsstPsf(kernel)
        fit = _fit(exposure, cur, S, badBitmask, fit["modelFull"], fitkw)
    return cur


def smoothPsfKernel(kernel, fwhm, coreFwhm=1.5, sigma=1.0, maxPattern=1.5):
    """Return a hole-free, smooth version of a PSF kernel.

    Inside ``coreFwhm`` FWHM the kernel is kept.  Outside it is rebuilt as
    ``radial profile(r) * exp(angular pattern)``: the profile is the median of each 1-px ring
    (so a few holed pixels cannot pull it down), forced to decrease outwards; the pattern is
    ``log(kernel / profile)`` median-filtered over 5x5 px (removes isolated holes and hot pixels),
    smoothed by ``sigma`` px and clipped to +-``maxPattern``.  Smooth structure (ellipticity,
    lobes, spikes wider than ~2 px) survives; pixel-scale defects do not.

    Parameters
    ----------
    kernel : `numpy.ndarray`
        PSF kernel (odd size, centred).
    fwhm : `float`
        PSF FWHM [px].
    coreFwhm : `float`, optional
        Radius kept unchanged [FWHM].
    sigma : `float`, optional
        Smoothing of the angular pattern [px].
    maxPattern : `float`, optional
        Clip of the log angular pattern.

    Returns
    -------
    kernel : `numpy.ndarray`
        Smoothed, non-negative kernel summing to 1.
    """
    from scipy.ndimage import median_filter, gaussian_filter
    k = np.asarray(kernel, float)
    H = k.shape[0] // 2
    u = np.arange(-H, H + 1)
    R = np.hypot(*np.meshgrid(u, u))
    floor = max(k.max(), 1e-30) * 1e-9
    rb = np.round(R).astype(int)
    nr = rb.max() + 1
    prof = np.array([np.median(k[rb == i]) if np.any(rb == i) else np.nan for i in range(nr)])
    r0 = coreFwhm * fwhm
    i0 = int(np.ceil(r0))
    for i in range(i0 + 1, nr):                    # monotone non-increasing outside the core
        if np.isfinite(prof[i - 1]):
            prof[i] = min(prof[i], prof[i - 1]) if np.isfinite(prof[i]) else prof[i - 1]
    prof = np.maximum(np.nan_to_num(prof, nan=floor), floor)
    base = np.exp(np.interp(R, np.arange(nr), np.log(prof)))
    pat = np.log(np.maximum(k, floor) / base)
    pat = gaussian_filter(median_filter(pat, size=5, mode="nearest"), sigma, mode="nearest")
    smooth = base * np.exp(np.clip(pat, -maxPattern, maxPattern))
    w = np.clip((r0 + 2.0 - R) / 2.0, 0.0, 1.0)
    out = w * k + (1.0 - w) * smooth
    out = np.clip(out, 0.0, None)
    return out / out.sum()


def regularizePsf(psf, point, P, verbose=True, label=""):
    """Return the PSF with its kernel passed through `smoothPsfKernel`.

    Parameters
    ----------
    psf : `lsst.afw.detection.Psf`
        PSF.
    point : `lsst.geom.Point2D`
        Position at which the kernel is taken.
    P : `dict`
        Pipeline options (``psfRegularize``, ``psf_reg_*``); no-op if regularisation is off.
    verbose : `bool`, optional
        Print how much flux moved.
    label : `str`, optional
        Suffix for the printed line.

    Returns
    -------
    psf : `lsst.afw.detection.Psf`
        The regularised (constant) PSF, or the input PSF.
    """
    if not P.get("psfRegularize"):
        return psf
    k = psf.computeKernelImage(point).array.astype(float)
    fwhm = L._fwhm(psf, point)
    new = smoothPsfKernel(k, fwhm, P["psfRegCoreFwhm"], P["psfRegSigma"])
    if verbose:
        H = k.shape[0] // 2
        u = np.arange(-H, H + 1)
        out = np.hypot(*np.meshgrid(u, u)) > P["psfRegCoreFwhm"] * fwhm
        ks = k / k.sum()
        print(f"  PSF regularised{label}: {100 * np.abs(new - ks)[out].sum():.2f}% of the flux "
              f"moved outside {P['psfRegCoreFwhm'] * fwhm:.1f}px (holes / pixel noise)")
    return L.kernelToLsstPsf(new)


def _radialProfile(img, ok, x, y, X0, Y0, edges):
    """Return the median of an image in annuli around (x, y); NaN where fewer than 5 good pixels."""
    R = edges[-1]
    ix0, iy0 = int(round(x)) - X0, int(round(y)) - Y0
    xa, xb = max(ix0 - int(R) - 1, 0), min(ix0 + int(R) + 2, img.shape[1])
    ya, yb = max(iy0 - int(R) - 1, 0), min(iy0 + int(R) + 2, img.shape[0])
    if xb <= xa or yb <= ya:
        return np.full(len(edges) - 1, np.nan)
    yy, xx = np.mgrid[Y0 + ya:Y0 + yb, X0 + xa:X0 + xb]
    r = np.hypot(xx - x, yy - y)
    v, m = img[ya:yb, xa:xb], ok[ya:yb, xa:xb]
    out = np.full(len(edges) - 1, np.nan)
    idx = np.digitize(r, edges) - 1
    for b in range(len(edges) - 1):
        sel = m & (idx == b)
        if sel.sum() >= 5:
            out[b] = np.median(v[sel])
    return out


def _azimuthalMean(k):
    """Return the image of a kernel's azimuthal mean (1-px annuli about the centre)."""
    h = k.shape[0] // 2
    u = np.arange(-h, h + 1)
    rb = np.round(np.hypot(*np.meshgrid(u, u))).astype(int)
    prof = np.bincount(rb.ravel(), k.ravel()) / np.maximum(np.bincount(rb.ravel()), 1)
    return prof[rb]


def _angularPattern(k, support):
    """Return a kernel's angular modulation, kernel / ring average (1 beyond ``support`` px)."""
    h = k.shape[0] // 2
    u = np.arange(-h, h + 1)
    R = np.hypot(*np.meshgrid(u, u))
    az = _azimuthalMean(k)
    pat = np.where(az > 1e-12 * max(k.max(), 1e-30), k / np.where(az > 0, az, 1.0), 1.0)
    pat = np.where(R <= support - 1, pat, 1.0)
    return np.clip(pat, 0.0, 50.0)


def extendPsfWingsKeep(psf, point, size=81, blend=(12.0, 16.0), keepAngular=False, verbose=True):
    """Extend a PSF with a Moffat tail, optionally keeping its angular structure.

    Like `fitLib.extendPsfWings`, but with ``keepAngular`` only the kernel's azimuthal mean is
    replaced by the Moffat beyond the blend radius, so spikes and other non-circular structure
    are kept.

    Parameters
    ----------
    psf : `lsst.afw.detection.Psf`
        PSF.
    point : `lsst.geom.Point2D`
        Position at which the kernel is taken.
    size : `int`, optional
        Size of the extended kernel [px].
    blend : `tuple` [`float`, `float`], optional
        Radial range of the cross-fade to the Moffat [px].
    keepAngular : `bool`, optional
        Keep the kernel's angular pattern in the tail.
    verbose : `bool`, optional
        Print the fitted tail.

    Returns
    -------
    psf : `lsst.afw.detection.Psf`
        The extended constant PSF (the input PSF if the fit fails).
    """
    kernel = psf.computeKernelImage(point).array.astype(float).copy()
    try:
        _, info = L.extendKernelWithMoffat(kernel, newSize=size)
    except Exception as e:
        print(f"  PSF wing extension failed ({type(e).__name__}: {e}) -- keeping current PSF")
        return psf
    A, a, b = info["A"], info["a"], info["beta"]
    H, h = size // 2, kernel.shape[0] // 2
    emb = np.zeros((2 * H + 1, 2 * H + 1))
    hh = min(h, H)
    emb[H - hh:H + hh + 1, H - hh:H + hh + 1] = kernel[h - hh:h + hh + 1, h - hh:h + hh + 1]
    v = np.arange(-H, H + 1)
    R = np.hypot(*np.meshgrid(v, v))
    b1 = min(blend[1], hh - 1.0)
    b0 = min(blend[0], b1 - 2.0)
    w = np.clip((b1 - R) / (b1 - b0), 0.0, 1.0)
    pat = _angularPattern(emb, hh) if keepAngular else 1.0
    new = w * emb + (1.0 - w) * A * (1.0 + (R / a) ** 2) ** (-b) * pat
    new /= new.sum()
    if verbose:
        print(f"  PSF wings extended{' (angular pattern kept)' if keepAngular else ''}: "
              f"{kernel.shape[0]}px -> {size}px; "
              f"Moffat a={a:.2f}px, beta={b:.2f}")
    return L.kernelToLsstPsf(new)


def wingProfilePsf(exposure, res, psf, P, verbose=True):
    """Return a PSF whose outer wings are measured from the stars.

    Neighbour-subtracted radial profiles of unsaturated bright stars (peak S/N >=
    ``wingMinSnr``), divided by their flux, give the normalised profile where they have signal;
    saturated stars are scaled to it where both are measured and extend it outwards.  A Moffat
    fitted to the combined profile replaces the kernel beyond ``wingBlendFwhm``.  Degenerate fits
    are rejected.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    res : `dict`
        Result of `runInnerAdaptive`.
    psf : `lsst.afw.detection.Psf`
        Current PSF.
    P : `dict`
        Pipeline options (``wing_*``, ``psfWingSize``).
    verbose : `bool`, optional
        Print the fitted profile.

    Returns
    -------
    psf : `lsst.afw.detection.Psf` or `None`
        The new PSF, or None if there are too few stars or the fit is degenerate.
    """
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)      # all-NaN annuli are expected
        return _wingProfilePsf(exposure, res, psf, P, verbose)


def _wingProfilePsf(exposure, res, psf, P, verbose=True):
    """Implement `wingProfilePsf` (warnings about empty annuli are silenced by the caller)."""
    from scipy.optimize import curve_fit
    bb = exposure.getBBox()
    c = bb.getCenter()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    kernel = psf.computeKernelImage(c).array.astype(float).copy()
    size = int(P["psfWingSize"] or kernel.shape[0]) | 1
    if kernel.shape[0] < size:                         # embed a smaller kernel
        big = np.zeros((size, size))
        o = (size - kernel.shape[0]) // 2
        big[o:o + kernel.shape[0], o:o + kernel.shape[0]] = kernel
        kernel = big
    H = kernel.shape[0] // 2
    fwhm = L._fwhm(psf, c)
    fit = res["snaps"][-1]["fit"]
    S = res["S"]
    amps, model, good = fit["amps"], fit["modelFull"], fit["good"]
    data = exposure.image.array
    psfPeak = kernel.max()
    noise = float(np.sqrt(np.nanmedian(exposure.variance.array)))
    coreLab = satCoreLabels(exposure, fwhm, P["satZoneMinPx"])
    core = _labelAt(S["xy"], coreLab, X0, Y0)
    snr = amps * psfPeak / noise
    edges = np.arange(0.0, H + 0.01, 1.0)
    rc = 0.5 * (edges[1:] + edges[:-1])
    U, Sat = [], []
    for j in np.flatnonzero(((snr >= P["wingMinSnr"]) | (core > 0)) & (S["ext"] == 0)):
        cj = fit["comps"][j]
        if cj is None or amps[j] <= 0:
            continue
        own = np.zeros_like(data)
        own[cj["ys"], cj["xs"]] = amps[j] * cj["stamp"]
        prof = _radialProfile(data - model + own, good, S["xy"][j, 0], S["xy"][j, 1], X0, Y0, edges)
        if core[j] > 0:
            Sat.append(prof)
        else:
            U.append(prof / amps[j])
    if len(U) < 2:
        if verbose:
            print(f"  wings from stars: only {len(U)} bright unsaturated star(s) -- skipped")
        return None
    U = np.array(U)
    u = np.nanmedian(U, axis=0)
    nu = np.sum(np.isfinite(U), axis=0)
    eu = 1.4826 * np.nanmedian(np.abs(U - u), axis=0) / np.sqrt(np.maximum(nu, 1))
    reliable = np.isfinite(u) & (u > 3 * eu) & (nu >= 2)
    comb = np.where(reliable, u, np.nan)
    nSatUsed = 0
    if Sat:
        profS = []
        for sp in Sat:
            m = reliable & np.isfinite(sp) & (rc >= 1.5 * fwhm)
            if m.sum() < 3:
                continue
            scale = np.sum(sp[m] * u[m]) / np.sum(u[m] ** 2)
            if scale > 0:
                profS.append(sp / scale)
        if profS:
            nSatUsed = len(profS)
            ps = np.nanmedian(np.array(profS), axis=0)
            beyond = ~reliable & np.isfinite(ps) & (ps > 0) & (rc > 1.5 * fwhm)
            comb = np.where(beyond, ps, comb)
    m = np.isfinite(comb) & (comb > 0) & (rc >= 1.5 * fwhm)
    if m.sum() < 4:
        if verbose:
            print("  wings from stars: profile not measurable -- skipped")
        return None
    moffat = lambda r, A, a, b: A * (1.0 + (r / a) ** 2) ** (-b)
    try:
        (A, a, b), _ = curve_fit(moffat, rc[m], comb[m], p0=[comb[m][0] * 2, fwhm, 2.5],
                                 sigma=comb[m], bounds=([0, 0.3, 1.05], [np.inf, 80.0, 12.0]),
                                 maxfev=10000)
    except Exception as e:
        if verbose:
            print(f"  wings from stars: Moffat fit failed ({e}) -- skipped")
        return None
    if a < 0.5 * fwhm or a > 79.0 or b < 1.5 or b > 11.9:
        if verbose:
            print(f"  wings from stars: Moffat fit degenerate (a={a:.2f}px, beta={b:.2f}; at a limit "
                  f"or unphysical) -- wings NOT changed")
        return None
    v = np.arange(-H, H + 1)
    Rg = np.hypot(*np.meshgrid(v, v))
    r0, r1 = P["wingBlendFwhm"][0] * fwhm, P["wingBlendFwhm"][1] * fwhm
    w = np.clip((r1 - Rg) / (r1 - r0), 0.0, 1.0)
    # the Moffat replaces only the RADIAL profile; the kernel's angular pattern (spikes, lobes)
    # relative to its ring average is kept
    pat = _angularPattern(kernel, H) if P["wingKeepAngular"] else 1.0
    new = w * kernel + (1.0 - w) * moffat(Rg, A, a, b) * pat
    new /= new.sum()
    tail = float(new[Rg > r1].sum())
    if verbose:
        print(f"  wings from stars: {len(U)} bright + {nSatUsed} saturated star(s); profile "
              f"measured to {rc[m].max():.0f}px; Moffat a={a:.2f}px beta={b:.2f}; "
              f"{100 * tail:.2f}% of flux beyond {r1:.0f}px (was "
              f"{100 * float(kernel[Rg > r1].sum() / kernel.sum()):.2f}%)")
    return L.kernelToLsstPsf(new)


_SPATIAL = {}      # id(psf) -> (psf, Kx, Ky): gradient terms of a spatially varying PSF


def spatialKernelPsf(K0, Kx, Ky, xc, yc, scale):
    """Build a linearly varying LSST PSF.

    ``K(x, y) = K0 + (x - xc)/scale * Kx + (y - yc)/scale * Ky``, as a LinearCombinationKernel
    with linear spatial functions (the same machinery as the PCA PSF).

    Parameters
    ----------
    K0 : `numpy.ndarray`
        Kernel at (xc, yc).
    Kx, Ky : `numpy.ndarray`
        Kernel change per ``scale`` px in x and y.
    xc, yc : `float`
        Reference position.
    scale : `float`
        Length scale [px].

    Returns
    -------
    psf : `lsst.meas.algorithms.KernelPsf`
        The spatially varying PSF.
    """
    import lsst.afw.math as afwMath
    import lsst.afw.image as afwImage
    kernels = [afwMath.FixedKernel(afwImage.makeImageFromArray(np.ascontiguousarray(k, dtype=np.float64)))
               for k in (K0, Kx, Ky)]
    funcs = []
    for prm in ([1.0, 0.0, 0.0], [-xc / scale, 1.0 / scale, 0.0], [-yc / scale, 0.0, 1.0 / scale]):
        f = afwMath.PolynomialFunction2D(1)
        f.setParameters(prm)
        funcs.append(f)
    try:
        lck = afwMath.LinearCombinationKernel(kernels, funcs)
    except TypeError:
        kl, fl = afwMath.KernelList(), []
        for k in kernels:
            kl.append(k)
        lck = afwMath.LinearCombinationKernel(kl, funcs)
    psf = L.KernelPsf(lck, L.geom.Point2D(xc, yc))
    _SPATIAL[id(psf)] = (psf, Kx, Ky)
    return psf


def _psfGradient(data, model, var, good, xy0, pos, amps, kernel, xc, yc, scale,
                 maxRadius=20, softSigma=1.0):
    """Fit the stars' normalised residual stamps per pixel as a + bx*dx + by*dy; return (bx, by, n)."""
    from scipy.ndimage import map_coordinates
    h = kernel.shape[0] // 2
    u = np.arange(-h, h + 1)
    V, U = np.meshgrid(u, u, indexing="ij")
    ny, nx = data.shape
    resid = data - model
    gf = good.astype(float)
    R, W, D = [], [], []
    for (x, y), a in zip(pos, amps):
        if a <= 0:
            continue
        yy, xx = (y - xy0[1]) + V, (x - xy0[0]) + U
        co = [yy.ravel(), xx.ravel()]
        inside = (yy >= 0) & (yy <= ny - 1) & (xx >= 0) & (xx <= nx - 1)
        r = map_coordinates(resid, co, order=3, mode="nearest").reshape(V.shape) / a
        gm = map_coordinates(gf, co, order=1, mode="nearest").reshape(V.shape) > 0.99
        v = map_coordinates(var, co, order=1, mode="nearest").reshape(V.shape)
        R.append(r)
        W.append(np.where(inside & gm & (v > 0), a * a / np.maximum(v, 1e-30), 0.0))
        D.append([1.0, (x - xc) / scale, (y - yc) / scale])
    R, W, D = np.array(R), np.array(W), np.array(D)
    if len(R) < 3:
        return None, None, len(R)
    M = np.einsum("jp,ja,jb->pab", W.reshape(len(R), -1), D, D)      # (npix, 3, 3)
    rhs = np.einsum("jp,jp,ja->pa", W.reshape(len(R), -1), R.reshape(len(R), -1), D)
    M += np.eye(3)[None] * 1e-12 * np.maximum(np.trace(M, axis1=1, axis2=2), 1e-30)[:, None, None]
    sol = np.linalg.solve(M, rhs[..., None])[..., 0]
    err = np.sqrt(np.clip(np.linalg.inv(M)[:, [1, 2], [1, 2]], 0, None))
    out = []
    rr = np.hypot(U, V)
    he = min(h, maxRadius)
    win = np.clip((he - rr) / max(0.4 * he, 1e-9), 0, 1)
    win = 0.5 - 0.5 * np.cos(np.pi * win)
    gy, gx = np.gradient(kernel)
    A = np.column_stack([gx.ravel(), gy.ravel()])
    for k in (1, 2):
        b = sol[:, k]
        e = err[:, k - 1]
        b = np.sign(b) * np.maximum(np.abs(b) - softSigma * e, 0.0)
        b = (b.reshape(kernel.shape)) * win
        coef, *_ = np.linalg.lstsq(A, b.ravel(), rcond=None)
        b = b - (A @ coef).reshape(b.shape)
        b -= (b.sum() / max(win.sum(), 1e-30)) * win          # zero total flux, no pedestal
        out.append(b)
    return out[0], out[1], len(R)


def refinePsfSpatial(exposure, psfOld, psfNew, res, badBitmask, P, fitkw, step=0.8,
                     verbose=True):
    """Add or update a linear spatial variation of the PSF, cross-validated.

    Measured from the clean stars' residuals in ``res`` (fitted with ``psfOld``); the gradient
    from half A is judged on half B and vice versa, and kept only if the held-out chi2 drops by at
    least ``psfSpatialMinGain`` compared with the constant ``psfNew``.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    psfOld : `lsst.afw.detection.Psf`
        PSF the residuals were fitted with.
    psfNew : `lsst.afw.detection.Psf`
        Constant PSF to add the variation to.
    res : `dict`
        Result of `runInnerAdaptive`.
    badBitmask : `int`
        Mask bits of pixels to ignore.
    P : `dict`
        Pipeline options (``psf_spatial_*``).
    fitkw : `dict`
        Keyword arguments of the joint fit.
    step : `float`, optional
        Fraction of the measured gradient applied.
    verbose : `bool`, optional
        Print the result.

    Returns
    -------
    psf : `lsst.afw.detection.Psf`
        The varying PSF, or ``psfNew`` if the variation was not kept.
    """
    from scipy.ndimage import gaussian_filter
    bb = exposure.getBBox()
    c = bb.getCenter()
    xc, yc = c.getX(), c.getY()
    scale = 0.5 * max(bb.getWidth(), bb.getHeight())
    fwhm = L._fwhm(psfNew, c)
    K0 = psfNew.computeKernelImage(c).array.astype(float).copy()
    noise = float(np.sqrt(np.nanmedian(exposure.variance.array)))
    idx = selectPsfStars(exposure, res, fwhm, K0.max(), noise, P, verbose=False)
    S, fit = res["S"], res["snaps"][-1]["fit"]
    pos = S["xy"][idx]
    if len(idx) < P["psfSpatialMinStars"] or np.std(pos[:, 0]) < 0.15 * scale \
            or np.std(pos[:, 1]) < 0.15 * scale:
        if verbose:
            print(f"  spatial PSF: {len(idx)} clean star(s), too few / too clustered -- constant PSF")
        return psfNew
    prev = _SPATIAL.get(id(psfOld))
    Kx0, Ky0 = (prev[1], prev[2]) if (prev is not None and prev[0] is psfOld
                                      and prev[1].shape == K0.shape) else (np.zeros_like(K0),) * 2
    data = exposure.image.array
    X0, Y0 = bb.getMinX(), bb.getMinY()
    mr = (P["stackKwargs"] or {}).get("maxRadius", 20)

    def grad(sel):
        bx, by, _ = _psfGradient(data, fit["modelFull"], fit["variance"], fit["good"], (X0, Y0),
                                 S["xy"][sel], fit["amps"][sel], K0, xc, yc, scale,
                                 maxRadius=mr, softSigma=P["psfSoftSigma"])
        if bx is None:
            return None
        if P["psfSmoothPx"]:
            bx, by = gaussian_filter(bx, P["psfSmoothPx"]), gaussian_filter(by, P["psfSmoothPx"])
        return spatialKernelPsf(K0, Kx0 + step * bx, Ky0 + step * by, xc, yc, scale)

    A, B = idx[0::2], idx[1::2]
    mA = _starMask(data.shape, X0, Y0, S["xy"][A], 2.5 * fwhm)
    mB = _starMask(data.shape, X0, Y0, S["xy"][B], 2.5 * fwhm)
    try:
        cA, cB = grad(A), grad(B)
    except Exception as e:
        print(f"  spatial PSF: could not build ({type(e).__name__}: {e}) -- constant PSF")
        return psfNew
    if cA is None or cB is None:
        return psfNew
    f0 = _fit(exposure, psfNew, S, badBitmask, fit["modelFull"], fitkw)
    fA = _fit(exposure, cA, S, badBitmask, fit["modelFull"], fitkw)
    fB = _fit(exposure, cB, S, badBitmask, fit["modelFull"], fitkw)
    before = _chi2On(exposure, f0, mB) + _chi2On(exposure, f0, mA)
    after = _chi2On(exposure, fA, mB) + _chi2On(exposure, fB, mA)
    gain = (before - after) / max(before, 1e-12)
    keep = gain >= P["psfSpatialMinGain"]
    if verbose:
        print(f"  spatial PSF from {len(idx)} clean stars (cross-validated): held-out chi2/px "
              f"{before / 2:.2f} (constant) -> {after / 2:.2f} (varying), {100 * gain:+.1f}% "
              f"-> {'KEPT' if keep else 'not kept'}")
    if not keep:
        return psfNew
    return grad(idx)


def fitBrightnessTerm(exposure, psf, res, badBitmask, P, fitkw, verbose=True):
    """Measure or update the brightness-dependent PSF term (brighter-fatter).

    From isolated, unsaturated stars with peak S/N >= ``bfMinSnr``: per 1-px ring, a weighted
    fit ``residual/flux = A(r) + B(r) * (log10 flux - L0)`` gives the slope B(r), smoothed,
    tapered and made flux-neutral.  Built from half of the stars and judged on the other half;
    kept only if the held-out chi2 improves.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    psf : `lsst.afw.detection.Psf`
        Current PSF.
    res : `dict`
        Result of `runInnerAdaptive`.
    badBitmask : `int`
        Mask bits of pixels to ignore.
    P : `dict`
        Pipeline options (``bf_*``).
    fitkw : `dict`
        Keyword arguments of the joint fit.
    verbose : `bool`, optional
        Print the result.

    Returns
    -------
    moved : `float`
        Fraction of the flux moved per dex of brightness (0 if nothing changed).
    """
    from scipy.ndimage import gaussian_filter1d
    bb = exposure.getBBox()
    c = bb.getCenter()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    kernel = psf.computeKernelImage(c).array.astype(float)
    fwhm = L._fwhm(psf, c)
    noise = float(np.sqrt(np.nanmedian(exposure.variance.array)))
    S = res["S"]
    fit = _fit(exposure, psf, S, badBitmask, res["snaps"][-1]["fit"]["modelFull"], fitkw)
    Pq = dict(P, psfStackMinSnr=P["bfMinSnr"], psfStackN=400)
    idx = selectPsfStars(exposure, dict(res, snaps=[dict(fit=fit)]), fwhm, kernel.max(),
                         noise, Pq, verbose=False)
    amps = fit["amps"]
    if len(idx) < 8 or np.ptp(np.log10(amps[idx])) < 0.5:
        if verbose:
            print(f"  brightness term: {len(idx)} isolated star(s) with S/N >= {P['bfMinSnr']:g} "
                  f"(need >= 8 spanning >= 0.5 dex) -- skipped")
        return 0.0
    idx = idx[np.argsort(amps[idx])]
    lf = np.log10(amps[idx])
    if _BF["L0"] is None:          # the PSF core is tuned on the S/N >= psfCoreMinSnr stars
        snr = amps[idx] * kernel.max() / noise
        coreSet = lf[snr >= P["psfCoreMinSnr"]]
        _BF["L0"] = float(np.median(coreSet if len(coreSet) >= 3 else lf))
    rmax = P["bfRadiusFwhm"] * fwhm
    nr = int(np.ceil(rmax))
    rc = np.arange(nr, dtype=float)                 # ring i = pixels with round(r) == i
    H = kernel.shape[0] // 2
    u = np.arange(-H, H + 1)
    rbk = np.floor(np.hypot(*np.meshgrid(u, u)) + 0.5).astype(int)
    area = np.array([np.sum(rbk == i) for i in range(nr)], float)
    Kring = np.array([kernel[rbk == i].mean() for i in range(nr)]) / kernel.sum()
    data = exposure.image.array
    resid = data - fit["modelFull"]
    var = fit["variance"]
    Yv = np.full((len(idx), nr), np.nan)
    Ys = np.full((len(idx), nr), np.inf)
    for n_, j in enumerate(idx):
        x, y = S["xy"][j]
        xa, xb = max(int(x - rmax) - 1 - X0, 0), min(int(x + rmax) + 2 - X0, data.shape[1])
        ya, yb = max(int(y - rmax) - 1 - Y0, 0), min(int(y + rmax) + 2 - Y0, data.shape[0])
        yy, xx = np.mgrid[ya + Y0:yb + Y0, xa + X0:xb + X0]
        rb = np.floor(np.hypot(xx - x, yy - y) + 0.5).astype(int)
        g = fit["good"][ya:yb, xa:xb]
        rr, vv = resid[ya:yb, xa:xb], var[ya:yb, xa:xb]
        for i in range(nr):
            m = g & (rb == i)
            if m.any():
                Yv[n_, i] = rr[m].mean() / amps[j]
                Ys[n_, i] = np.sqrt(vv[m].sum()) / m.sum() / amps[j]
    t = np.clip((rmax - rc) / (0.4 * rmax), 0.0, 1.0)
    taper = 0.5 * (1 - np.cos(np.pi * t))

    def slope(sel):
        B = np.zeros(nr)
        X = np.column_stack([np.ones(len(sel)), lf[sel] - _BF["L0"]])
        for i in range(nr):
            yv, ys = Yv[sel, i], Ys[sel, i]
            ok = np.isfinite(yv) & np.isfinite(ys)
            if ok.sum() < 4:
                continue
            w = 1.0 / (ys[ok] ** 2 + (0.02 * Kring[i]) ** 2)
            XtW = X[ok].T * w
            try:
                B[i] = np.linalg.solve(XtW @ X[ok], XtW @ yv[ok])[1]
            except np.linalg.LinAlgError:
                pass
        B = gaussian_filter1d(B, 0.7, mode="nearest") * taper
        Kt = Kring * taper
        B -= Kt * np.sum(B * area) / max(np.sum(Kt * area), 1e-30)       # flux-neutral
        return B

    oldR, oldP = _BF["r"], _BF["prof"]
    base = np.interp(rc, oldR, oldP, right=0.0) if oldP is not None else np.zeros(nr)

    def use(r, prof):
        _BF.update(r=r, prof=prof)
        _setBfFluxes(S["xy"], amps)
        _STAMP_CACHE.clear()

    iA, iB = np.arange(0, len(idx), 2), np.arange(1, len(idx), 2)
    mA = _starMask(data.shape, X0, Y0, S["xy"][idx[iA]], rmax)
    mB = _starMask(data.shape, X0, Y0, S["xy"][idx[iB]], rmax)
    gain, before, after = -1.0, np.nan, np.nan
    try:
        use(rc, base + slope(iA))
        fA = _fit(exposure, psf, S, badBitmask, fit["modelFull"], fitkw)
        use(rc, base + slope(iB))
        fB = _fit(exposure, psf, S, badBitmask, fit["modelFull"], fitkw)
        before = _chi2On(exposure, fit, mB) + _chi2On(exposure, fit, mA)
        after = _chi2On(exposure, fA, mB) + _chi2On(exposure, fB, mA)
        gain = (before - after) / max(before, 1e-12)
    except ValueError as e:
        print(f"  brightness term: {e}")
    if gain < P["psfCvMinGain"]:
        if oldP is None:
            _BF.update(r=None, prof=None)
            _STAMP_CACHE.clear()
        else:
            use(oldR, oldP)
        if verbose:
            print(f"  brightness term: {len(idx)} stars over {np.ptp(lf):.1f} dex, held-out chi2 "
                  f"{100 * gain:+.1f}% -> not changed")
        return 0.0
    dB = slope(np.arange(len(idx)))
    use(rc, base + dB)
    moved = float(0.5 * np.sum(np.abs(dB) * area))
    if verbose:
        core = 100 * _BF["prof"][0] / Kring[0]
        print(f"  brightness term: {len(idx)} stars over {np.ptp(lf):.1f} dex (L0 = 10^{_BF['L0']:.2f}), "
              f"held-out chi2/px {before / 2:.2f} -> {after / 2:.2f} ({100 * gain:+.1f}%) -> KEPT; "
              f"PSF peak changes {core:+.1f}% per dex of flux")
    return moved


BG_STATS_MASK = ("BAD", "SAT", "EDGE", "NO_DATA", "DETECTED", "CR", "INTRP")


def refitBackgroundLowfit(exposure, mode="plane", statsMask=BG_STATS_MASK, clip=3.0, nIter=5,
                          verbose=True):
    """Fit a constant or tilted-plane background to the unmasked pixels and subtract it in place.

    The robust alternative to `fitLib.refitBackground` for small cutouts of an image that is
    already background-subtracted: a spline / Chebyshev surface on 32 px bins is barely
    constrained where most pixels are masked and can swing by hundreds of counts in a corner.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure to background-subtract in place (sources flagged DETECTED).
    mode : `str`, optional
        'plane' (a + b x + c y), 'constant' or 'none'.
    statsMask : `tuple` [`str`], optional
        Mask planes whose pixels are ignored.
    clip : `float`, optional
        Sigma clipping threshold.
    nIter : `int`, optional
        Clipping iterations.
    verbose : `bool`, optional
        Print the fitted background.

    Returns
    -------
    exposure : `lsst.afw.image.Exposure`
        The background-subtracted exposure (same object).
    background : `lsst.afw.image.ImageF`
        The subtracted background image.

    Raises
    ------
    ValueError
        Raised if ``mode`` is not 'plane', 'constant' or 'none'.
    """
    if mode not in ("plane", "constant", "none"):
        raise ValueError(f"bgMode must be 'plane', 'constant', 'spline' or 'none', not {mode!r}")
    img = exposure.image.array
    ny, nx = img.shape
    bits = 0
    for name in statsMask:
        try:
            bits |= exposure.mask.getPlaneBitMask(name)
        except Exception:
            pass
    ok = ((exposure.mask.array & bits) == 0) & np.isfinite(img)
    yy, xx = np.mgrid[:ny, :nx]
    xs, ys = (xx - nx / 2) / max(nx, 1), (yy - ny / 2) / max(ny, 1)
    coef = np.zeros(3)
    if mode != "none" and ok.sum() >= 10:
        use = ok.copy()
        for _ in range(nIter):
            z = img[use]
            if mode == "plane":
                A = np.column_stack([np.ones(z.size), xs[use], ys[use]])
                coef = np.linalg.lstsq(A, z, rcond=None)[0]
            else:
                coef = np.array([np.median(z), 0.0, 0.0])
            r = img - (coef[0] + coef[1] * xs + coef[2] * ys)
            s = 1.4826 * np.median(np.abs(r[use] - np.median(r[use])))
            new = ok & (np.abs(r) < clip * max(s, 1e-12))
            if new.sum() < 10 or np.array_equal(new, use):
                break
            use = new
    bg = (coef[0] + coef[1] * xs + coef[2] * ys).astype(np.float32)
    img -= bg
    if verbose:
        print(f"  background refit ({mode}): level {coef[0]:.4g}, slope x {coef[1]:.3g}, "
              f"y {coef[2]:.3g} per cutout width, from {int(ok.sum())} px "
              f"({100 * ok.mean():.0f}% unmasked)")
    try:
        bgImg = L.afwImage.ImageF(bg)
        bgImg.setXY0(exposure.getBBox().getMin())
        if not isinstance(getattr(bgImg, "array", None), np.ndarray):
            raise TypeError("no afw image")
    except Exception:                                   # without the stack (tests)
        bgImg = types.SimpleNamespace(array=bg)
    return exposure, bgImg


def _kernelChange(psfOld, psfNew, point):
    """Return the relative L2 change between two PSFs' kernels at a point."""
    a = psfOld.computeKernelImage(point).array
    b = psfNew.computeKernelImage(point).array
    h = min(a.shape[0], b.shape[0]) // 2
    ca_, cb_ = a.shape[0] // 2, b.shape[0] // 2
    a = a[ca_ - h:ca_ + h + 1, ca_ - h:ca_ + h + 1]
    b = b[cb_ - h:cb_ + h + 1, cb_ - h:cb_ + h + 1]
    return float(np.linalg.norm(b - a) / max(np.linalg.norm(a), 1e-30))


def printTiming(timer=None):
    """Print the time spent per pipeline step.

    Parameters
    ----------
    timer : `lsst.meas.crowd.fitLib.StepTimer`, optional
        Timer (default: fitLib's ``TIMER``).
    """
    timer = timer or L.TIMER
    tot, cnt = defaultdict(float), defaultdict(int)
    for r in timer.records:
        if r["kind"] == "compute":
            tot[r["label"]] += r["seconds"]
            cnt[r["label"]] += 1
    print(f"\n{'step':40s} {'calls':>6s} {'total [s]':>10s}")
    for k, v in sorted(tot.items(), key=lambda kv: -kv[1])[:20]:
        print(f"{k:40s} {cnt[k]:6d} {v:10.2f}")
    print(f"wall clock: {time.perf_counter() - timer.tStart:.1f}s")


def runAdaptivePipeline(exposureCrop, psf, badBitmask, gaiaXy=None, gaiaRows=None,
                        starCat=None, bgImage=None, viz=False, gaiaGmag=None, **cfg):
    """Run the adaptive pipeline on a whole cutout, starting from Gaia positions.

    The replacement for `fitLib.runCutoutPipeline`: outer passes of `runInnerAdaptive`, with
    background refit and PSF refinement in between, until the source list and the PSF stop
    changing.

    Parameters
    ----------
    exposureCrop : `lsst.afw.image.Exposure`
        Cutout before background subtraction.
    psf : `lsst.afw.detection.Psf`
        Initial PSF.
    badBitmask : `int`
        Mask bits of pixels to ignore.
    gaiaXy : `numpy.ndarray`, optional
        Gaia starting positions (detector px), shape (N, 2); None or empty = no Gaia.
    gaiaRows : `numpy.ndarray`, optional
        Gaia catalogue row of each starting position.
    starCat : `lsst.afw.table.SourceCatalog`, optional
        Detector-level PSF stars.
    bgImage : `lsst.afw.image.Image`, optional
        Detector-level background, subtracted in the first pass.
    viz : `bool`, optional
        Draw the per-pass overview.
    gaiaGmag : `numpy.ndarray`, optional
        Gaia G magnitude of each starting position.
    **cfg
        Configuration: `fitLib.runCutoutPipeline` keys and ``PIPELINE_DEFAULTS`` keys.

    Returns
    -------
    outerHistory : `list` [`dict`]
        One dict per outer pass (as `fitLib.runCutoutPipeline`'s, plus ``psfChange``,
        ``modelErr``, ``bf`` and ``result['adaptive']``).
    """
    P = _params(cfg)
    _FAST["on"] = bool(P["fastStamps"])
    _FAST["cache"].clear()
    _bfReset()
    print(f"meas_crowd {__version__}: brightness term {'ON' if P['bfTerm'] else 'off'}, "
          f"PSF core median {'ON' if P['psfCoreMedian'] else 'off'}, "
          f"PSF regularise {'ON' if P['psfRegularize'] else 'off'}")
    enableSpeedups()
    L.TIMER.reset()
    cutout = exposureCrop.clone()
    cutout.setPsf(psf)
    cb = cutout.getBBox()
    c = cb.getCenter()
    if P["satGrowPx"]:
        from scipy.ndimage import binary_dilation
        sat = (cutout.mask.array & cutout.mask.getPlaneBitMask("SAT")) != 0
        if sat.any():
            grown = binary_dilation(sat, iterations=int(P["satGrowPx"])) & ~sat
            cutout.mask.addMaskPlane("SAT_GROWN")
            bit = cutout.mask.getPlaneBitMask("SAT_GROWN")
            cutout.mask.array[grown] |= bit
            badBitmask = badBitmask | bit
            print(f"excluding {int(grown.sum())} px within {P['satGrowPx']}px of saturated pixels")
    if P["nearSatFrac"]:
        lvl, how = saturationLevel(cutout, P)
        if lvl is None:
            print(f"near-saturated cores: {how}")
        else:
            from scipy.ndimage import binary_dilation, label
            dat = cutout.image.array
            sat = (cutout.mask.array & cutout.mask.getPlaneBitMask("SAT")) != 0
            ns = np.isfinite(dat) & (dat > P["nearSatFrac"] * lvl) & ~sat
            if ns.any():
                if P["nearSatGrowPx"]:
                    ns = binary_dilation(ns, iterations=int(P["nearSatGrowPx"]))
                ns &= (cutout.mask.array & badBitmask) == 0
                cutout.mask.addMaskPlane("NEARSAT")
                bit = cutout.mask.getPlaneBitMask("NEARSAT")
                cutout.mask.array[ns] |= bit
                badBitmask = badBitmask | bit
                nObj = label(ns)[1]
            else:
                nObj = 0
            print(f"near-saturated cores: saturation level {lvl:.4g} ({how}); excluding "
                  f"{int(ns.sum())} px above {P['nearSatFrac']:g} x that level in {nObj} "
                  f"star core(s) -- fitted from their wings")
    P["_satZones"], P["_satZoneMask"] = satZones(cutout, L._fwhm(psf, c), P)
    P["_haloRejects"] = np.zeros((0, 2))
    if P["psfWingSize"]:
        psf = extendPsfWingsKeep(psf, c, P["psfWingSize"], blend=P["psfExtBlend"],
                                 keepAngular=P["wingKeepAngular"])
    psf = regularizePsf(psf, c, P, label=" (initial)")
    raw = cutout.clone()
    expW = cutout.clone()
    if P["subtractInitialBg"] and bgImage is not None:
        try:
            expW.image.array[:] -= bgImage[cb].array
            print("pass 1 uses the detector-level background (subtracted)")
        except Exception as e:
            print(f"(could not subtract initial background: {e})")
    expW.setPsf(psf)

    S = _state(gaiaXy if gaiaXy is not None else None, fromGaia=True, gaia=gaiaRows,
               gmag=gaiaGmag)
    if len(S["xy"]) and P["fitGaiaOffset"]:
        noise = float(np.sqrt(np.nanmedian(expW.variance.array)))
        dx, dy, n = estimateGaiaOffset(expW, S, L._fwhm(psf, c), badBitmask, noise)
        if n:
            S["xy"] += (dx, dy)
            S["ref"] += (dx, dy)
    print(f"{len(S['xy'])} Gaia starting position(s) for cutout {cb}")

    starpos = None
    if starCat is not None and len(starCat):
        starpos = [(r.getX(), r.getY()) for r in starCat
                   if cb.contains(L.geom.Point2I(int(round(r.getX())), int(round(r.getY()))))]
    psfW, bgW, psfChange = psf, None, np.nan
    hist = []
    for k in range(P["maxOuter"]):
        clearStampCache()
        shapeNow = P["useShapeCut"] and (k > 0 or P["shapeInFirstPass"])
        thr0 = (P["gaiaThresholdStart"] if len(S["xy"]) else P["thresholdStart"]) if k == 0 \
            else P["thresholdMin"]
        print(f"\n{'=' * 70}\n===== Outer pass {k + 1} (cap {P['maxOuter']}), PSF FWHM "
              f"{L._fwhm(psfW, c):.3f}px, shape cut {'ON' if shapeNow else 'off'} =====")
        t0 = time.perf_counter()
        with L.TIMER.step("INNER LOOP (total)", outer=k, kind="wall"):
            res = runInnerAdaptive(expW, psfW, badBitmask, S, P, outerIdx=k,
                                   thresholdStart=thr0, shapeNow=shapeNow)
        S = res["S"]
        res["adaptive"]["seconds"] = time.perf_counter() - t0
        hist.append(dict(psf=psfW, exposure=expW, result=res, version=__version__,
                         bf=None if _BF["prof"] is None else
                         dict(r=_BF["r"].copy(), prof=_BF["prof"].copy(), L0=_BF["L0"]),
                         residStd=res["residStdHistory"][-1], nPeaks=len(res["peaks"]),
                         bgImage=bgW, psfChange=psfChange, badBitmask=badBitmask,
                         rejects=np.asarray(P["_haloRejects"], float).reshape(-1, 2).copy(),
                         poissonK=P["poissonK"]))
        if P["autoPoisson"] and res["snaps"]:
            try:
                import contextlib
                import io
                import types
                fake = types.SimpleNamespace(outerHistory=[hist[-1]], badBitmask=badBitmask,
                                             config=cfg)
                with contextlib.redirect_stdout(io.StringIO()):
                    dg = diagnoseResidual(fake, plot=False)
                share = dg.get("noiseShare", 0.0)
                hfR = dg["highfreq"]["ratio"]
                kNew = float(dg["highfreq"]["k"])
                sz = dg["bright-star disks"]["sigmaZ"]
                ratio = share
                if share >= 0.5 and len(hfR) and hfR[-1] > 1.15 and np.isfinite(kNew) \
                        and kNew > 0:
                    P["poissonK"] = max(P["poissonK"], kNew)
                    hist[-1]["poissonK"] = P["poissonK"]          # used by sigma='effective'
                    print(f"  noise check: bright-star excess is pixel-scale noise ({100 * share:.0f}%) "
                          f"-> variance += {P['poissonK']:.4g} x model from now on")
                else:
                    print(f"  noise check: bright-star disks sigmaZ {sz:.2f}; missing photon noise "
                          f"explains {100 * share:.0f}% of the excess -> variance unchanged")
            except Exception as e:
                print(f"  noise check failed ({type(e).__name__}: {e})")
        if P["modelErrMeasure"] and res["snaps"]:
            try:
                eps, _ = measureModelError(expW, res, P)
                hist[-1]["modelErr"] = eps
                if P["fitSysFromData"] and eps > 0:
                    P["sysFrac"] = max(eps, 0.002)
                    print(f"  -> fit weights use sysFrac = {P['sysFrac']:.3g} from now on")
            except Exception as e:
                print(f"  model-error measurement failed ({type(e).__name__}: {e})")
        ad = res["adaptive"]
        nNow = len(res["peaks"])
        growth = (abs(nNow - hist[-2]["nPeaks"]) / max(hist[-2]["nPeaks"], 1)
                  if len(hist) > 1 else np.inf)
        psfSettled = (not P["refitPsf"]) or (np.isfinite(psfChange) and psfChange < P["psfTol"])
        print(f"\nouter pass {k + 1}: {nNow} sources ({int(res['isGaia'].sum())} Gaia), "
              f"{ad['nInner']} inner iteration(s), stop: {ad['stop']}, residual "
              f"{'has EXCESS' if ad['missing'] else 'clean'}, source change "
              f"{100 * growth:.1f}%, PSF change into this pass {100 * psfChange:.2f}%, "
              f"{ad['seconds']:.0f}s")
        if viz and res["snaps"]:
            try:
                L.plotIterationOverview(expW, res["snaps"][-1],
                                        dict(residStd=res["residStdHistory"],
                                             nPeaks=res["nPeaksHistory"],
                                             nZero=res["nZeroAmpHistory"]),
                                        k, L._dataNorm(expW.image.array))
            except Exception as e:
                print(f"(overview plot failed: {e})")
        if k + 1 >= P["maxOuter"]:
            if ad["missing"]:
                print(f"WARNING: maxOuter={P['maxOuter']} reached with excess remaining")
            break
        if k + 1 >= P["minOuter"] and growth < P["outerGrowthTol"] and psfSettled:
            print("outer loop converged: " + ("residual clean, " if not ad["missing"] else
                  "excess remains but nothing changes any more, ") + "sources and PSF stable")
            break
        if not res["snaps"]:
            print("no sources -- stopping")
            break

        # ---- refit background + PSF from this pass's model ---------------------
        expN = raw.clone()
        with L.TIMER.step("mask + refit background", outer=k):
            srcMask = L.maskSourcesFromModel(expN, res["modelFull"])
            if P["bgSatMaskFwhm"]:
                from scipy.ndimage import distance_transform_edt
                satm = (expN.mask.array & expN.mask.getPlaneBitMask("SAT")) != 0
                if satm.any():
                    halo = distance_transform_edt(~satm) <= P["bgSatMaskFwhm"] * L._fwhm(psfW, c)
                    expN.mask.array[halo] |= expN.mask.getPlaneBitMask("DETECTED")
                    print(f"  background refit also ignores saturated-star halos "
                          f"({100 * halo.mean():.1f}% of the cutout)")
            try:
                if P["bgMode"] == "spline":
                    expN, bgW = L.refitBackground(expN, useApprox=P["bgUseApprox"],
                                                  binSize=P["bgBinSize"])
                else:
                    expN, bgW = refitBackgroundLowfit(expN, P["bgMode"])
            except Exception as e:
                print(f"  background refit failed ({type(e).__name__}: {e}) -- keeping image")
        psfNew = psfW
        expN.setPsf(psfW)
        if P["refitPsf"] == "stack":
            with L.TIMER.step("PSF refit (residual stack)", outer=k):
                try:
                    psfNew = refinePsfClean(
                        expW, psfW, res, badBitmask, P,
                        dict(supportFrac=P["supportFrac"], sysFrac=P["sysFrac"],
                             poissonK=P["poissonK"]),
                        nCycles=P["psfPolishCycles"])
                except Exception as e:
                    print(f"  residual-stack PSF refit failed ({type(e).__name__}: {e})")
                    psfNew = psfW
        elif P["refitPsf"]:
            ks = P["psfKernelSizes"] or tuple([P["psfKernelSize"]] +
                                              [s for s in (27, 21, 15, 11) if s < P["psfKernelSize"]])
            with L.TIMER.step("redetect + PSF refit", outer=k):
                try:
                    p_, _, _, cell, _ = L.redetectAndRefitPsfAdaptive(
                        expN, kernelSizes=ks, extraPeakPositions=res["peaks"],
                        forceIncludeStarPositions=starpos)
                    if cell is not None:
                        psfNew = p_
                except Exception as e:
                    print(f"  PSF redetermination failed ({type(e).__name__}: {e})")
        if P["psfWingSize"]:
            wp = None
            if P["wingsFromStars"]:
                try:
                    wp = wingProfilePsf(expW, res, psfNew, P)
                except Exception as e:
                    print(f"  wings from stars failed ({type(e).__name__}: {e})")
            if wp is not None:
                psfNew = wp
            elif P["wingsFromStars"] and psfNew is not psfW and P["refitPsf"] == "stack":
                pass                       # stack-refined PSF keeps its current (embedded) wings
            elif psfNew is not psfW:
                psfNew = extendPsfWingsKeep(psfNew, c, P["psfWingSize"], blend=P["psfExtBlend"],
                                            keepAngular=P["wingKeepAngular"])
        if P["psfCoreMedian"]:
            with L.TIMER.step("PSF core (median stack)", outer=k):
                try:
                    psfNew = refinePsfMedian(
                        expW, psfNew, res, badBitmask, P,
                        dict(supportFrac=P["supportFrac"], sysFrac=P["sysFrac"],
                             poissonK=P["poissonK"]),
                        fit0=res["snaps"][-1]["fit"] if psfNew is psfW else None)
                except Exception as e:
                    print(f"  PSF core (median stack) failed ({type(e).__name__}: {e})")
        if P["psfSpatial"] and P["refitPsf"] == "stack":
            try:
                psfNew = refinePsfSpatial(
                    expW, psfW, psfNew, res, badBitmask, P,
                    dict(supportFrac=P["supportFrac"], sysFrac=P["sysFrac"],
                         poissonK=P["poissonK"]))
            except Exception as e:
                print(f"  spatial PSF failed ({type(e).__name__}: {e}) -- constant PSF")
        if psfNew is not psfW:
            psfNew = regularizePsf(psfNew, c, P)
        psfChange = _kernelChange(psfW, psfNew, c) if psfNew is not psfW else 0.0
        if P["bfTerm"]:
            with L.TIMER.step("brightness-dependent PSF", outer=k):
                try:
                    for _ in range(2):
                        moved = fitBrightnessTerm(
                            expW, psfNew, res, badBitmask, P,
                            dict(supportFrac=P["supportFrac"], sysFrac=P["sysFrac"],
                                 poissonK=P["poissonK"]))
                        psfChange = max(psfChange, moved)     # a change -> not settled yet
                        if moved < P["psfTol"]:
                            break
                except Exception as e:
                    print(f"  brightness term failed ({type(e).__name__}: {e})")
        expN.setPsf(psfNew)
        expW, psfW = expN, psfNew
    printTiming()
    return hist


def dropOffDetectorGaia(exposure, xy, rows):
    """Drop Gaia starting positions that lie beyond the detector edge.

    A star outside the cutout whose nearest cutout pixel is flagged EDGE or NO_DATA lies off the
    detector: no pixel constrains its flux, so the joint fit can give it an absurd value whose
    wings then spoil the edge of the image.  Stars outside a cutout edge that is not a detector
    edge are kept (their wings fall on real pixels and their neighbours need them).

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The cutout.
    xy : `numpy.ndarray`
        Gaia positions (detector px), shape (N, 2).
    rows : `numpy.ndarray`
        Their rows in the Gaia table.

    Returns
    -------
    xy, rows : `numpy.ndarray`
        The positions and rows that are kept.
    """
    xy = np.asarray(xy, float).reshape(-1, 2)
    if not len(xy):
        return xy, rows
    bb = exposure.getBBox()
    mk = exposure.mask
    bits = 0
    for name in ("EDGE", "NO_DATA"):
        try:
            bits |= mk.getPlaneBitMask(name)
        except Exception:
            pass
    if not bits:
        return xy, rows
    ix = np.round(xy[:, 0]).astype(int) - bb.getMinX()
    iy = np.round(xy[:, 1]).astype(int) - bb.getMinY()
    ny, nx = mk.array.shape
    outside = (ix < 0) | (ix >= nx) | (iy < 0) | (iy >= ny)
    nearest = mk.array[np.clip(iy, 0, ny - 1), np.clip(ix, 0, nx - 1)]
    keep = ~(outside & ((nearest & bits) != 0))
    if (~keep).any():
        print(f"  dropped {int((~keep).sum())} Gaia star(s) beyond the detector edge")
    return xy[keep], np.asarray(rows)[keep]


def runCutoutGaia(rawCrop, psf, badBitmask, gaiaTbl=None, config=None, interiorBbox=None,
                  starCat=None, bgImage=None, viz=False, name="cutout", gmagMax=None,
                  **overrides):
    """Run the adaptive pipeline on one cropped exposure, starting from Gaia positions.

    Parameters
    ----------
    rawCrop : `lsst.afw.image.Exposure`
        Cutout before background subtraction.
    psf : `lsst.afw.detection.Psf`
        Initial PSF.
    badBitmask : `int`
        Mask bits of pixels to ignore.
    gaiaTbl : `astropy.table.Table`, optional
        Gaia catalogue (see `lsst.meas.crowd.gaiaCatalog.loadGaiaCatalog`); None = no Gaia.
    config : `dict`, optional
        Pipeline configuration.
    interiorBbox : `lsst.geom.Box2I`, optional
        The requested box without padding (default: the whole cutout).
    starCat : `lsst.afw.table.SourceCatalog`, optional
        Detector-level PSF stars.
    bgImage : `lsst.afw.image.Image`, optional
        Detector-level background.
    viz : `bool`, optional
        Draw the per-pass overview.
    name : `str`, optional
        Cutout name.
    gmagMax : `float`, optional
        Only start from Gaia stars brighter than this G magnitude.
    **overrides
        Configuration overrides.

    Returns
    -------
    run : `lsst.meas.crowd.session.CutoutRun`
        The result; its catalogue also has ``fromGaia``, ``gaiaSourceId``, ``gaiaG``,
        ``gaiaDx`` and ``gaiaDy``.
    """
    from . import session as cn
    cfg = cn.makeConfig(config, **overrides)
    P = _params(cfg)
    cb = rawCrop.getBBox()
    interiorBbox = interiorBbox or cb
    gaiaXy, rows = selectGaiaSources(gaiaTbl, cb, L._fwhm(psf, cb.getCenter()),
                                     P["gaiaPadFwhm"], gmagMax, P["gaiaMinSepPx"])
    if P["gaiaDropOffDetector"]:
        gaiaXy, rows = dropOffDetectorGaia(rawCrop, gaiaXy, rows)
    t0 = time.perf_counter()
    gm = (np.asarray(gaiaTbl["phot_g_mean_mag"], float)[rows]
          if gaiaTbl is not None and len(rows) and "phot_g_mean_mag" in gaiaTbl.colnames else None)
    hist = runAdaptivePipeline(rawCrop, psf, badBitmask, gaiaXy=gaiaXy, gaiaRows=rows,
                               starCat=starCat, bgImage=bgImage, viz=viz, gaiaGmag=gm, **cfg)
    res = hist[-1]["result"]
    S = res["S"]
    cat = cn.buildCatalog(hist, interiorBbox)
    cat["fromGaia"] = S["fromGaia"]
    gi = S["gaia"]
    has = gi >= 0
    sid = np.full(len(gi), -1, np.int64)
    gG = np.full(len(gi), np.nan)
    if gaiaTbl is not None and has.any():
        sid[has] = np.asarray(gaiaTbl["source_id"], np.int64)[gi[has]]
        if "phot_g_mean_mag" in gaiaTbl.colnames:
            gG[has] = np.asarray(gaiaTbl["phot_g_mean_mag"], float)[gi[has]]
    cat["gaiaSourceId"], cat["gaiaG"] = sid, gG
    d = S["xy"] - S["ref"]
    cat["gaiaDx"], cat["gaiaDy"] = d[:, 0], d[:, 1]
    cat["extended"] = S["ext"] > 0
    cat["extSigmaPx"] = S["ext"]
    summary = dict(name=name, status="ok", nSources=len(res["peaks"]),
                   nGaiaInput=int(len(gaiaXy)), nFromGaia=int(S["fromGaia"].sum()),
                   nInInterior=int(np.sum(cat["inInterior"])),
                   bboxFitted=[cb.getMinX(), cb.getMaxX(), cb.getMinY(), cb.getMaxY()],
                   interiorBbox=[interiorBbox.getMinX(), interiorBbox.getMaxX(),
                                 interiorBbox.getMinY(), interiorBbox.getMaxY()],
                   residStdPerOuter=[h["residStd"] for h in hist],
                   nPeaksPerOuter=[h["nPeaks"] for h in hist],
                   psfFwhmPxPerOuter=[L._fwhm(h["psf"], h["exposure"].getBBox().getCenter())
                                      for h in hist],
                   psfChangePerOuter=[h["psfChange"] for h in hist],
                   adaptivePerOuter=[{k: v for k, v in h["result"]["adaptive"].items()
                                      if k != "final"} for h in hist],
                   residualClean=not res["adaptive"]["missing"],
                   wallSeconds=time.perf_counter() - t0, pipelineKwargs=cfg)
    run = cn.CutoutRun(name, hist, cat, summary, rawCrop, interiorBbox,
                       hist[0].get("badBitmask", badBitmask), cfg, gaia=gaiaTbl)
    try:                                 # the photometry catalogue, with the Gaia comparison
        from .photometry import measurePhotometry
        run.photometry = measurePhotometry(run)
    except Exception as e:
        print(f"photometry failed ({type(e).__name__}: {e}); run.photometry is None")
    return run


_PAR = {}


def _parWorker(name):
    """Fit one cutout in a forked worker process (inherits the session, nothing is pickled)."""
    import pickle
    from . import session as cn
    session, config, outDir, gmagMax, overrides = _PAR["args"]
    job = session.jobs[name]
    run = runCutoutGaia(job["raw"], session.psf, session.badBitmask, session.gaia, config,
                        interiorBbox=job["interior"], starCat=session.starCat,
                        bgImage=session.bgImage, viz=False, name=name, gmagMax=gmagMax,
                        **overrides)
    path = None
    if outDir:
        path = os.path.join(outDir, f"visit{session.visit}_det{session.detector}", name)
        cn.saveRun(run, path)
    try:
        pickle.dumps(run)
        return name, run, run.summary, path
    except Exception as e:                     # LSST objects that do not pickle
        print(f"[{name}] result not transferable to the notebook ({e}); saved to {path}")
        return name, None, run.summary, path


def checkBrightNeighbours(run, brightSnr=300.0, withinFwhm=4.0, fluxRatio=5.0,
                          nControls=7, plot=True, maxShow=16, verbose=True, similarMin=3):
    """Test whether sources on the edges of bright stars are real.

    For every source within ``withinFwhm`` FWHM of a star at least ``fluxRatio`` times brighter
    (with peak S/N >= ``brightSnr``), its PSF-weighted flux is compared with the residual's flux
    at control positions at the same distance from the bright star, rotated in 45-degree steps.
    A mismatch pattern of the bright star (rings, spikes, lobes) puts similar flux at several
    angles; a real star is significantly brighter than its controls.  Verdicts: 'Gaia' (real by
    construction), 'real' (significance >= 5 and at most one similar angle), 'probable'
    (significance >= 3, or two similar angles), 'doubtful' (significance < 3, or similar flux at
    >= ``similarMin`` angles).

    Parameters
    ----------
    run : `lsst.meas.crowd.session.CutoutRun`
        A finished run.
    brightSnr : `float`, optional
        Minimum peak S/N of a bright star.
    withinFwhm : `float`, optional
        Neighbour radius [FWHM].
    fluxRatio : `float`, optional
        Minimum brightness ratio of the bright star to its neighbour.
    nControls : `int`, optional
        Number of control angles.
    plot : `bool`, optional
        Show thumbnails of the non-Gaia neighbours.
    maxShow : `int`, optional
        Maximum number of thumbnails.
    verbose : `bool`, optional
        Print the table.
    similarMin : `int`, optional
        Number of similar control angles that makes a source doubtful.

    Returns
    -------
    table : `astropy.table.Table`
        One row per neighbour: position, flux, significance, number of similar angles, verdict.
    """
    final = run.outerHistory[-1]
    exp, res = final["exposure"], final["result"]
    S, fit = res["S"], res["snaps"][-1]["fit"]
    bb = exp.getBBox()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    data = exp.image.array
    model = fit["modelFull"]          # PSF-only joint model: bright-star corrections would hide
    var, good = exp.variance.array, fit["good"]   # the repeating halo pattern at the controls
    psf = exp.getPsf()
    c = bb.getCenter()
    fwhm = L._fwhm(psf, c)
    psfPeak = psf.computeKernelImage(c).array.max()
    noise = float(np.sqrt(np.nanmedian(var)))
    amps = fit["amps"]
    snr = amps * psfPeak / noise
    xy = S["xy"]
    bright = np.flatnonzero(snr >= brightSnr)
    rows = []
    if not len(bright):
        print("no bright stars")
        return Table()
    btree = cKDTree(xy[bright])
    alltree = cKDTree(xy)
    resid = data - model
    R2 = (1.5 * fwhm) ** 2

    def psfFlux(img, x, y):
        """PSF-weighted amplitude (and its error) of img at (x, y) over good pixels."""
        st, sb = L.psfStamp(psf, x, y, bb, taper=L.WING_TAPER)
        if st is None:
            return np.nan, np.nan
        ys, xs = L.bboxToSlices(sb, bb)
        yy, xx = np.mgrid[sb.getMinY():sb.getMaxY() + 1, sb.getMinX():sb.getMaxX() + 1]
        m = good[ys, xs] & (((xx - x) ** 2 + (yy - y) ** 2) < R2) & (var[ys, xs] > 0)
        if m.sum() < 8:
            return np.nan, np.nan
        w = 1.0 / var[ys, xs][m]
        p = st[m]
        den = np.sum(p * p * w)
        return float(np.sum(p * img[ys, xs][m] * w) / den), float(1.0 / np.sqrt(den))

    for j in range(len(xy)):
        if snr[j] >= brightSnr:
            continue
        for k in btree.query_ball_point(xy[j], withinFwhm * fwhm):
            b = bright[k]
            if amps[b] < fluxRatio * max(amps[j], 1e-9):
                continue
            cj = fit["comps"][j]
            own = np.zeros_like(data)
            if cj is not None:
                own[cj["ys"], cj["xs"]] = amps[j] * cj["stamp"]
            aJ, eJ = psfFlux(resid + own, *xy[j])
            dx, dy = xy[j] - xy[b]
            d = float(np.hypot(dx, dy))
            th0 = np.arctan2(dy, dx)
            ctrl = []
            for q in range(1, nControls + 1):
                th = th0 + q * 2 * np.pi / (nControls + 1)
                cx, cy = xy[b, 0] + d * np.cos(th), xy[b, 1] + d * np.sin(th)
                if not (X0 <= cx <= bb.getMaxX() and Y0 <= cy <= bb.getMaxY()):
                    continue
                # a control on top of another (non-negligible) source measures THAT source's
                # mismatch, not the bright star's halo pattern -> skip it
                near = [k2 for k2 in alltree.query_ball_point((cx, cy), 1.5 * fwhm)
                        if k2 not in (j, b) and amps[k2] > 0.3 * max(amps[j], 1e-9)]
                if near:
                    continue
                aC, _ = psfFlux(resid, cx, cy)
                if np.isfinite(aC):
                    ctrl.append(aC)
            ctrl = np.array(ctrl)
            nSim = 0
            if len(ctrl) >= 3 and np.isfinite(aJ):
                med = float(np.median(ctrl))
                spread = float(1.4826 * np.median(np.abs(ctrl - med)))
                sig = (aJ - med) / max(eJ, spread, 1e-9)
                nSim = int((ctrl >= 0.5 * aJ).sum()) if aJ > 0 else len(ctrl)
            else:
                med, spread, sig = np.nan, np.nan, np.nan
            if S["fromGaia"][j]:
                verdict = "Gaia"
            elif not np.isfinite(sig):
                verdict = "n/a"
            elif nSim >= similarMin:   # the same "flux" repeats at other angles: halo pattern
                verdict = "doubtful"
            elif sig >= 5 and nSim < 2:
                verdict = "real"
            elif sig >= 3:              # (two similar angles can be chance in a crowded field)
                verdict = "probable"
            else:
                verdict = "doubtful"
            rows.append((j, xy[j, 0], xy[j, 1], int(b), d / fwhm, amps[j] / amps[b], snr[j],
                         aJ, med, spread, sig, len(ctrl), nSim, verdict))
            break
    t = Table(rows=rows, names=("id", "x", "y", "brightId", "distFwhm", "fluxRatio",
                                "peakSnr", "flux", "controlMedian", "controlScatter",
                                "significance", "nControls", "nSimilar", "verdict")) if rows else Table()
    if not len(t):
        print("no sources next to bright stars")
        return t
    if verbose:
        cnt = {v: int((t["verdict"] == v).sum()) for v in ("Gaia", "real", "probable", "doubtful", "n/a")}
        print(f"{len(t)} source(s) within {withinFwhm:g} FWHM of a >= {fluxRatio:g}x brighter star "
              f"(peak S/N >= {brightSnr:g}): " + ", ".join(f"{v} {n}" for v, n in cnt.items() if n))
        print(f"{'id':>5s} {'x':>8s} {'y':>8s} {'d[FWHM]':>7s} {'f/f_b':>7s} {'S/N':>6s} "
              f"{'signif':>6s} {'similar':>7s}  verdict   (similar = control angles with >= half "
              f"its flux)")
        for r in sorted(t, key=lambda r: (r["verdict"] != "doubtful", np.nan_to_num(r["significance"]))):
            print(f"{r['id']:5d} {r['x']:8.1f} {r['y']:8.1f} {r['distFwhm']:7.2f} {r['fluxRatio']:7.3f} "
                  f"{r['peakSnr']:6.1f} {r['significance']:6.1f} {r['nSimilar']:4d}/{r['nControls']:<2d}"
                  f"  {r['verdict']}")
    if plot:
        import matplotlib.pyplot as plt
        from astropy.visualization import simple_norm
        show = [r for r in sorted(t, key=lambda r: np.nan_to_num(r["significance"]))
                if r["verdict"] != "Gaia"][:maxShow]
        if show:
            h = int(round(3 * fwhm))
            fig, axes = plt.subplots(len(show), 3, figsize=(10, 3.2 * len(show)), squeeze=False)
            for axRow, r in zip(axes, show):
                j, b = int(r["id"]), int(r["brightId"])
                cj = fit["comps"][j]
                own = np.zeros_like(data)
                if cj is not None:
                    own[cj["ys"], cj["xs"]] = amps[j] * cj["stamp"]
                ix, iy = int(round(r["x"])) - X0, int(round(r["y"])) - Y0
                xa, xb_ = max(ix - h, 0), min(ix + h + 1, data.shape[1])
                ya, yb_ = max(iy - h, 0), min(iy + h + 1, data.shape[0])
                ext = (X0 + xa - .5, X0 + xb_ - .5, Y0 + ya - .5, Y0 + yb_ - .5)
                dwin = data[ya:yb_, xa:xb_]
                rw = (resid + own)[ya:yb_, xa:xb_]
                lim = max(np.nanpercentile(np.abs(rw), 99), 1e-9)
                for ax, img, kw, ttl in (
                        (axRow[0], dwin, dict(cmap=L.CMAP, norm=simple_norm(dwin, stretch="asinh",
                                                                            percent=99.5)), "data"),
                        (axRow[1], rw, dict(cmap="RdBu_r", vmin=-lim, vmax=lim),
                         "residual WITHOUT this source"),
                        (axRow[2], resid[ya:yb_, xa:xb_], dict(cmap="RdBu_r", vmin=-lim, vmax=lim),
                         "residual WITH it")):
                    ax.imshow(img, extent=ext, origin="lower", **kw)
                    ax.plot(r["x"], r["y"], "+", color="k", ms=10, mew=1.5)
                    ax.set_title(ttl, fontsize=9)
                    ax.set_xticks([]); ax.set_yticks([])
                    ax.annotate("", xy=(xy[b, 0], xy[b, 1]), xytext=(r["x"], r["y"]),
                                arrowprops=dict(arrowstyle="->", color="0.4", lw=0.8),
                                annotation_clip=True)
                axRow[0].set_ylabel(f"#{j}  {r['verdict']}\nsignif {r['significance']:.1f}",
                                    fontsize=9)
            fig.suptitle("Sources next to bright stars (least significant first; arrow points to "
                         "the bright star)", fontsize=11)
            plt.tight_layout()
            plt.show()
    return t


def diagnoseResidual(run, brightSnr=300.0, diskFwhm=3.0, plot=True):
    """Decide whether what is left around bright stars is noise or model error.

    For empty sky, bright-star disks (within ``diskFwhm`` FWHM of stars with peak S/N >=
    ``brightSnr``) and other source pixels: ``sigmaZ`` (robust std of residual/sigma; 1 =
    matches the variance plane) and ``coherence`` (robust std of the 3x3-averaged residual/sigma,
    x3, over sigmaZ, squared; white noise gives 1).  Also the pixel-to-pixel scatter against
    model flux, which measures missing photon noise.

    Parameters
    ----------
    run : `lsst.meas.crowd.session.CutoutRun`
        A finished run.
    brightSnr : `float`, optional
        Minimum peak S/N of a bright star.
    diskFwhm : `float`, optional
        Radius of the bright-star disks [FWHM].
    plot : `bool`, optional
        Show residual/sigma, its 3x3 average and its matched-filter S/N.

    Returns
    -------
    result : `dict`
        Per region ``sigmaZ``, ``coherence`` and pixel count; ``highfreq`` (pixel-to-pixel
        scatter vs model flux, with ``k``) and ``noiseShare`` (fraction of the bright-disk excess
        explained by missing photon noise).
    """
    from scipy.ndimage import uniform_filter
    final = run.outerHistory[-1]
    exp, res = final["exposure"], final["result"]
    S = res["S"]
    fit = res["snaps"][-1]["fit"]
    bb = exp.getBBox()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    data, model, var = exp.image.array, res["modelFull"], exp.variance.array
    bad = ((exp.mask.array & run.badBitmask) != 0) | ~np.isfinite(data - model) | ~(var > 0)
    fwhm = L._fwhm(exp.getPsf(), bb.getCenter())
    resid = data - model
    excl = unmodeledNearMask(resid, var, bad, fwhm)
    try:
        P = _params(run.config)
        excl |= satZones(exp, fwhm, P)[1] if P["satZoneFwhm"] else False
    except Exception:
        pass
    good = ~bad & ~excl
    z = np.where(good, resid / np.sqrt(np.where(var > 0, var, 1.0)), 0.0)
    noise = float(np.sqrt(np.nanmedian(var[good])))
    psfPeak = exp.getPsf().computeKernelImage(bb.getCenter()).array.max()
    snr = fit["amps"] * psfPeak / noise
    disks = _starMask(data.shape, X0, Y0, S["xy"][snr >= brightSnr], diskFwhm * fwhm) & good
    sky = good & (model < 0.3 * noise)
    other = good & ~sky & ~disks
    wgt = uniform_filter(good.astype(float), 3)
    zb = np.where(wgt > 0.999, uniform_filter(z, 3) / np.maximum(wgt, 1e-9), np.nan)

    def rs(x):
        x = x[np.isfinite(x)]
        if x.size < 20:
            return np.nan
        q = np.percentile(x, [25, 75])
        return (q[1] - q[0]) / 1.349

    out = {}
    print(f"{'region':22s} {'pixels':>8s} {'sigmaZ':>8s} {'coherence':>10s}")
    def cs(x):                       # clipped plain std: structure lives in the tails too
        x = x[np.isfinite(x)]
        x = x[np.abs(x) < 20]
        return float(np.std(x)) if x.size >= 20 else np.nan

    for name, m in (("empty sky", sky), ("bright-star disks", disks), ("other source pixels", other)):
        sz = rs(z[m])
        sp, sb = cs(z[m]), cs(zb[m])
        coh = (3 * sb / sp) ** 2 if np.isfinite(sp) and np.isfinite(sb) and sp > 0 else np.nan
        out[name] = dict(n=int(m.sum()), sigmaZ=sz, coherence=coh)
        print(f"{name:22s} {int(m.sum()):8d} {sz:8.3f} {coh:10.2f}")
    cSky = out["empty sky"]["coherence"]
    cB, sB = out["bright-star disks"]["coherence"], out["bright-star disks"]["sigmaZ"]
    print("\npixel-to-pixel scatter vs model flux:")
    hf = L.highfreqNoiseCheckArrays(resid, var, model, good, verbose=True)
    out["highfreq"] = hf
    # Decision: how much of the disks' excess variance (sigmaZ^2 - 1) is explained by
    # pixel-to-pixel scatter above the variance plane at the brightest model fluxes?
    hfBright = float(hf["ratio"][-1]) if len(hf["ratio"]) else np.nan
    if np.isfinite(sB) and np.isfinite(hfBright) and sB > 1.02:
        share = max(0.0, hfBright - 1.0) / max(sB ** 2 - 1.0, 1e-9)
        out["noiseShare"] = share
        if share >= 0.5:
            verdict = ("bright-star disks: excess is PIXEL-SCALE NOISE (pixel-to-pixel scatter "
                       "explains it) -- photon noise missing from the variance plane; use "
                       "poissonK=k.")
        elif share >= 0.25:
            verdict = "bright-star disks: MIXED -- part missing photon noise, part model structure."
        else:
            verdict = ("bright-star disks: MODEL STRUCTURE -- the variance plane is right "
                       "(pixel-to-pixel scatter ~ variance), the excess is coherent residual "
                       "the model misses.  See showBrightStarStack(run).")
        print(f"\nbright disks: sigmaZ^2 = {sB ** 2:.2f}; pixel-to-pixel scatter^2/variance at the "
              f"brightest fluxes = {hfBright:.2f}  ->  noise explains {100 * share:.0f}% of the "
              f"excess  ->  {verdict}")
        out["verdict"] = verdict
    if plot:
        import matplotlib.pyplot as plt
        snrmap = ca._matchedSnr(np.where(good, resid, 0.0), np.where(var > 0, var, 1.0), good,
                                fwhm, sky)[0]
        ext = (X0 - .5, X0 + data.shape[1] - .5, Y0 - .5, Y0 + data.shape[0] - .5)
        cm = plt.get_cmap("RdBu_r").copy()
        cm.set_bad(EXCLUDED_COLOR)
        fig, ax = plt.subplots(1, 3, figsize=(21, 6.5), sharex=True, sharey=True)
        panels = ((np.ma.masked_where(~good, z), 4, "residual / σ (per pixel)"),
                  (np.ma.masked_where(~np.isfinite(zb), 3 * zb), 4,
                   "3x3-averaged residual / σ  (x3: noise keeps the same scale)"),
                  (np.ma.masked_where(~good, snrmap), 5,
                   "matched-filter S/N of the residual (only coherent structure survives)"))
        for a, (img, lim, ttl) in zip(ax, panels):
            im = a.imshow(img, extent=ext, origin="lower", cmap=cm, vmin=-lim, vmax=lim)
            plt.colorbar(im, ax=a, fraction=0.046, pad=0.04)
            a.contour(X0 + np.arange(data.shape[1]), Y0 + np.arange(data.shape[0]),
                      disks.astype(float), levels=[0.5], colors="k", linewidths=0.6)
            a.set_title(ttl, fontsize=10)
        fig.suptitle("Noise or structure? (black contours = bright-star disks; grey = excluded)",
                     fontsize=11)
        plt.tight_layout()
        plt.show()
    return out


BAD_PIXEL_COLOR = "#ff00ff"     # magenta: in neither viridis/mako nor RdBu_r


NEARSAT_COLOR = "#ff9900"       # orange: near-saturated cores excluded by option B


EXCLUDED_COLOR_GRAY = "#3399ff"  # blue: same areas when the residual is shown in greyscale


EXCLUDED_COLOR = "#8c8c8c"      # grey: unmodeled flux next to bad pixels (not in either map)


def modeledNear(modeled, bad, fwhm):
    """Return the pixels closer to a modelled (near-saturated) core than to any other masked pixel.

    Residual blobs there are not grown into the excluded area by `unmodeledNearMask`.

    Parameters
    ----------
    modeled : `numpy.ndarray` [`bool`]
        Masked pixels whose star is in the model.
    bad : `numpy.ndarray` [`bool`]
        All masked pixels.
    fwhm : `float`
        PSF FWHM [px] (unused; kept for a uniform signature).

    Returns
    -------
    near : `numpy.ndarray` [`bool`]
        Pixels nearer to a modelled core than to other masked pixels.
    """
    from scipy.ndimage import distance_transform_edt
    if not modeled.any():
        return np.zeros(bad.shape, bool)
    hard = bad & ~modeled
    dm = distance_transform_edt(~modeled)
    dh = distance_transform_edt(~hard) if hard.any() else np.full(bad.shape, np.inf)
    return dm < dh


def unmodeledNearMask(resid, var, bad, fwhm, padFwhm=1.5, growSnr=2.5, reachFwhm=4.0,
                      modeled=None):
    """Find the pixels whose residual cannot be judged because a source there is partly masked.

    Everything within ``padFwhm`` FWHM of a bad pixel, plus any positive residual blob
    (matched-filter S/N > ``growSnr``) touching that pad, out to ``reachFwhm`` FWHM from the
    mask.  Masked cores whose star is in the model (``modeled``, near-saturated cores fitted from
    their wings) get no pad: the residual next to them is real information.

    Parameters
    ----------
    resid : `numpy.ndarray`
        Residual image.
    var : `numpy.ndarray`
        Variance image.
    bad : `numpy.ndarray` [`bool`]
        Bad pixels.
    fwhm : `float`
        PSF FWHM [px].
    padFwhm : `float`, optional
        Pad around bad pixels [FWHM].
    growSnr : `float`, optional
        Matched-filter S/N of residual blobs grown into the excluded area.
    reachFwhm : `float`, optional
        Maximum distance of the excluded area from the mask [FWHM].
    modeled : `numpy.ndarray` [`bool`], optional
        Masked pixels whose star is in the model.

    Returns
    -------
    excluded : `numpy.ndarray` [`bool`]
        Pixels to exclude (not including the bad pixels themselves).
    """
    from scipy.ndimage import distance_transform_edt, label
    if not bad.any():
        return np.zeros(bad.shape, bool)
    modeled = (modeled & bad) if modeled is not None else np.zeros(bad.shape, bool)
    hard = bad & ~modeled
    dist = distance_transform_edt(~bad)
    pad = ~bad & (distance_transform_edt(~hard) <= padFwhm * fwhm) if hard.any() \
        else np.zeros(bad.shape, bool)

    usable = ~bad & np.isfinite(resid) & np.isfinite(var) & (var > 0)
    snr = ca._matchedSnr(np.where(usable, resid, 0.0), np.where(usable, var, 1.0), usable, fwhm,
                         usable & ~pad)[0]
    lab, _ = label((snr > growSnr) & usable & ~modeledNear(modeled, bad, fwhm)
                   & (dist <= reachFwhm * fwhm))
    touch = np.unique(lab[pad & (lab > 0)])
    return pad | (np.isin(lab, touch[touch > 0]) & (lab > 0))


def runSessionGaia(session, config=None, names=None, outDir=None, viz=False, plotFinal=True,
                   gmagMax=None, nJobs=1, plotKwargs=None, **overrides):
    """Run the adaptive pipeline on every (or the named) cutout of a session.

    Calls `lsst.meas.crowd.gaiaCatalog.attachGaia` first if the session has no Gaia catalogue.

    Parameters
    ----------
    session : `lsst.meas.crowd.session.Session`
        The session.
    config : `dict`, optional
        Pipeline configuration.
    names : `list` [`str`], optional
        Cutouts to run (default: all).
    outDir : `str`, optional
        Save each run under ``<outDir>/visit<V>_det<D>/<name>/``; needed with ``nJobs > 1``,
        since results that cannot be sent back are only on disk.
    viz : `bool`, optional
        Draw the per-pass overview.
    plotFinal : `bool`, optional
        Plot the final comparison of each cutout.
    gmagMax : `float`, optional
        Only start from Gaia stars brighter than this G magnitude.
    nJobs : `int`, optional
        Fit cutouts in this many parallel worker processes (Linux fork; one core each).
    plotKwargs : `dict`, optional
        Keyword arguments of `lsst.meas.crowd.plotting.plotFinalComparison`.
    **overrides
        Configuration overrides (e.g. ``maxOuter=3``).

    Returns
    -------
    runs : `dict` [`str`, `lsst.meas.crowd.session.CutoutRun`]
        The results by cutout name.
    """
    from . import session as cn
    import matplotlib.pyplot as plt
    from .plotting import plotFinalComparison
    if getattr(session, "gaia", None) is None:
        try:
            attachGaia(session)
        except Exception as e:
            print(f"Gaia unavailable ({e}) -- running WITHOUT Gaia starting positions (adaptive loops still on)")
            session.gaia = None
    runs = {}
    todo = names or list(session.jobs)
    if nJobs and nJobs > 1 and len(todo) > 1:
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor, as_completed
        if not outDir:
            print("NOTE: set outDir with nJobs>1 so results are on disk even if they cannot "
                  "be sent back to the notebook")
        _PAR["args"] = (session, config, outDir, gmagMax, overrides)
        t0 = time.perf_counter()
        with ProcessPoolExecutor(min(nJobs, len(todo)), mp_context=mp.get_context("fork")) as ex:
            futs = {ex.submit(_parWorker, nm): nm for nm in todo}
            for f in as_completed(futs):
                nm = futs[f]
                try:
                    nm, run, summ, path = f.result()
                except Exception:
                    print(f"{nm} FAILED:")
                    traceback.print_exc()
                    continue
                print(f"{nm}: {summ['nSources']} sources ({summ['nFromGaia']} Gaia), residual "
                      f"{'clean' if summ['residualClean'] else 'still shows excess'}, "
                      f"{summ['wallSeconds']:.0f}s")
                if run is not None:
                    runs[nm] = run
                    if plotFinal:
                        plotFinalComparison(run, **(plotKwargs or {}))
                        plt.show()
                elif plotFinal and path:
                    cn.plotFinalComparisonFromFits(path)
        print(f"all {len(todo)} cutouts: {time.perf_counter() - t0:.0f}s wall")
        return runs
    for k, name in enumerate(todo, 1):
        job = session.jobs[name]
        print(f"{'#' * 78}\n[{k}/{len(todo)}] {name}  interior={job['interior']}")
        try:
            run = runCutoutGaia(job["raw"], session.psf, session.badBitmask, session.gaia,
                                config, interiorBbox=job["interior"], starCat=session.starCat,
                                bgImage=session.bgImage, viz=viz, name=name, gmagMax=gmagMax,
                                **overrides)
            runs[name] = run
            s = run.summary
            print(f"{name}: {s['nSources']} sources ({s['nFromGaia']} Gaia, "
                  f"{s['nInInterior']} in interior), residual "
                  f"{'clean' if s['residualClean'] else 'still shows excess'}, "
                  f"{s['wallSeconds']:.0f}s")
            if outDir:
                cn.saveRun(run, os.path.join(outDir, f"visit{session.visit}_det{session.detector}",
                                             name))
            if plotFinal:
                plotFinalComparison(run, **(plotKwargs or {}))
                plt.show()
        except Exception:
            traceback.print_exc()
    return runs


def report(run):
    """Print what the adaptive loops did, pass by pass.

    Parameters
    ----------
    run : `lsst.meas.crowd.session.CutoutRun`
        A finished run.
    """
    hist = getattr(run, "outerHistory", run)
    for i, h in enumerate(hist, 1):
        ad = h["result"]["adaptive"]
        res = h["result"]
        print(f"outer {i}: {h['nPeaks']} sources ({int(res['isGaia'].sum())} Gaia), "
              f"{ad['nInner']} inner it., PSF change into pass "
              f"{100 * h['psfChange']:.2f}%, {ad.get('seconds', np.nan):.0f}s, stop: {ad['stop']}")
        for r in ad["log"]:
            print(f"    it {r['it'] + 1:2d}  thr {r['threshold']:4.2f}σ  N={r['n']:5d}  "
                  f"+{r['added']:<4d} peaks +{r['pos']}/-{r['neg']} (excess {r['excess']:+d})  "
                  f"{r['action']}")
        if ad["final"] is not None:
            ca._printDiag(f"  [outer {i}]", ad["final"])
    last = hist[-1]["result"]["adaptive"]
    if last["missing"] and last["final"] is not None:
        print("top remaining positive residual peaks (x, y, S/N, px to nearest fitted):")
        for x, y, s, d in [c for c in last["final"]["candidates"] if c[2] >= last["target"]][:10]:
            print(f"    ({x:7.1f}, {y:7.1f})  S/N {s:5.1f}  {d:5.1f}")


# the names used before the "Gaia starting positions" naming, so existing notebooks keep working
SEEDED_DEFAULTS = PIPELINE_DEFAULTS
runSeededPipeline = runAdaptivePipeline
runInnerSeeded = runInnerAdaptive


# old snake_case names still work (with a FutureWarning)
from ._compat import moduleGetattr  # noqa: E402
__getattr__ = moduleGetattr(globals())
