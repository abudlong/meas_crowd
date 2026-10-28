# This file is part of meas_crowd.
#
# Developed for the LSST Data Management System.
# See the COPYRIGHT file at the top-level directory of this distribution.

"""Gaia-seeded, adaptive, iterative joint PSF-fitting pipeline (inner / outer loops, saturated and
bright-star handling, PSF refinement, extended sources, residual diagnostics).
"""

import math
import os
import time
import traceback
from collections import OrderedDict, defaultdict
import numpy as np
from scipy.linalg import cho_factor, solve_triangular
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from astropy.table import Table
from . import fitLib as L
from . import adaptive as ca
from .gaiaSeeds import (
    attach_gaia,
    select_seeds,
)


# Options of this driver on top of the run_cutout_pipeline config keys.
SEEDED_DEFAULTS = dict(
    # ---- adaptive inner loop
    max_inner=30,              # safety cap on inner iterations per outer pass
    seed_threshold_start=5.0,  # first detection threshold when Gaia seeds are used
    threshold_start=None,      # without seeds (None: inner_kwargs['threshold_start'])
    threshold_step=None,       # None: inner_kwargs['threshold_step']
    threshold_min=None,        # None: inner_kwargs['threshold_min']
    grow_frac=0.02,            # "still finding sources": added >= grow_frac * N ...
    grow_min=3,                # ... or >= grow_min sources -> hold the threshold
    patience=2,                # stop after this many iterations without progress
    detect_gain=1.0,           # detection runs on data - gain*model (1 = full joint model)
    target_snr=None,           # residual must be clean down to this S/N (None: min_snr)
    t_floor=3.0,               # never detect below this S/N
    fdr=0.25,                  # tolerated N_neg/N_pos when going deeper
    z_sig=2.0,                 # significant excess: >= z_sig*sqrt(N_pos+N_neg) ...
    min_excess=3,              # ... and >= this many sources
    edge_fwhm=2.0, halo_fwhm=2.5, min_sep_fwhm=0.75,   # residual-check clean region
    centroid_rounds_after_first=1,
    hopeless_frac=0.5,         # eliminate_weak: drop S/N < hopeless_frac*min_snr at once
    elim_drop_frac=0.25,
    # ---- adaptive outer loop
    min_outer=2,               # at least this many passes (pass 2 refines PSF + shape cut)
    max_outer=6,               # safety cap
    outer_growth_tol=0.01,     # source count changed by < 1% ...
    psf_tol=0.01,              # ... and PSF kernel changed by < 1% (L2)  -> converged
    # ---- seeds
    protect_seeds=True,
    seed_pad_fwhm=1.5,         # keep Gaia seeds up to this far OUTSIDE the cutout
    seed_max_drift_px=1.5,     # max distance a seed may move from its Gaia position
    seed_min_sep_px=0.5,       # Gaia duplicates closer than this are merged
    fit_gaia_offset=True,      # remove a global WCS-vs-data offset first
    subtract_initial_bg=True,  # subtract session.bgImage in pass 1
    keep_all_fits=False,       # keep full fits (all PSF stamps) in every snapshot
    # ---- saturated / bright stars
    sat_grow_px=2,             # also exclude pixels within this many px of SAT (non-linear,
                               # charge spill); 0 = off
    bright_corr=True,          # per-bright-star shift + width + flux correction terms
    bright_corr_sigma=100.0,   # ... for stars with amp*psf_peak/noise above this (or SAT)
    bright_corr_radius_fwhm=4.0,
    bright_corr_radial=True,   # ... plus radial-profile terms (core/wing shape mismatch)
    bright_core_rings_px=(0.0, 1.0, 2.0, 3.0, 4.0),  # bright (unsaturated) stars also get NARROW
    bright_core_ring_sigma=0.6,  #   radial terms (a central spike + thin rings at these radii [px],
                               #   Gaussian width this many px): the core of bright stars is sharper
                               #   than the PSF on a 1-2 px scale (spike / dip / ring residual)
    bright_core_modes=2,       # ... and each narrow term may be lopsided / elongated (angular modes
                               #   m <= this): the sharp core can sit a fraction of a pixel off the
                               #   rest of the star (core dipoles)
    bright_angular_modes=4,    # unsaturated bright stars: rings (<= 4 FWHM) may vary with
                               #   angle (cos/sin m*theta, m = 1..this) -- removes the lopsided
                               #   'dipole' residuals a shift term cannot.  0 = off
    sat_halo_modes=4,          # saturated stars: the halo rings may also vary with angle
                               #   (cos/sin m*theta, m = 1..this) -- spikes, bleed charge and
                               #   lopsided halos are not circular.  0 = radial only
    # ---- saturated-star zones
    sat_zone_fwhm=2.0,         # zone = pixels within sat_zone_fwhm*FWHM + sat_zone_scale*(half-
    sat_zone_scale=1.0,        #   width of the saturated region) of the SAT mask -- a band that
                               #   follows the mask's shape (bleed trails stay narrow)
    sat_zone_snr=10.0,         # non-Gaia sources inside a zone need S/N >= this (instead of
                               #   min_snr): real stars survive, halo/spike speckle does not.
                               #   The residual check ignores the zones. sat_zone_fwhm=0 -> off
    sat_zone_min_px=4,         # ignore saturated blobs smaller than this (hot pixels, CRs)
    psf_spatial=False,         # let the PSF vary LINEARLY across the cutout (kernel gradient in x
                               #   and y, measured from the clean stars); kept only if it lowers
                               #   chi2 on those stars by >= psf_spatial_min_gain
    psf_spatial_min_gain=0.02,
    psf_spatial_min_stars=10,
    psf_cv_min_gain=0.005,     # PSF updates are CROSS-VALIDATED: built from half the clean stars,
                               #   judged on the other half (both ways); kept only if the held-out
                               #   chi2 drops by >= this fraction (stops the PSF fitting noise)
    psf_smooth_px=0.6,         # Gaussian smoothing [px] of every PSF correction (a real PSF
                               #   change is smooth on the PSF scale; pixel noise is not)
    psf_soft_sigma=2.0,        # corrections below this many sigma per pixel are ignored
    psf_stack_n=25,            # PSF update between outer passes: up to this many CLEAN stars
    psf_stack_min_snr=50.0,    #   (peak S/N >= this, no saturated pixel within 3 FWHM, no
    psf_stack_iso_fwhm=2.5,    #   brighter-than-10% neighbour within psf_stack_iso_fwhm*FWHM,
                               #   not extended).  Saturated stars used to dominate the stack.
    mask_pad_fwhm=1.5,         # residual check ignores pixels within this many FWHM of bad
                               #   pixels (flux of sources whose centre is masked)
    sat_refine_wings=True,     # saturated stars without a Gaia position: refine the position
                               #   from the unsaturated wings (the core cannot be used)
    wing_keep_angular=False,   # wing extensions keep the kernel's angular pattern (spikes) beyond
                               #   the blend radius.  OFF: the visit kernel's outer part is mostly
                               #   noise and the pattern imprinted noisy rings on bright stars
    psf_core_median=False,     # (experimental, off: did not help on real data) after each outer pass: correct the PSF CORE with the MEDIAN residual
    psf_core_min_snr=100.0,    #   of the bright clean stars (peak S/N >= this; what
    psf_core_radius_fwhm=2.5,  #   show_bright_star_stack shows), within this many FWHM, tapered,
    psf_core_cycles=2,         #   cross-validated on two halves of the stars
    bf_term=False,             # (experimental, off: did not help on real data) BRIGHTNESS-DEPENDENT PSF (brighter-fatter): every star's PSF is
    bf_min_snr=30.0,           #   K + (log10 flux - L0) * B(r), with the radial profile B measured by
    bf_radius_fwhm=3.0,        #   regressing isolated stars' residual profiles on log flux (S/N >=
                               #   bf_min_snr, within bf_radius_fwhm*FWHM); flux-neutral, smooth,
                               #   cross-validated on two halves of the stars
    model_err_measure=True,    # (A) after each outer pass, MEASURE the PSF model error eps: the
    model_err_min_snr=50.0,    #   residual variance around stars (S/N >= this) beyond the pixel
                               #   variance, as (eps * model)^2.  plot_final_comparison then shows
                               #   residual / sqrt(var + (eps*model)^2) (sigma="model", default)
    fit_sys_from_data=False,   # ... and also use eps as sys_frac in the FIT weights / residual check
                               #   (default False: they keep sys_frac=0.03, which is more tolerant)
    near_sat_frac=0.0,         # (B, off by default; e.g. 0.7 to use) pixels above this fraction of the saturation level (not flagged
    near_sat_level=None,       #   SAT) are excluded like saturated ones: those stars are fitted
    near_sat_grow_px=1,        #   from their wings.  near_sat_level: saturation level in image
                               #   units (None = estimated from the SAT pixels).  0/None frac = off
    psf_regularize=True,       # after EVERY PSF update: beyond psf_reg_core_fwhm*FWHM the kernel is
    psf_reg_core_fwhm=1.5,     #   rebuilt as a smooth, monotone radial profile (ring medians) times
    psf_reg_sigma=1.0,         #   a smoothed angular pattern -- no holes or isolated hot pixels
    psf_ext_blend=(12.0, 16.0),  # initial 41->81 px wing extension: the kernel is replaced by a
                               #   fitted Moffat between these radii [px]; e.g. (17, 19.5) keeps
                               #   the PSF model's own profile almost to the kernel edge
    wing_blend_fwhm=(2.0, 3.0),  # measured wing profile replaces the kernel between these radii
                               #   [FWHM]; try (3.5, 5.0) to keep more of the PSF model's own wing
    wings_from_stars=True,     # measure the PSF's outer wings from the bright + saturated stars'
    wing_min_snr=300.0,        #   radial profiles (unsaturated stars with peak S/N >= this set
                               #   the shape; saturated stars extend it outwards) instead of
                               #   extrapolating the kernel's noisy 6-14 px profile
    bg_sat_mask_fwhm=6.0,      # background refit ignores pixels within this many FWHM of
                               #   saturated pixels (their halos were absorbed into the sky)
    drop_halo_mimics=True,     # end of each outer pass: remove non-Gaia sources next to bright
                               #   stars that check_bright_neighbours() calls 'doubtful', and
                               #   do not re-detect them in later passes
    halo_check_snr=300.0,      # ... bright = peak S/N >= this
    halo_check_fwhm=4.0,       # ... neighbours within this many FWHM
    halo_similar_min=3,        # ... 'doubtful' if >= this many control angles show >= half
                               #   its flux (a repeating halo pattern), or significance < 3
    auto_poisson=True,         # after each outer pass: if the residual on bright stars is
                               #   PIXEL-SCALE noise (diagnose_residual), measure the photon noise
                               #   the variance plane is missing (k per unit model flux) and use
                               #   var + k*model from then on (fit weights, S/N, residual check)
    one_source_per_sat=True,   # keep only the brightest source on each saturated core (its
                               #   masked centre cannot separate 2-3 sources; they split the
                               #   flux and leave blotches) -- also applies to Gaia entries
    # ---- extended sources / shape cut
    keep_extended=True,        # a source the shape cut finds too BROAD is kept with a
                               # Gaussian-broadened PSF (small galaxy) instead of deleted,
                               # so its flux does not stay in the residual
    shape_snr_test=False,      # the shape cut's own local-S/N test (eliminate_weak already
                               # enforces min_snr; the two disagreeing caused churn)
    later_shift_cap=0.5,       # centroid step cap after the first iteration of a pass
    drop_ring_mimics=True,     # end of every inner loop: a faint NON-Gaia source near a bright star
    ring_check_snr=150.0,      #   (peak S/N >= ring_check_snr, within ring_check_fwhm*FWHM, flux <
    ring_check_fwhm=4.0,       #   ring_max_ratio of the bright star) must stand out from the
    ring_max_ratio=0.05,       #   AZIMUTHAL MEDIAN of the residual at its distance from the star
    ring_snr=5.0,              #   with S/N >= ring_snr -- otherwise it is a fake absorbing the
                               #   star's wing excess: removed and not re-detected
    ring_bg_modes=4,           #   the wing level at each radius includes angular modes m <= this,
                               #   not just the azimuthal median
    ring_sym=True,             #   ... and the source must also stand out from the same spot rotated
                               #   by 90/180/270 deg about the star (4-fold diffraction pattern):
                               #   a real star is in one place, a spike/lobe is in all four
    split_test=True,           # end of every inner loop: a source (S/N split_min_snr .. bright_corr
    split_min_snr=20.0,        #   _sigma) with a positive residual peak (matched-filter S/N >=
    split_peak_snr=5.0,        #   split_peak_snr) 0.5-1.5 FWHM away gets a trial second source
    split_dchi2=25.0,          #   there; kept if chi2 within 2 FWHM drops by >= split_dchi2 and its
                               #   own S/N >= min_snr (unresolved pairs -> dipoles)
    ellip_test=True,           # end of every inner loop: a source (Gaia too) with S/N ellip_min_snr ..
    ellip_min_snr=15.0,        #   bright_corr_sigma and local chi2/px >= ellip_chi2_min is refitted,
    ellip_chi2_min=1.3,        #   together with the non-Gaia points within ellip_group_fwhm*FWHM, as
    ellip_group_fwhm=1.5,      #   ONE elliptical extended source (PSF x elongated Gaussian); kept if
    ellip_dchi2=10.0,          #   the local chi2 drops by >= ellip_dchi2 (+5 per extra fitted number:
                               #   a LONE point source -> ellipse adds 3, so it needs >= 25)
    final_polish=True,         # end of every inner loop (after all removals): re-optimise ALL
    final_polish_rounds=4,     #   positions jointly (Gauss-Newton steps capped at final_polish_cap
    final_polish_cap=0.3,      #   px, Gaia drift limit kept), kept only while the total chi2
                               #   drops -- removes dipoles from positions left stale when
                               #   neighbours were added / removed
    remember_rejects=True,     # do not re-add candidates already rejected in this pass
    # ---- speed
    fast_stamps=True,          # render PSF stamps / derivatives with the C++ kernels (NumPy if not
                               #   built) for spatially CONSTANT PSFs; each PSF is first checked
                               #   against psf.computeImage and the fast path is only used if they
                               #   agree to 2e-3 of the peak
    stamp_floor_sigma=0.05,    # crop each source's PSF stamp where its wing drops below this
                               # many sigma per pixel (faint sources get small stamps);
                               # None = always the full psf_wing_size stamp
)


# run_cutout_pipeline defaults for the keys this driver uses
_LIB_DEFAULTS = dict(
    inner_kwargs=None, support_frac=0.01, refine_centroids=True, centroid_min_amp=50.0,
    n_centroid_rounds=4, centroid_shift_tol=0.03, centroid_shift_schedule=None,
    coarse_centroid=True, coarse_centroid_hw_sigma=2.5, merge_separation_sigma=1.5,
    sys_frac=0.03, min_snr=5.0, prune_sep_sigma=1.5, reject_fwhm=0.75,
    bright_peak_sigma=100.0, freeze_fwhm=1.0, poisson_k=0.0, near_bright_fwhm=2.0,
    near_bright_snr=10.0, edge_margin_fwhm=0.5, use_eliminate_weak=True,
    use_shape_cut=True, shape_kwargs=None, shape_in_first_pass=False,
    refit_psf=True, stack_kwargs=None, psf_polish_cycles=3, psf_kernel_size=41,
    psf_kernel_sizes=None, bg_bin_size=32, bg_use_approx=True, psf_wing_size=81,
)


_IGNORED = {"n_outer_iterations", "n_inner_iterations", "outer_convergence_frac",
            "plot_every", "source_view", "max_pages_per_iter", "sources_per_page",
            "view_half", "x_range", "y_range", "viz"}


def _params(cfg):
    P = dict(_LIB_DEFAULTS)
    P.update(SEEDED_DEFAULTS)
    unknown = set(cfg) - set(P) - _IGNORED
    if unknown:
        raise ValueError(f"unknown option(s) {sorted(unknown)}")
    P.update({k: v for k, v in cfg.items() if k in P})
    ik = dict(n_iterations=10, threshold_start=8.0, threshold_step=-1.5, threshold_min=4.5,
              min_pixels=1)
    ik.update(P["inner_kwargs"] or {})
    for k in ("threshold_start", "threshold_step", "threshold_min"):
        if P[k] is None:
            P[k] = float(ik[k])
    P["min_pixels"] = int(ik.get("min_pixels", 1))
    if P["threshold_step"] >= 0:
        P["threshold_step"] = -abs(P["threshold_step"]) or -1.0
    P["t_floor"] = min(P["t_floor"], P["threshold_min"])
    return P


# =====================================================================
# Speedups (patched into cutout_fit_lib)
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


def _set_bf_fluxes(xy, amps):
    """Tell the stamp builder each source's current flux (for the brightness term)."""
    if _BF["prof"] is None or len(xy) == 0:
        return
    a = np.asarray(amps, float)
    ok = a > 0
    _BF.update(tree=cKDTree(np.asarray(xy)[ok]) if ok.any() else None, lf=np.log10(a[ok]))


def _bf_reset():
    _BF.update(r=None, prof=None, L0=None, tree=None, lf=None)
    _ELLIP.clear()
    _STAMP_CACHE.clear()


_ELLIP = {}       # source uid -> (axis ratio q <= 1, position angle [rad]) of an elliptical extended source


def _ellip_kernel(sa, q, th, prof="gauss"):
    """Normalised elliptical profile: 'gauss' (major sigma sa) or 'exp' (exponential disk, major
    scale length sa); minor axis q*sa, major axis at angle th.  The exponential is integrated over
    5x5 sub-pixels (its centre is a cusp)."""
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


def _set_ext_rule(S):
    m = S["ext"] > 0
    ell = [tuple(_ELLIP.get(int(u), (1.0, 0.0, "gauss"))) + ("gauss",) for u in S["uid"][m]]
    qs = np.array([e[0] for e in ell])
    ths = np.array([e[1] for e in ell])
    prs = [e[2] for e in ell]
    _EXT_RULE.update(tree=cKDTree(S["xy"][m]) if m.any() else None, sig=S["ext"][m], q=qs, th=ths,
                     prof=prs)


def clear_stamp_cache():
    _STAMP_CACHE.clear()
    _SIZE_RULE["tree"] = None


def _crop_stamp(arr, fx0, fy0, x, y, half):
    """Zero the stamp beyond radius `half` around (x, y) with a cosine taper over the outer
    quarter.  The array keeps its shape (callers rely on it); zeros never enter the fit."""
    yy, xx = np.mgrid[fy0:fy0 + arr.shape[0], fx0:fx0 + arr.shape[1]]
    w = np.clip((half - np.hypot(xx - x, yy - y)) / (0.25 * half), 0.0, 1.0)
    return arr * (0.5 - 0.5 * np.cos(np.pi * w))


def set_stamp_sizes(xy, amps, psf, point, noise, fwhm, floor_sigma=0.05, tol=1.5):
    """Choose a stamp half-width per source: the smallest of a few tiers beyond which
    amp * PSF stays below floor_sigma * noise.  Sources not in this list (new detections)
    get the full stamp.  Returns the half-widths (-1 = full)."""
    xy = np.asarray(xy, float).reshape(-1, 2)
    if floor_sigma is None or not len(xy):
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
            half[s * env[t] < floor_sigma] = t
    _SIZE_RULE.update(tree=cKDTree(xy), half=half, tol=tol)
    return half


def _psf_stamp_cached(psf, x, y, clip_bbox, taper=None):
    """Same result as cutout_fit_lib.psf_stamp, cached on the exact (x, y)."""
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
        arr, bb = _raw_psf_stamp(psf, x, y, taper)
        if eps:
            yy, xx = np.mgrid[bb.getMinY():bb.getMinY() + arr.shape[0],
                              bb.getMinX():bb.getMinX() + arr.shape[1]]
            add = eps * np.interp(np.hypot(xx - x, yy - y), _BF["r"], _BF["prof"], right=0.0)
            arr = np.clip(arr + add, 0.0, None)
        if half > 0:
            arr = _crop_stamp(arr, bb.getMinX(), bb.getMinY(), x, y, half)
        if sig > 0:
            tot = arr.sum()
            if eq < 0.999 or epr != "gauss":
                from scipy.signal import fftconvolve
                arr = fftconvolve(arr, _ellip_kernel(sig, eq, eth, epr), mode="same")
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
    bx0, by0 = max(fx0, clip_bbox.getMinX()), max(fy0, clip_bbox.getMinY())
    bx1, by1 = min(fx1, clip_bbox.getMaxX()), min(fy1, clip_bbox.getMaxY())
    if bx0 > bx1 or by0 > by1:
        return None, L.geom.Box2I()
    bb = L.geom.Box2I(L.geom.Point2I(bx0, by0), L.geom.Point2I(bx1, by1))
    return arr[by0 - fy0:by1 - fy0 + 1, bx0 - fx0:bx1 - fx0 + 1], bb


def _inv_diag(M):
    """diag(M^-1) of a symmetric positive (semi)definite matrix."""
    try:
        c, low = cho_factor(M, lower=True, check_finite=False)
        Linv = solve_triangular(np.tril(c), np.eye(len(M)), lower=True, check_finite=False)
        return np.sum(Linv ** 2, axis=0)
    except (np.linalg.LinAlgError, ValueError):
        return np.diag(np.linalg.pinv(M))


def amp_snr_full_fast(fit, noise_scale=1.0):
    """Drop-in for cutout_fit_lib.amp_snr_full: S/N = amp / sqrt(diag((A^T W A)^-1)),
    computed block by block (independent groups of overlapping sources) with Cholesky
    instead of one dense pinv of the whole matrix."""
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
    var_amp = np.full(n, np.inf)
    single = (counts[lab] == 1) & (d > 0)
    var_amp[single] = 1.0 / d[single]
    for c in np.flatnonzero(counts > 1):
        ii = order[starts[c]:starts[c] + counts[c]]
        var_amp[ii] = _inv_diag(N[ii][:, ii].toarray())
    sig = np.sqrt(np.clip(var_amp, 1e-30, None)) * noise_scale
    return np.where(np.isfinite(sig), fit["amps"] / sig, 0.0)


NNLS_MIN_SOURCES = 300            # below this lsq_linear is already fast (often 0 iterations)


def _nnls_normal(A, b):
    """Exact bounded (x >= 0) least squares via the normal equations: with N = A^T A = R^T R,
    ||Ax - b||^2 = ||Rx - c||^2 + const (c = R^-T A^T b), so NNLS runs on an n x n problem
    instead of the (pixels x n) one."""
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


def _lsq_linear_fast(A, b, bounds=(-np.inf, np.inf), **kw):
    """Patched into cutout_fit_lib: large non-negative joint fits use _nnls_normal."""
    from scipy.optimize import OptimizeResult
    from scipy.sparse import issparse
    try:
        is_nonneg = tuple(bounds) == (0, np.inf)
    except Exception:
        is_nonneg = False
    if A.shape[1] < NNLS_MIN_SOURCES or not issparse(A) or not is_nonneg:
        return _ORIG["lsq_linear"](A, b, bounds=bounds, **kw)
    try:
        x = _nnls_normal(A, b)
    except Exception:
        return _ORIG["lsq_linear"](A, b, bounds=bounds, **kw)
    res = A @ x - b
    return OptimizeResult(x=x, status=1, nit=0, cost=0.5 * float(res @ res), success=True)


def profile(fn, *args, top=30, **kw):
    """Run fn(*args, **kw) under cProfile and print where the time goes, e.g.
        runs = cg.profile(cg.run_session_gaia, sess, cfg, max_outer=2)
    Returns fn's result."""
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


def _fast_for(psf):
    """FastPsf for a spatially constant PSF that matches psf.computeImage, else None (cached)."""
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
            arr, bb = _ORIG["psf_stamp"](psf, x, y, _HUGE, taper=None)
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
    return L.geom.Box2I(L.geom.Point2I(int(b[0]), int(b[1])), L.geom.Point2I(int(b[2]), int(b[3])))


def _raw_psf_stamp(psf, x, y, taper):
    """Full (unclipped) PSF stamp at (x, y): fast path for constant PSFs, else psf.computeImage."""
    fast = _fast_for(psf)
    if fast is not None:
        arr, b = fast.stamp(x, y, (-10 ** 8, -10 ** 8, 10 ** 8, 10 ** 8), taper)
        return arr, _box(b)
    return _ORIG["psf_stamp"](psf, x, y, _HUGE, taper=taper)


def _psf_derivatives_fast(psf, x, y, bbox_clip, eps=0.25):
    """Drop-in for fitLib.compute_psf_derivatives using the fast path when possible."""
    fast = _fast_for(psf)
    if fast is None:
        return _ORIG["compute_psf_derivatives"](psf, x, y, bbox_clip, eps)
    clip = (bbox_clip.getMinX(), bbox_clip.getMinY(), bbox_clip.getMaxX(), bbox_clip.getMaxY())
    dx, dy, b = fast.derivatives(x, y, clip, eps)
    if dx is None:
        return None, None, L.geom.Box2I()
    return dx, dy, _box(b)


def enable_speedups():
    """Patch cutout_fit_lib with the cached psf_stamp and the fast honest S/N (idempotent)."""
    global _HUGE
    if "psf_stamp" not in _ORIG:
        _ORIG["psf_stamp"] = L.psf_stamp
        _ORIG["amp_snr_full"] = L.amp_snr_full
        _ORIG["lsq_linear"] = L.lsq_linear
        _ORIG["compute_psf_derivatives"] = L.compute_psf_derivatives
        _HUGE = L.geom.Box2I(L.geom.Point2I(-10 ** 8, -10 ** 8), L.geom.Extent2I(2 * 10 ** 8, 2 * 10 ** 8))
    L.psf_stamp = _psf_stamp_cached
    L.amp_snr_full = amp_snr_full_fast
    L.lsq_linear = _lsq_linear_fast
    L.compute_psf_derivatives = _psf_derivatives_fast


def disable_speedups():
    if "psf_stamp" in _ORIG:
        L.psf_stamp = _ORIG["psf_stamp"]
        L.amp_snr_full = _ORIG["amp_snr_full"]
        L.lsq_linear = _ORIG["lsq_linear"]
        L.compute_psf_derivatives = _ORIG["compute_psf_derivatives"]
    _FAST["cache"].clear()
    clear_stamp_cache()


# =====================================================================
# Source state: positions + seed bookkeeping, kept aligned through every filter
# =====================================================================
def _state(xy=None, seed=False, gaia=None, gmag=None):
    xy = np.asarray(xy if xy is not None else np.zeros((0, 2)), float).reshape(-1, 2)
    n = len(xy)
    g = np.full(n, np.nan) if gmag is None else np.asarray(gmag, float).reshape(n)
    return dict(xy=xy.copy(), seed=np.full(n, bool(seed)), ext=np.zeros(n), g=g,
                ref=xy.copy() if seed else np.full((n, 2), np.nan),
                gaia=(np.asarray(gaia, np.int64) if gaia is not None else np.full(n, -1, np.int64)),
                uid=np.arange(n, dtype=np.int64), next_uid=n)


def _take(S, keep):
    keep = np.asarray(keep)
    out = {k: S[k][keep] for k in ("xy", "seed", "ext", "g", "ref", "gaia", "uid")}
    out["next_uid"] = S["next_uid"]
    return out


def _add(S, xy_new):
    xy_new = np.asarray(xy_new, float).reshape(-1, 2)
    m = len(xy_new)
    out = dict(xy=np.vstack([S["xy"], xy_new]), seed=np.concatenate([S["seed"], np.zeros(m, bool)]),
               ext=np.concatenate([S["ext"], np.zeros(m)]),
               g=np.concatenate([S["g"], np.full(m, np.nan)]),
               ref=np.vstack([S["ref"], np.full((m, 2), np.nan)]),
               gaia=np.concatenate([S["gaia"], np.full(m, -1, np.int64)]),
               uid=np.concatenate([S["uid"], S["next_uid"] + np.arange(m)]),
               next_uid=S["next_uid"] + m)
    return out


def _xy(S):
    return [tuple(p) for p in S["xy"]]


def _clamp_seeds(S, max_drift):
    if not max_drift or not S["seed"].any():
        return
    s = S["seed"]
    d = S["xy"][s] - S["ref"][s]
    r = np.hypot(d[:, 0], d[:, 1])
    over = r > max_drift
    d[over] *= (max_drift / r[over])[:, None]
    S["xy"][s] = S["ref"][s] + d


def estimate_gaia_offset(exposure, S, fwhm, bad_bitmask, noise, min_peak_snr=20.0,
                         min_stars=5, verbose=True):
    """Median (data centroid - Gaia position) of bright, unsaturated, isolated-ish seeds.
    Uses a background-subtracted windowed centroid. Returns (dx, dy, n_used)."""
    data = exposure.image.array
    bad = (exposure.mask.array & bad_bitmask) != 0
    bb = exposure.getBBox()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    ny, nx = data.shape
    hw = max(3, int(round(1.25 * fwhm)))
    d = []
    others = cKDTree(S["xy"]) if len(S["xy"]) > 1 else None
    for x, y in S["xy"][S["seed"]]:
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
        if w.max() < min_peak_snr * noise or w.sum() <= 0:
            continue
        yy, xx = np.mgrid[-hw:hw + 1, -hw:hw + 1]
        cx = (xx * w).sum() / w.sum() + ix + X0
        cy = (yy * w).sum() / w.sum() + iy + Y0
        d.append((cx - x, cy - y))
    if len(d) < min_stars:
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
# Index-preserving versions of the cleaning steps (seeds protected)
# =====================================================================
def _fit(exposure, psf, S, bad_bitmask, var_model, fitkw):
    _set_ext_rule(S)
    return L.joint_psf_fit_diag(exposure, psf, _xy(S), bad_bitmask, var_model=var_model, **fitkw)


def _refine_subset(exposure, psf, S, fit, idx, max_shift, min_amp, freeze_sep, psf_sigma):
    """One Gauss-Newton position step for the sources `idx`, against the CURRENT joint
    residual (no model rebuild, unlike refine_peak_positions).  Sources fainter than
    min_amp, extended ones, and ones with a brighter neighbour within freeze_sep stay put."""
    eb = exposure.getBBox()
    data, var, good, amps = exposure.image.array, fit["variance"], fit["good"], fit["amps"]
    resid = data - fit["model_full"]
    xy = S["xy"]
    new = xy.copy()
    tree = cKDTree(xy) if len(xy) > 1 else None
    R2 = (2.5 * psf_sigma) ** 2
    for j in idx:
        a = amps[j]
        if a < min_amp or S["ext"][j] > 0:
            continue
        x, y = xy[j]
        if tree is not None and freeze_sep:
            if any(k != j and amps[k] > a for k in tree.query_ball_point(xy[j], freeze_sep)):
                continue
        dpdx, dpdy, bb = L.compute_psf_derivatives(psf, x, y, eb)
        if dpdx is None:
            continue
        ys, xs = L.bbox_to_slices(bb, eb)
        yy, xx = np.mgrid[bb.getMinY():bb.getMaxY() + 1, bb.getMinX():bb.getMaxX() + 1]
        g = good[ys, xs] & (((xx - x) ** 2 + (yy - y) ** 2) < R2)
        if g.sum() < 10:
            continue
        Pm, _ = L.psf_stamp(psf, x, y, bb, taper=L.WING_TAPER)
        if Pm is None or Pm.shape != dpdx.shape:
            continue
        w = 1.0 / np.sqrt(var[ys, xs][g])
        J = np.column_stack([Pm[g], a * dpdx[g], a * dpdy[g]]) * w[:, None]
        try:
            delta, *_ = np.linalg.lstsq(J, resid[ys, xs][g] * w, rcond=None)
        except np.linalg.LinAlgError:
            continue
        new[j] = (x + np.clip(delta[1], -max_shift, max_shift),
                  y + np.clip(delta[2], -max_shift, max_shift))
    return new


def _prune_merge(fit, S, min_snr, min_sep, bright, near_sep, near_snr, protect,
                 in_zone=None, zone_snr=None):
    """prune_and_merge that keeps the bookkeeping: returns (alive mask, new positions).
    Protected sources are never pruned; a non-protected source within min_sep of a
    protected one is absorbed by it (the protected position is kept)."""
    var, amps, good = fit["variance"], fit["amps"], fit["good"]
    xy = S["xy"].copy()
    n = len(xy)
    snr = np.zeros(n)
    for j, c in enumerate(fit["comps"]):
        if c is None:
            continue
        g = good[c["ys"], c["xs"]]
        snr[j] = amps[j] * np.sqrt(np.sum(c["stamp"][g] ** 2 / var[c["ys"], c["xs"]][g]))
    thresh = np.full(n, float(min_snr))
    if near_sep and bright is not None and np.any(bright):
        dd, _ = cKDTree(xy[bright]).query(xy)
        near = (dd > 0.5) & (dd < near_sep)
        thresh[near] = np.maximum(thresh[near], near_snr)
    n_zone_low = 0
    if in_zone is not None and zone_snr and in_zone.any():
        thresh[in_zone] = np.maximum(thresh[in_zone], zone_snr)
        n_zone_low = int((in_zone & (snr < thresh) & ~S["seed"]).sum())
    prot = S["seed"] if protect else np.zeros(n, bool)
    alive = (snr >= thresh) | prot
    n_low = int((~alive).sum())
    tree = cKDTree(S["xy"]) if n else None
    for a in np.lexsort((-amps, ~prot)):          # protected first, then brightest first
        if not alive[a]:
            continue
        for b in tree.query_ball_point(S["xy"][a], min_sep):
            if b == a or not alive[b] or prot[b]:
                continue
            if not prot[a] and amps[a] + amps[b] > 0:
                xy[a] = (amps[a] * xy[a] + amps[b] * xy[b]) / (amps[a] + amps[b])
            alive[b] = False
    print(f"  prune/merge: {n} -> {int(alive.sum())} ({n_low} low S/N"
          + (f", of which {n_zone_low} in saturated-star zones" if n_zone_low else "")
          + f", {n - n_low - int(alive.sum())} merged)")
    return alive, xy


def _eliminate_weak(exposure, psf, S, fit, bad_bitmask, P, fitkw, fwhm=None):
    n0 = len(S["xy"])
    rounds = 0
    scale = L.empirical_noise_scale(fit, exposure)
    min_snr = P["min_snr"]
    for _ in range(40):
        n = len(S["xy"])
        if n <= 1:
            break
        snr = amp_snr_full_fast(fit, scale)
        cand = ~S["seed"] if P["protect_seeds"] else np.ones(n, bool)
        weak = np.flatnonzero(cand & (snr < min_snr))
        if len(weak) == 0:
            break
        rounds += 1
        hopeless = snr[weak] < P["hopeless_frac"] * min_snr
        isolated = np.ones(len(weak), bool)          # no other weak source nearby: removing
        if fwhm and len(weak) > 1:                   # it barely changes anyone else's S/N
            cnt = np.array([len(l) for l in cKDTree(S["xy"][weak]).query_ball_point(
                S["xy"][weak], 2.0 * fwhm)])
            isolated = cnt <= 1
        rest = weak[~(hopeless | isolated)]
        rest = rest[np.argsort(snr[rest])][:int(np.ceil(P["elim_drop_frac"] * len(rest)))]
        drop = np.concatenate([weak[hopeless | isolated], rest])[:n - 1]
        keep = np.ones(n, bool)
        keep[drop] = False
        S = _take(S, keep)
        fit = _fit(exposure, psf, S, bad_bitmask, fit["model_full"], fitkw)
    print(f"  eliminate_weak: noise_scale={scale:.2f}, {n0} -> {len(S['xy'])} sources in "
          f"{rounds} round(s) (honest S/N>={min_snr:g}{', Gaia protected' if P['protect_seeds'] else ''})")
    return S, fit


def _shape_cut(exposure, psf, S, fit, bad_bitmask, P, fitkw, drop_frac=0.25, max_rounds=15,
               judge_uids=None, psf_sigma=None):
    """judge_uids: only these sources are judged (None = all non-Gaia point sources).
    Width failures that are too BROAD become extended sources if keep_extended."""
    kw = dict(P["shape_kwargs"] or {})
    mkw = {k: kw.pop(k) for k in ("bg_sigma", "width_scales", "min_pix") if k in kw}
    if not P["shape_snr_test"]:
        kw["min_snr_local"] = 1e-3
    n0, rejected, n_ext = len(S["xy"]), [], 0
    why = dict(snr=0, chi2=0, width=0, contrast=0)
    scale = L.empirical_noise_scale(fit, exposure)
    pscale = L.empirical_noise_scale(fit, exposure, matched=False)
    for _ in range(max_rounds):
        n = len(S["xy"])
        ok = (S["ext"] == 0) & (~S["seed"] if P["protect_seeds"] else True)
        if judge_uids is not None:
            ok &= np.isin(S["uid"], list(judge_uids))
        cand = np.flatnonzero(ok)
        if n <= 1 or len(cand) == 0:
            break
        sub = dict(fit, comps=[fit["comps"][j] for j in cand])
        m = L.shape_metrics(exposure, psf, [tuple(S["xy"][j]) for j in cand], sub,
                            noise_scale=scale, chi2_scale=pscale, **mkw)
        fail, score, flags = L.judge_shape(m, **kw)
        bad = np.flatnonzero(fail)
        if len(bad) == 0:
            break
        broad = (flags["width"] & ~flags["chi2"] & ~flags["contrast"] &
                 (np.nan_to_num(m["best_scale"]) > 1.0))
        if P["keep_extended"] and psf_sigma and broad[bad].any():
            conv = bad[broad[bad]]
            S["ext"][cand[conv]] = psf_sigma * np.sqrt(np.clip(m["best_scale"][conv] ** 2 - 1, 0.05, None))
            n_ext += len(conv)
            bad = bad[~broad[bad]]
            if len(bad) == 0:
                fit = _fit(exposure, psf, S, bad_bitmask, fit["model_full"], fitkw)
                continue
        order = bad[np.argsort(-score[bad])]
        sel = order[:min(max(1, int(np.ceil(drop_frac * len(bad)))), n - 1)]
        for i in sel:
            for k in why:
                why[k] += int(flags[k][i])
        drop = cand[sel]
        rejected += [(S["xy"][j][0], S["xy"][j][1], float(score[i])) for j, i in zip(drop, sel)]
        keep = np.ones(n, bool)
        keep[drop] = False
        S = _take(S, keep)
        fit = _fit(exposure, psf, S, bad_bitmask, fit["model_full"], fitkw)
    print(f"  shape cut: {n0} -> {len(S['xy'])} ({len(rejected)} removed, {n_ext} kept as extended; "
          f"removed for chi2 {why['chi2']}, width {why['width']}, contrast {why['contrast']}"
          + (f", S/N {why['snr']}" if P["shape_snr_test"] else "") + ")")
    return S, fit, rejected


def _satlike(exposure):
    """SAT | NEARSAT pixels (the latter if option B created that plane)."""
    m = exposure.mask
    out = (m.array & m.getPlaneBitMask("SAT")) != 0
    try:
        out |= (m.array & m.getPlaneBitMask("NEARSAT")) != 0
    except Exception:
        pass
    return out


def saturation_level(exposure, P):
    """(level, how) -- the saturation level in image units: near_sat_level if given, else the
    median of the SAT-flagged pixels, unless those look interpolated (below the brightest
    unflagged pixel), in which case the brightest unflagged pixel (a lower bound)."""
    if P.get("near_sat_level"):
        return float(P["near_sat_level"]), "given (near_sat_level)"
    data = exposure.image.array
    sat = (exposure.mask.array & exposure.mask.getPlaneBitMask("SAT")) != 0
    fin = np.isfinite(data)
    if (sat & fin).sum() < 5:
        return None, "no saturated pixels in this cutout -- set near_sat_level to use option B"
    lvl = float(np.nanmedian(data[sat & fin]))
    top = float(np.nanmax(data[~sat & fin]))
    if lvl >= 0.9 * top:
        return lvl, "median of the SAT-flagged pixels"
    return top, "brightest unflagged pixel (SAT pixel values look interpolated; a lower bound)"


def measure_model_error(exposure, res, P, verbose=True):
    """PSF model error eps (fraction of the model per pixel): around unsaturated point sources
    with peak S/N >= model_err_min_snr, the residual variance beyond the pixel variance is
    fitted as (eps * model)^2, from the robust (MAD) variance of residual/sigma in bins of
    model value (expected 1 + eps^2 * model^2 / var).
    Returns (eps, table) with rows (model, n_px, robust var of residual/sigma, model^2/var)."""
    fit = res["snaps"][-1]["fit"]
    S = res["S"]
    model = np.clip(res["model_full"], 0, None)
    data = exposure.image.array
    var = exposure.variance.array + P["poisson_k"] * model
    bb = exposure.getBBox()
    c = bb.getCenter()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    psf = exposure.getPsf()
    kern = psf.computeKernelImage(c).array
    fwhm = L._fwhm(psf, c)
    noise = float(np.sqrt(np.nanmedian(exposure.variance.array)))
    amps = fit["amps"]
    snr = amps * kern.max() / noise
    sel = np.flatnonzero((snr >= P["model_err_min_snr"]) & (S["ext"] == 0))
    satl = _satlike(exposure)
    if satl.any() and len(sel):
        from scipy.ndimage import binary_dilation
        near = binary_dilation(satl, iterations=int(np.ceil(3 * fwhm)))
        ix = np.clip(np.round(S["xy"][sel, 0]).astype(int) - X0, 0, satl.shape[1] - 1)
        iy = np.clip(np.round(S["xy"][sel, 1]).astype(int) - Y0, 0, satl.shape[0] - 1)
        sel = sel[~near[iy, ix]]
    if len(sel) < 3:
        return 0.0, []
    use = _star_mask(data.shape, X0, Y0, S["xy"][sel], 2.5 * fwhm) & fit["good"] & ~satl
    use &= np.isfinite(data) & (var > 0)
    m, r, v = model[use], (data - res["model_full"])[use], var[use]
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
        print(f"  PSF model error: {msg} ({len(sel)} stars, S/N >= {P['model_err_min_snr']:g}); "
              f"fit weights use sys_frac = {P['sys_frac']:.3g}")
    return eps, rows


def sat_zones(exposure, fwhm, P):
    """Saturated-star zones.  Returns (centres, mask): centres = (N, 3) array of
    (x, y, halo radius) per saturated blob, mask = bool image of the zone band, which
    follows the SAT mask's shape (pixels within w of the blob,
    w = sat_zone_fwhm*FWHM + sat_zone_scale*inscribed radius of the blob)."""
    from scipy.ndimage import label, center_of_mass, distance_transform_edt
    shape = exposure.image.array.shape
    empty = (np.zeros((0, 3)), np.zeros(shape, bool))
    if not P["sat_zone_fwhm"]:
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
        if area[l] < P["sat_zone_min_px"]:
            continue
        blob = lab == l
        r_in = float(inner[blob].max())
        w = P["sat_zone_fwhm"] * fwhm + P["sat_zone_scale"] * r_in
        mask |= distance_transform_edt(~blob) <= w
        cy, cx = center_of_mass(blob)
        cen.append((cx + bb.getMinX(), cy + bb.getMinY(), 2.0 * r_in + w))
        print(f"  saturated region at ({cen[-1][0]:7.1f}, {cen[-1][1]:7.1f}): {int(area[l])} px, "
              f"half-width {r_in:4.1f}px -> zone band {w:4.1f}px around it")
    if cen:
        print(f"  saturated-star zones cover {100 * mask.mean():.1f}% of the cutout; non-Gaia "
              f"sources there need S/N >= {P['sat_zone_snr']:g}")
    return np.array(cen).reshape(-1, 3), mask


def sat_core_labels(exposure, fwhm, min_px=4):
    """Label image of saturated cores, each grown by max(3 px, FWHM/2); 0 = none."""
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
        if area[l] < min_px:
            continue
        k += 1
        lab[binary_dilation(l0 == l, iterations=grow) & (lab == 0)] = k
    return lab


def _label_at(xy, lab, X0, Y0):
    xy = np.asarray(xy, float).reshape(-1, 2)
    if not len(xy) or not lab.any():
        return np.zeros(len(xy), np.int32)
    ix = np.round(xy[:, 0]).astype(int) - X0
    iy = np.round(xy[:, 1]).astype(int) - Y0
    ok = (ix >= 0) & (ix < lab.shape[1]) & (iy >= 0) & (iy < lab.shape[0])
    out = np.zeros(len(xy), np.int32)
    out[ok] = lab[iy[ok], ix[ok]]
    return out


def _one_per_core(S, amps, lab, X0, Y0):
    """Boolean keep-mask: one source per saturated core.  If Gaia has a star there, keep
    the BRIGHTEST Gaia star (G mag; its position is good to milliarcseconds, and the faint
    Gaia entries around bright stars are usually spurious).  Otherwise keep the source
    closest to the core's centre -- never pick by fitted flux, which is degenerate when
    the centre is masked."""
    keep = np.ones(len(S["xy"]), bool)
    core = _label_at(S["xy"], lab, X0, Y0)
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


def _refine_saturated(exposure, psf, S, fit, core_lab, X0, Y0, fwhm, n_iter=6, max_step=2.0):
    """Gauss-Newton position fit of saturated stars on their UNSATURATED wings (annulus out
    to ~3 FWHM beyond the core), holding all other sources fixed.  Skips stars with a
    Gaia position.  Returns new positions for all sources and the number moved."""
    eb = exposure.getBBox()
    data, var, good = exposure.image.array, fit["variance"], fit["good"]
    amps = fit["amps"]
    core = _label_at(S["xy"], core_lab, X0, Y0)
    xy = S["xy"].copy()
    moved = 0
    for j in np.flatnonzero(core > 0):
        if S["seed"][j] and np.isfinite(S["g"][j]):
            continue
        a = amps[j]
        if a <= 0 or fit["comps"][j] is None:
            continue
        c = fit["comps"][j]
        own = np.zeros_like(data)
        own[c["ys"], c["xs"]] = a * c["stamp"]
        others = fit["model_full"] - own
        rcore = np.sqrt((core_lab == core[j]).sum() / np.pi)
        x, y = xy[j]
        for _ in range(n_iter):
            dpdx, dpdy, bb = L.compute_psf_derivatives(psf, x, y, eb)
            if dpdx is None:
                break
            ys, xs = L.bbox_to_slices(bb, eb)
            Pm, _ = L.psf_stamp(psf, x, y, bb, taper=L.WING_TAPER)
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
            step = np.clip(d[1:], -max_step, max_step)
            x, y = x + step[0], y + step[1]
            if np.hypot(*step) < 0.02:
                break
        if np.hypot(x - xy[j, 0], y - xy[j, 1]) > 0.05:
            moved += 1
        xy[j] = (x, y)
    return xy, moved


def _zone_flags(xy, zone_mask, X0, Y0):
    xy = np.asarray(xy, float).reshape(-1, 2)
    if not len(xy) or not zone_mask.any():
        return np.zeros(len(xy), bool)
    ix = np.clip(np.round(xy[:, 0]).astype(int) - X0, 0, zone_mask.shape[1] - 1)
    iy = np.clip(np.round(xy[:, 1]).astype(int) - Y0, 0, zone_mask.shape[0] - 1)
    return zone_mask[iy, ix]


def bright_star_corrections(exposure, S, fit, model, psf_peak, noise, fwhm, bad_bitmask,
                            sigma_min=100.0, radius_fwhm=4.0, sweeps=2, radial=True,
                            zones=None, halo_modes=0, bright_modes=0, core_rings=(),
                            core_sigma=0.6, core_modes=0):
    """Local linear correction for every bright (or saturated) star, on top of the joint
    model, fitted on good pixels within radius_fwhm*FWHM:
        residual ~ amp * (a*P + bx*dP/dx + by*dP/dy + cxx*d2P/dx2 + cyy*d2P/dy2 + cxy*d2P/dxdy)
    i.e. flux, position (dipoles) and second-moment shape -- size AND ellipticity -- of
    each bright star may differ slightly from the shared PSF.  Weighted by the PIXEL
    variance (not the sys_frac-inflated one) so the core, where the mismatch is, counts.
    Returns (correction image, dict j -> dict(x, y, snr, n_sat, flux, dx, dy, chi2_before,
    chi2_after)).  Two sweeps, brightest first, for overlapping bright stars."""
    data, good = exposure.image.array, fit["good"]
    var = exposure.variance.array
    bb = exposure.getBBox()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    amps = fit["amps"]
    sat = (exposure.mask.array & exposure.mask.getPlaneBitMask("SAT")) != 0
    corr = np.zeros_like(model)
    snr = amps * psf_peak / noise
    picks = set(np.flatnonzero(snr > sigma_min).tolist())
    if sat.any():
        for j, (x, y) in enumerate(S["xy"]):
            ix, iy = int(round(x)) - X0, int(round(y)) - Y0
            if 0 <= iy < sat.shape[0] and 0 <= ix < sat.shape[1] and sat[iy, ix] and amps[j] > 0:
                picks.add(j)
    order = sorted(picks, key=lambda j: -amps[j])
    R = radius_fwhm * fwhm
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
            Rj, rks, is_sat = R, (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5), False
            if zones is not None and len(zones):           # saturated star: model its halo
                dz = np.hypot(zones[:, 0] - S["xy"][j, 0], zones[:, 1] - S["xy"][j, 1])
                k = int(np.argmin(dz))
                if dz[k] < 0.5 * zones[k, 2]:
                    Rj = max(R, zones[k, 2])
                    rks = tuple(np.arange(0.5, Rj / fwhm + 0.01, 0.75))
                    is_sat = True
            okpix = good[ys, xs] & (var[ys, xs] > 0)
            terms = [Pm, gx, gy, gxx, gyy, gyx]
            if radial:                  # smooth rings: core/halo profile differs from the PSF
                for rk in rks:
                    ring = np.exp(-0.5 * ((rr - rk * fwhm) / (0.5 * fwhm)) ** 2)
                    if (okpix & (ring > 0.5)).sum() >= 10:     # skip rings inside the mask
                        terms.append(Pm.max() * ring)
                        nm = (halo_modes if is_sat else (bright_modes if rk <= 4.0 else 0))
                        if nm:                                 # non-circular structure
                            th = np.arctan2(yy - py, xx - px)
                            for mm in range(1, int(nm) + 1):
                                terms += [Pm.max() * ring * np.cos(mm * th),
                                          Pm.max() * ring * np.sin(mm * th)]
            if core_rings and not is_sat:   # narrow core structure (1-2 px scale)
                thc = np.arctan2(yy - py, xx - px)
                for rk in core_rings:
                    ring = np.exp(-0.5 * ((rr - rk) / core_sigma) ** 2)
                    if (okpix & (ring > 0.5)).sum() >= 1:
                        terms.append(Pm.max() * ring)
                        base = ring * (rr / core_sigma if rk == 0 else 1.0)  # spike: its own dipole
                        for mm in range(1, int(core_modes) + 1):
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
                               n_sat=int((sat[ys, xs] & ring).sum()), flux=float(coef[0] * amps[j]),
                               dflux_pct=float(100 * coef[0]),
                               dx=float(-coef[1]), dy=float(-coef[2]),
                               chi2_before=float(np.mean(r ** 2)), chi2_after=float(np.mean(r2 ** 2)))
    return corr, info


def print_bright_corrections(info, max_rows=15):
    if not info:
        print("  no bright stars corrected")
        return
    print(f"  bright-star corrections ({len(info)} stars; chi2/px inside the fit radius with "
          f"the pixel variance, ~1 = residual at the noise level):")
    print(f"    {'x':>8s} {'y':>8s} {'peak S/N':>9s} {'SAT px':>6s} {'dx':>6s} {'dy':>6s} "
          f"{'dflux%':>7s} {'chi2 before':>11s} {'after':>7s}")
    for d in sorted(info.values(), key=lambda d: -d["snr"])[:max_rows]:
        print(f"    {d['x']:8.1f} {d['y']:8.1f} {d['snr']:9.0f} {d['n_sat']:6d} {d['dx']:+6.2f} "
              f"{d['dy']:+6.2f} {d['dflux_pct']:+7.2f} {d['chi2_before']:11.2f} {d['chi2_after']:7.2f}")


def _detect(working_residual, threshold, min_pixels):
    cfg = L.SourceDetectionTask.ConfigClass()
    cfg.thresholdValue = threshold
    cfg.thresholdType = "pixel_stdev"
    cfg.minPixels = min_pixels
    cfg.doTempLocalBackground = True
    cfg.tempLocalBackground.binSize = 16
    cfg.reEstimateBackground = False
    schema = L.afwTable.SourceTable.makeMinimalSchema()
    task = L.SourceDetectionTask(schema=schema, config=cfg)
    try:
        return task.run(L.afwTable.SourceTable.make(schema), working_residual).positive.getFootprints()
    except L.TooManyMaskedPixelsError:
        return None


def _slim(fit):
    return {k: v for k, v in fit.items() if k not in ("comps", "chi2map", "group_map")}


# =====================================================================
# Adaptive inner loop
# =====================================================================
def ring_mimics(exposure, S, fit, model_full, P, fwhm, psf_peak, noise):
    """Faint non-Gaia sources around bright stars that do NOT stand out from the azimuthal median
    of the residual at their distance from the star (fakes absorbing a smooth wing excess).
    Returns (indices to drop, n checked)."""
    data = exposure.image.array
    var, good, amps = fit["variance"], fit["good"], fit["amps"]
    xy = S["xy"]
    bb = exposure.getBBox()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    snr = amps * psf_peak / noise
    bright = np.flatnonzero((snr >= P["ring_check_snr"]) & (S["ext"] == 0))
    if not len(bright) or len(xy) < 2:
        return [], 0
    tree = cKDTree(xy)
    rmax = P["ring_check_fwhm"] * fwhm
    drop, n_checked = set(), 0
    for b in bright:
        cand = [k for k in tree.query_ball_point(xy[b], rmax)
                if k != b and not S["seed"][k] and S["ext"][k] == 0 and amps[k] > 0
                and amps[k] < P["ring_max_ratio"] * amps[b]
                and np.hypot(*(xy[k] - xy[b])) >= fwhm and fit["comps"][k] is not None]
        if not cand:
            continue
        x, y = xy[b]
        R = rmax + 2 * fwhm
        xa, xb = max(int(x - R) - X0, 0), min(int(x + R) + 2 - X0, data.shape[1])
        ya, yb = max(int(y - R) - Y0, 0), min(int(y + R) + 2 - Y0, data.shape[0])
        sl = (slice(ya, yb), slice(xa, xb))
        base = (data - model_full)[sl]
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
        M = int(P.get("ring_bg_modes", 0) or 0)
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
            n_checked += 1
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
            snr_k = f * np.sqrt(den)
            if P.get("ring_sym", False):
                from scipy.ndimage import map_coordinates
                dx, dy = xx[ap] - x, yy[ap] - y
                rot = []
                for c_, s_ in ((0, 1), (-1, 0), (0, -1)):        # 90, 180, 270 deg
                    xr, yr = x + dx * c_ - dy * s_, y + dx * s_ + dy * c_
                    rot.append(map_coordinates(res, [yr - (ya + Y0), xr - (xa + X0)], order=1,
                                               mode="nearest"))
                bsym = np.median(np.array(rot), axis=0)
                fs = np.sum(w[ap] * (rk[ap] - bsym) / v) / den
                snr_k = min(snr_k, fs * np.sqrt(den))
            if snr_k < P["ring_snr"]:
                drop.add(k)
    return sorted(drop), n_checked


def split_test(exposure, psf, S, fit, model_full, P, fwhm, psf_peak, noise, bad_bitmask, fitkw):
    """Unresolved pairs: for each source with split_min_snr <= peak S/N < bright_corr_sigma, look
    for a positive residual peak (matched-filter S/N >= split_peak_snr) 0.5-1.5 FWHM away and not
    near another source; add a trial source there.  All trials are fitted together; a trial is
    kept if the chi2 within 2 FWHM of it drops by >= split_dchi2 and its S/N >= min_snr.
    Returns (S, fit, n_tried, n_kept); fit is the input object if nothing was kept."""
    data = exposure.image.array
    var, good, amps = fit["variance"], fit["good"], fit["amps"]
    bb = exposure.getBBox()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    xy = S["xy"]
    resid = data - model_full
    usable = good & np.isfinite(resid) & (var > 0)
    zmap = ca._matched_snr(np.where(usable, resid, 0.0), np.where(usable, var, 1.0), usable, fwhm,
                           usable)[0]
    snr = amps * psf_peak / noise
    hosts = np.flatnonzero((snr >= P["split_min_snr"]) & (snr < P["bright_corr_sigma"])
                           & (S["ext"] == 0))
    tree = cKDTree(xy)
    trials, host_of = [], []
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
        if not np.isfinite(z[k]) or z[k] < P["split_peak_snr"]:
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
        host_of.append(j)
    if not trials:
        return S, fit, 0, 0
    S2 = _add(S, trials)
    f2 = _fit(exposure, psf, S2, bad_bitmask, fit["model_full"], fitkw)
    n0 = len(xy)
    keep = np.ones(len(S2["xy"]), bool)
    yy, xx = np.mgrid[Y0:Y0 + data.shape[0], X0:X0 + data.shape[1]]
    for t, (px, py) in enumerate(trials):
        jn = n0 + t
        x0_, x1_ = int(px - 2 * fwhm) - X0, int(px + 2 * fwhm) + 2 - X0
        y0_, y1_ = int(py - 2 * fwhm) - Y0, int(py + 2 * fwhm) + 2 - Y0
        sl = (slice(max(y0_, 0), max(y1_, 0)), slice(max(x0_, 0), max(x1_, 0)))
        m = good[sl] & (np.hypot(xx[sl] - px, yy[sl] - py) < 2 * fwhm) & (var[sl] > 0)
        c0 = np.sum((data[sl] - fit["model_full"][sl])[m] ** 2 / var[sl][m])
        c1 = np.sum((data[sl] - f2["model_full"][sl])[m] ** 2 / f2["variance"][sl][m])
        cj = f2["comps"][jn]
        if cj is None or f2["amps"][jn] <= 0:
            keep[jn] = False
            continue
        g = f2["good"][cj["ys"], cj["xs"]]
        s_new = f2["amps"][jn] * np.sqrt(np.sum(cj["stamp"][g] ** 2 / f2["variance"][cj["ys"], cj["xs"]][g]))
        if (c0 - c1) < P["split_dchi2"] or s_new < P["min_snr"]:
            keep[jn] = False
    n_kept = int(keep[n0:].sum())
    if not n_kept:
        return S, fit, len(trials), 0
    S3 = _take(S2, keep)
    f3 = _fit(exposure, psf, S3, bad_bitmask, fit["model_full"], fitkw)
    return S3, f3, len(trials), n_kept


def _ellip_model(psf, x, y, sa, q, th, shape, X0, Y0, prof="gauss"):
    """Image (full-array shape) of a unit-flux PSF convolved with an elliptical profile
    ('gauss' or 'exp'; major-axis scale sa [px], axis ratio q, angle th)."""
    from scipy.signal import fftconvolve
    arr, bb = _raw_psf_stamp(psf, x, y, None)
    arr = fftconvolve(arr, _ellip_kernel(sa, q, th, prof), mode="same")
    arr = arr / max(arr.sum(), 1e-30)
    out = np.zeros(shape)
    ya, xa = bb.getMinY() - Y0, bb.getMinX() - X0
    y0_, x0_ = max(ya, 0), max(xa, 0)
    y1_, x1_ = min(ya + arr.shape[0], shape[0]), min(xa + arr.shape[1], shape[1])
    if y1_ > y0_ and x1_ > x0_:
        out[y0_:y1_, x0_:x1_] = arr[y0_ - ya:y1_ - ya, x0_ - xa:x1_ - xa]
    return out


def ellip_test(exposure, psf, S, fit, model_full, P, fwhm, psf_peak, noise, bad_bitmask, fitkw,
               psf_sigma, verbose_ellip=True, only=None, explain=False):
    """Small galaxies.  A badly fitted source (local chi2/px >= ellip_chi2_min, S/N ellip_min_snr ..
    bright_corr_sigma) and the points within ellip_group_fwhm FWHM (at most one Gaia source) are
    refitted with alternative models:
        'one'      : ONE extended source replaces all of them
        'star+gal' : the Gaia source stays a point, ONE extended source replaces the others
    each with a 'gauss' or 'exp' (exponential disk) elliptical profile convolved with the PSF.
    The best model by chi2 + 5 per fitted number must beat the point sources by >= ellip_dchi2.
    Returns (S, fit, n_tried, n_kept)."""
    from scipy.optimize import least_squares
    data = exposure.image.array
    var, good, amps = fit["variance"], fit["good"], fit["amps"]
    bb = exposure.getBBox()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    xy = S["xy"]
    resid = data - model_full
    snr = amps * psf_peak / noise
    tree = cKDTree(xy)
    yy, xx = np.mgrid[Y0:Y0 + data.shape[0], X0:X0 + data.shape[1]]
    hosts = np.flatnonzero((snr >= P["ellip_min_snr"]) & (snr < P["bright_corr_sigma"])
                           & (S["ext"] == 0) & (amps > 0))
    if only is not None:
        for k in only:
            if explain and k not in hosts:
                print(f"    source {k} at ({xy[k, 0]:.1f}, {xy[k, 1]:.1f}) is not a candidate: peak S/N "
                      f"{snr[k]:.1f} (needs {P['ellip_min_snr']:g}..{P['bright_corr_sigma']:g}), "
                      f"extended={S['ext'][k] > 0}, flux {amps[k]:.3g}")
        hosts = np.array([k for k in only if k in hosts], int)
    used, accepted, n_try = set(), [], 0
    for k0 in hosts[np.argsort(-amps[hosts])]:
        j = int(k0)
        if j in used:
            continue
        x, y = xy[j]
        m2 = good & (np.hypot(xx - x, yy - y) < 2 * fwhm) & (var > 0)
        c2_loc = np.mean(resid[m2] ** 2 / var[m2]) if m2.sum() else np.nan
        if m2.sum() < 20 or c2_loc < P["ellip_chi2_min"]:
            if explain:
                print(f"    ({x:.1f}, {y:.1f}): local chi2/px {c2_loc:.2f} < {P['ellip_chi2_min']:g} "
                      f"(fit variance incl. sys_frac) -> not tried")
            continue
        near = [k for k in tree.query_ball_point(xy[j], P["ellip_group_fwhm"] * fwhm)
                if k != j and k not in used and S["ext"][k] == 0 and amps[k] > 0
                and fit["comps"][k] is not None]
        gaia_near = [k for k in near if S["seed"][k]]
        if (S["seed"][j] and gaia_near) or len(gaia_near) > 1:
            if explain:
                print(f"    ({x:.1f}, {y:.1f}): another Gaia source within {P['ellip_group_fwhm']:g} "
                      f"FWHM -> not merged (two Gaia stars)")
            continue
        gsrc = j if S["seed"][j] else (gaia_near[0] if gaia_near else None)
        grp = list(dict.fromkeys(([gsrc] if gsrc is not None else []) +
                                 [k for k in [j] + near if not S["seed"][k]]))
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
        chi_pt = float(np.sum((Tr - own[reg]) ** 2 * w))
        p_pt = 3 * len(grp)
        n_try += 1
        others = [k for k in grp if k != gsrc]
        variants = [("one", None, grp)]
        if gsrc is not None and others:
            variants.append(("star+gal", gsrc, others))
        best = None
        for vname, keep_pt, repl in variants:
            fr = amps[repl]
            gx = np.sum(fr * xy[repl, 0]) / fr.sum()
            gy = np.sum(fr * xy[repl, 1]) / fr.sum()
            fixed = stamps[keep_pt][reg] if keep_pt is not None else None
            for prof in ("gauss", "exp"):
                def comps_(p):
                    m = _ellip_model(psf, p[0], p[1], np.exp(p[2]), 0.1 + 0.9 / (1 + np.exp(-p[3])),
                                     p[4], data.shape, X0, Y0, prof)[reg]
                    A = [m] if fixed is None else [m, fixed]
                    A = np.array(A)
                    M = (A * w) @ A.T
                    try:
                        fl = np.linalg.solve(M + 1e-12 * np.eye(len(A)), (A * w) @ Tr)
                    except np.linalg.LinAlgError:
                        fl = np.zeros(len(A))
                    return A, fl

                def res_fn(p):
                    A, fl = comps_(p)
                    return (Tr - fl @ A) * np.sqrt(w)

                s0 = np.log(0.7 * psf_sigma if prof == "gauss" else 0.5 * psf_sigma)
                try:
                    sol = least_squares(res_fn, [gx, gy, s0, 1.0, 0.0],
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
                p_m = 6 + (3 if keep_pt is not None else 0)
                need = P["ellip_dchi2"] + 5.0 * max(0, p_m - p_pt)
                ok_ = (np.all(fl > 0) and chi_pt - chi >= need and sa >= 0.3 * psf_sigma
                       and q >= 0.15)
                if verbose_ellip:
                    print(f"    try {vname:8s} {prof:5s} at ({sol.x[0]:.1f}, {sol.x[1]:.1f}): "
                          f"{len(grp)} point(s) chi2 {chi_pt:.0f} -> {chi:.0f} (needs -{need:.0f}); "
                          f"scale {sa:.2f} px, q {q:.2f}, PA {np.degrees(th) % 180:.0f} deg"
                          f"{'' if np.all(fl > 0) else ', negative flux'} -> {'ok' if ok_ else 'no'}")
                score = chi - chi_pt + need
                if ok_ and (best is None or score < best[0]):
                    best = (score, vname, keep_pt, repl, sol.x[0], sol.x[1], sa, q, th, prof, chi)
        if best is not None:
            _, vname, keep_pt, repl, ex, ey, sa, q, th, prof, chi = best
            accepted.append((vname, keep_pt, repl, ex, ey, sa, q, th, prof, chi_pt, chi))
            used.update(grp)
    if not accepted:
        return S, fit, n_try, 0
    S = dict(S, xy=S["xy"].copy(), ext=S["ext"].copy(), ref=S["ref"].copy())
    keep = np.ones(len(S["xy"]), bool)
    for vname, keep_pt, repl, ex, ey, sa, q, th, prof, c0, c1 in accepted:
        host = repl[0] if keep_pt is not None else (next((k for k in repl if S["seed"][k]), repl[0]))
        S["xy"][host] = (ex, ey)
        if S["seed"][host]:            # a galaxy's Gaia centroid is unreliable: no drift limit
            S["ref"][host] = (ex, ey)
        S["ext"][host] = sa
        _ELLIP[int(S["uid"][host])] = (float(q), float(th), prof)
        for k in repl:
            if k != host:
                keep[k] = False
        print(f"    extended source ({vname}, {prof}) at ({ex:.1f}, {ey:.1f}) replaces {len(repl)} point "
              f"source(s){' next to the Gaia star' if keep_pt is not None else ''}: scale {sa:.2f} px, "
              f"axis ratio {q:.2f}, PA {np.degrees(th) % 180:.0f} deg; local chi2 {c0:.0f} -> {c1:.0f}")
    S = _take(S, keep)
    _clamp_seeds(S, P["seed_max_drift_px"])
    f2 = _fit(exposure, psf, S, bad_bitmask, fit["model_full"], fitkw)
    return S, f2, n_try, len(accepted)


def try_elliptical(run, xy, chi2_min=None, dchi2=None):
    """Debug / what-if on a FINISHED run: try the elliptical-source test on the source nearest
    to xy=(x, y) and print every step (why it is or is not tried, the fit, the chi2 change,
    any error).  Does not change the run.  chi2_min / dchi2 override the trigger / bar.
    Returns dict(kept, S, fit) -- S/fit with the elliptical source if it was kept."""
    if not _ORIG or _HUGE is None:
        enable_speedups()
    final = run.outer_history[-1]
    exp, res = final["exposure"], final["result"]
    S, fit = res["S"], res["snaps"][-1]["fit"]
    if fit.get("comps") is None:
        print("this run's last fit has no per-source stamps -- cannot test")
        return None
    P = dict(_params(run.config))
    P["poisson_k"] = float(final.get("poisson_k", 0.0) or 0.0)
    if chi2_min is not None:
        P["ellip_chi2_min"] = chi2_min
    if dchi2 is not None:
        P["ellip_dchi2"] = dchi2
    psf = exp.getPsf()
    c = exp.getBBox().getCenter()
    fwhm = L._fwhm(psf, c)
    psf_sigma = fwhm / 2.355
    psf_peak = psf.computeKernelImage(c).array.max()
    noise = float(np.sqrt(np.nanmedian(exp.variance.array)))
    j = int(cKDTree(S["xy"]).query(xy)[1])
    print(f"nearest source: index {j} at ({S['xy'][j, 0]:.1f}, {S['xy'][j, 1]:.1f}), "
          f"{'Gaia' if S['seed'][j] else 'non-Gaia'}, flux {fit['amps'][j]:.3g}, "
          f"peak S/N {fit['amps'][j] * psf_peak / noise:.1f}; FWHM {fwhm:.2f} px")
    S = {k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in S.items()}
    fitkw = dict(support_frac=P["support_frac"], sys_frac=P["sys_frac"], poisson_k=P["poisson_k"])
    n_before = len(_ELLIP)
    S2, f2, n_try, n_keep = ellip_test(exp, psf, S, fit, res["model_full"], P, fwhm, psf_peak, noise,
                                       run.bad_bitmask, fitkw, psf_sigma, only=[j], explain=True)
    if n_keep:                       # do not leave the what-if in the global table
        for u in list(_ELLIP)[n_before:]:
            _ELLIP.pop(u, None)
    print(f"tried {n_try}, kept {n_keep}")
    return dict(kept=bool(n_keep), S=S2, fit=f2)


def final_position_polish(exposure, psf, S, fit, bad_bitmask, P, fitkw, fwhm, psf_sigma):
    """Re-optimise ALL source positions jointly: Gauss-Newton steps (capped at final_polish_cap px,
    Gaia drift limit kept) + refit, repeated while the total chi2 drops.  Fixes dipoles left by
    positions that went stale when neighbours were added or removed.
    Returns (S, fit, info) -- fit is the input object if nothing improved."""
    data = exposure.image.array

    def chi_all(f):
        g = f["good"]
        return float(np.sum((data - f["model_full"])[g] ** 2 / f["variance"][g]))

    S = dict(S)
    c0 = c_prev = chi_all(fit)
    best = fit
    allidx = np.arange(len(S["xy"]))
    moved, n_r = 0, 0
    for r in range(1, P["final_polish_rounds"] + 1):
        before = S["xy"].copy()
        S["xy"] = _refine_subset(exposure, psf, S, best, allidx, P["final_polish_cap"],
                                 P["centroid_min_amp"], P["freeze_fwhm"] * fwhm, psf_sigma)
        _clamp_seeds(S, P["seed_max_drift_px"])
        f_new = _fit(exposure, psf, S, bad_bitmask, best["model_full"], fitkw)
        c_new = chi_all(f_new)
        if c_new >= c_prev:
            S["xy"] = before
            break
        d = np.hypot(*(S["xy"] - before).T)
        moved = max(moved, int((d > 0.01).sum()))
        best, c_prev, n_r = f_new, c_new, r
        if d.max() < P["centroid_shift_tol"]:
            break
    return S, best, dict(rounds=n_r, moved=moved, chi0=c0, chi1=c_prev)


def run_inner_seeded(exposure, psf, bad_bitmask, S, P, outer_idx=0, threshold_start=None,
                     shape_now=False, timer=None):
    """One outer pass: warm fit of S (Gaia seeds or the previous pass), then
    detect -> add -> joint fit -> refine -> clean, repeated until the residual check says
    stop.  Returns a result dict compatible with run_inner_loop_verbose's (plus 'S',
    'adaptive', 'is_seed', 'gaia_index')."""
    timer = timer or L.TIMER
    bbox = exposure.getBBox()
    ny, nx = bbox.getHeight(), bbox.getWidth()
    X0, Y0 = bbox.getMinX(), bbox.getMinY()
    data = exposure.image.array
    bad_mask = (exposure.mask.array & bad_bitmask) != 0
    sh = psf.computeShape(bbox.getCenter())
    psf_sigma = math.sqrt(0.5 * (sh.getIxx() + sh.getIyy()))
    fwhm = 2.355 * psf_sigma
    psf_peak = psf.computeKernelImage(bbox.getCenter()).array.max()
    noise = float(np.sqrt(np.nanmedian(exposure.variance.array)))
    fitkw = dict(support_frac=P["support_frac"], sys_frac=P["sys_frac"], poisson_k=P["poisson_k"])
    target = float(P["target_snr"] if P["target_snr"] is not None else P["min_snr"])
    t_floor, thr_min, step = P["t_floor"], P["threshold_min"], P["threshold_step"]
    thr = float(threshold_start if threshold_start is not None else P["threshold_start"])
    edge_px = (P["edge_margin_fwhm"] or 0.0) * fwhm
    pad_px = P["seed_pad_fwhm"] * fwhm
    merge_sep = (P["merge_separation_sigma"] or 2.0) * psf_sigma
    sched = list(P["centroid_shift_schedule"] or [1.5, 0.8, 0.4, 0.2, 0.1])
    sched += [sched[-1]] * max(0, P["n_centroid_rounds"] - len(sched))
    ts = sorted({3.0, 3.5, 4.0, 4.5, 5.0, 6.0, target, t_floor})

    core_lab = (sat_core_labels(exposure, fwhm, P["sat_zone_min_px"]) if P["one_source_per_sat"]
                else np.zeros((ny, nx), np.int32))
    okpix = ~bad_mask & np.isfinite(data)

    def rstd(img):                  # residual std over usable pixels only
        return float(np.nanstd(img[okpix]))

    zones = P.get("_sat_zones", np.zeros((0, 3)))
    zone_mask = P.get("_sat_zone_mask")
    if zone_mask is None or zone_mask.shape != (ny, nx):
        zone_mask = np.zeros((ny, nx), bool)
    if P["mask_pad_fwhm"] and bad_mask.any():           # flux of sources centred in the mask
        from scipy.ndimage import distance_transform_edt
        zone_mask = zone_mask | (distance_transform_edt(~bad_mask) <= P["mask_pad_fwhm"] * fwhm)

    def check(fit, model):
        resid = data - model
        return ca.residual_excess(resid, fit["variance"], fit["good"] & ~zone_mask, S["xy"],
                                  fit["amps"],
                                  X0, Y0, fwhm, psf_peak, noise, thresholds=ts,
                                  edge_fwhm=P["edge_fwhm"], bright_sigma=P["bright_peak_sigma"],
                                  halo_fwhm=P["halo_fwhm"], min_sep_fwhm=P["min_sep_fwhm"])

    fit, model_full = None, np.zeros((ny, nx))
    hist = dict(resid_std=[rstd(data)], n_peaks=[0], n_zero=[0])
    snaps, first_seen, log = [], {}, []
    bcorr = {}
    stamp_pt = bbox.getCenter()

    def corrected(fit):
        """model incl. bright-star correction terms; also updates the stamp-size rule"""
        nonlocal bcorr
        set_stamp_sizes(S["xy"], fit["amps"], psf, stamp_pt, noise, fwhm, P["stamp_floor_sigma"])
        _set_bf_fluxes(S["xy"], fit["amps"])
        if not P["bright_corr"]:
            bcorr = {}
            return fit["model_full"]
        corr, bcorr = bright_star_corrections(exposure, S, fit, fit["model_full"], psf_peak,
                                              noise, fwhm, bad_bitmask, P["bright_corr_sigma"],
                                              P["bright_corr_radius_fwhm"],
                                              radial=P["bright_corr_radial"], zones=zones,
                                              halo_modes=P["sat_halo_modes"],
                                              bright_modes=P["bright_angular_modes"],
                                              core_rings=P["bright_core_rings_px"],
                                              core_sigma=P["bright_core_ring_sigma"],
                                              core_modes=P["bright_core_modes"])
        return fit["model_full"] + corr
    prev_amp_uid = {}
    diag, stop, no_prog = None, "max_inner", 0
    prev_N, prev_e = None, None
    rejected = np.asarray(P.get("_halo_rejects", np.zeros((0, 2))), float).reshape(-1, 2).copy()
    shape_judged_once = False

    if len(S["xy"]):
        with timer.step("warm-start fit", outer=outer_idx):
            fit = _fit(exposure, psf, S, bad_bitmask, None, fitkw)
            fit = _fit(exposure, psf, S, bad_bitmask, fit["model_full"], fitkw)
        model_full = corrected(fit)
        prev_amp_uid = dict(zip(S["uid"].tolist(), fit["amps"]))
        hist = dict(resid_std=[rstd(data - model_full)], n_peaks=[len(S["xy"])],
                    n_zero=[fit["n_at_zero"]])
        print(f"  warm start: {len(S['xy'])} sources ({int(S['seed'].sum())} Gaia), residual "
              f"std {hist['resid_std'][0]:.3f}")

    for it in range(P["max_inner"]):
        times = {}

        def T(label):
            return timer.step(label, outer=outer_idx, inner=it)

        print(f"\n--- outer {outer_idx + 1} | inner {it + 1} (cap {P['max_inner']}) | "
              f"threshold={thr:.2f}σ | sources={len(S['xy'])} ---")
        det_in = data - P["detect_gain"] * model_full
        wr = exposure.clone()
        wr.image.array[:] = det_in
        if fit is not None:
            wr.variance.array[:] = fit["variance"]
        wr.setPsf(psf)
        with T("detect on residual") as rec:
            fps = _detect(wr, thr, P["min_pixels"])
        times["detect"] = rec["seconds"]
        if fps is None:
            stop = "detection failed (background estimation)"
            break
        raw_peaks = [(pk.getFx(), pk.getFy()) for fp in fps for pk in fp.getPeaks()]

        # ---- filter the new candidates
        with T("filter new peaks") as rec:
            new, _ = L.collapse_saturated_peaks(raw_peaks, exposure)
            if len(S["xy"]) and len(new) and fit is not None:
                is_bright = fit["amps"] * psf_peak > P["bright_peak_sigma"] * noise
                new, n_rej = L.reject_near_bright(new, _xy(S), is_bright, P["reject_fwhm"] * fwhm)
            if P["coarse_centroid"] and new:
                hw = max(int(round(P["coarse_centroid_hw_sigma"] * psf_sigma)), 3)
                new = [L.windowed_centroid(det_in, (X0, Y0), x, y, hw) for x, y in new]
            new, _ = L.drop_masked_peaks(new, exposure)
            new = [(x, y) for x, y in new
                   if X0 + edge_px <= x <= bbox.getMaxX() - edge_px
                   and Y0 + edge_px <= y <= bbox.getMaxY() - edge_px]
            if len(new) > 1:
                new = L.merge_close_peaks(new, exposure, min_separation=merge_sep)
            if new and len(S["xy"]):
                dd, _ = cKDTree(S["xy"]).query(np.asarray(new))
                new = [p for p, di in zip(new, dd) if di >= merge_sep]
            n_rem = 0
            if new and core_lab.any() and len(S["xy"]):     # core already has its source
                taken = set(_label_at(S["xy"], core_lab, X0, Y0).tolist()) - {0}
                lc = _label_at(new, core_lab, X0, Y0)
                new = [p for p, l in zip(new, lc) if l not in taken]
            if new and P["remember_rejects"] and len(rejected):
                dd, _ = cKDTree(rejected).query(np.asarray(new))
                n_rem = int((dd < merge_sep).sum())
                new = [p for p, di in zip(new, dd) if di >= merge_sep]
        times["filter"] = rec["seconds"]
        print(f"  {len(raw_peaks)} raw peaks -> {len(new)} new candidate(s)"
              + (f" ({n_rem} skipped: already rejected this pass)" if n_rem else ""))

        n_added = 0
        if (new or snaps == []) and (len(S["xy"]) + len(new)) > 0:
            first_new_uid = S["next_uid"]
            S = _add(S, new)
            with T("joint fit") as rec:
                fit = _fit(exposure, psf, S, bad_bitmask, model_full if fit is not None else None,
                           fitkw)
            times["fit"] = rec["seconds"]

            if P["refine_centroids"]:
                first = not snaps
                n_rounds = P["n_centroid_rounds"] if first else P["centroid_rounds_after_first"]
                if first:
                    idx = np.arange(len(S["xy"]))
                else:                          # new sources + their neighbours only
                    isnew = S["uid"] >= first_new_uid
                    idx = np.arange(len(S["xy"]))
                    if isnew.any():
                        dd, _ = cKDTree(S["xy"][isnew]).query(S["xy"])
                        idx = np.flatnonzero(dd < 2.0 * fwhm)
                    else:
                        idx = idx[:0]
                r, shift = -1, 0.0
                with T("refine centroids") as rec:
                    for r in range(n_rounds if len(idx) else 0):
                        before = S["xy"].copy()
                        cap = sched[r] if first else P["later_shift_cap"]
                        S["xy"] = _refine_subset(exposure, psf, S, fit, idx, cap,
                                                 P["centroid_min_amp"], P["freeze_fwhm"] * fwhm,
                                                 psf_sigma)
                        _clamp_seeds(S, P["seed_max_drift_px"])
                        fit = _fit(exposure, psf, S, bad_bitmask, fit["model_full"], fitkw)
                        shift = float(np.max(np.hypot(*(S["xy"] - before).T))) if len(before) else 0.0
                        if shift < P["centroid_shift_tol"]:
                            break
                if P["sat_refine_wings"] and core_lab.any():
                    S["xy"], n_mv = _refine_saturated(exposure, psf, S, fit, core_lab, X0, Y0, fwhm)
                    if n_mv:
                        fit = _fit(exposure, psf, S, bad_bitmask, fit["model_full"], fitkw)
                        print(f"    saturated stars: {n_mv} position(s) refined from the wings")
                times["centroids"] = rec["seconds"]
                print(f"    centroids: {len(idx)} source(s), {r + 1} round(s), last max shift "
                      f"{shift:.3f}px")

            with T("prune + merge + bounds") as rec:
                bright = fit["amps"] * psf_peak > P["bright_peak_sigma"] * noise
                alive, S["xy"] = _prune_merge(fit, S, P["min_snr"], P["prune_sep_sigma"] * psf_sigma,
                                              bright, P["near_bright_fwhm"] * fwhm,
                                              P["near_bright_snr"], P["protect_seeds"],
                                              in_zone=_zone_flags(S["xy"], zone_mask, X0, Y0),
                                              zone_snr=P["sat_zone_snr"])
                x, y = S["xy"][:, 0], S["xy"][:, 1]
                m = np.where(S["seed"], pad_px, -edge_px)        # seeds may sit just outside
                inside = ((x >= X0 - m) & (x <= bbox.getMaxX() + m) &
                          (y >= Y0 - m) & (y <= bbox.getMaxY() + m))
                n_edge = int((alive & ~inside).sum())
                S = _take(S, alive & inside)
                if core_lab.any() and len(S["xy"]):
                    amps_now = fit["amps"][np.flatnonzero(alive & inside)]
                    one = _one_per_core(S, amps_now, core_lab, X0, Y0)
                    if not one.all():
                        print(f"  saturated cores: removed {int((~one).sum())} extra source(s) "
                              f"({int((~one & S['seed']).sum())} of them Gaia entries)")
                        S = _take(S, one)
                if n_edge:
                    print(f"  dropped {n_edge} source(s) outside / within the edge margin")
                fit = _fit(exposure, psf, S, bad_bitmask, fit["model_full"], fitkw) if len(S["xy"]) else None
            times["prune"] = rec["seconds"]
            if fit is None:
                stop = "everything pruned"
                break

            if P["use_eliminate_weak"]:
                with T("eliminate weak") as rec:
                    S, fit = _eliminate_weak(exposure, psf, S, fit, bad_bitmask, P, fitkw, fwhm)
                times["eliminate weak"] = rec["seconds"]
            shape_rejected = []
            if shape_now:
                with T("shape cut") as rec:
                    ju = None
                    if shape_judged_once:          # only new sources and their neighbours
                        isnew = S["uid"] >= first_new_uid
                        ju = set()
                        if isnew.any():
                            dd, _ = cKDTree(S["xy"][isnew]).query(S["xy"])
                            ju = set(S["uid"][dd < 2.0 * fwhm].tolist())
                    S, fit, shape_rejected = _shape_cut(exposure, psf, S, fit, bad_bitmask, P, fitkw,
                                                        judge_uids=ju, psf_sigma=psf_sigma)
                    shape_judged_once = True
                times["shape cut"] = rec["seconds"]

            model_full = corrected(fit)
            resid = data - model_full
            n_added = int((S["uid"] >= first_new_uid).sum())
            if new:
                kept = set(S["uid"].tolist())
                lost = [p for k, p in enumerate(new) if first_new_uid + k not in kept]
                if lost:
                    rejected = np.vstack([rejected, np.asarray(lost, float).reshape(-1, 2)])
            is_new = S["uid"] >= first_new_uid
            prev_amp = np.array([prev_amp_uid.get(u, np.nan) for u in S["uid"].tolist()])
            prev_amp_uid = dict(zip(S["uid"].tolist(), fit["amps"]))
            for x, y in S["xy"]:
                first_seen.setdefault((round(x, 1), round(y, 1)), it)
            old_std, new_std = hist["resid_std"][-1], rstd(resid)
            hist["resid_std"].append(new_std)
            hist["n_peaks"].append(len(S["xy"]))
            hist["n_zero"].append(fit["n_at_zero"])
            if snaps and not P["keep_all_fits"]:
                snaps[-1]["fit"] = _slim(snaps[-1]["fit"])
            snaps.append(dict(it=it, threshold=thr, det_in=det_in, raw_peaks=raw_peaks,
                              collapsed_peaks=new, merged_peaks=_xy(S), merged_arr=S["xy"].copy(),
                              is_new=is_new, prev_amp=prev_amp,
                              footprint_mask=L._footprints_to_mask(fps, exposure), fit=fit,
                              residual=resid.copy(), std_before=old_std, std_after=new_std,
                              times=times, bad_mask=bad_mask, shape_rejected=shape_rejected,
                              is_seed=S["seed"].copy()))
            with T("residual check"):
                diag = check(fit, model_full)
            print(f"  {len(S['xy'])} sources ({int(S['seed'].sum())} Gaia, "
                  f"{int((S['ext'] > 0).sum())} extended, {n_added} new kept), "
                  f"residual std {old_std:.3f} -> {new_std:.3f}; time "
                  + ", ".join(f"{k}={v:.1f}s" for k, v in times.items()))
        elif not len(S["xy"]):
            print("  no sources yet")

        # ---- adaptive decision ---------------------------------------------------
        N = len(S["xy"])
        if diag is not None:
            missing, e, err = ca.is_missing(diag, target, P["z_sig"], P["min_excess"])
            p_, n_ = ca.counts_at(diag, target)
        else:
            missing, e, err, p_, n_ = True, 0, 0.0, 0, 0
        growing = n_added >= max(P["grow_min"], P["grow_frac"] * N)
        progress = n_added > 0 and (prev_N is None or N > prev_N or (prev_e is not None and e < prev_e))
        at_floor = thr <= thr_min + 1e-9
        thr_old = thr
        print(f"  residual check at >={target:g}σ: +{p_}/-{n_} peaks, excess {e:+d} ± {err:.1f} "
              f"-> {'MISSING sources' if missing else 'consistent with noise'}")
        if growing and not at_floor:
            action = f"still finding sources ({n_added} added) -> hold threshold {thr:.2f}σ"
        elif not at_floor:
            thr = max(thr + step, thr_min)
            action = f"lower threshold {thr_old:.2f} -> {thr:.2f}σ"
        elif growing and no_prog < P["patience"]:
            action = f"at threshold floor, still adding ({n_added}) -> continue"
        elif missing and thr > t_floor + 1e-9 and diag is not None:
            t_new = ca.pick_threshold(diag, P["fdr"], t_floor, thr)
            if t_new < thr - 1e-9:
                thr = t_new
                action = f"excess remains -> go deeper: threshold {thr_old:.2f} -> {thr:.2f}σ"
            else:
                stop, action = "excess remains, cannot go deeper (false-positive limit)", "STOP"
        else:
            stop = "residual clean" if not missing else "stalled (no progress)"
            action = "STOP"
        no_prog = 0 if (progress or thr != thr_old) else no_prog + 1
        if action != "STOP" and at_floor and no_prog >= P["patience"] and thr == thr_old:
            stop, action = f"no progress for {no_prog} iterations", "STOP"
        log.append(dict(it=it, threshold=thr_old, n=N, added=n_added, pos=p_, neg=n_, excess=e,
                        missing=missing, action=action))
        print(f"  => {action}" + (f" ({stop})" if action == "STOP" else ""))
        prev_N, prev_e = N, e
        if action == "STOP":
            break
    else:
        print(f"  WARNING: hit max_inner={P['max_inner']} while still working -- raise max_inner")

    if P["drop_halo_mimics"] and fit is not None and snaps and len(S["xy"]):
        import contextlib
        import io
        import types
        fake = types.SimpleNamespace(outer_history=[dict(
            exposure=exposure, result=dict(S=S, snaps=[dict(fit=fit)], model_full=model_full))])
        with contextlib.redirect_stdout(io.StringIO()):
            tb = check_bright_neighbours(fake, bright_snr=P["halo_check_snr"],
                                         within_fwhm=P["halo_check_fwhm"], plot=False,
                                         verbose=False, similar_min=P["halo_similar_min"])
        bad_ids = ([int(r["id"]) for r in tb if r["verdict"] == "doubtful"] if len(tb) else [])
        if bad_ids:
            keep = np.ones(len(S["xy"]), bool)
            keep[bad_ids] = False
            P["_halo_rejects"] = np.vstack([np.asarray(P.get("_halo_rejects", np.zeros((0, 2)))
                                                       ).reshape(-1, 2), S["xy"][~keep]])
            n_tot = len(tb)
            S = _take(S, keep)
            fit = _fit(exposure, psf, S, bad_bitmask, fit["model_full"], fitkw)
            model_full = corrected(fit)
            diag = check(fit, model_full)
            last = snaps[-1]
            last.update(fit=fit, merged_arr=S["xy"].copy(), merged_peaks=_xy(S),
                        residual=(data - model_full).copy(), is_seed=S["seed"].copy(),
                        is_new=np.asarray(last["is_new"])[keep] if len(last["is_new"]) == len(keep)
                        else np.zeros(len(S["xy"]), bool),
                        prev_amp=np.asarray(last["prev_amp"])[keep] if len(last["prev_amp"]) == len(keep)
                        else np.full(len(S["xy"]), np.nan))
            print(f"  halo check: removed {len(bad_ids)} of {n_tot} non-Gaia source(s) next to bright "
                  f"stars as halo mimics (repeat at >= {P['halo_similar_min']} control angles or "
                  f"significance < 3); they will not be re-detected")
    if P["drop_ring_mimics"] and fit is not None and snaps and len(S["xy"]):
        bad_ids, n_chk = ring_mimics(exposure, S, fit, model_full, P, fwhm, psf_peak, noise)
        if bad_ids:
            keep = np.ones(len(S["xy"]), bool)
            keep[bad_ids] = False
            P["_halo_rejects"] = np.vstack([np.asarray(P.get("_halo_rejects", np.zeros((0, 2)))
                                                       ).reshape(-1, 2), S["xy"][~keep]])
            S = _take(S, keep)
            fit = _fit(exposure, psf, S, bad_bitmask, fit["model_full"], fitkw)
            model_full = corrected(fit)
            diag = check(fit, model_full)
            last = snaps[-1]
            last.update(fit=fit, merged_arr=S["xy"].copy(), merged_peaks=_xy(S),
                        residual=(data - model_full).copy(), is_seed=S["seed"].copy(),
                        is_new=np.asarray(last["is_new"])[keep] if len(last["is_new"]) == len(keep)
                        else np.zeros(len(S["xy"]), bool),
                        prev_amp=np.asarray(last["prev_amp"])[keep] if len(last["prev_amp"]) == len(keep)
                        else np.full(len(S["xy"]), np.nan))
        print(f"  ring check: removed {len(bad_ids)} of {n_chk} faint non-Gaia source(s) around bright "
              f"stars that do not stand out from the star's wing at their radius (S/N < "
              f"{P['ring_snr']:g}); they will not be re-detected")
    def _polish(label):
        nonlocal S, fit, model_full, diag
        S, f_new, info = final_position_polish(exposure, psf, S, fit, bad_bitmask, P, fitkw, fwhm,
                                               psf_sigma)
        if f_new is not fit:
            fit = f_new
            model_full = corrected(fit)
            diag = check(fit, model_full)
            snaps[-1].update(fit=fit, merged_arr=S["xy"].copy(), merged_peaks=_xy(S),
                             residual=(data - model_full).copy(), is_seed=S["seed"].copy())
        print(f"  {label}: {info['rounds']} round(s), up to {info['moved']} source(s) "
              f"moved, total chi2 {info['chi0']:.4g} -> {info['chi1']:.4g} "
              f"({100 * (info['chi0'] - info['chi1']) / max(info['chi0'], 1e-12):+.2f}%)")

    if P["final_polish"] and fit is not None and snaps and len(S["xy"]):
        _polish("final position polish")
    if P["split_test"] and fit is not None and snaps and len(S["xy"]):
        S_new, f_new, n_try, n_keep = split_test(exposure, psf, S, fit, model_full, P, fwhm,
                                                 psf_peak, noise, bad_bitmask, fitkw)
        print(f"  split test: {n_try} source(s) with a lopsided residual tried as pairs, "
              f"{n_keep} second source(s) kept")
        if n_keep:
            S, fit = S_new, f_new
            model_full = corrected(fit)
            diag = check(fit, model_full)
            snaps[-1].update(fit=fit, merged_arr=S["xy"].copy(), merged_peaks=_xy(S),
                             residual=(data - model_full).copy(), is_seed=S["seed"].copy(),
                             is_new=np.zeros(len(S["xy"]), bool),
                             prev_amp=np.full(len(S["xy"]), np.nan))
            if P["final_polish"]:
                _polish("position polish after split")
    if P["ellip_test"] and fit is not None and snaps and len(S["xy"]):
        S_new, f_new, n_try, n_keep = ellip_test(exposure, psf, S, fit, model_full, P, fwhm,
                                                 psf_peak, noise, bad_bitmask, fitkw, psf_sigma)
        print(f"  elliptical test: {n_try} badly fitted source(s) tried as elliptical extended "
              f"sources, {n_keep} kept")
        if n_keep:
            S, fit = S_new, f_new
            model_full = corrected(fit)
            diag = check(fit, model_full)
            snaps[-1].update(fit=fit, merged_arr=S["xy"].copy(), merged_peaks=_xy(S),
                             residual=(data - model_full).copy(), is_seed=S["seed"].copy(),
                             is_new=np.zeros(len(S["xy"]), bool),
                             prev_amp=np.full(len(S["xy"]), np.nan))
    final_missing = bool(ca.is_missing(diag, target, P["z_sig"], P["min_excess"])[0]) if diag else True
    amps_out = np.array([]) if fit is None else fit["amps"].copy()
    for j, d in bcorr.items():
        if j < len(amps_out):
            amps_out[j] += d["flux"]
    if P["bright_corr"]:
        print_bright_corrections(bcorr)
    return dict(model_full=model_full, residual_full=data - model_full, peaks=_xy(S),
                amps=amps_out, bright_corr=bcorr,
                resid_std_history=hist["resid_std"], n_peaks_history=hist["n_peaks"],
                n_zero_amp_history=hist["n_zero"], first_seen=first_seen, snaps=snaps,
                S=S, is_seed=S["seed"].copy(), gaia_index=S["gaia"].copy(),
                adaptive=dict(final=diag, target=target, missing=final_missing, stop=stop,
                              log=log, n_inner=len(log)))


# =====================================================================
# Adaptive outer loop
# =====================================================================
def select_psf_stars(exposure, res, fwhm, psf_peak, noise, P, verbose=True):
    """Indices of clean, bright, isolated, unsaturated point sources for the PSF stack."""
    S = res["S"]
    fit = res["snaps"][-1]["fit"]
    amps = fit["amps"]
    xy = S["xy"]
    bb = exposure.getBBox()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    snr = amps * psf_peak / noise
    ok = (snr >= P["psf_stack_min_snr"]) & (S["ext"] == 0)
    sat = _satlike(exposure)
    if sat.any():
        from scipy.ndimage import binary_dilation
        near_sat = binary_dilation(sat, iterations=int(np.ceil(3 * fwhm)))
        ix = np.clip(np.round(xy[:, 0]).astype(int) - X0, 0, sat.shape[1] - 1)
        iy = np.clip(np.round(xy[:, 1]).astype(int) - Y0, 0, sat.shape[0] - 1)
        ok &= ~near_sat[iy, ix]
    h = 0.5 * fwhm * 4                                   # keep the stamp core inside the image
    ok &= ((xy[:, 0] > X0 + h) & (xy[:, 0] < bb.getMaxX() - h) &
           (xy[:, 1] > Y0 + h) & (xy[:, 1] < bb.getMaxY() - h))
    if len(xy) > 1:
        tree = cKDTree(xy)
        for j in np.flatnonzero(ok):
            nb = [k for k in tree.query_ball_point(xy[j], P["psf_stack_iso_fwhm"] * fwhm) if k != j]
            if any(amps[k] > 0.1 * amps[j] for k in nb):
                ok[j] = False
    rej = np.asarray(P.get("_halo_rejects", np.zeros((0, 2))), float).reshape(-1, 2)
    if len(rej):                    # their light is back in the residual -> would bias the stack
        dd, _ = cKDTree(rej).query(xy)
        ok &= dd > 3.0 * fwhm
    idx = np.flatnonzero(ok)
    idx = idx[np.argsort(-amps[idx])][:P["psf_stack_n"]]
    if verbose:
        print(f"  PSF stack: {len(idx)} clean star(s) (of {int((snr >= P['psf_stack_min_snr']).sum())} "
              f"with peak S/N >= {P['psf_stack_min_snr']:g}; saturated / crowded / extended excluded)")
    return idx


def _star_mask(shape, X0, Y0, pos, r):
    yy, xx = np.mgrid[Y0:Y0 + shape[0], X0:X0 + shape[1]]
    m = np.zeros(shape, bool)
    for x, y in pos:
        m |= (xx - x) ** 2 + (yy - y) ** 2 < r * r
    return m


def _chi2_on(exposure, f, m):
    m = m & f["good"]
    return float(np.mean((exposure.image.array - f["model_full"])[m] ** 2
                         / exposure.variance.array[m])) if m.any() else np.nan


def refine_psf_clean(exposure, psf, res, bad_bitmask, P, fitkw, n_cycles=3, verbose=True):
    """Residual-stack PSF update from clean stars, CROSS-VALIDATED.

    The clean stars are split into two halves A and B (alternating in brightness).  Each
    cycle builds a correction from A and judges it on B's pixels, and vice versa.  Only if
    the held-out chi2 improves (>= psf_cv_min_gain) is the correction from ALL stars applied.
    This stops the PSF absorbing the pixel noise of the stars it was built from (which
    always looks like an improvement on those same stars and shows up as speckle around
    every bright star).  Corrections are also smoothed (psf_smooth_px) and soft-thresholded
    at psf_soft_sigma."""
    from scipy.ndimage import gaussian_filter
    bb = exposure.getBBox()
    c = bb.getCenter()
    fwhm = L._fwhm(psf, c)
    psf_peak = psf.computeKernelImage(c).array.max()
    noise = float(np.sqrt(np.nanmedian(exposure.variance.array)))
    idx = select_psf_stars(exposure, res, fwhm, psf_peak, noise, P, verbose)
    if len(idx) < 6:
        print(f"  PSF stack: {len(idx)} clean star(s) -- too few to cross-validate, PSF not updated")
        return psf
    S = res["S"]
    fit = res["snaps"][-1]["fit"]
    data = exposure.image.array
    X0, Y0 = bb.getMinX(), bb.getMinY()
    A, B = idx[0::2], idx[1::2]
    mA = _star_mask(data.shape, X0, Y0, S["xy"][A], 2.5 * fwhm)
    mB = _star_mask(data.shape, X0, Y0, S["xy"][B], 2.5 * fwhm)
    skw = dict(P["stack_kwargs"] or {})
    skw.setdefault("soft_sigma", P["psf_soft_sigma"])
    sc = 0.5 * max(bb.getWidth(), bb.getHeight())

    def corrected_kernel(kernel, sel, f):
        new, _ = L.estimate_psf_correction(data, f["model_full"], f["variance"], f["good"],
                                           (X0, Y0), S["xy"][sel], f["amps"][sel], kernel,
                                           **dict(skw, n_use=len(sel)))
        d = gaussian_filter(new - kernel, P["psf_smooth_px"]) if P["psf_smooth_px"] else new - kernel
        out = kernel + d
        return out / out.sum()

    def as_psf(k, like):
        prev = _SPATIAL.get(id(like))
        if prev is not None and prev[0] is like and prev[1].shape == k.shape:
            return spatial_kernel_psf(k, prev[1], prev[2], c.getX(), c.getY(), sc)
        return L.kernel_to_lsst_psf(k)

    cur = psf
    for cyc in range(max(1, n_cycles)):
        kernel = cur.computeKernelImage(c).array.copy()
        try:
            fA = _fit(exposure, as_psf(corrected_kernel(kernel, A, fit), cur), S, bad_bitmask,
                      fit["model_full"], fitkw)
            fB = _fit(exposure, as_psf(corrected_kernel(kernel, B, fit), cur), S, bad_bitmask,
                      fit["model_full"], fitkw)
        except ValueError as e:
            print(f"    PSF cycle {cyc + 1}: {e}")
            break
        before = _chi2_on(exposure, fit, mB) + _chi2_on(exposure, fit, mA)
        after = _chi2_on(exposure, fA, mB) + _chi2_on(exposure, fB, mA)
        gain = (before - after) / max(before, 1e-12)
        keep = gain >= P["psf_cv_min_gain"]
        if verbose:
            print(f"    PSF cycle {cyc + 1}/{n_cycles}: held-out chi2/px on {len(idx)} PSF stars "
                  f"{before / 2:.2f} -> {after / 2:.2f} ({100 * gain:+.1f}%) "
                  f"({'kept' if keep else 'rejected, stopping'})")
        if not keep:
            break
        cand = as_psf(corrected_kernel(kernel, idx, fit), cur)
        fit = _fit(exposure, cand, S, bad_bitmask, fit["model_full"], fitkw)
        cur = cand
    return cur


def _residual_stack(resid, good, xy, amps, idx, X0, Y0, half, with_err=False):
    """Median over stars idx of resid / flux, sampled at each star's sub-pixel position on a
    (2*half+1)^2 grid centred on it (same grid as the PSF kernel image).  NaN where < 3 stars."""
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
    if not with_err:
        return med
    # error of the median: 1.25 * robust sigma / sqrt(n)
    err = 1.2533 * 1.4826 * np.nanmedian(np.abs(st - med), axis=0) / np.sqrt(np.maximum(nn, 1))
    return med, err


def refine_psf_median(exposure, psf, res, bad_bitmask, P, fitkw, fit0=None, verbose=True):
    """Correct the PSF CORE with the median residual/flux of the bright clean stars.

    The median is robust to neighbours, cosmic rays and per-star noise, so what it keeps is
    the shape error SHARED by the bright stars (the pattern show_bright_star_stack shows) --
    which the mean-based clean-star stack (soft threshold, fainter stars) can miss.
    Applied within psf_core_radius_fwhm*FWHM (cosine taper), smoothed by psf_smooth_px, and
    cross-validated: built from half the stars, judged on the other half's pixels."""
    from scipy.ndimage import gaussian_filter
    bb = exposure.getBBox()
    c = bb.getCenter()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    kernel = psf.computeKernelImage(c).array.astype(float).copy()
    fwhm = L._fwhm(psf, c)
    noise = float(np.sqrt(np.nanmedian(exposure.variance.array)))
    S = res["S"]
    fit = fit0 if fit0 is not None else _fit(exposure, psf, S, bad_bitmask,
                                             res["snaps"][-1]["fit"]["model_full"], fitkw)
    Pq = dict(P, psf_stack_min_snr=P["psf_core_min_snr"], psf_stack_n=60)
    idx = select_psf_stars(exposure, dict(res, snaps=[dict(fit=fit)]), fwhm, kernel.max(),
                           noise, Pq, verbose=False)
    if len(idx) < 6:
        if verbose:
            print(f"  PSF core: {len(idx)} bright clean star(s) (S/N >= {P['psf_core_min_snr']:g}) "
                  f"-- too few, not corrected")
        return psf
    rmax = P["psf_core_radius_fwhm"] * fwhm
    half = int(np.ceil(rmax)) + 1
    H = kernel.shape[0] // 2
    half = min(half, H)
    u = np.arange(-half, half + 1)
    R = np.hypot(*np.meshgrid(u, u))
    taper = np.where(R <= 0.6 * rmax, 1.0,
                     np.where(R >= rmax, 0.0, 0.5 * (1 + np.cos(np.pi * (R - 0.6 * rmax) / (0.4 * rmax)))))
    data = exposure.image.array
    A_, B_ = idx[0::2], idx[1::2]
    mA = _star_mask(data.shape, X0, Y0, S["xy"][A_], rmax)
    mB = _star_mask(data.shape, X0, Y0, S["xy"][B_], rmax)

    def corrected(f, sel, kbase):
        d, e = _residual_stack(data - f["model_full"], f["good"], S["xy"], f["amps"], sel, X0, Y0,
                               half, with_err=True)
        d, e = np.nan_to_num(d), np.nan_to_num(e, nan=np.inf)
        # keep only significant pixels (soft threshold): the outer, noise-dominated part of the
        # median punched holes into the kernel
        d = np.sign(d) * np.clip(np.abs(d) - P["psf_soft_sigma"] * e, 0.0, None) * taper
        if P["psf_smooth_px"]:
            d = gaussian_filter(d, P["psf_smooth_px"])
        k = kbase.copy()
        sub = k[H - half:H + half + 1, H - half:H + half + 1]
        # never remove more than 70% of the PSF at a pixel (no holes)
        k[H - half:H + half + 1, H - half:H + half + 1] = np.maximum(sub + d, 0.3 * sub)
        if P.get("psf_regularize"):
            k = smooth_psf_kernel(k, fwhm, P["psf_reg_core_fwhm"], P["psf_reg_sigma"])
        return k / k.sum()

    cur = psf
    for cyc in range(max(1, P["psf_core_cycles"])):
        try:
            fA = _fit(exposure, L.kernel_to_lsst_psf(corrected(fit, A_, kernel)), S, bad_bitmask,
                      fit["model_full"], fitkw)
            fB = _fit(exposure, L.kernel_to_lsst_psf(corrected(fit, B_, kernel)), S, bad_bitmask,
                      fit["model_full"], fitkw)
        except ValueError as e:
            print(f"    PSF core cycle {cyc + 1}: {e}")
            break
        before = _chi2_on(exposure, fit, mB) + _chi2_on(exposure, fit, mA)
        after = _chi2_on(exposure, fA, mB) + _chi2_on(exposure, fB, mA)
        gain = (before - after) / max(before, 1e-12)
        keep = gain >= P["psf_cv_min_gain"]
        if verbose:
            print(f"    PSF core cycle {cyc + 1}: median of {len(idx)} bright stars, held-out "
                  f"chi2/px {before / 2:.2f} -> {after / 2:.2f} ({100 * gain:+.1f}%) "
                  f"({'kept' if keep else 'rejected'})")
        if not keep:
            break
        kernel = corrected(fit, idx, kernel)
        cur = L.kernel_to_lsst_psf(kernel)
        fit = _fit(exposure, cur, S, bad_bitmask, fit["model_full"], fitkw)
    return cur


def smooth_psf_kernel(kernel, fwhm, core_fwhm=1.5, sigma=1.0, max_pattern=1.5):
    """Hole-free, smooth version of a PSF kernel.

    Inside core_fwhm*FWHM the kernel is kept (well measured).  Outside, it is rebuilt as
        radial profile(r) * exp(angular pattern)
    with the profile = 1-px ring MEDIANS (a few holed pixels cannot pull it down), forced to
    decrease monotonically outwards, and the pattern = log(kernel / profile) median-filtered
    (5x5: removes isolated holes / hot pixels) and Gaussian-smoothed (sigma px), clipped to
    +-max_pattern in log (factor ~4.5).  Smooth, real structure (ellipticity, lobes, spikes
    wider than ~2 px) survives; pixel-scale defects do not.  Blended over 2 px, renormalised."""
    from scipy.ndimage import median_filter, gaussian_filter
    k = np.asarray(kernel, float)
    H = k.shape[0] // 2
    u = np.arange(-H, H + 1)
    R = np.hypot(*np.meshgrid(u, u))
    floor = max(k.max(), 1e-30) * 1e-9
    rb = np.round(R).astype(int)
    nr = rb.max() + 1
    prof = np.array([np.median(k[rb == i]) if np.any(rb == i) else np.nan for i in range(nr)])
    r0 = core_fwhm * fwhm
    i0 = int(np.ceil(r0))
    for i in range(i0 + 1, nr):                    # monotone non-increasing outside the core
        if np.isfinite(prof[i - 1]):
            prof[i] = min(prof[i], prof[i - 1]) if np.isfinite(prof[i]) else prof[i - 1]
    prof = np.maximum(np.nan_to_num(prof, nan=floor), floor)
    base = np.exp(np.interp(R, np.arange(nr), np.log(prof)))
    pat = np.log(np.maximum(k, floor) / base)
    pat = gaussian_filter(median_filter(pat, size=5, mode="nearest"), sigma, mode="nearest")
    smooth = base * np.exp(np.clip(pat, -max_pattern, max_pattern))
    w = np.clip((r0 + 2.0 - R) / 2.0, 0.0, 1.0)
    out = w * k + (1.0 - w) * smooth
    out = np.clip(out, 0.0, None)
    return out / out.sum()


def regularize_psf(psf, point, P, verbose=True, label=""):
    """LSST PSF with its kernel passed through smooth_psf_kernel (no-op if psf_regularize off)."""
    if not P.get("psf_regularize"):
        return psf
    k = psf.computeKernelImage(point).array.astype(float)
    fwhm = L._fwhm(psf, point)
    new = smooth_psf_kernel(k, fwhm, P["psf_reg_core_fwhm"], P["psf_reg_sigma"])
    if verbose:
        H = k.shape[0] // 2
        u = np.arange(-H, H + 1)
        out = np.hypot(*np.meshgrid(u, u)) > P["psf_reg_core_fwhm"] * fwhm
        ks = k / k.sum()
        print(f"  PSF regularised{label}: {100 * np.abs(new - ks)[out].sum():.2f}% of the flux "
              f"moved outside {P['psf_reg_core_fwhm'] * fwhm:.1f}px (holes / pixel noise)")
    return L.kernel_to_lsst_psf(new)


def _radial_profile(img, ok, x, y, X0, Y0, edges):
    """Median of img in annuli around (x, y) (detector px); NaN where < 5 good pixels."""
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


def _azimuthal_mean(k):
    """Image of the azimuthal mean of k (1-px annuli about the centre)."""
    h = k.shape[0] // 2
    u = np.arange(-h, h + 1)
    rb = np.round(np.hypot(*np.meshgrid(u, u))).astype(int)
    prof = np.bincount(rb.ravel(), k.ravel()) / np.maximum(np.bincount(rb.ravel()), 1)
    return prof[rb]


def _angular_pattern(k, support):
    """k / (its ring average), i.e. the angular modulation (1 = circular); 1 beyond `support`
    px or where the ring average is not positive.  Clipped to [0, 50]."""
    h = k.shape[0] // 2
    u = np.arange(-h, h + 1)
    R = np.hypot(*np.meshgrid(u, u))
    az = _azimuthal_mean(k)
    pat = np.where(az > 1e-12 * max(k.max(), 1e-30), k / np.where(az > 0, az, 1.0), 1.0)
    pat = np.where(R <= support - 1, pat, 1.0)
    return np.clip(pat, 0.0, 50.0)


def extend_psf_wings_keep(psf, point, size=81, blend=(12.0, 16.0), keep_angular=False, verbose=True):
    """Like cutout_fit_lib.extend_psf_wings (Moffat tail out to size x size), but only the
    kernel's azimuthal mean is replaced by the Moffat beyond the blend radius -- diffraction
    spikes and other non-circular structure in the kernel are kept."""
    kernel = psf.computeKernelImage(point).array.astype(float).copy()
    try:
        _, info = L.extend_kernel_with_moffat(kernel, new_size=size)
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
    pat = _angular_pattern(emb, hh) if keep_angular else 1.0
    new = w * emb + (1.0 - w) * A * (1.0 + (R / a) ** 2) ** (-b) * pat
    new /= new.sum()
    if verbose:
        print(f"  PSF wings extended{' (angular pattern kept)' if keep_angular else ''}: "
              f"{kernel.shape[0]}px -> {size}px; "
              f"Moffat a={a:.2f}px, beta={b:.2f}")
    return L.kernel_to_lsst_psf(new)


def wing_profile_psf(exposure, res, psf, P, verbose=True):
    """PSF whose outer wings are MEASURED from the stars.

    Neighbour-subtracted radial profiles: unsaturated bright stars (peak S/N >= wing_min_snr),
    divided by their fitted flux, give the normalised profile where they have signal;
    saturated stars (unknown flux) are scaled to that profile on the radii where both are
    measured and extend it outwards.  A Moffat is fitted to the combined profile beyond
    1.5 FWHM and replaces the kernel beyond ~2.5 FWHM (core kept), then renormalised.
    Returns a new PSF, or None if there are too few stars."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)      # all-NaN annuli are expected
        return _wing_profile_psf(exposure, res, psf, P, verbose)


def _wing_profile_psf(exposure, res, psf, P, verbose=True):
    from scipy.optimize import curve_fit
    bb = exposure.getBBox()
    c = bb.getCenter()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    kernel = psf.computeKernelImage(c).array.astype(float).copy()
    size = int(P["psf_wing_size"] or kernel.shape[0]) | 1
    if kernel.shape[0] < size:                         # embed a smaller kernel
        big = np.zeros((size, size))
        o = (size - kernel.shape[0]) // 2
        big[o:o + kernel.shape[0], o:o + kernel.shape[0]] = kernel
        kernel = big
    H = kernel.shape[0] // 2
    fwhm = L._fwhm(psf, c)
    fit = res["snaps"][-1]["fit"]
    S = res["S"]
    amps, model, good = fit["amps"], fit["model_full"], fit["good"]
    data = exposure.image.array
    psf_peak = kernel.max()
    noise = float(np.sqrt(np.nanmedian(exposure.variance.array)))
    core_lab = sat_core_labels(exposure, fwhm, P["sat_zone_min_px"])
    core = _label_at(S["xy"], core_lab, X0, Y0)
    snr = amps * psf_peak / noise
    edges = np.arange(0.0, H + 0.01, 1.0)
    rc = 0.5 * (edges[1:] + edges[:-1])
    U, Sat = [], []
    for j in np.flatnonzero(((snr >= P["wing_min_snr"]) | (core > 0)) & (S["ext"] == 0)):
        cj = fit["comps"][j]
        if cj is None or amps[j] <= 0:
            continue
        own = np.zeros_like(data)
        own[cj["ys"], cj["xs"]] = amps[j] * cj["stamp"]
        prof = _radial_profile(data - model + own, good, S["xy"][j, 0], S["xy"][j, 1], X0, Y0, edges)
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
    n_sat_used = 0
    if Sat:
        prof_s = []
        for sp in Sat:
            m = reliable & np.isfinite(sp) & (rc >= 1.5 * fwhm)
            if m.sum() < 3:
                continue
            scale = np.sum(sp[m] * u[m]) / np.sum(u[m] ** 2)
            if scale > 0:
                prof_s.append(sp / scale)
        if prof_s:
            n_sat_used = len(prof_s)
            ps = np.nanmedian(np.array(prof_s), axis=0)
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
    r0, r1 = P["wing_blend_fwhm"][0] * fwhm, P["wing_blend_fwhm"][1] * fwhm
    w = np.clip((r1 - Rg) / (r1 - r0), 0.0, 1.0)
    # the Moffat replaces only the RADIAL profile; the kernel's angular pattern (spikes, lobes)
    # relative to its ring average is kept
    pat = _angular_pattern(kernel, H) if P["wing_keep_angular"] else 1.0
    new = w * kernel + (1.0 - w) * moffat(Rg, A, a, b) * pat
    new /= new.sum()
    tail = float(new[Rg > r1].sum())
    if verbose:
        print(f"  wings from stars: {len(U)} bright + {n_sat_used} saturated star(s); profile "
              f"measured to {rc[m].max():.0f}px; Moffat a={a:.2f}px beta={b:.2f}; "
              f"{100 * tail:.2f}% of flux beyond {r1:.0f}px (was "
              f"{100 * float(kernel[Rg > r1].sum() / kernel.sum()):.2f}%)")
    return L.kernel_to_lsst_psf(new)


_SPATIAL = {}      # id(psf) -> (psf, Kx, Ky): gradient terms of a spatially varying PSF


def spatial_kernel_psf(K0, Kx, Ky, xc, yc, scale):
    """LSST PSF  K(x, y) = K0 + (x-xc)/scale * Kx + (y-yc)/scale * Ky  (a LinearCombinationKernel
    with linear spatial functions -- the same machinery as the PCA PSF)."""
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


def _psf_gradient(data, model, var, good, xy0, pos, amps, kernel, xc, yc, scale,
                  max_radius=20, soft_sigma=1.0):
    """Per-pixel weighted fit of the stars' normalised residual stamps:
       R_j(u) = a(u) + bx(u) * dx_j + by(u) * dy_j,   dx_j = (x_j - xc)/scale, ...
    Returns (bx, by, n_stars): the PSF change per unit (scaled) x and y, soft-thresholded,
    tapered beyond max_radius, zero-sum, with position-degenerate parts removed."""
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
    he = min(h, max_radius)
    win = np.clip((he - rr) / max(0.4 * he, 1e-9), 0, 1)
    win = 0.5 - 0.5 * np.cos(np.pi * win)
    gy, gx = np.gradient(kernel)
    A = np.column_stack([gx.ravel(), gy.ravel()])
    for k in (1, 2):
        b = sol[:, k]
        e = err[:, k - 1]
        b = np.sign(b) * np.maximum(np.abs(b) - soft_sigma * e, 0.0)
        b = (b.reshape(kernel.shape)) * win
        coef, *_ = np.linalg.lstsq(A, b.ravel(), rcond=None)
        b = b - (A @ coef).reshape(b.shape)
        b -= (b.sum() / max(win.sum(), 1e-30)) * win          # zero total flux, no pedestal
        out.append(b)
    return out[0], out[1], len(R)


def refine_psf_spatial(exposure, psf_old, psf_new, res, bad_bitmask, P, fitkw, step=0.8,
                       verbose=True):
    """Add (or update) a linear spatial variation to psf_new, measured from the clean stars'
    residuals in `res` (fitted with psf_old).  CROSS-VALIDATED like refine_psf_clean: the
    gradient from half A is judged on half B and vice versa; kept only if the held-out chi2
    drops by >= psf_spatial_min_gain compared with the constant psf_new."""
    from scipy.ndimage import gaussian_filter
    bb = exposure.getBBox()
    c = bb.getCenter()
    xc, yc = c.getX(), c.getY()
    scale = 0.5 * max(bb.getWidth(), bb.getHeight())
    fwhm = L._fwhm(psf_new, c)
    K0 = psf_new.computeKernelImage(c).array.astype(float).copy()
    noise = float(np.sqrt(np.nanmedian(exposure.variance.array)))
    idx = select_psf_stars(exposure, res, fwhm, K0.max(), noise, P, verbose=False)
    S, fit = res["S"], res["snaps"][-1]["fit"]
    pos = S["xy"][idx]
    if len(idx) < P["psf_spatial_min_stars"] or np.std(pos[:, 0]) < 0.15 * scale \
            or np.std(pos[:, 1]) < 0.15 * scale:
        if verbose:
            print(f"  spatial PSF: {len(idx)} clean star(s), too few / too clustered -- constant PSF")
        return psf_new
    prev = _SPATIAL.get(id(psf_old))
    Kx0, Ky0 = (prev[1], prev[2]) if (prev is not None and prev[0] is psf_old
                                      and prev[1].shape == K0.shape) else (np.zeros_like(K0),) * 2
    data = exposure.image.array
    X0, Y0 = bb.getMinX(), bb.getMinY()
    mr = (P["stack_kwargs"] or {}).get("max_radius", 20)

    def grad(sel):
        bx, by, _ = _psf_gradient(data, fit["model_full"], fit["variance"], fit["good"], (X0, Y0),
                                  S["xy"][sel], fit["amps"][sel], K0, xc, yc, scale,
                                  max_radius=mr, soft_sigma=P["psf_soft_sigma"])
        if bx is None:
            return None
        if P["psf_smooth_px"]:
            bx, by = gaussian_filter(bx, P["psf_smooth_px"]), gaussian_filter(by, P["psf_smooth_px"])
        return spatial_kernel_psf(K0, Kx0 + step * bx, Ky0 + step * by, xc, yc, scale)

    A, B = idx[0::2], idx[1::2]
    mA = _star_mask(data.shape, X0, Y0, S["xy"][A], 2.5 * fwhm)
    mB = _star_mask(data.shape, X0, Y0, S["xy"][B], 2.5 * fwhm)
    try:
        cA, cB = grad(A), grad(B)
    except Exception as e:
        print(f"  spatial PSF: could not build ({type(e).__name__}: {e}) -- constant PSF")
        return psf_new
    if cA is None or cB is None:
        return psf_new
    f0 = _fit(exposure, psf_new, S, bad_bitmask, fit["model_full"], fitkw)
    fA = _fit(exposure, cA, S, bad_bitmask, fit["model_full"], fitkw)
    fB = _fit(exposure, cB, S, bad_bitmask, fit["model_full"], fitkw)
    before = _chi2_on(exposure, f0, mB) + _chi2_on(exposure, f0, mA)
    after = _chi2_on(exposure, fA, mB) + _chi2_on(exposure, fB, mA)
    gain = (before - after) / max(before, 1e-12)
    keep = gain >= P["psf_spatial_min_gain"]
    if verbose:
        print(f"  spatial PSF from {len(idx)} clean stars (cross-validated): held-out chi2/px "
              f"{before / 2:.2f} (constant) -> {after / 2:.2f} (varying), {100 * gain:+.1f}% "
              f"-> {'KEPT' if keep else 'not kept'}")
    if not keep:
        return psf_new
    return grad(idx)


def fit_brightness_term(exposure, psf, res, bad_bitmask, P, fitkw, verbose=True):
    """Measure / update the brightness-dependent PSF term (brighter-fatter).

    Isolated, unsaturated stars with peak S/N >= bf_min_snr: 1-px ring means of their
    residual / flux out to bf_radius_fwhm*FWHM.  Per ring, a weighted fit
        residual/flux = A(r) + B(r) * (log10 flux - L0)
    gives the slope B(r): how the profile changes per dex of brightness.  B is smoothed,
    tapered to 0 at the outer radius and made flux-neutral (fluxes stay fluxes); the
    intercept A is left to the PSF itself.  The update is built from half of the stars and
    judged on the other half (both ways); kept only if the held-out chi2 improves.
    Returns the fraction of flux moved per dex (0 if nothing changed)."""
    from scipy.ndimage import gaussian_filter1d
    bb = exposure.getBBox()
    c = bb.getCenter()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    kernel = psf.computeKernelImage(c).array.astype(float)
    fwhm = L._fwhm(psf, c)
    noise = float(np.sqrt(np.nanmedian(exposure.variance.array)))
    S = res["S"]
    fit = _fit(exposure, psf, S, bad_bitmask, res["snaps"][-1]["fit"]["model_full"], fitkw)
    Pq = dict(P, psf_stack_min_snr=P["bf_min_snr"], psf_stack_n=400)
    idx = select_psf_stars(exposure, dict(res, snaps=[dict(fit=fit)]), fwhm, kernel.max(),
                           noise, Pq, verbose=False)
    amps = fit["amps"]
    if len(idx) < 8 or np.ptp(np.log10(amps[idx])) < 0.5:
        if verbose:
            print(f"  brightness term: {len(idx)} isolated star(s) with S/N >= {P['bf_min_snr']:g} "
                  f"(need >= 8 spanning >= 0.5 dex) -- skipped")
        return 0.0
    idx = idx[np.argsort(amps[idx])]
    lf = np.log10(amps[idx])
    if _BF["L0"] is None:          # the PSF core is tuned on the S/N >= psf_core_min_snr stars
        snr = amps[idx] * kernel.max() / noise
        core_set = lf[snr >= P["psf_core_min_snr"]]
        _BF["L0"] = float(np.median(core_set if len(core_set) >= 3 else lf))
    rmax = P["bf_radius_fwhm"] * fwhm
    nr = int(np.ceil(rmax))
    rc = np.arange(nr, dtype=float)                 # ring i = pixels with round(r) == i
    H = kernel.shape[0] // 2
    u = np.arange(-H, H + 1)
    rbk = np.floor(np.hypot(*np.meshgrid(u, u)) + 0.5).astype(int)
    area = np.array([np.sum(rbk == i) for i in range(nr)], float)
    Kring = np.array([kernel[rbk == i].mean() for i in range(nr)]) / kernel.sum()
    data = exposure.image.array
    resid = data - fit["model_full"]
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

    old_r, old_p = _BF["r"], _BF["prof"]
    base = np.interp(rc, old_r, old_p, right=0.0) if old_p is not None else np.zeros(nr)

    def use(r, prof):
        _BF.update(r=r, prof=prof)
        _set_bf_fluxes(S["xy"], amps)
        _STAMP_CACHE.clear()

    iA, iB = np.arange(0, len(idx), 2), np.arange(1, len(idx), 2)
    mA = _star_mask(data.shape, X0, Y0, S["xy"][idx[iA]], rmax)
    mB = _star_mask(data.shape, X0, Y0, S["xy"][idx[iB]], rmax)
    gain, before, after = -1.0, np.nan, np.nan
    try:
        use(rc, base + slope(iA))
        fA = _fit(exposure, psf, S, bad_bitmask, fit["model_full"], fitkw)
        use(rc, base + slope(iB))
        fB = _fit(exposure, psf, S, bad_bitmask, fit["model_full"], fitkw)
        before = _chi2_on(exposure, fit, mB) + _chi2_on(exposure, fit, mA)
        after = _chi2_on(exposure, fA, mB) + _chi2_on(exposure, fB, mA)
        gain = (before - after) / max(before, 1e-12)
    except ValueError as e:
        print(f"  brightness term: {e}")
    if gain < P["psf_cv_min_gain"]:
        if old_p is None:
            _BF.update(r=None, prof=None)
            _STAMP_CACHE.clear()
        else:
            use(old_r, old_p)
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


def _kernel_change(psf_old, psf_new, point):
    a = psf_old.computeKernelImage(point).array
    b = psf_new.computeKernelImage(point).array
    h = min(a.shape[0], b.shape[0]) // 2
    ca_, cb_ = a.shape[0] // 2, b.shape[0] // 2
    a = a[ca_ - h:ca_ + h + 1, ca_ - h:ca_ + h + 1]
    b = b[cb_ - h:cb_ + h + 1, cb_ - h:cb_ + h + 1]
    return float(np.linalg.norm(b - a) / max(np.linalg.norm(a), 1e-30))


def print_timing(timer=None):
    timer = timer or L.TIMER
    tot, cnt = defaultdict(float), defaultdict(int)
    for r in timer.records:
        if r["kind"] == "compute":
            tot[r["label"]] += r["seconds"]
            cnt[r["label"]] += 1
    print(f"\n{'step':40s} {'calls':>6s} {'total [s]':>10s}")
    for k, v in sorted(tot.items(), key=lambda kv: -kv[1])[:20]:
        print(f"{k:40s} {cnt[k]:6d} {v:10.2f}")
    print(f"wall clock: {time.perf_counter() - timer.t_start:.1f}s")


def run_seeded_pipeline(exposure_crop, psf, bad_bitmask, seeds=None, seed_gaia_index=None,
                        starCat=None, bgImage=None, viz=False, seed_gmag=None, **cfg):
    """Gaia-seeded, adaptive replacement for cutout_fit_lib.run_cutout_pipeline on the
    WHOLE exposure_crop.  seeds: (N, 2) detector-pixel positions (None/empty = no seeds).
    cfg: run_cutout_pipeline keys + SEEDED_DEFAULTS keys.  Returns outer_history (same
    layout as run_cutout_pipeline, plus 'psf_change' and result['adaptive'])."""
    P = _params(cfg)
    _FAST["on"] = bool(P["fast_stamps"])
    _FAST["cache"].clear()
    _bf_reset()
    print(f"cutout_gaia {__version__}: brightness term {'ON' if P['bf_term'] else 'off'}, "
          f"PSF core median {'ON' if P['psf_core_median'] else 'off'}, "
          f"PSF regularise {'ON' if P['psf_regularize'] else 'off'}")
    enable_speedups()
    L.TIMER.reset()
    cutout = exposure_crop.clone()
    cutout.setPsf(psf)
    cb = cutout.getBBox()
    c = cb.getCenter()
    if P["sat_grow_px"]:
        from scipy.ndimage import binary_dilation
        sat = (cutout.mask.array & cutout.mask.getPlaneBitMask("SAT")) != 0
        if sat.any():
            grown = binary_dilation(sat, iterations=int(P["sat_grow_px"])) & ~sat
            cutout.mask.addMaskPlane("SAT_GROWN")
            bit = cutout.mask.getPlaneBitMask("SAT_GROWN")
            cutout.mask.array[grown] |= bit
            bad_bitmask = bad_bitmask | bit
            print(f"excluding {int(grown.sum())} px within {P['sat_grow_px']}px of saturated pixels")
    if P["near_sat_frac"]:
        lvl, how = saturation_level(cutout, P)
        if lvl is None:
            print(f"near-saturated cores: {how}")
        else:
            from scipy.ndimage import binary_dilation, label
            dat = cutout.image.array
            sat = (cutout.mask.array & cutout.mask.getPlaneBitMask("SAT")) != 0
            ns = np.isfinite(dat) & (dat > P["near_sat_frac"] * lvl) & ~sat
            if ns.any():
                if P["near_sat_grow_px"]:
                    ns = binary_dilation(ns, iterations=int(P["near_sat_grow_px"]))
                ns &= (cutout.mask.array & bad_bitmask) == 0
                cutout.mask.addMaskPlane("NEARSAT")
                bit = cutout.mask.getPlaneBitMask("NEARSAT")
                cutout.mask.array[ns] |= bit
                bad_bitmask = bad_bitmask | bit
                n_obj = label(ns)[1]
            else:
                n_obj = 0
            print(f"near-saturated cores: saturation level {lvl:.4g} ({how}); excluding "
                  f"{int(ns.sum())} px above {P['near_sat_frac']:g} x that level in {n_obj} "
                  f"star core(s) -- fitted from their wings")
    P["_sat_zones"], P["_sat_zone_mask"] = sat_zones(cutout, L._fwhm(psf, c), P)
    P["_halo_rejects"] = np.zeros((0, 2))
    if P["psf_wing_size"]:
        psf = extend_psf_wings_keep(psf, c, P["psf_wing_size"], blend=P["psf_ext_blend"],
                                    keep_angular=P["wing_keep_angular"])
    psf = regularize_psf(psf, c, P, label=" (initial)")
    raw = cutout.clone()
    exp_w = cutout.clone()
    if P["subtract_initial_bg"] and bgImage is not None:
        try:
            exp_w.image.array[:] -= bgImage[cb].array
            print("pass 1 uses the detector-level background (subtracted)")
        except Exception as e:
            print(f"(could not subtract initial background: {e})")
    exp_w.setPsf(psf)

    S = _state(seeds if seeds is not None else None, seed=True, gaia=seed_gaia_index,
               gmag=seed_gmag)
    if len(S["xy"]) and P["fit_gaia_offset"]:
        noise = float(np.sqrt(np.nanmedian(exp_w.variance.array)))
        dx, dy, n = estimate_gaia_offset(exp_w, S, L._fwhm(psf, c), bad_bitmask, noise)
        if n:
            S["xy"] += (dx, dy)
            S["ref"] += (dx, dy)
    print(f"{len(S['xy'])} Gaia seed(s) for cutout {cb}")

    starpos = None
    if starCat is not None and len(starCat):
        starpos = [(r.getX(), r.getY()) for r in starCat
                   if cb.contains(L.geom.Point2I(int(round(r.getX())), int(round(r.getY()))))]
    psf_w, bg_w, psf_change = psf, None, np.nan
    hist = []
    for k in range(P["max_outer"]):
        clear_stamp_cache()
        shape_now = P["use_shape_cut"] and (k > 0 or P["shape_in_first_pass"])
        thr0 = (P["seed_threshold_start"] if len(S["xy"]) else P["threshold_start"]) if k == 0 \
            else P["threshold_min"]
        print(f"\n{'=' * 70}\n===== Outer pass {k + 1} (cap {P['max_outer']}), PSF FWHM "
              f"{L._fwhm(psf_w, c):.3f}px, shape cut {'ON' if shape_now else 'off'} =====")
        t0 = time.perf_counter()
        with L.TIMER.step("INNER LOOP (total)", outer=k, kind="wall"):
            res = run_inner_seeded(exp_w, psf_w, bad_bitmask, S, P, outer_idx=k,
                                   threshold_start=thr0, shape_now=shape_now)
        S = res["S"]
        res["adaptive"]["seconds"] = time.perf_counter() - t0
        hist.append(dict(psf=psf_w, exposure=exp_w, result=res, version=__version__,
                         bf=None if _BF["prof"] is None else
                         dict(r=_BF["r"].copy(), prof=_BF["prof"].copy(), L0=_BF["L0"]),
                         resid_std=res["resid_std_history"][-1], n_peaks=len(res["peaks"]),
                         bg_image=bg_w, psf_change=psf_change, bad_bitmask=bad_bitmask,
                         poisson_k=P["poisson_k"]))
        if P["auto_poisson"] and res["snaps"]:
            try:
                import contextlib
                import io
                import types
                fake = types.SimpleNamespace(outer_history=[hist[-1]], bad_bitmask=bad_bitmask,
                                             config=cfg)
                with contextlib.redirect_stdout(io.StringIO()):
                    dg = diagnose_residual(fake, plot=False)
                share = dg.get("noise_share", 0.0)
                hf_r = dg["highfreq"]["ratio"]
                k_new = float(dg["highfreq"]["k"])
                sz = dg["bright-star disks"]["sigma_z"]
                ratio = share
                if share >= 0.5 and len(hf_r) and hf_r[-1] > 1.15 and np.isfinite(k_new) \
                        and k_new > 0:
                    P["poisson_k"] = max(P["poisson_k"], k_new)
                    hist[-1]["poisson_k"] = P["poisson_k"]          # used by sigma='effective'
                    print(f"  noise check: bright-star excess is pixel-scale noise ({100 * share:.0f}%) "
                          f"-> variance += {P['poisson_k']:.4g} x model from now on")
                else:
                    print(f"  noise check: bright-star disks sigma_z {sz:.2f}; missing photon noise "
                          f"explains {100 * share:.0f}% of the excess -> variance unchanged")
            except Exception as e:
                print(f"  noise check failed ({type(e).__name__}: {e})")
        if P["model_err_measure"] and res["snaps"]:
            try:
                eps, _ = measure_model_error(exp_w, res, P)
                hist[-1]["model_err"] = eps
                if P["fit_sys_from_data"] and eps > 0:
                    P["sys_frac"] = max(eps, 0.002)
                    print(f"  -> fit weights use sys_frac = {P['sys_frac']:.3g} from now on")
            except Exception as e:
                print(f"  model-error measurement failed ({type(e).__name__}: {e})")
        ad = res["adaptive"]
        n_now = len(res["peaks"])
        growth = (abs(n_now - hist[-2]["n_peaks"]) / max(hist[-2]["n_peaks"], 1)
                  if len(hist) > 1 else np.inf)
        psf_settled = (not P["refit_psf"]) or (np.isfinite(psf_change) and psf_change < P["psf_tol"])
        print(f"\nouter pass {k + 1}: {n_now} sources ({int(res['is_seed'].sum())} Gaia), "
              f"{ad['n_inner']} inner iteration(s), stop: {ad['stop']}, residual "
              f"{'has EXCESS' if ad['missing'] else 'clean'}, source change "
              f"{100 * growth:.1f}%, PSF change into this pass {100 * psf_change:.2f}%, "
              f"{ad['seconds']:.0f}s")
        if viz and res["snaps"]:
            try:
                L.plot_iteration_overview(exp_w, res["snaps"][-1],
                                          dict(resid_std=res["resid_std_history"],
                                               n_peaks=res["n_peaks_history"],
                                               n_zero=res["n_zero_amp_history"]),
                                          k, L._data_norm(exp_w.image.array))
            except Exception as e:
                print(f"(overview plot failed: {e})")
        if k + 1 >= P["max_outer"]:
            if ad["missing"]:
                print(f"WARNING: max_outer={P['max_outer']} reached with excess remaining")
            break
        if k + 1 >= P["min_outer"] and growth < P["outer_growth_tol"] and psf_settled:
            print("outer loop converged: " + ("residual clean, " if not ad["missing"] else
                  "excess remains but nothing changes any more, ") + "sources and PSF stable")
            break
        if not res["snaps"]:
            print("no sources -- stopping")
            break

        # ---- refit background + PSF from this pass's model ---------------------
        exp_n = raw.clone()
        with L.TIMER.step("mask + refit background", outer=k):
            src_mask = L.mask_sources_from_model(exp_n, res["model_full"])
            if P["bg_sat_mask_fwhm"]:
                from scipy.ndimage import distance_transform_edt
                satm = (exp_n.mask.array & exp_n.mask.getPlaneBitMask("SAT")) != 0
                if satm.any():
                    halo = distance_transform_edt(~satm) <= P["bg_sat_mask_fwhm"] * L._fwhm(psf_w, c)
                    exp_n.mask.array[halo] |= exp_n.mask.getPlaneBitMask("DETECTED")
                    print(f"  background refit also ignores saturated-star halos "
                          f"({100 * halo.mean():.1f}% of the cutout)")
            try:
                exp_n, bg_w = L.refit_background(exp_n, use_approx=P["bg_use_approx"],
                                                 bin_size=P["bg_bin_size"])
            except Exception as e:
                print(f"  background refit failed ({type(e).__name__}: {e}) -- keeping image")
        psf_new = psf_w
        exp_n.setPsf(psf_w)
        if P["refit_psf"] == "stack":
            with L.TIMER.step("PSF refit (residual stack)", outer=k):
                try:
                    psf_new = refine_psf_clean(
                        exp_w, psf_w, res, bad_bitmask, P,
                        dict(support_frac=P["support_frac"], sys_frac=P["sys_frac"],
                             poisson_k=P["poisson_k"]),
                        n_cycles=P["psf_polish_cycles"])
                except Exception as e:
                    print(f"  residual-stack PSF refit failed ({type(e).__name__}: {e})")
                    psf_new = psf_w
        elif P["refit_psf"]:
            ks = P["psf_kernel_sizes"] or tuple([P["psf_kernel_size"]] +
                                                 [s for s in (27, 21, 15, 11) if s < P["psf_kernel_size"]])
            with L.TIMER.step("redetect + PSF refit", outer=k):
                try:
                    p_, _, _, cell, _ = L.redetect_and_refit_psf_adaptive(
                        exp_n, kernel_sizes=ks, extra_peak_positions=res["peaks"],
                        force_include_star_positions=starpos)
                    if cell is not None:
                        psf_new = p_
                except Exception as e:
                    print(f"  PSF redetermination failed ({type(e).__name__}: {e})")
        if P["psf_wing_size"]:
            wp = None
            if P["wings_from_stars"]:
                try:
                    wp = wing_profile_psf(exp_w, res, psf_new, P)
                except Exception as e:
                    print(f"  wings from stars failed ({type(e).__name__}: {e})")
            if wp is not None:
                psf_new = wp
            elif P["wings_from_stars"] and psf_new is not psf_w and P["refit_psf"] == "stack":
                pass                       # stack-refined PSF keeps its current (embedded) wings
            elif psf_new is not psf_w:
                psf_new = extend_psf_wings_keep(psf_new, c, P["psf_wing_size"], blend=P["psf_ext_blend"],
                                               keep_angular=P["wing_keep_angular"])
        if P["psf_core_median"]:
            with L.TIMER.step("PSF core (median stack)", outer=k):
                try:
                    psf_new = refine_psf_median(
                        exp_w, psf_new, res, bad_bitmask, P,
                        dict(support_frac=P["support_frac"], sys_frac=P["sys_frac"],
                             poisson_k=P["poisson_k"]),
                        fit0=res["snaps"][-1]["fit"] if psf_new is psf_w else None)
                except Exception as e:
                    print(f"  PSF core (median stack) failed ({type(e).__name__}: {e})")
        if P["psf_spatial"] and P["refit_psf"] == "stack":
            try:
                psf_new = refine_psf_spatial(
                    exp_w, psf_w, psf_new, res, bad_bitmask, P,
                    dict(support_frac=P["support_frac"], sys_frac=P["sys_frac"],
                         poisson_k=P["poisson_k"]))
            except Exception as e:
                print(f"  spatial PSF failed ({type(e).__name__}: {e}) -- constant PSF")
        if psf_new is not psf_w:
            psf_new = regularize_psf(psf_new, c, P)
        psf_change = _kernel_change(psf_w, psf_new, c) if psf_new is not psf_w else 0.0
        if P["bf_term"]:
            with L.TIMER.step("brightness-dependent PSF", outer=k):
                try:
                    for _ in range(2):
                        moved = fit_brightness_term(
                            exp_w, psf_new, res, bad_bitmask, P,
                            dict(support_frac=P["support_frac"], sys_frac=P["sys_frac"],
                                 poisson_k=P["poisson_k"]))
                        psf_change = max(psf_change, moved)     # a change -> not settled yet
                        if moved < P["psf_tol"]:
                            break
                except Exception as e:
                    print(f"  brightness term failed ({type(e).__name__}: {e})")
        exp_n.setPsf(psf_new)
        exp_w, psf_w = exp_n, psf_new
    print_timing()
    return hist


def run_cutout_gaia(raw_crop, psf, bad_bitmask, gaia_tbl=None, config=None, interior_bbox=None,
                    starCat=None, bgImage=None, viz=False, name="cutout", gmag_max=None,
                    **overrides):
    """Seeded adaptive run on ONE cropped exposure. Returns a cutout_notebook.CutoutRun
    whose catalog also has from_gaia / gaia_source_id / gaia_G / gaia_dx / gaia_dy."""
    from . import session as cn
    cfg = cn.make_config(config, **overrides)
    P = _params(cfg)
    cb = raw_crop.getBBox()
    interior_bbox = interior_bbox or cb
    seeds, rows = select_seeds(gaia_tbl, cb, L._fwhm(psf, cb.getCenter()),
                               P["seed_pad_fwhm"], gmag_max, P["seed_min_sep_px"])
    t0 = time.perf_counter()
    gm = (np.asarray(gaia_tbl["phot_g_mean_mag"], float)[rows]
          if gaia_tbl is not None and len(rows) and "phot_g_mean_mag" in gaia_tbl.colnames else None)
    hist = run_seeded_pipeline(raw_crop, psf, bad_bitmask, seeds=seeds, seed_gaia_index=rows,
                               starCat=starCat, bgImage=bgImage, viz=viz, seed_gmag=gm, **cfg)
    res = hist[-1]["result"]
    S = res["S"]
    cat = cn.build_catalog(hist, interior_bbox)
    cat["from_gaia"] = S["seed"]
    gi = S["gaia"]
    has = gi >= 0
    sid = np.full(len(gi), -1, np.int64)
    gG = np.full(len(gi), np.nan)
    if gaia_tbl is not None and has.any():
        sid[has] = np.asarray(gaia_tbl["source_id"], np.int64)[gi[has]]
        if "phot_g_mean_mag" in gaia_tbl.colnames:
            gG[has] = np.asarray(gaia_tbl["phot_g_mean_mag"], float)[gi[has]]
    cat["gaia_source_id"], cat["gaia_G"] = sid, gG
    d = S["xy"] - S["ref"]
    cat["gaia_dx"], cat["gaia_dy"] = d[:, 0], d[:, 1]
    cat["extended"] = S["ext"] > 0
    cat["ext_sigma_px"] = S["ext"]
    summary = dict(name=name, status="ok", n_sources=len(res["peaks"]),
                   n_gaia_seeds=int(len(seeds)), n_from_gaia=int(S["seed"].sum()),
                   n_in_interior=int(np.sum(cat["in_interior"])),
                   bbox_fitted=[cb.getMinX(), cb.getMaxX(), cb.getMinY(), cb.getMaxY()],
                   interior_bbox=[interior_bbox.getMinX(), interior_bbox.getMaxX(),
                                  interior_bbox.getMinY(), interior_bbox.getMaxY()],
                   resid_std_per_outer=[h["resid_std"] for h in hist],
                   n_peaks_per_outer=[h["n_peaks"] for h in hist],
                   psf_fwhm_px_per_outer=[L._fwhm(h["psf"], h["exposure"].getBBox().getCenter())
                                          for h in hist],
                   psf_change_per_outer=[h["psf_change"] for h in hist],
                   adaptive_per_outer=[{k: v for k, v in h["result"]["adaptive"].items()
                                        if k != "final"} for h in hist],
                   residual_clean=not res["adaptive"]["missing"],
                   wall_seconds=time.perf_counter() - t0, pipeline_kwargs=cfg)
    return cn.CutoutRun(name, hist, cat, summary, raw_crop, interior_bbox,
                        hist[0].get("bad_bitmask", bad_bitmask), cfg)


_PAR = {}


def _par_worker(name):
    """Runs in a forked child: inherits the session from _PAR, no pickling of inputs."""
    import pickle
    from . import session as cn
    session, config, out_dir, gmag_max, overrides = _PAR["args"]
    job = session.jobs[name]
    run = run_cutout_gaia(job["raw"], session.psf, session.bad_bitmask, session.gaia, config,
                          interior_bbox=job["interior"], starCat=session.starCat,
                          bgImage=session.bgImage, viz=False, name=name, gmag_max=gmag_max,
                          **overrides)
    path = None
    if out_dir:
        path = os.path.join(out_dir, f"visit{session.visit}_det{session.detector}", name)
        cn.save_run(run, path)
    try:
        pickle.dumps(run)
        return name, run, run.summary, path
    except Exception as e:                     # LSST objects that do not pickle
        print(f"[{name}] result not transferable to the notebook ({e}); saved to {path}")
        return name, None, run.summary, path


def check_bright_neighbours(run, bright_snr=300.0, within_fwhm=4.0, flux_ratio=5.0,
                            n_controls=7, plot=True, max_show=16, verbose=True, similar_min=3):
    """Are the sources sitting on the edges of bright stars real?

    For every source within `within_fwhm` FWHM of a star at least `flux_ratio` times
    brighter (and with peak S/N >= bright_snr), compare its PSF-weighted flux with the flux
    the residual has at CONTROL positions: the same distance from the bright star, rotated
    by 45-degree steps.  A mismatch pattern of the bright star (rings, spikes, lobes) puts
    similar "flux" at several angles; a real star is significantly brighter than its
    controls.  significance = (flux - median(control flux)) / max(flux error, control scatter)

    verdict: 'Gaia'       Gaia star (real by construction)
             'doubtful'   comparable flux (>= half) at >= similar_min (3) control angles:
                          repeating halo pattern
             'real'       significance >= 5 and at most 1 similar angle
             'probable'   3 <= significance (or 2 similar angles -- can be chance)
             'probable'   3 <= significance < 5
             'doubtful'   significance < 3  -- indistinguishable from the bright star's halo
    Returns an astropy Table (also shows thumbnails of the non-Gaia ones if plot=True)."""
    final = run.outer_history[-1]
    exp, res = final["exposure"], final["result"]
    S, fit = res["S"], res["snaps"][-1]["fit"]
    bb = exp.getBBox()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    data = exp.image.array
    model = fit["model_full"]          # PSF-only joint model: bright-star corrections would hide
    var, good = exp.variance.array, fit["good"]   # the repeating halo pattern at the controls
    psf = exp.getPsf()
    c = bb.getCenter()
    fwhm = L._fwhm(psf, c)
    psf_peak = psf.computeKernelImage(c).array.max()
    noise = float(np.sqrt(np.nanmedian(var)))
    amps = fit["amps"]
    snr = amps * psf_peak / noise
    xy = S["xy"]
    bright = np.flatnonzero(snr >= bright_snr)
    rows = []
    if not len(bright):
        print("no bright stars")
        return Table()
    btree = cKDTree(xy[bright])
    alltree = cKDTree(xy)
    resid = data - model
    R2 = (1.5 * fwhm) ** 2

    def psf_flux(img, x, y):
        """PSF-weighted amplitude (and its error) of img at (x, y) over good pixels."""
        st, sb = L.psf_stamp(psf, x, y, bb, taper=L.WING_TAPER)
        if st is None:
            return np.nan, np.nan
        ys, xs = L.bbox_to_slices(sb, bb)
        yy, xx = np.mgrid[sb.getMinY():sb.getMaxY() + 1, sb.getMinX():sb.getMaxX() + 1]
        m = good[ys, xs] & (((xx - x) ** 2 + (yy - y) ** 2) < R2) & (var[ys, xs] > 0)
        if m.sum() < 8:
            return np.nan, np.nan
        w = 1.0 / var[ys, xs][m]
        p = st[m]
        den = np.sum(p * p * w)
        return float(np.sum(p * img[ys, xs][m] * w) / den), float(1.0 / np.sqrt(den))

    for j in range(len(xy)):
        if snr[j] >= bright_snr:
            continue
        for k in btree.query_ball_point(xy[j], within_fwhm * fwhm):
            b = bright[k]
            if amps[b] < flux_ratio * max(amps[j], 1e-9):
                continue
            cj = fit["comps"][j]
            own = np.zeros_like(data)
            if cj is not None:
                own[cj["ys"], cj["xs"]] = amps[j] * cj["stamp"]
            a_j, e_j = psf_flux(resid + own, *xy[j])
            dx, dy = xy[j] - xy[b]
            d = float(np.hypot(dx, dy))
            th0 = np.arctan2(dy, dx)
            ctrl = []
            for q in range(1, n_controls + 1):
                th = th0 + q * 2 * np.pi / (n_controls + 1)
                cx, cy = xy[b, 0] + d * np.cos(th), xy[b, 1] + d * np.sin(th)
                if not (X0 <= cx <= bb.getMaxX() and Y0 <= cy <= bb.getMaxY()):
                    continue
                # a control on top of another (non-negligible) source measures THAT source's
                # mismatch, not the bright star's halo pattern -> skip it
                near = [k2 for k2 in alltree.query_ball_point((cx, cy), 1.5 * fwhm)
                        if k2 not in (j, b) and amps[k2] > 0.3 * max(amps[j], 1e-9)]
                if near:
                    continue
                a_c, _ = psf_flux(resid, cx, cy)
                if np.isfinite(a_c):
                    ctrl.append(a_c)
            ctrl = np.array(ctrl)
            n_sim = 0
            if len(ctrl) >= 3 and np.isfinite(a_j):
                med = float(np.median(ctrl))
                spread = float(1.4826 * np.median(np.abs(ctrl - med)))
                sig = (a_j - med) / max(e_j, spread, 1e-9)
                n_sim = int((ctrl >= 0.5 * a_j).sum()) if a_j > 0 else len(ctrl)
            else:
                med, spread, sig = np.nan, np.nan, np.nan
            if S["seed"][j]:
                verdict = "Gaia"
            elif not np.isfinite(sig):
                verdict = "n/a"
            elif n_sim >= similar_min:   # the same "flux" repeats at other angles: halo pattern
                verdict = "doubtful"
            elif sig >= 5 and n_sim < 2:
                verdict = "real"
            elif sig >= 3:              # (two similar angles can be chance in a crowded field)
                verdict = "probable"
            else:
                verdict = "doubtful"
            rows.append((j, xy[j, 0], xy[j, 1], int(b), d / fwhm, amps[j] / amps[b], snr[j],
                         a_j, med, spread, sig, len(ctrl), n_sim, verdict))
            break
    t = Table(rows=rows, names=("id", "x", "y", "bright_id", "dist_fwhm", "flux_ratio",
                                "peak_snr", "flux", "control_median", "control_scatter",
                                "significance", "n_controls", "n_similar", "verdict")) if rows else Table()
    if not len(t):
        print("no sources next to bright stars")
        return t
    if verbose:
        cnt = {v: int((t["verdict"] == v).sum()) for v in ("Gaia", "real", "probable", "doubtful", "n/a")}
        print(f"{len(t)} source(s) within {within_fwhm:g} FWHM of a >= {flux_ratio:g}x brighter star "
              f"(peak S/N >= {bright_snr:g}): " + ", ".join(f"{v} {n}" for v, n in cnt.items() if n))
        print(f"{'id':>5s} {'x':>8s} {'y':>8s} {'d[FWHM]':>7s} {'f/f_b':>7s} {'S/N':>6s} "
              f"{'signif':>6s} {'similar':>7s}  verdict   (similar = control angles with >= half "
              f"its flux)")
        for r in sorted(t, key=lambda r: (r["verdict"] != "doubtful", np.nan_to_num(r["significance"]))):
            print(f"{r['id']:5d} {r['x']:8.1f} {r['y']:8.1f} {r['dist_fwhm']:7.2f} {r['flux_ratio']:7.3f} "
                  f"{r['peak_snr']:6.1f} {r['significance']:6.1f} {r['n_similar']:4d}/{r['n_controls']:<2d}"
                  f"  {r['verdict']}")
    if plot:
        import matplotlib.pyplot as plt
        from astropy.visualization import simple_norm
        show = [r for r in sorted(t, key=lambda r: np.nan_to_num(r["significance"]))
                if r["verdict"] != "Gaia"][:max_show]
        if show:
            h = int(round(3 * fwhm))
            fig, axes = plt.subplots(len(show), 3, figsize=(10, 3.2 * len(show)), squeeze=False)
            for ax_row, r in zip(axes, show):
                j, b = int(r["id"]), int(r["bright_id"])
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
                        (ax_row[0], dwin, dict(cmap=L.CMAP, norm=simple_norm(dwin, stretch="asinh",
                                                                               percent=99.5)), "data"),
                        (ax_row[1], rw, dict(cmap="RdBu_r", vmin=-lim, vmax=lim),
                         "residual WITHOUT this source"),
                        (ax_row[2], resid[ya:yb_, xa:xb_], dict(cmap="RdBu_r", vmin=-lim, vmax=lim),
                         "residual WITH it")):
                    ax.imshow(img, extent=ext, origin="lower", **kw)
                    ax.plot(r["x"], r["y"], "+", color="k", ms=10, mew=1.5)
                    ax.set_title(ttl, fontsize=9)
                    ax.set_xticks([]); ax.set_yticks([])
                    ax.annotate("", xy=(xy[b, 0], xy[b, 1]), xytext=(r["x"], r["y"]),
                                arrowprops=dict(arrowstyle="->", color="0.4", lw=0.8),
                                annotation_clip=True)
                ax_row[0].set_ylabel(f"#{j}  {r['verdict']}\nsignif {r['significance']:.1f}",
                                     fontsize=9)
            fig.suptitle("Sources next to bright stars (least significant first; arrow points to "
                         "the bright star)", fontsize=11)
            plt.tight_layout()
            plt.show()
    return t


def diagnose_residual(run, bright_snr=300.0, disk_fwhm=3.0, plot=True):
    """Is what is left around bright stars NOISE or MODEL ERROR?

    For three regions -- empty sky, bright-star disks (within disk_fwhm FWHM of stars with
    peak S/N >= bright_snr; bad / grey / saturated-zone pixels excluded) and everything else
    on sources -- it reports
      sigma_z    robust std of residual/sigma (1 = matches the variance plane)
      coherence  robust std of the 3x3-averaged residual/sigma, x3, divided by sigma_z, squared.
                 White (pixel-to-pixel) noise gives 1; correlated noise a bit more (see the
                 sky value); smooth structure the model misses gives clearly more.
    Bright disks with sigma_z > 1 but coherence ~ the sky value = pixel-scale noise (the
    star's own photon noise, larger than the variance plane says): no model can remove it.
    Coherence clearly above the sky value = real structure left to model.
    Also prints the pixel-to-pixel scatter vs model flux (k = extra variance per unit flux).
    Returns a dict of the numbers."""
    from scipy.ndimage import uniform_filter
    final = run.outer_history[-1]
    exp, res = final["exposure"], final["result"]
    S = res["S"]
    fit = res["snaps"][-1]["fit"]
    bb = exp.getBBox()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    data, model, var = exp.image.array, res["model_full"], exp.variance.array
    bad = ((exp.mask.array & run.bad_bitmask) != 0) | ~np.isfinite(data - model) | ~(var > 0)
    fwhm = L._fwhm(exp.getPsf(), bb.getCenter())
    resid = data - model
    excl = unmodeled_near_mask(resid, var, bad, fwhm)
    try:
        P = _params(run.config)
        excl |= sat_zones(exp, fwhm, P)[1] if P["sat_zone_fwhm"] else False
    except Exception:
        pass
    good = ~bad & ~excl
    z = np.where(good, resid / np.sqrt(np.where(var > 0, var, 1.0)), 0.0)
    noise = float(np.sqrt(np.nanmedian(var[good])))
    psf_peak = exp.getPsf().computeKernelImage(bb.getCenter()).array.max()
    snr = fit["amps"] * psf_peak / noise
    disks = _star_mask(data.shape, X0, Y0, S["xy"][snr >= bright_snr], disk_fwhm * fwhm) & good
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
    print(f"{'region':22s} {'pixels':>8s} {'sigma_z':>8s} {'coherence':>10s}")
    def cs(x):                       # clipped plain std: structure lives in the tails too
        x = x[np.isfinite(x)]
        x = x[np.abs(x) < 20]
        return float(np.std(x)) if x.size >= 20 else np.nan

    for name, m in (("empty sky", sky), ("bright-star disks", disks), ("other source pixels", other)):
        sz = rs(z[m])
        sp, sb = cs(z[m]), cs(zb[m])
        coh = (3 * sb / sp) ** 2 if np.isfinite(sp) and np.isfinite(sb) and sp > 0 else np.nan
        out[name] = dict(n=int(m.sum()), sigma_z=sz, coherence=coh)
        print(f"{name:22s} {int(m.sum()):8d} {sz:8.3f} {coh:10.2f}")
    c_sky = out["empty sky"]["coherence"]
    c_b, s_b = out["bright-star disks"]["coherence"], out["bright-star disks"]["sigma_z"]
    print("\npixel-to-pixel scatter vs model flux:")
    hf = L.highfreq_noise_check_arrays(resid, var, model, good, verbose=True)
    out["highfreq"] = hf
    # Decision: how much of the disks' excess variance (sigma_z^2 - 1) is explained by
    # pixel-to-pixel scatter above the variance plane at the brightest model fluxes?
    hf_bright = float(hf["ratio"][-1]) if len(hf["ratio"]) else np.nan
    if np.isfinite(s_b) and np.isfinite(hf_bright) and s_b > 1.02:
        share = max(0.0, hf_bright - 1.0) / max(s_b ** 2 - 1.0, 1e-9)
        out["noise_share"] = share
        if share >= 0.5:
            verdict = ("bright-star disks: excess is PIXEL-SCALE NOISE (pixel-to-pixel scatter "
                       "explains it) -- photon noise missing from the variance plane; use "
                       "poisson_k=k.")
        elif share >= 0.25:
            verdict = "bright-star disks: MIXED -- part missing photon noise, part model structure."
        else:
            verdict = ("bright-star disks: MODEL STRUCTURE -- the variance plane is right "
                       "(pixel-to-pixel scatter ~ variance), the excess is coherent residual "
                       "the model misses.  See show_bright_star_stack(run).")
        print(f"\nbright disks: sigma_z^2 = {s_b ** 2:.2f}; pixel-to-pixel scatter^2/variance at the "
              f"brightest fluxes = {hf_bright:.2f}  ->  noise explains {100 * share:.0f}% of the "
              f"excess  ->  {verdict}")
        out["verdict"] = verdict
    if plot:
        import matplotlib.pyplot as plt
        snrmap = ca._matched_snr(np.where(good, resid, 0.0), np.where(var > 0, var, 1.0), good,
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


def modeled_near(modeled, bad, fwhm):
    """Pixels closer to a modelled (near-saturated) core than to any other masked pixel, within
    reach -- residual blobs there are not grown into the grey area."""
    from scipy.ndimage import distance_transform_edt
    if not modeled.any():
        return np.zeros(bad.shape, bool)
    hard = bad & ~modeled
    dm = distance_transform_edt(~modeled)
    dh = distance_transform_edt(~hard) if hard.any() else np.full(bad.shape, np.inf)
    return dm < dh


def unmodeled_near_mask(resid, var, bad, fwhm, pad_fwhm=1.5, grow_snr=2.5, reach_fwhm=4.0,
                        modeled=None):
    """Pixels whose residual cannot be judged because a source there is (partly) under the
    bad-pixel mask: everything within pad_fwhm*FWHM of a bad pixel, plus any positive
    residual blob (matched-filter S/N > grow_snr) that touches that pad, limited to
    reach_fwhm*FWHM from the mask.
    modeled: masked pixels whose star IS in the model (near-saturated cores fitted from their
    wings) -- nothing is greyed around those: the residual next to them is real information
    (wing-fit quality, missed neighbours) and is shown and counted."""
    from scipy.ndimage import distance_transform_edt, label
    if not bad.any():
        return np.zeros(bad.shape, bool)
    modeled = (modeled & bad) if modeled is not None else np.zeros(bad.shape, bool)
    hard = bad & ~modeled
    dist = distance_transform_edt(~bad)
    pad = ~bad & (distance_transform_edt(~hard) <= pad_fwhm * fwhm) if hard.any() \
        else np.zeros(bad.shape, bool)

    usable = ~bad & np.isfinite(resid) & np.isfinite(var) & (var > 0)
    snr = ca._matched_snr(np.where(usable, resid, 0.0), np.where(usable, var, 1.0), usable, fwhm,
                          usable & ~pad)[0]
    lab, _ = label((snr > grow_snr) & usable & ~modeled_near(modeled, bad, fwhm)
                   & (dist <= reach_fwhm * fwhm))
    touch = np.unique(lab[pad & (lab > 0)])
    return pad | (np.isin(lab, touch[touch > 0]) & (lab > 0))


def run_session_gaia(session, config=None, names=None, out_dir=None, viz=False, plot_final=True,
                     gmag_max=None, n_jobs=1, plot_kwargs=None, **overrides):
    """Run every (or the named) cutout of a Session with Gaia seeds. Calls attach_gaia()
    first if session.gaia is missing. Returns {name: CutoutRun}.
    n_jobs > 1 runs cutouts in parallel worker processes (Linux fork; one core each).
    Results that cannot be sent back are still saved under out_dir (set it!).
    plot_kwargs: passed to plot_final_comparison, e.g. dict(circles="none")."""
    from . import session as cn
    import matplotlib.pyplot as plt
    from .plotting import plot_final_comparison
    if getattr(session, "gaia", None) is None:
        try:
            attach_gaia(session)
        except Exception as e:
            print(f"Gaia unavailable ({e}) -- running WITHOUT seeds (adaptive loops still on)")
            session.gaia = None
    runs = {}
    todo = names or list(session.jobs)
    if n_jobs and n_jobs > 1 and len(todo) > 1:
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor, as_completed
        if not out_dir:
            print("NOTE: set out_dir with n_jobs>1 so results are on disk even if they cannot "
                  "be sent back to the notebook")
        _PAR["args"] = (session, config, out_dir, gmag_max, overrides)
        t0 = time.perf_counter()
        with ProcessPoolExecutor(min(n_jobs, len(todo)), mp_context=mp.get_context("fork")) as ex:
            futs = {ex.submit(_par_worker, nm): nm for nm in todo}
            for f in as_completed(futs):
                nm = futs[f]
                try:
                    nm, run, summ, path = f.result()
                except Exception:
                    print(f"{nm} FAILED:")
                    traceback.print_exc()
                    continue
                print(f"{nm}: {summ['n_sources']} sources ({summ['n_from_gaia']} Gaia), residual "
                      f"{'clean' if summ['residual_clean'] else 'still shows excess'}, "
                      f"{summ['wall_seconds']:.0f}s")
                if run is not None:
                    runs[nm] = run
                    if plot_final:
                        plot_final_comparison(run, **(plot_kwargs or {}))
                        plt.show()
                elif plot_final and path:
                    cn.plot_final_comparison_from_fits(path)
        print(f"all {len(todo)} cutouts: {time.perf_counter() - t0:.0f}s wall")
        return runs
    for k, name in enumerate(todo, 1):
        job = session.jobs[name]
        print(f"{'#' * 78}\n[{k}/{len(todo)}] {name}  interior={job['interior']}")
        try:
            run = run_cutout_gaia(job["raw"], session.psf, session.bad_bitmask, session.gaia,
                                  config, interior_bbox=job["interior"], starCat=session.starCat,
                                  bgImage=session.bgImage, viz=viz, name=name, gmag_max=gmag_max,
                                  **overrides)
            runs[name] = run
            s = run.summary
            print(f"{name}: {s['n_sources']} sources ({s['n_from_gaia']} Gaia, "
                  f"{s['n_in_interior']} in interior), residual "
                  f"{'clean' if s['residual_clean'] else 'still shows excess'}, "
                  f"{s['wall_seconds']:.0f}s")
            if out_dir:
                cn.save_run(run, os.path.join(out_dir, f"visit{session.visit}_det{session.detector}",
                                              name))
            if plot_final:
                plot_final_comparison(run, **(plot_kwargs or {}))
                plt.show()
        except Exception:
            traceback.print_exc()
    return runs


def report(run):
    """What the adaptive loops did, pass by pass."""
    hist = getattr(run, "outer_history", run)
    for i, h in enumerate(hist, 1):
        ad = h["result"]["adaptive"]
        res = h["result"]
        print(f"outer {i}: {h['n_peaks']} sources ({int(res['is_seed'].sum())} Gaia), "
              f"{ad['n_inner']} inner it., PSF change into pass "
              f"{100 * h['psf_change']:.2f}%, {ad.get('seconds', np.nan):.0f}s, stop: {ad['stop']}")
        for r in ad["log"]:
            print(f"    it {r['it'] + 1:2d}  thr {r['threshold']:4.2f}σ  N={r['n']:5d}  "
                  f"+{r['added']:<4d} peaks +{r['pos']}/-{r['neg']} (excess {r['excess']:+d})  "
                  f"{r['action']}")
        if ad["final"] is not None:
            ca._print_diag(f"  [outer {i}]", ad["final"])
    last = hist[-1]["result"]["adaptive"]
    if last["missing"] and last["final"] is not None:
        print("top remaining positive residual peaks (x, y, S/N, px to nearest fitted):")
        for x, y, s, d in [c for c in last["final"]["candidates"] if c[2] >= last["target"]][:10]:
            print(f"    ({x:7.1f}, {y:7.1f})  S/N {s:5.1f}  {d:5.1f}")
