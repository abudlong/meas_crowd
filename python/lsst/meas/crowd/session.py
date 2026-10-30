# This file is part of meas_crowd.
#
# Developed for the LSST Data Management System.
# See the COPYRIGHT file at the top-level directory of this distribution.
#
# Vendored from the original notebook module cutout_notebook.py (behaviour unchanged; only the
# imports of the sibling modules were made package-relative).

"""Notebook-friendly front end for the pipeline (the original cutout_notebook.py).

Notebook-friendly front end for cutout_fit_lib.py (the run_cutouts.py / plot_from_fits.py
logic turned into plain functions; nothing here touches matplotlib's backend, so figures
show inline in Jupyter).

Quick start
-----------
    import cutout_notebook as cn

    cfg = cn.make_config("crowded_cfg_v5_testing.json",      # optional JSON / dict
                         n_outer=3, min_snr=4.0,              # any run_cutout_pipeline kwarg
                         shape_kwargs=dict(chi2_max=2.5))     # (see `deep` note below)
    cn.show_config(cfg)

    sess = cn.prepare_session(visit=2025062100433, detector=80,
                              cutouts=["3750,3750,250", (1000, 2000, 300, 300)],
                              pad=0, psf_source="refit")

    runs = cn.run_session(sess, cfg, viz=False, out_dir="results")   # dict name -> CutoutRun
    run = runs[next(iter(runs))]
    stats = cn.plot_final_comparison_robust(run)                      # new residual plot

    # later, without re-fitting (needs only numpy/matplotlib/astropy):
    cn.plot_final_comparison_from_fits("results/visit2025062100433_det80/<cutout name>")

Config merging
--------------
make_config() starts from DEFAULT_CONFIG (the run_cutouts.py defaults). A nested dict you
pass (e.g. shape_kwargs, stack_kwargs, inner_kwargs) REPLACES the whole default dict by
default -- exactly what run_cutouts.py --config does, so results match the CLI. Pass
deep=True to merge nested dicts key-by-key instead (then unspecified keys keep defaults).
"""
import copy
import json
import os
import time
import traceback
from dataclasses import dataclass, field

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy.table import Table
from astropy.visualization import simple_norm

from . import fitLib as L

# =====================================================================
# Config
# =====================================================================
DEFAULT_CONFIG = dict(
    n_outer_iterations=4,
    n_inner_iterations=10,
    plot_every=3,
    source_view="changed",
    max_pages_per_iter=2,
    sources_per_page=4,
    refine_centroids=True,
    centroid_min_amp=50.0,
    n_centroid_rounds=3,
    refit_psf="stack",
    stack_kwargs=dict(n_use=12, step=0.8, window=(0.85, 1.0), max_radius=20),
    psf_polish_cycles=3,
    psf_kernel_sizes=(41, 27, 21, 15, 11),
    merge_separation_sigma=1.5,
    sys_frac=0.03,
    poisson_k=0.0,
    min_snr=5.0,
    prune_sep_sigma=1.5,
    reject_fwhm=0.75,
    bright_peak_sigma=100.0,
    freeze_fwhm=1.0,
    near_bright_fwhm=2.0,
    near_bright_snr=10.0,
    edge_margin_fwhm=0.5,
    use_eliminate_weak=True,
    use_shape_cut=True,
    shape_kwargs=dict(min_snr_local=5.0, chi2_max=2.0, width_range=(0.75, 1.35),
                      dchi2_width=9.0, trust_snr=15.0, min_contrast=20.0,
                      weak_snr=10.0, dchi2_width_weak=6.5),
    shape_in_first_pass=False,
    psf_wing_size=81,
    inner_kwargs=dict(n_iterations=10, threshold_start=8.0, threshold_step=-1.5,
                      threshold_min=4.5, convergence_frac=0.0),
)


# The recommended settings for crowded fields (deeper detection, looser shape test, closer
# merging than DEFAULT_CONFIG); originally crowded_cfg_v5_testing.json, built in so no file is
# needed:
#   cn.make_config(CROWDED_FIELD_CONFIG, min_snr=4.0, save_to="my_run.json")
# With the default deep=False, stack_kwargs below REPLACES the default one (n_use, step,
# window fall back to the library defaults), exactly as the JSON did via the CLI.
CROWDED_FIELD_CONFIG = dict(
    min_snr=3.0, near_bright_snr=15.0, reject_fwhm=0.5, near_bright_fwhm=2.0,
    prune_sep_sigma=1.0, merge_separation_sigma=1.0, outer_convergence_frac=0,
    stack_kwargs=dict(max_radius=20),
    shape_kwargs=dict(min_snr_local=3.0, chi2_max=3.0, width_range=(0.7, 1.5),
                      dchi2_width=12.0, trust_snr=20.0, min_contrast=1.0,
                      weak_snr=8.0, dchi2_width_weak=12.0),
    inner_kwargs=dict(n_iterations=8, threshold_start=8.0, threshold_step=-1.0,
                      threshold_min=3.0, convergence_frac=0.0),
)

# the old name, kept so existing notebooks keep working
V5_TESTING_CONFIG = CROWDED_FIELD_CONFIG


def _merge(dst, src, deep):
    """Merge ``src`` into ``dst`` (key by key for nested dicts when ``deep``) and return ``dst``."""
    for k, v in src.items():
        if deep and isinstance(v, dict) and isinstance(dst.get(k), dict):
            _merge(dst[k], v, True)
        else:
            dst[k] = copy.deepcopy(v)
    return dst


def _normalize(cfg):
    """Turn the lists JSON gives back into the tuples the pipeline expects."""
    if isinstance(cfg.get("stack_kwargs"), dict) and "window" in cfg["stack_kwargs"]:
        cfg["stack_kwargs"]["window"] = tuple(cfg["stack_kwargs"]["window"])
    if isinstance(cfg.get("shape_kwargs"), dict) and "width_range" in cfg["shape_kwargs"]:
        cfg["shape_kwargs"]["width_range"] = tuple(cfg["shape_kwargs"]["width_range"])
    if "psf_kernel_sizes" in cfg and cfg["psf_kernel_sizes"] is not None:
        cfg["psf_kernel_sizes"] = tuple(cfg["psf_kernel_sizes"])
    return cfg


def save_config(cfg, path):
    """Write a config dict to a JSON file.

    Tuples become lists; `make_config` turns them back into tuples when the file is loaded.

    Parameters
    ----------
    cfg : `dict`
        Pipeline configuration.
    path : `str`
        Output file.

    Returns
    -------
    path : `str`
        The file written.
    """
    path = os.fspath(path)
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(cfg, fh, indent=2, default=_jsonable)
    print(f"saved config to {path}")
    return path


def make_config(base=None, deep=False, n_outer=None, n_inner=None, save_to=None,
                **overrides):
    """Build a pipeline configuration (keyword arguments of `fitLib.run_cutout_pipeline`).

    Parameters
    ----------
    base : `str` or `dict`, optional
        A JSON file or a dict, applied on top of ``DEFAULT_CONFIG``.
    deep : `bool`, optional
        Merge nested dicts key by key (default: replace them).
    n_outer : `int`, optional
        Shortcut for ``n_outer_iterations``.
    n_inner : `int`, optional
        Shortcut for ``n_inner_iterations`` and ``inner_kwargs['n_iterations']``.
    save_to : `str`, optional
        Also write the configuration to this JSON file.
    **overrides
        Any `fitLib.run_cutout_pipeline` keyword, applied last.

    Returns
    -------
    cfg : `dict`
        The configuration.
    """
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    base_path = base if isinstance(base, (str, os.PathLike)) else None
    if base is not None:
        if isinstance(base, (str, os.PathLike)):
            with open(base) as fh:
                base = json.load(fh)
        _merge(cfg, base, deep)
    _merge(cfg, overrides, deep)
    if n_outer:
        cfg["n_outer_iterations"] = n_outer
    if n_inner:
        cfg["n_inner_iterations"] = n_inner
        cfg["inner_kwargs"]["n_iterations"] = n_inner
    for k in ("x_range", "y_range", "viz", "save_to"):   # not pipeline kwargs
        cfg.pop(k, None)
    cfg = _normalize(cfg)
    if save_to:
        save_to = os.fspath(save_to)
        # a bare filename goes next to the base JSON (else the current working directory)
        if (not os.path.dirname(save_to) and base_path is not None):
            save_to = os.path.join(os.path.dirname(os.path.abspath(base_path)), save_to)
        save_config(cfg, os.path.abspath(save_to))
    return cfg


def check_config(cfg):
    """Check that every key of a configuration is a `fitLib.run_cutout_pipeline` option.

    Parameters
    ----------
    cfg : `dict`
        Pipeline configuration.

    Raises
    ------
    ValueError
        Raised if ``cfg`` has keys the pipeline does not accept.
    """
    import inspect
    ok = set(inspect.signature(L.run_cutout_pipeline).parameters)
    ok -= {"exposure_crop", "psf", "bad_bitmask", "starCat", "bgImage",
           "x_range", "y_range", "viz"}
    bad = sorted(set(cfg) - ok)
    if bad:
        raise ValueError(f"unknown run_cutout_pipeline option(s) in config: {bad}. "
                         f"Valid: {sorted(ok)}")


def show_config(cfg, diff_only=False):
    """Print a configuration.

    Parameters
    ----------
    cfg : `dict`
        Pipeline configuration.
    diff_only : `bool`, optional
        Only print keys that differ from ``DEFAULT_CONFIG``.
    """
    out = {k: v for k, v in cfg.items()
           if not diff_only or _normalize(copy.deepcopy(DEFAULT_CONFIG)).get(k) != v}
    print(json.dumps(out, indent=2, default=list))


# =====================================================================
# Cutout specs / session
# =====================================================================
def parse_cutout(spec):
    """Parse a cutout specification.

    Parameters
    ----------
    spec : `str`, `tuple` or `dict`
        'x0,y0,size', 'x0,y0,w,h', 'x0:x1,y0:y1', a tuple (x0, y0, size) or (x0, y0, w, h), or
        dict(x0=, y0=, w=, h=, name=); absolute detector pixels.

    Returns
    -------
    x0, y0 : `int`
        Corner of the box.
    w, h : `int`
        Size of the box.
    name : `str` or `None`
        Name of the cutout, if given.

    Raises
    ------
    ValueError
        Raised if ``spec`` cannot be parsed.
    """
    if isinstance(spec, dict):
        w = int(spec.get("w") or spec["size"])
        return int(spec["x0"]), int(spec["y0"]), w, int(spec.get("h") or w), spec.get("name")
    if isinstance(spec, (tuple, list)):
        v = [int(t) for t in spec]
    else:
        s = spec.replace(" ", "")
        if ":" in s:
            xs, ys = s.split(",")
            x0, x1 = (int(t) for t in xs.split(":"))
            y0, y1 = (int(t) for t in ys.split(":"))
            return x0, y0, x1 - x0, y1 - y0, None
        v = [int(float(t)) for t in s.split(",")]
    if len(v) == 3:
        return v[0], v[1], v[2], v[2], None
    if len(v) == 4:
        return v[0], v[1], v[2], v[3], None
    raise ValueError(f"bad cutout spec {spec!r}")


@dataclass
class Session:
    """Everything needed to fit cutouts of one visit and detector.

    Parameters
    ----------
    visit : `int`
        Visit id.
    detector : `int`
        Detector id.
    bad_bitmask : `int`
        Mask bits of pixels the fit ignores.
    psf : `lsst.afw.detection.Psf`
        Detector-level PSF.
    starCat : `lsst.afw.table.SourceCatalog` or `None`
        PSF stars of the detector-level PSF determination.
    bgImage : `lsst.afw.image.Image` or `None`
        Detector-level background.
    jobs : `dict` [`str`, `dict`]
        Per cutout: ``interior`` (`lsst.geom.Box2I`, the requested box) and ``raw`` (the padded
        cutout before background subtraction).
    wcs : `lsst.afw.geom.SkyWcs`, optional
        The visit image's WCS.
    psf_used_xy : `numpy.ndarray`, optional
        Positions (detector px, shape (N, 2)) of the stars the detector-level PSF was built from
        (``psf_source='refit'`` only).
    """

    visit: int
    detector: int
    bad_bitmask: int
    psf: object
    starCat: object
    bgImage: object
    jobs: dict = field(default_factory=dict)     # name -> dict(interior=Box2I, raw=Exposure)
    wcs: object = None
    psf_used_xy: object = None


def build_detector_psf(visit_image, psf_source="refit", return_used=False):
    """Set up the detector-level PSF.

    With ``psf_source='visit'`` the visit image's own PSF is used.  With 'refit' the background
    is subtracted, sources are detected and measured, stars selected and a PCA PSF determined
    (the visit image is modified in place).

    Parameters
    ----------
    visit_image : `lsst.afw.image.Exposure`
        The detector's visit image.
    psf_source : `str`, optional
        'refit' or 'visit'.
    return_used : `bool`, optional
        Also return the positions of the stars the PCA PSF was built from.

    Returns
    -------
    psf : `lsst.afw.detection.Psf`
        Detector-level PSF.
    starCat : `lsst.afw.table.SourceCatalog` or `None`
        Selected PSF stars (None for 'visit').
    bgImage : `lsst.afw.image.Image` or `None`
        Subtracted background (None for 'visit').
    used_xy : `numpy.ndarray` or `None`
        Positions of the stars used, shape (N, 2) (only with ``return_used``; None for 'visit').
    """
    if psf_source == "visit":
        print("PSF source: the visit_image's own PSF (no detector-level refit)")
        return (visit_image.getPsf(), None, None) + ((None,) if return_used else ())

    print("Detector-level setup: background subtraction -> detection -> PSF determination")
    t0 = time.perf_counter()
    bg_task = L.SubtractBackgroundTask()
    bg_task.config.algorithm = "AKIMA_SPLINE"
    bg_task.config.useApprox = False
    bgImage = bg_task.run(visit_image).background.getImage()

    L.InstallGaussianPsfTask().run(visit_image)
    schema = L.afwTable.SourceTable.makeMinimalSchema()
    det_cfg = L.SourceDetectionTask.ConfigClass()
    det_cfg.thresholdValue = 3.0
    det_cfg.thresholdType = "stdev"
    det_cfg.minPixels = 2
    det_cfg.doTempLocalBackground = True
    det_cfg.tempLocalBackground.binSize = 16
    det_cfg.background.binSize = 32
    detection = L.SourceDetectionTask(schema=schema, config=det_cfg)
    deblend = L.SourceDeblendTask(schema=schema)
    measure = L.SingleFrameMeasurementTask(schema=schema)
    table = L.afwTable.SourceTable.make(schema)
    sources = detection.run(table, visit_image).sources
    deblend.run(visit_image, sources)
    measure.run(measCat=sources, exposure=visit_image)
    sources = sources.copy(deep=True)
    print(f"  {len(sources)} detections at 3 sigma")

    starCat = L.ObjectSizeStarSelectorTask().run(sourceCat=sources,
                                                 exposure=visit_image).sourceCat
    print(f"  {len(starCat)} PSF-star candidates")
    cands = L.MakePsfCandidatesTask().run(starCat, exposure=visit_image).psfCandidates
    cfg = L.PcaPsfDeterminerTask.ConfigClass()
    cfg.spatialOrder = 1
    cfg.sizeCellX = cfg.sizeCellY = 512
    cfg.maxCandidates = 1500
    cfg.doRejectBlends = True
    cfg.nStarPerCell = 10
    psf, cellSet = L.PcaPsfDeterminerTask(config=cfg).determinePsf(
        exposure=visit_image, psfCandidateList=cands)
    visit_image.setPsf(psf)
    used = [c for cell in cellSet.getCellList() for c in cell if c.getStatus() == L.GOOD_STATUS]
    used_xy = np.array([(c.getXCenter(), c.getYCenter()) for c in used], float).reshape(-1, 2)
    print(f"  PSF determined from {len(used)} stars  ({time.perf_counter() - t0:.1f}s)")
    return (psf, starCat, bgImage) + ((used_xy,) if return_used else ())


def prepare_session(visit, detector, cutouts, pad=0, psf_source="refit",
                    instrument="LSSTCam", repo="dp2_prep",
                    collection="LSSTCam/runs/DRP/DP2",
                    visit_image_collections="LSSTCam/runs/*", visit_image=None):
    """Load a visit image, cut out the regions and set up the detector-level PSF.

    The cutouts are taken before the detector-wide background subtraction, as in the original
    notebook.

    Parameters
    ----------
    visit : `int`
        Visit id.
    detector : `int`
        Detector id.
    cutouts : `list`
        Cutout specifications (see `parse_cutout`); absolute detector pixels.
    pad : `int`, optional
        Fit a region grown by this many px; the interior keeps the requested box.
    psf_source : `str`, optional
        'refit' or 'visit' (see `build_detector_psf`).
    instrument : `str`, optional
        Instrument name.
    repo : `str`, optional
        Butler repository path or alias.
    collection : `str`, optional
        Default Butler collection.
    visit_image_collections : `str`, optional
        Collections searched for the visit image.
    visit_image : `lsst.afw.image.Exposure`, optional
        An already-loaded visit image (skips the Butler).

    Returns
    -------
    session : `Session`
        The session.

    Raises
    ------
    RuntimeError
        Raised if no visit image is found.
    """
    if visit_image is None:
        from lsst.daf.butler import Butler
        butler = Butler(repo, collections=collection)
        refs = butler.query_datasets(
            "visit_image", collections=visit_image_collections,
            where=f"instrument='{instrument}' AND visit={visit} AND detector={detector}",
            find_first=False)
        if not refs:
            raise RuntimeError("no visit_image found for that visit/detector")
        visit_image = butler.get("visit_image", collections=visit_image_collections,
                                 dataId=refs[0].dataId)
    det_bbox = visit_image.getBBox()
    print(f"Loaded visit_image: bbox {det_bbox}")
    mk = visit_image.mask
    bad_bitmask = (mk.getPlaneBitMask("SAT") | mk.getPlaneBitMask("CR") |
                   mk.getPlaneBitMask("BAD") | mk.getPlaneBitMask("INTRP") |
                   mk.getPlaneBitMask("EDGE"))

    G = L.geom
    jobs = {}
    for spec in cutouts:
        x0, y0, w, h, nm = parse_cutout(spec)
        req = G.Box2I(G.Point2I(x0, y0), G.Extent2I(w, h))
        req.clip(det_bbox)
        if req.isEmpty():
            print(f"!! cutout {x0},{y0},{w}x{h} lies outside the detector -- skipped")
            continue
        fit_bbox = G.Box2I(req)
        fit_bbox.grow(pad)
        fit_bbox.clip(det_bbox)
        name = nm or f"cutout_x{req.getMinX()}_y{req.getMinY()}_{req.getWidth()}x{req.getHeight()}"
        jobs[name] = dict(interior=req, raw=visit_image[fit_bbox].clone())

    psf, starCat, bgImage, used_xy = build_detector_psf(visit_image, psf_source, return_used=True)
    return Session(visit, detector, bad_bitmask, psf, starCat, bgImage, jobs,
                   wcs=visit_image.getWcs(), psf_used_xy=used_xy)


# =====================================================================
# Gaia (optional)
# =====================================================================
def query_gaia(wcs, bbox, tap_url="https://gea.esac.esa.int/tap-server/tap"):
    """Query Gaia DR3 at the ESA archive for the sources on a box.

    Parameters
    ----------
    wcs : `lsst.afw.geom.SkyWcs`
        WCS of the image.
    bbox : `lsst.geom.Box2I`
        Pixel box.
    tap_url : `str`, optional
        TAP service URL.

    Returns
    -------
    gaia : `astropy.table.Table`
        Gaia sources inside the box, with their pixel positions ``x``, ``y``.
    """
    import pyvo
    corners = [wcs.pixelToSky(x, y) for x, y in
               [(bbox.getMinX(), bbox.getMinY()), (bbox.getMinX(), bbox.getMaxY()),
                (bbox.getMaxX(), bbox.getMaxY()), (bbox.getMaxX(), bbox.getMinY())]]
    center = wcs.pixelToSky(bbox.getCenter())
    radius = max(center.separation(c).asDegrees() for c in corners) * 1.05
    ra, dec = center.getRa().asDegrees(), center.getDec().asDegrees()
    q = f"""SELECT source_id, ra, dec, parallax, pmra, pmdec,
                   phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag
            FROM gaiadr3.gaia_source
            WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, {radius}))"""
    gaia = pyvo.dal.TAPService(tap_url).search(q).to_table()
    x, y = wcs.skyToPixelArray(np.asarray(gaia["ra"], float), np.asarray(gaia["dec"], float),
                               degrees=True)
    inside = ((x >= bbox.getMinX()) & (x <= bbox.getMaxX()) &
              (y >= bbox.getMinY()) & (y <= bbox.getMaxY()))
    gaia = gaia[inside].copy()
    gaia["x"], gaia["y"] = x[inside], y[inside]
    return gaia


# =====================================================================
# Running
# =====================================================================
@dataclass
class CutoutRun:
    """Result of running the pipeline on one cutout.

    Parameters
    ----------
    name : `str`
        Cutout name.
    outer_history : `list` [`dict`]
        One dict per outer pass (``psf``, ``exposure``, ``result``, ...).
    catalog : `astropy.table.Table`
        Source catalogue (see `build_catalog`).
    summary : `dict`
        Run summary (timings, counts, residual statistics).
    raw_crop : `lsst.afw.image.Exposure`
        The input cutout before background subtraction.
    interior_bbox : `lsst.geom.Box2I`
        The requested box (without padding).
    bad_bitmask : `int`
        Mask bits of pixels the fit ignored.
    config : `dict`
        Configuration used.
    """

    name: str
    outer_history: list
    catalog: Table
    summary: dict
    raw_crop: object
    interior_bbox: object
    bad_bitmask: int
    config: dict


def build_catalog(outer_history, interior_bbox, gaia_tbl=None, match_radius=2.0):
    """Build the source catalogue of a finished run.

    Parameters
    ----------
    outer_history : `list` [`dict`]
        Outer-loop history.
    interior_bbox : `lsst.geom.Box2I`
        The requested box (``in_interior`` column).
    gaia_tbl : `astropy.table.Table`, optional
        Gaia sources (with ``x``, ``y``) to match.
    match_radius : `float`, optional
        Gaia match radius [px].

    Returns
    -------
    cat : `astropy.table.Table`
        Positions (pixel and sky), PSF flux, AB magnitude (assuming nJy), source metrics,
        ``in_interior`` and, with ``gaia_tbl``, the nearest Gaia source.
    """
    final = outer_history[-1]
    exp, res = final["exposure"], final["result"]
    peaks = np.asarray(res["peaks"], float).reshape(-1, 2)
    n = len(peaks)
    cat = Table()
    cat["id"] = np.arange(n)
    cat["x"], cat["y"] = peaks[:, 0], peaks[:, 1]
    xy0 = exp.getXY0()
    cat["x_cutout"], cat["y_cutout"] = peaks[:, 0] - xy0.getX(), peaks[:, 1] - xy0.getY()
    flux = np.asarray(res["amps"], float) if n else np.zeros(0)
    cat["psfFlux_nJy"] = flux
    with np.errstate(divide="ignore", invalid="ignore"):
        cat["mag_AB"] = np.where(flux > 0, -2.5 * np.log10(flux) + 31.4, np.nan)
    try:
        ra, dec = exp.getWcs().pixelToSkyArray(peaks[:, 0], peaks[:, 1], degrees=True)
        cat["ra"], cat["dec"] = np.asarray(ra), np.asarray(dec)
    except Exception as e:
        print(f"  (no sky coords in catalog: {e})")
    if n and res.get("snaps"):
        try:
            m = L.source_metrics(outer_history)
            cat["snr"], cat["chi2_core"] = m["snr"], m["chi2"]
            cat["nn_pix"], cat["edge_pix"], cat["group_size"] = m["nn"], m["edge"], m["group"]
        except Exception as e:
            print(f"  (source_metrics failed: {e})")
    ib = interior_bbox
    cat["in_interior"] = ((peaks[:, 0] >= ib.getMinX() - .5) & (peaks[:, 0] <= ib.getMaxX() + .5) &
                          (peaks[:, 1] >= ib.getMinY() - .5) & (peaks[:, 1] <= ib.getMaxY() + .5))
    if gaia_tbl is not None and n:
        from scipy.spatial import cKDTree
        d, i = cKDTree(np.column_stack([gaia_tbl["x"], gaia_tbl["y"]])).query(peaks)
        cat["gaia_source_id"] = np.asarray(gaia_tbl["source_id"])[i]
        cat["gaia_G"] = np.asarray(gaia_tbl["phot_g_mean_mag"], float)[i]
        cat["gaia_sep_pix"] = d
        cat["gaia_match"] = d <= match_radius
    return cat


def run_cutout(raw_crop, psf, bad_bitmask, config=None, interior_bbox=None, starCat=None,
               bgImage=None, viz=False, diagnostics=True, gaia_tbl=None, match_radius=2.0,
               name="cutout", **overrides):
    """Run the pipeline on one cropped exposure.

    Parameters
    ----------
    raw_crop : `lsst.afw.image.Exposure`
        The cutout (before background subtraction).
    psf : `lsst.afw.detection.Psf`
        Initial PSF.
    bad_bitmask : `int`
        Mask bits of pixels to ignore.
    config : `dict`, optional
        Pipeline configuration (default ``DEFAULT_CONFIG``).
    interior_bbox : `lsst.geom.Box2I`, optional
        The requested box (default: the whole cutout).
    starCat : `lsst.afw.table.SourceCatalog`, optional
        PSF stars forced into every PSF determination.
    bgImage : `lsst.afw.image.Image`, optional
        Detector background subtracted before the first pass.
    viz : `bool`, optional
        Draw the per-iteration diagnostics.
    diagnostics : `bool`, optional
        Run the end-of-run diagnostics.
    gaia_tbl : `astropy.table.Table`, optional
        Gaia sources to match in the catalogue.
    match_radius : `float`, optional
        Gaia match radius [px].
    name : `str`, optional
        Cutout name.
    **overrides
        Configuration overrides.

    Returns
    -------
    run : `CutoutRun`
        The result.
    """
    cfg = make_config(config, **overrides)   # idempotent for an already-complete config
    check_config(cfg)
    cb = raw_crop.getBBox()
    interior_bbox = interior_bbox or cb
    t0 = time.perf_counter()
    hist = L.run_cutout_pipeline(raw_crop, psf, bad_bitmask, starCat=starCat, bgImage=bgImage,
                                 x_range=(0, cb.getWidth()), y_range=(0, cb.getHeight()),
                                 viz=viz, **cfg)
    res = hist[-1]["result"]
    summary = dict(name=name, status="ok",
                   bbox_fitted=[cb.getMinX(), cb.getMaxX(), cb.getMinY(), cb.getMaxY()],
                   interior_bbox=[interior_bbox.getMinX(), interior_bbox.getMaxX(),
                                  interior_bbox.getMinY(), interior_bbox.getMaxY()],
                   n_sources=len(res["peaks"]),
                   resid_std_per_outer=[h["resid_std"] for h in hist],
                   n_peaks_per_outer=[h["n_peaks"] for h in hist],
                   psf_fwhm_px_per_outer=[L._fwhm(h["psf"], h["exposure"].getBBox().getCenter())
                                          for h in hist],
                   pipeline_kwargs=cfg)
    if diagnostics and len(res["peaks"]) and res.get("snaps"):
        try:
            summary["highfreq_noise_k"] = L.highfreq_noise_check(hist)["k"]
        except Exception as e:
            print(f"  highfreq_noise_check failed: {e}")
        try:
            snr_map = L.negative_source_check(hist, plot=viz)
            summary["negative_check"] = {str(t): dict(pos=int((snr_map >= t).sum()),
                                                      neg=int((snr_map <= -t).sum()))
                                         for t in (3, 4, 5)}
        except Exception as e:
            print(f"  negative_source_check failed: {e}")
    cat = build_catalog(hist, interior_bbox, gaia_tbl, match_radius)
    summary["n_in_interior"] = int(np.sum(cat["in_interior"]))
    summary["wall_seconds"] = time.perf_counter() - t0
    return CutoutRun(name, hist, cat, summary, raw_crop, interior_bbox, bad_bitmask, cfg)


def run_session(session, config=None, names=None, out_dir=None, viz=False, gaia=False,
                match_radius=2.0, plot_final=True, **overrides):
    """Run the pipeline on every (or the named) cutout of a session.

    A failure in one cutout is printed and skipped.

    Parameters
    ----------
    session : `Session`
        The session.
    config : `dict`, optional
        Pipeline configuration.
    names : `list` [`str`], optional
        Cutouts to run (default: all).
    out_dir : `str`, optional
        Save each run under ``<out_dir>/visit<V>_det<D>/<name>/`` (see `save_run`).
    viz : `bool`, optional
        Draw the per-iteration diagnostics.
    gaia : `bool`, optional
        Match the catalogues to Gaia (queried at the ESA archive).
    match_radius : `float`, optional
        Gaia match radius [px].
    plot_final : `bool`, optional
        Plot the final comparison of each cutout.
    **overrides
        Configuration overrides.

    Returns
    -------
    runs : `dict` [`str`, `CutoutRun`]
        The results by cutout name.
    """
    cfg = make_config(config, **overrides)
    runs = {}
    todo = names or list(session.jobs)
    for k, name in enumerate(todo, 1):
        job = session.jobs[name]
        print(f"{'#' * 78}\n[{k}/{len(todo)}] {name}  interior={job['interior']}")
        try:
            gaia_tbl = None
            if gaia:
                try:
                    gaia_tbl = query_gaia(session.wcs, job["raw"].getBBox())
                except Exception as e:
                    print(f"Gaia query failed ({type(e).__name__}: {e}) -- continuing without")
            run = run_cutout(job["raw"], session.psf, session.bad_bitmask, cfg,
                             interior_bbox=job["interior"], starCat=session.starCat,
                             bgImage=session.bgImage, viz=viz, gaia_tbl=gaia_tbl,
                             match_radius=match_radius, name=name)
            runs[name] = run
            print(f"{name}: {run.summary['n_sources']} sources "
                  f"({run.summary['n_in_interior']} in interior), "
                  f"{run.summary['wall_seconds']:.0f}s")
            if out_dir:
                save_run(run, os.path.join(out_dir,
                                           f"visit{session.visit}_det{session.detector}", name))
            if plot_final:
                plot_final_comparison_robust(run)
                plt.show()
        except Exception:
            traceback.print_exc()
    return runs


# =====================================================================
# Saving
# =====================================================================
def _jsonable(o):
    """Convert numpy scalars / arrays and other objects into JSON-serialisable values."""
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (set, tuple)):
        return list(o)
    return str(o)


def _wcs_header(exposure):
    """Return the exposure's WCS as FITS header cards (empty if it has none)."""
    hdr = {}
    try:
        meta = exposure.getWcs().getFitsMetadata()
        for k in meta.names():
            hdr[k] = meta.get(k)
        xy0 = exposure.getXY0()
        if "CRPIX1" in hdr:
            hdr["CRPIX1"] = float(hdr["CRPIX1"]) - xy0.getX()
            hdr["CRPIX2"] = float(hdr["CRPIX2"]) - xy0.getY()
    except Exception as e:
        print(f"  (no FITS WCS written: {e})")
    return hdr


def save_run(run, outdir):
    """Write a run's catalogue, images, PSF kernels and summary to a directory.

    Writes catalog.csv / catalog.fits, images.fits, psf_kernels.fits and summary.json, the layout
    `plot_final_comparison_from_fits` reads.

    Parameters
    ----------
    run : `CutoutRun`
        The run.
    outdir : `str`
        Output directory (created if needed).
    """
    from astropy.io import fits
    os.makedirs(outdir, exist_ok=True)
    final = run.outer_history[-1]
    exp, res = final["exposure"], final["result"]
    xy0 = exp.getXY0()
    hdr = fits.Header()
    for k, v in _wcs_header(exp).items():
        try:
            hdr[k] = v
        except Exception:
            pass
    hdr["X0"], hdr["Y0"] = xy0.getX(), xy0.getY()
    model, data = res["model_full"], exp.image.array
    hdus = [fits.PrimaryHDU(header=hdr)]
    for nm, arr in (("DATA_USED", data), ("MODEL", model), ("RESIDUAL", data - model),
                    ("VARIANCE", exp.variance.array), ("RAW_CUTOUT", run.raw_crop.image.array)):
        hdus.append(fits.ImageHDU(np.asarray(arr, dtype=np.float32), header=hdr, name=nm))
    hdus.append(fits.ImageHDU(np.asarray(exp.mask.array, dtype=np.int32), header=hdr, name="MASK"))
    hdus.append(fits.ImageHDU(((exp.mask.array & run.bad_bitmask) != 0).astype(np.uint8),
                              header=hdr, name="BAD_PIXELS"))
    fits.HDUList(hdus).writeto(os.path.join(outdir, "images.fits"), overwrite=True)

    k_hdus = [fits.PrimaryHDU()]
    for i, h in enumerate(run.outer_history):
        k = h["psf"].computeKernelImage(h["exposure"].getBBox().getCenter()).array
        k_hdus.append(fits.ImageHDU(np.asarray(k, dtype=np.float32), name=f"PSF_OUTER{i + 1}"))
    fits.HDUList(k_hdus).writeto(os.path.join(outdir, "psf_kernels.fits"), overwrite=True)

    run.catalog.write(os.path.join(outdir, "catalog.csv"), format="ascii.csv", overwrite=True)
    run.catalog.write(os.path.join(outdir, "catalog.fits"), overwrite=True)
    with open(os.path.join(outdir, "summary.json"), "w") as fh:
        json.dump(run.summary, fh, indent=2, default=_jsonable)
    with open(os.path.join(outdir, "config.json"), "w") as fh:     # reload with make_config()
        json.dump(run.config, fh, indent=2, default=_jsonable)
    print(f"saved run to {outdir}")


# =====================================================================
# Final comparison plot with robust residual statistics
# =====================================================================
EDGE_COLOR = "#ff00ff"
IQR_TO_SIGMA = 1.349          # sigma_robust = (Q75 - Q25) / 1.349 for a Gaussian


def robust_stats(values, valid=None):
    """Return the median and IQR-based robust sigma of finite (and valid) values.

    Parameters
    ----------
    values : `numpy.ndarray`
        Values.
    valid : `numpy.ndarray` [`bool`], optional
        Values to use.

    Returns
    -------
    stats : `dict`
        ``median``, ``sigma_iqr`` ((Q75 - Q25) / 1.349), ``q25``, ``q75`` and ``n``.
    """
    v = np.asarray(values, float)
    m = np.isfinite(v) if valid is None else (np.asarray(valid, bool) & np.isfinite(v))
    x = v[m]
    if x.size == 0:
        return dict(median=np.nan, sigma_iqr=np.nan, q25=np.nan, q75=np.nan, n=0)
    q25, med, q75 = np.percentile(x, [25, 50, 75])
    return dict(median=float(med), sigma_iqr=float((q75 - q25) / IQR_TO_SIGMA),
                q25=float(q25), q75=float(q75), n=int(x.size))


def _draw_edge_margin(ax, ext, margin_px, shade=True):
    """Outline and shade the band within ``margin_px`` of the image edge."""
    if not margin_px or margin_px <= 0:
        return
    xa, xb, ya, yb = ext
    if shade:
        for (a, b, c, d) in ((xa, xb, ya, ya + margin_px), (xa, xb, yb - margin_px, yb),
                             (xa, xa + margin_px, ya, yb), (xb - margin_px, xb, ya, yb)):
            ax.add_patch(patches.Rectangle((a, c), b - a, d - c, fc=EDGE_COLOR, ec="none",
                                           alpha=0.18, zorder=3))
    ax.add_patch(patches.Rectangle((xa + margin_px, ya + margin_px),
                                   (xb - xa) - 2 * margin_px, (yb - ya) - 2 * margin_px,
                                   fill=False, ec=EDGE_COLOR, lw=1.6, ls="--", zorder=4))


def _draw_final_comparison(raw, data, model, var, bad, x0, y0, peaks, radius, margin_px,
                           view=None, stats_in_view=True, show_peaks=True,
                           figsize=(26, 7), resid_clip_pct=99.0, z_clip=5.0):
    """Draw original | model | residual | residual/sigma from plain arrays; return (fig, stats).

    ``bad`` marks pixels that must not be used: they are drawn black in the residual panels and
    excluded from the statistics.
    """
    ny, nx = data.shape
    ext = (x0 - .5, x0 + nx - .5, y0 - .5, y0 + ny - .5)
    resid = data - model
    varok = np.isfinite(var) & (var > 0)
    sigma = np.sqrt(np.where(varok, var, np.nan))
    z = resid / sigma

    use = ~bad & np.isfinite(resid)                     # pixels allowed in the residual
    use_z = use & varok & np.isfinite(z)
    region = np.ones(data.shape, bool)
    if view is not None and stats_in_view:
        X = x0 + np.arange(nx)[None, :]
        Y = y0 + np.arange(ny)[:, None]
        region = (X >= view[0]) & (X <= view[1]) & (Y >= view[2]) & (Y <= view[3])
    st_r = robust_stats(resid, use & region)
    st_z = robust_stats(z, use_z & region)

    cm_res = plt.get_cmap("RdBu_r").copy()
    cm_res.set_bad("black")
    res_m = np.ma.masked_where(~use, resid)
    z_m = np.ma.masked_where(~use_z, z)
    sel = resid[use & region]
    rs = np.nanpercentile(np.abs(sel), resid_clip_pct) if sel.size else 1.0

    norm = simple_norm(raw, vmin=np.nanmin(raw), vmax=np.nanpercentile(raw, 99),
                       stretch="log", log_a=1000)
    cm_img = L.CMAP

    def sub(st, unit=""):
        return f"median {st['median']:+.2f}{unit},  robust σ {st['sigma_iqr']:.2f}{unit}  (N={st['n']})"

    panels = [
        ("Original visit image (cutout)", raw, dict(cmap=cm_img, norm=norm)),
        (f"Final model ({len(peaks)} sources)", model, dict(cmap=cm_img, norm=norm)),
        ("Residual = data − model\n" + sub(st_r), res_m, dict(cmap=cm_res, vmin=-rs, vmax=rs)),
        ("Residual / σ\n" + sub(st_z), z_m, dict(cmap=cm_res, vmin=-z_clip, vmax=z_clip)),
    ]
    fig, axes = plt.subplots(1, 4, figsize=figsize, sharex=True, sharey=True)
    for i, (ax, (title, img, kw)) in enumerate(zip(axes, panels)):
        im = ax.imshow(img, extent=ext, origin="lower", **kw)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        if i < 2 and bad.any():                              # red overlay on data / model
            ax.imshow(np.ma.masked_where(~bad, np.ones(bad.shape)), extent=ext,
                      origin="lower", cmap=matplotlib_listed("red"), alpha=0.35,
                      interpolation="nearest")
        _draw_edge_margin(ax, ext, margin_px)
        if show_peaks and len(peaks):
            for px, py in peaks:
                ax.add_patch(patches.Circle((px, py), radius, fill=False, ec="w",
                                            lw=0.5, alpha=0.6))
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("x [px]")
        ax.grid(False)
        if i >= 2:                                           # stats box on the residual panels
            st = st_r if i == 2 else st_z
            ax.text(0.02, 0.98, f"median = {st['median']:+.3f}\nσ_IQR = {st['sigma_iqr']:.3f}",
                    transform=ax.transAxes, va="top", ha="left", fontsize=10,
                    bbox=dict(fc="white", ec="k", alpha=0.85, pad=3), zorder=10)
    axes[0].set_ylabel("y [px]")
    if view is not None:
        axes[0].set_xlim(view[0], view[1])
        axes[0].set_ylim(view[2], view[3])
    else:
        axes[0].set_xlim(ext[0], ext[1])
        axes[0].set_ylim(ext[2], ext[3])
    handles = [patches.Patch(fc="red", alpha=0.35, label="masked (bad_bitmask)")]
    if margin_px:
        handles.insert(0, patches.Patch(fc=EDGE_COLOR, alpha=0.3, ec=EDGE_COLOR, ls="--",
                                        label=f"edge drop zone ({margin_px:.1f}px)"))
    axes[0].legend(handles=handles, loc="upper left", fontsize=8)
    fig.text(0.5, -0.01,
             "Residual stats: bad/non-finite pixels excluded (black in residual panels); "
             "median and robust σ = (Q75−Q25)/1.349"
             + (" over the plotted view" if (view is not None and stats_in_view) else
                " over the full array"),
             ha="center", fontsize=9)
    plt.tight_layout()
    return fig, dict(residual=st_r, residual_over_sigma=st_z)


def matplotlib_listed(color):
    """Return a one-colour matplotlib colour map.

    Parameters
    ----------
    color : `str`
        The colour.

    Returns
    -------
    cmap : `matplotlib.colors.ListedColormap`
        The colour map.
    """
    from matplotlib.colors import ListedColormap
    return ListedColormap([color])


def _print_stats(stats):
    """Print residual statistics."""
    for key, st in stats.items():
        print(f"{key:20s} median={st['median']:+9.3f}  robust σ(IQR/1.349)={st['sigma_iqr']:8.3f}  "
              f"Q25={st['q25']:+.3f}  Q75={st['q75']:+.3f}  N={st['n']}")


def plot_final_comparison_robust(run_or_history, exposure_crop=None, bad_bitmask=None,
                                 edge_margin_fwhm=None, view="interior", stats_in_view=True,
                                 show_peaks=True, figsize=(26, 7), verbose=True):
    """Plot original | final model | residual | residual/sigma with robust statistics.

    Bad pixels are black in the residual panels; median and IQR-robust sigma (good pixels only)
    are shown in the titles.

    Parameters
    ----------
    run_or_history : `CutoutRun` or `list` [`dict`]
        A run (preferred) or its outer-loop history.
    exposure_crop : `lsst.afw.image.Exposure`, optional
        Raw crop (only needed with a bare history).
    bad_bitmask : `int`, optional
        Mask bits of bad pixels (needed with a bare history).
    edge_margin_fwhm : `float`, optional
        Draw the edge drop-zone of this width [FWHM].
    view : `str` or `tuple`, optional
        'interior' (the run's interior box), None (the whole fitted array) or
        (xmin, xmax, ymin, ymax) in detector pixels.
    stats_in_view : `bool`, optional
        Restrict the statistics to the displayed view.
    show_peaks : `bool`, optional
        Mark the sources.
    figsize : `tuple` [`float`, `float`], optional
        Figure size [inches].
    verbose : `bool`, optional
        Print the statistics.

    Returns
    -------
    stats : `dict`
        ``residual`` and ``residual_over_sigma``: median, sigma_iqr, q25, q75, n.

    Raises
    ------
    ValueError
        Raised if a bare history is passed without ``bad_bitmask``.
    """
    if isinstance(run_or_history, CutoutRun):
        run = run_or_history
        hist, raw_exp = run.outer_history, run.raw_crop
        bad_bitmask = run.bad_bitmask if bad_bitmask is None else bad_bitmask
        if edge_margin_fwhm is None:
            edge_margin_fwhm = run.config.get("edge_margin_fwhm", 0.5)
        if isinstance(view, str) and view == "interior":
            b = run.interior_bbox
            view = (b.getMinX() - .5, b.getMaxX() + .5, b.getMinY() - .5, b.getMaxY() + .5)
    else:
        hist, raw_exp = run_or_history, exposure_crop
        if bad_bitmask is None:
            raise ValueError("bad_bitmask is required when passing a bare outer_history")
        if edge_margin_fwhm is None:
            edge_margin_fwhm = 0.5
        if isinstance(view, str):
            view = None
    final = hist[-1]
    exp, res = final["exposure"], final["result"]
    bbox = exp.getBBox()
    model = res["model_full"]
    data = exp.image.array
    raw = raw_exp[bbox].image.array if raw_exp is not None else data
    bad = (exp.mask.array & bad_bitmask) != 0
    peaks = np.asarray(res["peaks"], float).reshape(-1, 2)
    fwhm = L._fwhm(exp.getPsf(), bbox.getCenter())
    fig, stats = _draw_final_comparison(
        raw, data, model, exp.variance.array, bad, bbox.getMinX(), bbox.getMinY(), peaks,
        0.75 * fwhm, edge_margin_fwhm * fwhm if edge_margin_fwhm else 0.0,
        view=view, stats_in_view=stats_in_view, show_peaks=show_peaks, figsize=figsize)
    if verbose:
        _print_stats(stats)
    return stats


def plot_final_comparison_from_fits(d, view=None, out=None, stats_in_view=True,
                                    show_sources=True, show=True, dpi=110, verbose=True):
    """Rebuild the final-comparison figure from files written by `save_run`.

    No LSST stack and no re-fit needed.

    Parameters
    ----------
    d : `str`
        Directory written by `save_run`.
    view : `tuple` [`float`, `float`, `float`, `float`], optional
        (xmin, xmax, ymin, ymax) in detector pixels (default: the interior box from
        summary.json, else the whole array).
    out : `str`, optional
        Also save the figure to this PNG.
    stats_in_view : `bool`, optional
        Restrict the statistics to the view.
    show_sources : `bool`, optional
        Mark the catalogue's sources.
    show : `bool`, optional
        Display the figure.
    dpi : `int`, optional
        Resolution of the PNG.
    verbose : `bool`, optional
        Print the statistics.

    Returns
    -------
    stats : `dict`
        See `plot_final_comparison_robust`.
    """
    from astropy.io import fits
    with fits.open(os.path.join(d, "images.fits")) as hl:
        names = [h.name for h in hl]
        data = hl["DATA_USED"].data.astype(float)
        model = hl["MODEL"].data.astype(float)
        var = hl["VARIANCE"].data.astype(float) if "VARIANCE" in names else np.ones_like(data)
        raw = hl["RAW_CUTOUT"].data.astype(float) if "RAW_CUTOUT" in names else data
        bad = (hl["BAD_PIXELS"].data.astype(bool) if "BAD_PIXELS" in names
               else np.zeros(data.shape, bool))
        x0, y0 = hl["DATA_USED"].header.get("X0", 0), hl["DATA_USED"].header.get("Y0", 0)
    if raw.shape != data.shape:
        raw = data
    summ = {}
    sp = os.path.join(d, "summary.json")
    if os.path.exists(sp):
        with open(sp) as fh:
            summ = json.load(fh)
    fwhm = (summ.get("psf_fwhm_px_per_outer") or [3.5])[-1]
    edge_margin = (summ.get("pipeline_kwargs") or {}).get("edge_margin_fwhm", 0.5)
    if view is None and "interior_bbox" in summ:
        b = summ["interior_bbox"]
        view = (b[0] - .5, b[1] + .5, b[2] - .5, b[3] + .5)
    peaks = np.zeros((0, 2))
    if show_sources:
        for f in ("catalog.fits", "catalog.csv"):
            if os.path.exists(os.path.join(d, f)):
                t = Table.read(os.path.join(d, f))
                peaks = np.column_stack([t["x"], t["y"]])
                break
    fig, stats = _draw_final_comparison(
        raw, data, model, var, bad, x0, y0, peaks, 0.75 * fwhm, edge_margin * fwhm,
        view=view, stats_in_view=stats_in_view, show_peaks=show_sources)
    out = out or os.path.join(d, "final_comparison_from_fits.png")
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    if verbose:
        print("wrote", out)
        _print_stats(stats)
    if show:
        plt.show()
    else:
        plt.close(fig)
    return stats

