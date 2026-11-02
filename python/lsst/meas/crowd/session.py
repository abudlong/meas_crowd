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

"""Sessions, configuration and saving of runs.

`prepareSession` loads a visit image, cuts out the requested regions and sets up the PSF;
`makeConfig` builds a configuration (``CROWDED_FIELD_CONFIG`` is the recommended starting
point); `CutoutRun` holds a finished run and `saveRun` writes it to disk.

``makeConfig`` starts from ``DEFAULT_CONFIG``.  A nested dict you pass (e.g. ``shapeKwargs``,
``stackKwargs``, ``innerKwargs``) replaces the whole default dict; pass ``deep=True`` to merge
nested dicts key by key instead (unspecified keys keep their defaults).
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

from ._compat import camelKeys
from . import fitLib as L

# =====================================================================
# Config
# =====================================================================
DEFAULT_CONFIG = dict(
    nOuterIterations=4,
    nInnerIterations=10,
    plotEvery=3,
    sourceView="changed",
    maxPagesPerIter=2,
    sourcesPerPage=4,
    refineCentroids=True,
    centroidMinAmp=50.0,
    nCentroidRounds=3,
    refitPsf="stack",
    stackKwargs=dict(nUse=12, step=0.8, window=(0.85, 1.0), maxRadius=20),
    psfPolishCycles=3,
    psfKernelSizes=(41, 27, 21, 15, 11),
    mergeSeparationSigma=1.5,
    sysFrac=0.03,
    poissonK=0.0,
    minSnr=5.0,
    pruneSepSigma=1.5,
    rejectFwhm=0.75,
    brightPeakSigma=100.0,
    freezeFwhm=1.0,
    nearBrightFwhm=2.0,
    nearBrightSnr=10.0,
    edgeMarginFwhm=0.5,
    useEliminateWeak=True,
    useShapeCut=True,
    shapeKwargs=dict(minSnrLocal=5.0, chi2Max=2.0, widthRange=(0.75, 1.35),
                      dchi2Width=9.0, trustSnr=15.0, minContrast=20.0,
                      weakSnr=10.0, dchi2WidthWeak=6.5),
    shapeInFirstPass=False,
    psfWingSize=81,
    innerKwargs=dict(nIterations=10, thresholdStart=8.0, thresholdStep=-1.5,
                      thresholdMin=4.5, convergenceFrac=0.0),
)


# The recommended settings for crowded fields (deeper detection, looser shape test, closer
# merging than DEFAULT_CONFIG):
#   makeConfig(CROWDED_FIELD_CONFIG, minSnr=4.0, saveTo="my_run.json")
# With the default deep=False, stackKwargs below REPLACES the default one (nUse, step,
# window fall back to the library defaults), exactly as the JSON did via the CLI.
CROWDED_FIELD_CONFIG = dict(
    minSnr=3.0, nearBrightSnr=15.0, rejectFwhm=0.5, nearBrightFwhm=2.0,
    pruneSepSigma=1.0, mergeSeparationSigma=1.0, outerConvergenceFrac=0,
    stackKwargs=dict(maxRadius=20),
    shapeKwargs=dict(minSnrLocal=3.0, chi2Max=3.0, widthRange=(0.7, 1.5),
                      dchi2Width=12.0, trustSnr=20.0, minContrast=1.0,
                      weakSnr=8.0, dchi2WidthWeak=12.0),
    innerKwargs=dict(nIterations=8, thresholdStart=8.0, thresholdStep=-1.0,
                      thresholdMin=3.0, convergenceFrac=0.0),
)

# the old name, kept so existing notebooks keep working
V5_TESTING_CONFIG = CROWDED_FIELD_CONFIG

# CROWDED_FIELD_CONFIG for dense stellar fields (Galactic bulge and plane), where almost every
# source is a star and few stars are isolated, validated on six tiles of LSSTCam detector 144
# (visit 2025062100433):
# - the PSF residual-stack refit uses stars with no neighbour within 1 FWHM (instead of 2.5),
#   at S/N >= 30, up to 80 of them: with the default rule, only a star or two qualify and the
#   PSF is never refined;
# - every source is fitted as a point source (no extended-model test), so extended models
#   cannot absorb PSF mismatch.
CROWDED_STELLAR_FIELD_CONFIG = dict(
    CROWDED_FIELD_CONFIG,
    psfStackIsoFwhm=1.0, psfStackMinSnr=30.0, psfStackN=80,
    ellipTest=False, keepExtended=False,
)


def _merge(dst, src, deep):
    """Merge ``src`` into ``dst`` (key by key for nested dicts when ``deep``) and return ``dst``."""
    for k, v in camelKeys(src).items():           # old snake_case option names still work
        if deep and isinstance(v, dict) and isinstance(dst.get(k), dict):
            _merge(dst[k], v, True)
        else:
            dst[k] = copy.deepcopy(v)
    return dst


def _normalize(cfg):
    """Turn the lists JSON gives back into the tuples the pipeline expects."""
    if isinstance(cfg.get("stackKwargs"), dict) and "window" in cfg["stackKwargs"]:
        cfg["stackKwargs"]["window"] = tuple(cfg["stackKwargs"]["window"])
    if isinstance(cfg.get("shapeKwargs"), dict) and "widthRange" in cfg["shapeKwargs"]:
        cfg["shapeKwargs"]["widthRange"] = tuple(cfg["shapeKwargs"]["widthRange"])
    if "psfKernelSizes" in cfg and cfg["psfKernelSizes"] is not None:
        cfg["psfKernelSizes"] = tuple(cfg["psfKernelSizes"])
    return cfg


def saveConfig(cfg, path):
    """Write a config dict to a JSON file.

    Tuples become lists; `makeConfig` turns them back into tuples when the file is loaded.

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


def makeConfig(base=None, deep=False, nOuter=None, nInner=None, saveTo=None,
                **overrides):
    """Build a pipeline configuration (keyword arguments of `fitLib.runCutoutPipeline`).

    Parameters
    ----------
    base : `str` or `dict`, optional
        A JSON file or a dict, applied on top of ``DEFAULT_CONFIG``.
    deep : `bool`, optional
        Merge nested dicts key by key (default: replace them).
    nOuter : `int`, optional
        Shortcut for ``nOuterIterations``.
    nInner : `int`, optional
        Shortcut for ``nInnerIterations`` and ``innerKwargs['nIterations']``.
    saveTo : `str`, optional
        Also write the configuration to this JSON file.
    **overrides
        Any `fitLib.runCutoutPipeline` keyword, applied last.

    Returns
    -------
    cfg : `dict`
        The configuration.
    """
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    basePath = base if isinstance(base, (str, os.PathLike)) else None
    if base is not None:
        if isinstance(base, (str, os.PathLike)):
            with open(base) as fh:
                base = json.load(fh)
        _merge(cfg, base, deep)
    _merge(cfg, overrides, deep)
    if nOuter:
        cfg["nOuterIterations"] = nOuter
    if nInner:
        cfg["nInnerIterations"] = nInner
        cfg["innerKwargs"]["nIterations"] = nInner
    for k in ("xRange", "yRange", "viz", "saveTo"):   # not pipeline kwargs
        cfg.pop(k, None)
    cfg = _normalize(cfg)
    if saveTo:
        saveTo = os.fspath(saveTo)
        # a bare filename goes next to the base JSON (else the current working directory)
        if (not os.path.dirname(saveTo) and basePath is not None):
            saveTo = os.path.join(os.path.dirname(os.path.abspath(basePath)), saveTo)
        saveConfig(cfg, os.path.abspath(saveTo))
    return cfg


def checkConfig(cfg):
    """Check that every key of a configuration is a `fitLib.runCutoutPipeline` option.

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
    ok = set(inspect.signature(L.runCutoutPipeline).parameters)
    ok -= {"exposureCrop", "psf", "badBitmask", "starCat", "bgImage",
           "xRange", "yRange", "viz"}
    bad = sorted(set(cfg) - ok)
    if bad:
        raise ValueError(f"unknown runCutoutPipeline option(s) in config: {bad}. "
                         f"Valid: {sorted(ok)}")


def showConfig(cfg, diffOnly=False):
    """Print a configuration.

    Parameters
    ----------
    cfg : `dict`
        Pipeline configuration.
    diffOnly : `bool`, optional
        Only print keys that differ from ``DEFAULT_CONFIG``.
    """
    out = {k: v for k, v in cfg.items()
           if not diffOnly or _normalize(copy.deepcopy(DEFAULT_CONFIG)).get(k) != v}
    print(json.dumps(out, indent=2, default=list))


# =====================================================================
# Cutout specs / session
# =====================================================================
def parseCutout(spec):
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
    badBitmask : `int`
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
    psfUsedXy : `numpy.ndarray`, optional
        Positions (detector px, shape (N, 2)) of the stars the detector-level PSF was built from
        (``psfSource='refit'`` only).
    """

    visit: int
    detector: int
    badBitmask: int
    psf: object
    starCat: object
    bgImage: object
    jobs: dict = field(default_factory=dict)     # name -> dict(interior=Box2I, raw=Exposure)
    wcs: object = None
    psfUsedXy: object = None


def buildDetectorPsf(visitImage, psfSource="refit", returnUsed=False):
    """Set up the detector-level PSF.

    With ``psfSource='visit'`` the visit image's own PSF is used.  With 'refit' the background
    is subtracted, sources are detected and measured, stars selected and a PCA PSF determined
    (the visit image is modified in place).

    Parameters
    ----------
    visitImage : `lsst.afw.image.Exposure`
        The detector's visit image.
    psfSource : `str`, optional
        'refit' or 'visit'.
    returnUsed : `bool`, optional
        Also return the positions of the stars the PCA PSF was built from.

    Returns
    -------
    psf : `lsst.afw.detection.Psf`
        Detector-level PSF.
    starCat : `lsst.afw.table.SourceCatalog` or `None`
        Selected PSF stars (None for 'visit').
    bgImage : `lsst.afw.image.Image` or `None`
        Subtracted background (None for 'visit').
    usedXy : `numpy.ndarray` or `None`
        Positions of the stars used, shape (N, 2) (only with ``returnUsed``; None for 'visit').
    """
    if psfSource == "visit":
        print("PSF source: the visit_image's own PSF (no detector-level refit)")
        return (visitImage.getPsf(), None, None) + ((None,) if returnUsed else ())

    print("Detector-level setup: background subtraction -> detection -> PSF determination")
    t0 = time.perf_counter()
    bgTask = L.SubtractBackgroundTask()
    bgTask.config.algorithm = "AKIMA_SPLINE"
    bgTask.config.useApprox = False
    bgImage = bgTask.run(visitImage).background.getImage()

    L.InstallGaussianPsfTask().run(visitImage)
    schema = L.afwTable.SourceTable.makeMinimalSchema()
    detCfg = L.SourceDetectionTask.ConfigClass()
    detCfg.thresholdValue = 3.0
    detCfg.thresholdType = "stdev"
    detCfg.minPixels = 2
    detCfg.doTempLocalBackground = True
    detCfg.tempLocalBackground.binSize = 16
    detCfg.background.binSize = 32
    detection = L.SourceDetectionTask(schema=schema, config=detCfg)
    deblend = L.SourceDeblendTask(schema=schema)
    measure = L.SingleFrameMeasurementTask(schema=schema)
    table = L.afwTable.SourceTable.make(schema)
    sources = detection.run(table, visitImage).sources
    deblend.run(visitImage, sources)
    measure.run(measCat=sources, exposure=visitImage)
    sources = sources.copy(deep=True)
    print(f"  {len(sources)} detections at 3 sigma")

    starCat = L.ObjectSizeStarSelectorTask().run(sourceCat=sources,
                                                 exposure=visitImage).sourceCat
    print(f"  {len(starCat)} PSF-star candidates")
    cands = L.MakePsfCandidatesTask().run(starCat, exposure=visitImage).psfCandidates
    cfg = L.PcaPsfDeterminerTask.ConfigClass()
    cfg.spatialOrder = 1
    cfg.sizeCellX = cfg.sizeCellY = 512
    cfg.maxCandidates = 1500
    cfg.doRejectBlends = True
    cfg.nStarPerCell = 10
    psf, cellSet = L.PcaPsfDeterminerTask(config=cfg).determinePsf(
        exposure=visitImage, psfCandidateList=cands)
    visitImage.setPsf(psf)
    used = [c for cell in cellSet.getCellList() for c in cell if c.getStatus() == L.GOOD_STATUS]
    usedXy = np.array([(c.getXCenter(), c.getYCenter()) for c in used], float).reshape(-1, 2)
    print(f"  PSF determined from {len(used)} stars  ({time.perf_counter() - t0:.1f}s)")
    return (psf, starCat, bgImage) + ((usedXy,) if returnUsed else ())


def prepareSession(visit, detector, cutouts, pad=0, psfSource="refit",
                    instrument="LSSTCam", repo="dp2_prep",
                    collection="LSSTCam/runs/DRP/DP2",
                    visitImageCollections="LSSTCam/runs/*", visitImage=None):
    """Load a visit image, cut out the regions and set up the detector-level PSF.

    The cutouts are taken before the detector-wide background subtraction.

    Parameters
    ----------
    visit : `int`
        Visit id.
    detector : `int`
        Detector id.
    cutouts : `list`
        Cutout specifications (see `parseCutout`); absolute detector pixels.
    pad : `int`, optional
        Fit a region grown by this many px; the interior keeps the requested box.
    psfSource : `str`, optional
        'refit' or 'visit' (see `buildDetectorPsf`).
    instrument : `str`, optional
        Instrument name.
    repo : `str`, optional
        Butler repository path or alias.
    collection : `str`, optional
        Default Butler collection.
    visitImageCollections : `str`, optional
        Collections searched for the visit image.
    visitImage : `lsst.afw.image.Exposure`, optional
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
    if visitImage is None:
        from lsst.daf.butler import Butler
        butler = Butler(repo, collections=collection)
        refs = butler.query_datasets(
            "visit_image", collections=visitImageCollections,
            where=f"instrument='{instrument}' AND visit={visit} AND detector={detector}",
            find_first=False)
        if not refs:
            raise RuntimeError("no visit_image found for that visit/detector")
        visitImage = butler.get("visit_image", collections=visitImageCollections,
                                 dataId=refs[0].dataId)
    detBbox = visitImage.getBBox()
    print(f"Loaded visit_image: bbox {detBbox}")
    mk = visitImage.mask
    badBitmask = (mk.getPlaneBitMask("SAT") | mk.getPlaneBitMask("CR") |
                   mk.getPlaneBitMask("BAD") | mk.getPlaneBitMask("INTRP") |
                   mk.getPlaneBitMask("EDGE"))

    G = L.geom
    jobs = {}
    for spec in cutouts:
        x0, y0, w, h, nm = parseCutout(spec)
        req = G.Box2I(G.Point2I(x0, y0), G.Extent2I(w, h))
        req.clip(detBbox)
        if req.isEmpty():
            print(f"!! cutout {x0},{y0},{w}x{h} lies outside the detector -- skipped")
            continue
        fitBbox = G.Box2I(req)
        fitBbox.grow(pad)
        fitBbox.clip(detBbox)
        name = nm or f"cutout_x{req.getMinX()}_y{req.getMinY()}_{req.getWidth()}x{req.getHeight()}"
        jobs[name] = dict(interior=req, raw=visitImage[fitBbox].clone())

    psf, starCat, bgImage, usedXy = buildDetectorPsf(visitImage, psfSource, returnUsed=True)
    return Session(visit, detector, badBitmask, psf, starCat, bgImage, jobs,
                   wcs=visitImage.getWcs(), psfUsedXy=usedXy)


# =====================================================================
# Gaia (optional)
# =====================================================================
def queryGaia(wcs, bbox, tapUrl="https://gea.esac.esa.int/tap-server/tap"):
    """Query Gaia DR3 at the ESA archive for the sources on a box.

    Parameters
    ----------
    wcs : `lsst.afw.geom.SkyWcs`
        WCS of the image.
    bbox : `lsst.geom.Box2I`
        Pixel box.
    tapUrl : `str`, optional
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
    gaia = pyvo.dal.TAPService(tapUrl).search(q).to_table()
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
    outerHistory : `list` [`dict`]
        One dict per outer pass (``psf``, ``exposure``, ``result``, ...).
    catalog : `astropy.table.Table`
        Source catalogue (see `buildCatalog`).
    summary : `dict`
        Run summary (timings, counts, residual statistics).
    rawCrop : `lsst.afw.image.Exposure`
        The input cutout before background subtraction.
    interiorBbox : `lsst.geom.Box2I`
        The requested box (without padding).
    badBitmask : `int`
        Mask bits of pixels the fit ignored.
    config : `dict`
        Configuration used.
    gaia : `astropy.table.Table`, optional
        The Gaia catalogue whose positions started the run (used for the Gaia comparison columns).
    photometry : `astropy.table.Table`, optional
        The photometry catalogue (see `lsst.meas.crowd.photometry.measurePhotometry`), made at
        the end of the run and written by `saveRun` as ``photometry.fits``.
    """

    name: str
    outerHistory: list
    catalog: Table
    summary: dict
    rawCrop: object
    interiorBbox: object
    badBitmask: int
    config: dict
    gaia: object = None
    photometry: object = None


def buildCatalog(outerHistory, interiorBbox, gaiaTbl=None, matchRadius=2.0):
    """Build the source catalogue of a finished run.

    Parameters
    ----------
    outerHistory : `list` [`dict`]
        Outer-loop history.
    interiorBbox : `lsst.geom.Box2I`
        The requested box (``inInterior`` column).
    gaiaTbl : `astropy.table.Table`, optional
        Gaia sources (with ``x``, ``y``) to match.
    matchRadius : `float`, optional
        Gaia match radius [px].

    Returns
    -------
    cat : `astropy.table.Table`
        Positions (pixel and sky), PSF flux, AB magnitude (assuming nJy), source metrics,
        ``inInterior`` and, with ``gaiaTbl``, the nearest Gaia source.
    """
    final = outerHistory[-1]
    exp, res = final["exposure"], final["result"]
    peaks = np.asarray(res["peaks"], float).reshape(-1, 2)
    n = len(peaks)
    cat = Table()
    cat["id"] = np.arange(n)
    cat["x"], cat["y"] = peaks[:, 0], peaks[:, 1]
    xy0 = exp.getXY0()
    cat["xCutout"], cat["yCutout"] = peaks[:, 0] - xy0.getX(), peaks[:, 1] - xy0.getY()
    flux = np.asarray(res["amps"], float) if n else np.zeros(0)
    cat["psfFlux_nJy"] = flux
    with np.errstate(divide="ignore", invalid="ignore"):
        cat["magAB"] = np.where(flux > 0, -2.5 * np.log10(flux) + 31.4, np.nan)
    try:
        ra, dec = exp.getWcs().pixelToSkyArray(peaks[:, 0], peaks[:, 1], degrees=True)
        cat["ra"], cat["dec"] = np.asarray(ra), np.asarray(dec)
    except Exception as e:
        print(f"  (no sky coords in catalog: {e})")
    if n and res.get("snaps"):
        try:
            m = L.sourceMetrics(outerHistory)
            cat["snr"], cat["chi2Core"] = m["snr"], m["chi2"]
            cat["nnPix"], cat["edgePix"], cat["groupSize"] = m["nn"], m["edge"], m["group"]
        except Exception as e:
            print(f"  (sourceMetrics failed: {e})")
    ib = interiorBbox
    cat["inInterior"] = ((peaks[:, 0] >= ib.getMinX() - .5) & (peaks[:, 0] <= ib.getMaxX() + .5) &
                          (peaks[:, 1] >= ib.getMinY() - .5) & (peaks[:, 1] <= ib.getMaxY() + .5))
    if gaiaTbl is not None and n:
        from scipy.spatial import cKDTree
        d, i = cKDTree(np.column_stack([gaiaTbl["x"], gaiaTbl["y"]])).query(peaks)
        cat["gaiaSourceId"] = np.asarray(gaiaTbl["source_id"])[i]
        cat["gaiaG"] = np.asarray(gaiaTbl["phot_g_mean_mag"], float)[i]
        cat["gaiaSepPix"] = d
        cat["gaiaMatch"] = d <= matchRadius
    return cat


def runCutout(rawCrop, psf, badBitmask, config=None, interiorBbox=None, starCat=None,
               bgImage=None, viz=False, diagnostics=True, gaiaTbl=None, matchRadius=2.0,
               name="cutout", **overrides):
    """Run the pipeline on one cropped exposure.

    Parameters
    ----------
    rawCrop : `lsst.afw.image.Exposure`
        The cutout (before background subtraction).
    psf : `lsst.afw.detection.Psf`
        Initial PSF.
    badBitmask : `int`
        Mask bits of pixels to ignore.
    config : `dict`, optional
        Pipeline configuration (default ``DEFAULT_CONFIG``).
    interiorBbox : `lsst.geom.Box2I`, optional
        The requested box (default: the whole cutout).
    starCat : `lsst.afw.table.SourceCatalog`, optional
        PSF stars forced into every PSF determination.
    bgImage : `lsst.afw.image.Image`, optional
        Detector background subtracted before the first pass.
    viz : `bool`, optional
        Draw the per-iteration diagnostics.
    diagnostics : `bool`, optional
        Run the end-of-run diagnostics.
    gaiaTbl : `astropy.table.Table`, optional
        Gaia sources to match in the catalogue.
    matchRadius : `float`, optional
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
    cfg = makeConfig(config, **overrides)   # idempotent for an already-complete config
    checkConfig(cfg)
    cb = rawCrop.getBBox()
    interiorBbox = interiorBbox or cb
    t0 = time.perf_counter()
    hist = L.runCutoutPipeline(rawCrop, psf, badBitmask, starCat=starCat, bgImage=bgImage,
                                 xRange=(0, cb.getWidth()), yRange=(0, cb.getHeight()),
                                 viz=viz, **cfg)
    res = hist[-1]["result"]
    summary = dict(name=name, status="ok",
                   bboxFitted=[cb.getMinX(), cb.getMaxX(), cb.getMinY(), cb.getMaxY()],
                   interiorBbox=[interiorBbox.getMinX(), interiorBbox.getMaxX(),
                                  interiorBbox.getMinY(), interiorBbox.getMaxY()],
                   nSources=len(res["peaks"]),
                   residStdPerOuter=[h["residStd"] for h in hist],
                   nPeaksPerOuter=[h["nPeaks"] for h in hist],
                   psfFwhmPxPerOuter=[L._fwhm(h["psf"], h["exposure"].getBBox().getCenter())
                                          for h in hist],
                   pipelineKwargs=cfg)
    if diagnostics and len(res["peaks"]) and res.get("snaps"):
        try:
            summary["highfreqNoiseK"] = L.highfreqNoiseCheck(hist)["k"]
        except Exception as e:
            print(f"  highfreqNoiseCheck failed: {e}")
        try:
            snrMap = L.negativeSourceCheck(hist, plot=viz)
            summary["negativeCheck"] = {str(t): dict(pos=int((snrMap >= t).sum()),
                                                      neg=int((snrMap <= -t).sum()))
                                         for t in (3, 4, 5)}
        except Exception as e:
            print(f"  negativeSourceCheck failed: {e}")
    cat = buildCatalog(hist, interiorBbox, gaiaTbl, matchRadius)
    summary["nInInterior"] = int(np.sum(cat["inInterior"]))
    summary["wallSeconds"] = time.perf_counter() - t0
    return CutoutRun(name, hist, cat, summary, rawCrop, interiorBbox, badBitmask, cfg)


def runSession(session, config=None, names=None, outDir=None, viz=False, gaia=False,
                matchRadius=2.0, plotFinal=True, **overrides):
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
    outDir : `str`, optional
        Save each run under ``<outDir>/visit<V>_det<D>/<name>/`` (see `saveRun`).
    viz : `bool`, optional
        Draw the per-iteration diagnostics.
    gaia : `bool`, optional
        Match the catalogues to Gaia (queried at the ESA archive).
    matchRadius : `float`, optional
        Gaia match radius [px].
    plotFinal : `bool`, optional
        Plot the final comparison of each cutout.
    **overrides
        Configuration overrides.

    Returns
    -------
    runs : `dict` [`str`, `CutoutRun`]
        The results by cutout name.
    """
    cfg = makeConfig(config, **overrides)
    runs = {}
    todo = names or list(session.jobs)
    for k, name in enumerate(todo, 1):
        job = session.jobs[name]
        print(f"{'#' * 78}\n[{k}/{len(todo)}] {name}  interior={job['interior']}")
        try:
            gaiaTbl = None
            if gaia:
                try:
                    gaiaTbl = queryGaia(session.wcs, job["raw"].getBBox())
                except Exception as e:
                    print(f"Gaia query failed ({type(e).__name__}: {e}) -- continuing without")
            run = runCutout(job["raw"], session.psf, session.badBitmask, cfg,
                             interiorBbox=job["interior"], starCat=session.starCat,
                             bgImage=session.bgImage, viz=viz, gaiaTbl=gaiaTbl,
                             matchRadius=matchRadius, name=name)
            runs[name] = run
            print(f"{name}: {run.summary['nSources']} sources "
                  f"({run.summary['nInInterior']} in interior), "
                  f"{run.summary['wallSeconds']:.0f}s")
            if outDir:
                saveRun(run, os.path.join(outDir,
                                           f"visit{session.visit}_det{session.detector}", name))
            if plotFinal:
                plotFinalComparisonRobust(run)
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


def _wcsHeader(exposure):
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


def writePhotometry(run, outdir, verbose=True):
    """Write a run's photometry catalogue as photometry.fits and photometry.csv.

    The catalogue is computed (with the Gaia comparison when the run has its Gaia table) if the
    run does not have one yet.  A failure is reported, not raised, so the other products are
    still written.

    Parameters
    ----------
    run : `CutoutRun`
        The run.
    outdir : `str`
        Output directory (created if needed).
    verbose : `bool`, optional
        Print where the catalogue went.

    Returns
    -------
    path : `str` or `None`
        The FITS file written (None if the photometry failed).
    """
    from .photometry import fitsMeta, measurePhotometry
    os.makedirs(outdir, exist_ok=True)
    try:
        phot = getattr(run, "photometry", None)
        if phot is None:
            phot = measurePhotometry(run, verbose=False)
            run.photometry = phot
        out = phot.copy(copy_data=False)
        out.meta = fitsMeta(phot.meta)
        path = os.path.join(outdir, "photometry.fits")
        out.write(path, overwrite=True)
        out.write(os.path.join(outdir, "photometry.csv"), format="ascii.csv", overwrite=True)
        if verbose:
            print(f"  photometry catalogue ({len(phot)} sources) -> {path}")
        return path
    except Exception as e:
        print(f"  photometry catalogue not written ({type(e).__name__}: {e})")
        return None


def saveRun(run, outdir):
    """Write a run's catalogues, images, PSF kernels and summary to a directory.

    Writes photometry.fits / photometry.csv (the photometry catalogue with the Gaia comparison,
    computed here if the run has none), catalog.csv / catalog.fits, images.fits,
    psf_kernels.fits and summary.json, the layout `plotFinalComparisonFromFits` reads.

    Parameters
    ----------
    run : `CutoutRun`
        The run.
    outdir : `str`
        Output directory (created if needed).
    """
    from astropy.io import fits
    os.makedirs(outdir, exist_ok=True)
    final = run.outerHistory[-1]
    exp, res = final["exposure"], final["result"]
    xy0 = exp.getXY0()
    hdr = fits.Header()
    for k, v in _wcsHeader(exp).items():
        try:
            hdr[k] = v
        except Exception:
            pass
    hdr["X0"], hdr["Y0"] = xy0.getX(), xy0.getY()
    model, data = res["modelFull"], exp.image.array
    hdus = [fits.PrimaryHDU(header=hdr)]
    for nm, arr in (("DATA_USED", data), ("MODEL", model), ("RESIDUAL", data - model),
                    ("VARIANCE", exp.variance.array), ("RAW_CUTOUT", run.rawCrop.image.array)):
        hdus.append(fits.ImageHDU(np.asarray(arr, dtype=np.float32), header=hdr, name=nm))
    hdus.append(fits.ImageHDU(np.asarray(exp.mask.array, dtype=np.int32), header=hdr, name="MASK"))
    hdus.append(fits.ImageHDU(((exp.mask.array & run.badBitmask) != 0).astype(np.uint8),
                              header=hdr, name="BAD_PIXELS"))
    fits.HDUList(hdus).writeto(os.path.join(outdir, "images.fits"), overwrite=True)

    kHdus = [fits.PrimaryHDU()]
    for i, h in enumerate(run.outerHistory):
        k = h["psf"].computeKernelImage(h["exposure"].getBBox().getCenter()).array
        kHdus.append(fits.ImageHDU(np.asarray(k, dtype=np.float32), name=f"PSF_OUTER{i + 1}"))
    fits.HDUList(kHdus).writeto(os.path.join(outdir, "psf_kernels.fits"), overwrite=True)

    writePhotometry(run, outdir)
    run.catalog.write(os.path.join(outdir, "catalog.csv"), format="ascii.csv", overwrite=True)
    run.catalog.write(os.path.join(outdir, "catalog.fits"), overwrite=True)
    with open(os.path.join(outdir, "summary.json"), "w") as fh:
        json.dump(run.summary, fh, indent=2, default=_jsonable)
    with open(os.path.join(outdir, "config.json"), "w") as fh:     # reload with makeConfig()
        json.dump(run.config, fh, indent=2, default=_jsonable)
    print(f"saved run to {outdir}")


# =====================================================================
# Final comparison plot with robust residual statistics
# =====================================================================
EDGE_COLOR = "#ff00ff"
IQR_TO_SIGMA = 1.349          # sigma_robust = (Q75 - Q25) / 1.349 for a Gaussian


def robustStats(values, valid=None):
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
        ``median``, ``sigmaIqr`` ((Q75 - Q25) / 1.349), ``q25``, ``q75`` and ``n``.
    """
    v = np.asarray(values, float)
    m = np.isfinite(v) if valid is None else (np.asarray(valid, bool) & np.isfinite(v))
    x = v[m]
    if x.size == 0:
        return dict(median=np.nan, sigmaIqr=np.nan, q25=np.nan, q75=np.nan, n=0)
    q25, med, q75 = np.percentile(x, [25, 50, 75])
    return dict(median=float(med), sigmaIqr=float((q75 - q25) / IQR_TO_SIGMA),
                q25=float(q25), q75=float(q75), n=int(x.size))


def _drawEdgeMargin(ax, ext, marginPx, shade=True):
    """Outline and shade the band within ``marginPx`` of the image edge."""
    if not marginPx or marginPx <= 0:
        return
    xa, xb, ya, yb = ext
    if shade:
        for (a, b, c, d) in ((xa, xb, ya, ya + marginPx), (xa, xb, yb - marginPx, yb),
                             (xa, xa + marginPx, ya, yb), (xb - marginPx, xb, ya, yb)):
            ax.add_patch(patches.Rectangle((a, c), b - a, d - c, fc=EDGE_COLOR, ec="none",
                                           alpha=0.18, zorder=3))
    ax.add_patch(patches.Rectangle((xa + marginPx, ya + marginPx),
                                   (xb - xa) - 2 * marginPx, (yb - ya) - 2 * marginPx,
                                   fill=False, ec=EDGE_COLOR, lw=1.6, ls="--", zorder=4))


def _drawFinalComparison(raw, data, model, var, bad, x0, y0, peaks, radius, marginPx,
                           view=None, statsInView=True, showPeaks=True,
                           figsize=(26, 7), residClipPct=99.0, zClip=5.0):
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
    useZ = use & varok & np.isfinite(z)
    region = np.ones(data.shape, bool)
    if view is not None and statsInView:
        X = x0 + np.arange(nx)[None, :]
        Y = y0 + np.arange(ny)[:, None]
        region = (X >= view[0]) & (X <= view[1]) & (Y >= view[2]) & (Y <= view[3])
    stR = robustStats(resid, use & region)
    stZ = robustStats(z, useZ & region)

    cmRes = plt.get_cmap("RdBu_r").copy()
    cmRes.set_bad("black")
    resM = np.ma.masked_where(~use, resid)
    zM = np.ma.masked_where(~useZ, z)
    sel = resid[use & region]
    rs = np.nanpercentile(np.abs(sel), residClipPct) if sel.size else 1.0

    norm = simple_norm(raw, vmin=np.nanmin(raw), vmax=np.nanpercentile(raw, 99),
                       stretch="log", log_a=1000)
    cmImg = L.CMAP

    def sub(st, unit=""):
        return f"median {st['median']:+.2f}{unit},  robust σ {st['sigmaIqr']:.2f}{unit}  (N={st['n']})"

    panels = [
        ("Original visit image (cutout)", raw, dict(cmap=cmImg, norm=norm)),
        (f"Final model ({len(peaks)} sources)", model, dict(cmap=cmImg, norm=norm)),
        ("Residual = data − model\n" + sub(stR), resM, dict(cmap=cmRes, vmin=-rs, vmax=rs)),
        ("Residual / σ\n" + sub(stZ), zM, dict(cmap=cmRes, vmin=-zClip, vmax=zClip)),
    ]
    fig, axes = plt.subplots(1, 4, figsize=figsize, sharex=True, sharey=True)
    for i, (ax, (title, img, kw)) in enumerate(zip(axes, panels)):
        im = ax.imshow(img, extent=ext, origin="lower", **kw)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        if i < 2 and bad.any():                              # red overlay on data / model
            ax.imshow(np.ma.masked_where(~bad, np.ones(bad.shape)), extent=ext,
                      origin="lower", cmap=matplotlibListed("red"), alpha=0.35,
                      interpolation="nearest")
        _drawEdgeMargin(ax, ext, marginPx)
        if showPeaks and len(peaks):
            for px, py in peaks:
                ax.add_patch(patches.Circle((px, py), radius, fill=False, ec="w",
                                            lw=0.5, alpha=0.6))
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("x [px]")
        ax.grid(False)
        if i >= 2:                                           # stats box on the residual panels
            st = stR if i == 2 else stZ
            ax.text(0.02, 0.98, f"median = {st['median']:+.3f}\nσ_IQR = {st['sigmaIqr']:.3f}",
                    transform=ax.transAxes, va="top", ha="left", fontsize=10,
                    bbox=dict(fc="white", ec="k", alpha=0.85, pad=3), zorder=10)
    axes[0].set_ylabel("y [px]")
    if view is not None:
        axes[0].set_xlim(view[0], view[1])
        axes[0].set_ylim(view[2], view[3])
    else:
        axes[0].set_xlim(ext[0], ext[1])
        axes[0].set_ylim(ext[2], ext[3])
    handles = [patches.Patch(fc="red", alpha=0.35, label="masked (badBitmask)")]
    if marginPx:
        handles.insert(0, patches.Patch(fc=EDGE_COLOR, alpha=0.3, ec=EDGE_COLOR, ls="--",
                                        label=f"edge drop zone ({marginPx:.1f}px)"))
    axes[0].legend(handles=handles, loc="upper left", fontsize=8)
    fig.text(0.5, -0.01,
             "Residual stats: bad/non-finite pixels excluded (black in residual panels); "
             "median and robust σ = (Q75−Q25)/1.349"
             + (" over the plotted view" if (view is not None and statsInView) else
                " over the full array"),
             ha="center", fontsize=9)
    plt.tight_layout()
    return fig, dict(residual=stR, residualOverSigma=stZ)


def matplotlibListed(color):
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


def _printStats(stats):
    """Print residual statistics."""
    for key, st in stats.items():
        print(f"{key:20s} median={st['median']:+9.3f}  robust σ(IQR/1.349)={st['sigmaIqr']:8.3f}  "
              f"Q25={st['q25']:+.3f}  Q75={st['q75']:+.3f}  N={st['n']}")


def plotFinalComparisonRobust(runOrHistory, exposureCrop=None, badBitmask=None,
                                 edgeMarginFwhm=None, view="interior", statsInView=True,
                                 showPeaks=True, figsize=(26, 7), verbose=True):
    """Plot original | final model | residual | residual/sigma with robust statistics.

    Bad pixels are black in the residual panels; median and IQR-robust sigma (good pixels only)
    are shown in the titles.

    Parameters
    ----------
    runOrHistory : `CutoutRun` or `list` [`dict`]
        A run (preferred) or its outer-loop history.
    exposureCrop : `lsst.afw.image.Exposure`, optional
        Raw crop (only needed with a bare history).
    badBitmask : `int`, optional
        Mask bits of bad pixels (needed with a bare history).
    edgeMarginFwhm : `float`, optional
        Draw the edge drop-zone of this width [FWHM].
    view : `str` or `tuple`, optional
        'interior' (the run's interior box), None (the whole fitted array) or
        (xmin, xmax, ymin, ymax) in detector pixels.
    statsInView : `bool`, optional
        Restrict the statistics to the displayed view.
    showPeaks : `bool`, optional
        Mark the sources.
    figsize : `tuple` [`float`, `float`], optional
        Figure size [inches].
    verbose : `bool`, optional
        Print the statistics.

    Returns
    -------
    stats : `dict`
        ``residual`` and ``residualOverSigma``: median, sigmaIqr, q25, q75, n.

    Raises
    ------
    ValueError
        Raised if a bare history is passed without ``badBitmask``.
    """
    if isinstance(runOrHistory, CutoutRun):
        run = runOrHistory
        hist, rawExp = run.outerHistory, run.rawCrop
        badBitmask = run.badBitmask if badBitmask is None else badBitmask
        if edgeMarginFwhm is None:
            edgeMarginFwhm = run.config.get("edgeMarginFwhm", 0.5)
        if isinstance(view, str) and view == "interior":
            b = run.interiorBbox
            view = (b.getMinX() - .5, b.getMaxX() + .5, b.getMinY() - .5, b.getMaxY() + .5)
    else:
        hist, rawExp = runOrHistory, exposureCrop
        if badBitmask is None:
            raise ValueError("badBitmask is required when passing a bare outerHistory")
        if edgeMarginFwhm is None:
            edgeMarginFwhm = 0.5
        if isinstance(view, str):
            view = None
    final = hist[-1]
    exp, res = final["exposure"], final["result"]
    bbox = exp.getBBox()
    model = res["modelFull"]
    data = exp.image.array
    raw = rawExp[bbox].image.array if rawExp is not None else data
    bad = (exp.mask.array & badBitmask) != 0
    peaks = np.asarray(res["peaks"], float).reshape(-1, 2)
    fwhm = L._fwhm(exp.getPsf(), bbox.getCenter())
    fig, stats = _drawFinalComparison(
        raw, data, model, exp.variance.array, bad, bbox.getMinX(), bbox.getMinY(), peaks,
        0.75 * fwhm, edgeMarginFwhm * fwhm if edgeMarginFwhm else 0.0,
        view=view, statsInView=statsInView, showPeaks=showPeaks, figsize=figsize)
    if verbose:
        _printStats(stats)
    return stats


def plotFinalComparisonFromFits(d, view=None, out=None, statsInView=True,
                                    showSources=True, show=True, dpi=110, verbose=True):
    """Rebuild the final-comparison figure from files written by `saveRun`.

    No LSST stack and no re-fit needed.

    Parameters
    ----------
    d : `str`
        Directory written by `saveRun`.
    view : `tuple` [`float`, `float`, `float`, `float`], optional
        (xmin, xmax, ymin, ymax) in detector pixels (default: the interior box from
        summary.json, else the whole array).
    out : `str`, optional
        Also save the figure to this PNG.
    statsInView : `bool`, optional
        Restrict the statistics to the view.
    showSources : `bool`, optional
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
        See `plotFinalComparisonRobust`.
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
    fwhm = (summ.get("psfFwhmPxPerOuter") or [3.5])[-1]
    edgeMargin = (summ.get("pipelineKwargs") or {}).get("edgeMarginFwhm", 0.5)
    if view is None and "interiorBbox" in summ:
        b = summ["interiorBbox"]
        view = (b[0] - .5, b[1] + .5, b[2] - .5, b[3] + .5)
    peaks = np.zeros((0, 2))
    if showSources:
        for f in ("catalog.fits", "catalog.csv"):
            if os.path.exists(os.path.join(d, f)):
                t = Table.read(os.path.join(d, f))
                peaks = np.column_stack([t["x"], t["y"]])
                break
    fig, stats = _drawFinalComparison(
        raw, data, model, var, bad, x0, y0, peaks, 0.75 * fwhm, edgeMargin * fwhm,
        view=view, statsInView=statsInView, showPeaks=showSources)
    out = out or os.path.join(d, "final_comparison_from_fits.png")
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    if verbose:
        print("wrote", out)
        _printStats(stats)
    if show:
        plt.show()
    else:
        plt.close(fig)
    return stats



# old snake_case names still work (with a FutureWarning)
from ._compat import moduleGetattr  # noqa: E402
__getattr__ = moduleGetattr(globals())
