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

"""The Gaia DR3 catalogue of a visit and the Gaia starting positions of each cutout.

Butler reference catalogue (or ESA TAP), proper motions propagated to the visit epoch, and
the selection of the Gaia stars that start the fit of a cutout.
"""

import os
import time
import numpy as np
from scipy.spatial import cKDTree
from astropy.table import Table
from astropy.time import Time
from . import fitLib as L


GAIA_EPOCH = 2016.0          # Gaia DR3 reference epoch (TCB Julian year)


# =====================================================================
# Gaia catalog for a visit image
# =====================================================================
def visitEpoch(exposure):
    """Return the observation epoch of an exposure.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure with a VisitInfo.

    Returns
    -------
    epoch : `astropy.time.Time` or `None`
        The visit date, or None if it is missing or implausible.
    """
    try:
        d = exposure.getInfo().getVisitInfo().getDate()
        try:
            t = d.toAstropy()
        except AttributeError:
            from lsst.daf.base import DateTime
            t = Time(d.get(system=DateTime.MJD, scale=DateTime.TAI), format="mjd", scale="tai")
        return t if np.isfinite(t.jyear) and t.jyear > 2000 else None
    except Exception:
        return None


def _col(tbl, name, fill=np.nan):
    """Return a table column as a float array, or ``fill`` if the column is missing."""
    c = tbl[name]
    if hasattr(c, "filled"):
        c = c.filled(fill)
    return np.asarray(c, float)


def _skyCorners(wcs, bbox):
    """Return the sky positions of a box's corners."""
    pts = [(bbox.getMinX(), bbox.getMinY()), (bbox.getMinX(), bbox.getMaxY()),
           (bbox.getMaxX(), bbox.getMaxY()), (bbox.getMaxX(), bbox.getMinY())]
    return [wcs.pixelToSky(float(x), float(y)) for x, y in pts]


def _gaiaFromButler(wcs, bbox, epoch, butler, refcat, collections):
    """Load Gaia from the Butler reference catalogue (proper motions applied by the loader)."""
    import lsst.sphgeom as sph
    from lsst.meas.algorithms import ReferenceObjectLoader
    poly = sph.ConvexPolygon.convexHull(
        [sph.UnitVector3d(sph.LonLat.fromDegrees(c.getRa().asDegrees(), c.getDec().asDegrees()))
         for c in _skyCorners(wcs, bbox)])
    shards = [i for b, e in sph.HtmPixelization(7).envelope(poly) for i in range(b, e)]
    where = f"htm7 IN ({','.join(map(str, shards))})"
    try:        # the loader needs each shard's region -> data IDs WITH dimension records
        refs = list(butler.query_datasets(refcat, collections=collections, where=where,
                                          with_dimension_records=True))
    except TypeError:
        refs = list(butler.registry.queryDatasets(refcat, collections=collections,
                                                  where=where).expanded())
    if not refs:
        raise RuntimeError(f"no {refcat} shards in {collections}")
    loader = ReferenceObjectLoader(dataIds=[r.dataId for r in refs],
                                   refCats=[butler.getDeferred(r) for r in refs], name=refcat)
    res = loader.loadPixelBox(bbox, wcs, filterName="phot_g_mean", epoch=epoch)
    cat = res.refCat.copy(deep=True)

    def mag(field):
        try:
            f = np.asarray(cat[field], float)
        except Exception:
            return np.full(len(cat), np.nan)
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(f > 0, -2.5 * np.log10(f) + 31.4, np.nan)

    t = Table(dict(source_id=np.asarray(cat["id"], np.int64),
                   ra=np.degrees(np.asarray(cat["coord_ra"], float)),
                   dec=np.degrees(np.asarray(cat["coord_dec"], float)),
                   phot_g_mean_mag=mag(res.fluxField),
                   phot_bp_mean_mag=mag("phot_bp_mean_flux"),
                   phot_rp_mean_mag=mag("phot_rp_mean_flux")))
    t.meta["gaiaMagSystem"] = "AB, from the reference catalogue fluxes [nJy]"
    return t


def _gaiaFromTap(wcs, bbox, epoch, tapUrl):
    """Query Gaia at the ESA archive and propagate the positions to ``epoch``."""
    import pyvo
    corners = _skyCorners(wcs, bbox)
    center = wcs.pixelToSky(bbox.getCenter())
    radius = max(center.separation(c).asDegrees() for c in corners) * 1.05
    ra0, dec0 = center.getRa().asDegrees(), center.getDec().asDegrees()
    q = f"""SELECT source_id, ra, dec, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag,
                   phot_rp_mean_mag, ruwe
            FROM gaiadr3.gaia_source
            WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra0}, {dec0}, {radius}))"""
    svc = pyvo.dal.TAPService(tapUrl)
    try:
        t = svc.run_async(q, maxrec=10 ** 7).to_table()
    except Exception:
        t = svc.search(q, maxrec=10 ** 7).to_table()
    ra, dec = _col(t, "ra"), _col(t, "dec")
    dt = (epoch.jyear - GAIA_EPOCH) if epoch is not None else 0.0
    pmra, pmdec = _col(t, "pmra", 0.0), _col(t, "pmdec", 0.0)     # mas/yr (pmra has cos dec)
    decE = dec + pmdec * dt / 3.6e6
    raE = ra + pmra * dt / 3.6e6 / np.cos(np.radians(dec))
    out = Table(dict(source_id=np.asarray(t["source_id"], np.int64), ra=raE, dec=decE,
                     phot_g_mean_mag=_col(t, "phot_g_mean_mag"),
                     phot_bp_mean_mag=_col(t, "phot_bp_mean_mag"),
                     phot_rp_mean_mag=_col(t, "phot_rp_mean_mag"), ruwe=_col(t, "ruwe")))
    out.meta["gaiaMagSystem"] = "Vega (Gaia archive)"
    return out


def loadGaiaCatalog(wcs, bbox, epoch=None, source="auto", butler=None, repo="dp2_prep",
                    refcat="gaia_dr3_20230707", refcatCollections="refcats",
                    tapUrl="https://gea.esac.esa.int/tap-server/tap", cachePath=None,
                    pad=0, verbose=True):
    """Load the Gaia DR3 stars overlapping a box, with positions at the observation epoch.

    Parameters
    ----------
    wcs : `lsst.afw.geom.SkyWcs`
        WCS of the image.
    bbox : `lsst.geom.Box2I`
        Pixel box (detector pixels of ``wcs``).
    epoch : `astropy.time.Time`, optional
        Epoch of the positions (default: the Gaia epoch, J2016.0).
    source : `str`, optional
        'butler' (LSST reference catalogue), 'tap' (ESA archive) or 'auto' (Butler, then TAP).
    butler : `lsst.daf.butler.Butler`, optional
        Butler to use (default: one for ``repo``).
    repo : `str`, optional
        Butler repository path or alias.
    refcat : `str`, optional
        Reference catalogue dataset type.
    refcatCollections : `str`, optional
        Collections of the reference catalogue.
    tapUrl : `str`, optional
        TAP service URL.
    cachePath : `str`, optional
        .ecsv file: read if it exists, written otherwise.
    pad : `int`, optional
        Grow the box by this many px.
    verbose : `bool`, optional
        Print what was loaded.

    Returns
    -------
    gaia : `astropy.table.Table`
        ``source_id``, ``ra``, ``dec``, ``phot_g_mean_mag``, pixel positions ``x``, ``y`` and
        further Gaia columns; ``meta['epochJyear']`` is the epoch used.

    Raises
    ------
    RuntimeError
        Raised if Gaia could not be loaded from any source.
    """
    if cachePath and os.path.exists(cachePath):
        tbl = Table.read(cachePath)
        if verbose:
            print(f"Gaia: {len(tbl)} stars read from cache {cachePath}")
        return tbl
    G = L.geom
    qbox = G.Box2I(bbox)
    qbox.grow(int(pad))
    errs, tbl = [], None
    for src in (["butler", "tap"] if source == "auto" else [source]):
        try:
            t0 = time.perf_counter()
            if src == "butler":
                if butler is None:
                    from lsst.daf.butler import Butler
                    butler = Butler(repo)
                tbl = _gaiaFromButler(wcs, qbox, epoch, butler, refcat, refcatCollections)
            else:
                tbl = _gaiaFromTap(wcs, qbox, epoch, tapUrl)
            if verbose:
                print(f"Gaia: {len(tbl)} stars from {src} ({time.perf_counter() - t0:.1f}s, "
                      f"epoch {epoch.jyear if epoch is not None else GAIA_EPOCH:.2f})")
            break
        except Exception as e:
            errs.append(f"{src}: {type(e).__name__}: {e}")
            if verbose:
                print(f"Gaia via {src} failed -- {errs[-1]}")
    if tbl is None:
        raise RuntimeError("could not load Gaia: " + " | ".join(errs))
    x, y = wcs.skyToPixelArray(np.asarray(tbl["ra"], float), np.asarray(tbl["dec"], float),
                               degrees=True)
    tbl["x"], tbl["y"] = np.asarray(x, float), np.asarray(y, float)
    inside = ((tbl["x"] >= qbox.getMinX() - .5) & (tbl["x"] <= qbox.getMaxX() + .5) &
              (tbl["y"] >= qbox.getMinY() - .5) & (tbl["y"] <= qbox.getMaxY() + .5))
    tbl = tbl[inside]
    tbl.meta["epochJyear"] = float(epoch.jyear) if epoch is not None else GAIA_EPOCH
    if cachePath:
        os.makedirs(os.path.dirname(os.path.abspath(cachePath)), exist_ok=True)
        tbl.write(cachePath, overwrite=True)
    return tbl


def attachGaia(session, visitImage=None, pad=64, cacheDir="gaia_cache", **kw):
    """Load Gaia once for a session and store it as ``session.gaia``.

    Covers the whole visit image or, without one, the union of the session's cutouts; the
    catalogue is cached on disk.

    Parameters
    ----------
    session : `lsst.meas.crowd.session.Session`
        The session (``session.gaia`` is set).
    visitImage : `lsst.afw.image.Exposure`, optional
        The full visit image.
    pad : `int`, optional
        Grow the box by this many px.
    cacheDir : `str`, optional
        Cache directory (None = no cache).
    **kw
        Keyword arguments of `loadGaiaCatalog` (e.g. ``repo``).

    Returns
    -------
    gaia : `astropy.table.Table`
        The Gaia catalogue (see `loadGaiaCatalog`).
    """
    G = L.geom
    if visitImage is not None:
        bbox = visitImage.getBBox()
        epoch = visitEpoch(visitImage)
    else:
        boxes = [j["raw"].getBBox() for j in session.jobs.values()]
        bbox = G.Box2I(boxes[0])
        for b in boxes[1:]:
            bbox.include(b)
        epoch = visitEpoch(next(iter(session.jobs.values()))["raw"])
    if epoch is None:
        print("WARNING: no visit date found -- Gaia positions are NOT proper-motion corrected")
    cache = None
    if cacheDir:
        cache = os.path.join(cacheDir, f"gaia_visit{session.visit}_det{session.detector}_"
                             f"{bbox.getMinX()}_{bbox.getMinY()}_{bbox.getWidth()}x"
                             f"{bbox.getHeight()}.ecsv")
    session.gaia = loadGaiaCatalog(session.wcs, bbox, epoch, cachePath=cache, pad=pad, **kw)
    return session.gaia


# =====================================================================
# Gaia starting positions of a cutout
# =====================================================================
def selectGaiaSources(gaiaTbl, bbox, fwhm, padFwhm=1.5, gmagMax=None, minSepPx=0.5):
    """Select the Gaia stars whose positions start the fit of a cutout.

    Stars inside the box grown by ``padFwhm * fwhm``, brightest first, with near-duplicates
    (closer than ``minSepPx``) removed.

    Parameters
    ----------
    gaiaTbl : `astropy.table.Table` or `None`
        Gaia catalogue with ``x``, ``y`` (and ``phot_g_mean_mag``).
    bbox : `lsst.geom.Box2I`
        Cutout box.
    fwhm : `float`
        PSF FWHM [px].
    padFwhm : `float`, optional
        Growth of the box [FWHM].
    gmagMax : `float`, optional
        Only stars brighter than this G magnitude.
    minSepPx : `float`, optional
        Drop fainter stars within this many px of a brighter one.

    Returns
    -------
    xy : `numpy.ndarray`
        Starting positions, shape (N, 2).
    rows : `numpy.ndarray` [`int`]
        Their rows in ``gaiaTbl``.
    """
    if gaiaTbl is None or len(gaiaTbl) == 0:
        return np.zeros((0, 2)), np.zeros(0, np.int64)
    x, y = np.asarray(gaiaTbl["x"], float), np.asarray(gaiaTbl["y"], float)
    g = np.asarray(gaiaTbl["phot_g_mean_mag"], float) if "phot_g_mean_mag" in gaiaTbl.colnames \
        else np.zeros(len(x))
    pad = padFwhm * fwhm
    ok = (np.isfinite(x) & np.isfinite(y) & (x >= bbox.getMinX() - pad) & (x <= bbox.getMaxX() + pad)
          & (y >= bbox.getMinY() - pad) & (y <= bbox.getMaxY() + pad))
    if gmagMax is not None:
        ok &= np.nan_to_num(g, nan=99) <= gmagMax
    idx = np.flatnonzero(ok)
    idx = idx[np.argsort(np.nan_to_num(g[idx], nan=99))]           # brightest first
    keep = []
    if len(idx):
        xy = np.column_stack([x[idx], y[idx]])
        alive = np.ones(len(idx), bool)
        tree = cKDTree(xy)
        for i in range(len(idx)):
            if alive[i]:
                for j in tree.query_ball_point(xy[i], minSepPx):
                    if j > i:
                        alive[j] = False
        keep = idx[alive]
    keep = np.asarray(keep, np.int64)
    return np.column_stack([x[keep], y[keep]]) if len(keep) else np.zeros((0, 2)), keep


selectSeeds = selectGaiaSources     # the old name, so existing code keeps working


# old snake_case names still work (with a FutureWarning)
from ._compat import moduleGetattr  # noqa: E402
__getattr__ = moduleGetattr(globals())
