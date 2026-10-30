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

"""Whole-detector runs: tile, fit, stitch.

A detector is too large for one joint fit, so it is covered by a grid of tiles whose INTERIORS
partition the detector exactly; each tile is fitted on its interior grown by ``pad`` pixels (so
sources near a tile edge are fitted together with their neighbours from the next tile) and only
its interior is kept:

1. `fitTiles` fits the tiles (all of them, a region, or one share of them for a Slurm array
   task) and writes one ``tile_<name>.fits`` per tile -- data, model, variance, bad-pixel mask
   and the photometry table.  Tiles already on disk are skipped, so a job can be restarted.
2. `stitchTiles` assembles the interiors into full-detector images (``detector_mosaic.fits``:
   DATA, MODEL, RESIDUAL, VARIANCE, BAD, TILE extensions -- opens in DS9) and merges the
   catalogues (``detector_catalog.fits``), keeping every source once: in the tile whose interior
   contains it.
3. `plotDetector` shows data | model | residual | residual/sigma for the whole detector or a
   zoom, with optional tile boundaries.

The image called DATA is the image the fit used (after the pipeline's background subtraction), so
RESIDUAL = DATA - MODEL.
"""

__all__ = ["tileSpecs", "tileName", "saveTileProducts", "fitTiles", "stitchTiles",
           "loadMosaic", "loadTile", "plotTile", "plotDetector"]

import glob
import os
import time

import numpy as np

TILE_PREFIX = "tile_"


# ---------------------------------------------------------------- tiling
def tileSpecs(width, height, tile=400, x0=0, y0=0, region=None):
    """Return the interiors of a tile grid that partitions an area exactly.

    ``ceil(size / tile)`` tiles per axis with sizes as equal as possible.

    Parameters
    ----------
    width, height : `int`
        Size of the area [px].
    tile : `int`, optional
        Target tile size [px].
    x0, y0 : `int`, optional
        Origin of the area.
    region : `tuple` [`int`, `int`, `int`, `int`], optional
        (xmin, xmax, ymin, ymax): keep only tiles overlapping it.

    Returns
    -------
    specs : `list` [`tuple` [`int`, `int`, `int`, `int`]]
        (x0, y0, w, h) per tile.
    """
    nx = max(1, int(np.ceil(width / tile)))
    ny = max(1, int(np.ceil(height / tile)))
    ex = np.round(np.linspace(0, width, nx + 1)).astype(int) + x0
    ey = np.round(np.linspace(0, height, ny + 1)).astype(int) + y0
    out = []
    for j in range(ny):
        for i in range(nx):
            t = (int(ex[i]), int(ey[j]), int(ex[i + 1] - ex[i]), int(ey[j + 1] - ey[j]))
            if region is not None:
                if (t[0] + t[2] - 1 < region[0] or t[0] > region[1]
                        or t[1] + t[3] - 1 < region[2] or t[1] > region[3]):
                    continue
            out.append(t)
    return out


def tileName(spec):
    """Return the name of a tile, e.g. 'x0400_y0000_400x400'.

    Parameters
    ----------
    spec : `tuple` [`int`, `int`, `int`, `int`]
        (x0, y0, w, h).

    Returns
    -------
    name : `str`
        The tile name.
    """
    x0, y0, w, h = spec
    return f"x{x0:04d}_y{y0:04d}_{w}x{h}"


# ---------------------------------------------------------------- per-tile products
def saveTileProducts(run, path, phot=None, meta=None):
    """Write one tile's results to a FITS file.

    DATA, MODEL, VARIANCE and BAD images over the tile's fit box (header keys X0/Y0 = fit-box
    origin, IX0/IY0/IX1/IY1 = interior, inclusive) and a CATALOG table (the photometry, with an
    ``owned`` column: source inside this tile's interior, half-open so boundary sources belong to
    exactly one tile).  The file is written atomically.

    Parameters
    ----------
    run : `lsst.meas.crowd.session.CutoutRun`
        The tile's finished run.
    path : `str`
        Output file.
    phot : `astropy.table.Table`, optional
        Photometry (computed if not given).
    meta : `dict`, optional
        Extra header keywords.

    Returns
    -------
    path : `str`
        The file written.
    """
    from astropy.io import fits
    from astropy.table import Table
    final = run.outerHistory[-1]
    exp, res = final["exposure"], final["result"]
    bb = exp.getBBox()
    ib = run.interiorBbox
    data = np.asarray(exp.image.array, np.float32)
    model = np.asarray(res["modelFull"], np.float32)
    var = np.asarray(exp.variance.array, np.float32)
    bad = (((exp.mask.array & run.badBitmask) != 0) | ~np.isfinite(data)).astype(np.uint8)
    hdr = fits.Header()
    hdr["X0"], hdr["Y0"] = bb.getMinX(), bb.getMinY()
    hdr["IX0"], hdr["IY0"] = ib.getMinX(), ib.getMinY()
    hdr["IX1"], hdr["IY1"] = ib.getMaxX(), ib.getMaxY()
    hdr["TILE"] = run.name
    for k, v in (meta or {}).items():
        hdr[k[:8].upper()] = v
    if phot is None:
        from .photometry import measurePhotometry
        phot = measurePhotometry(run, verbose=False)
    cat = Table(phot, copy=True)
    x, y = np.asarray(cat["x"], float), np.asarray(cat["y"], float)
    # half-open ownership: a source on a tile boundary belongs to exactly one tile
    cat["owned"] = ((x >= ib.getMinX() - 0.5) & (x < ib.getMaxX() + 0.5) &
                    (y >= ib.getMinY() - 0.5) & (y < ib.getMaxY() + 0.5))
    cat["tile"] = run.name
    from .photometry import fitsMeta
    cat.meta = fitsMeta(cat.meta)                         # FITS keywords <= 8 characters
    hdus = [fits.PrimaryHDU(header=hdr)]
    for nm, arr in (("DATA", data), ("MODEL", model), ("VARIANCE", var), ("BAD", bad)):
        hdus.append(fits.ImageHDU(arr, header=hdr, name=nm))
    hdus.append(fits.table_to_hdu(cat))
    hdus[-1].name = "CATALOG"
    tmp = path + ".tmp"
    fits.HDUList(hdus).writeto(tmp, overwrite=True)
    os.replace(tmp, path)                  # atomic: a killed job never leaves half a file
    return path


# ---------------------------------------------------------------- fitting
def _loadVisitImage(visit, detector, repo, collection, instrument, visitImageCollections):
    """Load a visit image through the Butler."""
    from lsst.daf.butler import Butler
    butler = Butler(repo, collections=collection)
    where = f"instrument='{instrument}' AND visit={visit} AND detector={detector}"
    refs = butler.query_datasets("visit_image", collections=visitImageCollections, where=where,
                                 find_first=False)
    if not refs:
        raise RuntimeError(f"no visit_image for visit {visit}, detector {detector}")
    return butler.get("visit_image", collections=visitImageCollections, dataId=refs[0].dataId)


def fitTiles(visit, detector, config=None, outDir="detector_run", tile=400, pad=30,
             taskIndex=0, nTasks=1, region=None, repo="dp2_prep",
             collection="LSSTCam/runs/DRP/DP2", instrument="LSSTCam",
             visitImageCollections="LSSTCam/runs/*", psfSource="visit", overwrite=False,
             visitImage=None, **overrides):
    """Fit this task's share of a detector's tiles and write one file per tile.

    Tiles already on disk are skipped, so the job can be restarted; a failing tile is reported
    and the task moves on.

    Parameters
    ----------
    visit, detector : `int`
        Visit and detector ids.
    config : `dict`, optional
        Pipeline configuration (default ``makeConfig(CROWDED_FIELD_CONFIG)``).
    outDir : `str`, optional
        Output directory.
    tile, pad : `int`, optional
        Tile interior size and fit padding [px].
    taskIndex, nTasks : `int`, optional
        This task fits tiles ``taskIndex``, ``taskIndex + nTasks``, ... (one Slurm array task
        each).
    region : `tuple` [`int`, `int`, `int`, `int`], optional
        (xmin, xmax, ymin, ymax): only tiles overlapping it (for tests).
    repo : `str`, optional
        Butler repository path or alias.
    collection : `str`, optional
        Default Butler collection.
    instrument : `str`, optional
        Instrument name.
    visitImageCollections : `str`, optional
        Collections searched for the visit image.
    psfSource : `str`, optional
        'visit' or 'refit' (see `lsst.meas.crowd.session.buildDetectorPsf`).
    overwrite : `bool`, optional
        Refit tiles that already have a file.
    visitImage : `lsst.afw.image.Exposure`, optional
        An already-loaded visit image (skips the Butler).
    **overrides
        Keyword arguments of ``runSessionGaia`` (e.g. ``maxOuter=3, maxInner=15``).

    Returns
    -------
    paths : `list` [`str`]
        Tile files of this task that exist.
    """
    import contextlib
    import io
    from . import session as cn
    from .gaiaCatalog import attachGaia
    from .photometry import measurePhotometry
    from .pipeline import runSessionGaia
    os.makedirs(outDir, exist_ok=True)
    if visitImage is None:
        visitImage = _loadVisitImage(visit, detector, repo, collection, instrument,
                                     visitImageCollections)
    bb = visitImage.getBBox()
    specs = tileSpecs(bb.getWidth(), bb.getHeight(), tile, bb.getMinX(), bb.getMinY(), region)
    mine = specs[taskIndex::nTasks]
    paths = {tileName(s): os.path.join(outDir, f"{TILE_PREFIX}{tileName(s)}.fits") for s in mine}
    todo = [s for s in mine if overwrite or not os.path.exists(paths[tileName(s)])]
    print(f"detector {detector} of visit {visit}: {len(specs)} tiles of ~{tile} px (pad {pad}); "
          f"task {taskIndex + 1}/{nTasks} has {len(mine)}, {len(todo)} still to fit")
    if not todo:
        return sorted(paths.values())
    cfg = config if config is not None else cn.makeConfig(cn.CROWDED_FIELD_CONFIG)
    specsTodo = [dict(x0=s[0], y0=s[1], w=s[2], h=s[3], name=tileName(s)) for s in todo]
    sess = cn.prepareSession(visit, detector, specsTodo, pad=pad, psfSource=psfSource,
                             visitImage=visitImage)
    attachGaia(sess, cacheDir=os.path.join(outDir, "gaia_cache"),
               **({"repo": repo} if repo else {}))
    kw = dict(maxOuter=3, maxInner=15, psfPolishCycles=3)
    kw.update(overrides)
    tAll = time.perf_counter()
    for k, s in enumerate(todo):
        name = tileName(s)
        t0 = time.perf_counter()
        print(f"\n######## tile {k + 1}/{len(todo)}: {name}")
        try:
            runs = runSessionGaia(sess, cfg, names=[name], plotFinal=False, **kw)
            run = runs[name]
            with contextlib.redirect_stdout(io.StringIO()):
                phot = (run.photometry if getattr(run, "photometry", None) is not None
                        else measurePhotometry(run, verbose=False, gaia=getattr(sess, "gaia", None)))
            saveTileProducts(run, paths[name], phot,
                             meta=dict(visit=visit, detector=detector,
                                       seconds=time.perf_counter() - t0))
            nIn = int(np.sum(phot["inInterior"])) if "inInterior" in phot.colnames else len(phot)
            print(f"######## tile {name}: {nIn} sources in its interior, "
                  f"{time.perf_counter() - t0:.0f} s -> {paths[name]}")
        except Exception as e:                       # keep going: one bad tile must not kill the job
            import traceback
            traceback.print_exc()
            print(f"######## tile {name} FAILED ({type(e).__name__}: {e}) -- rerun to retry")
        finally:
            sess.jobs.pop(name, None)                # free the cutout
    print(f"\ntask {taskIndex + 1}/{nTasks}: {len(todo)} tile(s) in "
          f"{time.perf_counter() - tAll:.0f} s")
    return sorted(p for p in paths.values() if os.path.exists(p))


# ---------------------------------------------------------------- stitching
def stitchTiles(outDir="detector_run", mosaicName="detector_mosaic.fits",
                catalogName="detector_catalog.fits", verbose=True):
    """Assemble the tile files into full-detector images and one catalogue.

    Writes the mosaic (DATA, MODEL, RESIDUAL, VARIANCE, BAD and TILE extensions) and the
    catalogue (every source once, from the tile whose interior contains it) into ``outDir``.

    Parameters
    ----------
    outDir : `str`, optional
        Directory with the tile files.
    mosaicName : `str`, optional
        Mosaic file name.
    catalogName : `str`, optional
        Catalogue file name.
    verbose : `bool`, optional
        Print a summary.

    Returns
    -------
    mosaic : `dict`
        ``data``, ``model``, ``resid``, ``var``, ``bad``, ``tile`` (tile index map, -1 = not
        fitted), ``x0``, ``y0``, ``tiles`` (names), ``interiors`` and ``catalog``.

    Raises
    ------
    FileNotFoundError
        Raised if ``outDir`` has no tile files.
    """
    from astropy.io import fits
    from astropy.table import Table, vstack
    files = sorted(glob.glob(os.path.join(outDir, f"{TILE_PREFIX}*.fits")))
    if not files:
        raise FileNotFoundError(f"no {TILE_PREFIX}*.fits in {outDir}")
    heads = [fits.getheader(f, "DATA") for f in files]
    gx0 = min(h["IX0"] for h in heads)
    gy0 = min(h["IY0"] for h in heads)
    gx1 = max(h["IX1"] for h in heads)
    gy1 = max(h["IY1"] for h in heads)
    shape = (gy1 - gy0 + 1, gx1 - gx0 + 1)
    out = {k: np.full(shape, np.nan, np.float32) for k in ("data", "model", "var")}
    out["bad"] = np.zeros(shape, np.uint8)
    out["tile"] = np.full(shape, -1, np.int32)
    cats, names, interiors = [], [], []
    for t, (f, h) in enumerate(zip(files, heads)):
        with fits.open(f) as hl:
            fx0, fy0 = h["X0"], h["Y0"]
            ix0, iy0, ix1, iy1 = h["IX0"], h["IY0"], h["IX1"], h["IY1"]
            src = (slice(iy0 - fy0, iy1 - fy0 + 1), slice(ix0 - fx0, ix1 - fx0 + 1))
            dst = (slice(iy0 - gy0, iy1 - gy0 + 1), slice(ix0 - gx0, ix1 - gx0 + 1))
            out["data"][dst] = hl["DATA"].data[src]
            out["model"][dst] = hl["MODEL"].data[src]
            out["var"][dst] = hl["VARIANCE"].data[src]
            out["bad"][dst] = hl["BAD"].data[src]
            out["tile"][dst] = t
            cat = Table(hl["CATALOG"].data)
        names.append(h.get("TILE", os.path.basename(f)))
        interiors.append((ix0, iy0, ix1, iy1))
        cats.append(cat[np.asarray(cat["owned"], bool)])
    out["resid"] = out["data"] - out["model"]
    catalog = vstack(cats, metadata_conflicts="silent") if cats else Table()
    if len(catalog):
        catalog["id"] = np.arange(len(catalog))
    out.update(x0=gx0, y0=gy0, tiles=names, interiors=interiors, catalog=catalog)
    hdr = fits.Header()
    hdr["X0"], hdr["Y0"] = gx0, gy0
    hdr["NTILES"] = len(files)
    hdus = [fits.PrimaryHDU(header=hdr)]
    for nm, key in (("DATA", "data"), ("MODEL", "model"), ("RESIDUAL", "resid"),
                    ("VARIANCE", "var"), ("BAD", "bad"), ("TILE", "tile")):
        hdus.append(fits.ImageHDU(out[key], header=hdr, name=nm))
    fits.HDUList(hdus).writeto(os.path.join(outDir, mosaicName), overwrite=True)
    catalog.write(os.path.join(outDir, catalogName), overwrite=True)
    if verbose:
        nMissing = int((out["tile"] < 0).sum())
        print(f"stitched {len(files)} tile(s) -> {shape[1]} x {shape[0]} px mosaic, "
              f"{len(catalog)} sources"
              + (f"; {nMissing} px not covered by any finished tile" if nMissing else "")
              + f"\n  {os.path.join(outDir, mosaicName)}\n  {os.path.join(outDir, catalogName)}")
    return out


def loadMosaic(outDir="detector_run", mosaicName="detector_mosaic.fits",
               catalogName="detector_catalog.fits"):
    """Read the files written by `stitchTiles`.

    Parameters
    ----------
    outDir : `str`, optional
        Directory.
    mosaicName : `str`, optional
        Mosaic file name.
    catalogName : `str`, optional
        Catalogue file name.

    Returns
    -------
    mosaic : `dict`
        As returned by `stitchTiles` (without ``tiles`` and ``interiors``).
    """
    from astropy.io import fits
    from astropy.table import Table
    with fits.open(os.path.join(outDir, mosaicName)) as hl:
        out = dict(data=hl["DATA"].data, model=hl["MODEL"].data, resid=hl["RESIDUAL"].data,
                   var=hl["VARIANCE"].data, bad=hl["BAD"].data, tile=hl["TILE"].data,
                   x0=hl[0].header["X0"], y0=hl[0].header["Y0"])
    p = os.path.join(outDir, catalogName)
    out["catalog"] = Table.read(p) if os.path.exists(p) else None
    return out


def loadTile(path):
    """Read one tile file written by `saveTileProducts`.

    Parameters
    ----------
    path : `str`
        The ``tile_<name>.fits`` file.

    Returns
    -------
    tile : `dict`
        In the form of `loadMosaic`'s result over the tile's fit box (``data``, ``model``,
        ``resid``, ``var``, ``bad``, ``tile``, ``x0``, ``y0``, ``catalog``), plus ``name`` and
        ``interiors`` (the tile's interior as one (x0, y0, x1, y1) tuple, inclusive).
    """
    from astropy.io import fits
    from astropy.table import Table
    with fits.open(path) as hl:
        h = hl[0].header
        data = hl["DATA"].data.astype(float)
        model = hl["MODEL"].data.astype(float)
        out = dict(data=data, model=model, resid=data - model, var=hl["VARIANCE"].data.astype(float),
                   bad=hl["BAD"].data, tile=np.zeros(data.shape, np.int32), x0=h["X0"], y0=h["Y0"],
                   name=h.get("TILE", os.path.basename(path)),
                   interiors=[(h["IX0"], h["IY0"], h["IX1"], h["IY1"])])
    out["catalog"] = Table.read(path, hdu="CATALOG")
    return out


def plotTile(tile, zClip=5.0, grayscale=False, sources=True, figsize=(24, 6.5), out=None):
    """Plot data | model | residual | residual/sigma of one tile over its whole fit box.

    The tile's interior (the part kept in the mosaic) is outlined; outside it is the padding,
    fitted only so that sources near the edge have their neighbours.

    Parameters
    ----------
    tile : `dict` or `str`
        Output of `loadTile`, or the tile file.
    zClip : `float`, optional
        Colour range of residual/sigma.
    grayscale : `bool`, optional
        Residual panels in greyscale.
    sources : `bool`, optional
        Circle the tile catalogue's sources.
    figsize : `tuple` [`float`, `float`], optional
        Figure size [inches].
    out : `str`, optional
        Also save the figure to this file.

    Returns
    -------
    result : `dict`
        As `plotDetector`'s.
    """
    if isinstance(tile, str):
        tile = loadTile(tile)
    return plotDetector(tile, binning=1, showTiles=True, zClip=zClip, grayscale=grayscale,
                        figsize=figsize, sources=sources, out=out)


# ---------------------------------------------------------------- plotting
def _block(a, b, func=np.nanmean):
    """Reduce an image by b x b blocks (edges trimmed), for display."""
    if b <= 1:
        return a
    ny, nx = (a.shape[0] // b) * b, (a.shape[1] // b) * b
    return func(a[:ny, :nx].reshape(ny // b, b, nx // b, b), axis=(1, 3))


def plotDetector(mosaic, view=None, binning=None, showTiles=True, zClip=5.0, grayscale=False,
                 figsize=(26, 7.5), sources=False, out=None):
    """Plot data | model | residual | residual/sigma of a stitched detector or a zoom.

    Parameters
    ----------
    mosaic : `dict` or `str`
        Output of `stitchTiles` / `loadMosaic`, or the directory to load.
    view : `tuple` [`float`, `float`, `float`, `float`], optional
        (xmin, xmax, ymin, ymax) zoom in detector pixels.
    binning : `int`, optional
        Display block size (default: about 1000 px across the view); residual/sigma is binned as
        sum / sqrt(n) so its noise stays 1.
    showTiles : `bool`, optional
        Draw the tile boundaries.
    zClip : `float`, optional
        Colour range of residual/sigma.
    grayscale : `bool`, optional
        Residual panels in greyscale.
    figsize : `tuple` [`float`, `float`], optional
        Figure size [inches].
    sources : `bool`, optional
        Circle the catalogue's sources (for zooms).
    out : `str`, optional
        Also save the figure to this file.

    Returns
    -------
    result : `dict`
        ``residualSigma`` (robust sigma of residual/sigma over the view) and ``binning``.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from astropy.visualization import simple_norm
    if isinstance(mosaic, str):
        mosaic = loadMosaic(mosaic)
    gx0, gy0 = mosaic["x0"], mosaic["y0"]
    ny, nx = mosaic["data"].shape
    if view is None:
        view = (gx0, gx0 + nx - 1, gy0, gy0 + ny - 1)
    xa, xb = max(int(view[0]) - gx0, 0), min(int(view[1]) - gx0 + 1, nx)
    ya, yb = max(int(view[2]) - gy0, 0), min(int(view[3]) - gy0 + 1, ny)
    sl = (slice(ya, yb), slice(xa, xb))
    b = binning or max(1, int(round(max(xb - xa, yb - ya) / 1000)))
    bad = (mosaic["bad"][sl] > 0) | (mosaic["tile"][sl] < 0) | ~np.isfinite(mosaic["data"][sl])
    d = np.where(bad, np.nan, mosaic["data"][sl])
    m = np.where(bad, np.nan, mosaic["model"][sl])
    r = d - m
    v = np.where(bad | ~(mosaic["var"][sl] > 0), np.nan, mosaic["var"][sl])
    z = r / np.sqrt(v)
    good = np.isfinite(z)
    zs = np.where(good, z, 0.0)
    nb = _block(good.astype(float), b, np.sum)
    zb = np.where(nb > 0, _block(zs, b, np.sum) / np.sqrt(np.maximum(nb, 1)), np.nan) if b > 1 else z
    db, mb, rb = _block(d, b), _block(m, b), _block(r, b)
    badb = _block(bad.astype(float), b) > 0.5
    ext = (gx0 + xa - .5, gx0 + xa + db.shape[1] * b - .5,
           gy0 + ya - .5, gy0 + ya + db.shape[0] * b - .5)
    fin = db[np.isfinite(db)]
    norm = simple_norm(fin, vmin=np.nanpercentile(fin, 1), vmax=np.nanpercentile(fin, 99.5),
                       stretch="asinh") if fin.size else None
    rs = np.nanpercentile(np.abs(rb), 99) if np.isfinite(rb).any() else 1.0
    zg = z[good]
    q1, q3 = (np.percentile(zg, [25, 75]) if zg.size else (np.nan, np.nan))
    sig = (q3 - q1) / 1.349
    cmRes = plt.get_cmap("gray" if grayscale else "RdBu_r")
    panels = [(f"data (fitted image){', binned %dx%d' % (b, b) if b > 1 else ''}", db,
               dict(cmap="viridis", norm=norm)),
              ("model", mb, dict(cmap="viridis", norm=norm)),
              ("residual = data − model", rb, dict(cmap=cmRes, vmin=-rs, vmax=rs)),
              (f"residual / σ  (robust σ {sig:.3f} per px)", zb,
               dict(cmap=cmRes, vmin=-zClip, vmax=zClip))]
    fig, axes = plt.subplots(1, 4, figsize=figsize, sharex=True, sharey=True)
    for ax, (ttl, img, kw) in zip(axes, panels):
        im = ax.imshow(img, origin="lower", extent=ext, interpolation="nearest", **kw)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        if badb.any():
            ax.imshow(np.ma.masked_where(~badb, np.ones(badb.shape)), origin="lower", extent=ext,
                      cmap=ListedColormap(["#ff00ff"]), interpolation="nearest")
        if showTiles and "interiors" in mosaic:
            for ix0, iy0, ix1, iy1 in mosaic["interiors"]:
                ax.add_patch(plt.Rectangle((ix0 - .5, iy0 - .5), ix1 - ix0 + 1, iy1 - iy0 + 1,
                                           fill=False, ec="w" if ax in axes[:2] else "k", lw=0.4,
                                           alpha=0.5))
        if sources and mosaic.get("catalog") is not None and len(mosaic["catalog"]):
            c = mosaic["catalog"]
            ax.scatter(c["x"], c["y"], s=12, facecolors="none",
                       edgecolors="w" if ax in axes[:2] else "k", linewidths=0.4)
        ax.set_title(ttl, fontsize=11)
        ax.set_xlabel("x [px]")
    axes[0].set_ylabel("y [px]")
    axes[0].set_xlim(view[0] - .5, view[1] + .5)
    axes[0].set_ylim(view[2] - .5, view[3] + .5)
    nSrc = len(mosaic["catalog"]) if mosaic.get("catalog") is not None else 0
    what = f"tile {mosaic['name']}" if "name" in mosaic else "stitched detector"
    fig.suptitle(f"{what}: {nSrc} sources; magenta = bad / not fitted"
                 + ("; lines = tile interiors" if showTiles else ""), fontsize=11)
    plt.tight_layout()
    if out:
        fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.show()
    return dict(residualSigma=float(sig), binning=b)


# old snake_case names still work (with a FutureWarning)
from ._compat import moduleGetattr  # noqa: E402
__getattr__ = moduleGetattr(globals())
