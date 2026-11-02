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


"""Write every diagnostic plot of a run (or of a detector) to PNG files.

`saveRunPlots` saves, for one fitted cutout or tile, the PSF stars and their stamps, the PSF
model of every pass, the background of every pass, data | model | residual per outer and inner
round, what changed between passes, the final comparison, per-source pages and cutouts, the
bright-star diagnostics and the Gaia comparison, numbered in the order of the pipeline, plus an
``index.txt`` that says what each file shows.  `saveDetectorPlots` does the same for the
detector-level PSF (the stars the visit PSF was built from).  `fitTiles(..., plots=...)` and
``crowdDetector.py fit --plots`` call them.
"""

__all__ = ["saveRunPlots", "saveDetectorPlots", "showPsfStarStamps", "sampleTileNames"]

import contextlib
import io
import os
import traceback

import numpy as np

# sample tiles: the ones containing these points (fractions of the detector width / height):
# the four corners, the centre and the middle of the bottom and left edges
SAMPLE_POINTS = ((0.02, 0.02), (0.98, 0.02), (0.02, 0.98), (0.98, 0.98), (0.5, 0.5),
                 (0.5, 0.02), (0.02, 0.5))


def sampleTileNames(specs, width, height, x0=0, y0=0, points=SAMPLE_POINTS):
    """Return the names of the tiles that contain the sample points.

    Parameters
    ----------
    specs : `list` [`tuple` [`int`, `int`, `int`, `int`]]
        Tiles as (x0, y0, width, height), from `lsst.meas.crowd.detector.tileSpecs`.
    width, height : `int`
        Detector size [px].
    x0, y0 : `int`, optional
        Detector origin [px].
    points : `tuple` [`tuple` [`float`, `float`]], optional
        Sample points as fractions of the detector width and height.

    Returns
    -------
    names : `list` [`str`]
        Tile names (`lsst.meas.crowd.detector.tileName`), each once, in the order of ``points``.
    """
    from .detector import tileName
    names = []
    for fx, fy in points:
        x, y = x0 + fx * (width - 1), y0 + fy * (height - 1)
        for s in specs:
            if s[0] <= x < s[0] + s[2] and s[1] <= y < s[1] + s[3]:
                if tileName(s) not in names:
                    names.append(tileName(s))
                break
    return names


def showPsfStarStamps(image, xy, half=15, maxStamps=25, psf=None, title="", out=None):
    """Plot a stamp around each PSF star, followed by the PSF kernel.

    Parameters
    ----------
    image : `lsst.afw.image.Exposure`
        Image to cut the stamps from (detector pixel coordinates).
    xy : `numpy.ndarray`
        Star positions [px], shape (N, 2); the first ``maxStamps`` are shown.
    half : `int`, optional
        Half-size of each stamp [px].
    maxStamps : `int`, optional
        Maximum number of stamps.
    psf : `lsst.afw.detection.Psf`, optional
        PSF whose kernel is shown in the last panel (default: the image's own PSF).
    title : `str`, optional
        Figure title.
    out : `str`, optional
        Also save the figure to this file.

    Returns
    -------
    nShown : `int`
        Number of stamps drawn.
    """
    import matplotlib.pyplot as plt
    from astropy.visualization import simple_norm
    from . import fitLib as L
    data = np.asarray(image.image.array, float)
    bb = image.getBBox()
    x0, y0 = bb.getMinX(), bb.getMinY()
    ny, nx = data.shape
    stamps, where = [], []
    for x, y in np.asarray(xy, float).reshape(-1, 2):
        ix, iy = int(round(x)) - x0, int(round(y)) - y0
        if half <= ix < nx - half and half <= iy < ny - half:
            stamps.append(data[iy - half:iy + half + 1, ix - half:ix + half + 1])
            where.append((x, y))
        if len(stamps) >= maxStamps:
            break
    psf = psf or image.getPsf()
    c = bb.getCenter()
    kernel = np.asarray(psf.computeKernelImage(c).array, float)
    nPanel = len(stamps) + 1
    nCol = min(nPanel, 6)
    nRow = int(np.ceil(nPanel / nCol))
    fig, axes = plt.subplots(nRow, nCol, figsize=(2.6 * nCol, 2.6 * nRow + 0.6), squeeze=False)
    for ax in axes.ravel():
        ax.set_axis_off()
    for k, (st, (x, y)) in enumerate(zip(stamps, where)):
        ax = axes.ravel()[k]
        fin = st[np.isfinite(st)]
        norm = (simple_norm(fin, stretch="asinh", vmin=np.percentile(fin, 1), vmax=fin.max())
                if fin.size and fin.max() > np.percentile(fin, 1) else None)
        ax.imshow(st, origin="lower", cmap="gray_r", norm=norm, interpolation="nearest")
        ax.set_title(f"({x:.0f}, {y:.0f})", fontsize=8)
    ax = axes.ravel()[len(stamps)]
    pos = kernel[kernel > 0]
    ax.imshow(kernel, origin="lower", cmap="magma",
              norm=simple_norm(kernel, stretch="log", vmin=pos.max() * 1e-4 if pos.size else None))
    ax.set_title(f"PSF kernel, FWHM {L._fwhm(psf, c):.2f} px", fontsize=8)
    fig.suptitle(title or f"{len(stamps)} PSF stars", fontsize=11)
    plt.tight_layout()
    if out:
        fig.savefig(out, dpi=110)
    return len(stamps)


class _Saver:
    """Run plot functions and save every figure they open, numbered in call order."""

    def __init__(self, outDir, dpi=100, verbose=True):
        os.makedirs(outDir, exist_ok=True)
        self.outDir, self.dpi, self.verbose = outDir, dpi, verbose
        self.n, self.files, self.index, self.failed = 0, [], [], []
        self._listed = False          # set by a plot function that lists its own files

    def __call__(self, name, what, fn, *args, **kwargs):
        """Call ``fn(*args, **kwargs)``; save each new figure as ``NN_<name>[_k].png``."""
        import matplotlib.pyplot as plt
        before = set(plt.get_fignums())
        result = None
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                result = fn(*args, **kwargs)
        except Exception as e:                    # one failing plot must not stop the others
            self.failed.append(name)
            self.index.append(f"(not written) {name}: {type(e).__name__}: {e}")
            if self.verbose:
                print(f"    plot {name} skipped: {type(e).__name__}: {e}")
                traceback.print_exc(limit=2)
        new = [n for n in plt.get_fignums() if n not in before]
        if not new and name not in self.failed and not self._listed:
            self.index.append(f"(no figure) {name}: {what} -- nothing to show for this run")
        self._listed = False
        for k, num in enumerate(new):
            self.n += 1
            suffix = f"_{k + 1}" if len(new) > 1 else ""
            path = os.path.join(self.outDir, f"{self.n:02d}_{name}{suffix}.png")
            plt.figure(num).savefig(path, dpi=self.dpi, bbox_inches="tight")
            plt.close(num)
            self.files.append(path)
            self.index.append(f"{os.path.basename(path)}: {what}")
        return result

    def text(self, name, what, fn, *args, **kwargs):
        """Write what ``fn`` prints to ``<name>.txt``."""
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                fn(*args, **kwargs)
        except Exception as e:
            buf.write(f"\n{type(e).__name__}: {e}\n")
        path = os.path.join(self.outDir, f"{name}.txt")
        with open(path, "w") as fh:
            fh.write(buf.getvalue())
        self.files.append(path)
        self.index.append(f"{os.path.basename(path)}: {what}")

    def finish(self, title):
        """Write ``index.txt`` and return the files written."""
        with open(os.path.join(self.outDir, "index.txt"), "w") as fh:
            fh.write(title + "\n\n" + "\n".join(self.index) + "\n")
        if self.verbose:
            print(f"    {len(self.files)} plot file(s) in {self.outDir}"
                  + (f" ({len(self.failed)} skipped: {', '.join(self.failed)})" if self.failed else ""))
        return list(self.files)


def _magnitudeSpread(phot, n=12):
    """Source indices spread evenly from the brightest to the faintest interior source."""
    mag = np.asarray(phot["psfMag"], float)
    ok = np.isfinite(mag)
    if "inInterior" in phot.colnames:
        ok &= np.asarray(phot["inInterior"], bool)
    ids = np.asarray(phot["id"] if "id" in phot.colnames else np.arange(len(phot)))[ok]
    ids = ids[np.argsort(mag[ok])]
    if not len(ids):
        return []
    return sorted({int(i) for i in ids[np.linspace(0, len(ids) - 1, min(n, len(ids))).round().astype(int)]},
                  key=lambda i: list(ids).index(i))


def saveRunPlots(run, outDir, phot=None, initialBackground=None, sourcePages=4,
                 allSourcePages=False, dpi=100, verbose=True):
    """Write every diagnostic plot of one finished run to PNG files.

    In pipeline order: PSF stars and their stamps for each PSF update, the PSF model and FWHM of
    each pass, the background of each pass, data | model | residual per outer pass and per inner
    iteration, the sources gained and lost between passes, the final comparison, per-source
    pages (brightest first, and spread in magnitude), labelled source cutouts, the bright-star
    diagnostics and the Gaia comparison.  ``index.txt`` describes every file and ``report.txt``
    holds the pass-by-pass report.  A plot that fails is listed in ``index.txt`` and skipped.

    Parameters
    ----------
    run : `lsst.meas.crowd.session.CutoutRun`
        A finished run.
    outDir : `str`
        Directory for the files (created if needed).
    phot : `astropy.table.Table`, optional
        Photometry of the run (default ``run.photometry`` if present); needed for the
        magnitude-labelled plots and the Gaia comparison.
    initialBackground : `lsst.afw.image.Image`, optional
        Background subtracted before pass 1 (``session.bgImage``).
    sourcePages : `int`, optional
        Number of per-source pages, brightest sources first.
    allSourcePages : `bool`, optional
        Write a page for every source instead (can be hundreds of files).
    dpi : `int`, optional
        Resolution of the PNG files.
    verbose : `bool`, optional
        Print where the files went and which plots were skipped.

    Returns
    -------
    files : `list` [`str`]
        The files written.
    """
    import matplotlib.pyplot as plt
    from . import fitLib as L
    from . import plotting as pt
    from . import sourcePlots as spl
    from . import stagePlots as sp
    from .pipeline import report
    phot = phot if phot is not None else getattr(run, "photometry", None)
    hist = run.outerHistory
    n = len(hist)
    save = _Saver(outDir, dpi=dpi, verbose=verbose)
    plt.ioff()

    # 1. PSF: the stars each residual-stack update used, their stamps, and the PSF of each pass
    for k in range(max(n - 1, 1)):
        nxt = hist[min(k + 1, n - 1)]["psf"]
        label = f"after pass {k + 1}" if n > 1 else "of the final pass"
        stars = save(f"psfStars_pass{k + 1}_positions",
                     f"PSF candidates (orange) and stars used (green) by the PSF update {label}, "
                     "with the resulting kernel and its radial profile",
                     _psfStars, run, k, nxt)
        if stars is not None and len(stars["used"]):
            save(f"psfStars_pass{k + 1}_stamps", f"stamps of the stars used by the PSF update {label}",
                 showPsfStarStamps, hist[k]["exposure"], stars["used"], psf=nxt,
                 title=f"PSF stars used by the update {label} ({len(stars['used'])})")
    save("psfModels", "PSF kernel used in each outer pass", L.plotPsfModelsByOuterPass, hist)
    save("psfFwhm", "PSF FWHM used in each outer pass", L.plotPsfFwhmByOuterPass, hist)

    # 2. background
    save("backgrounds", "background entering each pass (top) and the change from the previous "
         "pass (bottom)", sp.showBackgrounds, run, initial=initialBackground)

    # 3. rounds
    save("rounds_outer", "data | model | residual | residual/sigma after each outer pass",
         sp.showRounds, run, level="outer")
    for k in range(n):
        save(f"rounds_inner_pass{k + 1}", f"data | model | residual | residual/sigma for each "
             f"inner iteration of pass {k + 1}", sp.showRounds, run, level="inner", outer=k)
    save("outerLoopSummary", "residual std and number of sources per outer pass",
         L.plotOuterLoopSummary, hist)
    for k in range(n - 1):
        save(f"comparePasses_{k + 1}to{k + 2}", f"sources lost and gained from pass {k + 1} to "
             f"{k + 2}, and why", sp.comparePasses, run, a=k, b=k + 1)

    # 4. final result
    save("final_comparison", "final data | model | residual | residual/sigma (tile interior)",
         pt.plotFinalComparison, run)
    save("final_comparison_full", "final data | model | residual | residual/sigma (whole fit "
         "box, with the padding)", pt.plotFinalComparison, run, view=None)
    save("final_residual_gray", "final residual in grayscale", pt.plotResidualGray, run)

    # 5. individual sources
    pagesDir = os.path.join(outDir, "source_pages")
    save("sourcePages_brightest", "per-source pages, brightest first: data | model | residual | "
         f"cross-sections (one overview + page per file in {os.path.basename(pagesDir)}/)",
         _linkedPages, run, pagesDir, None if allSourcePages else sourcePages, save)
    if phot is not None:
        pick = _magnitudeSpread(phot)
        if pick:
            save("sourcePages_magnitudeSpread", "per-source pages for 12 sources spread from the "
                 "brightest to the faintest", sp.showLinkedSourcePages, run, sources=pick,
                 sourcesPerPage=6, maxPages=None)
        for how in ("brightest", "worst", "faintest"):
            save(f"sourceCutouts_{how}", f"magnitude-labelled cutouts of the {how} sources",
                 spl.plotSourceCutouts, run, how=how, n=8, phot=phot)
    save("starZooms", "zooms on stars: data | model | residual", pt.showStarZooms, run)
    save("crossSections", "cross-sections through stars: data and model",
         pt.showCrossSections, run)
    save("brightStarStack", "stacked residual of the bright stars", pt.showBrightStarStack, run)
    save("coreExcess", "core-excess diagnostic of the bright stars", pt.coreExcessDiagnostic, run)
    save("starResidualTable", "residual flux around each star", pt.starResidualTable, run)

    # 6. photometry
    if phot is not None and "gaiaMatch" in phot.colnames:
        save("gaiaComparison", "PSF magnitudes against Gaia", pt.plotGaiaComparison, phot)
    save.text("report", "what each outer and inner pass did", report, run)
    return save.finish(f"Plots of {getattr(run, 'name', 'run')}")


def _psfStars(run, k, psf):
    """showPsfStars for the PSF update after pass ``k``; returns the stars."""
    from .stagePlots import psfStarPositions, showPsfStars
    stars = psfStarPositions(run=run, outer=k)
    showPsfStars(run.outerHistory[k]["exposure"], stars, psf=psf)
    return stars


def _linkedPages(run, pagesDir, maxPages, save):
    """Per-source pages into their own folder (one PNG each), listed in the index."""
    from .stagePlots import showLinkedSourcePages
    out = showLinkedSourcePages(run, sourceView="all", sourcesPerPage=6, maxPages=maxPages,
                                outDir=pagesDir)
    files = out.get("files", [])
    save.files.extend(files)
    rel = os.path.basename(pagesDir)
    save.index.extend(f"{rel}/{os.path.basename(f)}: "
                      + ("overview of the cutout with one box per source on the matching page"
                         if os.path.basename(f).startswith("overview") else
                         "one row per source: data | model | residual | x and y cross-sections")
                      for f in files)
    save._listed = True
    return out


def saveDetectorPlots(visitImage, outDir, catalog=None, session=None, dpi=100, verbose=True):
    """Write the detector-level PSF plots: the PSF stars over the detector, their stamps, the PSF.

    Parameters
    ----------
    visitImage : `lsst.afw.image.Exposure`
        The visit image.
    outDir : `str`
        Directory for the files (created if needed).
    catalog : `lsst.afw.table.SourceCatalog` or `astropy.table.Table`, optional
        The visit's source catalogue with ``calib_psf_candidate`` / ``calib_psf_used`` (the stars
        of the visit PSF).
    session : `lsst.meas.crowd.session.Session`, optional
        A session prepared with ``psfSource='refit'`` (the stars of the refitted PSF), used when
        no catalogue is given.
    dpi : `int`, optional
        Resolution of the PNG files.
    verbose : `bool`, optional
        Print where the files went.

    Returns
    -------
    files : `list` [`str`]
        The files written.
    """
    from .stagePlots import psfStarPositions, showPsfStars
    save = _Saver(outDir, dpi=dpi, verbose=verbose)
    stars = None
    if catalog is not None or (session is not None and session.starCat is not None):
        stars = save("psfStars_detector", "detector-level PSF candidates (orange) and stars used "
                     "(green) over the whole detector, with the kernel and its profile",
                     _detectorStars, visitImage, catalog, session, psfStarPositions, showPsfStars)
    else:
        save("psfKernel_detector", "the visit PSF kernel and its profile", showPsfStars,
             visitImage, None)
    if stars is not None:
        xy = stars["used"] if len(stars["used"]) else stars["candidates"]
        if len(xy):
            save("psfStars_detector_stamps", "stamps of the detector-level PSF stars",
                 showPsfStarStamps, visitImage, xy, half=20, maxStamps=30,
                 psf=session.psf if session is not None and catalog is None else None,
                 title=f"detector-level PSF stars ({stars['source']})")
    bb = visitImage.getBBox()
    for fx, fy in ((0.1, 0.1), (0.5, 0.5), (0.9, 0.9)):
        at = (bb.getMinX() + fx * bb.getWidth(), bb.getMinY() + fy * bb.getHeight())
        save(f"psfKernel_at_{at[0]:.0f}_{at[1]:.0f}", f"visit PSF kernel at ({at[0]:.0f}, {at[1]:.0f})",
             showPsfStars, visitImage, None, at=at)
    return save.finish("Detector-level PSF plots")


def _detectorStars(visitImage, catalog, session, psfStarPositions, showPsfStars):
    stars = (psfStarPositions(catalog=catalog) if catalog is not None
             else psfStarPositions(session=session))
    showPsfStars(visitImage, stars, psf=None if catalog is not None else session.psf)
    return stars
