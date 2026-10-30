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

"""Crowded-field iterative joint PSF fitting: the core library.

The joint non-negative PSF fit, source detection and filtering, centroid refinement, PSF and
background re-determination, the iterative inner / outer loop of `runCutoutPipeline`, and
the plots used to inspect each step.  The Gaia-started adaptive pipeline in
`lsst.meas.crowd.pipeline` builds on these functions.
"""
import math
import time
import warnings
from collections import Counter, OrderedDict, defaultdict
from contextlib import contextmanager

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as patches
import matplotlib.patheffects as pe
from matplotlib.colors import ListedColormap, TwoSlopeNorm
from matplotlib.backends.backend_pdf import PdfPages

from astropy.table import Table
from astropy.visualization import simple_norm, ZScaleInterval, ImageNormalize

from scipy.sparse import lil_matrix, csr_matrix, coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.optimize import lsq_linear, curve_fit
from scipy.spatial import cKDTree
from scipy.ndimage import label, binary_dilation, center_of_mass, map_coordinates
from scipy.signal import fftconvolve

import lsst.geom as geom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.afw.detection as afwDetect
import lsst.afw.geom as afwGeom

from lsst.meas.algorithms import (
    SubtractBackgroundTask, SubtractBackgroundConfig,
    SourceDetectionTask, SourceDetectionConfig,
    ObjectSizeStarSelectorTask, ObjectSizeStarSelectorConfig,
    KernelPsf,
)
from lsst.meas.algorithms.subtractBackground import TooManyMaskedPixelsError
from lsst.meas.algorithms.installGaussianPsf import InstallGaussianPsfTask, InstallGaussianPsfConfig
from lsst.meas.algorithms.makePsfCandidates import MakePsfCandidatesTask, MakePsfCandidatesConfig
from lsst.meas.algorithms.pcaPsfDeterminer import PcaPsfDeterminerTask, PcaPsfDeterminerConfig
from lsst.meas.deblender import SourceDeblendTask, SourceDeblendConfig
from lsst.meas.base import SingleFrameMeasurementTask, SingleFrameMeasurementConfig

warnings.filterwarnings("ignore")

from lsst.utils.threads import disable_implicit_threading
disable_implicit_threading()



# =====================================================================
# Inner-loop helpers: joint fit, detection, filtering, centroids
# =====================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy.visualization import simple_norm
from astropy.table import Table
import lsst.geom as geom
import lsst.afw.table as afwTable
from lsst.meas.algorithms import SourceDetectionTask
from lsst.meas.algorithms.subtractBackground import TooManyMaskedPixelsError
from scipy.sparse import lil_matrix, csr_matrix
from scipy.optimize import lsq_linear
from scipy.spatial import cKDTree
from scipy.ndimage import label, binary_dilation, center_of_mass


# =====================================================================
# Helpers
# =====================================================================

def bboxToSlices(bbox, parentBbox):
    """Convert a box into array slices relative to a parent box.

    Parameters
    ----------
    bbox : `lsst.geom.Box2I`
        Box to convert.
    parentBbox : `lsst.geom.Box2I`
        Box of the array the slices index (its origin is pixel [0, 0]).

    Returns
    -------
    slices : `tuple` [`slice`, `slice`]
        (y slice, x slice) selecting ``bbox`` from an array covering ``parentBbox``.
    """
    x0 = bbox.getMinX() - parentBbox.getMinX()
    y0 = bbox.getMinY() - parentBbox.getMinY()
    x1 = x0 + bbox.getWidth()
    y1 = y0 + bbox.getHeight()
    return slice(y0, y1), slice(x0, x1)


def collapseSaturatedPeaks(peakPositions, exposure, satBitmask=None, dilateIters=2):
    """Collapse every peak on or near a saturated blob into one peak per blob.

    The representative position is a flux-weighted centroid from the unsaturated ring just
    outside the saturated core (the flat top carries no usable position or flux information).
    Without this, noise on the saturated flat top produces many near-duplicate peaks whose PSF
    stamps are almost identical; the joint fit then splits the true amplitude between them
    arbitrarily, which shows up as a blocky, high-frequency patch in the model image.

    Parameters
    ----------
    peakPositions : `list` [`tuple` [`float`, `float`]]
        Peak (x, y) positions in parent pixel coordinates.
    exposure : `lsst.afw.image.Exposure`
        Exposure the peaks were detected on.
    satBitmask : `int`, optional
        Mask bits of saturated pixels (default: the SAT plane).
    dilateIters : `int`, optional
        Number of dilations of the saturated region when assigning peaks to a blob.

    Returns
    -------
    peaks : `list` [`tuple` [`float`, `float`]]
        Peaks with every saturated blob reduced to one representative peak.
    nRemoved : `int`
        Number of peaks removed.
    """
    if satBitmask is None:
        satBitmask = exposure.mask.getPlaneBitMask("SAT")

    mask = exposure.mask.array
    data = exposure.image.array
    xy0 = exposure.getXY0()
    x0, y0 = xy0.getX(), xy0.getY()

    sat = (mask & satBitmask) != 0
    if not sat.any():
        return list(peakPositions), 0

    labels, nLabels = label(sat)
    dilated = binary_dilation(sat, iterations=dilateIters)
    ring = dilated & ~sat

    blobCentroids = {}
    for lab in range(1, nLabels + 1):
        blobMask = labels == lab
        blobRing = binary_dilation(blobMask, iterations=dilateIters) & ring
        if blobRing.sum() == 0:
            cy, cx = center_of_mass(blobMask)
        else:
            weights = np.clip(data[blobRing], 0, None)
            ys, xs = np.nonzero(blobRing)
            if weights.sum() <= 0:
                cy, cx = center_of_mass(blobMask)
            else:
                cy = np.average(ys, weights=weights)
                cx = np.average(xs, weights=weights)
        blobCentroids[lab] = (cx + x0, cy + y0)

    collapsed = []
    nRemoved = 0
    seenBlobs = set()
    for (x, y) in peakPositions:
        px = int(round(x - x0))
        py = int(round(y - y0))
        inBounds = 0 <= py < sat.shape[0] and 0 <= px < sat.shape[1]
        lab = labels[py, px] if inBounds and sat[py, px] else 0

        if lab == 0:
            collapsed.append((x, y))
            continue

        nRemoved += 1
        if lab not in seenBlobs:
            seenBlobs.add(lab)
            collapsed.append(blobCentroids[lab])

    print(f"  collapsed {nRemoved} saturated-core peaks -> "
          f"{len(seenBlobs)} representative peak(s) across {nLabels} saturated blob(s)")
    return collapsed, nRemoved


def mergeClosePeaks(peakPositions, exposure, minSeparation=None):
    """Merge peaks closer than ``minSeparation``, keeping the brightest pixel's peak.

    Prevents near-duplicate columns in the design matrix (mostly from re-detecting the same
    wings) from destabilising the joint solve.

    Parameters
    ----------
    peakPositions : `list` [`tuple` [`float`, `float`]]
        Peak (x, y) positions in parent pixel coordinates.
    exposure : `lsst.afw.image.Exposure`
        Exposure whose pixel values rank the peaks.
    minSeparation : `float`, optional
        Merge radius [px] (default: twice the PSF sigma).

    Returns
    -------
    mergedPositions : `list` [`tuple` [`float`, `float`]]
        Surviving peak positions.
    """
    if len(peakPositions) == 0:
        return []
    if minSeparation is None:
        shape = exposure.getPsf().computeShape(exposure.getBBox().getCenter())
        psfSigma = np.sqrt(0.5 * (shape.getIxx() + shape.getIyy()))
        minSeparation = 2 * psfSigma

    positions = np.array(peakPositions)
    data = exposure.image.array
    xy0 = exposure.getXY0()

    peakVals = np.array([
        data[int(y - xy0.getY()), int(x - xy0.getX())]
        if 0 <= int(y - xy0.getY()) < data.shape[0] and 0 <= int(x - xy0.getX()) < data.shape[1]
        else -np.inf
        for x, y in positions
    ])

    tree = cKDTree(positions)
    pairs = tree.query_pairs(r=minSeparation)

    parent = list(range(len(positions)))
    def find(i):
        while parent[i] != i:
            i = parent[i]
        return i
    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
    for i, j in pairs:
        union(i, j)

    clusters = {}
    for i in range(len(positions)):
        clusters.setdefault(find(i), []).append(i)

    mergedPositions = []
    nMerged = 0
    for members in clusters.values():
        if len(members) > 1:
            nMerged += len(members) - 1
        best = members[np.argmax(peakVals[members])]
        mergedPositions.append(tuple(positions[best]))

    print(f"  merged {nMerged} duplicate/spurious peaks -- "
          f"{len(peakPositions)} -> {len(mergedPositions)}")
    return mergedPositions


def jointPsfFit(exposure, psf, peakPositions, badBitmask, ampZeroTol=1e-6):
    """Fit all PSF amplitudes simultaneously with a non-negative least-squares solve.

    One sparse, weighted, bound-constrained solve against the original exposure.  Unlike an
    unconstrained solve followed by clipping, ``lsq_linear`` with bounds (0, inf) lets the other
    amplitudes re-adjust when one source is pushed to zero, so overlapping sources get a real
    non-negative amplitude instead of one of them being dropped from the model.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Image to fit (image, variance and mask planes are used).
    psf : `lsst.afw.detection.Psf`
        PSF model.
    peakPositions : `list` [`tuple` [`float`, `float`]]
        Source (x, y) positions in parent pixel coordinates.
    badBitmask : `int`
        Mask bits of pixels to ignore.
    ampZeroTol : `float`, optional
        Amplitudes at or below this count as zero.

    Returns
    -------
    amps : `numpy.ndarray`
        Fitted amplitudes (fluxes), one per position.
    modelFull : `numpy.ndarray`
        Model image over the exposure's box.
    nAtZeroBound : `int`
        Number of sources that converged to (numerically) zero amplitude -- sources that still
        fail to fit rather than being absent.
    """
    expBbox = exposure.getBBox()
    ny, nx = expBbox.getHeight(), expBbox.getWidth()
    data = exposure.image.array
    variance = exposure.variance.array
    mask = exposure.mask.array

    goodGlobal = (mask & badBitmask) == 0
    goodGlobal &= np.isfinite(data) & np.isfinite(variance) & (variance > 0)

    pixelIndex = -np.ones((ny, nx), dtype=np.int64)
    flatPositions = np.where(goodGlobal.ravel())[0]
    pixelIndex.ravel()[flatPositions] = np.arange(len(flatPositions))
    nPixelsUsed = len(flatPositions)

    dVec = data.ravel()[flatPositions]
    wVec = 1.0 / variance.ravel()[flatPositions]

    A = lil_matrix((nPixelsUsed, len(peakPositions)))
    for j, (x, y) in enumerate(peakPositions):
        point = geom.Point2D(x, y)
        psfImage = psf.computeImage(point)
        bbox = psfImage.getBBox()
        bbox.clip(expBbox)
        if bbox.isEmpty():
            continue
        modelStamp = np.clip(psfImage[bbox].array, 0, None)

        y0 = bbox.getMinY() - expBbox.getMinY()
        x0 = bbox.getMinX() - expBbox.getMinX()
        nzSy, nzSx = np.nonzero(modelStamp)
        for sy, sx in zip(nzSy, nzSx):
            gy, gx = y0 + sy, x0 + sx
            if 0 <= gy < ny and 0 <= gx < nx:
                idx = pixelIndex[gy, gx]
                if idx >= 0:
                    A[idx, j] = modelStamp[sy, sx]

    A = csr_matrix(A)
    sqrtW = np.sqrt(wVec)
    A_weighted = A.multiply(sqrtW[:, None])
    dWeighted = dVec * sqrtW

    result = lsq_linear(A_weighted, dWeighted, bounds=(0, np.inf),
                         method='trf', lsmr_tol='auto', max_iter=2000)
    amps = result.x

    nAtZeroBound = int(np.sum(amps <= ampZeroTol))

    modelFull = np.zeros((ny, nx))
    for j, (x, y) in enumerate(peakPositions):
        if amps[j] <= ampZeroTol:
            continue
        point = geom.Point2D(x, y)
        psfImage = psf.computeImage(point)
        bbox = psfImage.getBBox()
        bbox.clip(expBbox)
        if bbox.isEmpty():
            continue
        modelStamp = np.clip(psfImage[bbox].array, 0, None)
        yslice, xslice = bboxToSlices(bbox, expBbox)
        modelFull[yslice, xslice] += amps[j] * modelStamp

    return amps, modelFull, nAtZeroBound


# =====================================================================
# Main iterative detection + joint fit loop
# =====================================================================

def iterativeJointDetection(
    exposure,
    psf,
    badBitmask,
    nIterations=10,
    thresholdStart=8.0,
    thresholdStep=-1.5,
    thresholdMin=1.5,
    minPixels=1,
    minNewPeaks=2,          # stop if a pass finds fewer new peaks than this
    convergenceFrac=0.02,    # stop if residual std improves by less than this fraction
    verbose=True,
):
    """Detect, fit and subtract iteratively until the residual stops yielding new sources.

    Each pass: (1) detect on the current residual with a decreasing threshold, so fainter
    sources are searched for once brighter ones are gone; (2) collapse multiple peaks on
    saturated cores; (3) merge the new peaks with all peaks found so far; (4) refit all of them
    jointly against the original exposure with a non-negative solve; (5) rebuild the residual.
    Every pass is a full re-solve, so errors from an early pass are not compounded.  The number
    of zero-amplitude sources is printed each pass: if it is large or not shrinking, real
    sources are being fitted down to nothing.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Image to fit.
    psf : `lsst.afw.detection.Psf`
        PSF model.
    badBitmask : `int`
        Mask bits of pixels to ignore.
    nIterations : `int`, optional
        Maximum number of passes.
    thresholdStart : `float`, optional
        Detection threshold [sigma] of the first pass.
    thresholdStep : `float`, optional
        Change of the threshold per pass (negative: decreasing).
    thresholdMin : `float`, optional
        Lowest detection threshold [sigma].
    minPixels : `int`, optional
        Minimum footprint size of a detection.
    minNewPeaks : `int`, optional
        Stop when a pass finds fewer new peaks than this.
    convergenceFrac : `float`, optional
        Stop when the residual standard deviation improves by less than this fraction.
    verbose : `bool`, optional
        Print progress.

    Returns
    -------
    result : `dict`
        ``modelFull``, ``residualFull``, ``peaks``, ``amps``, the per-pass histories
        (``residStdHistory``, ``nPeaksHistory``, ``nZeroAmpHistory``) and ``firstSeen``
        (pass index at which each peak, keyed by its position rounded to 0.1 px, was found).
    """
    expBbox = exposure.getBBox()

    workingResidual = exposure.clone()
    workingResidual.setPsf(psf)

    allPeaks = []
    residStdHistory = [np.nanstd(workingResidual.image.array)]
    nPeaksHistory = [0]
    nZeroAmpHistory = [0]
    firstSeen = {}   # (rounded x, rounded y) -> iteration index

    threshold = thresholdStart
    modelFull = np.zeros((expBbox.getHeight(), expBbox.getWidth()))
    ampsFinal = np.array([])

    for it in range(nIterations):
        if verbose:
            print(f"--- Pass {it+1}/{nIterations}  threshold={threshold:.2f} sigma  "
                  f"residual std={residStdHistory[-1]:.3f} nJy  "
                  f"(peaks so far={len(allPeaks)}) ---")

        detConfig = SourceDetectionTask.ConfigClass()
        detConfig.thresholdValue = threshold
        detConfig.thresholdType = "pixel_stdev"
        detConfig.minPixels = minPixels
        detConfig.doTempLocalBackground = True
        detConfig.tempLocalBackground.binSize = 16
        detConfig.reEstimateBackground = False  # residual has no real sky background left to fit

        schemaIt = afwTable.SourceTable.makeMinimalSchema()
        detectionTask_it = SourceDetectionTask(schema=schemaIt, config=detConfig)
        tableIt = afwTable.SourceTable.make(schemaIt)

        try:
            detRes_it = detectionTask_it.run(tableIt, workingResidual)
        except TooManyMaskedPixelsError:
            if verbose:
                print("  background estimation failed on residual -- stopping")
            break

        newPeaks = [
            (peak.getFx(), peak.getFy())
            for fp in detRes_it.positive.getFootprints()
            for peak in fp.getPeaks()
        ]
        if verbose:
            print(f"  {len(newPeaks)} raw peaks detected on residual this pass")

        # collapse spurious multi-peak detections on saturated cores
        newPeaks, nCollapsed = collapseSaturatedPeaks(newPeaks, exposure)

        if len(newPeaks) < minNewPeaks and it > 0:
            if verbose:
                print("  too few new peaks -- treating as converged, stopping")
            break

        mergedPeaks = mergeClosePeaks(allPeaks + newPeaks, exposure)

        amps, modelFull, nAtZeroBound = jointPsfFit(exposure, psf, mergedPeaks, badBitmask)

        # record which pass first introduced each surviving peak position
        for (x, y) in mergedPeaks:
            key = (round(x, 1), round(y, 1))
            if key not in firstSeen:
                firstSeen[key] = it

        workingResidual = exposure.clone()
        workingResidual.image.array[:] = exposure.image.array - modelFull
        workingResidual.setPsf(psf)

        allPeaks = mergedPeaks
        ampsFinal = amps
        nPeaksHistory.append(len(allPeaks))
        nZeroAmpHistory.append(nAtZeroBound)

        newResidStd = np.nanstd(workingResidual.image.array)
        oldResidStd = residStdHistory[-1]
        residStdHistory.append(newResidStd)

        if verbose:
            print(f"  refit {len(allPeaks)} peaks total ({nAtZeroBound} at zero amplitude) -> "
                  f"residual std {oldResidStd:.3f} -> {newResidStd:.3f}")
            if nAtZeroBound > 0:
                print(f"  NOTE: {nAtZeroBound} peak(s) fit to zero amplitude this pass -- "
                      f"their flux is still sitting in the residual, unmodeled")

        if oldResidStd > 0 and it > 0:
            improvement = abs(oldResidStd - newResidStd) / oldResidStd
            if improvement < convergenceFrac:
                if verbose:
                    print(f"  residual std improved by only {improvement*100:.2f}% "
                          f"(< {convergenceFrac*100:.0f}%) -- converged, stopping")
                break

        threshold = max(threshold + thresholdStep, thresholdMin)

    return {
        'modelFull': modelFull,
        'residualFull': workingResidual.image.array.copy(),
        'peaks': allPeaks,
        'amps': ampsFinal,
        'residStdHistory': residStdHistory,
        'nPeaksHistory': nPeaksHistory,
        'nZeroAmpHistory': nZeroAmpHistory,
        'firstSeen': firstSeen,
    }


# =====================================================================
# Tiled data / model / residual comparison
# =====================================================================

def plotTiledComparison(
    dataFull,
    modelFull,
    residualFull,
    nSources=None,
    nPasses=None,
    tileSize=100,
    cmapImg=None,
    globalScale=True,
    skipEmpty=True,
    emptyStdThresh=1e-6,
    figsizePerPanel=5,
    showRowProfile=True,
    rowSelect='maxFlux',
    stretch='log'
):
    """Plot data, model and residual (and a row profile) for every tile of an image.

    The data and model panels share one colour scale when ``globalScale`` is True; the
    residual panel always scales to its own tile.  With ``showRowProfile`` a fourth panel
    shows data and model along one row of the tile (marked with a dashed cyan line), so the
    model can be compared with the data through a crowded section.

    Parameters
    ----------
    dataFull, modelFull, residualFull : `numpy.ndarray`
        Data, model and residual images (same shape).
    nSources : `int`, optional
        Number of sources, for the title.
    nPasses : `int`, optional
        Number of detection passes, for the title.
    tileSize : `int`, optional
        Tile size [px].
    cmapImg : `str`, optional
        Colour map of the data and model panels (default: 'mako' if available, else 'viridis').
    globalScale : `bool`, optional
        One colour scale for all data / model panels instead of one per tile.
    skipEmpty : `bool`, optional
        Skip tiles whose data standard deviation is below ``emptyStdThresh``.
    emptyStdThresh : `float`, optional
        See ``skipEmpty``.
    figsizePerPanel : `float`, optional
        Size of one panel [inches].
    showRowProfile : `bool`, optional
        Add the row-profile panel.
    rowSelect : `str`, `int` or callable, optional
        Row of each tile to profile: 'maxFlux' (the row with the most flux), a row index local
        to the tile, or a function ``rowSelect(tileData) -> row index``.
    stretch : `str`, optional
        Stretch of the data / model colour scale.
    """
    if cmapImg is None:
        cmapImg = 'mako' if 'mako' in plt.colormaps() else 'viridis'

    ny, nx = dataFull.shape

    if globalScale:
        vmaxGlobal = np.nanpercentile(dataFull, 99.5)
        normGlobal = simple_norm(dataFull, vmin=np.nanmin(dataFull),
                                   vmax=vmaxGlobal, stretch=stretch)

    yEdges = list(range(0, ny, tileSize))
    xEdges = list(range(0, nx, tileSize))

    titleBits = []
    if nSources is not None:
        titleBits.append(f"{nSources} sources")
    if nPasses is not None:
        titleBits.append(f"{nPasses} passes")
    suffix = f" ({', '.join(titleBits)})" if titleBits else ""

    nPanels = 4 if showRowProfile else 3

    for y0 in yEdges:
        for x0 in xEdges:
            y1 = min(y0 + tileSize, ny)
            x1 = min(x0 + tileSize, nx)

            tileData = dataFull[y0:y1, x0:x1]

            if skipEmpty and np.nanstd(tileData) < emptyStdThresh:
                continue

            tileModel = modelFull[y0:y1, x0:x1]
            tileResid = residualFull[y0:y1, x0:x1]

            if globalScale:
                norm = normGlobal
            else:
                vmax = np.nanpercentile(tileData, 99.5)
                norm = simple_norm(tileData, vmin=np.nanmin(tileData),
                                    vmax=vmax, stretch='log')

            # residual scale: ALWAYS local to this tile
            resScale = np.nanpercentile(np.abs(tileResid), 99)
            if not np.isfinite(resScale) or resScale == 0:
                resScale = 1.0

            fig, axes = plt.subplots(1, nPanels,
                                      figsize=(figsizePerPanel * nPanels, figsizePerPanel))
            fig.patch.set_facecolor('white')

            im0 = axes[0].imshow(tileData, origin='lower', norm=norm, cmap=cmapImg)
            axes[0].set_title('Data')
            fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

            im1 = axes[1].imshow(tileModel, origin='lower', norm=norm, cmap=cmapImg)
            axes[1].set_title('Joint model')
            fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

            im2 = axes[2].imshow(tileResid, origin='lower', cmap='RdBu_r',
                                  vmin=-resScale, vmax=resScale)
            axes[2].set_title('Residual')
            fig.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04)

            for ax in axes[:3]:
                ax.set_xticks([]); ax.set_yticks([])
                ax.grid(False)

            if showRowProfile:
                if rowSelect == 'maxFlux':
                    rowIdxLocal = int(np.nanargmax(np.nansum(tileData, axis=1)))
                elif callable(rowSelect):
                    rowIdxLocal = int(rowSelect(tileData))
                else:
                    rowIdxLocal = int(rowSelect)
                rowIdxLocal = int(np.clip(rowIdxLocal, 0, tileData.shape[0] - 1))
                rowIdxGlobal = y0 + rowIdxLocal

                # mark the sliced row on the data/model panels for context
                axes[0].axhline(rowIdxLocal, color='cyan', lw=0.8, ls='--')
                axes[1].axhline(rowIdxLocal, color='cyan', lw=0.8, ls='--')

                xPix = np.arange(x0, x1)
                ax3 = axes[3]
                ax3.plot(xPix, tileData[rowIdxLocal, :], color='k', lw=1.0, label='data')
                ax3.plot(xPix, tileModel[rowIdxLocal, :], color='#DD8452', lw=1.3, label='model')
                ax3.plot(xPix, tileResid[rowIdxLocal, :], color='#4C72B0', lw=0.8,
                         alpha=0.7, label='residual')
                ax3.axhline(0, color='gray', lw=0.5, ls=':')
                ax3.set_title(f'Row {rowIdxGlobal} cross-section')
                ax3.set_xlabel('x [px]')
                ax3.legend(fontsize=7, loc='upper right')

            fig.suptitle(f'Iterative joint fit{suffix} — tile x[{x0}:{x1}] y[{y0}:{y1}]',
                         fontsize=12, fontweight='medium', y=1.03)
            plt.tight_layout()
            plt.show()


# =====================================================================
# Per-iteration detections overlay
# =====================================================================

def plotDetectionsByIteration(exposureCrop, resultIter):
    """Plot every fitted peak over the image, coded by the pass that first found it.

    Peaks whose final amplitude is (numerically) zero -- they hit the non-negativity bound --
    are marked as failed fits with a red 'x'; they are the ones to inspect if sources seem to
    be missing from the model.

    Parameters
    ----------
    exposureCrop : `lsst.afw.image.Exposure`
        Image the pipeline ran on.
    resultIter : `dict`
        Output of `iterativeJointDetection`.

    Returns
    -------
    sourcesTbl : `astropy.table.Table`
        Positions, amplitudes and detection pass of the sources.
    """
    xy0 = exposureCrop.getXY0()  # (x0, y0) of this crop in the original exposure's frame
    x0, y0 = xy0.getX(), xy0.getY()

    # --- assemble a sources table from the pipeline result ---
    positions = np.array(resultIter['peaks'])
    amps = resultIter['amps']
    firstSeen = resultIter['firstSeen']

    iterations = np.array([
        firstSeen.get((round(x, 1), round(y, 1)), -1)
        for x, y in positions
    ])

    sourcesTbl = Table()
    sourcesTbl['slot_Centroid_x'] = positions[:, 0]
    sourcesTbl['slot_Centroid_y'] = positions[:, 1]
    sourcesTbl['ampFit'] = amps
    sourcesTbl['iteration'] = iterations

    nItersRun = len(resultIter['residStdHistory']) - 1
    iterColors = cm.plasma(np.linspace(0, 1, max(nItersRun, 1)))
    markers = ['o', 's', '^', 'D', 'v', 'P', '*', 'X']

    # --- plot ---
    fig, ax = plt.subplots(figsize=(8, 8))
    dataCrop = exposureCrop.image.array
    vmax = np.nanpercentile(dataCrop, 99.5)
    norm = simple_norm(dataCrop, vmin=np.nanmin(dataCrop), vmax=vmax, stretch='log')
    ax.imshow(dataCrop, origin='lower', norm=norm, cmap='gray')

    for it in range(nItersRun):
        maskIter = sourcesTbl['iteration'] == it
        injected = maskIter & (sourcesTbl['ampFit'] > 0)
        failed = maskIter & ~(sourcesTbl['ampFit'] > 0)

        # shift into array-relative coords by subtracting the crop's XY0
        xInj = sourcesTbl['slot_Centroid_x'][injected] - x0
        yInj = sourcesTbl['slot_Centroid_y'][injected] - y0

        ax.scatter(xInj, yInj,
                   s=25, facecolors='none', edgecolors=[iterColors[it]],
                   marker=markers[it % len(markers)], linewidths=1.3,
                   label=f'pass {it+1} (n={injected.sum()})')

        if failed.any():
            xFail = sourcesTbl['slot_Centroid_x'][failed] - x0
            yFail = sourcesTbl['slot_Centroid_y'][failed] - y0
            ax.scatter(xFail, yFail,
                       s=25, color='red', marker='x', linewidths=1.3,
                       label=f'pass {it+1} failed fit (n={failed.sum()})' if it == 0 else None)

    ax.set_xlim(0, dataCrop.shape[1])
    ax.set_ylim(0, dataCrop.shape[0])
    ax.legend(fontsize=8, loc='upper left', bbox_to_anchor=(1.02, 1.0), framealpha=0.9)
    ax.set_title('Detections by iteration pass')
    ax.set_xticks([]); ax.set_yticks([])
    plt.tight_layout()
    plt.show()

    return sourcesTbl


# =====================================================================
# Outer loop: background and PSF re-determination
# =====================================================================

def maskSourcesFromModel(exposure, modelFull, noiseSigma=None,
                             nSigma=4.0, dilateIters=1, verbose=True):
    """Flag pixels where the joint-fit model has significant flux as DETECTED.

    The next background fit then treats them as source, not sky.  The joint model contains all
    sources found so far, including ones too faint for the first background pass.  A loose
    threshold flags most of the image and starves the background fit; if the printed masked
    fraction climbs much above ~50 %, tighten ``nSigma`` or ``dilateIters``.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure whose mask is modified in place.
    modelFull : `numpy.ndarray`
        Model image over the exposure's box.
    noiseSigma : `float`, optional
        Pixel noise (default: sqrt of the median variance).
    nSigma : `float`, optional
        Model pixels above ``nSigma * noiseSigma`` are flagged.
    dilateIters : `int`, optional
        Number of dilations of the flagged region.
    verbose : `bool`, optional
        Print the masked fraction.

    Returns
    -------
    srcMask : `numpy.ndarray` [`bool`]
        The flagged pixels.
    """
    if noiseSigma is None:
        noiseSigma = np.sqrt(np.nanmedian(exposure.variance.array))

    srcMask = modelFull > (nSigma * noiseSigma)
    if dilateIters > 0:
        srcMask = binary_dilation(srcMask, iterations=dilateIters)

    detectedBit = exposure.mask.getPlaneBitMask("DETECTED")
    exposure.mask.array[srcMask] |= detectedBit

    nMasked = int(srcMask.sum())
    frac = nMasked / srcMask.size
    if verbose:
        print(f"  masked {nMasked} px ({100 * frac:.1f}%) "
              f"as source for background refit")
    if frac > 0.5:
        print(f"  WARNING: {100 * frac:.1f}% of the image is masked as source -- "
              f"background refit and PSF redetermination will be starved of "
              f"sky pixels/star candidates. Consider raising nSigma or "
              f"lowering dilateIters.")
    return srcMask


def refitBackground(exposure, algorithm='AKIMA_SPLINE', binSize=32,
                      useApprox=True,
                      statsMask=("BAD", "SAT", "EDGE", "NO_DATA", "DETECTED",
                                  "CR", "INTRP")):
    """Re-derive the background ignoring masked pixels and subtract it in place.

    Pass a fresh clone of the raw (never background-subtracted) exposure, already flagged by
    `maskSourcesFromModel`, or the background is subtracted twice.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure to background-subtract in place.
    algorithm : `str`, optional
        Interpolation algorithm (only used when ``useApprox`` is False).
    binSize : `int`, optional
        Background bin size [px].
    useApprox : `bool`, optional
        Fit a smooth Chebyshev surface to the binned statistics instead of interpolating;
        more robust when many pixels are masked.
    statsMask : `tuple` [`str`], optional
        Mask planes ignored when computing the background statistics.

    Returns
    -------
    exposure : `lsst.afw.image.Exposure`
        The background-subtracted exposure (same object).
    background : `lsst.afw.image.Image`
        The subtracted background image.
    """
    bgConfig = SubtractBackgroundConfig()
    bgConfig.algorithm = algorithm
    bgConfig.useApprox = useApprox
    bgConfig.binSize = binSize
    bgConfig.ignoredPixelMask = list(statsMask)

    bgTask = SubtractBackgroundTask(config=bgConfig)
    bgResult = bgTask.run(exposure)
    return exposure, bgResult.background.getImage()


# =====================================================================
# Step B0: turn raw (x, y) peak positions into a measured SourceCatalog,
# so a previously-converged peak list (or an externally-supplied set of
# "known good" star positions) can feed the star selector directly.
# =====================================================================

def catalogFromPeaks(exposure, peakPositions, schema, measureTask,
                        stampRadius=15, excludePositions=None,
                        excludeRadius=5.0, verbose=True):
    """Build a measured SourceCatalog from raw (x, y) peak positions.

    Each peak gets a circular footprint and one Peak record, and the catalogue is measured with
    the same ``measureTask`` as fresh detections, so the fields are comparable and usable by the
    star selector.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure to measure on.
    peakPositions : `list` [`tuple` [`float`, `float`]]
        Peak positions in parent pixel coordinates.
    schema : `lsst.afw.table.Schema`
        The schema used to build ``measureTask``.
    measureTask : `lsst.meas.base.SingleFrameMeasurementTask`
        Measurement task.
    stampRadius : `int`, optional
        Footprint radius [px].
    excludePositions : `list` [`tuple` [`float`, `float`]], optional
        Skip peaks within ``excludeRadius`` of these positions (e.g. existing detections).
    excludeRadius : `float`, optional
        See ``excludePositions`` [px].
    verbose : `bool`, optional
        Print how many peaks were skipped.

    Returns
    -------
    catalog : `lsst.afw.table.SourceCatalog`
        The measured catalogue.
    """
    table = afwTable.SourceTable.make(schema)
    catalog = afwTable.SourceCatalog(table)
    expBbox = exposure.getBBox()

    excludeTree = None
    if excludePositions is not None and len(excludePositions) > 0:
        excludeTree = cKDTree(np.asarray(excludePositions))

    nSkippedDup = 0
    nSkippedBounds = 0
    for (x, y) in peakPositions:
        if excludeTree is not None:
            dist, _ = excludeTree.query([x, y])
            if dist < excludeRadius:
                nSkippedDup += 1
                continue

        center = geom.Point2I(int(round(x)), int(round(y)))
        if not expBbox.contains(center):
            nSkippedBounds += 1
            continue

        # NOTE: SpanSet/Stencil location can differ between stack versions --
        # see module docstring if this errors on your stack.
        spans = afwGeom.SpanSet.fromShape(stampRadius, afwGeom.Stencil.CIRCLE, offset=center)
        spans = spans.clippedTo(expBbox)
        if spans.getArea() == 0:
            nSkippedBounds += 1
            continue

        footprint = afwDetect.Footprint(spans, expBbox)
        footprint.addPeak(x, y, 0.0)
        record = catalog.addNew()
        record.setFootprint(footprint)

    if len(catalog) > 0:
        measureTask.run(measCat=catalog, exposure=exposure)
    catalog = catalog.copy(deep=True)

    if verbose:
        print(f"  catalogFromPeaks: measured {len(catalog)} peak(s) "
              f"({nSkippedDup} skipped as duplicates of existing sources, "
              f"{nSkippedBounds} skipped out of bounds)")
    return catalog


import math

import math

# =====================================================================
# Step B-helper 1: version-safe PSF-determiner kernel-size setter
# =====================================================================

def _setPsfDeterminerKernelSize(psfDeterminerConfig, kernelSize, verbose=True):
    """Set the size fields that exist on this stack's PcaPsfDeterminerConfig.

    ``stampSize`` always; ``kernelSizeMin`` / ``kernelSizeMax`` only where present (they are
    absent on some stack versions).

    Parameters
    ----------
    psfDeterminerConfig : `lsst.meas.algorithms.PcaPsfDeterminerConfig`
        Config to modify in place.
    kernelSize : `int`
        Kernel size [px].
    verbose : `bool`, optional
        Print the extra fields that were set.

    Returns
    -------
    psfDeterminerConfig : `lsst.meas.algorithms.PcaPsfDeterminerConfig`
        The same config.
    """
    psfDeterminerConfig.stampSize = kernelSize

    extraFields = ["kernelSizeMin", "kernelSizeMax"]
    setFields = []
    for field in extraFields:
        if hasattr(psfDeterminerConfig, field):
            setattr(psfDeterminerConfig, field, kernelSize)
            setFields.append(field)

    if verbose and setFields:
        print(f"  PSF determiner: also set {setFields} = {kernelSize}")
    return psfDeterminerConfig


# =====================================================================
# Step B-helper 2: drop PSF candidates whose stamp extraction silently failed
# =====================================================================

def _filterValidPsfCandidates(psfCandidateList, verbose=True):
    """Drop PSF candidates whose postage stamp is unusable (non-finite or empty).

    chi2 / amplitude / status are not checked: they keep their pre-fit defaults until the PSF
    determiner runs.

    Parameters
    ----------
    psfCandidateList : `list`
        PSF candidates from MakePsfCandidatesTask.
    verbose : `bool`, optional
        Print the dropped candidates.

    Returns
    -------
    good : `list`
        Candidates with a usable stamp.
    """
    good = []
    nBad = 0
    for cand in psfCandidateList:
        try:
            stamp = cand.getMaskedImage().getImage().array
        except Exception as e:
            nBad += 1
            if verbose:
                print(f"  dropping candidate id={cand.getSource().getId()}: "
                      f"stamp re-fetch failed ({e})")
            continue
        if not np.isfinite(stamp).all() or stamp.sum() <= 0:
            nBad += 1
            continue
        good.append(cand)

    if verbose and nBad > 0:
        print(f"  filtered out {nBad} PSF candidate(s) with unusable stamps "
              f"-> {len(good)} of {len(psfCandidateList)} remain")
    return good


# =====================================================================
# Step B: redetect + reselect stars + redetermine PSF on the freshly
# background-subtracted exposure.
# =====================================================================

def redetectAndRefitPsf(exposure, detThreshold=3.0, detMinPixels=2,
                            detLocalBgBin=16, detBgBin=32,
                            starSelectorConfig=None,
                            psfDeterminerConfig=None,
                            kernelSize=41,
                            extraPeakPositions=None,
                            forceIncludeStarPositions=None,
                            peakStampRadius=15,
                            verbose=True):
    """Detect, deblend, measure, select stars and determine a new PSF.

    Run in each outer pass, so it works on the newly cleaned background.  Extra peak positions (e.g. the previous pass's converged catalogue) are added
    before star selection; forced star positions are appended after selection, bypassing its
    cuts.  Candidates with failed stamps are filtered out before the PSF determiner runs (even
    one can make it reject the whole set).  If too few candidates survive, the exposure's
    current PSF is kept and a "NOT updated" note is printed.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Background-subtracted exposure.
    detThreshold : `float`, optional
        Detection threshold [sigma].
    detMinPixels : `int`, optional
        Minimum footprint size.
    detLocalBgBin : `int`, optional
        Bin size of the temporary local background [px].
    detBgBin : `int`, optional
        Bin size of the detection background [px].
    starSelectorConfig : `lsst.meas.algorithms.ObjectSizeStarSelectorConfig`, optional
        Star selector configuration.
    psfDeterminerConfig : `lsst.meas.algorithms.PcaPsfDeterminerConfig`, optional
        PSF determiner configuration.
    kernelSize : `int`, optional
        PSF kernel size [px], used for both the candidates and the determiner.
    extraPeakPositions : `list` [`tuple` [`float`, `float`]], optional
        Additional peaks measured and offered to the star selector.
    forceIncludeStarPositions : `list` [`tuple` [`float`, `float`]], optional
        Star positions appended to the selected stars regardless of the selector's cuts.
    peakStampRadius : `int`, optional
        Footprint radius of the extra / forced peaks [px].
    verbose : `bool`, optional
        Print progress.

    Returns
    -------
    psfNew : `lsst.afw.detection.Psf`
        The new PSF (or the current one if it could not be determined).
    sources : `lsst.afw.table.SourceCatalog`
        Measured detections.
    starCat : `lsst.afw.table.SourceCatalog`
        Selected PSF stars.
    cellSet : `lsst.afw.math.SpatialCellSet`
        PSF candidates used by the determiner.
    """
    installPsfTask = InstallGaussianPsfTask()
    installPsfTask.run(exposure)

    schema = afwTable.SourceTable.makeMinimalSchema()

    detConfig = SourceDetectionConfig()
    detConfig.thresholdValue = detThreshold
    detConfig.thresholdType = "stdev"
    detConfig.minPixels = detMinPixels
    detConfig.doTempLocalBackground = True
    detConfig.tempLocalBackground.binSize = detLocalBgBin
    detConfig.background.binSize = detBgBin

    detectionTask = SourceDetectionTask(schema=schema, config=detConfig)
    deblendTask = SourceDeblendTask(schema=schema)
    measureTask = SingleFrameMeasurementTask(schema=schema)

    table = afwTable.SourceTable.make(schema)
    detRes = detectionTask.run(table, exposure)
    sources = detRes.sources

    deblendTask.run(exposure, sources)
    measureTask.run(measCat=sources, exposure=exposure)
    sources = sources.copy(deep=True)

    # ---- fold in the full previously-converged peak catalog, if given ----
    if extraPeakPositions:
        freshXy = (np.array([(r.getX(), r.getY()) for r in sources])
                    if len(sources) else None)
        extraCat = catalogFromPeaks(
            exposure, extraPeakPositions, schema, measureTask,
            stampRadius=peakStampRadius, excludePositions=freshXy,
            excludeRadius=peakStampRadius, verbose=verbose)
        if len(extraCat) > 0:
            combined = afwTable.SourceCatalog(sources.table)
            combined.extend(sources, deep=True)
            combined.extend(extraCat, deep=True)
            # extend() does not guarantee contiguous storage, and
            # ObjectSizeStarSelectorTask (like several afw table consumers)
            # requires a contiguous catalog -- without this copy it raises
            # "Input catalogs for source selection must be contiguous."
            sources = combined.copy(deep=True)
            if verbose:
                print(f"  combined catalog for star selection: {len(sources)} sources "
                      f"({len(extraCat)} pulled in from the prior converged peak catalog)")

    starSelector = (ObjectSizeStarSelectorTask(config=starSelectorConfig)
                     if starSelectorConfig else ObjectSizeStarSelectorTask())
    selection = starSelector.run(sourceCat=sources, exposure=exposure)
    starCat = selection.sourceCat

    # ---- force the ORIGINAL star candidates back in, bypassing the selector ----
    if forceIncludeStarPositions:
        selectedXy = (np.array([(r.getX(), r.getY()) for r in starCat])
                       if len(starCat) else None)
        forcedCat = catalogFromPeaks(
            exposure, forceIncludeStarPositions, schema, measureTask,
            stampRadius=peakStampRadius, excludePositions=selectedXy,
            excludeRadius=peakStampRadius, verbose=verbose)
        if len(forcedCat) > 0:
            combinedStars = afwTable.SourceCatalog(starCat.table)
            combinedStars.extend(starCat, deep=True)
            combinedStars.extend(forcedCat, deep=True)
            nBefore = len(starCat)
            # same contiguity requirement as above -- MakePsfCandidatesTask
            # also expects/iterates a contiguous catalog.
            starCat = combinedStars.copy(deep=True)
            if verbose:
                print(f"  forced back {len(forcedCat)} original PSF-star "
                      f"candidate(s) -> {nBefore} -> {len(starCat)} total candidates")

    if verbose:
        print(f"  {len(starCat)} PSF-star candidates out of {len(sources)} detections")

    makeCandidatesConfig = MakePsfCandidatesTask.ConfigClass()
    makeCandidatesConfig.kernelSize = kernelSize
    makeCandidatesTask = MakePsfCandidatesTask(config=makeCandidatesConfig)
    candidatesResult = makeCandidatesTask.run(starCat, exposure=exposure)
    psfCandidateList = candidatesResult.psfCandidates

    # drop candidates whose stamp extraction silently failed -- this is
    # the fix for the "No viable PSF candidates survive" failures
    psfCandidateList = _filterValidPsfCandidates(psfCandidateList, verbose=verbose)

    if psfDeterminerConfig is None:
        psfDeterminerConfig = PcaPsfDeterminerTask.ConfigClass()
        psfDeterminerConfig.spatialOrder = 1
        psfDeterminerConfig.sizeCellX = 512
        psfDeterminerConfig.sizeCellY = 512
        psfDeterminerConfig.maxCandidates = 1500
        psfDeterminerConfig.doRejectBlends = False # changed from True
        psfDeterminerConfig.nStarPerCell = 10

    _setPsfDeterminerKernelSize(psfDeterminerConfig, kernelSize, verbose=verbose)

    psfPrev = exposure.getPsf()
    point = exposure.getBBox().getCenter()

    def _fwhm(psfObj):
        shape = psfObj.computeShape(point)
        return 2.355 * np.sqrt(0.5 * (shape.getIxx() + shape.getIyy()))

    fwhmPrev = _fwhm(psfPrev)

    psfDeterminer = PcaPsfDeterminerTask(config=psfDeterminerConfig)
    try:
        if len(psfCandidateList) == 0:
            raise IndexError("No viable PSF candidates survive (all candidates "
                              "filtered out before determinePsf was even called)")
        psfNew, cellSet_new = psfDeterminer.determinePsf(
            exposure=exposure, psfCandidateList=psfCandidateList)
        exposure.setPsf(psfNew)
        if verbose:
            fwhmNew = _fwhm(psfNew)
            print(f"  PSF updated: FWHM {fwhmPrev:.3f} -> {fwhmNew:.3f} px (crop center)")
    except IndexError as e:
        if verbose:
            print(f"  WARNING: PSF redetermination failed ({e}) with "
                  f"{len(psfCandidateList)} valid candidates from {len(starCat)} "
                  f"star-selected sources -- PSF NOT updated this pass "
                  f"(kept previous, FWHM {fwhmPrev:.3f} px)")
        psfNew = psfPrev
        cellSet_new = None

    return psfNew, sources, starCat, cellSet_new


# =====================================================================
# Step C: the outer loop itself
# =====================================================================

def outerLoopPsfBackgroundRefinement(
    rawExposureCrop,
    psfInitial,
    badBitmask,
    initialStarCat=None,
    nOuterIterations=3,
    nInnerIterations=10,
    innerKwargs=None,
    outerConvergenceFrac=0.02,
    bgUseApprox=True,
    psfKernelSize=41,
    verbose=True,
):
    """Alternate inner-loop fits with background and PSF re-derivation.

    Pass 0 runs `iterativeJointDetection` with the current background and PSF.  Pass k > 0
    masks every known source (from pass k-1's model) on a fresh clone of the raw image,
    re-derives the background, re-determines the PSF (offering pass k-1's catalogue and forcing
    in ``initialStarCat``) and reruns the inner loop.

    Parameters
    ----------
    rawExposureCrop : `lsst.afw.image.Exposure`
        Crop before any background subtraction.
    psfInitial : `lsst.afw.detection.Psf`
        PSF of the first pass.
    badBitmask : `int`
        Mask bits of pixels to ignore.
    initialStarCat : `lsst.afw.table.SourceCatalog`, optional
        Star selector output whose positions are forced into every PSF determination.
    nOuterIterations : `int`, optional
        Maximum number of outer passes.
    nInnerIterations : `int`, optional
        Maximum number of inner passes (when ``innerKwargs`` is not given).
    innerKwargs : `dict`, optional
        Keyword arguments of `iterativeJointDetection`.
    outerConvergenceFrac : `float`, optional
        Stop when the residual standard deviation improves by less than this fraction.
    bgUseApprox : `bool`, optional
        See `refitBackground`.
    psfKernelSize : `int`, optional
        PSF kernel size [px].
    verbose : `bool`, optional
        Print progress.

    Returns
    -------
    outerHistory : `list` [`dict`]
        One dict per outer pass: ``psf``, ``exposure``, ``result`` (inner-loop output),
        ``residStd``, ``nPeaks``.
    """
    if innerKwargs is None:
        innerKwargs = dict(nIterations=nInnerIterations, thresholdStart=8.0,
                             thresholdStep=-1.5, thresholdMin=1.0,
                             minPixels=1, minNewPeaks=1, convergenceFrac=0.02)

    originalStarPositions = None
    if initialStarCat is not None and len(initialStarCat) > 0:
        originalStarPositions = [
            (r.getCentroid().getX(), r.getCentroid().getY()) for r in initialStarCat
        ]
        if verbose:
            print(f"  {len(originalStarPositions)} original PSF-star position(s) "
                  f"will be forced into every subsequent outer pass")

    raw = rawExposureCrop.clone()  # kept pristine for every re-background pass
    exposureWorking = rawExposureCrop.clone()
    exposureWorking.setPsf(psfInitial)
    psfWorking = psfInitial
    bgImageWorking = None  # no separately-tracked bg image for outer pass 0

    outerHistory = []

    for outerIt in range(nOuterIterations):
        if verbose:
            print(f"\n===== Outer pass {outerIt + 1}/{nOuterIterations} =====")

        result = iterativeJointDetection(
            exposureWorking, psfWorking, badBitmask, **innerKwargs)
        finalResidStd = result['residStdHistory'][-1]

        outerHistory.append({
            'psf': psfWorking,
            'exposure': exposureWorking,
            'result': result,
            'residStd': finalResidStd,
            'nPeaks': len(result['peaks']),
            'bgImage': bgImageWorking,
        })

        if verbose:
            print(f"  outer pass {outerIt + 1} converged: "
                  f"{len(result['peaks'])} peaks, residual std={finalResidStd:.4f}")

        if outerIt > 0:
            prevStd = outerHistory[-2]['residStd']
            improvement = (abs(prevStd - finalResidStd) / prevStd
                            if prevStd > 0 else 0)
            if improvement < outerConvergenceFrac:
                if verbose:
                    print(f"  outer residual std improved by only "
                          f"{improvement * 100:.2f}% -- outer loop converged, stopping")
                break

        if outerIt == nOuterIterations - 1:
            break  # no need to re-derive background/PSF after the final pass

        # ---- re-derive background + PSF from this pass's converged model ----
        exposureNext = raw.clone()
        maskSourcesFromModel(exposureNext, result['modelFull'], verbose=verbose)
        exposureNext, bgImage_next = refitBackground(exposureNext, useApprox=bgUseApprox)

        psfNew, sourcesNew, starCat_new, cellSet_new = redetectAndRefitPsf(
            exposureNext, kernelSize=psfKernelSize,
            extraPeakPositions=result['peaks'],
            forceIncludeStarPositions=originalStarPositions,
            verbose=verbose)

        exposureWorking = exposureNext
        psfWorking = psfNew
        bgImageWorking = bgImage_next

    return outerHistory


# =====================================================================
# Diagnostics
# =====================================================================

def _extractPeakXy(peaks):
    """Extract (x, y) arrays from peaks of unknown type.

    Tries afw-style ``getX()``/``getY()``, ``.x``/``.y`` attributes, dicts with 'x'/'y' keys and
    finally indexable (x, y) pairs.

    Parameters
    ----------
    peaks : iterable
        The peaks.

    Returns
    -------
    xs, ys : `numpy.ndarray`
        Peak coordinates.
    """
    xs, ys = [], []
    for p in peaks:
        if hasattr(p, "getX") and hasattr(p, "getY"):
            xs.append(p.getX())
            ys.append(p.getY())
        elif hasattr(p, "x") and hasattr(p, "y"):
            xs.append(p.x)
            ys.append(p.y)
        elif isinstance(p, dict):
            xs.append(p["x"])
            ys.append(p["y"])
        else:
            xs.append(p[0])
            ys.append(p[1])
    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)


def plotPsfModelsByOuterPass(outerHistory, exposureCenter=None):
    """Plot the PSF kernel image used in each outer pass.

    Parameters
    ----------
    outerHistory : `list` [`dict`]
        Outer-loop history (each entry has ``psf`` and ``exposure``).
    exposureCenter : `lsst.geom.Point2D`, optional
        Position at which to render the kernel (default: the exposure centre).
    """
    n = len(outerHistory)
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 4))
    if n == 1:
        axes = [axes]

    for i, (ax, h) in enumerate(zip(axes, outerHistory)):
        point = exposureCenter or h["exposure"].getBBox().getCenter()
        psfImg = h["psf"].computeKernelImage(geom.Point2D(point)).array
        im = ax.imshow(psfImg, origin="lower", cmap="inferno")
        ax.set_title(f"Outer pass {i + 1} PSF")
        ax.set_xticks([])
        ax.set_yticks([])
        plt.colorbar(im, ax=ax, fraction=0.046)

    plt.tight_layout()
    plt.show()


def plotBackgroundByOuterPass(outerHistory):
    """Plot the background re-derived going into each outer pass.

    The first pass uses the background subtracted before the pipeline ran and is skipped.

    Parameters
    ----------
    outerHistory : `list` [`dict`]
        Outer-loop history (entries may have ``bgImage``).
    """
    passesWithBg = [(i, h) for i, h in enumerate(outerHistory)
                       if h.get("bgImage") is not None]
    if not passesWithBg:
        print("No re-derived background images to show yet -- these are only "
              "captured starting from outer pass 2 onward.")
        return

    fig, axes = plt.subplots(1, len(passesWithBg), figsize=(4 * len(passesWithBg), 4))
    if len(passesWithBg) == 1:
        axes = [axes]

    for ax, (i, h) in zip(axes, passesWithBg):
        im = ax.imshow(h["bgImage"].array, origin="lower", cmap="viridis")
        ax.set_title(f"Background entering pass {i + 1}")
        ax.set_xticks([])
        ax.set_yticks([])
        plt.colorbar(im, ax=ax, fraction=0.046)

    plt.tight_layout()
    plt.show()


def plotSourcesOverlaidByPass(outerHistory, markerSize=25, vminPct=5, vmaxPct=99):
    """Plot the data of each outer pass with its converged peak positions.

    Parameters
    ----------
    outerHistory : `list` [`dict`]
        Outer-loop history.
    markerSize : `float`, optional
        Marker size.
    vminPct, vmaxPct : `float`, optional
        Percentiles of the display range.
    """
    n = len(outerHistory)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 5))
    if n == 1:
        axes = [axes]

    for i, (ax, h) in enumerate(zip(axes, outerHistory)):
        exposure = h["exposure"]
        data = exposure.image.array
        vmin, vmax = np.nanpercentile(data, [vminPct, vmaxPct])
        ax.imshow(data, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)

        # FIX: `_extractPeakXy` returns positions in the exposure's
        # absolute (parent-frame) pixel coordinates (e.g. afw's
        # getFx()/getFy()). `data = exposure.image.array` above, plotted
        # with imshow, is ARRAY-relative -- pixel (0, 0) of the array is
        # NOT pixel (0, 0) of the parent frame unless the exposure's XY0
        # happens to be (0, 0). The previous version scattered the raw
        # (unshifted) positions on top of the array-relative image, which
        # is exactly why the markers showed up as one solid offset block
        # instead of tracking the stars: every peak was shifted by the
        # same (x0, y0) and clipped hard by the axis limits.
        xy0 = exposure.getXY0()
        x0, y0 = xy0.getX(), xy0.getY()

        xs, ys = _extractPeakXy(h["result"]["peaks"])
        ax.scatter(xs - x0, ys - y0, s=markerSize, facecolors="none",
                   edgecolors="#00FF00", linewidths=0.8)
        ax.set_xlim(0, data.shape[1])
        ax.set_ylim(0, data.shape[0])
        ax.set_title(f"Outer pass {i + 1}: {len(xs)} sources")
        ax.set_xticks([])
        ax.set_yticks([])

    plt.tight_layout()
    plt.show()

def plotOuterLoopSummary(outerHistory):
    """Plot the residual standard deviation and the source count against outer pass.

    Parameters
    ----------
    outerHistory : `list` [`dict`]
        Outer-loop history (entries have ``residStd`` and ``nPeaks``).
    """
    nOuter = len(outerHistory)
    passes = np.arange(1, nOuter + 1)
    residStds = [h['residStd'] for h in outerHistory]
    nPeaks = [h['nPeaks'] for h in outerHistory]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    fig.patch.set_facecolor('white')

    axes[0].plot(passes, residStds, 'o-', color='#DD8452')
    axes[0].set_xlabel('outer pass')
    axes[0].set_ylabel('final residual std [nJy]')
    axes[0].set_title('Residual std vs outer pass')
    axes[0].set_xticks(passes)

    axes[1].plot(passes, nPeaks, 'o-', color='#4C72B0')
    axes[1].set_xlabel('outer pass')
    axes[1].set_ylabel('total converged peaks')
    axes[1].set_title('Peak count vs outer pass')
    axes[1].set_xticks(passes)

    plt.tight_layout()
    plt.show()


def plotPsfFwhmByOuterPass(outerHistory, exposureCenter=None):
    """Plot the FWHM of the PSF used in each outer pass.

    Parameters
    ----------
    outerHistory : `list` [`dict`]
        Outer-loop history.
    exposureCenter : `lsst.geom.Point2D`, optional
        Position at which to evaluate the PSF (default: the exposure centre).

    Returns
    -------
    fwhms : `list` [`float`]
        FWHM [px] per outer pass.
    """
    fwhms = []
    for h in outerHistory:
        point = exposureCenter or h['exposure'].getBBox().getCenter()
        shape = h['psf'].computeShape(geom.Point2D(point))
        fwhm = 2.355 * np.sqrt(0.5 * (shape.getIxx() + shape.getIyy()))
        fwhms.append(fwhm)

    passes = np.arange(1, len(outerHistory) + 1)
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(passes, fwhms, 'o-', color='#55A868')
    ax.set_xlabel('outer pass')
    ax.set_ylabel('PSF FWHM [px] (crop center)')
    ax.set_title('PSF FWHM vs outer pass')
    ax.set_xticks(passes)
    plt.tight_layout()
    plt.show()
    return fwhms


# =====================================================================
# Cutout analysis: shape cut, S/N, diagnostics, the cutout pipeline
# =====================================================================

import time
from contextlib import contextmanager
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.patheffects as pe
from matplotlib.colors import ListedColormap
from astropy.visualization import simple_norm
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.optimize import lsq_linear
from scipy.spatial import cKDTree
from scipy.ndimage import map_coordinates
from scipy.optimize import curve_fit
from scipy.signal import fftconvolve

import lsst.geom as geom
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
from lsst.meas.algorithms import SourceDetectionTask
from lsst.meas.algorithms.subtractBackground import TooManyMaskedPixelsError

CMAP = 'mako' if 'mako' in plt.colormaps() else 'viridis'
GOOD_STATUS = getattr(afwMath.SpatialCellCandidate, 'GOOD', 1)

# no red (masked pixels) and no yellow/gold (support contours)
SRC_PALETTE = ['#00e5ff', '#ff00ff', '#76ff03', '#448aff',
               '#ff80ab', '#b388ff', '#ffffff', '#1de9b6']
_OUTLINE = [pe.withStroke(linewidth=3.5, foreground='k')]
_TEXT_OUTLINE = [pe.withStroke(linewidth=2.5, foreground='k')]


# =====================================================================
# Timing
# =====================================================================

class StepTimer:
    """Collect timing records of pipeline steps.

    Each record is a dict with ``label``, ``seconds``, ``outer``, ``inner`` and ``kind``
    ('compute' or 'viz').
    """

    def __init__(self):
        """Start an empty timer."""
        self.records = []
        self.tStart = time.perf_counter()

    @contextmanager
    def step(self, label, outer=None, inner=None, kind='compute', verbose=False):
        """Time a block of code (context manager).

        Parameters
        ----------
        label : `str`
            Step name.
        outer : `int`, optional
            Outer-pass index.
        inner : `int`, optional
            Inner-iteration index.
        kind : `str`, optional
            'compute' or 'viz'.
        verbose : `bool`, optional
            Print the elapsed time.

        Yields
        ------
        rec : `dict`
            The record; ``seconds`` is filled in when the block exits.
        """
        rec = dict(label=label, seconds=0.0, outer=outer, inner=inner, kind=kind)
        t0 = time.perf_counter()
        try:
            yield rec
        finally:
            rec['seconds'] = time.perf_counter() - t0
            self.records.append(rec)
            if verbose:
                print(f"    [time] {label}: {rec['seconds']:.3f}s")

    def add(self, label, seconds, outer=None, inner=None, kind='compute'):
        """Add a record with a known duration.

        Parameters
        ----------
        label : `str`
            Step name.
        seconds : `float`
            Duration [s].
        outer : `int`, optional
            Outer-pass index.
        inner : `int`, optional
            Inner-iteration index.
        kind : `str`, optional
            'compute' or 'viz'.
        """
        self.records.append(dict(label=label, seconds=seconds, outer=outer,
                                 inner=inner, kind=kind))

    def reset(self):
        """Remove all records and restart the clock."""
        self.records.clear()
        self.tStart = time.perf_counter()


TIMER = StepTimer()


def plotTimingSummary(timer=TIMER):
    """Plot the total time per step and the time per inner iteration.

    Parameters
    ----------
    timer : `StepTimer`, optional
        Timer to summarise (default: the module-level ``TIMER``).
    """
    comp = [r for r in timer.records if r['kind'] == 'compute']
    viz = [r for r in timer.records if r['kind'] == 'viz']
    if not comp:
        print("No timing records.")
        return

    totals = defaultdict(float)
    counts = defaultdict(int)
    for r in comp:
        totals[r['label']] += r['seconds']
        counts[r['label']] += 1

    print("\n================ TIMING SUMMARY ================")
    print(f"{'step':45s} {'calls':>6s} {'total [s]':>10s} {'mean [s]':>10s}")
    for lab, tot in sorted(totals.items(), key=lambda kv: -kv[1]):
        print(f"{lab:45s} {counts[lab]:6d} {tot:10.3f} {tot / counts[lab]:10.3f}")
    print(f"{'-' * 75}")
    print(f"compute total : {sum(totals.values()):.2f} s")
    print(f"plotting total: {sum(r['seconds'] for r in viz):.2f} s (not included above)")
    print(f"wall clock    : {time.perf_counter() - timer.tStart:.2f} s")

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(17, 6),
                                   gridspec_kw={'width_ratios': [1, 1.5]})
    labs = sorted(totals, key=lambda k: totals[k])
    ax0.barh(labs, [totals[k] for k in labs], color='#4C72B0')
    ax0.set_xlabel('total seconds')
    ax0.set_title('Total compute time by step')
    ax0.grid(axis='x', alpha=0.3)

    inner = [r for r in comp if r['inner'] is not None and r['outer'] is not None]
    keys = sorted({(r['outer'], r['inner']) for r in inner})
    innerLabels = sorted({r['label'] for r in inner})
    cmap = plt.cm.tab20(np.linspace(0, 1, max(len(innerLabels), 1)))
    bottoms = np.zeros(len(keys))
    for lab, col in zip(innerLabels, cmap):
        vals = np.array([sum(r['seconds'] for r in inner
                             if r['label'] == lab and (r['outer'], r['inner']) == k)
                         for k in keys])
        ax1.bar(range(len(keys)), vals, bottom=bottoms, color=col, label=lab)
        bottoms += vals
    ax1.set_xticks(range(len(keys)))
    ax1.set_xticklabels([f"o{o + 1}.i{i + 1}" for o, i in keys], rotation=60, fontsize=8)
    ax1.set_ylabel('seconds')
    ax1.set_title('Time per inner iteration (outer.inner)')
    ax1.legend(fontsize=7, ncol=2)
    ax1.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.show()
    plt.close(fig)


# =====================================================================
# Small geometry / plotting helpers
# =====================================================================

def makeCutout(exposureCrop, xRange=(0, 100), yRange=(0, 100)):
    """Cut a sub-exposure out of a crop, in coordinates relative to the crop's corner.

    Parameters
    ----------
    exposureCrop : `lsst.afw.image.Exposure`
        The crop.
    xRange, yRange : `tuple` [`int`, `int`], optional
        Pixel ranges relative to the crop's origin (end exclusive).

    Returns
    -------
    cutout : `lsst.afw.image.Exposure`
        A deep copy of the sub-exposure.

    Raises
    ------
    TypeError
        Raised if ``exposureCrop`` is not an Exposure.
    """
    if not hasattr(exposureCrop, 'getXY0'):
        raise TypeError(
            "makeCutout() got a "
            f"{type(exposureCrop).__name__}, not an lsst.afw.image Exposure. "
            "The `exposureCrop` variable was probably reassigned after "
            "`exposureCrop = exposure[...].clone()`. Re-run the cell that builds it "
            "(with .clone() and .setPsf(psf)) before calling runCutoutPipeline.")
    xy0 = exposureCrop.getXY0()
    bbox = geom.Box2I(
        geom.Point2I(xy0.getX() + xRange[0], xy0.getY() + yRange[0]),
        geom.Extent2I(xRange[1] - xRange[0], yRange[1] - yRange[0]))
    return exposureCrop[bbox].clone()


def _extent(bbox):
    """Return the matplotlib extent (pixel edges) of a box."""
    return (bbox.getMinX() - 0.5, bbox.getMaxX() + 0.5,
            bbox.getMinY() - 0.5, bbox.getMaxY() + 0.5)


def _coords(bbox):
    """Return the pixel x and y coordinates covered by a box."""
    xs = bbox.getMinX() + np.arange(bbox.getWidth())
    ys = bbox.getMinY() + np.arange(bbox.getHeight())
    return xs, ys


def _dataNorm(arr):
    """Return a log-stretch display normalisation for an image."""
    return simple_norm(arr, vmin=np.nanmin(arr), vmax=np.nanpercentile(arr, 99.5),
                       stretch='log')


def _overlayBad(ax, bad, extent):
    """Overlay bad pixels in translucent red."""
    if bad is None or not bad.any():
        return
    ax.imshow(np.ma.masked_where(~bad, np.ones(bad.shape)), extent=extent,
              origin='lower', cmap=ListedColormap(['red']), alpha=0.35,
              interpolation='nearest')


def _cleanAxes(ax):
    """Turn the axis grid off."""
    ax.grid(False)


def _peaksArray(peaks):
    """Return peaks as an (N, 2) float array."""
    return np.asarray(peaks, dtype=float).reshape(-1, 2)


def _fwhm(psfObj, point):
    """Return the FWHM [px] of a PSF at a point, from its second moments."""
    shape = psfObj.computeShape(point)
    return 2.355 * np.sqrt(0.5 * (shape.getIxx() + shape.getIyy()))


def _footprintsToMask(footprints, exposure):
    """Rasterise detection footprints into a boolean mask over an exposure."""
    bbox = exposure.getBBox()
    m = np.zeros((bbox.getHeight(), bbox.getWidth()), bool)
    try:
        for fp in footprints:
            fb = fp.getBBox()
            arr = np.asarray(fp.spans.asArray()).astype(bool)
            y0 = fb.getMinY() - bbox.getMinY()
            x0 = fb.getMinX() - bbox.getMinX()
            sub = m[y0:y0 + arr.shape[0], x0:x0 + arr.shape[1]]
            sub |= arr[:sub.shape[0], :sub.shape[1]]
    except Exception as e:  # SpanSet API differs between stack versions
        print(f"  (could not rasterize detection footprints: {e})")
    return m


def _matchPrevAmps(peaks, prevPeaks, prevAmps, tol=1.5):
    """Return each peak's amplitude from the previous iteration (NaN if it is new)."""
    out = np.full(len(peaks), np.nan)
    if len(peaks) == 0 or len(prevPeaks) == 0:
        return out
    d, i = cKDTree(prevPeaks).query(peaks)
    ok = d < tol
    out[ok] = np.asarray(prevAmps)[i[ok]]
    return out


# ---- PSF stamp with tamed wings / peak filtering -----------------------

# (r_start, r_end) as fractions of the stamp half-width. Inside r_start the PSF is
# untouched; between them it is cosine-tapered to exactly 0. None = raw kernel.
WING_TAPER = (0.8, 1.0)


def psfStamp(psf, x, y, clipBbox, taper=None):
    """Render the PSF at (x, y), clipped to a box, with its far wings tamed.

    With ``taper`` the median of the outermost pixel ring (edge pedestal) is subtracted and the
    stamp is tapered radially to zero, so each source ends smoothly in a disc instead of a
    hard-edged square.

    Parameters
    ----------
    psf : `lsst.afw.detection.Psf`
        PSF model.
    x, y : `float`
        Position in parent pixel coordinates.
    clipBbox : `lsst.geom.Box2I`
        Box to clip the stamp to.
    taper : `tuple` [`float`, `float`], optional
        Inner and outer taper radii as fractions of the stamp half-size.

    Returns
    -------
    stamp : `numpy.ndarray` or `None`
        The clipped stamp, or None if it misses ``clipBbox``.
    bbox : `lsst.geom.Box2I`
        Box of the clipped stamp.
    """
    img = psf.computeImage(geom.Point2D(x, y))
    full = img.getBBox()
    arr = np.array(img.array, dtype=float)
    if taper is not None:
        ring = np.concatenate([arr[0, :], arr[-1, :], arr[:, 0], arr[:, -1]])
        arr = arr - np.median(ring)
        h = 0.5 * (min(arr.shape) - 1)
        yy, xx = np.mgrid[full.getMinY():full.getMaxY() + 1,
                          full.getMinX():full.getMaxX() + 1]
        r = np.hypot(xx - x, yy - y)
        r0, r1 = taper[0] * h, taper[1] * h
        w = np.clip((r1 - r) / max(r1 - r0, 1e-9), 0.0, 1.0)
        arr = arr * (0.5 - 0.5 * np.cos(np.pi * w))
    bb = geom.Box2I(full)
    bb.clip(clipBbox)
    if bb.isEmpty():
        return None, bb
    sl = (slice(bb.getMinY() - full.getMinY(), bb.getMaxY() - full.getMinY() + 1),
          slice(bb.getMinX() - full.getMinX(), bb.getMaxX() - full.getMinX() + 1))
    return arr[sl], bb


def dropMaskedPeaks(peaks, exposure, planes=('EDGE', 'NO_DATA')):
    """Drop peaks on EDGE / NO_DATA pixels.

    SAT is deliberately not in the default list, so saturated stars are kept.

    Parameters
    ----------
    peaks : `list` [`tuple` [`float`, `float`]]
        Peak positions in parent pixel coordinates.
    exposure : `lsst.afw.image.Exposure`
        Exposure providing the mask.
    planes : `tuple` [`str`], optional
        Mask planes that reject a peak.

    Returns
    -------
    keep : `list` [`tuple` [`float`, `float`]]
        Remaining peaks.
    nDropped : `int`
        Number of peaks dropped.
    """
    mk = exposure.mask
    planesPresent = mk.getMaskPlaneDict()
    bits = 0
    for p in planes:
        if p in planesPresent:
            bits |= mk.getPlaneBitMask(p)
    if bits == 0:
        return list(peaks), 0
    bbox, arr = exposure.getBBox(), mk.array
    keep = []
    for x, y in peaks:
        ix, iy = int(round(x)) - bbox.getMinX(), int(round(y)) - bbox.getMinY()
        if 0 <= ix < arr.shape[1] and 0 <= iy < arr.shape[0] and (arr[iy, ix] & bits) == 0:
            keep.append((x, y))
    return keep, len(peaks) - len(keep)


def dropOutOfBounds(peaks, exposure, margin=0.0):
    """Drop peaks outside the exposure's box grown by ``margin`` px.

    Parameters
    ----------
    peaks : `list` [`tuple` [`float`, `float`]]
        Peak positions in parent pixel coordinates.
    exposure : `lsst.afw.image.Exposure`
        Exposure providing the box.
    margin : `float`, optional
        Growth of the box [px]; a negative margin drops peaks within ``|margin|`` px of the edge.

    Returns
    -------
    keep : `list` [`tuple` [`float`, `float`]]
        Remaining peaks.
    nDropped : `int`
        Number of peaks dropped.
    """
    bb = exposure.getBBox()
    x0, x1 = bb.getMinX() - margin, bb.getMaxX() + margin
    y0, y1 = bb.getMinY() - margin, bb.getMaxY() + margin
    keep = [(x, y) for x, y in peaks if x0 <= x <= x1 and y0 <= y <= y1]
    return keep, len(peaks) - len(keep)


# ---- over-counting guards ---------------------------------------------

def effectiveVariance(exposure, model=None, frac=0.03, poissonK=0.0):
    """Return the variance used to weight the fit.

    ``variance + (frac * model)**2 + poissonK * model``: the pixel variance plus a fractional
    model-error term and missing photon noise.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure providing the variance plane.
    model : `numpy.ndarray`, optional
        Current model image (no extra terms if None).
    frac : `float`, optional
        Fractional model error.
    poissonK : `float`, optional
        Extra variance per unit model flux.

    Returns
    -------
    var : `numpy.ndarray`
        The effective variance.
    """
    var = exposure.variance.array.copy()
    if model is not None:
        m = np.clip(model, 0, None)
        var += (frac * m) ** 2 + poissonK * m
    return var


def rejectNearBright(newPeaks, oldPeaks, oldIsBright, minSep):
    """Drop new peaks within ``minSep`` px of an existing bright source.

    Parameters
    ----------
    newPeaks : `list` [`tuple` [`float`, `float`]]
        Candidate peaks.
    oldPeaks : array-like
        Existing source positions, shape (N, 2).
    oldIsBright : array-like [`bool`]
        Which existing sources are bright.
    minSep : `float`
        Rejection radius [px].

    Returns
    -------
    keep : `list` [`tuple` [`float`, `float`]]
        Remaining candidates.
    n_rejected : `int`
        Number of candidates rejected.
    """
    old = np.asarray(oldPeaks, float).reshape(-1, 2)[np.asarray(oldIsBright, bool)]
    if len(old) == 0 or len(newPeaks) == 0:
        return list(newPeaks), 0
    d, _ = cKDTree(old).query(np.asarray(newPeaks, float))
    keep = [p for p, di in zip(newPeaks, d) if di > minSep]
    return keep, len(newPeaks) - len(keep)


def pruneAndMerge(fit, peaks, minSnr=5.0, minSep=3.0,
                    brightFlags=None, nearBrightSep=None, nearBrightSnr=8.0):
    """Drop low-S/N sources, then merge pairs closer than ``minSep`` px.

    The brighter source of a pair survives at the amplitude-weighted position.  Sources within
    ``nearBrightSep`` px of a bright source must reach ``nearBrightSnr`` instead of
    ``minSnr``.  Uses the cheap diagonal-only S/N; the honest S/N is applied by
    `eliminateWeak`.

    Parameters
    ----------
    fit : `dict`
        Joint fit (from `jointPsfFitDiag`) of ``peaks``.
    peaks : array-like
        Source positions, shape (N, 2).
    minSnr : `float`, optional
        Minimum S/N.
    minSep : `float`, optional
        Merge radius [px].
    brightFlags : array-like [`bool`], optional
        Which sources are bright.
    nearBrightSep : `float`, optional
        Radius around bright sources with the stricter S/N cut [px].
    nearBrightSnr : `float`, optional
        S/N required near bright sources.

    Returns
    -------
    peaks : `list` [`tuple` [`float`, `float`]]
        Surviving source positions.
    """
    var, amps = fit['variance'], fit['amps']
    peaks = np.asarray(peaks, float).copy()
    snr = np.zeros(len(peaks))
    for j, c in enumerate(fit['comps']):
        if c is None:
            continue
        g = fit['good'][c['ys'], c['xs']]
        v = var[c['ys'], c['xs']][g]
        snr[j] = amps[j] * np.sqrt(np.sum(c['stamp'][g] ** 2 / v))

    thresh = np.full(len(peaks), float(minSnr))
    nStrict = 0
    if nearBrightSep and brightFlags is not None:
        bf = np.asarray(brightFlags, bool)
        if bf.any() and len(bf) == len(peaks):
            d, _ = cKDTree(peaks[bf]).query(peaks)
            near = (d > 0.5) & (d < nearBrightSep)
            thresh[near] = np.maximum(thresh[near], nearBrightSnr)
            nStrict = int(near.sum())
    keep = snr >= thresh
    alive = keep.copy()
    for a in np.argsort(-amps):
        if not alive[a]:
            continue
        d = np.hypot(*(peaks - peaks[a]).T)
        for b in np.flatnonzero((d < minSep) & alive):
            if b == a:
                continue
            w = amps[a] + amps[b]
            peaks[a] = (amps[a] * peaks[a] + amps[b] * peaks[b]) / w
            alive[b] = False
    print(f"  prune/merge: {len(peaks)} -> {alive.sum()} "
          f"({(~keep).sum()} low S/N, {(keep & ~alive).sum()} merged; "
          f"{nStrict} held to S/N>={nearBrightSnr:g} near bright stars)")
    return [tuple(p) for p in peaks[alive]]


# ---- honest S/N + backward elimination ---------------------------------

def _pixelNoiseScale(fit, exposure, modelFrac=0.01):
    """Return the robust standard deviation of residual / sigma on sky pixels (at least 1)."""
    resid = exposure.image.array - fit['modelFull']
    z = resid / np.sqrt(fit['variance'])
    mmax = fit['modelFull'].max()
    sky = fit['good'] & (fit['modelFull'] < modelFrac * mmax)
    zz = z[sky]
    if zz.size < 50:
        zz = z[fit['good']]
    return max(1.0, 1.4826 * np.median(np.abs(zz - np.median(zz))))


def _matchedFilterNoiseScale(fit, exposure, coreHalf=10, floorSigma=0.25):
    """Return the robust standard deviation of the residual's matched-filter S/N on sky pixels.

    Unlike the pixel scale this includes pixel-to-pixel noise correlation that the variance
    plane ignores.  Returns 1.0 if there are too few clean pixels or no PSF.
    """
    psf = exposure.getPsf()
    if psf is None:
        return 1.0
    k = psf.computeKernelImage(exposure.getBBox().getCenter()).array
    h = k.shape[0] // 2
    c = min(coreHalf, h)
    k = k[h - c:h + c + 1, h - c:h + c + 1]      # core only: tails / ringing don't matter
    good, var, model = fit['good'], fit['variance'], fit['modelFull']
    resid = exposure.image.array - model
    inv = np.where(good, 1.0 / var, 0.0)
    num = fftconvolve(np.where(good, resid, 0.0) * inv, k[::-1, ::-1], mode='same')
    den = np.sqrt(np.clip(fftconvolve(inv, (k ** 2)[::-1, ::-1], mode='same'), 1e-30, None))
    snr = num / den
    sig = np.sqrt(np.nanmedian(var[good]))
    sky = good & (model < floorSigma * sig) & (den > 0.7 * np.median(den[good]))
    if sky.sum() < 200:
        return 1.0
    zz = snr[sky]
    return 1.4826 * np.median(np.abs(zz - np.median(zz)))


# How far S/N values move from the pixel-level scale toward the matched-filter scale.
#   0.0 = pixel scale only (original behaviour, most sensitive, most fakes)
#   1.0 = full matched-filter scale (strictest; can drop real faint sources)
#   0.5 = halfway
MATCHED_NOISE_WEIGHT = 0.5


def empiricalNoiseScale(fit, exposure, modelFrac=0.01, matched=True,
                          weight=None):
    """Return the empirical noise scale of a fit's residual.

    Parameters
    ----------
    fit : `dict`
        Joint fit (from `jointPsfFitDiag`).
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    modelFrac : `float`, optional
        Pixels with model below this fraction of the model maximum count as sky.
    matched : `bool`, optional
        True: scale for S/N values -- the pixel scale moved ``weight`` of the way towards the
        matched-filter scale when that is larger.  False: the pixel scale only (for chi2).
    weight : `float`, optional
        See ``matched`` (default ``MATCHED_NOISE_WEIGHT``).

    Returns
    -------
    scale : `float`
        Noise scale (at least 1).
    """
    s = _pixelNoiseScale(fit, exposure, modelFrac)
    if matched:
        w = MATCHED_NOISE_WEIGHT if weight is None else weight
        sm = _matchedFilterNoiseScale(fit, exposure)
        if sm > s:
            s = s + w * (sm - s)
    return s


def ampSnrFull(fit, noiseScale=1.0):
    """Return honest amplitude S/N values, including the covariance of overlapping sources.

    ``S/N = amp / sqrt(diag(N^-1))`` with ``N = A^T W A`` over all sources, divided by the
    empirical noise scale.

    Parameters
    ----------
    fit : `dict`
        Joint fit (from `jointPsfFitDiag`).
    noiseScale : `float`, optional
        Empirical noise scale dividing the S/N.

    Returns
    -------
    snr : `numpy.ndarray`
        S/N per source.
    """
    comps, var, good = fit['comps'], fit['variance'], fit['good']
    n = len(comps)
    pix = -np.ones(good.shape, dtype=np.int64)
    idx = np.flatnonzero(good.ravel())
    pix.ravel()[idx] = np.arange(len(idx))
    R, C, V = [], [], []
    for j, c in enumerate(comps):
        if c is None:
            continue
        p = pix[c['ys'], c['xs']]
        ok = p >= 0
        w = 1.0 / np.sqrt(var[c['ys'], c['xs']][ok])
        R.append(p[ok]); C.append(np.full(int(ok.sum()), j)); V.append(c['stamp'][ok] * w)
    if not R:
        return np.zeros(n)
    A = coo_matrix((np.concatenate(V), (np.concatenate(R), np.concatenate(C))),
                   shape=(len(idx), n)).tocsr()
    cov = np.linalg.pinv((A.T @ A).toarray())
    sig = np.sqrt(np.clip(np.diag(cov), 1e-30, None)) * noiseScale
    return fit['amps'] / sig


def eliminateWeak(exposure, psf, peaks, fit, badBitmask, minSnr=5.0,
                   dropFrac=0.25, maxRounds=30, **fitkw):
    """Remove sources below ``minSnr`` by backward elimination.

    Each round drops the weakest ~``dropFrac`` of the sub-threshold sources (at least one) and
    re-solves, until every source has an honest S/N of at least ``minSnr``.  The last source is
    never removed.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    psf : `lsst.afw.detection.Psf`
        PSF model.
    peaks : `list` [`tuple` [`float`, `float`]]
        Source positions; ``fit`` must correspond to them.
    fit : `dict`
        Joint fit of ``peaks``.
    badBitmask : `int`
        Mask bits of pixels to ignore.
    minSnr : `float`, optional
        Minimum honest S/N.
    dropFrac : `float`, optional
        Fraction of the failing sources dropped per round.
    maxRounds : `int`, optional
        Maximum number of rounds.
    **fitkw
        Keyword arguments of `jointPsfFitDiag`.

    Returns
    -------
    peaks : `list` [`tuple` [`float`, `float`]]
        Surviving source positions.
    fit : `dict`
        Joint fit of the survivors.
    """
    peaks = [tuple(p) for p in peaks]
    n0 = len(peaks)
    scale = empiricalNoiseScale(fit, exposure)
    for _ in range(maxRounds):
        if len(peaks) <= 1:
            break
        scale = empiricalNoiseScale(fit, exposure)
        snr = ampSnrFull(fit, scale)
        weak = np.flatnonzero(snr < minSnr)
        if len(weak) == 0:
            break
        order = weak[np.argsort(snr[weak])]
        nDrop = max(1, int(np.ceil(dropFrac * len(weak))))
        drop = set(order[:min(nDrop, len(peaks) - 1)].tolist())
        peaks = [p for i, p in enumerate(peaks) if i not in drop]
        fit = jointPsfFitDiag(exposure, psf, peaks, badBitmask,
                                 varModel=fit['modelFull'], **fitkw)
    print(f"  eliminateWeak: noiseScale={scale:.2f}, {n0} -> {len(peaks)} sources "
          f"(honest S/N>={minSnr:g})")
    return peaks, fit


# ---- PSF-shape constraint ----------------------------------------------

def _fitWithPlane(T, D, w, dx, dy):
    """Fit ``D ~ a*T + c + gx*dx + gy*dy`` by weighted least squares; return (a, a/sigma_a, chi2)."""
    A = np.column_stack([T, np.ones_like(T), dx, dy]) * w[:, None]
    b = D * w
    try:
        coef, *_ = np.linalg.lstsq(A, b, rcond=None)
        cov = np.linalg.inv(A.T @ A)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, np.nan
    r = b - A @ coef
    return float(coef[0]), float(coef[0] / np.sqrt(cov[0, 0])), float(r @ r)


def shapeMetrics(exposure, psf, peaks, fit, noiseScale=1.0, chi2Scale=None,
                  bgSigma=4.0, widthScales=(0.6, 0.75, 1.0, 1.3, 1.7), minPix=15):
    """Measure how PSF-like each source is.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    psf : `lsst.afw.detection.Psf`
        PSF model.
    peaks : `list` [`tuple` [`float`, `float`]]
        Source positions; ``fit`` must correspond to them (same order).
    fit : `dict`
        Joint fit of ``peaks``.
    noiseScale : `float`, optional
        Divides the S/N values (matched-filter calibrated).
    chi2Scale : `float`, optional
        Divides chi2 (pixel-level calibrated; default ``noiseScale``).
    bgSigma : `float`, optional
        Fit radius in units of the PSF sigma.
    widthScales : `tuple` [`float`], optional
        Trial PSF width scale factors.
    minPix : `int`, optional
        Minimum number of good pixels to judge a source.

    Returns
    -------
    metrics : `dict` [`str`, `numpy.ndarray`]
        ``snr`` (local, noise-scaled), ``chi2red``, ``bestScale``, ``dchi2Width`` (> 0: the
        data prefer ``bestScale`` over the PSF width) and ``contrast``; NaN = not judged.
    """
    cs = noiseScale if chi2Scale is None else chi2Scale
    bbox = exposure.getBBox()
    X0, Y0 = bbox.getMinX(), bbox.getMinY()
    data, model = exposure.image.array, fit['modelFull']
    var, good = fit['variance'], fit['good']
    ny, nx = data.shape
    shp = psf.computeShape(bbox.getCenter())
    sig = np.sqrt(0.5 * (shp.getIxx() + shp.getIyy()))
    R = bgSigma * sig
    h = int(np.ceil(R))
    n = len(peaks)
    out = dict(snr=np.full(n, np.nan), chi2red=np.full(n, np.nan),
               bestScale=np.full(n, np.nan), dchi2Width=np.full(n, np.nan),
               contrast=np.full(n, np.nan))
    buf = np.zeros((ny, nx))
    scales = np.asarray(widthScales, float)
    i1 = int(np.argmin(np.abs(scales - 1.0)))

    for j, (x, y) in enumerate(peaks):
        c = fit['comps'][j]
        if c is None:
            continue
        buf[c['ys'], c['xs']] = c['stamp']                       # unit-amplitude template
        ix, iy = int(round(x)) - X0, int(round(y)) - Y0
        xa, xb = max(ix - h, 0), min(ix + h + 1, nx)
        ya, yb = max(iy - h, 0), min(iy + h + 1, ny)
        sl = (slice(ya, yb), slice(xa, xb))
        P = buf[sl].copy()
        buf[c['ys'], c['xs']] = 0.0
        if xb <= xa or yb <= ya:
            continue
        D = data[sl] - (model[sl] - c['amp'] * P)                # neighbours removed
        # contrast = this source's peak / the NEIGHBOURS' model flux at its position.
        # PSF-wing error is a fraction of that neighbour flux, so a candidate that does
        # not stand well above it is indistinguishable from wing mismatch.
        jy, jx = iy - ya, ix - xa
        if 0 <= jy < P.shape[0] and 0 <= jx < P.shape[1]:
            contam = model[sl][jy, jx] - c['amp'] * P[jy, jx]
            out['contrast'][j] = c['amp'] * P.max() / max(contam, 1e-9)
        yy, xx = np.mgrid[Y0 + ya:Y0 + yb, X0 + xa:X0 + xb]
        dx, dy = xx - x, yy - y
        sel = good[sl] & (np.hypot(dx, dy) < R)
        if sel.sum() < minPix:
            continue
        w = 1.0 / np.sqrt(var[sl][sel])
        Dv, dxv, dyv = D[sel], dx[sel], dy[sel]

        _, sPsf, chi2 = _fitWithPlane(P[sel], Dv, w, dxv, dyv)
        out['snr'][j] = sPsf / noiseScale
        out['chi2red'][j] = chi2 / max(sel.sum() - 4, 1) / cs ** 2     # pixel-level scale

        r2 = dxv ** 2 + dyv ** 2
        snrW = np.array([_fitWithPlane(np.exp(-0.5 * r2 / (s * sig) ** 2),
                                          Dv, w, dxv, dyv)[1]
                          for s in scales]) / noiseScale
        if np.all(np.isfinite(snrW)):
            k = int(np.argmax(snrW))
            out['bestScale'][j] = scales[k]
            out['dchi2Width'][j] = snrW[k] * abs(snrW[k]) - snrW[i1] * abs(snrW[i1])
    return out


def judgeShape(m, minSnrLocal=5.0, chi2Max=2.0, widthRange=(0.75, 1.35),
                dchi2Width=9.0, trustSnr=30.0, minContrast=20.0,
                weakSnr=10.0, dchi2WidthWeak=4.0):
    """Decide which sources fail the PSF-shape tests.

    Sources with local S/N of at least ``trustSnr`` (including saturated ones) skip the chi2,
    width and contrast tests.  Untrusted sources with S/N below ``weakSnr`` use the looser
    width threshold ``dchi2WidthWeak``.

    Parameters
    ----------
    m : `dict` [`str`, `numpy.ndarray`]
        Output of `shapeMetrics`.
    minSnrLocal : `float`, optional
        Minimum local S/N.
    chi2Max : `float`, optional
        Maximum reduced chi2.
    widthRange : `tuple` [`float`, `float`], optional
        Allowed range of the best-fitting width scale.
    dchi2Width : `float`, optional
        Improvement in chi2 needed to reject on width.
    trustSnr : `float`, optional
        S/N above which a source is trusted.
    minContrast : `float`, optional
        Minimum ratio of the source's peak to the neighbours' model flux at its position (None
        or 0 = off); 20 ~ (5 % wing error) x (4x margin).
    weakSnr : `float`, optional
        S/N below which ``dchi2WidthWeak`` applies.
    dchi2WidthWeak : `float`, optional
        Width threshold for weak sources.

    Returns
    -------
    fail : `numpy.ndarray` [`bool`]
        Sources failing any test.
    score : `numpy.ndarray`
        Badness score (higher = worse); orders the elimination.
    flags : `dict` [`str`, `numpy.ndarray`]
        Per-test failures: ``snr``, ``chi2``, ``width``, ``contrast``.
    """
    snr, chi2 = m['snr'], m['chi2red']
    bs, dw = m['bestScale'], m['dchi2Width']
    con = m['contrast']
    valid = np.isfinite(snr)
    trusted = valid & (snr >= trustSnr)
    weak = valid & ~trusted & (snr < weakSnr)
    dwThr = np.where(weak, dchi2WidthWeak, dchi2Width)

    fSnr = valid & (snr < minSnrLocal)
    fChi = valid & ~trusted & (chi2 > chi2Max)
    fWid = (valid & ~trusted & np.isfinite(bs) &
             ((bs < widthRange[0]) | (bs > widthRange[1])) & (dw > dwThr))
    fCon = (valid & ~trusted & np.isfinite(con) & (con < minContrast)
             if minContrast else np.zeros_like(valid))
    with np.errstate(divide='ignore', invalid='ignore'):
        terms = [minSnrLocal / np.clip(np.nan_to_num(snr, nan=1e3), 1e-3, None),
                 np.nan_to_num(chi2 / chi2Max),
                 np.where(fWid, np.nan_to_num(dw) / dwThr, 0.0)]
        if minContrast:
            terms.append(minContrast / np.clip(np.nan_to_num(con, nan=1e6), 1e-3, None))
        score = np.maximum.reduce(terms)
    flags = dict(snr=fSnr, chi2=fChi, width=fWid, contrast=fCon)
    return fSnr | fChi | fWid | fCon, score, flags


def eliminateNonpsf(exposure, psf, peaks, fit, badBitmask, shapeKwargs=None,
                     dropFrac=0.25, maxRounds=15, verbose=True, **fitkw):
    """Remove sources that fail the PSF-shape test by backward elimination.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    psf : `lsst.afw.detection.Psf`
        PSF model.
    peaks : `list` [`tuple` [`float`, `float`]]
        Source positions; ``fit`` must correspond to them.
    fit : `dict`
        Joint fit of ``peaks``.
    badBitmask : `int`
        Mask bits of pixels to ignore.
    shapeKwargs : `dict`, optional
        Keyword arguments of `shapeMetrics` and `judgeShape`.
    dropFrac : `float`, optional
        Fraction of the failing sources dropped per round.
    maxRounds : `int`, optional
        Maximum number of rounds.
    verbose : `bool`, optional
        Print the rejected sources.
    **fitkw
        Keyword arguments of `jointPsfFitDiag`.

    Returns
    -------
    peaks : `list` [`tuple` [`float`, `float`]]
        Surviving source positions.
    fit : `dict`
        Joint fit of the survivors.
    rejected : `list` [`tuple` [`float`, `float`, `float`]]
        (x, y, score) of the removed sources.
    """
    kw = dict(shapeKwargs or {})
    mkw = {k: kw.pop(k) for k in ('bgSigma', 'widthScales', 'minPix') if k in kw}
    peaks = [tuple(p) for p in peaks]
    n0, rejected = len(peaks), []
    why = dict(snr=0, chi2=0, width=0, contrast=0)       # tally of tests that flagged dropped sources
    scale = pscale = 1.0
    rows = []                      # per-rejected-source metrics, printed at the end
    for _ in range(maxRounds):
        if len(peaks) <= 1:
            break
        scale = empiricalNoiseScale(fit, exposure)                  # for S/N values
        pscale = empiricalNoiseScale(fit, exposure, matched=False)  # for chi2
        m = shapeMetrics(exposure, psf, peaks, fit, noiseScale=scale,
                          chi2Scale=pscale, **mkw)
        fail, score, flags = judgeShape(m, **kw)
        idx = np.flatnonzero(fail)
        if len(idx) == 0:
            break
        order = idx[np.argsort(-score[idx])]
        nDrop = min(max(1, int(np.ceil(dropFrac * len(idx)))), len(peaks) - 1)
        drop = set(order[:nDrop].tolist())
        for i in drop:
            for key in why:
                why[key] += int(flags[key][i])
            rows.append((peaks[i][0], peaks[i][1], m['snr'][i], m['chi2red'][i],
                         m['bestScale'][i], m['dchi2Width'][i], m['contrast'][i],
                         '+'.join(k for k in why if flags[k][i]) or '-'))
        rejected += [(peaks[i][0], peaks[i][1], float(score[i])) for i in drop]
        peaks = [p for i, p in enumerate(peaks) if i not in drop]
        fit = jointPsfFitDiag(exposure, psf, peaks, badBitmask,
                                 varModel=fit['modelFull'], **fitkw)
    if verbose and rows:
        rows.sort(key=lambda r: -(r[2] if np.isfinite(r[2]) else -1e9))
        print("  shape-cut rejects, highest local S/N first (these are the ones most "
              "likely to be REAL):")
        print(f"    {'x':>8s} {'y':>8s} {'S/N':>6s} {'chi2':>6s} {'wscale':>6s} "
              f"{'dchi2w':>7s} {'contrast':>8s}  failed")
        for x, y, s, c, bs, dw, con, whyS in rows[:25]:
            print(f"    {x:8.1f} {y:8.1f} {s:6.1f} {c:6.2f} {bs:6.2f} {dw:7.1f} "
                  f"{con:8.1f}  {whyS}")
        if len(rows) > 25:
            print(f"    ... and {len(rows) - 25} more (all lower S/N)")
    if verbose:
        print(f"  shape cut: {n0} -> {len(peaks)} sources ({len(rejected)} removed; "
              f"flagged by: S/N {why['snr']}, chi2 {why['chi2']}, width {why['width']}, "
              f"contrast {why['contrast']} -- a source can trip several); "
              f"noise scales: S/N {scale:.2f}, pixel {pscale:.2f}")
    return peaks, fit, rejected


def plotShapeRejects(exposure, rejectedByIter, peaks=None,
                       title='sources removed by shape cut'):
    """Plot the sources removed by the shape cut over the data.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    rejectedByIter : `list` [`list` [`tuple`]]
        Per iteration, the (x, y, score) of rejected sources (magenta x).
    peaks : array-like, optional
        Kept source positions (cyan).
    title : `str`, optional
        Plot title.
    """
    bbox = exposure.getBBox()
    data = exposure.image.array
    rej = np.array([(x, y) for r in rejectedByIter for (x, y, _) in r]).reshape(-1, 2)
    fig, ax = plt.subplots(figsize=(8, 7.5))
    ax.imshow(data, extent=_extent(bbox), origin='lower', cmap='gray',
              norm=simple_norm(data, stretch='asinh', min_percent=1, max_percent=99.5))
    if peaks is not None and len(peaks):
        pk = _peaksArray(peaks)
        ax.scatter(pk[:, 0], pk[:, 1], facecolors='none', edgecolors='cyan', s=60,
                   label=f'kept ({len(pk)})')
    if len(rej):
        ax.scatter(rej[:, 0], rej[:, 1], marker='x', c='#ff00ff', s=50, lw=1.5,
                   label=f'rejected ({len(rej)})')
    ax.legend(fontsize=8)
    ax.set_title(title)
    plt.tight_layout()
    plt.show()


# ---- colour-coding / circle helpers ---------------------------------

def selectSourcePages(snap, sourceView='all', sourcesPerPage=6, maxPages=None):
    """Choose the sources that get per-source pages.

    Parameters
    ----------
    snap : `dict`
        Inner-loop snapshot (``fit``, ``prevAmp``, ``isNew``).
    sourceView : `str`, optional
        'all', 'new' (added this iteration) or 'changed' (new, or amplitude changed by more
        than 10 %).
    sourcesPerPage : `int`, optional
        Sources per page.
    maxPages : `int`, optional
        Maximum number of pages (None = no limit).

    Returns
    -------
    pages : `list` [`list` [`int`]]
        Source indices per page, brightest first.
    """
    amps, prevAmp, isNew = snap['fit']['amps'], snap['prevAmp'], snap['isNew']
    idx = np.arange(len(amps))
    if sourceView == 'new':
        idx = idx[isNew]
    elif sourceView == 'changed':
        with np.errstate(divide='ignore', invalid='ignore'):
            rel = np.abs(amps - prevAmp) / np.abs(prevAmp)
        idx = idx[isNew | (rel > 0.10) | ~np.isfinite(rel)]
    idx = idx[np.argsort(-amps[idx])]
    pages = [[int(j) for j in idx[i:i + sourcesPerPage]]
             for i in range(0, len(idx), sourcesPerPage)]
    return pages[:maxPages] if maxPages is not None else pages


def sourceColors(pages):
    """Assign each source on the pages a colour from ``SRC_PALETTE``.

    Parameters
    ----------
    pages : `list` [`list` [`int`]]
        Output of `selectSourcePages`.

    Returns
    -------
    colors : `dict` [`int`, `str`]
        Colour per source index.
    """
    flat = [j for page in pages for j in page]
    return {j: SRC_PALETTE[k % len(SRC_PALETTE)] for k, j in enumerate(flat)}


def _viewWindow(bbox, px, py, viewHalf):
    """Return the (unclipped) array index of a source and the clipped window around it."""
    X0, Y0 = bbox.getMinX(), bbox.getMinY()
    ix, iy = int(round(px)) - X0, int(round(py)) - Y0
    xa, xb = max(ix - viewHalf, 0), min(ix + viewHalf + 1, bbox.getWidth())
    ya, yb = max(iy - viewHalf, 0), min(iy + viewHalf + 1, bbox.getHeight())
    return ix, iy, xa, xb, ya, yb


def _circleRadius(exposure, circleRadius=None):
    """Return the marker radius: the given one, else 0.75 FWHM, else 5 px."""
    if circleRadius is not None:
        return circleRadius
    try:
        return 0.75 * _fwhm(exposure.getPsf(), exposure.getBBox().getCenter())
    except Exception:
        return 5.0


def _markSource(ax, x, y, color, radius, label=None):
    """Draw a circle (and optional label) at a source."""
    ax.add_patch(patches.Circle((x, y), radius, fill=False, ec=color, lw=1.8,
                                path_effects=_OUTLINE, zorder=6))
    if label is not None:
        ax.text(x + 0.75 * radius, y + 0.75 * radius, f"#{label}", color=color,
                fontsize=8, fontweight='bold', path_effects=_TEXT_OUTLINE, zorder=7)


def _frameAxes(ax, color, lw=3.0):
    """Colour the frame of an axes."""
    for s in ax.spines.values():
        s.set_edgecolor(color)
        s.set_linewidth(lw)


# =====================================================================
# Diagnostic joint fit
# =====================================================================

def jointPsfFitDiag(exposure, psf, peakPositions, badBitmask,
                       ampZeroTol=1e-6, supportFrac=0.01,
                       varModel=None, sysFrac=0.03, poissonK=0.0):
    """Fit all PSF amplitudes in one non-negative solve, keeping what is needed for diagnostics.

    Weighted, non-negative, one sparse solve against the original exposure.  PSF stamps are not
    clipped at zero (clipping would add a positive pedestal).

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Image to fit.
    psf : `lsst.afw.detection.Psf`
        PSF model.
    peakPositions : `list` [`tuple` [`float`, `float`]]
        Source positions in parent pixel coordinates.
    badBitmask : `int`
        Mask bits of pixels to ignore.
    ampZeroTol : `float`, optional
        Amplitudes at or below this count as zero.
    supportFrac : `float`, optional
        A source's significant support is where its stamp exceeds this fraction of its peak.
    varModel : `numpy.ndarray`, optional
        Model image used for the model-dependent variance terms (see `effectiveVariance`).
    sysFrac : `float`, optional
        Fractional model error added to the variance.
    poissonK : `float`, optional
        Extra variance per unit model flux.

    Returns
    -------
    fit : `dict`
        ``amps``, ``modelFull``, ``nAtZero``; ``comps`` (per source a dict with ``stamp``,
        ``sig`` (significant support), ``ys``/``xs`` (slices of its fit region), ``amp``,
        ``bboxParent``, or None if the stamp misses the image); coupled-source groups
        (``groupLabels``, ``nGroups``, ``groupMap``, ``groupChi2``); ``good`` (pixels
        used), ``chi2map``, ``variance`` (the effective variance), ``timings`` and the solver
        status (``lsqStatus``, ``lsqNit``, ``lsqCost``).
    """
    expBbox = exposure.getBBox()
    ny, nx = expBbox.getHeight(), expBbox.getWidth()
    data = exposure.image.array
    variance = effectiveVariance(exposure, varModel, sysFrac, poissonK)
    mask = exposure.mask.array

    good = ((mask & badBitmask) == 0) & np.isfinite(data) & np.isfinite(variance) & (variance > 0)
    pixelIndex = -np.ones((ny, nx), dtype=np.int64)
    flat = np.flatnonzero(good.ravel())
    pixelIndex.ravel()[flat] = np.arange(len(flat))
    nPix = len(flat)
    dVec = data.ravel()[flat]
    wVec = 1.0 / variance.ravel()[flat]
    nSrc = len(peakPositions)

    timings = {}
    t0 = time.perf_counter()
    comps = []
    rL, cL, vL, brL, bcL = [], [], [], [], []
    for j, (x, y) in enumerate(peakPositions):
        stamp, bbox = psfStamp(psf, x, y, expBbox, taper=WING_TAPER)
        if stamp is None:
            comps.append(None)
            continue
        ys_, xs_ = bboxToSlices(bbox, expBbox)
        subIdx = pixelIndex[ys_, xs_]
        smax = stamp.max()
        sig = stamp > supportFrac * smax if smax > 0 else np.zeros_like(stamp, bool)
        ok = (subIdx >= 0) & (stamp != 0)
        okb = ok & sig
        rL.append(subIdx[ok]); cL.append(np.full(int(ok.sum()), j)); vL.append(stamp[ok])
        brL.append(subIdx[okb]); bcL.append(np.full(int(okb.sum()), j))
        comps.append(dict(j=j, amp=0.0, stamp=stamp, sig=sig, ys=ys_, xs=xs_,
                          bboxParent=(bbox.getMinX(), bbox.getMaxX(),
                                       bbox.getMinY(), bbox.getMaxY())))

    cat = lambda l, dt: np.concatenate(l) if l else np.array([], dtype=dt)
    rows, cols, vals = cat(rL, np.int64), cat(cL, np.int64), cat(vL, float)
    A = coo_matrix((vals, (rows, cols)), shape=(nPix, nSrc)).tocsr()
    sqrtW = np.sqrt(wVec)
    Aw = csr_matrix(A.multiply(sqrtW[:, None]))
    dw = dVec * sqrtW
    timings['joint fit: build A'] = time.perf_counter() - t0

    t0 = time.perf_counter()
    res = lsq_linear(Aw, dw, bounds=(0, np.inf), method='trf', lsmr_tol='auto', max_iter=2000)
    amps = res.x
    timings['joint fit: lsq_linear'] = time.perf_counter() - t0

    t0 = time.perf_counter()
    modelFull = np.zeros((ny, nx))
    for j, c in enumerate(comps):
        if c is None:
            continue
        c['amp'] = float(amps[j])
        if amps[j] > ampZeroTol:
            modelFull[c['ys'], c['xs']] += amps[j] * c['stamp']

    br, bc = cat(brL, np.int64), cat(bcL, np.int64)
    B = coo_matrix((np.ones(len(br)), (br, bc)), shape=(nPix, nSrc)).tocsr()
    nGroups, labels = connected_components(B.T @ B, directed=False)
    groupMap = -np.ones((ny, nx), dtype=np.int64)
    for j, c in enumerate(comps):
        if c is None:
            continue
        sub = groupMap[c['ys'], c['xs']]
        sub[c['sig']] = labels[j]

    resid = data - modelFull
    chi2map = np.zeros((ny, nx))
    chi2map[good] = resid[good] ** 2 / variance[good]
    groupChi2 = {}
    for g in range(nGroups):
        m = (groupMap == g) & good
        groupChi2[g] = (chi2map[m].sum() / m.sum()) if m.sum() > 0 else np.nan
    timings['joint fit: model + groups'] = time.perf_counter() - t0

    return dict(amps=amps, modelFull=modelFull, nAtZero=int(np.sum(amps <= ampZeroTol)),
                comps=comps, groupLabels=labels, nGroups=nGroups, groupMap=groupMap,
                groupChi2=groupChi2, good=good, chi2map=chi2map, timings=timings,
                lsqStatus=res.status, lsqNit=res.nit, lsqCost=res.cost,
                variance=variance)


# =====================================================================
# Centroid refinement (fixes centroid-offset "dipole" residuals)
# =====================================================================

def computePsfDerivatives(psf, x, y, bboxClip, eps=0.25):
    """Compute central-difference derivative stamps dP/dx and dP/dy.

    ``computeImage`` centres its box on the pixel nearest the requested point, so the stamps at
    (x +- eps, y) and (x, y +- eps) need not share a box; all five boxes (and the clip box) are
    intersected first.

    Parameters
    ----------
    psf : `lsst.afw.detection.Psf`
        PSF model.
    x, y : `float`
        Position in parent pixel coordinates.
    bboxClip : `lsst.geom.Box2I`
        Box to clip to.
    eps : `float`, optional
        Finite-difference step [px].

    Returns
    -------
    dpdx, dpdy : `numpy.ndarray` or `None`
        Derivative stamps over ``bbox``, or None if the boxes do not overlap.
    bbox : `lsst.geom.Box2I`
        The common box.
    """
    imgC = psf.computeImage(geom.Point2D(x, y))
    imgXp = psf.computeImage(geom.Point2D(x + eps, y))
    imgXm = psf.computeImage(geom.Point2D(x - eps, y))
    imgYp = psf.computeImage(geom.Point2D(x, y + eps))
    imgYm = psf.computeImage(geom.Point2D(x, y - eps))

    common = geom.Box2I(imgC.getBBox())
    for im in (imgXp, imgXm, imgYp, imgYm):
        common.clip(im.getBBox())
    common.clip(bboxClip)
    if common.isEmpty():
        return None, None, common

    dpdx = (imgXp[common].array - imgXm[common].array) / (2 * eps)
    dpdy = (imgYp[common].array - imgYm[common].array) / (2 * eps)
    return dpdx, dpdy, common


def windowedCentroid(imageArray, xy0, x, y, halfWidth, maxFracOfWindow=0.9):
    """Return the flux-weighted centroid in a box around (x, y).

    Negative pixels are clipped to zero.  Falls back to (x, y) if there is too little positive
    flux or the centroid lands too close to the window edge (dominated by a neighbour).

    Parameters
    ----------
    imageArray : `numpy.ndarray`
        Image.
    xy0 : `tuple` [`int`, `int`]
        Parent pixel coordinates of ``imageArray[0, 0]``.
    x, y : `float`
        Starting position in parent pixel coordinates.
    halfWidth : `int`
        Half-size of the box [px].
    maxFracOfWindow : `float`, optional
        Reject centroids further than this fraction of ``halfWidth`` from the start.

    Returns
    -------
    cx, cy : `float`
        Centroid in parent pixel coordinates.
    """
    x0, y0 = xy0
    ny, nx = imageArray.shape
    ix, iy = int(round(x - x0)), int(round(y - y0))
    xa, xb = max(ix - halfWidth, 0), min(ix + halfWidth + 1, nx)
    ya, yb = max(iy - halfWidth, 0), min(iy + halfWidth + 1, ny)
    if xb <= xa or yb <= ya:
        return x, y
    stamp = np.clip(imageArray[ya:yb, xa:xb], 0, None)
    total = stamp.sum()
    if total <= 0:
        return x, y
    yy, xx = np.mgrid[ya:yb, xa:xb]
    cx = float((xx * stamp).sum() / total) + x0
    cy = float((yy * stamp).sum() / total) + y0
    if abs(cx - x) > maxFracOfWindow * halfWidth or \
       abs(cy - y) > maxFracOfWindow * halfWidth:
        return x, y
    return cx, cy


def refinePeakPositions(exposure, psf, peaks, amps, badBitmask,
                          minAmpForRefine=50.0, maxShift=0.6, nSubiters=3,
                          minGoodPixels=10, coreRadiusSigma=2.5,
                          variance=None, freezeSep=None):
    """Move each well-detected source towards the position that minimises its residual.

    Linearised Gauss-Newton with all other sources held fixed: over the source's core pixels,
    weighted by 1/sigma, ``(data - model) ~ dA * P + a * dx * dP/dx + a * dy * dP/dy``.  The dA
    column absorbs amplitude errors so they do not leak into the position; only the position
    step is applied and the caller re-solves the amplitudes.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    psf : `lsst.afw.detection.Psf`
        PSF model.
    peaks : `list` [`tuple` [`float`, `float`]]
        Source positions.
    amps : array-like
        Source amplitudes.
    badBitmask : `int`
        Mask bits of pixels to ignore.
    minAmpForRefine : `float`, optional
        Fainter sources keep their position.
    maxShift : `float`, optional
        Maximum shift per sub-step [px].
    nSubiters : `int`, optional
        Number of Gauss-Newton sub-steps.
    minGoodPixels : `int`, optional
        Minimum number of good core pixels.
    coreRadiusSigma : `float`, optional
        Core radius in units of the PSF sigma.
    variance : `numpy.ndarray`, optional
        Effective variance (default: the exposure's variance plane).
    freezeSep : `float`, optional
        Sources with a brighter neighbour within this many px keep their position (it is
        degenerate with the neighbour's amplitude).

    Returns
    -------
    peaks : `list` [`tuple` [`float`, `float`]]
        Refined positions.
    """
    expBbox = exposure.getBBox()
    data, mask = exposure.image.array, exposure.mask.array
    if variance is None:
        variance = exposure.variance.array
    good = ((mask & badBitmask) == 0) & np.isfinite(data) & np.isfinite(variance) & (variance > 0)
    peaks = [tuple(p) for p in peaks]
    amps = np.asarray(amps, float).copy()

    frozen = np.zeros(len(peaks), bool)
    if freezeSep is not None and len(peaks) > 1:
        P0 = np.asarray(peaks)
        for j, nb in enumerate(cKDTree(P0).query_ball_point(P0, freezeSep)):
            frozen[j] = any(k != j and amps[k] > amps[j] for k in nb)

    shp = psf.computeShape(expBbox.getCenter())
    R2 = (coreRadiusSigma * np.sqrt(0.5 * (shp.getIxx() + shp.getIyy()))) ** 2

    for _ in range(nSubiters):
        model = np.zeros_like(data)
        for (x, y), a in zip(peaks, amps):
            if a <= 0:
                continue
            stamp, bb = psfStamp(psf, x, y, expBbox, taper=WING_TAPER)
            if stamp is None:
                continue
            ys, xs = bboxToSlices(bb, expBbox)
            model[ys, xs] += a * stamp
        residAll = data - model

        moved = 0.0
        for j, ((x, y), a) in enumerate(zip(peaks, amps)):
            if a < minAmpForRefine or frozen[j]:
                continue
            dpdx, dpdy, bb = computePsfDerivatives(psf, x, y, expBbox)
            if dpdx is None:
                continue
            ys, xs = bboxToSlices(bb, expBbox)
            yy, xx = np.mgrid[bb.getMinY():bb.getMaxY() + 1, bb.getMinX():bb.getMaxX() + 1]
            g = good[ys, xs] & (((xx - x) ** 2 + (yy - y) ** 2) < R2)
            if g.sum() < minGoodPixels:
                continue
            P, _ = psfStamp(psf, x, y, bb, taper=WING_TAPER)
            w = 1.0 / np.sqrt(variance[ys, xs][g])
            J = np.column_stack([P[g], a * dpdx[g], a * dpdy[g]]) * w[:, None]
            r = residAll[ys, xs][g] * w
            try:
                delta, *_ = np.linalg.lstsq(J, r, rcond=None)
            except np.linalg.LinAlgError:
                continue
            dx, dy = np.clip(delta[1:], -maxShift, maxShift)
            peaks[j] = (x + dx, y + dy)
            moved = max(moved, abs(dx), abs(dy))
        if moved < 0.02:
            break
    return peaks


def redetectAndRefitPsfAdaptive(exposure, kernelSizes=(41, 27, 21, 15, 11), **kwargs):
    """Run `redetectAndRefitPsf`, trying kernel sizes from largest to smallest.

    MakePsfCandidatesTask drops every candidate whose stamp touches masked or edge pixels; on a
    small cutout a large kernel can drop all of them.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Background-subtracted exposure.
    kernelSizes : `tuple` [`int`], optional
        Kernel sizes to try [px].
    **kwargs
        Keyword arguments of `redetectAndRefitPsf`.

    Returns
    -------
    psfNew : `lsst.afw.detection.Psf`
        The new PSF (the previous one if every attempt failed).
    sources : `lsst.afw.table.SourceCatalog` or `None`
        Measured detections.
    starCat : `lsst.afw.table.SourceCatalog` or `None`
        Selected PSF stars.
    cellSet : `lsst.afw.math.SpatialCellSet` or `None`
        PSF candidates.
    kernelSize : `int` or `None`
        Kernel size that worked (None if all failed).
    """
    psfPrev = exposure.getPsf()
    fwhmPrev = _fwhm(psfPrev, exposure.getBBox().getCenter())
    for ks in kernelSizes:
        print(f"  [adaptive PSF refit] trying kernelSize={ks}...")
        try:
            psfNew, sourcesNew, starCat_new, cellSet_new = redetectAndRefitPsf(
                exposure, kernelSize=ks, **kwargs)
        except Exception as e:
            print(f"    kernelSize={ks} raised {type(e).__name__}: {e} -- trying smaller")
            continue
        if cellSet_new is not None:
            fwhmNew = _fwhm(psfNew, exposure.getBBox().getCenter())
            nUsed = sum(1 for cell in cellSet_new.getCellList() for c in cell
                         if c.getStatus() == GOOD_STATUS)
            print(f"    SUCCESS at kernelSize={ks}: {nUsed} candidates used, "
                  f"FWHM {fwhmPrev:.3f} -> {fwhmNew:.3f} px")
            exposure.setPsf(psfNew)
            return psfNew, sourcesNew, starCat_new, cellSet_new, ks
        print(f"    kernelSize={ks} failed (no viable candidates) -- trying smaller")

    print(f"  [adaptive PSF refit] FAILED at every kernelSize tried {kernelSizes} -- "
          f"keeping previous PSF (FWHM {fwhmPrev:.3f} px). Consider a larger cutout for "
          f"the PSF-refit step, or refitPsf=False / 'stack'.")
    exposure.setPsf(psfPrev)
    return psfPrev, None, None, None, None


# =====================================================================
# Per-iteration figures
# =====================================================================

def plotIterationOverview(exposure, snap, hist, outerIdx, normData,
                            colors=None, viewHalf=20, circleRadius=None):
    """Plot the overview of one inner iteration.

    Detection input, detections, fitted sources, model, residual and the iteration histories.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    snap : `dict`
        Inner-loop snapshot of the iteration.
    hist : `dict`
        Histories (``residStd``, ``nPeaks``, ``nZero``).
    outerIdx : `int`
        Outer-pass index.
    normData : `astropy.visualization.ImageNormalize`
        Display normalisation of the data.
    colors : `dict` [`int`, `str`], optional
        Colours of the sources with per-source pages.
    viewHalf : `int`, optional
        Half-size of the per-source windows [px].
    circleRadius : `float`, optional
        Marker radius [px].
    """
    colors = colors or {}
    radius = _circleRadius(exposure, circleRadius)
    bbox = exposure.getBBox()
    X0, Y0 = bbox.getMinX(), bbox.getMinY()
    ext = _extent(bbox)
    xs, ys = _coords(bbox)
    data = exposure.image.array
    fit = snap['fit']
    model = fit['modelFull']
    resid = snap['residual']
    peaks = snap['mergedArr']
    isNew = snap['isNew']
    bad = snap['badMask']
    it = snap['it']

    fig = plt.figure(figsize=(21, 10.5))
    gs = fig.add_gridspec(2, 4, hspace=0.32, wspace=0.28)

    # 1. detection input
    ax = fig.add_subplot(gs[0, 0])
    detIn = snap['detIn']
    im = ax.imshow(detIn, extent=ext, origin='lower', cmap=CMAP,
                   norm=simple_norm(detIn, stretch='asinh', min_percent=1, max_percent=99.5))
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    if snap['footprintMask'].any():
        ax.contour(xs, ys, snap['footprintMask'].astype(float), levels=[0.5],
                   colors='white', linewidths=0.6)
    raw = _peaksArray(snap['rawPeaks'])
    col = _peaksArray(snap['collapsedPeaks'])
    if len(raw):
        ax.scatter(raw[:, 0], raw[:, 1], marker='+', c='w', s=28,
                   label=f'raw peaks ({len(raw)})')
    if len(col):
        ax.scatter(col[:, 0], col[:, 1], facecolors='none', edgecolors='yellow', s=70,
                   lw=1.0, label=f'after sat-collapse + rejection ({len(col)})')
    _overlayBad(ax, bad, ext)
    ax.legend(fontsize=7, loc='upper right')
    what = 'data (bg-subtracted)' if it == 0 else 'residual from previous fit'
    ax.set_title(f"1. Detection input: {what}\nthreshold = {snap['threshold']:.2f}σ  (white = footprints)",
                 fontsize=10)
    _cleanAxes(ax)

    # 2. fit regions & the cutouts shown on the per-source pages
    ax = fig.add_subplot(gs[0, 1])
    im = ax.imshow(data, extent=ext, origin='lower', norm=normData, cmap=CMAP)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    for c in fit['comps']:
        if c is None or c['j'] in colors:
            continue
        x0b, x1b, y0b, y1b = c['bboxParent']
        ax.add_patch(patches.Rectangle((x0b - 0.5, y0b - 0.5), x1b - x0b + 1, y1b - y0b + 1,
                                       fill=False, ec='0.75', lw=0.4, alpha=0.5))
    labels = fit['groupLabels']
    sizes = np.bincount(labels, minlength=fit['nGroups'])
    for g in range(fit['nGroups']):
        gm = fit['groupMap'] == g
        if gm.any():
            ax.contour(xs, ys, gm.astype(float), levels=[0.5], colors='w',
                       linewidths=0.8, linestyles='--')
    if len(peaks):
        ax.scatter(peaks[~isNew, 0], peaks[~isNew, 1], facecolors='none', edgecolors='w',
                   s=25, lw=0.7, label='existing peak')
        ax.scatter(peaks[isNew, 0], peaks[isNew, 1], marker='*', c='red', s=70,
                   edgecolors='w', lw=0.4, label='NEW this iteration')
    rej = _peaksArray([(x, y) for (x, y, _) in snap.get('shapeRejected', [])])
    if len(rej):
        ax.scatter(rej[:, 0], rej[:, 1], marker='x', c='#ff00ff', s=45, lw=1.4,
                   label=f'removed by shape cut ({len(rej)})')
    for j, colr in colors.items():
        px, py = peaks[j]
        ix, iy, xa, xb, ya, yb = _viewWindow(bbox, px, py, viewHalf)
        ax.add_patch(patches.Rectangle((X0 + xa - 0.5, Y0 + ya - 0.5), xb - xa, yb - ya,
                                       fill=False, ec=colr, lw=1.8, path_effects=_OUTLINE))
        _markSource(ax, px, py, colr, radius, label=j)
    _overlayBad(ax, bad, ext)
    ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
    ax.legend(fontsize=7, loc='upper right')
    ax.set_title(f"2. Coloured box + circle = cutout on the per-source pages (#j)\n"
                 f"gray = other PSF stamps, white dashed = coupled groups ({fit['nGroups']}, "
                 f"max {sizes.max() if len(sizes) else 0} src)", fontsize=10)
    _cleanAxes(ax)

    # 3. model
    ax = fig.add_subplot(gs[0, 2])
    im = ax.imshow(model, extent=ext, origin='lower', norm=normData, cmap=CMAP)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_title(f"3. Joint PSF model ({len(peaks)} sources, "
                 f"{fit['nAtZero']} at zero amp)", fontsize=10)
    _cleanAxes(ax)

    # 4. residual
    ax = fig.add_subplot(gs[0, 3])
    rs = np.nanpercentile(np.abs(resid), 99)
    rs = rs if np.isfinite(rs) and rs > 0 else 1.0
    im = ax.imshow(resid, extent=ext, origin='lower', cmap='RdBu_r', vmin=-rs, vmax=rs)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    for g in range(fit['nGroups']):
        gm = fit['groupMap'] == g
        if gm.any():
            ax.contour(xs, ys, gm.astype(float), levels=[0.5], colors='k',
                       linewidths=0.5, alpha=0.5)
    ax.set_title(f"4. Residual = data - model\nstd {snap['stdBefore']:.3f} → {snap['stdAfter']:.3f}",
                 fontsize=10)
    _cleanAxes(ax)

    # 5. normalized residual histogram (EFFECTIVE variance)
    ax = fig.add_subplot(gs[1, 0])
    good = fit['good']
    varEff = fit['variance']
    z = (resid / np.sqrt(np.where(varEff > 0, varEff, np.inf)))[good]
    ax.hist(z, bins=np.linspace(-6, 6, 61), density=True, color='#4C72B0', alpha=0.7)
    g_ = np.linspace(-6, 6, 200)
    ax.plot(g_, np.exp(-0.5 * g_ ** 2) / np.sqrt(2 * np.pi), 'k--', lw=1, label='N(0,1)')
    ax.set_yscale('log'); ax.set_ylim(1e-4, 1)
    ax.set_xlabel('residual / σ_eff  (good pixels)')
    ax.set_title(f"5. Normalized residuals: mean={np.nanmean(z):.2f}, std={np.nanstd(z):.2f}", fontsize=10)
    ax.legend(fontsize=8)

    # 6. amplitudes
    ax = fig.add_subplot(gs[1, 1])
    amps = fit['amps']
    order = np.argsort(-amps)
    colorsBar = np.where(isNew[order], '#d62728', '#7f7f7f')
    ax.bar(np.arange(len(amps)), amps[order], color=colorsBar)
    zero = amps[order] <= 1e-6
    if zero.any():
        ax.scatter(np.arange(len(amps))[zero], np.zeros(zero.sum()), marker='x', c='k', s=30,
                   label='fit hit 0 bound', zorder=5)
        ax.legend(fontsize=8)
    ax.set_yscale('symlog', linthresh=1.0)
    ax.set_xlabel('source rank'); ax.set_ylabel('fitted amplitude (flux)')
    ax.set_title("6. Fitted amplitudes (red = new)", fontsize=10)

    # 7. history
    ax = fig.add_subplot(gs[1, 2])
    xsH = np.arange(len(hist['residStd']))
    ax.plot(xsH, hist['residStd'], 'o-', color='#DD8452')
    ax.set_xlabel('iteration (0 = before any fit)'); ax.set_ylabel('residual std', color='#DD8452')
    ax2 = ax.twinx()
    ax2.plot(xsH, hist['nPeaks'], 's--', color='#4C72B0')
    ax2.plot(xsH, hist['nZero'], '^:', color='k')
    ax2.set_ylabel('# peaks (blue) / # at zero (black)')
    ax2.grid(False)
    ax.set_title("7. Convergence history", fontsize=10)

    # 8. timing
    ax = fig.add_subplot(gs[1, 3])
    labs = list(snap['times'].keys())
    vals = [snap['times'][k] for k in labs]
    ax.barh(labs[::-1], vals[::-1], color='#55A868')
    ax.set_xlabel('seconds')
    ax.set_title(f"8. Timing this iteration (total {sum(vals):.2f}s)\n"
                 f"lsq_linear: {fit['lsqNit']} its, status {fit['lsqStatus']}", fontsize=10)
    ax.grid(axis='x', alpha=0.3)

    fig.suptitle(f"Outer pass {outerIdx + 1}  |  inner iteration {it + 1}  |  "
                 f"{len(peaks)} sources ({int(isNew.sum())} new)",
                 fontsize=14, fontweight='medium', y=0.995)
    plt.show()
    plt.close(fig)


# =====================================================================
# Instrumented inner loop
# =====================================================================

def runInnerLoopVerbose(exposure, psf, badBitmask, outerIdx=0,
                           nIterations=10, thresholdStart=8.0, thresholdStep=-1.5,
                           thresholdMin=4.5, minPixels=1, minNewPeaks=1,
                           convergenceFrac=0.02, viz=True, plotEvery=1,
                           sourceView='all', sourcesPerPage=6, maxPagesPerIter=None,
                           viewHalf=20, supportFrac=0.01, refineCentroids=True,
                           centroidMinAmp=50.0, nCentroidRounds=4,
                           centroidShiftTol=0.03, centroidShiftSchedule=None,
                           coarseCentroid=True, coarseCentroidHwSigma=2.5,
                           mergeSeparationSigma=1.5,
                           sysFrac=0.03, minSnr=5.0, pruneSepSigma=1.5,
                           rejectFwhm=0.75, brightPeakSigma=100.0, freezeFwhm=1.0, poissonK=0.0,
                           nearBrightFwhm=2.0, nearBrightSnr=10.0,
                           edgeMarginFwhm=0.5, useEliminateWeak=True,
                           useShapeCut=True, shapeKwargs=None,
                           timer=TIMER, verbose=True, initialPeaks=None):
    """Run the instrumented inner loop (detect, fit, refine, prune, repeat).

    An instrumented version of `iterativeJointDetection` with centroid refinement and guards
    against over-counting sources.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Image to fit.
    psf : `lsst.afw.detection.Psf`
        PSF model.
    badBitmask : `int`
        Mask bits of pixels to ignore.
    outerIdx : `int`, optional
        Outer-pass index (for labels).
    nIterations : `int`, optional
        Maximum number of iterations.
    thresholdStart, thresholdStep, thresholdMin : `float`, optional
        Detection threshold schedule [sigma].
    minPixels : `int`, optional
        Minimum footprint size.
    minNewPeaks : `int`, optional
        Stop when an iteration adds fewer peaks than this.
    convergenceFrac : `float`, optional
        Stop when the residual standard deviation improves by less than this fraction.
    viz : `bool`, optional
        Draw diagnostics.
    plotEvery : `int`, optional
        Draw every n-th iteration (the last always).
    sourceView : `str`, optional
        See `selectSourcePages`.
    sourcesPerPage : `int`, optional
        See `selectSourcePages`.
    maxPagesPerIter : `int`, optional
        Maximum per-source pages per iteration.
    viewHalf : `int`, optional
        Half-size of the per-source windows [px].
    supportFrac : `float`, optional
        See `jointPsfFitDiag`.
    refineCentroids : `bool`, optional
        Refine positions (`refinePeakPositions`).
    centroidMinAmp : `float`, optional
        Fainter sources keep their position.
    nCentroidRounds : `int`, optional
        Maximum centroid rounds per iteration.
    centroidShiftTol : `float`, optional
        Stop centroiding when the largest shift is below this [px].
    centroidShiftSchedule : `list` [`float`], optional
        Maximum shift per round [px].
    coarseCentroid : `bool`, optional
        Start new peaks from a windowed centroid.
    coarseCentroidHwSigma : `float`, optional
        Half-width of that window in PSF sigma.
    mergeSeparationSigma : `float`, optional
        Merge new peaks closer than this many PSF sigma.
    sysFrac : `float`, optional
        Fractional PSF-shape error added in quadrature to the variance.
    minSnr : `float`, optional
        S/N floor (cheap version in `pruneAndMerge`, honest one in `eliminateWeak`).
    pruneSepSigma : `float`, optional
        Pairs closer than this many PSF sigma are merged after refinement.
    rejectFwhm : `float`, optional
        New peaks within this many FWHM of a bright source are discarded.
    brightPeakSigma : `float`, optional
        A source is bright if its peak exceeds this many noise sigma.
    freezeFwhm : `float`, optional
        Sources with a brighter neighbour within this many FWHM keep their position.
    poissonK : `float`, optional
        Extra variance per unit model flux.
    nearBrightFwhm : `float`, optional
        Radius around bright stars with the stricter S/N cut [FWHM].
    nearBrightSnr : `float`, optional
        S/N required there.
    edgeMarginFwhm : `float`, optional
        Drop sources closer than this many FWHM to the edge.
    useEliminateWeak : `bool`, optional
        Backward elimination on the honest S/N.
    useShapeCut : `bool`, optional
        Backward elimination on PSF-likeness (`eliminateNonpsf`).
    shapeKwargs : `dict`, optional
        Keyword arguments of `shapeMetrics` / `judgeShape`.
    timer : `StepTimer`, optional
        Timer collecting the step durations.
    verbose : `bool`, optional
        Print progress.
    initialPeaks : `list` [`tuple` [`float`, `float`]], optional
        Warm-start source positions.

    Returns
    -------
    result : `dict`
        ``modelFull``, ``residualFull``, ``peaks``, ``amps``, the histories
        (``residStdHistory``, ``nPeaksHistory``, ``nZeroAmpHistory``), ``firstSeen`` and
        ``snaps`` (one dict per iteration with everything needed for plotting).
    """
    expBbox = exposure.getBBox()
    ny, nx = expBbox.getHeight(), expBbox.getWidth()
    normData = _dataNorm(exposure.image.array)
    badMask = (exposure.mask.array & badBitmask) != 0

    _sh = psf.computeShape(expBbox.getCenter())
    psfSigma = np.sqrt(0.5 * (_sh.getIxx() + _sh.getIyy()))
    fwhmPx = 2.355 * psfSigma
    psfPeak = psf.computeKernelImage(expBbox.getCenter()).array.max()
    noise = np.sqrt(np.nanmedian(exposure.variance.array))

    workingResidual = exposure.clone()
    workingResidual.setPsf(psf)

    allPeaks, prevAmps = [], np.array([])
    hist = dict(residStd=[np.nanstd(workingResidual.image.array)], nPeaks=[0], nZero=[0])
    firstSeen, snaps = {}, []
    threshold = thresholdStart
    modelFull = np.zeros((ny, nx))
    ampsFinal = np.array([])

    if initialPeaks is not None and len(initialPeaks):
        # ---- warm start (adaptive controller): fit the given peaks, then keep detecting ----
        allPeaks = [tuple(p) for p in initialPeaks]
        _f0 = jointPsfFitDiag(exposure, psf, allPeaks, badBitmask,
                                 supportFrac=supportFrac, varModel=None,
                                 sysFrac=sysFrac, poissonK=poissonK)
        _f0 = jointPsfFitDiag(exposure, psf, allPeaks, badBitmask,
                                 supportFrac=supportFrac, varModel=_f0['modelFull'],
                                 sysFrac=sysFrac, poissonK=poissonK)
        modelFull, ampsFinal, prevAmps = _f0['modelFull'], _f0['amps'], _f0['amps']
        workingResidual.image.array[:] = exposure.image.array - modelFull
        workingResidual.variance.array[:] = _f0['variance']
        hist = dict(residStd=[np.nanstd(workingResidual.image.array)],
                    nPeaks=[len(allPeaks)], nZero=[_f0['nAtZero']])
        print(f"  warm start from {len(allPeaks)} peaks "
              f"(residual std {hist['residStd'][0]:.3f})")

    for it in range(nIterations):
        print(f"\n--- outer {outerIdx + 1} | inner {it + 1}/{nIterations} | "
              f"threshold={threshold:.2f}σ | residual std={hist['residStd'][-1]:.3f} "
              f"| peaks so far={len(allPeaks)} ---")
        times = {}

        def S(label):
            return timer.step(label, outer=outerIdx, inner=it)

        detIn = workingResidual.image.array.copy()

        # --- detection on the current residual
        with S('detect on residual') as rec:
            cfg = SourceDetectionTask.ConfigClass()
            cfg.thresholdValue = threshold
            cfg.thresholdType = "pixel_stdev"
            cfg.minPixels = minPixels
            cfg.doTempLocalBackground = True
            cfg.tempLocalBackground.binSize = 16
            cfg.reEstimateBackground = False
            schemaIt = afwTable.SourceTable.makeMinimalSchema()
            taskIt = SourceDetectionTask(schema=schemaIt, config=cfg)
            tableIt = afwTable.SourceTable.make(schemaIt)
            try:
                det = taskIt.run(tableIt, workingResidual)
            except TooManyMaskedPixelsError:
                print("  background estimation failed on residual -- stopping")
                det = None
        if det is None:
            break
        times['detect on residual'] = rec['seconds']
        footprints = det.positive.getFootprints()
        rawPeaks = [(pk.getFx(), pk.getFy()) for fp in footprints for pk in fp.getPeaks()]
        fpMask = _footprintsToMask(footprints, exposure)
        print(f"  {len(rawPeaks)} raw peaks detected")

        with S('collapse saturated peaks') as rec:
            newPeaks, _ = collapseSaturatedPeaks(rawPeaks, exposure)
        times['collapse saturated peaks'] = rec['seconds']

        # reject new peaks sitting in the PSF-mismatch zone of a bright source
        if len(allPeaks) and len(newPeaks):
            with S('reject near bright') as rec:
                isBright = np.asarray(prevAmps) * psfPeak > brightPeakSigma * noise
                newPeaks, nRej = rejectNearBright(newPeaks, allPeaks, isBright,
                                                      rejectFwhm * fwhmPx)
            times['reject near bright'] = rec['seconds']
            if nRej:
                print(f"  rejected {nRej} new peak(s) within {rejectFwhm * fwhmPx:.1f}px "
                      f"of a bright source")

        if len(newPeaks) < minNewPeaks and it > 0:
            print("  too few new peaks -- treating as converged, stopping")
            break

        if coarseCentroid and len(newPeaks) > 0:
            with S('coarse centroid (new peaks)') as rec:
                hw = max(int(round(coarseCentroidHwSigma * psfSigma)), 3)
                xy0Local = (expBbox.getMinX(), expBbox.getMinY())
                newPeaks = [windowedCentroid(detIn, xy0Local, x, y, hw)
                             for (x, y) in newPeaks]
            times['coarse centroid (new peaks)'] = rec['seconds']

        with S('drop peaks on masked pixels') as rec:
            newPeaks, nDropped = dropMaskedPeaks(newPeaks, exposure)
        times['drop masked peaks'] = rec['seconds']
        if nDropped:
            print(f"  dropped {nDropped} new peak(s) sitting on EDGE/NO_DATA pixels")

        with S('merge close peaks') as rec:
            mergeMinSep = None
            if mergeSeparationSigma is not None:
                mergeMinSep = mergeSeparationSigma * psfSigma
            merged = mergeClosePeaks(allPeaks + newPeaks, exposure,
                                       minSeparation=mergeMinSep)
        times['merge close peaks'] = rec['seconds']
        nIn = len(allPeaks) + len(newPeaks)
        if mergeMinSep is not None:
            print(f"  mergeClosePeaks: {nIn} -> {len(merged)} peaks "
                  f"(min separation {mergeMinSep:.1f} px)")
        if len(merged) == 0:
            print("  no peaks at all -- stopping")
            break

        # varModel=modelFull (previous model) => systematic variance floor
        fit = jointPsfFitDiag(exposure, psf, merged, badBitmask,
                                 supportFrac=supportFrac,
                                 varModel=modelFull, sysFrac=sysFrac,
                                 poissonK=poissonK)
        for lab, sec in fit['timings'].items():
            timer.add(lab, sec, outer=outerIdx, inner=it)
            times[lab] = sec

        if refineCentroids:
            # alternate (refine position -> re-solve amplitude); step cap SHRINKS per round
            schedule = list(centroidShiftSchedule) if centroidShiftSchedule else \
                [1.5, 0.8, 0.4, 0.2, 0.1]
            if len(schedule) < nCentroidRounds:
                schedule += [schedule[-1]] * (nCentroidRounds - len(schedule))

            prevPositions = _peaksArray(merged)
            for roundI in range(nCentroidRounds):
                stepCap = schedule[roundI]
                with S(f'refine centroids (round {roundI + 1})') as rec:
                    merged = refinePeakPositions(
                        exposure, psf, merged, fit['amps'], badBitmask,
                        minAmpForRefine=centroidMinAmp, maxShift=stepCap,
                        variance=fit['variance'], freezeSep=freezeFwhm * fwhmPx)
                times[f'refine centroids r{roundI + 1}'] = rec['seconds']
                with S(f'joint fit: re-solve (round {roundI + 1})') as rec:
                    fit = jointPsfFitDiag(exposure, psf, merged, badBitmask,
                                             supportFrac=supportFrac,
                                             varModel=fit['modelFull'], sysFrac=sysFrac,
                                             poissonK=poissonK)
                for lab, sec in fit['timings'].items():
                    timer.add(lab, sec, outer=outerIdx, inner=it)
                times[f'joint fit: re-solve r{roundI + 1}'] = rec['seconds']

                newPositions = _peaksArray(merged)
                observedShift = (np.max(np.linalg.norm(newPositions - prevPositions, axis=1))
                                  if len(newPositions) else 0.0)
                prevPositions = newPositions
                print(f"    centroid round {roundI + 1} (step cap {stepCap:.2f}px): "
                      f"max position shift = {observedShift:.3f} px")
                if observedShift < centroidShiftTol:
                    break

        # prune low-S/N sources + merge pairs that collapsed together during refinement
        with S('prune + merge') as rec:
            brightFlags = fit['amps'] * psfPeak > brightPeakSigma * noise
            merged = pruneAndMerge(fit, merged, minSnr=minSnr,
                                     minSep=pruneSepSigma * psfSigma,
                                     brightFlags=brightFlags,
                                     nearBrightSep=nearBrightFwhm * fwhmPx,
                                     nearBrightSnr=nearBrightSnr)
        times['prune + merge'] = rec['seconds']

        # Must come AFTER pruneAndMerge (which indexes fit['comps'] by the pre-drop list).
        merged, nOob = dropOutOfBounds(merged, exposure)
        if nOob:
            print(f"  dropped {nOob} source(s) that drifted outside the image")

        if edgeMarginFwhm:
            merged, nEdge = dropOutOfBounds(merged, exposure,
                                                margin=-edgeMarginFwhm * fwhmPx)
            if nEdge:
                print(f"  dropped {nEdge} source(s) within {edgeMarginFwhm} FWHM "
                      f"({edgeMarginFwhm * fwhmPx:.1f}px) of the edge")

        if len(merged) == 0:
            print("  everything pruned/dropped -- stopping")
            break

        # fresh solve on the pruned list (fit must match `merged` before elimination)
        with S('joint fit: re-solve (post-prune)') as rec:
            fit = jointPsfFitDiag(exposure, psf, merged, badBitmask,
                                     supportFrac=supportFrac,
                                     varModel=fit['modelFull'], sysFrac=sysFrac,
                                     poissonK=poissonK)
        for lab, sec in fit['timings'].items():
            timer.add(lab, sec, outer=outerIdx, inner=it)
        times['joint fit: re-solve (post-prune)'] = rec['seconds']

        # backward elimination on the honest S/N
        if useEliminateWeak:
            with S('eliminate weak (honest S/N)') as rec:
                merged, fit = eliminateWeak(exposure, psf, merged, fit, badBitmask,
                                             minSnr=minSnr, supportFrac=supportFrac,
                                             sysFrac=sysFrac, poissonK=poissonK)
            times['eliminate weak'] = rec['seconds']

        # PSF-shape constraint (local-plane S/N, chi2, width test)
        shapeRejected = []
        if useShapeCut:
            with S('shape cut (PSF-likeness)') as rec:
                merged, fit, shapeRejected = eliminateNonpsf(
                    exposure, psf, merged, fit, badBitmask,
                    shapeKwargs=shapeKwargs, supportFrac=supportFrac,
                    sysFrac=sysFrac, poissonK=poissonK)
            times['shape cut'] = rec['seconds']

        modelFull, amps = fit['modelFull'], fit['amps']

        with S('rebuild residual') as rec:
            workingResidual = exposure.clone()
            workingResidual.image.array[:] = exposure.image.array - modelFull
            # detector thresholds on the EFFECTIVE variance
            workingResidual.variance.array[:] = fit['variance']
            workingResidual.setPsf(psf)
        times['rebuild residual'] = rec['seconds']

        mergedArr = _peaksArray(merged)
        priorArr = _peaksArray(allPeaks)
        if len(priorArr):
            d, _ = cKDTree(priorArr).query(mergedArr)
            isNew = d > 1.0
        else:
            isNew = np.ones(len(mergedArr), bool)
        prevAmp = _matchPrevAmps(mergedArr, priorArr, prevAmps)

        for (x, y) in merged:
            firstSeen.setdefault((round(x, 1), round(y, 1)), it)

        oldStd = hist['residStd'][-1]
        newStd = np.nanstd(workingResidual.image.array)
        allPeaks, prevAmps, ampsFinal = merged, amps, amps
        hist['residStd'].append(newStd)
        hist['nPeaks'].append(len(merged))
        hist['nZero'].append(fit['nAtZero'])

        sizes = np.bincount(fit['groupLabels'])
        print(f"  refit {len(merged)} peaks ({int(isNew.sum())} new, {fit['nAtZero']} at zero amp) "
              f"-> residual std {oldStd:.3f} -> {newStd:.3f}")
        print(f"  {fit['nGroups']} coupled groups; sizes: max={sizes.max()}, "
              f"median={np.median(sizes):.0f}, singletons={(sizes == 1).sum()}")
        print(f"  time: " + ", ".join(f"{k}={v:.2f}s" for k, v in times.items()))

        snap = dict(it=it, threshold=threshold, detIn=detIn, rawPeaks=rawPeaks,
                    collapsedPeaks=newPeaks, mergedPeaks=merged, mergedArr=mergedArr,
                    isNew=isNew, prevAmp=prevAmp, footprintMask=fpMask, fit=fit,
                    residual=workingResidual.image.array.copy(), stdBefore=oldStd,
                    stdAfter=newStd, times=times, badMask=badMask,
                    shapeRejected=shapeRejected)
        snaps.append(snap)

        isLastPlanned = (it == nIterations - 1)
        shouldPlot = viz and (it % max(plotEvery, 1) == 0 or isLastPlanned)
        if shouldPlot:
            with timer.step('plotting', outer=outerIdx, inner=it, kind='viz'):
                pages = selectSourcePages(snap, sourceView, sourcesPerPage,
                                            maxPagesPerIter)
                colors = sourceColors(pages)
                plotIterationOverview(exposure, snap, hist, outerIdx, normData,
                                        colors=colors, viewHalf=viewHalf)
                plotSourcePages(exposure, snap, outerIdx, pages=pages, colors=colors,
                                  sourceView=sourceView, viewHalf=viewHalf)

        if oldStd > 0 and it > 0:
            impr = abs(oldStd - newStd) / oldStd
            if impr < convergenceFrac:
                print(f"  residual std improved by only {impr * 100:.2f}% -- converged, stopping")
                break
        threshold = max(threshold + thresholdStep, thresholdMin)

    return dict(modelFull=modelFull, residualFull=workingResidual.image.array.copy(),
                peaks=allPeaks, amps=ampsFinal, residStdHistory=hist['residStd'],
                nPeaksHistory=hist['nPeaks'], nZeroAmpHistory=hist['nZero'],
                firstSeen=firstSeen, snaps=snaps)


# =====================================================================
# Cross-iteration summaries
# =====================================================================

def plotAmplitudeEvolution(result, outerIdx=0):
    """Plot the amplitude of every final source against inner iteration.

    Parameters
    ----------
    result : `dict`
        Inner-loop result with ``snaps``.
    outerIdx : `int`, optional
        Outer-pass index (for the title).
    """
    snaps = result['snaps']
    if not snaps:
        return
    ref = snaps[-1]['mergedArr']
    ampMat = np.full((len(snaps), len(ref)), np.nan)
    for i, s in enumerate(snaps):
        if len(s['mergedArr']) == 0:
            continue
        d, k = cKDTree(s['mergedArr']).query(ref)
        ok = d < 2.0
        ampMat[i, ok] = s['fit']['amps'][k[ok]]

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(16, 5.5))
    final = np.nan_to_num(ampMat[-1], nan=0.0)
    cols = plt.cm.viridis(np.argsort(np.argsort(final)) / max(len(final) - 1, 1))
    for k in range(len(ref)):
        ax0.plot(np.arange(1, len(snaps) + 1), ampMat[:, k], 'o-', color=cols[k], ms=3, lw=1)
    ax0.set_yscale('symlog', linthresh=1.0)
    ax0.set_xlabel('inner iteration'); ax0.set_ylabel('fitted amplitude')
    ax0.set_title(f'Outer {outerIdx + 1}: amplitude of each final source vs iteration')

    with np.errstate(divide='ignore', invalid='ignore'):
        frac = (ampMat[1:] - ampMat[:-1]) / np.abs(ampMat[:-1]) if len(snaps) > 1 else None
    if frac is not None:
        for k in range(len(ref)):
            ax1.plot(np.arange(2, len(snaps) + 1), 100 * frac[:, k], '.-', color=cols[k], lw=0.8)
        ax1.axhline(0, color='k', lw=0.6)
        ax1.set_ylim(-100, 100)
    ax1.set_xlabel('inner iteration'); ax1.set_ylabel('% change vs previous iteration')
    ax1.set_title('How much each source moved when refit')
    plt.tight_layout()
    plt.show()
    plt.close(fig)


def plotSourceEvolution(exposure, result, outerIdx=0, viewHalf=20):
    """Plot model and residual cross-sections of three sources at every iteration.

    The brightest, median and faintest final sources: the x cross-section of the model (top)
    and of the residual (bottom).

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    result : `dict`
        Inner-loop result with ``snaps``.
    outerIdx : `int`, optional
        Outer-pass index (for the title).
    viewHalf : `int`, optional
        Half-length of the cross-sections [px].
    """
    snaps = result['snaps']
    if not snaps:
        return
    bbox = exposure.getBBox()
    X0, Y0 = bbox.getMinX(), bbox.getMinY()
    data = exposure.image.array
    final = snaps[-1]
    peaks, amps = final['mergedArr'], final['fit']['amps']
    order = np.argsort(amps)
    picks = {'dimmest': order[0], 'median': order[len(order) // 2], 'brightest': order[-1]}
    cols = plt.cm.viridis(np.linspace(0, 1, len(snaps)))

    fig, axes = plt.subplots(2, 3, figsize=(18, 8), sharex='col')
    for c, (name, j) in enumerate(picks.items()):
        px, py = peaks[j]
        iy = int(np.clip(int(round(py)) - Y0, 0, data.shape[0] - 1))
        ix = int(np.clip(int(round(px)) - X0, 0, data.shape[1] - 1))
        xa, xb = max(ix - viewHalf, 0), min(ix + viewHalf + 1, data.shape[1])
        xx = X0 + np.arange(xa, xb)
        axes[0, c].plot(xx, data[iy, xa:xb], 'k', lw=2, label='data')
        for i, s in enumerate(snaps):
            axes[0, c].plot(xx, s['fit']['modelFull'][iy, xa:xb], color=cols[i], lw=1.2,
                            label=f'iter {i + 1}')
            axes[1, c].plot(xx, data[iy, xa:xb] - s['fit']['modelFull'][iy, xa:xb],
                            color=cols[i], lw=1.0)
        axes[1, c].axhline(0, color='gray', lw=0.6)
        for a in axes[:, c]:
            a.axvline(px, color='red', lw=0.7, alpha=0.6)
        axes[0, c].set_title(f'{name}: ({px:.1f},{py:.1f}) final amp={amps[j]:.1f}')
        axes[1, c].set_xlabel('x [px]')
    axes[0, 0].set_ylabel('data / model'); axes[1, 0].set_ylabel('residual')
    axes[0, 0].legend(fontsize=6, ncol=2)
    fig.suptitle(f'Outer {outerIdx + 1}: model & residual through the source row, all iterations',
                 fontsize=13)
    plt.tight_layout()
    plt.show()
    plt.close(fig)


# =====================================================================
# Outer-loop step figures
# =====================================================================

def plotCutoutSetup(cutout, psf, badBitmask, starPositions=None):
    """Plot the cutout before fitting: data, sigma map, mask planes and PSF.

    Parameters
    ----------
    cutout : `lsst.afw.image.Exposure`
        The cutout.
    psf : `lsst.afw.detection.Psf`
        PSF model.
    badBitmask : `int`
        Mask bits of bad pixels (shown in red).
    starPositions : array-like, optional
        Original PSF star positions to mark.
    """
    bbox = cutout.getBBox()
    ext = _extent(bbox)
    data = cutout.image.array
    bad = (cutout.mask.array & badBitmask) != 0
    fig, axes = plt.subplots(1, 4, figsize=(21, 5))

    im = axes[0].imshow(data, extent=ext, origin='lower', norm=_dataNorm(data), cmap=CMAP)
    plt.colorbar(im, ax=axes[0], fraction=0.046, pad=0.04)
    _overlayBad(axes[0], bad, ext)
    if starPositions is not None and len(starPositions):
        sp = _peaksArray(starPositions)
        axes[0].scatter(sp[:, 0], sp[:, 1], facecolors='none', edgecolors='lime', s=80,
                        label='original PSF stars')
        axes[0].legend(fontsize=8)
    axes[0].set_title(f'Cutout data ({bbox.getWidth()}x{bbox.getHeight()} px), red = masked')

    im = axes[1].imshow(np.sqrt(cutout.variance.array), extent=ext, origin='lower', cmap='magma')
    plt.colorbar(im, ax=axes[1], fraction=0.046, pad=0.04)
    axes[1].set_title('σ map (sqrt variance)')

    planes = cutout.mask.getMaskPlaneDict()
    txt = []
    for name in ['SAT', 'CR', 'BAD', 'INTRP', 'EDGE', 'NO_DATA', 'DETECTED']:
        if name in planes:
            n = int(((cutout.mask.array & cutout.mask.getPlaneBitMask(name)) != 0).sum())
            txt.append(f'{name}: {n}')
    axes[2].imshow(bad, extent=ext, origin='lower', cmap='gray_r')
    axes[2].set_title('Pixels excluded from fit\n' + ', '.join(txt), fontsize=9)

    c = bbox.getCenter()
    k = psf.computeKernelImage(c).array
    im = axes[3].imshow(k, origin='lower', cmap='inferno')
    plt.colorbar(im, ax=axes[3], fraction=0.046, pad=0.04)
    axes[3].set_title(f'Starting PSF kernel (FWHM≈{_fwhm(psf, c):.2f}px, {k.shape[0]}x{k.shape[1]})')
    for a in axes:
        _cleanAxes(a)
    plt.tight_layout()
    plt.show()
    plt.close(fig)


def plotOuterStepBackground(rawExp, expNext, srcMask, modelFull, bgNext, bgPrev,
                               outerIdx):
    """Plot the background refit of one outer step.

    Input image, model and source mask, new and previous background, and the result.

    Parameters
    ----------
    rawExp : `lsst.afw.image.Exposure`
        Input to the background refit.
    expNext : `lsst.afw.image.Exposure`
        Background-subtracted result.
    srcMask : `numpy.ndarray` [`bool`]
        Source mask used.
    modelFull : `numpy.ndarray`
        Joint-fit model.
    bgNext : `lsst.afw.image.Image` or `None`
        New background.
    bgPrev : `lsst.afw.image.Image` or `None`
        Previous background.
    outerIdx : `int`
        Outer-pass index.

    Returns
    -------
    norm : `astropy.visualization.ImageNormalize`
        Display normalisation used for the images.
    """
    bbox = rawExp.getBBox()
    ext = _extent(bbox)
    xs, ys = _coords(bbox)
    before, after = rawExp.image.array, expNext.image.array
    fig, axes = plt.subplots(2, 4, figsize=(21, 9.5))

    def asinh(a):
        return simple_norm(a, stretch='asinh', min_percent=1, max_percent=99.5)

    im = axes[0, 0].imshow(before, extent=ext, origin='lower', cmap=CMAP, norm=asinh(before))
    plt.colorbar(im, ax=axes[0, 0], fraction=0.046, pad=0.04)
    axes[0, 0].set_title('Input to background refit')

    im = axes[0, 1].imshow(modelFull, extent=ext, origin='lower', cmap=CMAP,
                           norm=simple_norm(modelFull, stretch='asinh', percent=99.5))
    axes[0, 1].contour(xs, ys, srcMask.astype(float), levels=[0.5], colors='red', linewidths=0.8)
    plt.colorbar(im, ax=axes[0, 1], fraction=0.046, pad=0.04)
    axes[0, 1].set_title(f'Joint model + source mask (red)\n{100 * srcMask.mean():.1f}% masked')

    bgarr = bgNext.array if bgNext is not None else np.zeros_like(before)
    im = axes[0, 2].imshow(bgarr, extent=ext, origin='lower', cmap='viridis')
    plt.colorbar(im, ax=axes[0, 2], fraction=0.046, pad=0.04)
    axes[0, 2].set_title('Newly fit background (subtracted)')

    im = axes[0, 3].imshow(after, extent=ext, origin='lower', cmap=CMAP, norm=asinh(after))
    plt.colorbar(im, ax=axes[0, 3], fraction=0.046, pad=0.04)
    axes[0, 3].set_title('After background refit')

    cy, cx = before.shape[0] // 2, before.shape[1] // 2
    if bgPrev is not None:
        axes[1, 0].plot(xs, bgPrev[cy, :], label='previous bg')
        axes[1, 1].plot(ys, bgPrev[:, cx], label='previous bg')
    axes[1, 0].plot(xs, bgarr[cy, :], label='new bg')
    axes[1, 1].plot(ys, bgarr[:, cx], label='new bg')
    axes[1, 0].set_title(f'Background x cross-section (row {cy})'); axes[1, 0].legend()
    axes[1, 1].set_title(f'Background y cross-section (col {cx})'); axes[1, 1].legend()
    axes[1, 0].set_xlabel('x [px]'); axes[1, 1].set_xlabel('y [px]')

    axes[1, 2].plot(xs, np.where(srcMask[cy, :], np.nan, before[cy, :]), label='before (sources masked)')
    axes[1, 2].plot(xs, np.where(srcMask[cy, :], np.nan, after[cy, :]), label='after (sources masked)')
    axes[1, 2].axhline(0, color='gray', lw=0.6)
    axes[1, 2].set_title('Sky-only row before/after'); axes[1, 2].legend(fontsize=8)
    axes[1, 2].set_xlabel('x [px]')

    sky = ~srcMask
    lo, hi = np.nanpercentile(before[sky], [1, 99]) if sky.any() else (-1, 1)
    bins = np.linspace(lo, hi, 60)
    axes[1, 3].hist(before[sky], bins=bins, alpha=0.5, label=f'before (med={np.median(before[sky]):.2f})')
    axes[1, 3].hist(after[sky], bins=bins, alpha=0.5, label=f'after (med={np.median(after[sky]):.2f})')
    axes[1, 3].set_title('Sky pixel values'); axes[1, 3].legend(fontsize=8)
    for a in axes[0]:
        _cleanAxes(a)
    fig.suptitle(f'Outer {outerIdx + 1} → {outerIdx + 2}: mask sources & refit background',
                 fontsize=13)
    plt.tight_layout()
    plt.show()
    plt.close(fig)


def plotPsfCandidateStamps(cellSet, maxStamps=24, title=''):
    """Plot the PSF candidate stamps, used ones first.

    Parameters
    ----------
    cellSet : `lsst.afw.math.SpatialCellSet` or `None`
        PSF candidates.
    maxStamps : `int`, optional
        Maximum number of stamps.
    title : `str`, optional
        Title suffix.
    """
    if cellSet is None:
        print("  (no cellSet -- PSF determination did not run, nothing to show)")
        return
    cands = [c for cell in cellSet.getCellList() for c in cell]
    if not cands:
        print("  (cellSet empty)")
        return
    cands = sorted(cands, key=lambda c: c.getStatus() != GOOD_STATUS)[:maxStamps]
    ncol = min(6, len(cands))
    nrow = int(np.ceil(len(cands) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.1 * ncol, 3.1 * nrow), squeeze=False)
    for a in axes.ravel():
        a.axis('off')
    for a, c in zip(axes.ravel(), cands):
        st = c.getMaskedImage().image.array
        a.imshow(st, origin='lower', cmap='plasma',
                 norm=simple_norm(np.nan_to_num(st), percent=99.5, stretch='log'))
        src = c.getSource()
        good = c.getStatus() == GOOD_STATUS
        a.set_title(f"({src.getX():.0f},{src.getY():.0f}) {'USED' if good else 'rejected'}",
                    fontsize=8, color='green' if good else 'red')
    fig.suptitle(f'PSF candidate stamps {title}', fontsize=12)
    plt.tight_layout()
    plt.show()
    plt.close(fig)


def plotOuterStepPsf(expNext, sources, starCat, forcedPositions, psfPrev, psfNew,
                        psfUpdated, outerIdx):
    """Plot the PSF refit of one outer step: star selection and old vs new PSF.

    Parameters
    ----------
    expNext : `lsst.afw.image.Exposure`
        Background-subtracted exposure.
    sources : `lsst.afw.table.SourceCatalog` or `None`
        Fresh detections.
    starCat : `lsst.afw.table.SourceCatalog` or `None`
        PSF candidates.
    forcedPositions : array-like or `None`
        Forced original star positions.
    psfPrev, psfNew : `lsst.afw.detection.Psf`
        Previous and new PSF.
    psfUpdated : `bool`
        Whether the PSF was updated.
    outerIdx : `int`
        Outer-pass index.
    """
    bbox = expNext.getBBox()
    ext = _extent(bbox)
    data = expNext.image.array
    c = bbox.getCenter()
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))

    axes[0, 0].imshow(data, extent=ext, origin='lower', cmap='gray',
                      norm=simple_norm(data, stretch='asinh', min_percent=1, max_percent=99.5))
    if sources is not None and len(sources):
        sp = np.array([(r.getX(), r.getY()) for r in sources])
        axes[0, 0].scatter(sp[:, 0], sp[:, 1], s=8, c='cyan', label=f'fresh detections ({len(sp)})')
    if starCat is not None and len(starCat):
        st = np.array([(r.getX(), r.getY()) for r in starCat])
        axes[0, 0].scatter(st[:, 0], st[:, 1], facecolors='none', edgecolors='orange', s=80,
                           label=f'PSF candidates ({len(st)})')
    if forcedPositions is not None and len(forcedPositions):
        fp = _peaksArray(forcedPositions)
        axes[0, 0].scatter(fp[:, 0], fp[:, 1], marker='s', facecolors='none', edgecolors='lime',
                           s=110, label=f'forced original stars ({len(fp)})')
    axes[0, 0].set_xlim(ext[0], ext[1]); axes[0, 0].set_ylim(ext[2], ext[3])
    axes[0, 0].legend(fontsize=8, loc='upper right')
    axes[0, 0].set_title('Star selection on refit-background image')
    _cleanAxes(axes[0, 0])

    kp = psfPrev.computeKernelImage(c).array
    kn = psfNew.computeKernelImage(c).array
    for a, k, lab, psfO in ((axes[0, 1], kp, 'previous', psfPrev), (axes[0, 2], kn, 'new', psfNew)):
        im = a.imshow(k, origin='lower', cmap='inferno')
        plt.colorbar(im, ax=a, fraction=0.046, pad=0.04)
        a.set_title(f'{lab} PSF  FWHM={_fwhm(psfO, c):.3f}px  ({k.shape[0]}x{k.shape[1]})')
        _cleanAxes(a)
    if not psfUpdated:
        axes[0, 2].set_title(axes[0, 2].get_title() + '\n(NOT updated -- same as previous)', color='red')

    if kp.shape == kn.shape:
        d = kn - kp
        s = np.nanmax(np.abs(d)) or 1.0
        im = axes[1, 0].imshow(d, origin='lower', cmap='RdBu_r', vmin=-s, vmax=s)
        plt.colorbar(im, ax=axes[1, 0], fraction=0.046, pad=0.04)
        axes[1, 0].set_title('new - previous')
        _cleanAxes(axes[1, 0])
    else:
        axes[1, 0].axis('off')
        axes[1, 0].text(0.5, 0.5, 'kernel sizes differ', ha='center', transform=axes[1, 0].transAxes)

    for a, sel, lab in ((axes[1, 1], 'x', 'x'), (axes[1, 2], 'y', 'y')):
        for k, name, colr in ((kp, 'previous', 'tab:blue'), (kn, 'new', 'tab:orange')):
            cy, cx = k.shape[0] // 2, k.shape[1] // 2
            prof = k[cy, :] if sel == 'x' else k[:, cx]
            a.plot(np.arange(len(prof)) - len(prof) // 2, prof, color=colr, label=name,
                   ls='-' if name == 'previous' else '--')
        a.set_yscale('log')
        a.set_xlabel(f'{lab} offset [px]')
        a.set_title(f'PSF {lab} cross-section (log)')
        a.legend()
    fig.suptitle(f'Outer {outerIdx + 1} → {outerIdx + 2}: PSF redetermination', fontsize=13)
    plt.tight_layout()
    plt.show()
    plt.close(fig)


# =====================================================================
# PSF refinement from the bright stars' residuals (+ noise diagnostic)
# =====================================================================
# After a run:
#     highfreqNoiseCheck(outerHistory)                  # is variance too small?
#     psfNew, info = refinePsfFromFit(outerHistory)   # one update step
#     plotPsfCorrection(info)
# or let the pipeline do it each outer pass with refitPsf='stack'.

def highfreqNoiseCheckArrays(resid, variance, model, good, nbins=6, verbose=True):
    """Compare the residual's pixel-to-pixel scatter with the variance plane.

    Binned by model flux: ratio ~ 1 everywhere means the variance is right; a ratio growing with
    model flux means the variance lacks source photon noise.

    Parameters
    ----------
    resid : `numpy.ndarray`
        Residual image.
    variance : `numpy.ndarray`
        Variance image.
    model : `numpy.ndarray`
        Model image.
    good : `numpy.ndarray` [`bool`]
        Usable pixels.
    nbins : `int`, optional
        Number of model-flux bins.
    verbose : `bool`, optional
        Print the table.

    Returns
    -------
    result : `dict`
        ``flux`` (bin centres), ``ratio`` (scatter^2 / variance) and ``k`` (extra variance per
        unit model flux: var_true ~ var + k * model).
    """
    d = (resid[:, 1:] - resid[:, :-1]) / np.sqrt(2)
    m = 0.5 * (model[:, 1:] + model[:, :-1])
    v = 0.5 * (variance[:, 1:] + variance[:, :-1])
    g = good[:, 1:] & good[:, :-1]
    d, m, v = d[g], m[g], v[g]
    pos = m > 0
    edges = np.unique(np.quantile(m[pos], np.linspace(0, 1, nbins + 1))) if pos.any() else []
    ratio, mid, dv = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (m >= lo) & (m <= hi)
        if s.sum() < 20:
            continue
        rob = (1.4826 * np.median(np.abs(d[s] - np.median(d[s])))) ** 2
        ratio.append(rob / np.mean(v[s]))
        mid.append(np.mean(m[s]))
        dv.append(rob - np.mean(v[s]))
    mid, dv, ratio = np.array(mid), np.array(dv), np.array(ratio)
    k = float(np.sum(dv[-2:] * mid[-2:]) / np.sum(mid[-2:] ** 2)) if len(mid) >= 2 else np.nan
    if verbose:
        print("model flux (bin mean)   scatter^2 / variance-plane")
        for a, b in zip(mid, ratio):
            print(f"  {a:12.1f}          {b:6.2f}")
        print(f"extra variance per unit model flux (bright bins): k = {k:.4g}"
              f"  -> use var + k*model in effectiveVariance if ratio grows with flux")
    return dict(flux=mid, ratio=ratio, k=k)


def highfreqNoiseCheck(outerHistory, outerIdx=-1, snapIdx=-1, **kw):
    """Run `highfreqNoiseCheckArrays` on one snapshot of a run.

    Parameters
    ----------
    outerHistory : `list` [`dict`]
        Outer-loop history.
    outerIdx : `int`, optional
        Outer pass (default: last).
    snapIdx : `int`, optional
        Inner-loop snapshot (default: last).
    **kw
        Keyword arguments of `highfreqNoiseCheckArrays`.

    Returns
    -------
    result : `dict`
        See `highfreqNoiseCheckArrays`.
    """
    h = outerHistory[outerIdx]
    snap = h['result']['snaps'][snapIdx]
    fit = snap['fit']
    resid = h['exposure'].image.array - fit['modelFull']
    return highfreqNoiseCheckArrays(resid, fit['variance'], fit['modelFull'], fit['good'], **kw)


def estimatePsfCorrection(data, model, variance, good, xy0, positions, amps, kernel,
                            nUse=8, step=0.8, clip=3.0, softSigma=1.0,
                            window=(0.6, 1.0), maxRadius=None):
    """Estimate an improved PSF kernel from neighbour-subtracted residuals.

    For each bright star, (data - model) / amp sampled on the kernel grid around its fitted
    position estimates (true PSF - current PSF).  The weighted, clipped mean over stars is
    soft-thresholded, windowed to zero at the kernel edge, stripped of the components degenerate
    with position (dP/dx, dP/dy), added to the kernel and renormalised.

    Parameters
    ----------
    data : `numpy.ndarray`
        Image.
    model : `numpy.ndarray`
        Joint-fit model.
    variance : `numpy.ndarray`
        Variance used for the weights.
    good : `numpy.ndarray` [`bool`]
        Usable pixels.
    xy0 : `tuple` [`int`, `int`]
        Parent pixel coordinates of ``data[0, 0]``.
    positions : array-like
        Star positions, shape (N, 2).
    amps : array-like
        Star amplitudes.
    kernel : `numpy.ndarray`
        Current PSF kernel.
    nUse : `int`, optional
        Number of brightest stars used.
    step : `float`, optional
        Fraction of the correction applied.
    clip : `float`, optional
        Sigma clipping of the stack.
    softSigma : `float`, optional
        Soft threshold of the correction [sigma].
    window : `tuple` [`float`, `float`], optional
        Inner and outer edge of the window as fractions of the kernel half-size.
    maxRadius : `float`, optional
        Do not change the kernel beyond this radius [px].

    Returns
    -------
    new : `numpy.ndarray`
        Corrected, normalised kernel.
    info : `dict`
        ``old``, ``new``, ``delta``, ``rawDelta``, ``err``, ``nStars``, ``used``.
    """
    positions = np.asarray(positions, float)
    amps = np.asarray(amps, float)
    ny, nx = data.shape
    resid = data - model
    goodF = good.astype(float)
    h = kernel.shape[0] // 2
    u = np.arange(-h, h + 1)
    V, U = np.meshgrid(u, u, indexing='ij')

    order = np.argsort(-amps)[:nUse]
    R, W = [], []
    for j in order:
        if amps[j] <= 0:
            continue
        xj, yj = positions[j, 0] - xy0[0], positions[j, 1] - xy0[1]
        yy, xx = yj + V, xj + U
        coords = [yy.ravel(), xx.ravel()]
        inside = (yy >= 0) & (yy <= ny - 1) & (xx >= 0) & (xx <= nx - 1)
        r = map_coordinates(resid, coords, order=3, mode='nearest').reshape(V.shape) / amps[j]
        gm = map_coordinates(goodF, coords, order=1, mode='nearest').reshape(V.shape) > 0.99
        v = map_coordinates(variance, coords, order=1, mode='nearest').reshape(V.shape)
        R.append(r)
        W.append(np.where(inside & gm & (v > 0), amps[j] ** 2 / np.maximum(v, 1e-30), 0.0))
    if not R:
        raise ValueError("no usable sources for PSF correction")
    R, W = np.array(R), np.array(W)

    if len(R) >= 4:
        med = np.median(R, axis=0)
        mad = 1.4826 * np.median(np.abs(R - med), axis=0)
        floor = 1.0 / np.sqrt(np.maximum(W.mean(axis=0), 1e-30))
        W = np.where(np.abs(R - med) > clip * np.maximum(mad, floor), 0.0, W)

    Wsum = W.sum(axis=0)
    with np.errstate(divide='ignore', invalid='ignore'):
        delta = np.where(Wsum > 0, (W * R).sum(axis=0) / Wsum, 0.0)
        err = np.where(Wsum > 0, 1.0 / np.sqrt(Wsum), np.inf)
    rawDelta = delta.copy()
    delta = np.sign(delta) * np.maximum(np.abs(delta) - softSigma * err, 0.0)

    # maxRadius: the correction is tapered to zero at this radius (px) instead of at the
    # kernel edge -- use it with an extended (Moffat-tailed) kernel so the noisy outer
    # residuals never modify the analytic tail, e.g. stackKwargs=dict(maxRadius=20).
    he = h if maxRadius is None else min(h, maxRadius)
    win = window if window else (0.0, 1.0)
    rr = np.hypot(U, V)
    w = np.clip((he - rr) / max(he * (win[1] - win[0]), 1e-9), 0, 1)
    if window or maxRadius is not None:
        delta = delta * (0.5 - 0.5 * np.cos(np.pi * w))

    gy, gx = np.gradient(kernel)
    A = np.column_stack([gx.ravel(), gy.ravel()])
    coef, *_ = np.linalg.lstsq(A, delta.ravel(), rcond=None)
    delta = delta - (A @ coef).reshape(delta.shape)

    new = kernel + step * delta
    new = new / new.sum()
    info = dict(old=kernel, new=new, delta=delta, rawDelta=rawDelta, err=err,
                nStars=len(R), used=[int(j) for j in order])
    return new, info


def kernelToLsstPsf(kernel):
    """Wrap a kernel array as a spatially constant LSST PSF.

    Parameters
    ----------
    kernel : `numpy.ndarray`
        Square kernel of odd size, centred at n//2, summing to 1.

    Returns
    -------
    psf : `lsst.meas.algorithms.KernelPsf`
        The PSF.
    """
    import lsst.afw.image as afwImage
    import lsst.afw.math as afwMath
    try:
        from lsst.meas.algorithms import KernelPsf
    except ImportError:
        from lsst.afw.detection import KernelPsf
    img = afwImage.makeImageFromArray(np.ascontiguousarray(kernel, dtype=np.float64))
    return KernelPsf(afwMath.FixedKernel(img))


def extendKernelWithMoffat(kernel, newSize=81, fitR=(6.0, 14.0), blend=(12.0, 16.0)):
    """Replace the noisy outer part of a PSF kernel by an analytic Moffat tail.

    The azimuthal-mean profile between ``fitR[0]`` and ``fitR[1]`` px is fitted with
    ``A * (1 + (r/a)^2)^-b`` (relative weighting, so faint outer bins count).  The new kernel
    equals the original inside ``blend[0]``, the Moffat beyond ``blend[1]``, with a linear
    cross-fade between, and is renormalised -- fitted amplitudes then become total fluxes
    including the tail.  Works on an already-extended kernel too.

    Parameters
    ----------
    kernel : `numpy.ndarray`
        PSF kernel.
    newSize : `int`, optional
        Size of the new kernel [px] (odd).
    fitR : `tuple` [`float`, `float`], optional
        Radial range of the Moffat fit [px].
    blend : `tuple` [`float`, `float`], optional
        Radial range of the cross-fade [px].

    Returns
    -------
    new_kernel : `numpy.ndarray`
        Extended, normalised kernel.
    info : `dict`
        Moffat ``A``, ``a``, ``beta`` and ``tailFrac`` (flux fraction beyond ``blend[1]``).

    Raises
    ------
    RuntimeError
        Raised if the Moffat fit fails.
    """
    n = kernel.shape[0]
    h = n // 2
    u = np.arange(-h, h + 1)
    r = np.hypot(*np.meshgrid(u, u))
    rHi = min(fitR[1], h - 1.0)
    edges = np.arange(fitR[0], rHi + 1.0)
    rc = 0.5 * (edges[1:] + edges[:-1])
    prof = np.array([kernel[(r >= lo) & (r < hi)].mean() for lo, hi in zip(edges[:-1], edges[1:])])
    moffat = lambda rr, A, a, b: A * (1.0 + (rr / a) ** 2) ** (-b)
    sigma = np.clip(prof, 1e-6 * np.nanmax(np.abs(prof)), None)      # relative weighting
    (A, a, b), _ = curve_fit(moffat, rc, prof, p0=[max(prof[0], 1e-12), 4.0, 3.0], sigma=sigma,
                             bounds=([0.0, 0.5, 1.05], [np.inf, 60.0, 12.0]), maxfev=5000)

    H = newSize // 2
    v = np.arange(-H, H + 1)
    R = np.hypot(*np.meshgrid(v, v))
    tail = moffat(R, A, a, b)
    hh = min(h, H)
    core = np.zeros_like(tail)
    core[H - hh:H + hh + 1, H - hh:H + hh + 1] = kernel[h - hh:h + hh + 1, h - hh:h + hh + 1]
    b1 = min(blend[1], hh - 1.0)
    b0 = min(blend[0], b1 - 2.0)
    w = np.clip((b1 - R) / (b1 - b0), 0.0, 1.0)                      # 1 inside b0, 0 beyond b1
    out = w * core + (1.0 - w) * tail
    out = out / out.sum()
    info = dict(A=float(A), a=float(a), beta=float(b), tailFrac=float(out[R > b1].sum()))
    return out, info


def extendPsfWings(psf, point, size=81, verbose=True, **kw):
    """Return a constant KernelPsf with a Moffat tail, taken from any PSF at a point.

    Parameters
    ----------
    psf : `lsst.afw.detection.Psf`
        Input PSF.
    point : `lsst.geom.Point2D`
        Position at which the kernel is taken.
    size : `int`, optional
        Size of the extended kernel [px].
    verbose : `bool`, optional
        Print the fitted tail.
    **kw
        Keyword arguments of `extendKernelWithMoffat`.

    Returns
    -------
    psf : `lsst.afw.detection.Psf`
        The extended PSF, or the input PSF if the fit fails.
    """
    kernel = psf.computeKernelImage(point).array.astype(float).copy()
    try:
        new, info = extendKernelWithMoffat(kernel, newSize=size, **kw)
    except Exception as e:
        print(f"  PSF wing extension failed ({type(e).__name__}: {e}) -- keeping current PSF")
        return psf
    if verbose:
        print(f"  PSF wings extended: {kernel.shape[0]}px -> {size}px kernel; Moffat "
              f"a={info['a']:.2f}px, beta={info['beta']:.2f}; "
              f"{100 * info['tailFrac']:.2f}% of flux in the tail beyond the blend radius")
    return kernelToLsstPsf(new)


def refinePsfFromFit(outerHistory, outerIdx=-1, snapIdx=-1, **kw):
    """Make one PSF update step from a finished run.

    Parameters
    ----------
    outerHistory : `list` [`dict`]
        Outer-loop history.
    outerIdx : `int`, optional
        Outer pass (default: last).
    snapIdx : `int`, optional
        Inner-loop snapshot (default: last).
    **kw
        Keyword arguments of `estimatePsfCorrection`.

    Returns
    -------
    psf : `lsst.meas.algorithms.KernelPsf`
        Updated PSF.
    info : `dict`
        See `estimatePsfCorrection`.
    """
    h = outerHistory[outerIdx]
    exp, psf = h['exposure'], h['psf']
    snap = h['result']['snaps'][snapIdx]
    fit = snap['fit']
    kernel = psf.computeKernelImage(exp.getBBox().getCenter()).array.copy()
    xy0 = (exp.getXY0().getX(), exp.getXY0().getY())
    new, info = estimatePsfCorrection(exp.image.array, fit['modelFull'], fit['variance'],
                                        fit['good'], xy0, snap['mergedArr'], fit['amps'],
                                        kernel, **kw)
    return kernelToLsstPsf(new), info


def polishPsfFixedSources(exposure, psf, result, badBitmask, nCycles=3,
                             stackKwargs=None, supportFrac=0.01, sysFrac=0.03,
                             poissonK=0.0, refineCentroids=True, centroidMinAmp=50.0,
                             freezeSep=None, maxShift=0.3, verbose=True):
    """Iterate PSF correction and re-solve without re-detecting sources.

    A cycle is kept only if the residual rms over good pixels drops.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    psf : `lsst.afw.detection.Psf`
        Current PSF.
    result : `dict`
        Inner-loop result.
    badBitmask : `int`
        Mask bits of pixels to ignore.
    nCycles : `int`, optional
        Maximum number of cycles.
    stackKwargs : `dict`, optional
        Keyword arguments of `estimatePsfCorrection`.
    supportFrac, sysFrac, poissonK : `float`, optional
        See `jointPsfFitDiag`.
    refineCentroids : `bool`, optional
        Refine positions after each PSF change.
    centroidMinAmp : `float`, optional
        Fainter sources keep their position.
    freezeSep : `float`, optional
        See `refinePeakPositions`.
    maxShift : `float`, optional
        Maximum position shift [px].
    verbose : `bool`, optional
        Print each cycle.

    Returns
    -------
    psf : `lsst.afw.detection.Psf`
        The best PSF.
    info : `dict` or `None`
        Information of the last kept correction (see `estimatePsfCorrection`).
    """
    snap = result['snaps'][-1]
    peaks = [tuple(p) for p in snap['mergedArr']]
    fit = snap['fit']
    data = exposure.image.array

    def rms(f):
        return float(np.sqrt(np.mean((data - f['modelFull'])[f['good']] ** 2)))

    cur, best, info = psf, rms(fit), None
    for c in range(nCycles):
        hist = [dict(exposure=exposure, psf=cur,
                     result=dict(snaps=[dict(fit=fit, mergedArr=_peaksArray(peaks))]))]
        try:
            cand, cinfo = refinePsfFromFit(hist, -1, -1, **(stackKwargs or {}))
        except ValueError as e:
            print(f"    polish cycle {c + 1}: {e} -- stopping")
            break
        info = cinfo
        fitC = jointPsfFitDiag(exposure, cand, peaks, badBitmask,
                                   supportFrac=supportFrac, varModel=fit['modelFull'],
                                   sysFrac=sysFrac, poissonK=poissonK)
        peaksC = peaks
        if refineCentroids:
            peaksC = refinePeakPositions(
                exposure, cand, peaks, fitC['amps'], badBitmask,
                minAmpForRefine=centroidMinAmp, maxShift=maxShift,
                variance=fitC['variance'], freezeSep=freezeSep)
            fitC = jointPsfFitDiag(exposure, cand, peaksC, badBitmask,
                                       supportFrac=supportFrac,
                                       varModel=fitC['modelFull'],
                                       sysFrac=sysFrac, poissonK=poissonK)
        r = rms(fitC)
        ok = r < best
        if verbose:
            print(f"    polish cycle {c + 1}/{nCycles}: residual rms {best:.3f} -> {r:.3f} "
                  f"({'kept' if ok else 'REJECTED, stopping'})")
        if not ok:
            break
        cur, fit, peaks, best = cand, fitC, peaksC, r
    return cur, info


def plotPsfCorrection(info):
    """Plot the old kernel, the proposed correction and the radial profiles.

    Parameters
    ----------
    info : `dict`
        Output of `estimatePsfCorrection`.
    """
    old, new, d, raw, err = info['old'], info['new'], info['delta'], info['rawDelta'], info['err']
    h = old.shape[0] // 2
    ext = (-h - .5, h + .5, -h - .5, h + .5)
    fig, ax = plt.subplots(1, 4, figsize=(20, 4.6))
    im = ax[0].imshow(old, origin='lower', extent=ext, cmap='inferno'); plt.colorbar(im, ax=ax[0], fraction=.046)
    ax[0].set_title('current kernel')
    s = np.nanmax(np.abs(raw[np.isfinite(raw)])) or 1
    im = ax[1].imshow(raw, origin='lower', extent=ext, cmap='RdBu_r', vmin=-s, vmax=s); plt.colorbar(im, ax=ax[1], fraction=.046)
    ax[1].set_title(f"stacked (data-model)/amp, {info['nStars']} stars\n(true PSF - current PSF)")
    s = np.nanmax(np.abs(d)) or 1
    im = ax[2].imshow(d, origin='lower', extent=ext, cmap='RdBu_r', vmin=-s, vmax=s); plt.colorbar(im, ax=ax[2], fraction=.046)
    ax[2].set_title('correction applied (thresholded, windowed)')
    r = np.hypot(*np.meshgrid(np.arange(-h, h + 1), np.arange(-h, h + 1)))
    edges = np.arange(0, h + 1)
    rc = .5 * (edges[1:] + edges[:-1])
    prof = lambda a: [a[(r >= lo) & (r < hi)].mean() for lo, hi in zip(edges[:-1], edges[1:])]
    ax[3].plot(rc, prof(old), 'k', label='old kernel')
    ax[3].plot(rc, prof(new), '--', color='#DD8452', label='new kernel')
    ax[3].plot(rc, prof(raw), color='gray', lw=1, label='raw stacked diff')
    ax[3].set_yscale('symlog', linthresh=1e-5); ax[3].set_xlabel('r [px]')
    ax[3].set_title('azimuthal profiles'); ax[3].legend(fontsize=8)
    plt.tight_layout(); plt.show()


# =====================================================================
# Full pipeline on the cutout
# =====================================================================

def _wantsMoreOuter(result, outerIdx):
    """Return False: never extend the outer loop (replaced by ``adaptive.install()``)."""
    return False


def runCutoutPipeline(exposureCrop, psf, badBitmask, starCat=None, bgImage=None,
                        xRange=(0, 100), yRange=(0, 100),
                        nOuterIterations=2, nInnerIterations=10, innerKwargs=None,
                        outerConvergenceFrac=0.02, refitPsf=True, psfKernelSize=41,
                        psfKernelSizes=None,
                        bgBinSize=32, bgUseApprox=True,
                        viz=True, plotEvery=1, sourceView='all', sourcesPerPage=6,
                        maxPagesPerIter=None, viewHalf=20, supportFrac=0.01,
                        refineCentroids=True, centroidMinAmp=50.0,
                        nCentroidRounds=4, centroidShiftTol=0.03,
                        centroidShiftSchedule=None, coarseCentroid=True,
                        coarseCentroidHwSigma=2.5, mergeSeparationSigma=1.5,
                        sysFrac=0.03, minSnr=5.0, pruneSepSigma=1.5,
                        rejectFwhm=0.75, brightPeakSigma=100.0, freezeFwhm=1.0, poissonK=0.0,
                        nearBrightFwhm=2.0, nearBrightSnr=10.0,
                        edgeMarginFwhm=0.5, useEliminateWeak=True,
                        useShapeCut=True, shapeKwargs=None, shapeInFirstPass=False,
                        stackKwargs=None, psfPolishCycles=3, psfWingSize=81):
    """Run the outer loop (inner fit, background and PSF refinement) on a cutout.

    Parameters
    ----------
    exposureCrop : `lsst.afw.image.Exposure`
        Region to process (before background subtraction).
    psf : `lsst.afw.detection.Psf`
        Initial PSF.
    badBitmask : `int`
        Mask bits of pixels to ignore.
    starCat : `lsst.afw.table.SourceCatalog`, optional
        Original PSF stars, forced into every PSF determination.
    bgImage : `lsst.afw.image.Image`, optional
        Detector background subtracted before the first pass.
    xRange, yRange : `tuple` [`int`, `int`], optional
        Sub-region of ``exposureCrop`` (see `makeCutout`).
    nOuterIterations : `int`, optional
        Number of outer passes.
    nInnerIterations : `int`, optional
        Maximum inner iterations per pass.
    innerKwargs : `dict`, optional
        Keyword arguments of `runInnerLoopVerbose`.
    outerConvergenceFrac : `float`, optional
        Stop when the residual standard deviation improves by less than this fraction.
    refitPsf : `bool` or `str`, optional
        False, True (PCA star-selector refit) or 'stack' (residual-stack refit from the bright
        stars, see `refinePsfFromFit`).
    psfKernelSize : `int`, optional
        PSF kernel size [px].
    psfKernelSizes : `tuple` [`int`], optional
        Kernel sizes tried by `redetectAndRefitPsfAdaptive`.
    bgBinSize : `int`, optional
        Background bin size [px].
    bgUseApprox : `bool`, optional
        See `refitBackground`.
    viz : `bool`, optional
        Draw diagnostics.
    plotEvery : `int`, optional
        Draw every n-th inner iteration (the last always).
    sourceView : `str`, optional
        'all', 'new' or 'changed': which sources get a per-source page.
    sourcesPerPage : `int`, optional
        Sources per page.
    maxPagesPerIter : `int`, optional
        Maximum per-source pages per iteration (None = no limit).
    viewHalf : `int`, optional
        Half-size of the per-source windows [px].
    supportFrac : `float`, optional
        See `jointPsfFitDiag`.
    refineCentroids : `bool`, optional
        Refine positions (`refinePeakPositions`).
    centroidMinAmp : `float`, optional
        Fainter sources keep their position.
    nCentroidRounds, centroidShiftTol, centroidShiftSchedule : optional
        Centroid refinement controls (see `runInnerLoopVerbose`).
    coarseCentroid, coarseCentroidHwSigma : optional
        Windowed starting centroid of new peaks (see `runInnerLoopVerbose`).
    mergeSeparationSigma : `float`, optional
        Merge new peaks closer than this many PSF sigma.
    sysFrac, minSnr, pruneSepSigma, rejectFwhm, brightPeakSigma : `float`, optional
        Over-counting guards (see `runInnerLoopVerbose`).
    freezeFwhm, poissonK, nearBrightFwhm, nearBrightSnr : `float`, optional
        Over-counting guards (see `runInnerLoopVerbose`).
    edgeMarginFwhm : `float`, optional
        Drop sources within this many FWHM of the cutout edge.
    useEliminateWeak : `bool`, optional
        Backward elimination on the honest S/N.
    useShapeCut : `bool`, optional
        PSF-shape constraint (`eliminateNonpsf`).
    shapeKwargs : `dict`, optional
        Keyword arguments of `shapeMetrics` / `judgeShape`.
    shapeInFirstPass : `bool`, optional
        Apply the shape cut in the first outer pass too (by default only once the PSF has been
        refined: a rough first PSF can make the cut reject real faint companions).
    stackKwargs : `dict`, optional
        Keyword arguments of `refinePsfFromFit` (``refitPsf='stack'``).
    psfPolishCycles : `int`, optional
        PSF correction / re-solve cycles per outer pass (``refitPsf='stack'``).
    psfWingSize : `int`, optional
        Convert every PSF to a constant kernel of this size with a Moffat tail (see
        `extendPsfWings`); None = off.  Stops the model ending at the stamp edge, which leaves
        a ring of fake sources half a kernel width from bright stars.

    Returns
    -------
    outerHistory : `list` [`dict`]
        One dict per outer pass: ``psf``, ``exposure``, ``result``, ``residStd``, ``nPeaks``,
        ``bgImage``.
    """
    TIMER.reset()
    if psfKernelSizes is None:
        fallbacks = [k for k in (27, 21, 15, 11) if k < psfKernelSize]
        psfKernelSizes = tuple([psfKernelSize] + fallbacks)
    innerKwargs = dict(nIterations=nInnerIterations, thresholdStart=8.0,
                        thresholdStep=-1.5, thresholdMin=4.5, minPixels=1,
                        minNewPeaks=1, convergenceFrac=0.02, **(innerKwargs or {})) \
        if innerKwargs is None or 'nIterations' not in innerKwargs else dict(innerKwargs)

    with TIMER.step('make cutout'):
        cutout = makeCutout(exposureCrop, xRange, yRange)
        cutout.setPsf(psf)
    cb = cutout.getBBox()
    if psfWingSize:
        psf = extendPsfWings(psf, cb.getCenter(), psfWingSize)
        cutout.setPsf(psf)
    print(f"Cutout parent-frame bbox: x[{cb.getMinX()}, {cb.getMaxX()}] "
          f"y[{cb.getMinY()}, {cb.getMaxY()}]")

    originalStarPositions = None
    if starCat is not None and len(starCat) > 0:
        allpos = [(r.getCentroid().getX(), r.getCentroid().getY()) for r in starCat]
        originalStarPositions = [(x, y) for x, y in allpos
                                   if cb.contains(geom.Point2I(int(round(x)), int(round(y))))]
        print(f"{len(originalStarPositions)} of {len(allpos)} original PSF stars fall inside the cutout")

    if viz:
        plotCutoutSetup(cutout, psf, badBitmask, originalStarPositions)

    bgPrev = None
    if bgImage is not None:
        try:
            bgPrev = bgImage[cb].array.copy()
        except Exception as e:
            print(f"(could not crop initial bgImage: {e})")

    raw = cutout.clone()
    exposureWorking = cutout.clone()
    exposureWorking.setPsf(psf)
    psfWorking = psf
    bgImageWorking = None
    outerHistory = []

    for outerIt in __import__("itertools").count():   # open-ended; see break below
        print(f"\n{'=' * 70}\n===== Outer pass {outerIt + 1}/{nOuterIterations} "
              f"(PSF FWHM={_fwhm(psfWorking, cb.getCenter()):.3f}px) =====\n{'=' * 70}")
        shapeNow = useShapeCut and (outerIt > 0 or shapeInFirstPass)
        print(f"  PSF-shape cut this pass: {'ON' if shapeNow else 'off'}")
        with TIMER.step('INNER LOOP (total, incl. plotting)', outer=outerIt, kind='wall'):
            result = runInnerLoopVerbose(
                exposureWorking, psfWorking, badBitmask, outerIdx=outerIt,
                viz=viz, plotEvery=plotEvery, sourceView=sourceView,
                sourcesPerPage=sourcesPerPage, maxPagesPerIter=maxPagesPerIter,
                viewHalf=viewHalf, supportFrac=supportFrac,
                refineCentroids=refineCentroids, centroidMinAmp=centroidMinAmp,
                nCentroidRounds=nCentroidRounds, centroidShiftTol=centroidShiftTol,
                centroidShiftSchedule=centroidShiftSchedule,
                coarseCentroid=coarseCentroid,
                coarseCentroidHwSigma=coarseCentroidHwSigma,
                mergeSeparationSigma=mergeSeparationSigma,
                sysFrac=sysFrac, minSnr=minSnr, pruneSepSigma=pruneSepSigma,
                rejectFwhm=rejectFwhm, brightPeakSigma=brightPeakSigma,
                freezeFwhm=freezeFwhm, poissonK=poissonK,
                nearBrightFwhm=nearBrightFwhm, nearBrightSnr=nearBrightSnr,
                edgeMarginFwhm=edgeMarginFwhm, useEliminateWeak=useEliminateWeak,
                useShapeCut=shapeNow, shapeKwargs=shapeKwargs,
                **innerKwargs)
        finalStd = result['residStdHistory'][-1]
        outerHistory.append(dict(psf=psfWorking, exposure=exposureWorking, result=result,
                                  residStd=finalStd, nPeaks=len(result['peaks']),
                                  bgImage=bgImageWorking))
        print(f"\nouter pass {outerIt + 1}: {len(result['peaks'])} peaks, residual std={finalStd:.4f}")

        if viz:
            with TIMER.step('plotting', outer=outerIt, kind='viz'):
                plotAmplitudeEvolution(result, outerIt)
                plotSourceEvolution(exposureWorking, result, outerIt, viewHalf=viewHalf)

        if outerIt > 0:
            prev = outerHistory[-2]['residStd']
            impr = abs(prev - finalStd) / prev if prev > 0 else 0
            if impr < outerConvergenceFrac:
                print(f"outer residual std improved by only {impr * 100:.2f}% -- outer loop converged")
                break
        if outerIt >= nOuterIterations - 1 and not _wantsMoreOuter(result, outerIt):
            break

        # ---- re-derive background + PSF ------------------------------------
        exposureNext = raw.clone()
        with TIMER.step('mask sources from model', outer=outerIt):
            srcMask = maskSourcesFromModel(exposureNext, result['modelFull'])
        bgNext = None
        with TIMER.step('refit background', outer=outerIt):
            try:
                exposureNext, bgNext = refitBackground(
                    exposureNext, useApprox=bgUseApprox, binSize=bgBinSize)
            except Exception as e:
                print(f"  background refit failed ({type(e).__name__}: {e}) -- "
                      f"keeping un-refit image")
        if viz:
            with TIMER.step('plotting', outer=outerIt, kind='viz'):
                plotOuterStepBackground(raw, exposureNext, srcMask, result['modelFull'],
                                           bgNext, bgPrev, outerIt)
        bgPrev = bgNext.array.copy() if bgNext is not None else bgPrev

        psfNew, sourcesNew, starCat_new, cellSet_new = psfWorking, None, None, None
        exposureNext.setPsf(psfWorking)
        if refitPsf == 'stack':
            with TIMER.step('PSF refit (residual stack)', outer=outerIt):
                try:
                    if psfPolishCycles > 1:
                        psfNew, pinfo = polishPsfFixedSources(
                            exposureWorking, psfWorking, result, badBitmask,
                            nCycles=psfPolishCycles, stackKwargs=stackKwargs,
                            supportFrac=supportFrac, sysFrac=sysFrac,
                            poissonK=poissonK, refineCentroids=refineCentroids,
                            centroidMinAmp=centroidMinAmp,
                            freezeSep=freezeFwhm * _fwhm(psfWorking, cb.getCenter()))
                    else:
                        psfNew, pinfo = refinePsfFromFit(outerHistory, -1,
                                                             **(stackKwargs or {}))
                    print(f"  PSF updated from residual stack of {pinfo['nStars']} stars: "
                          f"FWHM {_fwhm(psfWorking, cb.getCenter()):.3f} -> "
                          f"{_fwhm(psfNew, cb.getCenter()):.3f} px")
                    if viz:
                        plotPsfCorrection(pinfo)
                except Exception as e:
                    print(f"  residual-stack PSF refit raised {type(e).__name__}: {e} -- "
                          f"keeping current PSF")
                    psfNew = psfWorking
            exposureNext.setPsf(psfNew)
        elif refitPsf:
            with TIMER.step('redetect + PSF refit', outer=outerIt):
                try:
                    psfNew, sourcesNew, starCat_new, cellSet_new, ksUsed = \
                        redetectAndRefitPsfAdaptive(
                            exposureNext, kernelSizes=psfKernelSizes,
                            extraPeakPositions=result['peaks'],
                            forceIncludeStarPositions=originalStarPositions)
                except Exception as e:
                    print(f"  PSF redetermination raised {type(e).__name__}: {e} -- "
                          f"keeping current PSF")
                    cellSet_new = None
            psfUpdated = cellSet_new is not None
            if not psfUpdated:
                psfNew = psfWorking
            exposureNext.setPsf(psfNew)
            if viz:
                with TIMER.step('plotting', outer=outerIt, kind='viz'):
                    plotOuterStepPsf(exposureNext, sourcesNew, starCat_new,
                                        originalStarPositions, psfWorking, psfNew,
                                        psfUpdated, outerIt)
                    plotPsfCandidateStamps(cellSet_new, title=f'(outer pass {outerIt + 2})')

        if psfWingSize and psfNew is not psfWorking:
            # PSF changed this pass: refresh the Moffat tail so it matches the new core
            psfNew = extendPsfWings(psfNew, cb.getCenter(), psfWingSize)
            exposureNext.setPsf(psfNew)

        exposureWorking, psfWorking, bgImageWorking = exposureNext, psfNew, bgNext

    plotTimingSummary(TIMER)
    return outerHistory


def plotWingDiagnostic(outerHistory, nStars=3, rmax=40):
    """Plot the azimuthal profiles of the brightest sources: data vs their own model.

    Data beyond the dashed line well above the model means a real halo the kernel cannot
    represent; a difference centred on zero without a step means only core / PSF-shape mismatch
    remains.

    Parameters
    ----------
    outerHistory : `list` [`dict`]
        Outer-loop history.
    nStars : `int`, optional
        Number of brightest sources shown.
    rmax : `float`, optional
        Maximum radius [px].
    """
    final = outerHistory[-1]
    exp, res = final['exposure'], final['result']
    snap = res['snaps'][-1]
    fit, peaks = snap['fit'], snap['mergedArr']
    bbox = exp.getBBox()
    ny, nx = bbox.getHeight(), bbox.getWidth()
    data, model, good = exp.image.array, fit['modelFull'], fit['good']
    yy, xx = np.mgrid[bbox.getMinY():bbox.getMaxY() + 1, bbox.getMinX():bbox.getMaxX() + 1]
    edges = np.arange(0, rmax + 1)
    rc = 0.5 * (edges[1:] + edges[:-1])

    picks = np.argsort(-fit['amps'])[:nStars]
    fig, axes = plt.subplots(1, len(picks), figsize=(6.2 * len(picks), 5), squeeze=False)
    for ax, j in zip(axes[0], picks):
        c = fit['comps'][j]
        if c is None:
            continue
        own = np.zeros((ny, nx))
        own[c['ys'], c['xs']] = c['amp'] * c['stamp']
        othersRemoved = data - (model - own)
        r = np.hypot(xx - peaks[j][0], yy - peaks[j][1])
        dProf, mProf = [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            sel = good & (r >= lo) & (r < hi)
            ok = sel.sum() > 3
            dProf.append(np.mean(othersRemoved[sel]) if ok else np.nan)
            mProf.append(np.mean(own[sel]) if ok else np.nan)
        dProf, mProf = np.array(dProf), np.array(mProf)
        ax.plot(rc, dProf, 'k', lw=1.5, label='data - neighbours')
        ax.plot(rc, mProf, '--', color='#DD8452', lw=1.5, label='own model')
        ax.plot(rc, dProf - mProf, color='gray', lw=1.2, label='difference')
        h = 0.5 * (c['stamp'].shape[0] - 1)
        ax.axvline(h, color='r', ls='--', lw=1, label=f'stamp half-width ({h:.0f}px)')
        if WING_TAPER is not None:
            ax.axvline(WING_TAPER[0] * h, color='r', ls=':', lw=1, label='taper starts')
        ax.axhline(0, color='gray', lw=0.5)
        ax.set_yscale('symlog', linthresh=10)
        ax.set_xlabel('radius [px]'); ax.set_ylabel('mean value in annulus')
        ax.set_title(f"source #{j} ({peaks[j][0]:.1f},{peaks[j][1]:.1f}) amp={fit['amps'][j]:.0f}")
        ax.legend(fontsize=7)
    plt.tight_layout()
    plt.show()
    plt.close(fig)


# =====================================================================
# Per-source cutout figures: data | model | residual | cross-section
# =====================================================================

def _window(bbox, px, py, half):
    """Return the (unclipped) array index of a source and the clipped window around it."""
    X0, Y0 = bbox.getMinX(), bbox.getMinY()
    ix, iy = int(round(px)) - X0, int(round(py)) - Y0
    xa, xb = max(ix - half, 0), min(ix + half + 1, bbox.getWidth())
    ya, yb = max(iy - half, 0), min(iy + half + 1, bbox.getHeight())
    return ix, iy, xa, xb, ya, yb


def _selectSources(peaks, amps, sourceIds, topN, minAmp):
    """Return the indices of the sources to plot (given ids, or the brightest above a flux)."""
    if sourceIds is not None:
        return [int(j) for j in sourceIds]
    order = np.argsort(-amps)
    order = [int(j) for j in order if amps[j] > minAmp]
    return order[:topN] if topN is not None else order


def plotSourceCutouts(outerHistory, outerIdx=-1, snapIdx=-1,
                        sourceIds=None, topN=None, minAmp=0.0,
                        half=15, crossAxis='x', showOwn=True,
                        sourcesPerFig=4, pdfPath=None, show=True,
                        cmap=None, markNeighbours=True):
    """Plot cutouts and cross-sections of individual sources of a run.

    Parameters
    ----------
    outerHistory : `list` [`dict`]
        Outer-loop history.
    outerIdx : `int`, optional
        Outer pass (default: last).
    snapIdx : `int`, optional
        Inner-loop snapshot (default: last).
    sourceIds : `list` [`int`], optional
        Sources to show; overrides ``topN`` / ``minAmp``.
    topN : `int`, optional
        Show the brightest N sources.
    minAmp : `float`, optional
        Skip sources with amplitude at or below this.
    half : `int`, optional
        Cutout half-width [px].
    crossAxis : `str`, optional
        'x', 'y' or 'both'.
    showOwn : `bool`, optional
        Draw the source's own PSF component in the cross-section.
    sourcesPerFig : `int`, optional
        Sources per figure.
    pdfPath : `str`, optional
        Also save every figure as a page of this PDF.
    show : `bool`, optional
        Display the figures.
    cmap : `str`, optional
        Colour map of the images.
    markNeighbours : `bool`, optional
        Mark the neighbouring sources.

    Returns
    -------
    figs : `list` [`matplotlib.figure.Figure`]
        The figures.
    """
    cmap = cmap or CMAP
    h = outerHistory[outerIdx]
    exposure = h['exposure']
    snap = h['result']['snaps'][snapIdx]
    fit = snap['fit']
    data = exposure.image.array
    model = fit['modelFull']
    resid = data - model
    peaks = np.asarray(snap['mergedArr'], float)
    amps = fit['amps']
    bad = snap['badMask']
    bbox = exposure.getBBox()
    X0, Y0 = bbox.getMinX(), bbox.getMinY()
    ny, nx = data.shape

    ids = _selectSources(peaks, amps, sourceIds, topN, minAmp)
    if not ids:
        print("No sources selected.")
        return []

    axesList = {'x': ['x'], 'y': ['y'], 'both': ['x', 'y']}[crossAxis]
    ncol = 3 + len(axesList)

    pdf = None
    if pdfPath is not None:
        from matplotlib.backends.backend_pdf import PdfPages
        pdf = PdfPages(pdfPath)

    figs = []
    pages = [ids[i:i + sourcesPerFig] for i in range(0, len(ids), sourcesPerFig)]
    for pi, page in enumerate(pages):
        fig, axes = plt.subplots(len(page), ncol, figsize=(4.6 * ncol, 3.6 * len(page)),
                                 squeeze=False)
        for r, j in enumerate(page):
            px, py = peaks[j]
            ix, iy, xa, xb, ya, yb = _window(bbox, px, py, half)

            if xb <= xa or yb <= ya:
                for a in axes[r]:
                    a.axis('off')
                axes[r, 0].text(0.5, 0.5, f"#{j} ({px:.1f}, {py:.1f}) outside image",
                                ha='center', transform=axes[r, 0].transAxes)
                continue

            iyC = int(np.clip(iy, 0, ny - 1))
            ixC = int(np.clip(ix, 0, nx - 1))

            ext = (X0 + xa - 0.5, X0 + xb - 0.5, Y0 + ya - 0.5, Y0 + yb - 0.5)
            d, m, rs = data[ya:yb, xa:xb], model[ya:yb, xa:xb], resid[ya:yb, xa:xb]

            lo = np.nanmin(np.concatenate([d.ravel(), m.ravel()]))
            hi = np.nanmax(np.concatenate([d.ravel(), m.ravel()]))
            norm = simple_norm(d, vmin=lo, vmax=hi, stretch='asinh')
            rlim = np.nanpercentile(np.abs(rs), 99)
            rlim = rlim if np.isfinite(rlim) and rlim > 0 else 1.0

            comp = fit['comps'][j]
            own = np.zeros((ny, nx))
            if comp is not None:
                own[comp['ys'], comp['xs']] = comp['amp'] * comp['stamp']

            panels = [(d, dict(norm=norm, cmap=cmap), 'data'),
                      (m, dict(norm=norm, cmap=cmap), 'model'),
                      (rs, dict(cmap='RdBu_r', vmin=-rlim, vmax=rlim), 'residual')]
            for c, (arr, kw, name) in enumerate(panels):
                ax = axes[r, c]
                im = ax.imshow(arr, extent=ext, origin='lower', **kw)
                plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                ax.add_patch(patches.Circle((px, py), 3.0, fill=False, ec='w', lw=1.0))
                ax.scatter([px], [py], marker='x', c='w', s=30, lw=1.2)
                if markNeighbours:
                    inw = ((peaks[:, 0] > ext[0]) & (peaks[:, 0] < ext[1]) &
                           (peaks[:, 1] > ext[2]) & (peaks[:, 1] < ext[3]))
                    inw[j] = False
                    if inw.any():
                        ax.scatter(peaks[inw, 0], peaks[inw, 1], marker='+', c='w', s=25, lw=0.8)
                if bad[ya:yb, xa:xb].any():
                    ax.imshow(np.ma.masked_where(~bad[ya:yb, xa:xb], np.ones_like(d)),
                              extent=ext, origin='lower', cmap='autumn', alpha=0.35,
                              interpolation='nearest')
                ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
                ax.grid(False)
                if name == 'residual':
                    ax.set_title(f"residual  (rms {np.nanstd(rs):.1f})", fontsize=9, loc='left')
                elif name == 'model':
                    ax.set_title(f"model  amp={amps[j]:.1f}", fontsize=9, loc='left')
                else:
                    ax.set_title(f"#{j}  ({px:.1f}, {py:.1f})", fontsize=9, loc='left')

            for k, axis in enumerate(axesList):
                ax = axes[r, 3 + k]
                if axis == 'x':
                    coord = X0 + np.arange(xa, xb)
                    sl, centre = (iyC, slice(xa, xb)), px
                else:
                    coord = Y0 + np.arange(ya, yb)
                    sl, centre = (slice(ya, yb), ixC), py
                ax.plot(coord, data[sl], 'k-', lw=1.6, label='data')
                ax.plot(coord, model[sl], '-', color='#DD8452', lw=1.6, alpha=0.9, label='model')
                if showOwn:
                    ax.plot(coord, own[sl], ':', color='#1f77b4', lw=1.3, label='this source')
                ax.plot(coord, resid[sl], '-', color='gray', lw=0.8, alpha=0.7, label='residual')
                ax.axhline(0, color='gray', lw=0.5)
                ax.axvline(centre, color='c', lw=1.0, ls='--')
                b = bad[sl]
                if b.any():
                    ax.plot(coord[b], np.zeros(int(b.sum())), '|', color='red', ms=8,
                            transform=ax.get_xaxis_transform())
                ax.set_xlim(coord[0], coord[-1])
                ax.set_xlabel(f'{axis} [px]')
                ax.set_title(f"{axis} cross-section", fontsize=9, loc='left')
                if r == 0 and k == 0:
                    ax.legend(fontsize=7, loc='upper right')

        fig.suptitle(f"Outer {outerIdx % len(outerHistory) + 1} "
                     f"final fit: per-source cutouts (page {pi + 1}/{len(pages)})",
                     fontsize=12, y=1.0)
        plt.tight_layout()
        if pdf is not None:
            pdf.savefig(fig, bbox_inches='tight')
        if show:
            plt.show()
        figs.append(fig)
        if not show:
            plt.close(fig)

    if pdf is not None:
        pdf.close()
        print(f"saved {len(pages)} page(s) to {pdfPath}")
    return figs


# =====================================================================
# False-positive diagnostics
# =====================================================================

def negativeSourceCheck(outerHistory, outerIdx=-1, snapIdx=-1, thresholds=(3, 4, 5),
                          plot=True):
    """Estimate the number of false detections from the residual itself.

    Counts pixels of the matched-filter S/N map of the final residual at >= +N sigma and at
    <= -N sigma.  Real sources are positive only, so the negative excursions approximate the
    positive noise excursions (false positives); PSF-mismatch dipoles inflate both.

    Parameters
    ----------
    outerHistory : `list` [`dict`]
        Outer-loop history.
    outerIdx : `int`, optional
        Outer pass (default: last).
    snapIdx : `int`, optional
        Inner-loop snapshot (default: last).
    thresholds : `tuple` [`float`], optional
        S/N thresholds.
    plot : `bool`, optional
        Plot the S/N map and histogram.

    Returns
    -------
    snrMap : `numpy.ndarray`
        Matched-filter S/N map of the residual.
    """
    h = outerHistory[outerIdx]
    exp = h['exposure']
    fit = h['result']['snaps'][snapIdx]['fit']
    k = h['psf'].computeKernelImage(exp.getBBox().getCenter()).array
    resid = exp.image.array - fit['modelFull']
    var = fit['variance']
    good = fit['good']
    inv = np.where(good, 1.0 / var, 0.0)
    num = fftconvolve(np.where(good, resid, 0.0) * inv, k[::-1, ::-1], mode='same')
    den = np.sqrt(np.clip(fftconvolve(inv, (k ** 2)[::-1, ::-1], mode='same'), 1e-30, None))
    scale = empiricalNoiseScale(fit, exp)
    pscale = empiricalNoiseScale(fit, exp, matched=False)
    snrMap = num / den / scale
    snrMap[~good] = 0.0
    print(f"empirical noise scale: S/N (matched) = {scale:.2f}, pixel = {pscale:.2f}")
    print(f"{'N sigma':>8s} {'>= +N':>8s} {'<= -N':>8s}")
    for t in thresholds:
        print(f"{t:8.1f} {(snrMap >= t).sum():8d} {(snrMap <= -t).sum():8d}")
    if plot:
        bbox = exp.getBBox()
        fig, ax = plt.subplots(1, 2, figsize=(14, 5.5))
        lim = max(5.0, float(np.nanpercentile(np.abs(snrMap), 99.9)))
        im = ax[0].imshow(snrMap, extent=_extent(bbox), origin='lower', cmap='RdBu_r',
                          vmin=-lim, vmax=lim)
        plt.colorbar(im, ax=ax[0], fraction=0.046, pad=0.04)
        pk = _peaksArray(h['result']['snaps'][snapIdx]['mergedArr'])
        if len(pk):
            ax[0].scatter(pk[:, 0], pk[:, 1], facecolors='none', edgecolors='k', s=30, lw=0.6)
        ax[0].set_title('Matched-filter S/N of residual (circles = fitted sources)')
        bins = np.linspace(-lim, lim, 61)
        ax[1].hist(snrMap[good], bins=bins, color='#4C72B0', alpha=0.8)
        ax[1].set_yscale('log')
        ax[1].set_xlabel('S/N of residual'); ax[1].set_title('Histogram (should be symmetric)')
        plt.tight_layout(); plt.show()
    return snrMap


def sourceSnrHistogram(outerHistory, outerIdx=-1, snapIdx=-1, minSnr=5.0):
    """Plot the histogram of the honest S/N of every final source.

    A real population shows a gap or knee; a smooth pile-up just above ``minSnr`` means noise.

    Parameters
    ----------
    outerHistory : `list` [`dict`]
        Outer-loop history.
    outerIdx : `int`, optional
        Outer pass (default: last).
    snapIdx : `int`, optional
        Inner-loop snapshot (default: last).
    minSnr : `float`, optional
        S/N threshold marked in the plot.

    Returns
    -------
    snr : `numpy.ndarray`
        S/N per source.
    """
    h = outerHistory[outerIdx]
    fit = h['result']['snaps'][snapIdx]['fit']
    scale = empiricalNoiseScale(fit, h['exposure'])
    snr = ampSnrFull(fit, scale)
    print(f"{len(snr)} sources; honest S/N: min={snr.min():.1f}, "
          f"median={np.median(snr):.1f}, max={snr.max():.0f}; "
          f"{(snr < minSnr + 2).sum()} within +2 of the cut")
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(np.log10(np.clip(snr, 1e-1, None)), bins=30, color='#55A868')
    ax.axvline(np.log10(minSnr), color='r', ls='--', label=f'minSnr={minSnr:g}')
    ax.set_xlabel('log10 honest S/N'); ax.set_ylabel('# sources'); ax.legend()
    ax.set_title('Final source S/N distribution')
    plt.tight_layout(); plt.show()
    return snr



# =====================================================================
# Final comparison plots
# =====================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import ListedColormap
from astropy.visualization import simple_norm

EDGE_COLOR = '#ff00ff'   # magenta: distinct from red (masked) and the RdBu / mako maps


def drawEdgeMargin(ax, bbox, marginPx, color=EDGE_COLOR, shade=True):
    """Outline and shade the band within ``marginPx`` of the image border.

    This is the region removed by ``dropOutOfBounds(..., margin=-marginPx)``.

    Parameters
    ----------
    ax : `matplotlib.axes.Axes`
        Axes to draw on.
    bbox : `lsst.geom.Box2I`
        Image box.
    marginPx : `float`
        Band width [px].
    color : `str`, optional
        Line colour.
    shade : `bool`, optional
        Shade the band.
    """
    if not marginPx or marginPx <= 0:
        return
    x0, x1 = bbox.getMinX(), bbox.getMaxX()
    y0, y1 = bbox.getMinY(), bbox.getMaxY()
    if shade:
        X, Y = np.meshgrid(np.arange(x0, x1 + 1), np.arange(y0, y1 + 1))
        inner = ((X >= x0 + marginPx) & (X <= x1 - marginPx) &
                 (Y >= y0 + marginPx) & (Y <= y1 - marginPx))
        band = ~inner
        ax.imshow(np.ma.masked_where(~band, np.ones(band.shape)),
                  extent=_extent(bbox), origin='lower',
                  cmap=ListedColormap([color]), alpha=0.18,
                  interpolation='nearest', zorder=3)
    ax.add_patch(patches.Rectangle((x0 + marginPx, y0 + marginPx),
                                   (x1 - x0) - 2 * marginPx,
                                   (y1 - y0) - 2 * marginPx,
                                   fill=False, ec=color, lw=1.6, ls='--', zorder=4))


def plotFinalComparison(outerHistory, exposureCrop, badBitmask=None,
                          edgeMarginFwhm=0.5, view=None,
                          showPeaks=True, circleRadius=None, figsize=(26, 7)):
    """Plot original | final model | residual | residual/sigma over the fitted region.

    Parameters
    ----------
    outerHistory : `list` [`dict`]
        Outer-loop history.
    exposureCrop : `lsst.afw.image.Exposure`
        The input crop (original image panel).
    badBitmask : `int`, optional
        Mask bits of bad pixels.
    edgeMarginFwhm : `float`, optional
        Edge margin used in the run; draws its boundary (None / 0 = do not draw).
    view : `tuple` [`float`, `float`, `float`, `float`], optional
        (xmin, xmax, ymin, ymax) zoom.
    showPeaks : `bool`, optional
        Mark the sources.
    circleRadius : `float`, optional
        Marker radius [px].
    figsize : `tuple` [`float`, `float`], optional
        Figure size [inches].
    """
    final = outerHistory[-1]
    exp = final['exposure']
    res = final['result']
    model = res['modelFull']
    peaks = _peaksArray(res['peaks'])
    bbox = exp.getBBox()

    origArr = exposureCrop[bbox].image.array
    resid = exp.image.array - model
    sigma = np.sqrt(np.where(exp.variance.array > 0, exp.variance.array, np.nan))
    zmap = resid / sigma
    bad = None if badBitmask is None else (exp.mask.array & badBitmask) != 0

    ext = _extent(bbox)
    norm = simple_norm(origArr, vmin=np.nanmin(origArr),
                       vmax=np.nanpercentile(origArr, 99),
                       stretch='log', log_a=1000)
    rsel = resid[~bad] if bad is not None and (~bad).any() else resid
    rs = np.nanpercentile(np.abs(rsel), 99)
    zs = 5.0
    radius = _circleRadius(exp, circleRadius)
    marginPx = (edgeMarginFwhm * _fwhm(exp.getPsf(), bbox.getCenter())
                 if edgeMarginFwhm else 0.0)

    fig, axes = plt.subplots(1, 4, figsize=figsize, sharex=True, sharey=True)
    panels = [
        (origArr, dict(norm=norm, cmap=CMAP), 'Original visit image (cutout)'),
        (model, dict(norm=norm, cmap=CMAP), f'Final model ({len(peaks)} sources)'),
        (resid, dict(cmap='RdBu_r', vmin=-rs, vmax=rs),
         f'Residual = data - model  (std {np.nanstd(resid):.2f})'),
        (zmap, dict(cmap='RdBu_r', vmin=-zs, vmax=zs),
         f'Residual / σ  (std {np.nanstd(zmap):.2f})'),
    ]
    for ax, (arr, kw, title) in zip(axes, panels):
        im = ax.imshow(arr, extent=ext, origin='lower', **kw)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        if bad is not None:
            _overlayBad(ax, bad, ext)
        drawEdgeMargin(ax, bbox, marginPx)
        if showPeaks and len(peaks):
            for x, y in peaks:
                ax.add_patch(patches.Circle((x, y), radius, fill=False, ec='w',
                                            lw=0.5, alpha=0.6))
        ax.set_title(title, fontsize=11)
        ax.set_xlabel('x [px]')
        _cleanAxes(ax)

    if view is not None:
        axes[0].set_xlim(view[0], view[1])
        axes[0].set_ylim(view[2], view[3])   # sharex/sharey propagates to all panels
    else:
        axes[0].set_xlim(ext[0], ext[1])
        axes[0].set_ylim(ext[2], ext[3])

    if marginPx:
        axes[0].legend(handles=[
            patches.Patch(fc=EDGE_COLOR, alpha=0.3, ec=EDGE_COLOR, ls='--',
                          label=f'edge drop zone ({edgeMarginFwhm} FWHM = {marginPx:.1f}px)'),
            patches.Patch(fc='red', alpha=0.35, label='masked (badBitmask)')],
            loc='upper left', fontsize=8)
    axes[0].set_ylabel('y [px]')
    plt.tight_layout()
    plt.show()
    plt.close(fig)



# =====================================================================
# Representative sources
# =====================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from collections import OrderedDict
from scipy.spatial import cKDTree
from matplotlib.backends.backend_pdf import PdfPages
from astropy.visualization import simple_norm


# =====================================================================
# plotSourcePages (with pdf= / show= / snr= / dpi=)
# =====================================================================

def plotSourcePages(exposure, snap, outerIdx, pages=None, colors=None,
                      sourceView='all', sourcesPerPage=6, maxPages=None,
                      viewHalf=20, circleRadius=None,
                      pdf=None, show=True, snr=None, dpi=80):
    """Plot one row per source: data | model | residual | x-section | y-section.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    snap : `dict`
        Inner-loop snapshot.
    outerIdx : `int`
        Outer-pass index.
    pages : `list` [`list` [`int`]], optional
        Source indices per page (default: `selectSourcePages`).
    colors : `dict` [`int`, `str`], optional
        Colour per source.
    sourceView : `str`, optional
        See `selectSourcePages`.
    sourcesPerPage : `int`, optional
        See `selectSourcePages`.
    maxPages : `int`, optional
        See `selectSourcePages`.
    viewHalf : `int`, optional
        Half-size of each window [px].
    circleRadius : `float`, optional
        Marker radius [px].
    pdf : `matplotlib.backends.backend_pdf.PdfPages`, optional
        Open PDF; every page is saved into it.
    show : `bool`, optional
        Also display each page (slow for hundreds of sources).
    snr : `numpy.ndarray`, optional
        Honest S/N per source, shown in the model title.
    dpi : `int`, optional
        Resolution of the PDF pages.
    """
    bbox = exposure.getBBox()
    X0, Y0 = bbox.getMinX(), bbox.getMinY()
    ny, nx = bbox.getHeight(), bbox.getWidth()
    data = exposure.image.array
    fit = snap['fit']
    model = fit['modelFull']
    resid = data - model
    peaks, amps = snap['mergedArr'], fit['amps']
    prevAmp, isNew, bad = snap['prevAmp'], snap['isNew'], snap['badMask']

    if pages is None:
        pages = selectSourcePages(snap, sourceView, sourcesPerPage, maxPages)
    if not pages:
        print("  (no sources selected for per-source pages)")
        return
    if colors is None:
        colors = sourceColors(pages)
    radius = _circleRadius(exposure, circleRadius)

    for p, page in enumerate(pages):
        n = len(page)
        fig, axes = plt.subplots(n, 5, figsize=(21, 3.3 * n), squeeze=False)
        for r, j in enumerate(page):
            col = colors[int(j)]
            px, py = peaks[j]
            comp = fit['comps'][j]
            ix, iy, xa, xb, ya, yb = _viewWindow(bbox, px, py, viewHalf)
            # clamp the cross-section row/column into the image
            ix = int(np.clip(ix, 0, nx - 1))
            iy = int(np.clip(iy, 0, ny - 1))
            if xb <= xa or yb <= ya:
                for a in axes[r]:
                    a.axis('off')
                axes[r, 0].text(0.5, 0.5, f"#{j} outside image", ha='center',
                                transform=axes[r, 0].transAxes)
                continue
            wext = (X0 + xa - 0.5, X0 + xb - 0.5, Y0 + ya - 0.5, Y0 + yb - 0.5)
            dwin, mwin, rwin = data[ya:yb, xa:xb], model[ya:yb, xa:xb], resid[ya:yb, xa:xb]
            nrm = simple_norm(dwin, vmin=np.nanmin(dwin), vmax=np.nanmax(dwin), stretch='asinh')
            rs = np.nanpercentile(np.abs(rwin), 99)
            rs = rs if np.isfinite(rs) and rs > 0 else 1.0

            own = np.zeros((ny, nx))
            ownSig = np.zeros((ny, nx), bool)
            if comp is not None:
                own[comp['ys'], comp['xs']] = comp['amp'] * comp['stamp']
                ownSig[comp['ys'], comp['xs']] = comp['sig']
            wxs, wys = X0 + np.arange(xa, xb), Y0 + np.arange(ya, yb)

            # --- image windows -------------------------------------------------
            for cI, (arr, kw) in enumerate([
                    (dwin, dict(norm=nrm, cmap=CMAP)),
                    (mwin, dict(norm=nrm, cmap=CMAP)),
                    (rwin, dict(cmap='RdBu_r', vmin=-rs, vmax=rs))]):
                ax = axes[r, cI]
                im = ax.imshow(arr, extent=wext, origin='lower', **kw)
                plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                if ownSig[ya:yb, xa:xb].any():
                    ax.contour(wxs, wys, ownSig[ya:yb, xa:xb].astype(float), levels=[0.5],
                               colors='gold', linewidths=1.2)
                ax.scatter([px], [py], marker='x', c=col, s=40, lw=1.5, zorder=6)
                _markSource(ax, px, py, col, radius, label=j if cI == 0 else None)
                inwin = ((peaks[:, 0] > wext[0]) & (peaks[:, 0] < wext[1]) &
                         (peaks[:, 1] > wext[2]) & (peaks[:, 1] < wext[3]))
                inwin[j] = False
                if inwin.any():
                    ax.scatter(peaks[inwin, 0], peaks[inwin, 1], marker='+', c='w', s=30, lw=0.9)
                _overlayBad(ax, bad[ya:yb, xa:xb], wext)
                if comp is not None:
                    x0b, x1b, y0b, y1b = comp['bboxParent']
                    ax.add_patch(patches.Rectangle((x0b - 0.5, y0b - 0.5), x1b - x0b + 1,
                                                   y1b - y0b + 1, fill=False, ec=col, lw=1.5,
                                                   ls='--', path_effects=_OUTLINE))
                ax.set_xlim(wext[0], wext[1]); ax.set_ylim(wext[2], wext[3])
                _cleanAxes(ax)
                _frameAxes(ax, col)

            g = fit['groupLabels'][j]
            tag = 'NEW' if isNew[j] else 'refit'
            dtxt = ''
            if np.isfinite(prevAmp[j]) and prevAmp[j] != 0:
                dtxt = f"  Δ={100 * (amps[j] - prevAmp[j]) / abs(prevAmp[j]):+.1f}%"
            stxt = f"  S/N={snr[j]:.1f}" if snr is not None else ''
            zflag = '  [AT ZERO BOUND]' if amps[j] <= 1e-6 else ''
            axes[r, 0].set_title(f"#{j} ({px:.1f},{py:.1f}) [{tag}] G{g}", fontsize=9, loc='left',
                                 bbox=dict(fc=col, ec='k', pad=2))
            axes[r, 1].set_title(f"model  amp={amps[j]:.1f}{dtxt}{stxt}{zflag}",
                                 fontsize=9, loc='left')
            axes[r, 2].set_title(f"residual  (group χ²ν={fit['groupChi2'].get(g, np.nan):.2f})",
                                 fontsize=9, loc='left')
            if r == 0:
                axes[r, 0].text(0.02, 0.02,
                                'circle = this source\ndashed box = PSF stamp\ngold = significant support',
                                transform=axes[r, 0].transAxes, fontsize=6, color='w',
                                bbox=dict(fc='k', alpha=0.5, lw=0))

            # --- cross-sections -----------------------------------------------
            for cI, axis in ((3, 'x'), (4, 'y')):
                ax = axes[r, cI]
                if axis == 'x':
                    coord, sl, centre = wxs, (iy, slice(xa, xb)), px
                    c0, c1 = ((comp['bboxParent'][0], comp['bboxParent'][1]) if comp else (None, None))
                else:
                    coord, sl, centre = wys, (slice(ya, yb), ix), py
                    c0, c1 = ((comp['bboxParent'][2], comp['bboxParent'][3]) if comp else (None, None))
                badline = bad[sl]
                ax.plot(coord, data[sl], color='k', lw=1.2, label='data')
                ax.plot(coord, model[sl], color='#DD8452', lw=1.5, ls='--', label='total model')
                ax.plot(coord, own[sl], color='#1f77b4', lw=1.3, label="this source's PSF")
                ax.plot(coord, model[sl] - own[sl], color='#2ca02c', lw=1.0, ls=':', label='neighbours')
                ax.plot(coord, resid[sl], color='gray', lw=0.9, alpha=0.8, label='residual')
                ax.axhline(0, color='gray', lw=0.5)
                ax.axvline(centre, color=col, lw=2.0, path_effects=_OUTLINE)
                for s in (-radius, radius):
                    ax.axvline(centre + s, color=col, lw=1.3, ls=':', path_effects=_OUTLINE)
                if c0 is not None:
                    ax.axvspan(c0 - 0.5, c1 + 0.5, color=col, alpha=0.10)
                ax.fill_between(coord, 0, 1, where=ownSig[sl], transform=ax.get_xaxis_transform(),
                                color='gold', alpha=0.15)
                if badline.any():
                    ax.plot(coord[badline], np.zeros(int(badline.sum())), '|', color='red',
                            transform=ax.get_xaxis_transform(), ms=8)
                ax.set_xlim(coord[0], coord[-1])
                ax.set_xlabel(f'{axis} [px]')
                ax.set_title(f"{axis}-cross-section through source #{j}", fontsize=9)
                _frameAxes(ax, col)
                if r == 0 and cI == 3:
                    ax.legend(fontsize=6, loc='upper right')

        fig.suptitle(f"Outer {outerIdx + 1} / inner iter {snap['it'] + 1}: per-source fits "
                     f"(page {p + 1}/{len(pages)}, view={sourceView}); "
                     f"colours match overview panel 2", fontsize=13, y=1.0)
        plt.tight_layout()
        if pdf is not None:
            pdf.savefig(fig, dpi=dpi)
        if show:
            plt.show()
        plt.close(fig)


# =====================================================================
# Per-source metrics
# =====================================================================

def sourceMetrics(outerHistory, outerIdx=-1, snapIdx=-1, coreFwhm=1.5):
    """Return per-source quantities used to stratify the sample.

    Parameters
    ----------
    outerHistory : `list` [`dict`]
        Outer-loop history.
    outerIdx : `int`, optional
        Outer pass (default: last).
    snapIdx : `int`, optional
        Inner-loop snapshot (default: last).
    coreFwhm : `float`, optional
        Radius of the chi2 region [FWHM].

    Returns
    -------
    metrics : `dict`
        ``peaks``, ``amps``, ``snr`` (honest S/N), ``nn`` (distance to the nearest other source
        [px]), ``edge`` (distance to the image border [px]), ``chi2`` (mean resid^2 / sigma^2
        within ``coreFwhm`` FWHM), ``group`` (size of the coupled group), ``fwhm``,
        ``outerPos`` and the exposure and snapshot used.
    """
    outerPos = outerIdx if outerIdx >= 0 else len(outerHistory) + outerIdx
    h = outerHistory[outerPos]
    exposure = h['exposure']
    snap = h['result']['snaps'][snapIdx]
    fit = snap['fit']
    peaks = np.asarray(snap['mergedArr'], float)
    bbox = exposure.getBBox()
    X0, Y0 = bbox.getMinX(), bbox.getMinY()
    ny, nx = bbox.getHeight(), bbox.getWidth()

    snr = ampSnrFull(fit, empiricalNoiseScale(fit, exposure))
    fwhm = _fwhm(exposure.getPsf(), bbox.getCenter())
    r = coreFwhm * fwhm

    if len(peaks) > 1:
        d, _ = cKDTree(peaks).query(peaks, k=2)
        nn = d[:, 1]
    else:
        nn = np.full(len(peaks), np.inf)

    edge = np.minimum.reduce([peaks[:, 0] - bbox.getMinX(), bbox.getMaxX() - peaks[:, 0],
                              peaks[:, 1] - bbox.getMinY(), bbox.getMaxY() - peaks[:, 1]])

    resid = exposure.image.array - fit['modelFull']
    var, good = fit['variance'], fit['good']
    chi2 = np.full(len(peaks), np.nan)
    rr = int(np.ceil(r))
    for j, (px, py) in enumerate(peaks):
        ix, iy = int(round(px)) - X0, int(round(py)) - Y0
        xa, xb = max(ix - rr, 0), min(ix + rr + 1, nx)
        ya, yb = max(iy - rr, 0), min(iy + rr + 1, ny)
        if xb <= xa or yb <= ya:
            continue
        yy, xx = np.mgrid[ya:yb, xa:xb]
        sel = good[ya:yb, xa:xb] & (((xx + X0 - px) ** 2 + (yy + Y0 - py) ** 2) <= r ** 2)
        if sel.sum() > 5:
            chi2[j] = np.mean(resid[ya:yb, xa:xb][sel] ** 2 / var[ya:yb, xa:xb][sel])

    group = np.bincount(fit['groupLabels'])[fit['groupLabels']]
    return dict(peaks=peaks, amps=fit['amps'], snr=snr, nn=nn, edge=edge, chi2=chi2,
                group=group, fwhm=fwhm, outerPos=outerPos, exposure=exposure,
                snap=snap, fit=fit)


# =====================================================================
# Selection
# =====================================================================

def selectRepresentativeSources(m, nPerCategory=4, seed=0, categories=None,
                                  extraIds=None):
    """Select a representative sample of sources per category.

    Parameters
    ----------
    m : `dict`
        Output of `sourceMetrics`.
    nPerCategory : `int` or `dict` [`str`, `int`], optional
        Sources drawn per category.
    seed : `int`, optional
        Random seed of the tercile draws.
    categories : `list` [`str`], optional
        Categories to keep, in display order.
    extraIds : `list` [`int`], optional
        Source ids appended as a 'picked' category.

    Returns
    -------
    selection : `collections.OrderedDict` [`str`, `list` [`int`]]
        Source ids per category.
    """
    rng = np.random.default_rng(seed)
    snr, amps, nn, edge, chi2 = m['snr'], m['amps'], m['nn'], m['edge'], m['chi2']
    nSrc = len(snr)
    used = set()

    def nOf(cat):
        return nPerCategory.get(cat, 4) if isinstance(nPerCategory, dict) else nPerCategory

    def take(cands, n):
        out = []
        for j in cands:
            j = int(j)
            if j in used:
                continue
            out.append(j)
            used.add(j)
            if len(out) >= n:
                break
        return out

    sel = OrderedDict()
    sel['brightest'] = take(np.argsort(-amps), nOf('brightest'))

    lo, hi = np.quantile(snr, [1 / 3, 2 / 3])
    terciles = OrderedDict([('bright', np.flatnonzero(snr >= hi)),
                            ('mid', np.flatnonzero((snr >= lo) & (snr < hi))),
                            ('faint', np.flatnonzero(snr < lo))])
    for name, idx in terciles.items():
        sel[name] = take(rng.permutation(idx), nOf(name))

    sel['marginal (lowest S/N)'] = take(np.argsort(snr), nOf('marginal (lowest S/N)'))
    sel['crowded (nearest neighbour)'] = take(np.argsort(nn), nOf('crowded (nearest neighbour)'))

    chiWorst = np.where(np.isfinite(chi2), chi2, -np.inf)
    sel['worst fit (core chi2)'] = take(np.argsort(-chiWorst), nOf('worst fit (core chi2)'))

    brightHalf = snr >= np.median(snr)
    chiBest = np.where(np.isfinite(chi2) & brightHalf, chi2, np.inf)
    sel['best fit (bright half)'] = take(np.argsort(chiBest), nOf('best fit (bright half)'))

    sel['near edge'] = take(np.argsort(edge), nOf('near edge'))

    if extraIds is not None:
        sel['picked'] = take(np.asarray(extraIds, int), nSrc)

    if categories is not None:
        sel = OrderedDict((c, sel[c]) for c in categories if c in sel)
    return OrderedDict((c, ids) for c, ids in sel.items() if ids)


# =====================================================================
# Plot
# =====================================================================

def plotRepresentativeSources(outerHistory, outerIdx=-1, snapIdx=-1,
                                nPerCategory=4, seed=0, categories=None,
                                extraIds=None, sourcesPerPage=4, viewHalf=20,
                                pdfPath=None, show=True, dpi=80):
    """Select a representative sample and plot it with `plotSourcePages`.

    Each category is split into pages of ``sourcesPerPage`` rows, so a page never mixes
    categories; the printed table maps every source id to its category.

    Parameters
    ----------
    outerHistory : `list` [`dict`]
        Outer-loop history.
    outerIdx : `int`, optional
        Outer pass (default: last).
    snapIdx : `int`, optional
        Inner-loop snapshot (default: last).
    nPerCategory : `int` or `dict`, optional
        See `selectRepresentativeSources`.
    seed : `int`, optional
        See `selectRepresentativeSources`.
    categories : `list` [`str`], optional
        See `selectRepresentativeSources`.
    extraIds : `list` [`int`], optional
        See `selectRepresentativeSources`.
    sourcesPerPage : `int`, optional
        Rows per page.
    viewHalf : `int`, optional
        Half-size of each window [px].
    pdfPath : `str`, optional
        Also write all pages to this PDF.
    show : `bool`, optional
        Display the pages.
    dpi : `int`, optional
        Resolution of the PDF pages.

    Returns
    -------
    selection : `collections.OrderedDict` [`str`, `list` [`int`]]
        Source ids per category.
    """
    m = sourceMetrics(outerHistory, outerIdx, snapIdx)
    sel = selectRepresentativeSources(m, nPerCategory, seed, categories, extraIds)

    # console table
    print(f"{'category':30s} {'id':>4s} {'x':>7s} {'y':>7s} {'amp':>10s} {'S/N':>7s} "
          f"{'nn[px]':>7s} {'chi2':>6s} {'grp':>4s}")
    for cat, ids in sel.items():
        for j in ids:
            x, y = m['peaks'][j]
            print(f"{cat:30s} {j:4d} {x:7.1f} {y:7.1f} {m['amps'][j]:10.1f} {m['snr'][j]:7.1f} "
                  f"{m['nn'][j]:7.2f} {m['chi2'][j]:6.2f} {m['group'][j]:4d}")

    pages, pageCat = [], []
    for cat, ids in sel.items():
        for i in range(0, len(ids), sourcesPerPage):
            pages.append(ids[i:i + sourcesPerPage])
            pageCat.append(cat)
    print("\npage -> category:")
    for p, cat in enumerate(pageCat, 1):
        print(f"  page {p:2d}: {cat}  (ids {pages[p - 1]})")
    if not pages:
        print("no sources selected")
        return sel

    colors = sourceColors(pages)
    kw = dict(pages=pages, colors=colors, viewHalf=viewHalf, show=show,
              snr=m['snr'], dpi=dpi)
    if pdfPath is not None:
        with PdfPages(pdfPath) as pdf:
            plotSourcePages(m['exposure'], m['snap'], m['outerPos'], pdf=pdf, **kw)
        print(f"saved {len(pages)} page(s) to {pdfPath}")
    else:
        plotSourcePages(m['exposure'], m['snap'], m['outerPos'], pdf=None, **kw)
    return sel



# old snake_case names still work (with a FutureWarning)
from ._compat import moduleGetattr  # noqa: E402
__getattr__ = moduleGetattr(globals())
