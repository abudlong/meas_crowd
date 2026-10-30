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
background re-determination, the iterative inner / outer loop of `run_cutout_pipeline`, and
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



# ######################################################################
# [from notebook cell: inner loop helpers]
# ######################################################################

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

def bbox_to_slices(bbox, parent_bbox):
    """Convert a box into array slices relative to a parent box.

    Parameters
    ----------
    bbox : `lsst.geom.Box2I`
        Box to convert.
    parent_bbox : `lsst.geom.Box2I`
        Box of the array the slices index (its origin is pixel [0, 0]).

    Returns
    -------
    slices : `tuple` [`slice`, `slice`]
        (y slice, x slice) selecting ``bbox`` from an array covering ``parent_bbox``.
    """
    x0 = bbox.getMinX() - parent_bbox.getMinX()
    y0 = bbox.getMinY() - parent_bbox.getMinY()
    x1 = x0 + bbox.getWidth()
    y1 = y0 + bbox.getHeight()
    return slice(y0, y1), slice(x0, x1)


def collapse_saturated_peaks(peak_positions, exposure, sat_bitmask=None, dilate_iters=2):
    """Collapse every peak on or near a saturated blob into one peak per blob.

    The representative position is a flux-weighted centroid from the unsaturated ring just
    outside the saturated core (the flat top carries no usable position or flux information).
    Without this, noise on the saturated flat top produces many near-duplicate peaks whose PSF
    stamps are almost identical; the joint fit then splits the true amplitude between them
    arbitrarily, which shows up as a blocky, high-frequency patch in the model image.

    Parameters
    ----------
    peak_positions : `list` [`tuple` [`float`, `float`]]
        Peak (x, y) positions in parent pixel coordinates.
    exposure : `lsst.afw.image.Exposure`
        Exposure the peaks were detected on.
    sat_bitmask : `int`, optional
        Mask bits of saturated pixels (default: the SAT plane).
    dilate_iters : `int`, optional
        Number of dilations of the saturated region when assigning peaks to a blob.

    Returns
    -------
    peaks : `list` [`tuple` [`float`, `float`]]
        Peaks with every saturated blob reduced to one representative peak.
    n_removed : `int`
        Number of peaks removed.
    """
    if sat_bitmask is None:
        sat_bitmask = exposure.mask.getPlaneBitMask("SAT")

    mask = exposure.mask.array
    data = exposure.image.array
    xy0 = exposure.getXY0()
    x0, y0 = xy0.getX(), xy0.getY()

    sat = (mask & sat_bitmask) != 0
    if not sat.any():
        return list(peak_positions), 0

    labels, n_labels = label(sat)
    dilated = binary_dilation(sat, iterations=dilate_iters)
    ring = dilated & ~sat

    blob_centroids = {}
    for lab in range(1, n_labels + 1):
        blob_mask = labels == lab
        blob_ring = binary_dilation(blob_mask, iterations=dilate_iters) & ring
        if blob_ring.sum() == 0:
            cy, cx = center_of_mass(blob_mask)
        else:
            weights = np.clip(data[blob_ring], 0, None)
            ys, xs = np.nonzero(blob_ring)
            if weights.sum() <= 0:
                cy, cx = center_of_mass(blob_mask)
            else:
                cy = np.average(ys, weights=weights)
                cx = np.average(xs, weights=weights)
        blob_centroids[lab] = (cx + x0, cy + y0)

    collapsed = []
    n_removed = 0
    seen_blobs = set()
    for (x, y) in peak_positions:
        px = int(round(x - x0))
        py = int(round(y - y0))
        in_bounds = 0 <= py < sat.shape[0] and 0 <= px < sat.shape[1]
        lab = labels[py, px] if in_bounds and sat[py, px] else 0

        if lab == 0:
            collapsed.append((x, y))
            continue

        n_removed += 1
        if lab not in seen_blobs:
            seen_blobs.add(lab)
            collapsed.append(blob_centroids[lab])

    print(f"  collapsed {n_removed} saturated-core peaks -> "
          f"{len(seen_blobs)} representative peak(s) across {n_labels} saturated blob(s)")
    return collapsed, n_removed


def merge_close_peaks(peak_positions, exposure, min_separation=None):
    """Merge peaks closer than ``min_separation``, keeping the brightest pixel's peak.

    Prevents near-duplicate columns in the design matrix (mostly from re-detecting the same
    wings) from destabilising the joint solve.

    Parameters
    ----------
    peak_positions : `list` [`tuple` [`float`, `float`]]
        Peak (x, y) positions in parent pixel coordinates.
    exposure : `lsst.afw.image.Exposure`
        Exposure whose pixel values rank the peaks.
    min_separation : `float`, optional
        Merge radius [px] (default: twice the PSF sigma).

    Returns
    -------
    merged_positions : `list` [`tuple` [`float`, `float`]]
        Surviving peak positions.
    """
    if len(peak_positions) == 0:
        return []
    if min_separation is None:
        shape = exposure.getPsf().computeShape(exposure.getBBox().getCenter())
        psf_sigma = np.sqrt(0.5 * (shape.getIxx() + shape.getIyy()))
        min_separation = 2 * psf_sigma

    positions = np.array(peak_positions)
    data = exposure.image.array
    xy0 = exposure.getXY0()

    peak_vals = np.array([
        data[int(y - xy0.getY()), int(x - xy0.getX())]
        if 0 <= int(y - xy0.getY()) < data.shape[0] and 0 <= int(x - xy0.getX()) < data.shape[1]
        else -np.inf
        for x, y in positions
    ])

    tree = cKDTree(positions)
    pairs = tree.query_pairs(r=min_separation)

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

    merged_positions = []
    n_merged = 0
    for members in clusters.values():
        if len(members) > 1:
            n_merged += len(members) - 1
        best = members[np.argmax(peak_vals[members])]
        merged_positions.append(tuple(positions[best]))

    print(f"  merged {n_merged} duplicate/spurious peaks -- "
          f"{len(peak_positions)} -> {len(merged_positions)}")
    return merged_positions


def joint_psf_fit(exposure, psf, peak_positions, bad_bitmask, amp_zero_tol=1e-6):
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
    peak_positions : `list` [`tuple` [`float`, `float`]]
        Source (x, y) positions in parent pixel coordinates.
    bad_bitmask : `int`
        Mask bits of pixels to ignore.
    amp_zero_tol : `float`, optional
        Amplitudes at or below this count as zero.

    Returns
    -------
    amps : `numpy.ndarray`
        Fitted amplitudes (fluxes), one per position.
    model_full : `numpy.ndarray`
        Model image over the exposure's box.
    n_at_zero_bound : `int`
        Number of sources that converged to (numerically) zero amplitude -- sources that still
        fail to fit rather than being absent.
    """
    exp_bbox = exposure.getBBox()
    ny, nx = exp_bbox.getHeight(), exp_bbox.getWidth()
    data = exposure.image.array
    variance = exposure.variance.array
    mask = exposure.mask.array

    good_global = (mask & bad_bitmask) == 0
    good_global &= np.isfinite(data) & np.isfinite(variance) & (variance > 0)

    pixel_index = -np.ones((ny, nx), dtype=np.int64)
    flat_positions = np.where(good_global.ravel())[0]
    pixel_index.ravel()[flat_positions] = np.arange(len(flat_positions))
    n_pixels_used = len(flat_positions)

    d_vec = data.ravel()[flat_positions]
    w_vec = 1.0 / variance.ravel()[flat_positions]

    A = lil_matrix((n_pixels_used, len(peak_positions)))
    for j, (x, y) in enumerate(peak_positions):
        point = geom.Point2D(x, y)
        psfImage = psf.computeImage(point)
        bbox = psfImage.getBBox()
        bbox.clip(exp_bbox)
        if bbox.isEmpty():
            continue
        model_stamp = np.clip(psfImage[bbox].array, 0, None)

        y0 = bbox.getMinY() - exp_bbox.getMinY()
        x0 = bbox.getMinX() - exp_bbox.getMinX()
        nz_sy, nz_sx = np.nonzero(model_stamp)
        for sy, sx in zip(nz_sy, nz_sx):
            gy, gx = y0 + sy, x0 + sx
            if 0 <= gy < ny and 0 <= gx < nx:
                idx = pixel_index[gy, gx]
                if idx >= 0:
                    A[idx, j] = model_stamp[sy, sx]

    A = csr_matrix(A)
    sqrt_w = np.sqrt(w_vec)
    A_weighted = A.multiply(sqrt_w[:, None])
    d_weighted = d_vec * sqrt_w

    result = lsq_linear(A_weighted, d_weighted, bounds=(0, np.inf),
                         method='trf', lsmr_tol='auto', max_iter=2000)
    amps = result.x

    n_at_zero_bound = int(np.sum(amps <= amp_zero_tol))

    model_full = np.zeros((ny, nx))
    for j, (x, y) in enumerate(peak_positions):
        if amps[j] <= amp_zero_tol:
            continue
        point = geom.Point2D(x, y)
        psfImage = psf.computeImage(point)
        bbox = psfImage.getBBox()
        bbox.clip(exp_bbox)
        if bbox.isEmpty():
            continue
        model_stamp = np.clip(psfImage[bbox].array, 0, None)
        yslice, xslice = bbox_to_slices(bbox, exp_bbox)
        model_full[yslice, xslice] += amps[j] * model_stamp

    return amps, model_full, n_at_zero_bound


# =====================================================================
# Main iterative detection + joint fit loop
# =====================================================================

def iterative_joint_detection(
    exposure,
    psf,
    bad_bitmask,
    n_iterations=10,
    threshold_start=8.0,
    threshold_step=-1.5,
    threshold_min=1.5,
    min_pixels=1,
    min_new_peaks=2,          # stop if a pass finds fewer new peaks than this
    convergence_frac=0.02,    # stop if residual std improves by less than this fraction
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
    bad_bitmask : `int`
        Mask bits of pixels to ignore.
    n_iterations : `int`, optional
        Maximum number of passes.
    threshold_start : `float`, optional
        Detection threshold [sigma] of the first pass.
    threshold_step : `float`, optional
        Change of the threshold per pass (negative: decreasing).
    threshold_min : `float`, optional
        Lowest detection threshold [sigma].
    min_pixels : `int`, optional
        Minimum footprint size of a detection.
    min_new_peaks : `int`, optional
        Stop when a pass finds fewer new peaks than this.
    convergence_frac : `float`, optional
        Stop when the residual standard deviation improves by less than this fraction.
    verbose : `bool`, optional
        Print progress.

    Returns
    -------
    result : `dict`
        ``model_full``, ``residual_full``, ``peaks``, ``amps``, the per-pass histories
        (``resid_std_history``, ``n_peaks_history``, ``n_zero_amp_history``) and ``first_seen``
        (pass index at which each peak, keyed by its position rounded to 0.1 px, was found).
    """
    exp_bbox = exposure.getBBox()

    working_residual = exposure.clone()
    working_residual.setPsf(psf)

    all_peaks = []
    resid_std_history = [np.nanstd(working_residual.image.array)]
    n_peaks_history = [0]
    n_zero_amp_history = [0]
    first_seen = {}   # (rounded x, rounded y) -> iteration index

    threshold = threshold_start
    model_full = np.zeros((exp_bbox.getHeight(), exp_bbox.getWidth()))
    amps_final = np.array([])

    for it in range(n_iterations):
        if verbose:
            print(f"--- Pass {it+1}/{n_iterations}  threshold={threshold:.2f} sigma  "
                  f"residual std={resid_std_history[-1]:.3f} nJy  "
                  f"(peaks so far={len(all_peaks)}) ---")

        det_config = SourceDetectionTask.ConfigClass()
        det_config.thresholdValue = threshold
        det_config.thresholdType = "pixel_stdev"
        det_config.minPixels = min_pixels
        det_config.doTempLocalBackground = True
        det_config.tempLocalBackground.binSize = 16
        det_config.reEstimateBackground = False  # residual has no real sky background left to fit

        schema_it = afwTable.SourceTable.makeMinimalSchema()
        detectionTask_it = SourceDetectionTask(schema=schema_it, config=det_config)
        table_it = afwTable.SourceTable.make(schema_it)

        try:
            detRes_it = detectionTask_it.run(table_it, working_residual)
        except TooManyMaskedPixelsError:
            if verbose:
                print("  background estimation failed on residual -- stopping")
            break

        new_peaks = [
            (peak.getFx(), peak.getFy())
            for fp in detRes_it.positive.getFootprints()
            for peak in fp.getPeaks()
        ]
        if verbose:
            print(f"  {len(new_peaks)} raw peaks detected on residual this pass")

        # collapse spurious multi-peak detections on saturated cores
        new_peaks, n_collapsed = collapse_saturated_peaks(new_peaks, exposure)

        if len(new_peaks) < min_new_peaks and it > 0:
            if verbose:
                print("  too few new peaks -- treating as converged, stopping")
            break

        merged_peaks = merge_close_peaks(all_peaks + new_peaks, exposure)

        amps, model_full, n_at_zero_bound = joint_psf_fit(exposure, psf, merged_peaks, bad_bitmask)

        # record which pass first introduced each surviving peak position
        for (x, y) in merged_peaks:
            key = (round(x, 1), round(y, 1))
            if key not in first_seen:
                first_seen[key] = it

        working_residual = exposure.clone()
        working_residual.image.array[:] = exposure.image.array - model_full
        working_residual.setPsf(psf)

        all_peaks = merged_peaks
        amps_final = amps
        n_peaks_history.append(len(all_peaks))
        n_zero_amp_history.append(n_at_zero_bound)

        new_resid_std = np.nanstd(working_residual.image.array)
        old_resid_std = resid_std_history[-1]
        resid_std_history.append(new_resid_std)

        if verbose:
            print(f"  refit {len(all_peaks)} peaks total ({n_at_zero_bound} at zero amplitude) -> "
                  f"residual std {old_resid_std:.3f} -> {new_resid_std:.3f}")
            if n_at_zero_bound > 0:
                print(f"  NOTE: {n_at_zero_bound} peak(s) fit to zero amplitude this pass -- "
                      f"their flux is still sitting in the residual, unmodeled")

        if old_resid_std > 0 and it > 0:
            improvement = abs(old_resid_std - new_resid_std) / old_resid_std
            if improvement < convergence_frac:
                if verbose:
                    print(f"  residual std improved by only {improvement*100:.2f}% "
                          f"(< {convergence_frac*100:.0f}%) -- converged, stopping")
                break

        threshold = max(threshold + threshold_step, threshold_min)

    return {
        'model_full': model_full,
        'residual_full': working_residual.image.array.copy(),
        'peaks': all_peaks,
        'amps': amps_final,
        'resid_std_history': resid_std_history,
        'n_peaks_history': n_peaks_history,
        'n_zero_amp_history': n_zero_amp_history,
        'first_seen': first_seen,
    }


# =====================================================================
# Tiled data / model / residual comparison
# =====================================================================

def plot_tiled_comparison(
    data_full,
    model_full,
    residual_full,
    n_sources=None,
    n_passes=None,
    tile_size=100,
    cmap_img=None,
    global_scale=True,
    skip_empty=True,
    empty_std_thresh=1e-6,
    figsize_per_panel=5,
    show_row_profile=True,
    row_select='max_flux',
    stretch='log'
):
    """Plot data, model and residual (and a row profile) for every tile of an image.

    The data and model panels share one colour scale when ``global_scale`` is True; the
    residual panel always scales to its own tile.  With ``show_row_profile`` a fourth panel
    shows data and model along one row of the tile (marked with a dashed cyan line), so the
    model can be compared with the data through a crowded section.

    Parameters
    ----------
    data_full, model_full, residual_full : `numpy.ndarray`
        Data, model and residual images (same shape).
    n_sources : `int`, optional
        Number of sources, for the title.
    n_passes : `int`, optional
        Number of detection passes, for the title.
    tile_size : `int`, optional
        Tile size [px].
    cmap_img : `str`, optional
        Colour map of the data and model panels (default: 'mako' if available, else 'viridis').
    global_scale : `bool`, optional
        One colour scale for all data / model panels instead of one per tile.
    skip_empty : `bool`, optional
        Skip tiles whose data standard deviation is below ``empty_std_thresh``.
    empty_std_thresh : `float`, optional
        See ``skip_empty``.
    figsize_per_panel : `float`, optional
        Size of one panel [inches].
    show_row_profile : `bool`, optional
        Add the row-profile panel.
    row_select : `str`, `int` or callable, optional
        Row of each tile to profile: 'max_flux' (the row with the most flux), a row index local
        to the tile, or a function ``row_select(tile_data) -> row index``.
    stretch : `str`, optional
        Stretch of the data / model colour scale.
    """
    if cmap_img is None:
        cmap_img = 'mako' if 'mako' in plt.colormaps() else 'viridis'

    ny, nx = data_full.shape

    if global_scale:
        vmax_global = np.nanpercentile(data_full, 99.5)
        norm_global = simple_norm(data_full, vmin=np.nanmin(data_full),
                                   vmax=vmax_global, stretch=stretch)

    y_edges = list(range(0, ny, tile_size))
    x_edges = list(range(0, nx, tile_size))

    title_bits = []
    if n_sources is not None:
        title_bits.append(f"{n_sources} sources")
    if n_passes is not None:
        title_bits.append(f"{n_passes} passes")
    suffix = f" ({', '.join(title_bits)})" if title_bits else ""

    n_panels = 4 if show_row_profile else 3

    for y0 in y_edges:
        for x0 in x_edges:
            y1 = min(y0 + tile_size, ny)
            x1 = min(x0 + tile_size, nx)

            tile_data = data_full[y0:y1, x0:x1]

            if skip_empty and np.nanstd(tile_data) < empty_std_thresh:
                continue

            tile_model = model_full[y0:y1, x0:x1]
            tile_resid = residual_full[y0:y1, x0:x1]

            if global_scale:
                norm = norm_global
            else:
                vmax = np.nanpercentile(tile_data, 99.5)
                norm = simple_norm(tile_data, vmin=np.nanmin(tile_data),
                                    vmax=vmax, stretch='log')

            # residual scale: ALWAYS local to this tile
            res_scale = np.nanpercentile(np.abs(tile_resid), 99)
            if not np.isfinite(res_scale) or res_scale == 0:
                res_scale = 1.0

            fig, axes = plt.subplots(1, n_panels,
                                      figsize=(figsize_per_panel * n_panels, figsize_per_panel))
            fig.patch.set_facecolor('white')

            im0 = axes[0].imshow(tile_data, origin='lower', norm=norm, cmap=cmap_img)
            axes[0].set_title('Data')
            fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

            im1 = axes[1].imshow(tile_model, origin='lower', norm=norm, cmap=cmap_img)
            axes[1].set_title('Joint model')
            fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

            im2 = axes[2].imshow(tile_resid, origin='lower', cmap='RdBu_r',
                                  vmin=-res_scale, vmax=res_scale)
            axes[2].set_title('Residual')
            fig.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04)

            for ax in axes[:3]:
                ax.set_xticks([]); ax.set_yticks([])
                ax.grid(False)

            if show_row_profile:
                if row_select == 'max_flux':
                    row_idx_local = int(np.nanargmax(np.nansum(tile_data, axis=1)))
                elif callable(row_select):
                    row_idx_local = int(row_select(tile_data))
                else:
                    row_idx_local = int(row_select)
                row_idx_local = int(np.clip(row_idx_local, 0, tile_data.shape[0] - 1))
                row_idx_global = y0 + row_idx_local

                # mark the sliced row on the data/model panels for context
                axes[0].axhline(row_idx_local, color='cyan', lw=0.8, ls='--')
                axes[1].axhline(row_idx_local, color='cyan', lw=0.8, ls='--')

                x_pix = np.arange(x0, x1)
                ax3 = axes[3]
                ax3.plot(x_pix, tile_data[row_idx_local, :], color='k', lw=1.0, label='data')
                ax3.plot(x_pix, tile_model[row_idx_local, :], color='#DD8452', lw=1.3, label='model')
                ax3.plot(x_pix, tile_resid[row_idx_local, :], color='#4C72B0', lw=0.8,
                         alpha=0.7, label='residual')
                ax3.axhline(0, color='gray', lw=0.5, ls=':')
                ax3.set_title(f'Row {row_idx_global} cross-section')
                ax3.set_xlabel('x [px]')
                ax3.legend(fontsize=7, loc='upper right')

            fig.suptitle(f'Iterative joint fit{suffix} — tile x[{x0}:{x1}] y[{y0}:{y1}]',
                         fontsize=12, fontweight='medium', y=1.03)
            plt.tight_layout()
            plt.show()


# =====================================================================
# Per-iteration detections overlay
# =====================================================================

def plot_detections_by_iteration(exposure_crop, result_iter):
    """Plot every fitted peak over the image, coded by the pass that first found it.

    Peaks whose final amplitude is (numerically) zero -- they hit the non-negativity bound --
    are marked as failed fits with a red 'x'; they are the ones to inspect if sources seem to
    be missing from the model.

    Parameters
    ----------
    exposure_crop : `lsst.afw.image.Exposure`
        Image the pipeline ran on.
    result_iter : `dict`
        Output of `iterative_joint_detection`.

    Returns
    -------
    sources_tbl : `astropy.table.Table`
        Positions, amplitudes and detection pass of the sources.
    """
    xy0 = exposure_crop.getXY0()  # (x0, y0) of this crop in the original exposure's frame
    x0, y0 = xy0.getX(), xy0.getY()

    # --- assemble a sources table from the pipeline result ---
    positions = np.array(result_iter['peaks'])
    amps = result_iter['amps']
    first_seen = result_iter['first_seen']

    iterations = np.array([
        first_seen.get((round(x, 1), round(y, 1)), -1)
        for x, y in positions
    ])

    sources_tbl = Table()
    sources_tbl['slot_Centroid_x'] = positions[:, 0]
    sources_tbl['slot_Centroid_y'] = positions[:, 1]
    sources_tbl['amp_fit'] = amps
    sources_tbl['iteration'] = iterations

    n_iters_run = len(result_iter['resid_std_history']) - 1
    iter_colors = cm.plasma(np.linspace(0, 1, max(n_iters_run, 1)))
    markers = ['o', 's', '^', 'D', 'v', 'P', '*', 'X']

    # --- plot ---
    fig, ax = plt.subplots(figsize=(8, 8))
    data_crop = exposure_crop.image.array
    vmax = np.nanpercentile(data_crop, 99.5)
    norm = simple_norm(data_crop, vmin=np.nanmin(data_crop), vmax=vmax, stretch='log')
    ax.imshow(data_crop, origin='lower', norm=norm, cmap='gray')

    for it in range(n_iters_run):
        mask_iter = sources_tbl['iteration'] == it
        injected = mask_iter & (sources_tbl['amp_fit'] > 0)
        failed = mask_iter & ~(sources_tbl['amp_fit'] > 0)

        # shift into array-relative coords by subtracting the crop's XY0
        x_inj = sources_tbl['slot_Centroid_x'][injected] - x0
        y_inj = sources_tbl['slot_Centroid_y'][injected] - y0

        ax.scatter(x_inj, y_inj,
                   s=25, facecolors='none', edgecolors=[iter_colors[it]],
                   marker=markers[it % len(markers)], linewidths=1.3,
                   label=f'pass {it+1} (n={injected.sum()})')

        if failed.any():
            x_fail = sources_tbl['slot_Centroid_x'][failed] - x0
            y_fail = sources_tbl['slot_Centroid_y'][failed] - y0
            ax.scatter(x_fail, y_fail,
                       s=25, color='red', marker='x', linewidths=1.3,
                       label=f'pass {it+1} failed fit (n={failed.sum()})' if it == 0 else None)

    ax.set_xlim(0, data_crop.shape[1])
    ax.set_ylim(0, data_crop.shape[0])
    ax.legend(fontsize=8, loc='upper left', bbox_to_anchor=(1.02, 1.0), framealpha=0.9)
    ax.set_title('Detections by iteration pass')
    ax.set_xticks([]); ax.set_yticks([])
    plt.tight_layout()
    plt.show()

    return sources_tbl


# ######################################################################
# [from notebook cell: outer loop]
# ######################################################################

def mask_sources_from_model(exposure, model_full, noise_sigma=None,
                             n_sigma=4.0, dilate_iters=1, verbose=True):
    """Flag pixels where the joint-fit model has significant flux as DETECTED.

    The next background fit then treats them as source, not sky.  The joint model contains all
    sources found so far, including ones too faint for the first background pass.  A loose
    threshold flags most of the image and starves the background fit; if the printed masked
    fraction climbs much above ~50 %, tighten ``n_sigma`` or ``dilate_iters``.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure whose mask is modified in place.
    model_full : `numpy.ndarray`
        Model image over the exposure's box.
    noise_sigma : `float`, optional
        Pixel noise (default: sqrt of the median variance).
    n_sigma : `float`, optional
        Model pixels above ``n_sigma * noise_sigma`` are flagged.
    dilate_iters : `int`, optional
        Number of dilations of the flagged region.
    verbose : `bool`, optional
        Print the masked fraction.

    Returns
    -------
    src_mask : `numpy.ndarray` [`bool`]
        The flagged pixels.
    """
    if noise_sigma is None:
        noise_sigma = np.sqrt(np.nanmedian(exposure.variance.array))

    src_mask = model_full > (n_sigma * noise_sigma)
    if dilate_iters > 0:
        src_mask = binary_dilation(src_mask, iterations=dilate_iters)

    detected_bit = exposure.mask.getPlaneBitMask("DETECTED")
    exposure.mask.array[src_mask] |= detected_bit

    n_masked = int(src_mask.sum())
    frac = n_masked / src_mask.size
    if verbose:
        print(f"  masked {n_masked} px ({100 * frac:.1f}%) "
              f"as source for background refit")
    if frac > 0.5:
        print(f"  WARNING: {100 * frac:.1f}% of the image is masked as source -- "
              f"background refit and PSF redetermination will be starved of "
              f"sky pixels/star candidates. Consider raising n_sigma or "
              f"lowering dilate_iters.")
    return src_mask


def refit_background(exposure, algorithm='AKIMA_SPLINE', bin_size=32,
                      use_approx=True,
                      stats_mask=("BAD", "SAT", "EDGE", "NO_DATA", "DETECTED",
                                  "CR", "INTRP")):
    """Re-derive the background ignoring masked pixels and subtract it in place.

    Pass a fresh clone of the raw (never background-subtracted) exposure, already flagged by
    `mask_sources_from_model`, or the background is subtracted twice.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure to background-subtract in place.
    algorithm : `str`, optional
        Interpolation algorithm (only used when ``use_approx`` is False).
    bin_size : `int`, optional
        Background bin size [px].
    use_approx : `bool`, optional
        Fit a smooth Chebyshev surface to the binned statistics instead of interpolating;
        more robust when many pixels are masked.
    stats_mask : `tuple` [`str`], optional
        Mask planes ignored when computing the background statistics.

    Returns
    -------
    exposure : `lsst.afw.image.Exposure`
        The background-subtracted exposure (same object).
    background : `lsst.afw.image.Image`
        The subtracted background image.
    """
    bg_config = SubtractBackgroundConfig()
    bg_config.algorithm = algorithm
    bg_config.useApprox = use_approx
    bg_config.binSize = bin_size
    bg_config.ignoredPixelMask = list(stats_mask)

    bgTask = SubtractBackgroundTask(config=bg_config)
    bgResult = bgTask.run(exposure)
    return exposure, bgResult.background.getImage()


# =====================================================================
# Step B0: turn raw (x, y) peak positions into a measured SourceCatalog,
# so a previously-converged peak list (or an externally-supplied set of
# "known good" star positions) can feed the star selector directly.
# =====================================================================

def catalog_from_peaks(exposure, peak_positions, schema, measureTask,
                        stamp_radius=15, exclude_positions=None,
                        exclude_radius=5.0, verbose=True):
    """Build a measured SourceCatalog from raw (x, y) peak positions.

    Each peak gets a circular footprint and one Peak record, and the catalogue is measured with
    the same ``measureTask`` as fresh detections, so the fields are comparable and usable by the
    star selector.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure to measure on.
    peak_positions : `list` [`tuple` [`float`, `float`]]
        Peak positions in parent pixel coordinates.
    schema : `lsst.afw.table.Schema`
        The schema used to build ``measureTask``.
    measureTask : `lsst.meas.base.SingleFrameMeasurementTask`
        Measurement task.
    stamp_radius : `int`, optional
        Footprint radius [px].
    exclude_positions : `list` [`tuple` [`float`, `float`]], optional
        Skip peaks within ``exclude_radius`` of these positions (e.g. existing detections).
    exclude_radius : `float`, optional
        See ``exclude_positions`` [px].
    verbose : `bool`, optional
        Print how many peaks were skipped.

    Returns
    -------
    catalog : `lsst.afw.table.SourceCatalog`
        The measured catalogue.
    """
    table = afwTable.SourceTable.make(schema)
    catalog = afwTable.SourceCatalog(table)
    exp_bbox = exposure.getBBox()

    exclude_tree = None
    if exclude_positions is not None and len(exclude_positions) > 0:
        exclude_tree = cKDTree(np.asarray(exclude_positions))

    n_skipped_dup = 0
    n_skipped_bounds = 0
    for (x, y) in peak_positions:
        if exclude_tree is not None:
            dist, _ = exclude_tree.query([x, y])
            if dist < exclude_radius:
                n_skipped_dup += 1
                continue

        center = geom.Point2I(int(round(x)), int(round(y)))
        if not exp_bbox.contains(center):
            n_skipped_bounds += 1
            continue

        # NOTE: SpanSet/Stencil location can differ between stack versions --
        # see module docstring if this errors on your stack.
        spans = afwGeom.SpanSet.fromShape(stamp_radius, afwGeom.Stencil.CIRCLE, offset=center)
        spans = spans.clippedTo(exp_bbox)
        if spans.getArea() == 0:
            n_skipped_bounds += 1
            continue

        footprint = afwDetect.Footprint(spans, exp_bbox)
        footprint.addPeak(x, y, 0.0)
        record = catalog.addNew()
        record.setFootprint(footprint)

    if len(catalog) > 0:
        measureTask.run(measCat=catalog, exposure=exposure)
    catalog = catalog.copy(deep=True)

    if verbose:
        print(f"  catalog_from_peaks: measured {len(catalog)} peak(s) "
              f"({n_skipped_dup} skipped as duplicates of existing sources, "
              f"{n_skipped_bounds} skipped out of bounds)")
    return catalog


import math

import math

# =====================================================================
# Step B-helper 1: version-safe PSF-determiner kernel-size setter
# =====================================================================

def _set_psf_determiner_kernel_size(psf_determiner_config, kernel_size, verbose=True):
    """Set the size fields that exist on this stack's PcaPsfDeterminerConfig.

    ``stampSize`` always; ``kernelSizeMin`` / ``kernelSizeMax`` only where present (they are
    absent on some stack versions).

    Parameters
    ----------
    psf_determiner_config : `lsst.meas.algorithms.PcaPsfDeterminerConfig`
        Config to modify in place.
    kernel_size : `int`
        Kernel size [px].
    verbose : `bool`, optional
        Print the extra fields that were set.

    Returns
    -------
    psf_determiner_config : `lsst.meas.algorithms.PcaPsfDeterminerConfig`
        The same config.
    """
    psf_determiner_config.stampSize = kernel_size

    extra_fields = ["kernelSizeMin", "kernelSizeMax"]
    set_fields = []
    for field in extra_fields:
        if hasattr(psf_determiner_config, field):
            setattr(psf_determiner_config, field, kernel_size)
            set_fields.append(field)

    if verbose and set_fields:
        print(f"  PSF determiner: also set {set_fields} = {kernel_size}")
    return psf_determiner_config


# =====================================================================
# Step B-helper 2: drop PSF candidates whose stamp extraction silently failed
# =====================================================================

def _filter_valid_psf_candidates(psfCandidateList, verbose=True):
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
    n_bad = 0
    for cand in psfCandidateList:
        try:
            stamp = cand.getMaskedImage().getImage().array
        except Exception as e:
            n_bad += 1
            if verbose:
                print(f"  dropping candidate id={cand.getSource().getId()}: "
                      f"stamp re-fetch failed ({e})")
            continue
        if not np.isfinite(stamp).all() or stamp.sum() <= 0:
            n_bad += 1
            continue
        good.append(cand)

    if verbose and n_bad > 0:
        print(f"  filtered out {n_bad} PSF candidate(s) with unusable stamps "
              f"-> {len(good)} of {len(psfCandidateList)} remain")
    return good


# =====================================================================
# Step B: redetect + reselect stars + redetermine PSF on the freshly
# background-subtracted exposure.
# =====================================================================

def redetect_and_refit_psf(exposure, det_threshold=3.0, det_min_pixels=2,
                            det_local_bg_bin=16, det_bg_bin=32,
                            star_selector_config=None,
                            psf_determiner_config=None,
                            kernel_size=41,
                            extra_peak_positions=None,
                            force_include_star_positions=None,
                            peak_stamp_radius=15,
                            verbose=True):
    """Detect, deblend, measure, select stars and determine a new PSF.

    Mirrors the original one-shot pipeline, so each outer pass works on the newly cleaned
    background.  Extra peak positions (e.g. the previous pass's converged catalogue) are added
    before star selection; forced star positions are appended after selection, bypassing its
    cuts.  Candidates with failed stamps are filtered out before the PSF determiner runs (even
    one can make it reject the whole set).  If too few candidates survive, the exposure's
    current PSF is kept and a "NOT updated" note is printed.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Background-subtracted exposure.
    det_threshold : `float`, optional
        Detection threshold [sigma].
    det_min_pixels : `int`, optional
        Minimum footprint size.
    det_local_bg_bin : `int`, optional
        Bin size of the temporary local background [px].
    det_bg_bin : `int`, optional
        Bin size of the detection background [px].
    star_selector_config : `lsst.meas.algorithms.ObjectSizeStarSelectorConfig`, optional
        Star selector configuration.
    psf_determiner_config : `lsst.meas.algorithms.PcaPsfDeterminerConfig`, optional
        PSF determiner configuration.
    kernel_size : `int`, optional
        PSF kernel size [px], used for both the candidates and the determiner.
    extra_peak_positions : `list` [`tuple` [`float`, `float`]], optional
        Additional peaks measured and offered to the star selector.
    force_include_star_positions : `list` [`tuple` [`float`, `float`]], optional
        Star positions appended to the selected stars regardless of the selector's cuts.
    peak_stamp_radius : `int`, optional
        Footprint radius of the extra / forced peaks [px].
    verbose : `bool`, optional
        Print progress.

    Returns
    -------
    psf_new : `lsst.afw.detection.Psf`
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

    det_config = SourceDetectionConfig()
    det_config.thresholdValue = det_threshold
    det_config.thresholdType = "stdev"
    det_config.minPixels = det_min_pixels
    det_config.doTempLocalBackground = True
    det_config.tempLocalBackground.binSize = det_local_bg_bin
    det_config.background.binSize = det_bg_bin

    detectionTask = SourceDetectionTask(schema=schema, config=det_config)
    deblendTask = SourceDeblendTask(schema=schema)
    measureTask = SingleFrameMeasurementTask(schema=schema)

    table = afwTable.SourceTable.make(schema)
    detRes = detectionTask.run(table, exposure)
    sources = detRes.sources

    deblendTask.run(exposure, sources)
    measureTask.run(measCat=sources, exposure=exposure)
    sources = sources.copy(deep=True)

    # ---- fold in the full previously-converged peak catalog, if given ----
    if extra_peak_positions:
        fresh_xy = (np.array([(r.getX(), r.getY()) for r in sources])
                    if len(sources) else None)
        extra_cat = catalog_from_peaks(
            exposure, extra_peak_positions, schema, measureTask,
            stamp_radius=peak_stamp_radius, exclude_positions=fresh_xy,
            exclude_radius=peak_stamp_radius, verbose=verbose)
        if len(extra_cat) > 0:
            combined = afwTable.SourceCatalog(sources.table)
            combined.extend(sources, deep=True)
            combined.extend(extra_cat, deep=True)
            # extend() does not guarantee contiguous storage, and
            # ObjectSizeStarSelectorTask (like several afw table consumers)
            # requires a contiguous catalog -- without this copy it raises
            # "Input catalogs for source selection must be contiguous."
            sources = combined.copy(deep=True)
            if verbose:
                print(f"  combined catalog for star selection: {len(sources)} sources "
                      f"({len(extra_cat)} pulled in from the prior converged peak catalog)")

    starSelector = (ObjectSizeStarSelectorTask(config=star_selector_config)
                     if star_selector_config else ObjectSizeStarSelectorTask())
    selection = starSelector.run(sourceCat=sources, exposure=exposure)
    starCat = selection.sourceCat

    # ---- force the ORIGINAL star candidates back in, bypassing the selector ----
    if force_include_star_positions:
        selected_xy = (np.array([(r.getX(), r.getY()) for r in starCat])
                       if len(starCat) else None)
        forced_cat = catalog_from_peaks(
            exposure, force_include_star_positions, schema, measureTask,
            stamp_radius=peak_stamp_radius, exclude_positions=selected_xy,
            exclude_radius=peak_stamp_radius, verbose=verbose)
        if len(forced_cat) > 0:
            combined_stars = afwTable.SourceCatalog(starCat.table)
            combined_stars.extend(starCat, deep=True)
            combined_stars.extend(forced_cat, deep=True)
            n_before = len(starCat)
            # same contiguity requirement as above -- MakePsfCandidatesTask
            # also expects/iterates a contiguous catalog.
            starCat = combined_stars.copy(deep=True)
            if verbose:
                print(f"  forced back {len(forced_cat)} original PSF-star "
                      f"candidate(s) -> {n_before} -> {len(starCat)} total candidates")

    if verbose:
        print(f"  {len(starCat)} PSF-star candidates out of {len(sources)} detections")

    makeCandidatesConfig = MakePsfCandidatesTask.ConfigClass()
    makeCandidatesConfig.kernelSize = kernel_size
    makeCandidatesTask = MakePsfCandidatesTask(config=makeCandidatesConfig)
    candidatesResult = makeCandidatesTask.run(starCat, exposure=exposure)
    psfCandidateList = candidatesResult.psfCandidates

    # drop candidates whose stamp extraction silently failed -- this is
    # the fix for the "No viable PSF candidates survive" failures
    psfCandidateList = _filter_valid_psf_candidates(psfCandidateList, verbose=verbose)

    if psf_determiner_config is None:
        psf_determiner_config = PcaPsfDeterminerTask.ConfigClass()
        psf_determiner_config.spatialOrder = 1
        psf_determiner_config.sizeCellX = 512
        psf_determiner_config.sizeCellY = 512
        psf_determiner_config.maxCandidates = 1500
        psf_determiner_config.doRejectBlends = False # changed from True
        psf_determiner_config.nStarPerCell = 10

    _set_psf_determiner_kernel_size(psf_determiner_config, kernel_size, verbose=verbose)

    psf_prev = exposure.getPsf()
    point = exposure.getBBox().getCenter()

    def _fwhm(psf_obj):
        shape = psf_obj.computeShape(point)
        return 2.355 * np.sqrt(0.5 * (shape.getIxx() + shape.getIyy()))

    fwhm_prev = _fwhm(psf_prev)

    psfDeterminer = PcaPsfDeterminerTask(config=psf_determiner_config)
    try:
        if len(psfCandidateList) == 0:
            raise IndexError("No viable PSF candidates survive (all candidates "
                              "filtered out before determinePsf was even called)")
        psf_new, cellSet_new = psfDeterminer.determinePsf(
            exposure=exposure, psfCandidateList=psfCandidateList)
        exposure.setPsf(psf_new)
        if verbose:
            fwhm_new = _fwhm(psf_new)
            print(f"  PSF updated: FWHM {fwhm_prev:.3f} -> {fwhm_new:.3f} px (crop center)")
    except IndexError as e:
        if verbose:
            print(f"  WARNING: PSF redetermination failed ({e}) with "
                  f"{len(psfCandidateList)} valid candidates from {len(starCat)} "
                  f"star-selected sources -- PSF NOT updated this pass "
                  f"(kept previous, FWHM {fwhm_prev:.3f} px)")
        psf_new = psf_prev
        cellSet_new = None

    return psf_new, sources, starCat, cellSet_new


# =====================================================================
# Step C: the outer loop itself
# =====================================================================

def outer_loop_psf_background_refinement(
    raw_exposure_crop,
    psf_initial,
    bad_bitmask,
    initial_star_cat=None,
    n_outer_iterations=3,
    n_inner_iterations=10,
    inner_kwargs=None,
    outer_convergence_frac=0.02,
    bg_use_approx=True,
    psf_kernel_size=41,
    verbose=True,
):
    """Alternate inner-loop fits with background and PSF re-derivation.

    Pass 0 runs `iterative_joint_detection` with the current background and PSF.  Pass k > 0
    masks every known source (from pass k-1's model) on a fresh clone of the raw image,
    re-derives the background, re-determines the PSF (offering pass k-1's catalogue and forcing
    in ``initial_star_cat``) and reruns the inner loop.

    Parameters
    ----------
    raw_exposure_crop : `lsst.afw.image.Exposure`
        Crop before any background subtraction.
    psf_initial : `lsst.afw.detection.Psf`
        PSF of the first pass.
    bad_bitmask : `int`
        Mask bits of pixels to ignore.
    initial_star_cat : `lsst.afw.table.SourceCatalog`, optional
        Star selector output whose positions are forced into every PSF determination.
    n_outer_iterations : `int`, optional
        Maximum number of outer passes.
    n_inner_iterations : `int`, optional
        Maximum number of inner passes (when ``inner_kwargs`` is not given).
    inner_kwargs : `dict`, optional
        Keyword arguments of `iterative_joint_detection`.
    outer_convergence_frac : `float`, optional
        Stop when the residual standard deviation improves by less than this fraction.
    bg_use_approx : `bool`, optional
        See `refit_background`.
    psf_kernel_size : `int`, optional
        PSF kernel size [px].
    verbose : `bool`, optional
        Print progress.

    Returns
    -------
    outer_history : `list` [`dict`]
        One dict per outer pass: ``psf``, ``exposure``, ``result`` (inner-loop output),
        ``resid_std``, ``n_peaks``.
    """
    if inner_kwargs is None:
        inner_kwargs = dict(n_iterations=n_inner_iterations, threshold_start=8.0,
                             threshold_step=-1.5, threshold_min=1.0,
                             min_pixels=1, min_new_peaks=1, convergence_frac=0.02)

    original_star_positions = None
    if initial_star_cat is not None and len(initial_star_cat) > 0:
        original_star_positions = [
            (r.getCentroid().getX(), r.getCentroid().getY()) for r in initial_star_cat
        ]
        if verbose:
            print(f"  {len(original_star_positions)} original PSF-star position(s) "
                  f"will be forced into every subsequent outer pass")

    raw = raw_exposure_crop.clone()  # kept pristine for every re-background pass
    exposure_working = raw_exposure_crop.clone()
    exposure_working.setPsf(psf_initial)
    psf_working = psf_initial
    bg_image_working = None  # no separately-tracked bg image for outer pass 0

    outer_history = []

    for outer_it in range(n_outer_iterations):
        if verbose:
            print(f"\n===== Outer pass {outer_it + 1}/{n_outer_iterations} =====")

        result = iterative_joint_detection(
            exposure_working, psf_working, bad_bitmask, **inner_kwargs)
        final_resid_std = result['resid_std_history'][-1]

        outer_history.append({
            'psf': psf_working,
            'exposure': exposure_working,
            'result': result,
            'resid_std': final_resid_std,
            'n_peaks': len(result['peaks']),
            'bg_image': bg_image_working,
        })

        if verbose:
            print(f"  outer pass {outer_it + 1} converged: "
                  f"{len(result['peaks'])} peaks, residual std={final_resid_std:.4f}")

        if outer_it > 0:
            prev_std = outer_history[-2]['resid_std']
            improvement = (abs(prev_std - final_resid_std) / prev_std
                            if prev_std > 0 else 0)
            if improvement < outer_convergence_frac:
                if verbose:
                    print(f"  outer residual std improved by only "
                          f"{improvement * 100:.2f}% -- outer loop converged, stopping")
                break

        if outer_it == n_outer_iterations - 1:
            break  # no need to re-derive background/PSF after the final pass

        # ---- re-derive background + PSF from this pass's converged model ----
        exposure_next = raw.clone()
        mask_sources_from_model(exposure_next, result['model_full'], verbose=verbose)
        exposure_next, bgImage_next = refit_background(exposure_next, use_approx=bg_use_approx)

        psf_new, sources_new, starCat_new, cellSet_new = redetect_and_refit_psf(
            exposure_next, kernel_size=psf_kernel_size,
            extra_peak_positions=result['peaks'],
            force_include_star_positions=original_star_positions,
            verbose=verbose)

        exposure_working = exposure_next
        psf_working = psf_new
        bg_image_working = bgImage_next

    return outer_history


# =====================================================================
# Diagnostics
# =====================================================================

def _extract_peak_xy(peaks):
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


def plot_psf_models_by_outer_pass(outer_history, exposure_center=None):
    """Plot the PSF kernel image used in each outer pass.

    Parameters
    ----------
    outer_history : `list` [`dict`]
        Outer-loop history (each entry has ``psf`` and ``exposure``).
    exposure_center : `lsst.geom.Point2D`, optional
        Position at which to render the kernel (default: the exposure centre).
    """
    n = len(outer_history)
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 4))
    if n == 1:
        axes = [axes]

    for i, (ax, h) in enumerate(zip(axes, outer_history)):
        point = exposure_center or h["exposure"].getBBox().getCenter()
        psf_img = h["psf"].computeKernelImage(geom.Point2D(point)).array
        im = ax.imshow(psf_img, origin="lower", cmap="inferno")
        ax.set_title(f"Outer pass {i + 1} PSF")
        ax.set_xticks([])
        ax.set_yticks([])
        plt.colorbar(im, ax=ax, fraction=0.046)

    plt.tight_layout()
    plt.show()


def plot_background_by_outer_pass(outer_history):
    """Plot the background re-derived going into each outer pass.

    The first pass uses the background subtracted before the pipeline ran and is skipped.

    Parameters
    ----------
    outer_history : `list` [`dict`]
        Outer-loop history (entries may have ``bg_image``).
    """
    passes_with_bg = [(i, h) for i, h in enumerate(outer_history)
                       if h.get("bg_image") is not None]
    if not passes_with_bg:
        print("No re-derived background images to show yet -- these are only "
              "captured starting from outer pass 2 onward.")
        return

    fig, axes = plt.subplots(1, len(passes_with_bg), figsize=(4 * len(passes_with_bg), 4))
    if len(passes_with_bg) == 1:
        axes = [axes]

    for ax, (i, h) in zip(axes, passes_with_bg):
        im = ax.imshow(h["bg_image"].array, origin="lower", cmap="viridis")
        ax.set_title(f"Background entering pass {i + 1}")
        ax.set_xticks([])
        ax.set_yticks([])
        plt.colorbar(im, ax=ax, fraction=0.046)

    plt.tight_layout()
    plt.show()


def plot_sources_overlaid_by_pass(outer_history, marker_size=25, vmin_pct=5, vmax_pct=99):
    """Plot the data of each outer pass with its converged peak positions.

    Parameters
    ----------
    outer_history : `list` [`dict`]
        Outer-loop history.
    marker_size : `float`, optional
        Marker size.
    vmin_pct, vmax_pct : `float`, optional
        Percentiles of the display range.
    """
    n = len(outer_history)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 5))
    if n == 1:
        axes = [axes]

    for i, (ax, h) in enumerate(zip(axes, outer_history)):
        exposure = h["exposure"]
        data = exposure.image.array
        vmin, vmax = np.nanpercentile(data, [vmin_pct, vmax_pct])
        ax.imshow(data, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)

        # FIX: `_extract_peak_xy` returns positions in the exposure's
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

        xs, ys = _extract_peak_xy(h["result"]["peaks"])
        ax.scatter(xs - x0, ys - y0, s=marker_size, facecolors="none",
                   edgecolors="#00FF00", linewidths=0.8)
        ax.set_xlim(0, data.shape[1])
        ax.set_ylim(0, data.shape[0])
        ax.set_title(f"Outer pass {i + 1}: {len(xs)} sources")
        ax.set_xticks([])
        ax.set_yticks([])

    plt.tight_layout()
    plt.show()

def plot_outer_loop_summary(outer_history):
    """Plot the residual standard deviation and the source count against outer pass.

    Parameters
    ----------
    outer_history : `list` [`dict`]
        Outer-loop history (entries have ``resid_std`` and ``n_peaks``).
    """
    n_outer = len(outer_history)
    passes = np.arange(1, n_outer + 1)
    resid_stds = [h['resid_std'] for h in outer_history]
    n_peaks = [h['n_peaks'] for h in outer_history]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    fig.patch.set_facecolor('white')

    axes[0].plot(passes, resid_stds, 'o-', color='#DD8452')
    axes[0].set_xlabel('outer pass')
    axes[0].set_ylabel('final residual std [nJy]')
    axes[0].set_title('Residual std vs outer pass')
    axes[0].set_xticks(passes)

    axes[1].plot(passes, n_peaks, 'o-', color='#4C72B0')
    axes[1].set_xlabel('outer pass')
    axes[1].set_ylabel('total converged peaks')
    axes[1].set_title('Peak count vs outer pass')
    axes[1].set_xticks(passes)

    plt.tight_layout()
    plt.show()


def plot_psf_fwhm_by_outer_pass(outer_history, exposure_center=None):
    """Plot the FWHM of the PSF used in each outer pass.

    Parameters
    ----------
    outer_history : `list` [`dict`]
        Outer-loop history.
    exposure_center : `lsst.geom.Point2D`, optional
        Position at which to evaluate the PSF (default: the exposure centre).

    Returns
    -------
    fwhms : `list` [`float`]
        FWHM [px] per outer pass.
    """
    fwhms = []
    for h in outer_history:
        point = exposure_center or h['exposure'].getBBox().getCenter()
        shape = h['psf'].computeShape(geom.Point2D(point))
        fwhm = 2.355 * np.sqrt(0.5 * (shape.getIxx() + shape.getIyy()))
        fwhms.append(fwhm)

    passes = np.arange(1, len(outer_history) + 1)
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(passes, fwhms, 'o-', color='#55A868')
    ax.set_xlabel('outer pass')
    ax.set_ylabel('PSF FWHM [px] (crop center)')
    ax.set_title('PSF FWHM vs outer pass')
    ax.set_xticks(passes)
    plt.tight_layout()
    plt.show()
    return fwhms


# ######################################################################
# [from notebook cell: cutout analysis]
# ######################################################################

"""
visualize_cutout_fit_v2.py
==========================
Iterative joint PSF fit + outer-loop PSF/background refinement on a cutout, with
VISUALIZATION and TIMING of every step. This is visualize_cutout_fit.py merged with the
new PSF-SHAPE CONSTRAINT (formerly shape_constraint.py).

Needs in the session (from your other scripts):
    exposure_crop, psf, bad_bitmask, starCat, bgImage (optional)
    collapse_saturated_peaks, merge_close_peaks, bbox_to_slices,
    mask_sources_from_model, refit_background, redetect_and_refit_psf

IMPORTANT: do NOT define the older compute_psf_derivatives / refine_peak_positions
(forward-difference, eps=0.5) in your session AFTER loading this file -- they would
silently override the fixed versions here.

NEW in this version: PSF-shape constraint (eliminate_nonpsf)
------------------------------------------------------------
Per source, on NEIGHBOUR-SUBTRACTED data  D_j = data - (model - own_j):
  1. local-plane S/N : fit D_j ~ a*P + c + gx*dx + gy*dy over a ~4 sigma disc. The
     constant + tilt absorb halo pedestals / gradients near bright stars, so a source
     must be a compact PSF-shaped bump on top of the local background.
  2. reduced chi2 of that fit (noise-scaled): rings / arcs / dipoles -> chi2 >> 1.
  3. width test: refit with Gaussians of sigma = s * psf_sigma. If the data
     significantly prefer a width outside width_range -> too broad (background blob)
     or too narrow (spike). Faint untrusted sources (local S/N < weak_snr) use the
     looser threshold dchi2_width_weak.
Sources with local S/N >= trust_snr (incl. saturated ones) skip tests 2-3. Failures are
removed by backward elimination (worst ~25% per round, re-solve, re-test).
Run order in the inner loop: ... eliminate_weak -> eliminate_nonpsf.
New knobs: use_shape_cut, shape_kwargs, shape_in_first_pass (run_cutout_pipeline).

NEW: S/N calibrated against the MATCHED-FILTER noise
----------------------------------------------------
empirical_noise_scale() now returns max(pixel-MAD scale, matched-filter scale). The
matched-filter scale is the robust std of the PSF-matched-filter S/N map of the residual
on source-free sky, so it includes pixel-to-pixel noise correlation (warped coadds,
covariance the variance plane ignores). Pixel MAD alone cannot see that, which inflates
every S/N and lets faint noise bumps through. chi2 stays on the pixel scale
(empirical_noise_scale(..., matched=False)). The scales are printed by eliminate_weak
and eliminate_nonpsf: if the S/N scale is clearly above the pixel scale, correlated
noise was the cause.

Earlier rounds (kept): empirical_noise_scale, amp_snr_full, eliminate_weak,
edge_margin_fwhm, negative_source_check, source_snr_histogram, faint-companion fixes,
stack-based PSF refit with polishing cycles, effective_variance (sys_frac),
prune_and_merge, reject_near_bright, centroid freezing near brighter neighbours,
Gauss-Newton centroiding with central-difference PSF derivatives, tapered PSF stamps.

NOTE on "groups": the actual solve is ONE joint bounded least-squares over all sources.
"Groups" are a diagnostic (overlap of PSF stamps above `support_frac` of peak).
"""

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
        self.t_start = time.perf_counter()

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
        self.t_start = time.perf_counter()


TIMER = StepTimer()


def plot_timing_summary(timer=TIMER):
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
    print(f"wall clock    : {time.perf_counter() - timer.t_start:.2f} s")

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(17, 6),
                                   gridspec_kw={'width_ratios': [1, 1.5]})
    labs = sorted(totals, key=lambda k: totals[k])
    ax0.barh(labs, [totals[k] for k in labs], color='#4C72B0')
    ax0.set_xlabel('total seconds')
    ax0.set_title('Total compute time by step')
    ax0.grid(axis='x', alpha=0.3)

    inner = [r for r in comp if r['inner'] is not None and r['outer'] is not None]
    keys = sorted({(r['outer'], r['inner']) for r in inner})
    inner_labels = sorted({r['label'] for r in inner})
    cmap = plt.cm.tab20(np.linspace(0, 1, max(len(inner_labels), 1)))
    bottoms = np.zeros(len(keys))
    for lab, col in zip(inner_labels, cmap):
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

def make_cutout(exposure_crop, x_range=(0, 100), y_range=(0, 100)):
    """Cut a sub-exposure out of a crop, in coordinates relative to the crop's corner.

    Parameters
    ----------
    exposure_crop : `lsst.afw.image.Exposure`
        The crop.
    x_range, y_range : `tuple` [`int`, `int`], optional
        Pixel ranges relative to the crop's origin (end exclusive).

    Returns
    -------
    cutout : `lsst.afw.image.Exposure`
        A deep copy of the sub-exposure.

    Raises
    ------
    TypeError
        Raised if ``exposure_crop`` is not an Exposure.
    """
    if not hasattr(exposure_crop, 'getXY0'):
        raise TypeError(
            "make_cutout() got a "
            f"{type(exposure_crop).__name__}, not an lsst.afw.image Exposure. "
            "The `exposure_crop` variable was probably reassigned after "
            "`exposure_crop = exposure[...].clone()`. Re-run the cell that builds it "
            "(with .clone() and .setPsf(psf)) before calling run_cutout_pipeline.")
    xy0 = exposure_crop.getXY0()
    bbox = geom.Box2I(
        geom.Point2I(xy0.getX() + x_range[0], xy0.getY() + y_range[0]),
        geom.Extent2I(x_range[1] - x_range[0], y_range[1] - y_range[0]))
    return exposure_crop[bbox].clone()


def _extent(bbox):
    """Return the matplotlib extent (pixel edges) of a box."""
    return (bbox.getMinX() - 0.5, bbox.getMaxX() + 0.5,
            bbox.getMinY() - 0.5, bbox.getMaxY() + 0.5)


def _coords(bbox):
    """Return the pixel x and y coordinates covered by a box."""
    xs = bbox.getMinX() + np.arange(bbox.getWidth())
    ys = bbox.getMinY() + np.arange(bbox.getHeight())
    return xs, ys


def _data_norm(arr):
    """Return a log-stretch display normalisation for an image."""
    return simple_norm(arr, vmin=np.nanmin(arr), vmax=np.nanpercentile(arr, 99.5),
                       stretch='log')


def _overlay_bad(ax, bad, extent):
    """Overlay bad pixels in translucent red."""
    if bad is None or not bad.any():
        return
    ax.imshow(np.ma.masked_where(~bad, np.ones(bad.shape)), extent=extent,
              origin='lower', cmap=ListedColormap(['red']), alpha=0.35,
              interpolation='nearest')


def _clean_axes(ax):
    """Turn the axis grid off."""
    ax.grid(False)


def _peaks_array(peaks):
    """Return peaks as an (N, 2) float array."""
    return np.asarray(peaks, dtype=float).reshape(-1, 2)


def _fwhm(psf_obj, point):
    """Return the FWHM [px] of a PSF at a point, from its second moments."""
    shape = psf_obj.computeShape(point)
    return 2.355 * np.sqrt(0.5 * (shape.getIxx() + shape.getIyy()))


def _footprints_to_mask(footprints, exposure):
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


def _match_prev_amps(peaks, prev_peaks, prev_amps, tol=1.5):
    """Return each peak's amplitude from the previous iteration (NaN if it is new)."""
    out = np.full(len(peaks), np.nan)
    if len(peaks) == 0 or len(prev_peaks) == 0:
        return out
    d, i = cKDTree(prev_peaks).query(peaks)
    ok = d < tol
    out[ok] = np.asarray(prev_amps)[i[ok]]
    return out


# ---- PSF stamp with tamed wings / peak filtering -----------------------

# (r_start, r_end) as fractions of the stamp half-width. Inside r_start the PSF is
# untouched; between them it is cosine-tapered to exactly 0. None = raw kernel.
WING_TAPER = (0.8, 1.0)


def psf_stamp(psf, x, y, clip_bbox, taper=None):
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
    clip_bbox : `lsst.geom.Box2I`
        Box to clip the stamp to.
    taper : `tuple` [`float`, `float`], optional
        Inner and outer taper radii as fractions of the stamp half-size.

    Returns
    -------
    stamp : `numpy.ndarray` or `None`
        The clipped stamp, or None if it misses ``clip_bbox``.
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
    bb.clip(clip_bbox)
    if bb.isEmpty():
        return None, bb
    sl = (slice(bb.getMinY() - full.getMinY(), bb.getMaxY() - full.getMinY() + 1),
          slice(bb.getMinX() - full.getMinX(), bb.getMaxX() - full.getMinX() + 1))
    return arr[sl], bb


def drop_masked_peaks(peaks, exposure, planes=('EDGE', 'NO_DATA')):
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
    n_dropped : `int`
        Number of peaks dropped.
    """
    mk = exposure.mask
    planes_present = mk.getMaskPlaneDict()
    bits = 0
    for p in planes:
        if p in planes_present:
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


def drop_out_of_bounds(peaks, exposure, margin=0.0):
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
    n_dropped : `int`
        Number of peaks dropped.
    """
    bb = exposure.getBBox()
    x0, x1 = bb.getMinX() - margin, bb.getMaxX() + margin
    y0, y1 = bb.getMinY() - margin, bb.getMaxY() + margin
    keep = [(x, y) for x, y in peaks if x0 <= x <= x1 and y0 <= y <= y1]
    return keep, len(peaks) - len(keep)


# ---- over-counting guards ---------------------------------------------

def effective_variance(exposure, model=None, frac=0.03, poisson_k=0.0):
    """Return the variance used to weight the fit.

    ``variance + (frac * model)**2 + poisson_k * model``: the pixel variance plus a fractional
    model-error term and missing photon noise.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure providing the variance plane.
    model : `numpy.ndarray`, optional
        Current model image (no extra terms if None).
    frac : `float`, optional
        Fractional model error.
    poisson_k : `float`, optional
        Extra variance per unit model flux.

    Returns
    -------
    var : `numpy.ndarray`
        The effective variance.
    """
    var = exposure.variance.array.copy()
    if model is not None:
        m = np.clip(model, 0, None)
        var += (frac * m) ** 2 + poisson_k * m
    return var


def reject_near_bright(new_peaks, old_peaks, old_is_bright, min_sep):
    """Drop new peaks within ``min_sep`` px of an existing bright source.

    Parameters
    ----------
    new_peaks : `list` [`tuple` [`float`, `float`]]
        Candidate peaks.
    old_peaks : array-like
        Existing source positions, shape (N, 2).
    old_is_bright : array-like [`bool`]
        Which existing sources are bright.
    min_sep : `float`
        Rejection radius [px].

    Returns
    -------
    keep : `list` [`tuple` [`float`, `float`]]
        Remaining candidates.
    n_rejected : `int`
        Number of candidates rejected.
    """
    old = np.asarray(old_peaks, float).reshape(-1, 2)[np.asarray(old_is_bright, bool)]
    if len(old) == 0 or len(new_peaks) == 0:
        return list(new_peaks), 0
    d, _ = cKDTree(old).query(np.asarray(new_peaks, float))
    keep = [p for p, di in zip(new_peaks, d) if di > min_sep]
    return keep, len(new_peaks) - len(keep)


def prune_and_merge(fit, peaks, min_snr=5.0, min_sep=3.0,
                    bright_flags=None, near_bright_sep=None, near_bright_snr=8.0):
    """Drop low-S/N sources, then merge pairs closer than ``min_sep`` px.

    The brighter source of a pair survives at the amplitude-weighted position.  Sources within
    ``near_bright_sep`` px of a bright source must reach ``near_bright_snr`` instead of
    ``min_snr``.  Uses the cheap diagonal-only S/N; the honest S/N is applied by
    `eliminate_weak`.

    Parameters
    ----------
    fit : `dict`
        Joint fit (from `joint_psf_fit_diag`) of ``peaks``.
    peaks : array-like
        Source positions, shape (N, 2).
    min_snr : `float`, optional
        Minimum S/N.
    min_sep : `float`, optional
        Merge radius [px].
    bright_flags : array-like [`bool`], optional
        Which sources are bright.
    near_bright_sep : `float`, optional
        Radius around bright sources with the stricter S/N cut [px].
    near_bright_snr : `float`, optional
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

    thresh = np.full(len(peaks), float(min_snr))
    n_strict = 0
    if near_bright_sep and bright_flags is not None:
        bf = np.asarray(bright_flags, bool)
        if bf.any() and len(bf) == len(peaks):
            d, _ = cKDTree(peaks[bf]).query(peaks)
            near = (d > 0.5) & (d < near_bright_sep)
            thresh[near] = np.maximum(thresh[near], near_bright_snr)
            n_strict = int(near.sum())
    keep = snr >= thresh
    alive = keep.copy()
    for a in np.argsort(-amps):
        if not alive[a]:
            continue
        d = np.hypot(*(peaks - peaks[a]).T)
        for b in np.flatnonzero((d < min_sep) & alive):
            if b == a:
                continue
            w = amps[a] + amps[b]
            peaks[a] = (amps[a] * peaks[a] + amps[b] * peaks[b]) / w
            alive[b] = False
    print(f"  prune/merge: {len(peaks)} -> {alive.sum()} "
          f"({(~keep).sum()} low S/N, {(keep & ~alive).sum()} merged; "
          f"{n_strict} held to S/N>={near_bright_snr:g} near bright stars)")
    return [tuple(p) for p in peaks[alive]]


# ---- honest S/N + backward elimination ---------------------------------

def _pixel_noise_scale(fit, exposure, model_frac=0.01):
    """Return the robust standard deviation of residual / sigma on sky pixels (at least 1)."""
    resid = exposure.image.array - fit['model_full']
    z = resid / np.sqrt(fit['variance'])
    mmax = fit['model_full'].max()
    sky = fit['good'] & (fit['model_full'] < model_frac * mmax)
    zz = z[sky]
    if zz.size < 50:
        zz = z[fit['good']]
    return max(1.0, 1.4826 * np.median(np.abs(zz - np.median(zz))))


def _matched_filter_noise_scale(fit, exposure, core_half=10, floor_sigma=0.25):
    """Return the robust standard deviation of the residual's matched-filter S/N on sky pixels.

    Unlike the pixel scale this includes pixel-to-pixel noise correlation that the variance
    plane ignores.  Returns 1.0 if there are too few clean pixels or no PSF.
    """
    psf = exposure.getPsf()
    if psf is None:
        return 1.0
    k = psf.computeKernelImage(exposure.getBBox().getCenter()).array
    h = k.shape[0] // 2
    c = min(core_half, h)
    k = k[h - c:h + c + 1, h - c:h + c + 1]      # core only: tails / ringing don't matter
    good, var, model = fit['good'], fit['variance'], fit['model_full']
    resid = exposure.image.array - model
    inv = np.where(good, 1.0 / var, 0.0)
    num = fftconvolve(np.where(good, resid, 0.0) * inv, k[::-1, ::-1], mode='same')
    den = np.sqrt(np.clip(fftconvolve(inv, (k ** 2)[::-1, ::-1], mode='same'), 1e-30, None))
    snr = num / den
    sig = np.sqrt(np.nanmedian(var[good]))
    sky = good & (model < floor_sigma * sig) & (den > 0.7 * np.median(den[good]))
    if sky.sum() < 200:
        return 1.0
    zz = snr[sky]
    return 1.4826 * np.median(np.abs(zz - np.median(zz)))


# How far S/N values move from the pixel-level scale toward the matched-filter scale.
#   0.0 = pixel scale only (original behaviour, most sensitive, most fakes)
#   1.0 = full matched-filter scale (strictest; can drop real faint sources)
#   0.5 = halfway
MATCHED_NOISE_WEIGHT = 0.5


def empirical_noise_scale(fit, exposure, model_frac=0.01, matched=True,
                          weight=None):
    """Return the empirical noise scale of a fit's residual.

    Parameters
    ----------
    fit : `dict`
        Joint fit (from `joint_psf_fit_diag`).
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    model_frac : `float`, optional
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
    s = _pixel_noise_scale(fit, exposure, model_frac)
    if matched:
        w = MATCHED_NOISE_WEIGHT if weight is None else weight
        sm = _matched_filter_noise_scale(fit, exposure)
        if sm > s:
            s = s + w * (sm - s)
    return s


def amp_snr_full(fit, noise_scale=1.0):
    """Return honest amplitude S/N values, including the covariance of overlapping sources.

    ``S/N = amp / sqrt(diag(N^-1))`` with ``N = A^T W A`` over all sources, divided by the
    empirical noise scale.

    Parameters
    ----------
    fit : `dict`
        Joint fit (from `joint_psf_fit_diag`).
    noise_scale : `float`, optional
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
    sig = np.sqrt(np.clip(np.diag(cov), 1e-30, None)) * noise_scale
    return fit['amps'] / sig


def eliminate_weak(exposure, psf, peaks, fit, bad_bitmask, min_snr=5.0,
                   drop_frac=0.25, max_rounds=30, **fitkw):
    """Remove sources below ``min_snr`` by backward elimination.

    Each round drops the weakest ~``drop_frac`` of the sub-threshold sources (at least one) and
    re-solves, until every source has an honest S/N of at least ``min_snr``.  The last source is
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
    bad_bitmask : `int`
        Mask bits of pixels to ignore.
    min_snr : `float`, optional
        Minimum honest S/N.
    drop_frac : `float`, optional
        Fraction of the failing sources dropped per round.
    max_rounds : `int`, optional
        Maximum number of rounds.
    **fitkw
        Keyword arguments of `joint_psf_fit_diag`.

    Returns
    -------
    peaks : `list` [`tuple` [`float`, `float`]]
        Surviving source positions.
    fit : `dict`
        Joint fit of the survivors.
    """
    peaks = [tuple(p) for p in peaks]
    n0 = len(peaks)
    scale = empirical_noise_scale(fit, exposure)
    for _ in range(max_rounds):
        if len(peaks) <= 1:
            break
        scale = empirical_noise_scale(fit, exposure)
        snr = amp_snr_full(fit, scale)
        weak = np.flatnonzero(snr < min_snr)
        if len(weak) == 0:
            break
        order = weak[np.argsort(snr[weak])]
        n_drop = max(1, int(np.ceil(drop_frac * len(weak))))
        drop = set(order[:min(n_drop, len(peaks) - 1)].tolist())
        peaks = [p for i, p in enumerate(peaks) if i not in drop]
        fit = joint_psf_fit_diag(exposure, psf, peaks, bad_bitmask,
                                 var_model=fit['model_full'], **fitkw)
    print(f"  eliminate_weak: noise_scale={scale:.2f}, {n0} -> {len(peaks)} sources "
          f"(honest S/N>={min_snr:g})")
    return peaks, fit


# ---- NEW: PSF-shape constraint -----------------------------------------

def _fit_with_plane(T, D, w, dx, dy):
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


def shape_metrics(exposure, psf, peaks, fit, noise_scale=1.0, chi2_scale=None,
                  bg_sigma=4.0, width_scales=(0.6, 0.75, 1.0, 1.3, 1.7), min_pix=15):
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
    noise_scale : `float`, optional
        Divides the S/N values (matched-filter calibrated).
    chi2_scale : `float`, optional
        Divides chi2 (pixel-level calibrated; default ``noise_scale``).
    bg_sigma : `float`, optional
        Fit radius in units of the PSF sigma.
    width_scales : `tuple` [`float`], optional
        Trial PSF width scale factors.
    min_pix : `int`, optional
        Minimum number of good pixels to judge a source.

    Returns
    -------
    metrics : `dict` [`str`, `numpy.ndarray`]
        ``snr`` (local, noise-scaled), ``chi2red``, ``best_scale``, ``dchi2_width`` (> 0: the
        data prefer ``best_scale`` over the PSF width) and ``contrast``; NaN = not judged.
    """
    cs = noise_scale if chi2_scale is None else chi2_scale
    bbox = exposure.getBBox()
    X0, Y0 = bbox.getMinX(), bbox.getMinY()
    data, model = exposure.image.array, fit['model_full']
    var, good = fit['variance'], fit['good']
    ny, nx = data.shape
    shp = psf.computeShape(bbox.getCenter())
    sig = np.sqrt(0.5 * (shp.getIxx() + shp.getIyy()))
    R = bg_sigma * sig
    h = int(np.ceil(R))
    n = len(peaks)
    out = dict(snr=np.full(n, np.nan), chi2red=np.full(n, np.nan),
               best_scale=np.full(n, np.nan), dchi2_width=np.full(n, np.nan),
               contrast=np.full(n, np.nan))
    buf = np.zeros((ny, nx))
    scales = np.asarray(width_scales, float)
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
        if sel.sum() < min_pix:
            continue
        w = 1.0 / np.sqrt(var[sl][sel])
        Dv, dxv, dyv = D[sel], dx[sel], dy[sel]

        _, s_psf, chi2 = _fit_with_plane(P[sel], Dv, w, dxv, dyv)
        out['snr'][j] = s_psf / noise_scale
        out['chi2red'][j] = chi2 / max(sel.sum() - 4, 1) / cs ** 2     # pixel-level scale

        r2 = dxv ** 2 + dyv ** 2
        snr_w = np.array([_fit_with_plane(np.exp(-0.5 * r2 / (s * sig) ** 2),
                                          Dv, w, dxv, dyv)[1]
                          for s in scales]) / noise_scale
        if np.all(np.isfinite(snr_w)):
            k = int(np.argmax(snr_w))
            out['best_scale'][j] = scales[k]
            out['dchi2_width'][j] = snr_w[k] * abs(snr_w[k]) - snr_w[i1] * abs(snr_w[i1])
    return out


def judge_shape(m, min_snr_local=5.0, chi2_max=2.0, width_range=(0.75, 1.35),
                dchi2_width=9.0, trust_snr=30.0, min_contrast=20.0,
                weak_snr=10.0, dchi2_width_weak=4.0):
    """Decide which sources fail the PSF-shape tests.

    Sources with local S/N of at least ``trust_snr`` (including saturated ones) skip the chi2,
    width and contrast tests.  Untrusted sources with S/N below ``weak_snr`` use the looser
    width threshold ``dchi2_width_weak``.

    Parameters
    ----------
    m : `dict` [`str`, `numpy.ndarray`]
        Output of `shape_metrics`.
    min_snr_local : `float`, optional
        Minimum local S/N.
    chi2_max : `float`, optional
        Maximum reduced chi2.
    width_range : `tuple` [`float`, `float`], optional
        Allowed range of the best-fitting width scale.
    dchi2_width : `float`, optional
        Improvement in chi2 needed to reject on width.
    trust_snr : `float`, optional
        S/N above which a source is trusted.
    min_contrast : `float`, optional
        Minimum ratio of the source's peak to the neighbours' model flux at its position (None
        or 0 = off); 20 ~ (5 % wing error) x (4x margin).
    weak_snr : `float`, optional
        S/N below which ``dchi2_width_weak`` applies.
    dchi2_width_weak : `float`, optional
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
    bs, dw = m['best_scale'], m['dchi2_width']
    con = m['contrast']
    valid = np.isfinite(snr)
    trusted = valid & (snr >= trust_snr)
    weak = valid & ~trusted & (snr < weak_snr)
    dw_thr = np.where(weak, dchi2_width_weak, dchi2_width)

    f_snr = valid & (snr < min_snr_local)
    f_chi = valid & ~trusted & (chi2 > chi2_max)
    f_wid = (valid & ~trusted & np.isfinite(bs) &
             ((bs < width_range[0]) | (bs > width_range[1])) & (dw > dw_thr))
    f_con = (valid & ~trusted & np.isfinite(con) & (con < min_contrast)
             if min_contrast else np.zeros_like(valid))
    with np.errstate(divide='ignore', invalid='ignore'):
        terms = [min_snr_local / np.clip(np.nan_to_num(snr, nan=1e3), 1e-3, None),
                 np.nan_to_num(chi2 / chi2_max),
                 np.where(f_wid, np.nan_to_num(dw) / dw_thr, 0.0)]
        if min_contrast:
            terms.append(min_contrast / np.clip(np.nan_to_num(con, nan=1e6), 1e-3, None))
        score = np.maximum.reduce(terms)
    flags = dict(snr=f_snr, chi2=f_chi, width=f_wid, contrast=f_con)
    return f_snr | f_chi | f_wid | f_con, score, flags


def eliminate_nonpsf(exposure, psf, peaks, fit, bad_bitmask, shape_kwargs=None,
                     drop_frac=0.25, max_rounds=15, verbose=True, **fitkw):
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
    bad_bitmask : `int`
        Mask bits of pixels to ignore.
    shape_kwargs : `dict`, optional
        Keyword arguments of `shape_metrics` and `judge_shape`.
    drop_frac : `float`, optional
        Fraction of the failing sources dropped per round.
    max_rounds : `int`, optional
        Maximum number of rounds.
    verbose : `bool`, optional
        Print the rejected sources.
    **fitkw
        Keyword arguments of `joint_psf_fit_diag`.

    Returns
    -------
    peaks : `list` [`tuple` [`float`, `float`]]
        Surviving source positions.
    fit : `dict`
        Joint fit of the survivors.
    rejected : `list` [`tuple` [`float`, `float`, `float`]]
        (x, y, score) of the removed sources.
    """
    kw = dict(shape_kwargs or {})
    mkw = {k: kw.pop(k) for k in ('bg_sigma', 'width_scales', 'min_pix') if k in kw}
    peaks = [tuple(p) for p in peaks]
    n0, rejected = len(peaks), []
    why = dict(snr=0, chi2=0, width=0, contrast=0)       # tally of tests that flagged dropped sources
    scale = pscale = 1.0
    rows = []                      # per-rejected-source metrics, printed at the end
    for _ in range(max_rounds):
        if len(peaks) <= 1:
            break
        scale = empirical_noise_scale(fit, exposure)                  # for S/N values
        pscale = empirical_noise_scale(fit, exposure, matched=False)  # for chi2
        m = shape_metrics(exposure, psf, peaks, fit, noise_scale=scale,
                          chi2_scale=pscale, **mkw)
        fail, score, flags = judge_shape(m, **kw)
        idx = np.flatnonzero(fail)
        if len(idx) == 0:
            break
        order = idx[np.argsort(-score[idx])]
        n_drop = min(max(1, int(np.ceil(drop_frac * len(idx)))), len(peaks) - 1)
        drop = set(order[:n_drop].tolist())
        for i in drop:
            for key in why:
                why[key] += int(flags[key][i])
            rows.append((peaks[i][0], peaks[i][1], m['snr'][i], m['chi2red'][i],
                         m['best_scale'][i], m['dchi2_width'][i], m['contrast'][i],
                         '+'.join(k for k in why if flags[k][i]) or '-'))
        rejected += [(peaks[i][0], peaks[i][1], float(score[i])) for i in drop]
        peaks = [p for i, p in enumerate(peaks) if i not in drop]
        fit = joint_psf_fit_diag(exposure, psf, peaks, bad_bitmask,
                                 var_model=fit['model_full'], **fitkw)
    if verbose and rows:
        rows.sort(key=lambda r: -(r[2] if np.isfinite(r[2]) else -1e9))
        print("  shape-cut rejects, highest local S/N first (these are the ones most "
              "likely to be REAL):")
        print(f"    {'x':>8s} {'y':>8s} {'S/N':>6s} {'chi2':>6s} {'wscale':>6s} "
              f"{'dchi2w':>7s} {'contrast':>8s}  failed")
        for x, y, s, c, bs, dw, con, why_s in rows[:25]:
            print(f"    {x:8.1f} {y:8.1f} {s:6.1f} {c:6.2f} {bs:6.2f} {dw:7.1f} "
                  f"{con:8.1f}  {why_s}")
        if len(rows) > 25:
            print(f"    ... and {len(rows) - 25} more (all lower S/N)")
    if verbose:
        print(f"  shape cut: {n0} -> {len(peaks)} sources ({len(rejected)} removed; "
              f"flagged by: S/N {why['snr']}, chi2 {why['chi2']}, width {why['width']}, "
              f"contrast {why['contrast']} -- a source can trip several); "
              f"noise scales: S/N {scale:.2f}, pixel {pscale:.2f}")
    return peaks, fit, rejected


def plot_shape_rejects(exposure, rejected_by_iter, peaks=None,
                       title='sources removed by shape cut'):
    """Plot the sources removed by the shape cut over the data.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    rejected_by_iter : `list` [`list` [`tuple`]]
        Per iteration, the (x, y, score) of rejected sources (magenta x).
    peaks : array-like, optional
        Kept source positions (cyan).
    title : `str`, optional
        Plot title.
    """
    bbox = exposure.getBBox()
    data = exposure.image.array
    rej = np.array([(x, y) for r in rejected_by_iter for (x, y, _) in r]).reshape(-1, 2)
    fig, ax = plt.subplots(figsize=(8, 7.5))
    ax.imshow(data, extent=_extent(bbox), origin='lower', cmap='gray',
              norm=simple_norm(data, stretch='asinh', min_percent=1, max_percent=99.5))
    if peaks is not None and len(peaks):
        pk = _peaks_array(peaks)
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

def select_source_pages(snap, source_view='all', sources_per_page=6, max_pages=None):
    """Choose the sources that get per-source pages.

    Parameters
    ----------
    snap : `dict`
        Inner-loop snapshot (``fit``, ``prev_amp``, ``is_new``).
    source_view : `str`, optional
        'all', 'new' (added this iteration) or 'changed' (new, or amplitude changed by more
        than 10 %).
    sources_per_page : `int`, optional
        Sources per page.
    max_pages : `int`, optional
        Maximum number of pages (None = no limit).

    Returns
    -------
    pages : `list` [`list` [`int`]]
        Source indices per page, brightest first.
    """
    amps, prev_amp, is_new = snap['fit']['amps'], snap['prev_amp'], snap['is_new']
    idx = np.arange(len(amps))
    if source_view == 'new':
        idx = idx[is_new]
    elif source_view == 'changed':
        with np.errstate(divide='ignore', invalid='ignore'):
            rel = np.abs(amps - prev_amp) / np.abs(prev_amp)
        idx = idx[is_new | (rel > 0.10) | ~np.isfinite(rel)]
    idx = idx[np.argsort(-amps[idx])]
    pages = [[int(j) for j in idx[i:i + sources_per_page]]
             for i in range(0, len(idx), sources_per_page)]
    return pages[:max_pages] if max_pages is not None else pages


def source_colors(pages):
    """Assign each source on the pages a colour from ``SRC_PALETTE``.

    Parameters
    ----------
    pages : `list` [`list` [`int`]]
        Output of `select_source_pages`.

    Returns
    -------
    colors : `dict` [`int`, `str`]
        Colour per source index.
    """
    flat = [j for page in pages for j in page]
    return {j: SRC_PALETTE[k % len(SRC_PALETTE)] for k, j in enumerate(flat)}


def _view_window(bbox, px, py, view_half):
    """Return the (unclipped) array index of a source and the clipped window around it."""
    X0, Y0 = bbox.getMinX(), bbox.getMinY()
    ix, iy = int(round(px)) - X0, int(round(py)) - Y0
    xa, xb = max(ix - view_half, 0), min(ix + view_half + 1, bbox.getWidth())
    ya, yb = max(iy - view_half, 0), min(iy + view_half + 1, bbox.getHeight())
    return ix, iy, xa, xb, ya, yb


def _circle_radius(exposure, circle_radius=None):
    """Return the marker radius: the given one, else 0.75 FWHM, else 5 px."""
    if circle_radius is not None:
        return circle_radius
    try:
        return 0.75 * _fwhm(exposure.getPsf(), exposure.getBBox().getCenter())
    except Exception:
        return 5.0


def _mark_source(ax, x, y, color, radius, label=None):
    """Draw a circle (and optional label) at a source."""
    ax.add_patch(patches.Circle((x, y), radius, fill=False, ec=color, lw=1.8,
                                path_effects=_OUTLINE, zorder=6))
    if label is not None:
        ax.text(x + 0.75 * radius, y + 0.75 * radius, f"#{label}", color=color,
                fontsize=8, fontweight='bold', path_effects=_TEXT_OUTLINE, zorder=7)


def _frame_axes(ax, color, lw=3.0):
    """Colour the frame of an axes."""
    for s in ax.spines.values():
        s.set_edgecolor(color)
        s.set_linewidth(lw)


# =====================================================================
# Diagnostic joint fit
# =====================================================================

def joint_psf_fit_diag(exposure, psf, peak_positions, bad_bitmask,
                       amp_zero_tol=1e-6, support_frac=0.01,
                       var_model=None, sys_frac=0.03, poisson_k=0.0):
    """Fit all PSF amplitudes in one non-negative solve, keeping what is needed for diagnostics.

    Weighted, non-negative, one sparse solve against the original exposure.  PSF stamps are not
    clipped at zero (clipping would add a positive pedestal).

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Image to fit.
    psf : `lsst.afw.detection.Psf`
        PSF model.
    peak_positions : `list` [`tuple` [`float`, `float`]]
        Source positions in parent pixel coordinates.
    bad_bitmask : `int`
        Mask bits of pixels to ignore.
    amp_zero_tol : `float`, optional
        Amplitudes at or below this count as zero.
    support_frac : `float`, optional
        A source's significant support is where its stamp exceeds this fraction of its peak.
    var_model : `numpy.ndarray`, optional
        Model image used for the model-dependent variance terms (see `effective_variance`).
    sys_frac : `float`, optional
        Fractional model error added to the variance.
    poisson_k : `float`, optional
        Extra variance per unit model flux.

    Returns
    -------
    fit : `dict`
        ``amps``, ``model_full``, ``n_at_zero``; ``comps`` (per source a dict with ``stamp``,
        ``sig`` (significant support), ``ys``/``xs`` (slices of its fit region), ``amp``,
        ``bbox_parent``, or None if the stamp misses the image); coupled-source groups
        (``group_labels``, ``n_groups``, ``group_map``, ``group_chi2``); ``good`` (pixels
        used), ``chi2map``, ``variance`` (the effective variance), ``timings`` and the solver
        status (``lsq_status``, ``lsq_nit``, ``lsq_cost``).
    """
    exp_bbox = exposure.getBBox()
    ny, nx = exp_bbox.getHeight(), exp_bbox.getWidth()
    data = exposure.image.array
    variance = effective_variance(exposure, var_model, sys_frac, poisson_k)
    mask = exposure.mask.array

    good = ((mask & bad_bitmask) == 0) & np.isfinite(data) & np.isfinite(variance) & (variance > 0)
    pixel_index = -np.ones((ny, nx), dtype=np.int64)
    flat = np.flatnonzero(good.ravel())
    pixel_index.ravel()[flat] = np.arange(len(flat))
    n_pix = len(flat)
    d_vec = data.ravel()[flat]
    w_vec = 1.0 / variance.ravel()[flat]
    n_src = len(peak_positions)

    timings = {}
    t0 = time.perf_counter()
    comps = []
    r_l, c_l, v_l, br_l, bc_l = [], [], [], [], []
    for j, (x, y) in enumerate(peak_positions):
        stamp, bbox = psf_stamp(psf, x, y, exp_bbox, taper=WING_TAPER)
        if stamp is None:
            comps.append(None)
            continue
        ys_, xs_ = bbox_to_slices(bbox, exp_bbox)
        sub_idx = pixel_index[ys_, xs_]
        smax = stamp.max()
        sig = stamp > support_frac * smax if smax > 0 else np.zeros_like(stamp, bool)
        ok = (sub_idx >= 0) & (stamp != 0)
        okb = ok & sig
        r_l.append(sub_idx[ok]); c_l.append(np.full(int(ok.sum()), j)); v_l.append(stamp[ok])
        br_l.append(sub_idx[okb]); bc_l.append(np.full(int(okb.sum()), j))
        comps.append(dict(j=j, amp=0.0, stamp=stamp, sig=sig, ys=ys_, xs=xs_,
                          bbox_parent=(bbox.getMinX(), bbox.getMaxX(),
                                       bbox.getMinY(), bbox.getMaxY())))

    cat = lambda l, dt: np.concatenate(l) if l else np.array([], dtype=dt)
    rows, cols, vals = cat(r_l, np.int64), cat(c_l, np.int64), cat(v_l, float)
    A = coo_matrix((vals, (rows, cols)), shape=(n_pix, n_src)).tocsr()
    sqrt_w = np.sqrt(w_vec)
    Aw = csr_matrix(A.multiply(sqrt_w[:, None]))
    dw = d_vec * sqrt_w
    timings['joint fit: build A'] = time.perf_counter() - t0

    t0 = time.perf_counter()
    res = lsq_linear(Aw, dw, bounds=(0, np.inf), method='trf', lsmr_tol='auto', max_iter=2000)
    amps = res.x
    timings['joint fit: lsq_linear'] = time.perf_counter() - t0

    t0 = time.perf_counter()
    model_full = np.zeros((ny, nx))
    for j, c in enumerate(comps):
        if c is None:
            continue
        c['amp'] = float(amps[j])
        if amps[j] > amp_zero_tol:
            model_full[c['ys'], c['xs']] += amps[j] * c['stamp']

    br, bc = cat(br_l, np.int64), cat(bc_l, np.int64)
    B = coo_matrix((np.ones(len(br)), (br, bc)), shape=(n_pix, n_src)).tocsr()
    n_groups, labels = connected_components(B.T @ B, directed=False)
    group_map = -np.ones((ny, nx), dtype=np.int64)
    for j, c in enumerate(comps):
        if c is None:
            continue
        sub = group_map[c['ys'], c['xs']]
        sub[c['sig']] = labels[j]

    resid = data - model_full
    chi2map = np.zeros((ny, nx))
    chi2map[good] = resid[good] ** 2 / variance[good]
    group_chi2 = {}
    for g in range(n_groups):
        m = (group_map == g) & good
        group_chi2[g] = (chi2map[m].sum() / m.sum()) if m.sum() > 0 else np.nan
    timings['joint fit: model + groups'] = time.perf_counter() - t0

    return dict(amps=amps, model_full=model_full, n_at_zero=int(np.sum(amps <= amp_zero_tol)),
                comps=comps, group_labels=labels, n_groups=n_groups, group_map=group_map,
                group_chi2=group_chi2, good=good, chi2map=chi2map, timings=timings,
                lsq_status=res.status, lsq_nit=res.nit, lsq_cost=res.cost,
                variance=variance)


# =====================================================================
# Centroid refinement (fixes centroid-offset "dipole" residuals)
# =====================================================================

def compute_psf_derivatives(psf, x, y, bbox_clip, eps=0.25):
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
    bbox_clip : `lsst.geom.Box2I`
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
    img_c = psf.computeImage(geom.Point2D(x, y))
    img_xp = psf.computeImage(geom.Point2D(x + eps, y))
    img_xm = psf.computeImage(geom.Point2D(x - eps, y))
    img_yp = psf.computeImage(geom.Point2D(x, y + eps))
    img_ym = psf.computeImage(geom.Point2D(x, y - eps))

    common = geom.Box2I(img_c.getBBox())
    for im in (img_xp, img_xm, img_yp, img_ym):
        common.clip(im.getBBox())
    common.clip(bbox_clip)
    if common.isEmpty():
        return None, None, common

    dpdx = (img_xp[common].array - img_xm[common].array) / (2 * eps)
    dpdy = (img_yp[common].array - img_ym[common].array) / (2 * eps)
    return dpdx, dpdy, common


def windowed_centroid(image_array, xy0, x, y, half_width, max_frac_of_window=0.9):
    """Return the flux-weighted centroid in a box around (x, y).

    Negative pixels are clipped to zero.  Falls back to (x, y) if there is too little positive
    flux or the centroid lands too close to the window edge (dominated by a neighbour).

    Parameters
    ----------
    image_array : `numpy.ndarray`
        Image.
    xy0 : `tuple` [`int`, `int`]
        Parent pixel coordinates of ``image_array[0, 0]``.
    x, y : `float`
        Starting position in parent pixel coordinates.
    half_width : `int`
        Half-size of the box [px].
    max_frac_of_window : `float`, optional
        Reject centroids further than this fraction of ``half_width`` from the start.

    Returns
    -------
    cx, cy : `float`
        Centroid in parent pixel coordinates.
    """
    x0, y0 = xy0
    ny, nx = image_array.shape
    ix, iy = int(round(x - x0)), int(round(y - y0))
    xa, xb = max(ix - half_width, 0), min(ix + half_width + 1, nx)
    ya, yb = max(iy - half_width, 0), min(iy + half_width + 1, ny)
    if xb <= xa or yb <= ya:
        return x, y
    stamp = np.clip(image_array[ya:yb, xa:xb], 0, None)
    total = stamp.sum()
    if total <= 0:
        return x, y
    yy, xx = np.mgrid[ya:yb, xa:xb]
    cx = float((xx * stamp).sum() / total) + x0
    cy = float((yy * stamp).sum() / total) + y0
    if abs(cx - x) > max_frac_of_window * half_width or \
       abs(cy - y) > max_frac_of_window * half_width:
        return x, y
    return cx, cy


def refine_peak_positions(exposure, psf, peaks, amps, bad_bitmask,
                          min_amp_for_refine=50.0, max_shift=0.6, n_subiters=3,
                          min_good_pixels=10, core_radius_sigma=2.5,
                          variance=None, freeze_sep=None):
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
    bad_bitmask : `int`
        Mask bits of pixels to ignore.
    min_amp_for_refine : `float`, optional
        Fainter sources keep their position.
    max_shift : `float`, optional
        Maximum shift per sub-step [px].
    n_subiters : `int`, optional
        Number of Gauss-Newton sub-steps.
    min_good_pixels : `int`, optional
        Minimum number of good core pixels.
    core_radius_sigma : `float`, optional
        Core radius in units of the PSF sigma.
    variance : `numpy.ndarray`, optional
        Effective variance (default: the exposure's variance plane).
    freeze_sep : `float`, optional
        Sources with a brighter neighbour within this many px keep their position (it is
        degenerate with the neighbour's amplitude).

    Returns
    -------
    peaks : `list` [`tuple` [`float`, `float`]]
        Refined positions.
    """
    exp_bbox = exposure.getBBox()
    data, mask = exposure.image.array, exposure.mask.array
    if variance is None:
        variance = exposure.variance.array
    good = ((mask & bad_bitmask) == 0) & np.isfinite(data) & np.isfinite(variance) & (variance > 0)
    peaks = [tuple(p) for p in peaks]
    amps = np.asarray(amps, float).copy()

    frozen = np.zeros(len(peaks), bool)
    if freeze_sep is not None and len(peaks) > 1:
        P0 = np.asarray(peaks)
        for j, nb in enumerate(cKDTree(P0).query_ball_point(P0, freeze_sep)):
            frozen[j] = any(k != j and amps[k] > amps[j] for k in nb)

    shp = psf.computeShape(exp_bbox.getCenter())
    R2 = (core_radius_sigma * np.sqrt(0.5 * (shp.getIxx() + shp.getIyy()))) ** 2

    for _ in range(n_subiters):
        model = np.zeros_like(data)
        for (x, y), a in zip(peaks, amps):
            if a <= 0:
                continue
            stamp, bb = psf_stamp(psf, x, y, exp_bbox, taper=WING_TAPER)
            if stamp is None:
                continue
            ys, xs = bbox_to_slices(bb, exp_bbox)
            model[ys, xs] += a * stamp
        resid_all = data - model

        moved = 0.0
        for j, ((x, y), a) in enumerate(zip(peaks, amps)):
            if a < min_amp_for_refine or frozen[j]:
                continue
            dpdx, dpdy, bb = compute_psf_derivatives(psf, x, y, exp_bbox)
            if dpdx is None:
                continue
            ys, xs = bbox_to_slices(bb, exp_bbox)
            yy, xx = np.mgrid[bb.getMinY():bb.getMaxY() + 1, bb.getMinX():bb.getMaxX() + 1]
            g = good[ys, xs] & (((xx - x) ** 2 + (yy - y) ** 2) < R2)
            if g.sum() < min_good_pixels:
                continue
            P, _ = psf_stamp(psf, x, y, bb, taper=WING_TAPER)
            w = 1.0 / np.sqrt(variance[ys, xs][g])
            J = np.column_stack([P[g], a * dpdx[g], a * dpdy[g]]) * w[:, None]
            r = resid_all[ys, xs][g] * w
            try:
                delta, *_ = np.linalg.lstsq(J, r, rcond=None)
            except np.linalg.LinAlgError:
                continue
            dx, dy = np.clip(delta[1:], -max_shift, max_shift)
            peaks[j] = (x + dx, y + dy)
            moved = max(moved, abs(dx), abs(dy))
        if moved < 0.02:
            break
    return peaks


def redetect_and_refit_psf_adaptive(exposure, kernel_sizes=(41, 27, 21, 15, 11), **kwargs):
    """Run `redetect_and_refit_psf`, trying kernel sizes from largest to smallest.

    MakePsfCandidatesTask drops every candidate whose stamp touches masked or edge pixels; on a
    small cutout a large kernel can drop all of them.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Background-subtracted exposure.
    kernel_sizes : `tuple` [`int`], optional
        Kernel sizes to try [px].
    **kwargs
        Keyword arguments of `redetect_and_refit_psf`.

    Returns
    -------
    psf_new : `lsst.afw.detection.Psf`
        The new PSF (the previous one if every attempt failed).
    sources : `lsst.afw.table.SourceCatalog` or `None`
        Measured detections.
    starCat : `lsst.afw.table.SourceCatalog` or `None`
        Selected PSF stars.
    cellSet : `lsst.afw.math.SpatialCellSet` or `None`
        PSF candidates.
    kernel_size : `int` or `None`
        Kernel size that worked (None if all failed).
    """
    psf_prev = exposure.getPsf()
    fwhm_prev = _fwhm(psf_prev, exposure.getBBox().getCenter())
    for ks in kernel_sizes:
        print(f"  [adaptive PSF refit] trying kernel_size={ks}...")
        try:
            psf_new, sources_new, starCat_new, cellSet_new = redetect_and_refit_psf(
                exposure, kernel_size=ks, **kwargs)
        except Exception as e:
            print(f"    kernel_size={ks} raised {type(e).__name__}: {e} -- trying smaller")
            continue
        if cellSet_new is not None:
            fwhm_new = _fwhm(psf_new, exposure.getBBox().getCenter())
            n_used = sum(1 for cell in cellSet_new.getCellList() for c in cell
                         if c.getStatus() == GOOD_STATUS)
            print(f"    SUCCESS at kernel_size={ks}: {n_used} candidates used, "
                  f"FWHM {fwhm_prev:.3f} -> {fwhm_new:.3f} px")
            exposure.setPsf(psf_new)
            return psf_new, sources_new, starCat_new, cellSet_new, ks
        print(f"    kernel_size={ks} failed (no viable candidates) -- trying smaller")

    print(f"  [adaptive PSF refit] FAILED at every kernel_size tried {kernel_sizes} -- "
          f"keeping previous PSF (FWHM {fwhm_prev:.3f} px). Consider a larger cutout for "
          f"the PSF-refit step, or refit_psf=False / 'stack'.")
    exposure.setPsf(psf_prev)
    return psf_prev, None, None, None, None


# =====================================================================
# Per-iteration figures
# =====================================================================

def plot_iteration_overview(exposure, snap, hist, outer_idx, norm_data,
                            colors=None, view_half=20, circle_radius=None):
    """Plot the overview of one inner iteration.

    Detection input, detections, fitted sources, model, residual and the iteration histories.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    snap : `dict`
        Inner-loop snapshot of the iteration.
    hist : `dict`
        Histories (``resid_std``, ``n_peaks``, ``n_zero``).
    outer_idx : `int`
        Outer-pass index.
    norm_data : `astropy.visualization.ImageNormalize`
        Display normalisation of the data.
    colors : `dict` [`int`, `str`], optional
        Colours of the sources with per-source pages.
    view_half : `int`, optional
        Half-size of the per-source windows [px].
    circle_radius : `float`, optional
        Marker radius [px].
    """
    colors = colors or {}
    radius = _circle_radius(exposure, circle_radius)
    bbox = exposure.getBBox()
    X0, Y0 = bbox.getMinX(), bbox.getMinY()
    ext = _extent(bbox)
    xs, ys = _coords(bbox)
    data = exposure.image.array
    fit = snap['fit']
    model = fit['model_full']
    resid = snap['residual']
    peaks = snap['merged_arr']
    is_new = snap['is_new']
    bad = snap['bad_mask']
    it = snap['it']

    fig = plt.figure(figsize=(21, 10.5))
    gs = fig.add_gridspec(2, 4, hspace=0.32, wspace=0.28)

    # 1. detection input
    ax = fig.add_subplot(gs[0, 0])
    det_in = snap['det_in']
    im = ax.imshow(det_in, extent=ext, origin='lower', cmap=CMAP,
                   norm=simple_norm(det_in, stretch='asinh', min_percent=1, max_percent=99.5))
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    if snap['footprint_mask'].any():
        ax.contour(xs, ys, snap['footprint_mask'].astype(float), levels=[0.5],
                   colors='white', linewidths=0.6)
    raw = _peaks_array(snap['raw_peaks'])
    col = _peaks_array(snap['collapsed_peaks'])
    if len(raw):
        ax.scatter(raw[:, 0], raw[:, 1], marker='+', c='w', s=28,
                   label=f'raw peaks ({len(raw)})')
    if len(col):
        ax.scatter(col[:, 0], col[:, 1], facecolors='none', edgecolors='yellow', s=70,
                   lw=1.0, label=f'after sat-collapse + rejection ({len(col)})')
    _overlay_bad(ax, bad, ext)
    ax.legend(fontsize=7, loc='upper right')
    what = 'data (bg-subtracted)' if it == 0 else 'residual from previous fit'
    ax.set_title(f"1. Detection input: {what}\nthreshold = {snap['threshold']:.2f}σ  (white = footprints)",
                 fontsize=10)
    _clean_axes(ax)

    # 2. fit regions & the cutouts shown on the per-source pages
    ax = fig.add_subplot(gs[0, 1])
    im = ax.imshow(data, extent=ext, origin='lower', norm=norm_data, cmap=CMAP)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    for c in fit['comps']:
        if c is None or c['j'] in colors:
            continue
        x0b, x1b, y0b, y1b = c['bbox_parent']
        ax.add_patch(patches.Rectangle((x0b - 0.5, y0b - 0.5), x1b - x0b + 1, y1b - y0b + 1,
                                       fill=False, ec='0.75', lw=0.4, alpha=0.5))
    labels = fit['group_labels']
    sizes = np.bincount(labels, minlength=fit['n_groups'])
    for g in range(fit['n_groups']):
        gm = fit['group_map'] == g
        if gm.any():
            ax.contour(xs, ys, gm.astype(float), levels=[0.5], colors='w',
                       linewidths=0.8, linestyles='--')
    if len(peaks):
        ax.scatter(peaks[~is_new, 0], peaks[~is_new, 1], facecolors='none', edgecolors='w',
                   s=25, lw=0.7, label='existing peak')
        ax.scatter(peaks[is_new, 0], peaks[is_new, 1], marker='*', c='red', s=70,
                   edgecolors='w', lw=0.4, label='NEW this iteration')
    rej = _peaks_array([(x, y) for (x, y, _) in snap.get('shape_rejected', [])])
    if len(rej):
        ax.scatter(rej[:, 0], rej[:, 1], marker='x', c='#ff00ff', s=45, lw=1.4,
                   label=f'removed by shape cut ({len(rej)})')
    for j, colr in colors.items():
        px, py = peaks[j]
        ix, iy, xa, xb, ya, yb = _view_window(bbox, px, py, view_half)
        ax.add_patch(patches.Rectangle((X0 + xa - 0.5, Y0 + ya - 0.5), xb - xa, yb - ya,
                                       fill=False, ec=colr, lw=1.8, path_effects=_OUTLINE))
        _mark_source(ax, px, py, colr, radius, label=j)
    _overlay_bad(ax, bad, ext)
    ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
    ax.legend(fontsize=7, loc='upper right')
    ax.set_title(f"2. Coloured box + circle = cutout on the per-source pages (#j)\n"
                 f"gray = other PSF stamps, white dashed = coupled groups ({fit['n_groups']}, "
                 f"max {sizes.max() if len(sizes) else 0} src)", fontsize=10)
    _clean_axes(ax)

    # 3. model
    ax = fig.add_subplot(gs[0, 2])
    im = ax.imshow(model, extent=ext, origin='lower', norm=norm_data, cmap=CMAP)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_title(f"3. Joint PSF model ({len(peaks)} sources, "
                 f"{fit['n_at_zero']} at zero amp)", fontsize=10)
    _clean_axes(ax)

    # 4. residual
    ax = fig.add_subplot(gs[0, 3])
    rs = np.nanpercentile(np.abs(resid), 99)
    rs = rs if np.isfinite(rs) and rs > 0 else 1.0
    im = ax.imshow(resid, extent=ext, origin='lower', cmap='RdBu_r', vmin=-rs, vmax=rs)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    for g in range(fit['n_groups']):
        gm = fit['group_map'] == g
        if gm.any():
            ax.contour(xs, ys, gm.astype(float), levels=[0.5], colors='k',
                       linewidths=0.5, alpha=0.5)
    ax.set_title(f"4. Residual = data - model\nstd {snap['std_before']:.3f} → {snap['std_after']:.3f}",
                 fontsize=10)
    _clean_axes(ax)

    # 5. normalized residual histogram (EFFECTIVE variance)
    ax = fig.add_subplot(gs[1, 0])
    good = fit['good']
    var_eff = fit['variance']
    z = (resid / np.sqrt(np.where(var_eff > 0, var_eff, np.inf)))[good]
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
    colors_bar = np.where(is_new[order], '#d62728', '#7f7f7f')
    ax.bar(np.arange(len(amps)), amps[order], color=colors_bar)
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
    xs_h = np.arange(len(hist['resid_std']))
    ax.plot(xs_h, hist['resid_std'], 'o-', color='#DD8452')
    ax.set_xlabel('iteration (0 = before any fit)'); ax.set_ylabel('residual std', color='#DD8452')
    ax2 = ax.twinx()
    ax2.plot(xs_h, hist['n_peaks'], 's--', color='#4C72B0')
    ax2.plot(xs_h, hist['n_zero'], '^:', color='k')
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
                 f"lsq_linear: {fit['lsq_nit']} its, status {fit['lsq_status']}", fontsize=10)
    ax.grid(axis='x', alpha=0.3)

    fig.suptitle(f"Outer pass {outer_idx + 1}  |  inner iteration {it + 1}  |  "
                 f"{len(peaks)} sources ({int(is_new.sum())} new)",
                 fontsize=14, fontweight='medium', y=0.995)
    plt.show()
    plt.close(fig)


def plot_source_pages(exposure, snap, outer_idx, pages=None, colors=None,
                      source_view='all', sources_per_page=6, max_pages=None,
                      view_half=20, circle_radius=None):
    """Plot one row per source: data | model | residual | x-section | y-section.

    Superseded by the second definition of this function further down (kept from the
    notebook).  Each source has its own colour (matching overview panel 2); the dashed box is
    the PSF stamp used by the solver, gold the significant support.
    """
    bbox = exposure.getBBox()
    X0, Y0 = bbox.getMinX(), bbox.getMinY()
    ny, nx = bbox.getHeight(), bbox.getWidth()
    data = exposure.image.array
    fit = snap['fit']
    model = fit['model_full']
    resid = data - model
    peaks, amps = snap['merged_arr'], fit['amps']
    prev_amp, is_new, bad = snap['prev_amp'], snap['is_new'], snap['bad_mask']

    if pages is None:
        pages = select_source_pages(snap, source_view, sources_per_page, max_pages)
    if not pages:
        print("  (no sources selected for per-source pages)")
        return
    if colors is None:
        colors = source_colors(pages)
    radius = _circle_radius(exposure, circle_radius)

    for p, page in enumerate(pages):
        n = len(page)
        fig, axes = plt.subplots(n, 5, figsize=(21, 3.3 * n), squeeze=False)
        for r, j in enumerate(page):
            col = colors[int(j)]
            px, py = peaks[j]
            comp = fit['comps'][j]
            ix, iy, xa, xb, ya, yb = _view_window(bbox, px, py, view_half)
            ix = int(np.clip(ix, 0, nx - 1)); iy = int(np.clip(iy, 0, ny - 1))
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
            own_sig = np.zeros((ny, nx), bool)
            if comp is not None:
                own[comp['ys'], comp['xs']] = comp['amp'] * comp['stamp']
                own_sig[comp['ys'], comp['xs']] = comp['sig']
            wxs, wys = X0 + np.arange(xa, xb), Y0 + np.arange(ya, yb)

            for c_i, (arr, kw) in enumerate([
                    (dwin, dict(norm=nrm, cmap=CMAP)),
                    (mwin, dict(norm=nrm, cmap=CMAP)),
                    (rwin, dict(cmap='RdBu_r', vmin=-rs, vmax=rs))]):
                ax = axes[r, c_i]
                im = ax.imshow(arr, extent=wext, origin='lower', **kw)
                plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                if own_sig[ya:yb, xa:xb].any():
                    ax.contour(wxs, wys, own_sig[ya:yb, xa:xb].astype(float), levels=[0.5],
                               colors='gold', linewidths=1.2)
                ax.scatter([px], [py], marker='x', c=col, s=40, lw=1.5, zorder=6)
                _mark_source(ax, px, py, col, radius, label=j if c_i == 0 else None)
                inwin = ((peaks[:, 0] > wext[0]) & (peaks[:, 0] < wext[1]) &
                         (peaks[:, 1] > wext[2]) & (peaks[:, 1] < wext[3]))
                inwin[j] = False
                if inwin.any():
                    ax.scatter(peaks[inwin, 0], peaks[inwin, 1], marker='+', c='w', s=30, lw=0.9)
                _overlay_bad(ax, bad[ya:yb, xa:xb], wext)
                if comp is not None:
                    x0b, x1b, y0b, y1b = comp['bbox_parent']
                    ax.add_patch(patches.Rectangle((x0b - 0.5, y0b - 0.5), x1b - x0b + 1,
                                                   y1b - y0b + 1, fill=False, ec=col, lw=1.5,
                                                   ls='--', path_effects=_OUTLINE))
                ax.set_xlim(wext[0], wext[1]); ax.set_ylim(wext[2], wext[3])
                _clean_axes(ax)
                _frame_axes(ax, col)

            g = fit['group_labels'][j]
            tag = 'NEW' if is_new[j] else 'refit'
            dtxt = ''
            if np.isfinite(prev_amp[j]) and prev_amp[j] != 0:
                dtxt = f"  Δ={100 * (amps[j] - prev_amp[j]) / abs(prev_amp[j]):+.1f}%"
            zflag = '  [AT ZERO BOUND]' if amps[j] <= 1e-6 else ''
            axes[r, 0].set_title(f"#{j} ({px:.1f},{py:.1f}) [{tag}] G{g}", fontsize=9, loc='left',
                                 bbox=dict(fc=col, ec='k', pad=2))
            axes[r, 1].set_title(f"model  amp={amps[j]:.1f}{dtxt}{zflag}", fontsize=9, loc='left')
            axes[r, 2].set_title(f"residual  (group χ²ν={fit['group_chi2'].get(g, np.nan):.2f})",
                                 fontsize=9, loc='left')
            if r == 0:
                axes[r, 0].text(0.02, 0.02,
                                'circle = this source\ndashed box = PSF stamp\ngold = significant support',
                                transform=axes[r, 0].transAxes, fontsize=6, color='w',
                                bbox=dict(fc='k', alpha=0.5, lw=0))

            for c_i, axis in ((3, 'x'), (4, 'y')):
                ax = axes[r, c_i]
                if axis == 'x':
                    coord, sl, centre = wxs, (iy, slice(xa, xb)), px
                    c0, c1 = ((comp['bbox_parent'][0], comp['bbox_parent'][1]) if comp else (None, None))
                else:
                    coord, sl, centre = wys, (slice(ya, yb), ix), py
                    c0, c1 = ((comp['bbox_parent'][2], comp['bbox_parent'][3]) if comp else (None, None))
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
                ax.fill_between(coord, 0, 1, where=own_sig[sl], transform=ax.get_xaxis_transform(),
                                color='gold', alpha=0.15)
                if badline.any():
                    ax.plot(coord[badline], np.zeros(int(badline.sum())), '|', color='red',
                            transform=ax.get_xaxis_transform(), ms=8)
                ax.set_xlim(coord[0], coord[-1])
                ax.set_xlabel(f'{axis} [px]')
                ax.set_title(f"{axis}-cross-section through source #{j}", fontsize=9)
                _frame_axes(ax, col)
                if r == 0 and c_i == 3:
                    ax.legend(fontsize=6, loc='upper right')
        fig.suptitle(f"Outer {outer_idx + 1} / inner iter {snap['it'] + 1}: per-source fits "
                     f"(page {p + 1}/{len(pages)}, sorted by amplitude, view={source_view}); "
                     f"colours match overview panel 2", fontsize=13, y=1.0)
        plt.tight_layout()
        plt.show()
        plt.close(fig)


# =====================================================================
# Instrumented inner loop
# =====================================================================

def run_inner_loop_verbose(exposure, psf, bad_bitmask, outer_idx=0,
                           n_iterations=10, threshold_start=8.0, threshold_step=-1.5,
                           threshold_min=4.5, min_pixels=1, min_new_peaks=1,
                           convergence_frac=0.02, viz=True, plot_every=1,
                           source_view='all', sources_per_page=6, max_pages_per_iter=None,
                           view_half=20, support_frac=0.01, refine_centroids=True,
                           centroid_min_amp=50.0, n_centroid_rounds=4,
                           centroid_shift_tol=0.03, centroid_shift_schedule=None,
                           coarse_centroid=True, coarse_centroid_hw_sigma=2.5,
                           merge_separation_sigma=1.5,
                           sys_frac=0.03, min_snr=5.0, prune_sep_sigma=1.5,
                           reject_fwhm=0.75, bright_peak_sigma=100.0, freeze_fwhm=1.0, poisson_k=0.0,
                           near_bright_fwhm=2.0, near_bright_snr=10.0,
                           edge_margin_fwhm=0.5, use_eliminate_weak=True,
                           use_shape_cut=True, shape_kwargs=None,
                           timer=TIMER, verbose=True, initial_peaks=None):
    """Run the instrumented inner loop (detect, fit, refine, prune, repeat).

    An instrumented version of `iterative_joint_detection` with centroid refinement and guards
    against over-counting sources.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Image to fit.
    psf : `lsst.afw.detection.Psf`
        PSF model.
    bad_bitmask : `int`
        Mask bits of pixels to ignore.
    outer_idx : `int`, optional
        Outer-pass index (for labels).
    n_iterations : `int`, optional
        Maximum number of iterations.
    threshold_start, threshold_step, threshold_min : `float`, optional
        Detection threshold schedule [sigma].
    min_pixels : `int`, optional
        Minimum footprint size.
    min_new_peaks : `int`, optional
        Stop when an iteration adds fewer peaks than this.
    convergence_frac : `float`, optional
        Stop when the residual standard deviation improves by less than this fraction.
    viz : `bool`, optional
        Draw diagnostics.
    plot_every : `int`, optional
        Draw every n-th iteration (the last always).
    source_view : `str`, optional
        See `select_source_pages`.
    sources_per_page : `int`, optional
        See `select_source_pages`.
    max_pages_per_iter : `int`, optional
        Maximum per-source pages per iteration.
    view_half : `int`, optional
        Half-size of the per-source windows [px].
    support_frac : `float`, optional
        See `joint_psf_fit_diag`.
    refine_centroids : `bool`, optional
        Refine positions (`refine_peak_positions`).
    centroid_min_amp : `float`, optional
        Fainter sources keep their position.
    n_centroid_rounds : `int`, optional
        Maximum centroid rounds per iteration.
    centroid_shift_tol : `float`, optional
        Stop centroiding when the largest shift is below this [px].
    centroid_shift_schedule : `list` [`float`], optional
        Maximum shift per round [px].
    coarse_centroid : `bool`, optional
        Start new peaks from a windowed centroid.
    coarse_centroid_hw_sigma : `float`, optional
        Half-width of that window in PSF sigma.
    merge_separation_sigma : `float`, optional
        Merge new peaks closer than this many PSF sigma.
    sys_frac : `float`, optional
        Fractional PSF-shape error added in quadrature to the variance.
    min_snr : `float`, optional
        S/N floor (cheap version in `prune_and_merge`, honest one in `eliminate_weak`).
    prune_sep_sigma : `float`, optional
        Pairs closer than this many PSF sigma are merged after refinement.
    reject_fwhm : `float`, optional
        New peaks within this many FWHM of a bright source are discarded.
    bright_peak_sigma : `float`, optional
        A source is bright if its peak exceeds this many noise sigma.
    freeze_fwhm : `float`, optional
        Sources with a brighter neighbour within this many FWHM keep their position.
    poisson_k : `float`, optional
        Extra variance per unit model flux.
    near_bright_fwhm : `float`, optional
        Radius around bright stars with the stricter S/N cut [FWHM].
    near_bright_snr : `float`, optional
        S/N required there.
    edge_margin_fwhm : `float`, optional
        Drop sources closer than this many FWHM to the edge.
    use_eliminate_weak : `bool`, optional
        Backward elimination on the honest S/N.
    use_shape_cut : `bool`, optional
        Backward elimination on PSF-likeness (`eliminate_nonpsf`).
    shape_kwargs : `dict`, optional
        Keyword arguments of `shape_metrics` / `judge_shape`.
    timer : `StepTimer`, optional
        Timer collecting the step durations.
    verbose : `bool`, optional
        Print progress.
    initial_peaks : `list` [`tuple` [`float`, `float`]], optional
        Warm-start source positions.

    Returns
    -------
    result : `dict`
        ``model_full``, ``residual_full``, ``peaks``, ``amps``, the histories
        (``resid_std_history``, ``n_peaks_history``, ``n_zero_amp_history``), ``first_seen`` and
        ``snaps`` (one dict per iteration with everything needed for plotting).
    """
    exp_bbox = exposure.getBBox()
    ny, nx = exp_bbox.getHeight(), exp_bbox.getWidth()
    norm_data = _data_norm(exposure.image.array)
    bad_mask = (exposure.mask.array & bad_bitmask) != 0

    _sh = psf.computeShape(exp_bbox.getCenter())
    psf_sigma = np.sqrt(0.5 * (_sh.getIxx() + _sh.getIyy()))
    fwhm_px = 2.355 * psf_sigma
    psf_peak = psf.computeKernelImage(exp_bbox.getCenter()).array.max()
    noise = np.sqrt(np.nanmedian(exposure.variance.array))

    working_residual = exposure.clone()
    working_residual.setPsf(psf)

    all_peaks, prev_amps = [], np.array([])
    hist = dict(resid_std=[np.nanstd(working_residual.image.array)], n_peaks=[0], n_zero=[0])
    first_seen, snaps = {}, []
    threshold = threshold_start
    model_full = np.zeros((ny, nx))
    amps_final = np.array([])

    if initial_peaks is not None and len(initial_peaks):
        # ---- warm start (cutout_adaptive): fit the given peaks, then keep detecting ----
        all_peaks = [tuple(p) for p in initial_peaks]
        _f0 = joint_psf_fit_diag(exposure, psf, all_peaks, bad_bitmask,
                                 support_frac=support_frac, var_model=None,
                                 sys_frac=sys_frac, poisson_k=poisson_k)
        _f0 = joint_psf_fit_diag(exposure, psf, all_peaks, bad_bitmask,
                                 support_frac=support_frac, var_model=_f0['model_full'],
                                 sys_frac=sys_frac, poisson_k=poisson_k)
        model_full, amps_final, prev_amps = _f0['model_full'], _f0['amps'], _f0['amps']
        working_residual.image.array[:] = exposure.image.array - model_full
        working_residual.variance.array[:] = _f0['variance']
        hist = dict(resid_std=[np.nanstd(working_residual.image.array)],
                    n_peaks=[len(all_peaks)], n_zero=[_f0['n_at_zero']])
        print(f"  warm start from {len(all_peaks)} peaks "
              f"(residual std {hist['resid_std'][0]:.3f})")

    for it in range(n_iterations):
        print(f"\n--- outer {outer_idx + 1} | inner {it + 1}/{n_iterations} | "
              f"threshold={threshold:.2f}σ | residual std={hist['resid_std'][-1]:.3f} "
              f"| peaks so far={len(all_peaks)} ---")
        times = {}

        def S(label):
            return timer.step(label, outer=outer_idx, inner=it)

        det_in = working_residual.image.array.copy()

        # --- detection on the current residual
        with S('detect on residual') as rec:
            cfg = SourceDetectionTask.ConfigClass()
            cfg.thresholdValue = threshold
            cfg.thresholdType = "pixel_stdev"
            cfg.minPixels = min_pixels
            cfg.doTempLocalBackground = True
            cfg.tempLocalBackground.binSize = 16
            cfg.reEstimateBackground = False
            schema_it = afwTable.SourceTable.makeMinimalSchema()
            task_it = SourceDetectionTask(schema=schema_it, config=cfg)
            table_it = afwTable.SourceTable.make(schema_it)
            try:
                det = task_it.run(table_it, working_residual)
            except TooManyMaskedPixelsError:
                print("  background estimation failed on residual -- stopping")
                det = None
        if det is None:
            break
        times['detect on residual'] = rec['seconds']
        footprints = det.positive.getFootprints()
        raw_peaks = [(pk.getFx(), pk.getFy()) for fp in footprints for pk in fp.getPeaks()]
        fp_mask = _footprints_to_mask(footprints, exposure)
        print(f"  {len(raw_peaks)} raw peaks detected")

        with S('collapse saturated peaks') as rec:
            new_peaks, _ = collapse_saturated_peaks(raw_peaks, exposure)
        times['collapse saturated peaks'] = rec['seconds']

        # reject new peaks sitting in the PSF-mismatch zone of a bright source
        if len(all_peaks) and len(new_peaks):
            with S('reject near bright') as rec:
                is_bright = np.asarray(prev_amps) * psf_peak > bright_peak_sigma * noise
                new_peaks, n_rej = reject_near_bright(new_peaks, all_peaks, is_bright,
                                                      reject_fwhm * fwhm_px)
            times['reject near bright'] = rec['seconds']
            if n_rej:
                print(f"  rejected {n_rej} new peak(s) within {reject_fwhm * fwhm_px:.1f}px "
                      f"of a bright source")

        if len(new_peaks) < min_new_peaks and it > 0:
            print("  too few new peaks -- treating as converged, stopping")
            break

        if coarse_centroid and len(new_peaks) > 0:
            with S('coarse centroid (new peaks)') as rec:
                hw = max(int(round(coarse_centroid_hw_sigma * psf_sigma)), 3)
                xy0_local = (exp_bbox.getMinX(), exp_bbox.getMinY())
                new_peaks = [windowed_centroid(det_in, xy0_local, x, y, hw)
                             for (x, y) in new_peaks]
            times['coarse centroid (new peaks)'] = rec['seconds']

        with S('drop peaks on masked pixels') as rec:
            new_peaks, n_dropped = drop_masked_peaks(new_peaks, exposure)
        times['drop masked peaks'] = rec['seconds']
        if n_dropped:
            print(f"  dropped {n_dropped} new peak(s) sitting on EDGE/NO_DATA pixels")

        with S('merge close peaks') as rec:
            merge_min_sep = None
            if merge_separation_sigma is not None:
                merge_min_sep = merge_separation_sigma * psf_sigma
            merged = merge_close_peaks(all_peaks + new_peaks, exposure,
                                       min_separation=merge_min_sep)
        times['merge close peaks'] = rec['seconds']
        n_in = len(all_peaks) + len(new_peaks)
        if merge_min_sep is not None:
            print(f"  merge_close_peaks: {n_in} -> {len(merged)} peaks "
                  f"(min separation {merge_min_sep:.1f} px)")
        if len(merged) == 0:
            print("  no peaks at all -- stopping")
            break

        # var_model=model_full (previous model) => systematic variance floor
        fit = joint_psf_fit_diag(exposure, psf, merged, bad_bitmask,
                                 support_frac=support_frac,
                                 var_model=model_full, sys_frac=sys_frac,
                                 poisson_k=poisson_k)
        for lab, sec in fit['timings'].items():
            timer.add(lab, sec, outer=outer_idx, inner=it)
            times[lab] = sec

        if refine_centroids:
            # alternate (refine position -> re-solve amplitude); step cap SHRINKS per round
            schedule = list(centroid_shift_schedule) if centroid_shift_schedule else \
                [1.5, 0.8, 0.4, 0.2, 0.1]
            if len(schedule) < n_centroid_rounds:
                schedule += [schedule[-1]] * (n_centroid_rounds - len(schedule))

            prev_positions = _peaks_array(merged)
            for round_i in range(n_centroid_rounds):
                step_cap = schedule[round_i]
                with S(f'refine centroids (round {round_i + 1})') as rec:
                    merged = refine_peak_positions(
                        exposure, psf, merged, fit['amps'], bad_bitmask,
                        min_amp_for_refine=centroid_min_amp, max_shift=step_cap,
                        variance=fit['variance'], freeze_sep=freeze_fwhm * fwhm_px)
                times[f'refine centroids r{round_i + 1}'] = rec['seconds']
                with S(f'joint fit: re-solve (round {round_i + 1})') as rec:
                    fit = joint_psf_fit_diag(exposure, psf, merged, bad_bitmask,
                                             support_frac=support_frac,
                                             var_model=fit['model_full'], sys_frac=sys_frac,
                                             poisson_k=poisson_k)
                for lab, sec in fit['timings'].items():
                    timer.add(lab, sec, outer=outer_idx, inner=it)
                times[f'joint fit: re-solve r{round_i + 1}'] = rec['seconds']

                new_positions = _peaks_array(merged)
                observed_shift = (np.max(np.linalg.norm(new_positions - prev_positions, axis=1))
                                  if len(new_positions) else 0.0)
                prev_positions = new_positions
                print(f"    centroid round {round_i + 1} (step cap {step_cap:.2f}px): "
                      f"max position shift = {observed_shift:.3f} px")
                if observed_shift < centroid_shift_tol:
                    break

        # prune low-S/N sources + merge pairs that collapsed together during refinement
        with S('prune + merge') as rec:
            bright_flags = fit['amps'] * psf_peak > bright_peak_sigma * noise
            merged = prune_and_merge(fit, merged, min_snr=min_snr,
                                     min_sep=prune_sep_sigma * psf_sigma,
                                     bright_flags=bright_flags,
                                     near_bright_sep=near_bright_fwhm * fwhm_px,
                                     near_bright_snr=near_bright_snr)
        times['prune + merge'] = rec['seconds']

        # Must come AFTER prune_and_merge (which indexes fit['comps'] by the pre-drop list).
        merged, n_oob = drop_out_of_bounds(merged, exposure)
        if n_oob:
            print(f"  dropped {n_oob} source(s) that drifted outside the image")

        if edge_margin_fwhm:
            merged, n_edge = drop_out_of_bounds(merged, exposure,
                                                margin=-edge_margin_fwhm * fwhm_px)
            if n_edge:
                print(f"  dropped {n_edge} source(s) within {edge_margin_fwhm} FWHM "
                      f"({edge_margin_fwhm * fwhm_px:.1f}px) of the edge")

        if len(merged) == 0:
            print("  everything pruned/dropped -- stopping")
            break

        # fresh solve on the pruned list (fit must match `merged` before elimination)
        with S('joint fit: re-solve (post-prune)') as rec:
            fit = joint_psf_fit_diag(exposure, psf, merged, bad_bitmask,
                                     support_frac=support_frac,
                                     var_model=fit['model_full'], sys_frac=sys_frac,
                                     poisson_k=poisson_k)
        for lab, sec in fit['timings'].items():
            timer.add(lab, sec, outer=outer_idx, inner=it)
        times['joint fit: re-solve (post-prune)'] = rec['seconds']

        # backward elimination on the honest S/N
        if use_eliminate_weak:
            with S('eliminate weak (honest S/N)') as rec:
                merged, fit = eliminate_weak(exposure, psf, merged, fit, bad_bitmask,
                                             min_snr=min_snr, support_frac=support_frac,
                                             sys_frac=sys_frac, poisson_k=poisson_k)
            times['eliminate weak'] = rec['seconds']

        # NEW: PSF-shape constraint (local-plane S/N, chi2, width test)
        shape_rejected = []
        if use_shape_cut:
            with S('shape cut (PSF-likeness)') as rec:
                merged, fit, shape_rejected = eliminate_nonpsf(
                    exposure, psf, merged, fit, bad_bitmask,
                    shape_kwargs=shape_kwargs, support_frac=support_frac,
                    sys_frac=sys_frac, poisson_k=poisson_k)
            times['shape cut'] = rec['seconds']

        model_full, amps = fit['model_full'], fit['amps']

        with S('rebuild residual') as rec:
            working_residual = exposure.clone()
            working_residual.image.array[:] = exposure.image.array - model_full
            # detector thresholds on the EFFECTIVE variance
            working_residual.variance.array[:] = fit['variance']
            working_residual.setPsf(psf)
        times['rebuild residual'] = rec['seconds']

        merged_arr = _peaks_array(merged)
        prior_arr = _peaks_array(all_peaks)
        if len(prior_arr):
            d, _ = cKDTree(prior_arr).query(merged_arr)
            is_new = d > 1.0
        else:
            is_new = np.ones(len(merged_arr), bool)
        prev_amp = _match_prev_amps(merged_arr, prior_arr, prev_amps)

        for (x, y) in merged:
            first_seen.setdefault((round(x, 1), round(y, 1)), it)

        old_std = hist['resid_std'][-1]
        new_std = np.nanstd(working_residual.image.array)
        all_peaks, prev_amps, amps_final = merged, amps, amps
        hist['resid_std'].append(new_std)
        hist['n_peaks'].append(len(merged))
        hist['n_zero'].append(fit['n_at_zero'])

        sizes = np.bincount(fit['group_labels'])
        print(f"  refit {len(merged)} peaks ({int(is_new.sum())} new, {fit['n_at_zero']} at zero amp) "
              f"-> residual std {old_std:.3f} -> {new_std:.3f}")
        print(f"  {fit['n_groups']} coupled groups; sizes: max={sizes.max()}, "
              f"median={np.median(sizes):.0f}, singletons={(sizes == 1).sum()}")
        print(f"  time: " + ", ".join(f"{k}={v:.2f}s" for k, v in times.items()))

        snap = dict(it=it, threshold=threshold, det_in=det_in, raw_peaks=raw_peaks,
                    collapsed_peaks=new_peaks, merged_peaks=merged, merged_arr=merged_arr,
                    is_new=is_new, prev_amp=prev_amp, footprint_mask=fp_mask, fit=fit,
                    residual=working_residual.image.array.copy(), std_before=old_std,
                    std_after=new_std, times=times, bad_mask=bad_mask,
                    shape_rejected=shape_rejected)
        snaps.append(snap)

        is_last_planned = (it == n_iterations - 1)
        should_plot = viz and (it % max(plot_every, 1) == 0 or is_last_planned)
        if should_plot:
            with timer.step('plotting', outer=outer_idx, inner=it, kind='viz'):
                pages = select_source_pages(snap, source_view, sources_per_page,
                                            max_pages_per_iter)
                colors = source_colors(pages)
                plot_iteration_overview(exposure, snap, hist, outer_idx, norm_data,
                                        colors=colors, view_half=view_half)
                plot_source_pages(exposure, snap, outer_idx, pages=pages, colors=colors,
                                  source_view=source_view, view_half=view_half)

        if old_std > 0 and it > 0:
            impr = abs(old_std - new_std) / old_std
            if impr < convergence_frac:
                print(f"  residual std improved by only {impr * 100:.2f}% -- converged, stopping")
                break
        threshold = max(threshold + threshold_step, threshold_min)

    return dict(model_full=model_full, residual_full=working_residual.image.array.copy(),
                peaks=all_peaks, amps=amps_final, resid_std_history=hist['resid_std'],
                n_peaks_history=hist['n_peaks'], n_zero_amp_history=hist['n_zero'],
                first_seen=first_seen, snaps=snaps)


# =====================================================================
# Cross-iteration summaries
# =====================================================================

def plot_amplitude_evolution(result, outer_idx=0):
    """Plot the amplitude of every final source against inner iteration.

    Parameters
    ----------
    result : `dict`
        Inner-loop result with ``snaps``.
    outer_idx : `int`, optional
        Outer-pass index (for the title).
    """
    snaps = result['snaps']
    if not snaps:
        return
    ref = snaps[-1]['merged_arr']
    amp_mat = np.full((len(snaps), len(ref)), np.nan)
    for i, s in enumerate(snaps):
        if len(s['merged_arr']) == 0:
            continue
        d, k = cKDTree(s['merged_arr']).query(ref)
        ok = d < 2.0
        amp_mat[i, ok] = s['fit']['amps'][k[ok]]

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(16, 5.5))
    final = np.nan_to_num(amp_mat[-1], nan=0.0)
    cols = plt.cm.viridis(np.argsort(np.argsort(final)) / max(len(final) - 1, 1))
    for k in range(len(ref)):
        ax0.plot(np.arange(1, len(snaps) + 1), amp_mat[:, k], 'o-', color=cols[k], ms=3, lw=1)
    ax0.set_yscale('symlog', linthresh=1.0)
    ax0.set_xlabel('inner iteration'); ax0.set_ylabel('fitted amplitude')
    ax0.set_title(f'Outer {outer_idx + 1}: amplitude of each final source vs iteration')

    with np.errstate(divide='ignore', invalid='ignore'):
        frac = (amp_mat[1:] - amp_mat[:-1]) / np.abs(amp_mat[:-1]) if len(snaps) > 1 else None
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


def plot_source_evolution(exposure, result, outer_idx=0, view_half=20):
    """Plot model and residual cross-sections of three sources at every iteration.

    The brightest, median and faintest final sources: the x cross-section of the model (top)
    and of the residual (bottom).

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    result : `dict`
        Inner-loop result with ``snaps``.
    outer_idx : `int`, optional
        Outer-pass index (for the title).
    view_half : `int`, optional
        Half-length of the cross-sections [px].
    """
    snaps = result['snaps']
    if not snaps:
        return
    bbox = exposure.getBBox()
    X0, Y0 = bbox.getMinX(), bbox.getMinY()
    data = exposure.image.array
    final = snaps[-1]
    peaks, amps = final['merged_arr'], final['fit']['amps']
    order = np.argsort(amps)
    picks = {'dimmest': order[0], 'median': order[len(order) // 2], 'brightest': order[-1]}
    cols = plt.cm.viridis(np.linspace(0, 1, len(snaps)))

    fig, axes = plt.subplots(2, 3, figsize=(18, 8), sharex='col')
    for c, (name, j) in enumerate(picks.items()):
        px, py = peaks[j]
        iy = int(np.clip(int(round(py)) - Y0, 0, data.shape[0] - 1))
        ix = int(np.clip(int(round(px)) - X0, 0, data.shape[1] - 1))
        xa, xb = max(ix - view_half, 0), min(ix + view_half + 1, data.shape[1])
        xx = X0 + np.arange(xa, xb)
        axes[0, c].plot(xx, data[iy, xa:xb], 'k', lw=2, label='data')
        for i, s in enumerate(snaps):
            axes[0, c].plot(xx, s['fit']['model_full'][iy, xa:xb], color=cols[i], lw=1.2,
                            label=f'iter {i + 1}')
            axes[1, c].plot(xx, data[iy, xa:xb] - s['fit']['model_full'][iy, xa:xb],
                            color=cols[i], lw=1.0)
        axes[1, c].axhline(0, color='gray', lw=0.6)
        for a in axes[:, c]:
            a.axvline(px, color='red', lw=0.7, alpha=0.6)
        axes[0, c].set_title(f'{name}: ({px:.1f},{py:.1f}) final amp={amps[j]:.1f}')
        axes[1, c].set_xlabel('x [px]')
    axes[0, 0].set_ylabel('data / model'); axes[1, 0].set_ylabel('residual')
    axes[0, 0].legend(fontsize=6, ncol=2)
    fig.suptitle(f'Outer {outer_idx + 1}: model & residual through the source row, all iterations',
                 fontsize=13)
    plt.tight_layout()
    plt.show()
    plt.close(fig)


# =====================================================================
# Outer-loop step figures
# =====================================================================

def plot_cutout_setup(cutout, psf, bad_bitmask, star_positions=None):
    """Plot the cutout before fitting: data, sigma map, mask planes and PSF.

    Parameters
    ----------
    cutout : `lsst.afw.image.Exposure`
        The cutout.
    psf : `lsst.afw.detection.Psf`
        PSF model.
    bad_bitmask : `int`
        Mask bits of bad pixels (shown in red).
    star_positions : array-like, optional
        Original PSF star positions to mark.
    """
    bbox = cutout.getBBox()
    ext = _extent(bbox)
    data = cutout.image.array
    bad = (cutout.mask.array & bad_bitmask) != 0
    fig, axes = plt.subplots(1, 4, figsize=(21, 5))

    im = axes[0].imshow(data, extent=ext, origin='lower', norm=_data_norm(data), cmap=CMAP)
    plt.colorbar(im, ax=axes[0], fraction=0.046, pad=0.04)
    _overlay_bad(axes[0], bad, ext)
    if star_positions is not None and len(star_positions):
        sp = _peaks_array(star_positions)
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
        _clean_axes(a)
    plt.tight_layout()
    plt.show()
    plt.close(fig)


def plot_outer_step_background(raw_exp, exp_next, src_mask, model_full, bg_next, bg_prev,
                               outer_idx):
    """Plot the background refit of one outer step.

    Input image, model and source mask, new and previous background, and the result.

    Parameters
    ----------
    raw_exp : `lsst.afw.image.Exposure`
        Input to the background refit.
    exp_next : `lsst.afw.image.Exposure`
        Background-subtracted result.
    src_mask : `numpy.ndarray` [`bool`]
        Source mask used.
    model_full : `numpy.ndarray`
        Joint-fit model.
    bg_next : `lsst.afw.image.Image` or `None`
        New background.
    bg_prev : `lsst.afw.image.Image` or `None`
        Previous background.
    outer_idx : `int`
        Outer-pass index.

    Returns
    -------
    norm : `astropy.visualization.ImageNormalize`
        Display normalisation used for the images.
    """
    bbox = raw_exp.getBBox()
    ext = _extent(bbox)
    xs, ys = _coords(bbox)
    before, after = raw_exp.image.array, exp_next.image.array
    fig, axes = plt.subplots(2, 4, figsize=(21, 9.5))

    def asinh(a):
        return simple_norm(a, stretch='asinh', min_percent=1, max_percent=99.5)

    im = axes[0, 0].imshow(before, extent=ext, origin='lower', cmap=CMAP, norm=asinh(before))
    plt.colorbar(im, ax=axes[0, 0], fraction=0.046, pad=0.04)
    axes[0, 0].set_title('Input to background refit')

    im = axes[0, 1].imshow(model_full, extent=ext, origin='lower', cmap=CMAP,
                           norm=simple_norm(model_full, stretch='asinh', percent=99.5))
    axes[0, 1].contour(xs, ys, src_mask.astype(float), levels=[0.5], colors='red', linewidths=0.8)
    plt.colorbar(im, ax=axes[0, 1], fraction=0.046, pad=0.04)
    axes[0, 1].set_title(f'Joint model + source mask (red)\n{100 * src_mask.mean():.1f}% masked')

    bgarr = bg_next.array if bg_next is not None else np.zeros_like(before)
    im = axes[0, 2].imshow(bgarr, extent=ext, origin='lower', cmap='viridis')
    plt.colorbar(im, ax=axes[0, 2], fraction=0.046, pad=0.04)
    axes[0, 2].set_title('Newly fit background (subtracted)')

    im = axes[0, 3].imshow(after, extent=ext, origin='lower', cmap=CMAP, norm=asinh(after))
    plt.colorbar(im, ax=axes[0, 3], fraction=0.046, pad=0.04)
    axes[0, 3].set_title('After background refit')

    cy, cx = before.shape[0] // 2, before.shape[1] // 2
    if bg_prev is not None:
        axes[1, 0].plot(xs, bg_prev[cy, :], label='previous bg')
        axes[1, 1].plot(ys, bg_prev[:, cx], label='previous bg')
    axes[1, 0].plot(xs, bgarr[cy, :], label='new bg')
    axes[1, 1].plot(ys, bgarr[:, cx], label='new bg')
    axes[1, 0].set_title(f'Background x cross-section (row {cy})'); axes[1, 0].legend()
    axes[1, 1].set_title(f'Background y cross-section (col {cx})'); axes[1, 1].legend()
    axes[1, 0].set_xlabel('x [px]'); axes[1, 1].set_xlabel('y [px]')

    axes[1, 2].plot(xs, np.where(src_mask[cy, :], np.nan, before[cy, :]), label='before (sources masked)')
    axes[1, 2].plot(xs, np.where(src_mask[cy, :], np.nan, after[cy, :]), label='after (sources masked)')
    axes[1, 2].axhline(0, color='gray', lw=0.6)
    axes[1, 2].set_title('Sky-only row before/after'); axes[1, 2].legend(fontsize=8)
    axes[1, 2].set_xlabel('x [px]')

    sky = ~src_mask
    lo, hi = np.nanpercentile(before[sky], [1, 99]) if sky.any() else (-1, 1)
    bins = np.linspace(lo, hi, 60)
    axes[1, 3].hist(before[sky], bins=bins, alpha=0.5, label=f'before (med={np.median(before[sky]):.2f})')
    axes[1, 3].hist(after[sky], bins=bins, alpha=0.5, label=f'after (med={np.median(after[sky]):.2f})')
    axes[1, 3].set_title('Sky pixel values'); axes[1, 3].legend(fontsize=8)
    for a in axes[0]:
        _clean_axes(a)
    fig.suptitle(f'Outer {outer_idx + 1} → {outer_idx + 2}: mask sources & refit background',
                 fontsize=13)
    plt.tight_layout()
    plt.show()
    plt.close(fig)


def plot_psf_candidate_stamps(cellSet, max_stamps=24, title=''):
    """Plot the PSF candidate stamps, used ones first.

    Parameters
    ----------
    cellSet : `lsst.afw.math.SpatialCellSet` or `None`
        PSF candidates.
    max_stamps : `int`, optional
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
    cands = sorted(cands, key=lambda c: c.getStatus() != GOOD_STATUS)[:max_stamps]
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


def plot_outer_step_psf(exp_next, sources, starCat, forced_positions, psf_prev, psf_new,
                        psf_updated, outer_idx):
    """Plot the PSF refit of one outer step: star selection and old vs new PSF.

    Parameters
    ----------
    exp_next : `lsst.afw.image.Exposure`
        Background-subtracted exposure.
    sources : `lsst.afw.table.SourceCatalog` or `None`
        Fresh detections.
    starCat : `lsst.afw.table.SourceCatalog` or `None`
        PSF candidates.
    forced_positions : array-like or `None`
        Forced original star positions.
    psf_prev, psf_new : `lsst.afw.detection.Psf`
        Previous and new PSF.
    psf_updated : `bool`
        Whether the PSF was updated.
    outer_idx : `int`
        Outer-pass index.
    """
    bbox = exp_next.getBBox()
    ext = _extent(bbox)
    data = exp_next.image.array
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
    if forced_positions is not None and len(forced_positions):
        fp = _peaks_array(forced_positions)
        axes[0, 0].scatter(fp[:, 0], fp[:, 1], marker='s', facecolors='none', edgecolors='lime',
                           s=110, label=f'forced original stars ({len(fp)})')
    axes[0, 0].set_xlim(ext[0], ext[1]); axes[0, 0].set_ylim(ext[2], ext[3])
    axes[0, 0].legend(fontsize=8, loc='upper right')
    axes[0, 0].set_title('Star selection on refit-background image')
    _clean_axes(axes[0, 0])

    kp = psf_prev.computeKernelImage(c).array
    kn = psf_new.computeKernelImage(c).array
    for a, k, lab, psf_o in ((axes[0, 1], kp, 'previous', psf_prev), (axes[0, 2], kn, 'new', psf_new)):
        im = a.imshow(k, origin='lower', cmap='inferno')
        plt.colorbar(im, ax=a, fraction=0.046, pad=0.04)
        a.set_title(f'{lab} PSF  FWHM={_fwhm(psf_o, c):.3f}px  ({k.shape[0]}x{k.shape[1]})')
        _clean_axes(a)
    if not psf_updated:
        axes[0, 2].set_title(axes[0, 2].get_title() + '\n(NOT updated -- same as previous)', color='red')

    if kp.shape == kn.shape:
        d = kn - kp
        s = np.nanmax(np.abs(d)) or 1.0
        im = axes[1, 0].imshow(d, origin='lower', cmap='RdBu_r', vmin=-s, vmax=s)
        plt.colorbar(im, ax=axes[1, 0], fraction=0.046, pad=0.04)
        axes[1, 0].set_title('new - previous')
        _clean_axes(axes[1, 0])
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
    fig.suptitle(f'Outer {outer_idx + 1} → {outer_idx + 2}: PSF redetermination', fontsize=13)
    plt.tight_layout()
    plt.show()
    plt.close(fig)


# =====================================================================
# PSF refinement from the bright stars' residuals (+ noise diagnostic)
# =====================================================================
# After a run:
#     highfreq_noise_check(outer_history)                  # is variance too small?
#     psf_new, info = refine_psf_from_fit(outer_history)   # one update step
#     plot_psf_correction(info)
# or let the pipeline do it each outer pass with refit_psf='stack'.

def highfreq_noise_check_arrays(resid, variance, model, good, nbins=6, verbose=True):
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
              f"  -> use var + k*model in effective_variance if ratio grows with flux")
    return dict(flux=mid, ratio=ratio, k=k)


def highfreq_noise_check(outer_history, outer_idx=-1, snap_idx=-1, **kw):
    """Run `highfreq_noise_check_arrays` on one snapshot of a run.

    Parameters
    ----------
    outer_history : `list` [`dict`]
        Outer-loop history.
    outer_idx : `int`, optional
        Outer pass (default: last).
    snap_idx : `int`, optional
        Inner-loop snapshot (default: last).
    **kw
        Keyword arguments of `highfreq_noise_check_arrays`.

    Returns
    -------
    result : `dict`
        See `highfreq_noise_check_arrays`.
    """
    h = outer_history[outer_idx]
    snap = h['result']['snaps'][snap_idx]
    fit = snap['fit']
    resid = h['exposure'].image.array - fit['model_full']
    return highfreq_noise_check_arrays(resid, fit['variance'], fit['model_full'], fit['good'], **kw)


def estimate_psf_correction(data, model, variance, good, xy0, positions, amps, kernel,
                            n_use=8, step=0.8, clip=3.0, soft_sigma=1.0,
                            window=(0.6, 1.0), max_radius=None):
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
    n_use : `int`, optional
        Number of brightest stars used.
    step : `float`, optional
        Fraction of the correction applied.
    clip : `float`, optional
        Sigma clipping of the stack.
    soft_sigma : `float`, optional
        Soft threshold of the correction [sigma].
    window : `tuple` [`float`, `float`], optional
        Inner and outer edge of the window as fractions of the kernel half-size.
    max_radius : `float`, optional
        Do not change the kernel beyond this radius [px].

    Returns
    -------
    new : `numpy.ndarray`
        Corrected, normalised kernel.
    info : `dict`
        ``old``, ``new``, ``delta``, ``raw_delta``, ``err``, ``n_stars``, ``used``.
    """
    positions = np.asarray(positions, float)
    amps = np.asarray(amps, float)
    ny, nx = data.shape
    resid = data - model
    good_f = good.astype(float)
    h = kernel.shape[0] // 2
    u = np.arange(-h, h + 1)
    V, U = np.meshgrid(u, u, indexing='ij')

    order = np.argsort(-amps)[:n_use]
    R, W = [], []
    for j in order:
        if amps[j] <= 0:
            continue
        xj, yj = positions[j, 0] - xy0[0], positions[j, 1] - xy0[1]
        yy, xx = yj + V, xj + U
        coords = [yy.ravel(), xx.ravel()]
        inside = (yy >= 0) & (yy <= ny - 1) & (xx >= 0) & (xx <= nx - 1)
        r = map_coordinates(resid, coords, order=3, mode='nearest').reshape(V.shape) / amps[j]
        gm = map_coordinates(good_f, coords, order=1, mode='nearest').reshape(V.shape) > 0.99
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
    raw_delta = delta.copy()
    delta = np.sign(delta) * np.maximum(np.abs(delta) - soft_sigma * err, 0.0)

    # max_radius: the correction is tapered to zero at this radius (px) instead of at the
    # kernel edge -- use it with an extended (Moffat-tailed) kernel so the noisy outer
    # residuals never modify the analytic tail, e.g. stack_kwargs=dict(max_radius=20).
    he = h if max_radius is None else min(h, max_radius)
    win = window if window else (0.0, 1.0)
    rr = np.hypot(U, V)
    w = np.clip((he - rr) / max(he * (win[1] - win[0]), 1e-9), 0, 1)
    if window or max_radius is not None:
        delta = delta * (0.5 - 0.5 * np.cos(np.pi * w))

    gy, gx = np.gradient(kernel)
    A = np.column_stack([gx.ravel(), gy.ravel()])
    coef, *_ = np.linalg.lstsq(A, delta.ravel(), rcond=None)
    delta = delta - (A @ coef).reshape(delta.shape)

    new = kernel + step * delta
    new = new / new.sum()
    info = dict(old=kernel, new=new, delta=delta, raw_delta=raw_delta, err=err,
                n_stars=len(R), used=[int(j) for j in order])
    return new, info


def kernel_to_lsst_psf(kernel):
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


def extend_kernel_with_moffat(kernel, new_size=81, fit_r=(6.0, 14.0), blend=(12.0, 16.0)):
    """Replace the noisy outer part of a PSF kernel by an analytic Moffat tail.

    The azimuthal-mean profile between ``fit_r[0]`` and ``fit_r[1]`` px is fitted with
    ``A * (1 + (r/a)^2)^-b`` (relative weighting, so faint outer bins count).  The new kernel
    equals the original inside ``blend[0]``, the Moffat beyond ``blend[1]``, with a linear
    cross-fade between, and is renormalised -- fitted amplitudes then become total fluxes
    including the tail.  Works on an already-extended kernel too.

    Parameters
    ----------
    kernel : `numpy.ndarray`
        PSF kernel.
    new_size : `int`, optional
        Size of the new kernel [px] (odd).
    fit_r : `tuple` [`float`, `float`], optional
        Radial range of the Moffat fit [px].
    blend : `tuple` [`float`, `float`], optional
        Radial range of the cross-fade [px].

    Returns
    -------
    new_kernel : `numpy.ndarray`
        Extended, normalised kernel.
    info : `dict`
        Moffat ``A``, ``a``, ``beta`` and ``tail_frac`` (flux fraction beyond ``blend[1]``).

    Raises
    ------
    RuntimeError
        Raised if the Moffat fit fails.
    """
    n = kernel.shape[0]
    h = n // 2
    u = np.arange(-h, h + 1)
    r = np.hypot(*np.meshgrid(u, u))
    r_hi = min(fit_r[1], h - 1.0)
    edges = np.arange(fit_r[0], r_hi + 1.0)
    rc = 0.5 * (edges[1:] + edges[:-1])
    prof = np.array([kernel[(r >= lo) & (r < hi)].mean() for lo, hi in zip(edges[:-1], edges[1:])])
    moffat = lambda rr, A, a, b: A * (1.0 + (rr / a) ** 2) ** (-b)
    sigma = np.clip(prof, 1e-6 * np.nanmax(np.abs(prof)), None)      # relative weighting
    (A, a, b), _ = curve_fit(moffat, rc, prof, p0=[max(prof[0], 1e-12), 4.0, 3.0], sigma=sigma,
                             bounds=([0.0, 0.5, 1.05], [np.inf, 60.0, 12.0]), maxfev=5000)

    H = new_size // 2
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
    info = dict(A=float(A), a=float(a), beta=float(b), tail_frac=float(out[R > b1].sum()))
    return out, info


def extend_psf_wings(psf, point, size=81, verbose=True, **kw):
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
        Keyword arguments of `extend_kernel_with_moffat`.

    Returns
    -------
    psf : `lsst.afw.detection.Psf`
        The extended PSF, or the input PSF if the fit fails.
    """
    kernel = psf.computeKernelImage(point).array.astype(float).copy()
    try:
        new, info = extend_kernel_with_moffat(kernel, new_size=size, **kw)
    except Exception as e:
        print(f"  PSF wing extension failed ({type(e).__name__}: {e}) -- keeping current PSF")
        return psf
    if verbose:
        print(f"  PSF wings extended: {kernel.shape[0]}px -> {size}px kernel; Moffat "
              f"a={info['a']:.2f}px, beta={info['beta']:.2f}; "
              f"{100 * info['tail_frac']:.2f}% of flux in the tail beyond the blend radius")
    return kernel_to_lsst_psf(new)


def refine_psf_from_fit(outer_history, outer_idx=-1, snap_idx=-1, **kw):
    """Make one PSF update step from a finished run.

    Parameters
    ----------
    outer_history : `list` [`dict`]
        Outer-loop history.
    outer_idx : `int`, optional
        Outer pass (default: last).
    snap_idx : `int`, optional
        Inner-loop snapshot (default: last).
    **kw
        Keyword arguments of `estimate_psf_correction`.

    Returns
    -------
    psf : `lsst.meas.algorithms.KernelPsf`
        Updated PSF.
    info : `dict`
        See `estimate_psf_correction`.
    """
    h = outer_history[outer_idx]
    exp, psf = h['exposure'], h['psf']
    snap = h['result']['snaps'][snap_idx]
    fit = snap['fit']
    kernel = psf.computeKernelImage(exp.getBBox().getCenter()).array.copy()
    xy0 = (exp.getXY0().getX(), exp.getXY0().getY())
    new, info = estimate_psf_correction(exp.image.array, fit['model_full'], fit['variance'],
                                        fit['good'], xy0, snap['merged_arr'], fit['amps'],
                                        kernel, **kw)
    return kernel_to_lsst_psf(new), info


def polish_psf_fixed_sources(exposure, psf, result, bad_bitmask, n_cycles=3,
                             stack_kwargs=None, support_frac=0.01, sys_frac=0.03,
                             poisson_k=0.0, refine_centroids=True, centroid_min_amp=50.0,
                             freeze_sep=None, max_shift=0.3, verbose=True):
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
    bad_bitmask : `int`
        Mask bits of pixels to ignore.
    n_cycles : `int`, optional
        Maximum number of cycles.
    stack_kwargs : `dict`, optional
        Keyword arguments of `estimate_psf_correction`.
    support_frac, sys_frac, poisson_k : `float`, optional
        See `joint_psf_fit_diag`.
    refine_centroids : `bool`, optional
        Refine positions after each PSF change.
    centroid_min_amp : `float`, optional
        Fainter sources keep their position.
    freeze_sep : `float`, optional
        See `refine_peak_positions`.
    max_shift : `float`, optional
        Maximum position shift [px].
    verbose : `bool`, optional
        Print each cycle.

    Returns
    -------
    psf : `lsst.afw.detection.Psf`
        The best PSF.
    info : `dict` or `None`
        Information of the last kept correction (see `estimate_psf_correction`).
    """
    snap = result['snaps'][-1]
    peaks = [tuple(p) for p in snap['merged_arr']]
    fit = snap['fit']
    data = exposure.image.array

    def rms(f):
        return float(np.sqrt(np.mean((data - f['model_full'])[f['good']] ** 2)))

    cur, best, info = psf, rms(fit), None
    for c in range(n_cycles):
        hist = [dict(exposure=exposure, psf=cur,
                     result=dict(snaps=[dict(fit=fit, merged_arr=_peaks_array(peaks))]))]
        try:
            cand, cinfo = refine_psf_from_fit(hist, -1, -1, **(stack_kwargs or {}))
        except ValueError as e:
            print(f"    polish cycle {c + 1}: {e} -- stopping")
            break
        info = cinfo
        fit_c = joint_psf_fit_diag(exposure, cand, peaks, bad_bitmask,
                                   support_frac=support_frac, var_model=fit['model_full'],
                                   sys_frac=sys_frac, poisson_k=poisson_k)
        peaks_c = peaks
        if refine_centroids:
            peaks_c = refine_peak_positions(
                exposure, cand, peaks, fit_c['amps'], bad_bitmask,
                min_amp_for_refine=centroid_min_amp, max_shift=max_shift,
                variance=fit_c['variance'], freeze_sep=freeze_sep)
            fit_c = joint_psf_fit_diag(exposure, cand, peaks_c, bad_bitmask,
                                       support_frac=support_frac,
                                       var_model=fit_c['model_full'],
                                       sys_frac=sys_frac, poisson_k=poisson_k)
        r = rms(fit_c)
        ok = r < best
        if verbose:
            print(f"    polish cycle {c + 1}/{n_cycles}: residual rms {best:.3f} -> {r:.3f} "
                  f"({'kept' if ok else 'REJECTED, stopping'})")
        if not ok:
            break
        cur, fit, peaks, best = cand, fit_c, peaks_c, r
    return cur, info


def plot_psf_correction(info):
    """Plot the old kernel, the proposed correction and the radial profiles.

    Parameters
    ----------
    info : `dict`
        Output of `estimate_psf_correction`.
    """
    old, new, d, raw, err = info['old'], info['new'], info['delta'], info['raw_delta'], info['err']
    h = old.shape[0] // 2
    ext = (-h - .5, h + .5, -h - .5, h + .5)
    fig, ax = plt.subplots(1, 4, figsize=(20, 4.6))
    im = ax[0].imshow(old, origin='lower', extent=ext, cmap='inferno'); plt.colorbar(im, ax=ax[0], fraction=.046)
    ax[0].set_title('current kernel')
    s = np.nanmax(np.abs(raw[np.isfinite(raw)])) or 1
    im = ax[1].imshow(raw, origin='lower', extent=ext, cmap='RdBu_r', vmin=-s, vmax=s); plt.colorbar(im, ax=ax[1], fraction=.046)
    ax[1].set_title(f"stacked (data-model)/amp, {info['n_stars']} stars\n(true PSF - current PSF)")
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

def _wants_more_outer(result, outer_idx):
    """Return False: never extend the outer loop (replaced by ``adaptive.install()``)."""
    return False


def run_cutout_pipeline(exposure_crop, psf, bad_bitmask, starCat=None, bgImage=None,
                        x_range=(0, 100), y_range=(0, 100),
                        n_outer_iterations=2, n_inner_iterations=10, inner_kwargs=None,
                        outer_convergence_frac=0.02, refit_psf=True, psf_kernel_size=41,
                        psf_kernel_sizes=None,
                        bg_bin_size=32, bg_use_approx=True,
                        viz=True, plot_every=1, source_view='all', sources_per_page=6,
                        max_pages_per_iter=None, view_half=20, support_frac=0.01,
                        refine_centroids=True, centroid_min_amp=50.0,
                        n_centroid_rounds=4, centroid_shift_tol=0.03,
                        centroid_shift_schedule=None, coarse_centroid=True,
                        coarse_centroid_hw_sigma=2.5, merge_separation_sigma=1.5,
                        sys_frac=0.03, min_snr=5.0, prune_sep_sigma=1.5,
                        reject_fwhm=0.75, bright_peak_sigma=100.0, freeze_fwhm=1.0, poisson_k=0.0,
                        near_bright_fwhm=2.0, near_bright_snr=10.0,
                        edge_margin_fwhm=0.5, use_eliminate_weak=True,
                        use_shape_cut=True, shape_kwargs=None, shape_in_first_pass=False,
                        stack_kwargs=None, psf_polish_cycles=3, psf_wing_size=81):
    """Run the outer loop (inner fit, background and PSF refinement) on a cutout.

    Parameters
    ----------
    exposure_crop : `lsst.afw.image.Exposure`
        Region to process (before background subtraction).
    psf : `lsst.afw.detection.Psf`
        Initial PSF.
    bad_bitmask : `int`
        Mask bits of pixels to ignore.
    starCat : `lsst.afw.table.SourceCatalog`, optional
        Original PSF stars, forced into every PSF determination.
    bgImage : `lsst.afw.image.Image`, optional
        Detector background subtracted before the first pass.
    x_range, y_range : `tuple` [`int`, `int`], optional
        Sub-region of ``exposure_crop`` (see `make_cutout`).
    n_outer_iterations : `int`, optional
        Number of outer passes.
    n_inner_iterations : `int`, optional
        Maximum inner iterations per pass.
    inner_kwargs : `dict`, optional
        Keyword arguments of `run_inner_loop_verbose`.
    outer_convergence_frac : `float`, optional
        Stop when the residual standard deviation improves by less than this fraction.
    refit_psf : `bool` or `str`, optional
        False, True (PCA star-selector refit) or 'stack' (residual-stack refit from the bright
        stars, see `refine_psf_from_fit`).
    psf_kernel_size : `int`, optional
        PSF kernel size [px].
    psf_kernel_sizes : `tuple` [`int`], optional
        Kernel sizes tried by `redetect_and_refit_psf_adaptive`.
    bg_bin_size : `int`, optional
        Background bin size [px].
    bg_use_approx : `bool`, optional
        See `refit_background`.
    viz : `bool`, optional
        Draw diagnostics.
    plot_every : `int`, optional
        Draw every n-th inner iteration (the last always).
    source_view : `str`, optional
        'all', 'new' or 'changed': which sources get a per-source page.
    sources_per_page : `int`, optional
        Sources per page.
    max_pages_per_iter : `int`, optional
        Maximum per-source pages per iteration (None = no limit).
    view_half : `int`, optional
        Half-size of the per-source windows [px].
    support_frac : `float`, optional
        See `joint_psf_fit_diag`.
    refine_centroids : `bool`, optional
        Refine positions (`refine_peak_positions`).
    centroid_min_amp : `float`, optional
        Fainter sources keep their position.
    n_centroid_rounds, centroid_shift_tol, centroid_shift_schedule : optional
        Centroid refinement controls (see `run_inner_loop_verbose`).
    coarse_centroid, coarse_centroid_hw_sigma : optional
        Windowed starting centroid of new peaks (see `run_inner_loop_verbose`).
    merge_separation_sigma : `float`, optional
        Merge new peaks closer than this many PSF sigma.
    sys_frac, min_snr, prune_sep_sigma, reject_fwhm, bright_peak_sigma : `float`, optional
        Over-counting guards (see `run_inner_loop_verbose`).
    freeze_fwhm, poisson_k, near_bright_fwhm, near_bright_snr : `float`, optional
        Over-counting guards (see `run_inner_loop_verbose`).
    edge_margin_fwhm : `float`, optional
        Drop sources within this many FWHM of the cutout edge.
    use_eliminate_weak : `bool`, optional
        Backward elimination on the honest S/N.
    use_shape_cut : `bool`, optional
        PSF-shape constraint (`eliminate_nonpsf`).
    shape_kwargs : `dict`, optional
        Keyword arguments of `shape_metrics` / `judge_shape`.
    shape_in_first_pass : `bool`, optional
        Apply the shape cut in the first outer pass too (by default only once the PSF has been
        refined: a rough first PSF can make the cut reject real faint companions).
    stack_kwargs : `dict`, optional
        Keyword arguments of `refine_psf_from_fit` (``refit_psf='stack'``).
    psf_polish_cycles : `int`, optional
        PSF correction / re-solve cycles per outer pass (``refit_psf='stack'``).
    psf_wing_size : `int`, optional
        Convert every PSF to a constant kernel of this size with a Moffat tail (see
        `extend_psf_wings`); None = off.  Stops the model ending at the stamp edge, which leaves
        a ring of fake sources half a kernel width from bright stars.

    Returns
    -------
    outer_history : `list` [`dict`]
        One dict per outer pass: ``psf``, ``exposure``, ``result``, ``resid_std``, ``n_peaks``,
        ``bg_image``.
    """
    TIMER.reset()
    if psf_kernel_sizes is None:
        fallbacks = [k for k in (27, 21, 15, 11) if k < psf_kernel_size]
        psf_kernel_sizes = tuple([psf_kernel_size] + fallbacks)
    inner_kwargs = dict(n_iterations=n_inner_iterations, threshold_start=8.0,
                        threshold_step=-1.5, threshold_min=4.5, min_pixels=1,
                        min_new_peaks=1, convergence_frac=0.02, **(inner_kwargs or {})) \
        if inner_kwargs is None or 'n_iterations' not in inner_kwargs else dict(inner_kwargs)

    with TIMER.step('make cutout'):
        cutout = make_cutout(exposure_crop, x_range, y_range)
        cutout.setPsf(psf)
    cb = cutout.getBBox()
    if psf_wing_size:
        psf = extend_psf_wings(psf, cb.getCenter(), psf_wing_size)
        cutout.setPsf(psf)
    print(f"Cutout parent-frame bbox: x[{cb.getMinX()}, {cb.getMaxX()}] "
          f"y[{cb.getMinY()}, {cb.getMaxY()}]")

    original_star_positions = None
    if starCat is not None and len(starCat) > 0:
        allpos = [(r.getCentroid().getX(), r.getCentroid().getY()) for r in starCat]
        original_star_positions = [(x, y) for x, y in allpos
                                   if cb.contains(geom.Point2I(int(round(x)), int(round(y))))]
        print(f"{len(original_star_positions)} of {len(allpos)} original PSF stars fall inside the cutout")

    if viz:
        plot_cutout_setup(cutout, psf, bad_bitmask, original_star_positions)

    bg_prev = None
    if bgImage is not None:
        try:
            bg_prev = bgImage[cb].array.copy()
        except Exception as e:
            print(f"(could not crop initial bgImage: {e})")

    raw = cutout.clone()
    exposure_working = cutout.clone()
    exposure_working.setPsf(psf)
    psf_working = psf
    bg_image_working = None
    outer_history = []

    for outer_it in __import__("itertools").count():   # open-ended; see break below
        print(f"\n{'=' * 70}\n===== Outer pass {outer_it + 1}/{n_outer_iterations} "
              f"(PSF FWHM={_fwhm(psf_working, cb.getCenter()):.3f}px) =====\n{'=' * 70}")
        shape_now = use_shape_cut and (outer_it > 0 or shape_in_first_pass)
        print(f"  PSF-shape cut this pass: {'ON' if shape_now else 'off'}")
        with TIMER.step('INNER LOOP (total, incl. plotting)', outer=outer_it, kind='wall'):
            result = run_inner_loop_verbose(
                exposure_working, psf_working, bad_bitmask, outer_idx=outer_it,
                viz=viz, plot_every=plot_every, source_view=source_view,
                sources_per_page=sources_per_page, max_pages_per_iter=max_pages_per_iter,
                view_half=view_half, support_frac=support_frac,
                refine_centroids=refine_centroids, centroid_min_amp=centroid_min_amp,
                n_centroid_rounds=n_centroid_rounds, centroid_shift_tol=centroid_shift_tol,
                centroid_shift_schedule=centroid_shift_schedule,
                coarse_centroid=coarse_centroid,
                coarse_centroid_hw_sigma=coarse_centroid_hw_sigma,
                merge_separation_sigma=merge_separation_sigma,
                sys_frac=sys_frac, min_snr=min_snr, prune_sep_sigma=prune_sep_sigma,
                reject_fwhm=reject_fwhm, bright_peak_sigma=bright_peak_sigma,
                freeze_fwhm=freeze_fwhm, poisson_k=poisson_k,
                near_bright_fwhm=near_bright_fwhm, near_bright_snr=near_bright_snr,
                edge_margin_fwhm=edge_margin_fwhm, use_eliminate_weak=use_eliminate_weak,
                use_shape_cut=shape_now, shape_kwargs=shape_kwargs,
                **inner_kwargs)
        final_std = result['resid_std_history'][-1]
        outer_history.append(dict(psf=psf_working, exposure=exposure_working, result=result,
                                  resid_std=final_std, n_peaks=len(result['peaks']),
                                  bg_image=bg_image_working))
        print(f"\nouter pass {outer_it + 1}: {len(result['peaks'])} peaks, residual std={final_std:.4f}")

        if viz:
            with TIMER.step('plotting', outer=outer_it, kind='viz'):
                plot_amplitude_evolution(result, outer_it)
                plot_source_evolution(exposure_working, result, outer_it, view_half=view_half)

        if outer_it > 0:
            prev = outer_history[-2]['resid_std']
            impr = abs(prev - final_std) / prev if prev > 0 else 0
            if impr < outer_convergence_frac:
                print(f"outer residual std improved by only {impr * 100:.2f}% -- outer loop converged")
                break
        if outer_it >= n_outer_iterations - 1 and not _wants_more_outer(result, outer_it):
            break

        # ---- re-derive background + PSF ------------------------------------
        exposure_next = raw.clone()
        with TIMER.step('mask sources from model', outer=outer_it):
            src_mask = mask_sources_from_model(exposure_next, result['model_full'])
        bg_next = None
        with TIMER.step('refit background', outer=outer_it):
            try:
                exposure_next, bg_next = refit_background(
                    exposure_next, use_approx=bg_use_approx, bin_size=bg_bin_size)
            except Exception as e:
                print(f"  background refit failed ({type(e).__name__}: {e}) -- "
                      f"keeping un-refit image")
        if viz:
            with TIMER.step('plotting', outer=outer_it, kind='viz'):
                plot_outer_step_background(raw, exposure_next, src_mask, result['model_full'],
                                           bg_next, bg_prev, outer_it)
        bg_prev = bg_next.array.copy() if bg_next is not None else bg_prev

        psf_new, sources_new, starCat_new, cellSet_new = psf_working, None, None, None
        exposure_next.setPsf(psf_working)
        if refit_psf == 'stack':
            with TIMER.step('PSF refit (residual stack)', outer=outer_it):
                try:
                    if psf_polish_cycles > 1:
                        psf_new, pinfo = polish_psf_fixed_sources(
                            exposure_working, psf_working, result, bad_bitmask,
                            n_cycles=psf_polish_cycles, stack_kwargs=stack_kwargs,
                            support_frac=support_frac, sys_frac=sys_frac,
                            poisson_k=poisson_k, refine_centroids=refine_centroids,
                            centroid_min_amp=centroid_min_amp,
                            freeze_sep=freeze_fwhm * _fwhm(psf_working, cb.getCenter()))
                    else:
                        psf_new, pinfo = refine_psf_from_fit(outer_history, -1,
                                                             **(stack_kwargs or {}))
                    print(f"  PSF updated from residual stack of {pinfo['n_stars']} stars: "
                          f"FWHM {_fwhm(psf_working, cb.getCenter()):.3f} -> "
                          f"{_fwhm(psf_new, cb.getCenter()):.3f} px")
                    if viz:
                        plot_psf_correction(pinfo)
                except Exception as e:
                    print(f"  residual-stack PSF refit raised {type(e).__name__}: {e} -- "
                          f"keeping current PSF")
                    psf_new = psf_working
            exposure_next.setPsf(psf_new)
        elif refit_psf:
            with TIMER.step('redetect + PSF refit', outer=outer_it):
                try:
                    psf_new, sources_new, starCat_new, cellSet_new, ks_used = \
                        redetect_and_refit_psf_adaptive(
                            exposure_next, kernel_sizes=psf_kernel_sizes,
                            extra_peak_positions=result['peaks'],
                            force_include_star_positions=original_star_positions)
                except Exception as e:
                    print(f"  PSF redetermination raised {type(e).__name__}: {e} -- "
                          f"keeping current PSF")
                    cellSet_new = None
            psf_updated = cellSet_new is not None
            if not psf_updated:
                psf_new = psf_working
            exposure_next.setPsf(psf_new)
            if viz:
                with TIMER.step('plotting', outer=outer_it, kind='viz'):
                    plot_outer_step_psf(exposure_next, sources_new, starCat_new,
                                        original_star_positions, psf_working, psf_new,
                                        psf_updated, outer_it)
                    plot_psf_candidate_stamps(cellSet_new, title=f'(outer pass {outer_it + 2})')

        if psf_wing_size and psf_new is not psf_working:
            # PSF changed this pass: refresh the Moffat tail so it matches the new core
            psf_new = extend_psf_wings(psf_new, cb.getCenter(), psf_wing_size)
            exposure_next.setPsf(psf_new)

        exposure_working, psf_working, bg_image_working = exposure_next, psf_new, bg_next

    plot_timing_summary(TIMER)
    return outer_history


def plot_wing_diagnostic(outer_history, n_stars=3, rmax=40):
    """Plot the azimuthal profiles of the brightest sources: data vs their own model.

    Data beyond the dashed line well above the model means a real halo the kernel cannot
    represent; a difference centred on zero without a step means only core / PSF-shape mismatch
    remains.

    Parameters
    ----------
    outer_history : `list` [`dict`]
        Outer-loop history.
    n_stars : `int`, optional
        Number of brightest sources shown.
    rmax : `float`, optional
        Maximum radius [px].
    """
    final = outer_history[-1]
    exp, res = final['exposure'], final['result']
    snap = res['snaps'][-1]
    fit, peaks = snap['fit'], snap['merged_arr']
    bbox = exp.getBBox()
    ny, nx = bbox.getHeight(), bbox.getWidth()
    data, model, good = exp.image.array, fit['model_full'], fit['good']
    yy, xx = np.mgrid[bbox.getMinY():bbox.getMaxY() + 1, bbox.getMinX():bbox.getMaxX() + 1]
    edges = np.arange(0, rmax + 1)
    rc = 0.5 * (edges[1:] + edges[:-1])

    picks = np.argsort(-fit['amps'])[:n_stars]
    fig, axes = plt.subplots(1, len(picks), figsize=(6.2 * len(picks), 5), squeeze=False)
    for ax, j in zip(axes[0], picks):
        c = fit['comps'][j]
        if c is None:
            continue
        own = np.zeros((ny, nx))
        own[c['ys'], c['xs']] = c['amp'] * c['stamp']
        others_removed = data - (model - own)
        r = np.hypot(xx - peaks[j][0], yy - peaks[j][1])
        d_prof, m_prof = [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            sel = good & (r >= lo) & (r < hi)
            ok = sel.sum() > 3
            d_prof.append(np.mean(others_removed[sel]) if ok else np.nan)
            m_prof.append(np.mean(own[sel]) if ok else np.nan)
        d_prof, m_prof = np.array(d_prof), np.array(m_prof)
        ax.plot(rc, d_prof, 'k', lw=1.5, label='data - neighbours')
        ax.plot(rc, m_prof, '--', color='#DD8452', lw=1.5, label='own model')
        ax.plot(rc, d_prof - m_prof, color='gray', lw=1.2, label='difference')
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


def _select_sources(peaks, amps, source_ids, top_n, min_amp):
    """Return the indices of the sources to plot (given ids, or the brightest above a flux)."""
    if source_ids is not None:
        return [int(j) for j in source_ids]
    order = np.argsort(-amps)
    order = [int(j) for j in order if amps[j] > min_amp]
    return order[:top_n] if top_n is not None else order


def plot_source_cutouts(outer_history, outer_idx=-1, snap_idx=-1,
                        source_ids=None, top_n=None, min_amp=0.0,
                        half=15, cross_axis='x', show_own=True,
                        sources_per_fig=4, pdf_path=None, show=True,
                        cmap=None, mark_neighbours=True):
    """Plot cutouts and cross-sections of individual sources of a run.

    Parameters
    ----------
    outer_history : `list` [`dict`]
        Outer-loop history.
    outer_idx : `int`, optional
        Outer pass (default: last).
    snap_idx : `int`, optional
        Inner-loop snapshot (default: last).
    source_ids : `list` [`int`], optional
        Sources to show; overrides ``top_n`` / ``min_amp``.
    top_n : `int`, optional
        Show the brightest N sources.
    min_amp : `float`, optional
        Skip sources with amplitude at or below this.
    half : `int`, optional
        Cutout half-width [px].
    cross_axis : `str`, optional
        'x', 'y' or 'both'.
    show_own : `bool`, optional
        Draw the source's own PSF component in the cross-section.
    sources_per_fig : `int`, optional
        Sources per figure.
    pdf_path : `str`, optional
        Also save every figure as a page of this PDF.
    show : `bool`, optional
        Display the figures.
    cmap : `str`, optional
        Colour map of the images.
    mark_neighbours : `bool`, optional
        Mark the neighbouring sources.

    Returns
    -------
    figs : `list` [`matplotlib.figure.Figure`]
        The figures.
    """
    cmap = cmap or CMAP
    h = outer_history[outer_idx]
    exposure = h['exposure']
    snap = h['result']['snaps'][snap_idx]
    fit = snap['fit']
    data = exposure.image.array
    model = fit['model_full']
    resid = data - model
    peaks = np.asarray(snap['merged_arr'], float)
    amps = fit['amps']
    bad = snap['bad_mask']
    bbox = exposure.getBBox()
    X0, Y0 = bbox.getMinX(), bbox.getMinY()
    ny, nx = data.shape

    ids = _select_sources(peaks, amps, source_ids, top_n, min_amp)
    if not ids:
        print("No sources selected.")
        return []

    axes_list = {'x': ['x'], 'y': ['y'], 'both': ['x', 'y']}[cross_axis]
    ncol = 3 + len(axes_list)

    pdf = None
    if pdf_path is not None:
        from matplotlib.backends.backend_pdf import PdfPages
        pdf = PdfPages(pdf_path)

    figs = []
    pages = [ids[i:i + sources_per_fig] for i in range(0, len(ids), sources_per_fig)]
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

            iy_c = int(np.clip(iy, 0, ny - 1))
            ix_c = int(np.clip(ix, 0, nx - 1))

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
                if mark_neighbours:
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

            for k, axis in enumerate(axes_list):
                ax = axes[r, 3 + k]
                if axis == 'x':
                    coord = X0 + np.arange(xa, xb)
                    sl, centre = (iy_c, slice(xa, xb)), px
                else:
                    coord = Y0 + np.arange(ya, yb)
                    sl, centre = (slice(ya, yb), ix_c), py
                ax.plot(coord, data[sl], 'k-', lw=1.6, label='data')
                ax.plot(coord, model[sl], '-', color='#DD8452', lw=1.6, alpha=0.9, label='model')
                if show_own:
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

        fig.suptitle(f"Outer {outer_idx % len(outer_history) + 1} "
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
        print(f"saved {len(pages)} page(s) to {pdf_path}")
    return figs


# =====================================================================
# False-positive diagnostics
# =====================================================================

def negative_source_check(outer_history, outer_idx=-1, snap_idx=-1, thresholds=(3, 4, 5),
                          plot=True):
    """Estimate the number of false detections from the residual itself.

    Counts pixels of the matched-filter S/N map of the final residual at >= +N sigma and at
    <= -N sigma.  Real sources are positive only, so the negative excursions approximate the
    positive noise excursions (false positives); PSF-mismatch dipoles inflate both.

    Parameters
    ----------
    outer_history : `list` [`dict`]
        Outer-loop history.
    outer_idx : `int`, optional
        Outer pass (default: last).
    snap_idx : `int`, optional
        Inner-loop snapshot (default: last).
    thresholds : `tuple` [`float`], optional
        S/N thresholds.
    plot : `bool`, optional
        Plot the S/N map and histogram.

    Returns
    -------
    snr_map : `numpy.ndarray`
        Matched-filter S/N map of the residual.
    """
    h = outer_history[outer_idx]
    exp = h['exposure']
    fit = h['result']['snaps'][snap_idx]['fit']
    k = h['psf'].computeKernelImage(exp.getBBox().getCenter()).array
    resid = exp.image.array - fit['model_full']
    var = fit['variance']
    good = fit['good']
    inv = np.where(good, 1.0 / var, 0.0)
    num = fftconvolve(np.where(good, resid, 0.0) * inv, k[::-1, ::-1], mode='same')
    den = np.sqrt(np.clip(fftconvolve(inv, (k ** 2)[::-1, ::-1], mode='same'), 1e-30, None))
    scale = empirical_noise_scale(fit, exp)
    pscale = empirical_noise_scale(fit, exp, matched=False)
    snr_map = num / den / scale
    snr_map[~good] = 0.0
    print(f"empirical noise scale: S/N (matched) = {scale:.2f}, pixel = {pscale:.2f}")
    print(f"{'N sigma':>8s} {'>= +N':>8s} {'<= -N':>8s}")
    for t in thresholds:
        print(f"{t:8.1f} {(snr_map >= t).sum():8d} {(snr_map <= -t).sum():8d}")
    if plot:
        bbox = exp.getBBox()
        fig, ax = plt.subplots(1, 2, figsize=(14, 5.5))
        lim = max(5.0, float(np.nanpercentile(np.abs(snr_map), 99.9)))
        im = ax[0].imshow(snr_map, extent=_extent(bbox), origin='lower', cmap='RdBu_r',
                          vmin=-lim, vmax=lim)
        plt.colorbar(im, ax=ax[0], fraction=0.046, pad=0.04)
        pk = _peaks_array(h['result']['snaps'][snap_idx]['merged_arr'])
        if len(pk):
            ax[0].scatter(pk[:, 0], pk[:, 1], facecolors='none', edgecolors='k', s=30, lw=0.6)
        ax[0].set_title('Matched-filter S/N of residual (circles = fitted sources)')
        bins = np.linspace(-lim, lim, 61)
        ax[1].hist(snr_map[good], bins=bins, color='#4C72B0', alpha=0.8)
        ax[1].set_yscale('log')
        ax[1].set_xlabel('S/N of residual'); ax[1].set_title('Histogram (should be symmetric)')
        plt.tight_layout(); plt.show()
    return snr_map


def source_snr_histogram(outer_history, outer_idx=-1, snap_idx=-1, min_snr=5.0):
    """Plot the histogram of the honest S/N of every final source.

    A real population shows a gap or knee; a smooth pile-up just above ``min_snr`` means noise.

    Parameters
    ----------
    outer_history : `list` [`dict`]
        Outer-loop history.
    outer_idx : `int`, optional
        Outer pass (default: last).
    snap_idx : `int`, optional
        Inner-loop snapshot (default: last).
    min_snr : `float`, optional
        S/N threshold marked in the plot.

    Returns
    -------
    snr : `numpy.ndarray`
        S/N per source.
    """
    h = outer_history[outer_idx]
    fit = h['result']['snaps'][snap_idx]['fit']
    scale = empirical_noise_scale(fit, h['exposure'])
    snr = amp_snr_full(fit, scale)
    print(f"{len(snr)} sources; honest S/N: min={snr.min():.1f}, "
          f"median={np.median(snr):.1f}, max={snr.max():.0f}; "
          f"{(snr < min_snr + 2).sum()} within +2 of the cut")
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(np.log10(np.clip(snr, 1e-1, None)), bins=30, color='#55A868')
    ax.axvline(np.log10(min_snr), color='r', ls='--', label=f'min_snr={min_snr:g}')
    ax.set_xlabel('log10 honest S/N'); ax.set_ylabel('# sources'); ax.legend()
    ax.set_title('Final source S/N distribution')
    plt.tight_layout(); plt.show()
    return snr



# ######################################################################
# [from notebook cell: final comparison]
# ######################################################################

"""
final_comparison_with_edge_margin.py
====================================
Paste AFTER visualize_cutout_fit.py has loaded. Replaces plot_final_comparison.

Adds a magenta boundary showing where sources get dropped by edge_margin_fwhm
(sources inside the shaded band were removed by the pipeline; sources centred
outside the image can never be detected at all).
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import ListedColormap
from astropy.visualization import simple_norm

EDGE_COLOR = '#ff00ff'   # magenta: distinct from red (masked) and the RdBu / mako maps


def draw_edge_margin(ax, bbox, margin_px, color=EDGE_COLOR, shade=True):
    """Outline and shade the band within ``margin_px`` of the image border.

    This is the region removed by ``drop_out_of_bounds(..., margin=-margin_px)``.

    Parameters
    ----------
    ax : `matplotlib.axes.Axes`
        Axes to draw on.
    bbox : `lsst.geom.Box2I`
        Image box.
    margin_px : `float`
        Band width [px].
    color : `str`, optional
        Line colour.
    shade : `bool`, optional
        Shade the band.
    """
    if not margin_px or margin_px <= 0:
        return
    x0, x1 = bbox.getMinX(), bbox.getMaxX()
    y0, y1 = bbox.getMinY(), bbox.getMaxY()
    if shade:
        X, Y = np.meshgrid(np.arange(x0, x1 + 1), np.arange(y0, y1 + 1))
        inner = ((X >= x0 + margin_px) & (X <= x1 - margin_px) &
                 (Y >= y0 + margin_px) & (Y <= y1 - margin_px))
        band = ~inner
        ax.imshow(np.ma.masked_where(~band, np.ones(band.shape)),
                  extent=_extent(bbox), origin='lower',
                  cmap=ListedColormap([color]), alpha=0.18,
                  interpolation='nearest', zorder=3)
    ax.add_patch(patches.Rectangle((x0 + margin_px, y0 + margin_px),
                                   (x1 - x0) - 2 * margin_px,
                                   (y1 - y0) - 2 * margin_px,
                                   fill=False, ec=color, lw=1.6, ls='--', zorder=4))


def plot_final_comparison(outer_history, exposure_crop, bad_bitmask=None,
                          edge_margin_fwhm=0.5, view=None,
                          show_peaks=True, circle_radius=None, figsize=(26, 7)):
    """Plot original | final model | residual | residual/sigma over the fitted region.

    Parameters
    ----------
    outer_history : `list` [`dict`]
        Outer-loop history.
    exposure_crop : `lsst.afw.image.Exposure`
        The input crop (original image panel).
    bad_bitmask : `int`, optional
        Mask bits of bad pixels.
    edge_margin_fwhm : `float`, optional
        Edge margin used in the run; draws its boundary (None / 0 = do not draw).
    view : `tuple` [`float`, `float`, `float`, `float`], optional
        (xmin, xmax, ymin, ymax) zoom.
    show_peaks : `bool`, optional
        Mark the sources.
    circle_radius : `float`, optional
        Marker radius [px].
    figsize : `tuple` [`float`, `float`], optional
        Figure size [inches].
    """
    final = outer_history[-1]
    exp = final['exposure']
    res = final['result']
    model = res['model_full']
    peaks = _peaks_array(res['peaks'])
    bbox = exp.getBBox()

    orig_arr = exposure_crop[bbox].image.array
    resid = exp.image.array - model
    sigma = np.sqrt(np.where(exp.variance.array > 0, exp.variance.array, np.nan))
    zmap = resid / sigma
    bad = None if bad_bitmask is None else (exp.mask.array & bad_bitmask) != 0

    ext = _extent(bbox)
    norm = simple_norm(orig_arr, vmin=np.nanmin(orig_arr),
                       vmax=np.nanpercentile(orig_arr, 99),
                       stretch='log', log_a=1000)
    rsel = resid[~bad] if bad is not None and (~bad).any() else resid
    rs = np.nanpercentile(np.abs(rsel), 99)
    zs = 5.0
    radius = _circle_radius(exp, circle_radius)
    margin_px = (edge_margin_fwhm * _fwhm(exp.getPsf(), bbox.getCenter())
                 if edge_margin_fwhm else 0.0)

    fig, axes = plt.subplots(1, 4, figsize=figsize, sharex=True, sharey=True)
    panels = [
        (orig_arr, dict(norm=norm, cmap=CMAP), 'Original visit image (cutout)'),
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
            _overlay_bad(ax, bad, ext)
        draw_edge_margin(ax, bbox, margin_px)
        if show_peaks and len(peaks):
            for x, y in peaks:
                ax.add_patch(patches.Circle((x, y), radius, fill=False, ec='w',
                                            lw=0.5, alpha=0.6))
        ax.set_title(title, fontsize=11)
        ax.set_xlabel('x [px]')
        _clean_axes(ax)

    if view is not None:
        axes[0].set_xlim(view[0], view[1])
        axes[0].set_ylim(view[2], view[3])   # sharex/sharey propagates to all panels
    else:
        axes[0].set_xlim(ext[0], ext[1])
        axes[0].set_ylim(ext[2], ext[3])

    if margin_px:
        axes[0].legend(handles=[
            patches.Patch(fc=EDGE_COLOR, alpha=0.3, ec=EDGE_COLOR, ls='--',
                          label=f'edge drop zone ({edge_margin_fwhm} FWHM = {margin_px:.1f}px)'),
            patches.Patch(fc='red', alpha=0.35, label='masked (bad_bitmask)')],
            loc='upper left', fontsize=8)
    axes[0].set_ylabel('y [px]')
    plt.tight_layout()
    plt.show()
    plt.close(fig)



# ######################################################################
# [from notebook cell: representative sources]
# ######################################################################

"""
representative_sources.py
=========================
Pick a representative sample of the fitted sources and show, for each one,
data | model | residual | x-section | y-section.

Paste AFTER visualize_cutout_fit.py (helpers, pipeline). This file carries its own
copy of the UPDATED plot_source_pages (pdf= / show= / snr= / dpi= arguments), so it
works even if your session still has the original one. Being defined here, it also
replaces the old one for anything you run afterwards.

Sampling strategy (each source appears at most once, in the first category that
claims it, so the sample is diverse rather than 12 copies of the brightest star):

  brightest      top-N by fitted amplitude
  bright / mid / faint   random draw from the top / middle / bottom third of honest S/N
  marginal       lowest honest S/N (closest to the min_snr cut: most likely false)
  crowded        smallest nearest-neighbour distance (blends)
  worst fit      highest local chi2 in the source core (PSF / model failures)
  best fit       lowest local chi2 among the brighter half (what "good" looks like)
  near edge      smallest distance to the image border (edge behaviour)

Each category becomes its own page(s) so you can see immediately which regime a
row belongs to; a table printed to the console lists id, position and metrics.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from collections import OrderedDict
from scipy.spatial import cKDTree
from matplotlib.backends.backend_pdf import PdfPages
from astropy.visualization import simple_norm


# =====================================================================
# plot_source_pages (updated: pdf= / show= / snr= / dpi=)
# =====================================================================

def plot_source_pages(exposure, snap, outer_idx, pages=None, colors=None,
                      source_view='all', sources_per_page=6, max_pages=None,
                      view_half=20, circle_radius=None,
                      pdf=None, show=True, snr=None, dpi=80):
    """Plot one row per source: data | model | residual | x-section | y-section.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        The fitted exposure.
    snap : `dict`
        Inner-loop snapshot.
    outer_idx : `int`
        Outer-pass index.
    pages : `list` [`list` [`int`]], optional
        Source indices per page (default: `select_source_pages`).
    colors : `dict` [`int`, `str`], optional
        Colour per source.
    source_view : `str`, optional
        See `select_source_pages`.
    sources_per_page : `int`, optional
        See `select_source_pages`.
    max_pages : `int`, optional
        See `select_source_pages`.
    view_half : `int`, optional
        Half-size of each window [px].
    circle_radius : `float`, optional
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
    model = fit['model_full']
    resid = data - model
    peaks, amps = snap['merged_arr'], fit['amps']
    prev_amp, is_new, bad = snap['prev_amp'], snap['is_new'], snap['bad_mask']

    if pages is None:
        pages = select_source_pages(snap, source_view, sources_per_page, max_pages)
    if not pages:
        print("  (no sources selected for per-source pages)")
        return
    if colors is None:
        colors = source_colors(pages)
    radius = _circle_radius(exposure, circle_radius)

    for p, page in enumerate(pages):
        n = len(page)
        fig, axes = plt.subplots(n, 5, figsize=(21, 3.3 * n), squeeze=False)
        for r, j in enumerate(page):
            col = colors[int(j)]
            px, py = peaks[j]
            comp = fit['comps'][j]
            ix, iy, xa, xb, ya, yb = _view_window(bbox, px, py, view_half)
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
            own_sig = np.zeros((ny, nx), bool)
            if comp is not None:
                own[comp['ys'], comp['xs']] = comp['amp'] * comp['stamp']
                own_sig[comp['ys'], comp['xs']] = comp['sig']
            wxs, wys = X0 + np.arange(xa, xb), Y0 + np.arange(ya, yb)

            # --- image windows -------------------------------------------------
            for c_i, (arr, kw) in enumerate([
                    (dwin, dict(norm=nrm, cmap=CMAP)),
                    (mwin, dict(norm=nrm, cmap=CMAP)),
                    (rwin, dict(cmap='RdBu_r', vmin=-rs, vmax=rs))]):
                ax = axes[r, c_i]
                im = ax.imshow(arr, extent=wext, origin='lower', **kw)
                plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                if own_sig[ya:yb, xa:xb].any():
                    ax.contour(wxs, wys, own_sig[ya:yb, xa:xb].astype(float), levels=[0.5],
                               colors='gold', linewidths=1.2)
                ax.scatter([px], [py], marker='x', c=col, s=40, lw=1.5, zorder=6)
                _mark_source(ax, px, py, col, radius, label=j if c_i == 0 else None)
                inwin = ((peaks[:, 0] > wext[0]) & (peaks[:, 0] < wext[1]) &
                         (peaks[:, 1] > wext[2]) & (peaks[:, 1] < wext[3]))
                inwin[j] = False
                if inwin.any():
                    ax.scatter(peaks[inwin, 0], peaks[inwin, 1], marker='+', c='w', s=30, lw=0.9)
                _overlay_bad(ax, bad[ya:yb, xa:xb], wext)
                if comp is not None:
                    x0b, x1b, y0b, y1b = comp['bbox_parent']
                    ax.add_patch(patches.Rectangle((x0b - 0.5, y0b - 0.5), x1b - x0b + 1,
                                                   y1b - y0b + 1, fill=False, ec=col, lw=1.5,
                                                   ls='--', path_effects=_OUTLINE))
                ax.set_xlim(wext[0], wext[1]); ax.set_ylim(wext[2], wext[3])
                _clean_axes(ax)
                _frame_axes(ax, col)

            g = fit['group_labels'][j]
            tag = 'NEW' if is_new[j] else 'refit'
            dtxt = ''
            if np.isfinite(prev_amp[j]) and prev_amp[j] != 0:
                dtxt = f"  Δ={100 * (amps[j] - prev_amp[j]) / abs(prev_amp[j]):+.1f}%"
            stxt = f"  S/N={snr[j]:.1f}" if snr is not None else ''
            zflag = '  [AT ZERO BOUND]' if amps[j] <= 1e-6 else ''
            axes[r, 0].set_title(f"#{j} ({px:.1f},{py:.1f}) [{tag}] G{g}", fontsize=9, loc='left',
                                 bbox=dict(fc=col, ec='k', pad=2))
            axes[r, 1].set_title(f"model  amp={amps[j]:.1f}{dtxt}{stxt}{zflag}",
                                 fontsize=9, loc='left')
            axes[r, 2].set_title(f"residual  (group χ²ν={fit['group_chi2'].get(g, np.nan):.2f})",
                                 fontsize=9, loc='left')
            if r == 0:
                axes[r, 0].text(0.02, 0.02,
                                'circle = this source\ndashed box = PSF stamp\ngold = significant support',
                                transform=axes[r, 0].transAxes, fontsize=6, color='w',
                                bbox=dict(fc='k', alpha=0.5, lw=0))

            # --- cross-sections -----------------------------------------------
            for c_i, axis in ((3, 'x'), (4, 'y')):
                ax = axes[r, c_i]
                if axis == 'x':
                    coord, sl, centre = wxs, (iy, slice(xa, xb)), px
                    c0, c1 = ((comp['bbox_parent'][0], comp['bbox_parent'][1]) if comp else (None, None))
                else:
                    coord, sl, centre = wys, (slice(ya, yb), ix), py
                    c0, c1 = ((comp['bbox_parent'][2], comp['bbox_parent'][3]) if comp else (None, None))
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
                ax.fill_between(coord, 0, 1, where=own_sig[sl], transform=ax.get_xaxis_transform(),
                                color='gold', alpha=0.15)
                if badline.any():
                    ax.plot(coord[badline], np.zeros(int(badline.sum())), '|', color='red',
                            transform=ax.get_xaxis_transform(), ms=8)
                ax.set_xlim(coord[0], coord[-1])
                ax.set_xlabel(f'{axis} [px]')
                ax.set_title(f"{axis}-cross-section through source #{j}", fontsize=9)
                _frame_axes(ax, col)
                if r == 0 and c_i == 3:
                    ax.legend(fontsize=6, loc='upper right')

        fig.suptitle(f"Outer {outer_idx + 1} / inner iter {snap['it'] + 1}: per-source fits "
                     f"(page {p + 1}/{len(pages)}, view={source_view}); "
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

def source_metrics(outer_history, outer_idx=-1, snap_idx=-1, core_fwhm=1.5):
    """Return per-source quantities used to stratify the sample.

    Parameters
    ----------
    outer_history : `list` [`dict`]
        Outer-loop history.
    outer_idx : `int`, optional
        Outer pass (default: last).
    snap_idx : `int`, optional
        Inner-loop snapshot (default: last).
    core_fwhm : `float`, optional
        Radius of the chi2 region [FWHM].

    Returns
    -------
    metrics : `dict`
        ``peaks``, ``amps``, ``snr`` (honest S/N), ``nn`` (distance to the nearest other source
        [px]), ``edge`` (distance to the image border [px]), ``chi2`` (mean resid^2 / sigma^2
        within ``core_fwhm`` FWHM), ``group`` (size of the coupled group), ``fwhm``,
        ``outer_pos`` and the exposure and snapshot used.
    """
    outer_pos = outer_idx if outer_idx >= 0 else len(outer_history) + outer_idx
    h = outer_history[outer_pos]
    exposure = h['exposure']
    snap = h['result']['snaps'][snap_idx]
    fit = snap['fit']
    peaks = np.asarray(snap['merged_arr'], float)
    bbox = exposure.getBBox()
    X0, Y0 = bbox.getMinX(), bbox.getMinY()
    ny, nx = bbox.getHeight(), bbox.getWidth()

    snr = amp_snr_full(fit, empirical_noise_scale(fit, exposure))
    fwhm = _fwhm(exposure.getPsf(), bbox.getCenter())
    r = core_fwhm * fwhm

    if len(peaks) > 1:
        d, _ = cKDTree(peaks).query(peaks, k=2)
        nn = d[:, 1]
    else:
        nn = np.full(len(peaks), np.inf)

    edge = np.minimum.reduce([peaks[:, 0] - bbox.getMinX(), bbox.getMaxX() - peaks[:, 0],
                              peaks[:, 1] - bbox.getMinY(), bbox.getMaxY() - peaks[:, 1]])

    resid = exposure.image.array - fit['model_full']
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

    group = np.bincount(fit['group_labels'])[fit['group_labels']]
    return dict(peaks=peaks, amps=fit['amps'], snr=snr, nn=nn, edge=edge, chi2=chi2,
                group=group, fwhm=fwhm, outer_pos=outer_pos, exposure=exposure,
                snap=snap, fit=fit)


# =====================================================================
# Selection
# =====================================================================

def select_representative_sources(m, n_per_category=4, seed=0, categories=None,
                                  extra_ids=None):
    """Select a representative sample of sources per category.

    Parameters
    ----------
    m : `dict`
        Output of `source_metrics`.
    n_per_category : `int` or `dict` [`str`, `int`], optional
        Sources drawn per category.
    seed : `int`, optional
        Random seed of the tercile draws.
    categories : `list` [`str`], optional
        Categories to keep, in display order.
    extra_ids : `list` [`int`], optional
        Source ids appended as a 'picked' category.

    Returns
    -------
    selection : `collections.OrderedDict` [`str`, `list` [`int`]]
        Source ids per category.
    """
    rng = np.random.default_rng(seed)
    snr, amps, nn, edge, chi2 = m['snr'], m['amps'], m['nn'], m['edge'], m['chi2']
    n_src = len(snr)
    used = set()

    def n_of(cat):
        return n_per_category.get(cat, 4) if isinstance(n_per_category, dict) else n_per_category

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
    sel['brightest'] = take(np.argsort(-amps), n_of('brightest'))

    lo, hi = np.quantile(snr, [1 / 3, 2 / 3])
    terciles = OrderedDict([('bright', np.flatnonzero(snr >= hi)),
                            ('mid', np.flatnonzero((snr >= lo) & (snr < hi))),
                            ('faint', np.flatnonzero(snr < lo))])
    for name, idx in terciles.items():
        sel[name] = take(rng.permutation(idx), n_of(name))

    sel['marginal (lowest S/N)'] = take(np.argsort(snr), n_of('marginal (lowest S/N)'))
    sel['crowded (nearest neighbour)'] = take(np.argsort(nn), n_of('crowded (nearest neighbour)'))

    chi_worst = np.where(np.isfinite(chi2), chi2, -np.inf)
    sel['worst fit (core chi2)'] = take(np.argsort(-chi_worst), n_of('worst fit (core chi2)'))

    bright_half = snr >= np.median(snr)
    chi_best = np.where(np.isfinite(chi2) & bright_half, chi2, np.inf)
    sel['best fit (bright half)'] = take(np.argsort(chi_best), n_of('best fit (bright half)'))

    sel['near edge'] = take(np.argsort(edge), n_of('near edge'))

    if extra_ids is not None:
        sel['picked'] = take(np.asarray(extra_ids, int), n_src)

    if categories is not None:
        sel = OrderedDict((c, sel[c]) for c in categories if c in sel)
    return OrderedDict((c, ids) for c, ids in sel.items() if ids)


# =====================================================================
# Plot
# =====================================================================

def plot_representative_sources(outer_history, outer_idx=-1, snap_idx=-1,
                                n_per_category=4, seed=0, categories=None,
                                extra_ids=None, sources_per_page=4, view_half=20,
                                pdf_path=None, show=True, dpi=80):
    """Select a representative sample and plot it with `plot_source_pages`.

    Each category is split into pages of ``sources_per_page`` rows, so a page never mixes
    categories; the printed table maps every source id to its category.

    Parameters
    ----------
    outer_history : `list` [`dict`]
        Outer-loop history.
    outer_idx : `int`, optional
        Outer pass (default: last).
    snap_idx : `int`, optional
        Inner-loop snapshot (default: last).
    n_per_category : `int` or `dict`, optional
        See `select_representative_sources`.
    seed : `int`, optional
        See `select_representative_sources`.
    categories : `list` [`str`], optional
        See `select_representative_sources`.
    extra_ids : `list` [`int`], optional
        See `select_representative_sources`.
    sources_per_page : `int`, optional
        Rows per page.
    view_half : `int`, optional
        Half-size of each window [px].
    pdf_path : `str`, optional
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
    m = source_metrics(outer_history, outer_idx, snap_idx)
    sel = select_representative_sources(m, n_per_category, seed, categories, extra_ids)

    # console table
    print(f"{'category':30s} {'id':>4s} {'x':>7s} {'y':>7s} {'amp':>10s} {'S/N':>7s} "
          f"{'nn[px]':>7s} {'chi2':>6s} {'grp':>4s}")
    for cat, ids in sel.items():
        for j in ids:
            x, y = m['peaks'][j]
            print(f"{cat:30s} {j:4d} {x:7.1f} {y:7.1f} {m['amps'][j]:10.1f} {m['snr'][j]:7.1f} "
                  f"{m['nn'][j]:7.2f} {m['chi2'][j]:6.2f} {m['group'][j]:4d}")

    pages, page_cat = [], []
    for cat, ids in sel.items():
        for i in range(0, len(ids), sources_per_page):
            pages.append(ids[i:i + sources_per_page])
            page_cat.append(cat)
    print("\npage -> category:")
    for p, cat in enumerate(page_cat, 1):
        print(f"  page {p:2d}: {cat}  (ids {pages[p - 1]})")
    if not pages:
        print("no sources selected")
        return sel

    colors = source_colors(pages)
    kw = dict(pages=pages, colors=colors, view_half=view_half, show=show,
              snr=m['snr'], dpi=dpi)
    if pdf_path is not None:
        with PdfPages(pdf_path) as pdf:
            plot_source_pages(m['exposure'], m['snap'], m['outer_pos'], pdf=pdf, **kw)
        print(f"saved {len(pages)} page(s) to {pdf_path}")
    else:
        plot_source_pages(m['exposure'], m['snap'], m['outer_pos'], pdf=None, **kw)
    return sel

