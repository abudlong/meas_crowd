# This file is part of meas_crowd.
#
# Developed for the LSST Data Management System.
# See the COPYRIGHT file at the top-level directory of this distribution.

"""Fast PSF stamps for a spatially CONSTANT PSF.

``lsst.afw`` renders a stamp with ``psf.computeImage(point)``: the kernel image, shifted by the
sub-pixel part of ``point`` with a Lanczos-5 kernel.  The pipeline asks for tens of thousands of
such stamps (five per source per centroid step for the derivatives alone), and the per-call
overhead dominates.  Here the kernel is taken ONCE and every shift is done in C++
(``_crowdLib``, pybind11) or, if the extension is not built, with an equivalent NumPy
implementation (separable Lanczos weight matrices).

`FastPsf.validate` compares against ``psf.computeImage`` and the pipeline only uses the fast
path when they agree (``fast_stamps`` option), so results do not change.
"""

__all__ = ["HAVE_CXX", "lanczos_weights", "shift_image", "taper_stamp", "render_stamps",
           "render_derivatives", "FastPsf"]

import numpy as np

try:
    from . import _crowdLib as _lib
    HAVE_CXX = True
except ImportError:                     # extension not built: NumPy fallback
    _lib = None
    HAVE_CXX = False


def lanczos_weights(frac, order=5):
    """(first, weights): normalised Lanczos weights so that
    ``out[i] = sum_j weights[j] * in[i + first + j]`` approximates ``in(i - frac)``."""
    first = int(np.floor(-frac)) - order + 1
    x = -frac - (first + np.arange(2 * order))
    w = np.sinc(x) * np.sinc(x / order)
    w[np.abs(x) >= order] = 0.0
    return first, w / w.sum()


def _shift_matrix(n, frac, order):
    """(n, n) matrix M with ``M @ v`` = v shifted by +frac (zeros outside)."""
    first, w = lanczos_weights(frac, order)
    M = np.zeros((n, n))
    i = np.arange(n)
    for j, wj in enumerate(w):
        m = i + first + j
        ok = (m >= 0) & (m < n)
        M[i[ok], m[ok]] += wj
    return M


def shift_image(image, dx, dy, order=5, use_cxx=None):
    """``image`` shifted by (dx, dy) px: out(x, y) = in(x - dx, y - dy) (afw offsetImage)."""
    image = np.ascontiguousarray(image, dtype=float)
    if (HAVE_CXX if use_cxx is None else use_cxx):
        return _lib.shiftImage(image, float(dx), float(dy), int(order))
    ny, nx = image.shape
    return _shift_matrix(ny, dy, order) @ image @ _shift_matrix(nx, dx, order).T


def taper_stamp(stamp, cx, cy, taper, use_cxx=None):
    """cutout_fit_lib.psf_stamp's wing taper: subtract the median of the outer pixel ring, then a
    raised cosine from taper[0]*h to taper[1]*h around (cx, cy) (stamp pixel coordinates)."""
    stamp = np.ascontiguousarray(stamp, dtype=float)
    if (HAVE_CXX if use_cxx is None else use_cxx):
        return _lib.taperStamp(stamp, float(cx), float(cy), float(taper[0]), float(taper[1]))
    ring = np.concatenate([stamp[0, :], stamp[-1, :], stamp[:, 0], stamp[:, -1]])
    arr = stamp - np.median(ring)
    h = 0.5 * (min(arr.shape) - 1)
    yy, xx = np.mgrid[:arr.shape[0], :arr.shape[1]]
    r = np.hypot(xx - cx, yy - cy)
    r0, r1 = taper[0] * h, taper[1] * h
    w = np.clip((r1 - r) / max(r1 - r0, 1e-9), 0.0, 1.0)
    return arr * (0.5 - 0.5 * np.cos(np.pi * w))


def _split(v):
    i = np.floor(np.asarray(v, float) + 0.5)
    return i.astype(np.int64), np.asarray(v, float) - i


def render_stamps(kernel, xs, ys, order=5, taper=None, use_cxx=None):
    """Stamps of a constant PSF at many positions.

    Returns (stamps[n, ny, nx], x0[n], y0[n]); stamp i has origin (x0[i], y0[i]) =
    (round(x) - nx//2, round(y) - ny//2), like ``psf.computeImage(Point2D(x, y))``."""
    kernel = np.ascontiguousarray(kernel, dtype=float)
    xs = np.ascontiguousarray(np.atleast_1d(xs), dtype=float)
    ys = np.ascontiguousarray(np.atleast_1d(ys), dtype=float)
    if (HAVE_CXX if use_cxx is None else use_cxx):
        t0, t1 = taper if taper is not None else (0.0, 1.0)
        return _lib.renderStamps(kernel, xs, ys, int(order), taper is not None, float(t0), float(t1))
    ny, nx = kernel.shape
    ix, fx = _split(xs)
    iy, fy = _split(ys)
    out = np.empty((len(xs), ny, nx))
    for k in range(len(xs)):
        out[k] = shift_image(kernel, fx[k], fy[k], order, use_cxx=False)
        if taper is not None:
            out[k] = taper_stamp(out[k], nx // 2 + fx[k], ny // 2 + fy[k], taper, use_cxx=False)
    return out, ix - nx // 2, iy - ny // 2


def render_derivatives(kernel, xs, ys, eps=0.25, order=5, use_cxx=None):
    """Central-difference dP/dx, dP/dy (step eps) on the centred stamps' origins.
    Returns (dx[n, ny, nx], dy[n, ny, nx], x0[n], y0[n])."""
    kernel = np.ascontiguousarray(kernel, dtype=float)
    xs = np.ascontiguousarray(np.atleast_1d(xs), dtype=float)
    ys = np.ascontiguousarray(np.atleast_1d(ys), dtype=float)
    if (HAVE_CXX if use_cxx is None else use_cxx):
        return _lib.renderDerivatives(kernel, xs, ys, int(order), float(eps))
    ny, nx = kernel.shape
    ix, fx = _split(xs)
    iy, fy = _split(ys)
    dx = np.empty((len(xs), ny, nx))
    dy = np.empty_like(dx)
    for k in range(len(xs)):
        dx[k] = (shift_image(kernel, fx[k] + eps, fy[k], order, False)
                 - shift_image(kernel, fx[k] - eps, fy[k], order, False)) / (2 * eps)
        dy[k] = (shift_image(kernel, fx[k], fy[k] + eps, order, False)
                 - shift_image(kernel, fx[k], fy[k] - eps, order, False)) / (2 * eps)
    return dx, dy, ix - nx // 2, iy - ny // 2


class FastPsf:
    """Stamps and derivatives of a spatially constant PSF, with the same return conventions as
    ``fitLib.psf_stamp`` and ``fitLib.compute_psf_derivatives``.

    Parameters
    ----------
    kernel : `numpy.ndarray`
        The PSF kernel image (odd size, centred), e.g. ``psf.computeKernelImage(point).array``.
    order : `int`
        Lanczos order (afw renders PSF images with lanczos5).
    """

    def __init__(self, kernel, order=5):
        self.kernel = np.ascontiguousarray(kernel, dtype=float)
        self.order = order

    @staticmethod
    def _clip(x0, y0, ny, nx, clip):
        """Intersection of the stamp box with clip=(xmin, ymin, xmax, ymax) (inclusive), as
        (slice_y, slice_x, (bx0, by0, bx1, by1)), or None if empty."""
        bx0, by0 = max(x0, clip[0]), max(y0, clip[1])
        bx1, by1 = min(x0 + nx - 1, clip[2]), min(y0 + ny - 1, clip[3])
        if bx0 > bx1 or by0 > by1:
            return None
        return (slice(by0 - y0, by1 - y0 + 1), slice(bx0 - x0, bx1 - x0 + 1), (bx0, by0, bx1, by1))

    def stamp(self, x, y, clip, taper=None):
        """(array, (xmin, ymin, xmax, ymax)) of the stamp at (x, y) clipped to ``clip``; the
        array is None if the stamp misses ``clip``."""
        st, x0, y0 = render_stamps(self.kernel, [x], [y], self.order, taper)
        ny, nx = self.kernel.shape
        c = self._clip(int(x0[0]), int(y0[0]), ny, nx, clip)
        if c is None:
            return None, None
        return st[0][c[0], c[1]], c[2]

    def derivatives(self, x, y, clip, eps=0.25):
        """(dP/dx, dP/dy, (xmin, ymin, xmax, ymax)) clipped to ``clip``; (None, None, None) if
        empty."""
        dx, dy, x0, y0 = render_derivatives(self.kernel, [x], [y], eps, self.order)
        ny, nx = self.kernel.shape
        c = self._clip(int(x0[0]), int(y0[0]), ny, nx, clip)
        if c is None:
            return None, None, None
        return dx[0][c[0], c[1]], dy[0][c[0], c[1]], c[2]

    def validate(self, reference, positions, rtol=2e-3):
        """Max |fast - reference| / peak over ``positions``; ``reference(x, y)`` must return
        (array, (xmin, ymin)) of the full, untapered stamp (e.g. from psf.computeImage).
        Returns (ok, worst)."""
        worst = 0.0
        ny, nx = self.kernel.shape
        for x, y in positions:
            ref, (rx0, ry0) = reference(x, y)
            st, x0, y0 = render_stamps(self.kernel, [x], [y], self.order)
            if (int(x0[0]), int(y0[0])) != (rx0, ry0) or ref.shape != (ny, nx):
                return False, np.inf
            worst = max(worst, float(np.max(np.abs(st[0] - ref)) / max(np.max(np.abs(ref)), 1e-30)))
        return worst <= rtol, worst
