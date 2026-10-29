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
    """Return normalised Lanczos weights for a sub-pixel shift.

    Parameters
    ----------
    frac : `float`
        Shift [px].
    order : `int`, optional
        Lanczos order.

    Returns
    -------
    first : `int`
        Offset of the first weight: ``out[i] = sum_j w[j] * in[i + first + j]`` approximates
        ``in(i - frac)``.
    weights : `numpy.ndarray`
        The ``2 * order`` weights (summing to 1).
    """
    first = int(np.floor(-frac)) - order + 1
    x = -frac - (first + np.arange(2 * order))
    w = np.sinc(x) * np.sinc(x / order)
    w[np.abs(x) >= order] = 0.0
    return first, w / w.sum()


def _shift_matrix(n, frac, order):
    """Return the (n, n) matrix that shifts a vector by ``frac`` (zeros outside)."""
    first, w = lanczos_weights(frac, order)
    M = np.zeros((n, n))
    i = np.arange(n)
    for j, wj in enumerate(w):
        m = i + first + j
        ok = (m >= 0) & (m < n)
        M[i[ok], m[ok]] += wj
    return M


def shift_image(image, dx, dy, order=5, use_cxx=None):
    """Shift an image by a sub-pixel offset with a separable Lanczos kernel.

    Same convention as ``lsst.afw.math.offsetImage``: ``out(x, y) = in(x - dx, y - dy)``; pixels
    outside the input count as zero.

    Parameters
    ----------
    image : `numpy.ndarray`
        Image.
    dx, dy : `float`
        Shift [px].
    order : `int`, optional
        Lanczos order.
    use_cxx : `bool`, optional
        Use the C++ kernels (default: if built).

    Returns
    -------
    shifted : `numpy.ndarray`
        The shifted image.
    """
    image = np.ascontiguousarray(image, dtype=float)
    if (HAVE_CXX if use_cxx is None else use_cxx):
        return _lib.shiftImage(image, float(dx), float(dy), int(order))
    ny, nx = image.shape
    return _shift_matrix(ny, dy, order) @ image @ _shift_matrix(nx, dx, order).T


def taper_stamp(stamp, cx, cy, taper, use_cxx=None):
    """Tame the wings of a PSF stamp like ``fitLib.psf_stamp(taper=...)``.

    Subtract the median of the outer pixel ring, then apply a raised cosine from ``taper[0] * h``
    to ``taper[1] * h`` around (cx, cy), h being half the stamp size.

    Parameters
    ----------
    stamp : `numpy.ndarray`
        PSF stamp.
    cx, cy : `float`
        Centre in stamp pixel coordinates.
    taper : `tuple` [`float`, `float`]
        Inner and outer taper radii as fractions of h.
    use_cxx : `bool`, optional
        Use the C++ kernels (default: if built).

    Returns
    -------
    stamp : `numpy.ndarray`
        The tapered stamp.
    """
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
    """Split positions into nearest pixel (floor(x + 0.5)) and sub-pixel offset."""
    i = np.floor(np.asarray(v, float) + 0.5)
    return i.astype(np.int64), np.asarray(v, float) - i


def render_stamps(kernel, xs, ys, order=5, taper=None, use_cxx=None):
    """Render stamps of a constant PSF at many positions.

    Stamp i has origin (round(x) - nx//2, round(y) - ny//2) and is the kernel shifted by the
    position's sub-pixel offset, like ``psf.computeImage(Point2D(x, y))``.

    Parameters
    ----------
    kernel : `numpy.ndarray`
        PSF kernel image (odd size, centred).
    xs, ys : array-like
        Positions in parent pixel coordinates.
    order : `int`, optional
        Lanczos order.
    taper : `tuple` [`float`, `float`], optional
        Wing taper (see `taper_stamp`); None = no taper.
    use_cxx : `bool`, optional
        Use the C++ kernels (default: if built).

    Returns
    -------
    stamps : `numpy.ndarray`
        Shape (n, ny, nx).
    x0, y0 : `numpy.ndarray` [`int`]
        Stamp origins.
    """
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
    """Render central-difference derivative stamps dP/dx, dP/dy of a constant PSF.

    Parameters
    ----------
    kernel : `numpy.ndarray`
        PSF kernel image (odd size, centred).
    xs, ys : array-like
        Positions in parent pixel coordinates.
    eps : `float`, optional
        Finite-difference step [px].
    order : `int`, optional
        Lanczos order.
    use_cxx : `bool`, optional
        Use the C++ kernels (default: if built).

    Returns
    -------
    dx, dy : `numpy.ndarray`
        Derivative stamps, shape (n, ny, nx), on the centred stamps' origins.
    x0, y0 : `numpy.ndarray` [`int`]
        Stamp origins.
    """
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
    """Stamps and derivatives of a spatially constant PSF.

    Return conventions match ``fitLib.psf_stamp`` and ``fitLib.compute_psf_derivatives``, with
    boxes given as (xmin, ymin, xmax, ymax) tuples.

    Parameters
    ----------
    kernel : `numpy.ndarray`
        PSF kernel image (odd size, centred), e.g. ``psf.computeKernelImage(point).array``.
    order : `int`, optional
        Lanczos order (afw renders PSF images with lanczos5).
    """

    def __init__(self, kernel, order=5):
        """Store the kernel as a contiguous float array."""
        self.kernel = np.ascontiguousarray(kernel, dtype=float)
        self.order = order

    @staticmethod
    def _clip(x0, y0, ny, nx, clip):
        """Intersect a stamp box with a clip box; return (slice_y, slice_x, box) or None."""
        bx0, by0 = max(x0, clip[0]), max(y0, clip[1])
        bx1, by1 = min(x0 + nx - 1, clip[2]), min(y0 + ny - 1, clip[3])
        if bx0 > bx1 or by0 > by1:
            return None
        return (slice(by0 - y0, by1 - y0 + 1), slice(bx0 - x0, bx1 - x0 + 1), (bx0, by0, bx1, by1))

    def stamp(self, x, y, clip, taper=None):
        """Render the stamp at (x, y), clipped to a box.

        Parameters
        ----------
        x, y : `float`
            Position in parent pixel coordinates.
        clip : `tuple` [`int`, `int`, `int`, `int`]
            Clip box (xmin, ymin, xmax, ymax), inclusive.
        taper : `tuple` [`float`, `float`], optional
            Wing taper (see `taper_stamp`).

        Returns
        -------
        stamp : `numpy.ndarray` or `None`
            The clipped stamp (None if it misses ``clip``).
        box : `tuple` [`int`, `int`, `int`, `int`] or `None`
            Its box (xmin, ymin, xmax, ymax).
        """
        st, x0, y0 = render_stamps(self.kernel, [x], [y], self.order, taper)
        ny, nx = self.kernel.shape
        c = self._clip(int(x0[0]), int(y0[0]), ny, nx, clip)
        if c is None:
            return None, None
        return st[0][c[0], c[1]], c[2]

    def derivatives(self, x, y, clip, eps=0.25):
        """Render the derivative stamps at (x, y), clipped to a box.

        Parameters
        ----------
        x, y : `float`
            Position in parent pixel coordinates.
        clip : `tuple` [`int`, `int`, `int`, `int`]
            Clip box (xmin, ymin, xmax, ymax), inclusive.
        eps : `float`, optional
            Finite-difference step [px].

        Returns
        -------
        dpdx, dpdy : `numpy.ndarray` or `None`
            Derivative stamps (None if they miss ``clip``).
        box : `tuple` [`int`, `int`, `int`, `int`] or `None`
            Their box (xmin, ymin, xmax, ymax).
        """
        dx, dy, x0, y0 = render_derivatives(self.kernel, [x], [y], eps, self.order)
        ny, nx = self.kernel.shape
        c = self._clip(int(x0[0]), int(y0[0]), ny, nx, clip)
        if c is None:
            return None, None, None
        return dx[0][c[0], c[1]], dy[0][c[0], c[1]], c[2]

    def validate(self, reference, positions, rtol=2e-3):
        """Compare the fast stamps with a reference renderer.

        Parameters
        ----------
        reference : callable
            ``reference(x, y) -> (array, (xmin, ymin))``: the full, untapered stamp (e.g. from
            ``psf.computeImage``).
        positions : iterable [`tuple` [`float`, `float`]]
            Positions to compare at.
        rtol : `float`, optional
            Maximum allowed difference relative to the peak.

        Returns
        -------
        ok : `bool`
            Whether all stamps agree (same origin and shape, difference within ``rtol``).
        worst : `float`
            Largest difference relative to the peak.
        """
        worst = 0.0
        ny, nx = self.kernel.shape
        for x, y in positions:
            ref, (rx0, ry0) = reference(x, y)
            st, x0, y0 = render_stamps(self.kernel, [x], [y], self.order)
            if (int(x0[0]), int(y0[0])) != (rx0, ry0) or ref.shape != (ny, nx):
                return False, np.inf
            worst = max(worst, float(np.max(np.abs(st[0] - ref)) / max(np.max(np.abs(ref)), 1e-30)))
        return worst <= rtol, worst
