# This file is part of meas_crowd.
#
# Developed for the LSST Data Management System.
# See the COPYRIGHT file at the top-level directory of this distribution.

"""Minimal stand-ins for the Rubin stack, so the unit tests also run where it is not installed
(laptops, GitHub CI).  `install()` does nothing when ``lsst.afw`` is importable.

Only ``lsst.geom``'s boxes and points are functional; every other stack module is a placeholder
whose attributes are inert dummies (enough to import the package, not to run the stack).
"""

import importlib
import sys
import types

import numpy as np

HAVE_STACK = False


class _Dummy:
    """Inert placeholder: callable, subscriptable, any attribute."""

    def __init__(self, *a, **k):
        pass

    def __call__(self, *a, **k):
        return _Dummy()

    def __getattr__(self, name):
        return _Dummy()

    def __getitem__(self, k):
        return _Dummy()


class _StubModule(types.ModuleType):
    def __getattr__(self, name):
        if name.startswith("__"):
            raise AttributeError(name)
        return type(name, (_Dummy,), {})


# ---------------------------------------------------------------- lsst.geom (functional)
class Point2I(tuple):
    def __new__(cls, x=0, y=0):
        return super().__new__(cls, (int(x), int(y)))

    def getX(self):
        return self[0]

    def getY(self):
        return self[1]


class Point2D(tuple):
    def __new__(cls, x=0.0, y=0.0):
        if isinstance(x, tuple):                    # Point2D(point), like lsst.geom
            x, y = x
        return super().__new__(cls, (float(x), float(y)))

    def getX(self):
        return self[0]

    def getY(self):
        return self[1]


class Extent2I(Point2I):
    pass


class Box2I:
    def __init__(self, a=None, b=None):
        if a is None:                       # empty box
            self.x0, self.y0, self.x1, self.y1 = 0, 0, -1, -1
        elif isinstance(a, Box2I):
            self.x0, self.y0, self.x1, self.y1 = a.x0, a.y0, a.x1, a.y1
        elif isinstance(b, Extent2I):
            self.x0, self.y0 = a[0], a[1]
            self.x1, self.y1 = a[0] + b[0] - 1, a[1] + b[1] - 1
        else:
            self.x0, self.y0, self.x1, self.y1 = a[0], a[1], b[0], b[1]

    def getMinX(self):
        return self.x0

    def getMinY(self):
        return self.y0

    def getMaxX(self):
        return self.x1

    def getMaxY(self):
        return self.y1

    def getWidth(self):
        return self.x1 - self.x0 + 1

    def getHeight(self):
        return self.y1 - self.y0 + 1

    def getCenter(self):
        return Point2D(0.5 * (self.x0 + self.x1), 0.5 * (self.y0 + self.y1))

    def isEmpty(self):
        return self.x1 < self.x0 or self.y1 < self.y0

    def clip(self, other):
        self.x0, self.y0 = max(self.x0, other.x0), max(self.y0, other.y0)
        self.x1, self.y1 = min(self.x1, other.x1), min(self.y1, other.y1)

    def contains(self, p):
        return self.x0 <= p[0] <= self.x1 and self.y0 <= p[1] <= self.y1


_STACK_MODULES = [
    "lsst.geom", "lsst.afw", "lsst.afw.image", "lsst.afw.math", "lsst.afw.table",
    "lsst.afw.detection", "lsst.afw.geom", "lsst.meas", "lsst.meas.algorithms",
    "lsst.meas.algorithms.subtractBackground", "lsst.meas.algorithms.installGaussianPsf",
    "lsst.meas.algorithms.makePsfCandidates", "lsst.meas.algorithms.pcaPsfDeterminer",
    "lsst.meas.deblender", "lsst.meas.base", "lsst.utils", "lsst.utils.threads",
    "lsst.daf", "lsst.daf.butler",
]


def install():
    """Make ``lsst.meas.crowd`` importable without the stack.  Returns True if the real stack is
    present (stubs not installed)."""
    global HAVE_STACK
    if isinstance(sys.modules.get("lsst.afw.image"), _StubModule):
        return False                        # stubs already installed by an earlier test module
    try:
        importlib.import_module("lsst.afw.image")
        HAVE_STACK = True
        return True
    except ImportError:
        pass
    import lsst  # our own namespace package (python/lsst)
    import lsst.meas  # noqa: F401
    for name in _STACK_MODULES:
        if name in ("lsst.meas",):
            continue
        if name not in sys.modules:
            mod = _StubModule(name)
            mod.__path__ = []
            sys.modules[name] = mod
            parent, _, child = name.rpartition(".")
            setattr(sys.modules[parent], child, mod)
    g = sys.modules["lsst.geom"]
    g.Box2I, g.Point2I, g.Point2D, g.Extent2I = Box2I, Point2I, Point2D, Extent2I
    sys.modules["lsst.utils.threads"].disable_implicit_threading = lambda: None
    return False


# ---------------------------------------------------------------- test fixtures
class FakeMask:
    planes = {"SAT": 1, "NEARSAT": 2, "BAD": 4, "SAT_GROWN": 8, "DETECTED": 16, "EDGE": 32,
              "NO_DATA": 64}

    def __init__(self, arr):
        self.array = arr

    def getPlaneBitMask(self, name):
        return self.planes[name]

    def addMaskPlane(self, name):
        self.planes.setdefault(name, 1 << len(self.planes))


class FakePsf:
    """Constant PSF from a kernel image (computeKernelImage / computeImage like a KernelPsf)."""

    def __init__(self, kernel, fwhm=None):
        self.kernel = np.asarray(kernel, float)
        self.fwhm = fwhm

    def computeKernelImage(self, point=None):
        return types.SimpleNamespace(array=self.kernel)

    def getAveragePosition(self):
        return Point2D(0.0, 0.0)

    def computeShape(self, point=None):
        k = self.kernel / self.kernel.sum()
        hy, hx = k.shape[0] // 2, k.shape[1] // 2
        yy, xx = np.mgrid[-hy:hy + 1, -hx:hx + 1]
        ixx, iyy = float(np.sum(k * xx ** 2)), float(np.sum(k * yy ** 2))
        return types.SimpleNamespace(getIxx=lambda: ixx, getIyy=lambda: iyy)


class FakeExposure:
    def __init__(self, image, variance, mask=None, psf=None, x0=0, y0=0, photo_scale=None):
        ny, nx = image.shape
        self.image = types.SimpleNamespace(array=image)
        self.variance = types.SimpleNamespace(array=variance)
        self.mask = FakeMask(np.zeros((ny, nx), int) if mask is None else mask)
        self._psf = psf
        self._bbox = Box2I(Point2I(x0, y0), Point2I(x0 + nx - 1, y0 + ny - 1))
        self._scale = photo_scale

    def getBBox(self):
        return self._bbox

    def getPsf(self):
        return self._psf

    def getPhotoCalib(self):
        if self._scale is None:
            return None
        return types.SimpleNamespace(instFluxToNanojansky=lambda f: f * self._scale)

    def getWcs(self):
        raise RuntimeError("no WCS in tests")


def gaussian_kernel(sigma, half):
    u = np.arange(-half, half + 1)
    X, Y = np.meshgrid(u, u)
    k = np.exp(-0.5 * (X ** 2 + Y ** 2) / sigma ** 2)
    return k / k.sum()


def gaussian_comp(shape, x, y, sigma, half, x0=0, y0=0):
    """Joint-fit component dict (stamp, ys, xs) of a unit Gaussian at (x, y)."""
    ny, nx = shape
    ix, iy = int(round(x)) - x0, int(round(y)) - y0
    ya, yb = max(iy - half, 0), min(iy + half + 1, ny)
    xa, xb = max(ix - half, 0), min(ix + half + 1, nx)
    yy, xx = np.mgrid[ya:yb, xa:xb]
    g = np.exp(-0.5 * ((xx + x0 - x) ** 2 + (yy + y0 - y) ** 2) / sigma ** 2) / (2 * np.pi * sigma ** 2)
    return dict(stamp=g, ys=slice(ya, yb), xs=slice(xa, xb))


def place(shape, comp, amp=1.0):
    out = np.zeros(shape)
    out[comp["ys"], comp["xs"]] += amp * comp["stamp"]
    return out
