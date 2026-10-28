"""Time PSF stamp rendering: lsst.afw psf.computeImage vs meas_crowd's fast kernels.

    setup -r /path/to/meas_crowd   (or PYTHONPATH=/path/to/meas_crowd/python)
    python examples/benchmark_stamps.py
"""
import time

import numpy as np

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.geom as geom
from lsst.meas.algorithms import KernelPsf

from lsst.meas.crowd import fastKernels as fk

h = 40                                              # 81 x 81 kernel, as in the pipeline
u = np.arange(-h, h + 1)
X, Y = np.meshgrid(u, u)
k = np.exp(-0.5 * (X ** 2 + Y ** 2) / 2.7 ** 2) + 0.003 * (1 + (np.hypot(X, Y) / 9) ** 2) ** -2.5
psf = KernelPsf(afwMath.FixedKernel(afwImage.ImageD(k / k.sum())), geom.Point2D(0, 0))
kern = psf.computeKernelImage(geom.Point2D(0, 0)).array
rng = np.random.default_rng(0)
xy = rng.uniform(0, 300, (400, 2))


def timeit(fn, label, n_calls):
    t = time.perf_counter()
    fn()
    dt = time.perf_counter() - t
    print(f"{label:52s} {1e3 * dt:8.1f} ms  ({1e6 * dt / n_calls:6.0f} us per stamp)")


timeit(lambda: [psf.computeImage(geom.Point2D(x, y)) for x, y in xy], "afw psf.computeImage, 400 stamps", 400)
timeit(lambda: fk.render_stamps(kern, xy[:, 0], xy[:, 1], use_cxx=False), "meas_crowd NumPy, 400 stamps", 400)
if fk.HAVE_CXX:
    timeit(lambda: fk.render_stamps(kern, xy[:, 0], xy[:, 1], use_cxx=True),
           "meas_crowd C++, 400 stamps", 400)
timeit(lambda: [[psf.computeImage(geom.Point2D(x + dx, y + dy))
                 for dx, dy in ((.25, 0), (-.25, 0), (0, .25), (0, -.25))] for x, y in xy],
       "afw derivatives (4 computeImage each), 400 sources", 400)
if fk.HAVE_CXX:
    timeit(lambda: fk.render_derivatives(kern, xy[:, 0], xy[:, 1], use_cxx=True),
           "meas_crowd C++ derivatives, 400 sources", 400)
