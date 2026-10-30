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

"""The fast stamps reproduce lsst.afw's psf.computeImage for a constant (Kernel) PSF.
Needs the Rubin stack; skipped otherwise."""

import unittest

import numpy as np

import lsstStubs

HAVE_STACK = lsstStubs.install()


@unittest.skipUnless(HAVE_STACK, "needs the Rubin stack (lsst.afw)")
class StackComparisonTestCase(unittest.TestCase):
    def setUp(self):
        import lsst.afw.image as afwImage
        import lsst.afw.math as afwMath
        import lsst.geom as geom
        from lsst.meas.algorithms import KernelPsf
        h = 20
        u = np.arange(-h, h + 1)
        X, Y = np.meshgrid(u, u)
        k = np.exp(-0.5 * ((X / 2.6) ** 2 + (Y / 2.2) ** 2)) + 0.002 * np.exp(-np.hypot(X, Y) / 6)
        img = afwImage.ImageD(k / k.sum())
        self.psf = KernelPsf(afwMath.FixedKernel(img), geom.Point2D(0, 0))
        self.geom = geom

    def testStampsMatchComputeImage(self):
        from lsst.meas.crowd.fastKernels import FastPsf
        fast = FastPsf(self.psf.computeKernelImage(self.geom.Point2D(0, 0)).array)
        rng = np.random.default_rng(0)
        for x, y in rng.uniform(-500, 500, (20, 2)):
            ref = self.psf.computeImage(self.geom.Point2D(x, y))
            arr, b = fast.stamp(x, y, (-10 ** 6, -10 ** 6, 10 ** 6, 10 ** 6))
            self.assertEqual(b[:2], (ref.getX0(), ref.getY0()))
            self.assertLess(np.abs(arr - ref.array).max() / ref.array.max(), 2e-3)

    def testDerivativesMatchFitLib(self):
        from lsst.meas.crowd import fitLib
        from lsst.meas.crowd.fastKernels import FastPsf
        fast = FastPsf(self.psf.computeKernelImage(self.geom.Point2D(0, 0)).array)
        clip = self.geom.Box2I(self.geom.Point2I(-1000, -1000), self.geom.Point2I(1000, 1000))
        for x, y in [(10.1, 20.2), (-33.4, 5.05), (100.49, -7.7)]:
            dx0, dy0, bb = fitLib.compute_psf_derivatives(self.psf, x, y, clip)
            dx1, dy1, b = fast.derivatives(x, y, (bb.getMinX(), bb.getMinY(), bb.getMaxX(), bb.getMaxY()))
            scale = np.abs(dx0).max()
            self.assertLess(np.abs(dx1 - dx0).max() / scale, 5e-3)
            self.assertLess(np.abs(dy1 - dy0).max() / scale, 5e-3)


if __name__ == "__main__":
    unittest.main()
