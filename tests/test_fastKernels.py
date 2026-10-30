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

"""Tests of the fast PSF-stamp kernels (C++ and NumPy implementations)."""

import unittest

import numpy as np

import lsstStubs

lsstStubs.install()

from lsst.meas.crowd import fastKernels as fk  # noqa: E402


def _kernel(half=20, sigma=2.4):
    u = np.arange(-half, half + 1)
    X, Y = np.meshgrid(u, u)
    k = np.exp(-0.5 * (X ** 2 + Y ** 2) / sigma ** 2) * (1 + 0.15 * np.cos(2 * np.arctan2(Y, X)))
    return k / k.sum()


def _gauss(half, sigma, cx, cy):
    u = np.arange(-half, half + 1)
    X, Y = np.meshgrid(u, u)
    g = np.exp(-0.5 * ((X - cx) ** 2 + (Y - cy) ** 2) / sigma ** 2)
    return g / g.sum()


class LanczosWeightsTestCase(unittest.TestCase):
    def testNormalised(self):
        for f in (-0.5, -0.31, 0.0, 0.2, 0.49):
            first, w = fk.lanczosWeights(f, 5)
            self.assertEqual(len(w), 10)
            self.assertAlmostEqual(w.sum(), 1.0, places=12)

    def testZeroShiftIsIdentity(self):
        first, w = fk.lanczosWeights(0.0, 5)
        self.assertAlmostEqual(w[-first], 1.0, places=12)
        self.assertAlmostEqual(np.abs(np.delete(w, -first)).max(), 0.0, places=12)


class ShiftTestCase(unittest.TestCase):
    def testShiftMatchesAnalyticGaussian(self):
        g = _gauss(20, 2.6, 0, 0)
        for dx, dy in ((0.3, -0.2), (-0.45, 0.1), (0.0, 0.5)):
            for use in (False, True) if fk.HAVE_CXX else (False,):
                out = fk.shiftImage(g, dx, dy, useCxx=use)
                ref = _gauss(20, 2.6, dx, dy)
                self.assertLess(np.abs(out - ref).max() / ref.max(), 3e-3)

    def testIntegerShift(self):
        g = _gauss(20, 2.0, 0, 0)
        out = fk.shiftImage(g, 1.0, 0.0, useCxx=False)
        np.testing.assert_allclose(out[:, 1:], g[:, :-1], atol=1e-12)

    @unittest.skipUnless(fk.HAVE_CXX, "C++ extension not built")
    def testCxxEqualsNumpy(self):
        k = _kernel()
        rng = np.random.default_rng(1)
        xs, ys = rng.uniform(10, 90, 50), rng.uniform(10, 90, 50)
        for taper in (None, (0.75, 1.0)):
            a = fk.renderStamps(k, xs, ys, taper=taper, useCxx=True)
            b = fk.renderStamps(k, xs, ys, taper=taper, useCxx=False)
            np.testing.assert_allclose(a[0], b[0], atol=1e-13 * k.max())
            np.testing.assert_array_equal(a[1], b[1])
            np.testing.assert_array_equal(a[2], b[2])
        a = fk.renderDerivatives(k, xs, ys, useCxx=True)
        b = fk.renderDerivatives(k, xs, ys, useCxx=False)
        np.testing.assert_allclose(a[0], b[0], atol=1e-12 * k.max())
        np.testing.assert_allclose(a[1], b[1], atol=1e-12 * k.max())


class RenderTestCase(unittest.TestCase):
    def testOriginConvention(self):
        k = _kernel(half=10)
        st, x0, y0 = fk.renderStamps(k, [50.4, 50.6, 50.5], [20.0, 19.5, 20.49], useCxx=False)
        # afw: origin = floor(x + 0.5) - half
        self.assertEqual(list(x0), [40, 41, 41])
        self.assertEqual(list(y0), [10, 10, 10])
        self.assertEqual(st.shape, (3, 21, 21))

    def testFluxConserved(self):
        k = _gauss(15, 2.0, 0, 0)
        st, _, _ = fk.renderStamps(k, [10.37], [12.81], useCxx=False)
        self.assertAlmostEqual(st[0].sum(), 1.0, places=5)

    def testDerivativesMatchAnalytic(self):
        """Central differences (step 0.25) of a shifted Gaussian = those of the exact Gaussian."""
        s, eps = 2.2, 0.25
        k = _gauss(15, s, 0, 0)
        x, y = 30.2, 40.35
        dx, dy, x0, y0 = fk.renderDerivatives(k, [x], [y], eps=eps)
        u = np.arange(31)
        X, Y = np.meshgrid(u + x0[0], u + y0[0])

        def g(cx, cy):
            v = np.exp(-0.5 * ((X - cx) ** 2 + (Y - cy) ** 2) / s ** 2)
            return v / (2 * np.pi * s ** 2)

        rx = (g(x + eps, y) - g(x - eps, y)) / (2 * eps)
        ry = (g(x, y + eps) - g(x, y - eps)) / (2 * eps)
        np.testing.assert_allclose(dx[0], rx, atol=1.5e-2 * np.abs(rx).max())  # Lanczos-5 accuracy, amplified by 1/(2 eps)
        np.testing.assert_allclose(dy[0], ry, atol=1.5e-2 * np.abs(ry).max())

    def testTaperMatchesReference(self):
        """Same as fitLib.psfStamp(taper=...): pedestal + raised cosine."""
        k = _kernel(half=12) + 1e-4
        st = fk.renderStamps(k, [5.3], [7.8], useCxx=False)[0][0]
        out = fk.taperStamp(st, 12.3, 11.8, (0.75, 1.0), useCxx=False)
        ring = np.concatenate([st[0, :], st[-1, :], st[:, 0], st[:, -1]])
        ref = st - np.median(ring)
        h = 0.5 * (min(st.shape) - 1)
        yy, xx = np.mgrid[:25, :25]
        r = np.hypot(xx - 12.3, yy - 11.8)
        w = np.clip((h - r) / (0.25 * h), 0, 1)
        np.testing.assert_allclose(out, ref * (0.5 - 0.5 * np.cos(np.pi * w)), atol=1e-15)
        if fk.HAVE_CXX:
            np.testing.assert_allclose(fk.taperStamp(st, 12.3, 11.8, (0.75, 1.0), useCxx=True),
                                       out, atol=1e-15)


class FastPsfTestCase(unittest.TestCase):
    def testStampClip(self):
        f = fk.FastPsf(_kernel(half=10))
        arr, b = f.stamp(5.2, 5.0, (0, 0, 99, 99))
        self.assertEqual(b, (0, 0, 15, 15))
        self.assertEqual(arr.shape, (16, 16))
        arr, b = f.stamp(500.0, 5.0, (0, 0, 99, 99))
        self.assertIsNone(arr)

    def testDerivativesClip(self):
        f = fk.FastPsf(_kernel(half=10))
        dx, dy, b = f.derivatives(50.1, 50.2, (45, 45, 99, 99))
        self.assertEqual(b, (45, 45, 60, 60))
        self.assertEqual(dx.shape, dy.shape)

    def testValidate(self):
        k = _kernel(half=10)
        f = fk.FastPsf(k)

        def good(x, y):
            st, x0, y0 = fk.renderStamps(k, [x], [y], useCxx=False)
            return st[0], (int(x0[0]), int(y0[0]))

        def wrong(x, y):
            st, x0, y0 = fk.renderStamps(k * 1.1, [x], [y], useCxx=False)
            return st[0], (int(x0[0]), int(y0[0]))

        self.assertTrue(f.validate(good, [(10.2, 11.7), (40.5, 3.3)])[0])
        self.assertFalse(f.validate(wrong, [(10.2, 11.7)])[0])


if __name__ == "__main__":
    unittest.main()
