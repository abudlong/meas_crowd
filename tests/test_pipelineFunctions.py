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

"""Tests of the pipeline's building blocks on synthetic images with a known answer: PSF model
regularisation and correction, bright-star correction terms, model-error measurement, saturation
level, ring / pair / galaxy handling, position polish and the fast-stamp switch."""

import types
import unittest
from unittest import mock

import numpy as np

import lsstStubs

lsstStubs.install()

from lsst.meas.crowd import pipeline as pl  # noqa: E402
from lsst.meas.crowd import fastKernels as fk  # noqa: E402

L = pl.L
geom = L.geom


def gaussImage(shape, x, y, sigma, x0=0, y0=0):
    yy, xx = np.mgrid[y0:y0 + shape[0], x0:x0 + shape[1]]
    return np.exp(-0.5 * ((xx - x) ** 2 + (yy - y) ** 2) / sigma ** 2) / (2 * np.pi * sigma ** 2)


def lstsqFit(data, var, compFn):
    """A stand-in for the joint fit: unconstrained weighted least squares of compFn stamps."""
    shape = data.shape

    def fit(exposure, psf, S, bad, varModel, fitkw):
        comps = [compFn(x, y) for x, y in S["xy"]]
        A = np.array([lsstStubs.place(shape, c).ravel() for c in comps]).T
        w = 1 / np.sqrt(var.ravel())
        a = np.linalg.lstsq(A * w[:, None], data.ravel() * w, rcond=None)[0]
        return dict(amps=a, modelFull=(A @ a).reshape(shape), good=np.ones(shape, bool),
                    variance=var, comps=comps)
    return fit


class PsfKernelTestCase(unittest.TestCase):
    def testSmoothFillsHoles(self):
        H = 40
        u = np.arange(-H, H + 1)
        X, Y = np.meshgrid(u, u)
        R = np.hypot(X, Y)
        true = (1 + (np.hypot(X / 1.1, Y / 0.9) / 3.0) ** 2) ** -3.0
        true /= true.sum()
        bad = true.copy()
        for x, y in [(10, 0), (-10, 1), (0, -11), (7, 8)]:
            bad[H + y, H + x] *= 0.01
        bad += np.random.default_rng(1).normal(0, 2e-6, bad.shape) * (R > 6)   # wing pixel noise
        bad /= bad.sum()
        new = pl.smoothPsfKernel(bad, 4.3)
        self.assertGreaterEqual(new.min(), 0.0)
        self.assertAlmostEqual(new.sum(), 1.0, places=10)
        ratio = new[R <= 1.5 * 4.3] / bad[R <= 1.5 * 4.3]     # core kept (up to renormalisation)
        self.assertLess(np.ptp(ratio), 1e-9)
        self.assertLess(abs(ratio[0] - 1), 2e-3)
        for x, y in [(10, 0), (7, 8)]:
            self.assertGreater(new[H + y, H + x], 0.5 * true[H + y, H + x])         # hole filled
        m = (R > 6) & (R < 25)
        self.assertLess(np.abs(new - true)[m].sum(), np.abs(bad - true)[m].sum())
        # ellipticity survives
        self.assertAlmostEqual(np.sum(new * X ** 2) / np.sum(new * Y ** 2),
                               np.sum(true * X ** 2) / np.sum(true * Y ** 2), delta=0.08)

    def _coreCase(self, truthFn):
        n, h = 420, 20
        u = np.arange(-h, h + 1)
        X, Y = np.meshgrid(u, u)

        def g(s, dx=0):
            k = np.exp(-0.5 * ((X - dx) ** 2 + Y ** 2) / s ** 2)
            return k / k.sum()

        truthk, modelk = truthFn(g), g(2.2)

        def placek(k, x, y):
            m = np.zeros((n, n))
            m[int(y) - h:int(y) + h + 1, int(x) - h:int(x) + h + 1] = k
            return m

        rng = np.random.default_rng(3)
        xy = np.array([(x, y) for x in range(40, 390, 50) for y in range(40, 390, 50)], float)
        amps = rng.uniform(1e5, 6e5, len(xy))
        data = sum(a * placek(truthk, *p) for a, p in zip(amps, xy)) + rng.normal(0, 10, (n, n))
        var = np.full((n, n), 100.0)

        def mkpsf(k):
            return types.SimpleNamespace(computeKernelImage=lambda pt, k=k: types.SimpleNamespace(array=k),
                                         k=k)

        def fakefit(exposure, psf, S, bad, vm, kw):
            A = np.array([placek(psf.k, *p).ravel() for p in S["xy"]]).T
            a = np.linalg.lstsq(A, data.ravel(), rcond=None)[0]
            return dict(amps=a, modelFull=(A @ a).reshape(n, n), good=np.ones((n, n), bool))

        exp = lsstStubs.FakeExposure(data, var)
        S = pl._state(xy, fromGaia=True)
        P = dict(pl._params({}))
        P["_haloRejects"] = np.zeros((0, 2))
        with mock.patch.object(pl, "_fit", fakefit), \
                mock.patch.object(L, "kernelToLsstPsf", mkpsf, create=True), \
                mock.patch.object(L, "_fwhm", lambda p, pt: 5.0, create=True):
            fit = fakefit(exp, mkpsf(modelk), S, 0, None, {})
            new = pl.refinePsfMedian(exp, mkpsf(modelk), dict(S=S, snaps=[dict(fit=fit)]), 0, P, {},
                                     fit0=fit, verbose=False)
        return truthk, modelk, new.k

    def testCoreCorrectionFixesWrongPsf(self):
        truthk, modelk, newk = self._coreCase(lambda g: 0.85 * g(1.9) + 0.15 * g(1.9, dx=1.2))
        self.assertLess(np.abs(newk - truthk).sum(), 0.25 * np.abs(modelk - truthk).sum())
        self.assertGreaterEqual(newk.min(), 0.0)

    def testCoreCorrectionLeavesCorrectPsf(self):
        truthk, modelk, newk = self._coreCase(lambda g: g(2.2))
        np.testing.assert_allclose(newk, modelk, atol=1e-12)


class BrightStarTermsTestCase(unittest.TestCase):
    def _star(self, offset):
        n, h, s0, F = 80, 25, 1.9, 2.2e6
        yy, xx = np.mgrid[:n, :n]

        def g(x, y, s):
            k = np.exp(-0.5 * ((xx - x) ** 2 + (yy - y) ** 2) / s ** 2)
            return k / k.sum()

        x0, y0 = 40.2, 39.7
        psf = g(x0, y0, s0)
        truth = F * (0.97 * g(x0, y0, s0) + 0.03 * g(x0 + offset, y0, 0.9))   # sharper core
        var = 550 + truth
        data = truth + np.random.default_rng(0).normal(0, 1, (n, n)) * np.sqrt(var)
        sy, sx = slice(int(y0) - h, int(y0) + h + 1), slice(int(x0) - h, int(x0) + h + 1)
        fit = dict(amps=np.array([F]), good=np.ones((n, n), bool),
                   comps=[dict(stamp=psf[sy, sx], ys=sy, xs=sx)])
        exp = lsstStubs.FakeExposure(data, var)
        S = pl._state(np.array([[x0, y0]]), fromGaia=True)
        return exp, S, fit, F * psf, psf.max(), 2.355 * s0, (int(round(y0)), int(round(x0)))

    def testCoreRings(self):
        exp, S, fit, model, peak, fwhm, c = self._star(0.0)
        out = {}
        for rings in ((), (0.0, 1.0, 2.0, 3.0, 4.0)):
            corr, info = pl.brightStarCorrections(exp, S, fit, model, peak, 23.5, fwhm, 1,
                                                  radial=True, brightModes=4, coreRings=rings)
            out[bool(rings)] = (info[0]["chi2After"], abs((exp.image.array - model - corr)[c]))
        self.assertGreater(out[False][0], 1.5)
        self.assertLess(out[True][0], 1.1)
        self.assertLess(out[True][1], 0.1 * out[False][1])

    def testLopsidedCore(self):
        exp, S, fit, model, peak, fwhm, c = self._star(0.3)
        chi = {}
        for modes in (0, 2):
            corr, info = pl.brightStarCorrections(exp, S, fit, model, peak, 23.5, fwhm, 1,
                                                  radial=True, brightModes=4,
                                                  coreRings=(0.0, 1.0, 2.0, 3.0, 4.0),
                                                  coreModes=modes)
            chi[modes] = info[0]["chi2After"]
        self.assertLess(chi[2], chi[0])
        self.assertLess(chi[2], 1.0)


class ModelErrorTestCase(unittest.TestCase):
    def _run(self, epsTrue):
        n, h, s0 = 360, 14, 1.9
        u = np.arange(-h, h + 1)
        K = np.exp(-0.5 * (u[None] ** 2 + u[:, None] ** 2) / s0 ** 2)
        K /= K.sum()
        psf = lsstStubs.FakePsf(K)
        rng = np.random.default_rng(3)
        xy = np.array([(x, y) for x in range(25, 340, 30) for y in range(25, 340, 30)], float)
        fl = 10 ** rng.uniform(4, 6.3, len(xy))
        model = np.zeros((n, n))
        for f, (x, y) in zip(fl, xy):
            model[int(y) - h:int(y) + h + 1, int(x) - h:int(x) + h + 1] += f * K
        truth = model * (1 + epsTrue * rng.normal(0, 1, (n, n)))
        data = truth + rng.normal(0, 1, (n, n)) * np.sqrt(100 + truth)
        exp = lsstStubs.FakeExposure(data, 100 + model, psf=psf)
        res = dict(S=pl._state(xy, fromGaia=True), modelFull=model,
                   snaps=[dict(fit=dict(amps=fl, good=np.ones((n, n), bool)))])
        with mock.patch.object(L, "_fwhm", lambda p, pt: 2.355 * s0, create=True):
            eps, _ = pl.measureModelError(exp, res, dict(pl._params({})), verbose=False)
        return eps

    def testRecoversModelError(self):
        self.assertAlmostEqual(self._run(0.005), 0.005, delta=0.0015)

    def testZeroForPerfectModel(self):
        self.assertEqual(self._run(0.0), 0.0)


class SaturationLevelTestCase(unittest.TestCase):
    def testLevels(self):
        d = np.random.default_rng(0).normal(0, 10, (100, 100))
        m = np.zeros((100, 100), int)
        d[40:45, 40:45] = 98000
        m[40:45, 40:45] = 1
        d[70, 70] = 60000
        exp = lsstStubs.FakeExposure(d, np.ones_like(d), mask=m)
        P = dict(pl._params({}))
        self.assertEqual(pl.saturationLevel(exp, P)[0], 98000.0)
        d[40:45, 40:45] = 500                      # interpolated SAT pixels
        lvl, how = pl.saturationLevel(exp, P)
        self.assertEqual(lvl, 60000.0)
        self.assertIn("lower bound", how)
        m[:] = 0
        self.assertIsNone(pl.saturationLevel(exp, P)[0])
        P["nearSatLevel"] = 1.2e5
        self.assertEqual(pl.saturationLevel(exp, P)[0], 1.2e5)


class RingMimicTestCase(unittest.TestCase):
    def _case(self, lobes):
        n, s0, h = 120, 1.9, 12
        fwhm = 2.355 * s0
        yy, xx = np.mgrid[:n, :n]

        def comp(x, y):
            return lsstStubs.gaussianComp((n, n), x, y, s0, h)

        FB = 2e6
        r = np.hypot(xx - 60, yy - 60)
        th = np.arctan2(yy - 60, xx - 60)
        shape = (0.3 + 1.4 * np.cos(2 * th) ** 8) if lobes else 1.0
        ring = FB * 2.2e-5 * np.exp(-0.5 * ((r - 10.5) / 1.5) ** 2) * shape
        real = [(60 + 11 * np.cos(0.8), 60 + 11 * np.sin(0.8), 1e4 if lobes else 2.5e4),
                (60 + 7 * np.cos(4.0), 60 + 7 * np.sin(4.0), 1.5e4)]
        truth = lsstStubs.place((n, n), comp(60, 60), FB) + ring + \
            sum(lsstStubs.place((n, n), comp(x, y), f) for x, y, f in real)
        var = 100 + truth
        data = truth + np.random.default_rng(0).normal(0, 1, (n, n)) * np.sqrt(var)
        ang = [0, np.pi / 2, np.pi, 1.5 * np.pi] if lobes else [0.3, 1.3, 2.6, 3.5, 4.6, 5.6]
        fakes = [(60 + 10.5 * np.cos(a), 60 + 10.5 * np.sin(a)) for a in ang]
        xy = np.array([(60.0, 60.0)] + [(x, y) for x, y, _ in real] + fakes)
        comps = [comp(x, y) for x, y in xy]
        A = np.array([lsstStubs.place((n, n), c).ravel() for c in comps]).T
        amps = np.linalg.lstsq(A, data.ravel(), rcond=None)[0]
        S = pl._state(xy)
        S["fromGaia"][0] = True
        fit = dict(amps=amps, variance=var, good=np.ones((n, n), bool), comps=comps)
        exp = lsstStubs.FakeExposure(data, var)
        drop, nchk = pl.ringMimics(exp, S, fit, (A @ amps).reshape(n, n), dict(pl._params({})),
                                   fwhm, 1 / (2 * np.pi * s0 ** 2), 10.0)
        return set(drop), set(range(3, 3 + len(fakes)))

    def testSmoothRing(self):
        drop, fakes = self._case(False)
        self.assertEqual(drop, fakes)

    def testDiffractionLobes(self):
        drop, fakes = self._case(True)
        self.assertEqual(drop, fakes)


class SplitAndPolishTestCase(unittest.TestCase):
    def setUp(self):
        self.n, self.s0, self.h = 240, 1.9, 10
        n, s0, h = self.n, self.s0, self.h
        rng = np.random.default_rng(3)
        self.singles = [(x, y) for x in (30, 90, 150, 210) for y in (30, 90)]
        self.pairs = [(x, y) for x in (30, 90, 150, 210) for y in (150, 210)][:5]
        truth = np.zeros((n, n))
        modelXy = []
        for x, y in self.singles:
            truth += lsstStubs.place((n, n), lsstStubs.gaussianComp((n, n), x + .3, y - .2, s0, h),
                                     rng.uniform(4e3, 1.2e4))
            modelXy.append((x + .3, y - .2))
        for x, y in self.pairs:
            a, b, f1, f2 = (x, y), (x + 2.5, y + 0.8), 8e3, 3e3
            for p, f in ((a, f1), (b, f2)):
                truth += lsstStubs.place((n, n), lsstStubs.gaussianComp((n, n), *p, s0, h), f)
            modelXy.append(((f1 * a[0] + f2 * b[0]) / (f1 + f2), (f1 * a[1] + f2 * b[1]) / (f1 + f2)))
        self.var = 100 + truth
        self.data = truth + rng.normal(0, 1, (n, n)) * np.sqrt(self.var)
        self.modelXy = np.array(modelXy)
        self.fakefit = lstsqFit(self.data, self.var,
                                lambda x, y: lsstStubs.gaussianComp((n, n), x, y, s0, h))
        self.exp = lsstStubs.FakeExposure(self.data, self.var)

    def testSplitFindsPairs(self):
        from scipy.ndimage import gaussian_filter

        def msnr(resid, var, usable, fwhm, noiseMask):
            s = fwhm / 2.355
            num = gaussian_filter(resid, s) * 4 * np.pi * s * s
            den = np.sqrt(gaussian_filter(var, s / np.sqrt(2)) * 4 * np.pi * s * s / 2)
            return (num / den / np.sqrt(2),)

        S = pl._state(self.modelXy, fromGaia=True)
        with mock.patch.object(pl, "_fit", self.fakefit), mock.patch.object(pl.ca, "_matchedSnr", msnr):
            f0 = self.fakefit(self.exp, None, S, 0, None, {})
            S2, f2, ntry, nkeep = pl.splitTest(self.exp, None, S, f0, f0["modelFull"],
                                               dict(pl._params({})), 2.355 * self.s0,
                                               1 / (2 * np.pi * self.s0 ** 2), 10.0, 0, {})
        self.assertEqual(nkeep, len(self.pairs))
        new = S2["xy"][len(self.modelXy):]
        for x, y in self.pairs:      # one new source next to each pair ...
            self.assertLess(np.min(np.hypot(new[:, 0] - x - 1.2, new[:, 1] - y - 0.4)), 6.0)
        for x, y in self.singles:    # ... and none at a single star
            self.assertGreater(np.min(np.hypot(new[:, 0] - x, new[:, 1] - y)), 10.0)

    def testPolishFixesStalePositions(self):
        n, s0, h = self.n, self.s0, self.h

        def stamp(x, y):
            ix, iy = int(round(x)), int(round(y))
            yy, xx = np.mgrid[iy - h:iy + h + 1, ix - h:ix + h + 1]
            g = np.exp(-0.5 * ((xx - x) ** 2 + (yy - y) ** 2) / s0 ** 2) / (2 * np.pi * s0 ** 2)
            return g, xx, yy, geom.Box2I(geom.Point2I(ix - h, iy - h), geom.Point2I(ix + h, iy + h))

        def deriv(psf, x, y, eb):
            g, xx, yy, b = stamp(x, y)
            return g * (xx - x) / s0 ** 2, g * (yy - y) / s0 ** 2, b

        rng = np.random.default_rng(4)
        truth = np.array([(x + rng.uniform(-.5, .5), y + rng.uniform(-.5, .5))
                          for x in range(20, 230, 30) for y in range(20, 230, 30)])
        fl = 10 ** rng.uniform(3.5, 5.0, len(truth))
        img = sum(f * lsstStubs.place((n, n), lsstStubs.gaussianComp((n, n), x, y, s0, h))
                  for f, (x, y) in zip(fl, truth))
        var = np.full((n, n), 100.0)
        data = img + rng.normal(0, 10, (n, n))
        exp = lsstStubs.FakeExposure(data, var)
        fakefit = lstsqFit(data, var, lambda x, y: lsstStubs.gaussianComp((n, n), x, y, s0, h))
        stale = truth + rng.normal(0, 0.2, truth.shape) * (rng.uniform(size=(len(truth), 1)) < 0.4)
        S = pl._state(stale, fromGaia=True)
        S["ref"] = truth.copy()
        P = dict(pl._params({}))
        P["finalPolishRounds"] = 8
        with mock.patch.object(pl, "_fit", fakefit), \
                mock.patch.object(L, "computePsfDerivatives", deriv, create=True), \
                mock.patch.object(L, "psfStamp", lambda psf, x, y, bb, taper=None: stamp(x, y)[::3],
                                  create=True), \
                mock.patch.object(L, "bboxToSlices", lambda b, eb: (slice(b.getMinY(), b.getMaxY() + 1),
                                                                    slice(b.getMinX(), b.getMaxX() + 1)),
                                  create=True), \
                mock.patch.object(L, "WING_TAPER", None, create=True):
            f0 = fakefit(exp, None, S, 0, None, {})
            S2, f1, info = pl.finalPositionPolish(exp, None, S, f0, 0, P, {}, 2.355 * s0, s0)
        e0, e1 = np.hypot(*(stale - truth).T), np.hypot(*(S2["xy"] - truth).T)
        self.assertLess(info["chi1"], info["chi0"])
        self.assertLess(e1.mean(), 0.3 * e0.mean())
        self.assertLess(info["chi1"] / n ** 2, 1.02)


class EllipticalSourceTestCase(unittest.TestCase):
    n, s0, h = 120, 1.9, 14

    def _run(self, truthFn, xy, gaiaIdx):
        from scipy.signal import fftconvolve
        n, s0, h = self.n, self.s0, self.h
        big = geom.Box2I(geom.Point2I(0, 0), geom.Point2I(n - 1, n - 1))

        def gstamp(x, y):
            ix, iy = int(round(x)), int(round(y))
            yy, xx = np.mgrid[iy - h:iy + h + 1, ix - h:ix + h + 1]
            k = np.exp(-0.5 * ((xx - x) ** 2 + (yy - y) ** 2) / s0 ** 2)
            return k / k.sum(), geom.Box2I(geom.Point2I(ix - h, iy - h), geom.Point2I(ix + h, iy + h))

        def place(arr, b):
            m = np.zeros((n, n))
            m[b.getMinY():b.getMaxY() + 1, b.getMinX():b.getMaxX() + 1] = arr
            return m

        truth = truthFn(gstamp, place, fftconvolve)
        var = 500 + truth
        data = truth + np.random.default_rng(1).normal(0, 1, (n, n)) * np.sqrt(var)

        def fakefit(exposure, psf, S, bad, vm, kw):
            pl._setExtRule(S)
            comps = []
            for x, y in S["xy"]:
                a, bb = pl._psfStampCached(psf, x, y, big)
                comps.append(dict(stamp=a, ys=slice(bb.getMinY(), bb.getMaxY() + 1),
                                  xs=slice(bb.getMinX(), bb.getMaxX() + 1)))
            A = np.array([place(c["stamp"], geom.Box2I(geom.Point2I(c["xs"].start, c["ys"].start),
                                                       geom.Point2I(c["xs"].stop - 1, c["ys"].stop - 1))).ravel()
                          for c in comps]).T
            a = np.linalg.lstsq(A, data.ravel(), rcond=None)[0]
            return dict(amps=a, modelFull=(A @ a).reshape(n, n), good=np.ones((n, n), bool),
                        variance=var, comps=comps)

        exp = lsstStubs.FakeExposure(data, var)
        S = pl._state(np.array(xy, float))
        S["fromGaia"][gaiaIdx] = True
        S["ref"] = S["xy"].copy()
        pl._bfReset()
        pl.clearStampCache()
        with mock.patch.object(pl, "_fit", fakefit), \
                mock.patch.object(pl, "_rawPsfStamp", lambda psf, x, y, taper: gstamp(x, y)), \
                mock.patch.object(pl, "_ORIG", {"psfStamp": lambda psf, x, y, bb, taper=None: gstamp(x, y)}), \
                mock.patch.object(pl, "_HUGE", big):
            f0 = fakefit(exp, None, S, 0, None, {})
            S2, f2, ntry, nkeep = pl.ellipTest(exp, None, S, f0, f0["modelFull"], dict(pl._params({})),
                                               2.355 * s0, gstamp(50, 50)[0].max(), 10.0, 0, {}, s0,
                                               verboseEllip=False)
            chi0 = float(np.mean((data - f0["modelFull"]) ** 2 / var))
            chi1 = float(np.mean((data - f2["modelFull"]) ** 2 / var))
        pl._bfReset()
        pl.clearStampCache()
        return S2, nkeep, chi0, chi1

    def testGalaxyReplacesPoints(self):
        def truth(gstamp, place, conv):
            st, b = gstamp(60.3, 58.7)
            gal = place(conv(st, pl._ellipKernel(1.6, 0.6, np.radians(25)), mode="same"), b) * 4.4e4
            st2, b2 = gstamp(85.2, 80.4)
            return gal + place(st2, b2) * 4e4

        S2, nkeep, chi0, chi1 = self._run(truth, [[59.0, 58.0], [62.5, 60.0], [85.2, 80.4]], [0, 2])
        self.assertEqual(nkeep, 1)
        self.assertEqual(len(S2["xy"]), 2)
        self.assertLess(chi1, chi0)
        self.assertLess(chi1, 1.01)
        self.assertGreater(S2["ext"][0], 1.0)              # the Gaia source became the galaxy
        self.assertEqual(S2["ext"][1], 0.0)                # the star stays a point

    def testStarPlusGalaxy(self):
        def truth(gstamp, place, conv):
            st, b = gstamp(64.5, 60.0)
            gal = place(conv(st, pl._ellipKernel(1.8, 0.45, np.radians(20), "exp"), mode="same"), b) * 3e4
            return gal + place(*gstamp(59.2, 58.4)) * 4e4 + place(*gstamp(85.2, 80.4)) * 4e4

        S2, nkeep, chi0, chi1 = self._run(truth, [[59.2, 58.4], [64.0, 59.8], [85.2, 80.4]], [0, 2])
        self.assertEqual(nkeep, 1)
        self.assertEqual(len(S2["xy"]), 3)
        self.assertEqual(S2["ext"][0], 0.0)                # Gaia star kept as a point
        self.assertGreater(S2["ext"][1], 1.0)              # its neighbour became the galaxy
        self.assertLess(chi1, 1.01)

    def testStarStaysPoint(self):
        def truth(gstamp, place, conv):
            return place(*gstamp(60.3, 58.7)) * 4.4e4 + place(*gstamp(85.2, 80.4)) * 4e4

        S2, nkeep, chi0, chi1 = self._run(truth, [[60.3, 58.7], [85.2, 80.4]], [0, 1])
        self.assertEqual(nkeep, 0)


class FastStampSwitchTestCase(unittest.TestCase):
    """The pipeline uses the fast stamps only when they agree with psf.computeImage."""

    def setUp(self):
        self.k = lsstStubs.gaussianKernel(2.3, 12)
        self.psf = lsstStubs.FakePsf(self.k)
        pl._FAST["cache"].clear()
        pl._FAST["on"] = True

    def tearDown(self):
        pl._FAST["cache"].clear()

    def _ref(self, kernel):
        def ref(psf, x, y, bb, taper=None):
            st, x0, y0 = fk.renderStamps(kernel, [x], [y], useCxx=False)
            return st[0], geom.Box2I(geom.Point2I(int(x0[0]), int(y0[0])),
                                     geom.Point2I(int(x0[0]) + 24, int(y0[0]) + 24))
        return ref

    def testAgreeingPsfUsesFastPath(self):
        with mock.patch.dict(pl._ORIG, {"psfStamp": self._ref(self.k)}), \
                mock.patch.object(pl, "_HUGE", geom.Box2I()):
            pl._FAST["verbose"] = False
            fast = pl._fastFor(self.psf)
            self.assertIsNotNone(fast)
            arr, bb = pl._rawPsfStamp(self.psf, 100.3, 50.6, None)
            self.assertEqual((bb.getMinX(), bb.getMinY()), (88, 39))
            ref = fk.renderStamps(self.k, [100.3], [50.6], useCxx=False)[0][0]
            np.testing.assert_allclose(arr, ref, atol=1e-14)

    def testVaryingPsfFallsBack(self):
        with mock.patch.dict(pl._ORIG, {"psfStamp": self._ref(self.k * 1.05)}), \
                mock.patch.object(pl, "_HUGE", geom.Box2I()):
            pl._FAST["verbose"] = False
            self.assertIsNone(pl._fastFor(self.psf))

    def testSwitchOff(self):
        pl._FAST["on"] = False
        self.assertIsNone(pl._fastFor(self.psf))
        pl._FAST["on"] = True


class BackgroundLowFitTestCase(unittest.TestCase):
    """The plane / constant background refit between outer passes."""

    def _exposure(self, seed=0):
        rng = np.random.default_rng(seed)
        ny, nx = 120, 150
        yy, xx = np.mgrid[:ny, :nx]
        sky = 4.0 + 6.0 * (xx - nx / 2) / nx - 3.0 * (yy - ny / 2) / ny
        img = sky + rng.normal(0, 2.0, (ny, nx))
        mask = np.zeros((ny, nx), int)
        img[40:60, 50:80] += 500.0                       # a star, flagged DETECTED
        mask[40:60, 50:80] |= lsstStubs.FakeMask.planes["DETECTED"]
        img[100:110, 10:20] += 80.0                      # an unflagged blob: must be clipped
        exp = lsstStubs.FakeExposure(img.copy(), np.full((ny, nx), 4.0), mask=mask, x0=300, y0=500)
        return exp, sky, img

    def testPlaneRecovered(self):
        exp, sky, img = self._exposure()
        exp, bg = pl.refitBackgroundLowfit(exp, "plane", verbose=False)
        np.testing.assert_allclose(bg.array, sky, atol=0.15)
        np.testing.assert_allclose(exp.image.array, img - bg.array)      # subtracted in place

    def testConstantAndNone(self):
        exp, sky, _ = self._exposure(1)
        _, bg = pl.refitBackgroundLowfit(exp, "constant", verbose=False)
        self.assertAlmostEqual(float(bg.array.mean()), 4.0, delta=0.3)
        self.assertEqual(np.ptp(bg.array), 0.0)
        exp, _, img = self._exposure(2)
        _, bg = pl.refitBackgroundLowfit(exp, "none", verbose=False)
        self.assertFalse(bg.array.any())
        np.testing.assert_array_equal(exp.image.array, img)
        with self.assertRaises(ValueError):
            pl.refitBackgroundLowfit(exp, "spline")

    def testDefaultIsPlane(self):
        self.assertEqual(pl.PIPELINE_DEFAULTS["bgMode"], "plane")


class OffDetectorGaiaTestCase(unittest.TestCase):
    def testDropsOnlyBeyondDetectorEdge(self):
        img = np.zeros((50, 60))
        mask = np.zeros((50, 60), int)
        mask[:, -3:] = lsstStubs.FakeMask.planes["EDGE"]          # right edge = detector edge
        exp = lsstStubs.FakeExposure(img, img + 1, mask=mask, x0=100, y0=200)
        xy = np.array([[130.0, 220.0],     # inside -> kept
                       [165.0, 220.0],     # beyond the right (detector) edge -> dropped
                       [95.0, 220.0],      # beyond the left edge, real pixels next to it -> kept
                       [130.0, 255.0]])    # beyond the top edge, real pixels -> kept
        kept, rows = pl.dropOffDetectorGaia(exp, xy, np.arange(4))
        self.assertEqual(list(rows), [0, 2, 3])
        np.testing.assert_array_equal(kept, xy[[0, 2, 3]])

    def testNoEdgeMaskKeepsEverything(self):
        img = np.zeros((20, 20))
        exp = lsstStubs.FakeExposure(img, img + 1)
        xy = np.array([[-5.0, 10.0], [10.0, 10.0]])
        kept, rows = pl.dropOffDetectorGaia(exp, xy, np.array([7, 8]))
        self.assertEqual(list(rows), [7, 8])
        kept, rows = pl.dropOffDetectorGaia(exp, np.zeros((0, 2)), np.zeros(0, int))
        self.assertEqual(len(kept), 0)


class StellarConfigTestCase(unittest.TestCase):
    """CROWDED_STELLAR_FIELD_CONFIG and the command line's --config names."""

    def testPresetIsCrowdedPlusStellarSettings(self):
        from lsst.meas.crowd import session
        stellar = dict(psfStackIsoFwhm=1.0, psfStackMinSnr=30.0, psfStackN=80,
                       ellipTest=False, keepExtended=False)
        self.assertEqual(session.CROWDED_STELLAR_FIELD_CONFIG,
                         dict(session.CROWDED_FIELD_CONFIG, **stellar))
        cfg = session.makeConfig(session.CROWDED_STELLAR_FIELD_CONFIG)
        for k, v in stellar.items():
            self.assertEqual(cfg[k], v)
        self.assertEqual(cfg["minSnr"], session.CROWDED_FIELD_CONFIG["minSnr"])
        self.assertTrue(pl.PIPELINE_DEFAULTS["ellipTest"])          # library defaults unchanged

    def testCommandLineConfigNames(self):
        import os
        import runpy
        from lsst.meas.crowd import session
        script = os.path.join(os.path.dirname(__file__), os.pardir, "bin.src", "crowdDetector.py")
        ns = runpy.run_path(script, run_name="crowdDetector")
        self.assertEqual(set(ns["CONFIGS"]), {"stellar", "crowded"})
        self.assertIs(getattr(session, ns["CONFIGS"]["stellar"]),
                      session.CROWDED_STELLAR_FIELD_CONFIG)
        self.assertIs(getattr(session, ns["CONFIGS"]["crowded"]), session.CROWDED_FIELD_CONFIG)


if __name__ == "__main__":
    unittest.main()
