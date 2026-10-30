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

"""Tests of the residual diagnostics and the run-level plots."""

import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import lsstStubs  # noqa: E402

lsstStubs.install()

from lsst.meas.crowd import pipeline as pl  # noqa: E402
from lsst.meas.crowd import plotting  # noqa: E402
import fakeRun  # noqa: E402

plt.show = lambda *a, **k: None


def _run(**kw):
    return fakeRun.makeRun(nSide=11, shape=(400, 400), half=15, sigma=1.9, sky=100.0,
                           fluxRange=(3e4, 2e6), seed=5, **kw)


class StarResidualTableTestCase(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def testPerfectModel(self):
        t = plotting.starResidualTable(_run(), plot=False, nWorst=0)
        self.assertTrue(any("PHOTON NOISE" in v for v in t["verdict"]))
        self.assertLess(np.median(t["chi2Final"]), 1.1)

    def testBrighterFatter(self):
        run = _run(sigmaFn=lambda f, x, y: 1.9 * (1 + 0.09 * (np.log10(f) - 5)))
        t = plotting.starResidualTable(run, plot=False, nWorst=0)
        self.assertTrue(any("depends on brightness" in v for v in t["verdict"]))
        self.assertLess(t["coreVsFlux"][0], -0.5)

    def testPositionDependence(self):
        run = _run(sigmaFn=lambda f, x, y: 1.9 * (1 + 0.06 * (x - 200) / 200))
        t = plotting.starResidualTable(run, plot=False, nWorst=0)
        self.assertTrue(any("POSITION" in v for v in t["verdict"]))

    def testPlot(self):
        plotting.starResidualTable(_run(), plot=True, nWorst=3)


class CoreExcessTestCase(unittest.TestCase):
    def testBroadPsf(self):
        out = plotting.coreExcessDiagnostic(_run(modelSigma=1.9 * 1.03), plot=False)
        self.assertTrue(any("too broad" in v for v in out["verdict"]))
        self.assertGreater(np.median(out["excess"]), 0.02)

    def testFluxDependent(self):
        run = _run(sigmaFn=lambda f, x, y: 1.9 * (1 - 0.03 * (np.log10(f) - 4.5)))
        out = plotting.coreExcessDiagnostic(run, plot=False)
        self.assertTrue(any("GROWS WITH FLUX" in v for v in out["verdict"]))

    def testPerfect(self):
        out = plotting.coreExcessDiagnostic(_run(), plot=False)
        self.assertEqual(out["verdict"], ["no significant central excess"])


class UnmodeledMaskTestCase(unittest.TestCase):
    def testNoPadAroundModeledCores(self):
        n, s0 = 160, 2.0
        fw = 2.355 * s0
        yy, xx = np.mgrid[:n, :n]
        star = lambda x, y, f: f * np.exp(-0.5 * ((xx - x) ** 2 + (yy - y) ** 2) / s0 ** 2) / (2 * np.pi * s0 ** 2)
        truth = star(50, 50, 6e5) + star(56, 53, 2e4) + star(110, 100, 3e5)
        model = star(50, 50, 6e5) + star(110, 100, 3e5)          # neighbour (56, 53) missing
        var = 100 + truth
        data = truth + np.random.default_rng(2).normal(0, 1, (n, n)) * np.sqrt(var)
        bad = np.zeros((n, n), bool)
        nearsat = bad.copy()
        nearsat[47:54, 47:54] = True
        bad |= nearsat
        bad[:, 140:143] = True                                    # a bad column
        old = pl.unmodeledNearMask(data - model, var, bad, fw)
        new = pl.unmodeledNearMask(data - model, var, bad, fw, modeled=nearsat)
        self.assertGreater(old[30:70, 30:70].sum(), 100)
        self.assertEqual(new[30:70, 30:70].sum(), 0)              # nothing greyed at the core
        self.assertTrue(new[:, 136:140].all())                    # the bad column keeps its pad
        self.assertFalse(pl.unmodeledNearMask(data, var, np.zeros((n, n), bool), fw).any())


class PlotSmokeTestCase(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def testFinalComparison(self):
        run = fakeRun.makeRun()
        for kw in (dict(), dict(grayscale=True), dict(panels="residuals", sigma="pixel"),
                   dict(circles="none", excludeNearMask=False)):
            st = plotting.plotFinalComparison(run, view=None, verbose=False, **kw)
            self.assertLess(abs(st["residualOverSigma"]["sigmaIqr"] - 1.0), 0.1)

    def testCrossSectionsAndZooms(self):
        run = fakeRun.makeRun()
        out = plotting.showCrossSections(run, stars=[(1040, 2040)], half=10)
        self.assertEqual(len(out), 1)
        plotting.showStarZooms(run, n=2, half=10)

    def testBrightStarStack(self):
        out = plotting.showBrightStarStack(_run(), minSnr=100, plot=True)
        self.assertGreater(out["n"], 5)
        self.assertLess(np.nanmax(np.abs(out["psfOnly"])), 1e-3)


if __name__ == "__main__":
    unittest.main()
