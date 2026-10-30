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

"""Tests of the per-source cutout data and plots."""

import os
import tempfile
import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import lsstStubs  # noqa: E402

lsstStubs.install()

from lsst.meas.crowd import sourcePlots as sp  # noqa: E402
import fakeRun  # noqa: E402


class CutoutDataTestCase(unittest.TestCase):
    def setUp(self):
        shape = (80, 90)
        self.c = lsstStubs.gaussianComp(shape, 140.3, 235.7, 2.0, 20, x0=100, y0=200)
        self.model = lsstStubs.place(shape, self.c, 5000.0)
        self.var = 25.0 + self.model
        self.bad = np.zeros(shape, bool)

    def testPerfectModel(self):
        d = sp.sourceCutoutData(self.model, self.model, self.var, self.bad, 140.3, 235.7,
                                x0=100, y0=200, half=10)
        self.assertEqual(d["window"], (130, 150, 226, 246))
        self.assertEqual(d["data"].shape, (21, 21))
        for key in ("xcut", "ycut", "radial"):
            self.assertAlmostEqual(np.nanmax(np.abs(d[key]["resid"])), 0.0)
            self.assertEqual(d["chi2"][key[0]], 0.0)
        # the x-cut peaks at the source
        self.assertLess(abs(d["xcut"]["offset"][np.nanargmax(d["xcut"]["data"])]), 0.5)

    def testNoiseOnlyChi2(self):
        rng = np.random.default_rng(0)
        data = self.model + rng.normal(0, 1, self.model.shape) * np.sqrt(self.var)
        d = sp.sourceCutoutData(data, self.model, self.var, self.bad, 140.3, 235.7, 100, 200,
                                half=15, width=3)
        self.assertLess(abs(d["chi2"]["x"] - 1), 0.8)
        self.assertLess(abs(d["chi2"]["r"] - 1), 0.8)

    def testBadPixelsIgnored(self):
        bad = self.bad.copy()
        bad[35, 30:50] = True                # the source's row
        data = self.model.copy()
        data[bad] = 1e9
        d = sp.sourceCutoutData(data, self.model, self.var, bad, 140.3, 235.7, 100, 200,
                                half=10, width=3)
        self.assertLess(np.nanmax(np.abs(d["xcut"]["resid"])), 1e-6)

    def testEdge(self):
        d = sp.sourceCutoutData(self.model, self.model, self.var, self.bad, 101.0, 201.0,
                                100, 200, half=10)
        self.assertEqual(d["window"][:3:2], (100, 200))


class SelectTestCase(unittest.TestCase):
    def testSelections(self):
        phot = {"psfSnr": np.array([5, 50, 100, 30, 200.0])}
        chi2 = np.array([9, 1.0, 3.0, 2.0, 1.5])
        self.assertEqual(sp.selectSources(phot, "worst", 2, chi2=chi2), [2, 3])  # S/N 5 excluded
        self.assertEqual(sp.selectSources(phot, "brightest", 2), [4, 2])
        self.assertEqual(sp.selectSources(phot, "faintest", 1), [3])
        self.assertEqual(len(sp.selectSources(phot, "random", 3)), 3)
        xy = np.array([[0, 0], [10, 10], [20, 20], [30, 30], [40, 40.0]])
        self.assertEqual(sp.selectSources(phot, positions=[(19, 21), (1, 1)], xy=xy), [2, 0])
        with self.assertRaises(ValueError):
            sp.selectSources(phot, "nonsense", 2)


class PlotSmokeTestCase(unittest.TestCase):
    def testPlotAndPdf(self):
        run = fakeRun.makeRun()
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "cutouts.pdf")
            rows = sp.plotSourceCutouts(run, how="brightest", n=3, perPage=2, pdfPath=path,
                                        show=False)
            self.assertTrue(os.path.getsize(path) > 1000)
        self.assertEqual(len(rows), 3)
        self.assertTrue(np.all(np.isfinite(rows["psfMag"])))
        plt.close("all")


if __name__ == "__main__":
    unittest.main()
