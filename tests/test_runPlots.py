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

"""Tests of the stage plots: PSF stars, backgrounds per pass, rounds, linked source pages."""

import os
import tempfile
import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from astropy.table import Table  # noqa: E402

import lsstStubs  # noqa: E402

lsstStubs.install()

import fakeRun  # noqa: E402
from lsst.meas.crowd import detector as dt  # noqa: E402
from lsst.meas.crowd import runPlots as rp  # noqa: E402
from test_stagePlots import _twoPassRun  # noqa: E402

plt.show = lambda *a, **k: None


class SampleTilesTestCase(unittest.TestCase):
    """Which tiles get the full set of plots."""

    def testCornersCentreAndEdges(self):
        specs = dt.tileSpecs(4072, 4000, 400)
        names = rp.sampleTileNames(specs, 4072, 4000)
        self.assertEqual(len(names), 7)
        self.assertEqual(names[0], dt.tileName(specs[0]))                 # bottom-left corner
        self.assertIn(dt.tileName(specs[-1]), names)                     # top-right corner
        self.assertEqual(len(set(names)), len(names))

    def testOneTileDetector(self):
        specs = dt.tileSpecs(300, 300, 400)
        self.assertEqual(rp.sampleTileNames(specs, 300, 300), [dt.tileName(specs[0])])


class SaveRunPlotsTestCase(unittest.TestCase):
    """saveRunPlots writes every plot of a run and an index describing them."""

    def testWritesEveryStage(self):
        run, bb, bg = _twoPassRun()
        with tempfile.TemporaryDirectory() as d:
            files = rp.saveRunPlots(run, d, verbose=False)
            self.assertTrue(files)
            for f in files:
                self.assertTrue(os.path.exists(f), f)
            names = sorted(os.listdir(d))
            for stage in ("psfModels", "backgrounds", "rounds_outer", "rounds_inner_pass2",
                          "final_comparison", "final_residual_gray"):
                self.assertTrue(any(stage in n for n in names), stage)
            self.assertIn("report.txt", names)
            with open(os.path.join(d, "index.txt")) as fh:
                index = fh.read()
            for n in names:
                if n.endswith(".png"):
                    self.assertIn(n, index)                 # every figure is described
            pngs = [n for n in names if n.endswith(".png")]
            self.assertEqual(pngs, sorted(pngs))             # numbered in pipeline order
        self.assertEqual(plt.get_fignums(), [])              # every figure closed

    def testFailingPlotIsSkipped(self):
        with tempfile.TemporaryDirectory() as d:
            save = rp._Saver(d, verbose=False)

            def good():
                plt.figure()
                plt.plot([0, 1])

            def bad():
                plt.figure()
                raise RuntimeError("boom")

            save("bad", "fails", bad)
            save("good", "works", good)
            files = save.finish("test")
            self.assertEqual(save.failed, ["bad"])
            self.assertTrue(any(f.endswith("good.png") for f in files))
            with open(os.path.join(d, "index.txt")) as fh:
                self.assertIn("(not written) bad: RuntimeError: boom", fh.read())


class DetectorPlotsTestCase(unittest.TestCase):
    """PSF-star stamps and the detector-level PSF plots."""

    def testStampsAndDetectorPlots(self):
        run = fakeRun.makeRun(seed=5)
        exp = run.outerHistory[0]["exposure"]
        xy = np.asarray(run.outerHistory[0]["result"]["S"]["xy"], float)
        self.assertEqual(rp.showPsfStarStamps(exp, xy, half=10, maxStamps=5), 5)
        self.assertEqual(rp.showPsfStarStamps(exp, np.array([[-50.0, -50.0]]), half=10), 0)
        plt.close("all")
        cat = Table(dict(x=xy[:, 0], y=xy[:, 1], calib_psf_candidate=np.ones(len(xy), bool),
                         calib_psf_used=np.arange(len(xy)) % 2 == 0))
        with tempfile.TemporaryDirectory() as d:
            files = rp.saveDetectorPlots(exp, d, catalog=cat, verbose=False)
            names = [os.path.basename(f) for f in files]
            self.assertTrue(any("psfStars_detector_stamps" in n for n in names))
            self.assertTrue(any("psfStars_detector.png" in n or "psfStars_detector_" in n
                                for n in names))
            with tempfile.TemporaryDirectory() as d2:          # no catalogue: kernels only
                files = rp.saveDetectorPlots(exp, d2, verbose=False)
                self.assertTrue(any("psfKernel_detector" in f for f in files))


class FitTilesPlotsTestCase(unittest.TestCase):
    """fitTiles(plots=...) writes the plots of the chosen tiles, refitting tiles without them."""

    def setUp(self):
        from unittest import mock
        from lsst.meas.crowd import gaiaCatalog, pipeline, session
        from test_detector import _Synthetic
        self.syn = _Synthetic()
        self.tmp = tempfile.TemporaryDirectory()
        self.out = self.tmp.name
        self.runPlots, self.detPlots = [], []

        def fakeRunPlots(run, outDir, **kw):
            os.makedirs(outDir, exist_ok=True)
            open(os.path.join(outDir, "index.txt"), "w").close()
            self.runPlots.append(os.path.basename(outDir))
            return []

        def fakeDetPlots(visitImage, outDir, catalog=None, session=None, **kw):
            os.makedirs(outDir, exist_ok=True)
            open(os.path.join(outDir, "index.txt"), "w").close()
            self.detPlots.append(outDir)
            return []

        self.patches = [mock.patch.object(session, "prepareSession", self.syn.prepareSession),
                        mock.patch.object(session, "makeConfig", lambda *a, **k: {}),
                        mock.patch.object(gaiaCatalog, "attachGaia", lambda *a, **k: None),
                        mock.patch.object(pipeline, "runSessionGaia", self.syn.runSessionGaia),
                        mock.patch.object(rp, "saveRunPlots", fakeRunPlots),
                        mock.patch.object(rp, "saveDetectorPlots", fakeDetPlots),
                        mock.patch.object(dt, "_loadPsfStarCatalog", lambda *a: None)]
        for p in self.patches:
            p.start()

    def tearDown(self):
        for p in self.patches:
            p.stop()
        self.tmp.cleanup()
        plt.close("all")

    def _fit(self, **kw):
        return dt.fitTiles(1, 2, outDir=self.out, tile=200, pad=20,
                           visitImage=self.syn.visitImage, **kw)

    def testSampleTilesAndDetectorPlots(self):
        self._fit(plots="sample")
        specs = dt.tileSpecs(600, 500, 200)
        expected = rp.sampleTileNames(specs, 600, 500)
        self.assertEqual(sorted(self.runPlots), sorted(expected))
        self.assertEqual(self.detPlots, [os.path.join(self.out, "plots", "detector")])
        nCalls = len(self.syn.calls)
        self._fit(plots="sample")                                  # all there: nothing redone
        self.assertEqual(len(self.syn.calls), nCalls)
        self.assertEqual(len(self.detPlots), 1)

    def testPlotsForAnAlreadyFittedTile(self):
        self._fit()                                                # no plots
        self.assertEqual(self.runPlots, [])
        self.assertEqual(self.detPlots, [])
        nCalls = len(self.syn.calls)
        name = dt.tileName(dt.tileSpecs(600, 500, 200)[4])
        self._fit(plots=[name])                                    # refit that tile only
        self.assertEqual(self.runPlots, [name])
        self.assertEqual(len(self.syn.calls), nCalls + 1)


if __name__ == "__main__":
    unittest.main()
