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
import types
import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from astropy.table import Table  # noqa: E402

import lsstStubs  # noqa: E402

lsstStubs.install()

import fakeRun  # noqa: E402
from lsst.meas.crowd import stagePlots as sp  # noqa: E402

plt.show = lambda *a, **k: None


def _twoPassRun():
    """A fake run with two outer passes, background images and three inner snapshots."""
    run = fakeRun.makeRun(seed=3)
    final = run.outerHistory[0]
    exp, res = final["exposure"], final["result"]
    data, model = exp.image.array, res["modelFull"]
    n = len(res["S"]["xy"])
    fit = res["snaps"][0]["fit"]
    fit.update(groupLabels=np.zeros(n, int), groupChi2={0: 1.0}, nGroups=1, nAtZero=0)
    x0, y0 = exp.getBBox().getMinX(), exp.getBBox().getMinY()
    for c, a in zip(fit["comps"], fit["amps"]):      # the fields the real joint fit adds
        c.update(amp=a, sig=c["stamp"] > 1e-3,
                 bboxParent=(x0 + c["xs"].start, x0 + c["xs"].stop - 1,
                             y0 + c["ys"].start, y0 + c["ys"].stop - 1))
    # inner snapshots: the source list growing (first half, then all, then all again)
    half = fit["amps"].copy()
    half[n // 2:] = 0.0
    comps = fit["comps"]
    mHalf = sum(lsstStubs.place(data.shape, c, a) for c, a in zip(comps, half))
    snaps = [dict(it=0, threshold=8.0, residual=data - mHalf, mergedArr=res["S"]["xy"][:n // 2],
                  fit=dict(amps=half[:n // 2])),
             dict(it=1, threshold=6.5, residual=data - model, mergedArr=res["S"]["xy"],
                  fit=dict(amps=fit["amps"])),
             dict(it=2, threshold=5.0, residual=data - model, mergedArr=res["S"]["xy"],
                  fit=fit, prevAmp=fit["amps"].copy(), isNew=np.zeros(n, bool),
                  badMask=np.zeros(data.shape, bool))]
    res["snaps"] = snaps
    bb = exp.getBBox()
    ny, nx = data.shape
    bg2 = np.full((ny, nx), 3.0) + np.linspace(0, 1, nx)[None, :]
    first = dict(final, bgImage=None)
    second = dict(final, bgImage=types.SimpleNamespace(array=bg2))
    run.outerHistory = [first, second]
    run.config = dict(satZoneFwhm=0)
    return run, bb, bg2


class PsfStarsTestCase(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def testCatalogFlags(self):
        cat = Table(dict(slot_Centroid_x=[1., 2., 3., 4.], slot_Centroid_y=[5., 6., 7., 8.],
                         calib_psf_candidate=[True, True, True, False],
                         calib_psf_used=[True, False, True, False]))
        st = sp.psfStarPositions(catalog=cat)
        self.assertEqual(len(st["candidates"]), 3)
        np.testing.assert_array_equal(st["used"], [[1., 5.], [3., 7.]])

    def testSession(self):
        sess = types.SimpleNamespace(starCat=Table(dict(x=[10., 20.], y=[30., 40.])),
                                     psfUsedXy=np.array([[10., 30.]]))
        st = sp.psfStarPositions(session=sess)
        self.assertEqual(len(st["candidates"]), 2)
        self.assertEqual(len(st["used"]), 1)
        with self.assertRaises(ValueError):           # psfSource='visit': no star catalogue
            sp.psfStarPositions(session=types.SimpleNamespace(starCat=None))
        with self.assertRaises(ValueError):
            sp.psfStarPositions()

    def testRunStarsAreBrightSources(self):
        run, _, _ = _twoPassRun()
        st = sp.psfStarPositions(run=run)
        xy = run.outerHistory[-1]["result"]["S"]["xy"]
        self.assertGreater(len(st["used"]), 0)
        self.assertLessEqual(len(st["used"]), len(st["candidates"]))
        for p in st["used"]:                          # every used star is one of the sources
            self.assertLess(np.min(np.hypot(*(xy - p).T)), 1e-9)

    def testPlotKernel(self):
        run, _, _ = _twoPassRun()
        exp = run.outerHistory[-1]["exposure"]
        with tempfile.TemporaryDirectory() as d:
            png = os.path.join(d, "psf.png")
            info = sp.showPsfStars(exp, sp.psfStarPositions(run=run), out=png)
            self.assertGreater(os.path.getsize(png), 1000)
        self.assertAlmostEqual(info["fwhm"], 2.355 * 2.0, delta=0.15)   # fake PSF: sigma 2 px
        self.assertAlmostEqual(info["kernel"].sum(), 1.0, places=6)


class BackgroundAndRoundsTestCase(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def testBackgrounds(self):
        run, bb, bg2 = _twoPassRun()
        ny, nx = bg2.shape
        initial = types.SimpleNamespace(array=np.full((ny, nx), 2.0))
        bgs = sp.showBackgrounds(run, initial=initial)
        self.assertEqual(len(bgs), 2)
        np.testing.assert_array_equal(bgs[1], bg2)
        bgs = sp.showBackgrounds(run)                  # visit image: nothing subtracted first
        self.assertIsNone(bgs[0])

    def testOuterRounds(self):
        run, _, _ = _twoPassRun()
        rows = sp.showRounds(run, level="outer")
        self.assertEqual(len(rows), 2)
        self.assertAlmostEqual(rows[-1]["chiStd"], 1.0, delta=0.1)      # perfect model: noise

    def testInnerRounds(self):
        run, _, _ = _twoPassRun()
        rows = sp.showRounds(run, level="inner", view=(1010, 1100, 2010, 2100))
        self.assertEqual([r["nSources"] for r in rows], [8, 16, 16])
        self.assertGreater(rows[0]["chiStd"], rows[-1]["chiStd"])      # missing sources show
        self.assertEqual(len(sp.showRounds(run, level="inner", maxRows=2)), 2)
        with self.assertRaises(ValueError):
            sp.showRounds(run, level="tiles")


class LinkedPagesTestCase(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def testColoursMatch(self):
        run, _, _ = _twoPassRun()
        xy = run.outerHistory[-1]["result"]["S"]["xy"]
        out = sp.showLinkedSourcePages(run, sources=[0, tuple(xy[5] + 0.3)], viewHalf=12)
        self.assertEqual(out["pages"], [[0, 5]])
        self.assertEqual(set(out["colors"]), {0, 5})
        # the overview boxes use exactly the page colours
        fig = plt.figure(1)
        edge = {matplotlib.colors.to_hex(p.get_edgecolor()) for a in fig.axes for p in a.patches
                if isinstance(p, matplotlib.patches.Rectangle)}
        self.assertEqual(edge, {matplotlib.colors.to_hex(c) for c in out["colors"].values()})

    def testDefaultSelectionAndSlimSnapshot(self):
        run, _, _ = _twoPassRun()
        out = sp.showLinkedSourcePages(run, sourcesPerPage=4, showPages=False)
        self.assertEqual(len(out["pages"]), 1)
        amps = run.outerHistory[-1]["result"]["snaps"][-1]["fit"]["amps"]
        self.assertEqual(out["pages"][0], list(np.argsort(-amps)[:4]))    # brightest first
        with self.assertRaises(ValueError):              # older snapshot: stamps dropped
            sp.showLinkedSourcePages(run, snap=0)

    def testOutDirWritesEveryPageAndClosesFigures(self):
        run, _, _ = _twoPassRun()
        with tempfile.TemporaryDirectory() as d:
            out = sp.showLinkedSourcePages(run, sourcesPerPage=4, maxPages=None, outDir=d)
            self.assertEqual(len(out["pages"]), 4)                     # 16 sources / 4
            self.assertEqual(len(out["files"]), 8)                     # overview + page each
            self.assertTrue(all(os.path.getsize(f) > 1000 for f in out["files"]))
        self.assertEqual(plt.get_fignums(), [])                        # nothing kept in memory
        # each overview only boxes its own page's sources
        self.assertEqual(sum(len(p) for p in out["pages"]), 16)


class ComparePassesTestCase(unittest.TestCase):
    """Sources lost / gained between two outer passes, with the reason."""

    def setUp(self):
        from lsst.meas.crowd import pipeline
        base = fakeRun.makeRun(seed=4)
        h = base.outerHistory[0]
        exp, fit = h["exposure"], h["result"]["snaps"][0]["fit"]
        xy, amps, comps = h["result"]["S"]["xy"], fit["amps"], fit["comps"]
        shape = exp.image.array.shape
        self.xy = xy
        self.fake = (xy[0] + xy[1]) / 2 + (7.0, 0.0)         # an empty spot between stars
        fakeComp = lsstStubs.gaussianComp(shape, *self.fake, 2.0, 15, 1000, 2000)

        def onePass(idx, extra=(), rejects=()):
            pxy = np.vstack([xy[idx]] + [np.atleast_2d(p) for p, _, _ in extra])
            pamp = np.concatenate([amps[idx], [a for _, a, _ in extra]])
            model = sum(lsstStubs.place(shape, comps[j], amps[j]) for j in idx)
            for _, a, c in extra:
                model = model + lsstStubs.place(shape, c, a)
            S = pipeline._state(pxy, fromGaia=False)
            res = dict(S=S, snaps=[dict(fit=dict(amps=pamp))], modelFull=model)
            return dict(exposure=exp, psf=h["psf"], result=res,
                        rejects=np.asarray(rejects, float).reshape(-1, 2))

        n = len(xy)
        aIdx = [j for j in range(n) if j != 7]                 # pass a misses star 7 ...
        bIdx = [j for j in range(n) if j not in (3, 9)]        # pass b drops stars 3 and 9
        self.run = types.SimpleNamespace(
            outerHistory=[onePass(aIdx, extra=[(self.fake, 3000.0, fakeComp)]),   # ... has a fake
                          onePass(bIdx, rejects=[xy[9]])],   # star 9: removed by the halo check
            badBitmask=4)

    def tearDown(self):
        plt.close("all")

    def _row(self, t, kind, p):
        d = np.hypot(t["x"] - p[0], t["y"] - p[1])
        sel = (t["kind"] == kind) & (d < 0.5)
        self.assertEqual(int(sel.sum()), 1, (kind, p))
        return t[sel][0]

    def testVerdicts(self):
        t = sp.comparePasses(self.run, plot=False)
        self.assertEqual(int(np.sum(t["kind"] == "lost")), 3)         # stars 3, 9 and the fake
        self.assertEqual(int(np.sum(t["kind"] == "gained")), 1)       # star 7
        r = self._row(t, "lost", self.xy[3])
        self.assertTrue(r["verdict"].startswith("LOST"))
        self.assertAlmostEqual(r["leftFrac"], 1.0, delta=0.15)       # all its flux is in the residual
        self.assertIn("blocked", self._row(t, "lost", self.xy[9])["verdict"])
        r = self._row(t, "lost", self.fake)
        self.assertIn("spurious", r["verdict"])
        self.assertLess(abs(r["leftFrac"]), 0.2)
        self.assertIn("missing from the residual", self._row(t, "gained", self.xy[7])["verdict"])

    def testPlotAndErrors(self):
        with tempfile.TemporaryDirectory() as d:
            png = os.path.join(d, "cmp.png")
            sp.comparePasses(self.run, out=png)
            self.assertGreater(os.path.getsize(png), 1000)
        with self.assertRaises(ValueError):
            sp.comparePasses(types.SimpleNamespace(outerHistory=self.run.outerHistory[:1]))


if __name__ == "__main__":
    unittest.main()
