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

"""Tests of the whole-detector tiling, per-tile products, stitching and plotting.  The per-tile
fit itself is replaced by a perfect synthetic model (it needs the stack); everything around it
-- task splitting, resuming, files, ownership, mosaic -- is the real code."""

import os
import tempfile
import types
import unittest
from unittest import mock

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import lsstStubs  # noqa: E402

lsstStubs.install()

from lsst.meas.crowd import detector as det  # noqa: E402
from lsst.meas.crowd import gaiaCatalog, pipeline, session  # noqa: E402

geom = pipeline.L.geom
plt.show = lambda *a, **k: None


class TileSpecTestCase(unittest.TestCase):
    def testPartition(self):
        for (w, h, t) in ((4072, 4000, 400), (600, 500, 200), (100, 90, 400), (1001, 999, 250)):
            specs = det.tileSpecs(w, h, t, x0=7, y0=3)
            cover = np.zeros((h, w), int)
            for x0, y0, tw, th in specs:
                self.assertLessEqual(tw, t + 1)
                cover[y0 - 3:y0 - 3 + th, x0 - 7:x0 - 7 + tw] += 1
            self.assertTrue(np.all(cover == 1), (w, h, t))      # every pixel in exactly one tile
        self.assertEqual(len(det.tileSpecs(4072, 4000, 400)), 11 * 10)

    def testRegion(self):
        specs = det.tileSpecs(600, 500, 200, region=(250, 330, 0, 50))
        self.assertEqual([s[:2] for s in specs], [(200, 0)])
        self.assertEqual(det.tileName((200, 0, 200, 167)), "x0200_y0000_200x167")


class _Synthetic:
    """A 'detector' with Gaussian stars and a perfect per-tile fit."""

    def __init__(self, shape=(500, 600), n=60, sigma=2.0, half=12, seed=1):
        rng = np.random.default_rng(seed)
        self.shape, self.sigma, self.half = shape, sigma, half
        self.xy = np.column_stack([rng.uniform(5, shape[1] - 5, n), rng.uniform(5, shape[0] - 5, n)])
        self.flux = 10 ** rng.uniform(3.5, 5, n)
        truth = np.zeros(shape)
        for (x, y), f in zip(self.xy, self.flux):
            truth += lsstStubs.place(shape, lsstStubs.gaussianComp(shape, x, y, sigma, half), f)
        self.truth = truth
        self.var = 100 + truth
        self.data = truth + rng.normal(0, 1, shape) * np.sqrt(self.var)
        self.psf = lsstStubs.FakePsf(lsstStubs.gaussianKernel(sigma, half))
        self.calls = []
        bb = geom.Box2I(geom.Point2I(0, 0), geom.Point2I(shape[1] - 1, shape[0] - 1))
        self.visitImage = types.SimpleNamespace(getBBox=lambda: bb)

    def prepareSession(self, visit, detector, cutouts, pad=0, psfSource=None, visitImage=None):
        jobs = {}
        for c in cutouts:
            ib = geom.Box2I(geom.Point2I(c["x0"], c["y0"]),
                            geom.Point2I(c["x0"] + c["w"] - 1, c["y0"] + c["h"] - 1))
            fx0, fy0 = max(c["x0"] - pad, 0), max(c["y0"] - pad, 0)
            fx1 = min(c["x0"] + c["w"] - 1 + pad, self.shape[1] - 1)
            fy1 = min(c["y0"] + c["h"] - 1 + pad, self.shape[0] - 1)
            sl = (slice(fy0, fy1 + 1), slice(fx0, fx1 + 1))
            exp = lsstStubs.FakeExposure(self.data[sl].copy(), self.var[sl].copy(), psf=self.psf,
                                         x0=fx0, y0=fy0)
            jobs[c["name"]] = dict(interior=ib, raw=exp)
        return types.SimpleNamespace(jobs=jobs, gaia=None)

    def runSessionGaia(self, sess, cfg, names=None, plotFinal=False, **kw):
        name = names[0]
        self.calls.append(name)
        job = sess.jobs[name]
        exp = job["raw"]
        bb = exp.getBBox()
        x0, y0 = bb.getMinX(), bb.getMinY()
        shape = exp.image.array.shape
        inside = ((self.xy[:, 0] > x0 - 3) & (self.xy[:, 0] < bb.getMaxX() + 3) &
                  (self.xy[:, 1] > y0 - 3) & (self.xy[:, 1] < bb.getMaxY() + 3))
        xy, fl = self.xy[inside], self.flux[inside]
        comps = [lsstStubs.gaussianComp(shape, x, y, self.sigma, self.half, x0, y0) for x, y in xy]
        model = self.truth[y0:y0 + shape[0], x0:x0 + shape[1]].copy()
        S = pipeline._state(xy, fromGaia=True)
        fit = dict(amps=fl, modelFull=model, good=np.ones(shape, bool), variance=exp.variance.array,
                   comps=comps)
        res = dict(S=S, snaps=[dict(fit=fit)], modelFull=model, amps=fl.copy(), brightCorr={},
                   peaks=[tuple(p) for p in xy])
        run = types.SimpleNamespace(name=name, outerHistory=[dict(exposure=exp, result=res)],
                                    badBitmask=4, interiorBbox=job["interior"], config={})
        return {name: run}


class FitStitchTestCase(unittest.TestCase):
    def setUp(self):
        self.syn = _Synthetic()
        self.tmp = tempfile.TemporaryDirectory()
        self.out = self.tmp.name
        self.patches = [mock.patch.object(session, "prepareSession", self.syn.prepareSession),
                        mock.patch.object(session, "makeConfig", lambda *a, **k: {}),
                        mock.patch.object(gaiaCatalog, "attachGaia", lambda *a, **k: None),
                        mock.patch.object(pipeline, "runSessionGaia", self.syn.runSessionGaia)]
        for p in self.patches:
            p.start()

    def tearDown(self):
        for p in self.patches:
            p.stop()
        self.tmp.cleanup()
        plt.close("all")

    def _fit(self, **kw):
        return det.fitTiles(1, 2, outDir=self.out, tile=200, pad=20,
                            visitImage=self.syn.visitImage, **kw)

    def testTasksResumeAndStitch(self):
        a = self._fit(taskIndex=0, nTasks=2)
        b = self._fit(taskIndex=1, nTasks=2)
        self.assertEqual(len(a) + len(b), 9)                        # 3 x 3 tiles, split in two
        self.assertEqual(len(set(a) & set(b)), 0)
        nCalls = len(self.syn.calls)
        self._fit()                                                 # everything done: no refits
        self.assertEqual(len(self.syn.calls), nCalls)
        m = det.stitchTiles(self.out, verbose=False)
        np.testing.assert_allclose(m["data"], self.syn.data.astype(np.float32), rtol=1e-6)
        np.testing.assert_allclose(m["model"], self.syn.truth.astype(np.float32), rtol=1e-5, atol=1e-3)
        self.assertTrue(np.all(m["tile"] >= 0))
        z = m["resid"] / np.sqrt(m["var"])
        self.assertLess(abs(np.std(z) - 1), 0.05)
        # every star exactly once, with the right flux
        cat = m["catalog"]
        self.assertEqual(len(cat), len(self.syn.xy))
        order = np.lexsort((cat["y"], cat["x"]))
        torder = np.lexsort((self.syn.xy[:, 1], self.syn.xy[:, 0]))
        np.testing.assert_allclose(np.asarray(cat["x"])[order], self.syn.xy[torder, 0])
        np.testing.assert_allclose(np.asarray(cat["psfFlux"])[order], self.syn.flux[torder])
        # files on disk
        self.assertTrue(os.path.exists(os.path.join(self.out, "detector_mosaic.fits")))
        m2 = det.loadMosaic(self.out)
        np.testing.assert_array_equal(m2["tile"], m["tile"])
        self.assertEqual(len(m2["catalog"]), len(cat))

    def testRegionAndGaps(self):
        files = self._fit(region=(0, 150, 0, 150))                  # only the corner tile
        self.assertEqual(len(files), 1)
        m = det.stitchTiles(self.out, verbose=False)
        self.assertEqual(m["data"].shape, (167, 200))

    def testFailingTileDoesNotStopTask(self):
        good = self.syn.runSessionGaia

        def flaky(sess, cfg, names=None, **kw):
            if names[0].startswith("x0200_y0167"):
                raise RuntimeError("boom")
            return good(sess, cfg, names=names, **kw)

        with mock.patch.object(pipeline, "runSessionGaia", flaky):
            files = self._fit()
        self.assertEqual(len(files), 8)
        m = det.stitchTiles(self.out, verbose=False)
        self.assertTrue(np.any(m["tile"] < 0))                      # the gap is marked

    def testPlot(self):
        self._fit()
        m = det.stitchTiles(self.out, verbose=False)
        png = os.path.join(self.out, "view.png")
        st = det.plotDetector(m, out=png)
        self.assertTrue(os.path.getsize(png) > 1000)
        self.assertLess(abs(st["residualSigma"] - 1), 0.05)
        st = det.plotDetector(self.out, view=(100, 300, 50, 250), binning=1, sources=True,
                              grayscale=True)
        self.assertEqual(st["binning"], 1)

    def testTile(self):
        files = self._fit(region=(0, 150, 0, 150))                  # the corner tile
        t = det.loadTile(files[0])
        ix0, iy0, ix1, iy1 = t["interiors"][0]
        self.assertEqual((ix0, iy0), (0, 0))
        self.assertEqual(t["data"].shape[1], ix1 - ix0 + 1 + 20)    # interior + pad on one side
        self.assertEqual(int(np.sum(t["catalog"]["owned"])),
                         int(np.sum((self.syn.xy[:, 0] < ix1 + .5) & (self.syn.xy[:, 1] < iy1 + .5))))
        png = os.path.join(self.out, "tile.png")
        st = det.plotTile(files[0], out=png)
        self.assertTrue(os.path.getsize(png) > 1000)
        self.assertLess(abs(st["residualSigma"] - 1), 0.1)

    def testResidualImages(self):
        from astropy.io import fits
        files = self._fit(region=(0, 150, 0, 150))                  # the corner tile
        with fits.open(files[0]) as hl:
            names = [h.name for h in hl]
            d, m, v = (hl[k].data.astype(float) for k in ("DATA", "MODEL", "VARIANCE"))
            r, z = hl["RESIDUAL"].data, hl["RESIDUAL_SIGMA"].data
        self.assertEqual(names[1:7], ["DATA", "MODEL", "RESIDUAL", "RESIDUAL_SIGMA", "VARIANCE",
                                      "BAD"])
        np.testing.assert_allclose(r, d - m, rtol=1e-5, atol=1e-3)
        ok = np.isfinite(z)
        np.testing.assert_allclose(z[ok], ((d - m) / np.sqrt(v))[ok], rtol=1e-4, atol=1e-4)
        det.stitchTiles(self.out, verbose=False)
        with fits.open(os.path.join(self.out, "detector_mosaic.fits")) as hl:
            self.assertEqual([h.name for h in hl][1:], ["DATA", "MODEL", "RESIDUAL",
                                                        "RESIDUAL_SIGMA", "VARIANCE", "BAD",
                                                        "TILE"])
            self.assertLess(abs(np.nanstd(hl["RESIDUAL_SIGMA"].data) - 1), 0.1)
        self.assertIn("residSigma", det.loadMosaic(self.out))
        self.assertIn("residSigma", det.loadTile(files[0]))

    def testOldTilesGetResiduals(self):
        from astropy.io import fits
        files = self._fit(region=(0, 150, 0, 150))
        with fits.open(files[0]) as hl:                             # a tile from an older version
            old = fits.HDUList([h.copy() for h in hl if h.name not in ("RESIDUAL", "RESIDUAL_SIGMA")])
        old.writeto(files[0], overwrite=True)
        self.assertEqual(det.addTileResiduals(self.out), files)
        self.assertEqual(det.addTileResiduals(self.out), [])         # already there: untouched
        with fits.open(files[0]) as hl:
            self.assertEqual([h.name for h in hl][1:8], ["DATA", "MODEL", "RESIDUAL",
                                                         "RESIDUAL_SIGMA", "VARIANCE", "BAD",
                                                         "CATALOG"])
        self.assertEqual(len(det.stitchTiles(self.out, verbose=False)["catalog"]),
                         int(np.sum((self.syn.xy[:, 0] < 199.5) & (self.syn.xy[:, 1] < 166.5))))

    def testCli(self):
        import importlib.util
        path = os.path.join(os.path.dirname(__file__), "..", "bin.src", "crowdDetector.py")
        spec = importlib.util.spec_from_file_location("crowdDetector", path)
        cli = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cli)
        self._fit()
        self.assertEqual(cli.main(["stitch", "--out", self.out]), 0)
        self.assertEqual(cli.main(["plot", "--out", self.out, "--view", "0", "200", "0", "200"]), 0)
        self.assertTrue(os.path.exists(os.path.join(self.out, "detector_view.png")))


if __name__ == "__main__":
    unittest.main()
