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

"""Tests of the photometry helpers (magnitudes, apertures, Fisher errors)."""

import types
import unittest

import numpy as np

import lsstStubs

lsstStubs.install()

from lsst.meas.crowd import photometry as ph  # noqa: E402


class MagnitudeTestCase(unittest.TestCase):
    def testFluxToMag(self):
        m = ph.fluxToMag([3631e9, 1.0, 0.0, -5.0])
        self.assertAlmostEqual(m[0], 0.0, places=2)          # 3631 Jy = AB 0
        self.assertAlmostEqual(m[1], 31.4, places=12)
        self.assertTrue(np.isnan(m[2]) and np.isnan(m[3]))

    def testMagError(self):
        m, e = ph.fluxToMag([100.0], [1.0])
        self.assertAlmostEqual(e[0], 2.5 / np.log(10) * 0.01, places=12)

    def testZeroPoint(self):
        m = ph.fluxToMag([10.0], zeroPoint=25.0)
        self.assertAlmostEqual(m[0], 22.5, places=12)

    def testZeroPointFromRun(self):
        img = np.zeros((5, 5))
        run = types.SimpleNamespace(outerHistory=[dict(exposure=lsstStubs.FakeExposure(img, img + 1))])
        zp, how = ph.zeroPointFromRun(run)
        self.assertEqual(zp, ph.AB_NJY_ZERO_POINT)
        run.outerHistory[0]["exposure"] = lsstStubs.FakeExposure(img, img + 1, photoScale=10.0)
        zp, how = ph.zeroPointFromRun(run)
        self.assertAlmostEqual(zp, 31.4 - 2.5, places=12)    # 1 count = 10 nJy
        self.assertEqual(how, "PhotoCalib")


class ApertureTestCase(unittest.TestCase):
    def testPixelFractionsArea(self):
        yy, xx = np.mgrid[-10:11, -10:11]
        for r in (2.0, 3.7, 6.25):
            f = ph._pixelFractions(xx + 0.3, yy - 0.2, r)
            self.assertAlmostEqual(f.sum(), np.pi * r * r, delta=0.02 * np.pi * r * r)

    def testApertureFluxOfGaussian(self):
        shape, s = (60, 60), 2.0
        c = lsstStubs.gaussianComp(shape, 30.3, 29.6, s, 25)
        img = lsstStubs.place(shape, c, 1000.0)
        var = np.ones(shape)
        good = np.ones(shape, bool)
        fl, err, cov = ph.apertureFlux(img, var, good, 30.3, 29.6, 15.0)
        self.assertAlmostEqual(fl, 1000.0, delta=0.5)
        self.assertAlmostEqual(cov, 1.0)
        self.assertAlmostEqual(err, np.sqrt(np.pi * 15.0 ** 2), delta=0.03 * err)
        # 1-sigma radius of a 2-d Gaussian holds 1 - exp(-1/2) = 39.3%
        fl, _, _ = ph.apertureFlux(img, var, good, 30.3, 29.6, s)
        self.assertAlmostEqual(fl / 1000.0, 1 - np.exp(-0.5), delta=0.01)

    def testCoverage(self):
        shape = (40, 40)
        good = np.ones(shape, bool)
        good[:, :20] = False
        fl, err, cov = ph.apertureFlux(np.ones(shape), np.ones(shape), good, 19.5, 20.0, 5.0)
        self.assertAlmostEqual(cov, 0.5, delta=0.05)

    def testOffImage(self):
        fl, err, cov = ph.apertureFlux(np.ones((10, 10)), np.ones((10, 10)),
                                       np.ones((10, 10), bool), 100.0, 100.0, 3.0)
        self.assertTrue(np.isnan(fl))
        self.assertEqual(cov, 0.0)

    def testApertureCorrection(self):
        k = lsstStubs.gaussianKernel(2.0, 20)
        self.assertAlmostEqual(ph.apertureCorrection(k, 2.0), 1 - np.exp(-0.5), delta=0.01)
        self.assertAlmostEqual(ph.apertureCorrection(k, 30.0), 1.0, places=6)


class FisherErrorTestCase(unittest.TestCase):
    def testIsolatedSource(self):
        """sigma(amp) = 1 / sqrt(sum P^2 / var) for one source."""
        shape = (50, 50)
        c = lsstStubs.gaussianComp(shape, 25.0, 25.0, 2.0, 15)
        var = np.full(shape, 4.0)
        err = ph.fisherFluxErrors([c], var, np.ones(shape, bool))
        self.assertAlmostEqual(err[0], 1 / np.sqrt(np.sum(c["stamp"] ** 2 / 4.0)), places=10)

    def testBlendHasLargerErrors(self):
        shape = (60, 60)
        var = np.ones(shape)
        good = np.ones(shape, bool)
        iso = ph.fisherFluxErrors([lsstStubs.gaussianComp(shape, 30, 30, 2.0, 15)], var, good)
        pair = ph.fisherFluxErrors([lsstStubs.gaussianComp(shape, 30, 30, 2.0, 15),
                                    lsstStubs.gaussianComp(shape, 32, 30, 2.0, 15)], var, good)
        far = ph.fisherFluxErrors([lsstStubs.gaussianComp(shape, 15, 15, 2.0, 12),
                                   lsstStubs.gaussianComp(shape, 45, 45, 2.0, 12)], var, good)
        self.assertGreater(pair[0], 1.5 * iso[0])
        self.assertAlmostEqual(far[0], iso[0], places=8)

    def testMissingComponent(self):
        shape = (30, 30)
        err = ph.fisherFluxErrors([None, lsstStubs.gaussianComp(shape, 15, 15, 2.0, 10)],
                                  np.ones(shape), np.ones(shape, bool))
        self.assertTrue(np.isinf(err[0]))
        self.assertTrue(np.isfinite(err[1]))

    def testErrorsMatchScatter(self):
        """The Fisher error is the actual scatter of least-squares amplitudes."""
        shape = (40, 40)
        comps = [lsstStubs.gaussianComp(shape, 20, 20, 2.0, 12),
                 lsstStubs.gaussianComp(shape, 22.5, 21, 2.0, 12)]
        A = np.array([lsstStubs.place(shape, c).ravel() for c in comps]).T
        truth = A @ np.array([500.0, 300.0])
        rng = np.random.default_rng(3)
        amps = np.array([np.linalg.lstsq(A, truth + rng.normal(0, 1, truth.size), rcond=None)[0]
                         for _ in range(400)])
        err = ph.fisherFluxErrors(comps, np.ones(shape), np.ones(shape, bool))
        np.testing.assert_allclose(amps.std(axis=0), err, rtol=0.12)


class MeasurePhotometryTestCase(unittest.TestCase):
    def testSyntheticRun(self):
        import fakeRun
        run = fakeRun.makeRun(nSide=5, seed=4)
        t = ph.measurePhotometry(run, verbose=False)
        self.assertEqual(len(t), 25)
        self.assertEqual(t.meta["zeroPointSource"], "image units assumed to be nJy")
        # PSF fluxes agree with the truth within their errors
        pull = (t["psfFlux"] - run.truthFlux) / t["psfFluxErr"]
        self.assertLess(abs(np.mean(pull)), 1.0)
        self.assertLess(np.std(pull), 1.6)
        # aperture photometry (corrected to total) agrees with PSF photometry for bright stars
        bright = t["psfSnr"] > 50
        self.assertLess(abs(np.median(t["apMinusPsf"][bright])), 0.01)
        np.testing.assert_allclose(t["psfMag"], 31.4 - 2.5 * np.log10(t["psfFlux"]))

    def testPhotoCalib(self):
        import fakeRun
        run = fakeRun.makeRun(nSide=2, photoScale=2.0)
        t = ph.measurePhotometry(run, verbose=False)
        np.testing.assert_allclose(t["psfMag"], 31.4 - 2.5 * np.log10(2.0 * t["psfFlux"]))


class GaiaComparisonTestCase(unittest.TestCase):
    """Matching to Gaia and the flux / position differences."""

    def _tables(self):
        from astropy.table import Table
        dec0 = 10.0
        gaia = Table(dict(source_id=np.arange(100, 106), x=[10., 20., 30., 40., 25., 500.],
                          y=[10., 20., 30., 40., 15., 500.], ra=np.full(6, 50.0),
                          dec=np.full(6, dec0), phot_g_mean_mag=[15., 16., 17., 18., 19., 20.],
                          phot_bp_mean_mag=[15.5, 16.5, 17.5, 18.5, 19.5, 20.5],
                          phot_rp_mean_mag=[14.8, 15.8, 16.8, 17.8, 18.8, 19.8]))
        gaia.meta["gaiaMagSystem"] = "Vega (Gaia archive)"
        G = np.array([15., 16., 17., 18., 21.])
        phot = Table(dict(x=[10.1, 20.1, 30.1, 40.4, 33.0], y=[10., 20., 30., 40.2, 12.0],
                          gaiaIndex=[0, 1, 2, -1, -1], psfMag=G + 0.3, psfSnr=np.full(5, 100.),
                          extended=np.zeros(5, bool),
                          ra=50.0 + np.array([2, 2, 2, 0, 0]) / 3.6e6 / np.cos(np.radians(dec0)),
                          dec=dec0 + np.array([1, 1, 1, 0, 0]) / 3.6e6))
        phot["psfFlux"] = 10 ** (-0.4 * (phot["psfMag"] - ph.AB_NJY_ZERO_POINT))
        phot.meta["zeroPoint"] = ph.AB_NJY_ZERO_POINT
        return phot, gaia

    def testMatchingAndDifferences(self):
        phot, gaia = self._tables()
        t = ph.compareWithGaia(phot, gaia, verbose=False)
        np.testing.assert_array_equal(t["gaiaMatch"], [True, True, True, True, False])
        self.assertEqual(list(t["gaiaMatchType"]), ["linked", "linked", "linked", "position", ""])
        np.testing.assert_array_equal(t["gaiaSourceId"], [100, 101, 102, 103, -1])
        np.testing.assert_allclose(t["gaiaDx"][:3], 0.1, atol=1e-9)
        np.testing.assert_allclose(t["gaiaDra"][:3], 2.0, atol=1e-3)       # mas, with cos(dec)
        np.testing.assert_allclose(t["gaiaDdec"][:3], 1.0, atol=1e-3)
        np.testing.assert_allclose(t["gaiaSep"][:3], np.sqrt(5), atol=1e-3)
        self.assertAlmostEqual(t.meta["GAIAOFF"], 0.3, places=9)            # psfMag - G
        np.testing.assert_allclose(t["gaiaDmagRel"][:4], 0.0, atol=1e-9)
        np.testing.assert_allclose(t["gaiaFluxRatio"][:4], 1.0, rtol=1e-9)
        np.testing.assert_allclose(t["gaiaBpRp"][:4], 0.7, atol=1e-9)
        self.assertEqual(t.meta["GAIAUNM"], 1)            # star at (25, 15) has no fitted source
        self.assertTrue(np.isnan(t["gaiaG"][4]))
        self.assertNotIn("gaiaMatch", phot.colnames)    # input not modified

    def testRadiusAndNoSky(self):
        phot, gaia = self._tables()
        del phot["ra"], phot["dec"]
        t = ph.compareWithGaia(phot, gaia, matchRadiusPx=0.3, verbose=False)
        self.assertFalse(t["gaiaMatch"][3])               # 0.45 px away
        self.assertNotIn("gaiaDra", t.colnames)

    def testFromSyntheticRun(self):
        import fakeRun
        from astropy.table import Table
        run = fakeRun.makeRun(seed=7)
        xy = run.outerHistory[-1]["result"]["S"]["xy"]
        G = ph.AB_NJY_ZERO_POINT - 2.5 * np.log10(run.truthFlux) - 0.2       # offset +0.2
        gaia = Table(dict(source_id=np.arange(len(xy)), x=xy[:, 0] + 0.05, y=xy[:, 1],
                          ra=np.zeros(len(xy)), dec=np.zeros(len(xy)), phot_g_mean_mag=G))
        t = ph.measurePhotometry(run, verbose=False, gaia=gaia)
        self.assertTrue(np.all(t["gaiaMatch"]))
        self.assertAlmostEqual(t.meta["GAIAOFF"], 0.2, delta=0.01)
        self.assertLess(np.std(t["gaiaFluxRatio"]), 0.02)

    def testRunGaiaUsedByDefault(self):
        import fakeRun
        from astropy.table import Table
        run = fakeRun.makeRun(seed=7)
        xy = run.outerHistory[-1]["result"]["S"]["xy"]
        run.gaia = Table(dict(source_id=np.arange(len(xy)), x=xy[:, 0], y=xy[:, 1],
                              ra=np.zeros(len(xy)), dec=np.zeros(len(xy)),
                              phot_g_mean_mag=np.full(len(xy), 18.0)))
        self.assertIn("gaiaDmag", ph.measurePhotometry(run, verbose=False).colnames)
        self.assertNotIn("gaiaDmag", ph.measurePhotometry(run, verbose=False, gaia=False).colnames)

    def testWritePhotometry(self):
        import os
        import tempfile
        import fakeRun
        from astropy.io import fits
        from astropy.table import Table
        from lsst.meas.crowd import session
        run = fakeRun.makeRun(seed=7)
        xy = run.outerHistory[-1]["result"]["S"]["xy"]
        run.gaia = Table(dict(source_id=np.arange(len(xy)), x=xy[:, 0], y=xy[:, 1],
                              ra=np.zeros(len(xy)), dec=np.zeros(len(xy)),
                              phot_g_mean_mag=np.full(len(xy), 18.0)))
        run.photometry = None
        with tempfile.TemporaryDirectory() as d:
            path = session.writePhotometry(run, d, verbose=False)
            self.assertTrue(os.path.exists(os.path.join(d, "photometry.csv")))
            back = Table.read(path)
            hdr = fits.getheader(path, 1)
        self.assertEqual(len(back), len(xy))
        self.assertIn("gaiaDmagRel", back.colnames)
        self.assertEqual(hdr["ZP"], ph.AB_NJY_ZERO_POINT)        # short FITS keywords
        self.assertIn("GAIAOFF", hdr)
        self.assertIsNotNone(run.photometry)                     # kept on the run

    def testFitsMeta(self):
        m = ph.fitsMeta(dict(zeroPoint=31.4, apertureRadiusPx=3.0, GAIAOFF=np.nan,
                             GAIASYS="Vega", obj=[1, 2]))
        self.assertEqual(m, dict(ZP=31.4, APRAD=3.0, GAIASYS="Vega"))

    def testPlot(self):
        import os
        import tempfile
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from lsst.meas.crowd import plotting
        phot, gaia = self._tables()
        t = ph.compareWithGaia(phot, gaia, verbose=False)
        with tempfile.TemporaryDirectory() as d:
            png = os.path.join(d, "gaia.png")
            out = plotting.plotGaiaComparison(t, out=png)
            self.assertGreater(os.path.getsize(png), 1000)
        self.assertEqual(out["n"], 4)
        self.assertAlmostEqual(out["dx"], 2.0 * 3 / 4 + 0.0 / 4, delta=1.0)
        with self.assertRaises(ValueError):
            plotting.plotGaiaComparison(phot)
        plt.close("all")


if __name__ == "__main__":
    unittest.main()
