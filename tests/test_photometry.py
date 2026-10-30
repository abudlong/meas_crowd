# This file is part of meas_crowd.
#
# Developed for the LSST Data Management System.
# See the COPYRIGHT file at the top-level directory of this distribution.

"""Tests of the photometry helpers (magnitudes, apertures, Fisher errors)."""

import types
import unittest

import numpy as np

import lsstStubs

lsstStubs.install()

from lsst.meas.crowd import photometry as ph  # noqa: E402


class MagnitudeTestCase(unittest.TestCase):
    def testFluxToMag(self):
        m = ph.flux_to_mag([3631e9, 1.0, 0.0, -5.0])
        self.assertAlmostEqual(m[0], 0.0, places=2)          # 3631 Jy = AB 0
        self.assertAlmostEqual(m[1], 31.4, places=12)
        self.assertTrue(np.isnan(m[2]) and np.isnan(m[3]))

    def testMagError(self):
        m, e = ph.flux_to_mag([100.0], [1.0])
        self.assertAlmostEqual(e[0], 2.5 / np.log(10) * 0.01, places=12)

    def testZeroPoint(self):
        m = ph.flux_to_mag([10.0], zero_point=25.0)
        self.assertAlmostEqual(m[0], 22.5, places=12)

    def testZeroPointFromRun(self):
        img = np.zeros((5, 5))
        run = types.SimpleNamespace(outer_history=[dict(exposure=lsstStubs.FakeExposure(img, img + 1))])
        zp, how = ph.zero_point_from_run(run)
        self.assertEqual(zp, ph.AB_NJY_ZERO_POINT)
        run.outer_history[0]["exposure"] = lsstStubs.FakeExposure(img, img + 1, photo_scale=10.0)
        zp, how = ph.zero_point_from_run(run)
        self.assertAlmostEqual(zp, 31.4 - 2.5, places=12)    # 1 count = 10 nJy
        self.assertEqual(how, "PhotoCalib")


class ApertureTestCase(unittest.TestCase):
    def testPixelFractionsArea(self):
        yy, xx = np.mgrid[-10:11, -10:11]
        for r in (2.0, 3.7, 6.25):
            f = ph._pixel_fractions(xx + 0.3, yy - 0.2, r)
            self.assertAlmostEqual(f.sum(), np.pi * r * r, delta=0.02 * np.pi * r * r)

    def testApertureFluxOfGaussian(self):
        shape, s = (60, 60), 2.0
        c = lsstStubs.gaussian_comp(shape, 30.3, 29.6, s, 25)
        img = lsstStubs.place(shape, c, 1000.0)
        var = np.ones(shape)
        good = np.ones(shape, bool)
        fl, err, cov = ph.aperture_flux(img, var, good, 30.3, 29.6, 15.0)
        self.assertAlmostEqual(fl, 1000.0, delta=0.5)
        self.assertAlmostEqual(cov, 1.0)
        self.assertAlmostEqual(err, np.sqrt(np.pi * 15.0 ** 2), delta=0.03 * err)
        # 1-sigma radius of a 2-d Gaussian holds 1 - exp(-1/2) = 39.3%
        fl, _, _ = ph.aperture_flux(img, var, good, 30.3, 29.6, s)
        self.assertAlmostEqual(fl / 1000.0, 1 - np.exp(-0.5), delta=0.01)

    def testCoverage(self):
        shape = (40, 40)
        good = np.ones(shape, bool)
        good[:, :20] = False
        fl, err, cov = ph.aperture_flux(np.ones(shape), np.ones(shape), good, 19.5, 20.0, 5.0)
        self.assertAlmostEqual(cov, 0.5, delta=0.05)

    def testOffImage(self):
        fl, err, cov = ph.aperture_flux(np.ones((10, 10)), np.ones((10, 10)),
                                        np.ones((10, 10), bool), 100.0, 100.0, 3.0)
        self.assertTrue(np.isnan(fl))
        self.assertEqual(cov, 0.0)

    def testApertureCorrection(self):
        k = lsstStubs.gaussian_kernel(2.0, 20)
        self.assertAlmostEqual(ph.aperture_correction(k, 2.0), 1 - np.exp(-0.5), delta=0.01)
        self.assertAlmostEqual(ph.aperture_correction(k, 30.0), 1.0, places=6)


class FisherErrorTestCase(unittest.TestCase):
    def testIsolatedSource(self):
        """sigma(amp) = 1 / sqrt(sum P^2 / var) for one source."""
        shape = (50, 50)
        c = lsstStubs.gaussian_comp(shape, 25.0, 25.0, 2.0, 15)
        var = np.full(shape, 4.0)
        err = ph.fisher_flux_errors([c], var, np.ones(shape, bool))
        self.assertAlmostEqual(err[0], 1 / np.sqrt(np.sum(c["stamp"] ** 2 / 4.0)), places=10)

    def testBlendHasLargerErrors(self):
        shape = (60, 60)
        var = np.ones(shape)
        good = np.ones(shape, bool)
        iso = ph.fisher_flux_errors([lsstStubs.gaussian_comp(shape, 30, 30, 2.0, 15)], var, good)
        pair = ph.fisher_flux_errors([lsstStubs.gaussian_comp(shape, 30, 30, 2.0, 15),
                                      lsstStubs.gaussian_comp(shape, 32, 30, 2.0, 15)], var, good)
        far = ph.fisher_flux_errors([lsstStubs.gaussian_comp(shape, 15, 15, 2.0, 12),
                                     lsstStubs.gaussian_comp(shape, 45, 45, 2.0, 12)], var, good)
        self.assertGreater(pair[0], 1.5 * iso[0])
        self.assertAlmostEqual(far[0], iso[0], places=8)

    def testMissingComponent(self):
        shape = (30, 30)
        err = ph.fisher_flux_errors([None, lsstStubs.gaussian_comp(shape, 15, 15, 2.0, 10)],
                                    np.ones(shape), np.ones(shape, bool))
        self.assertTrue(np.isinf(err[0]))
        self.assertTrue(np.isfinite(err[1]))

    def testErrorsMatchScatter(self):
        """The Fisher error is the actual scatter of least-squares amplitudes."""
        shape = (40, 40)
        comps = [lsstStubs.gaussian_comp(shape, 20, 20, 2.0, 12),
                 lsstStubs.gaussian_comp(shape, 22.5, 21, 2.0, 12)]
        A = np.array([lsstStubs.place(shape, c).ravel() for c in comps]).T
        truth = A @ np.array([500.0, 300.0])
        rng = np.random.default_rng(3)
        amps = np.array([np.linalg.lstsq(A, truth + rng.normal(0, 1, truth.size), rcond=None)[0]
                         for _ in range(400)])
        err = ph.fisher_flux_errors(comps, np.ones(shape), np.ones(shape, bool))
        np.testing.assert_allclose(amps.std(axis=0), err, rtol=0.12)


class MeasurePhotometryTestCase(unittest.TestCase):
    def testSyntheticRun(self):
        import fakeRun
        run = fakeRun.make_run(n_side=5, seed=4)
        t = ph.measure_photometry(run, verbose=False)
        self.assertEqual(len(t), 25)
        self.assertEqual(t.meta["zero_point_source"], "image units assumed to be nJy")
        # PSF fluxes agree with the truth within their errors
        pull = (t["psf_flux"] - run.truth_flux) / t["psf_flux_err"]
        self.assertLess(abs(np.mean(pull)), 1.0)
        self.assertLess(np.std(pull), 1.6)
        # aperture photometry (corrected to total) agrees with PSF photometry for bright stars
        bright = t["psf_snr"] > 50
        self.assertLess(abs(np.median(t["ap_minus_psf"][bright])), 0.01)
        np.testing.assert_allclose(t["psf_mag"], 31.4 - 2.5 * np.log10(t["psf_flux"]))

    def testPhotoCalib(self):
        import fakeRun
        run = fakeRun.make_run(n_side=2, photo_scale=2.0)
        t = ph.measure_photometry(run, verbose=False)
        np.testing.assert_allclose(t["psf_mag"], 31.4 - 2.5 * np.log10(2.0 * t["psf_flux"]))


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
        gaia.meta["gaia_mag_system"] = "Vega (Gaia archive)"
        G = np.array([15., 16., 17., 18., 21.])
        phot = Table(dict(x=[10.1, 20.1, 30.1, 40.4, 33.0], y=[10., 20., 30., 40.2, 12.0],
                          gaia_index=[0, 1, 2, -1, -1], psf_mag=G + 0.3, psf_snr=np.full(5, 100.),
                          extended=np.zeros(5, bool),
                          ra=50.0 + np.array([2, 2, 2, 0, 0]) / 3.6e6 / np.cos(np.radians(dec0)),
                          dec=dec0 + np.array([1, 1, 1, 0, 0]) / 3.6e6))
        phot["psf_flux"] = 10 ** (-0.4 * (phot["psf_mag"] - ph.AB_NJY_ZERO_POINT))
        phot.meta["zero_point"] = ph.AB_NJY_ZERO_POINT
        return phot, gaia

    def testMatchingAndDifferences(self):
        phot, gaia = self._tables()
        t = ph.compare_with_gaia(phot, gaia, verbose=False)
        np.testing.assert_array_equal(t["gaia_match"], [True, True, True, True, False])
        self.assertEqual(list(t["gaia_match_type"]), ["seed", "seed", "seed", "position", ""])
        np.testing.assert_array_equal(t["gaia_source_id"], [100, 101, 102, 103, -1])
        np.testing.assert_allclose(t["gaia_dx"][:3], 0.1, atol=1e-9)
        np.testing.assert_allclose(t["gaia_dra"][:3], 2.0, atol=1e-3)       # mas, with cos(dec)
        np.testing.assert_allclose(t["gaia_ddec"][:3], 1.0, atol=1e-3)
        np.testing.assert_allclose(t["gaia_sep"][:3], np.sqrt(5), atol=1e-3)
        self.assertAlmostEqual(t.meta["GAIAOFF"], 0.3, places=9)            # psf_mag - G
        np.testing.assert_allclose(t["gaia_dmag_rel"][:4], 0.0, atol=1e-9)
        np.testing.assert_allclose(t["gaia_flux_ratio"][:4], 1.0, rtol=1e-9)
        np.testing.assert_allclose(t["gaia_bp_rp"][:4], 0.7, atol=1e-9)
        self.assertEqual(t.meta["GAIAUNM"], 1)            # star at (25, 15) has no fitted source
        self.assertTrue(np.isnan(t["gaia_G"][4]))
        self.assertNotIn("gaia_match", phot.colnames)    # input not modified

    def testRadiusAndNoSky(self):
        phot, gaia = self._tables()
        del phot["ra"], phot["dec"]
        t = ph.compare_with_gaia(phot, gaia, match_radius_px=0.3, verbose=False)
        self.assertFalse(t["gaia_match"][3])               # 0.45 px away
        self.assertNotIn("gaia_dra", t.colnames)

    def testFromSyntheticRun(self):
        import fakeRun
        from astropy.table import Table
        run = fakeRun.make_run(seed=7)
        xy = run.outer_history[-1]["result"]["S"]["xy"]
        G = ph.AB_NJY_ZERO_POINT - 2.5 * np.log10(run.truth_flux) - 0.2       # offset +0.2
        gaia = Table(dict(source_id=np.arange(len(xy)), x=xy[:, 0] + 0.05, y=xy[:, 1],
                          ra=np.zeros(len(xy)), dec=np.zeros(len(xy)), phot_g_mean_mag=G))
        t = ph.measure_photometry(run, verbose=False, gaia=gaia)
        self.assertTrue(np.all(t["gaia_match"]))
        self.assertAlmostEqual(t.meta["GAIAOFF"], 0.2, delta=0.01)
        self.assertLess(np.std(t["gaia_flux_ratio"]), 0.02)

    def testRunGaiaUsedByDefault(self):
        import fakeRun
        from astropy.table import Table
        run = fakeRun.make_run(seed=7)
        xy = run.outer_history[-1]["result"]["S"]["xy"]
        run.gaia = Table(dict(source_id=np.arange(len(xy)), x=xy[:, 0], y=xy[:, 1],
                              ra=np.zeros(len(xy)), dec=np.zeros(len(xy)),
                              phot_g_mean_mag=np.full(len(xy), 18.0)))
        self.assertIn("gaia_dmag", ph.measure_photometry(run, verbose=False).colnames)
        self.assertNotIn("gaia_dmag", ph.measure_photometry(run, verbose=False, gaia=False).colnames)

    def testWritePhotometry(self):
        import os
        import tempfile
        import fakeRun
        from astropy.io import fits
        from astropy.table import Table
        from lsst.meas.crowd import session
        run = fakeRun.make_run(seed=7)
        xy = run.outer_history[-1]["result"]["S"]["xy"]
        run.gaia = Table(dict(source_id=np.arange(len(xy)), x=xy[:, 0], y=xy[:, 1],
                              ra=np.zeros(len(xy)), dec=np.zeros(len(xy)),
                              phot_g_mean_mag=np.full(len(xy), 18.0)))
        run.photometry = None
        with tempfile.TemporaryDirectory() as d:
            path = session.write_photometry(run, d, verbose=False)
            self.assertTrue(os.path.exists(os.path.join(d, "photometry.csv")))
            back = Table.read(path)
            hdr = fits.getheader(path, 1)
        self.assertEqual(len(back), len(xy))
        self.assertIn("gaia_dmag_rel", back.colnames)
        self.assertEqual(hdr["ZP"], ph.AB_NJY_ZERO_POINT)        # short FITS keywords
        self.assertIn("GAIAOFF", hdr)
        self.assertIsNotNone(run.photometry)                     # kept on the run

    def testFitsMeta(self):
        m = ph.fits_meta(dict(zero_point=31.4, aperture_radius_px=3.0, GAIAOFF=np.nan,
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
        t = ph.compare_with_gaia(phot, gaia, verbose=False)
        with tempfile.TemporaryDirectory() as d:
            png = os.path.join(d, "gaia.png")
            out = plotting.plot_gaia_comparison(t, out=png)
            self.assertGreater(os.path.getsize(png), 1000)
        self.assertEqual(out["n"], 4)
        self.assertAlmostEqual(out["dx"], 2.0 * 3 / 4 + 0.0 / 4, delta=1.0)
        with self.assertRaises(ValueError):
            plotting.plot_gaia_comparison(phot)
        plt.close("all")


if __name__ == "__main__":
    unittest.main()
