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

"""A small synthetic finished run (same structure as run_session_gaia's CutoutRun) for tests of
the plotting / photometry / diagnostic functions."""

import types

import numpy as np

import lsstStubs


def make_run(n_side=4, sigma=2.0, half=15, shape=(120, 130), x0=1000, y0=2000, sky=100.0,
             seed=0, flux_range=(3e3, 3e5), model_sigma=None, photo_scale=None, sigma_fn=None):
    """Grid of Gaussian stars with Poisson noise; the 'fit' uses Gaussians of model_sigma
    (default: the true sigma, i.e. a perfect model).  sigma_fn(flux, x, y) -> true sigma of each
    star (e.g. brighter-fatter).  Returns a CutoutRun-like namespace."""
    rng = np.random.default_rng(seed)
    ny, nx = shape
    gx = np.linspace(x0 + 20, x0 + nx - 20, n_side)
    gy = np.linspace(y0 + 20, y0 + ny - 20, n_side)
    xy = np.array([(x + rng.uniform(-.5, .5), y + rng.uniform(-.5, .5)) for x in gx for y in gy])
    flux = 10 ** rng.uniform(np.log10(flux_range[0]), np.log10(flux_range[1]), len(xy))
    truth = np.zeros(shape)
    for (x, y), f in zip(xy, flux):
        st = sigma if sigma_fn is None else sigma_fn(f, x - x0, y - y0)
        truth += lsstStubs.place(shape, lsstStubs.gaussian_comp(shape, x, y, st, half, x0, y0), f)
    var = sky + truth
    data = truth + rng.normal(0, 1, shape) * np.sqrt(var)
    ms = sigma if model_sigma is None else model_sigma
    comps = [lsstStubs.gaussian_comp(shape, x, y, ms, half, x0, y0) for x, y in xy]
    A = np.array([lsstStubs.place(shape, c).ravel() for c in comps]).T
    amps = np.linalg.lstsq(A / np.sqrt(var.ravel())[:, None], data.ravel() / np.sqrt(var.ravel()),
                           rcond=None)[0]
    model = (A @ amps).reshape(shape)
    psf = lsstStubs.FakePsf(lsstStubs.gaussian_kernel(ms, half))
    exp = lsstStubs.FakeExposure(data, var, psf=psf, x0=x0, y0=y0, photo_scale=photo_scale)
    from lsst.meas.crowd import pipeline
    S = pipeline._state(xy, from_gaia=True)
    S["g"][:] = 20.0
    good = np.ones(shape, bool)
    fit = dict(amps=amps, model_full=model, good=good, variance=var, comps=comps)
    res = dict(S=S, snaps=[dict(fit=fit)], model_full=model, amps=amps.copy(),
               peaks=[tuple(p) for p in xy], bright_corr={}, is_gaia=S["from_gaia"].copy())
    final = dict(exposure=exp, result=res, psf=psf, poisson_k=0.0, model_err=0.0)
    return types.SimpleNamespace(name="synthetic", outer_history=[final], bad_bitmask=4,
                                 interior_bbox=exp.getBBox(), config=dict(sat_zone_fwhm=0),
                                 truth_flux=flux, raw_crop=_RawCrop(data))


class _RawCrop:
    """run.raw_crop[bbox].image.array"""

    def __init__(self, data):
        self.data = data

    def __getitem__(self, bbox):
        return types.SimpleNamespace(image=types.SimpleNamespace(array=self.data))
