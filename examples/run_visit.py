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

"""Fit two 300x300 cutouts of one visit and look at the results.

Run cell by cell in a notebook (each "# %%" block), with the Rubin stack and meas_crowd set up::

    setup -r /path/to/meas_crowd      # after `scons` has built it
"""
# %% setup
import time

import lsst.meas.crowd as crowd

print("meas_crowd driver version:", crowd.pipeline.__version__)
print("C++ kernels available:", crowd.HAVE_CXX)

cfg = crowd.make_config(crowd.CROWDED_FIELD_CONFIG)
sess = crowd.prepare_session(
    visit=2025062100433, detector=144,
    cutouts=["1000,1000,300", "2500,2000,300"],     # x0, y0, size (detector pixels)
    pad=20,
    psf_source="visit",
)
crowd.attach_gaia(sess)                              # Gaia DR3 catalogue (cached after the first time)

# %% fit
t0 = time.perf_counter()
runs = crowd.run_session_gaia(sess, cfg, out_dir="results",
                              max_outer=3, max_inner=15, psf_polish_cycles=3, plot_final=False)
print(f"fit: {time.perf_counter() - t0:.0f} s")

# %% data | model | residual
for name, run in runs.items():
    crowd.plot_final_comparison(run)
    crowd.plot_residual_gray(run)

# %% photometry + per-source cutouts with cross-sections and magnitudes
for name, run in runs.items():
    phot = crowd.measure_photometry(run)
    phot.write(f"results/{name}_photometry.fits", overwrite=True)
    crowd.plot_source_cutouts(run, how="worst", n=8, phot=phot,
                              pdf_path=f"results/{name}_worst_sources.pdf")
    crowd.plot_source_cutouts(run, how="brightest", n=4, phot=phot)

# %% diagnostics
for name, run in runs.items():
    crowd.diagnose_residual(run)
    crowd.star_residual_table(run)
