"""Fit two 300x300 cutouts of one visit and look at the results.

Run cell by cell in a notebook (each "# %%" block), with the Rubin stack and meas_crowd set up::

    setup -r /path/to/meas_crowd      # after `scons` has built it
"""
# %% setup
import time

import lsst.meas.crowd as crowd

print("meas_crowd driver version:", crowd.pipeline.__version__)
print("C++ kernels available:", crowd.HAVE_CXX)

cfg = crowd.make_config(crowd.V5_TESTING_CONFIG)
sess = crowd.prepare_session(
    visit=2025062100433, detector=144,
    cutouts=["1000,1000,300", "2500,2000,300"],     # x0, y0, size (detector pixels)
    pad=20,
    psf_source="visit",
)
crowd.attach_gaia(sess)                              # Gaia DR3 seeds (cached after the first time)

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
