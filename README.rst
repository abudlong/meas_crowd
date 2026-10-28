##########
meas_crowd
##########

Gaia-seeded, adaptive, iterative joint PSF photometry for crowded Rubin/LSST fields.

The fit starts from the Gaia DR3 sources that fall on a visit cutout, fits all sources jointly
(non-negative least squares), subtracts the model, detects new sources in the residual and
repeats.  Inner and outer loops stop on their own when the residual is consistent with noise.
Between outer passes the PSF is refined from the cutout's own stars and the background is refit.
Saturated and bright stars get local correction terms, close pairs are split, small galaxies are
modelled as elliptical extended sources, and the residual is checked for fake detections around
bright stars.

Layout
======

The package follows the Rubin DM layout (``sconsUtils`` + ``eups``)::

    SConstruct, ups/              build / eups configuration
    include/lsst/meas/crowd/      C++ headers
    src/                          C++ sources (lib/libmeas_crowd)
    python/lsst/meas/crowd/
        __init__.py               public API (submodules load on first use)
        pipeline.py               driver: run_session_gaia, inner / outer loops, PSF refinement
        gaiaSeeds.py              Gaia catalogue (Butler refcat or ESA TAP), seed selection
        photometry.py             PSF and aperture magnitudes with errors
        plotting.py               data | model | residual, residual diagnostics, star zooms
        sourcePlots.py            per-source cutouts + x / y / radial cross-sections + magnitudes
        fastKernels.py            fast PSF stamps (C++ via pybind11, NumPy fallback)
        _crowdLib.cc              pybind11 wrapper of src/
        fitLib.py, adaptive.py, session.py
                                  the original cutout_fit_lib / cutout_adaptive / cutout_notebook
                                  modules, vendored unchanged apart from package-relative imports
    tests/                        unit tests (run with or without the stack)
    examples/                     run_visit.py (notebook cells), benchmark_stamps.py
    doc/lsst.meas.crowd/          package documentation

Building
========

With the Rubin stack::

    setup -r .
    scons            # builds lib/, the pybind11 module, runs the tests

Without the stack (laptop / CI; only the C++ extension is built, the pipeline itself needs the
stack)::

    pip install pybind11 numpy scipy matplotlib astropy pytest
    python setup.py build_ext --inplace
    export PYTHONPATH=$PWD/python

If the extension is not built, ``fastKernels`` falls back to an equivalent NumPy implementation
(identical results, about 6x slower).

Running
=======

``examples/run_visit.py`` has the full workflow as notebook cells.  In short::

    import lsst.meas.crowd as crowd

    cfg = crowd.make_config(crowd.V5_TESTING_CONFIG)
    sess = crowd.prepare_session(visit=2025062100433, detector=144,
                                 cutouts=["1000,1000,300", "2500,2000,300"], pad=20,
                                 psf_source="visit")
    crowd.attach_gaia(sess)
    runs = crowd.run_session_gaia(sess, cfg, out_dir="results",
                                  max_outer=3, max_inner=15, psf_polish_cycles=3)

    run = runs["cutout_x1000_y1000_300x300"]
    crowd.plot_final_comparison(run)                  # data | model | residual | residual/sigma
    phot = crowd.measure_photometry(run)              # astropy Table, PSF + aperture magnitudes
    crowd.plot_source_cutouts(run, how="worst", n=8, phot=phot, pdf_path="worst.pdf")

All function names of the old ``cutout_gaia.py`` are unchanged, so ``import lsst.meas.crowd as cg``
works with existing notebooks.

Plots
-----

``plot_final_comparison``   data | model | residual | residual/sigma (``grayscale=True``,
                            ``panels="residuals"``, ``circles="none"``)
``plot_residual_gray``      the two residual panels, large, greyscale
``plot_source_cutouts``     per source: data | model | residual | residual/sigma | x-cut | y-cut |
                            radial profile, labelled with PSF and aperture magnitudes, S/N,
                            Gaia G and local chi2; selections ``worst`` / ``brightest`` /
                            ``faintest`` / ``random`` or explicit positions; optional PDF
``show_star_zooms``, ``show_cross_sections``, ``show_bright_star_stack``,
``star_residual_table``, ``diagnose_residual``, ``core_excess_diagnostic``, ``compare_runs``

Photometry
----------

``measure_photometry(run)`` returns one row per source: position (pixel and sky), Gaia flag and G,
extended flag, PSF flux with its Fisher error (overlapping sources inverted together, so blends
get honest errors), PSF magnitude, aperture flux on the neighbour-subtracted image corrected to
total with the PSF's curve of growth, aperture magnitude, and ap - PSF.  Magnitudes use the
exposure's PhotoCalib when present, else the nJy AB zero point (31.4).

Speed
=====

Most of the run time goes into rendering PSF stamps: ``psf.computeImage`` for every source and
centroid step, five times per source for the position derivatives.  The PSF is constant over a
cutout once the pipeline has extended its wings, so ``fastKernels`` takes the kernel image once and
does the Lanczos-5 sub-pixel shifts itself, in C++:

=====================================  ==========  =========
400 stamps of 81 x 81 px                NumPy       C++
=====================================  ==========  =========
stamps (with wing taper)                300 ms      51 ms
position derivatives (4 shifts each)    548 ms      109 ms
=====================================  ==========  =========

Run ``python examples/benchmark_stamps.py`` with the stack for the comparison with
``psf.computeImage`` on your machine.  The fast path is on by default (``fast_stamps=True``) but is
only used after checking, for each PSF, that it reproduces ``psf.computeImage`` to 2e-3 of the
peak; otherwise the pipeline falls back and says so in the log.  ``n_jobs=`` in
``run_session_gaia`` fits several cutouts in parallel.

Tests
=====

::

    PYTHONPATH=python:tests MPLBACKEND=Agg pytest tests

The tests run without the Rubin stack (``tests/lsstStubs.py`` stands in for it) on synthetic
images with a known answer.  ``tests/test_fastStampsStack.py`` compares the fast stamps with
``lsst.afw`` and only runs when the stack is set up.

Still to do
===========

* choose a licence (Rubin packages use GPLv3) and add ``LICENSE``;
* wrap the driver as an ``lsst.pipe.base`` Task with a ``pex_config`` Config (options are a
  plain dict today, ``SEEDED_DEFAULTS``);
* Rubin naming conventions (camelCase) for the public functions;
* end-to-end tests on a small real cutout (needs data access in CI).
