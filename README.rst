##########
meas_crowd
##########

Adaptive, iterative joint PSF photometry for crowded Rubin/LSST fields, starting from Gaia positions.

The fit starts from the Gaia DR3 sources that fall on a visit cutout, fits all sources jointly
(non-negative least squares), subtracts the model, detects new sources in the residual and
repeats.  Inner and outer loops stop on their own when the residual is consistent with noise.
Between outer passes the PSF is refined from the cutout's own stars and the background is refit.
Saturated and bright stars get local correction terms, close pairs are split, small galaxies are
modelled as elliptical extended sources, and the residual is checked for fake detections around
bright stars.  A whole detector is fitted as a grid of overlapping tiles that are stitched back
into one mosaic and one catalogue.

**Documentation:**

* `doc/pipeline_guide.md <doc/pipeline_guide.md>`_: the pipeline flow and the key parameters, a
  recipe for every plot, and how to run and check the unit tests;
* `doc/algorithm.md <doc/algorithm.md>`_: how every step works, why it is there, and the exact
  decision rules and thresholds;
* `doc/catalog_columns.md <doc/catalog_columns.md>`_: every column of the output catalogue and
  its units.
* `doc/function_reference.md <doc/function_reference.md>`_: every public function and class, its
  parameters and what it returns (generated from the docstrings by
  ``doc/make_function_reference.py``);
* `doc/unit_tests.md <doc/unit_tests.md>`_: every unit test, what it checks, its input and what it
  must return to pass.

Pipeline at a glance
====================

1. **Setup** (``prepare_session``): load the ``visit_image``, cut out each box grown by ``pad`` px,
   build the bad-pixel mask and the detector-level PSF (``psf_source="visit"`` or ``"refit"``).
2. **Gaia catalogue** (``attach_gaia``): Gaia DR3 at the visit epoch, cached on disk; its stars give
   the fit its starting positions.
3. **Preparation**: grow the saturated mask, extend the PSF wings to 81 px, remove a global
   Gaia-vs-image offset.
4. **Outer pass**: the adaptive inner loop (detect, add, joint fit, refine, clean; the threshold
   steps down until the residual is clean), then bright-star corrections, fake-source removal,
   position polish, pair and galaxy tests.
5. **Converged?** Source count and PSF each changed by less than 1 % -> stop.
6. **Refit** the background (raw cutout, model masked) and the PSF (residual stack of clean
   stars) and go back to 4.
7. **Outputs**: a ``CutoutRun`` with every pass, the catalogue and, on request, photometry.

PSF source
----------

``psf_source="visit"`` (default for detector runs) uses the PSF stored with the visit image; the
image is already background-subtracted.  ``psf_source="refit"`` subtracts a background, detects
stars, selects them by size and fits a PCA PSF; the background (``session.bgImage``) is then
subtracted in pass 1, and the candidate and used stars are kept (``session.starCat``,
``session.psf_used_xy``).  Either way the PSF is refined inside each cutout between passes
(``refit_psf="stack"``).

Layout
======

The package follows the Rubin DM layout (``sconsUtils`` + ``eups``)::

    SConstruct, ups/              build / eups configuration
    include/lsst/meas/crowd/      C++ headers
    src/                          C++ sources (lib/libmeas_crowd)
    python/lsst/meas/crowd/
        __init__.py               public API (submodules load on first use)
        pipeline.py               driver: run_session_gaia, inner / outer loops, PSF refinement
        gaiaCatalog.py            Gaia catalogue (Butler refcat or ESA TAP), starting positions
        photometry.py             PSF and aperture magnitudes with errors
        plotting.py               data | model | residual, residual diagnostics, star zooms
        sourcePlots.py            per-source cutouts + x / y / radial cross-sections + magnitudes
        stagePlots.py             PSF stars + model PSF, backgrounds, rounds, linked source pages
        detector.py               whole-detector tiling, stitching, mosaic and tile plots
        fastKernels.py            fast PSF stamps (C++ via pybind11, NumPy fallback)
        _crowdLib.cc              pybind11 wrapper of src/
        fitLib.py, adaptive.py, session.py
                                  the original cutout_fit_lib / cutout_adaptive / cutout_notebook
                                  modules (package-relative imports, numpydoc docstrings)
    bin.src/crowdDetector.py      command line: fit / stitch / plot a whole detector
    tests/                        unit tests (run with or without the stack)
    examples/                     run_visit.py (notebook cells), run_visit_batch.py (terminal),
                                  usdf_job.sh, usdf_detector_array.sh (Slurm), benchmark_stamps.py
    doc/                          user guide, algorithm, function reference, catalogue, unit tests

Building
========

With the Rubin stack (e.g. at USDF)::

    source /sdf/group/rubin/sw/w_latest/loadLSST.bash
    setup lsst_distrib
    cd meas_crowd
    setup -k -r .      # set up this package (pybind11 comes with the conda environment)
    scons              # builds lib/ and the pybind11 module, runs the tests

Without the stack (laptop / CI; only the C++ extension is built, the pipeline itself needs the
stack)::

    pip install pybind11 numpy scipy matplotlib astropy pytest
    python setup.py build_ext --inplace
    export PYTHONPATH=$PWD/python

If the extension is not built, ``fastKernels`` falls back to an equivalent NumPy implementation
(identical results, about 6x slower).

Running
=======

Cutouts (notebook or Python)
----------------------------

``examples/run_visit.py`` has the full workflow as notebook cells.  In short::

    import lsst.meas.crowd as crowd

    cfg = crowd.make_config(crowd.CROWDED_FIELD_CONFIG)
    sess = crowd.prepare_session(visit=2025062100433, detector=144,
                                 cutouts=["1000,1000,300", "2500,2000,300"], pad=20,
                                 psf_source="visit")
    crowd.attach_gaia(sess)
    runs = crowd.run_session_gaia(sess, cfg, out_dir="results", max_outer=3, max_inner=15)

    run = runs["cutout_x1000_y1000_300x300"]
    crowd.report(run)                                 # what each pass did
    crowd.plot_final_comparison(run)                  # data | model | residual | residual/sigma
    phot = crowd.measure_photometry(run)              # astropy Table, PSF + aperture magnitudes

Options go to ``make_config`` or as keyword overrides to ``run_session_gaia``;
``crowd.PIPELINE_DEFAULTS`` lists them all with comments, and an unknown name raises an error.
All function names of the old ``cutout_gaia.py`` are unchanged, so ``import lsst.meas.crowd as cg``
works with existing notebooks.

Cutouts from a terminal (USDF)
------------------------------

::

    python examples/run_visit_batch.py --visit 2025062100433 --detector 144 \
        --cutouts 1000,1000,300 2500,2000,300 --out results_term     # figures saved as PNG
    sbatch examples/usdf_job.sh                                      # the same as a Slurm job

Whole detector
--------------

The detector is covered by 400 px tiles whose interiors partition it exactly; each tile is
fitted with 30 px padding and written to ``tile_<name>.fits``.  Tasks skip tiles already on disk,
so a job can be restarted.  ``stitch_tiles`` writes ``detector_mosaic.fits`` (DATA, MODEL,
RESIDUAL, VARIANCE, BAD, TILE; opens in DS9) and ``detector_catalog.fits`` (every source once)::

    VISIT=2025062100433 DETECTOR=144 bash examples/usdf_detector_array.sh   # Slurm array + stitch

    # or step by step (--region limits a test run to a few tiles)
    python bin.src/crowdDetector.py fit 2025062100433 144 --out $OUT --region 1000 1800 1000 1800
    python bin.src/crowdDetector.py stitch --out $OUT
    python bin.src/crowdDetector.py plot --out $OUT

Plots
=====

Final results:

============================  ==================================================================
``plot_final_comparison``     data | model | residual | residual/sigma (``grayscale=True``,
                              ``panels="residuals"``, ``circles="none"``)
``plot_residual_gray``        the two residual panels, large, greyscale
``plot_source_cutouts``       per source: data | model | residual | residual/sigma | x-cut |
                              y-cut | radial profile, labelled with PSF and aperture magnitudes,
                              S/N, Gaia G and local chi2; optional PDF
``show_star_zooms``           zooms of the worst (or chosen) stars with their neighbours
``show_cross_sections``       x, y and radial cuts with residual and +-1 sigma band
``diagnose_residual``,        what the leftover residual is made of, star by star
``star_residual_table``,
``core_excess_diagnostic``,
``show_bright_star_stack``
``compare_runs``              side-by-side numbers for several runs of the same cutouts
``plot_gaia_comparison``      magnitude difference with Gaia against G and BP - RP, and the
                              position offsets (scatter and across the image)
============================  ==================================================================

Intermediate stages (``stagePlots``):

============================  ==================================================================
``psf_star_positions``,       candidate and used PSF stars (visit catalogue flags, the
``show_psf_stars``            detector-level refit, or the in-cutout stack) over the image, with
                              the PSF kernel and its radial profile
``show_backgrounds``          background entering each outer pass and the change between passes
``show_rounds``               data | model | residual | residual/sigma per outer pass
                              (``level="outer"``) or per inner iteration (``level="inner"``)
``show_linked_source_pages``  per-source pages (data | model | residual | x / y cuts) framed in
                              the colour of their box on the full data, model and residual;
                              ``out_dir=`` writes one PNG per page (memory stays flat)
``compare_passes``            which sources one outer pass lost or gained relative to another,
                              and why (flux left in the residual, merged, halo / ring check,
                              moved, spurious)
============================  ==================================================================

Whole detector (``detector``):

============================  ==================================================================
``plot_tile``, ``load_tile``  one tile over its fit box, interior outlined
``plot_detector``             the stitched mosaic, binned for display, or a zoom (``view=``)
============================  ==================================================================

Photometry
==========

``measure_photometry(run)`` returns one row per source: position (pixel and sky), Gaia flag and G,
extended flag, PSF flux with its Fisher error (overlapping sources inverted together, so blends
get honest errors), PSF magnitude, aperture flux on the neighbour-subtracted image corrected to
total with the PSF's curve of growth, aperture magnitude, and ap - PSF.  Magnitudes use the
exposure's PhotoCalib when present, else the nJy AB zero point (31.4).

Every cutout fit makes this catalogue at the end (``run.photometry``), and
``run_session_gaia(..., out_dir=...)`` saves it as ``photometry.fits`` / ``photometry.csv`` next
to the images; whole-detector runs put it in each tile file and ``detector_catalog.fits``.

The same table carries the comparison with Gaia (the run keeps the Gaia catalogue it started
from as ``run.gaia``; ``compare_with_gaia(phot, gaia)`` adds the columns to any table): a match flag (the Gaia star the source started from, or the nearest unclaimed one within
1 px), the Gaia source id, G and BP - RP, position offsets in px and mas, ``psf_mag - G`` with and
without the median offset, and the flux ratio to Gaia.  ``doc/catalog_columns.md`` lists every
column with its unit.

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
``run_session_gaia`` fits several cutouts in parallel, and a detector run spreads its tiles over a
Slurm array.

Tests
=====

::

    cd tests
    PYTHONPATH=../python MPLBACKEND=Agg python -m pytest -q .        # all (82 tests)
    PYTHONPATH=../python MPLBACKEND=Agg python -m pytest -v test_photometry.py::ApertureTestCase

With the stack set up, ``scons tests`` runs them too.  The tests run without the Rubin stack
(``tests/lsstStubs.py`` stands in for it) on synthetic images with a known answer;
``tests/fakeRun.py`` builds a small finished run for the plotting and photometry tests.
``test_fastStampsStack.py`` compares the fast stamps with ``lsst.afw`` and only runs when the stack
is set up (the 2 skips otherwise).  ``test_docstrings.py`` checks that every public function has a
numpydoc docstring; CI (``.github/workflows/build.yaml``) also runs flake8 and pydocstyle.

Licence
=======

GNU General Public License v3 or later, like other Rubin packages: see ``LICENSE``; copyright
holders are listed in ``COPYRIGHT``.  Every source file carries the standard LSST licence header.

Still to do
===========

* wrap the driver as an ``lsst.pipe.base`` Task with a ``pex_config`` Config (options are a
  plain dict today, ``PIPELINE_DEFAULTS``);
* Rubin naming conventions (camelCase) for the public functions;
* end-to-end tests on a small real cutout (needs data access in CI).
