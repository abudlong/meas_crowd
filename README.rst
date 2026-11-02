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
  must return to pass;
* `doc/command_line.md <doc/command_line.md>`_: a whole-detector run from the command line at
  USDF (test, Slurm array, following the jobs, outputs, troubleshooting).

Pipeline at a glance
====================

1. **Setup** (``prepareSession``): load the ``visit_image``, cut out each box grown by ``pad`` px,
   build the bad-pixel mask and the detector-level PSF (``psfSource="visit"`` or ``"refit"``).
2. **Gaia catalogue** (``attachGaia``): Gaia DR3 at the visit epoch, cached on disk; its stars give
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

``psfSource="visit"`` (default for detector runs) uses the PSF stored with the visit image; the
image is already background-subtracted.  ``psfSource="refit"`` subtracts a background, detects
stars, selects them by size and fits a PCA PSF; the background (``session.bgImage``) is then
subtracted in pass 1, and the candidate and used stars are kept (``session.starCat``,
``session.psfUsedXy``).  Either way the PSF is refined inside each cutout between passes
(``refitPsf="stack"``).

Layout
======

The package follows the Rubin DM layout (``sconsUtils`` + ``eups``)::

    SConstruct, ups/              build / eups configuration
    include/lsst/meas/crowd/      C++ headers
    src/                          C++ sources (lib/libmeas_crowd)
    python/lsst/meas/crowd/
        __init__.py               public API (submodules load on first use)
        pipeline.py               driver: runSessionGaia, inner / outer loops, PSF refinement
        gaiaCatalog.py            Gaia catalogue (Butler refcat or ESA TAP), starting positions
        photometry.py             PSF and aperture magnitudes with errors
        plotting.py               data | model | residual, residual diagnostics, star zooms
        sourcePlots.py            per-source cutouts + x / y / radial cross-sections + magnitudes
        stagePlots.py             PSF stars + model PSF, backgrounds, rounds, linked source pages
        detector.py               whole-detector tiling, stitching, mosaic and tile plots
        runPlots.py               every plot of a run / the detector PSF / chosen tiles, to PNG
        fastKernels.py            fast PSF stamps (C++ via pybind11, NumPy fallback)
        _crowdLib.cc              pybind11 wrapper of src/
        fitLib.py                 core library: joint fit, detection, PSF / background tools
        adaptive.py               residual check deciding whether sources are still missing
        session.py                sessions, configuration, saving runs
    bin.src/crowdDetector.py      command line: fit / stitch / plot a whole detector
    tests/                        unit tests (run with or without the stack)
    examples/                     run_visit.py (notebook cells), run_visit_batch.py (terminal),
                                  usdf_job.sh, usdf_detector_array.sh (Slurm), benchmark_stamps.py
    doc/                          user guide, algorithm, function reference, catalogue, unit tests,
                                  command-line runs

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

    cfg = crowd.makeConfig(crowd.CROWDED_FIELD_CONFIG)
    sess = crowd.prepareSession(visit=2025062100433, detector=144,
                                 cutouts=["1000,1000,300", "2500,2000,300"], pad=20,
                                 psfSource="visit")
    crowd.attachGaia(sess)
    runs = crowd.runSessionGaia(sess, cfg, outDir="results", maxOuter=3, maxInner=15)

    run = runs["cutout_x1000_y1000_300x300"]
    crowd.report(run)                                 # what each pass did
    crowd.plotFinalComparison(run)                  # data | model | residual | residual/sigma
    phot = crowd.measurePhotometry(run)              # astropy Table, PSF + aperture magnitudes

Options go to ``makeConfig`` or as keyword overrides to ``runSessionGaia``;
``crowd.PIPELINE_DEFAULTS`` lists them all with comments, and an unknown name raises an error.
Every public function of ``pipeline``, ``plotting`` and ``gaiaCatalog`` is also available as
``crowd.<name>``.

Names
-----

The code follows the Rubin DM naming conventions: functions, arguments, configuration options
and catalogue columns are camelCase (``runSessionGaia``, ``maxOuter``, ``psfMag``), classes are
CapWords and constants UPPER_CASE.  Older snake_case names (``run_session_gaia``,
``max_outer``, ``stack_kwargs=dict(n_use=...)``, configuration JSON files saved by older
versions) still work and raise a ``FutureWarning`` naming the new spelling.  Output file names
(``photometry.fits``, ``detector_catalog.fits``, ...) are unchanged; catalogues written before
the switch have snake_case column names.

Cutouts from a terminal (USDF)
------------------------------

::

    python examples/run_visit_batch.py --visit 2025062100433 --detector 144 \
        --cutouts 1000,1000,300 2500,2000,300 --out results_term     # figures saved as PNG
    sbatch examples/usdf_job.sh                                      # the same as a Slurm job

Whole detector
--------------

The detector is covered by 400 px tiles whose interiors partition it exactly; each tile is
fitted with 30 px padding and written to ``tile_<name>.fits`` (DATA, MODEL, RESIDUAL,
RESIDUAL_SIGMA, VARIANCE, BAD and the catalogue).  Tasks skip tiles already on disk, so a job can
be restarted.  ``stitchTiles`` writes ``detector_mosaic.fits`` (the same images plus TILE; opens
in DS9) and ``detector_catalog.fits`` (every source once).  ``doc/command_line.md`` walks through
a run at USDF step by step::

    VISIT=2025062100433 DETECTOR=144 bash examples/usdf_detector_array.sh   # Slurm array + stitch

    # or step by step (--region limits a test run to a few tiles)
    python bin.src/crowdDetector.py fit 2025062100433 144 --out $OUT --region 1000 1800 1000 1800
    python bin.src/crowdDetector.py stitch --out $OUT
    python bin.src/crowdDetector.py plot --out $OUT

``--out`` is optional: without it the results go to ``$MEAS_CROWD_OUT`` if that is set, else to
``DEFAULT_OUT`` at the top of ``bin.src/crowdDetector.py``.
``fit`` uses ``CROWDED_STELLAR_FIELD_CONFIG`` (dense stellar fields) unless given
``--config crowded`` (``CROWDED_FIELD_CONFIG``) or ``--config my_config.json``; the Slurm script
takes ``CONFIG=``.

Plots
=====

Final results:

============================  ==================================================================
``plotFinalComparison``       data | model | residual | residual/sigma (``grayscale=True``,
                              ``panels="residuals"``, ``circles="none"``)
``plotResidualGray``          the two residual panels, large, greyscale
``plotSourceCutouts``         per source: data | model | residual | residual/sigma | x-cut |
                              y-cut | radial profile, labelled with PSF and aperture magnitudes,
                              S/N, Gaia G and local chi2; optional PDF
``showStarZooms``             zooms of the worst (or chosen) stars with their neighbours
``showCrossSections``         x, y and radial cuts with residual and +-1 sigma band
``diagnoseResidual``,         what the leftover residual is made of, star by star
``starResidualTable``,
``coreExcessDiagnostic``,
``showBrightStarStack``
``compareRuns``               side-by-side numbers for several runs of the same cutouts
``plotGaiaComparison``        magnitude difference with Gaia against G and BP - RP, and the
                              position offsets (scatter and across the image)
============================  ==================================================================

Intermediate stages (``stagePlots``):

============================  ==================================================================
``psfStarPositions``,         candidate and used PSF stars (visit catalogue flags, the
``showPsfStars``              detector-level refit, or the in-cutout stack) over the image, with
                              the PSF kernel and its radial profile
``showBackgrounds``           background entering each outer pass and the change between passes
``showRounds``                data | model | residual | residual/sigma per outer pass
                              (``level="outer"``) or per inner iteration (``level="inner"``)
``showLinkedSourcePages``     per-source pages (data | model | residual | x / y cuts) framed in
                              the colour of their box on the full data, model and residual;
                              ``outDir=`` writes one PNG per page (memory stays flat)
``comparePasses``             which sources one outer pass lost or gained relative to another,
                              and why (flux left in the residual, merged, halo / ring check,
                              moved, spurious)
============================  ==================================================================

Whole detector (``detector``):

============================  ==================================================================
``plotTile``, ``loadTile``    one tile over its fit box, interior outlined
``plotDetector``              the stitched mosaic, binned for display, or a zoom (``view=``)
============================  ==================================================================

Everything at once (``runPlots``):

============================  ==================================================================
``saveRunPlots``              every plot of one run (PSF stars, stamps and model per pass,
                              backgrounds, rounds, source pages and cutouts, final comparison,
                              bright-star diagnostics, Gaia) as numbered PNGs + ``index.txt``
``saveDetectorPlots``         the detector-level PSF stars, their stamps and the PSF kernel
``fitTiles(plots=...)``       the same for chosen tiles of a detector run
                              (``crowdDetector.py fit --plots sample|all``)
============================  ==================================================================

Photometry
==========

``measurePhotometry(run)`` returns one row per source: position (pixel and sky), Gaia flag and G,
extended flag, PSF flux with its Fisher error (overlapping sources inverted together, so blends
get honest errors), PSF magnitude, aperture flux on the neighbour-subtracted image corrected to
total with the PSF's curve of growth, aperture magnitude, and ap - PSF.  Magnitudes use the
exposure's PhotoCalib when present, else the nJy AB zero point (31.4).

Every cutout fit makes this catalogue at the end (``run.photometry``), and
``runSessionGaia(..., outDir=...)`` saves it as ``photometry.fits`` / ``photometry.csv`` next
to the images; whole-detector runs put it in each tile file and ``detector_catalog.fits``.

The same table carries the comparison with Gaia (the run keeps the Gaia catalogue it started
from as ``run.gaia``; ``compareWithGaia(phot, gaia)`` adds the columns to any table): a match flag (the Gaia star the source started from, or the nearest unclaimed one within
1 px), the Gaia source id, G and BP - RP, position offsets in px and mas, ``psfMag - G`` with and
without the median offset, and the flux ratio to Gaia.  ``doc/catalog_columns.md`` lists every
column with its unit.

Speed
=====

Most of the run time goes into rendering PSF stamps: ``psf.computeImage`` for every source and
centroid step, five times per source for the position derivatives.  The PSF is constant over a
cutout once the pipeline has extended its wings, so ``fastKernels`` takes the kernel image once and
does the Lanczos-5 sub-pixel shifts itself, in C++:

=====================================  ==========  =========
400 stamps of 81 x 81 px               NumPy       C++
=====================================  ==========  =========
stamps (with wing taper)               300 ms      51 ms
position derivatives (4 shifts each)   548 ms      109 ms
=====================================  ==========  =========

Run ``python examples/benchmark_stamps.py`` with the stack for the comparison with
``psf.computeImage`` on your machine.  The fast path is on by default (``fastStamps=True``) but is
only used after checking, for each PSF, that it reproduces ``psf.computeImage`` to 2e-3 of the
peak; otherwise the pipeline falls back and says so in the log.  ``nJobs=`` in
``runSessionGaia`` fits several cutouts in parallel, and a detector run spreads its tiles over a
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
