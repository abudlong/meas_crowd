.. py:currentmodule:: lsst.meas.crowd

.. _lsst.meas.crowd:

###############
lsst.meas.crowd
###############

Adaptive, iterative joint PSF photometry in crowded fields, starting from Gaia positions.

The fit starts from the Gaia DR3 sources in the visit, fits all sources jointly (non-negative
least squares), subtracts the model, finds new sources in the residual and repeats; the inner and
outer loops stop when the residual is consistent with noise.  The PSF is refined from the cutout's
own stars between outer passes.  Saturated and bright stars get extra correction terms; close
pairs and small galaxies are handled by dedicated tests.

User guide
==========

``doc/algorithm.md``: how every step works and its decision rules.
``doc/pipeline_guide.md``: the pipeline flow, key parameters, how-to recipes for every plot and
how to run and check the unit tests.  ``doc/catalog_columns.md``: the output catalogue's columns
and units.  ``doc/unit_tests.md``: every unit test, its input and what it checks.
``doc/function_reference.md``: every public function, its parameters and return values.

Modules
=======

``pipeline``      the driver: ``runSessionGaia``, ``runCutoutGaia``, inner / outer loops
``gaiaCatalog``   Gaia catalogue (Butler refcat or ESA TAP) and starting positions
``photometry``    PSF and aperture magnitudes with errors (``measurePhotometry``)
``plotting``      run-level plots and residual diagnostics
``sourcePlots``   per-source cutouts, cross-sections and magnitudes (``plotSourceCutouts``)
``stagePlots``    PSF stars, backgrounds, rounds and colour-linked source pages
``detector``      whole-detector tiling, stitching and mosaic / tile plots
``fastKernels``   C++ (pybind11) / NumPy rendering of PSF stamps for a constant PSF
``fitLib``        core library: joint fit, detection, PSF and background tools
``adaptive``      residual check deciding whether sources are still missing
``session``       sessions, configuration, saving runs
