.. py:currentmodule:: lsst.meas.crowd

.. _lsst.meas.crowd:

###############
lsst.meas.crowd
###############

Gaia-seeded, adaptive, iterative joint PSF photometry in crowded fields.

The fit starts from the Gaia DR3 sources in the visit, fits all sources jointly (non-negative
least squares), subtracts the model, finds new sources in the residual and repeats; the inner and
outer loops stop when the residual is consistent with noise.  The PSF is refined from the cutout's
own stars between outer passes.  Saturated and bright stars get extra correction terms; close
pairs and small galaxies are handled by dedicated tests.

User guide
==========

``doc/pipeline_guide.md``: the pipeline flow, key parameters, how-to recipes for every plot and
how to run and check the unit tests.  ``doc/catalog_columns.md``: the output catalogue's columns
and units.

Modules
=======

``pipeline``      the driver: ``run_session_gaia``, ``run_cutout_gaia``, inner / outer loops
``gaiaSeeds``     Gaia catalogue (Butler refcat or ESA TAP) and seed selection
``photometry``    PSF and aperture magnitudes with errors (``measure_photometry``)
``plotting``      run-level plots and residual diagnostics
``sourcePlots``   per-source cutouts, cross-sections and magnitudes (``plot_source_cutouts``)
``stagePlots``    PSF stars, backgrounds, rounds and colour-linked source pages
``detector``      whole-detector tiling, stitching and mosaic / tile plots
``fastKernels``   C++ (pybind11) / NumPy rendering of PSF stamps for a constant PSF
``fitLib``, ``adaptive``, ``session``  the original notebook modules (vendored)
