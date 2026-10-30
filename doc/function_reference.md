# meas_crowd: function reference

Every public function and class of `lsst.meas.crowd`: what it does, its parameters (with types and defaults) and what it returns.  Generated from the docstrings by `doc/make_function_reference.py`; do not edit by hand (run `python doc/make_function_reference.py` after changing a docstring).  How the pipeline fits together is in [pipeline_guide.md](pipeline_guide.md); the output catalogue in [catalog_columns.md](catalog_columns.md).

175 functions and 4 classes in 11 modules.  Everything listed here can be used as `crowd.<name>` after `import lsst.meas.crowd as crowd`, or from its module (`from lsst.meas.crowd import photometry`).  Names starting with `_` are internal helpers and are not listed.

## Start here

| Function | What it is for |
| --- | --- |
| [`session.make_config`](#sessionmake_config) | Build a pipeline configuration (keyword arguments of `fitLib.run_cutout_pipeline`). |
| [`session.prepare_session`](#sessionprepare_session) | Load a visit image, cut out the regions and set up the detector-level PSF. |
| [`gaiaCatalog.attach_gaia`](#gaiacatalogattach_gaia) | Load Gaia once for a session and store it as ``session.gaia``. |
| [`pipeline.run_session_gaia`](#pipelinerun_session_gaia) | Run the seeded pipeline on every (or the named) cutout of a session. |
| [`pipeline.report`](#pipelinereport) | Print what the adaptive loops did, pass by pass. |
| [`photometry.measure_photometry`](#photometrymeasure_photometry) | Measure PSF and aperture magnitudes of a run's sources. |
| [`plotting.plot_final_comparison`](#plottingplot_final_comparison) | Plot original \| final model \| residual \| residual/sigma of a run. |
| [`plotting.plot_gaia_comparison`](#plottingplot_gaia_comparison) | Plot how the fitted photometry and positions compare with Gaia. |
| [`stagePlots.show_rounds`](#stageplotsshow_rounds) | Plot data \| model \| residual \| residual/sigma for each round of the fit. |
| [`stagePlots.show_linked_source_pages`](#stageplotsshow_linked_source_pages) | Plot per-source pages whose coloured frames match boxes on the data, model and residual. |
| [`stagePlots.compare_passes`](#stageplotscompare_passes) | Find the sources that one outer pass lost or gained relative to another, and why. |
| [`detector.fit_tiles`](#detectorfit_tiles) | Fit this task's share of a detector's tiles and write one file per tile. |
| [`detector.stitch_tiles`](#detectorstitch_tiles) | Assemble the tile files into full-detector images and one catalogue. |
| [`detector.plot_detector`](#detectorplot_detector) | Plot data \| model \| residual \| residual/sigma of a stitched detector or a zoom. |

## Modules

| Module | Role | Public functions / classes |
| --- | --- | --- |
| [`pipeline`](#pipeline) | the driver: run the Gaia-seeded pipeline on cutouts, inner / outer loops, PSF refinement, diagnostics | 39 |
| [`session`](#session) | load a visit, cut out regions, set up the PSF; configuration; saving runs | 19 |
| [`gaiaCatalog`](#gaiacatalog) | the Gaia catalogue and the seeds of each cutout | 4 |
| [`photometry`](#photometry) | PSF and aperture photometry, errors, the Gaia comparison | 9 |
| [`plotting`](#plotting) | final data / model / residual plots and residual diagnostics | 9 |
| [`stagePlots`](#stageplots) | plots of the intermediate stages and pass-to-pass comparison | 6 |
| [`sourcePlots`](#sourceplots) | per-source cutouts, cross-sections and magnitudes | 3 |
| [`detector`](#detector) | whole-detector runs: tiles, stitching, mosaic and tile plots | 9 |
| [`fastKernels`](#fastkernels) | fast PSF stamps (C++ with a NumPy fallback) | 6 |
| [`fitLib`](#fitlib) | the original cutout_fit_lib: joint fit, detection, PSF and background tools, plots | 68 |
| [`adaptive`](#adaptive) | the original cutout_adaptive: the adaptive inner-loop decisions | 7 |

## pipeline

`lsst.meas.crowd.pipeline`: the driver: run the Gaia-seeded pipeline on cutouts, inner / outer loops, PSF refinement, diagnostics.

Gaia-seeded, adaptive, iterative joint PSF-fitting pipeline.

| Name | Summary |
| --- | --- |
| [`clear_stamp_cache`](#pipelineclear_stamp_cache) | Empty the PSF stamp cache and the stamp-size rule. |
| [`set_stamp_sizes`](#pipelineset_stamp_sizes) | Choose a PSF stamp half-width per source. |
| [`amp_snr_full_fast`](#pipelineamp_snr_full_fast) | Return honest amplitude S/N values, computed block by block. |
| [`profile`](#pipelineprofile) | Run a function under cProfile and print where the time goes. |
| [`enable_speedups`](#pipelineenable_speedups) | Patch fitLib with the cached / fast PSF stamps, fast S/N and fast NNLS (idempotent). |
| [`disable_speedups`](#pipelinedisable_speedups) | Undo `enable_speedups` and empty the caches. |
| [`estimate_gaia_offset`](#pipelineestimate_gaia_offset) | Measure the global offset between the image and the Gaia positions. |
| [`saturation_level`](#pipelinesaturation_level) | Estimate the saturation level of an image. |
| [`measure_model_error`](#pipelinemeasure_model_error) | Measure the PSF model error as a fraction of the model per pixel. |
| [`sat_zones`](#pipelinesat_zones) | Find the zones around saturated stars. |
| [`sat_core_labels`](#pipelinesat_core_labels) | Label the saturated cores, each grown by max(3 px, FWHM/2). |
| [`bright_star_corrections`](#pipelinebright_star_corrections) | Fit local correction terms for every bright or saturated star. |
| [`print_bright_corrections`](#pipelineprint_bright_corrections) | Print the table of bright-star corrections. |
| [`ring_mimics`](#pipelinering_mimics) | Find faint non-Gaia sources that are fakes absorbing a bright star's wing excess. |
| [`split_test`](#pipelinesplit_test) | Split unresolved pairs. |
| [`ellip_test`](#pipelineellip_test) | Model badly fitted sources as small galaxies. |
| [`try_elliptical`](#pipelinetry_elliptical) | Try the elliptical-source test on one source of a finished run, printing every step. |
| [`final_position_polish`](#pipelinefinal_position_polish) | Re-optimise all source positions jointly. |
| [`run_inner_seeded`](#pipelinerun_inner_seeded) | Run one outer pass of the seeded pipeline. |
| [`select_psf_stars`](#pipelineselect_psf_stars) | Select clean, bright, isolated, unsaturated point sources for the PSF stack. |
| [`refine_psf_clean`](#pipelinerefine_psf_clean) | Update the PSF from the residual stack of clean stars, cross-validated. |
| [`refine_psf_median`](#pipelinerefine_psf_median) | Correct the PSF core with the median residual/flux of the bright clean stars. |
| [`smooth_psf_kernel`](#pipelinesmooth_psf_kernel) | Return a hole-free, smooth version of a PSF kernel. |
| [`regularize_psf`](#pipelineregularize_psf) | Return the PSF with its kernel passed through `smooth_psf_kernel`. |
| [`extend_psf_wings_keep`](#pipelineextend_psf_wings_keep) | Extend a PSF with a Moffat tail, optionally keeping its angular structure. |
| [`wing_profile_psf`](#pipelinewing_profile_psf) | Return a PSF whose outer wings are measured from the stars. |
| [`spatial_kernel_psf`](#pipelinespatial_kernel_psf) | Build a linearly varying LSST PSF. |
| [`refine_psf_spatial`](#pipelinerefine_psf_spatial) | Add or update a linear spatial variation of the PSF, cross-validated. |
| [`fit_brightness_term`](#pipelinefit_brightness_term) | Measure or update the brightness-dependent PSF term (brighter-fatter). |
| [`refit_background_lowfit`](#pipelinerefit_background_lowfit) | Fit a constant or tilted-plane background to the unmasked pixels and subtract it in place. |
| [`print_timing`](#pipelineprint_timing) | Print the time spent per pipeline step. |
| [`run_seeded_pipeline`](#pipelinerun_seeded_pipeline) | Run the Gaia-seeded adaptive pipeline on a whole cutout. |
| [`run_cutout_gaia`](#pipelinerun_cutout_gaia) | Run the seeded adaptive pipeline on one cropped exposure. |
| [`check_bright_neighbours`](#pipelinecheck_bright_neighbours) | Test whether sources on the edges of bright stars are real. |
| [`diagnose_residual`](#pipelinediagnose_residual) | Decide whether what is left around bright stars is noise or model error. |
| [`modeled_near`](#pipelinemodeled_near) | Return the pixels closer to a modelled (near-saturated) core than to any other masked pixel. |
| [`unmodeled_near_mask`](#pipelineunmodeled_near_mask) | Find the pixels whose residual cannot be judged because a source there is partly masked. |
| [`run_session_gaia`](#pipelinerun_session_gaia) | Run the seeded pipeline on every (or the named) cutout of a session. |
| [`report`](#pipelinereport) | Print what the adaptive loops did, pass by pass. |

### pipeline.clear_stamp_cache

`clear_stamp_cache()`

Empty the PSF stamp cache and the stamp-size rule.

Call after changing the PSF in place; a new PSF object gets its own cache automatically.

### pipeline.set_stamp_sizes

`set_stamp_sizes(xy, amps, psf, point, noise, fwhm, floor_sigma=0.05, tol=1.5)`

Choose a PSF stamp half-width per source.

The smallest of a few tiers beyond which ``amp * PSF`` stays below ``floor_sigma * noise``, so faint sources get small stamps.  Sources not in the list (new detections) get the full stamp.

| Parameter | Type | Description |
| --- | --- | --- |
| `xy` | `numpy.ndarray` | Source positions, shape (N, 2). |
| `amps` | `numpy.ndarray` | Source amplitudes. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `point` | `lsst.geom.Point2D` | Position at which the PSF profile is evaluated. |
| `noise` | `float` | Typical pixel noise. |
| `fwhm` | `float` | PSF FWHM [px]. |
| `floor_sigma` | `float`, optional | Crop where the source's wing drops below this many noise sigma. |
| `tol` | `float`, optional | Positions within this many px of a listed source use its stamp size. |

| Returns | Type | Description |
| --- | --- | --- |
| `half` | `numpy.ndarray` [`int`] | Half-width per source (-1 = full stamp). |

### pipeline.amp_snr_full_fast

`amp_snr_full_fast(fit, noise_scale=1.0)`

Return honest amplitude S/N values, computed block by block.

Drop-in for `fitLib.amp_snr_full`: ``S/N = amp / sqrt(diag((A^T W A)^-1))``, inverted per independent group of overlapping sources with Cholesky instead of one dense pseudo-inverse.

| Parameter | Type | Description |
| --- | --- | --- |
| `fit` | `dict` | Joint fit (from `fitLib.joint_psf_fit_diag`). |
| `noise_scale` | `float`, optional | Empirical noise scale dividing the S/N. |

| Returns | Type | Description |
| --- | --- | --- |
| `snr` | `numpy.ndarray` | S/N per source. |

### pipeline.profile

`profile(fn, *args, top=30, **kw)`

Run a function under cProfile and print where the time goes.

For example ``runs = crowd.profile(crowd.run_session_gaia, sess, cfg, max_outer=2)``.

| Parameter | Type | Description |
| --- | --- | --- |
| `fn` | callable | Function to profile. |
| `*args` | – | Positional arguments of ``fn``. |
| `top` | `int`, optional | Number of entries printed. |
| `**kw` | – | Keyword arguments of ``fn``. |

| Returns | Type | Description |
| --- | --- | --- |
| `result` | `object` | What ``fn`` returned. |

### pipeline.enable_speedups

`enable_speedups()`

Patch fitLib with the cached / fast PSF stamps, fast S/N and fast NNLS (idempotent).

Called by the pipeline itself; also speeds up fitLib's own `fitLib.run_cutout_pipeline`.

### pipeline.disable_speedups

`disable_speedups()`

Undo `enable_speedups` and empty the caches.

### pipeline.estimate_gaia_offset

`estimate_gaia_offset(exposure, S, fwhm, bad_bitmask, noise, min_peak_snr=20.0, min_stars=5, verbose=True)`

Measure the global offset between the image and the Gaia positions.

The median (data centroid - Gaia position) of bright, unsaturated, fairly isolated seeds, from background-subtracted windowed centroids.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The cutout. |
| `S` | `dict` | Source state with the Gaia seeds (see ``_state``). |
| `fwhm` | `float` | PSF FWHM [px]. |
| `bad_bitmask` | `int` | Mask bits of pixels to ignore. |
| `noise` | `float` | Typical pixel noise. |
| `min_peak_snr` | `float`, optional | Minimum peak S/N of a star used. |
| `min_stars` | `int`, optional | Minimum number of stars; otherwise no offset is applied. |
| `verbose` | `bool`, optional | Print the result. |

| Returns | Type | Description |
| --- | --- | --- |
| `dx, dy` | `float` | Offset (data - Gaia) [px]; 0 if too few stars. |
| `n_used` | `int` | Number of stars used. |

### pipeline.saturation_level

`saturation_level(exposure, P)`

Estimate the saturation level of an image.

``near_sat_level`` if given; else the median of the SAT-flagged pixels, unless those look interpolated (below the brightest unflagged pixel), in which case the brightest unflagged pixel, a lower bound.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The image. |
| `P` | `dict` | Pipeline options (uses ``near_sat_level``). |

| Returns | Type | Description |
| --- | --- | --- |
| `level` | `float` or `None` | Saturation level in image units (None if it cannot be estimated). |
| `how` | `str` | How it was obtained. |

### pipeline.measure_model_error

`measure_model_error(exposure, res, P, verbose=True)`

Measure the PSF model error as a fraction of the model per pixel.

Around unsaturated point sources with peak S/N >= ``model_err_min_snr``, the robust (MAD) variance of residual/sigma in bins of model value is fitted as ``1 + eps^2 * model^2 / var``. Values not significant at 2 sigma give eps = 0.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `res` | `dict` | Result of `run_inner_seeded`. |
| `P` | `dict` | Pipeline options. |
| `verbose` | `bool`, optional | Print the result. |

| Returns | Type | Description |
| --- | --- | --- |
| `eps` | `float` | Model error (0 if not significant). |
| `table` | `list` [`tuple`] | Per model bin: (model, number of pixels, robust variance of residual/sigma, model^2/var). |

### pipeline.sat_zones

`sat_zones(exposure, fwhm, P)`

Find the zones around saturated stars.

The zone band follows the SAT mask's shape: pixels within ``w = sat_zone_fwhm * FWHM + sat_zone_scale * inscribed radius`` of each blob.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The cutout. |
| `fwhm` | `float` | PSF FWHM [px]. |
| `P` | `dict` | Pipeline options (``sat_zone_*``). |

| Returns | Type | Description |
| --- | --- | --- |
| `centres` | `numpy.ndarray` | (x, y, halo radius) per saturated blob, shape (N, 3). |
| `mask` | `numpy.ndarray` [`bool`] | The zone band. |

### pipeline.sat_core_labels

`sat_core_labels(exposure, fwhm, min_px=4)`

Label the saturated cores, each grown by max(3 px, FWHM/2).

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The cutout. |
| `fwhm` | `float` | PSF FWHM [px]. |
| `min_px` | `int`, optional | Ignore saturated blobs smaller than this (hot pixels, cosmic rays). |

| Returns | Type | Description |
| --- | --- | --- |
| `labels` | `numpy.ndarray` [`int`] | Core label per pixel (0 = none). |

### pipeline.bright_star_corrections

`bright_star_corrections(exposure, S, fit, model, psf_peak, noise, fwhm, bad_bitmask, sigma_min=100.0, radius_fwhm=4.0, sweeps=2, radial=True, zones=None, halo_modes=0, bright_modes=0, core_rings=(), core_sigma=0.6, core_modes=0)`

Fit local correction terms for every bright or saturated star.

On top of the joint model, within ``radius_fwhm`` FWHM, each bright star's residual is fitted with ``amp * (a*P + bx*dP/dx + by*dP/dy + cxx*d2P/dx2 + cyy*d2P/dy2 + cxy*d2P/dxdy)`` (flux, position, size and ellipticity), optional radial rings with angular modes (halo and wing structure) and optional narrow core terms (core sharper than the PSF, possibly off-centre). Weighted by the pixel variance (not the sys_frac-inflated one) so the core, where the mismatch is, counts.  Two sweeps, brightest first, for overlapping bright stars.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `S` | `dict` | Source state. |
| `fit` | `dict` | Joint fit. |
| `model` | `numpy.ndarray` | Joint model. |
| `psf_peak` | `float` | Peak value of the unit-flux PSF. |
| `noise` | `float` | Typical pixel noise. |
| `fwhm` | `float` | PSF FWHM [px]. |
| `bad_bitmask` | `int` | Mask bits of pixels to ignore. |
| `sigma_min` | `float`, optional | Stars with peak S/N above this (or saturated) are corrected. |
| `radius_fwhm` | `float`, optional | Fit radius [FWHM]. |
| `sweeps` | `int`, optional | Number of sweeps. |
| `radial` | `bool`, optional | Add radial ring terms. |
| `zones` | `numpy.ndarray`, optional | Saturated-star zones (see `sat_zones`); stars in them get halo terms out to the zone. |
| `halo_modes` | `int`, optional | Angular modes of the saturated stars' rings. |
| `bright_modes` | `int`, optional | Angular modes of the unsaturated stars' rings. |
| `core_rings` | `tuple` [`float`], optional | Radii [px] of the narrow core terms (unsaturated stars). |
| `core_sigma` | `float`, optional | Width [px] of the narrow core terms. |
| `core_modes` | `int`, optional | Angular modes of the narrow core terms. |

| Returns | Type | Description |
| --- | --- | --- |
| `corr` | `numpy.ndarray` | Correction image (add to the joint model). |
| `info` | `dict` [`int`, `dict`] | Per corrected star: ``x``, ``y``, ``snr``, ``n_sat``, ``flux``, ``dflux_pct``, ``dx``, ``dy``, ``chi2_before``, ``chi2_after``. |

### pipeline.print_bright_corrections

`print_bright_corrections(info, max_rows=15)`

Print the table of bright-star corrections.

| Parameter | Type | Description |
| --- | --- | --- |
| `info` | `dict` [`int`, `dict`] | Second output of `bright_star_corrections`. |
| `max_rows` | `int`, optional | Number of (brightest) stars printed. |

### pipeline.ring_mimics

`ring_mimics(exposure, S, fit, model_full, P, fwhm, psf_peak, noise)`

Find faint non-Gaia sources that are fakes absorbing a bright star's wing excess.

A faint source near a bright star must stand out from the residual's wing level at its distance from the star (including angular modes) and from the same spot rotated by 90, 180 and 270 degrees (four-fold diffraction structure).

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `S` | `dict` | Source state. |
| `fit` | `dict` | Joint fit. |
| `model_full` | `numpy.ndarray` | Model including the bright-star corrections. |
| `P` | `dict` | Pipeline options (``ring_*``). |
| `fwhm` | `float` | PSF FWHM [px]. |
| `psf_peak` | `float` | Peak value of the unit-flux PSF. |
| `noise` | `float` | Typical pixel noise. |

| Returns | Type | Description |
| --- | --- | --- |
| `drop` | `list` [`int`] | Indices of the sources to remove. |
| `n_checked` | `int` | Number of sources checked. |

### pipeline.split_test

`split_test(exposure, psf, S, fit, model_full, P, fwhm, psf_peak, noise, bad_bitmask, fitkw)`

Split unresolved pairs.

For each source with ``split_min_snr <= peak S/N < bright_corr_sigma``, a positive residual peak (matched-filter S/N >= ``split_peak_snr``) 0.5-1.5 FWHM away and away from other sources gets a trial source.  All trials are fitted together; a trial is kept if the chi2 within 2 FWHM drops by at least ``split_dchi2`` and its S/N is at least ``min_snr``.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `S` | `dict` | Source state. |
| `fit` | `dict` | Joint fit. |
| `model_full` | `numpy.ndarray` | Model including the bright-star corrections. |
| `P` | `dict` | Pipeline options (``split_*``). |
| `fwhm` | `float` | PSF FWHM [px]. |
| `psf_peak` | `float` | Peak value of the unit-flux PSF. |
| `noise` | `float` | Typical pixel noise. |
| `bad_bitmask` | `int` | Mask bits of pixels to ignore. |
| `fitkw` | `dict` | Keyword arguments of the joint fit. |

| Returns | Type | Description |
| --- | --- | --- |
| `S` | `dict` | Source state with the kept second sources. |
| `fit` | `dict` | Joint fit (the input object if nothing was kept). |
| `n_tried` | `int` | Number of trial sources. |
| `n_kept` | `int` | Number kept. |

### pipeline.ellip_test

`ellip_test(exposure, psf, S, fit, model_full, P, fwhm, psf_peak, noise, bad_bitmask, fitkw, psf_sigma, verbose_ellip=True, only=None, explain=False)`

Model badly fitted sources as small galaxies.

A badly fitted source (local chi2/px >= ``ellip_chi2_min``, S/N between ``ellip_min_snr`` and ``bright_corr_sigma``) and the points within ``ellip_group_fwhm`` FWHM (at most one Gaia source) are refitted with alternative models: 'one' (one extended source replaces all of them) and 'star+gal' (the Gaia source stays a point, one extended source replaces the others), each with a Gaussian or exponential elliptical profile convolved with the PSF.  The best model, penalised by 5 per extra fitted number, must beat the point sources by at least ``ellip_dchi2``.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `S` | `dict` | Source state. |
| `fit` | `dict` | Joint fit. |
| `model_full` | `numpy.ndarray` | Model including the bright-star corrections. |
| `P` | `dict` | Pipeline options (``ellip_*``). |
| `fwhm` | `float` | PSF FWHM [px]. |
| `psf_peak` | `float` | Peak value of the unit-flux PSF. |
| `noise` | `float` | Typical pixel noise. |
| `bad_bitmask` | `int` | Mask bits of pixels to ignore. |
| `fitkw` | `dict` | Keyword arguments of the joint fit. |
| `psf_sigma` | `float` | PSF sigma [px]. |
| `verbose_ellip` | `bool`, optional | Print every model tried. |
| `only` | `list` [`int`], optional | Only try these source indices. |
| `explain` | `bool`, optional | Print why sources are not tried, and tracebacks of failed fits. |

| Returns | Type | Description |
| --- | --- | --- |
| `S` | `dict` | Source state with the extended sources. |
| `fit` | `dict` | Joint fit. |
| `n_tried` | `int` | Number of sources tried. |
| `n_kept` | `int` | Number of extended sources kept. |

### pipeline.try_elliptical

`try_elliptical(run, xy, chi2_min=None, dchi2=None)`

Try the elliptical-source test on one source of a finished run, printing every step.

Shows why the source is or is not tried, the fitted models, the chi2 change and any error. The run is not changed.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `xy` | `tuple` [`float`, `float`] | Position; the nearest source is tested. |
| `chi2_min` | `float`, optional | Override the trigger (``ellip_chi2_min``). |
| `dchi2` | `float`, optional | Override the acceptance bar (``ellip_dchi2``). |

| Returns | Type | Description |
| --- | --- | --- |
| `result` | `dict` or `None` | ``kept``, ``S`` and ``fit`` (with the extended source if it was kept); None if the run's last fit has no per-source stamps. |

### pipeline.final_position_polish

`final_position_polish(exposure, psf, S, fit, bad_bitmask, P, fitkw, fwhm, psf_sigma)`

Re-optimise all source positions jointly.

Gauss-Newton steps (capped at ``final_polish_cap`` px, Gaia drift limit kept) and refit, repeated while the total chi2 drops.  Fixes dipoles from positions that went stale when neighbours were added or removed.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `S` | `dict` | Source state. |
| `fit` | `dict` | Joint fit of ``S``. |
| `bad_bitmask` | `int` | Mask bits of pixels to ignore. |
| `P` | `dict` | Pipeline options (``final_polish_*``). |
| `fitkw` | `dict` | Keyword arguments of the joint fit. |
| `fwhm` | `float` | PSF FWHM [px]. |
| `psf_sigma` | `float` | PSF sigma [px]. |

| Returns | Type | Description |
| --- | --- | --- |
| `S` | `dict` | Source state with the polished positions. |
| `fit` | `dict` | Joint fit (the input object if nothing improved). |
| `info` | `dict` | ``rounds``, ``moved``, ``chi0`` and ``chi1`` (total chi2 before / after). |

### pipeline.run_inner_seeded

`run_inner_seeded(exposure, psf, bad_bitmask, S, P, outer_idx=0, threshold_start=None, shape_now=False, timer=None)`

Run one outer pass of the seeded pipeline.

A warm fit of the current sources (Gaia seeds or the previous pass), then detect -> add -> joint fit -> refine -> clean, repeated until the residual check says stop; then the halo, ring, position-polish, pair-split and galaxy tests.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Background-subtracted cutout. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `bad_bitmask` | `int` | Mask bits of pixels to ignore. |
| `S` | `dict` | Starting source state. |
| `P` | `dict` | Pipeline options. |
| `outer_idx` | `int`, optional | Outer-pass index. |
| `threshold_start` | `float`, optional | First detection threshold [sigma]. |
| `shape_now` | `bool`, optional | Apply the shape cut in this pass. |
| `timer` | `lsst.meas.crowd.fitLib.StepTimer`, optional | Timer collecting the step durations. |

| Returns | Type | Description |
| --- | --- | --- |
| `result` | `dict` | Compatible with `fitLib.run_inner_loop_verbose`'s result (``model_full``, ``residual_full``, ``peaks``, ``amps``, histories, ``snaps``), plus ``S``, ``bright_corr``, ``adaptive`` (the loop decisions), ``is_seed`` and ``gaia_index``. |

### pipeline.select_psf_stars

`select_psf_stars(exposure, res, fwhm, psf_peak, noise, P, verbose=True)`

Select clean, bright, isolated, unsaturated point sources for the PSF stack.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `res` | `dict` | Result of `run_inner_seeded`. |
| `fwhm` | `float` | PSF FWHM [px]. |
| `psf_peak` | `float` | Peak value of the unit-flux PSF. |
| `noise` | `float` | Typical pixel noise. |
| `P` | `dict` | Pipeline options (``psf_stack_*``). |
| `verbose` | `bool`, optional | Print the selection. |

| Returns | Type | Description |
| --- | --- | --- |
| `idx` | `numpy.ndarray` [`int`] | Indices of the selected stars, brightest first. |

### pipeline.refine_psf_clean

`refine_psf_clean(exposure, psf, res, bad_bitmask, P, fitkw, n_cycles=3, verbose=True)`

Update the PSF from the residual stack of clean stars, cross-validated.

The clean stars are split into two halves A and B (alternating in brightness).  Each cycle builds a correction from A and judges it on B's pixels, and vice versa; only if the held-out chi2 improves by at least ``psf_cv_min_gain`` is the correction from all stars applied.  This stops the PSF absorbing the pixel noise of its own stars, which shows up as speckle around every bright star.  Corrections are smoothed (``psf_smooth_px``) and soft-thresholded (``psf_soft_sigma``).

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psf` | `lsst.afw.detection.Psf` | Current PSF. |
| `res` | `dict` | Result of `run_inner_seeded`. |
| `bad_bitmask` | `int` | Mask bits of pixels to ignore. |
| `P` | `dict` | Pipeline options. |
| `fitkw` | `dict` | Keyword arguments of the joint fit. |
| `n_cycles` | `int`, optional | Maximum number of correction cycles. |
| `verbose` | `bool`, optional | Print each cycle. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | The updated PSF (the input PSF if no cycle was kept). |

### pipeline.refine_psf_median

`refine_psf_median(exposure, psf, res, bad_bitmask, P, fitkw, fit0=None, verbose=True)`

Correct the PSF core with the median residual/flux of the bright clean stars.

The median is robust to neighbours, cosmic rays and per-star noise, so it keeps the shape error shared by the bright stars.  Applied within ``psf_core_radius_fwhm`` FWHM (cosine taper), soft-thresholded, smoothed and cross-validated on two halves of the stars.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psf` | `lsst.afw.detection.Psf` | Current PSF. |
| `res` | `dict` | Result of `run_inner_seeded`. |
| `bad_bitmask` | `int` | Mask bits of pixels to ignore. |
| `P` | `dict` | Pipeline options (``psf_core_*``). |
| `fitkw` | `dict` | Keyword arguments of the joint fit. |
| `fit0` | `dict`, optional | Joint fit with ``psf`` (computed if not given). |
| `verbose` | `bool`, optional | Print each cycle. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | The corrected PSF (the input PSF if no correction was kept). |

### pipeline.smooth_psf_kernel

`smooth_psf_kernel(kernel, fwhm, core_fwhm=1.5, sigma=1.0, max_pattern=1.5)`

Return a hole-free, smooth version of a PSF kernel.

Inside ``core_fwhm`` FWHM the kernel is kept.  Outside it is rebuilt as ``radial profile(r) * exp(angular pattern)``: the profile is the median of each 1-px ring (so a few holed pixels cannot pull it down), forced to decrease outwards; the pattern is ``log(kernel / profile)`` median-filtered over 5x5 px (removes isolated holes and hot pixels), smoothed by ``sigma`` px and clipped to +-``max_pattern``.  Smooth structure (ellipticity, lobes, spikes wider than ~2 px) survives; pixel-scale defects do not.

| Parameter | Type | Description |
| --- | --- | --- |
| `kernel` | `numpy.ndarray` | PSF kernel (odd size, centred). |
| `fwhm` | `float` | PSF FWHM [px]. |
| `core_fwhm` | `float`, optional | Radius kept unchanged [FWHM]. |
| `sigma` | `float`, optional | Smoothing of the angular pattern [px]. |
| `max_pattern` | `float`, optional | Clip of the log angular pattern. |

| Returns | Type | Description |
| --- | --- | --- |
| `kernel` | `numpy.ndarray` | Smoothed, non-negative kernel summing to 1. |

### pipeline.regularize_psf

`regularize_psf(psf, point, P, verbose=True, label='')`

Return the PSF with its kernel passed through `smooth_psf_kernel`.

| Parameter | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | PSF. |
| `point` | `lsst.geom.Point2D` | Position at which the kernel is taken. |
| `P` | `dict` | Pipeline options (``psf_regularize``, ``psf_reg_*``); no-op if regularisation is off. |
| `verbose` | `bool`, optional | Print how much flux moved. |
| `label` | `str`, optional | Suffix for the printed line. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | The regularised (constant) PSF, or the input PSF. |

### pipeline.extend_psf_wings_keep

`extend_psf_wings_keep(psf, point, size=81, blend=(12.0, 16.0), keep_angular=False, verbose=True)`

Extend a PSF with a Moffat tail, optionally keeping its angular structure.

Like `fitLib.extend_psf_wings`, but with ``keep_angular`` only the kernel's azimuthal mean is replaced by the Moffat beyond the blend radius, so spikes and other non-circular structure are kept.

| Parameter | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | PSF. |
| `point` | `lsst.geom.Point2D` | Position at which the kernel is taken. |
| `size` | `int`, optional | Size of the extended kernel [px]. |
| `blend` | `tuple` [`float`, `float`], optional | Radial range of the cross-fade to the Moffat [px]. |
| `keep_angular` | `bool`, optional | Keep the kernel's angular pattern in the tail. |
| `verbose` | `bool`, optional | Print the fitted tail. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | The extended constant PSF (the input PSF if the fit fails). |

### pipeline.wing_profile_psf

`wing_profile_psf(exposure, res, psf, P, verbose=True)`

Return a PSF whose outer wings are measured from the stars.

Neighbour-subtracted radial profiles of unsaturated bright stars (peak S/N >= ``wing_min_snr``), divided by their flux, give the normalised profile where they have signal; saturated stars are scaled to it where both are measured and extend it outwards.  A Moffat fitted to the combined profile replaces the kernel beyond ``wing_blend_fwhm``.  Degenerate fits are rejected.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `res` | `dict` | Result of `run_inner_seeded`. |
| `psf` | `lsst.afw.detection.Psf` | Current PSF. |
| `P` | `dict` | Pipeline options (``wing_*``, ``psf_wing_size``). |
| `verbose` | `bool`, optional | Print the fitted profile. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` or `None` | The new PSF, or None if there are too few stars or the fit is degenerate. |

### pipeline.spatial_kernel_psf

`spatial_kernel_psf(K0, Kx, Ky, xc, yc, scale)`

Build a linearly varying LSST PSF.

``K(x, y) = K0 + (x - xc)/scale * Kx + (y - yc)/scale * Ky``, as a LinearCombinationKernel with linear spatial functions (the same machinery as the PCA PSF).

| Parameter | Type | Description |
| --- | --- | --- |
| `K0` | `numpy.ndarray` | Kernel at (xc, yc). |
| `Kx, Ky` | `numpy.ndarray` | Kernel change per ``scale`` px in x and y. |
| `xc, yc` | `float` | Reference position. |
| `scale` | `float` | Length scale [px]. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.meas.algorithms.KernelPsf` | The spatially varying PSF. |

### pipeline.refine_psf_spatial

`refine_psf_spatial(exposure, psf_old, psf_new, res, bad_bitmask, P, fitkw, step=0.8, verbose=True)`

Add or update a linear spatial variation of the PSF, cross-validated.

Measured from the clean stars' residuals in ``res`` (fitted with ``psf_old``); the gradient from half A is judged on half B and vice versa, and kept only if the held-out chi2 drops by at least ``psf_spatial_min_gain`` compared with the constant ``psf_new``.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psf_old` | `lsst.afw.detection.Psf` | PSF the residuals were fitted with. |
| `psf_new` | `lsst.afw.detection.Psf` | Constant PSF to add the variation to. |
| `res` | `dict` | Result of `run_inner_seeded`. |
| `bad_bitmask` | `int` | Mask bits of pixels to ignore. |
| `P` | `dict` | Pipeline options (``psf_spatial_*``). |
| `fitkw` | `dict` | Keyword arguments of the joint fit. |
| `step` | `float`, optional | Fraction of the measured gradient applied. |
| `verbose` | `bool`, optional | Print the result. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | The varying PSF, or ``psf_new`` if the variation was not kept. |

### pipeline.fit_brightness_term

`fit_brightness_term(exposure, psf, res, bad_bitmask, P, fitkw, verbose=True)`

Measure or update the brightness-dependent PSF term (brighter-fatter).

From isolated, unsaturated stars with peak S/N >= ``bf_min_snr``: per 1-px ring, a weighted fit ``residual/flux = A(r) + B(r) * (log10 flux - L0)`` gives the slope B(r), smoothed, tapered and made flux-neutral.  Built from half of the stars and judged on the other half; kept only if the held-out chi2 improves.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psf` | `lsst.afw.detection.Psf` | Current PSF. |
| `res` | `dict` | Result of `run_inner_seeded`. |
| `bad_bitmask` | `int` | Mask bits of pixels to ignore. |
| `P` | `dict` | Pipeline options (``bf_*``). |
| `fitkw` | `dict` | Keyword arguments of the joint fit. |
| `verbose` | `bool`, optional | Print the result. |

| Returns | Type | Description |
| --- | --- | --- |
| `moved` | `float` | Fraction of the flux moved per dex of brightness (0 if nothing changed). |

### pipeline.refit_background_lowfit

`refit_background_lowfit(exposure, mode='plane', stats_mask=('BAD', 'SAT', 'EDGE', 'NO_DATA', 'DETECTED', 'CR', 'INTRP'), clip=3.0, n_iter=5, verbose=True)`

Fit a constant or tilted-plane background to the unmasked pixels and subtract it in place.

The robust alternative to `fitLib.refit_background` for small cutouts of an image that is already background-subtracted: a spline / Chebyshev surface on 32 px bins is barely constrained where most pixels are masked and can swing by hundreds of counts in a corner.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Exposure to background-subtract in place (sources flagged DETECTED). |
| `mode` | `str`, optional | 'plane' (a + b x + c y), 'constant' or 'none'. |
| `stats_mask` | `tuple` [`str`], optional | Mask planes whose pixels are ignored. |
| `clip` | `float`, optional | Sigma clipping threshold. |
| `n_iter` | `int`, optional | Clipping iterations. |
| `verbose` | `bool`, optional | Print the fitted background. |

| Returns | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The background-subtracted exposure (same object). |
| `background` | `lsst.afw.image.ImageF` | The subtracted background image. |

Raises: `ValueError` (Raised if ``mode`` is not 'plane', 'constant' or 'none'.)

### pipeline.print_timing

`print_timing(timer=None)`

Print the time spent per pipeline step.

| Parameter | Type | Description |
| --- | --- | --- |
| `timer` | `lsst.meas.crowd.fitLib.StepTimer`, optional | Timer (default: fitLib's ``TIMER``). |

### pipeline.run_seeded_pipeline

`run_seeded_pipeline(exposure_crop, psf, bad_bitmask, seeds=None, seed_gaia_index=None, starCat=None, bgImage=None, viz=False, seed_gmag=None, **cfg)`

Run the Gaia-seeded adaptive pipeline on a whole cutout.

The replacement for `fitLib.run_cutout_pipeline`: outer passes of `run_inner_seeded`, with background refit and PSF refinement in between, until the source list and the PSF stop changing.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure_crop` | `lsst.afw.image.Exposure` | Cutout before background subtraction. |
| `psf` | `lsst.afw.detection.Psf` | Initial PSF. |
| `bad_bitmask` | `int` | Mask bits of pixels to ignore. |
| `seeds` | `numpy.ndarray`, optional | Seed positions (detector px), shape (N, 2); None or empty = no seeds. |
| `seed_gaia_index` | `numpy.ndarray`, optional | Gaia catalogue row of each seed. |
| `starCat` | `lsst.afw.table.SourceCatalog`, optional | Detector-level PSF stars. |
| `bgImage` | `lsst.afw.image.Image`, optional | Detector-level background, subtracted in the first pass. |
| `viz` | `bool`, optional | Draw the per-pass overview. |
| `seed_gmag` | `numpy.ndarray`, optional | Gaia G magnitude of each seed. |
| `**cfg` | – | Configuration: `fitLib.run_cutout_pipeline` keys and ``SEEDED_DEFAULTS`` keys. |

| Returns | Type | Description |
| --- | --- | --- |
| `outer_history` | `list` [`dict`] | One dict per outer pass (as `fitLib.run_cutout_pipeline`'s, plus ``psf_change``, ``model_err``, ``bf`` and ``result['adaptive']``). |

### pipeline.run_cutout_gaia

`run_cutout_gaia(raw_crop, psf, bad_bitmask, gaia_tbl=None, config=None, interior_bbox=None, starCat=None, bgImage=None, viz=False, name='cutout', gmag_max=None, **overrides)`

Run the seeded adaptive pipeline on one cropped exposure.

| Parameter | Type | Description |
| --- | --- | --- |
| `raw_crop` | `lsst.afw.image.Exposure` | Cutout before background subtraction. |
| `psf` | `lsst.afw.detection.Psf` | Initial PSF. |
| `bad_bitmask` | `int` | Mask bits of pixels to ignore. |
| `gaia_tbl` | `astropy.table.Table`, optional | Gaia catalogue (see `lsst.meas.crowd.gaiaCatalog.gaia_catalog`); None = no seeds. |
| `config` | `dict`, optional | Pipeline configuration. |
| `interior_bbox` | `lsst.geom.Box2I`, optional | The requested box without padding (default: the whole cutout). |
| `starCat` | `lsst.afw.table.SourceCatalog`, optional | Detector-level PSF stars. |
| `bgImage` | `lsst.afw.image.Image`, optional | Detector-level background. |
| `viz` | `bool`, optional | Draw the per-pass overview. |
| `name` | `str`, optional | Cutout name. |
| `gmag_max` | `float`, optional | Only seed Gaia sources brighter than this G magnitude. |
| `**overrides` | – | Configuration overrides. |

| Returns | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | The result; its catalogue also has ``from_gaia``, ``gaia_source_id``, ``gaia_G``, ``gaia_dx`` and ``gaia_dy``. |

### pipeline.check_bright_neighbours

`check_bright_neighbours(run, bright_snr=300.0, within_fwhm=4.0, flux_ratio=5.0, n_controls=7, plot=True, max_show=16, verbose=True, similar_min=3)`

Test whether sources on the edges of bright stars are real.

For every source within ``within_fwhm`` FWHM of a star at least ``flux_ratio`` times brighter (with peak S/N >= ``bright_snr``), its PSF-weighted flux is compared with the residual's flux at control positions at the same distance from the bright star, rotated in 45-degree steps. A mismatch pattern of the bright star (rings, spikes, lobes) puts similar flux at several angles; a real star is significantly brighter than its controls.  Verdicts: 'Gaia' (real by construction), 'real' (significance >= 5 and at most one similar angle), 'probable' (significance >= 3, or two similar angles), 'doubtful' (significance < 3, or similar flux at >= ``similar_min`` angles).

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `bright_snr` | `float`, optional | Minimum peak S/N of a bright star. |
| `within_fwhm` | `float`, optional | Neighbour radius [FWHM]. |
| `flux_ratio` | `float`, optional | Minimum brightness ratio of the bright star to its neighbour. |
| `n_controls` | `int`, optional | Number of control angles. |
| `plot` | `bool`, optional | Show thumbnails of the non-Gaia neighbours. |
| `max_show` | `int`, optional | Maximum number of thumbnails. |
| `verbose` | `bool`, optional | Print the table. |
| `similar_min` | `int`, optional | Number of similar control angles that makes a source doubtful. |

| Returns | Type | Description |
| --- | --- | --- |
| `table` | `astropy.table.Table` | One row per neighbour: position, flux, significance, number of similar angles, verdict. |

### pipeline.diagnose_residual

`diagnose_residual(run, bright_snr=300.0, disk_fwhm=3.0, plot=True)`

Decide whether what is left around bright stars is noise or model error.

For empty sky, bright-star disks (within ``disk_fwhm`` FWHM of stars with peak S/N >= ``bright_snr``) and other source pixels: ``sigma_z`` (robust std of residual/sigma; 1 = matches the variance plane) and ``coherence`` (robust std of the 3x3-averaged residual/sigma, x3, over sigma_z, squared; white noise gives 1).  Also the pixel-to-pixel scatter against model flux, which measures missing photon noise.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `bright_snr` | `float`, optional | Minimum peak S/N of a bright star. |
| `disk_fwhm` | `float`, optional | Radius of the bright-star disks [FWHM]. |
| `plot` | `bool`, optional | Show residual/sigma, its 3x3 average and its matched-filter S/N. |

| Returns | Type | Description |
| --- | --- | --- |
| `result` | `dict` | Per region ``sigma_z``, ``coherence`` and pixel count; ``highfreq`` (pixel-to-pixel scatter vs model flux, with ``k``) and ``noise_share`` (fraction of the bright-disk excess explained by missing photon noise). |

### pipeline.modeled_near

`modeled_near(modeled, bad, fwhm)`

Return the pixels closer to a modelled (near-saturated) core than to any other masked pixel.

Residual blobs there are not grown into the excluded area by `unmodeled_near_mask`.

| Parameter | Type | Description |
| --- | --- | --- |
| `modeled` | `numpy.ndarray` [`bool`] | Masked pixels whose star is in the model. |
| `bad` | `numpy.ndarray` [`bool`] | All masked pixels. |
| `fwhm` | `float` | PSF FWHM [px] (unused; kept for a uniform signature). |

| Returns | Type | Description |
| --- | --- | --- |
| `near` | `numpy.ndarray` [`bool`] | Pixels nearer to a modelled core than to other masked pixels. |

### pipeline.unmodeled_near_mask

`unmodeled_near_mask(resid, var, bad, fwhm, pad_fwhm=1.5, grow_snr=2.5, reach_fwhm=4.0, modeled=None)`

Find the pixels whose residual cannot be judged because a source there is partly masked.

Everything within ``pad_fwhm`` FWHM of a bad pixel, plus any positive residual blob (matched-filter S/N > ``grow_snr``) touching that pad, out to ``reach_fwhm`` FWHM from the mask.  Masked cores whose star is in the model (``modeled``, near-saturated cores fitted from their wings) get no pad: the residual next to them is real information.

| Parameter | Type | Description |
| --- | --- | --- |
| `resid` | `numpy.ndarray` | Residual image. |
| `var` | `numpy.ndarray` | Variance image. |
| `bad` | `numpy.ndarray` [`bool`] | Bad pixels. |
| `fwhm` | `float` | PSF FWHM [px]. |
| `pad_fwhm` | `float`, optional | Pad around bad pixels [FWHM]. |
| `grow_snr` | `float`, optional | Matched-filter S/N of residual blobs grown into the excluded area. |
| `reach_fwhm` | `float`, optional | Maximum distance of the excluded area from the mask [FWHM]. |
| `modeled` | `numpy.ndarray` [`bool`], optional | Masked pixels whose star is in the model. |

| Returns | Type | Description |
| --- | --- | --- |
| `excluded` | `numpy.ndarray` [`bool`] | Pixels to exclude (not including the bad pixels themselves). |

### pipeline.run_session_gaia

`run_session_gaia(session, config=None, names=None, out_dir=None, viz=False, plot_final=True, gmag_max=None, n_jobs=1, plot_kwargs=None, **overrides)`

Run the seeded pipeline on every (or the named) cutout of a session.

Calls `lsst.meas.crowd.gaiaCatalog.attach_gaia` first if the session has no Gaia catalogue.

| Parameter | Type | Description |
| --- | --- | --- |
| `session` | `lsst.meas.crowd.session.Session` | The session. |
| `config` | `dict`, optional | Pipeline configuration. |
| `names` | `list` [`str`], optional | Cutouts to run (default: all). |
| `out_dir` | `str`, optional | Save each run under ``<out_dir>/visit<V>_det<D>/<name>/``; needed with ``n_jobs > 1``, since results that cannot be sent back are only on disk. |
| `viz` | `bool`, optional | Draw the per-pass overview. |
| `plot_final` | `bool`, optional | Plot the final comparison of each cutout. |
| `gmag_max` | `float`, optional | Only seed Gaia sources brighter than this G magnitude. |
| `n_jobs` | `int`, optional | Fit cutouts in this many parallel worker processes (Linux fork; one core each). |
| `plot_kwargs` | `dict`, optional | Keyword arguments of `lsst.meas.crowd.plotting.plot_final_comparison`. |
| `**overrides` | – | Configuration overrides (e.g. ``max_outer=3``). |

| Returns | Type | Description |
| --- | --- | --- |
| `runs` | `dict` [`str`, `lsst.meas.crowd.session.CutoutRun`] | The results by cutout name. |

### pipeline.report

`report(run)`

Print what the adaptive loops did, pass by pass.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |

## session

`lsst.meas.crowd.session`: load a visit, cut out regions, set up the PSF; configuration; saving runs.

Notebook-friendly front end for the pipeline (the original cutout_notebook.py).

| Name | Summary |
| --- | --- |
| [`save_config`](#sessionsave_config) | Write a config dict to a JSON file. |
| [`make_config`](#sessionmake_config) | Build a pipeline configuration (keyword arguments of `fitLib.run_cutout_pipeline`). |
| [`check_config`](#sessioncheck_config) | Check that every key of a configuration is a `fitLib.run_cutout_pipeline` option. |
| [`show_config`](#sessionshow_config) | Print a configuration. |
| [`parse_cutout`](#sessionparse_cutout) | Parse a cutout specification. |
| [`Session`](#sessionsession) (class) | Everything needed to fit cutouts of one visit and detector. |
| [`build_detector_psf`](#sessionbuild_detector_psf) | Set up the detector-level PSF. |
| [`prepare_session`](#sessionprepare_session) | Load a visit image, cut out the regions and set up the detector-level PSF. |
| [`query_gaia`](#sessionquery_gaia) | Query Gaia DR3 at the ESA archive for the sources on a box. |
| [`CutoutRun`](#sessioncutoutrun) (class) | Result of running the pipeline on one cutout. |
| [`build_catalog`](#sessionbuild_catalog) | Build the source catalogue of a finished run. |
| [`run_cutout`](#sessionrun_cutout) | Run the pipeline on one cropped exposure. |
| [`run_session`](#sessionrun_session) | Run the pipeline on every (or the named) cutout of a session. |
| [`write_photometry`](#sessionwrite_photometry) | Write a run's photometry catalogue as photometry.fits and photometry.csv. |
| [`save_run`](#sessionsave_run) | Write a run's catalogues, images, PSF kernels and summary to a directory. |
| [`robust_stats`](#sessionrobust_stats) | Return the median and IQR-based robust sigma of finite (and valid) values. |
| [`matplotlib_listed`](#sessionmatplotlib_listed) | Return a one-colour matplotlib colour map. |
| [`plot_final_comparison_robust`](#sessionplot_final_comparison_robust) | Plot original \| final model \| residual \| residual/sigma with robust statistics. |
| [`plot_final_comparison_from_fits`](#sessionplot_final_comparison_from_fits) | Rebuild the final-comparison figure from files written by `save_run`. |

### session.save_config

`save_config(cfg, path)`

Write a config dict to a JSON file.

Tuples become lists; `make_config` turns them back into tuples when the file is loaded.

| Parameter | Type | Description |
| --- | --- | --- |
| `cfg` | `dict` | Pipeline configuration. |
| `path` | `str` | Output file. |

| Returns | Type | Description |
| --- | --- | --- |
| `path` | `str` | The file written. |

### session.make_config

`make_config(base=None, deep=False, n_outer=None, n_inner=None, save_to=None, **overrides)`

Build a pipeline configuration (keyword arguments of `fitLib.run_cutout_pipeline`).

| Parameter | Type | Description |
| --- | --- | --- |
| `base` | `str` or `dict`, optional | A JSON file or a dict, applied on top of ``DEFAULT_CONFIG``. |
| `deep` | `bool`, optional | Merge nested dicts key by key (default: replace them). |
| `n_outer` | `int`, optional | Shortcut for ``n_outer_iterations``. |
| `n_inner` | `int`, optional | Shortcut for ``n_inner_iterations`` and ``inner_kwargs['n_iterations']``. |
| `save_to` | `str`, optional | Also write the configuration to this JSON file. |
| `**overrides` | – | Any `fitLib.run_cutout_pipeline` keyword, applied last. |

| Returns | Type | Description |
| --- | --- | --- |
| `cfg` | `dict` | The configuration. |

### session.check_config

`check_config(cfg)`

Check that every key of a configuration is a `fitLib.run_cutout_pipeline` option.

| Parameter | Type | Description |
| --- | --- | --- |
| `cfg` | `dict` | Pipeline configuration. |

Raises: `ValueError` (Raised if ``cfg`` has keys the pipeline does not accept.)

### session.show_config

`show_config(cfg, diff_only=False)`

Print a configuration.

| Parameter | Type | Description |
| --- | --- | --- |
| `cfg` | `dict` | Pipeline configuration. |
| `diff_only` | `bool`, optional | Only print keys that differ from ``DEFAULT_CONFIG``. |

### session.parse_cutout

`parse_cutout(spec)`

Parse a cutout specification.

| Parameter | Type | Description |
| --- | --- | --- |
| `spec` | `str`, `tuple` or `dict` | 'x0,y0,size', 'x0,y0,w,h', 'x0:x1,y0:y1', a tuple (x0, y0, size) or (x0, y0, w, h), or dict(x0=, y0=, w=, h=, name=); absolute detector pixels. |

| Returns | Type | Description |
| --- | --- | --- |
| `x0, y0` | `int` | Corner of the box. |
| `w, h` | `int` | Size of the box. |
| `name` | `str` or `None` | Name of the cutout, if given. |

Raises: `ValueError` (Raised if ``spec`` cannot be parsed.)

### session.Session

Class.  `Session(visit, detector, bad_bitmask, psf, starCat, bgImage, jobs=<factory>, wcs=None, psf_used_xy=None)`

Everything needed to fit cutouts of one visit and detector.

| Parameter | Type | Description |
| --- | --- | --- |
| `visit` | `int` | Visit id. |
| `detector` | `int` | Detector id. |
| `bad_bitmask` | `int` | Mask bits of pixels the fit ignores. |
| `psf` | `lsst.afw.detection.Psf` | Detector-level PSF. |
| `starCat` | `lsst.afw.table.SourceCatalog` or `None` | PSF stars of the detector-level PSF determination. |
| `bgImage` | `lsst.afw.image.Image` or `None` | Detector-level background. |
| `jobs` | `dict` [`str`, `dict`] | Per cutout: ``interior`` (`lsst.geom.Box2I`, the requested box) and ``raw`` (the padded cutout before background subtraction). |
| `wcs` | `lsst.afw.geom.SkyWcs`, optional | The visit image's WCS. |
| `psf_used_xy` | `numpy.ndarray`, optional | Positions (detector px, shape (N, 2)) of the stars the detector-level PSF was built from (``psf_source='refit'`` only). |

### session.build_detector_psf

`build_detector_psf(visit_image, psf_source='refit', return_used=False)`

Set up the detector-level PSF.

With ``psf_source='visit'`` the visit image's own PSF is used.  With 'refit' the background is subtracted, sources are detected and measured, stars selected and a PCA PSF determined (the visit image is modified in place).

| Parameter | Type | Description |
| --- | --- | --- |
| `visit_image` | `lsst.afw.image.Exposure` | The detector's visit image. |
| `psf_source` | `str`, optional | 'refit' or 'visit'. |
| `return_used` | `bool`, optional | Also return the positions of the stars the PCA PSF was built from. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | Detector-level PSF. |
| `starCat` | `lsst.afw.table.SourceCatalog` or `None` | Selected PSF stars (None for 'visit'). |
| `bgImage` | `lsst.afw.image.Image` or `None` | Subtracted background (None for 'visit'). |
| `used_xy` | `numpy.ndarray` or `None` | Positions of the stars used, shape (N, 2) (only with ``return_used``; None for 'visit'). |

### session.prepare_session

`prepare_session(visit, detector, cutouts, pad=0, psf_source='refit', instrument='LSSTCam', repo='dp2_prep', collection='LSSTCam/runs/DRP/DP2', visit_image_collections='LSSTCam/runs/*', visit_image=None)`

Load a visit image, cut out the regions and set up the detector-level PSF.

The cutouts are taken before the detector-wide background subtraction, as in the original notebook.

| Parameter | Type | Description |
| --- | --- | --- |
| `visit` | `int` | Visit id. |
| `detector` | `int` | Detector id. |
| `cutouts` | `list` | Cutout specifications (see `parse_cutout`); absolute detector pixels. |
| `pad` | `int`, optional | Fit a region grown by this many px; the interior keeps the requested box. |
| `psf_source` | `str`, optional | 'refit' or 'visit' (see `build_detector_psf`). |
| `instrument` | `str`, optional | Instrument name. |
| `repo` | `str`, optional | Butler repository path or alias. |
| `collection` | `str`, optional | Default Butler collection. |
| `visit_image_collections` | `str`, optional | Collections searched for the visit image. |
| `visit_image` | `lsst.afw.image.Exposure`, optional | An already-loaded visit image (skips the Butler). |

| Returns | Type | Description |
| --- | --- | --- |
| `session` | `Session` | The session. |

Raises: `RuntimeError` (Raised if no visit image is found.)

### session.query_gaia

`query_gaia(wcs, bbox, tap_url='https://gea.esac.esa.int/tap-server/tap')`

Query Gaia DR3 at the ESA archive for the sources on a box.

| Parameter | Type | Description |
| --- | --- | --- |
| `wcs` | `lsst.afw.geom.SkyWcs` | WCS of the image. |
| `bbox` | `lsst.geom.Box2I` | Pixel box. |
| `tap_url` | `str`, optional | TAP service URL. |

| Returns | Type | Description |
| --- | --- | --- |
| `gaia` | `astropy.table.Table` | Gaia sources inside the box, with their pixel positions ``x``, ``y``. |

### session.CutoutRun

Class.  `CutoutRun(name, outer_history, catalog, summary, raw_crop, interior_bbox, bad_bitmask, config, gaia=None, photometry=None)`

Result of running the pipeline on one cutout.

| Parameter | Type | Description |
| --- | --- | --- |
| `name` | `str` | Cutout name. |
| `outer_history` | `list` [`dict`] | One dict per outer pass (``psf``, ``exposure``, ``result``, ...). |
| `catalog` | `astropy.table.Table` | Source catalogue (see `build_catalog`). |
| `summary` | `dict` | Run summary (timings, counts, residual statistics). |
| `raw_crop` | `lsst.afw.image.Exposure` | The input cutout before background subtraction. |
| `interior_bbox` | `lsst.geom.Box2I` | The requested box (without padding). |
| `bad_bitmask` | `int` | Mask bits of pixels the fit ignored. |
| `config` | `dict` | Configuration used. |
| `gaia` | `astropy.table.Table`, optional | The Gaia catalogue the run was seeded from (used for the Gaia comparison columns). |
| `photometry` | `astropy.table.Table`, optional | The photometry catalogue (see `lsst.meas.crowd.photometry.measure_photometry`), made at the end of the run and written by `save_run` as ``photometry.fits``. |

### session.build_catalog

`build_catalog(outer_history, interior_bbox, gaia_tbl=None, match_radius=2.0)`

Build the source catalogue of a finished run.

| Parameter | Type | Description |
| --- | --- | --- |
| `outer_history` | `list` [`dict`] | Outer-loop history. |
| `interior_bbox` | `lsst.geom.Box2I` | The requested box (``in_interior`` column). |
| `gaia_tbl` | `astropy.table.Table`, optional | Gaia sources (with ``x``, ``y``) to match. |
| `match_radius` | `float`, optional | Gaia match radius [px]. |

| Returns | Type | Description |
| --- | --- | --- |
| `cat` | `astropy.table.Table` | Positions (pixel and sky), PSF flux, AB magnitude (assuming nJy), source metrics, ``in_interior`` and, with ``gaia_tbl``, the nearest Gaia source. |

### session.run_cutout

`run_cutout(raw_crop, psf, bad_bitmask, config=None, interior_bbox=None, starCat=None, bgImage=None, viz=False, diagnostics=True, gaia_tbl=None, match_radius=2.0, name='cutout', **overrides)`

Run the pipeline on one cropped exposure.

| Parameter | Type | Description |
| --- | --- | --- |
| `raw_crop` | `lsst.afw.image.Exposure` | The cutout (before background subtraction). |
| `psf` | `lsst.afw.detection.Psf` | Initial PSF. |
| `bad_bitmask` | `int` | Mask bits of pixels to ignore. |
| `config` | `dict`, optional | Pipeline configuration (default ``DEFAULT_CONFIG``). |
| `interior_bbox` | `lsst.geom.Box2I`, optional | The requested box (default: the whole cutout). |
| `starCat` | `lsst.afw.table.SourceCatalog`, optional | PSF stars forced into every PSF determination. |
| `bgImage` | `lsst.afw.image.Image`, optional | Detector background subtracted before the first pass. |
| `viz` | `bool`, optional | Draw the per-iteration diagnostics. |
| `diagnostics` | `bool`, optional | Run the end-of-run diagnostics. |
| `gaia_tbl` | `astropy.table.Table`, optional | Gaia sources to match in the catalogue. |
| `match_radius` | `float`, optional | Gaia match radius [px]. |
| `name` | `str`, optional | Cutout name. |
| `**overrides` | – | Configuration overrides. |

| Returns | Type | Description |
| --- | --- | --- |
| `run` | `CutoutRun` | The result. |

### session.run_session

`run_session(session, config=None, names=None, out_dir=None, viz=False, gaia=False, match_radius=2.0, plot_final=True, **overrides)`

Run the pipeline on every (or the named) cutout of a session.

A failure in one cutout is printed and skipped.

| Parameter | Type | Description |
| --- | --- | --- |
| `session` | `Session` | The session. |
| `config` | `dict`, optional | Pipeline configuration. |
| `names` | `list` [`str`], optional | Cutouts to run (default: all). |
| `out_dir` | `str`, optional | Save each run under ``<out_dir>/visit<V>_det<D>/<name>/`` (see `save_run`). |
| `viz` | `bool`, optional | Draw the per-iteration diagnostics. |
| `gaia` | `bool`, optional | Match the catalogues to Gaia (queried at the ESA archive). |
| `match_radius` | `float`, optional | Gaia match radius [px]. |
| `plot_final` | `bool`, optional | Plot the final comparison of each cutout. |
| `**overrides` | – | Configuration overrides. |

| Returns | Type | Description |
| --- | --- | --- |
| `runs` | `dict` [`str`, `CutoutRun`] | The results by cutout name. |

### session.write_photometry

`write_photometry(run, outdir, verbose=True)`

Write a run's photometry catalogue as photometry.fits and photometry.csv.

The catalogue is computed (with the Gaia comparison when the run has its Gaia table) if the run does not have one yet.  A failure is reported, not raised, so the other products are still written.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `CutoutRun` | The run. |
| `outdir` | `str` | Output directory (created if needed). |
| `verbose` | `bool`, optional | Print where the catalogue went. |

| Returns | Type | Description |
| --- | --- | --- |
| `path` | `str` or `None` | The FITS file written (None if the photometry failed). |

### session.save_run

`save_run(run, outdir)`

Write a run's catalogues, images, PSF kernels and summary to a directory.

Writes photometry.fits / photometry.csv (the photometry catalogue with the Gaia comparison, computed here if the run has none), catalog.csv / catalog.fits, images.fits, psf_kernels.fits and summary.json, the layout `plot_final_comparison_from_fits` reads.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `CutoutRun` | The run. |
| `outdir` | `str` | Output directory (created if needed). |

### session.robust_stats

`robust_stats(values, valid=None)`

Return the median and IQR-based robust sigma of finite (and valid) values.

| Parameter | Type | Description |
| --- | --- | --- |
| `values` | `numpy.ndarray` | Values. |
| `valid` | `numpy.ndarray` [`bool`], optional | Values to use. |

| Returns | Type | Description |
| --- | --- | --- |
| `stats` | `dict` | ``median``, ``sigma_iqr`` ((Q75 - Q25) / 1.349), ``q25``, ``q75`` and ``n``. |

### session.matplotlib_listed

`matplotlib_listed(color)`

Return a one-colour matplotlib colour map.

| Parameter | Type | Description |
| --- | --- | --- |
| `color` | `str` | The colour. |

| Returns | Type | Description |
| --- | --- | --- |
| `cmap` | `matplotlib.colors.ListedColormap` | The colour map. |

### session.plot_final_comparison_robust

`plot_final_comparison_robust(run_or_history, exposure_crop=None, bad_bitmask=None, edge_margin_fwhm=None, view='interior', stats_in_view=True, show_peaks=True, figsize=(26, 7), verbose=True)`

Plot original | final model | residual | residual/sigma with robust statistics.

Bad pixels are black in the residual panels; median and IQR-robust sigma (good pixels only) are shown in the titles.

| Parameter | Type | Description |
| --- | --- | --- |
| `run_or_history` | `CutoutRun` or `list` [`dict`] | A run (preferred) or its outer-loop history. |
| `exposure_crop` | `lsst.afw.image.Exposure`, optional | Raw crop (only needed with a bare history). |
| `bad_bitmask` | `int`, optional | Mask bits of bad pixels (needed with a bare history). |
| `edge_margin_fwhm` | `float`, optional | Draw the edge drop-zone of this width [FWHM]. |
| `view` | `str` or `tuple`, optional | 'interior' (the run's interior box), None (the whole fitted array) or (xmin, xmax, ymin, ymax) in detector pixels. |
| `stats_in_view` | `bool`, optional | Restrict the statistics to the displayed view. |
| `show_peaks` | `bool`, optional | Mark the sources. |
| `figsize` | `tuple` [`float`, `float`], optional | Figure size [inches]. |
| `verbose` | `bool`, optional | Print the statistics. |

| Returns | Type | Description |
| --- | --- | --- |
| `stats` | `dict` | ``residual`` and ``residual_over_sigma``: median, sigma_iqr, q25, q75, n. |

Raises: `ValueError` (Raised if a bare history is passed without ``bad_bitmask``.)

### session.plot_final_comparison_from_fits

`plot_final_comparison_from_fits(d, view=None, out=None, stats_in_view=True, show_sources=True, show=True, dpi=110, verbose=True)`

Rebuild the final-comparison figure from files written by `save_run`.

No LSST stack and no re-fit needed.

| Parameter | Type | Description |
| --- | --- | --- |
| `d` | `str` | Directory written by `save_run`. |
| `view` | `tuple` [`float`, `float`, `float`, `float`], optional | (xmin, xmax, ymin, ymax) in detector pixels (default: the interior box from summary.json, else the whole array). |
| `out` | `str`, optional | Also save the figure to this PNG. |
| `stats_in_view` | `bool`, optional | Restrict the statistics to the view. |
| `show_sources` | `bool`, optional | Mark the catalogue's sources. |
| `show` | `bool`, optional | Display the figure. |
| `dpi` | `int`, optional | Resolution of the PNG. |
| `verbose` | `bool`, optional | Print the statistics. |

| Returns | Type | Description |
| --- | --- | --- |
| `stats` | `dict` | See `plot_final_comparison_robust`. |

## gaiaCatalog

`lsst.meas.crowd.gaiaCatalog`: the Gaia catalogue and the seeds of each cutout.

Gaia DR3 seed catalogue for a visit.

| Name | Summary |
| --- | --- |
| [`visit_epoch`](#gaiacatalogvisit_epoch) | Return the observation epoch of an exposure. |
| [`gaia_catalog`](#gaiacataloggaia_catalog) | Load the Gaia DR3 stars overlapping a box, with positions at the observation epoch. |
| [`attach_gaia`](#gaiacatalogattach_gaia) | Load Gaia once for a session and store it as ``session.gaia``. |
| [`select_seeds`](#gaiacatalogselect_seeds) | Select the Gaia stars that seed a cutout. |

### gaiaCatalog.visit_epoch

`visit_epoch(exposure)`

Return the observation epoch of an exposure.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Exposure with a VisitInfo. |

| Returns | Type | Description |
| --- | --- | --- |
| `epoch` | `astropy.time.Time` or `None` | The visit date, or None if it is missing or implausible. |

### gaiaCatalog.gaia_catalog

`gaia_catalog(wcs, bbox, epoch=None, source='auto', butler=None, repo='dp2_prep', refcat='gaia_dr3_20230707', refcat_collections='refcats', tap_url='https://gea.esac.esa.int/tap-server/tap', cache_path=None, pad=0, verbose=True)`

Load the Gaia DR3 stars overlapping a box, with positions at the observation epoch.

| Parameter | Type | Description |
| --- | --- | --- |
| `wcs` | `lsst.afw.geom.SkyWcs` | WCS of the image. |
| `bbox` | `lsst.geom.Box2I` | Pixel box (detector pixels of ``wcs``). |
| `epoch` | `astropy.time.Time`, optional | Epoch of the positions (default: the Gaia epoch, J2016.0). |
| `source` | `str`, optional | 'butler' (LSST reference catalogue), 'tap' (ESA archive) or 'auto' (Butler, then TAP). |
| `butler` | `lsst.daf.butler.Butler`, optional | Butler to use (default: one for ``repo``). |
| `repo` | `str`, optional | Butler repository path or alias. |
| `refcat` | `str`, optional | Reference catalogue dataset type. |
| `refcat_collections` | `str`, optional | Collections of the reference catalogue. |
| `tap_url` | `str`, optional | TAP service URL. |
| `cache_path` | `str`, optional | .ecsv file: read if it exists, written otherwise. |
| `pad` | `int`, optional | Grow the box by this many px. |
| `verbose` | `bool`, optional | Print what was loaded. |

| Returns | Type | Description |
| --- | --- | --- |
| `gaia` | `astropy.table.Table` | ``source_id``, ``ra``, ``dec``, ``phot_g_mean_mag``, pixel positions ``x``, ``y`` and further Gaia columns; ``meta['epoch_jyear']`` is the epoch used. |

Raises: `RuntimeError` (Raised if Gaia could not be loaded from any source.)

### gaiaCatalog.attach_gaia

`attach_gaia(session, visit_image=None, pad=64, cache_dir='gaia_cache', **kw)`

Load Gaia once for a session and store it as ``session.gaia``.

Covers the whole visit image or, without one, the union of the session's cutouts; the catalogue is cached on disk.

| Parameter | Type | Description |
| --- | --- | --- |
| `session` | `lsst.meas.crowd.session.Session` | The session (``session.gaia`` is set). |
| `visit_image` | `lsst.afw.image.Exposure`, optional | The full visit image. |
| `pad` | `int`, optional | Grow the box by this many px. |
| `cache_dir` | `str`, optional | Cache directory (None = no cache). |
| `**kw` | – | Keyword arguments of `gaia_catalog` (e.g. ``repo``). |

| Returns | Type | Description |
| --- | --- | --- |
| `gaia` | `astropy.table.Table` | The Gaia catalogue (see `gaia_catalog`). |

### gaiaCatalog.select_seeds

`select_seeds(gaia_tbl, bbox, fwhm, pad_fwhm=1.5, gmag_max=None, min_sep_px=0.5)`

Select the Gaia stars that seed a cutout.

Stars inside the box grown by ``pad_fwhm * fwhm``, brightest first, with near-duplicates (closer than ``min_sep_px``) removed.

| Parameter | Type | Description |
| --- | --- | --- |
| `gaia_tbl` | `astropy.table.Table` or `None` | Gaia catalogue with ``x``, ``y`` (and ``phot_g_mean_mag``). |
| `bbox` | `lsst.geom.Box2I` | Cutout box. |
| `fwhm` | `float` | PSF FWHM [px]. |
| `pad_fwhm` | `float`, optional | Growth of the box [FWHM]. |
| `gmag_max` | `float`, optional | Only stars brighter than this G magnitude. |
| `min_sep_px` | `float`, optional | Drop fainter stars within this many px of a brighter one. |

| Returns | Type | Description |
| --- | --- | --- |
| `xy` | `numpy.ndarray` | Seed positions, shape (N, 2). |
| `rows` | `numpy.ndarray` [`int`] | Their rows in ``gaia_tbl``. |

## photometry

`lsst.meas.crowd.photometry`: PSF and aperture photometry, errors, the Gaia comparison.

Magnitudes of the fitted sources.

| Name | Summary |
| --- | --- |
| [`flux_to_mag`](#photometryflux_to_mag) | Convert fluxes (and errors) to magnitudes. |
| [`aperture_flux`](#photometryaperture_flux) | Sum the flux in a circular aperture. |
| [`aperture_correction`](#photometryaperture_correction) | Return the fraction of a PSF kernel's flux inside a radius. |
| [`fisher_flux_errors`](#photometryfisher_flux_errors) | Return the 1-sigma errors of jointly fitted amplitudes. |
| [`psf_flux_errors`](#photometrypsf_flux_errors) | Return the Fisher errors of the PSF fluxes of a run's final joint fit. |
| [`zero_point_from_run`](#photometryzero_point_from_run) | Return the magnitude zero point of a run's image units. |
| [`compare_with_gaia`](#photometrycompare_with_gaia) | Add Gaia comparison columns to a photometry table. |
| [`fits_meta`](#photometryfits_meta) | Return catalogue metadata with FITS-compatible keywords and values. |
| [`measure_photometry`](#photometrymeasure_photometry) | Measure PSF and aperture magnitudes of a run's sources. |

### photometry.flux_to_mag

`flux_to_mag(flux, flux_err=None, zero_point=31.4)`

Convert fluxes (and errors) to magnitudes.

| Parameter | Type | Description |
| --- | --- | --- |
| `flux` | array-like | Fluxes, in units where ``zero_point`` applies. |
| `flux_err` | array-like, optional | Flux errors. |
| `zero_point` | `float`, optional | Magnitude of unit flux (default: AB for nJy). |

| Returns | Type | Description |
| --- | --- | --- |
| `mag` | `numpy.ndarray` | Magnitudes (NaN for non-positive fluxes). |
| `mag_err` | `numpy.ndarray` | Magnitude errors (only if ``flux_err`` is given). |

### photometry.aperture_flux

`aperture_flux(image, variance, good, x, y, radius, x0=0, y0=0)`

Sum the flux in a circular aperture.

| Parameter | Type | Description |
| --- | --- | --- |
| `image` | `numpy.ndarray` | Image. |
| `variance` | `numpy.ndarray` | Variance image. |
| `good` | `numpy.ndarray` [`bool`] | Usable pixels (others are left out of the sum). |
| `x, y` | `float` | Centre in parent pixel coordinates. |
| `radius` | `float` | Aperture radius [px]. |
| `x0, y0` | `int`, optional | Parent pixel coordinates of ``image[0, 0]``. |

| Returns | Type | Description |
| --- | --- | --- |
| `flux` | `float` | Aperture flux (NaN if the aperture misses the image). |
| `flux_err` | `float` | Its error. |
| `coverage` | `float` | Fraction of the aperture area on good pixels. |

### photometry.aperture_correction

`aperture_correction(kernel, radius)`

Return the fraction of a PSF kernel's flux inside a radius.

| Parameter | Type | Description |
| --- | --- | --- |
| `kernel` | `numpy.ndarray` | PSF kernel (odd size, centred). |
| `radius` | `float` | Aperture radius [px]. |

| Returns | Type | Description |
| --- | --- | --- |
| `fraction` | `float` | Enclosed flux fraction. |

### photometry.fisher_flux_errors

`fisher_flux_errors(comps, variance, good, n=None)`

Return the 1-sigma errors of jointly fitted amplitudes.

``sqrt(diag((A^T W A)^-1))``, inverted per group of overlapping sources, so blends get larger, honest errors.

| Parameter | Type | Description |
| --- | --- | --- |
| `comps` | `list` [`dict` or `None`] | The joint fit's per-source components (``stamp``, ``ys``, ``xs``). |
| `variance` | `numpy.ndarray` | Variance used as weights. |
| `good` | `numpy.ndarray` [`bool`] | Pixels used in the fit. |
| `n` | `int`, optional | Number of sources (default: ``len(comps)``). |

| Returns | Type | Description |
| --- | --- | --- |
| `err` | `numpy.ndarray` | Error per source (inf where a source has no usable pixels). |

### photometry.psf_flux_errors

`psf_flux_errors(run)`

Return the Fisher errors of the PSF fluxes of a run's final joint fit.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |

| Returns | Type | Description |
| --- | --- | --- |
| `err` | `numpy.ndarray` | Error per source (fit variance, including the sys_frac term). |

### photometry.zero_point_from_run

`zero_point_from_run(run)`

Return the magnitude zero point of a run's image units.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |

| Returns | Type | Description |
| --- | --- | --- |
| `zero_point` | `float` | Zero point (from the exposure's PhotoCalib, else 31.4 for nJy). |
| `how` | `str` | Where it came from. |

### photometry.compare_with_gaia

`compare_with_gaia(phot, gaia, match_radius_px=1.0, offset_min_snr=50.0, verbose=True)`

Add Gaia comparison columns to a photometry table.

Every fitted source is paired with a Gaia star: a Gaia-seeded source with the star it was seeded from (``gaia_index``), any other source with the nearest unclaimed Gaia star within ``match_radius_px``.  Positions are compared with Gaia's at the visit epoch (as placed on the detector by the visit WCS); magnitudes with Gaia G.  Gaia G is a much broader band than any LSST filter, so ``gaia_dmag`` contains a zero-point offset and a colour term: the median offset of bright, point-like matches is removed in ``gaia_dmag_rel``, and ``gaia_bp_rp`` is given to fit the colour term.

| Parameter | Type | Description |
| --- | --- | --- |
| `phot` | `astropy.table.Table` | Output of `measure_photometry` (``x``, ``y``, ``psf_flux``, ``psf_mag``, ``gaia_index``; ``ra``, ``dec`` if the exposure had a WCS). |
| `gaia` | `astropy.table.Table` | The Gaia catalogue the run was seeded from (``session.gaia``: ``source_id``, ``ra``, ``dec``, ``x``, ``y``, ``phot_g_mean_mag``, optionally ``phot_bp_mean_mag`` and ``phot_rp_mean_mag``). |
| `match_radius_px` | `float`, optional | Match radius for sources that were not Gaia-seeded [px]. |
| `offset_min_snr` | `float`, optional | Minimum PSF S/N of the matches that set the median magnitude offset. |
| `verbose` | `bool`, optional | Print a summary. |

| Returns | Type | Description |
| --- | --- | --- |
| `table` | `astropy.table.Table` | A copy of ``phot`` with ``gaia_match``, ``gaia_match_type``, ``gaia_source_id``, ``gaia_ra``, ``gaia_dec``, ``gaia_x``, ``gaia_y``, ``gaia_dx``, ``gaia_dy``, ``gaia_dra``, ``gaia_ddec``, ``gaia_sep``, ``gaia_G``, ``gaia_bp_rp``, ``gaia_dmag``, ``gaia_dmag_rel``, ``gaia_flux`` and ``gaia_flux_ratio`` (see doc/catalog_columns.md); ``meta`` gets ``GAIAOFF`` (the median offset), ``GAIANOFF`` (stars used), ``GAIAUNM`` (Gaia stars in the fitted area without a fitted source) and ``GAIASYS``. |

### photometry.fits_meta

`fits_meta(meta)`

Return catalogue metadata with FITS-compatible keywords and values.

Long names are shortened (`FITS_META_KEYS`, else the first 8 characters, upper case); values that are not numbers, strings or booleans, and non-finite numbers, are dropped.

| Parameter | Type | Description |
| --- | --- | --- |
| `meta` | `dict` | Table metadata. |

| Returns | Type | Description |
| --- | --- | --- |
| `meta` | `dict` | The metadata to write. |

### photometry.measure_photometry

`measure_photometry(run, aperture_fwhm=1.5, zero_point=None, verbose=True, gaia=None, gaia_match_px=1.0)`

Measure PSF and aperture magnitudes of a run's sources.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run (from ``run_session_gaia``). |
| `aperture_fwhm` | `float`, optional | Aperture radius [FWHM]; aperture fluxes are corrected to total. |
| `zero_point` | `float`, optional | Override the zero point (default: PhotoCalib, else 31.4 for nJy). |
| `verbose` | `bool`, optional | Print a summary. |
| `gaia` | `astropy.table.Table` or `bool`, optional | The Gaia catalogue for the Gaia comparison columns of `compare_with_gaia`; None (default) uses the run's own (``run.gaia``, the catalogue it was seeded from) if it has one, False leaves the columns out. |
| `gaia_match_px` | `float`, optional | Gaia match radius for sources that were not Gaia-seeded [px]. |

| Returns | Type | Description |
| --- | --- | --- |
| `table` | `astropy.table.Table` | ``id``, ``x``, ``y``, ``ra``, ``dec``, ``gaia``, ``gaia_index``, ``gaia_G``, ``extended``, ``in_interior``, ``psf_flux``, ``psf_flux_err``, ``psf_snr``, ``psf_mag``, ``psf_mag_err``, ``ap_flux``, ``ap_flux_err``, ``ap_coverage``, ``ap_mag``, ``ap_mag_err``, ``ap_minus_psf`` (and, with ``gaia``, the `compare_with_gaia` columns); ``meta`` has the zero point and aperture settings.  The columns and their units are described in ``doc/catalog_columns.md``. |

## plotting

`lsst.meas.crowd.plotting`: final data / model / residual plots and residual diagnostics.

Diagnostic plots of a finished run.

| Name | Summary |
| --- | --- |
| [`show_bright_star_stack`](#plottingshow_bright_star_stack) | Stack the bright stars' residuals to show the pattern they all share. |
| [`show_star_zooms`](#plottingshow_star_zooms) | Show zooms of individual stars: data \| model \| residual \| residual/sigma. |
| [`show_cross_sections`](#plottingshow_cross_sections) | Plot cuts through stars along x, along y and the radial profile. |
| [`compare_runs`](#plottingcompare_runs) | Print side-by-side numbers for several runs of the same cutouts. |
| [`core_excess_diagnostic`](#plottingcore_excess_diagnostic) | Find out why bright stars' cores are sharper than the model. |
| [`plot_residual_gray`](#plottingplot_residual_gray) | Plot the two residual panels of `plot_final_comparison`, larger and in greyscale. |
| [`star_residual_table`](#plottingstar_residual_table) | Find which stars carry the leftover residual, and what they have in common. |
| [`plot_final_comparison`](#plottingplot_final_comparison) | Plot original \| final model \| residual \| residual/sigma of a run. |
| [`plot_gaia_comparison`](#plottingplot_gaia_comparison) | Plot how the fitted photometry and positions compare with Gaia. |

### plotting.show_bright_star_stack

`show_bright_star_stack(run, min_snr=200.0, half=20, plot=True)`

Stack the bright stars' residuals to show the pattern they all share.

The median of the bright, clean, unsaturated stars' residuals, each divided by the star's flux and centred on its sub-pixel position, once with the PSF-only joint model and once with the final model (including the per-star corrections).  A pattern that survives the median is common to all bright stars, i.e. a feature missing from the PSF model; noise and star-specific errors average out.  Also prints the core error of the brighter and fainter half of the stars (a brighter-fatter test) and shows the PSF kernel.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `min_snr` | `float`, optional | Minimum peak S/N of the stacked stars. |
| `half` | `int`, optional | Half-size of the stack [px]. |
| `plot` | `bool`, optional | Show the kernel and the two stacks. |

| Returns | Type | Description |
| --- | --- | --- |
| `stack` | `dict` | ``psf_only`` and ``final`` (stacked residual / flux images), ``n`` (number of stars) and, with at least six stars, ``by_brightness`` (core and ring residual of the brighter and fainter half). |

### plotting.show_star_zooms

`show_star_zooms(run, stars=None, n=8, half=12, sigma='model', grayscale=False, min_snr=30.0, cuts=False)`

Show zooms of individual stars: data | model | residual | residual/sigma.

Every source in the window is marked (o = Gaia-seeded, x = found in the residual), with a line from each Gaia source's catalogue position (+) to its fitted position.  Prints each star's neighbours (offset, flux ratio, Gaia or not, distance from its Gaia position).

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `stars` | `list` [`tuple` [`float`, `float`]], optional | Detector positions; None = the ``n`` worst stars of `star_residual_table`. |
| `n` | `int`, optional | Number of worst stars shown. |
| `half` | `int`, optional | Half-size of each zoom [px]. |
| `sigma` | `str`, optional | 'model', 'effective' or 'pixel' (see `plot_final_comparison`). |
| `grayscale` | `bool`, optional | Residual panels in greyscale. |
| `min_snr` | `float`, optional | Minimum peak S/N when choosing the worst stars. |
| `cuts` | `bool`, optional | Also draw `show_cross_sections` for the same stars. |

| Returns | Type | Description |
| --- | --- | --- |
| `indices` | `list` [`int`] or `None` | Indices of the stars shown (None if there were none). |

### plotting.show_cross_sections

`show_cross_sections(run, stars=None, n=4, half=15, width=1, sigma='model', min_snr=30.0, log=True)`

Plot cuts through stars along x, along y and the radial profile.

Top: data (points) and model (line); bottom: residual with its +-1 sigma band, so a ring shows as a bump in the radial residual at its radius.  Other sources in a cut are marked (dotted: Gaia, dashed: found in the residual).

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `stars` | `list` [`tuple` [`float`, `float`]], optional | Detector positions; None = the ``n`` worst stars of `star_residual_table`. |
| `n` | `int`, optional | Number of worst stars shown. |
| `half` | `int`, optional | Half-length of the cuts [px]. |
| `width` | `int`, optional | Rows / columns averaged for the x / y cuts (odd; 1 = single pixel). |
| `sigma` | `str`, optional | 'model', 'effective' or 'pixel' (see `plot_final_comparison`). |
| `min_snr` | `float`, optional | Minimum peak S/N when choosing the worst stars. |
| `log` | `bool`, optional | Symmetric-log axis for data and model (wings visible). |

| Returns | Type | Description |
| --- | --- | --- |
| `cuts` | `list` [`dict`] or `None` | Per star ``index`` and ``cuts`` (the x, y and radial profiles); None if no stars. |

### plotting.compare_runs

`compare_runs(variants, min_snr=200.0)`

Print side-by-side numbers for several runs of the same cutouts.

For each cutout and variant: robust residual/sigma, number of sources and the stacked PSF-only residual of the bright stars (median residual / flux x 1e5) in radial bins -- a ring shows as a large value in one bin; values closer to 0 mean a better PSF.

| Parameter | Type | Description |
| --- | --- | --- |
| `variants` | `dict` [`str`, `dict`] | Variant name -> runs (as returned by ``run_session_gaia``), e.g. ``{"default": runs, "wider wing join": runs_b}``. |
| `min_snr` | `float`, optional | Minimum peak S/N of the stacked stars. |

### plotting.core_excess_diagnostic

`core_excess_diagnostic(run, min_snr=100.0, plot=True)`

Find out why bright stars' cores are sharper than the model.

For each unsaturated, fairly isolated star with peak S/N >= ``min_snr``, the central excess (data - model summed over the central 3x3 px, divided by the model there) is correlated with flux (sensor effects or a flux-dependent PSF), peak / saturation level (non-linearity) and sub-pixel phase (stamps blurred by sub-pixel shifting).  No trend but a significant excess means the PSF core itself is too broad.  A trend only counts if it is large compared with the excess itself.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `min_snr` | `float`, optional | Minimum peak S/N of the stars used. |
| `plot` | `bool`, optional | Plot the excess against each quantity. |

| Returns | Type | Description |
| --- | --- | --- |
| `result` | `dict` or `None` | Per-star arrays (``idx``, ``flux``, ``snr``, ``peak``, ``excess``, ``excess_psf_only``, ``centre_px``, ``phase``, ``err``), ``correlations`` and ``verdict``; None if fewer than six stars. |

### plotting.plot_residual_gray

`plot_residual_gray(run, **kw)`

Plot the two residual panels of `plot_final_comparison`, larger and in greyscale.

Black = negative, white = positive; areas that cannot be modelled are drawn in blue.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `**kw` | – | Keyword arguments of `plot_final_comparison`. |

| Returns | Type | Description |
| --- | --- | --- |
| `stats` | `dict` | See `plot_final_comparison`. |

### plotting.star_residual_table

`star_residual_table(run, min_snr=30.0, r_fwhm=2.5, plot=True, n_worst=10)`

Find which stars carry the leftover residual, and what they have in common.

For every unsaturated point source with peak S/N >= ``min_snr``: chi2 per pixel within ``r_fwhm`` FWHM for the final and the PSF-only model, the PSF-only core residual (flux within 1 FWHM / star flux; < 0 = star broader than the PSF) and the mismatch (residual flux beyond the noise / star flux, independent of brightness).  The mismatch is rank-correlated with flux, position and crowding, and outlier stars are counted.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `min_snr` | `float`, optional | Minimum peak S/N. |
| `r_fwhm` | `float`, optional | Radius of the chi2 region [FWHM]. |
| `plot` | `bool`, optional | Plot chi2 and core residual against flux and the stars' positions. |
| `n_worst` | `int`, optional | Number of worst stars listed. |

| Returns | Type | Description |
| --- | --- | --- |
| `table` | `dict` or `None` | Per-star arrays (``idx``, ``flux``, ``snr``, ``x``, ``y``, ``dist``, ``n_nb``, ``chi2_final``, ``chi2_psf``, ``core_psf``, ``npix``, ``mismatch``), ``correlations``, ``core_vs_flux``, ``outliers`` and ``verdict``; None if fewer than eight stars. |

### plotting.plot_final_comparison

`plot_final_comparison(run, view='interior', bad_color='#ff00ff', circles='images', show_zones=True, exclude_near_mask=True, excluded_color='#8c8c8c', mask_pad_fwhm=1.5, z_clip=5.0, resid_clip_pct=99.0, figsize=(26, 7), sigma='model', nearsat_color='#ff9900', grayscale=False, panels='all', verbose=True)`

Plot original | final model | residual | residual/sigma of a run.

Bad pixels are drawn opaque in one colour (default magenta) in all four panels and excluded from the statistics.  Saturated-star zones have a dashed outline, extended sources dashed circles.  Pixels whose flux belongs to sources under the mask (which the model cannot fit) are greyed out and excluded too.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `view` | `str` or `tuple`, optional | 'interior' (the requested box), None (the whole fitted array) or (xmin, xmax, ymin, ymax) in detector pixels. |
| `bad_color` | `str`, optional | Colour of bad pixels. |
| `circles` | `str`, optional | Where to draw the source circles: 'images' (data and model), 'all' or 'none'. |
| `show_zones` | `bool`, optional | Outline the saturated-star zones. |
| `exclude_near_mask` | `bool`, optional | Grey out (and exclude) pixels next to bad pixels that cannot be modelled. |
| `excluded_color` | `str`, optional | Colour of those pixels. |
| `mask_pad_fwhm` | `float`, optional | Pad around bad pixels [FWHM]. |
| `z_clip` | `float`, optional | Colour range of residual/sigma. |
| `resid_clip_pct` | `float`, optional | Percentile setting the colour range of the residual. |
| `figsize` | `tuple` [`float`, `float`], optional | Figure size [inches]. |
| `sigma` | `str`, optional | 'model' = sqrt(variance + (eps*model)^2 + k*model) with the measured PSF model error eps; 'effective' = without eps; 'pixel' = the variance plane only. |
| `nearsat_color` | `str`, optional | Colour of near-saturated cores fitted from their wings. |
| `grayscale` | `bool`, optional | Residual panels in greyscale (excluded areas then in blue). |
| `panels` | `str`, optional | 'all' (four panels) or 'residuals' (the two residual panels, larger). |
| `verbose` | `bool`, optional | Print the statistics. |

| Returns | Type | Description |
| --- | --- | --- |
| `stats` | `dict` | ``residual`` and ``residual_over_sigma`` (median, robust sigma, quartiles, number of pixels) and ``excluded_px``. |

### plotting.plot_gaia_comparison

`plot_gaia_comparison(phot, min_snr=20.0, figsize=(22, 5.5), out=None)`

Plot how the fitted photometry and positions compare with Gaia.

Four panels for the Gaia-matched point sources with PSF S/N >= ``min_snr``: magnitude difference (offset removed) against Gaia G and against Gaia BP - RP (the colour term between G and the LSST band), the position offsets, and the offsets as arrows across the image (a pattern there is a WCS distortion; scatter is centroid noise).

| Parameter | Type | Description |
| --- | --- | --- |
| `phot` | `astropy.table.Table` | Output of `lsst.meas.crowd.photometry.compare_with_gaia` (or of ``measure_photometry`` with ``gaia=``). |
| `min_snr` | `float`, optional | Minimum PSF S/N of the plotted sources. |
| `figsize` | `tuple` [`float`, `float`], optional | Figure size [inches]. |
| `out` | `str`, optional | Also save the figure to this file. |

| Returns | Type | Description |
| --- | --- | --- |
| `summary` | `dict` | ``n`` (sources plotted), ``dmag_scatter`` (robust std of ``gaia_dmag_rel`` [mag]), ``offset`` (the removed median offset [mag]) and ``dx``, ``dy`` (median position offsets [mas] or, without sky coordinates, [px]). |

Raises: `ValueError` (Raised if the table has no Gaia comparison columns.)

## stagePlots

`lsst.meas.crowd.stagePlots`: plots of the intermediate stages and pass-to-pass comparison.

Plots of the pipeline's intermediate stages.

| Name | Summary |
| --- | --- |
| [`psf_star_positions`](#stageplotspsf_star_positions) | Return the candidate and used PSF stars of one of the pipeline's PSF determinations. |
| [`show_psf_stars`](#stageplotsshow_psf_stars) | Plot the PSF stars over an image, the PSF kernel and its radial profile. |
| [`show_backgrounds`](#stageplotsshow_backgrounds) | Plot the background entering every outer pass, and the change from pass to pass. |
| [`show_rounds`](#stageplotsshow_rounds) | Plot data \| model \| residual \| residual/sigma for each round of the fit. |
| [`compare_passes`](#stageplotscompare_passes) | Find the sources that one outer pass lost or gained relative to another, and why. |
| [`show_linked_source_pages`](#stageplotsshow_linked_source_pages) | Plot per-source pages whose coloured frames match boxes on the data, model and residual. |

### stagePlots.psf_star_positions

`psf_star_positions(session=None, catalog=None, run=None, outer=-1)`

Return the candidate and used PSF stars of one of the pipeline's PSF determinations.

Exactly one source is used, in this order: a catalogue with ``calib_psf_candidate`` / ``calib_psf_used`` flags (the visit PSF's own stars), a session (the detector-level ``psf_source='refit'`` PSF), or a run (the clean stars the residual-stack PSF update after pass ``outer`` used).

| Parameter | Type | Description |
| --- | --- | --- |
| `session` | `lsst.meas.crowd.session.Session`, optional | A session prepared with ``psf_source='refit'``. |
| `catalog` | `lsst.afw.table.SourceCatalog` or `astropy.table.Table`, optional | A source catalogue with PSF flags, e.g. the visit's ``single_visit_star_unstandardized``. |
| `run` | `lsst.meas.crowd.session.CutoutRun`, optional | A finished run. |
| `outer` | `int`, optional | Outer pass of ``run`` whose result the PSF update used. |

| Returns | Type | Description |
| --- | --- | --- |
| `stars` | `dict` | ``candidates`` and ``used`` (positions in detector px, shape (N, 2)) and ``source`` (a description). |

Raises: `ValueError` (Raised if no source of PSF stars is given.)

### stagePlots.show_psf_stars

`show_psf_stars(image, stars=None, psf=None, at=None, view=None, figsize=(19, 6), out=None)`

Plot the PSF stars over an image, the PSF kernel and its radial profile.

| Parameter | Type | Description |
| --- | --- | --- |
| `image` | `lsst.afw.image.Exposure` | The visit image (or a cutout) to draw the stars on. |
| `stars` | `dict`, optional | Output of `psf_star_positions`; None = draw no stars. |
| `psf` | `lsst.afw.detection.Psf`, optional | PSF to show (default: the image's own PSF). |
| `at` | `tuple` [`float`, `float`], optional | Detector position at which to render the kernel (default: the image centre). |
| `view` | `tuple` [`float`, `float`, `float`, `float`], optional | (xmin, xmax, ymin, ymax) zoom in detector pixels. |
| `figsize` | `tuple` [`float`, `float`], optional | Figure size [inches]. |
| `out` | `str`, optional | Also save the figure to this file. |

| Returns | Type | Description |
| --- | --- | --- |
| `info` | `dict` | ``kernel`` (the rendered kernel image), ``fwhm`` [px], ``n_candidates`` and ``n_used``. |

### stagePlots.show_backgrounds

`show_backgrounds(run, initial=None, figsize_per_panel=4.2, out=None)`

Plot the background entering every outer pass, and the change from pass to pass.

Pass 1 uses ``initial`` (the detector-level background of ``psf_source='refit'``; a visit image is already background-subtracted, so there is none).  Before each later pass the background is refitted on the raw cutout with the previous pass's sources masked.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `initial` | `lsst.afw.image.Image`, optional | The background subtracted before pass 1 (``session.bgImage``). |
| `figsize_per_panel` | `float`, optional | Size of each panel [inches]. |
| `out` | `str`, optional | Also save the figure to this file. |

| Returns | Type | Description |
| --- | --- | --- |
| `backgrounds` | `list` [`numpy.ndarray` or `None`] | Background entering each pass (None = nothing was subtracted). |

### stagePlots.show_rounds

`show_rounds(run, level='outer', outer=-1, rounds=None, max_rows=8, view=None, z_clip=5.0, figsize_per_panel=3.6, out=None)`

Plot data | model | residual | residual/sigma for each round of the fit.

``level='outer'``: one row per outer pass (each with its own background and PSF). ``level='inner'``: one row per inner iteration of outer pass ``outer`` (the source list growing as the detection threshold steps down).  Data and model share one colour scale, and the residual one scale for all rows, so the rows can be compared.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `level` | `str`, optional | 'outer' or 'inner'. |
| `outer` | `int`, optional | Outer pass whose inner iterations are shown (``level='inner'``). |
| `rounds` | `list` [`int`], optional | Rounds to show (default: up to ``max_rows``, evenly spaced, the last included). |
| `max_rows` | `int`, optional | Maximum number of rows when ``rounds`` is not given. |
| `view` | `tuple` [`float`, `float`, `float`, `float`], optional | (xmin, xmax, ymin, ymax) zoom in detector pixels. |
| `z_clip` | `float`, optional | Colour range of residual/sigma. |
| `figsize_per_panel` | `float`, optional | Size of each panel [inches]. |
| `out` | `str`, optional | Also save the figure to this file. |

| Returns | Type | Description |
| --- | --- | --- |
| `rows` | `list` [`dict`] | Per row ``label``, ``n_sources``, ``resid_std`` (robust std of the residual) and ``chi_std`` (robust std of residual/sigma) over the good pixels in view. |

Raises: `ValueError` (Raised if ``level`` is not 'outer' or 'inner'.)

### stagePlots.compare_passes

`compare_passes(run, a=-2, b=-1, match_px=1.0, n_show=8, half=12, plot=True, out=None)`

Find the sources that one outer pass lost or gained relative to another, and why.

Sources are matched by position (within ``match_px``).  For every source of pass ``a`` without a match in pass ``b`` the residual of pass ``b`` is measured at its position with a PSF-weighted (matched-filter) flux: if most of its flux is still there, a real source was lost; if a neighbour within 1.5 FWHM gained that flux, it was merged; if a halo or ring check removed it, it is also blocked from re-detection in later passes.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run with at least two outer passes. |
| `a, b` | `int`, optional | The two outer passes to compare (default: the last two). |
| `match_px` | `float`, optional | Matching radius [px]. |
| `n_show` | `int`, optional | Number of lost sources drawn (brightest first). |
| `half` | `int`, optional | Half-size of each zoom [px]. |
| `plot` | `bool`, optional | Draw the overview and the zooms. |
| `out` | `str`, optional | Also save the figure to this file. |

| Returns | Type | Description |
| --- | --- | --- |
| `table` | `astropy.table.Table` | One row per lost or gained source: ``kind`` ('lost' / 'gained'), ``x``, ``y``, ``flux``, ``snr`` (peak S/N in the pass it exists in), ``gaia``, ``flux_left`` and ``left_snr`` (matched-filter flux and S/N in the residual of the pass without it), ``left_frac``, ``nearest_px`` (distance to the nearest source of the other pass), ``nearest_dflux`` (that source's flux change) and ``verdict``. |

Raises: `ValueError` (Raised if the run has fewer than two outer passes.)

### stagePlots.show_linked_source_pages

`show_linked_source_pages(run, outer=-1, snap=-1, sources=None, source_view='all', sources_per_page=6, max_pages=1, view_half=20, show_pages=True, out_dir=None, pdf=None, dpi=80, figsize=(21, 7.5), out=None)`

Plot per-source pages whose coloured frames match boxes on the data, model and residual.

For each page: an overview (data | model | residual of the whole cutout) with one coloured box per source on that page -- the window its row shows -- then the page itself: one row per source with data | model | residual | x cross-section | y cross-section, framed in the same colour.

For many pages use ``out_dir``: every figure is written to its own PNG file and closed at once, so memory stays flat.  A `PdfPages` (``pdf``) keeps every image in memory until it is closed (matplotlib), which can exhaust a notebook's memory for hundreds of sources.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `outer` | `int`, optional | Outer pass. |
| `snap` | `int`, optional | Inner-iteration snapshot of that pass; only the last one keeps the per-source stamps unless the run used ``keep_all_fits=True``. |
| `sources` | `list`, optional | Source indices, or (x, y) detector positions (the nearest source is used); None = choose with ``source_view``. |
| `source_view` | `str`, optional | 'all' (brightest first), 'new' or 'changed' (see `fitLib.select_source_pages`). |
| `sources_per_page` | `int`, optional | Rows per page. |
| `max_pages` | `int`, optional | Maximum number of pages (None = all). |
| `view_half` | `int`, optional | Half-size of each source window [px] (the size of the boxes). |
| `show_pages` | `bool`, optional | Draw the pages (False = only the overviews). |
| `out_dir` | `str`, optional | Write ``overview_NNN.png`` and ``page_NNN.png`` here instead of showing them. |
| `pdf` | `matplotlib.backends.backend_pdf.PdfPages`, optional | Open PDF; the overviews and pages are saved into it instead of being shown. |
| `dpi` | `int`, optional | Resolution of the saved pages. |
| `figsize` | `tuple` [`float`, `float`], optional | Size of the overview [inches]. |
| `out` | `str`, optional | Also save the first overview to this file. |

| Returns | Type | Description |
| --- | --- | --- |
| `linked` | `dict` | ``pages`` (source indices per page), ``colors`` (colour per source index) and ``files`` (the PNG files written with ``out_dir``). |

Raises: `ValueError` (Raised if the snapshot no longer has the per-source stamps.)

## sourcePlots

`lsst.meas.crowd.sourcePlots`: per-source cutouts, cross-sections and magnitudes.

Per-source cutouts, cross-sections and magnitudes.

| Name | Summary |
| --- | --- |
| [`source_cutout_data`](#sourceplotssource_cutout_data) | Compute the arrays one cutout row shows for a source. |
| [`select_sources`](#sourceplotsselect_sources) | Choose the sources to show. |
| [`plot_source_cutouts`](#sourceplotsplot_source_cutouts) | Plot cutouts and cross-sections of individual sources, labelled with their magnitudes. |

### sourcePlots.source_cutout_data

`source_cutout_data(data, model, var, bad, x, y, x0=0, y0=0, half=15, width=1)`

Compute the arrays one cutout row shows for a source.

| Parameter | Type | Description |
| --- | --- | --- |
| `data` | `numpy.ndarray` | Image. |
| `model` | `numpy.ndarray` | Model image. |
| `var` | `numpy.ndarray` | Variance image. |
| `bad` | `numpy.ndarray` [`bool`] | Bad pixels. |
| `x, y` | `float` | Source position in parent pixel coordinates. |
| `x0, y0` | `int`, optional | Parent pixel coordinates of ``data[0, 0]``. |
| `half` | `int`, optional | Half-size of the window [px]. |
| `width` | `int`, optional | Rows / columns averaged in the x / y cuts. |

| Returns | Type | Description |
| --- | --- | --- |
| `cutout` | `dict` | ``window`` (xmin, xmax, ymin, ymax), the ``data``, ``model``, ``resid``, ``z`` (residual / sigma) and ``bad`` arrays, ``xcut`` / ``ycut`` / ``radial`` (each a dict with the coordinate, ``data``, ``model``, ``resid`` and ``sigma``) and ``chi2`` (per point of each cut). |

### sourcePlots.select_sources

`select_sources(phot, how='worst', n=6, positions=None, xy=None, chi2=None, min_snr=20.0, seed=0)`

Choose the sources to show.

| Parameter | Type | Description |
| --- | --- | --- |
| `phot` | `astropy.table.Table` or `dict` | Photometry with ``psf_snr``. |
| `how` | `str`, optional | 'worst' (highest ``chi2``), 'brightest', 'faintest' or 'random', among sources with S/N >= ``min_snr``. |
| `n` | `int`, optional | Number of sources. |
| `positions` | `list` [`tuple` [`float`, `float`]], optional | Show the nearest source to each position instead. |
| `xy` | `numpy.ndarray`, optional | Source positions (needed with ``positions``). |
| `chi2` | `numpy.ndarray`, optional | Local chi2 per source (needed for 'worst'). |
| `min_snr` | `float`, optional | Minimum S/N. |
| `seed` | `int`, optional | Random seed for 'random'. |

| Returns | Type | Description |
| --- | --- | --- |
| `indices` | `list` [`int`] | Selected source indices. |

Raises: `ValueError` (Raised if ``how`` is not a known selection.)

### sourcePlots.plot_source_cutouts

`plot_source_cutouts(run, how='worst', n=6, positions=None, half=15, width=1, sigma='model', phot=None, min_snr=20.0, log=True, per_page=4, pdf_path=None, show=True, grayscale=False)`

Plot cutouts and cross-sections of individual sources, labelled with their magnitudes.

One row per source: data | model | residual | residual/sigma | x-cut | y-cut | radial profile.  The cuts show data (points) and model (line), and the residual with its +-1 sigma band; other sources in the window are marked (o Gaia, x found in the residual).

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `how` | `str`, optional | 'worst', 'brightest', 'faintest' or 'random' (see `select_sources`). |
| `n` | `int`, optional | Number of sources. |
| `positions` | `list` [`tuple` [`float`, `float`]], optional | Show the sources nearest these detector positions. |
| `half` | `int`, optional | Half-size of each window [px]. |
| `width` | `int`, optional | Rows / columns averaged in the cuts. |
| `sigma` | `str`, optional | 'model' (variance plus measured model error), 'effective' or 'pixel'. |
| `phot` | `astropy.table.Table`, optional | Output of `lsst.meas.crowd.photometry.measure_photometry` (computed if not given). |
| `min_snr` | `float`, optional | Minimum S/N of the sources selected. |
| `log` | `bool`, optional | Symmetric-log axis for data and model. |
| `per_page` | `int`, optional | Sources per figure. |
| `pdf_path` | `str`, optional | Also write all pages to this PDF. |
| `show` | `bool`, optional | Display the figures. |
| `grayscale` | `bool`, optional | Residual panels in greyscale. |

| Returns | Type | Description |
| --- | --- | --- |
| `rows` | `astropy.table.Table` | Photometry rows of the sources shown. |

## detector

`lsst.meas.crowd.detector`: whole-detector runs: tiles, stitching, mosaic and tile plots.

Whole-detector runs: tile, fit, stitch.

| Name | Summary |
| --- | --- |
| [`tile_specs`](#detectortile_specs) | Return the interiors of a tile grid that partitions an area exactly. |
| [`tile_name`](#detectortile_name) | Return the name of a tile, e.g. 'x0400_y0000_400x400'. |
| [`save_tile_products`](#detectorsave_tile_products) | Write one tile's results to a FITS file. |
| [`fit_tiles`](#detectorfit_tiles) | Fit this task's share of a detector's tiles and write one file per tile. |
| [`stitch_tiles`](#detectorstitch_tiles) | Assemble the tile files into full-detector images and one catalogue. |
| [`load_mosaic`](#detectorload_mosaic) | Read the files written by `stitch_tiles`. |
| [`load_tile`](#detectorload_tile) | Read one tile file written by `save_tile_products`. |
| [`plot_tile`](#detectorplot_tile) | Plot data \| model \| residual \| residual/sigma of one tile over its whole fit box. |
| [`plot_detector`](#detectorplot_detector) | Plot data \| model \| residual \| residual/sigma of a stitched detector or a zoom. |

### detector.tile_specs

`tile_specs(width, height, tile=400, x0=0, y0=0, region=None)`

Return the interiors of a tile grid that partitions an area exactly.

``ceil(size / tile)`` tiles per axis with sizes as equal as possible.

| Parameter | Type | Description |
| --- | --- | --- |
| `width, height` | `int` | Size of the area [px]. |
| `tile` | `int`, optional | Target tile size [px]. |
| `x0, y0` | `int`, optional | Origin of the area. |
| `region` | `tuple` [`int`, `int`, `int`, `int`], optional | (xmin, xmax, ymin, ymax): keep only tiles overlapping it. |

| Returns | Type | Description |
| --- | --- | --- |
| `specs` | `list` [`tuple` [`int`, `int`, `int`, `int`]] | (x0, y0, w, h) per tile. |

### detector.tile_name

`tile_name(spec)`

Return the name of a tile, e.g. 'x0400_y0000_400x400'.

| Parameter | Type | Description |
| --- | --- | --- |
| `spec` | `tuple` [`int`, `int`, `int`, `int`] | (x0, y0, w, h). |

| Returns | Type | Description |
| --- | --- | --- |
| `name` | `str` | The tile name. |

### detector.save_tile_products

`save_tile_products(run, path, phot=None, meta=None)`

Write one tile's results to a FITS file.

DATA, MODEL, VARIANCE and BAD images over the tile's fit box (header keys X0/Y0 = fit-box origin, IX0/IY0/IX1/IY1 = interior, inclusive) and a CATALOG table (the photometry, with an ``owned`` column: source inside this tile's interior, half-open so boundary sources belong to exactly one tile).  The file is written atomically.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | The tile's finished run. |
| `path` | `str` | Output file. |
| `phot` | `astropy.table.Table`, optional | Photometry (computed if not given). |
| `meta` | `dict`, optional | Extra header keywords. |

| Returns | Type | Description |
| --- | --- | --- |
| `path` | `str` | The file written. |

### detector.fit_tiles

`fit_tiles(visit, detector, config=None, out_dir='detector_run', tile=400, pad=30, task_index=0, n_tasks=1, region=None, repo='dp2_prep', collection='LSSTCam/runs/DRP/DP2', instrument='LSSTCam', visit_image_collections='LSSTCam/runs/*', psf_source='visit', overwrite=False, visit_image=None, **overrides)`

Fit this task's share of a detector's tiles and write one file per tile.

Tiles already on disk are skipped, so the job can be restarted; a failing tile is reported and the task moves on.

| Parameter | Type | Description |
| --- | --- | --- |
| `visit, detector` | `int` | Visit and detector ids. |
| `config` | `dict`, optional | Pipeline configuration (default ``make_config(CROWDED_FIELD_CONFIG)``). |
| `out_dir` | `str`, optional | Output directory. |
| `tile, pad` | `int`, optional | Tile interior size and fit padding [px]. |
| `task_index, n_tasks` | `int`, optional | This task fits tiles ``task_index``, ``task_index + n_tasks``, ... (one Slurm array task each). |
| `region` | `tuple` [`int`, `int`, `int`, `int`], optional | (xmin, xmax, ymin, ymax): only tiles overlapping it (for tests). |
| `repo` | `str`, optional | Butler repository path or alias. |
| `collection` | `str`, optional | Default Butler collection. |
| `instrument` | `str`, optional | Instrument name. |
| `visit_image_collections` | `str`, optional | Collections searched for the visit image. |
| `psf_source` | `str`, optional | 'visit' or 'refit' (see `lsst.meas.crowd.session.build_detector_psf`). |
| `overwrite` | `bool`, optional | Refit tiles that already have a file. |
| `visit_image` | `lsst.afw.image.Exposure`, optional | An already-loaded visit image (skips the Butler). |
| `**overrides` | – | Keyword arguments of ``run_session_gaia`` (e.g. ``max_outer=3, max_inner=15``). |

| Returns | Type | Description |
| --- | --- | --- |
| `paths` | `list` [`str`] | Tile files of this task that exist. |

### detector.stitch_tiles

`stitch_tiles(out_dir='detector_run', mosaic_name='detector_mosaic.fits', catalog_name='detector_catalog.fits', verbose=True)`

Assemble the tile files into full-detector images and one catalogue.

Writes the mosaic (DATA, MODEL, RESIDUAL, VARIANCE, BAD and TILE extensions) and the catalogue (every source once, from the tile whose interior contains it) into ``out_dir``.

| Parameter | Type | Description |
| --- | --- | --- |
| `out_dir` | `str`, optional | Directory with the tile files. |
| `mosaic_name` | `str`, optional | Mosaic file name. |
| `catalog_name` | `str`, optional | Catalogue file name. |
| `verbose` | `bool`, optional | Print a summary. |

| Returns | Type | Description |
| --- | --- | --- |
| `mosaic` | `dict` | ``data``, ``model``, ``resid``, ``var``, ``bad``, ``tile`` (tile index map, -1 = not fitted), ``x0``, ``y0``, ``tiles`` (names), ``interiors`` and ``catalog``. |

Raises: `FileNotFoundError` (Raised if ``out_dir`` has no tile files.)

### detector.load_mosaic

`load_mosaic(out_dir='detector_run', mosaic_name='detector_mosaic.fits', catalog_name='detector_catalog.fits')`

Read the files written by `stitch_tiles`.

| Parameter | Type | Description |
| --- | --- | --- |
| `out_dir` | `str`, optional | Directory. |
| `mosaic_name` | `str`, optional | Mosaic file name. |
| `catalog_name` | `str`, optional | Catalogue file name. |

| Returns | Type | Description |
| --- | --- | --- |
| `mosaic` | `dict` | As returned by `stitch_tiles` (without ``tiles`` and ``interiors``). |

### detector.load_tile

`load_tile(path)`

Read one tile file written by `save_tile_products`.

| Parameter | Type | Description |
| --- | --- | --- |
| `path` | `str` | The ``tile_<name>.fits`` file. |

| Returns | Type | Description |
| --- | --- | --- |
| `tile` | `dict` | In the form of `load_mosaic`'s result over the tile's fit box (``data``, ``model``, ``resid``, ``var``, ``bad``, ``tile``, ``x0``, ``y0``, ``catalog``), plus ``name`` and ``interiors`` (the tile's interior as one (x0, y0, x1, y1) tuple, inclusive). |

### detector.plot_tile

`plot_tile(tile, z_clip=5.0, grayscale=False, sources=True, figsize=(24, 6.5), out=None)`

Plot data | model | residual | residual/sigma of one tile over its whole fit box.

The tile's interior (the part kept in the mosaic) is outlined; outside it is the padding, fitted only so that sources near the edge have their neighbours.

| Parameter | Type | Description |
| --- | --- | --- |
| `tile` | `dict` or `str` | Output of `load_tile`, or the tile file. |
| `z_clip` | `float`, optional | Colour range of residual/sigma. |
| `grayscale` | `bool`, optional | Residual panels in greyscale. |
| `sources` | `bool`, optional | Circle the tile catalogue's sources. |
| `figsize` | `tuple` [`float`, `float`], optional | Figure size [inches]. |
| `out` | `str`, optional | Also save the figure to this file. |

| Returns | Type | Description |
| --- | --- | --- |
| `result` | `dict` | As `plot_detector`'s. |

### detector.plot_detector

`plot_detector(mosaic, view=None, binning=None, show_tiles=True, z_clip=5.0, grayscale=False, figsize=(26, 7.5), sources=False, out=None)`

Plot data | model | residual | residual/sigma of a stitched detector or a zoom.

| Parameter | Type | Description |
| --- | --- | --- |
| `mosaic` | `dict` or `str` | Output of `stitch_tiles` / `load_mosaic`, or the directory to load. |
| `view` | `tuple` [`float`, `float`, `float`, `float`], optional | (xmin, xmax, ymin, ymax) zoom in detector pixels. |
| `binning` | `int`, optional | Display block size (default: about 1000 px across the view); residual/sigma is binned as sum / sqrt(n) so its noise stays 1. |
| `show_tiles` | `bool`, optional | Draw the tile boundaries. |
| `z_clip` | `float`, optional | Colour range of residual/sigma. |
| `grayscale` | `bool`, optional | Residual panels in greyscale. |
| `figsize` | `tuple` [`float`, `float`], optional | Figure size [inches]. |
| `sources` | `bool`, optional | Circle the catalogue's sources (for zooms). |
| `out` | `str`, optional | Also save the figure to this file. |

| Returns | Type | Description |
| --- | --- | --- |
| `result` | `dict` | ``residual_sigma`` (robust sigma of residual/sigma over the view) and ``binning``. |

## fastKernels

`lsst.meas.crowd.fastKernels`: fast PSF stamps (C++ with a NumPy fallback).

Fast PSF stamps for a spatially CONSTANT PSF.

| Name | Summary |
| --- | --- |
| [`lanczos_weights`](#fastkernelslanczos_weights) | Return normalised Lanczos weights for a sub-pixel shift. |
| [`shift_image`](#fastkernelsshift_image) | Shift an image by a sub-pixel offset with a separable Lanczos kernel. |
| [`taper_stamp`](#fastkernelstaper_stamp) | Tame the wings of a PSF stamp like ``fitLib.psf_stamp(taper=...)``. |
| [`render_stamps`](#fastkernelsrender_stamps) | Render stamps of a constant PSF at many positions. |
| [`render_derivatives`](#fastkernelsrender_derivatives) | Render central-difference derivative stamps dP/dx, dP/dy of a constant PSF. |
| [`FastPsf`](#fastkernelsfastpsf) (class) | Stamps and derivatives of a spatially constant PSF. |

### fastKernels.lanczos_weights

`lanczos_weights(frac, order=5)`

Return normalised Lanczos weights for a sub-pixel shift.

| Parameter | Type | Description |
| --- | --- | --- |
| `frac` | `float` | Shift [px]. |
| `order` | `int`, optional | Lanczos order. |

| Returns | Type | Description |
| --- | --- | --- |
| `first` | `int` | Offset of the first weight: ``out[i] = sum_j w[j] * in[i + first + j]`` approximates ``in(i - frac)``. |
| `weights` | `numpy.ndarray` | The ``2 * order`` weights (summing to 1). |

### fastKernels.shift_image

`shift_image(image, dx, dy, order=5, use_cxx=None)`

Shift an image by a sub-pixel offset with a separable Lanczos kernel.

Same convention as ``lsst.afw.math.offsetImage``: ``out(x, y) = in(x - dx, y - dy)``; pixels outside the input count as zero.

| Parameter | Type | Description |
| --- | --- | --- |
| `image` | `numpy.ndarray` | Image. |
| `dx, dy` | `float` | Shift [px]. |
| `order` | `int`, optional | Lanczos order. |
| `use_cxx` | `bool`, optional | Use the C++ kernels (default: if built). |

| Returns | Type | Description |
| --- | --- | --- |
| `shifted` | `numpy.ndarray` | The shifted image. |

### fastKernels.taper_stamp

`taper_stamp(stamp, cx, cy, taper, use_cxx=None)`

Tame the wings of a PSF stamp like ``fitLib.psf_stamp(taper=...)``.

Subtract the median of the outer pixel ring, then apply a raised cosine from ``taper[0] * h`` to ``taper[1] * h`` around (cx, cy), h being half the stamp size.

| Parameter | Type | Description |
| --- | --- | --- |
| `stamp` | `numpy.ndarray` | PSF stamp. |
| `cx, cy` | `float` | Centre in stamp pixel coordinates. |
| `taper` | `tuple` [`float`, `float`] | Inner and outer taper radii as fractions of h. |
| `use_cxx` | `bool`, optional | Use the C++ kernels (default: if built). |

| Returns | Type | Description |
| --- | --- | --- |
| `stamp` | `numpy.ndarray` | The tapered stamp. |

### fastKernels.render_stamps

`render_stamps(kernel, xs, ys, order=5, taper=None, use_cxx=None)`

Render stamps of a constant PSF at many positions.

Stamp i has origin (round(x) - nx//2, round(y) - ny//2) and is the kernel shifted by the position's sub-pixel offset, like ``psf.computeImage(Point2D(x, y))``.

| Parameter | Type | Description |
| --- | --- | --- |
| `kernel` | `numpy.ndarray` | PSF kernel image (odd size, centred). |
| `xs, ys` | array-like | Positions in parent pixel coordinates. |
| `order` | `int`, optional | Lanczos order. |
| `taper` | `tuple` [`float`, `float`], optional | Wing taper (see `taper_stamp`); None = no taper. |
| `use_cxx` | `bool`, optional | Use the C++ kernels (default: if built). |

| Returns | Type | Description |
| --- | --- | --- |
| `stamps` | `numpy.ndarray` | Shape (n, ny, nx). |
| `x0, y0` | `numpy.ndarray` [`int`] | Stamp origins. |

### fastKernels.render_derivatives

`render_derivatives(kernel, xs, ys, eps=0.25, order=5, use_cxx=None)`

Render central-difference derivative stamps dP/dx, dP/dy of a constant PSF.

| Parameter | Type | Description |
| --- | --- | --- |
| `kernel` | `numpy.ndarray` | PSF kernel image (odd size, centred). |
| `xs, ys` | array-like | Positions in parent pixel coordinates. |
| `eps` | `float`, optional | Finite-difference step [px]. |
| `order` | `int`, optional | Lanczos order. |
| `use_cxx` | `bool`, optional | Use the C++ kernels (default: if built). |

| Returns | Type | Description |
| --- | --- | --- |
| `dx, dy` | `numpy.ndarray` | Derivative stamps, shape (n, ny, nx), on the centred stamps' origins. |
| `x0, y0` | `numpy.ndarray` [`int`] | Stamp origins. |

### fastKernels.FastPsf

Class.  `FastPsf(kernel, order=5)`

Stamps and derivatives of a spatially constant PSF.

Return conventions match ``fitLib.psf_stamp`` and ``fitLib.compute_psf_derivatives``, with boxes given as (xmin, ymin, xmax, ymax) tuples.

| Parameter | Type | Description |
| --- | --- | --- |
| `kernel` | `numpy.ndarray` | PSF kernel image (odd size, centred), e.g. ``psf.computeKernelImage(point).array``. |
| `order` | `int`, optional | Lanczos order (afw renders PSF images with lanczos5). |

#### fastKernels.FastPsf.stamp

`stamp(self, x, y, clip, taper=None)`

Render the stamp at (x, y), clipped to a box.

| Parameter | Type | Description |
| --- | --- | --- |
| `x, y` | `float` | Position in parent pixel coordinates. |
| `clip` | `tuple` [`int`, `int`, `int`, `int`] | Clip box (xmin, ymin, xmax, ymax), inclusive. |
| `taper` | `tuple` [`float`, `float`], optional | Wing taper (see `taper_stamp`). |

| Returns | Type | Description |
| --- | --- | --- |
| `stamp` | `numpy.ndarray` or `None` | The clipped stamp (None if it misses ``clip``). |
| `box` | `tuple` [`int`, `int`, `int`, `int`] or `None` | Its box (xmin, ymin, xmax, ymax). |

#### fastKernels.FastPsf.derivatives

`derivatives(self, x, y, clip, eps=0.25)`

Render the derivative stamps at (x, y), clipped to a box.

| Parameter | Type | Description |
| --- | --- | --- |
| `x, y` | `float` | Position in parent pixel coordinates. |
| `clip` | `tuple` [`int`, `int`, `int`, `int`] | Clip box (xmin, ymin, xmax, ymax), inclusive. |
| `eps` | `float`, optional | Finite-difference step [px]. |

| Returns | Type | Description |
| --- | --- | --- |
| `dpdx, dpdy` | `numpy.ndarray` or `None` | Derivative stamps (None if they miss ``clip``). |
| `box` | `tuple` [`int`, `int`, `int`, `int`] or `None` | Their box (xmin, ymin, xmax, ymax). |

#### fastKernels.FastPsf.validate

`validate(self, reference, positions, rtol=0.002)`

Compare the fast stamps with a reference renderer.

| Parameter | Type | Description |
| --- | --- | --- |
| `reference` | callable | ``reference(x, y) -> (array, (xmin, ymin))``: the full, untapered stamp (e.g. from ``psf.computeImage``). |
| `positions` | iterable [`tuple` [`float`, `float`]] | Positions to compare at. |
| `rtol` | `float`, optional | Maximum allowed difference relative to the peak. |

| Returns | Type | Description |
| --- | --- | --- |
| `ok` | `bool` | Whether all stamps agree (same origin and shape, difference within ``rtol``). |
| `worst` | `float` | Largest difference relative to the peak. |

## fitLib

`lsst.meas.crowd.fitLib`: the original cutout_fit_lib: joint fit, detection, PSF and background tools, plots.

Crowded-field iterative joint PSF fitting (the original cutout_fit_lib.py).

| Name | Summary |
| --- | --- |
| [`bbox_to_slices`](#fitlibbbox_to_slices) | Convert a box into array slices relative to a parent box. |
| [`collapse_saturated_peaks`](#fitlibcollapse_saturated_peaks) | Collapse every peak on or near a saturated blob into one peak per blob. |
| [`merge_close_peaks`](#fitlibmerge_close_peaks) | Merge peaks closer than ``min_separation``, keeping the brightest pixel's peak. |
| [`joint_psf_fit`](#fitlibjoint_psf_fit) | Fit all PSF amplitudes simultaneously with a non-negative least-squares solve. |
| [`iterative_joint_detection`](#fitlibiterative_joint_detection) | Detect, fit and subtract iteratively until the residual stops yielding new sources. |
| [`plot_tiled_comparison`](#fitlibplot_tiled_comparison) | Plot data, model and residual (and a row profile) for every tile of an image. |
| [`plot_detections_by_iteration`](#fitlibplot_detections_by_iteration) | Plot every fitted peak over the image, coded by the pass that first found it. |
| [`mask_sources_from_model`](#fitlibmask_sources_from_model) | Flag pixels where the joint-fit model has significant flux as DETECTED. |
| [`refit_background`](#fitlibrefit_background) | Re-derive the background ignoring masked pixels and subtract it in place. |
| [`catalog_from_peaks`](#fitlibcatalog_from_peaks) | Build a measured SourceCatalog from raw (x, y) peak positions. |
| [`redetect_and_refit_psf`](#fitlibredetect_and_refit_psf) | Detect, deblend, measure, select stars and determine a new PSF. |
| [`outer_loop_psf_background_refinement`](#fitlibouter_loop_psf_background_refinement) | Alternate inner-loop fits with background and PSF re-derivation. |
| [`plot_psf_models_by_outer_pass`](#fitlibplot_psf_models_by_outer_pass) | Plot the PSF kernel image used in each outer pass. |
| [`plot_background_by_outer_pass`](#fitlibplot_background_by_outer_pass) | Plot the background re-derived going into each outer pass. |
| [`plot_sources_overlaid_by_pass`](#fitlibplot_sources_overlaid_by_pass) | Plot the data of each outer pass with its converged peak positions. |
| [`plot_outer_loop_summary`](#fitlibplot_outer_loop_summary) | Plot the residual standard deviation and the source count against outer pass. |
| [`plot_psf_fwhm_by_outer_pass`](#fitlibplot_psf_fwhm_by_outer_pass) | Plot the FWHM of the PSF used in each outer pass. |
| [`StepTimer`](#fitlibsteptimer) (class) | Collect timing records of pipeline steps. |
| [`plot_timing_summary`](#fitlibplot_timing_summary) | Plot the total time per step and the time per inner iteration. |
| [`make_cutout`](#fitlibmake_cutout) | Cut a sub-exposure out of a crop, in coordinates relative to the crop's corner. |
| [`psf_stamp`](#fitlibpsf_stamp) | Render the PSF at (x, y), clipped to a box, with its far wings tamed. |
| [`drop_masked_peaks`](#fitlibdrop_masked_peaks) | Drop peaks on EDGE / NO_DATA pixels. |
| [`drop_out_of_bounds`](#fitlibdrop_out_of_bounds) | Drop peaks outside the exposure's box grown by ``margin`` px. |
| [`effective_variance`](#fitlibeffective_variance) | Return the variance used to weight the fit. |
| [`reject_near_bright`](#fitlibreject_near_bright) | Drop new peaks within ``min_sep`` px of an existing bright source. |
| [`prune_and_merge`](#fitlibprune_and_merge) | Drop low-S/N sources, then merge pairs closer than ``min_sep`` px. |
| [`empirical_noise_scale`](#fitlibempirical_noise_scale) | Return the empirical noise scale of a fit's residual. |
| [`amp_snr_full`](#fitlibamp_snr_full) | Return honest amplitude S/N values, including the covariance of overlapping sources. |
| [`eliminate_weak`](#fitlibeliminate_weak) | Remove sources below ``min_snr`` by backward elimination. |
| [`shape_metrics`](#fitlibshape_metrics) | Measure how PSF-like each source is. |
| [`judge_shape`](#fitlibjudge_shape) | Decide which sources fail the PSF-shape tests. |
| [`eliminate_nonpsf`](#fitlibeliminate_nonpsf) | Remove sources that fail the PSF-shape test by backward elimination. |
| [`plot_shape_rejects`](#fitlibplot_shape_rejects) | Plot the sources removed by the shape cut over the data. |
| [`select_source_pages`](#fitlibselect_source_pages) | Choose the sources that get per-source pages. |
| [`source_colors`](#fitlibsource_colors) | Assign each source on the pages a colour from ``SRC_PALETTE``. |
| [`joint_psf_fit_diag`](#fitlibjoint_psf_fit_diag) | Fit all PSF amplitudes in one non-negative solve, keeping what is needed for diagnostics. |
| [`compute_psf_derivatives`](#fitlibcompute_psf_derivatives) | Compute central-difference derivative stamps dP/dx and dP/dy. |
| [`windowed_centroid`](#fitlibwindowed_centroid) | Return the flux-weighted centroid in a box around (x, y). |
| [`refine_peak_positions`](#fitlibrefine_peak_positions) | Move each well-detected source towards the position that minimises its residual. |
| [`redetect_and_refit_psf_adaptive`](#fitlibredetect_and_refit_psf_adaptive) | Run `redetect_and_refit_psf`, trying kernel sizes from largest to smallest. |
| [`plot_iteration_overview`](#fitlibplot_iteration_overview) | Plot the overview of one inner iteration. |
| [`run_inner_loop_verbose`](#fitlibrun_inner_loop_verbose) | Run the instrumented inner loop (detect, fit, refine, prune, repeat). |
| [`plot_amplitude_evolution`](#fitlibplot_amplitude_evolution) | Plot the amplitude of every final source against inner iteration. |
| [`plot_source_evolution`](#fitlibplot_source_evolution) | Plot model and residual cross-sections of three sources at every iteration. |
| [`plot_cutout_setup`](#fitlibplot_cutout_setup) | Plot the cutout before fitting: data, sigma map, mask planes and PSF. |
| [`plot_outer_step_background`](#fitlibplot_outer_step_background) | Plot the background refit of one outer step. |
| [`plot_psf_candidate_stamps`](#fitlibplot_psf_candidate_stamps) | Plot the PSF candidate stamps, used ones first. |
| [`plot_outer_step_psf`](#fitlibplot_outer_step_psf) | Plot the PSF refit of one outer step: star selection and old vs new PSF. |
| [`highfreq_noise_check_arrays`](#fitlibhighfreq_noise_check_arrays) | Compare the residual's pixel-to-pixel scatter with the variance plane. |
| [`highfreq_noise_check`](#fitlibhighfreq_noise_check) | Run `highfreq_noise_check_arrays` on one snapshot of a run. |
| [`estimate_psf_correction`](#fitlibestimate_psf_correction) | Estimate an improved PSF kernel from neighbour-subtracted residuals. |
| [`kernel_to_lsst_psf`](#fitlibkernel_to_lsst_psf) | Wrap a kernel array as a spatially constant LSST PSF. |
| [`extend_kernel_with_moffat`](#fitlibextend_kernel_with_moffat) | Replace the noisy outer part of a PSF kernel by an analytic Moffat tail. |
| [`extend_psf_wings`](#fitlibextend_psf_wings) | Return a constant KernelPsf with a Moffat tail, taken from any PSF at a point. |
| [`refine_psf_from_fit`](#fitlibrefine_psf_from_fit) | Make one PSF update step from a finished run. |
| [`polish_psf_fixed_sources`](#fitlibpolish_psf_fixed_sources) | Iterate PSF correction and re-solve without re-detecting sources. |
| [`plot_psf_correction`](#fitlibplot_psf_correction) | Plot the old kernel, the proposed correction and the radial profiles. |
| [`run_cutout_pipeline`](#fitlibrun_cutout_pipeline) | Run the outer loop (inner fit, background and PSF refinement) on a cutout. |
| [`plot_wing_diagnostic`](#fitlibplot_wing_diagnostic) | Plot the azimuthal profiles of the brightest sources: data vs their own model. |
| [`plot_source_cutouts`](#fitlibplot_source_cutouts) | Plot cutouts and cross-sections of individual sources of a run. |
| [`negative_source_check`](#fitlibnegative_source_check) | Estimate the number of false detections from the residual itself. |
| [`source_snr_histogram`](#fitlibsource_snr_histogram) | Plot the histogram of the honest S/N of every final source. |
| [`draw_edge_margin`](#fitlibdraw_edge_margin) | Outline and shade the band within ``margin_px`` of the image border. |
| [`plot_final_comparison`](#fitlibplot_final_comparison) | Plot original \| final model \| residual \| residual/sigma over the fitted region. |
| [`plot_source_pages`](#fitlibplot_source_pages) | Plot one row per source: data \| model \| residual \| x-section \| y-section. |
| [`source_metrics`](#fitlibsource_metrics) | Return per-source quantities used to stratify the sample. |
| [`select_representative_sources`](#fitlibselect_representative_sources) | Select a representative sample of sources per category. |
| [`plot_representative_sources`](#fitlibplot_representative_sources) | Select a representative sample and plot it with `plot_source_pages`. |

### fitLib.bbox_to_slices

`bbox_to_slices(bbox, parent_bbox)`

Convert a box into array slices relative to a parent box.

| Parameter | Type | Description |
| --- | --- | --- |
| `bbox` | `lsst.geom.Box2I` | Box to convert. |
| `parent_bbox` | `lsst.geom.Box2I` | Box of the array the slices index (its origin is pixel [0, 0]). |

| Returns | Type | Description |
| --- | --- | --- |
| `slices` | `tuple` [`slice`, `slice`] | (y slice, x slice) selecting ``bbox`` from an array covering ``parent_bbox``. |

### fitLib.collapse_saturated_peaks

`collapse_saturated_peaks(peak_positions, exposure, sat_bitmask=None, dilate_iters=2)`

Collapse every peak on or near a saturated blob into one peak per blob.

The representative position is a flux-weighted centroid from the unsaturated ring just outside the saturated core (the flat top carries no usable position or flux information). Without this, noise on the saturated flat top produces many near-duplicate peaks whose PSF stamps are almost identical; the joint fit then splits the true amplitude between them arbitrarily, which shows up as a blocky, high-frequency patch in the model image.

| Parameter | Type | Description |
| --- | --- | --- |
| `peak_positions` | `list` [`tuple` [`float`, `float`]] | Peak (x, y) positions in parent pixel coordinates. |
| `exposure` | `lsst.afw.image.Exposure` | Exposure the peaks were detected on. |
| `sat_bitmask` | `int`, optional | Mask bits of saturated pixels (default: the SAT plane). |
| `dilate_iters` | `int`, optional | Number of dilations of the saturated region when assigning peaks to a blob. |

| Returns | Type | Description |
| --- | --- | --- |
| `peaks` | `list` [`tuple` [`float`, `float`]] | Peaks with every saturated blob reduced to one representative peak. |
| `n_removed` | `int` | Number of peaks removed. |

### fitLib.merge_close_peaks

`merge_close_peaks(peak_positions, exposure, min_separation=None)`

Merge peaks closer than ``min_separation``, keeping the brightest pixel's peak.

Prevents near-duplicate columns in the design matrix (mostly from re-detecting the same wings) from destabilising the joint solve.

| Parameter | Type | Description |
| --- | --- | --- |
| `peak_positions` | `list` [`tuple` [`float`, `float`]] | Peak (x, y) positions in parent pixel coordinates. |
| `exposure` | `lsst.afw.image.Exposure` | Exposure whose pixel values rank the peaks. |
| `min_separation` | `float`, optional | Merge radius [px] (default: twice the PSF sigma). |

| Returns | Type | Description |
| --- | --- | --- |
| `merged_positions` | `list` [`tuple` [`float`, `float`]] | Surviving peak positions. |

### fitLib.joint_psf_fit

`joint_psf_fit(exposure, psf, peak_positions, bad_bitmask, amp_zero_tol=1e-06)`

Fit all PSF amplitudes simultaneously with a non-negative least-squares solve.

One sparse, weighted, bound-constrained solve against the original exposure.  Unlike an unconstrained solve followed by clipping, ``lsq_linear`` with bounds (0, inf) lets the other amplitudes re-adjust when one source is pushed to zero, so overlapping sources get a real non-negative amplitude instead of one of them being dropped from the model.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Image to fit (image, variance and mask planes are used). |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `peak_positions` | `list` [`tuple` [`float`, `float`]] | Source (x, y) positions in parent pixel coordinates. |
| `bad_bitmask` | `int` | Mask bits of pixels to ignore. |
| `amp_zero_tol` | `float`, optional | Amplitudes at or below this count as zero. |

| Returns | Type | Description |
| --- | --- | --- |
| `amps` | `numpy.ndarray` | Fitted amplitudes (fluxes), one per position. |
| `model_full` | `numpy.ndarray` | Model image over the exposure's box. |
| `n_at_zero_bound` | `int` | Number of sources that converged to (numerically) zero amplitude -- sources that still fail to fit rather than being absent. |

### fitLib.iterative_joint_detection

`iterative_joint_detection(exposure, psf, bad_bitmask, n_iterations=10, threshold_start=8.0, threshold_step=-1.5, threshold_min=1.5, min_pixels=1, min_new_peaks=2, convergence_frac=0.02, verbose=True)`

Detect, fit and subtract iteratively until the residual stops yielding new sources.

Each pass: (1) detect on the current residual with a decreasing threshold, so fainter sources are searched for once brighter ones are gone; (2) collapse multiple peaks on saturated cores; (3) merge the new peaks with all peaks found so far; (4) refit all of them jointly against the original exposure with a non-negative solve; (5) rebuild the residual. Every pass is a full re-solve, so errors from an early pass are not compounded.  The number of zero-amplitude sources is printed each pass: if it is large or not shrinking, real sources are being fitted down to nothing.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Image to fit. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `bad_bitmask` | `int` | Mask bits of pixels to ignore. |
| `n_iterations` | `int`, optional | Maximum number of passes. |
| `threshold_start` | `float`, optional | Detection threshold [sigma] of the first pass. |
| `threshold_step` | `float`, optional | Change of the threshold per pass (negative: decreasing). |
| `threshold_min` | `float`, optional | Lowest detection threshold [sigma]. |
| `min_pixels` | `int`, optional | Minimum footprint size of a detection. |
| `min_new_peaks` | `int`, optional | Stop when a pass finds fewer new peaks than this. |
| `convergence_frac` | `float`, optional | Stop when the residual standard deviation improves by less than this fraction. |
| `verbose` | `bool`, optional | Print progress. |

| Returns | Type | Description |
| --- | --- | --- |
| `result` | `dict` | ``model_full``, ``residual_full``, ``peaks``, ``amps``, the per-pass histories (``resid_std_history``, ``n_peaks_history``, ``n_zero_amp_history``) and ``first_seen`` (pass index at which each peak, keyed by its position rounded to 0.1 px, was found). |

### fitLib.plot_tiled_comparison

`plot_tiled_comparison(data_full, model_full, residual_full, n_sources=None, n_passes=None, tile_size=100, cmap_img=None, global_scale=True, skip_empty=True, empty_std_thresh=1e-06, figsize_per_panel=5, show_row_profile=True, row_select='max_flux', stretch='log')`

Plot data, model and residual (and a row profile) for every tile of an image.

The data and model panels share one colour scale when ``global_scale`` is True; the residual panel always scales to its own tile.  With ``show_row_profile`` a fourth panel shows data and model along one row of the tile (marked with a dashed cyan line), so the model can be compared with the data through a crowded section.

| Parameter | Type | Description |
| --- | --- | --- |
| `data_full, model_full, residual_full` | `numpy.ndarray` | Data, model and residual images (same shape). |
| `n_sources` | `int`, optional | Number of sources, for the title. |
| `n_passes` | `int`, optional | Number of detection passes, for the title. |
| `tile_size` | `int`, optional | Tile size [px]. |
| `cmap_img` | `str`, optional | Colour map of the data and model panels (default: 'mako' if available, else 'viridis'). |
| `global_scale` | `bool`, optional | One colour scale for all data / model panels instead of one per tile. |
| `skip_empty` | `bool`, optional | Skip tiles whose data standard deviation is below ``empty_std_thresh``. |
| `empty_std_thresh` | `float`, optional | See ``skip_empty``. |
| `figsize_per_panel` | `float`, optional | Size of one panel [inches]. |
| `show_row_profile` | `bool`, optional | Add the row-profile panel. |
| `row_select` | `str`, `int` or callable, optional | Row of each tile to profile: 'max_flux' (the row with the most flux), a row index local to the tile, or a function ``row_select(tile_data) -> row index``. |
| `stretch` | `str`, optional | Stretch of the data / model colour scale. |

### fitLib.plot_detections_by_iteration

`plot_detections_by_iteration(exposure_crop, result_iter)`

Plot every fitted peak over the image, coded by the pass that first found it.

Peaks whose final amplitude is (numerically) zero -- they hit the non-negativity bound -- are marked as failed fits with a red 'x'; they are the ones to inspect if sources seem to be missing from the model.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure_crop` | `lsst.afw.image.Exposure` | Image the pipeline ran on. |
| `result_iter` | `dict` | Output of `iterative_joint_detection`. |

| Returns | Type | Description |
| --- | --- | --- |
| `sources_tbl` | `astropy.table.Table` | Positions, amplitudes and detection pass of the sources. |

### fitLib.mask_sources_from_model

`mask_sources_from_model(exposure, model_full, noise_sigma=None, n_sigma=4.0, dilate_iters=1, verbose=True)`

Flag pixels where the joint-fit model has significant flux as DETECTED.

The next background fit then treats them as source, not sky.  The joint model contains all sources found so far, including ones too faint for the first background pass.  A loose threshold flags most of the image and starves the background fit; if the printed masked fraction climbs much above ~50 %, tighten ``n_sigma`` or ``dilate_iters``.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Exposure whose mask is modified in place. |
| `model_full` | `numpy.ndarray` | Model image over the exposure's box. |
| `noise_sigma` | `float`, optional | Pixel noise (default: sqrt of the median variance). |
| `n_sigma` | `float`, optional | Model pixels above ``n_sigma * noise_sigma`` are flagged. |
| `dilate_iters` | `int`, optional | Number of dilations of the flagged region. |
| `verbose` | `bool`, optional | Print the masked fraction. |

| Returns | Type | Description |
| --- | --- | --- |
| `src_mask` | `numpy.ndarray` [`bool`] | The flagged pixels. |

### fitLib.refit_background

`refit_background(exposure, algorithm='AKIMA_SPLINE', bin_size=32, use_approx=True, stats_mask=('BAD', 'SAT', 'EDGE', 'NO_DATA', 'DETECTED', 'CR', 'INTRP'))`

Re-derive the background ignoring masked pixels and subtract it in place.

Pass a fresh clone of the raw (never background-subtracted) exposure, already flagged by `mask_sources_from_model`, or the background is subtracted twice.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Exposure to background-subtract in place. |
| `algorithm` | `str`, optional | Interpolation algorithm (only used when ``use_approx`` is False). |
| `bin_size` | `int`, optional | Background bin size [px]. |
| `use_approx` | `bool`, optional | Fit a smooth Chebyshev surface to the binned statistics instead of interpolating; more robust when many pixels are masked. |
| `stats_mask` | `tuple` [`str`], optional | Mask planes ignored when computing the background statistics. |

| Returns | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The background-subtracted exposure (same object). |
| `background` | `lsst.afw.image.Image` | The subtracted background image. |

### fitLib.catalog_from_peaks

`catalog_from_peaks(exposure, peak_positions, schema, measureTask, stamp_radius=15, exclude_positions=None, exclude_radius=5.0, verbose=True)`

Build a measured SourceCatalog from raw (x, y) peak positions.

Each peak gets a circular footprint and one Peak record, and the catalogue is measured with the same ``measureTask`` as fresh detections, so the fields are comparable and usable by the star selector.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Exposure to measure on. |
| `peak_positions` | `list` [`tuple` [`float`, `float`]] | Peak positions in parent pixel coordinates. |
| `schema` | `lsst.afw.table.Schema` | The schema used to build ``measureTask``. |
| `measureTask` | `lsst.meas.base.SingleFrameMeasurementTask` | Measurement task. |
| `stamp_radius` | `int`, optional | Footprint radius [px]. |
| `exclude_positions` | `list` [`tuple` [`float`, `float`]], optional | Skip peaks within ``exclude_radius`` of these positions (e.g. existing detections). |
| `exclude_radius` | `float`, optional | See ``exclude_positions`` [px]. |
| `verbose` | `bool`, optional | Print how many peaks were skipped. |

| Returns | Type | Description |
| --- | --- | --- |
| `catalog` | `lsst.afw.table.SourceCatalog` | The measured catalogue. |

### fitLib.redetect_and_refit_psf

`redetect_and_refit_psf(exposure, det_threshold=3.0, det_min_pixels=2, det_local_bg_bin=16, det_bg_bin=32, star_selector_config=None, psf_determiner_config=None, kernel_size=41, extra_peak_positions=None, force_include_star_positions=None, peak_stamp_radius=15, verbose=True)`

Detect, deblend, measure, select stars and determine a new PSF.

Mirrors the original one-shot pipeline, so each outer pass works on the newly cleaned background.  Extra peak positions (e.g. the previous pass's converged catalogue) are added before star selection; forced star positions are appended after selection, bypassing its cuts.  Candidates with failed stamps are filtered out before the PSF determiner runs (even one can make it reject the whole set).  If too few candidates survive, the exposure's current PSF is kept and a "NOT updated" note is printed.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Background-subtracted exposure. |
| `det_threshold` | `float`, optional | Detection threshold [sigma]. |
| `det_min_pixels` | `int`, optional | Minimum footprint size. |
| `det_local_bg_bin` | `int`, optional | Bin size of the temporary local background [px]. |
| `det_bg_bin` | `int`, optional | Bin size of the detection background [px]. |
| `star_selector_config` | `lsst.meas.algorithms.ObjectSizeStarSelectorConfig`, optional | Star selector configuration. |
| `psf_determiner_config` | `lsst.meas.algorithms.PcaPsfDeterminerConfig`, optional | PSF determiner configuration. |
| `kernel_size` | `int`, optional | PSF kernel size [px], used for both the candidates and the determiner. |
| `extra_peak_positions` | `list` [`tuple` [`float`, `float`]], optional | Additional peaks measured and offered to the star selector. |
| `force_include_star_positions` | `list` [`tuple` [`float`, `float`]], optional | Star positions appended to the selected stars regardless of the selector's cuts. |
| `peak_stamp_radius` | `int`, optional | Footprint radius of the extra / forced peaks [px]. |
| `verbose` | `bool`, optional | Print progress. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf_new` | `lsst.afw.detection.Psf` | The new PSF (or the current one if it could not be determined). |
| `sources` | `lsst.afw.table.SourceCatalog` | Measured detections. |
| `starCat` | `lsst.afw.table.SourceCatalog` | Selected PSF stars. |
| `cellSet` | `lsst.afw.math.SpatialCellSet` | PSF candidates used by the determiner. |

### fitLib.outer_loop_psf_background_refinement

`outer_loop_psf_background_refinement(raw_exposure_crop, psf_initial, bad_bitmask, initial_star_cat=None, n_outer_iterations=3, n_inner_iterations=10, inner_kwargs=None, outer_convergence_frac=0.02, bg_use_approx=True, psf_kernel_size=41, verbose=True)`

Alternate inner-loop fits with background and PSF re-derivation.

Pass 0 runs `iterative_joint_detection` with the current background and PSF.  Pass k > 0 masks every known source (from pass k-1's model) on a fresh clone of the raw image, re-derives the background, re-determines the PSF (offering pass k-1's catalogue and forcing in ``initial_star_cat``) and reruns the inner loop.

| Parameter | Type | Description |
| --- | --- | --- |
| `raw_exposure_crop` | `lsst.afw.image.Exposure` | Crop before any background subtraction. |
| `psf_initial` | `lsst.afw.detection.Psf` | PSF of the first pass. |
| `bad_bitmask` | `int` | Mask bits of pixels to ignore. |
| `initial_star_cat` | `lsst.afw.table.SourceCatalog`, optional | Star selector output whose positions are forced into every PSF determination. |
| `n_outer_iterations` | `int`, optional | Maximum number of outer passes. |
| `n_inner_iterations` | `int`, optional | Maximum number of inner passes (when ``inner_kwargs`` is not given). |
| `inner_kwargs` | `dict`, optional | Keyword arguments of `iterative_joint_detection`. |
| `outer_convergence_frac` | `float`, optional | Stop when the residual standard deviation improves by less than this fraction. |
| `bg_use_approx` | `bool`, optional | See `refit_background`. |
| `psf_kernel_size` | `int`, optional | PSF kernel size [px]. |
| `verbose` | `bool`, optional | Print progress. |

| Returns | Type | Description |
| --- | --- | --- |
| `outer_history` | `list` [`dict`] | One dict per outer pass: ``psf``, ``exposure``, ``result`` (inner-loop output), ``resid_std``, ``n_peaks``. |

### fitLib.plot_psf_models_by_outer_pass

`plot_psf_models_by_outer_pass(outer_history, exposure_center=None)`

Plot the PSF kernel image used in each outer pass.

| Parameter | Type | Description |
| --- | --- | --- |
| `outer_history` | `list` [`dict`] | Outer-loop history (each entry has ``psf`` and ``exposure``). |
| `exposure_center` | `lsst.geom.Point2D`, optional | Position at which to render the kernel (default: the exposure centre). |

### fitLib.plot_background_by_outer_pass

`plot_background_by_outer_pass(outer_history)`

Plot the background re-derived going into each outer pass.

The first pass uses the background subtracted before the pipeline ran and is skipped.

| Parameter | Type | Description |
| --- | --- | --- |
| `outer_history` | `list` [`dict`] | Outer-loop history (entries may have ``bg_image``). |

### fitLib.plot_sources_overlaid_by_pass

`plot_sources_overlaid_by_pass(outer_history, marker_size=25, vmin_pct=5, vmax_pct=99)`

Plot the data of each outer pass with its converged peak positions.

| Parameter | Type | Description |
| --- | --- | --- |
| `outer_history` | `list` [`dict`] | Outer-loop history. |
| `marker_size` | `float`, optional | Marker size. |
| `vmin_pct, vmax_pct` | `float`, optional | Percentiles of the display range. |

### fitLib.plot_outer_loop_summary

`plot_outer_loop_summary(outer_history)`

Plot the residual standard deviation and the source count against outer pass.

| Parameter | Type | Description |
| --- | --- | --- |
| `outer_history` | `list` [`dict`] | Outer-loop history (entries have ``resid_std`` and ``n_peaks``). |

### fitLib.plot_psf_fwhm_by_outer_pass

`plot_psf_fwhm_by_outer_pass(outer_history, exposure_center=None)`

Plot the FWHM of the PSF used in each outer pass.

| Parameter | Type | Description |
| --- | --- | --- |
| `outer_history` | `list` [`dict`] | Outer-loop history. |
| `exposure_center` | `lsst.geom.Point2D`, optional | Position at which to evaluate the PSF (default: the exposure centre). |

| Returns | Type | Description |
| --- | --- | --- |
| `fwhms` | `list` [`float`] | FWHM [px] per outer pass. |

### fitLib.StepTimer

Class.  `StepTimer()`

Collect timing records of pipeline steps.

Each record is a dict with ``label``, ``seconds``, ``outer``, ``inner`` and ``kind`` ('compute' or 'viz').

#### fitLib.StepTimer.step

`step(self, label, outer=None, inner=None, kind='compute', verbose=False)`

Time a block of code (context manager).

| Parameter | Type | Description |
| --- | --- | --- |
| `label` | `str` | Step name. |
| `outer` | `int`, optional | Outer-pass index. |
| `inner` | `int`, optional | Inner-iteration index. |
| `kind` | `str`, optional | 'compute' or 'viz'. |
| `verbose` | `bool`, optional | Print the elapsed time. |

| Yields | Type | Description |
| --- | --- | --- |
| `rec` | `dict` | The record; ``seconds`` is filled in when the block exits. |

#### fitLib.StepTimer.add

`add(self, label, seconds, outer=None, inner=None, kind='compute')`

Add a record with a known duration.

| Parameter | Type | Description |
| --- | --- | --- |
| `label` | `str` | Step name. |
| `seconds` | `float` | Duration [s]. |
| `outer` | `int`, optional | Outer-pass index. |
| `inner` | `int`, optional | Inner-iteration index. |
| `kind` | `str`, optional | 'compute' or 'viz'. |

#### fitLib.StepTimer.reset

`reset(self)`

Remove all records and restart the clock.

### fitLib.plot_timing_summary

`plot_timing_summary(timer=TIMER)`

Plot the total time per step and the time per inner iteration.

| Parameter | Type | Description |
| --- | --- | --- |
| `timer` | `StepTimer`, optional | Timer to summarise (default: the module-level ``TIMER``). |

### fitLib.make_cutout

`make_cutout(exposure_crop, x_range=(0, 100), y_range=(0, 100))`

Cut a sub-exposure out of a crop, in coordinates relative to the crop's corner.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure_crop` | `lsst.afw.image.Exposure` | The crop. |
| `x_range, y_range` | `tuple` [`int`, `int`], optional | Pixel ranges relative to the crop's origin (end exclusive). |

| Returns | Type | Description |
| --- | --- | --- |
| `cutout` | `lsst.afw.image.Exposure` | A deep copy of the sub-exposure. |

Raises: `TypeError` (Raised if ``exposure_crop`` is not an Exposure.)

### fitLib.psf_stamp

`psf_stamp(psf, x, y, clip_bbox, taper=None)`

Render the PSF at (x, y), clipped to a box, with its far wings tamed.

With ``taper`` the median of the outermost pixel ring (edge pedestal) is subtracted and the stamp is tapered radially to zero, so each source ends smoothly in a disc instead of a hard-edged square.

| Parameter | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `x, y` | `float` | Position in parent pixel coordinates. |
| `clip_bbox` | `lsst.geom.Box2I` | Box to clip the stamp to. |
| `taper` | `tuple` [`float`, `float`], optional | Inner and outer taper radii as fractions of the stamp half-size. |

| Returns | Type | Description |
| --- | --- | --- |
| `stamp` | `numpy.ndarray` or `None` | The clipped stamp, or None if it misses ``clip_bbox``. |
| `bbox` | `lsst.geom.Box2I` | Box of the clipped stamp. |

### fitLib.drop_masked_peaks

`drop_masked_peaks(peaks, exposure, planes=('EDGE', 'NO_DATA'))`

Drop peaks on EDGE / NO_DATA pixels.

SAT is deliberately not in the default list, so saturated stars are kept.

| Parameter | Type | Description |
| --- | --- | --- |
| `peaks` | `list` [`tuple` [`float`, `float`]] | Peak positions in parent pixel coordinates. |
| `exposure` | `lsst.afw.image.Exposure` | Exposure providing the mask. |
| `planes` | `tuple` [`str`], optional | Mask planes that reject a peak. |

| Returns | Type | Description |
| --- | --- | --- |
| `keep` | `list` [`tuple` [`float`, `float`]] | Remaining peaks. |
| `n_dropped` | `int` | Number of peaks dropped. |

### fitLib.drop_out_of_bounds

`drop_out_of_bounds(peaks, exposure, margin=0.0)`

Drop peaks outside the exposure's box grown by ``margin`` px.

| Parameter | Type | Description |
| --- | --- | --- |
| `peaks` | `list` [`tuple` [`float`, `float`]] | Peak positions in parent pixel coordinates. |
| `exposure` | `lsst.afw.image.Exposure` | Exposure providing the box. |
| `margin` | `float`, optional | Growth of the box [px]; a negative margin drops peaks within ``\|margin\|`` px of the edge. |

| Returns | Type | Description |
| --- | --- | --- |
| `keep` | `list` [`tuple` [`float`, `float`]] | Remaining peaks. |
| `n_dropped` | `int` | Number of peaks dropped. |

### fitLib.effective_variance

`effective_variance(exposure, model=None, frac=0.03, poisson_k=0.0)`

Return the variance used to weight the fit.

``variance + (frac * model)**2 + poisson_k * model``: the pixel variance plus a fractional model-error term and missing photon noise.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Exposure providing the variance plane. |
| `model` | `numpy.ndarray`, optional | Current model image (no extra terms if None). |
| `frac` | `float`, optional | Fractional model error. |
| `poisson_k` | `float`, optional | Extra variance per unit model flux. |

| Returns | Type | Description |
| --- | --- | --- |
| `var` | `numpy.ndarray` | The effective variance. |

### fitLib.reject_near_bright

`reject_near_bright(new_peaks, old_peaks, old_is_bright, min_sep)`

Drop new peaks within ``min_sep`` px of an existing bright source.

| Parameter | Type | Description |
| --- | --- | --- |
| `new_peaks` | `list` [`tuple` [`float`, `float`]] | Candidate peaks. |
| `old_peaks` | array-like | Existing source positions, shape (N, 2). |
| `old_is_bright` | array-like [`bool`] | Which existing sources are bright. |
| `min_sep` | `float` | Rejection radius [px]. |

| Returns | Type | Description |
| --- | --- | --- |
| `keep` | `list` [`tuple` [`float`, `float`]] | Remaining candidates. |
| `n_rejected` | `int` | Number of candidates rejected. |

### fitLib.prune_and_merge

`prune_and_merge(fit, peaks, min_snr=5.0, min_sep=3.0, bright_flags=None, near_bright_sep=None, near_bright_snr=8.0)`

Drop low-S/N sources, then merge pairs closer than ``min_sep`` px.

The brighter source of a pair survives at the amplitude-weighted position.  Sources within ``near_bright_sep`` px of a bright source must reach ``near_bright_snr`` instead of ``min_snr``.  Uses the cheap diagonal-only S/N; the honest S/N is applied by `eliminate_weak`.

| Parameter | Type | Description |
| --- | --- | --- |
| `fit` | `dict` | Joint fit (from `joint_psf_fit_diag`) of ``peaks``. |
| `peaks` | array-like | Source positions, shape (N, 2). |
| `min_snr` | `float`, optional | Minimum S/N. |
| `min_sep` | `float`, optional | Merge radius [px]. |
| `bright_flags` | array-like [`bool`], optional | Which sources are bright. |
| `near_bright_sep` | `float`, optional | Radius around bright sources with the stricter S/N cut [px]. |
| `near_bright_snr` | `float`, optional | S/N required near bright sources. |

| Returns | Type | Description |
| --- | --- | --- |
| `peaks` | `list` [`tuple` [`float`, `float`]] | Surviving source positions. |

### fitLib.empirical_noise_scale

`empirical_noise_scale(fit, exposure, model_frac=0.01, matched=True, weight=None)`

Return the empirical noise scale of a fit's residual.

| Parameter | Type | Description |
| --- | --- | --- |
| `fit` | `dict` | Joint fit (from `joint_psf_fit_diag`). |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `model_frac` | `float`, optional | Pixels with model below this fraction of the model maximum count as sky. |
| `matched` | `bool`, optional | True: scale for S/N values -- the pixel scale moved ``weight`` of the way towards the matched-filter scale when that is larger.  False: the pixel scale only (for chi2). |
| `weight` | `float`, optional | See ``matched`` (default ``MATCHED_NOISE_WEIGHT``). |

| Returns | Type | Description |
| --- | --- | --- |
| `scale` | `float` | Noise scale (at least 1). |

### fitLib.amp_snr_full

`amp_snr_full(fit, noise_scale=1.0)`

Return honest amplitude S/N values, including the covariance of overlapping sources.

``S/N = amp / sqrt(diag(N^-1))`` with ``N = A^T W A`` over all sources, divided by the empirical noise scale.

| Parameter | Type | Description |
| --- | --- | --- |
| `fit` | `dict` | Joint fit (from `joint_psf_fit_diag`). |
| `noise_scale` | `float`, optional | Empirical noise scale dividing the S/N. |

| Returns | Type | Description |
| --- | --- | --- |
| `snr` | `numpy.ndarray` | S/N per source. |

### fitLib.eliminate_weak

`eliminate_weak(exposure, psf, peaks, fit, bad_bitmask, min_snr=5.0, drop_frac=0.25, max_rounds=30, **fitkw)`

Remove sources below ``min_snr`` by backward elimination.

Each round drops the weakest ~``drop_frac`` of the sub-threshold sources (at least one) and re-solves, until every source has an honest S/N of at least ``min_snr``.  The last source is never removed.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `peaks` | `list` [`tuple` [`float`, `float`]] | Source positions; ``fit`` must correspond to them. |
| `fit` | `dict` | Joint fit of ``peaks``. |
| `bad_bitmask` | `int` | Mask bits of pixels to ignore. |
| `min_snr` | `float`, optional | Minimum honest S/N. |
| `drop_frac` | `float`, optional | Fraction of the failing sources dropped per round. |
| `max_rounds` | `int`, optional | Maximum number of rounds. |
| `**fitkw` | – | Keyword arguments of `joint_psf_fit_diag`. |

| Returns | Type | Description |
| --- | --- | --- |
| `peaks` | `list` [`tuple` [`float`, `float`]] | Surviving source positions. |
| `fit` | `dict` | Joint fit of the survivors. |

### fitLib.shape_metrics

`shape_metrics(exposure, psf, peaks, fit, noise_scale=1.0, chi2_scale=None, bg_sigma=4.0, width_scales=(0.6, 0.75, 1.0, 1.3, 1.7), min_pix=15)`

Measure how PSF-like each source is.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `peaks` | `list` [`tuple` [`float`, `float`]] | Source positions; ``fit`` must correspond to them (same order). |
| `fit` | `dict` | Joint fit of ``peaks``. |
| `noise_scale` | `float`, optional | Divides the S/N values (matched-filter calibrated). |
| `chi2_scale` | `float`, optional | Divides chi2 (pixel-level calibrated; default ``noise_scale``). |
| `bg_sigma` | `float`, optional | Fit radius in units of the PSF sigma. |
| `width_scales` | `tuple` [`float`], optional | Trial PSF width scale factors. |
| `min_pix` | `int`, optional | Minimum number of good pixels to judge a source. |

| Returns | Type | Description |
| --- | --- | --- |
| `metrics` | `dict` [`str`, `numpy.ndarray`] | ``snr`` (local, noise-scaled), ``chi2red``, ``best_scale``, ``dchi2_width`` (> 0: the data prefer ``best_scale`` over the PSF width) and ``contrast``; NaN = not judged. |

### fitLib.judge_shape

`judge_shape(m, min_snr_local=5.0, chi2_max=2.0, width_range=(0.75, 1.35), dchi2_width=9.0, trust_snr=30.0, min_contrast=20.0, weak_snr=10.0, dchi2_width_weak=4.0)`

Decide which sources fail the PSF-shape tests.

Sources with local S/N of at least ``trust_snr`` (including saturated ones) skip the chi2, width and contrast tests.  Untrusted sources with S/N below ``weak_snr`` use the looser width threshold ``dchi2_width_weak``.

| Parameter | Type | Description |
| --- | --- | --- |
| `m` | `dict` [`str`, `numpy.ndarray`] | Output of `shape_metrics`. |
| `min_snr_local` | `float`, optional | Minimum local S/N. |
| `chi2_max` | `float`, optional | Maximum reduced chi2. |
| `width_range` | `tuple` [`float`, `float`], optional | Allowed range of the best-fitting width scale. |
| `dchi2_width` | `float`, optional | Improvement in chi2 needed to reject on width. |
| `trust_snr` | `float`, optional | S/N above which a source is trusted. |
| `min_contrast` | `float`, optional | Minimum ratio of the source's peak to the neighbours' model flux at its position (None or 0 = off); 20 ~ (5 % wing error) x (4x margin). |
| `weak_snr` | `float`, optional | S/N below which ``dchi2_width_weak`` applies. |
| `dchi2_width_weak` | `float`, optional | Width threshold for weak sources. |

| Returns | Type | Description |
| --- | --- | --- |
| `fail` | `numpy.ndarray` [`bool`] | Sources failing any test. |
| `score` | `numpy.ndarray` | Badness score (higher = worse); orders the elimination. |
| `flags` | `dict` [`str`, `numpy.ndarray`] | Per-test failures: ``snr``, ``chi2``, ``width``, ``contrast``. |

### fitLib.eliminate_nonpsf

`eliminate_nonpsf(exposure, psf, peaks, fit, bad_bitmask, shape_kwargs=None, drop_frac=0.25, max_rounds=15, verbose=True, **fitkw)`

Remove sources that fail the PSF-shape test by backward elimination.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `peaks` | `list` [`tuple` [`float`, `float`]] | Source positions; ``fit`` must correspond to them. |
| `fit` | `dict` | Joint fit of ``peaks``. |
| `bad_bitmask` | `int` | Mask bits of pixels to ignore. |
| `shape_kwargs` | `dict`, optional | Keyword arguments of `shape_metrics` and `judge_shape`. |
| `drop_frac` | `float`, optional | Fraction of the failing sources dropped per round. |
| `max_rounds` | `int`, optional | Maximum number of rounds. |
| `verbose` | `bool`, optional | Print the rejected sources. |
| `**fitkw` | – | Keyword arguments of `joint_psf_fit_diag`. |

| Returns | Type | Description |
| --- | --- | --- |
| `peaks` | `list` [`tuple` [`float`, `float`]] | Surviving source positions. |
| `fit` | `dict` | Joint fit of the survivors. |
| `rejected` | `list` [`tuple` [`float`, `float`, `float`]] | (x, y, score) of the removed sources. |

### fitLib.plot_shape_rejects

`plot_shape_rejects(exposure, rejected_by_iter, peaks=None, title='sources removed by shape cut')`

Plot the sources removed by the shape cut over the data.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `rejected_by_iter` | `list` [`list` [`tuple`]] | Per iteration, the (x, y, score) of rejected sources (magenta x). |
| `peaks` | array-like, optional | Kept source positions (cyan). |
| `title` | `str`, optional | Plot title. |

### fitLib.select_source_pages

`select_source_pages(snap, source_view='all', sources_per_page=6, max_pages=None)`

Choose the sources that get per-source pages.

| Parameter | Type | Description |
| --- | --- | --- |
| `snap` | `dict` | Inner-loop snapshot (``fit``, ``prev_amp``, ``is_new``). |
| `source_view` | `str`, optional | 'all', 'new' (added this iteration) or 'changed' (new, or amplitude changed by more than 10 %). |
| `sources_per_page` | `int`, optional | Sources per page. |
| `max_pages` | `int`, optional | Maximum number of pages (None = no limit). |

| Returns | Type | Description |
| --- | --- | --- |
| `pages` | `list` [`list` [`int`]] | Source indices per page, brightest first. |

### fitLib.source_colors

`source_colors(pages)`

Assign each source on the pages a colour from ``SRC_PALETTE``.

| Parameter | Type | Description |
| --- | --- | --- |
| `pages` | `list` [`list` [`int`]] | Output of `select_source_pages`. |

| Returns | Type | Description |
| --- | --- | --- |
| `colors` | `dict` [`int`, `str`] | Colour per source index. |

### fitLib.joint_psf_fit_diag

`joint_psf_fit_diag(exposure, psf, peak_positions, bad_bitmask, amp_zero_tol=1e-06, support_frac=0.01, var_model=None, sys_frac=0.03, poisson_k=0.0)`

Fit all PSF amplitudes in one non-negative solve, keeping what is needed for diagnostics.

Weighted, non-negative, one sparse solve against the original exposure.  PSF stamps are not clipped at zero (clipping would add a positive pedestal).

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Image to fit. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `peak_positions` | `list` [`tuple` [`float`, `float`]] | Source positions in parent pixel coordinates. |
| `bad_bitmask` | `int` | Mask bits of pixels to ignore. |
| `amp_zero_tol` | `float`, optional | Amplitudes at or below this count as zero. |
| `support_frac` | `float`, optional | A source's significant support is where its stamp exceeds this fraction of its peak. |
| `var_model` | `numpy.ndarray`, optional | Model image used for the model-dependent variance terms (see `effective_variance`). |
| `sys_frac` | `float`, optional | Fractional model error added to the variance. |
| `poisson_k` | `float`, optional | Extra variance per unit model flux. |

| Returns | Type | Description |
| --- | --- | --- |
| `fit` | `dict` | ``amps``, ``model_full``, ``n_at_zero``; ``comps`` (per source a dict with ``stamp``, ``sig`` (significant support), ``ys``/``xs`` (slices of its fit region), ``amp``, ``bbox_parent``, or None if the stamp misses the image); coupled-source groups (``group_labels``, ``n_groups``, ``group_map``, ``group_chi2``); ``good`` (pixels used), ``chi2map``, ``variance`` (the effective variance), ``timings`` and the solver status (``lsq_status``, ``lsq_nit``, ``lsq_cost``). |

### fitLib.compute_psf_derivatives

`compute_psf_derivatives(psf, x, y, bbox_clip, eps=0.25)`

Compute central-difference derivative stamps dP/dx and dP/dy.

``computeImage`` centres its box on the pixel nearest the requested point, so the stamps at (x +- eps, y) and (x, y +- eps) need not share a box; all five boxes (and the clip box) are intersected first.

| Parameter | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `x, y` | `float` | Position in parent pixel coordinates. |
| `bbox_clip` | `lsst.geom.Box2I` | Box to clip to. |
| `eps` | `float`, optional | Finite-difference step [px]. |

| Returns | Type | Description |
| --- | --- | --- |
| `dpdx, dpdy` | `numpy.ndarray` or `None` | Derivative stamps over ``bbox``, or None if the boxes do not overlap. |
| `bbox` | `lsst.geom.Box2I` | The common box. |

### fitLib.windowed_centroid

`windowed_centroid(image_array, xy0, x, y, half_width, max_frac_of_window=0.9)`

Return the flux-weighted centroid in a box around (x, y).

Negative pixels are clipped to zero.  Falls back to (x, y) if there is too little positive flux or the centroid lands too close to the window edge (dominated by a neighbour).

| Parameter | Type | Description |
| --- | --- | --- |
| `image_array` | `numpy.ndarray` | Image. |
| `xy0` | `tuple` [`int`, `int`] | Parent pixel coordinates of ``image_array[0, 0]``. |
| `x, y` | `float` | Starting position in parent pixel coordinates. |
| `half_width` | `int` | Half-size of the box [px]. |
| `max_frac_of_window` | `float`, optional | Reject centroids further than this fraction of ``half_width`` from the start. |

| Returns | Type | Description |
| --- | --- | --- |
| `cx, cy` | `float` | Centroid in parent pixel coordinates. |

### fitLib.refine_peak_positions

`refine_peak_positions(exposure, psf, peaks, amps, bad_bitmask, min_amp_for_refine=50.0, max_shift=0.6, n_subiters=3, min_good_pixels=10, core_radius_sigma=2.5, variance=None, freeze_sep=None)`

Move each well-detected source towards the position that minimises its residual.

Linearised Gauss-Newton with all other sources held fixed: over the source's core pixels, weighted by 1/sigma, ``(data - model) ~ dA * P + a * dx * dP/dx + a * dy * dP/dy``.  The dA column absorbs amplitude errors so they do not leak into the position; only the position step is applied and the caller re-solves the amplitudes.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `peaks` | `list` [`tuple` [`float`, `float`]] | Source positions. |
| `amps` | array-like | Source amplitudes. |
| `bad_bitmask` | `int` | Mask bits of pixels to ignore. |
| `min_amp_for_refine` | `float`, optional | Fainter sources keep their position. |
| `max_shift` | `float`, optional | Maximum shift per sub-step [px]. |
| `n_subiters` | `int`, optional | Number of Gauss-Newton sub-steps. |
| `min_good_pixels` | `int`, optional | Minimum number of good core pixels. |
| `core_radius_sigma` | `float`, optional | Core radius in units of the PSF sigma. |
| `variance` | `numpy.ndarray`, optional | Effective variance (default: the exposure's variance plane). |
| `freeze_sep` | `float`, optional | Sources with a brighter neighbour within this many px keep their position (it is degenerate with the neighbour's amplitude). |

| Returns | Type | Description |
| --- | --- | --- |
| `peaks` | `list` [`tuple` [`float`, `float`]] | Refined positions. |

### fitLib.redetect_and_refit_psf_adaptive

`redetect_and_refit_psf_adaptive(exposure, kernel_sizes=(41, 27, 21, 15, 11), **kwargs)`

Run `redetect_and_refit_psf`, trying kernel sizes from largest to smallest.

MakePsfCandidatesTask drops every candidate whose stamp touches masked or edge pixels; on a small cutout a large kernel can drop all of them.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Background-subtracted exposure. |
| `kernel_sizes` | `tuple` [`int`], optional | Kernel sizes to try [px]. |
| `**kwargs` | – | Keyword arguments of `redetect_and_refit_psf`. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf_new` | `lsst.afw.detection.Psf` | The new PSF (the previous one if every attempt failed). |
| `sources` | `lsst.afw.table.SourceCatalog` or `None` | Measured detections. |
| `starCat` | `lsst.afw.table.SourceCatalog` or `None` | Selected PSF stars. |
| `cellSet` | `lsst.afw.math.SpatialCellSet` or `None` | PSF candidates. |
| `kernel_size` | `int` or `None` | Kernel size that worked (None if all failed). |

### fitLib.plot_iteration_overview

`plot_iteration_overview(exposure, snap, hist, outer_idx, norm_data, colors=None, view_half=20, circle_radius=None)`

Plot the overview of one inner iteration.

Detection input, detections, fitted sources, model, residual and the iteration histories.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `snap` | `dict` | Inner-loop snapshot of the iteration. |
| `hist` | `dict` | Histories (``resid_std``, ``n_peaks``, ``n_zero``). |
| `outer_idx` | `int` | Outer-pass index. |
| `norm_data` | `astropy.visualization.ImageNormalize` | Display normalisation of the data. |
| `colors` | `dict` [`int`, `str`], optional | Colours of the sources with per-source pages. |
| `view_half` | `int`, optional | Half-size of the per-source windows [px]. |
| `circle_radius` | `float`, optional | Marker radius [px]. |

### fitLib.run_inner_loop_verbose

`run_inner_loop_verbose(exposure, psf, bad_bitmask, outer_idx=0, n_iterations=10, threshold_start=8.0, threshold_step=-1.5, threshold_min=4.5, min_pixels=1, min_new_peaks=1, convergence_frac=0.02, viz=True, plot_every=1, source_view='all', sources_per_page=6, max_pages_per_iter=None, view_half=20, support_frac=0.01, refine_centroids=True, centroid_min_amp=50.0, n_centroid_rounds=4, centroid_shift_tol=0.03, centroid_shift_schedule=None, coarse_centroid=True, coarse_centroid_hw_sigma=2.5, merge_separation_sigma=1.5, sys_frac=0.03, min_snr=5.0, prune_sep_sigma=1.5, reject_fwhm=0.75, bright_peak_sigma=100.0, freeze_fwhm=1.0, poisson_k=0.0, near_bright_fwhm=2.0, near_bright_snr=10.0, edge_margin_fwhm=0.5, use_eliminate_weak=True, use_shape_cut=True, shape_kwargs=None, timer=TIMER, verbose=True, initial_peaks=None)`

Run the instrumented inner loop (detect, fit, refine, prune, repeat).

An instrumented version of `iterative_joint_detection` with centroid refinement and guards against over-counting sources.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Image to fit. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `bad_bitmask` | `int` | Mask bits of pixels to ignore. |
| `outer_idx` | `int`, optional | Outer-pass index (for labels). |
| `n_iterations` | `int`, optional | Maximum number of iterations. |
| `threshold_start, threshold_step, threshold_min` | `float`, optional | Detection threshold schedule [sigma]. |
| `min_pixels` | `int`, optional | Minimum footprint size. |
| `min_new_peaks` | `int`, optional | Stop when an iteration adds fewer peaks than this. |
| `convergence_frac` | `float`, optional | Stop when the residual standard deviation improves by less than this fraction. |
| `viz` | `bool`, optional | Draw diagnostics. |
| `plot_every` | `int`, optional | Draw every n-th iteration (the last always). |
| `source_view` | `str`, optional | See `select_source_pages`. |
| `sources_per_page` | `int`, optional | See `select_source_pages`. |
| `max_pages_per_iter` | `int`, optional | Maximum per-source pages per iteration. |
| `view_half` | `int`, optional | Half-size of the per-source windows [px]. |
| `support_frac` | `float`, optional | See `joint_psf_fit_diag`. |
| `refine_centroids` | `bool`, optional | Refine positions (`refine_peak_positions`). |
| `centroid_min_amp` | `float`, optional | Fainter sources keep their position. |
| `n_centroid_rounds` | `int`, optional | Maximum centroid rounds per iteration. |
| `centroid_shift_tol` | `float`, optional | Stop centroiding when the largest shift is below this [px]. |
| `centroid_shift_schedule` | `list` [`float`], optional | Maximum shift per round [px]. |
| `coarse_centroid` | `bool`, optional | Start new peaks from a windowed centroid. |
| `coarse_centroid_hw_sigma` | `float`, optional | Half-width of that window in PSF sigma. |
| `merge_separation_sigma` | `float`, optional | Merge new peaks closer than this many PSF sigma. |
| `sys_frac` | `float`, optional | Fractional PSF-shape error added in quadrature to the variance. |
| `min_snr` | `float`, optional | S/N floor (cheap version in `prune_and_merge`, honest one in `eliminate_weak`). |
| `prune_sep_sigma` | `float`, optional | Pairs closer than this many PSF sigma are merged after refinement. |
| `reject_fwhm` | `float`, optional | New peaks within this many FWHM of a bright source are discarded. |
| `bright_peak_sigma` | `float`, optional | A source is bright if its peak exceeds this many noise sigma. |
| `freeze_fwhm` | `float`, optional | Sources with a brighter neighbour within this many FWHM keep their position. |
| `poisson_k` | `float`, optional | Extra variance per unit model flux. |
| `near_bright_fwhm` | `float`, optional | Radius around bright stars with the stricter S/N cut [FWHM]. |
| `near_bright_snr` | `float`, optional | S/N required there. |
| `edge_margin_fwhm` | `float`, optional | Drop sources closer than this many FWHM to the edge. |
| `use_eliminate_weak` | `bool`, optional | Backward elimination on the honest S/N. |
| `use_shape_cut` | `bool`, optional | Backward elimination on PSF-likeness (`eliminate_nonpsf`). |
| `shape_kwargs` | `dict`, optional | Keyword arguments of `shape_metrics` / `judge_shape`. |
| `timer` | `StepTimer`, optional | Timer collecting the step durations. |
| `verbose` | `bool`, optional | Print progress. |
| `initial_peaks` | `list` [`tuple` [`float`, `float`]], optional | Warm-start source positions. |

| Returns | Type | Description |
| --- | --- | --- |
| `result` | `dict` | ``model_full``, ``residual_full``, ``peaks``, ``amps``, the histories (``resid_std_history``, ``n_peaks_history``, ``n_zero_amp_history``), ``first_seen`` and ``snaps`` (one dict per iteration with everything needed for plotting). |

### fitLib.plot_amplitude_evolution

`plot_amplitude_evolution(result, outer_idx=0)`

Plot the amplitude of every final source against inner iteration.

| Parameter | Type | Description |
| --- | --- | --- |
| `result` | `dict` | Inner-loop result with ``snaps``. |
| `outer_idx` | `int`, optional | Outer-pass index (for the title). |

### fitLib.plot_source_evolution

`plot_source_evolution(exposure, result, outer_idx=0, view_half=20)`

Plot model and residual cross-sections of three sources at every iteration.

The brightest, median and faintest final sources: the x cross-section of the model (top) and of the residual (bottom).

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `result` | `dict` | Inner-loop result with ``snaps``. |
| `outer_idx` | `int`, optional | Outer-pass index (for the title). |
| `view_half` | `int`, optional | Half-length of the cross-sections [px]. |

### fitLib.plot_cutout_setup

`plot_cutout_setup(cutout, psf, bad_bitmask, star_positions=None)`

Plot the cutout before fitting: data, sigma map, mask planes and PSF.

| Parameter | Type | Description |
| --- | --- | --- |
| `cutout` | `lsst.afw.image.Exposure` | The cutout. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `bad_bitmask` | `int` | Mask bits of bad pixels (shown in red). |
| `star_positions` | array-like, optional | Original PSF star positions to mark. |

### fitLib.plot_outer_step_background

`plot_outer_step_background(raw_exp, exp_next, src_mask, model_full, bg_next, bg_prev, outer_idx)`

Plot the background refit of one outer step.

Input image, model and source mask, new and previous background, and the result.

| Parameter | Type | Description |
| --- | --- | --- |
| `raw_exp` | `lsst.afw.image.Exposure` | Input to the background refit. |
| `exp_next` | `lsst.afw.image.Exposure` | Background-subtracted result. |
| `src_mask` | `numpy.ndarray` [`bool`] | Source mask used. |
| `model_full` | `numpy.ndarray` | Joint-fit model. |
| `bg_next` | `lsst.afw.image.Image` or `None` | New background. |
| `bg_prev` | `lsst.afw.image.Image` or `None` | Previous background. |
| `outer_idx` | `int` | Outer-pass index. |

| Returns | Type | Description |
| --- | --- | --- |
| `norm` | `astropy.visualization.ImageNormalize` | Display normalisation used for the images. |

### fitLib.plot_psf_candidate_stamps

`plot_psf_candidate_stamps(cellSet, max_stamps=24, title='')`

Plot the PSF candidate stamps, used ones first.

| Parameter | Type | Description |
| --- | --- | --- |
| `cellSet` | `lsst.afw.math.SpatialCellSet` or `None` | PSF candidates. |
| `max_stamps` | `int`, optional | Maximum number of stamps. |
| `title` | `str`, optional | Title suffix. |

### fitLib.plot_outer_step_psf

`plot_outer_step_psf(exp_next, sources, starCat, forced_positions, psf_prev, psf_new, psf_updated, outer_idx)`

Plot the PSF refit of one outer step: star selection and old vs new PSF.

| Parameter | Type | Description |
| --- | --- | --- |
| `exp_next` | `lsst.afw.image.Exposure` | Background-subtracted exposure. |
| `sources` | `lsst.afw.table.SourceCatalog` or `None` | Fresh detections. |
| `starCat` | `lsst.afw.table.SourceCatalog` or `None` | PSF candidates. |
| `forced_positions` | array-like or `None` | Forced original star positions. |
| `psf_prev, psf_new` | `lsst.afw.detection.Psf` | Previous and new PSF. |
| `psf_updated` | `bool` | Whether the PSF was updated. |
| `outer_idx` | `int` | Outer-pass index. |

### fitLib.highfreq_noise_check_arrays

`highfreq_noise_check_arrays(resid, variance, model, good, nbins=6, verbose=True)`

Compare the residual's pixel-to-pixel scatter with the variance plane.

Binned by model flux: ratio ~ 1 everywhere means the variance is right; a ratio growing with model flux means the variance lacks source photon noise.

| Parameter | Type | Description |
| --- | --- | --- |
| `resid` | `numpy.ndarray` | Residual image. |
| `variance` | `numpy.ndarray` | Variance image. |
| `model` | `numpy.ndarray` | Model image. |
| `good` | `numpy.ndarray` [`bool`] | Usable pixels. |
| `nbins` | `int`, optional | Number of model-flux bins. |
| `verbose` | `bool`, optional | Print the table. |

| Returns | Type | Description |
| --- | --- | --- |
| `result` | `dict` | ``flux`` (bin centres), ``ratio`` (scatter^2 / variance) and ``k`` (extra variance per unit model flux: var_true ~ var + k * model). |

### fitLib.highfreq_noise_check

`highfreq_noise_check(outer_history, outer_idx=-1, snap_idx=-1, **kw)`

Run `highfreq_noise_check_arrays` on one snapshot of a run.

| Parameter | Type | Description |
| --- | --- | --- |
| `outer_history` | `list` [`dict`] | Outer-loop history. |
| `outer_idx` | `int`, optional | Outer pass (default: last). |
| `snap_idx` | `int`, optional | Inner-loop snapshot (default: last). |
| `**kw` | – | Keyword arguments of `highfreq_noise_check_arrays`. |

| Returns | Type | Description |
| --- | --- | --- |
| `result` | `dict` | See `highfreq_noise_check_arrays`. |

### fitLib.estimate_psf_correction

`estimate_psf_correction(data, model, variance, good, xy0, positions, amps, kernel, n_use=8, step=0.8, clip=3.0, soft_sigma=1.0, window=(0.6, 1.0), max_radius=None)`

Estimate an improved PSF kernel from neighbour-subtracted residuals.

For each bright star, (data - model) / amp sampled on the kernel grid around its fitted position estimates (true PSF - current PSF).  The weighted, clipped mean over stars is soft-thresholded, windowed to zero at the kernel edge, stripped of the components degenerate with position (dP/dx, dP/dy), added to the kernel and renormalised.

| Parameter | Type | Description |
| --- | --- | --- |
| `data` | `numpy.ndarray` | Image. |
| `model` | `numpy.ndarray` | Joint-fit model. |
| `variance` | `numpy.ndarray` | Variance used for the weights. |
| `good` | `numpy.ndarray` [`bool`] | Usable pixels. |
| `xy0` | `tuple` [`int`, `int`] | Parent pixel coordinates of ``data[0, 0]``. |
| `positions` | array-like | Star positions, shape (N, 2). |
| `amps` | array-like | Star amplitudes. |
| `kernel` | `numpy.ndarray` | Current PSF kernel. |
| `n_use` | `int`, optional | Number of brightest stars used. |
| `step` | `float`, optional | Fraction of the correction applied. |
| `clip` | `float`, optional | Sigma clipping of the stack. |
| `soft_sigma` | `float`, optional | Soft threshold of the correction [sigma]. |
| `window` | `tuple` [`float`, `float`], optional | Inner and outer edge of the window as fractions of the kernel half-size. |
| `max_radius` | `float`, optional | Do not change the kernel beyond this radius [px]. |

| Returns | Type | Description |
| --- | --- | --- |
| `new` | `numpy.ndarray` | Corrected, normalised kernel. |
| `info` | `dict` | ``old``, ``new``, ``delta``, ``raw_delta``, ``err``, ``n_stars``, ``used``. |

### fitLib.kernel_to_lsst_psf

`kernel_to_lsst_psf(kernel)`

Wrap a kernel array as a spatially constant LSST PSF.

| Parameter | Type | Description |
| --- | --- | --- |
| `kernel` | `numpy.ndarray` | Square kernel of odd size, centred at n//2, summing to 1. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.meas.algorithms.KernelPsf` | The PSF. |

### fitLib.extend_kernel_with_moffat

`extend_kernel_with_moffat(kernel, new_size=81, fit_r=(6.0, 14.0), blend=(12.0, 16.0))`

Replace the noisy outer part of a PSF kernel by an analytic Moffat tail.

The azimuthal-mean profile between ``fit_r[0]`` and ``fit_r[1]`` px is fitted with ``A * (1 + (r/a)^2)^-b`` (relative weighting, so faint outer bins count).  The new kernel equals the original inside ``blend[0]``, the Moffat beyond ``blend[1]``, with a linear cross-fade between, and is renormalised -- fitted amplitudes then become total fluxes including the tail.  Works on an already-extended kernel too.

| Parameter | Type | Description |
| --- | --- | --- |
| `kernel` | `numpy.ndarray` | PSF kernel. |
| `new_size` | `int`, optional | Size of the new kernel [px] (odd). |
| `fit_r` | `tuple` [`float`, `float`], optional | Radial range of the Moffat fit [px]. |
| `blend` | `tuple` [`float`, `float`], optional | Radial range of the cross-fade [px]. |

| Returns | Type | Description |
| --- | --- | --- |
| `new_kernel` | `numpy.ndarray` | Extended, normalised kernel. |
| `info` | `dict` | Moffat ``A``, ``a``, ``beta`` and ``tail_frac`` (flux fraction beyond ``blend[1]``). |

Raises: `RuntimeError` (Raised if the Moffat fit fails.)

### fitLib.extend_psf_wings

`extend_psf_wings(psf, point, size=81, verbose=True, **kw)`

Return a constant KernelPsf with a Moffat tail, taken from any PSF at a point.

| Parameter | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | Input PSF. |
| `point` | `lsst.geom.Point2D` | Position at which the kernel is taken. |
| `size` | `int`, optional | Size of the extended kernel [px]. |
| `verbose` | `bool`, optional | Print the fitted tail. |
| `**kw` | – | Keyword arguments of `extend_kernel_with_moffat`. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | The extended PSF, or the input PSF if the fit fails. |

### fitLib.refine_psf_from_fit

`refine_psf_from_fit(outer_history, outer_idx=-1, snap_idx=-1, **kw)`

Make one PSF update step from a finished run.

| Parameter | Type | Description |
| --- | --- | --- |
| `outer_history` | `list` [`dict`] | Outer-loop history. |
| `outer_idx` | `int`, optional | Outer pass (default: last). |
| `snap_idx` | `int`, optional | Inner-loop snapshot (default: last). |
| `**kw` | – | Keyword arguments of `estimate_psf_correction`. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.meas.algorithms.KernelPsf` | Updated PSF. |
| `info` | `dict` | See `estimate_psf_correction`. |

### fitLib.polish_psf_fixed_sources

`polish_psf_fixed_sources(exposure, psf, result, bad_bitmask, n_cycles=3, stack_kwargs=None, support_frac=0.01, sys_frac=0.03, poisson_k=0.0, refine_centroids=True, centroid_min_amp=50.0, freeze_sep=None, max_shift=0.3, verbose=True)`

Iterate PSF correction and re-solve without re-detecting sources.

A cycle is kept only if the residual rms over good pixels drops.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psf` | `lsst.afw.detection.Psf` | Current PSF. |
| `result` | `dict` | Inner-loop result. |
| `bad_bitmask` | `int` | Mask bits of pixels to ignore. |
| `n_cycles` | `int`, optional | Maximum number of cycles. |
| `stack_kwargs` | `dict`, optional | Keyword arguments of `estimate_psf_correction`. |
| `support_frac, sys_frac, poisson_k` | `float`, optional | See `joint_psf_fit_diag`. |
| `refine_centroids` | `bool`, optional | Refine positions after each PSF change. |
| `centroid_min_amp` | `float`, optional | Fainter sources keep their position. |
| `freeze_sep` | `float`, optional | See `refine_peak_positions`. |
| `max_shift` | `float`, optional | Maximum position shift [px]. |
| `verbose` | `bool`, optional | Print each cycle. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | The best PSF. |
| `info` | `dict` or `None` | Information of the last kept correction (see `estimate_psf_correction`). |

### fitLib.plot_psf_correction

`plot_psf_correction(info)`

Plot the old kernel, the proposed correction and the radial profiles.

| Parameter | Type | Description |
| --- | --- | --- |
| `info` | `dict` | Output of `estimate_psf_correction`. |

### fitLib.run_cutout_pipeline

`run_cutout_pipeline(exposure_crop, psf, bad_bitmask, starCat=None, bgImage=None, x_range=(0, 100), y_range=(0, 100), n_outer_iterations=2, n_inner_iterations=10, inner_kwargs=None, outer_convergence_frac=0.02, refit_psf=True, psf_kernel_size=41, psf_kernel_sizes=None, bg_bin_size=32, bg_use_approx=True, viz=True, plot_every=1, source_view='all', sources_per_page=6, max_pages_per_iter=None, view_half=20, support_frac=0.01, refine_centroids=True, centroid_min_amp=50.0, n_centroid_rounds=4, centroid_shift_tol=0.03, centroid_shift_schedule=None, coarse_centroid=True, coarse_centroid_hw_sigma=2.5, merge_separation_sigma=1.5, sys_frac=0.03, min_snr=5.0, prune_sep_sigma=1.5, reject_fwhm=0.75, bright_peak_sigma=100.0, freeze_fwhm=1.0, poisson_k=0.0, near_bright_fwhm=2.0, near_bright_snr=10.0, edge_margin_fwhm=0.5, use_eliminate_weak=True, use_shape_cut=True, shape_kwargs=None, shape_in_first_pass=False, stack_kwargs=None, psf_polish_cycles=3, psf_wing_size=81)`

Run the outer loop (inner fit, background and PSF refinement) on a cutout.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure_crop` | `lsst.afw.image.Exposure` | Region to process (before background subtraction). |
| `psf` | `lsst.afw.detection.Psf` | Initial PSF. |
| `bad_bitmask` | `int` | Mask bits of pixels to ignore. |
| `starCat` | `lsst.afw.table.SourceCatalog`, optional | Original PSF stars, forced into every PSF determination. |
| `bgImage` | `lsst.afw.image.Image`, optional | Detector background subtracted before the first pass. |
| `x_range, y_range` | `tuple` [`int`, `int`], optional | Sub-region of ``exposure_crop`` (see `make_cutout`). |
| `n_outer_iterations` | `int`, optional | Number of outer passes. |
| `n_inner_iterations` | `int`, optional | Maximum inner iterations per pass. |
| `inner_kwargs` | `dict`, optional | Keyword arguments of `run_inner_loop_verbose`. |
| `outer_convergence_frac` | `float`, optional | Stop when the residual standard deviation improves by less than this fraction. |
| `refit_psf` | `bool` or `str`, optional | False, True (PCA star-selector refit) or 'stack' (residual-stack refit from the bright stars, see `refine_psf_from_fit`). |
| `psf_kernel_size` | `int`, optional | PSF kernel size [px]. |
| `psf_kernel_sizes` | `tuple` [`int`], optional | Kernel sizes tried by `redetect_and_refit_psf_adaptive`. |
| `bg_bin_size` | `int`, optional | Background bin size [px]. |
| `bg_use_approx` | `bool`, optional | See `refit_background`. |
| `viz` | `bool`, optional | Draw diagnostics. |
| `plot_every` | `int`, optional | Draw every n-th inner iteration (the last always). |
| `source_view` | `str`, optional | 'all', 'new' or 'changed': which sources get a per-source page. |
| `sources_per_page` | `int`, optional | Sources per page. |
| `max_pages_per_iter` | `int`, optional | Maximum per-source pages per iteration (None = no limit). |
| `view_half` | `int`, optional | Half-size of the per-source windows [px]. |
| `support_frac` | `float`, optional | See `joint_psf_fit_diag`. |
| `refine_centroids` | `bool`, optional | Refine positions (`refine_peak_positions`). |
| `centroid_min_amp` | `float`, optional | Fainter sources keep their position. |
| `n_centroid_rounds, centroid_shift_tol, centroid_shift_schedule` | optional | Centroid refinement controls (see `run_inner_loop_verbose`). |
| `coarse_centroid, coarse_centroid_hw_sigma` | optional | Windowed starting centroid of new peaks (see `run_inner_loop_verbose`). |
| `merge_separation_sigma` | `float`, optional | Merge new peaks closer than this many PSF sigma. |
| `sys_frac, min_snr, prune_sep_sigma, reject_fwhm, bright_peak_sigma` | `float`, optional | Over-counting guards (see `run_inner_loop_verbose`). |
| `freeze_fwhm, poisson_k, near_bright_fwhm, near_bright_snr` | `float`, optional | Over-counting guards (see `run_inner_loop_verbose`). |
| `edge_margin_fwhm` | `float`, optional | Drop sources within this many FWHM of the cutout edge. |
| `use_eliminate_weak` | `bool`, optional | Backward elimination on the honest S/N. |
| `use_shape_cut` | `bool`, optional | PSF-shape constraint (`eliminate_nonpsf`). |
| `shape_kwargs` | `dict`, optional | Keyword arguments of `shape_metrics` / `judge_shape`. |
| `shape_in_first_pass` | `bool`, optional | Apply the shape cut in the first outer pass too (by default only once the PSF has been refined: a rough first PSF can make the cut reject real faint companions). |
| `stack_kwargs` | `dict`, optional | Keyword arguments of `refine_psf_from_fit` (``refit_psf='stack'``). |
| `psf_polish_cycles` | `int`, optional | PSF correction / re-solve cycles per outer pass (``refit_psf='stack'``). |
| `psf_wing_size` | `int`, optional | Convert every PSF to a constant kernel of this size with a Moffat tail (see `extend_psf_wings`); None = off.  Stops the model ending at the stamp edge, which leaves a ring of fake sources half a kernel width from bright stars. |

| Returns | Type | Description |
| --- | --- | --- |
| `outer_history` | `list` [`dict`] | One dict per outer pass: ``psf``, ``exposure``, ``result``, ``resid_std``, ``n_peaks``, ``bg_image``. |

### fitLib.plot_wing_diagnostic

`plot_wing_diagnostic(outer_history, n_stars=3, rmax=40)`

Plot the azimuthal profiles of the brightest sources: data vs their own model.

Data beyond the dashed line well above the model means a real halo the kernel cannot represent; a difference centred on zero without a step means only core / PSF-shape mismatch remains.

| Parameter | Type | Description |
| --- | --- | --- |
| `outer_history` | `list` [`dict`] | Outer-loop history. |
| `n_stars` | `int`, optional | Number of brightest sources shown. |
| `rmax` | `float`, optional | Maximum radius [px]. |

### fitLib.plot_source_cutouts

`plot_source_cutouts(outer_history, outer_idx=-1, snap_idx=-1, source_ids=None, top_n=None, min_amp=0.0, half=15, cross_axis='x', show_own=True, sources_per_fig=4, pdf_path=None, show=True, cmap=None, mark_neighbours=True)`

Plot cutouts and cross-sections of individual sources of a run.

| Parameter | Type | Description |
| --- | --- | --- |
| `outer_history` | `list` [`dict`] | Outer-loop history. |
| `outer_idx` | `int`, optional | Outer pass (default: last). |
| `snap_idx` | `int`, optional | Inner-loop snapshot (default: last). |
| `source_ids` | `list` [`int`], optional | Sources to show; overrides ``top_n`` / ``min_amp``. |
| `top_n` | `int`, optional | Show the brightest N sources. |
| `min_amp` | `float`, optional | Skip sources with amplitude at or below this. |
| `half` | `int`, optional | Cutout half-width [px]. |
| `cross_axis` | `str`, optional | 'x', 'y' or 'both'. |
| `show_own` | `bool`, optional | Draw the source's own PSF component in the cross-section. |
| `sources_per_fig` | `int`, optional | Sources per figure. |
| `pdf_path` | `str`, optional | Also save every figure as a page of this PDF. |
| `show` | `bool`, optional | Display the figures. |
| `cmap` | `str`, optional | Colour map of the images. |
| `mark_neighbours` | `bool`, optional | Mark the neighbouring sources. |

| Returns | Type | Description |
| --- | --- | --- |
| `figs` | `list` [`matplotlib.figure.Figure`] | The figures. |

### fitLib.negative_source_check

`negative_source_check(outer_history, outer_idx=-1, snap_idx=-1, thresholds=(3, 4, 5), plot=True)`

Estimate the number of false detections from the residual itself.

Counts pixels of the matched-filter S/N map of the final residual at >= +N sigma and at <= -N sigma.  Real sources are positive only, so the negative excursions approximate the positive noise excursions (false positives); PSF-mismatch dipoles inflate both.

| Parameter | Type | Description |
| --- | --- | --- |
| `outer_history` | `list` [`dict`] | Outer-loop history. |
| `outer_idx` | `int`, optional | Outer pass (default: last). |
| `snap_idx` | `int`, optional | Inner-loop snapshot (default: last). |
| `thresholds` | `tuple` [`float`], optional | S/N thresholds. |
| `plot` | `bool`, optional | Plot the S/N map and histogram. |

| Returns | Type | Description |
| --- | --- | --- |
| `snr_map` | `numpy.ndarray` | Matched-filter S/N map of the residual. |

### fitLib.source_snr_histogram

`source_snr_histogram(outer_history, outer_idx=-1, snap_idx=-1, min_snr=5.0)`

Plot the histogram of the honest S/N of every final source.

A real population shows a gap or knee; a smooth pile-up just above ``min_snr`` means noise.

| Parameter | Type | Description |
| --- | --- | --- |
| `outer_history` | `list` [`dict`] | Outer-loop history. |
| `outer_idx` | `int`, optional | Outer pass (default: last). |
| `snap_idx` | `int`, optional | Inner-loop snapshot (default: last). |
| `min_snr` | `float`, optional | S/N threshold marked in the plot. |

| Returns | Type | Description |
| --- | --- | --- |
| `snr` | `numpy.ndarray` | S/N per source. |

### fitLib.draw_edge_margin

`draw_edge_margin(ax, bbox, margin_px, color='#ff00ff', shade=True)`

Outline and shade the band within ``margin_px`` of the image border.

This is the region removed by ``drop_out_of_bounds(..., margin=-margin_px)``.

| Parameter | Type | Description |
| --- | --- | --- |
| `ax` | `matplotlib.axes.Axes` | Axes to draw on. |
| `bbox` | `lsst.geom.Box2I` | Image box. |
| `margin_px` | `float` | Band width [px]. |
| `color` | `str`, optional | Line colour. |
| `shade` | `bool`, optional | Shade the band. |

### fitLib.plot_final_comparison

`plot_final_comparison(outer_history, exposure_crop, bad_bitmask=None, edge_margin_fwhm=0.5, view=None, show_peaks=True, circle_radius=None, figsize=(26, 7))`

Plot original | final model | residual | residual/sigma over the fitted region.

| Parameter | Type | Description |
| --- | --- | --- |
| `outer_history` | `list` [`dict`] | Outer-loop history. |
| `exposure_crop` | `lsst.afw.image.Exposure` | The input crop (original image panel). |
| `bad_bitmask` | `int`, optional | Mask bits of bad pixels. |
| `edge_margin_fwhm` | `float`, optional | Edge margin used in the run; draws its boundary (None / 0 = do not draw). |
| `view` | `tuple` [`float`, `float`, `float`, `float`], optional | (xmin, xmax, ymin, ymax) zoom. |
| `show_peaks` | `bool`, optional | Mark the sources. |
| `circle_radius` | `float`, optional | Marker radius [px]. |
| `figsize` | `tuple` [`float`, `float`], optional | Figure size [inches]. |

### fitLib.plot_source_pages

`plot_source_pages(exposure, snap, outer_idx, pages=None, colors=None, source_view='all', sources_per_page=6, max_pages=None, view_half=20, circle_radius=None, pdf=None, show=True, snr=None, dpi=80)`

Plot one row per source: data | model | residual | x-section | y-section.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `snap` | `dict` | Inner-loop snapshot. |
| `outer_idx` | `int` | Outer-pass index. |
| `pages` | `list` [`list` [`int`]], optional | Source indices per page (default: `select_source_pages`). |
| `colors` | `dict` [`int`, `str`], optional | Colour per source. |
| `source_view` | `str`, optional | See `select_source_pages`. |
| `sources_per_page` | `int`, optional | See `select_source_pages`. |
| `max_pages` | `int`, optional | See `select_source_pages`. |
| `view_half` | `int`, optional | Half-size of each window [px]. |
| `circle_radius` | `float`, optional | Marker radius [px]. |
| `pdf` | `matplotlib.backends.backend_pdf.PdfPages`, optional | Open PDF; every page is saved into it. |
| `show` | `bool`, optional | Also display each page (slow for hundreds of sources). |
| `snr` | `numpy.ndarray`, optional | Honest S/N per source, shown in the model title. |
| `dpi` | `int`, optional | Resolution of the PDF pages. |

### fitLib.source_metrics

`source_metrics(outer_history, outer_idx=-1, snap_idx=-1, core_fwhm=1.5)`

Return per-source quantities used to stratify the sample.

| Parameter | Type | Description |
| --- | --- | --- |
| `outer_history` | `list` [`dict`] | Outer-loop history. |
| `outer_idx` | `int`, optional | Outer pass (default: last). |
| `snap_idx` | `int`, optional | Inner-loop snapshot (default: last). |
| `core_fwhm` | `float`, optional | Radius of the chi2 region [FWHM]. |

| Returns | Type | Description |
| --- | --- | --- |
| `metrics` | `dict` | ``peaks``, ``amps``, ``snr`` (honest S/N), ``nn`` (distance to the nearest other source [px]), ``edge`` (distance to the image border [px]), ``chi2`` (mean resid^2 / sigma^2 within ``core_fwhm`` FWHM), ``group`` (size of the coupled group), ``fwhm``, ``outer_pos`` and the exposure and snapshot used. |

### fitLib.select_representative_sources

`select_representative_sources(m, n_per_category=4, seed=0, categories=None, extra_ids=None)`

Select a representative sample of sources per category.

| Parameter | Type | Description |
| --- | --- | --- |
| `m` | `dict` | Output of `source_metrics`. |
| `n_per_category` | `int` or `dict` [`str`, `int`], optional | Sources drawn per category. |
| `seed` | `int`, optional | Random seed of the tercile draws. |
| `categories` | `list` [`str`], optional | Categories to keep, in display order. |
| `extra_ids` | `list` [`int`], optional | Source ids appended as a 'picked' category. |

| Returns | Type | Description |
| --- | --- | --- |
| `selection` | `collections.OrderedDict` [`str`, `list` [`int`]] | Source ids per category. |

### fitLib.plot_representative_sources

`plot_representative_sources(outer_history, outer_idx=-1, snap_idx=-1, n_per_category=4, seed=0, categories=None, extra_ids=None, sources_per_page=4, view_half=20, pdf_path=None, show=True, dpi=80)`

Select a representative sample and plot it with `plot_source_pages`.

Each category is split into pages of ``sources_per_page`` rows, so a page never mixes categories; the printed table maps every source id to its category.

| Parameter | Type | Description |
| --- | --- | --- |
| `outer_history` | `list` [`dict`] | Outer-loop history. |
| `outer_idx` | `int`, optional | Outer pass (default: last). |
| `snap_idx` | `int`, optional | Inner-loop snapshot (default: last). |
| `n_per_category` | `int` or `dict`, optional | See `select_representative_sources`. |
| `seed` | `int`, optional | See `select_representative_sources`. |
| `categories` | `list` [`str`], optional | See `select_representative_sources`. |
| `extra_ids` | `list` [`int`], optional | See `select_representative_sources`. |
| `sources_per_page` | `int`, optional | Rows per page. |
| `view_half` | `int`, optional | Half-size of each window [px]. |
| `pdf_path` | `str`, optional | Also write all pages to this PDF. |
| `show` | `bool`, optional | Display the pages. |
| `dpi` | `int`, optional | Resolution of the PDF pages. |

| Returns | Type | Description |
| --- | --- | --- |
| `selection` | `collections.OrderedDict` [`str`, `list` [`int`]] | Source ids per category. |

## adaptive

`lsst.meas.crowd.adaptive`: the original cutout_adaptive: the adaptive inner-loop decisions.

Residual-driven adaptive control for the cutout pipeline (the original cutout_adaptive.py).

| Name | Summary |
| --- | --- |
| [`residual_excess`](#adaptiveresidual_excess) | Count positive and negative peaks of the residual's matched-filter S/N. |
| [`counts_at`](#adaptivecounts_at) | Return the positive and negative peak counts at an S/N threshold. |
| [`is_missing`](#adaptiveis_missing) | Decide whether the residual still contains missing sources at an S/N threshold. |
| [`pick_threshold`](#adaptivepick_threshold) | Return the lowest S/N threshold the field can support. |
| [`install`](#adaptiveinstall) | Route fitLib's inner loop through the adaptive controller. |
| [`uninstall`](#adaptiveuninstall) | Restore fitLib's original inner and outer loops. |
| [`report`](#adaptivereport) | Print what the adaptive controller did in every outer pass of a finished run. |

### adaptive.residual_excess

`residual_excess(resid, var, good, peaks_xy, amps, x0, y0, fwhm, psf_peak, noise, thresholds=(3.0, 3.5, 4.0, 4.5, 5.0, 6.0), edge_fwhm=2.0, bright_sigma=100.0, halo_fwhm=2.5, min_sep_fwhm=0.75)`

Count positive and negative peaks of the residual's matched-filter S/N.

Only local maxima in the clean region count: away from the cutout edge, from bright sources (halo radius growing with flux) and from fitted sources themselves.  Real missing sources are positive only, so the positive excess over the negative count estimates how many are missing.

| Parameter | Type | Description |
| --- | --- | --- |
| `resid` | `numpy.ndarray` | Residual image. |
| `var` | `numpy.ndarray` | Variance image. |
| `good` | `numpy.ndarray` [`bool`] | Usable pixels. |
| `peaks_xy` | array-like | Fitted source positions (parent pixel coordinates), shape (N, 2). |
| `amps` | array-like | Fitted source amplitudes. |
| `x0, y0` | `int` | Parent pixel coordinates of ``resid[0, 0]``. |
| `fwhm` | `float` | PSF FWHM [px]. |
| `psf_peak` | `float` | Peak value of the unit-flux PSF. |
| `noise` | `float` | Typical pixel noise. |
| `thresholds` | `tuple` [`float`], optional | S/N thresholds at which to count. |
| `edge_fwhm` | `float`, optional | Ignore peaks within this many FWHM of the edge. |
| `bright_sigma` | `float`, optional | A source is bright if ``amp * psf_peak / noise`` exceeds this. |
| `halo_fwhm` | `float`, optional | Ignore peaks within this many FWHM of a bright source (grows with flux). |
| `min_sep_fwhm` | `float`, optional | Peaks closer than this to a fitted source are its own mismatch residual. |

| Returns | Type | Description |
| --- | --- | --- |
| `result` | `dict` | ``t`` (thresholds, ascending), ``pos`` / ``neg`` (number of peaks with \|S/N\| >= t), ``scale``, ``offset`` (S/N calibration), ``area_frac`` (fraction of the image in the clean region) and ``candidates`` ((x, y, snr, distance to the nearest fitted source) of the positive peaks). |

### adaptive.counts_at

`counts_at(d, t)`

Return the positive and negative peak counts at an S/N threshold.

| Parameter | Type | Description |
| --- | --- | --- |
| `d` | `dict` | Output of `residual_excess`. |
| `t` | `float` | S/N threshold (the nearest one in ``d['t']`` is used). |

| Returns | Type | Description |
| --- | --- | --- |
| `n_pos` | `int` | Number of positive peaks. |
| `n_neg` | `int` | Number of negative peaks. |

### adaptive.is_missing

`is_missing(d, t, z_sig=2.0, min_excess=3)`

Decide whether the residual still contains missing sources at an S/N threshold.

| Parameter | Type | Description |
| --- | --- | --- |
| `d` | `dict` | Output of `residual_excess`. |
| `t` | `float` | S/N threshold. |
| `z_sig` | `float`, optional | The excess must exceed ``z_sig * sqrt(N_pos + N_neg)`` ... |
| `min_excess` | `int`, optional | ... and be at least this many sources. |

| Returns | Type | Description |
| --- | --- | --- |
| `missing` | `bool` | Whether sources are missing. |
| `excess` | `int` | Positive minus negative peak count. |
| `error` | `float` | Its uncertainty, sqrt(N_pos + N_neg). |

### adaptive.pick_threshold

`pick_threshold(d, fdr, t_floor, t_cap)`

Return the lowest S/N threshold the field can support.

The lowest S/N in [``t_floor``, ``t_cap``] whose estimated false fraction (N_neg + 1) / N_pos is at most ``fdr`` (with N_pos >= 3).

| Parameter | Type | Description |
| --- | --- | --- |
| `d` | `dict` | Output of `residual_excess`. |
| `fdr` | `float` | Tolerated false fraction. |
| `t_floor` | `float` | Lowest threshold allowed. |
| `t_cap` | `float` | Highest threshold; returned if no lower one qualifies. |

| Returns | Type | Description |
| --- | --- | --- |
| `threshold` | `float` | The chosen S/N threshold. |

### adaptive.install

`install(**overrides)`

Route fitLib's inner loop through the adaptive controller.

Also enables the adaptive outer loop.  Requires the hooks added by ``patch_cutout_fit_lib.py``.

| Parameter | Type | Description |
| --- | --- | --- |
| `**overrides` | – | Settings that replace entries of ``ADAPT``. |

Raises: `ValueError` (Raised if an override is not an ``ADAPT`` option.); `RuntimeError` (Raised if fitLib has not been patched with the hooks.)

### adaptive.uninstall

`uninstall()`

Restore fitLib's original inner and outer loops.

### adaptive.report

`report(run_or_history)`

Print what the adaptive controller did in every outer pass of a finished run.

| Parameter | Type | Description |
| --- | --- | --- |
| `run_or_history` | `CutoutRun` or `list` [`dict`] | A finished run or its outer-loop history. |
