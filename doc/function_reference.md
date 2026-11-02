# meas_crowd: function reference

Every public function and class of `lsst.meas.crowd`: what it does, its parameters (with types and defaults) and what it returns.  Generated from the docstrings by `doc/make_function_reference.py`; do not edit by hand (run `python doc/make_function_reference.py` after changing a docstring).  How the pipeline fits together is in [pipeline_guide.md](pipeline_guide.md); the output catalogue in [catalog_columns.md](catalog_columns.md).

181 functions and 4 classes in 12 modules.  Everything listed here can be used as `crowd.<name>` after `import lsst.meas.crowd as crowd`, or from its module (`from lsst.meas.crowd import photometry`).  Names starting with `_` are internal helpers and are not listed.

## Start here

| Function | What it is for |
| --- | --- |
| [`session.makeConfig`](#sessionmakeconfig) | Build a pipeline configuration (keyword arguments of `fitLib.runCutoutPipeline`). |
| [`session.prepareSession`](#sessionpreparesession) | Load a visit image, cut out the regions and set up the detector-level PSF. |
| [`gaiaCatalog.attachGaia`](#gaiacatalogattachgaia) | Load Gaia once for a session and store it as ``session.gaia``. |
| [`pipeline.runSessionGaia`](#pipelinerunsessiongaia) | Run the adaptive pipeline on every (or the named) cutout of a session. |
| [`pipeline.report`](#pipelinereport) | Print what the adaptive loops did, pass by pass. |
| [`photometry.measurePhotometry`](#photometrymeasurephotometry) | Measure PSF and aperture magnitudes of a run's sources. |
| [`plotting.plotFinalComparison`](#plottingplotfinalcomparison) | Plot original \| final model \| residual \| residual/sigma of a run. |
| [`plotting.plotGaiaComparison`](#plottingplotgaiacomparison) | Plot how the fitted photometry and positions compare with Gaia. |
| [`stagePlots.showRounds`](#stageplotsshowrounds) | Plot data \| model \| residual \| residual/sigma for each round of the fit. |
| [`stagePlots.showLinkedSourcePages`](#stageplotsshowlinkedsourcepages) | Plot per-source pages whose coloured frames match boxes on the data, model and residual. |
| [`stagePlots.comparePasses`](#stageplotscomparepasses) | Find the sources that one outer pass lost or gained relative to another, and why. |
| [`detector.fitTiles`](#detectorfittiles) | Fit this task's share of a detector's tiles and write one file per tile. |
| [`detector.stitchTiles`](#detectorstitchtiles) | Assemble the tile files into full-detector images and one catalogue. |
| [`detector.plotDetector`](#detectorplotdetector) | Plot data \| model \| residual \| residual/sigma of a stitched detector or a zoom. |
| [`runPlots.saveRunPlots`](#runplotssaverunplots) | Write every diagnostic plot of one finished run to PNG files. |

## Modules

| Module | Role | Public functions / classes |
| --- | --- | --- |
| [`pipeline`](#pipeline) | the driver: run the pipeline on cutouts from Gaia starting positions, inner / outer loops, PSF refinement, diagnostics | 40 |
| [`session`](#session) | load a visit, cut out regions, set up the PSF; configuration; saving runs | 19 |
| [`gaiaCatalog`](#gaiacatalog) | the Gaia catalogue and the starting positions of each cutout | 4 |
| [`photometry`](#photometry) | PSF and aperture photometry, errors, the Gaia comparison | 9 |
| [`plotting`](#plotting) | final data / model / residual plots and residual diagnostics | 9 |
| [`stagePlots`](#stageplots) | plots of the intermediate stages and pass-to-pass comparison | 6 |
| [`sourcePlots`](#sourceplots) | per-source cutouts, cross-sections and magnitudes | 3 |
| [`detector`](#detector) | whole-detector runs: tiles, stitching, mosaic and tile plots | 10 |
| [`runPlots`](#runplots) | write every plot of a run, of the detector-level PSF, or of chosen tiles | 4 |
| [`fastKernels`](#fastkernels) | fast PSF stamps (C++ with a NumPy fallback) | 6 |
| [`fitLib`](#fitlib) | the core library: joint fit, detection, PSF and background tools, plots | 68 |
| [`adaptive`](#adaptive) | the residual check deciding whether sources are still missing | 7 |

## pipeline

`lsst.meas.crowd.pipeline`: the driver: run the pipeline on cutouts from Gaia starting positions, inner / outer loops, PSF refinement, diagnostics.

Adaptive, iterative joint PSF-fitting pipeline starting from Gaia positions.

| Name | Summary |
| --- | --- |
| [`clearStampCache`](#pipelineclearstampcache) | Empty the PSF stamp cache and the stamp-size rule. |
| [`setStampSizes`](#pipelinesetstampsizes) | Choose a PSF stamp half-width per source. |
| [`ampSnrFullFast`](#pipelineampsnrfullfast) | Return honest amplitude S/N values, computed block by block. |
| [`profile`](#pipelineprofile) | Run a function under cProfile and print where the time goes. |
| [`enableSpeedups`](#pipelineenablespeedups) | Patch fitLib with the cached / fast PSF stamps, fast S/N and fast NNLS (idempotent). |
| [`disableSpeedups`](#pipelinedisablespeedups) | Undo `enableSpeedups` and empty the caches. |
| [`estimateGaiaOffset`](#pipelineestimategaiaoffset) | Measure the global offset between the image and the Gaia positions. |
| [`saturationLevel`](#pipelinesaturationlevel) | Estimate the saturation level of an image. |
| [`measureModelError`](#pipelinemeasuremodelerror) | Measure the PSF model error as a fraction of the model per pixel. |
| [`satZones`](#pipelinesatzones) | Find the zones around saturated stars. |
| [`satCoreLabels`](#pipelinesatcorelabels) | Label the saturated cores, each grown by max(3 px, FWHM/2). |
| [`brightStarCorrections`](#pipelinebrightstarcorrections) | Fit local correction terms for every bright or saturated star. |
| [`printBrightCorrections`](#pipelineprintbrightcorrections) | Print the table of bright-star corrections. |
| [`ringMimics`](#pipelineringmimics) | Find faint non-Gaia sources that are fakes absorbing a bright star's wing excess. |
| [`splitTest`](#pipelinesplittest) | Split unresolved pairs. |
| [`ellipTest`](#pipelineelliptest) | Model badly fitted sources as small galaxies. |
| [`tryElliptical`](#pipelinetryelliptical) | Try the elliptical-source test on one source of a finished run, printing every step. |
| [`finalPositionPolish`](#pipelinefinalpositionpolish) | Re-optimise all source positions jointly. |
| [`runInnerAdaptive`](#pipelineruninneradaptive) | Run one outer pass of the adaptive pipeline. |
| [`selectPsfStars`](#pipelineselectpsfstars) | Select clean, bright, isolated, unsaturated point sources for the PSF stack. |
| [`refinePsfClean`](#pipelinerefinepsfclean) | Update the PSF from the residual stack of clean stars, cross-validated. |
| [`refinePsfMedian`](#pipelinerefinepsfmedian) | Correct the PSF core with the median residual/flux of the bright clean stars. |
| [`smoothPsfKernel`](#pipelinesmoothpsfkernel) | Return a hole-free, smooth version of a PSF kernel. |
| [`regularizePsf`](#pipelineregularizepsf) | Return the PSF with its kernel passed through `smoothPsfKernel`. |
| [`extendPsfWingsKeep`](#pipelineextendpsfwingskeep) | Extend a PSF with a Moffat tail, optionally keeping its angular structure. |
| [`wingProfilePsf`](#pipelinewingprofilepsf) | Return a PSF whose outer wings are measured from the stars. |
| [`spatialKernelPsf`](#pipelinespatialkernelpsf) | Build a linearly varying LSST PSF. |
| [`refinePsfSpatial`](#pipelinerefinepsfspatial) | Add or update a linear spatial variation of the PSF, cross-validated. |
| [`fitBrightnessTerm`](#pipelinefitbrightnessterm) | Measure or update the brightness-dependent PSF term (brighter-fatter). |
| [`refitBackgroundLowfit`](#pipelinerefitbackgroundlowfit) | Fit a constant or tilted-plane background to the unmasked pixels and subtract it in place. |
| [`printTiming`](#pipelineprinttiming) | Print the time spent per pipeline step. |
| [`runAdaptivePipeline`](#pipelinerunadaptivepipeline) | Run the adaptive pipeline on a whole cutout, starting from Gaia positions. |
| [`dropOffDetectorGaia`](#pipelinedropoffdetectorgaia) | Drop Gaia starting positions that lie beyond the detector edge. |
| [`runCutoutGaia`](#pipelineruncutoutgaia) | Run the adaptive pipeline on one cropped exposure, starting from Gaia positions. |
| [`checkBrightNeighbours`](#pipelinecheckbrightneighbours) | Test whether sources on the edges of bright stars are real. |
| [`diagnoseResidual`](#pipelinediagnoseresidual) | Decide whether what is left around bright stars is noise or model error. |
| [`modeledNear`](#pipelinemodelednear) | Return the pixels closer to a modelled (near-saturated) core than to any other masked pixel. |
| [`unmodeledNearMask`](#pipelineunmodelednearmask) | Find the pixels whose residual cannot be judged because a source there is partly masked. |
| [`runSessionGaia`](#pipelinerunsessiongaia) | Run the adaptive pipeline on every (or the named) cutout of a session. |
| [`report`](#pipelinereport) | Print what the adaptive loops did, pass by pass. |

### pipeline.clearStampCache

`clearStampCache()`

Empty the PSF stamp cache and the stamp-size rule.

Call after changing the PSF in place; a new PSF object gets its own cache automatically.

### pipeline.setStampSizes

`setStampSizes(xy, amps, psf, point, noise, fwhm, floorSigma=0.05, tol=1.5)`

Choose a PSF stamp half-width per source.

The smallest of a few tiers beyond which ``amp * PSF`` stays below ``floorSigma * noise``, so faint sources get small stamps.  Sources not in the list (new detections) get the full stamp.

| Parameter | Type | Description |
| --- | --- | --- |
| `xy` | `numpy.ndarray` | Source positions, shape (N, 2). |
| `amps` | `numpy.ndarray` | Source amplitudes. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `point` | `lsst.geom.Point2D` | Position at which the PSF profile is evaluated. |
| `noise` | `float` | Typical pixel noise. |
| `fwhm` | `float` | PSF FWHM [px]. |
| `floorSigma` | `float`, optional | Crop where the source's wing drops below this many noise sigma. |
| `tol` | `float`, optional | Positions within this many px of a listed source use its stamp size. |

| Returns | Type | Description |
| --- | --- | --- |
| `half` | `numpy.ndarray` [`int`] | Half-width per source (-1 = full stamp). |

### pipeline.ampSnrFullFast

`ampSnrFullFast(fit, noiseScale=1.0)`

Return honest amplitude S/N values, computed block by block.

Drop-in for `fitLib.ampSnrFull`: ``S/N = amp / sqrt(diag((A^T W A)^-1))``, inverted per independent group of overlapping sources with Cholesky instead of one dense pseudo-inverse.

| Parameter | Type | Description |
| --- | --- | --- |
| `fit` | `dict` | Joint fit (from `fitLib.jointPsfFitDiag`). |
| `noiseScale` | `float`, optional | Empirical noise scale dividing the S/N. |

| Returns | Type | Description |
| --- | --- | --- |
| `snr` | `numpy.ndarray` | S/N per source. |

### pipeline.profile

`profile(fn, *args, top=30, **kw)`

Run a function under cProfile and print where the time goes.

For example ``runs = crowd.profile(crowd.runSessionGaia, sess, cfg, maxOuter=2)``.

| Parameter | Type | Description |
| --- | --- | --- |
| `fn` | callable | Function to profile. |
| `*args` | – | Positional arguments of ``fn``. |
| `top` | `int`, optional | Number of entries printed. |
| `**kw` | – | Keyword arguments of ``fn``. |

| Returns | Type | Description |
| --- | --- | --- |
| `result` | `object` | What ``fn`` returned. |

### pipeline.enableSpeedups

`enableSpeedups()`

Patch fitLib with the cached / fast PSF stamps, fast S/N and fast NNLS (idempotent).

Called by the pipeline itself; also speeds up fitLib's own `fitLib.runCutoutPipeline`.

### pipeline.disableSpeedups

`disableSpeedups()`

Undo `enableSpeedups` and empty the caches.

### pipeline.estimateGaiaOffset

`estimateGaiaOffset(exposure, S, fwhm, badBitmask, noise, minPeakSnr=20.0, minStars=5, verbose=True)`

Measure the global offset between the image and the Gaia positions.

The median (data centroid - Gaia position) of bright, unsaturated, fairly isolated Gaia sources, from background-subtracted windowed centroids.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The cutout. |
| `S` | `dict` | Source state with the Gaia starting positions (see ``_state``). |
| `fwhm` | `float` | PSF FWHM [px]. |
| `badBitmask` | `int` | Mask bits of pixels to ignore. |
| `noise` | `float` | Typical pixel noise. |
| `minPeakSnr` | `float`, optional | Minimum peak S/N of a star used. |
| `minStars` | `int`, optional | Minimum number of stars; otherwise no offset is applied. |
| `verbose` | `bool`, optional | Print the result. |

| Returns | Type | Description |
| --- | --- | --- |
| `dx, dy` | `float` | Offset (data - Gaia) [px]; 0 if too few stars. |
| `nUsed` | `int` | Number of stars used. |

### pipeline.saturationLevel

`saturationLevel(exposure, P)`

Estimate the saturation level of an image.

``nearSatLevel`` if given; else the median of the SAT-flagged pixels, unless those look interpolated (below the brightest unflagged pixel), in which case the brightest unflagged pixel, a lower bound.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The image. |
| `P` | `dict` | Pipeline options (uses ``nearSatLevel``). |

| Returns | Type | Description |
| --- | --- | --- |
| `level` | `float` or `None` | Saturation level in image units (None if it cannot be estimated). |
| `how` | `str` | How it was obtained. |

### pipeline.measureModelError

`measureModelError(exposure, res, P, verbose=True)`

Measure the PSF model error as a fraction of the model per pixel.

Around unsaturated point sources with peak S/N >= ``modelErrMinSnr``, the robust (MAD) variance of residual/sigma in bins of model value is fitted as ``1 + eps^2 * model^2 / var``. Values not significant at 2 sigma give eps = 0.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `res` | `dict` | Result of `runInnerAdaptive`. |
| `P` | `dict` | Pipeline options. |
| `verbose` | `bool`, optional | Print the result. |

| Returns | Type | Description |
| --- | --- | --- |
| `eps` | `float` | Model error (0 if not significant). |
| `table` | `list` [`tuple`] | Per model bin: (model, number of pixels, robust variance of residual/sigma, model^2/var). |

### pipeline.satZones

`satZones(exposure, fwhm, P)`

Find the zones around saturated stars.

The zone band follows the SAT mask's shape: pixels within ``w = satZoneFwhm * FWHM + satZoneScale * inscribed radius`` of each blob.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The cutout. |
| `fwhm` | `float` | PSF FWHM [px]. |
| `P` | `dict` | Pipeline options (``sat_zone_*``). |

| Returns | Type | Description |
| --- | --- | --- |
| `centres` | `numpy.ndarray` | (x, y, halo radius) per saturated blob, shape (N, 3). |
| `mask` | `numpy.ndarray` [`bool`] | The zone band. |

### pipeline.satCoreLabels

`satCoreLabels(exposure, fwhm, minPx=4)`

Label the saturated cores, each grown by max(3 px, FWHM/2).

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The cutout. |
| `fwhm` | `float` | PSF FWHM [px]. |
| `minPx` | `int`, optional | Ignore saturated blobs smaller than this (hot pixels, cosmic rays). |

| Returns | Type | Description |
| --- | --- | --- |
| `labels` | `numpy.ndarray` [`int`] | Core label per pixel (0 = none). |

### pipeline.brightStarCorrections

`brightStarCorrections(exposure, S, fit, model, psfPeak, noise, fwhm, badBitmask, sigmaMin=100.0, radiusFwhm=4.0, sweeps=2, radial=True, zones=None, haloModes=0, brightModes=0, coreRings=(), coreSigma=0.6, coreModes=0)`

Fit local correction terms for every bright or saturated star.

On top of the joint model, within ``radiusFwhm`` FWHM, each bright star's residual is fitted with ``amp * (a*P + bx*dP/dx + by*dP/dy + cxx*d2P/dx2 + cyy*d2P/dy2 + cxy*d2P/dxdy)`` (flux, position, size and ellipticity), optional radial rings with angular modes (halo and wing structure) and optional narrow core terms (core sharper than the PSF, possibly off-centre). Weighted by the pixel variance (not the sysFrac-inflated one) so the core, where the mismatch is, counts.  Two sweeps, brightest first, for overlapping bright stars.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `S` | `dict` | Source state. |
| `fit` | `dict` | Joint fit. |
| `model` | `numpy.ndarray` | Joint model. |
| `psfPeak` | `float` | Peak value of the unit-flux PSF. |
| `noise` | `float` | Typical pixel noise. |
| `fwhm` | `float` | PSF FWHM [px]. |
| `badBitmask` | `int` | Mask bits of pixels to ignore. |
| `sigmaMin` | `float`, optional | Stars with peak S/N above this (or saturated) are corrected. |
| `radiusFwhm` | `float`, optional | Fit radius [FWHM]. |
| `sweeps` | `int`, optional | Number of sweeps. |
| `radial` | `bool`, optional | Add radial ring terms. |
| `zones` | `numpy.ndarray`, optional | Saturated-star zones (see `satZones`); stars in them get halo terms out to the zone. |
| `haloModes` | `int`, optional | Angular modes of the saturated stars' rings. |
| `brightModes` | `int`, optional | Angular modes of the unsaturated stars' rings. |
| `coreRings` | `tuple` [`float`], optional | Radii [px] of the narrow core terms (unsaturated stars). |
| `coreSigma` | `float`, optional | Width [px] of the narrow core terms. |
| `coreModes` | `int`, optional | Angular modes of the narrow core terms. |

| Returns | Type | Description |
| --- | --- | --- |
| `corr` | `numpy.ndarray` | Correction image (add to the joint model). |
| `info` | `dict` [`int`, `dict`] | Per corrected star: ``x``, ``y``, ``snr``, ``nSat``, ``flux``, ``dfluxPct``, ``dx``, ``dy``, ``chi2Before``, ``chi2After``. |

### pipeline.printBrightCorrections

`printBrightCorrections(info, maxRows=15)`

Print the table of bright-star corrections.

| Parameter | Type | Description |
| --- | --- | --- |
| `info` | `dict` [`int`, `dict`] | Second output of `brightStarCorrections`. |
| `maxRows` | `int`, optional | Number of (brightest) stars printed. |

### pipeline.ringMimics

`ringMimics(exposure, S, fit, modelFull, P, fwhm, psfPeak, noise)`

Find faint non-Gaia sources that are fakes absorbing a bright star's wing excess.

A faint source near a bright star must stand out from the residual's wing level at its distance from the star (including angular modes) and from the same spot rotated by 90, 180 and 270 degrees (four-fold diffraction structure).

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `S` | `dict` | Source state. |
| `fit` | `dict` | Joint fit. |
| `modelFull` | `numpy.ndarray` | Model including the bright-star corrections. |
| `P` | `dict` | Pipeline options (``ring_*``). |
| `fwhm` | `float` | PSF FWHM [px]. |
| `psfPeak` | `float` | Peak value of the unit-flux PSF. |
| `noise` | `float` | Typical pixel noise. |

| Returns | Type | Description |
| --- | --- | --- |
| `drop` | `list` [`int`] | Indices of the sources to remove. |
| `nChecked` | `int` | Number of sources checked. |

### pipeline.splitTest

`splitTest(exposure, psf, S, fit, modelFull, P, fwhm, psfPeak, noise, badBitmask, fitkw)`

Split unresolved pairs.

For each source with ``splitMinSnr <= peak S/N < brightCorrSigma``, a positive residual peak (matched-filter S/N >= ``splitPeakSnr``) 0.5-1.5 FWHM away and away from other sources gets a trial source.  All trials are fitted together; a trial is kept if the chi2 within 2 FWHM drops by at least ``splitDchi2`` and its S/N is at least ``minSnr``.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `S` | `dict` | Source state. |
| `fit` | `dict` | Joint fit. |
| `modelFull` | `numpy.ndarray` | Model including the bright-star corrections. |
| `P` | `dict` | Pipeline options (``split_*``). |
| `fwhm` | `float` | PSF FWHM [px]. |
| `psfPeak` | `float` | Peak value of the unit-flux PSF. |
| `noise` | `float` | Typical pixel noise. |
| `badBitmask` | `int` | Mask bits of pixels to ignore. |
| `fitkw` | `dict` | Keyword arguments of the joint fit. |

| Returns | Type | Description |
| --- | --- | --- |
| `S` | `dict` | Source state with the kept second sources. |
| `fit` | `dict` | Joint fit (the input object if nothing was kept). |
| `n_tried` | `int` | Number of trial sources. |
| `nKept` | `int` | Number kept. |

### pipeline.ellipTest

`ellipTest(exposure, psf, S, fit, modelFull, P, fwhm, psfPeak, noise, badBitmask, fitkw, psfSigma, verboseEllip=True, only=None, explain=False)`

Model badly fitted sources as small galaxies.

A badly fitted source (local chi2/px >= ``ellipChi2Min``, S/N between ``ellipMinSnr`` and ``brightCorrSigma``) and the points within ``ellipGroupFwhm`` FWHM (at most one Gaia source) are refitted with alternative models: 'one' (one extended source replaces all of them) and 'star+gal' (the Gaia source stays a point, one extended source replaces the others), each with a Gaussian or exponential elliptical profile convolved with the PSF.  The best model, penalised by 5 per extra fitted number, must beat the point sources by at least ``ellipDchi2``.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `S` | `dict` | Source state. |
| `fit` | `dict` | Joint fit. |
| `modelFull` | `numpy.ndarray` | Model including the bright-star corrections. |
| `P` | `dict` | Pipeline options (``ellip_*``). |
| `fwhm` | `float` | PSF FWHM [px]. |
| `psfPeak` | `float` | Peak value of the unit-flux PSF. |
| `noise` | `float` | Typical pixel noise. |
| `badBitmask` | `int` | Mask bits of pixels to ignore. |
| `fitkw` | `dict` | Keyword arguments of the joint fit. |
| `psfSigma` | `float` | PSF sigma [px]. |
| `verboseEllip` | `bool`, optional | Print every model tried. |
| `only` | `list` [`int`], optional | Only try these source indices. |
| `explain` | `bool`, optional | Print why sources are not tried, and tracebacks of failed fits. |

| Returns | Type | Description |
| --- | --- | --- |
| `S` | `dict` | Source state with the extended sources. |
| `fit` | `dict` | Joint fit. |
| `n_tried` | `int` | Number of sources tried. |
| `nKept` | `int` | Number of extended sources kept. |

### pipeline.tryElliptical

`tryElliptical(run, xy, chi2Min=None, dchi2=None)`

Try the elliptical-source test on one source of a finished run, printing every step.

Shows why the source is or is not tried, the fitted models, the chi2 change and any error. The run is not changed.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `xy` | `tuple` [`float`, `float`] | Position; the nearest source is tested. |
| `chi2Min` | `float`, optional | Override the trigger (``ellipChi2Min``). |
| `dchi2` | `float`, optional | Override the acceptance bar (``ellipDchi2``). |

| Returns | Type | Description |
| --- | --- | --- |
| `result` | `dict` or `None` | ``kept``, ``S`` and ``fit`` (with the extended source if it was kept); None if the run's last fit has no per-source stamps. |

### pipeline.finalPositionPolish

`finalPositionPolish(exposure, psf, S, fit, badBitmask, P, fitkw, fwhm, psfSigma)`

Re-optimise all source positions jointly.

Gauss-Newton steps (capped at ``finalPolishCap`` px, Gaia drift limit kept) and refit, repeated while the total chi2 drops.  Fixes dipoles from positions that went stale when neighbours were added or removed.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `S` | `dict` | Source state. |
| `fit` | `dict` | Joint fit of ``S``. |
| `badBitmask` | `int` | Mask bits of pixels to ignore. |
| `P` | `dict` | Pipeline options (``final_polish_*``). |
| `fitkw` | `dict` | Keyword arguments of the joint fit. |
| `fwhm` | `float` | PSF FWHM [px]. |
| `psfSigma` | `float` | PSF sigma [px]. |

| Returns | Type | Description |
| --- | --- | --- |
| `S` | `dict` | Source state with the polished positions. |
| `fit` | `dict` | Joint fit (the input object if nothing improved). |
| `info` | `dict` | ``rounds``, ``moved``, ``chi0`` and ``chi1`` (total chi2 before / after). |

### pipeline.runInnerAdaptive

`runInnerAdaptive(exposure, psf, badBitmask, S, P, outerIdx=0, thresholdStart=None, shapeNow=False, timer=None)`

Run one outer pass of the adaptive pipeline.

A warm fit of the current sources (Gaia starting positions or the previous pass), then detect -> add -> joint fit -> refine -> clean, repeated until the residual check says stop; then the halo, ring, position-polish, pair-split and galaxy tests.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Background-subtracted cutout. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `badBitmask` | `int` | Mask bits of pixels to ignore. |
| `S` | `dict` | Starting source state. |
| `P` | `dict` | Pipeline options. |
| `outerIdx` | `int`, optional | Outer-pass index. |
| `thresholdStart` | `float`, optional | First detection threshold [sigma]. |
| `shapeNow` | `bool`, optional | Apply the shape cut in this pass. |
| `timer` | `lsst.meas.crowd.fitLib.StepTimer`, optional | Timer collecting the step durations. |

| Returns | Type | Description |
| --- | --- | --- |
| `result` | `dict` | Compatible with `fitLib.runInnerLoopVerbose`'s result (``modelFull``, ``residualFull``, ``peaks``, ``amps``, histories, ``snaps``), plus ``S``, ``brightCorr``, ``adaptive`` (the loop decisions), ``isGaia`` and ``gaiaIndex``. |

### pipeline.selectPsfStars

`selectPsfStars(exposure, res, fwhm, psfPeak, noise, P, verbose=True)`

Select clean, bright, isolated, unsaturated point sources for the PSF stack.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `res` | `dict` | Result of `runInnerAdaptive`. |
| `fwhm` | `float` | PSF FWHM [px]. |
| `psfPeak` | `float` | Peak value of the unit-flux PSF. |
| `noise` | `float` | Typical pixel noise. |
| `P` | `dict` | Pipeline options (``psf_stack_*``). |
| `verbose` | `bool`, optional | Print the selection. |

| Returns | Type | Description |
| --- | --- | --- |
| `idx` | `numpy.ndarray` [`int`] | Indices of the selected stars, brightest first. |

### pipeline.refinePsfClean

`refinePsfClean(exposure, psf, res, badBitmask, P, fitkw, nCycles=3, verbose=True)`

Update the PSF from the residual stack of clean stars, cross-validated.

The clean stars are split into two halves A and B (alternating in brightness).  Each cycle builds a correction from A and judges it on B's pixels, and vice versa; only if the held-out chi2 improves by at least ``psfCvMinGain`` is the correction from all stars applied.  This stops the PSF absorbing the pixel noise of its own stars, which shows up as speckle around every bright star.  Corrections are smoothed (``psfSmoothPx``) and soft-thresholded (``psfSoftSigma``).

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psf` | `lsst.afw.detection.Psf` | Current PSF. |
| `res` | `dict` | Result of `runInnerAdaptive`. |
| `badBitmask` | `int` | Mask bits of pixels to ignore. |
| `P` | `dict` | Pipeline options. |
| `fitkw` | `dict` | Keyword arguments of the joint fit. |
| `nCycles` | `int`, optional | Maximum number of correction cycles. |
| `verbose` | `bool`, optional | Print each cycle. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | The updated PSF (the input PSF if no cycle was kept). |

### pipeline.refinePsfMedian

`refinePsfMedian(exposure, psf, res, badBitmask, P, fitkw, fit0=None, verbose=True)`

Correct the PSF core with the median residual/flux of the bright clean stars.

The median is robust to neighbours, cosmic rays and per-star noise, so it keeps the shape error shared by the bright stars.  Applied within ``psfCoreRadiusFwhm`` FWHM (cosine taper), soft-thresholded, smoothed and cross-validated on two halves of the stars.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psf` | `lsst.afw.detection.Psf` | Current PSF. |
| `res` | `dict` | Result of `runInnerAdaptive`. |
| `badBitmask` | `int` | Mask bits of pixels to ignore. |
| `P` | `dict` | Pipeline options (``psf_core_*``). |
| `fitkw` | `dict` | Keyword arguments of the joint fit. |
| `fit0` | `dict`, optional | Joint fit with ``psf`` (computed if not given). |
| `verbose` | `bool`, optional | Print each cycle. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | The corrected PSF (the input PSF if no correction was kept). |

### pipeline.smoothPsfKernel

`smoothPsfKernel(kernel, fwhm, coreFwhm=1.5, sigma=1.0, maxPattern=1.5)`

Return a hole-free, smooth version of a PSF kernel.

Inside ``coreFwhm`` FWHM the kernel is kept.  Outside it is rebuilt as ``radial profile(r) * exp(angular pattern)``: the profile is the median of each 1-px ring (so a few holed pixels cannot pull it down), forced to decrease outwards; the pattern is ``log(kernel / profile)`` median-filtered over 5x5 px (removes isolated holes and hot pixels), smoothed by ``sigma`` px and clipped to +-``maxPattern``.  Smooth structure (ellipticity, lobes, spikes wider than ~2 px) survives; pixel-scale defects do not.

| Parameter | Type | Description |
| --- | --- | --- |
| `kernel` | `numpy.ndarray` | PSF kernel (odd size, centred). |
| `fwhm` | `float` | PSF FWHM [px]. |
| `coreFwhm` | `float`, optional | Radius kept unchanged [FWHM]. |
| `sigma` | `float`, optional | Smoothing of the angular pattern [px]. |
| `maxPattern` | `float`, optional | Clip of the log angular pattern. |

| Returns | Type | Description |
| --- | --- | --- |
| `kernel` | `numpy.ndarray` | Smoothed, non-negative kernel summing to 1. |

### pipeline.regularizePsf

`regularizePsf(psf, point, P, verbose=True, label='')`

Return the PSF with its kernel passed through `smoothPsfKernel`.

| Parameter | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | PSF. |
| `point` | `lsst.geom.Point2D` | Position at which the kernel is taken. |
| `P` | `dict` | Pipeline options (``psfRegularize``, ``psf_reg_*``); no-op if regularisation is off. |
| `verbose` | `bool`, optional | Print how much flux moved. |
| `label` | `str`, optional | Suffix for the printed line. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | The regularised (constant) PSF, or the input PSF. |

### pipeline.extendPsfWingsKeep

`extendPsfWingsKeep(psf, point, size=81, blend=(12.0, 16.0), keepAngular=False, verbose=True)`

Extend a PSF with a Moffat tail, optionally keeping its angular structure.

Like `fitLib.extendPsfWings`, but with ``keepAngular`` only the kernel's azimuthal mean is replaced by the Moffat beyond the blend radius, so spikes and other non-circular structure are kept.

| Parameter | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | PSF. |
| `point` | `lsst.geom.Point2D` | Position at which the kernel is taken. |
| `size` | `int`, optional | Size of the extended kernel [px]. |
| `blend` | `tuple` [`float`, `float`], optional | Radial range of the cross-fade to the Moffat [px]. |
| `keepAngular` | `bool`, optional | Keep the kernel's angular pattern in the tail. |
| `verbose` | `bool`, optional | Print the fitted tail. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | The extended constant PSF (the input PSF if the fit fails). |

### pipeline.wingProfilePsf

`wingProfilePsf(exposure, res, psf, P, verbose=True)`

Return a PSF whose outer wings are measured from the stars.

Neighbour-subtracted radial profiles of unsaturated bright stars (peak S/N >= ``wingMinSnr``), divided by their flux, give the normalised profile where they have signal; saturated stars are scaled to it where both are measured and extend it outwards.  A Moffat fitted to the combined profile replaces the kernel beyond ``wingBlendFwhm``.  Degenerate fits are rejected.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `res` | `dict` | Result of `runInnerAdaptive`. |
| `psf` | `lsst.afw.detection.Psf` | Current PSF. |
| `P` | `dict` | Pipeline options (``wing_*``, ``psfWingSize``). |
| `verbose` | `bool`, optional | Print the fitted profile. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` or `None` | The new PSF, or None if there are too few stars or the fit is degenerate. |

### pipeline.spatialKernelPsf

`spatialKernelPsf(K0, Kx, Ky, xc, yc, scale)`

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

### pipeline.refinePsfSpatial

`refinePsfSpatial(exposure, psfOld, psfNew, res, badBitmask, P, fitkw, step=0.8, verbose=True)`

Add or update a linear spatial variation of the PSF, cross-validated.

Measured from the clean stars' residuals in ``res`` (fitted with ``psfOld``); the gradient from half A is judged on half B and vice versa, and kept only if the held-out chi2 drops by at least ``psfSpatialMinGain`` compared with the constant ``psfNew``.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psfOld` | `lsst.afw.detection.Psf` | PSF the residuals were fitted with. |
| `psfNew` | `lsst.afw.detection.Psf` | Constant PSF to add the variation to. |
| `res` | `dict` | Result of `runInnerAdaptive`. |
| `badBitmask` | `int` | Mask bits of pixels to ignore. |
| `P` | `dict` | Pipeline options (``psf_spatial_*``). |
| `fitkw` | `dict` | Keyword arguments of the joint fit. |
| `step` | `float`, optional | Fraction of the measured gradient applied. |
| `verbose` | `bool`, optional | Print the result. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | The varying PSF, or ``psfNew`` if the variation was not kept. |

### pipeline.fitBrightnessTerm

`fitBrightnessTerm(exposure, psf, res, badBitmask, P, fitkw, verbose=True)`

Measure or update the brightness-dependent PSF term (brighter-fatter).

From isolated, unsaturated stars with peak S/N >= ``bfMinSnr``: per 1-px ring, a weighted fit ``residual/flux = A(r) + B(r) * (log10 flux - L0)`` gives the slope B(r), smoothed, tapered and made flux-neutral.  Built from half of the stars and judged on the other half; kept only if the held-out chi2 improves.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psf` | `lsst.afw.detection.Psf` | Current PSF. |
| `res` | `dict` | Result of `runInnerAdaptive`. |
| `badBitmask` | `int` | Mask bits of pixels to ignore. |
| `P` | `dict` | Pipeline options (``bf_*``). |
| `fitkw` | `dict` | Keyword arguments of the joint fit. |
| `verbose` | `bool`, optional | Print the result. |

| Returns | Type | Description |
| --- | --- | --- |
| `moved` | `float` | Fraction of the flux moved per dex of brightness (0 if nothing changed). |

### pipeline.refitBackgroundLowfit

`refitBackgroundLowfit(exposure, mode='plane', statsMask=('BAD', 'SAT', 'EDGE', 'NO_DATA', 'DETECTED', 'CR', 'INTRP'), clip=3.0, nIter=5, verbose=True)`

Fit a constant or tilted-plane background to the unmasked pixels and subtract it in place.

The robust alternative to `fitLib.refitBackground` for small cutouts of an image that is already background-subtracted: a spline / Chebyshev surface on 32 px bins is barely constrained where most pixels are masked and can swing by hundreds of counts in a corner.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Exposure to background-subtract in place (sources flagged DETECTED). |
| `mode` | `str`, optional | 'plane' (a + b x + c y), 'constant' or 'none'. |
| `statsMask` | `tuple` [`str`], optional | Mask planes whose pixels are ignored. |
| `clip` | `float`, optional | Sigma clipping threshold. |
| `nIter` | `int`, optional | Clipping iterations. |
| `verbose` | `bool`, optional | Print the fitted background. |

| Returns | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The background-subtracted exposure (same object). |
| `background` | `lsst.afw.image.ImageF` | The subtracted background image. |

Raises: `ValueError` (Raised if ``mode`` is not 'plane', 'constant' or 'none'.)

### pipeline.printTiming

`printTiming(timer=None)`

Print the time spent per pipeline step.

| Parameter | Type | Description |
| --- | --- | --- |
| `timer` | `lsst.meas.crowd.fitLib.StepTimer`, optional | Timer (default: fitLib's ``TIMER``). |

### pipeline.runAdaptivePipeline

`runAdaptivePipeline(exposureCrop, psf, badBitmask, gaiaXy=None, gaiaRows=None, starCat=None, bgImage=None, viz=False, gaiaGmag=None, **cfg)`

Run the adaptive pipeline on a whole cutout, starting from Gaia positions.

The replacement for `fitLib.runCutoutPipeline`: outer passes of `runInnerAdaptive`, with background refit and PSF refinement in between, until the source list and the PSF stop changing.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposureCrop` | `lsst.afw.image.Exposure` | Cutout before background subtraction. |
| `psf` | `lsst.afw.detection.Psf` | Initial PSF. |
| `badBitmask` | `int` | Mask bits of pixels to ignore. |
| `gaiaXy` | `numpy.ndarray`, optional | Gaia starting positions (detector px), shape (N, 2); None or empty = no Gaia. |
| `gaiaRows` | `numpy.ndarray`, optional | Gaia catalogue row of each starting position. |
| `starCat` | `lsst.afw.table.SourceCatalog`, optional | Detector-level PSF stars. |
| `bgImage` | `lsst.afw.image.Image`, optional | Detector-level background, subtracted in the first pass. |
| `viz` | `bool`, optional | Draw the per-pass overview. |
| `gaiaGmag` | `numpy.ndarray`, optional | Gaia G magnitude of each starting position. |
| `**cfg` | – | Configuration: `fitLib.runCutoutPipeline` keys and ``PIPELINE_DEFAULTS`` keys. |

| Returns | Type | Description |
| --- | --- | --- |
| `outerHistory` | `list` [`dict`] | One dict per outer pass (as `fitLib.runCutoutPipeline`'s, plus ``psfChange``, ``modelErr``, ``bf`` and ``result['adaptive']``). |

### pipeline.dropOffDetectorGaia

`dropOffDetectorGaia(exposure, xy, rows)`

Drop Gaia starting positions that lie beyond the detector edge.

A star outside the cutout whose nearest cutout pixel is flagged EDGE or NO_DATA lies off the detector: no pixel constrains its flux, so the joint fit can give it an absurd value whose wings then spoil the edge of the image.  Stars outside a cutout edge that is not a detector edge are kept (their wings fall on real pixels and their neighbours need them).

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The cutout. |
| `xy` | `numpy.ndarray` | Gaia positions (detector px), shape (N, 2). |
| `rows` | `numpy.ndarray` | Their rows in the Gaia table. |

| Returns | Type | Description |
| --- | --- | --- |
| `xy, rows` | `numpy.ndarray` | The positions and rows that are kept. |

### pipeline.runCutoutGaia

`runCutoutGaia(rawCrop, psf, badBitmask, gaiaTbl=None, config=None, interiorBbox=None, starCat=None, bgImage=None, viz=False, name='cutout', gmagMax=None, **overrides)`

Run the adaptive pipeline on one cropped exposure, starting from Gaia positions.

| Parameter | Type | Description |
| --- | --- | --- |
| `rawCrop` | `lsst.afw.image.Exposure` | Cutout before background subtraction. |
| `psf` | `lsst.afw.detection.Psf` | Initial PSF. |
| `badBitmask` | `int` | Mask bits of pixels to ignore. |
| `gaiaTbl` | `astropy.table.Table`, optional | Gaia catalogue (see `lsst.meas.crowd.gaiaCatalog.loadGaiaCatalog`); None = no Gaia. |
| `config` | `dict`, optional | Pipeline configuration. |
| `interiorBbox` | `lsst.geom.Box2I`, optional | The requested box without padding (default: the whole cutout). |
| `starCat` | `lsst.afw.table.SourceCatalog`, optional | Detector-level PSF stars. |
| `bgImage` | `lsst.afw.image.Image`, optional | Detector-level background. |
| `viz` | `bool`, optional | Draw the per-pass overview. |
| `name` | `str`, optional | Cutout name. |
| `gmagMax` | `float`, optional | Only start from Gaia stars brighter than this G magnitude. |
| `**overrides` | – | Configuration overrides. |

| Returns | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | The result; its catalogue also has ``fromGaia``, ``gaiaSourceId``, ``gaiaG``, ``gaiaDx`` and ``gaiaDy``. |

### pipeline.checkBrightNeighbours

`checkBrightNeighbours(run, brightSnr=300.0, withinFwhm=4.0, fluxRatio=5.0, nControls=7, plot=True, maxShow=16, verbose=True, similarMin=3)`

Test whether sources on the edges of bright stars are real.

For every source within ``withinFwhm`` FWHM of a star at least ``fluxRatio`` times brighter (with peak S/N >= ``brightSnr``), its PSF-weighted flux is compared with the residual's flux at control positions at the same distance from the bright star, rotated in 45-degree steps. A mismatch pattern of the bright star (rings, spikes, lobes) puts similar flux at several angles; a real star is significantly brighter than its controls.  Verdicts: 'Gaia' (real by construction), 'real' (significance >= 5 and at most one similar angle), 'probable' (significance >= 3, or two similar angles), 'doubtful' (significance < 3, or similar flux at >= ``similarMin`` angles).

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `brightSnr` | `float`, optional | Minimum peak S/N of a bright star. |
| `withinFwhm` | `float`, optional | Neighbour radius [FWHM]. |
| `fluxRatio` | `float`, optional | Minimum brightness ratio of the bright star to its neighbour. |
| `nControls` | `int`, optional | Number of control angles. |
| `plot` | `bool`, optional | Show thumbnails of the non-Gaia neighbours. |
| `maxShow` | `int`, optional | Maximum number of thumbnails. |
| `verbose` | `bool`, optional | Print the table. |
| `similarMin` | `int`, optional | Number of similar control angles that makes a source doubtful. |

| Returns | Type | Description |
| --- | --- | --- |
| `table` | `astropy.table.Table` | One row per neighbour: position, flux, significance, number of similar angles, verdict. |

### pipeline.diagnoseResidual

`diagnoseResidual(run, brightSnr=300.0, diskFwhm=3.0, plot=True)`

Decide whether what is left around bright stars is noise or model error.

For empty sky, bright-star disks (within ``diskFwhm`` FWHM of stars with peak S/N >= ``brightSnr``) and other source pixels: ``sigmaZ`` (robust std of residual/sigma; 1 = matches the variance plane) and ``coherence`` (robust std of the 3x3-averaged residual/sigma, x3, over sigmaZ, squared; white noise gives 1).  Also the pixel-to-pixel scatter against model flux, which measures missing photon noise.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `brightSnr` | `float`, optional | Minimum peak S/N of a bright star. |
| `diskFwhm` | `float`, optional | Radius of the bright-star disks [FWHM]. |
| `plot` | `bool`, optional | Show residual/sigma, its 3x3 average and its matched-filter S/N. |

| Returns | Type | Description |
| --- | --- | --- |
| `result` | `dict` | Per region ``sigmaZ``, ``coherence`` and pixel count; ``highfreq`` (pixel-to-pixel scatter vs model flux, with ``k``) and ``noiseShare`` (fraction of the bright-disk excess explained by missing photon noise). |

### pipeline.modeledNear

`modeledNear(modeled, bad, fwhm)`

Return the pixels closer to a modelled (near-saturated) core than to any other masked pixel.

Residual blobs there are not grown into the excluded area by `unmodeledNearMask`.

| Parameter | Type | Description |
| --- | --- | --- |
| `modeled` | `numpy.ndarray` [`bool`] | Masked pixels whose star is in the model. |
| `bad` | `numpy.ndarray` [`bool`] | All masked pixels. |
| `fwhm` | `float` | PSF FWHM [px] (unused; kept for a uniform signature). |

| Returns | Type | Description |
| --- | --- | --- |
| `near` | `numpy.ndarray` [`bool`] | Pixels nearer to a modelled core than to other masked pixels. |

### pipeline.unmodeledNearMask

`unmodeledNearMask(resid, var, bad, fwhm, padFwhm=1.5, growSnr=2.5, reachFwhm=4.0, modeled=None)`

Find the pixels whose residual cannot be judged because a source there is partly masked.

Everything within ``padFwhm`` FWHM of a bad pixel, plus any positive residual blob (matched-filter S/N > ``growSnr``) touching that pad, out to ``reachFwhm`` FWHM from the mask.  Masked cores whose star is in the model (``modeled``, near-saturated cores fitted from their wings) get no pad: the residual next to them is real information.

| Parameter | Type | Description |
| --- | --- | --- |
| `resid` | `numpy.ndarray` | Residual image. |
| `var` | `numpy.ndarray` | Variance image. |
| `bad` | `numpy.ndarray` [`bool`] | Bad pixels. |
| `fwhm` | `float` | PSF FWHM [px]. |
| `padFwhm` | `float`, optional | Pad around bad pixels [FWHM]. |
| `growSnr` | `float`, optional | Matched-filter S/N of residual blobs grown into the excluded area. |
| `reachFwhm` | `float`, optional | Maximum distance of the excluded area from the mask [FWHM]. |
| `modeled` | `numpy.ndarray` [`bool`], optional | Masked pixels whose star is in the model. |

| Returns | Type | Description |
| --- | --- | --- |
| `excluded` | `numpy.ndarray` [`bool`] | Pixels to exclude (not including the bad pixels themselves). |

### pipeline.runSessionGaia

`runSessionGaia(session, config=None, names=None, outDir=None, viz=False, plotFinal=True, gmagMax=None, nJobs=1, plotKwargs=None, **overrides)`

Run the adaptive pipeline on every (or the named) cutout of a session.

Calls `lsst.meas.crowd.gaiaCatalog.attachGaia` first if the session has no Gaia catalogue.

| Parameter | Type | Description |
| --- | --- | --- |
| `session` | `lsst.meas.crowd.session.Session` | The session. |
| `config` | `dict`, optional | Pipeline configuration. |
| `names` | `list` [`str`], optional | Cutouts to run (default: all). |
| `outDir` | `str`, optional | Save each run under ``<outDir>/visit<V>_det<D>/<name>/``; needed with ``nJobs > 1``, since results that cannot be sent back are only on disk. |
| `viz` | `bool`, optional | Draw the per-pass overview. |
| `plotFinal` | `bool`, optional | Plot the final comparison of each cutout. |
| `gmagMax` | `float`, optional | Only start from Gaia stars brighter than this G magnitude. |
| `nJobs` | `int`, optional | Fit cutouts in this many parallel worker processes (Linux fork; one core each). |
| `plotKwargs` | `dict`, optional | Keyword arguments of `lsst.meas.crowd.plotting.plotFinalComparison`. |
| `**overrides` | – | Configuration overrides (e.g. ``maxOuter=3``). |

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

Sessions, configuration and saving of runs.

| Name | Summary |
| --- | --- |
| [`saveConfig`](#sessionsaveconfig) | Write a config dict to a JSON file. |
| [`makeConfig`](#sessionmakeconfig) | Build a pipeline configuration (keyword arguments of `fitLib.runCutoutPipeline`). |
| [`checkConfig`](#sessioncheckconfig) | Check that every key of a configuration is a `fitLib.runCutoutPipeline` option. |
| [`showConfig`](#sessionshowconfig) | Print a configuration. |
| [`parseCutout`](#sessionparsecutout) | Parse a cutout specification. |
| [`Session`](#sessionsession) (class) | Everything needed to fit cutouts of one visit and detector. |
| [`buildDetectorPsf`](#sessionbuilddetectorpsf) | Set up the detector-level PSF. |
| [`prepareSession`](#sessionpreparesession) | Load a visit image, cut out the regions and set up the detector-level PSF. |
| [`queryGaia`](#sessionquerygaia) | Query Gaia DR3 at the ESA archive for the sources on a box. |
| [`CutoutRun`](#sessioncutoutrun) (class) | Result of running the pipeline on one cutout. |
| [`buildCatalog`](#sessionbuildcatalog) | Build the source catalogue of a finished run. |
| [`runCutout`](#sessionruncutout) | Run the pipeline on one cropped exposure. |
| [`runSession`](#sessionrunsession) | Run the pipeline on every (or the named) cutout of a session. |
| [`writePhotometry`](#sessionwritephotometry) | Write a run's photometry catalogue as photometry.fits and photometry.csv. |
| [`saveRun`](#sessionsaverun) | Write a run's catalogues, images, PSF kernels and summary to a directory. |
| [`robustStats`](#sessionrobuststats) | Return the median and IQR-based robust sigma of finite (and valid) values. |
| [`matplotlibListed`](#sessionmatplotliblisted) | Return a one-colour matplotlib colour map. |
| [`plotFinalComparisonRobust`](#sessionplotfinalcomparisonrobust) | Plot original \| final model \| residual \| residual/sigma with robust statistics. |
| [`plotFinalComparisonFromFits`](#sessionplotfinalcomparisonfromfits) | Rebuild the final-comparison figure from files written by `saveRun`. |

### session.saveConfig

`saveConfig(cfg, path)`

Write a config dict to a JSON file.

Tuples become lists; `makeConfig` turns them back into tuples when the file is loaded.

| Parameter | Type | Description |
| --- | --- | --- |
| `cfg` | `dict` | Pipeline configuration. |
| `path` | `str` | Output file. |

| Returns | Type | Description |
| --- | --- | --- |
| `path` | `str` | The file written. |

### session.makeConfig

`makeConfig(base=None, deep=False, nOuter=None, nInner=None, saveTo=None, **overrides)`

Build a pipeline configuration (keyword arguments of `fitLib.runCutoutPipeline`).

| Parameter | Type | Description |
| --- | --- | --- |
| `base` | `str` or `dict`, optional | A JSON file or a dict, applied on top of ``DEFAULT_CONFIG``. |
| `deep` | `bool`, optional | Merge nested dicts key by key (default: replace them). |
| `nOuter` | `int`, optional | Shortcut for ``nOuterIterations``. |
| `nInner` | `int`, optional | Shortcut for ``nInnerIterations`` and ``innerKwargs['nIterations']``. |
| `saveTo` | `str`, optional | Also write the configuration to this JSON file. |
| `**overrides` | – | Any `fitLib.runCutoutPipeline` keyword, applied last. |

| Returns | Type | Description |
| --- | --- | --- |
| `cfg` | `dict` | The configuration. |

### session.checkConfig

`checkConfig(cfg)`

Check that every key of a configuration is a `fitLib.runCutoutPipeline` option.

| Parameter | Type | Description |
| --- | --- | --- |
| `cfg` | `dict` | Pipeline configuration. |

Raises: `ValueError` (Raised if ``cfg`` has keys the pipeline does not accept.)

### session.showConfig

`showConfig(cfg, diffOnly=False)`

Print a configuration.

| Parameter | Type | Description |
| --- | --- | --- |
| `cfg` | `dict` | Pipeline configuration. |
| `diffOnly` | `bool`, optional | Only print keys that differ from ``DEFAULT_CONFIG``. |

### session.parseCutout

`parseCutout(spec)`

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

Class.  `Session(visit, detector, badBitmask, psf, starCat, bgImage, jobs=<factory>, wcs=None, psfUsedXy=None)`

Everything needed to fit cutouts of one visit and detector.

| Parameter | Type | Description |
| --- | --- | --- |
| `visit` | `int` | Visit id. |
| `detector` | `int` | Detector id. |
| `badBitmask` | `int` | Mask bits of pixels the fit ignores. |
| `psf` | `lsst.afw.detection.Psf` | Detector-level PSF. |
| `starCat` | `lsst.afw.table.SourceCatalog` or `None` | PSF stars of the detector-level PSF determination. |
| `bgImage` | `lsst.afw.image.Image` or `None` | Detector-level background. |
| `jobs` | `dict` [`str`, `dict`] | Per cutout: ``interior`` (`lsst.geom.Box2I`, the requested box) and ``raw`` (the padded cutout before background subtraction). |
| `wcs` | `lsst.afw.geom.SkyWcs`, optional | The visit image's WCS. |
| `psfUsedXy` | `numpy.ndarray`, optional | Positions (detector px, shape (N, 2)) of the stars the detector-level PSF was built from (``psfSource='refit'`` only). |

### session.buildDetectorPsf

`buildDetectorPsf(visitImage, psfSource='refit', returnUsed=False)`

Set up the detector-level PSF.

With ``psfSource='visit'`` the visit image's own PSF is used.  With 'refit' the background is subtracted, sources are detected and measured, stars selected and a PCA PSF determined (the visit image is modified in place).

| Parameter | Type | Description |
| --- | --- | --- |
| `visitImage` | `lsst.afw.image.Exposure` | The detector's visit image. |
| `psfSource` | `str`, optional | 'refit' or 'visit'. |
| `returnUsed` | `bool`, optional | Also return the positions of the stars the PCA PSF was built from. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | Detector-level PSF. |
| `starCat` | `lsst.afw.table.SourceCatalog` or `None` | Selected PSF stars (None for 'visit'). |
| `bgImage` | `lsst.afw.image.Image` or `None` | Subtracted background (None for 'visit'). |
| `usedXy` | `numpy.ndarray` or `None` | Positions of the stars used, shape (N, 2) (only with ``returnUsed``; None for 'visit'). |

### session.prepareSession

`prepareSession(visit, detector, cutouts, pad=0, psfSource='refit', instrument='LSSTCam', repo='dp2_prep', collection='LSSTCam/runs/DRP/DP2', visitImageCollections='LSSTCam/runs/*', visitImage=None)`

Load a visit image, cut out the regions and set up the detector-level PSF.

The cutouts are taken before the detector-wide background subtraction.

| Parameter | Type | Description |
| --- | --- | --- |
| `visit` | `int` | Visit id. |
| `detector` | `int` | Detector id. |
| `cutouts` | `list` | Cutout specifications (see `parseCutout`); absolute detector pixels. |
| `pad` | `int`, optional | Fit a region grown by this many px; the interior keeps the requested box. |
| `psfSource` | `str`, optional | 'refit' or 'visit' (see `buildDetectorPsf`). |
| `instrument` | `str`, optional | Instrument name. |
| `repo` | `str`, optional | Butler repository path or alias. |
| `collection` | `str`, optional | Default Butler collection. |
| `visitImageCollections` | `str`, optional | Collections searched for the visit image. |
| `visitImage` | `lsst.afw.image.Exposure`, optional | An already-loaded visit image (skips the Butler). |

| Returns | Type | Description |
| --- | --- | --- |
| `session` | `Session` | The session. |

Raises: `RuntimeError` (Raised if no visit image is found.)

### session.queryGaia

`queryGaia(wcs, bbox, tapUrl='https://gea.esac.esa.int/tap-server/tap')`

Query Gaia DR3 at the ESA archive for the sources on a box.

| Parameter | Type | Description |
| --- | --- | --- |
| `wcs` | `lsst.afw.geom.SkyWcs` | WCS of the image. |
| `bbox` | `lsst.geom.Box2I` | Pixel box. |
| `tapUrl` | `str`, optional | TAP service URL. |

| Returns | Type | Description |
| --- | --- | --- |
| `gaia` | `astropy.table.Table` | Gaia sources inside the box, with their pixel positions ``x``, ``y``. |

### session.CutoutRun

Class.  `CutoutRun(name, outerHistory, catalog, summary, rawCrop, interiorBbox, badBitmask, config, gaia=None, photometry=None)`

Result of running the pipeline on one cutout.

| Parameter | Type | Description |
| --- | --- | --- |
| `name` | `str` | Cutout name. |
| `outerHistory` | `list` [`dict`] | One dict per outer pass (``psf``, ``exposure``, ``result``, ...). |
| `catalog` | `astropy.table.Table` | Source catalogue (see `buildCatalog`). |
| `summary` | `dict` | Run summary (timings, counts, residual statistics). |
| `rawCrop` | `lsst.afw.image.Exposure` | The input cutout before background subtraction. |
| `interiorBbox` | `lsst.geom.Box2I` | The requested box (without padding). |
| `badBitmask` | `int` | Mask bits of pixels the fit ignored. |
| `config` | `dict` | Configuration used. |
| `gaia` | `astropy.table.Table`, optional | The Gaia catalogue whose positions started the run (used for the Gaia comparison columns). |
| `photometry` | `astropy.table.Table`, optional | The photometry catalogue (see `lsst.meas.crowd.photometry.measurePhotometry`), made at the end of the run and written by `saveRun` as ``photometry.fits``. |

### session.buildCatalog

`buildCatalog(outerHistory, interiorBbox, gaiaTbl=None, matchRadius=2.0)`

Build the source catalogue of a finished run.

| Parameter | Type | Description |
| --- | --- | --- |
| `outerHistory` | `list` [`dict`] | Outer-loop history. |
| `interiorBbox` | `lsst.geom.Box2I` | The requested box (``inInterior`` column). |
| `gaiaTbl` | `astropy.table.Table`, optional | Gaia sources (with ``x``, ``y``) to match. |
| `matchRadius` | `float`, optional | Gaia match radius [px]. |

| Returns | Type | Description |
| --- | --- | --- |
| `cat` | `astropy.table.Table` | Positions (pixel and sky), PSF flux, AB magnitude (assuming nJy), source metrics, ``inInterior`` and, with ``gaiaTbl``, the nearest Gaia source. |

### session.runCutout

`runCutout(rawCrop, psf, badBitmask, config=None, interiorBbox=None, starCat=None, bgImage=None, viz=False, diagnostics=True, gaiaTbl=None, matchRadius=2.0, name='cutout', **overrides)`

Run the pipeline on one cropped exposure.

| Parameter | Type | Description |
| --- | --- | --- |
| `rawCrop` | `lsst.afw.image.Exposure` | The cutout (before background subtraction). |
| `psf` | `lsst.afw.detection.Psf` | Initial PSF. |
| `badBitmask` | `int` | Mask bits of pixels to ignore. |
| `config` | `dict`, optional | Pipeline configuration (default ``DEFAULT_CONFIG``). |
| `interiorBbox` | `lsst.geom.Box2I`, optional | The requested box (default: the whole cutout). |
| `starCat` | `lsst.afw.table.SourceCatalog`, optional | PSF stars forced into every PSF determination. |
| `bgImage` | `lsst.afw.image.Image`, optional | Detector background subtracted before the first pass. |
| `viz` | `bool`, optional | Draw the per-iteration diagnostics. |
| `diagnostics` | `bool`, optional | Run the end-of-run diagnostics. |
| `gaiaTbl` | `astropy.table.Table`, optional | Gaia sources to match in the catalogue. |
| `matchRadius` | `float`, optional | Gaia match radius [px]. |
| `name` | `str`, optional | Cutout name. |
| `**overrides` | – | Configuration overrides. |

| Returns | Type | Description |
| --- | --- | --- |
| `run` | `CutoutRun` | The result. |

### session.runSession

`runSession(session, config=None, names=None, outDir=None, viz=False, gaia=False, matchRadius=2.0, plotFinal=True, **overrides)`

Run the pipeline on every (or the named) cutout of a session.

A failure in one cutout is printed and skipped.

| Parameter | Type | Description |
| --- | --- | --- |
| `session` | `Session` | The session. |
| `config` | `dict`, optional | Pipeline configuration. |
| `names` | `list` [`str`], optional | Cutouts to run (default: all). |
| `outDir` | `str`, optional | Save each run under ``<outDir>/visit<V>_det<D>/<name>/`` (see `saveRun`). |
| `viz` | `bool`, optional | Draw the per-iteration diagnostics. |
| `gaia` | `bool`, optional | Match the catalogues to Gaia (queried at the ESA archive). |
| `matchRadius` | `float`, optional | Gaia match radius [px]. |
| `plotFinal` | `bool`, optional | Plot the final comparison of each cutout. |
| `**overrides` | – | Configuration overrides. |

| Returns | Type | Description |
| --- | --- | --- |
| `runs` | `dict` [`str`, `CutoutRun`] | The results by cutout name. |

### session.writePhotometry

`writePhotometry(run, outdir, verbose=True)`

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

### session.saveRun

`saveRun(run, outdir)`

Write a run's catalogues, images, PSF kernels and summary to a directory.

Writes photometry.fits / photometry.csv (the photometry catalogue with the Gaia comparison, computed here if the run has none), catalog.csv / catalog.fits, images.fits, psf_kernels.fits and summary.json, the layout `plotFinalComparisonFromFits` reads.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `CutoutRun` | The run. |
| `outdir` | `str` | Output directory (created if needed). |

### session.robustStats

`robustStats(values, valid=None)`

Return the median and IQR-based robust sigma of finite (and valid) values.

| Parameter | Type | Description |
| --- | --- | --- |
| `values` | `numpy.ndarray` | Values. |
| `valid` | `numpy.ndarray` [`bool`], optional | Values to use. |

| Returns | Type | Description |
| --- | --- | --- |
| `stats` | `dict` | ``median``, ``sigmaIqr`` ((Q75 - Q25) / 1.349), ``q25``, ``q75`` and ``n``. |

### session.matplotlibListed

`matplotlibListed(color)`

Return a one-colour matplotlib colour map.

| Parameter | Type | Description |
| --- | --- | --- |
| `color` | `str` | The colour. |

| Returns | Type | Description |
| --- | --- | --- |
| `cmap` | `matplotlib.colors.ListedColormap` | The colour map. |

### session.plotFinalComparisonRobust

`plotFinalComparisonRobust(runOrHistory, exposureCrop=None, badBitmask=None, edgeMarginFwhm=None, view='interior', statsInView=True, showPeaks=True, figsize=(26, 7), verbose=True)`

Plot original | final model | residual | residual/sigma with robust statistics.

Bad pixels are black in the residual panels; median and IQR-robust sigma (good pixels only) are shown in the titles.

| Parameter | Type | Description |
| --- | --- | --- |
| `runOrHistory` | `CutoutRun` or `list` [`dict`] | A run (preferred) or its outer-loop history. |
| `exposureCrop` | `lsst.afw.image.Exposure`, optional | Raw crop (only needed with a bare history). |
| `badBitmask` | `int`, optional | Mask bits of bad pixels (needed with a bare history). |
| `edgeMarginFwhm` | `float`, optional | Draw the edge drop-zone of this width [FWHM]. |
| `view` | `str` or `tuple`, optional | 'interior' (the run's interior box), None (the whole fitted array) or (xmin, xmax, ymin, ymax) in detector pixels. |
| `statsInView` | `bool`, optional | Restrict the statistics to the displayed view. |
| `showPeaks` | `bool`, optional | Mark the sources. |
| `figsize` | `tuple` [`float`, `float`], optional | Figure size [inches]. |
| `verbose` | `bool`, optional | Print the statistics. |

| Returns | Type | Description |
| --- | --- | --- |
| `stats` | `dict` | ``residual`` and ``residualOverSigma``: median, sigmaIqr, q25, q75, n. |

Raises: `ValueError` (Raised if a bare history is passed without ``badBitmask``.)

### session.plotFinalComparisonFromFits

`plotFinalComparisonFromFits(d, view=None, out=None, statsInView=True, showSources=True, show=True, dpi=110, verbose=True)`

Rebuild the final-comparison figure from files written by `saveRun`.

No LSST stack and no re-fit needed.

| Parameter | Type | Description |
| --- | --- | --- |
| `d` | `str` | Directory written by `saveRun`. |
| `view` | `tuple` [`float`, `float`, `float`, `float`], optional | (xmin, xmax, ymin, ymax) in detector pixels (default: the interior box from summary.json, else the whole array). |
| `out` | `str`, optional | Also save the figure to this PNG. |
| `statsInView` | `bool`, optional | Restrict the statistics to the view. |
| `showSources` | `bool`, optional | Mark the catalogue's sources. |
| `show` | `bool`, optional | Display the figure. |
| `dpi` | `int`, optional | Resolution of the PNG. |
| `verbose` | `bool`, optional | Print the statistics. |

| Returns | Type | Description |
| --- | --- | --- |
| `stats` | `dict` | See `plotFinalComparisonRobust`. |

## gaiaCatalog

`lsst.meas.crowd.gaiaCatalog`: the Gaia catalogue and the starting positions of each cutout.

The Gaia DR3 catalogue of a visit and the Gaia starting positions of each cutout.

| Name | Summary |
| --- | --- |
| [`visitEpoch`](#gaiacatalogvisitepoch) | Return the observation epoch of an exposure. |
| [`loadGaiaCatalog`](#gaiacatalogloadgaiacatalog) | Load the Gaia DR3 stars overlapping a box, with positions at the observation epoch. |
| [`attachGaia`](#gaiacatalogattachgaia) | Load Gaia once for a session and store it as ``session.gaia``. |
| [`selectGaiaSources`](#gaiacatalogselectgaiasources) | Select the Gaia stars whose positions start the fit of a cutout. |

### gaiaCatalog.visitEpoch

`visitEpoch(exposure)`

Return the observation epoch of an exposure.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Exposure with a VisitInfo. |

| Returns | Type | Description |
| --- | --- | --- |
| `epoch` | `astropy.time.Time` or `None` | The visit date, or None if it is missing or implausible. |

### gaiaCatalog.loadGaiaCatalog

`loadGaiaCatalog(wcs, bbox, epoch=None, source='auto', butler=None, repo='dp2_prep', refcat='gaia_dr3_20230707', refcatCollections='refcats', tapUrl='https://gea.esac.esa.int/tap-server/tap', cachePath=None, pad=0, verbose=True)`

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
| `refcatCollections` | `str`, optional | Collections of the reference catalogue. |
| `tapUrl` | `str`, optional | TAP service URL. |
| `cachePath` | `str`, optional | .ecsv file: read if it exists, written otherwise. |
| `pad` | `int`, optional | Grow the box by this many px. |
| `verbose` | `bool`, optional | Print what was loaded. |

| Returns | Type | Description |
| --- | --- | --- |
| `gaia` | `astropy.table.Table` | ``source_id``, ``ra``, ``dec``, ``phot_g_mean_mag``, pixel positions ``x``, ``y`` and further Gaia columns; ``meta['epochJyear']`` is the epoch used. |

Raises: `RuntimeError` (Raised if Gaia could not be loaded from any source.)

### gaiaCatalog.attachGaia

`attachGaia(session, visitImage=None, pad=64, cacheDir='gaia_cache', **kw)`

Load Gaia once for a session and store it as ``session.gaia``.

Covers the whole visit image or, without one, the union of the session's cutouts; the catalogue is cached on disk.

| Parameter | Type | Description |
| --- | --- | --- |
| `session` | `lsst.meas.crowd.session.Session` | The session (``session.gaia`` is set). |
| `visitImage` | `lsst.afw.image.Exposure`, optional | The full visit image. |
| `pad` | `int`, optional | Grow the box by this many px. |
| `cacheDir` | `str`, optional | Cache directory (None = no cache). |
| `**kw` | – | Keyword arguments of `loadGaiaCatalog` (e.g. ``repo``). |

| Returns | Type | Description |
| --- | --- | --- |
| `gaia` | `astropy.table.Table` | The Gaia catalogue (see `loadGaiaCatalog`). |

### gaiaCatalog.selectGaiaSources

`selectGaiaSources(gaiaTbl, bbox, fwhm, padFwhm=1.5, gmagMax=None, minSepPx=0.5)`

Select the Gaia stars whose positions start the fit of a cutout.

Stars inside the box grown by ``padFwhm * fwhm``, brightest first, with near-duplicates (closer than ``minSepPx``) removed.

| Parameter | Type | Description |
| --- | --- | --- |
| `gaiaTbl` | `astropy.table.Table` or `None` | Gaia catalogue with ``x``, ``y`` (and ``phot_g_mean_mag``). |
| `bbox` | `lsst.geom.Box2I` | Cutout box. |
| `fwhm` | `float` | PSF FWHM [px]. |
| `padFwhm` | `float`, optional | Growth of the box [FWHM]. |
| `gmagMax` | `float`, optional | Only stars brighter than this G magnitude. |
| `minSepPx` | `float`, optional | Drop fainter stars within this many px of a brighter one. |

| Returns | Type | Description |
| --- | --- | --- |
| `xy` | `numpy.ndarray` | Starting positions, shape (N, 2). |
| `rows` | `numpy.ndarray` [`int`] | Their rows in ``gaiaTbl``. |

## photometry

`lsst.meas.crowd.photometry`: PSF and aperture photometry, errors, the Gaia comparison.

Magnitudes of the fitted sources.

| Name | Summary |
| --- | --- |
| [`fluxToMag`](#photometryfluxtomag) | Convert fluxes (and errors) to magnitudes. |
| [`apertureFlux`](#photometryapertureflux) | Sum the flux in a circular aperture. |
| [`apertureCorrection`](#photometryaperturecorrection) | Return the fraction of a PSF kernel's flux inside a radius. |
| [`fisherFluxErrors`](#photometryfisherfluxerrors) | Return the 1-sigma errors of jointly fitted amplitudes. |
| [`psfFluxErrors`](#photometrypsffluxerrors) | Return the Fisher errors of the PSF fluxes of a run's final joint fit. |
| [`zeroPointFromRun`](#photometryzeropointfromrun) | Return the magnitude zero point of a run's image units. |
| [`compareWithGaia`](#photometrycomparewithgaia) | Add Gaia comparison columns to a photometry table. |
| [`fitsMeta`](#photometryfitsmeta) | Return catalogue metadata with FITS-compatible keywords and values. |
| [`measurePhotometry`](#photometrymeasurephotometry) | Measure PSF and aperture magnitudes of a run's sources. |

### photometry.fluxToMag

`fluxToMag(flux, fluxErr=None, zeroPoint=31.4)`

Convert fluxes (and errors) to magnitudes.

| Parameter | Type | Description |
| --- | --- | --- |
| `flux` | array-like | Fluxes, in units where ``zeroPoint`` applies. |
| `fluxErr` | array-like, optional | Flux errors. |
| `zeroPoint` | `float`, optional | Magnitude of unit flux (default: AB for nJy). |

| Returns | Type | Description |
| --- | --- | --- |
| `mag` | `numpy.ndarray` | Magnitudes (NaN for non-positive fluxes). |
| `mag_err` | `numpy.ndarray` | Magnitude errors (only if ``fluxErr`` is given). |

### photometry.apertureFlux

`apertureFlux(image, variance, good, x, y, radius, x0=0, y0=0)`

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
| `fluxErr` | `float` | Its error. |
| `coverage` | `float` | Fraction of the aperture area on good pixels. |

### photometry.apertureCorrection

`apertureCorrection(kernel, radius)`

Return the fraction of a PSF kernel's flux inside a radius.

| Parameter | Type | Description |
| --- | --- | --- |
| `kernel` | `numpy.ndarray` | PSF kernel (odd size, centred). |
| `radius` | `float` | Aperture radius [px]. |

| Returns | Type | Description |
| --- | --- | --- |
| `fraction` | `float` | Enclosed flux fraction. |

### photometry.fisherFluxErrors

`fisherFluxErrors(comps, variance, good, n=None)`

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

### photometry.psfFluxErrors

`psfFluxErrors(run)`

Return the Fisher errors of the PSF fluxes of a run's final joint fit.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |

| Returns | Type | Description |
| --- | --- | --- |
| `err` | `numpy.ndarray` | Error per source (fit variance, including the sysFrac term). |

### photometry.zeroPointFromRun

`zeroPointFromRun(run)`

Return the magnitude zero point of a run's image units.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |

| Returns | Type | Description |
| --- | --- | --- |
| `zeroPoint` | `float` | Zero point (from the exposure's PhotoCalib, else 31.4 for nJy). |
| `how` | `str` | Where it came from. |

### photometry.compareWithGaia

`compareWithGaia(phot, gaia, matchRadiusPx=1.0, offsetMinSnr=50.0, verbose=True)`

Add Gaia comparison columns to a photometry table.

Every fitted source is paired with a Gaia star: a source started from a Gaia position with that star (``gaiaIndex``), any other source with the nearest unclaimed Gaia star within ``matchRadiusPx``.  Positions are compared with Gaia's at the visit epoch (as placed on the detector by the visit WCS); magnitudes with Gaia G.  Gaia G is a much broader band than any LSST filter, so ``gaiaDmag`` contains a zero-point offset and a colour term: the median offset of bright, point-like matches is removed in ``gaiaDmagRel``, and ``gaiaBpRp`` is given to fit the colour term.

| Parameter | Type | Description |
| --- | --- | --- |
| `phot` | `astropy.table.Table` | Output of `measurePhotometry` (``x``, ``y``, ``psfFlux``, ``psfMag``, ``gaiaIndex``; ``ra``, ``dec`` if the exposure had a WCS). |
| `gaia` | `astropy.table.Table` | The Gaia catalogue the run started from (``session.gaia``: ``source_id``, ``ra``, ``dec``, ``x``, ``y``, ``phot_g_mean_mag``, optionally ``phot_bp_mean_mag`` and ``phot_rp_mean_mag``). |
| `matchRadiusPx` | `float`, optional | Match radius for sources that did not start from a Gaia position [px]. |
| `offsetMinSnr` | `float`, optional | Minimum PSF S/N of the matches that set the median magnitude offset. |
| `verbose` | `bool`, optional | Print a summary. |

| Returns | Type | Description |
| --- | --- | --- |
| `table` | `astropy.table.Table` | A copy of ``phot`` with ``gaiaMatch``, ``gaiaMatchType``, ``gaiaSourceId``, ``gaiaRa``, ``gaiaDec``, ``gaiaX``, ``gaiaY``, ``gaiaDx``, ``gaiaDy``, ``gaiaDra``, ``gaiaDdec``, ``gaiaSep``, ``gaiaG``, ``gaiaBpRp``, ``gaiaDmag``, ``gaiaDmagRel``, ``gaiaFlux`` and ``gaiaFluxRatio`` (see doc/catalog_columns.md); ``meta`` gets ``GAIAOFF`` (the median offset), ``GAIANOFF`` (stars used), ``GAIAUNM`` (Gaia stars in the fitted area without a fitted source) and ``GAIASYS``. |

### photometry.fitsMeta

`fitsMeta(meta)`

Return catalogue metadata with FITS-compatible keywords and values.

Long names are shortened (`FITS_META_KEYS`, else the first 8 characters, upper case); values that are not numbers, strings or booleans, and non-finite numbers, are dropped.

| Parameter | Type | Description |
| --- | --- | --- |
| `meta` | `dict` | Table metadata. |

| Returns | Type | Description |
| --- | --- | --- |
| `meta` | `dict` | The metadata to write. |

### photometry.measurePhotometry

`measurePhotometry(run, apertureFwhm=1.5, zeroPoint=None, verbose=True, gaia=None, gaiaMatchPx=1.0)`

Measure PSF and aperture magnitudes of a run's sources.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run (from ``runSessionGaia``). |
| `apertureFwhm` | `float`, optional | Aperture radius [FWHM]; aperture fluxes are corrected to total. |
| `zeroPoint` | `float`, optional | Override the zero point (default: PhotoCalib, else 31.4 for nJy). |
| `verbose` | `bool`, optional | Print a summary. |
| `gaia` | `astropy.table.Table` or `bool`, optional | The Gaia catalogue for the Gaia comparison columns of `compareWithGaia`; None (default) uses the run's own (``run.gaia``, the catalogue it started from) if it has one, False leaves the columns out. |
| `gaiaMatchPx` | `float`, optional | Gaia match radius for sources that did not start from a Gaia position [px]. |

| Returns | Type | Description |
| --- | --- | --- |
| `table` | `astropy.table.Table` | ``id``, ``x``, ``y``, ``ra``, ``dec``, ``gaia``, ``gaiaIndex``, ``gaiaG``, ``extended``, ``inInterior``, ``psfFlux``, ``psfFluxErr``, ``psfSnr``, ``psfMag``, ``psfMagErr``, ``apFlux``, ``apFluxErr``, ``apCoverage``, ``apMag``, ``apMagErr``, ``apMinusPsf`` (and, with ``gaia``, the `compareWithGaia` columns); ``meta`` has the zero point and aperture settings.  The columns and their units are described in ``doc/catalog_columns.md``. |

## plotting

`lsst.meas.crowd.plotting`: final data / model / residual plots and residual diagnostics.

Diagnostic plots of a finished run.

| Name | Summary |
| --- | --- |
| [`showBrightStarStack`](#plottingshowbrightstarstack) | Stack the bright stars' residuals to show the pattern they all share. |
| [`showStarZooms`](#plottingshowstarzooms) | Show zooms of individual stars: data \| model \| residual \| residual/sigma. |
| [`showCrossSections`](#plottingshowcrosssections) | Plot cuts through stars along x, along y and the radial profile. |
| [`compareRuns`](#plottingcompareruns) | Print side-by-side numbers for several runs of the same cutouts. |
| [`coreExcessDiagnostic`](#plottingcoreexcessdiagnostic) | Find out why bright stars' cores are sharper than the model. |
| [`plotResidualGray`](#plottingplotresidualgray) | Plot the two residual panels of `plotFinalComparison`, larger and in greyscale. |
| [`starResidualTable`](#plottingstarresidualtable) | Find which stars carry the leftover residual, and what they have in common. |
| [`plotFinalComparison`](#plottingplotfinalcomparison) | Plot original \| final model \| residual \| residual/sigma of a run. |
| [`plotGaiaComparison`](#plottingplotgaiacomparison) | Plot how the fitted photometry and positions compare with Gaia. |

### plotting.showBrightStarStack

`showBrightStarStack(run, minSnr=200.0, half=20, plot=True)`

Stack the bright stars' residuals to show the pattern they all share.

The median of the bright, clean, unsaturated stars' residuals, each divided by the star's flux and centred on its sub-pixel position, once with the PSF-only joint model and once with the final model (including the per-star corrections).  A pattern that survives the median is common to all bright stars, i.e. a feature missing from the PSF model; noise and star-specific errors average out.  Also prints the core error of the brighter and fainter half of the stars (a brighter-fatter test) and shows the PSF kernel.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `minSnr` | `float`, optional | Minimum peak S/N of the stacked stars. |
| `half` | `int`, optional | Half-size of the stack [px]. |
| `plot` | `bool`, optional | Show the kernel and the two stacks. |

| Returns | Type | Description |
| --- | --- | --- |
| `stack` | `dict` | ``psfOnly`` and ``final`` (stacked residual / flux images), ``n`` (number of stars) and, with at least six stars, ``byBrightness`` (core and ring residual of the brighter and fainter half). |

### plotting.showStarZooms

`showStarZooms(run, stars=None, n=8, half=12, sigma='model', grayscale=False, minSnr=30.0, cuts=False)`

Show zooms of individual stars: data | model | residual | residual/sigma.

Every source in the window is marked (o = started from Gaia, x = found in the residual), with a line from each Gaia source's catalogue position (+) to its fitted position.  Prints each star's neighbours (offset, flux ratio, Gaia or not, distance from its Gaia position).

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `stars` | `list` [`tuple` [`float`, `float`]], optional | Detector positions; None = the ``n`` worst stars of `starResidualTable`. |
| `n` | `int`, optional | Number of worst stars shown. |
| `half` | `int`, optional | Half-size of each zoom [px]. |
| `sigma` | `str`, optional | 'model', 'effective' or 'pixel' (see `plotFinalComparison`). |
| `grayscale` | `bool`, optional | Residual panels in greyscale. |
| `minSnr` | `float`, optional | Minimum peak S/N when choosing the worst stars. |
| `cuts` | `bool`, optional | Also draw `showCrossSections` for the same stars. |

| Returns | Type | Description |
| --- | --- | --- |
| `indices` | `list` [`int`] or `None` | Indices of the stars shown (None if there were none). |

### plotting.showCrossSections

`showCrossSections(run, stars=None, n=4, half=15, width=1, sigma='model', minSnr=30.0, log=True)`

Plot cuts through stars along x, along y and the radial profile.

Top: data (points) and model (line); bottom: residual with its +-1 sigma band, so a ring shows as a bump in the radial residual at its radius.  Other sources in a cut are marked (dotted: Gaia, dashed: found in the residual).

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `stars` | `list` [`tuple` [`float`, `float`]], optional | Detector positions; None = the ``n`` worst stars of `starResidualTable`. |
| `n` | `int`, optional | Number of worst stars shown. |
| `half` | `int`, optional | Half-length of the cuts [px]. |
| `width` | `int`, optional | Rows / columns averaged for the x / y cuts (odd; 1 = single pixel). |
| `sigma` | `str`, optional | 'model', 'effective' or 'pixel' (see `plotFinalComparison`). |
| `minSnr` | `float`, optional | Minimum peak S/N when choosing the worst stars. |
| `log` | `bool`, optional | Symmetric-log axis for data and model (wings visible). |

| Returns | Type | Description |
| --- | --- | --- |
| `cuts` | `list` [`dict`] or `None` | Per star ``index`` and ``cuts`` (the x, y and radial profiles); None if no stars. |

### plotting.compareRuns

`compareRuns(variants, minSnr=200.0)`

Print side-by-side numbers for several runs of the same cutouts.

For each cutout and variant: robust residual/sigma, number of sources and the stacked PSF-only residual of the bright stars (median residual / flux x 1e5) in radial bins -- a ring shows as a large value in one bin; values closer to 0 mean a better PSF.

| Parameter | Type | Description |
| --- | --- | --- |
| `variants` | `dict` [`str`, `dict`] | Variant name -> runs (as returned by ``runSessionGaia``), e.g. ``{"default": runs, "wider wing join": runs_b}``. |
| `minSnr` | `float`, optional | Minimum peak S/N of the stacked stars. |

### plotting.coreExcessDiagnostic

`coreExcessDiagnostic(run, minSnr=100.0, plot=True)`

Find out why bright stars' cores are sharper than the model.

For each unsaturated, fairly isolated star with peak S/N >= ``minSnr``, the central excess (data - model summed over the central 3x3 px, divided by the model there) is correlated with flux (sensor effects or a flux-dependent PSF), peak / saturation level (non-linearity) and sub-pixel phase (stamps blurred by sub-pixel shifting).  No trend but a significant excess means the PSF core itself is too broad.  A trend only counts if it is large compared with the excess itself.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `minSnr` | `float`, optional | Minimum peak S/N of the stars used. |
| `plot` | `bool`, optional | Plot the excess against each quantity. |

| Returns | Type | Description |
| --- | --- | --- |
| `result` | `dict` or `None` | Per-star arrays (``idx``, ``flux``, ``snr``, ``peak``, ``excess``, ``excessPsfOnly``, ``centrePx``, ``phase``, ``err``), ``correlations`` and ``verdict``; None if fewer than six stars. |

### plotting.plotResidualGray

`plotResidualGray(run, **kw)`

Plot the two residual panels of `plotFinalComparison`, larger and in greyscale.

Black = negative, white = positive; areas that cannot be modelled are drawn in blue.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `**kw` | – | Keyword arguments of `plotFinalComparison`. |

| Returns | Type | Description |
| --- | --- | --- |
| `stats` | `dict` | See `plotFinalComparison`. |

### plotting.starResidualTable

`starResidualTable(run, minSnr=30.0, rFwhm=2.5, plot=True, nWorst=10)`

Find which stars carry the leftover residual, and what they have in common.

For every unsaturated point source with peak S/N >= ``minSnr``: chi2 per pixel within ``rFwhm`` FWHM for the final and the PSF-only model, the PSF-only core residual (flux within 1 FWHM / star flux; < 0 = star broader than the PSF) and the mismatch (residual flux beyond the noise / star flux, independent of brightness).  The mismatch is rank-correlated with flux, position and crowding, and outlier stars are counted.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `minSnr` | `float`, optional | Minimum peak S/N. |
| `rFwhm` | `float`, optional | Radius of the chi2 region [FWHM]. |
| `plot` | `bool`, optional | Plot chi2 and core residual against flux and the stars' positions. |
| `nWorst` | `int`, optional | Number of worst stars listed. |

| Returns | Type | Description |
| --- | --- | --- |
| `table` | `dict` or `None` | Per-star arrays (``idx``, ``flux``, ``snr``, ``x``, ``y``, ``dist``, ``nNb``, ``chi2Final``, ``chi2Psf``, ``corePsf``, ``npix``, ``mismatch``), ``correlations``, ``coreVsFlux``, ``outliers`` and ``verdict``; None if fewer than eight stars. |

### plotting.plotFinalComparison

`plotFinalComparison(run, view='interior', badColor='#ff00ff', circles='images', showZones=True, excludeNearMask=True, excludedColor='#8c8c8c', maskPadFwhm=1.5, zClip=5.0, residClipPct=99.0, figsize=(26, 7), sigma='model', nearsatColor='#ff9900', grayscale=False, panels='all', verbose=True)`

Plot original | final model | residual | residual/sigma of a run.

Bad pixels are drawn opaque in one colour (default magenta) in all four panels and excluded from the statistics.  Saturated-star zones have a dashed outline, extended sources dashed circles.  Pixels whose flux belongs to sources under the mask (which the model cannot fit) are greyed out and excluded too.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `view` | `str` or `tuple`, optional | 'interior' (the requested box), None (the whole fitted array) or (xmin, xmax, ymin, ymax) in detector pixels. |
| `badColor` | `str`, optional | Colour of bad pixels. |
| `circles` | `str`, optional | Where to draw the source circles: 'images' (data and model), 'all' or 'none'. |
| `showZones` | `bool`, optional | Outline the saturated-star zones. |
| `excludeNearMask` | `bool`, optional | Grey out (and exclude) pixels next to bad pixels that cannot be modelled. |
| `excludedColor` | `str`, optional | Colour of those pixels. |
| `maskPadFwhm` | `float`, optional | Pad around bad pixels [FWHM]. |
| `zClip` | `float`, optional | Colour range of residual/sigma. |
| `residClipPct` | `float`, optional | Percentile setting the colour range of the residual. |
| `figsize` | `tuple` [`float`, `float`], optional | Figure size [inches]. |
| `sigma` | `str`, optional | 'model' = sqrt(variance + (eps*model)^2 + k*model) with the measured PSF model error eps; 'effective' = without eps; 'pixel' = the variance plane only. |
| `nearsatColor` | `str`, optional | Colour of near-saturated cores fitted from their wings. |
| `grayscale` | `bool`, optional | Residual panels in greyscale (excluded areas then in blue). |
| `panels` | `str`, optional | 'all' (four panels) or 'residuals' (the two residual panels, larger). |
| `verbose` | `bool`, optional | Print the statistics. |

| Returns | Type | Description |
| --- | --- | --- |
| `stats` | `dict` | ``residual`` and ``residualOverSigma`` (median, robust sigma, quartiles, number of pixels) and ``excludedPx``. |

### plotting.plotGaiaComparison

`plotGaiaComparison(phot, minSnr=20.0, figsize=(22, 5.5), out=None)`

Plot how the fitted photometry and positions compare with Gaia.

Four panels for the Gaia-matched point sources with PSF S/N >= ``minSnr``: magnitude difference (offset removed) against Gaia G and against Gaia BP - RP (the colour term between G and the LSST band), the position offsets, and the offsets as arrows across the image (a pattern there is a WCS distortion; scatter is centroid noise).

| Parameter | Type | Description |
| --- | --- | --- |
| `phot` | `astropy.table.Table` | Output of `lsst.meas.crowd.photometry.compareWithGaia` (or of ``measurePhotometry`` with ``gaia=``). |
| `minSnr` | `float`, optional | Minimum PSF S/N of the plotted sources. |
| `figsize` | `tuple` [`float`, `float`], optional | Figure size [inches]. |
| `out` | `str`, optional | Also save the figure to this file. |

| Returns | Type | Description |
| --- | --- | --- |
| `summary` | `dict` | ``n`` (sources plotted), ``dmagScatter`` (robust std of ``gaiaDmagRel`` [mag]), ``offset`` (the removed median offset [mag]) and ``dx``, ``dy`` (median position offsets [mas] or, without sky coordinates, [px]). |

Raises: `ValueError` (Raised if the table has no Gaia comparison columns.)

## stagePlots

`lsst.meas.crowd.stagePlots`: plots of the intermediate stages and pass-to-pass comparison.

Plots of the pipeline's intermediate stages.

| Name | Summary |
| --- | --- |
| [`psfStarPositions`](#stageplotspsfstarpositions) | Return the candidate and used PSF stars of one of the pipeline's PSF determinations. |
| [`showPsfStars`](#stageplotsshowpsfstars) | Plot the PSF stars over an image, the PSF kernel and its radial profile. |
| [`showBackgrounds`](#stageplotsshowbackgrounds) | Plot the background entering every outer pass, and the change from pass to pass. |
| [`showRounds`](#stageplotsshowrounds) | Plot data \| model \| residual \| residual/sigma for each round of the fit. |
| [`comparePasses`](#stageplotscomparepasses) | Find the sources that one outer pass lost or gained relative to another, and why. |
| [`showLinkedSourcePages`](#stageplotsshowlinkedsourcepages) | Plot per-source pages whose coloured frames match boxes on the data, model and residual. |

### stagePlots.psfStarPositions

`psfStarPositions(session=None, catalog=None, run=None, outer=-1)`

Return the candidate and used PSF stars of one of the pipeline's PSF determinations.

Exactly one source is used, in this order: a catalogue with ``calib_psf_candidate`` / ``calib_psf_used`` flags (the visit PSF's own stars), a session (the detector-level ``psfSource='refit'`` PSF), or a run (the clean stars the residual-stack PSF update after pass ``outer`` used).

| Parameter | Type | Description |
| --- | --- | --- |
| `session` | `lsst.meas.crowd.session.Session`, optional | A session prepared with ``psfSource='refit'``. |
| `catalog` | `lsst.afw.table.SourceCatalog` or `astropy.table.Table`, optional | A source catalogue with PSF flags, e.g. the visit's ``single_visit_star_unstandardized``. |
| `run` | `lsst.meas.crowd.session.CutoutRun`, optional | A finished run. |
| `outer` | `int`, optional | Outer pass of ``run`` whose result the PSF update used. |

| Returns | Type | Description |
| --- | --- | --- |
| `stars` | `dict` | ``candidates`` and ``used`` (positions in detector px, shape (N, 2)) and ``source`` (a description). |

Raises: `ValueError` (Raised if no source of PSF stars is given.)

### stagePlots.showPsfStars

`showPsfStars(image, stars=None, psf=None, at=None, view=None, figsize=(19, 6), out=None)`

Plot the PSF stars over an image, the PSF kernel and its radial profile.

| Parameter | Type | Description |
| --- | --- | --- |
| `image` | `lsst.afw.image.Exposure` | The visit image (or a cutout) to draw the stars on. |
| `stars` | `dict`, optional | Output of `psfStarPositions`; None = draw no stars. |
| `psf` | `lsst.afw.detection.Psf`, optional | PSF to show (default: the image's own PSF). |
| `at` | `tuple` [`float`, `float`], optional | Detector position at which to render the kernel (default: the image centre). |
| `view` | `tuple` [`float`, `float`, `float`, `float`], optional | (xmin, xmax, ymin, ymax) zoom in detector pixels. |
| `figsize` | `tuple` [`float`, `float`], optional | Figure size [inches]. |
| `out` | `str`, optional | Also save the figure to this file. |

| Returns | Type | Description |
| --- | --- | --- |
| `info` | `dict` | ``kernel`` (the rendered kernel image), ``fwhm`` [px], ``nCandidates`` and ``nUsed``. |

### stagePlots.showBackgrounds

`showBackgrounds(run, initial=None, figsizePerPanel=4.2, out=None)`

Plot the background entering every outer pass, and the change from pass to pass.

Pass 1 uses ``initial`` (the detector-level background of ``psfSource='refit'``; a visit image is already background-subtracted, so there is none).  Before each later pass the background is refitted on the raw cutout with the previous pass's sources masked.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `initial` | `lsst.afw.image.Image`, optional | The background subtracted before pass 1 (``session.bgImage``). |
| `figsizePerPanel` | `float`, optional | Size of each panel [inches]. |
| `out` | `str`, optional | Also save the figure to this file. |

| Returns | Type | Description |
| --- | --- | --- |
| `backgrounds` | `list` [`numpy.ndarray` or `None`] | Background entering each pass (None = nothing was subtracted). |

### stagePlots.showRounds

`showRounds(run, level='outer', outer=-1, rounds=None, maxRows=8, view=None, zClip=5.0, figsizePerPanel=3.6, out=None)`

Plot data | model | residual | residual/sigma for each round of the fit.

``level='outer'``: one row per outer pass (each with its own background and PSF). ``level='inner'``: one row per inner iteration of outer pass ``outer`` (the source list growing as the detection threshold steps down).  Data and model share one colour scale, and the residual one scale for all rows, so the rows can be compared.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `level` | `str`, optional | 'outer' or 'inner'. |
| `outer` | `int`, optional | Outer pass whose inner iterations are shown (``level='inner'``). |
| `rounds` | `list` [`int`], optional | Rounds to show (default: up to ``maxRows``, evenly spaced, the last included). |
| `maxRows` | `int`, optional | Maximum number of rows when ``rounds`` is not given. |
| `view` | `tuple` [`float`, `float`, `float`, `float`], optional | (xmin, xmax, ymin, ymax) zoom in detector pixels. |
| `zClip` | `float`, optional | Colour range of residual/sigma. |
| `figsizePerPanel` | `float`, optional | Size of each panel [inches]. |
| `out` | `str`, optional | Also save the figure to this file. |

| Returns | Type | Description |
| --- | --- | --- |
| `rows` | `list` [`dict`] | Per row ``label``, ``nSources``, ``residStd`` (robust std of the residual) and ``chiStd`` (robust std of residual/sigma) over the good pixels in view. |

Raises: `ValueError` (Raised if ``level`` is not 'outer' or 'inner'.)

### stagePlots.comparePasses

`comparePasses(run, a=-2, b=-1, matchPx=1.0, nShow=8, half=12, plot=True, out=None)`

Find the sources that one outer pass lost or gained relative to another, and why.

Sources are matched by position (within ``matchPx``).  For every source of pass ``a`` without a match in pass ``b`` the residual of pass ``b`` is measured at its position with a PSF-weighted (matched-filter) flux: if most of its flux is still there, a real source was lost; if a neighbour within 1.5 FWHM gained that flux, it was merged; if a halo or ring check removed it, it is also blocked from re-detection in later passes.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run with at least two outer passes. |
| `a, b` | `int`, optional | The two outer passes to compare (default: the last two). |
| `matchPx` | `float`, optional | Matching radius [px]. |
| `nShow` | `int`, optional | Number of lost sources drawn (brightest first). |
| `half` | `int`, optional | Half-size of each zoom [px]. |
| `plot` | `bool`, optional | Draw the overview and the zooms. |
| `out` | `str`, optional | Also save the figure to this file. |

| Returns | Type | Description |
| --- | --- | --- |
| `table` | `astropy.table.Table` | One row per lost or gained source: ``kind`` ('lost' / 'gained'), ``x``, ``y``, ``flux``, ``snr`` (peak S/N in the pass it exists in), ``gaia``, ``fluxLeft`` and ``leftSnr`` (matched-filter flux and S/N in the residual of the pass without it), ``leftFrac``, ``nearestPx`` (distance to the nearest source of the other pass), ``nearestDflux`` (that source's flux change) and ``verdict``. |

Raises: `ValueError` (Raised if the run has fewer than two outer passes.)

### stagePlots.showLinkedSourcePages

`showLinkedSourcePages(run, outer=-1, snap=-1, sources=None, sourceView='all', sourcesPerPage=6, maxPages=1, viewHalf=20, showPages=True, outDir=None, pdf=None, dpi=80, figsize=(21, 7.5), out=None)`

Plot per-source pages whose coloured frames match boxes on the data, model and residual.

For each page: an overview (data | model | residual of the whole cutout) with one coloured box per source on that page -- the window its row shows -- then the page itself: one row per source with data | model | residual | x cross-section | y cross-section, framed in the same colour.

For many pages use ``outDir``: every figure is written to its own PNG file and closed at once, so memory stays flat.  A `PdfPages` (``pdf``) keeps every image in memory until it is closed (matplotlib), which can exhaust a notebook's memory for hundreds of sources.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `outer` | `int`, optional | Outer pass. |
| `snap` | `int`, optional | Inner-iteration snapshot of that pass; only the last one keeps the per-source stamps unless the run used ``keepAllFits=True``. |
| `sources` | `list`, optional | Source indices, or (x, y) detector positions (the nearest source is used); None = choose with ``sourceView``. |
| `sourceView` | `str`, optional | 'all' (brightest first), 'new' or 'changed' (see `fitLib.selectSourcePages`). |
| `sourcesPerPage` | `int`, optional | Rows per page. |
| `maxPages` | `int`, optional | Maximum number of pages (None = all). |
| `viewHalf` | `int`, optional | Half-size of each source window [px] (the size of the boxes). |
| `showPages` | `bool`, optional | Draw the pages (False = only the overviews). |
| `outDir` | `str`, optional | Write ``overview_NNN.png`` and ``page_NNN.png`` here instead of showing them. |
| `pdf` | `matplotlib.backends.backend_pdf.PdfPages`, optional | Open PDF; the overviews and pages are saved into it instead of being shown. |
| `dpi` | `int`, optional | Resolution of the saved pages. |
| `figsize` | `tuple` [`float`, `float`], optional | Size of the overview [inches]. |
| `out` | `str`, optional | Also save the first overview to this file. |

| Returns | Type | Description |
| --- | --- | --- |
| `linked` | `dict` | ``pages`` (source indices per page), ``colors`` (colour per source index) and ``files`` (the PNG files written with ``outDir``). |

Raises: `ValueError` (Raised if the snapshot no longer has the per-source stamps.)

## sourcePlots

`lsst.meas.crowd.sourcePlots`: per-source cutouts, cross-sections and magnitudes.

Per-source cutouts, cross-sections and magnitudes.

| Name | Summary |
| --- | --- |
| [`sourceCutoutData`](#sourceplotssourcecutoutdata) | Compute the arrays one cutout row shows for a source. |
| [`selectSources`](#sourceplotsselectsources) | Choose the sources to show. |
| [`plotSourceCutouts`](#sourceplotsplotsourcecutouts) | Plot cutouts and cross-sections of individual sources, labelled with their magnitudes. |

### sourcePlots.sourceCutoutData

`sourceCutoutData(data, model, var, bad, x, y, x0=0, y0=0, half=15, width=1)`

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

### sourcePlots.selectSources

`selectSources(phot, how='worst', n=6, positions=None, xy=None, chi2=None, minSnr=20.0, seed=0)`

Choose the sources to show.

| Parameter | Type | Description |
| --- | --- | --- |
| `phot` | `astropy.table.Table` or `dict` | Photometry with ``psfSnr``. |
| `how` | `str`, optional | 'worst' (highest ``chi2``), 'brightest', 'faintest' or 'random', among sources with S/N >= ``minSnr``. |
| `n` | `int`, optional | Number of sources. |
| `positions` | `list` [`tuple` [`float`, `float`]], optional | Show the nearest source to each position instead. |
| `xy` | `numpy.ndarray`, optional | Source positions (needed with ``positions``). |
| `chi2` | `numpy.ndarray`, optional | Local chi2 per source (needed for 'worst'). |
| `minSnr` | `float`, optional | Minimum S/N. |
| `seed` | `int`, optional | Random seed for 'random'. |

| Returns | Type | Description |
| --- | --- | --- |
| `indices` | `list` [`int`] | Selected source indices. |

Raises: `ValueError` (Raised if ``how`` is not a known selection.)

### sourcePlots.plotSourceCutouts

`plotSourceCutouts(run, how='worst', n=6, positions=None, half=15, width=1, sigma='model', phot=None, minSnr=20.0, log=True, perPage=4, pdfPath=None, show=True, grayscale=False)`

Plot cutouts and cross-sections of individual sources, labelled with their magnitudes.

One row per source: data | model | residual | residual/sigma | x-cut | y-cut | radial profile.  The cuts show data (points) and model (line), and the residual with its +-1 sigma band; other sources in the window are marked (o Gaia, x found in the residual).

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `how` | `str`, optional | 'worst', 'brightest', 'faintest' or 'random' (see `selectSources`). |
| `n` | `int`, optional | Number of sources. |
| `positions` | `list` [`tuple` [`float`, `float`]], optional | Show the sources nearest these detector positions. |
| `half` | `int`, optional | Half-size of each window [px]. |
| `width` | `int`, optional | Rows / columns averaged in the cuts. |
| `sigma` | `str`, optional | 'model' (variance plus measured model error), 'effective' or 'pixel'. |
| `phot` | `astropy.table.Table`, optional | Output of `lsst.meas.crowd.photometry.measurePhotometry` (computed if not given). |
| `minSnr` | `float`, optional | Minimum S/N of the sources selected. |
| `log` | `bool`, optional | Symmetric-log axis for data and model. |
| `perPage` | `int`, optional | Sources per figure. |
| `pdfPath` | `str`, optional | Also write all pages to this PDF. |
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
| [`tileSpecs`](#detectortilespecs) | Return the interiors of a tile grid that partitions an area exactly. |
| [`tileName`](#detectortilename) | Return the name of a tile, e.g. 'x0400_y0000_400x400'. |
| [`addTileResiduals`](#detectoraddtileresiduals) | Add the RESIDUAL and RESIDUAL_SIGMA images to tile files written without them. |
| [`saveTileProducts`](#detectorsavetileproducts) | Write one tile's results to a FITS file. |
| [`fitTiles`](#detectorfittiles) | Fit this task's share of a detector's tiles and write one file per tile. |
| [`stitchTiles`](#detectorstitchtiles) | Assemble the tile files into full-detector images and one catalogue. |
| [`loadMosaic`](#detectorloadmosaic) | Read the files written by `stitchTiles`. |
| [`loadTile`](#detectorloadtile) | Read one tile file written by `saveTileProducts`. |
| [`plotTile`](#detectorplottile) | Plot data \| model \| residual \| residual/sigma of one tile over its whole fit box. |
| [`plotDetector`](#detectorplotdetector) | Plot data \| model \| residual \| residual/sigma of a stitched detector or a zoom. |

### detector.tileSpecs

`tileSpecs(width, height, tile=400, x0=0, y0=0, region=None)`

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

### detector.tileName

`tileName(spec)`

Return the name of a tile, e.g. 'x0400_y0000_400x400'.

| Parameter | Type | Description |
| --- | --- | --- |
| `spec` | `tuple` [`int`, `int`, `int`, `int`] | (x0, y0, w, h). |

| Returns | Type | Description |
| --- | --- | --- |
| `name` | `str` | The tile name. |

### detector.addTileResiduals

`addTileResiduals(outDir='detector_run')`

Add the RESIDUAL and RESIDUAL_SIGMA images to tile files written without them.

Tiles fitted by older versions have only DATA, MODEL, VARIANCE and BAD; this adds the two residual images (nothing is refitted).  Files that already have them are left alone.

| Parameter | Type | Description |
| --- | --- | --- |
| `outDir` | `str`, optional | Directory with the tile files. |

| Returns | Type | Description |
| --- | --- | --- |
| `updated` | `list` [`str`] | The files that were rewritten. |

### detector.saveTileProducts

`saveTileProducts(run, path, phot=None, meta=None)`

Write one tile's results to a FITS file.

DATA, MODEL, RESIDUAL, RESIDUAL_SIGMA, VARIANCE and BAD images over the tile's fit box (RESIDUAL = DATA - MODEL, RESIDUAL_SIGMA = RESIDUAL / sqrt(VARIANCE); header keys X0/Y0 = fit-box origin, IX0/IY0/IX1/IY1 = interior, inclusive) and a CATALOG table (the photometry, with an ``owned`` column: source inside this tile's interior, half-open so boundary sources belong to exactly one tile).  The file is written atomically.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | The tile's finished run. |
| `path` | `str` | Output file. |
| `phot` | `astropy.table.Table`, optional | Photometry (computed if not given). |
| `meta` | `dict`, optional | Extra header keywords. |

| Returns | Type | Description |
| --- | --- | --- |
| `path` | `str` | The file written. |

### detector.fitTiles

`fitTiles(visit, detector, config=None, outDir='detector_run', tile=400, pad=30, taskIndex=0, nTasks=1, region=None, repo='dp2_prep', collection='LSSTCam/runs/DRP/DP2', instrument='LSSTCam', visitImageCollections='LSSTCam/runs/*', psfSource='visit', overwrite=False, visitImage=None, plots=None, plotDir=None, allSourcePages=False, **overrides)`

Fit this task's share of a detector's tiles and write one file per tile.

Tiles already on disk are skipped, so the job can be restarted; a failing tile is reported and the task moves on.  With ``plots``, every diagnostic plot of the chosen tiles (`lsst.meas.crowd.runPlots.saveRunPlots`) goes to ``<plotDir>/<tile name>/``, and task 0 also writes the detector-level PSF plots to ``<plotDir>/detector/``; a chosen tile whose plots are missing is refitted even if its tile file exists.

| Parameter | Type | Description |
| --- | --- | --- |
| `visit, detector` | `int` | Visit and detector ids. |
| `config` | `dict`, optional | Pipeline configuration (default ``makeConfig(CROWDED_FIELD_CONFIG)``). |
| `outDir` | `str`, optional | Output directory. |
| `tile, pad` | `int`, optional | Tile interior size and fit padding [px]. |
| `taskIndex, nTasks` | `int`, optional | This task fits tiles ``taskIndex``, ``taskIndex + nTasks``, ... (one Slurm array task each). |
| `region` | `tuple` [`int`, `int`, `int`, `int`], optional | (xmin, xmax, ymin, ymax): only tiles overlapping it (for tests). |
| `repo` | `str`, optional | Butler repository path or alias. |
| `collection` | `str`, optional | Default Butler collection. |
| `instrument` | `str`, optional | Instrument name. |
| `visitImageCollections` | `str`, optional | Collections searched for the visit image. |
| `psfSource` | `str`, optional | 'visit' or 'refit' (see `lsst.meas.crowd.session.buildDetectorPsf`). |
| `overwrite` | `bool`, optional | Refit tiles that already have a file. |
| `visitImage` | `lsst.afw.image.Exposure`, optional | An already-loaded visit image (skips the Butler). |
| `plots` | `str` or `list` [`str`], optional | Tiles to write all plots for: None (no plots), 'sample' (the tiles at the corners, the centre and two edge midpoints, `lsst.meas.crowd.runPlots.sampleTileNames`), 'all', or a list of tile names. |
| `plotDir` | `str`, optional | Directory for the plots (default ``<outDir>/plots``). |
| `allSourcePages` | `bool`, optional | Write a per-source page for every source of a plotted tile (default: 4 pages). |
| `**overrides` | – | Keyword arguments of ``runSessionGaia`` (e.g. ``maxOuter=3, maxInner=15``). |

| Returns | Type | Description |
| --- | --- | --- |
| `paths` | `list` [`str`] | Tile files of this task that exist. |

### detector.stitchTiles

`stitchTiles(outDir='detector_run', mosaicName='detector_mosaic.fits', catalogName='detector_catalog.fits', verbose=True)`

Assemble the tile files into full-detector images and one catalogue.

Writes the mosaic (DATA, MODEL, RESIDUAL, RESIDUAL_SIGMA, VARIANCE, BAD and TILE extensions) and the catalogue (every source once, from the tile whose interior contains it) into ``outDir``.  Tile files from older versions without the residual images get them added (`addTileResiduals`).

| Parameter | Type | Description |
| --- | --- | --- |
| `outDir` | `str`, optional | Directory with the tile files. |
| `mosaicName` | `str`, optional | Mosaic file name. |
| `catalogName` | `str`, optional | Catalogue file name. |
| `verbose` | `bool`, optional | Print a summary. |

| Returns | Type | Description |
| --- | --- | --- |
| `mosaic` | `dict` | ``data``, ``model``, ``resid``, ``residSigma``, ``var``, ``bad``, ``tile`` (tile index map, -1 = not fitted), ``x0``, ``y0``, ``tiles`` (names), ``interiors`` and ``catalog``. |

Raises: `FileNotFoundError` (Raised if ``outDir`` has no tile files.)

### detector.loadMosaic

`loadMosaic(outDir='detector_run', mosaicName='detector_mosaic.fits', catalogName='detector_catalog.fits')`

Read the files written by `stitchTiles`.

| Parameter | Type | Description |
| --- | --- | --- |
| `outDir` | `str`, optional | Directory. |
| `mosaicName` | `str`, optional | Mosaic file name. |
| `catalogName` | `str`, optional | Catalogue file name. |

| Returns | Type | Description |
| --- | --- | --- |
| `mosaic` | `dict` | As returned by `stitchTiles` (without ``tiles`` and ``interiors``). |

### detector.loadTile

`loadTile(path)`

Read one tile file written by `saveTileProducts`.

| Parameter | Type | Description |
| --- | --- | --- |
| `path` | `str` | The ``tile_<name>.fits`` file. |

| Returns | Type | Description |
| --- | --- | --- |
| `tile` | `dict` | In the form of `loadMosaic`'s result over the tile's fit box (``data``, ``model``, ``resid``, ``residSigma``, ``var``, ``bad``, ``tile``, ``x0``, ``y0``, ``catalog``), plus ``name`` and ``interiors`` (the tile's interior as one (x0, y0, x1, y1) tuple, inclusive). |

### detector.plotTile

`plotTile(tile, zClip=5.0, grayscale=False, sources=True, figsize=(24, 6.5), out=None)`

Plot data | model | residual | residual/sigma of one tile over its whole fit box.

The tile's interior (the part kept in the mosaic) is outlined; outside it is the padding, fitted only so that sources near the edge have their neighbours.

| Parameter | Type | Description |
| --- | --- | --- |
| `tile` | `dict` or `str` | Output of `loadTile`, or the tile file. |
| `zClip` | `float`, optional | Colour range of residual/sigma. |
| `grayscale` | `bool`, optional | Residual panels in greyscale. |
| `sources` | `bool`, optional | Circle the tile catalogue's sources. |
| `figsize` | `tuple` [`float`, `float`], optional | Figure size [inches]. |
| `out` | `str`, optional | Also save the figure to this file. |

| Returns | Type | Description |
| --- | --- | --- |
| `result` | `dict` | As `plotDetector`'s. |

### detector.plotDetector

`plotDetector(mosaic, view=None, binning=None, showTiles=True, zClip=5.0, grayscale=False, figsize=(26, 7.5), sources=False, out=None)`

Plot data | model | residual | residual/sigma of a stitched detector or a zoom.

| Parameter | Type | Description |
| --- | --- | --- |
| `mosaic` | `dict` or `str` | Output of `stitchTiles` / `loadMosaic`, or the directory to load. |
| `view` | `tuple` [`float`, `float`, `float`, `float`], optional | (xmin, xmax, ymin, ymax) zoom in detector pixels. |
| `binning` | `int`, optional | Display block size (default: about 1000 px across the view); residual/sigma is binned as sum / sqrt(n) so its noise stays 1. |
| `showTiles` | `bool`, optional | Draw the tile boundaries. |
| `zClip` | `float`, optional | Colour range of residual/sigma. |
| `grayscale` | `bool`, optional | Residual panels in greyscale. |
| `figsize` | `tuple` [`float`, `float`], optional | Figure size [inches]. |
| `sources` | `bool`, optional | Circle the catalogue's sources (for zooms). |
| `out` | `str`, optional | Also save the figure to this file. |

| Returns | Type | Description |
| --- | --- | --- |
| `result` | `dict` | ``residualSigma`` (robust sigma of residual/sigma over the view) and ``binning``. |

## runPlots

`lsst.meas.crowd.runPlots`: write every plot of a run, of the detector-level PSF, or of chosen tiles.

Write every diagnostic plot of a run (or of a detector) to PNG files.

| Name | Summary |
| --- | --- |
| [`sampleTileNames`](#runplotssampletilenames) | Return the names of the tiles that contain the sample points. |
| [`showPsfStarStamps`](#runplotsshowpsfstarstamps) | Plot a stamp around each PSF star, followed by the PSF kernel. |
| [`saveRunPlots`](#runplotssaverunplots) | Write every diagnostic plot of one finished run to PNG files. |
| [`saveDetectorPlots`](#runplotssavedetectorplots) | Write the detector-level PSF plots: the PSF stars over the detector, their stamps, the PSF. |

### runPlots.sampleTileNames

`sampleTileNames(specs, width, height, x0=0, y0=0, points=((0.02, 0.02), (0.98, 0.02), (0.02, 0.98), (0.98, 0.98), (0.5, 0.5), (0.5, 0.02), (0.02, 0.5)))`

Return the names of the tiles that contain the sample points.

| Parameter | Type | Description |
| --- | --- | --- |
| `specs` | `list` [`tuple` [`int`, `int`, `int`, `int`]] | Tiles as (x0, y0, width, height), from `lsst.meas.crowd.detector.tileSpecs`. |
| `width, height` | `int` | Detector size [px]. |
| `x0, y0` | `int`, optional | Detector origin [px]. |
| `points` | `tuple` [`tuple` [`float`, `float`]], optional | Sample points as fractions of the detector width and height. |

| Returns | Type | Description |
| --- | --- | --- |
| `names` | `list` [`str`] | Tile names (`lsst.meas.crowd.detector.tileName`), each once, in the order of ``points``. |

### runPlots.showPsfStarStamps

`showPsfStarStamps(image, xy, half=15, maxStamps=25, psf=None, title='', out=None)`

Plot a stamp around each PSF star, followed by the PSF kernel.

| Parameter | Type | Description |
| --- | --- | --- |
| `image` | `lsst.afw.image.Exposure` | Image to cut the stamps from (detector pixel coordinates). |
| `xy` | `numpy.ndarray` | Star positions [px], shape (N, 2); the first ``maxStamps`` are shown. |
| `half` | `int`, optional | Half-size of each stamp [px]. |
| `maxStamps` | `int`, optional | Maximum number of stamps. |
| `psf` | `lsst.afw.detection.Psf`, optional | PSF whose kernel is shown in the last panel (default: the image's own PSF). |
| `title` | `str`, optional | Figure title. |
| `out` | `str`, optional | Also save the figure to this file. |

| Returns | Type | Description |
| --- | --- | --- |
| `nShown` | `int` | Number of stamps drawn. |

### runPlots.saveRunPlots

`saveRunPlots(run, outDir, phot=None, initialBackground=None, sourcePages=4, allSourcePages=False, dpi=100, verbose=True)`

Write every diagnostic plot of one finished run to PNG files.

In pipeline order: PSF stars and their stamps for each PSF update, the PSF model and FWHM of each pass, the background of each pass, data | model | residual per outer pass and per inner iteration, the sources gained and lost between passes, the final comparison, per-source pages (brightest first, and spread in magnitude), labelled source cutouts, the bright-star diagnostics and the Gaia comparison.  ``index.txt`` describes every file and ``report.txt`` holds the pass-by-pass report.  A plot that fails is listed in ``index.txt`` and skipped.

| Parameter | Type | Description |
| --- | --- | --- |
| `run` | `lsst.meas.crowd.session.CutoutRun` | A finished run. |
| `outDir` | `str` | Directory for the files (created if needed). |
| `phot` | `astropy.table.Table`, optional | Photometry of the run (default ``run.photometry`` if present); needed for the magnitude-labelled plots and the Gaia comparison. |
| `initialBackground` | `lsst.afw.image.Image`, optional | Background subtracted before pass 1 (``session.bgImage``). |
| `sourcePages` | `int`, optional | Number of per-source pages, brightest sources first. |
| `allSourcePages` | `bool`, optional | Write a page for every source instead (can be hundreds of files). |
| `dpi` | `int`, optional | Resolution of the PNG files. |
| `verbose` | `bool`, optional | Print where the files went and which plots were skipped. |

| Returns | Type | Description |
| --- | --- | --- |
| `files` | `list` [`str`] | The files written. |

### runPlots.saveDetectorPlots

`saveDetectorPlots(visitImage, outDir, catalog=None, session=None, dpi=100, verbose=True)`

Write the detector-level PSF plots: the PSF stars over the detector, their stamps, the PSF.

| Parameter | Type | Description |
| --- | --- | --- |
| `visitImage` | `lsst.afw.image.Exposure` | The visit image. |
| `outDir` | `str` | Directory for the files (created if needed). |
| `catalog` | `lsst.afw.table.SourceCatalog` or `astropy.table.Table`, optional | The visit's source catalogue with ``calib_psf_candidate`` / ``calib_psf_used`` (the stars of the visit PSF). |
| `session` | `lsst.meas.crowd.session.Session`, optional | A session prepared with ``psfSource='refit'`` (the stars of the refitted PSF), used when no catalogue is given. |
| `dpi` | `int`, optional | Resolution of the PNG files. |
| `verbose` | `bool`, optional | Print where the files went. |

| Returns | Type | Description |
| --- | --- | --- |
| `files` | `list` [`str`] | The files written. |

## fastKernels

`lsst.meas.crowd.fastKernels`: fast PSF stamps (C++ with a NumPy fallback).

Fast PSF stamps for a spatially CONSTANT PSF.

| Name | Summary |
| --- | --- |
| [`lanczosWeights`](#fastkernelslanczosweights) | Return normalised Lanczos weights for a sub-pixel shift. |
| [`shiftImage`](#fastkernelsshiftimage) | Shift an image by a sub-pixel offset with a separable Lanczos kernel. |
| [`taperStamp`](#fastkernelstaperstamp) | Tame the wings of a PSF stamp like ``fitLib.psfStamp(taper=...)``. |
| [`renderStamps`](#fastkernelsrenderstamps) | Render stamps of a constant PSF at many positions. |
| [`renderDerivatives`](#fastkernelsrenderderivatives) | Render central-difference derivative stamps dP/dx, dP/dy of a constant PSF. |
| [`FastPsf`](#fastkernelsfastpsf) (class) | Stamps and derivatives of a spatially constant PSF. |

### fastKernels.lanczosWeights

`lanczosWeights(frac, order=5)`

Return normalised Lanczos weights for a sub-pixel shift.

| Parameter | Type | Description |
| --- | --- | --- |
| `frac` | `float` | Shift [px]. |
| `order` | `int`, optional | Lanczos order. |

| Returns | Type | Description |
| --- | --- | --- |
| `first` | `int` | Offset of the first weight: ``out[i] = sum_j w[j] * in[i + first + j]`` approximates ``in(i - frac)``. |
| `weights` | `numpy.ndarray` | The ``2 * order`` weights (summing to 1). |

### fastKernels.shiftImage

`shiftImage(image, dx, dy, order=5, useCxx=None)`

Shift an image by a sub-pixel offset with a separable Lanczos kernel.

Same convention as ``lsst.afw.math.offsetImage``: ``out(x, y) = in(x - dx, y - dy)``; pixels outside the input count as zero.

| Parameter | Type | Description |
| --- | --- | --- |
| `image` | `numpy.ndarray` | Image. |
| `dx, dy` | `float` | Shift [px]. |
| `order` | `int`, optional | Lanczos order. |
| `useCxx` | `bool`, optional | Use the C++ kernels (default: if built). |

| Returns | Type | Description |
| --- | --- | --- |
| `shifted` | `numpy.ndarray` | The shifted image. |

### fastKernels.taperStamp

`taperStamp(stamp, cx, cy, taper, useCxx=None)`

Tame the wings of a PSF stamp like ``fitLib.psfStamp(taper=...)``.

Subtract the median of the outer pixel ring, then apply a raised cosine from ``taper[0] * h`` to ``taper[1] * h`` around (cx, cy), h being half the stamp size.

| Parameter | Type | Description |
| --- | --- | --- |
| `stamp` | `numpy.ndarray` | PSF stamp. |
| `cx, cy` | `float` | Centre in stamp pixel coordinates. |
| `taper` | `tuple` [`float`, `float`] | Inner and outer taper radii as fractions of h. |
| `useCxx` | `bool`, optional | Use the C++ kernels (default: if built). |

| Returns | Type | Description |
| --- | --- | --- |
| `stamp` | `numpy.ndarray` | The tapered stamp. |

### fastKernels.renderStamps

`renderStamps(kernel, xs, ys, order=5, taper=None, useCxx=None)`

Render stamps of a constant PSF at many positions.

Stamp i has origin (round(x) - nx//2, round(y) - ny//2) and is the kernel shifted by the position's sub-pixel offset, like ``psf.computeImage(Point2D(x, y))``.

| Parameter | Type | Description |
| --- | --- | --- |
| `kernel` | `numpy.ndarray` | PSF kernel image (odd size, centred). |
| `xs, ys` | array-like | Positions in parent pixel coordinates. |
| `order` | `int`, optional | Lanczos order. |
| `taper` | `tuple` [`float`, `float`], optional | Wing taper (see `taperStamp`); None = no taper. |
| `useCxx` | `bool`, optional | Use the C++ kernels (default: if built). |

| Returns | Type | Description |
| --- | --- | --- |
| `stamps` | `numpy.ndarray` | Shape (n, ny, nx). |
| `x0, y0` | `numpy.ndarray` [`int`] | Stamp origins. |

### fastKernels.renderDerivatives

`renderDerivatives(kernel, xs, ys, eps=0.25, order=5, useCxx=None)`

Render central-difference derivative stamps dP/dx, dP/dy of a constant PSF.

| Parameter | Type | Description |
| --- | --- | --- |
| `kernel` | `numpy.ndarray` | PSF kernel image (odd size, centred). |
| `xs, ys` | array-like | Positions in parent pixel coordinates. |
| `eps` | `float`, optional | Finite-difference step [px]. |
| `order` | `int`, optional | Lanczos order. |
| `useCxx` | `bool`, optional | Use the C++ kernels (default: if built). |

| Returns | Type | Description |
| --- | --- | --- |
| `dx, dy` | `numpy.ndarray` | Derivative stamps, shape (n, ny, nx), on the centred stamps' origins. |
| `x0, y0` | `numpy.ndarray` [`int`] | Stamp origins. |

### fastKernels.FastPsf

Class.  `FastPsf(kernel, order=5)`

Stamps and derivatives of a spatially constant PSF.

Return conventions match ``fitLib.psfStamp`` and ``fitLib.computePsfDerivatives``, with boxes given as (xmin, ymin, xmax, ymax) tuples.

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
| `taper` | `tuple` [`float`, `float`], optional | Wing taper (see `taperStamp`). |

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

`lsst.meas.crowd.fitLib`: the core library: joint fit, detection, PSF and background tools, plots.

Crowded-field iterative joint PSF fitting: the core library.

| Name | Summary |
| --- | --- |
| [`bboxToSlices`](#fitlibbboxtoslices) | Convert a box into array slices relative to a parent box. |
| [`collapseSaturatedPeaks`](#fitlibcollapsesaturatedpeaks) | Collapse every peak on or near a saturated blob into one peak per blob. |
| [`mergeClosePeaks`](#fitlibmergeclosepeaks) | Merge peaks closer than ``minSeparation``, keeping the brightest pixel's peak. |
| [`jointPsfFit`](#fitlibjointpsffit) | Fit all PSF amplitudes simultaneously with a non-negative least-squares solve. |
| [`iterativeJointDetection`](#fitlibiterativejointdetection) | Detect, fit and subtract iteratively until the residual stops yielding new sources. |
| [`plotTiledComparison`](#fitlibplottiledcomparison) | Plot data, model and residual (and a row profile) for every tile of an image. |
| [`plotDetectionsByIteration`](#fitlibplotdetectionsbyiteration) | Plot every fitted peak over the image, coded by the pass that first found it. |
| [`maskSourcesFromModel`](#fitlibmasksourcesfrommodel) | Flag pixels where the joint-fit model has significant flux as DETECTED. |
| [`refitBackground`](#fitlibrefitbackground) | Re-derive the background ignoring masked pixels and subtract it in place. |
| [`catalogFromPeaks`](#fitlibcatalogfrompeaks) | Build a measured SourceCatalog from raw (x, y) peak positions. |
| [`redetectAndRefitPsf`](#fitlibredetectandrefitpsf) | Detect, deblend, measure, select stars and determine a new PSF. |
| [`outerLoopPsfBackgroundRefinement`](#fitlibouterlooppsfbackgroundrefinement) | Alternate inner-loop fits with background and PSF re-derivation. |
| [`plotPsfModelsByOuterPass`](#fitlibplotpsfmodelsbyouterpass) | Plot the PSF kernel image used in each outer pass. |
| [`plotBackgroundByOuterPass`](#fitlibplotbackgroundbyouterpass) | Plot the background re-derived going into each outer pass. |
| [`plotSourcesOverlaidByPass`](#fitlibplotsourcesoverlaidbypass) | Plot the data of each outer pass with its converged peak positions. |
| [`plotOuterLoopSummary`](#fitlibplotouterloopsummary) | Plot the residual standard deviation and the source count against outer pass. |
| [`plotPsfFwhmByOuterPass`](#fitlibplotpsffwhmbyouterpass) | Plot the FWHM of the PSF used in each outer pass. |
| [`StepTimer`](#fitlibsteptimer) (class) | Collect timing records of pipeline steps. |
| [`plotTimingSummary`](#fitlibplottimingsummary) | Plot the total time per step and the time per inner iteration. |
| [`makeCutout`](#fitlibmakecutout) | Cut a sub-exposure out of a crop, in coordinates relative to the crop's corner. |
| [`psfStamp`](#fitlibpsfstamp) | Render the PSF at (x, y), clipped to a box, with its far wings tamed. |
| [`dropMaskedPeaks`](#fitlibdropmaskedpeaks) | Drop peaks on EDGE / NO_DATA pixels. |
| [`dropOutOfBounds`](#fitlibdropoutofbounds) | Drop peaks outside the exposure's box grown by ``margin`` px. |
| [`effectiveVariance`](#fitlibeffectivevariance) | Return the variance used to weight the fit. |
| [`rejectNearBright`](#fitlibrejectnearbright) | Drop new peaks within ``minSep`` px of an existing bright source. |
| [`pruneAndMerge`](#fitlibpruneandmerge) | Drop low-S/N sources, then merge pairs closer than ``minSep`` px. |
| [`empiricalNoiseScale`](#fitlibempiricalnoisescale) | Return the empirical noise scale of a fit's residual. |
| [`ampSnrFull`](#fitlibampsnrfull) | Return honest amplitude S/N values, including the covariance of overlapping sources. |
| [`eliminateWeak`](#fitlibeliminateweak) | Remove sources below ``minSnr`` by backward elimination. |
| [`shapeMetrics`](#fitlibshapemetrics) | Measure how PSF-like each source is. |
| [`judgeShape`](#fitlibjudgeshape) | Decide which sources fail the PSF-shape tests. |
| [`eliminateNonpsf`](#fitlibeliminatenonpsf) | Remove sources that fail the PSF-shape test by backward elimination. |
| [`plotShapeRejects`](#fitlibplotshaperejects) | Plot the sources removed by the shape cut over the data. |
| [`selectSourcePages`](#fitlibselectsourcepages) | Choose the sources that get per-source pages. |
| [`sourceColors`](#fitlibsourcecolors) | Assign each source on the pages a colour from ``SRC_PALETTE``. |
| [`jointPsfFitDiag`](#fitlibjointpsffitdiag) | Fit all PSF amplitudes in one non-negative solve, keeping what is needed for diagnostics. |
| [`computePsfDerivatives`](#fitlibcomputepsfderivatives) | Compute central-difference derivative stamps dP/dx and dP/dy. |
| [`windowedCentroid`](#fitlibwindowedcentroid) | Return the flux-weighted centroid in a box around (x, y). |
| [`refinePeakPositions`](#fitlibrefinepeakpositions) | Move each well-detected source towards the position that minimises its residual. |
| [`redetectAndRefitPsfAdaptive`](#fitlibredetectandrefitpsfadaptive) | Run `redetectAndRefitPsf`, trying kernel sizes from largest to smallest. |
| [`plotIterationOverview`](#fitlibplotiterationoverview) | Plot the overview of one inner iteration. |
| [`runInnerLoopVerbose`](#fitlibruninnerloopverbose) | Run the instrumented inner loop (detect, fit, refine, prune, repeat). |
| [`plotAmplitudeEvolution`](#fitlibplotamplitudeevolution) | Plot the amplitude of every final source against inner iteration. |
| [`plotSourceEvolution`](#fitlibplotsourceevolution) | Plot model and residual cross-sections of three sources at every iteration. |
| [`plotCutoutSetup`](#fitlibplotcutoutsetup) | Plot the cutout before fitting: data, sigma map, mask planes and PSF. |
| [`plotOuterStepBackground`](#fitlibplotouterstepbackground) | Plot the background refit of one outer step. |
| [`plotPsfCandidateStamps`](#fitlibplotpsfcandidatestamps) | Plot the PSF candidate stamps, used ones first. |
| [`plotOuterStepPsf`](#fitlibplotoutersteppsf) | Plot the PSF refit of one outer step: star selection and old vs new PSF. |
| [`highfreqNoiseCheckArrays`](#fitlibhighfreqnoisecheckarrays) | Compare the residual's pixel-to-pixel scatter with the variance plane. |
| [`highfreqNoiseCheck`](#fitlibhighfreqnoisecheck) | Run `highfreqNoiseCheckArrays` on one snapshot of a run. |
| [`estimatePsfCorrection`](#fitlibestimatepsfcorrection) | Estimate an improved PSF kernel from neighbour-subtracted residuals. |
| [`kernelToLsstPsf`](#fitlibkerneltolsstpsf) | Wrap a kernel array as a spatially constant LSST PSF. |
| [`extendKernelWithMoffat`](#fitlibextendkernelwithmoffat) | Replace the noisy outer part of a PSF kernel by an analytic Moffat tail. |
| [`extendPsfWings`](#fitlibextendpsfwings) | Return a constant KernelPsf with a Moffat tail, taken from any PSF at a point. |
| [`refinePsfFromFit`](#fitlibrefinepsffromfit) | Make one PSF update step from a finished run. |
| [`polishPsfFixedSources`](#fitlibpolishpsffixedsources) | Iterate PSF correction and re-solve without re-detecting sources. |
| [`plotPsfCorrection`](#fitlibplotpsfcorrection) | Plot the old kernel, the proposed correction and the radial profiles. |
| [`runCutoutPipeline`](#fitlibruncutoutpipeline) | Run the outer loop (inner fit, background and PSF refinement) on a cutout. |
| [`plotWingDiagnostic`](#fitlibplotwingdiagnostic) | Plot the azimuthal profiles of the brightest sources: data vs their own model. |
| [`plotSourceCutouts`](#fitlibplotsourcecutouts) | Plot cutouts and cross-sections of individual sources of a run. |
| [`negativeSourceCheck`](#fitlibnegativesourcecheck) | Estimate the number of false detections from the residual itself. |
| [`sourceSnrHistogram`](#fitlibsourcesnrhistogram) | Plot the histogram of the honest S/N of every final source. |
| [`drawEdgeMargin`](#fitlibdrawedgemargin) | Outline and shade the band within ``marginPx`` of the image border. |
| [`plotFinalComparison`](#fitlibplotfinalcomparison) | Plot original \| final model \| residual \| residual/sigma over the fitted region. |
| [`plotSourcePages`](#fitlibplotsourcepages) | Plot one row per source: data \| model \| residual \| x-section \| y-section. |
| [`sourceMetrics`](#fitlibsourcemetrics) | Return per-source quantities used to stratify the sample. |
| [`selectRepresentativeSources`](#fitlibselectrepresentativesources) | Select a representative sample of sources per category. |
| [`plotRepresentativeSources`](#fitlibplotrepresentativesources) | Select a representative sample and plot it with `plotSourcePages`. |

### fitLib.bboxToSlices

`bboxToSlices(bbox, parentBbox)`

Convert a box into array slices relative to a parent box.

| Parameter | Type | Description |
| --- | --- | --- |
| `bbox` | `lsst.geom.Box2I` | Box to convert. |
| `parentBbox` | `lsst.geom.Box2I` | Box of the array the slices index (its origin is pixel [0, 0]). |

| Returns | Type | Description |
| --- | --- | --- |
| `slices` | `tuple` [`slice`, `slice`] | (y slice, x slice) selecting ``bbox`` from an array covering ``parentBbox``. |

### fitLib.collapseSaturatedPeaks

`collapseSaturatedPeaks(peakPositions, exposure, satBitmask=None, dilateIters=2)`

Collapse every peak on or near a saturated blob into one peak per blob.

The representative position is a flux-weighted centroid from the unsaturated ring just outside the saturated core (the flat top carries no usable position or flux information). Without this, noise on the saturated flat top produces many near-duplicate peaks whose PSF stamps are almost identical; the joint fit then splits the true amplitude between them arbitrarily, which shows up as a blocky, high-frequency patch in the model image.

| Parameter | Type | Description |
| --- | --- | --- |
| `peakPositions` | `list` [`tuple` [`float`, `float`]] | Peak (x, y) positions in parent pixel coordinates. |
| `exposure` | `lsst.afw.image.Exposure` | Exposure the peaks were detected on. |
| `satBitmask` | `int`, optional | Mask bits of saturated pixels (default: the SAT plane). |
| `dilateIters` | `int`, optional | Number of dilations of the saturated region when assigning peaks to a blob. |

| Returns | Type | Description |
| --- | --- | --- |
| `peaks` | `list` [`tuple` [`float`, `float`]] | Peaks with every saturated blob reduced to one representative peak. |
| `nRemoved` | `int` | Number of peaks removed. |

### fitLib.mergeClosePeaks

`mergeClosePeaks(peakPositions, exposure, minSeparation=None)`

Merge peaks closer than ``minSeparation``, keeping the brightest pixel's peak.

Prevents near-duplicate columns in the design matrix (mostly from re-detecting the same wings) from destabilising the joint solve.

| Parameter | Type | Description |
| --- | --- | --- |
| `peakPositions` | `list` [`tuple` [`float`, `float`]] | Peak (x, y) positions in parent pixel coordinates. |
| `exposure` | `lsst.afw.image.Exposure` | Exposure whose pixel values rank the peaks. |
| `minSeparation` | `float`, optional | Merge radius [px] (default: twice the PSF sigma). |

| Returns | Type | Description |
| --- | --- | --- |
| `mergedPositions` | `list` [`tuple` [`float`, `float`]] | Surviving peak positions. |

### fitLib.jointPsfFit

`jointPsfFit(exposure, psf, peakPositions, badBitmask, ampZeroTol=1e-06)`

Fit all PSF amplitudes simultaneously with a non-negative least-squares solve.

One sparse, weighted, bound-constrained solve against the original exposure.  Unlike an unconstrained solve followed by clipping, ``lsq_linear`` with bounds (0, inf) lets the other amplitudes re-adjust when one source is pushed to zero, so overlapping sources get a real non-negative amplitude instead of one of them being dropped from the model.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Image to fit (image, variance and mask planes are used). |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `peakPositions` | `list` [`tuple` [`float`, `float`]] | Source (x, y) positions in parent pixel coordinates. |
| `badBitmask` | `int` | Mask bits of pixels to ignore. |
| `ampZeroTol` | `float`, optional | Amplitudes at or below this count as zero. |

| Returns | Type | Description |
| --- | --- | --- |
| `amps` | `numpy.ndarray` | Fitted amplitudes (fluxes), one per position. |
| `modelFull` | `numpy.ndarray` | Model image over the exposure's box. |
| `nAtZeroBound` | `int` | Number of sources that converged to (numerically) zero amplitude -- sources that still fail to fit rather than being absent. |

### fitLib.iterativeJointDetection

`iterativeJointDetection(exposure, psf, badBitmask, nIterations=10, thresholdStart=8.0, thresholdStep=-1.5, thresholdMin=1.5, minPixels=1, minNewPeaks=2, convergenceFrac=0.02, verbose=True)`

Detect, fit and subtract iteratively until the residual stops yielding new sources.

Each pass: (1) detect on the current residual with a decreasing threshold, so fainter sources are searched for once brighter ones are gone; (2) collapse multiple peaks on saturated cores; (3) merge the new peaks with all peaks found so far; (4) refit all of them jointly against the original exposure with a non-negative solve; (5) rebuild the residual. Every pass is a full re-solve, so errors from an early pass are not compounded.  The number of zero-amplitude sources is printed each pass: if it is large or not shrinking, real sources are being fitted down to nothing.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Image to fit. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `badBitmask` | `int` | Mask bits of pixels to ignore. |
| `nIterations` | `int`, optional | Maximum number of passes. |
| `thresholdStart` | `float`, optional | Detection threshold [sigma] of the first pass. |
| `thresholdStep` | `float`, optional | Change of the threshold per pass (negative: decreasing). |
| `thresholdMin` | `float`, optional | Lowest detection threshold [sigma]. |
| `minPixels` | `int`, optional | Minimum footprint size of a detection. |
| `minNewPeaks` | `int`, optional | Stop when a pass finds fewer new peaks than this. |
| `convergenceFrac` | `float`, optional | Stop when the residual standard deviation improves by less than this fraction. |
| `verbose` | `bool`, optional | Print progress. |

| Returns | Type | Description |
| --- | --- | --- |
| `result` | `dict` | ``modelFull``, ``residualFull``, ``peaks``, ``amps``, the per-pass histories (``residStdHistory``, ``nPeaksHistory``, ``nZeroAmpHistory``) and ``firstSeen`` (pass index at which each peak, keyed by its position rounded to 0.1 px, was found). |

### fitLib.plotTiledComparison

`plotTiledComparison(dataFull, modelFull, residualFull, nSources=None, nPasses=None, tileSize=100, cmapImg=None, globalScale=True, skipEmpty=True, emptyStdThresh=1e-06, figsizePerPanel=5, showRowProfile=True, rowSelect='maxFlux', stretch='log')`

Plot data, model and residual (and a row profile) for every tile of an image.

The data and model panels share one colour scale when ``globalScale`` is True; the residual panel always scales to its own tile.  With ``showRowProfile`` a fourth panel shows data and model along one row of the tile (marked with a dashed cyan line), so the model can be compared with the data through a crowded section.

| Parameter | Type | Description |
| --- | --- | --- |
| `dataFull, modelFull, residualFull` | `numpy.ndarray` | Data, model and residual images (same shape). |
| `nSources` | `int`, optional | Number of sources, for the title. |
| `nPasses` | `int`, optional | Number of detection passes, for the title. |
| `tileSize` | `int`, optional | Tile size [px]. |
| `cmapImg` | `str`, optional | Colour map of the data and model panels (default: 'mako' if available, else 'viridis'). |
| `globalScale` | `bool`, optional | One colour scale for all data / model panels instead of one per tile. |
| `skipEmpty` | `bool`, optional | Skip tiles whose data standard deviation is below ``emptyStdThresh``. |
| `emptyStdThresh` | `float`, optional | See ``skipEmpty``. |
| `figsizePerPanel` | `float`, optional | Size of one panel [inches]. |
| `showRowProfile` | `bool`, optional | Add the row-profile panel. |
| `rowSelect` | `str`, `int` or callable, optional | Row of each tile to profile: 'maxFlux' (the row with the most flux), a row index local to the tile, or a function ``rowSelect(tileData) -> row index``. |
| `stretch` | `str`, optional | Stretch of the data / model colour scale. |

### fitLib.plotDetectionsByIteration

`plotDetectionsByIteration(exposureCrop, resultIter)`

Plot every fitted peak over the image, coded by the pass that first found it.

Peaks whose final amplitude is (numerically) zero -- they hit the non-negativity bound -- are marked as failed fits with a red 'x'; they are the ones to inspect if sources seem to be missing from the model.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposureCrop` | `lsst.afw.image.Exposure` | Image the pipeline ran on. |
| `resultIter` | `dict` | Output of `iterativeJointDetection`. |

| Returns | Type | Description |
| --- | --- | --- |
| `sourcesTbl` | `astropy.table.Table` | Positions, amplitudes and detection pass of the sources. |

### fitLib.maskSourcesFromModel

`maskSourcesFromModel(exposure, modelFull, noiseSigma=None, nSigma=4.0, dilateIters=1, verbose=True)`

Flag pixels where the joint-fit model has significant flux as DETECTED.

The next background fit then treats them as source, not sky.  The joint model contains all sources found so far, including ones too faint for the first background pass.  A loose threshold flags most of the image and starves the background fit; if the printed masked fraction climbs much above ~50 %, tighten ``nSigma`` or ``dilateIters``.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Exposure whose mask is modified in place. |
| `modelFull` | `numpy.ndarray` | Model image over the exposure's box. |
| `noiseSigma` | `float`, optional | Pixel noise (default: sqrt of the median variance). |
| `nSigma` | `float`, optional | Model pixels above ``nSigma * noiseSigma`` are flagged. |
| `dilateIters` | `int`, optional | Number of dilations of the flagged region. |
| `verbose` | `bool`, optional | Print the masked fraction. |

| Returns | Type | Description |
| --- | --- | --- |
| `srcMask` | `numpy.ndarray` [`bool`] | The flagged pixels. |

### fitLib.refitBackground

`refitBackground(exposure, algorithm='AKIMA_SPLINE', binSize=32, useApprox=True, statsMask=('BAD', 'SAT', 'EDGE', 'NO_DATA', 'DETECTED', 'CR', 'INTRP'))`

Re-derive the background ignoring masked pixels and subtract it in place.

Pass a fresh clone of the raw (never background-subtracted) exposure, already flagged by `maskSourcesFromModel`, or the background is subtracted twice.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Exposure to background-subtract in place. |
| `algorithm` | `str`, optional | Interpolation algorithm (only used when ``useApprox`` is False). |
| `binSize` | `int`, optional | Background bin size [px]. |
| `useApprox` | `bool`, optional | Fit a smooth Chebyshev surface to the binned statistics instead of interpolating; more robust when many pixels are masked. |
| `statsMask` | `tuple` [`str`], optional | Mask planes ignored when computing the background statistics. |

| Returns | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The background-subtracted exposure (same object). |
| `background` | `lsst.afw.image.Image` | The subtracted background image. |

### fitLib.catalogFromPeaks

`catalogFromPeaks(exposure, peakPositions, schema, measureTask, stampRadius=15, excludePositions=None, excludeRadius=5.0, verbose=True)`

Build a measured SourceCatalog from raw (x, y) peak positions.

Each peak gets a circular footprint and one Peak record, and the catalogue is measured with the same ``measureTask`` as fresh detections, so the fields are comparable and usable by the star selector.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Exposure to measure on. |
| `peakPositions` | `list` [`tuple` [`float`, `float`]] | Peak positions in parent pixel coordinates. |
| `schema` | `lsst.afw.table.Schema` | The schema used to build ``measureTask``. |
| `measureTask` | `lsst.meas.base.SingleFrameMeasurementTask` | Measurement task. |
| `stampRadius` | `int`, optional | Footprint radius [px]. |
| `excludePositions` | `list` [`tuple` [`float`, `float`]], optional | Skip peaks within ``excludeRadius`` of these positions (e.g. existing detections). |
| `excludeRadius` | `float`, optional | See ``excludePositions`` [px]. |
| `verbose` | `bool`, optional | Print how many peaks were skipped. |

| Returns | Type | Description |
| --- | --- | --- |
| `catalog` | `lsst.afw.table.SourceCatalog` | The measured catalogue. |

### fitLib.redetectAndRefitPsf

`redetectAndRefitPsf(exposure, detThreshold=3.0, detMinPixels=2, detLocalBgBin=16, detBgBin=32, starSelectorConfig=None, psfDeterminerConfig=None, kernelSize=41, extraPeakPositions=None, forceIncludeStarPositions=None, peakStampRadius=15, verbose=True)`

Detect, deblend, measure, select stars and determine a new PSF.

Run in each outer pass, so it works on the newly cleaned background.  Extra peak positions (e.g. the previous pass's converged catalogue) are added before star selection; forced star positions are appended after selection, bypassing its cuts.  Candidates with failed stamps are filtered out before the PSF determiner runs (even one can make it reject the whole set).  If too few candidates survive, the exposure's current PSF is kept and a "NOT updated" note is printed.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Background-subtracted exposure. |
| `detThreshold` | `float`, optional | Detection threshold [sigma]. |
| `detMinPixels` | `int`, optional | Minimum footprint size. |
| `detLocalBgBin` | `int`, optional | Bin size of the temporary local background [px]. |
| `detBgBin` | `int`, optional | Bin size of the detection background [px]. |
| `starSelectorConfig` | `lsst.meas.algorithms.ObjectSizeStarSelectorConfig`, optional | Star selector configuration. |
| `psfDeterminerConfig` | `lsst.meas.algorithms.PcaPsfDeterminerConfig`, optional | PSF determiner configuration. |
| `kernelSize` | `int`, optional | PSF kernel size [px], used for both the candidates and the determiner. |
| `extraPeakPositions` | `list` [`tuple` [`float`, `float`]], optional | Additional peaks measured and offered to the star selector. |
| `forceIncludeStarPositions` | `list` [`tuple` [`float`, `float`]], optional | Star positions appended to the selected stars regardless of the selector's cuts. |
| `peakStampRadius` | `int`, optional | Footprint radius of the extra / forced peaks [px]. |
| `verbose` | `bool`, optional | Print progress. |

| Returns | Type | Description |
| --- | --- | --- |
| `psfNew` | `lsst.afw.detection.Psf` | The new PSF (or the current one if it could not be determined). |
| `sources` | `lsst.afw.table.SourceCatalog` | Measured detections. |
| `starCat` | `lsst.afw.table.SourceCatalog` | Selected PSF stars. |
| `cellSet` | `lsst.afw.math.SpatialCellSet` | PSF candidates used by the determiner. |

### fitLib.outerLoopPsfBackgroundRefinement

`outerLoopPsfBackgroundRefinement(rawExposureCrop, psfInitial, badBitmask, initialStarCat=None, nOuterIterations=3, nInnerIterations=10, innerKwargs=None, outerConvergenceFrac=0.02, bgUseApprox=True, psfKernelSize=41, verbose=True)`

Alternate inner-loop fits with background and PSF re-derivation.

Pass 0 runs `iterativeJointDetection` with the current background and PSF.  Pass k > 0 masks every known source (from pass k-1's model) on a fresh clone of the raw image, re-derives the background, re-determines the PSF (offering pass k-1's catalogue and forcing in ``initialStarCat``) and reruns the inner loop.

| Parameter | Type | Description |
| --- | --- | --- |
| `rawExposureCrop` | `lsst.afw.image.Exposure` | Crop before any background subtraction. |
| `psfInitial` | `lsst.afw.detection.Psf` | PSF of the first pass. |
| `badBitmask` | `int` | Mask bits of pixels to ignore. |
| `initialStarCat` | `lsst.afw.table.SourceCatalog`, optional | Star selector output whose positions are forced into every PSF determination. |
| `nOuterIterations` | `int`, optional | Maximum number of outer passes. |
| `nInnerIterations` | `int`, optional | Maximum number of inner passes (when ``innerKwargs`` is not given). |
| `innerKwargs` | `dict`, optional | Keyword arguments of `iterativeJointDetection`. |
| `outerConvergenceFrac` | `float`, optional | Stop when the residual standard deviation improves by less than this fraction. |
| `bgUseApprox` | `bool`, optional | See `refitBackground`. |
| `psfKernelSize` | `int`, optional | PSF kernel size [px]. |
| `verbose` | `bool`, optional | Print progress. |

| Returns | Type | Description |
| --- | --- | --- |
| `outerHistory` | `list` [`dict`] | One dict per outer pass: ``psf``, ``exposure``, ``result`` (inner-loop output), ``residStd``, ``nPeaks``. |

### fitLib.plotPsfModelsByOuterPass

`plotPsfModelsByOuterPass(outerHistory, exposureCenter=None)`

Plot the PSF kernel image used in each outer pass.

| Parameter | Type | Description |
| --- | --- | --- |
| `outerHistory` | `list` [`dict`] | Outer-loop history (each entry has ``psf`` and ``exposure``). |
| `exposureCenter` | `lsst.geom.Point2D`, optional | Position at which to render the kernel (default: the exposure centre). |

### fitLib.plotBackgroundByOuterPass

`plotBackgroundByOuterPass(outerHistory)`

Plot the background re-derived going into each outer pass.

The first pass uses the background subtracted before the pipeline ran and is skipped.

| Parameter | Type | Description |
| --- | --- | --- |
| `outerHistory` | `list` [`dict`] | Outer-loop history (entries may have ``bgImage``). |

### fitLib.plotSourcesOverlaidByPass

`plotSourcesOverlaidByPass(outerHistory, markerSize=25, vminPct=5, vmaxPct=99)`

Plot the data of each outer pass with its converged peak positions.

| Parameter | Type | Description |
| --- | --- | --- |
| `outerHistory` | `list` [`dict`] | Outer-loop history. |
| `markerSize` | `float`, optional | Marker size. |
| `vminPct, vmaxPct` | `float`, optional | Percentiles of the display range. |

### fitLib.plotOuterLoopSummary

`plotOuterLoopSummary(outerHistory)`

Plot the residual standard deviation and the source count against outer pass.

| Parameter | Type | Description |
| --- | --- | --- |
| `outerHistory` | `list` [`dict`] | Outer-loop history (entries have ``residStd`` and ``nPeaks``). |

### fitLib.plotPsfFwhmByOuterPass

`plotPsfFwhmByOuterPass(outerHistory, exposureCenter=None)`

Plot the FWHM of the PSF used in each outer pass.

| Parameter | Type | Description |
| --- | --- | --- |
| `outerHistory` | `list` [`dict`] | Outer-loop history. |
| `exposureCenter` | `lsst.geom.Point2D`, optional | Position at which to evaluate the PSF (default: the exposure centre). |

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

### fitLib.plotTimingSummary

`plotTimingSummary(timer=TIMER)`

Plot the total time per step and the time per inner iteration.

| Parameter | Type | Description |
| --- | --- | --- |
| `timer` | `StepTimer`, optional | Timer to summarise (default: the module-level ``TIMER``). |

### fitLib.makeCutout

`makeCutout(exposureCrop, xRange=(0, 100), yRange=(0, 100))`

Cut a sub-exposure out of a crop, in coordinates relative to the crop's corner.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposureCrop` | `lsst.afw.image.Exposure` | The crop. |
| `xRange, yRange` | `tuple` [`int`, `int`], optional | Pixel ranges relative to the crop's origin (end exclusive). |

| Returns | Type | Description |
| --- | --- | --- |
| `cutout` | `lsst.afw.image.Exposure` | A deep copy of the sub-exposure. |

Raises: `TypeError` (Raised if ``exposureCrop`` is not an Exposure.)

### fitLib.psfStamp

`psfStamp(psf, x, y, clipBbox, taper=None)`

Render the PSF at (x, y), clipped to a box, with its far wings tamed.

With ``taper`` the median of the outermost pixel ring (edge pedestal) is subtracted and the stamp is tapered radially to zero, so each source ends smoothly in a disc instead of a hard-edged square.

| Parameter | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `x, y` | `float` | Position in parent pixel coordinates. |
| `clipBbox` | `lsst.geom.Box2I` | Box to clip the stamp to. |
| `taper` | `tuple` [`float`, `float`], optional | Inner and outer taper radii as fractions of the stamp half-size. |

| Returns | Type | Description |
| --- | --- | --- |
| `stamp` | `numpy.ndarray` or `None` | The clipped stamp, or None if it misses ``clipBbox``. |
| `bbox` | `lsst.geom.Box2I` | Box of the clipped stamp. |

### fitLib.dropMaskedPeaks

`dropMaskedPeaks(peaks, exposure, planes=('EDGE', 'NO_DATA'))`

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
| `nDropped` | `int` | Number of peaks dropped. |

### fitLib.dropOutOfBounds

`dropOutOfBounds(peaks, exposure, margin=0.0)`

Drop peaks outside the exposure's box grown by ``margin`` px.

| Parameter | Type | Description |
| --- | --- | --- |
| `peaks` | `list` [`tuple` [`float`, `float`]] | Peak positions in parent pixel coordinates. |
| `exposure` | `lsst.afw.image.Exposure` | Exposure providing the box. |
| `margin` | `float`, optional | Growth of the box [px]; a negative margin drops peaks within ``\|margin\|`` px of the edge. |

| Returns | Type | Description |
| --- | --- | --- |
| `keep` | `list` [`tuple` [`float`, `float`]] | Remaining peaks. |
| `nDropped` | `int` | Number of peaks dropped. |

### fitLib.effectiveVariance

`effectiveVariance(exposure, model=None, frac=0.03, poissonK=0.0)`

Return the variance used to weight the fit.

``variance + (frac * model)**2 + poissonK * model``: the pixel variance plus a fractional model-error term and missing photon noise.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Exposure providing the variance plane. |
| `model` | `numpy.ndarray`, optional | Current model image (no extra terms if None). |
| `frac` | `float`, optional | Fractional model error. |
| `poissonK` | `float`, optional | Extra variance per unit model flux. |

| Returns | Type | Description |
| --- | --- | --- |
| `var` | `numpy.ndarray` | The effective variance. |

### fitLib.rejectNearBright

`rejectNearBright(newPeaks, oldPeaks, oldIsBright, minSep)`

Drop new peaks within ``minSep`` px of an existing bright source.

| Parameter | Type | Description |
| --- | --- | --- |
| `newPeaks` | `list` [`tuple` [`float`, `float`]] | Candidate peaks. |
| `oldPeaks` | array-like | Existing source positions, shape (N, 2). |
| `oldIsBright` | array-like [`bool`] | Which existing sources are bright. |
| `minSep` | `float` | Rejection radius [px]. |

| Returns | Type | Description |
| --- | --- | --- |
| `keep` | `list` [`tuple` [`float`, `float`]] | Remaining candidates. |
| `n_rejected` | `int` | Number of candidates rejected. |

### fitLib.pruneAndMerge

`pruneAndMerge(fit, peaks, minSnr=5.0, minSep=3.0, brightFlags=None, nearBrightSep=None, nearBrightSnr=8.0)`

Drop low-S/N sources, then merge pairs closer than ``minSep`` px.

The brighter source of a pair survives at the amplitude-weighted position.  Sources within ``nearBrightSep`` px of a bright source must reach ``nearBrightSnr`` instead of ``minSnr``.  Uses the cheap diagonal-only S/N; the honest S/N is applied by `eliminateWeak`.

| Parameter | Type | Description |
| --- | --- | --- |
| `fit` | `dict` | Joint fit (from `jointPsfFitDiag`) of ``peaks``. |
| `peaks` | array-like | Source positions, shape (N, 2). |
| `minSnr` | `float`, optional | Minimum S/N. |
| `minSep` | `float`, optional | Merge radius [px]. |
| `brightFlags` | array-like [`bool`], optional | Which sources are bright. |
| `nearBrightSep` | `float`, optional | Radius around bright sources with the stricter S/N cut [px]. |
| `nearBrightSnr` | `float`, optional | S/N required near bright sources. |

| Returns | Type | Description |
| --- | --- | --- |
| `peaks` | `list` [`tuple` [`float`, `float`]] | Surviving source positions. |

### fitLib.empiricalNoiseScale

`empiricalNoiseScale(fit, exposure, modelFrac=0.01, matched=True, weight=None)`

Return the empirical noise scale of a fit's residual.

| Parameter | Type | Description |
| --- | --- | --- |
| `fit` | `dict` | Joint fit (from `jointPsfFitDiag`). |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `modelFrac` | `float`, optional | Pixels with model below this fraction of the model maximum count as sky. |
| `matched` | `bool`, optional | True: scale for S/N values -- the pixel scale moved ``weight`` of the way towards the matched-filter scale when that is larger.  False: the pixel scale only (for chi2). |
| `weight` | `float`, optional | See ``matched`` (default ``MATCHED_NOISE_WEIGHT``). |

| Returns | Type | Description |
| --- | --- | --- |
| `scale` | `float` | Noise scale (at least 1). |

### fitLib.ampSnrFull

`ampSnrFull(fit, noiseScale=1.0)`

Return honest amplitude S/N values, including the covariance of overlapping sources.

``S/N = amp / sqrt(diag(N^-1))`` with ``N = A^T W A`` over all sources, divided by the empirical noise scale.

| Parameter | Type | Description |
| --- | --- | --- |
| `fit` | `dict` | Joint fit (from `jointPsfFitDiag`). |
| `noiseScale` | `float`, optional | Empirical noise scale dividing the S/N. |

| Returns | Type | Description |
| --- | --- | --- |
| `snr` | `numpy.ndarray` | S/N per source. |

### fitLib.eliminateWeak

`eliminateWeak(exposure, psf, peaks, fit, badBitmask, minSnr=5.0, dropFrac=0.25, maxRounds=30, **fitkw)`

Remove sources below ``minSnr`` by backward elimination.

Each round drops the weakest ~``dropFrac`` of the sub-threshold sources (at least one) and re-solves, until every source has an honest S/N of at least ``minSnr``.  The last source is never removed.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `peaks` | `list` [`tuple` [`float`, `float`]] | Source positions; ``fit`` must correspond to them. |
| `fit` | `dict` | Joint fit of ``peaks``. |
| `badBitmask` | `int` | Mask bits of pixels to ignore. |
| `minSnr` | `float`, optional | Minimum honest S/N. |
| `dropFrac` | `float`, optional | Fraction of the failing sources dropped per round. |
| `maxRounds` | `int`, optional | Maximum number of rounds. |
| `**fitkw` | – | Keyword arguments of `jointPsfFitDiag`. |

| Returns | Type | Description |
| --- | --- | --- |
| `peaks` | `list` [`tuple` [`float`, `float`]] | Surviving source positions. |
| `fit` | `dict` | Joint fit of the survivors. |

### fitLib.shapeMetrics

`shapeMetrics(exposure, psf, peaks, fit, noiseScale=1.0, chi2Scale=None, bgSigma=4.0, widthScales=(0.6, 0.75, 1.0, 1.3, 1.7), minPix=15)`

Measure how PSF-like each source is.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `peaks` | `list` [`tuple` [`float`, `float`]] | Source positions; ``fit`` must correspond to them (same order). |
| `fit` | `dict` | Joint fit of ``peaks``. |
| `noiseScale` | `float`, optional | Divides the S/N values (matched-filter calibrated). |
| `chi2Scale` | `float`, optional | Divides chi2 (pixel-level calibrated; default ``noiseScale``). |
| `bgSigma` | `float`, optional | Fit radius in units of the PSF sigma. |
| `widthScales` | `tuple` [`float`], optional | Trial PSF width scale factors. |
| `minPix` | `int`, optional | Minimum number of good pixels to judge a source. |

| Returns | Type | Description |
| --- | --- | --- |
| `metrics` | `dict` [`str`, `numpy.ndarray`] | ``snr`` (local, noise-scaled), ``chi2red``, ``bestScale``, ``dchi2Width`` (> 0: the data prefer ``bestScale`` over the PSF width) and ``contrast``; NaN = not judged. |

### fitLib.judgeShape

`judgeShape(m, minSnrLocal=5.0, chi2Max=2.0, widthRange=(0.75, 1.35), dchi2Width=9.0, trustSnr=30.0, minContrast=20.0, weakSnr=10.0, dchi2WidthWeak=4.0)`

Decide which sources fail the PSF-shape tests.

Sources with local S/N of at least ``trustSnr`` (including saturated ones) skip the chi2, width and contrast tests.  Untrusted sources with S/N below ``weakSnr`` use the looser width threshold ``dchi2WidthWeak``.

| Parameter | Type | Description |
| --- | --- | --- |
| `m` | `dict` [`str`, `numpy.ndarray`] | Output of `shapeMetrics`. |
| `minSnrLocal` | `float`, optional | Minimum local S/N. |
| `chi2Max` | `float`, optional | Maximum reduced chi2. |
| `widthRange` | `tuple` [`float`, `float`], optional | Allowed range of the best-fitting width scale. |
| `dchi2Width` | `float`, optional | Improvement in chi2 needed to reject on width. |
| `trustSnr` | `float`, optional | S/N above which a source is trusted. |
| `minContrast` | `float`, optional | Minimum ratio of the source's peak to the neighbours' model flux at its position (None or 0 = off); 20 ~ (5 % wing error) x (4x margin). |
| `weakSnr` | `float`, optional | S/N below which ``dchi2WidthWeak`` applies. |
| `dchi2WidthWeak` | `float`, optional | Width threshold for weak sources. |

| Returns | Type | Description |
| --- | --- | --- |
| `fail` | `numpy.ndarray` [`bool`] | Sources failing any test. |
| `score` | `numpy.ndarray` | Badness score (higher = worse); orders the elimination. |
| `flags` | `dict` [`str`, `numpy.ndarray`] | Per-test failures: ``snr``, ``chi2``, ``width``, ``contrast``. |

### fitLib.eliminateNonpsf

`eliminateNonpsf(exposure, psf, peaks, fit, badBitmask, shapeKwargs=None, dropFrac=0.25, maxRounds=15, verbose=True, **fitkw)`

Remove sources that fail the PSF-shape test by backward elimination.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `peaks` | `list` [`tuple` [`float`, `float`]] | Source positions; ``fit`` must correspond to them. |
| `fit` | `dict` | Joint fit of ``peaks``. |
| `badBitmask` | `int` | Mask bits of pixels to ignore. |
| `shapeKwargs` | `dict`, optional | Keyword arguments of `shapeMetrics` and `judgeShape`. |
| `dropFrac` | `float`, optional | Fraction of the failing sources dropped per round. |
| `maxRounds` | `int`, optional | Maximum number of rounds. |
| `verbose` | `bool`, optional | Print the rejected sources. |
| `**fitkw` | – | Keyword arguments of `jointPsfFitDiag`. |

| Returns | Type | Description |
| --- | --- | --- |
| `peaks` | `list` [`tuple` [`float`, `float`]] | Surviving source positions. |
| `fit` | `dict` | Joint fit of the survivors. |
| `rejected` | `list` [`tuple` [`float`, `float`, `float`]] | (x, y, score) of the removed sources. |

### fitLib.plotShapeRejects

`plotShapeRejects(exposure, rejectedByIter, peaks=None, title='sources removed by shape cut')`

Plot the sources removed by the shape cut over the data.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `rejectedByIter` | `list` [`list` [`tuple`]] | Per iteration, the (x, y, score) of rejected sources (magenta x). |
| `peaks` | array-like, optional | Kept source positions (cyan). |
| `title` | `str`, optional | Plot title. |

### fitLib.selectSourcePages

`selectSourcePages(snap, sourceView='all', sourcesPerPage=6, maxPages=None)`

Choose the sources that get per-source pages.

| Parameter | Type | Description |
| --- | --- | --- |
| `snap` | `dict` | Inner-loop snapshot (``fit``, ``prevAmp``, ``isNew``). |
| `sourceView` | `str`, optional | 'all', 'new' (added this iteration) or 'changed' (new, or amplitude changed by more than 10 %). |
| `sourcesPerPage` | `int`, optional | Sources per page. |
| `maxPages` | `int`, optional | Maximum number of pages (None = no limit). |

| Returns | Type | Description |
| --- | --- | --- |
| `pages` | `list` [`list` [`int`]] | Source indices per page, brightest first. |

### fitLib.sourceColors

`sourceColors(pages)`

Assign each source on the pages a colour from ``SRC_PALETTE``.

| Parameter | Type | Description |
| --- | --- | --- |
| `pages` | `list` [`list` [`int`]] | Output of `selectSourcePages`. |

| Returns | Type | Description |
| --- | --- | --- |
| `colors` | `dict` [`int`, `str`] | Colour per source index. |

### fitLib.jointPsfFitDiag

`jointPsfFitDiag(exposure, psf, peakPositions, badBitmask, ampZeroTol=1e-06, supportFrac=0.01, varModel=None, sysFrac=0.03, poissonK=0.0)`

Fit all PSF amplitudes in one non-negative solve, keeping what is needed for diagnostics.

Weighted, non-negative, one sparse solve against the original exposure.  PSF stamps are not clipped at zero (clipping would add a positive pedestal).

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Image to fit. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `peakPositions` | `list` [`tuple` [`float`, `float`]] | Source positions in parent pixel coordinates. |
| `badBitmask` | `int` | Mask bits of pixels to ignore. |
| `ampZeroTol` | `float`, optional | Amplitudes at or below this count as zero. |
| `supportFrac` | `float`, optional | A source's significant support is where its stamp exceeds this fraction of its peak. |
| `varModel` | `numpy.ndarray`, optional | Model image used for the model-dependent variance terms (see `effectiveVariance`). |
| `sysFrac` | `float`, optional | Fractional model error added to the variance. |
| `poissonK` | `float`, optional | Extra variance per unit model flux. |

| Returns | Type | Description |
| --- | --- | --- |
| `fit` | `dict` | ``amps``, ``modelFull``, ``nAtZero``; ``comps`` (per source a dict with ``stamp``, ``sig`` (significant support), ``ys``/``xs`` (slices of its fit region), ``amp``, ``bboxParent``, or None if the stamp misses the image); coupled-source groups (``groupLabels``, ``nGroups``, ``groupMap``, ``groupChi2``); ``good`` (pixels used), ``chi2map``, ``variance`` (the effective variance), ``timings`` and the solver status (``lsqStatus``, ``lsqNit``, ``lsqCost``). |

### fitLib.computePsfDerivatives

`computePsfDerivatives(psf, x, y, bboxClip, eps=0.25)`

Compute central-difference derivative stamps dP/dx and dP/dy.

``computeImage`` centres its box on the pixel nearest the requested point, so the stamps at (x +- eps, y) and (x, y +- eps) need not share a box; all five boxes (and the clip box) are intersected first.

| Parameter | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `x, y` | `float` | Position in parent pixel coordinates. |
| `bboxClip` | `lsst.geom.Box2I` | Box to clip to. |
| `eps` | `float`, optional | Finite-difference step [px]. |

| Returns | Type | Description |
| --- | --- | --- |
| `dpdx, dpdy` | `numpy.ndarray` or `None` | Derivative stamps over ``bbox``, or None if the boxes do not overlap. |
| `bbox` | `lsst.geom.Box2I` | The common box. |

### fitLib.windowedCentroid

`windowedCentroid(imageArray, xy0, x, y, halfWidth, maxFracOfWindow=0.9)`

Return the flux-weighted centroid in a box around (x, y).

Negative pixels are clipped to zero.  Falls back to (x, y) if there is too little positive flux or the centroid lands too close to the window edge (dominated by a neighbour).

| Parameter | Type | Description |
| --- | --- | --- |
| `imageArray` | `numpy.ndarray` | Image. |
| `xy0` | `tuple` [`int`, `int`] | Parent pixel coordinates of ``imageArray[0, 0]``. |
| `x, y` | `float` | Starting position in parent pixel coordinates. |
| `halfWidth` | `int` | Half-size of the box [px]. |
| `maxFracOfWindow` | `float`, optional | Reject centroids further than this fraction of ``halfWidth`` from the start. |

| Returns | Type | Description |
| --- | --- | --- |
| `cx, cy` | `float` | Centroid in parent pixel coordinates. |

### fitLib.refinePeakPositions

`refinePeakPositions(exposure, psf, peaks, amps, badBitmask, minAmpForRefine=50.0, maxShift=0.6, nSubiters=3, minGoodPixels=10, coreRadiusSigma=2.5, variance=None, freezeSep=None)`

Move each well-detected source towards the position that minimises its residual.

Linearised Gauss-Newton with all other sources held fixed: over the source's core pixels, weighted by 1/sigma, ``(data - model) ~ dA * P + a * dx * dP/dx + a * dy * dP/dy``.  The dA column absorbs amplitude errors so they do not leak into the position; only the position step is applied and the caller re-solves the amplitudes.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `peaks` | `list` [`tuple` [`float`, `float`]] | Source positions. |
| `amps` | array-like | Source amplitudes. |
| `badBitmask` | `int` | Mask bits of pixels to ignore. |
| `minAmpForRefine` | `float`, optional | Fainter sources keep their position. |
| `maxShift` | `float`, optional | Maximum shift per sub-step [px]. |
| `nSubiters` | `int`, optional | Number of Gauss-Newton sub-steps. |
| `minGoodPixels` | `int`, optional | Minimum number of good core pixels. |
| `coreRadiusSigma` | `float`, optional | Core radius in units of the PSF sigma. |
| `variance` | `numpy.ndarray`, optional | Effective variance (default: the exposure's variance plane). |
| `freezeSep` | `float`, optional | Sources with a brighter neighbour within this many px keep their position (it is degenerate with the neighbour's amplitude). |

| Returns | Type | Description |
| --- | --- | --- |
| `peaks` | `list` [`tuple` [`float`, `float`]] | Refined positions. |

### fitLib.redetectAndRefitPsfAdaptive

`redetectAndRefitPsfAdaptive(exposure, kernelSizes=(41, 27, 21, 15, 11), **kwargs)`

Run `redetectAndRefitPsf`, trying kernel sizes from largest to smallest.

MakePsfCandidatesTask drops every candidate whose stamp touches masked or edge pixels; on a small cutout a large kernel can drop all of them.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Background-subtracted exposure. |
| `kernelSizes` | `tuple` [`int`], optional | Kernel sizes to try [px]. |
| `**kwargs` | – | Keyword arguments of `redetectAndRefitPsf`. |

| Returns | Type | Description |
| --- | --- | --- |
| `psfNew` | `lsst.afw.detection.Psf` | The new PSF (the previous one if every attempt failed). |
| `sources` | `lsst.afw.table.SourceCatalog` or `None` | Measured detections. |
| `starCat` | `lsst.afw.table.SourceCatalog` or `None` | Selected PSF stars. |
| `cellSet` | `lsst.afw.math.SpatialCellSet` or `None` | PSF candidates. |
| `kernelSize` | `int` or `None` | Kernel size that worked (None if all failed). |

### fitLib.plotIterationOverview

`plotIterationOverview(exposure, snap, hist, outerIdx, normData, colors=None, viewHalf=20, circleRadius=None)`

Plot the overview of one inner iteration.

Detection input, detections, fitted sources, model, residual and the iteration histories.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `snap` | `dict` | Inner-loop snapshot of the iteration. |
| `hist` | `dict` | Histories (``residStd``, ``nPeaks``, ``nZero``). |
| `outerIdx` | `int` | Outer-pass index. |
| `normData` | `astropy.visualization.ImageNormalize` | Display normalisation of the data. |
| `colors` | `dict` [`int`, `str`], optional | Colours of the sources with per-source pages. |
| `viewHalf` | `int`, optional | Half-size of the per-source windows [px]. |
| `circleRadius` | `float`, optional | Marker radius [px]. |

### fitLib.runInnerLoopVerbose

`runInnerLoopVerbose(exposure, psf, badBitmask, outerIdx=0, nIterations=10, thresholdStart=8.0, thresholdStep=-1.5, thresholdMin=4.5, minPixels=1, minNewPeaks=1, convergenceFrac=0.02, viz=True, plotEvery=1, sourceView='all', sourcesPerPage=6, maxPagesPerIter=None, viewHalf=20, supportFrac=0.01, refineCentroids=True, centroidMinAmp=50.0, nCentroidRounds=4, centroidShiftTol=0.03, centroidShiftSchedule=None, coarseCentroid=True, coarseCentroidHwSigma=2.5, mergeSeparationSigma=1.5, sysFrac=0.03, minSnr=5.0, pruneSepSigma=1.5, rejectFwhm=0.75, brightPeakSigma=100.0, freezeFwhm=1.0, poissonK=0.0, nearBrightFwhm=2.0, nearBrightSnr=10.0, edgeMarginFwhm=0.5, useEliminateWeak=True, useShapeCut=True, shapeKwargs=None, timer=TIMER, verbose=True, initialPeaks=None)`

Run the instrumented inner loop (detect, fit, refine, prune, repeat).

An instrumented version of `iterativeJointDetection` with centroid refinement and guards against over-counting sources.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | Image to fit. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `badBitmask` | `int` | Mask bits of pixels to ignore. |
| `outerIdx` | `int`, optional | Outer-pass index (for labels). |
| `nIterations` | `int`, optional | Maximum number of iterations. |
| `thresholdStart, thresholdStep, thresholdMin` | `float`, optional | Detection threshold schedule [sigma]. |
| `minPixels` | `int`, optional | Minimum footprint size. |
| `minNewPeaks` | `int`, optional | Stop when an iteration adds fewer peaks than this. |
| `convergenceFrac` | `float`, optional | Stop when the residual standard deviation improves by less than this fraction. |
| `viz` | `bool`, optional | Draw diagnostics. |
| `plotEvery` | `int`, optional | Draw every n-th iteration (the last always). |
| `sourceView` | `str`, optional | See `selectSourcePages`. |
| `sourcesPerPage` | `int`, optional | See `selectSourcePages`. |
| `maxPagesPerIter` | `int`, optional | Maximum per-source pages per iteration. |
| `viewHalf` | `int`, optional | Half-size of the per-source windows [px]. |
| `supportFrac` | `float`, optional | See `jointPsfFitDiag`. |
| `refineCentroids` | `bool`, optional | Refine positions (`refinePeakPositions`). |
| `centroidMinAmp` | `float`, optional | Fainter sources keep their position. |
| `nCentroidRounds` | `int`, optional | Maximum centroid rounds per iteration. |
| `centroidShiftTol` | `float`, optional | Stop centroiding when the largest shift is below this [px]. |
| `centroidShiftSchedule` | `list` [`float`], optional | Maximum shift per round [px]. |
| `coarseCentroid` | `bool`, optional | Start new peaks from a windowed centroid. |
| `coarseCentroidHwSigma` | `float`, optional | Half-width of that window in PSF sigma. |
| `mergeSeparationSigma` | `float`, optional | Merge new peaks closer than this many PSF sigma. |
| `sysFrac` | `float`, optional | Fractional PSF-shape error added in quadrature to the variance. |
| `minSnr` | `float`, optional | S/N floor (cheap version in `pruneAndMerge`, honest one in `eliminateWeak`). |
| `pruneSepSigma` | `float`, optional | Pairs closer than this many PSF sigma are merged after refinement. |
| `rejectFwhm` | `float`, optional | New peaks within this many FWHM of a bright source are discarded. |
| `brightPeakSigma` | `float`, optional | A source is bright if its peak exceeds this many noise sigma. |
| `freezeFwhm` | `float`, optional | Sources with a brighter neighbour within this many FWHM keep their position. |
| `poissonK` | `float`, optional | Extra variance per unit model flux. |
| `nearBrightFwhm` | `float`, optional | Radius around bright stars with the stricter S/N cut [FWHM]. |
| `nearBrightSnr` | `float`, optional | S/N required there. |
| `edgeMarginFwhm` | `float`, optional | Drop sources closer than this many FWHM to the edge. |
| `useEliminateWeak` | `bool`, optional | Backward elimination on the honest S/N. |
| `useShapeCut` | `bool`, optional | Backward elimination on PSF-likeness (`eliminateNonpsf`). |
| `shapeKwargs` | `dict`, optional | Keyword arguments of `shapeMetrics` / `judgeShape`. |
| `timer` | `StepTimer`, optional | Timer collecting the step durations. |
| `verbose` | `bool`, optional | Print progress. |
| `initialPeaks` | `list` [`tuple` [`float`, `float`]], optional | Warm-start source positions. |

| Returns | Type | Description |
| --- | --- | --- |
| `result` | `dict` | ``modelFull``, ``residualFull``, ``peaks``, ``amps``, the histories (``residStdHistory``, ``nPeaksHistory``, ``nZeroAmpHistory``), ``firstSeen`` and ``snaps`` (one dict per iteration with everything needed for plotting). |

### fitLib.plotAmplitudeEvolution

`plotAmplitudeEvolution(result, outerIdx=0)`

Plot the amplitude of every final source against inner iteration.

| Parameter | Type | Description |
| --- | --- | --- |
| `result` | `dict` | Inner-loop result with ``snaps``. |
| `outerIdx` | `int`, optional | Outer-pass index (for the title). |

### fitLib.plotSourceEvolution

`plotSourceEvolution(exposure, result, outerIdx=0, viewHalf=20)`

Plot model and residual cross-sections of three sources at every iteration.

The brightest, median and faintest final sources: the x cross-section of the model (top) and of the residual (bottom).

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `result` | `dict` | Inner-loop result with ``snaps``. |
| `outerIdx` | `int`, optional | Outer-pass index (for the title). |
| `viewHalf` | `int`, optional | Half-length of the cross-sections [px]. |

### fitLib.plotCutoutSetup

`plotCutoutSetup(cutout, psf, badBitmask, starPositions=None)`

Plot the cutout before fitting: data, sigma map, mask planes and PSF.

| Parameter | Type | Description |
| --- | --- | --- |
| `cutout` | `lsst.afw.image.Exposure` | The cutout. |
| `psf` | `lsst.afw.detection.Psf` | PSF model. |
| `badBitmask` | `int` | Mask bits of bad pixels (shown in red). |
| `starPositions` | array-like, optional | Original PSF star positions to mark. |

### fitLib.plotOuterStepBackground

`plotOuterStepBackground(rawExp, expNext, srcMask, modelFull, bgNext, bgPrev, outerIdx)`

Plot the background refit of one outer step.

Input image, model and source mask, new and previous background, and the result.

| Parameter | Type | Description |
| --- | --- | --- |
| `rawExp` | `lsst.afw.image.Exposure` | Input to the background refit. |
| `expNext` | `lsst.afw.image.Exposure` | Background-subtracted result. |
| `srcMask` | `numpy.ndarray` [`bool`] | Source mask used. |
| `modelFull` | `numpy.ndarray` | Joint-fit model. |
| `bgNext` | `lsst.afw.image.Image` or `None` | New background. |
| `bgPrev` | `lsst.afw.image.Image` or `None` | Previous background. |
| `outerIdx` | `int` | Outer-pass index. |

| Returns | Type | Description |
| --- | --- | --- |
| `norm` | `astropy.visualization.ImageNormalize` | Display normalisation used for the images. |

### fitLib.plotPsfCandidateStamps

`plotPsfCandidateStamps(cellSet, maxStamps=24, title='')`

Plot the PSF candidate stamps, used ones first.

| Parameter | Type | Description |
| --- | --- | --- |
| `cellSet` | `lsst.afw.math.SpatialCellSet` or `None` | PSF candidates. |
| `maxStamps` | `int`, optional | Maximum number of stamps. |
| `title` | `str`, optional | Title suffix. |

### fitLib.plotOuterStepPsf

`plotOuterStepPsf(expNext, sources, starCat, forcedPositions, psfPrev, psfNew, psfUpdated, outerIdx)`

Plot the PSF refit of one outer step: star selection and old vs new PSF.

| Parameter | Type | Description |
| --- | --- | --- |
| `expNext` | `lsst.afw.image.Exposure` | Background-subtracted exposure. |
| `sources` | `lsst.afw.table.SourceCatalog` or `None` | Fresh detections. |
| `starCat` | `lsst.afw.table.SourceCatalog` or `None` | PSF candidates. |
| `forcedPositions` | array-like or `None` | Forced original star positions. |
| `psfPrev, psfNew` | `lsst.afw.detection.Psf` | Previous and new PSF. |
| `psfUpdated` | `bool` | Whether the PSF was updated. |
| `outerIdx` | `int` | Outer-pass index. |

### fitLib.highfreqNoiseCheckArrays

`highfreqNoiseCheckArrays(resid, variance, model, good, nbins=6, verbose=True)`

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

### fitLib.highfreqNoiseCheck

`highfreqNoiseCheck(outerHistory, outerIdx=-1, snapIdx=-1, **kw)`

Run `highfreqNoiseCheckArrays` on one snapshot of a run.

| Parameter | Type | Description |
| --- | --- | --- |
| `outerHistory` | `list` [`dict`] | Outer-loop history. |
| `outerIdx` | `int`, optional | Outer pass (default: last). |
| `snapIdx` | `int`, optional | Inner-loop snapshot (default: last). |
| `**kw` | – | Keyword arguments of `highfreqNoiseCheckArrays`. |

| Returns | Type | Description |
| --- | --- | --- |
| `result` | `dict` | See `highfreqNoiseCheckArrays`. |

### fitLib.estimatePsfCorrection

`estimatePsfCorrection(data, model, variance, good, xy0, positions, amps, kernel, nUse=8, step=0.8, clip=3.0, softSigma=1.0, window=(0.6, 1.0), maxRadius=None)`

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
| `nUse` | `int`, optional | Number of brightest stars used. |
| `step` | `float`, optional | Fraction of the correction applied. |
| `clip` | `float`, optional | Sigma clipping of the stack. |
| `softSigma` | `float`, optional | Soft threshold of the correction [sigma]. |
| `window` | `tuple` [`float`, `float`], optional | Inner and outer edge of the window as fractions of the kernel half-size. |
| `maxRadius` | `float`, optional | Do not change the kernel beyond this radius [px]. |

| Returns | Type | Description |
| --- | --- | --- |
| `new` | `numpy.ndarray` | Corrected, normalised kernel. |
| `info` | `dict` | ``old``, ``new``, ``delta``, ``rawDelta``, ``err``, ``nStars``, ``used``. |

### fitLib.kernelToLsstPsf

`kernelToLsstPsf(kernel)`

Wrap a kernel array as a spatially constant LSST PSF.

| Parameter | Type | Description |
| --- | --- | --- |
| `kernel` | `numpy.ndarray` | Square kernel of odd size, centred at n//2, summing to 1. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.meas.algorithms.KernelPsf` | The PSF. |

### fitLib.extendKernelWithMoffat

`extendKernelWithMoffat(kernel, newSize=81, fitR=(6.0, 14.0), blend=(12.0, 16.0))`

Replace the noisy outer part of a PSF kernel by an analytic Moffat tail.

The azimuthal-mean profile between ``fitR[0]`` and ``fitR[1]`` px is fitted with ``A * (1 + (r/a)^2)^-b`` (relative weighting, so faint outer bins count).  The new kernel equals the original inside ``blend[0]``, the Moffat beyond ``blend[1]``, with a linear cross-fade between, and is renormalised -- fitted amplitudes then become total fluxes including the tail.  Works on an already-extended kernel too.

| Parameter | Type | Description |
| --- | --- | --- |
| `kernel` | `numpy.ndarray` | PSF kernel. |
| `newSize` | `int`, optional | Size of the new kernel [px] (odd). |
| `fitR` | `tuple` [`float`, `float`], optional | Radial range of the Moffat fit [px]. |
| `blend` | `tuple` [`float`, `float`], optional | Radial range of the cross-fade [px]. |

| Returns | Type | Description |
| --- | --- | --- |
| `new_kernel` | `numpy.ndarray` | Extended, normalised kernel. |
| `info` | `dict` | Moffat ``A``, ``a``, ``beta`` and ``tailFrac`` (flux fraction beyond ``blend[1]``). |

Raises: `RuntimeError` (Raised if the Moffat fit fails.)

### fitLib.extendPsfWings

`extendPsfWings(psf, point, size=81, verbose=True, **kw)`

Return a constant KernelPsf with a Moffat tail, taken from any PSF at a point.

| Parameter | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | Input PSF. |
| `point` | `lsst.geom.Point2D` | Position at which the kernel is taken. |
| `size` | `int`, optional | Size of the extended kernel [px]. |
| `verbose` | `bool`, optional | Print the fitted tail. |
| `**kw` | – | Keyword arguments of `extendKernelWithMoffat`. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | The extended PSF, or the input PSF if the fit fails. |

### fitLib.refinePsfFromFit

`refinePsfFromFit(outerHistory, outerIdx=-1, snapIdx=-1, **kw)`

Make one PSF update step from a finished run.

| Parameter | Type | Description |
| --- | --- | --- |
| `outerHistory` | `list` [`dict`] | Outer-loop history. |
| `outerIdx` | `int`, optional | Outer pass (default: last). |
| `snapIdx` | `int`, optional | Inner-loop snapshot (default: last). |
| `**kw` | – | Keyword arguments of `estimatePsfCorrection`. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.meas.algorithms.KernelPsf` | Updated PSF. |
| `info` | `dict` | See `estimatePsfCorrection`. |

### fitLib.polishPsfFixedSources

`polishPsfFixedSources(exposure, psf, result, badBitmask, nCycles=3, stackKwargs=None, supportFrac=0.01, sysFrac=0.03, poissonK=0.0, refineCentroids=True, centroidMinAmp=50.0, freezeSep=None, maxShift=0.3, verbose=True)`

Iterate PSF correction and re-solve without re-detecting sources.

A cycle is kept only if the residual rms over good pixels drops.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `psf` | `lsst.afw.detection.Psf` | Current PSF. |
| `result` | `dict` | Inner-loop result. |
| `badBitmask` | `int` | Mask bits of pixels to ignore. |
| `nCycles` | `int`, optional | Maximum number of cycles. |
| `stackKwargs` | `dict`, optional | Keyword arguments of `estimatePsfCorrection`. |
| `supportFrac, sysFrac, poissonK` | `float`, optional | See `jointPsfFitDiag`. |
| `refineCentroids` | `bool`, optional | Refine positions after each PSF change. |
| `centroidMinAmp` | `float`, optional | Fainter sources keep their position. |
| `freezeSep` | `float`, optional | See `refinePeakPositions`. |
| `maxShift` | `float`, optional | Maximum position shift [px]. |
| `verbose` | `bool`, optional | Print each cycle. |

| Returns | Type | Description |
| --- | --- | --- |
| `psf` | `lsst.afw.detection.Psf` | The best PSF. |
| `info` | `dict` or `None` | Information of the last kept correction (see `estimatePsfCorrection`). |

### fitLib.plotPsfCorrection

`plotPsfCorrection(info)`

Plot the old kernel, the proposed correction and the radial profiles.

| Parameter | Type | Description |
| --- | --- | --- |
| `info` | `dict` | Output of `estimatePsfCorrection`. |

### fitLib.runCutoutPipeline

`runCutoutPipeline(exposureCrop, psf, badBitmask, starCat=None, bgImage=None, xRange=(0, 100), yRange=(0, 100), nOuterIterations=2, nInnerIterations=10, innerKwargs=None, outerConvergenceFrac=0.02, refitPsf=True, psfKernelSize=41, psfKernelSizes=None, bgBinSize=32, bgUseApprox=True, viz=True, plotEvery=1, sourceView='all', sourcesPerPage=6, maxPagesPerIter=None, viewHalf=20, supportFrac=0.01, refineCentroids=True, centroidMinAmp=50.0, nCentroidRounds=4, centroidShiftTol=0.03, centroidShiftSchedule=None, coarseCentroid=True, coarseCentroidHwSigma=2.5, mergeSeparationSigma=1.5, sysFrac=0.03, minSnr=5.0, pruneSepSigma=1.5, rejectFwhm=0.75, brightPeakSigma=100.0, freezeFwhm=1.0, poissonK=0.0, nearBrightFwhm=2.0, nearBrightSnr=10.0, edgeMarginFwhm=0.5, useEliminateWeak=True, useShapeCut=True, shapeKwargs=None, shapeInFirstPass=False, stackKwargs=None, psfPolishCycles=3, psfWingSize=81)`

Run the outer loop (inner fit, background and PSF refinement) on a cutout.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposureCrop` | `lsst.afw.image.Exposure` | Region to process (before background subtraction). |
| `psf` | `lsst.afw.detection.Psf` | Initial PSF. |
| `badBitmask` | `int` | Mask bits of pixels to ignore. |
| `starCat` | `lsst.afw.table.SourceCatalog`, optional | Original PSF stars, forced into every PSF determination. |
| `bgImage` | `lsst.afw.image.Image`, optional | Detector background subtracted before the first pass. |
| `xRange, yRange` | `tuple` [`int`, `int`], optional | Sub-region of ``exposureCrop`` (see `makeCutout`). |
| `nOuterIterations` | `int`, optional | Number of outer passes. |
| `nInnerIterations` | `int`, optional | Maximum inner iterations per pass. |
| `innerKwargs` | `dict`, optional | Keyword arguments of `runInnerLoopVerbose`. |
| `outerConvergenceFrac` | `float`, optional | Stop when the residual standard deviation improves by less than this fraction. |
| `refitPsf` | `bool` or `str`, optional | False, True (PCA star-selector refit) or 'stack' (residual-stack refit from the bright stars, see `refinePsfFromFit`). |
| `psfKernelSize` | `int`, optional | PSF kernel size [px]. |
| `psfKernelSizes` | `tuple` [`int`], optional | Kernel sizes tried by `redetectAndRefitPsfAdaptive`. |
| `bgBinSize` | `int`, optional | Background bin size [px]. |
| `bgUseApprox` | `bool`, optional | See `refitBackground`. |
| `viz` | `bool`, optional | Draw diagnostics. |
| `plotEvery` | `int`, optional | Draw every n-th inner iteration (the last always). |
| `sourceView` | `str`, optional | 'all', 'new' or 'changed': which sources get a per-source page. |
| `sourcesPerPage` | `int`, optional | Sources per page. |
| `maxPagesPerIter` | `int`, optional | Maximum per-source pages per iteration (None = no limit). |
| `viewHalf` | `int`, optional | Half-size of the per-source windows [px]. |
| `supportFrac` | `float`, optional | See `jointPsfFitDiag`. |
| `refineCentroids` | `bool`, optional | Refine positions (`refinePeakPositions`). |
| `centroidMinAmp` | `float`, optional | Fainter sources keep their position. |
| `nCentroidRounds, centroidShiftTol, centroidShiftSchedule` | optional | Centroid refinement controls (see `runInnerLoopVerbose`). |
| `coarseCentroid, coarseCentroidHwSigma` | optional | Windowed starting centroid of new peaks (see `runInnerLoopVerbose`). |
| `mergeSeparationSigma` | `float`, optional | Merge new peaks closer than this many PSF sigma. |
| `sysFrac, minSnr, pruneSepSigma, rejectFwhm, brightPeakSigma` | `float`, optional | Over-counting guards (see `runInnerLoopVerbose`). |
| `freezeFwhm, poissonK, nearBrightFwhm, nearBrightSnr` | `float`, optional | Over-counting guards (see `runInnerLoopVerbose`). |
| `edgeMarginFwhm` | `float`, optional | Drop sources within this many FWHM of the cutout edge. |
| `useEliminateWeak` | `bool`, optional | Backward elimination on the honest S/N. |
| `useShapeCut` | `bool`, optional | PSF-shape constraint (`eliminateNonpsf`). |
| `shapeKwargs` | `dict`, optional | Keyword arguments of `shapeMetrics` / `judgeShape`. |
| `shapeInFirstPass` | `bool`, optional | Apply the shape cut in the first outer pass too (by default only once the PSF has been refined: a rough first PSF can make the cut reject real faint companions). |
| `stackKwargs` | `dict`, optional | Keyword arguments of `refinePsfFromFit` (``refitPsf='stack'``). |
| `psfPolishCycles` | `int`, optional | PSF correction / re-solve cycles per outer pass (``refitPsf='stack'``). |
| `psfWingSize` | `int`, optional | Convert every PSF to a constant kernel of this size with a Moffat tail (see `extendPsfWings`); None = off.  Stops the model ending at the stamp edge, which leaves a ring of fake sources half a kernel width from bright stars. |

| Returns | Type | Description |
| --- | --- | --- |
| `outerHistory` | `list` [`dict`] | One dict per outer pass: ``psf``, ``exposure``, ``result``, ``residStd``, ``nPeaks``, ``bgImage``. |

### fitLib.plotWingDiagnostic

`plotWingDiagnostic(outerHistory, nStars=3, rmax=40)`

Plot the azimuthal profiles of the brightest sources: data vs their own model.

Data beyond the dashed line well above the model means a real halo the kernel cannot represent; a difference centred on zero without a step means only core / PSF-shape mismatch remains.

| Parameter | Type | Description |
| --- | --- | --- |
| `outerHistory` | `list` [`dict`] | Outer-loop history. |
| `nStars` | `int`, optional | Number of brightest sources shown. |
| `rmax` | `float`, optional | Maximum radius [px]. |

### fitLib.plotSourceCutouts

`plotSourceCutouts(outerHistory, outerIdx=-1, snapIdx=-1, sourceIds=None, topN=None, minAmp=0.0, half=15, crossAxis='x', showOwn=True, sourcesPerFig=4, pdfPath=None, show=True, cmap=None, markNeighbours=True)`

Plot cutouts and cross-sections of individual sources of a run.

| Parameter | Type | Description |
| --- | --- | --- |
| `outerHistory` | `list` [`dict`] | Outer-loop history. |
| `outerIdx` | `int`, optional | Outer pass (default: last). |
| `snapIdx` | `int`, optional | Inner-loop snapshot (default: last). |
| `sourceIds` | `list` [`int`], optional | Sources to show; overrides ``topN`` / ``minAmp``. |
| `topN` | `int`, optional | Show the brightest N sources. |
| `minAmp` | `float`, optional | Skip sources with amplitude at or below this. |
| `half` | `int`, optional | Cutout half-width [px]. |
| `crossAxis` | `str`, optional | 'x', 'y' or 'both'. |
| `showOwn` | `bool`, optional | Draw the source's own PSF component in the cross-section. |
| `sourcesPerFig` | `int`, optional | Sources per figure. |
| `pdfPath` | `str`, optional | Also save every figure as a page of this PDF. |
| `show` | `bool`, optional | Display the figures. |
| `cmap` | `str`, optional | Colour map of the images. |
| `markNeighbours` | `bool`, optional | Mark the neighbouring sources. |

| Returns | Type | Description |
| --- | --- | --- |
| `figs` | `list` [`matplotlib.figure.Figure`] | The figures. |

### fitLib.negativeSourceCheck

`negativeSourceCheck(outerHistory, outerIdx=-1, snapIdx=-1, thresholds=(3, 4, 5), plot=True)`

Estimate the number of false detections from the residual itself.

Counts pixels of the matched-filter S/N map of the final residual at >= +N sigma and at <= -N sigma.  Real sources are positive only, so the negative excursions approximate the positive noise excursions (false positives); PSF-mismatch dipoles inflate both.

| Parameter | Type | Description |
| --- | --- | --- |
| `outerHistory` | `list` [`dict`] | Outer-loop history. |
| `outerIdx` | `int`, optional | Outer pass (default: last). |
| `snapIdx` | `int`, optional | Inner-loop snapshot (default: last). |
| `thresholds` | `tuple` [`float`], optional | S/N thresholds. |
| `plot` | `bool`, optional | Plot the S/N map and histogram. |

| Returns | Type | Description |
| --- | --- | --- |
| `snrMap` | `numpy.ndarray` | Matched-filter S/N map of the residual. |

### fitLib.sourceSnrHistogram

`sourceSnrHistogram(outerHistory, outerIdx=-1, snapIdx=-1, minSnr=5.0)`

Plot the histogram of the honest S/N of every final source.

A real population shows a gap or knee; a smooth pile-up just above ``minSnr`` means noise.

| Parameter | Type | Description |
| --- | --- | --- |
| `outerHistory` | `list` [`dict`] | Outer-loop history. |
| `outerIdx` | `int`, optional | Outer pass (default: last). |
| `snapIdx` | `int`, optional | Inner-loop snapshot (default: last). |
| `minSnr` | `float`, optional | S/N threshold marked in the plot. |

| Returns | Type | Description |
| --- | --- | --- |
| `snr` | `numpy.ndarray` | S/N per source. |

### fitLib.drawEdgeMargin

`drawEdgeMargin(ax, bbox, marginPx, color='#ff00ff', shade=True)`

Outline and shade the band within ``marginPx`` of the image border.

This is the region removed by ``dropOutOfBounds(..., margin=-marginPx)``.

| Parameter | Type | Description |
| --- | --- | --- |
| `ax` | `matplotlib.axes.Axes` | Axes to draw on. |
| `bbox` | `lsst.geom.Box2I` | Image box. |
| `marginPx` | `float` | Band width [px]. |
| `color` | `str`, optional | Line colour. |
| `shade` | `bool`, optional | Shade the band. |

### fitLib.plotFinalComparison

`plotFinalComparison(outerHistory, exposureCrop, badBitmask=None, edgeMarginFwhm=0.5, view=None, showPeaks=True, circleRadius=None, figsize=(26, 7))`

Plot original | final model | residual | residual/sigma over the fitted region.

| Parameter | Type | Description |
| --- | --- | --- |
| `outerHistory` | `list` [`dict`] | Outer-loop history. |
| `exposureCrop` | `lsst.afw.image.Exposure` | The input crop (original image panel). |
| `badBitmask` | `int`, optional | Mask bits of bad pixels. |
| `edgeMarginFwhm` | `float`, optional | Edge margin used in the run; draws its boundary (None / 0 = do not draw). |
| `view` | `tuple` [`float`, `float`, `float`, `float`], optional | (xmin, xmax, ymin, ymax) zoom. |
| `showPeaks` | `bool`, optional | Mark the sources. |
| `circleRadius` | `float`, optional | Marker radius [px]. |
| `figsize` | `tuple` [`float`, `float`], optional | Figure size [inches]. |

### fitLib.plotSourcePages

`plotSourcePages(exposure, snap, outerIdx, pages=None, colors=None, sourceView='all', sourcesPerPage=6, maxPages=None, viewHalf=20, circleRadius=None, pdf=None, show=True, snr=None, dpi=80)`

Plot one row per source: data | model | residual | x-section | y-section.

| Parameter | Type | Description |
| --- | --- | --- |
| `exposure` | `lsst.afw.image.Exposure` | The fitted exposure. |
| `snap` | `dict` | Inner-loop snapshot. |
| `outerIdx` | `int` | Outer-pass index. |
| `pages` | `list` [`list` [`int`]], optional | Source indices per page (default: `selectSourcePages`). |
| `colors` | `dict` [`int`, `str`], optional | Colour per source. |
| `sourceView` | `str`, optional | See `selectSourcePages`. |
| `sourcesPerPage` | `int`, optional | See `selectSourcePages`. |
| `maxPages` | `int`, optional | See `selectSourcePages`. |
| `viewHalf` | `int`, optional | Half-size of each window [px]. |
| `circleRadius` | `float`, optional | Marker radius [px]. |
| `pdf` | `matplotlib.backends.backend_pdf.PdfPages`, optional | Open PDF; every page is saved into it. |
| `show` | `bool`, optional | Also display each page (slow for hundreds of sources). |
| `snr` | `numpy.ndarray`, optional | Honest S/N per source, shown in the model title. |
| `dpi` | `int`, optional | Resolution of the PDF pages. |

### fitLib.sourceMetrics

`sourceMetrics(outerHistory, outerIdx=-1, snapIdx=-1, coreFwhm=1.5)`

Return per-source quantities used to stratify the sample.

| Parameter | Type | Description |
| --- | --- | --- |
| `outerHistory` | `list` [`dict`] | Outer-loop history. |
| `outerIdx` | `int`, optional | Outer pass (default: last). |
| `snapIdx` | `int`, optional | Inner-loop snapshot (default: last). |
| `coreFwhm` | `float`, optional | Radius of the chi2 region [FWHM]. |

| Returns | Type | Description |
| --- | --- | --- |
| `metrics` | `dict` | ``peaks``, ``amps``, ``snr`` (honest S/N), ``nn`` (distance to the nearest other source [px]), ``edge`` (distance to the image border [px]), ``chi2`` (mean resid^2 / sigma^2 within ``coreFwhm`` FWHM), ``group`` (size of the coupled group), ``fwhm``, ``outerPos`` and the exposure and snapshot used. |

### fitLib.selectRepresentativeSources

`selectRepresentativeSources(m, nPerCategory=4, seed=0, categories=None, extraIds=None)`

Select a representative sample of sources per category.

| Parameter | Type | Description |
| --- | --- | --- |
| `m` | `dict` | Output of `sourceMetrics`. |
| `nPerCategory` | `int` or `dict` [`str`, `int`], optional | Sources drawn per category. |
| `seed` | `int`, optional | Random seed of the tercile draws. |
| `categories` | `list` [`str`], optional | Categories to keep, in display order. |
| `extraIds` | `list` [`int`], optional | Source ids appended as a 'picked' category. |

| Returns | Type | Description |
| --- | --- | --- |
| `selection` | `collections.OrderedDict` [`str`, `list` [`int`]] | Source ids per category. |

### fitLib.plotRepresentativeSources

`plotRepresentativeSources(outerHistory, outerIdx=-1, snapIdx=-1, nPerCategory=4, seed=0, categories=None, extraIds=None, sourcesPerPage=4, viewHalf=20, pdfPath=None, show=True, dpi=80)`

Select a representative sample and plot it with `plotSourcePages`.

Each category is split into pages of ``sourcesPerPage`` rows, so a page never mixes categories; the printed table maps every source id to its category.

| Parameter | Type | Description |
| --- | --- | --- |
| `outerHistory` | `list` [`dict`] | Outer-loop history. |
| `outerIdx` | `int`, optional | Outer pass (default: last). |
| `snapIdx` | `int`, optional | Inner-loop snapshot (default: last). |
| `nPerCategory` | `int` or `dict`, optional | See `selectRepresentativeSources`. |
| `seed` | `int`, optional | See `selectRepresentativeSources`. |
| `categories` | `list` [`str`], optional | See `selectRepresentativeSources`. |
| `extraIds` | `list` [`int`], optional | See `selectRepresentativeSources`. |
| `sourcesPerPage` | `int`, optional | Rows per page. |
| `viewHalf` | `int`, optional | Half-size of each window [px]. |
| `pdfPath` | `str`, optional | Also write all pages to this PDF. |
| `show` | `bool`, optional | Display the pages. |
| `dpi` | `int`, optional | Resolution of the PDF pages. |

| Returns | Type | Description |
| --- | --- | --- |
| `selection` | `collections.OrderedDict` [`str`, `list` [`int`]] | Source ids per category. |

## adaptive

`lsst.meas.crowd.adaptive`: the residual check deciding whether sources are still missing.

Residual check that decides whether sources are still missing.

| Name | Summary |
| --- | --- |
| [`residualExcess`](#adaptiveresidualexcess) | Count positive and negative peaks of the residual's matched-filter S/N. |
| [`countsAt`](#adaptivecountsat) | Return the positive and negative peak counts at an S/N threshold. |
| [`isMissing`](#adaptiveismissing) | Decide whether the residual still contains missing sources at an S/N threshold. |
| [`pickThreshold`](#adaptivepickthreshold) | Return the lowest S/N threshold the field can support. |
| [`install`](#adaptiveinstall) | Route fitLib's inner loop through the adaptive controller. |
| [`uninstall`](#adaptiveuninstall) | Restore fitLib's original inner and outer loops. |
| [`report`](#adaptivereport) | Print what the adaptive controller did in every outer pass of a finished run. |

### adaptive.residualExcess

`residualExcess(resid, var, good, peaksXy, amps, x0, y0, fwhm, psfPeak, noise, thresholds=(3.0, 3.5, 4.0, 4.5, 5.0, 6.0), edgeFwhm=2.0, brightSigma=100.0, haloFwhm=2.5, minSepFwhm=0.75)`

Count positive and negative peaks of the residual's matched-filter S/N.

Only local maxima in the clean region count: away from the cutout edge, from bright sources (halo radius growing with flux) and from fitted sources themselves.  Real missing sources are positive only, so the positive excess over the negative count estimates how many are missing.

| Parameter | Type | Description |
| --- | --- | --- |
| `resid` | `numpy.ndarray` | Residual image. |
| `var` | `numpy.ndarray` | Variance image. |
| `good` | `numpy.ndarray` [`bool`] | Usable pixels. |
| `peaksXy` | array-like | Fitted source positions (parent pixel coordinates), shape (N, 2). |
| `amps` | array-like | Fitted source amplitudes. |
| `x0, y0` | `int` | Parent pixel coordinates of ``resid[0, 0]``. |
| `fwhm` | `float` | PSF FWHM [px]. |
| `psfPeak` | `float` | Peak value of the unit-flux PSF. |
| `noise` | `float` | Typical pixel noise. |
| `thresholds` | `tuple` [`float`], optional | S/N thresholds at which to count. |
| `edgeFwhm` | `float`, optional | Ignore peaks within this many FWHM of the edge. |
| `brightSigma` | `float`, optional | A source is bright if ``amp * psfPeak / noise`` exceeds this. |
| `haloFwhm` | `float`, optional | Ignore peaks within this many FWHM of a bright source (grows with flux). |
| `minSepFwhm` | `float`, optional | Peaks closer than this to a fitted source are its own mismatch residual. |

| Returns | Type | Description |
| --- | --- | --- |
| `result` | `dict` | ``t`` (thresholds, ascending), ``pos`` / ``neg`` (number of peaks with \|S/N\| >= t), ``scale``, ``offset`` (S/N calibration), ``areaFrac`` (fraction of the image in the clean region) and ``candidates`` ((x, y, snr, distance to the nearest fitted source) of the positive peaks). |

### adaptive.countsAt

`countsAt(d, t)`

Return the positive and negative peak counts at an S/N threshold.

| Parameter | Type | Description |
| --- | --- | --- |
| `d` | `dict` | Output of `residualExcess`. |
| `t` | `float` | S/N threshold (the nearest one in ``d['t']`` is used). |

| Returns | Type | Description |
| --- | --- | --- |
| `n_pos` | `int` | Number of positive peaks. |
| `n_neg` | `int` | Number of negative peaks. |

### adaptive.isMissing

`isMissing(d, t, zSig=2.0, minExcess=3)`

Decide whether the residual still contains missing sources at an S/N threshold.

| Parameter | Type | Description |
| --- | --- | --- |
| `d` | `dict` | Output of `residualExcess`. |
| `t` | `float` | S/N threshold. |
| `zSig` | `float`, optional | The excess must exceed ``zSig * sqrt(N_pos + N_neg)`` ... |
| `minExcess` | `int`, optional | ... and be at least this many sources. |

| Returns | Type | Description |
| --- | --- | --- |
| `missing` | `bool` | Whether sources are missing. |
| `excess` | `int` | Positive minus negative peak count. |
| `error` | `float` | Its uncertainty, sqrt(N_pos + N_neg). |

### adaptive.pickThreshold

`pickThreshold(d, fdr, tFloor, tCap)`

Return the lowest S/N threshold the field can support.

The lowest S/N in [``tFloor``, ``tCap``] whose estimated false fraction (N_neg + 1) / N_pos is at most ``fdr`` (with N_pos >= 3).

| Parameter | Type | Description |
| --- | --- | --- |
| `d` | `dict` | Output of `residualExcess`. |
| `fdr` | `float` | Tolerated false fraction. |
| `tFloor` | `float` | Lowest threshold allowed. |
| `tCap` | `float` | Highest threshold; returned if no lower one qualifies. |

| Returns | Type | Description |
| --- | --- | --- |
| `threshold` | `float` | The chosen S/N threshold. |

### adaptive.install

`install(**overrides)`

Route fitLib's inner loop through the adaptive controller.

Also enables the adaptive outer loop (through fitLib's ``_wantsMoreOuter`` hook).

| Parameter | Type | Description |
| --- | --- | --- |
| `**overrides` | – | Settings that replace entries of ``ADAPT``. |

Raises: `ValueError` (Raised if an override is not an ``ADAPT`` option.); `RuntimeError` (Raised if fitLib's inner loop has no ``initialPeaks`` argument.)

### adaptive.uninstall

`uninstall()`

Restore fitLib's original inner and outer loops.

### adaptive.report

`report(runOrHistory)`

Print what the adaptive controller did in every outer pass of a finished run.

| Parameter | Type | Description |
| --- | --- | --- |
| `runOrHistory` | `CutoutRun` or `list` [`dict`] | A finished run or its outer-loop history. |
