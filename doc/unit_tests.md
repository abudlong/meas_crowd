# meas_crowd: unit tests

Every unit test in `tests/`: what it checks, what it feeds in, and what it must get back to
pass.  How to run them and read the results is in the
[pipeline guide](pipeline_guide.md#unit-tests-running-them-and-checking-the-results).

The suite has **97 tests in 9 files**.  Without the Rubin stack 95 run and the 2 in
`test_fastStampsStack.py` are skipped; `testCxxEqualsNumpy` is also skipped if the C++ extension
was not built.  `test_docstrings.py` checks that this page lists every test and that the
[function reference](function_reference.md) is up to date.

- [How the tests work](#how-the-tests-work)
- [test_fastKernels.py](#test_fastkernelspy): fast PSF stamps (C++ and NumPy)
- [test_fastStampsStack.py](#test_faststampsstackpy): fast stamps against the stack
- [test_photometry.py](#test_photometrypy): magnitudes, apertures, errors, Gaia comparison
- [test_pipelineFunctions.py](#test_pipelinefunctionspy): the pipeline's building blocks
- [test_diagnostics.py](#test_diagnosticspy): residual diagnostics and run-level plots
- [test_sourcePlots.py](#test_sourceplotspy): per-source cutouts
- [test_stagePlots.py](#test_stageplotspy): stage plots and pass-to-pass comparison
- [test_detector.py](#test_detectorpy): whole-detector tiling, stitching and plots
- [test_docstrings.py](#test_docstringspy): documentation checks

## How the tests work

Each test builds a small **synthetic** image or table whose correct answer is known (stars of
known flux and position, a known PSF, known noise), runs one function on it and compares the
result with that answer.  No test needs data access.

Shared helpers:

| Helper | What it provides |
| --- | --- |
| `tests/lsstStubs.py` | Stand-ins for the Rubin stack so the tests run without it: working `lsst.geom` boxes and points, `FakeExposure` (image, variance, mask, PSF, optional PhotoCalib), `FakePsf` (a fixed kernel image), `gaussian_kernel`, `gaussian_comp` (a unit Gaussian stamp at a sub-pixel position, in the joint fit's component format) and `place` (add a stamp into an image).  `install()` does nothing when the real stack is set up |
| `tests/fakeRun.py` | `make_run(...)`: a finished run with the same structure as `run_session_gaia`'s result: a grid of Gaussian stars (default 4 × 4 in a 120 × 130 px image at detector position (1000, 2000), fluxes 3×10³–3×10⁵, sky 100) with Poisson noise, fitted by least squares with the true PSF (a perfect model) unless `model_sigma` or `sigma_fn` make the model wrong on purpose.  `run.truth_flux` holds the true fluxes |

Plotting tests use the non-interactive `Agg` backend; they check that the figure or file is made
and that returned numbers are right, not how the picture looks.

## test_fastKernels.py

The fast PSF-stamp kernels in `fastKernels.py` (Lanczos-5 sub-pixel shifts, stamps, derivatives,
taper), in both the C++ and the NumPy implementation.

| Test | What it checks | Input | Passes when |
| --- | --- | --- | --- |
| `LanczosWeightsTestCase.testNormalised` | Lanczos-5 interpolation weights conserve flux | sub-pixel shifts −0.5, −0.31, 0, 0.2, 0.49 | 10 weights each, summing to 1 (to 10⁻¹²) |
| `LanczosWeightsTestCase.testZeroShiftIsIdentity` | A zero shift leaves the image unchanged | shift 0 | the weight on the centre pixel is 1, all others 0 |
| `ShiftTestCase.testShiftMatchesAnalyticGaussian` | `shift_image` moves an image by a fraction of a pixel correctly | a Gaussian (σ = 2.6 px), shifts (0.3, −0.2), (−0.45, 0.1), (0, 0.5); NumPy and C++ | result matches the exactly shifted Gaussian to 0.3 % of the peak |
| `ShiftTestCase.testIntegerShift` | A whole-pixel shift is an exact copy | Gaussian, shift (1, 0) | columns shifted by one, to 10⁻¹² |
| `ShiftTestCase.testCxxEqualsNumpy` | The C++ and NumPy versions give the same numbers | an elliptical kernel, 50 random positions, with and without taper; stamps and derivatives | identical to 10⁻¹² of the peak (skipped if C++ is not built) |
| `RenderTestCase.testOriginConvention` | Stamps sit where `lsst.afw` puts them | positions x = 50.4, 50.6, 50.5 | stamp origin = floor(x + 0.5) − half, as afw's `computeImage` |
| `RenderTestCase.testFluxConserved` | A rendered stamp keeps unit flux | Gaussian kernel at (10.37, 12.81) | stamp sums to 1 (to 10⁻⁵) |
| `RenderTestCase.testDerivativesMatchAnalytic` | Position derivatives used by the centroid fit | Gaussian (σ = 2.2 px) at (30.2, 40.35), step 0.25 px | match the exact Gaussian's central differences to 1.5 % of the peak |
| `RenderTestCase.testTaperMatchesReference` | The wing taper equals the original library's (pedestal removed, raised-cosine edge) | elliptical kernel + pedestal, taper (0.75, 1.0) | equal to the formula to 10⁻¹⁵ (NumPy and C++) |
| `FastPsfTestCase.testStampClip` | Stamps are clipped to the image box | `FastPsf` stamp near the corner, and one far outside | box (0, 0, 15, 15), shape 16 × 16; outside → None |
| `FastPsfTestCase.testDerivativesClip` | Derivative stamps are clipped too | derivatives at (50.1, 50.2) with box starting at (45, 45) | box (45, 45, 60, 60), both derivative arrays the same shape |
| `FastPsfTestCase.testValidate` | The check against the reference PSF accepts a match and rejects a mismatch | a reference renderer with the same kernel, and one 10 % too bright | first validates True, second False |

## test_fastStampsStack.py

Runs only with the Rubin stack (skipped otherwise).

| Test | What it checks | Input | Passes when |
| --- | --- | --- | --- |
| `StackComparisonTestCase.testStampsMatchComputeImage` | Fast stamps reproduce afw's `psf.computeImage` | a `KernelPsf` (elliptical core + faint exponential wing), 20 random positions | same origin; pixels agree to 0.2 % of the peak |
| `StackComparisonTestCase.testDerivativesMatchFitLib` | Fast derivatives reproduce the original `fitLib.compute_psf_derivatives` | the same PSF at 3 positions | agree to 0.5 % of the largest derivative |

## test_photometry.py

Magnitudes, aperture photometry, flux errors, the photometry catalogue and the Gaia comparison
(`photometry.py`).

| Test | What it checks | Input | Passes when |
| --- | --- | --- | --- |
| `MagnitudeTestCase.testFluxToMag` | Flux → AB magnitude with the nJy zero point | 3631 Jy, 1 nJy, 0, −5 | 0.00, 31.4, NaN, NaN |
| `MagnitudeTestCase.testMagError` | Magnitude error | flux 100 ± 1 | 1.0857 × 0.01 mag |
| `MagnitudeTestCase.testZeroPoint` | A given zero point is used | flux 10, ZP 25 | 22.5 |
| `MagnitudeTestCase.testZeroPointFromRun` | Zero point from the exposure | exposure without PhotoCalib; with 1 count = 10 nJy | 31.4 ("nJy assumed"); 28.9 ("PhotoCalib") |
| `ApertureTestCase.testPixelFractionsArea` | Partial-pixel aperture weights | circles r = 2, 3.7, 6.25 px off pixel centre | weights sum to πr² (within 2 %) |
| `ApertureTestCase.testApertureFluxOfGaussian` | Aperture flux, error and coverage | Gaussian of flux 1000 (σ = 2 px), unit variance; r = 15 px and r = 1σ | 1000 ± 0.5 with error √(πr²) and coverage 1; 39.3 % of the flux inside 1σ |
| `ApertureTestCase.testCoverage` | Fraction of the aperture on good pixels | left half of the image masked, aperture on the boundary | coverage 0.5 |
| `ApertureTestCase.testOffImage` | An aperture off the image | centre at (100, 100) on a 10 × 10 image | flux NaN, coverage 0 |
| `ApertureTestCase.testApertureCorrection` | PSF curve of growth | Gaussian kernel σ = 2; r = 2 and r = 30 | 0.393; 1.0 |
| `FisherErrorTestCase.testIsolatedSource` | Flux error of one source | one Gaussian, variance 4 | 1 / √(Σ P² / var) |
| `FisherErrorTestCase.testBlendHasLargerErrors` | Blended sources get larger errors, separate ones do not | one source; a pair 2 px apart; a pair far apart | pair error > 1.5 × single; far pair = single |
| `FisherErrorTestCase.testMissingComponent` | A source without a stamp | components [None, Gaussian] | its error is infinite, the other finite |
| `FisherErrorTestCase.testErrorsMatchScatter` | The Fisher errors are the real scatter | two overlapping sources, 400 noise realisations fitted by least squares | standard deviation of the fitted fluxes = predicted error (within 12 %) |
| `MeasurePhotometryTestCase.testSyntheticRun` | The full photometry catalogue | `fakeRun` with 25 stars | 25 rows; (PSF flux − truth) / error has mean < 1 and std < 1.6; ap − PSF median < 0.01 mag for S/N > 50; `psf_mag` = 31.4 − 2.5 log10(flux) |
| `MeasurePhotometryTestCase.testPhotoCalib` | Magnitudes use the PhotoCalib | `fakeRun` with 1 count = 2 nJy | `psf_mag` = 31.4 − 2.5 log10(2 × flux) |
| `GaiaComparisonTestCase.testMatchingAndDifferences` | Gaia matching and every difference column | 5 sources (3 Gaia-seeded, 1 within 0.45 px of a Gaia star, 1 unmatched), 6 Gaia stars; known offsets 0.1 px, 2 mas in RA, 1 mas in Dec, +0.3 mag | match flags and types (seed / position / none), source ids, `gaia_dx` = 0.1, `gaia_dra` = 2, `gaia_ddec` = 1, `gaia_sep` = √5 mas, offset 0.3, `gaia_dmag_rel` = 0, flux ratio 1, BP − RP 0.7, 1 unmatched Gaia star; input table unchanged |
| `GaiaComparisonTestCase.testRadiusAndNoSky` | Match radius; table without sky coordinates | the same, radius 0.3 px, no `ra`/`dec` | the 0.45 px source is not matched; no mas columns |
| `GaiaComparisonTestCase.testFromSyntheticRun` | Gaia comparison through `measure_photometry` | `fakeRun` + Gaia table with G = true magnitude − 0.2 | every source matched, offset 0.2 ± 0.01 mag, flux-ratio scatter < 2 % |
| `GaiaComparisonTestCase.testRunGaiaUsedByDefault` | The run's own Gaia table is used | `fakeRun` with `run.gaia` set | Gaia columns present by default, absent with `gaia=False` |
| `GaiaComparisonTestCase.testWritePhotometry` | The photometry catalogue is saved | `fakeRun` with Gaia, `write_photometry` | `photometry.fits` and `.csv` written; all rows and Gaia columns read back; FITS header has `ZP` and `GAIAOFF`; `run.photometry` kept |
| `GaiaComparisonTestCase.testFitsMeta` | Header keywords made FITS-safe | metadata with long names, a NaN and a list | shortened to `ZP`, `APRAD`; NaN and list dropped |
| `GaiaComparisonTestCase.testPlot` | `plot_gaia_comparison` | the 5-source table | PNG written; 4 sources plotted; ValueError without Gaia columns |

## test_pipelineFunctions.py

The pipeline's building blocks (`pipeline.py`), each on a synthetic image with a known answer.
The joint fit itself is replaced by a plain least-squares fit where the test needs one.

| Test | What it checks | Input | Passes when |
| --- | --- | --- | --- |
| `PsfKernelTestCase.testSmoothFillsHoles` | PSF kernel regularisation | a Moffat kernel with 4 wing pixels at 1 % and wing noise | non-negative, sums to 1, core unchanged, holes filled (> 50 % of truth), closer to the truth in the wings, ellipticity kept |
| `PsfKernelTestCase.testCoreCorrectionFixesWrongPsf` | PSF core correction from the median star residual | 49 bright stars with a lopsided core, fitted with a plain Gaussian PSF | corrected kernel 4× closer to the truth, non-negative |
| `PsfKernelTestCase.testCoreCorrectionLeavesCorrectPsf` | ... and leaves a correct PSF alone | the same stars with the correct PSF | kernel unchanged (to 10⁻¹²) |
| `BrightStarTermsTestCase.testCoreRings` | Per-star correction terms for a sharper core | a bright star (flux 2.2×10⁶) with 3 % of its light in a core sharper than the PSF | χ² > 1.5 without the narrow core rings, < 1.1 with; central residual 10× smaller |
| `BrightStarTermsTestCase.testLopsidedCore` | Lopsided (off-centre) core terms | the same star with the sharp core 0.3 px off | χ² lower with core modes, and < 1 |
| `ModelErrorTestCase.testRecoversModelError` | Measuring the PSF model error ε | 121 stars whose truth differs from the model by 0.5 % per pixel | ε = 0.005 ± 0.0015 |
| `ModelErrorTestCase.testZeroForPerfectModel` | ... and zero when the model is right | the same with ε = 0 | ε = 0 |
| `SaturationLevelTestCase.testLevels` | Saturation level from the SAT pixels | SAT pixels at 98000; then interpolated SAT pixels and one bright pixel at 60000; then no SAT; then a given level | 98000; 60000 ("lower bound"); None; the given 1.2×10⁵ |
| `RingMimicTestCase.testSmoothRing` | Removing fake sources on a bright star's ring | a bright star with a smooth ring, 6 fake sources on the ring, 2 real faint stars nearby | exactly the 6 fakes removed |
| `RingMimicTestCase.testDiffractionLobes` | ... with 4-fold diffraction lobes | the same with lobes and 4 fakes on them | exactly the 4 fakes removed |
| `SplitAndPolishTestCase.testSplitFindsPairs` | Splitting unresolved pairs | 8 single stars and 5 pairs 2.6 px apart, each pair fitted as one source | 5 sources added, one next to each pair, none at a single star |
| `SplitAndPolishTestCase.testPolishFixesStalePositions` | Joint position polish | 49 stars, 40 % of them with positions off by 0.2 px | χ² drops, mean position error 3× smaller, χ²/pixel < 1.02 |
| `EllipticalSourceTestCase.testGalaxyReplacesPoints` | The elliptical galaxy model | an elliptical galaxy fitted as 2 points (one Gaia), plus a star | the 2 points become 1 galaxy (on the Gaia source), the star stays a point, χ² < 1.01 |
| `EllipticalSourceTestCase.testStarPlusGalaxy` | Star next to a galaxy | a Gaia star with a galaxy 5 px away | the star stays a point, its neighbour becomes the galaxy, χ² < 1.01 |
| `EllipticalSourceTestCase.testStarStaysPoint` | No false galaxies | two point sources | nothing changed |
| `FastStampSwitchTestCase.testAgreeingPsfUsesFastPath` | The fast stamps are used when they agree with the PSF | a PSF whose reference stamps match the fast ones | fast path on; stamp origin (88, 39); pixels equal to 10⁻¹⁴ |
| `FastStampSwitchTestCase.testVaryingPsfFallsBack` | ... and not when they disagree | reference stamps 5 % different | fast path off (falls back to `computeImage`) |
| `FastStampSwitchTestCase.testSwitchOff` | `fast_stamps=False` | switch off | fast path off |
| `BackgroundLowFitTestCase.testPlaneRecovered` | Plane background refit | a tilted sky plane with noise, a masked star and an unmasked 80-count blob | plane recovered to 0.15 counts; subtracted in place |
| `BackgroundLowFitTestCase.testConstantAndNone` | Constant and no refit; bad mode | the same image | constant = 4.0 ± 0.3 and flat; "none" subtracts nothing; "spline" here raises ValueError |
| `BackgroundLowFitTestCase.testDefaultIsPlane` | The default | `SEEDED_DEFAULTS` | `bg_mode` is "plane" |

## test_diagnostics.py

Residual diagnostics and run-level plots (`plotting.py`, `pipeline.unmodeled_near_mask`), on
`fakeRun` runs of 121 stars (400 × 400 px) unless noted.

| Test | What it checks | Input | Passes when |
| --- | --- | --- | --- |
| `StarResidualTableTestCase.testPerfectModel` | `star_residual_table` on a perfect model | perfect model | verdict "photon noise"; median χ² < 1.1 |
| `StarResidualTableTestCase.testBrighterFatter` | Detects a brightness-dependent PSF | true width grows 9 % per decade of flux | verdict "depends on brightness"; core-vs-flux correlation < −0.5 |
| `StarResidualTableTestCase.testPositionDependence` | Detects a position-dependent PSF | true width varies by ±6 % from left to right | verdict mentions "position" |
| `StarResidualTableTestCase.testPlot` | Its plots | perfect model, plot on | runs without error |
| `CoreExcessTestCase.testBroadPsf` | `core_excess_diagnostic` finds a too-broad PSF | model 3 % broader than the truth | verdict "too broad"; median core excess > 2 % |
| `CoreExcessTestCase.testFluxDependent` | ... and a flux-dependent core | true width shrinks with flux | verdict "grows with flux" |
| `CoreExcessTestCase.testPerfect` | ... and nothing for a perfect model | perfect model | "no significant central excess" |
| `UnmodeledMaskTestCase.testNoPadAroundModeledCores` | The "cannot be modelled" mask | a near-saturated core (modelled from its wings), a missing neighbour, a bad column | > 100 px greyed at the core without the `modeled` hint, 0 with it; the bad column keeps its margin; nothing greyed with no bad pixels |
| `PlotSmokeTestCase.testFinalComparison` | `plot_final_comparison` options | default `fakeRun`; 4 option sets | runs; residual/σ spread = 1 ± 0.1 |
| `PlotSmokeTestCase.testCrossSectionsAndZooms` | `show_cross_sections`, `show_star_zooms` | one star; the 2 worst stars | one set of cuts returned; zooms run |
| `PlotSmokeTestCase.testBrightStarStack` | `show_bright_star_stack` | bright stars (S/N ≥ 100), perfect model | > 5 stars stacked; stacked residual < 10⁻³ of the flux |

## test_sourcePlots.py

Per-source cutouts, cross-sections and their selection (`sourcePlots.py`).

| Test | What it checks | Input | Passes when |
| --- | --- | --- | --- |
| `CutoutDataTestCase.testPerfectModel` | `source_cutout_data` windows and cuts | a Gaussian star (flux 5000) at (140.3, 235.7), data = model | window (130, 150, 226, 246), 21 × 21; all residuals and χ² 0; x-cut peaks within 0.5 px of the star |
| `CutoutDataTestCase.testNoiseOnlyChi2` | χ² of the cuts with noise only | the same star + Poisson noise | χ² per cut ≈ 1 (within 0.8) |
| `CutoutDataTestCase.testBadPixelsIgnored` | Bad pixels are left out of the cuts | the star's row masked and set to 10⁹ | x-cut residual stays 0 |
| `CutoutDataTestCase.testEdge` | A source at the image corner | star at (101, 201) | window clipped at the image edge (100, 200) |
| `SelectTestCase.testSelections` | `select_sources` | 5 sources with S/N and χ² | worst = [2, 3] (S/N 5 excluded), brightest = [4, 2], faintest = [3], random 3, by position = [2, 0]; unknown mode raises ValueError |
| `PlotSmokeTestCase.testPlotAndPdf` | `plot_source_cutouts` with a PDF | default `fakeRun`, 3 brightest, 2 per page | PDF written; 3 rows with finite magnitudes |

## test_stagePlots.py

The stage plots and the pass-to-pass comparison (`stagePlots.py`).  Most use a `fakeRun` turned
into two outer passes with three inner snapshots (the first missing half the stars) and a tilted
background in pass 2.

| Test | What it checks | Input | Passes when |
| --- | --- | --- | --- |
| `PsfStarsTestCase.testCatalogFlags` | PSF stars from a catalogue's flags | table with `slot_Centroid_x/y`, `calib_psf_candidate`, `calib_psf_used` | 3 candidates; the 2 used stars at the right positions |
| `PsfStarsTestCase.testSession` | PSF stars of a refit session | session with `starCat` and `psf_used_xy` | 2 candidates, 1 used; ValueError for a visit-PSF session or no input |
| `PsfStarsTestCase.testRunStarsAreBrightSources` | PSF stars of the in-cutout stack | two-pass run | some stars used, no more than the candidates, each one of the run's sources |
| `PsfStarsTestCase.testPlotKernel` | `show_psf_stars` | two-pass run | PNG written; FWHM = 2.355 × 2 ± 0.15 px; kernel sums to 1 |
| `BackgroundAndRoundsTestCase.testBackgrounds` | `show_backgrounds` | initial background 2.0; pass-2 tilted plane; then no initial | 2 backgrounds, pass 2 equal to the input; None for pass 1 without an initial one |
| `BackgroundAndRoundsTestCase.testOuterRounds` | `show_rounds(level="outer")` | two-pass run, perfect final model | 2 rows; last residual/σ spread 1 ± 0.1 |
| `BackgroundAndRoundsTestCase.testInnerRounds` | `show_rounds(level="inner")` | 3 snapshots, zoomed view | source counts 8, 16, 16; residual/σ larger while stars are missing; `max_rows=2` gives 2 rows; bad `level` raises ValueError |
| `LinkedPagesTestCase.testColoursMatch` | Overview boxes match the page colours | sources by index and by position | pages [[0, 5]]; box colours = page colours |
| `LinkedPagesTestCase.testDefaultSelectionAndSlimSnapshot` | Default selection; old snapshots | 4 per page | the 4 brightest; ValueError for a snapshot without stamps |
| `LinkedPagesTestCase.testOutDirWritesEveryPageAndClosesFigures` | `out_dir` writes one PNG per page and frees memory | all 16 sources, 4 per page | 4 pages, 8 files (overview + page each), no figures left open |
| `ComparePassesTestCase.testVerdicts` | `compare_passes` finds and explains changes | pass a: 15 stars + 1 fake source; pass b: drops 2 real stars (one marked as removed by the halo check), gains 1 | 3 lost, 1 gained; real star "LOST" with ~100 % of its flux left; halo-check star "blocked"; fake "spurious" with < 20 % left; gained star "missing from the residual" |
| `ComparePassesTestCase.testPlotAndErrors` | Its plot; a one-pass run | the same | PNG written; ValueError with one pass |

## test_detector.py

Whole-detector tiling, per-tile files, stitching, plots and the command line (`detector.py`,
`bin.src/crowdDetector.py`).  The per-tile fit (which needs the stack) is replaced by a perfect
synthetic one: a 600 × 500 px "detector" with 60 Gaussian stars; everything around it is the
real code.

| Test | What it checks | Input | Passes when |
| --- | --- | --- | --- |
| `TileSpecTestCase.testPartition` | Tile interiors cover the detector exactly once | 4 detector / tile sizes incl. LSSTCam 4072 × 4000 with 400 px tiles | every pixel in exactly one tile; 110 tiles for LSSTCam |
| `TileSpecTestCase.testRegion` | Only tiles overlapping a region; tile names | region in a 600 × 500 detector | one tile at (200, 0); name "x0200_y0000_200x167" |
| `FitStitchTestCase.testTasksResumeAndStitch` | Slurm-style task split, restart, stitch | 2 tasks of 9 tiles, then a rerun | tasks share the 9 tiles without overlap; rerun refits nothing; mosaic data and model equal the inputs; residual/σ spread ≈ 1; every star once with the right flux; files reload |
| `FitStitchTestCase.testRegionAndGaps` | A partial run | only the corner tile | 1 file; mosaic 200 × 167 px |
| `FitStitchTestCase.testFailingTileDoesNotStopTask` | One failing tile | one tile raises an error | the other 8 written; the gap marked in the TILE map |
| `FitStitchTestCase.testPlot` | `plot_detector` | stitched mosaic; a zoom, greyscale, sources circled | PNG written; residual/σ spread ≈ 1; binning 1 for the zoom |
| `FitStitchTestCase.testTile` | `load_tile`, `plot_tile` | the corner tile | interior at (0, 0), width = interior + 20 px pad; `owned` count right; PNG written; residual/σ ≈ 1 |
| `FitStitchTestCase.testCli` | `crowdDetector.py stitch` and `plot` | fitted tiles | both return 0; `detector_view.png` written |

## test_docstrings.py

| Test | What it checks | Input | Passes when |
| --- | --- | --- | --- |
| `DocstringTestCase.testNumpydoc` | Every public function, class and method has an LSST (numpydoc) docstring | all 11 modules | numpydoc's validator reports no problems (Rubin's usual exemptions; skipped if numpydoc is not installed) |
| `FunctionReferenceTestCase.testFunctionReferenceUpToDate` | The [function reference](function_reference.md) matches the docstrings | `doc/make_function_reference.py` run in memory | the generated page equals `doc/function_reference.md` (skipped if numpydoc is not installed) |
| `UnitTestDocTestCase.testUnitTestDocListsEveryTest` | This page is complete | the test methods of every `tests/test_*.py` | every `Class.method` appears on this page |
