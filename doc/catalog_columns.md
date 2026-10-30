# meas_crowd: output catalogue

What each column of the source catalogue means and its units.  The pipeline itself is described
in [pipeline_guide.md](pipeline_guide.md).

- [Which file is which](#which-file-is-which)
- [Conventions](#conventions)
- [Photometry catalogue](#photometry-catalogue)
    - [Measurement columns](#measurement-columns)
    - [Gaia comparison columns](#gaia-comparison-columns)
- [Whole-detector columns](#whole-detector-columns)
- [Header keywords (table metadata)](#header-keywords-table-metadata)
- [The run catalogue (`catalog.fits` from `saveRun`)](#the-run-catalogue-catalogfits-from-saverun)
- [Selecting a clean sample](#selecting-a-clean-sample)

## Which file is which

| Catalogue | Made by | Where | Columns |
| --- | --- | --- | --- |
| Photometry catalogue (the main output) | automatically at the end of every cutout fit (`runSessionGaia`, `runCutoutGaia`); again on demand with `crowd.measurePhotometry(run)` | `run.photometry` in memory; **saved** as `<outDir>/visit<V>_det<D>/<cutout>/photometry.fits` and `photometry.csv` when `runSessionGaia` has an `outDir` (or by `saveRun` / `writePhotometry`) | [measurement](#measurement-columns) + [Gaia comparison](#gaia-comparison-columns) |
| Tile catalogue | `fitTiles` | `CATALOG` extension of each `tile_<name>.fits` | the photometry catalogue + [whole-detector](#whole-detector-columns) |
| Detector catalogue | `stitchTiles` | `detector_catalog.fits` | the tile catalogues merged, each source once (`owned` rows only) |
| Run catalogue | `runSessionGaia(..., outDir=...)` / `saveRun` | `<outDir>/visit<V>_det<D>/<cutout>/catalog.fits` and `.csv` | [a smaller, older set](#the-run-catalogue-catalogfits-from-saverun) |

One row is one fitted source.  `id` is the source's index in the final fit: the same number as
`#j` on the source pages (`showLinkedSourcePages`) and in `comparePasses`.

## Conventions

- **Pixel positions** (`x`, `y`, `gaiaX`, `gaiaY`) are detector coordinates in the LSST
  convention: the same frame as the `visit_image` (not relative to the cutout), with integer
  values at pixel **centres** (the pixel `(1000, 2000)` spans 999.5–1000.5, 1999.5–2000.5).
- **Sky positions** (`ra`, `dec`, `gaiaRa`, `gaiaDec`) are ICRS degrees.  The fitted ones are
  converted from `x`, `y` with the visit's WCS; Gaia's are at the **visit epoch** (proper motions
  applied).
- **Fluxes** are in the **image's units**.  LSST `visit_image`s are calibrated to **nJy**, so
  `psfFlux` etc. are in nJy for normal runs.  `meta["zeroPoint"]` (header `ZP`) converts them:
  `mag = ZP − 2.5 log10(flux)`; ZP = 31.4 for nJy.
- **Magnitudes** are **AB**.  Gaia G is not: see [Gaia comparison columns](#gaia-comparison-columns).
- Missing values are `NaN` (floats), `-1` (integer ids / indices) or `""` (strings).
- Differences are always **fitted − Gaia**.

## Photometry catalogue

One table holds both the measurements and the comparison with Gaia.  The Gaia columns are
included whenever the run started from Gaia positions (the run keeps its Gaia table as `run.gaia`, and
`measurePhotometry` uses it by default); `measurePhotometry(run, gaia=False)` leaves them out.

### Measurement columns

Columns written by `measurePhotometry`, in order.  `ra`, `dec` are present only when the
exposure has a WCS (always for a `visit_image`).

| Column | Type | Unit | Description |
| --- | --- | --- | --- |
| `id` | int | – | Source index in the final fit (0 … N−1) |
| `x`, `y` | float | px | Fitted centre, detector coordinates (see [Conventions](#conventions)) |
| `ra`, `dec` | float | deg | Sky position of the fitted centre (ICRS, from the visit WCS) |
| `gaia` | bool | – | The source **started from a Gaia position** (and kept that link through the fit) |
| `gaiaIndex` | int | – | Row of that Gaia star in `session.gaia` (−1 = did not start from Gaia) |
| `gaiaG` | float | mag | Gaia G of the matched star (sources started from Gaia, and position matches) |
| `extended` | bool | – | Fitted as an extended source (PSF broadened by a Gaussian, or the elliptical galaxy model) |
| `inInterior` | bool | – | Centre lies inside the requested box (False = in the `pad` margin; fitted only so edge stars have their neighbours) |
| `psfFlux` | float | image units (nJy) | Total flux from the joint PSF fit (the PSF is normalised to unit sum, so the amplitude is the total flux).  For bright stars it includes the flux of their per-star correction terms |
| `psfFluxErr` | float | image units (nJy) | 1σ error of `psfFlux` from the joint fit's Fisher matrix, inverted per group of overlapping sources, so blended sources get larger, honest errors |
| `psfSnr` | float | – | `psfFlux / psfFluxErr` |
| `psfMag` | float | mag (AB) | `ZP − 2.5 log10(psfFlux)`; NaN for flux ≤ 0 |
| `psfMagErr` | float | mag | `1.0857 × psfFluxErr / psfFlux` |
| `apFlux` | float | image units (nJy) | Aperture flux on the image with all **other** sources' models subtracted, in a circle of radius `APRAD` (1.5 FWHM), corrected to total with the PSF's curve of growth (`APCORR`) |
| `apFluxErr` | float | image units (nJy) | 1σ error of `apFlux` from the variance plane (aperture-corrected) |
| `apCoverage` | float | 0–1 | Fraction of the aperture area on good (unmasked, finite) pixels; < 1 near masks and edges, where `apFlux` is less reliable |
| `apMag`, `apMagErr` | float | mag (AB), mag | Magnitude of `apFlux` and its error |
| `apMinusPsf` | float | mag | `apMag − psfMag`; ≈ 0 for well-fitted point sources, < 0 for extended sources or a PSF that is too narrow |

### Gaia comparison columns

Part of the photometry catalogue, after the measurement columns (made by `compareWithGaia`,
which `measurePhotometry` calls; `crowd.compareWithGaia(phot, gaia)` adds them to a table
made without them).
A source is matched to the Gaia star whose position it started from; any other source is matched
to the nearest Gaia star within `gaiaMatchPx` (default 1 px) that no other source has claimed.

| Column | Type | Unit | Description |
| --- | --- | --- | --- |
| `gaiaMatch` | bool | – | The source has a Gaia counterpart |
| `gaiaMatchType` | str | – | `"linked"` (started from this Gaia star's position), `"position"` (nearest Gaia star within the radius) or `""` |
| `gaiaSourceId` | int | – | Gaia DR3 `source_id` (−1 = no match) |
| `gaiaRa`, `gaiaDec` | float | deg | Gaia position at the visit epoch (ICRS) |
| `gaiaX`, `gaiaY` | float | px | Gaia position on the detector, through the visit WCS |
| `gaiaDx`, `gaiaDy` | float | px | `x − gaiaX`, `y − gaiaY` |
| `gaiaDra` | float | mas | `(ra − gaiaRa) × cos(gaiaDec)`: offset along RA, on the sky |
| `gaiaDdec` | float | mas | `dec − gaiaDec` |
| `gaiaSep` | float | mas | Total separation, `sqrt(gaiaDra² + gaiaDdec²)` |
| `gaiaG` | float | mag | Gaia G of the counterpart (see the note on systems below) |
| `gaiaBpRp` | float | mag | Gaia BP − RP colour (NaN if the Gaia catalogue has no BP / RP, e.g. an older cache) |
| `gaiaDmag` | float | mag | `psfMag − gaiaG` |
| `gaiaDmagRel` | float | mag | `gaiaDmag − GAIAOFF`: the difference with the median offset removed (≈ 0 for most stars; what remains is noise, the colour term and real errors) |
| `gaiaFlux` | float | image units (nJy) | Flux predicted from Gaia: `10^(−0.4 (gaiaG + GAIAOFF − ZP))` |
| `gaiaFluxRatio` | float | – | `psfFlux / gaiaFlux` (≈ 1; > 1 = brighter than Gaia predicts) |

Notes:

- **Why an offset.**  Gaia G is a very broad band (≈ 330–1050 nm); an LSST band is much
  narrower, so `psfMag − G` differs from star to star with colour, and between filters by a
  constant.  `GAIAOFF` is the 3σ-clipped median of `gaiaDmag` over matched point sources with
  `psfSnr ≥ 50`.  The colour term is left in `gaiaDmagRel`; fit it against `gaiaBpRp` if
  needed (`plotGaiaComparison` shows it).
- **G magnitude system.**  Gaia from the ESA archive (TAP) gives Gaia's own (Vega) G; Gaia from
  the Butler reference catalogue gives `−2.5 log10(flux) + 31.4` of the catalogue's G flux in
  nJy.  The header keyword `GAIASYS` says which.  The offset absorbs the difference.
- **Position offsets.**  A common offset of all stars is a WCS offset (the pipeline removes one
  from the Gaia starting positions before fitting, but the catalogue compares with Gaia as placed by the visit WCS);
  a pattern across the image is a WCS distortion; the scatter is the centroid noise.
- **Gaia stars without a source.**  `GAIAUNM` counts Gaia stars inside the fitted area that no
  source matched (too faint, masked, merged into a neighbour, or missing from the fit).

## Whole-detector columns

Added by `fitTiles` (tile catalogues) and kept in `detector_catalog.fits`.

| Column | Type | Unit | Description |
| --- | --- | --- | --- |
| `owned` | bool | – | The centre lies in this tile's interior (half-open: a source exactly on a boundary belongs to one tile).  `stitchTiles` keeps only owned rows, so every source appears once |
| `tile` | str | – | Tile name, `x<x0>_y<y0>_<w>x<h>` (interior origin and size, px) |

Each tile is fitted on its own, so `GAIAOFF` and `gaiaDmagRel` are per tile; for detector-wide
work recompute the offset from `gaiaDmag` over the whole catalogue.

## Header keywords (table metadata)

In memory these are `phot.meta[...]`; in the tile and detector FITS files the long names are
shortened to 8-character keywords.

| In memory | FITS keyword | Unit | Meaning |
| --- | --- | --- | --- |
| `zeroPoint` | `ZP` | mag | Zero point: `mag = ZP − 2.5 log10(flux)` |
| `zeroPointSource` | `ZPSRC` | – | Where it came from: `PhotoCalib`, `given`, or nJy assumed |
| `apertureRadiusPx` | `APRAD` | px | Aperture radius of `apFlux` |
| `apertureCorrection` | `APCORR` | – | Fraction of a point source's flux inside the aperture (`apFlux` is divided by it) |
| `fwhmPx` | `FWHM` | px | PSF FWHM at the cutout centre |
| `GAIAOFF` | `GAIAOFF` | mag | Median `psfMag − G` removed in `gaiaDmagRel` (absent from FITS headers when no star qualified) |
| `GAIANOFF` | `GAIANOFF` | – | Number of stars that set `GAIAOFF` |
| `GAIAUNM` | `GAIAUNM` | – | Gaia stars in the fitted area without a matched source |
| `GAIASYS` | `GAIASYS` | – | System of `gaiaG` (Vega from the Gaia archive, or AB from the reference catalogue) |

The tile files' primary header also has `X0`, `Y0` (fit-box origin), `IX0`, `IY0`, `IX1`, `IY1`
(interior, inclusive), `TILE`, `VISIT`, `DETECTOR` and `SECONDS` (fit time).

## The run catalogue (`catalog.fits` from `saveRun`)

`saveRun` (called by `runSessionGaia(..., outDir=...)`) also writes `run.catalog`, a smaller
table made during the fit, next to `photometry.fits`.  Prefer the photometry catalogue: it has
errors, aperture photometry and the Gaia comparison.

| Column | Unit | Description |
| --- | --- | --- |
| `id`, `x`, `y`, `ra`, `dec`, `inInterior`, `extended` | as above | as above |
| `xCutout`, `yCutout` | px | Position relative to the fitted cutout's lower-left pixel |
| `psfFlux_nJy` | nJy (image units) | Same as `psfFlux` |
| `magAB` | mag | `31.4 − 2.5 log10(psfFlux_nJy)` (assumes nJy) |
| `snr` | – | S/N of the fitted amplitude (from the joint fit) |
| `chi2Core` | – | Mean (residual/σ)² within 1.5 FWHM of the source |
| `nnPix` | px | Distance to the nearest other source |
| `edgePix` | px | Distance to the cutout edge |
| `groupSize` | – | Number of sources fitted jointly with it (overlapping group) |
| `fromGaia` | – | Same as `gaia` |
| `gaiaSourceId`, `gaiaG` | –, mag | The Gaia star the source started from |
| `gaiaDx`, `gaiaDy` | px | Fitted position minus the Gaia **starting** position after the pipeline's global offset correction (so not the same as the comparison columns above, which use Gaia as placed by the WCS) |
| `extSigmaPx` | px | Size of the Gaussian broadening of an extended source (0 = point source) |

## Selecting a clean sample

```python
import numpy as np
good = (phot["inInterior"] & ~phot["extended"] & (phot["psfSnr"] > 20)
        & (phot["apCoverage"] > 0.95) & (np.abs(phot["apMinusPsf"]) < 0.1))
stars = phot[good]
gaia_stars = stars[stars["gaiaMatch"]]          # with a Gaia counterpart
```
