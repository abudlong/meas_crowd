# meas_crowd: output catalogue

What each column of the source catalogue means and its units.  The pipeline itself is described
in [pipeline_guide.md](pipeline_guide.md).

- [Which file is which](#which-file-is-which)
- [Conventions](#conventions)
- [Photometry catalogue](#photometry-catalogue)
- [Gaia comparison columns](#gaia-comparison-columns)
- [Whole-detector columns](#whole-detector-columns)
- [Header keywords (table metadata)](#header-keywords-table-metadata)
- [The run catalogue (`catalog.fits` from `save_run`)](#the-run-catalogue-catalogfits-from-save_run)
- [Selecting a clean sample](#selecting-a-clean-sample)

## Which file is which

| Catalogue | Made by | Where | Columns |
| --- | --- | --- | --- |
| Photometry catalogue (the main output) | `crowd.measure_photometry(run, gaia=sess.gaia)` | in memory (an `astropy.table.Table`); `phot.write("x.fits")` to save | [Photometry](#photometry-catalogue) + [Gaia](#gaia-comparison-columns) |
| Tile catalogue | `fit_tiles` | `CATALOG` extension of each `tile_<name>.fits` | Photometry + Gaia + [whole-detector](#whole-detector-columns) |
| Detector catalogue | `stitch_tiles` | `detector_catalog.fits` | the tile catalogues merged, each source once (`owned` rows only) |
| Run catalogue | `run_session_gaia(..., out_dir=...)` / `save_run` | `<out_dir>/visit<V>_det<D>/<cutout>/catalog.fits` and `.csv` | [a smaller, older set](#the-run-catalogue-catalogfits-from-save_run) |

One row is one fitted source.  `id` is the source's index in the final fit: the same number as
`#j` on the source pages (`show_linked_source_pages`) and in `compare_passes`.

## Conventions

- **Pixel positions** (`x`, `y`, `gaia_x`, `gaia_y`) are detector coordinates in the LSST
  convention: the same frame as the `visit_image` (not relative to the cutout), with integer
  values at pixel **centres** (the pixel `(1000, 2000)` spans 999.5–1000.5, 1999.5–2000.5).
- **Sky positions** (`ra`, `dec`, `gaia_ra`, `gaia_dec`) are ICRS degrees.  The fitted ones are
  converted from `x`, `y` with the visit's WCS; Gaia's are at the **visit epoch** (proper motions
  applied).
- **Fluxes** are in the **image's units**.  LSST `visit_image`s are calibrated to **nJy**, so
  `psf_flux` etc. are in nJy for normal runs.  `meta["zero_point"]` (header `ZP`) converts them:
  `mag = ZP − 2.5 log10(flux)`; ZP = 31.4 for nJy.
- **Magnitudes** are **AB**.  Gaia G is not: see [`gaia_G`](#gaia-comparison-columns).
- Missing values are `NaN` (floats), `-1` (integer ids / indices) or `""` (strings).
- Differences are always **fitted − Gaia**.

## Photometry catalogue

Columns written by `measure_photometry`, in order.  `ra`, `dec` are present only when the
exposure has a WCS (always for a `visit_image`).

| Column | Type | Unit | Description |
| --- | --- | --- | --- |
| `id` | int | – | Source index in the final fit (0 … N−1) |
| `x`, `y` | float | px | Fitted centre, detector coordinates (see [Conventions](#conventions)) |
| `ra`, `dec` | float | deg | Sky position of the fitted centre (ICRS, from the visit WCS) |
| `gaia` | bool | – | The source was **seeded** from a Gaia star (and kept that link through the fit) |
| `gaia_index` | int | – | Row of that Gaia star in `session.gaia` (−1 = not seeded) |
| `gaia_G` | float | mag | Gaia G of the matched star (seeded sources; with `gaia=`, also position matches) |
| `extended` | bool | – | Fitted as an extended source (PSF broadened by a Gaussian, or the elliptical galaxy model) |
| `in_interior` | bool | – | Centre lies inside the requested box (False = in the `pad` margin; fitted only so edge stars have their neighbours) |
| `psf_flux` | float | image units (nJy) | Total flux from the joint PSF fit (the PSF is normalised to unit sum, so the amplitude is the total flux).  For bright stars it includes the flux of their per-star correction terms |
| `psf_flux_err` | float | image units (nJy) | 1σ error of `psf_flux` from the joint fit's Fisher matrix, inverted per group of overlapping sources, so blended sources get larger, honest errors |
| `psf_snr` | float | – | `psf_flux / psf_flux_err` |
| `psf_mag` | float | mag (AB) | `ZP − 2.5 log10(psf_flux)`; NaN for flux ≤ 0 |
| `psf_mag_err` | float | mag | `1.0857 × psf_flux_err / psf_flux` |
| `ap_flux` | float | image units (nJy) | Aperture flux on the image with all **other** sources' models subtracted, in a circle of radius `APRAD` (1.5 FWHM), corrected to total with the PSF's curve of growth (`APCORR`) |
| `ap_flux_err` | float | image units (nJy) | 1σ error of `ap_flux` from the variance plane (aperture-corrected) |
| `ap_coverage` | float | 0–1 | Fraction of the aperture area on good (unmasked, finite) pixels; < 1 near masks and edges, where `ap_flux` is less reliable |
| `ap_mag`, `ap_mag_err` | float | mag (AB), mag | Magnitude of `ap_flux` and its error |
| `ap_minus_psf` | float | mag | `ap_mag − psf_mag`; ≈ 0 for well-fitted point sources, < 0 for extended sources or a PSF that is too narrow |

## Gaia comparison columns

Added by `measure_photometry(run, gaia=sess.gaia)` or `crowd.compare_with_gaia(phot, gaia)`.
A source is matched to the Gaia star it was seeded from; a source that was not seeded is matched
to the nearest Gaia star within `gaia_match_px` (default 1 px) that no other source has claimed.

| Column | Type | Unit | Description |
| --- | --- | --- | --- |
| `gaia_match` | bool | – | The source has a Gaia counterpart |
| `gaia_match_type` | str | – | `"seed"` (seeded from it), `"position"` (nearest Gaia star within the radius) or `""` |
| `gaia_source_id` | int | – | Gaia DR3 `source_id` (−1 = no match) |
| `gaia_ra`, `gaia_dec` | float | deg | Gaia position at the visit epoch (ICRS) |
| `gaia_x`, `gaia_y` | float | px | Gaia position on the detector, through the visit WCS |
| `gaia_dx`, `gaia_dy` | float | px | `x − gaia_x`, `y − gaia_y` |
| `gaia_dra` | float | mas | `(ra − gaia_ra) × cos(gaia_dec)`: offset along RA, on the sky |
| `gaia_ddec` | float | mas | `dec − gaia_dec` |
| `gaia_sep` | float | mas | Total separation, `sqrt(gaia_dra² + gaia_ddec²)` |
| `gaia_G` | float | mag | Gaia G of the counterpart (see the note on systems below) |
| `gaia_bp_rp` | float | mag | Gaia BP − RP colour (NaN if the Gaia catalogue has no BP / RP, e.g. an older cache) |
| `gaia_dmag` | float | mag | `psf_mag − gaia_G` |
| `gaia_dmag_rel` | float | mag | `gaia_dmag − GAIAOFF`: the difference with the median offset removed (≈ 0 for most stars; what remains is noise, the colour term and real errors) |
| `gaia_flux` | float | image units (nJy) | Flux predicted from Gaia: `10^(−0.4 (gaia_G + GAIAOFF − ZP))` |
| `gaia_flux_ratio` | float | – | `psf_flux / gaia_flux` (≈ 1; > 1 = brighter than Gaia predicts) |

Notes:

- **Why an offset.**  Gaia G is a very broad band (≈ 330–1050 nm); an LSST band is much
  narrower, so `psf_mag − G` differs from star to star with colour, and between filters by a
  constant.  `GAIAOFF` is the 3σ-clipped median of `gaia_dmag` over matched point sources with
  `psf_snr ≥ 50`.  The colour term is left in `gaia_dmag_rel`; fit it against `gaia_bp_rp` if
  needed (`plot_gaia_comparison` shows it).
- **G magnitude system.**  Gaia from the ESA archive (TAP) gives Gaia's own (Vega) G; Gaia from
  the Butler reference catalogue gives `−2.5 log10(flux) + 31.4` of the catalogue's G flux in
  nJy.  The header keyword `GAIASYS` says which.  The offset absorbs the difference.
- **Position offsets.**  A common offset of all stars is a WCS offset (the pipeline removes one
  from the seeds before fitting, but the catalogue compares with Gaia as placed by the visit WCS);
  a pattern across the image is a WCS distortion; the scatter is the centroid noise.
- **Gaia stars without a source.**  `GAIAUNM` counts Gaia stars inside the fitted area that no
  source matched (too faint, masked, merged into a neighbour, or missing from the fit).

## Whole-detector columns

Added by `fit_tiles` (tile catalogues) and kept in `detector_catalog.fits`.

| Column | Type | Unit | Description |
| --- | --- | --- | --- |
| `owned` | bool | – | The centre lies in this tile's interior (half-open: a source exactly on a boundary belongs to one tile).  `stitch_tiles` keeps only owned rows, so every source appears once |
| `tile` | str | – | Tile name, `x<x0>_y<y0>_<w>x<h>` (interior origin and size, px) |

Each tile is fitted on its own, so `GAIAOFF` and `gaia_dmag_rel` are per tile; for detector-wide
work recompute the offset from `gaia_dmag` over the whole catalogue.

## Header keywords (table metadata)

In memory these are `phot.meta[...]`; in the tile and detector FITS files the long names are
shortened to 8-character keywords.

| In memory | FITS keyword | Unit | Meaning |
| --- | --- | --- | --- |
| `zero_point` | `ZP` | mag | Zero point: `mag = ZP − 2.5 log10(flux)` |
| `zero_point_source` | `ZPSRC` | – | Where it came from: `PhotoCalib`, `given`, or nJy assumed |
| `aperture_radius_px` | `APRAD` | px | Aperture radius of `ap_flux` |
| `aperture_correction` | `APCORR` | – | Fraction of a point source's flux inside the aperture (`ap_flux` is divided by it) |
| `fwhm_px` | `FWHM` | px | PSF FWHM at the cutout centre |
| `GAIAOFF` | `GAIAOFF` | mag | Median `psf_mag − G` removed in `gaia_dmag_rel` (absent from FITS headers when no star qualified) |
| `GAIANOFF` | `GAIANOFF` | – | Number of stars that set `GAIAOFF` |
| `GAIAUNM` | `GAIAUNM` | – | Gaia stars in the fitted area without a matched source |
| `GAIASYS` | `GAIASYS` | – | System of `gaia_G` (Vega from the Gaia archive, or AB from the reference catalogue) |

The tile files' primary header also has `X0`, `Y0` (fit-box origin), `IX0`, `IY0`, `IX1`, `IY1`
(interior, inclusive), `TILE`, `VISIT`, `DETECTOR` and `SECONDS` (fit time).

## The run catalogue (`catalog.fits` from `save_run`)

`save_run` (called by `run_session_gaia(..., out_dir=...)`) writes `run.catalog`, a smaller
table made during the fit.  Prefer the photometry catalogue: it has errors, aperture photometry
and the Gaia comparison.

| Column | Unit | Description |
| --- | --- | --- |
| `id`, `x`, `y`, `ra`, `dec`, `in_interior`, `extended` | as above | as above |
| `x_cutout`, `y_cutout` | px | Position relative to the fitted cutout's lower-left pixel |
| `psfFlux_nJy` | nJy (image units) | Same as `psf_flux` |
| `mag_AB` | mag | `31.4 − 2.5 log10(psfFlux_nJy)` (assumes nJy) |
| `snr` | – | S/N of the fitted amplitude (from the joint fit) |
| `chi2_core` | – | Mean (residual/σ)² within 1.5 FWHM of the source |
| `nn_pix` | px | Distance to the nearest other source |
| `edge_pix` | px | Distance to the cutout edge |
| `group_size` | – | Number of sources fitted jointly with it (overlapping group) |
| `from_gaia` | – | Same as `gaia` |
| `gaia_source_id`, `gaia_G` | –, mag | The seeding Gaia star |
| `gaia_dx`, `gaia_dy` | px | Fitted position minus the Gaia **seed** position after the pipeline's global offset correction (so not the same as the comparison columns above, which use Gaia as placed by the WCS) |
| `ext_sigma_px` | px | Size of the Gaussian broadening of an extended source (0 = point source) |

## Selecting a clean sample

```python
import numpy as np
good = (phot["in_interior"] & ~phot["extended"] & (phot["psf_snr"] > 20)
        & (phot["ap_coverage"] > 0.95) & (np.abs(phot["ap_minus_psf"]) < 0.1))
stars = phot[good]
gaia_stars = stars[stars["gaia_match"]]          # with a Gaia counterpart
```
