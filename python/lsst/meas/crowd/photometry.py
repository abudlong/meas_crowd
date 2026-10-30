# This file is part of meas_crowd.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Magnitudes of the fitted sources.

For every source of a finished run:

* **PSF flux** -- the joint-fit amplitude (plus the bright-star correction terms' flux), with
  its error from the inverse Fisher matrix of the joint fit (overlapping sources are inverted
  together, so blends get honest, larger errors);
* **aperture flux** -- summed in a circle on the image with all *other* sources subtracted,
  corrected to total flux with the PSF's own curve of growth;
* AB **magnitudes** for both, through the exposure's PhotoCalib when present, else a zero point.

The numerical pieces (`flux_to_mag`, `aperture_flux`, `aperture_correction`,
`fisher_flux_errors`) work on plain arrays and are unit tested without the Rubin stack.
"""

__all__ = ["AB_NJY_ZERO_POINT", "flux_to_mag", "aperture_flux", "aperture_correction",
           "fisher_flux_errors", "psf_flux_errors", "zero_point_from_run", "compare_with_gaia",
           "measure_photometry", "fits_meta", "FITS_META_KEYS"]

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

AB_NJY_ZERO_POINT = 31.4      # m_AB = -2.5 log10(flux / nJy) + 31.4


def flux_to_mag(flux, flux_err=None, zero_point=AB_NJY_ZERO_POINT):
    """Convert fluxes (and errors) to magnitudes.

    Parameters
    ----------
    flux : array-like
        Fluxes, in units where ``zero_point`` applies.
    flux_err : array-like, optional
        Flux errors.
    zero_point : `float`, optional
        Magnitude of unit flux (default: AB for nJy).

    Returns
    -------
    mag : `numpy.ndarray`
        Magnitudes (NaN for non-positive fluxes).
    mag_err : `numpy.ndarray`
        Magnitude errors (only if ``flux_err`` is given).
    """
    flux = np.asarray(flux, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        mag = np.where(flux > 0, zero_point - 2.5 * np.log10(np.where(flux > 0, flux, 1.0)), np.nan)
        if flux_err is None:
            return mag
        err = np.where(flux > 0, 2.5 / np.log(10) * np.asarray(flux_err, float) / flux, np.nan)
    return mag, err


def _pixel_fractions(dx, dy, radius, sub=5):
    """Return the fraction of each pixel inside a circle (5x5 sub-sampling at the edge)."""
    d = np.hypot(dx, dy)
    frac = (d <= radius - 0.7072).astype(float)
    edge = (d > radius - 0.7072) & (d < radius + 0.7072)
    if edge.any():
        o = (np.arange(sub) + 0.5) / sub - 0.5
        ox, oy = np.meshgrid(o, o)
        ex, ey = dx[edge][:, None, None] + ox, dy[edge][:, None, None] + oy
        frac[edge] = np.mean(np.hypot(ex, ey) <= radius, axis=(1, 2))
    return frac


def aperture_flux(image, variance, good, x, y, radius, x0=0, y0=0):
    """Sum the flux in a circular aperture.

    Parameters
    ----------
    image : `numpy.ndarray`
        Image.
    variance : `numpy.ndarray`
        Variance image.
    good : `numpy.ndarray` [`bool`]
        Usable pixels (others are left out of the sum).
    x, y : `float`
        Centre in parent pixel coordinates.
    radius : `float`
        Aperture radius [px].
    x0, y0 : `int`, optional
        Parent pixel coordinates of ``image[0, 0]``.

    Returns
    -------
    flux : `float`
        Aperture flux (NaN if the aperture misses the image).
    flux_err : `float`
        Its error.
    coverage : `float`
        Fraction of the aperture area on good pixels.
    """
    ny, nx = image.shape
    r = int(np.ceil(radius)) + 1
    ix, iy = int(round(x)) - x0, int(round(y)) - y0
    xa, xb = max(ix - r, 0), min(ix + r + 1, nx)
    ya, yb = max(iy - r, 0), min(iy + r + 1, ny)
    if xb <= xa or yb <= ya:
        return np.nan, np.nan, 0.0
    yy, xx = np.mgrid[ya:yb, xa:xb]
    f = _pixel_fractions(xx + x0 - x, yy + y0 - y, radius)
    g = good[ya:yb, xa:xb] & np.isfinite(image[ya:yb, xa:xb]) & (variance[ya:yb, xa:xb] > 0)
    area = f.sum()
    if area <= 0:
        return np.nan, np.nan, 0.0
    fl = float(np.sum((f * image[ya:yb, xa:xb])[g]))
    err = float(np.sqrt(np.sum((f ** 2 * variance[ya:yb, xa:xb])[g])))
    return fl, err, float(f[g].sum() / area)


def aperture_correction(kernel, radius):
    """Return the fraction of a PSF kernel's flux inside a radius.

    Parameters
    ----------
    kernel : `numpy.ndarray`
        PSF kernel (odd size, centred).
    radius : `float`
        Aperture radius [px].

    Returns
    -------
    fraction : `float`
        Enclosed flux fraction.
    """
    k = np.asarray(kernel, float)
    h = k.shape[0] // 2, k.shape[1] // 2
    yy, xx = np.mgrid[-h[0]:h[0] + 1, -h[1]:h[1] + 1]
    f = _pixel_fractions(xx.astype(float), yy.astype(float), radius)
    return float(np.sum(f * k) / k.sum())


def fisher_flux_errors(comps, variance, good, n=None):
    """Return the 1-sigma errors of jointly fitted amplitudes.

    ``sqrt(diag((A^T W A)^-1))``, inverted per group of overlapping sources, so blends get
    larger, honest errors.

    Parameters
    ----------
    comps : `list` [`dict` or `None`]
        The joint fit's per-source components (``stamp``, ``ys``, ``xs``).
    variance : `numpy.ndarray`
        Variance used as weights.
    good : `numpy.ndarray` [`bool`]
        Pixels used in the fit.
    n : `int`, optional
        Number of sources (default: ``len(comps)``).

    Returns
    -------
    err : `numpy.ndarray`
        Error per source (inf where a source has no usable pixels).
    """
    n = len(comps) if n is None else n
    pix = -np.ones(good.shape, dtype=np.int64)
    idx = np.flatnonzero(good.ravel())
    pix.ravel()[idx] = np.arange(len(idx))
    R, C, V = [], [], []
    for j, c in enumerate(comps):
        if c is None:
            continue
        p = pix[c["ys"], c["xs"]]
        ok = (p >= 0) & (c["stamp"] != 0)
        w = 1.0 / np.sqrt(variance[c["ys"], c["xs"]][ok])
        R.append(p[ok])
        C.append(np.full(int(ok.sum()), j))
        V.append(c["stamp"][ok] * w)
    err = np.full(n, np.inf)
    if not R:
        return err
    A = coo_matrix((np.concatenate(V), (np.concatenate(R), np.concatenate(C))),
                   shape=(len(idx), n)).tocsr()
    F = (A.T @ A).tocsr()
    _, lab = connected_components(F, directed=False)
    for g in np.unique(lab):
        ii = np.flatnonzero(lab == g)
        block = F[ii][:, ii].toarray()
        if not np.any(np.diag(block) > 0):
            continue
        keep = np.diag(block) > 0
        sub = block[np.ix_(keep, keep)]
        try:
            cov = np.linalg.inv(sub)
        except np.linalg.LinAlgError:
            cov = np.linalg.pinv(sub)
        err[ii[keep]] = np.sqrt(np.clip(np.diag(cov), 0, None))
    return err


def psf_flux_errors(run):
    """Return the Fisher errors of the PSF fluxes of a run's final joint fit.

    Parameters
    ----------
    run : `lsst.meas.crowd.session.CutoutRun`
        A finished run.

    Returns
    -------
    err : `numpy.ndarray`
        Error per source (fit variance, including the sys_frac term).
    """
    res = run.outer_history[-1]["result"]
    fit = res["snaps"][-1]["fit"]
    return fisher_flux_errors(fit["comps"], fit["variance"], fit["good"], len(fit["amps"]))


def zero_point_from_run(run):
    """Return the magnitude zero point of a run's image units.

    Parameters
    ----------
    run : `lsst.meas.crowd.session.CutoutRun`
        A finished run.

    Returns
    -------
    zero_point : `float`
        Zero point (from the exposure's PhotoCalib, else 31.4 for nJy).
    how : `str`
        Where it came from.
    """
    exp = run.outer_history[-1]["exposure"]
    try:
        pc = exp.getPhotoCalib()
    except Exception:
        pc = None
    if pc is not None:
        try:
            s = float(pc.instFluxToNanojansky(1.0))
            if np.isfinite(s) and s > 0:
                return AB_NJY_ZERO_POINT - 2.5 * np.log10(s), "PhotoCalib"
        except Exception:
            pass
    return AB_NJY_ZERO_POINT, "image units assumed to be nJy"


def compare_with_gaia(phot, gaia, match_radius_px=1.0, offset_min_snr=50.0, verbose=True):
    """Add Gaia comparison columns to a photometry table.

    Every fitted source is paired with a Gaia star: a source started from a Gaia position with
    that star (``gaia_index``), any other source with the nearest unclaimed Gaia star within
    ``match_radius_px``.  Positions are compared with Gaia's at the visit epoch (as placed on the
    detector by the visit WCS); magnitudes with Gaia G.  Gaia G is a much broader band than any
    LSST filter, so ``gaia_dmag`` contains a zero-point offset and a colour term: the median
    offset of bright, point-like matches is removed in ``gaia_dmag_rel``, and ``gaia_bp_rp`` is
    given to fit the colour term.

    Parameters
    ----------
    phot : `astropy.table.Table`
        Output of `measure_photometry` (``x``, ``y``, ``psf_flux``, ``psf_mag``, ``gaia_index``;
        ``ra``, ``dec`` if the exposure had a WCS).
    gaia : `astropy.table.Table`
        The Gaia catalogue the run started from (``session.gaia``: ``source_id``, ``ra``,
        ``dec``, ``x``, ``y``, ``phot_g_mean_mag``, optionally ``phot_bp_mean_mag`` and
        ``phot_rp_mean_mag``).
    match_radius_px : `float`, optional
        Match radius for sources that did not start from a Gaia position [px].
    offset_min_snr : `float`, optional
        Minimum PSF S/N of the matches that set the median magnitude offset.
    verbose : `bool`, optional
        Print a summary.

    Returns
    -------
    table : `astropy.table.Table`
        A copy of ``phot`` with ``gaia_match``, ``gaia_match_type``, ``gaia_source_id``,
        ``gaia_ra``, ``gaia_dec``, ``gaia_x``, ``gaia_y``, ``gaia_dx``, ``gaia_dy``,
        ``gaia_dra``, ``gaia_ddec``, ``gaia_sep``, ``gaia_G``, ``gaia_bp_rp``, ``gaia_dmag``,
        ``gaia_dmag_rel``, ``gaia_flux`` and ``gaia_flux_ratio`` (see doc/catalog_columns.md);
        ``meta`` gets ``GAIAOFF`` (the median offset), ``GAIANOFF`` (stars used),
        ``GAIAUNM`` (Gaia stars in the fitted area without a fitted source) and ``GAIASYS``.
    """
    from astropy.table import Table
    from scipy.spatial import cKDTree
    t = Table(phot, copy=True)
    n = len(t)
    gx, gy = np.asarray(gaia["x"], float), np.asarray(gaia["y"], float)
    x, y = np.asarray(t["x"], float), np.asarray(t["y"], float)
    row = np.full(n, -1, np.int64)
    kind = np.array([""] * n, dtype="U8")
    gi = np.asarray(t["gaia_index"], np.int64) if "gaia_index" in t.colnames else np.full(n, -1)
    linked = (gi >= 0) & (gi < len(gaia))
    row[linked], kind[linked] = gi[linked], "linked"
    claimed = set(row[linked].tolist())
    ok = np.isfinite(gx) & np.isfinite(gy)
    if ok.any() and (~linked).any():
        idx_ok = np.flatnonzero(ok)
        tree = cKDTree(np.column_stack([gx[ok], gy[ok]]))
        rest = np.flatnonzero(~linked)
        d, k = tree.query(np.column_stack([x[rest], y[rest]]))
        for j in np.argsort(d):                        # closest pairs first
            r = int(idx_ok[k[j]])
            if d[j] <= match_radius_px and r not in claimed:
                row[rest[j]], kind[rest[j]] = r, "position"
                claimed.add(r)
    m = row >= 0

    def col(name, dtype=float, fill=np.nan):
        out = np.full(n, fill, dtype)
        if name in gaia.colnames:
            out[m] = np.asarray(gaia[name], dtype)[row[m]]
        return out

    t["gaia_match"] = m
    t["gaia_match_type"] = kind
    t["gaia_source_id"] = col("source_id", np.int64, -1)
    t["gaia_ra"], t["gaia_dec"] = col("ra"), col("dec")
    t["gaia_x"], t["gaia_y"] = col("x"), col("y")
    t["gaia_dx"], t["gaia_dy"] = x - t["gaia_x"], y - t["gaia_y"]
    if "ra" in t.colnames and "dec" in t.colnames:
        cosd = np.cos(np.radians(np.asarray(t["gaia_dec"], float)))
        dra = (np.asarray(t["ra"], float) - t["gaia_ra"] + 180.0) % 360.0 - 180.0
        t["gaia_dra"] = dra * cosd * 3.6e6
        t["gaia_ddec"] = (np.asarray(t["dec"], float) - t["gaia_dec"]) * 3.6e6
        t["gaia_sep"] = np.hypot(t["gaia_dra"], t["gaia_ddec"])
    t["gaia_G"] = col("phot_g_mean_mag")
    t["gaia_bp_rp"] = col("phot_bp_mean_mag") - col("phot_rp_mean_mag")
    t["gaia_dmag"] = np.asarray(t["psf_mag"], float) - t["gaia_G"]
    use = (m & np.isfinite(t["gaia_dmag"]) & (np.asarray(t["psf_snr"], float) >= offset_min_snr))
    if "extended" in t.colnames:
        use &= ~np.asarray(t["extended"], bool)
    off, n_off = np.nan, int(use.sum())
    if n_off:
        v = np.asarray(t["gaia_dmag"], float)[use]
        for _ in range(5):                             # 3-sigma clipped median
            med = np.median(v)
            s = 1.4826 * np.median(np.abs(v - med)) or 1e-3
            keep = np.abs(v - med) < 3 * s
            if keep.all() or keep.sum() < 3:
                break
            v = v[keep]
        off, n_off = float(np.median(v)), len(v)
    t["gaia_dmag_rel"] = t["gaia_dmag"] - off
    zp = float(t.meta.get("zero_point", AB_NJY_ZERO_POINT))
    t["gaia_flux"] = 10 ** (-0.4 * (t["gaia_G"] + off - zp))
    with np.errstate(divide="ignore", invalid="ignore"):
        t["gaia_flux_ratio"] = np.asarray(t["psf_flux"], float) / t["gaia_flux"]
    inside = ok & (gx >= x.min() - 0.5) & (gx <= x.max() + 0.5) & (gy >= y.min() - 0.5) \
        & (gy <= y.max() + 0.5) if n else ok & False
    n_unm = int(np.sum(inside & ~np.isin(np.arange(len(gaia)), row[m])))
    t.meta.update(GAIAOFF=off, GAIANOFF=n_off, GAIAUNM=n_unm,
                  GAIASYS=str(gaia.meta.get("gaia_mag_system", "unknown")))
    if verbose:
        rel = np.asarray(t["gaia_dmag_rel"], float)[use]
        sc = 1.4826 * np.nanmedian(np.abs(rel)) if n_off else np.nan
        sep = (np.nanmedian(np.asarray(t["gaia_sep"])[m]) if "gaia_sep" in t.colnames and m.any()
               else np.nan)
        print(f"Gaia comparison: {int(m.sum())} of {n} sources matched "
              f"({int(np.sum(kind == 'linked'))} linked, {int(np.sum(kind == 'position'))} by position); "
              f"{n_unm} Gaia stars in the area without a fitted source; psf_mag - G median "
              f"{off:+.3f} from {n_off} stars (S/N >= {offset_min_snr:g}), scatter {sc:.3f} mag; "
              f"median separation {sep:.1f} mas")
    return t


# FITS header keywords (<= 8 characters) of the photometry catalogue's metadata
FITS_META_KEYS = dict(zero_point="ZP", zero_point_source="ZPSRC", aperture_radius_px="APRAD",
                      aperture_correction="APCORR", fwhm_px="FWHM")


def fits_meta(meta):
    """Return catalogue metadata with FITS-compatible keywords and values.

    Long names are shortened (`FITS_META_KEYS`, else the first 8 characters, upper case);
    values that are not numbers, strings or booleans, and non-finite numbers, are dropped.

    Parameters
    ----------
    meta : `dict`
        Table metadata.

    Returns
    -------
    meta : `dict`
        The metadata to write.
    """
    out = {}
    for k, v in meta.items():
        if isinstance(v, (str, bool, np.bool_)) or (
                isinstance(v, (int, float, np.floating, np.integer)) and np.isfinite(v)):
            out[FITS_META_KEYS.get(k, k[:8].upper())] = v
    return out


def measure_photometry(run, aperture_fwhm=1.5, zero_point=None, verbose=True, gaia=None,
                       gaia_match_px=1.0):
    """Measure PSF and aperture magnitudes of a run's sources.

    Parameters
    ----------
    run : `lsst.meas.crowd.session.CutoutRun`
        A finished run (from ``run_session_gaia``).
    aperture_fwhm : `float`, optional
        Aperture radius [FWHM]; aperture fluxes are corrected to total.
    zero_point : `float`, optional
        Override the zero point (default: PhotoCalib, else 31.4 for nJy).
    verbose : `bool`, optional
        Print a summary.
    gaia : `astropy.table.Table` or `bool`, optional
        The Gaia catalogue for the Gaia comparison columns of `compare_with_gaia`; None (default)
        uses the run's own (``run.gaia``, the catalogue it started from) if it has one, False
        leaves the columns out.
    gaia_match_px : `float`, optional
        Gaia match radius for sources that did not start from a Gaia position [px].

    Returns
    -------
    table : `astropy.table.Table`
        ``id``, ``x``, ``y``, ``ra``, ``dec``, ``gaia``, ``gaia_index``, ``gaia_G``,
        ``extended``, ``in_interior``, ``psf_flux``, ``psf_flux_err``, ``psf_snr``, ``psf_mag``,
        ``psf_mag_err``, ``ap_flux``, ``ap_flux_err``, ``ap_coverage``, ``ap_mag``,
        ``ap_mag_err``, ``ap_minus_psf`` (and, with ``gaia``, the `compare_with_gaia`
        columns); ``meta`` has the zero point and aperture settings.  The columns and their
        units are described in ``doc/catalog_columns.md``.
    """
    from astropy.table import Table
    from . import fitLib as L
    final = run.outer_history[-1]
    exp, res = final["exposure"], final["result"]
    S, fit = res["S"], res["snaps"][-1]["fit"]
    bb = exp.getBBox()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    data, model = exp.image.array, res["model_full"]
    var = exp.variance.array
    good = fit["good"]
    xy = S["xy"]
    n = len(xy)
    flux = np.asarray(res["amps"], float) if len(res["amps"]) == n else np.asarray(fit["amps"], float)
    ferr = fisher_flux_errors(fit["comps"], fit["variance"], good, n)
    zp, how = (zero_point, "given") if zero_point is not None else zero_point_from_run(run)
    psf = exp.getPsf()
    c = bb.getCenter()
    fwhm = L._fwhm(psf, c)
    radius = aperture_fwhm * fwhm
    apcorr = aperture_correction(psf.computeKernelImage(c).array, radius)
    resid = data - model
    ap = np.full(n, np.nan)
    ape = np.full(n, np.nan)
    cov = np.zeros(n)
    for j in range(n):
        cj = fit["comps"][j]
        own = np.zeros_like(data)
        if cj is not None:
            own[cj["ys"], cj["xs"]] = fit["amps"][j] * cj["stamp"]
        # this source's light alone: residual + its own model (neighbours stay subtracted).
        # A bright star's correction terms are part of its light too; they are not stored per
        # star, so the correction light within its stamp is used (exact unless two corrected
        # bright stars overlap)
        img = resid + own
        if j in res.get("bright_corr", {}):
            img = img + (res["model_full"] - fit["model_full"]) * (own != 0)
        fl, er, cv = aperture_flux(img, var, good, xy[j, 0], xy[j, 1], radius, X0, Y0)
        ap[j], ape[j], cov[j] = fl / apcorr, er / apcorr, cv
    pm, pme = flux_to_mag(flux, ferr, zp)
    am, ame = flux_to_mag(ap, ape, zp)
    t = Table()
    t["id"] = np.arange(n)
    t["x"], t["y"] = xy[:, 0], xy[:, 1]
    try:
        ra, dec = exp.getWcs().pixelToSkyArray(xy[:, 0], xy[:, 1], degrees=True)
        t["ra"], t["dec"] = np.asarray(ra), np.asarray(dec)
    except Exception:
        pass
    t["gaia"] = np.asarray(S["from_gaia"], bool)
    t["gaia_index"] = np.asarray(S["gaia"])
    t["gaia_G"] = np.asarray(S["g"], float)
    t["extended"] = np.asarray(S["ext"]) > 0
    ib = getattr(run, "interior_bbox", None)
    if ib is not None:
        t["in_interior"] = ((xy[:, 0] >= ib.getMinX() - .5) & (xy[:, 0] <= ib.getMaxX() + .5) &
                            (xy[:, 1] >= ib.getMinY() - .5) & (xy[:, 1] <= ib.getMaxY() + .5))
    t["psf_flux"], t["psf_flux_err"] = flux, ferr
    with np.errstate(divide="ignore", invalid="ignore"):
        t["psf_snr"] = np.where(ferr > 0, flux / ferr, 0.0)
    t["psf_mag"], t["psf_mag_err"] = pm, pme
    t["ap_flux"], t["ap_flux_err"], t["ap_coverage"] = ap, ape, cov
    t["ap_mag"], t["ap_mag_err"] = am, ame
    t["ap_minus_psf"] = am - pm
    t.meta.update(zero_point=zp, zero_point_source=how, aperture_radius_px=radius,
                  aperture_correction=apcorr, fwhm_px=fwhm)
    if verbose:
        ok = np.isfinite(t["ap_minus_psf"]) & (t["psf_snr"] > 20) & (cov > 0.95) & ~t["extended"]
        med = np.nanmedian(t["ap_minus_psf"][ok]) if ok.any() else np.nan
        print(f"photometry: {n} sources; zero point {zp:.3f} ({how}); aperture r = {radius:.1f} px "
              f"({aperture_fwhm:g} FWHM, aperture correction {apcorr:.3f}); median ap - PSF "
              f"(S/N > 20, point sources) = {med:+.3f} mag")
    if gaia is None:
        gaia = getattr(run, "gaia", None)
    if gaia is not None and gaia is not False and len(gaia):
        t = compare_with_gaia(t, gaia, match_radius_px=gaia_match_px, verbose=verbose)
    return t
