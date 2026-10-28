# This file is part of meas_crowd.
#
# Developed for the LSST Data Management System.
# See the COPYRIGHT file at the top-level directory of this distribution.

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
           "fisher_flux_errors", "psf_flux_errors", "zero_point_from_run", "measure_photometry"]

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

AB_NJY_ZERO_POINT = 31.4      # m_AB = -2.5 log10(flux / nJy) + 31.4


def flux_to_mag(flux, flux_err=None, zero_point=AB_NJY_ZERO_POINT):
    """Magnitude (and error) of flux in units where ``zero_point`` applies.
    Non-positive fluxes give NaN magnitudes.  Returns mag, or (mag, mag_err)."""
    flux = np.asarray(flux, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        mag = np.where(flux > 0, zero_point - 2.5 * np.log10(np.where(flux > 0, flux, 1.0)), np.nan)
        if flux_err is None:
            return mag
        err = np.where(flux > 0, 2.5 / np.log(10) * np.asarray(flux_err, float) / flux, np.nan)
    return mag, err


def _pixel_fractions(dx, dy, radius, sub=5):
    """Fraction of each pixel (centres at offsets dx, dy from the aperture centre) inside a circle
    of the given radius, by sub x sub sub-sampling (exact for pixels fully in or out)."""
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
    """Flux and error in a circular aperture of ``radius`` px around (x, y) (parent coordinates;
    the arrays start at (x0, y0)).  Bad pixels are replaced by nothing: the returned
    ``coverage`` (fraction of the aperture area on good pixels) tells how complete the sum is.
    Returns (flux, flux_err, coverage)."""
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
    """Fraction of a (centred, odd-sized) PSF kernel's flux inside ``radius`` px."""
    k = np.asarray(kernel, float)
    h = k.shape[0] // 2, k.shape[1] // 2
    yy, xx = np.mgrid[-h[0]:h[0] + 1, -h[1]:h[1] + 1]
    f = _pixel_fractions(xx.astype(float), yy.astype(float), radius)
    return float(np.sum(f * k) / k.sum())


def fisher_flux_errors(comps, variance, good, n=None):
    """1-sigma errors of the jointly fitted amplitudes: sqrt(diag((A^T W A)^-1)), inverted per
    group of overlapping sources.  ``comps`` is the joint fit's list of per-source dicts
    (stamp, ys, xs) or None.  Returns an array (inf where a source has no usable pixels)."""
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
    """Fisher errors of the PSF fluxes of a run's final joint fit (fit variance, i.e. including
    the sys_frac term used for the weights)."""
    res = run.outer_history[-1]["result"]
    fit = res["snaps"][-1]["fit"]
    return fisher_flux_errors(fit["comps"], fit["variance"], fit["good"], len(fit["amps"]))


def zero_point_from_run(run):
    """(zero_point, how): magnitude zero point of the run's image units -- from the exposure's
    PhotoCalib if it has one, else the nJy AB zero point."""
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


def measure_photometry(run, aperture_fwhm=1.5, zero_point=None, verbose=True):
    """Catalogue of the run's sources with PSF and aperture magnitudes.

    Parameters
    ----------
    run : `lsst.meas.crowd.session.CutoutRun`
        A finished run (from ``run_session_gaia``).
    aperture_fwhm : `float`
        Aperture radius in units of the PSF FWHM (aperture fluxes are corrected to total).
    zero_point : `float`, optional
        Override the zero point (default: PhotoCalib, else 31.4 for nJy).

    Returns
    -------
    table : `astropy.table.Table`
        id, x, y, ra, dec, gaia (bool), gaia_index, gaia_G, extended, in_interior,
        psf_flux, psf_flux_err, psf_snr, psf_mag, psf_mag_err,
        ap_flux, ap_flux_err, ap_coverage, ap_mag, ap_mag_err, ap_minus_psf.
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
    t["gaia"] = np.asarray(S["seed"], bool)
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
    return t
