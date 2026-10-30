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

The numerical pieces (`fluxToMag`, `apertureFlux`, `apertureCorrection`,
`fisherFluxErrors`) work on plain arrays and are unit tested without the Rubin stack.
"""

__all__ = ["AB_NJY_ZERO_POINT", "fluxToMag", "apertureFlux", "apertureCorrection",
           "fisherFluxErrors", "psfFluxErrors", "zeroPointFromRun", "compareWithGaia",
           "measurePhotometry", "fitsMeta", "FITS_META_KEYS"]

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

AB_NJY_ZERO_POINT = 31.4      # m_AB = -2.5 log10(flux / nJy) + 31.4


def fluxToMag(flux, fluxErr=None, zeroPoint=AB_NJY_ZERO_POINT):
    """Convert fluxes (and errors) to magnitudes.

    Parameters
    ----------
    flux : array-like
        Fluxes, in units where ``zeroPoint`` applies.
    fluxErr : array-like, optional
        Flux errors.
    zeroPoint : `float`, optional
        Magnitude of unit flux (default: AB for nJy).

    Returns
    -------
    mag : `numpy.ndarray`
        Magnitudes (NaN for non-positive fluxes).
    mag_err : `numpy.ndarray`
        Magnitude errors (only if ``fluxErr`` is given).
    """
    flux = np.asarray(flux, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        mag = np.where(flux > 0, zeroPoint - 2.5 * np.log10(np.where(flux > 0, flux, 1.0)), np.nan)
        if fluxErr is None:
            return mag
        err = np.where(flux > 0, 2.5 / np.log(10) * np.asarray(fluxErr, float) / flux, np.nan)
    return mag, err


def _pixelFractions(dx, dy, radius, sub=5):
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


def apertureFlux(image, variance, good, x, y, radius, x0=0, y0=0):
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
    fluxErr : `float`
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
    f = _pixelFractions(xx + x0 - x, yy + y0 - y, radius)
    g = good[ya:yb, xa:xb] & np.isfinite(image[ya:yb, xa:xb]) & (variance[ya:yb, xa:xb] > 0)
    area = f.sum()
    if area <= 0:
        return np.nan, np.nan, 0.0
    fl = float(np.sum((f * image[ya:yb, xa:xb])[g]))
    err = float(np.sqrt(np.sum((f ** 2 * variance[ya:yb, xa:xb])[g])))
    return fl, err, float(f[g].sum() / area)


def apertureCorrection(kernel, radius):
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
    f = _pixelFractions(xx.astype(float), yy.astype(float), radius)
    return float(np.sum(f * k) / k.sum())


def fisherFluxErrors(comps, variance, good, n=None):
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


def psfFluxErrors(run):
    """Return the Fisher errors of the PSF fluxes of a run's final joint fit.

    Parameters
    ----------
    run : `lsst.meas.crowd.session.CutoutRun`
        A finished run.

    Returns
    -------
    err : `numpy.ndarray`
        Error per source (fit variance, including the sysFrac term).
    """
    res = run.outerHistory[-1]["result"]
    fit = res["snaps"][-1]["fit"]
    return fisherFluxErrors(fit["comps"], fit["variance"], fit["good"], len(fit["amps"]))


def zeroPointFromRun(run):
    """Return the magnitude zero point of a run's image units.

    Parameters
    ----------
    run : `lsst.meas.crowd.session.CutoutRun`
        A finished run.

    Returns
    -------
    zeroPoint : `float`
        Zero point (from the exposure's PhotoCalib, else 31.4 for nJy).
    how : `str`
        Where it came from.
    """
    exp = run.outerHistory[-1]["exposure"]
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


def compareWithGaia(phot, gaia, matchRadiusPx=1.0, offsetMinSnr=50.0, verbose=True):
    """Add Gaia comparison columns to a photometry table.

    Every fitted source is paired with a Gaia star: a source started from a Gaia position with
    that star (``gaiaIndex``), any other source with the nearest unclaimed Gaia star within
    ``matchRadiusPx``.  Positions are compared with Gaia's at the visit epoch (as placed on the
    detector by the visit WCS); magnitudes with Gaia G.  Gaia G is a much broader band than any
    LSST filter, so ``gaiaDmag`` contains a zero-point offset and a colour term: the median
    offset of bright, point-like matches is removed in ``gaiaDmagRel``, and ``gaiaBpRp`` is
    given to fit the colour term.

    Parameters
    ----------
    phot : `astropy.table.Table`
        Output of `measurePhotometry` (``x``, ``y``, ``psfFlux``, ``psfMag``, ``gaiaIndex``;
        ``ra``, ``dec`` if the exposure had a WCS).
    gaia : `astropy.table.Table`
        The Gaia catalogue the run started from (``session.gaia``: ``source_id``, ``ra``,
        ``dec``, ``x``, ``y``, ``phot_g_mean_mag``, optionally ``phot_bp_mean_mag`` and
        ``phot_rp_mean_mag``).
    matchRadiusPx : `float`, optional
        Match radius for sources that did not start from a Gaia position [px].
    offsetMinSnr : `float`, optional
        Minimum PSF S/N of the matches that set the median magnitude offset.
    verbose : `bool`, optional
        Print a summary.

    Returns
    -------
    table : `astropy.table.Table`
        A copy of ``phot`` with ``gaiaMatch``, ``gaiaMatchType``, ``gaiaSourceId``,
        ``gaiaRa``, ``gaiaDec``, ``gaiaX``, ``gaiaY``, ``gaiaDx``, ``gaiaDy``,
        ``gaiaDra``, ``gaiaDdec``, ``gaiaSep``, ``gaiaG``, ``gaiaBpRp``, ``gaiaDmag``,
        ``gaiaDmagRel``, ``gaiaFlux`` and ``gaiaFluxRatio`` (see doc/catalog_columns.md);
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
    gi = np.asarray(t["gaiaIndex"], np.int64) if "gaiaIndex" in t.colnames else np.full(n, -1)
    linked = (gi >= 0) & (gi < len(gaia))
    row[linked], kind[linked] = gi[linked], "linked"
    claimed = set(row[linked].tolist())
    ok = np.isfinite(gx) & np.isfinite(gy)
    if ok.any() and (~linked).any():
        idxOk = np.flatnonzero(ok)
        tree = cKDTree(np.column_stack([gx[ok], gy[ok]]))
        rest = np.flatnonzero(~linked)
        d, k = tree.query(np.column_stack([x[rest], y[rest]]))
        for j in np.argsort(d):                        # closest pairs first
            r = int(idxOk[k[j]])
            if d[j] <= matchRadiusPx and r not in claimed:
                row[rest[j]], kind[rest[j]] = r, "position"
                claimed.add(r)
    m = row >= 0

    def col(name, dtype=float, fill=np.nan):
        out = np.full(n, fill, dtype)
        if name in gaia.colnames:
            out[m] = np.asarray(gaia[name], dtype)[row[m]]
        return out

    t["gaiaMatch"] = m
    t["gaiaMatchType"] = kind
    t["gaiaSourceId"] = col("source_id", np.int64, -1)
    t["gaiaRa"], t["gaiaDec"] = col("ra"), col("dec")
    t["gaiaX"], t["gaiaY"] = col("x"), col("y")
    t["gaiaDx"], t["gaiaDy"] = x - t["gaiaX"], y - t["gaiaY"]
    if "ra" in t.colnames and "dec" in t.colnames:
        cosd = np.cos(np.radians(np.asarray(t["gaiaDec"], float)))
        dra = (np.asarray(t["ra"], float) - t["gaiaRa"] + 180.0) % 360.0 - 180.0
        t["gaiaDra"] = dra * cosd * 3.6e6
        t["gaiaDdec"] = (np.asarray(t["dec"], float) - t["gaiaDec"]) * 3.6e6
        t["gaiaSep"] = np.hypot(t["gaiaDra"], t["gaiaDdec"])
    t["gaiaG"] = col("phot_g_mean_mag")
    t["gaiaBpRp"] = col("phot_bp_mean_mag") - col("phot_rp_mean_mag")
    t["gaiaDmag"] = np.asarray(t["psfMag"], float) - t["gaiaG"]
    use = (m & np.isfinite(t["gaiaDmag"]) & (np.asarray(t["psfSnr"], float) >= offsetMinSnr))
    if "extended" in t.colnames:
        use &= ~np.asarray(t["extended"], bool)
    off, nOff = np.nan, int(use.sum())
    if nOff:
        v = np.asarray(t["gaiaDmag"], float)[use]
        for _ in range(5):                             # 3-sigma clipped median
            med = np.median(v)
            s = 1.4826 * np.median(np.abs(v - med)) or 1e-3
            keep = np.abs(v - med) < 3 * s
            if keep.all() or keep.sum() < 3:
                break
            v = v[keep]
        off, nOff = float(np.median(v)), len(v)
    t["gaiaDmagRel"] = t["gaiaDmag"] - off
    zp = float(t.meta.get("zeroPoint", AB_NJY_ZERO_POINT))
    t["gaiaFlux"] = 10 ** (-0.4 * (t["gaiaG"] + off - zp))
    with np.errstate(divide="ignore", invalid="ignore"):
        t["gaiaFluxRatio"] = np.asarray(t["psfFlux"], float) / t["gaiaFlux"]
    inside = ok & (gx >= x.min() - 0.5) & (gx <= x.max() + 0.5) & (gy >= y.min() - 0.5) \
        & (gy <= y.max() + 0.5) if n else ok & False
    nUnm = int(np.sum(inside & ~np.isin(np.arange(len(gaia)), row[m])))
    t.meta.update(GAIAOFF=off, GAIANOFF=nOff, GAIAUNM=nUnm,
                  GAIASYS=str(gaia.meta.get("gaiaMagSystem", "unknown")))
    if verbose:
        rel = np.asarray(t["gaiaDmagRel"], float)[use]
        sc = 1.4826 * np.nanmedian(np.abs(rel)) if nOff else np.nan
        sep = (np.nanmedian(np.asarray(t["gaiaSep"])[m]) if "gaiaSep" in t.colnames and m.any()
               else np.nan)
        print(f"Gaia comparison: {int(m.sum())} of {n} sources matched "
              f"({int(np.sum(kind == 'linked'))} linked, {int(np.sum(kind == 'position'))} by position); "
              f"{nUnm} Gaia stars in the area without a fitted source; psfMag - G median "
              f"{off:+.3f} from {nOff} stars (S/N >= {offsetMinSnr:g}), scatter {sc:.3f} mag; "
              f"median separation {sep:.1f} mas")
    return t


# FITS header keywords (<= 8 characters) of the photometry catalogue's metadata
FITS_META_KEYS = dict(zeroPoint="ZP", zeroPointSource="ZPSRC", apertureRadiusPx="APRAD",
                      apertureCorrection="APCORR", fwhmPx="FWHM")


def fitsMeta(meta):
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


def measurePhotometry(run, apertureFwhm=1.5, zeroPoint=None, verbose=True, gaia=None,
                      gaiaMatchPx=1.0):
    """Measure PSF and aperture magnitudes of a run's sources.

    Parameters
    ----------
    run : `lsst.meas.crowd.session.CutoutRun`
        A finished run (from ``runSessionGaia``).
    apertureFwhm : `float`, optional
        Aperture radius [FWHM]; aperture fluxes are corrected to total.
    zeroPoint : `float`, optional
        Override the zero point (default: PhotoCalib, else 31.4 for nJy).
    verbose : `bool`, optional
        Print a summary.
    gaia : `astropy.table.Table` or `bool`, optional
        The Gaia catalogue for the Gaia comparison columns of `compareWithGaia`; None (default)
        uses the run's own (``run.gaia``, the catalogue it started from) if it has one, False
        leaves the columns out.
    gaiaMatchPx : `float`, optional
        Gaia match radius for sources that did not start from a Gaia position [px].

    Returns
    -------
    table : `astropy.table.Table`
        ``id``, ``x``, ``y``, ``ra``, ``dec``, ``gaia``, ``gaiaIndex``, ``gaiaG``,
        ``extended``, ``inInterior``, ``psfFlux``, ``psfFluxErr``, ``psfSnr``, ``psfMag``,
        ``psfMagErr``, ``apFlux``, ``apFluxErr``, ``apCoverage``, ``apMag``,
        ``apMagErr``, ``apMinusPsf`` (and, with ``gaia``, the `compareWithGaia`
        columns); ``meta`` has the zero point and aperture settings.  The columns and their
        units are described in ``doc/catalog_columns.md``.
    """
    from astropy.table import Table
    from . import fitLib as L
    final = run.outerHistory[-1]
    exp, res = final["exposure"], final["result"]
    S, fit = res["S"], res["snaps"][-1]["fit"]
    bb = exp.getBBox()
    X0, Y0 = bb.getMinX(), bb.getMinY()
    data, model = exp.image.array, res["modelFull"]
    var = exp.variance.array
    good = fit["good"]
    xy = S["xy"]
    n = len(xy)
    flux = np.asarray(res["amps"], float) if len(res["amps"]) == n else np.asarray(fit["amps"], float)
    ferr = fisherFluxErrors(fit["comps"], fit["variance"], good, n)
    zp, how = (zeroPoint, "given") if zeroPoint is not None else zeroPointFromRun(run)
    psf = exp.getPsf()
    c = bb.getCenter()
    fwhm = L._fwhm(psf, c)
    radius = apertureFwhm * fwhm
    apcorr = apertureCorrection(psf.computeKernelImage(c).array, radius)
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
        if j in res.get("brightCorr", {}):
            img = img + (res["modelFull"] - fit["modelFull"]) * (own != 0)
        fl, er, cv = apertureFlux(img, var, good, xy[j, 0], xy[j, 1], radius, X0, Y0)
        ap[j], ape[j], cov[j] = fl / apcorr, er / apcorr, cv
    pm, pme = fluxToMag(flux, ferr, zp)
    am, ame = fluxToMag(ap, ape, zp)
    t = Table()
    t["id"] = np.arange(n)
    t["x"], t["y"] = xy[:, 0], xy[:, 1]
    try:
        ra, dec = exp.getWcs().pixelToSkyArray(xy[:, 0], xy[:, 1], degrees=True)
        t["ra"], t["dec"] = np.asarray(ra), np.asarray(dec)
    except Exception:
        pass
    t["gaia"] = np.asarray(S["fromGaia"], bool)
    t["gaiaIndex"] = np.asarray(S["gaia"])
    t["gaiaG"] = np.asarray(S["g"], float)
    t["extended"] = np.asarray(S["ext"]) > 0
    ib = getattr(run, "interiorBbox", None)
    if ib is not None:
        t["inInterior"] = ((xy[:, 0] >= ib.getMinX() - .5) & (xy[:, 0] <= ib.getMaxX() + .5) &
                           (xy[:, 1] >= ib.getMinY() - .5) & (xy[:, 1] <= ib.getMaxY() + .5))
    t["psfFlux"], t["psfFluxErr"] = flux, ferr
    with np.errstate(divide="ignore", invalid="ignore"):
        t["psfSnr"] = np.where(ferr > 0, flux / ferr, 0.0)
    t["psfMag"], t["psfMagErr"] = pm, pme
    t["apFlux"], t["apFluxErr"], t["apCoverage"] = ap, ape, cov
    t["apMag"], t["apMagErr"] = am, ame
    t["apMinusPsf"] = am - pm
    t.meta.update(zeroPoint=zp, zeroPointSource=how, apertureRadiusPx=radius,
                  apertureCorrection=apcorr, fwhmPx=fwhm)
    if verbose:
        ok = np.isfinite(t["apMinusPsf"]) & (t["psfSnr"] > 20) & (cov > 0.95) & ~t["extended"]
        med = np.nanmedian(t["apMinusPsf"][ok]) if ok.any() else np.nan
        print(f"photometry: {n} sources; zero point {zp:.3f} ({how}); aperture r = {radius:.1f} px "
              f"({apertureFwhm:g} FWHM, aperture correction {apcorr:.3f}); median ap - PSF "
              f"(S/N > 20, point sources) = {med:+.3f} mag")
    if gaia is None:
        gaia = getattr(run, "gaia", None)
    if gaia is not None and gaia is not False and len(gaia):
        t = compareWithGaia(t, gaia, matchRadiusPx=gaiaMatchPx, verbose=verbose)
    return t


# old snake_case names still work (with a FutureWarning)
from ._compat import moduleGetattr  # noqa: E402
__getattr__ = moduleGetattr(globals())
