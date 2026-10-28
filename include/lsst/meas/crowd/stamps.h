// -*- LSST-C++ -*-
/*
 * This file is part of meas_crowd.
 *
 * Developed for the LSST Data Management System.
 * See the COPYRIGHT file at the top-level directory of this distribution.
 */
#ifndef LSST_MEAS_CROWD_STAMPS_H
#define LSST_MEAS_CROWD_STAMPS_H

#include <cstddef>
#include <vector>

namespace lsst {
namespace meas {
namespace crowd {

/**
 * Normalised Lanczos weights for interpolating a sampled function at position p = m0 - frac,
 * i.e. a shift of the samples by +frac pixels.
 *
 * @param[in]  frac   Sub-pixel shift (any value; |frac| <= 0.5 is typical).
 * @param[in]  order  Lanczos order n (afw uses 5 for PSF images).
 * @param[out] first  Offset of weights[0] relative to the output pixel.
 * @param[out] weights  2n weights, summing to 1.
 */
void lanczosWeights(double frac, int order, int &first, std::vector<double> &weights);

/**
 * Shift a (ny, nx) row-major image by (dx, dy) pixels with a separable Lanczos kernel;
 * pixels outside the input are treated as zero.  Same convention as
 * lsst::afw::math::offsetImage: out(x, y) = in(x - dx, y - dy).
 */
void shiftImage(double const *in, std::size_t ny, std::size_t nx, double dx, double dy, int order,
                double *out);

/**
 * Tame the wings of a PSF stamp in place, exactly like cutout_fit_lib.psf_stamp(taper=...):
 * subtract the median of the outermost pixel ring, then multiply by a raised-cosine taper that
 * is 1 inside r0 = t0 * h and 0 outside r1 = t1 * h (h = half the stamp size, r measured from
 * (cx, cy) in stamp pixel coordinates).
 */
void taperStamp(double *arr, std::size_t ny, std::size_t nx, double cx, double cy, double t0,
                double t1);

/**
 * Render PSF stamps for many sources from one (constant) kernel image.
 *
 * Source i is centred at (xs[i], ys[i]) (parent pixel coordinates).  Its stamp has the kernel's
 * size, origin (round(x) - hx, round(y) - hy), and is the kernel shifted by the fractional part
 * of the position -- what lsst::afw::detection::Psf::computeImage returns for a KernelPsf.
 *
 * @param[out] out     n * ny * nx values.
 * @param[out] x0, y0  stamp origins (n each).
 */
void renderStamps(double const *kernel, std::size_t ny, std::size_t nx, double const *xs,
                  double const *ys, std::size_t n, int order, bool taper, double t0, double t1,
                  double *out, long *x0, long *y0);

/**
 * Central-difference derivatives dP/dx, dP/dy (step eps) of the stamps of renderStamps (no
 * taper), all on the same origin as the centred stamp.
 */
void renderDerivatives(double const *kernel, std::size_t ny, std::size_t nx, double const *xs,
                       double const *ys, std::size_t n, int order, double eps, double *dx,
                       double *dy, long *x0, long *y0);

}  // namespace crowd
}  // namespace meas
}  // namespace lsst

#endif  // LSST_MEAS_CROWD_STAMPS_H
