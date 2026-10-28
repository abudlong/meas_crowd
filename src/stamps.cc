// -*- LSST-C++ -*-
/*
 * This file is part of meas_crowd.
 *
 * Developed for the LSST Data Management System.
 * See the COPYRIGHT file at the top-level directory of this distribution.
 */
#include "lsst/meas/crowd/stamps.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lsst {
namespace meas {
namespace crowd {

namespace {

constexpr double PI = 3.14159265358979323846;

double sinc(double x) {
    if (std::abs(x) < 1e-12) return 1.0;
    double const px = PI * x;
    return std::sin(px) / px;
}

double lanczos(double x, int order) {
    if (std::abs(x) >= order) return 0.0;
    return sinc(x) * sinc(x / order);
}

// 1-d pass: out[i] = sum_j w[j] * in[(i + first + j) * stride]; zero outside [0, n)
void pass1d(double const *in, std::size_t n, std::size_t stride, int first,
            std::vector<double> const &w, double *out, std::size_t ostride) {
    long const nl = static_cast<long>(n);
    long const nw = static_cast<long>(w.size());
    double const *wp = w.data();
    // interior: every tap inside the input, no bounds checks
    long const lo = std::max(0L, -static_cast<long>(first));
    long const hi = std::min(nl, nl - first - nw + 1);
    for (long i = 0; i < nl; ++i) {
        double s = 0.0;
        if (i >= lo && i < hi) {
            double const *src = in + (i + first) * static_cast<long>(stride);
            for (long j = 0; j < nw; ++j) s += wp[j] * src[j * static_cast<long>(stride)];
        } else {
            for (long j = 0; j < nw; ++j) {
                long const m = i + first + j;
                if (m >= 0 && m < nl) s += wp[j] * in[m * static_cast<long>(stride)];
            }
        }
        out[i * ostride] = s;
    }
}

}  // namespace

void lanczosWeights(double frac, int order, int &first, std::vector<double> &weights) {
    if (order < 1) throw std::invalid_argument("Lanczos order must be >= 1");
    first = static_cast<int>(std::floor(-frac)) - order + 1;
    weights.assign(2 * order, 0.0);
    double sum = 0.0;
    for (int j = 0; j < 2 * order; ++j) {
        weights[j] = lanczos(-frac - (first + j), order);
        sum += weights[j];
    }
    if (sum != 0.0) {
        for (double &w : weights) w /= sum;
    }
}

void shiftImage(double const *in, std::size_t ny, std::size_t nx, double dx, double dy, int order,
                double *out) {
    int fx, fy;
    std::vector<double> wx, wy;
    lanczosWeights(dx, order, fx, wx);
    lanczosWeights(dy, order, fy, wy);
    std::vector<double> tmp(ny * nx);
    for (std::size_t r = 0; r < ny; ++r) {  // along x
        pass1d(in + r * nx, nx, 1, fx, wx, tmp.data() + r * nx, 1);
    }
    long const nyl = static_cast<long>(ny);
    for (long r = 0; r < nyl; ++r) {  // along y, row by row (contiguous)
        double *o = out + r * nx;
        std::fill(o, o + nx, 0.0);
        for (std::size_t j = 0; j < wy.size(); ++j) {
            long const m = r + fy + static_cast<long>(j);
            if (m < 0 || m >= nyl) continue;
            double const w = wy[j];
            double const *t = tmp.data() + m * nx;
            for (std::size_t c = 0; c < nx; ++c) o[c] += w * t[c];
        }
    }
}

void taperStamp(double *arr, std::size_t ny, std::size_t nx, double cx, double cy, double t0,
                double t1) {
    if (ny == 0 || nx == 0) return;
    std::vector<double> ring;
    ring.reserve(2 * (nx + ny));
    for (std::size_t c = 0; c < nx; ++c) {
        ring.push_back(arr[c]);
        ring.push_back(arr[(ny - 1) * nx + c]);
    }
    for (std::size_t r = 0; r < ny; ++r) {
        ring.push_back(arr[r * nx]);
        ring.push_back(arr[r * nx + nx - 1]);
    }
    // numpy-compatible median (mean of the two middle values for an even count)
    std::size_t const n = ring.size();
    std::nth_element(ring.begin(), ring.begin() + n / 2, ring.end());
    double med = ring[n / 2];
    if (n % 2 == 0) {
        double const lo = *std::max_element(ring.begin(), ring.begin() + n / 2);
        med = 0.5 * (med + lo);
    }
    double const h = 0.5 * (std::min(nx, ny) - 1.0);
    double const r0 = t0 * h, r1 = t1 * h;
    double const span = std::max(r1 - r0, 1e-9);
    for (std::size_t r = 0; r < ny; ++r) {
        for (std::size_t c = 0; c < nx; ++c) {
            double const dx = c - cx, dy = r - cy;
            double const rad2 = dx * dx + dy * dy;
            double f;
            if (rad2 <= r0 * r0 && r0 < r1) {
                f = 1.0;
            } else if (rad2 >= r1 * r1) {
                f = 0.0;
            } else {
                double const w = std::min(std::max((r1 - std::sqrt(rad2)) / span, 0.0), 1.0);
                f = 0.5 - 0.5 * std::cos(PI * w);
            }
            arr[r * nx + c] = (arr[r * nx + c] - med) * f;
        }
    }
}

void renderStamps(double const *kernel, std::size_t ny, std::size_t nx, double const *xs,
                  double const *ys, std::size_t n, int order, bool taper, double t0, double t1,
                  double *out, long *x0, long *y0) {
    long const hx = static_cast<long>(nx / 2), hy = static_cast<long>(ny / 2);
    for (std::size_t i = 0; i < n; ++i) {
        double const ix = std::floor(xs[i] + 0.5), iy = std::floor(ys[i] + 0.5);
        double const fx = xs[i] - ix, fy = ys[i] - iy;
        double *o = out + i * ny * nx;
        shiftImage(kernel, ny, nx, fx, fy, order, o);
        if (taper) taperStamp(o, ny, nx, hx + fx, hy + fy, t0, t1);
        x0[i] = static_cast<long>(ix) - hx;
        y0[i] = static_cast<long>(iy) - hy;
    }
}

void renderDerivatives(double const *kernel, std::size_t ny, std::size_t nx, double const *xs,
                       double const *ys, std::size_t n, int order, double eps, double *dx,
                       double *dy, long *x0, long *y0) {
    long const hx = static_cast<long>(nx / 2), hy = static_cast<long>(ny / 2);
    std::size_t const sz = ny * nx;
    std::vector<double> p(sz), m(sz);
    for (std::size_t i = 0; i < n; ++i) {
        double const ix = std::floor(xs[i] + 0.5), iy = std::floor(ys[i] + 0.5);
        double const fx = xs[i] - ix, fy = ys[i] - iy;
        shiftImage(kernel, ny, nx, fx + eps, fy, order, p.data());
        shiftImage(kernel, ny, nx, fx - eps, fy, order, m.data());
        for (std::size_t k = 0; k < sz; ++k) dx[i * sz + k] = (p[k] - m[k]) / (2 * eps);
        shiftImage(kernel, ny, nx, fx, fy + eps, order, p.data());
        shiftImage(kernel, ny, nx, fx, fy - eps, order, m.data());
        for (std::size_t k = 0; k < sz; ++k) dy[i * sz + k] = (p[k] - m[k]) / (2 * eps);
        x0[i] = static_cast<long>(ix) - hx;
        y0[i] = static_cast<long>(iy) - hy;
    }
}

}  // namespace crowd
}  // namespace meas
}  // namespace lsst
