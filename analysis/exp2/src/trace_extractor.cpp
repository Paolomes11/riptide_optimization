/* Copyright 2026 Giulio Mesini — Apache License 2.0 */

#include "trace_extractor.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace {

static double deg_to_rad(double deg) {
  return deg * (M_PI / 180.0);
}

static double rad_to_deg(double rad) {
  return rad * (180.0 / M_PI);
}

static double normalize_axial_deg(double deg) {
  double a = deg;
  while (a <= -90.0)
    a += 180.0;
  while (a > 90.0)
    a -= 180.0;
  return a;
}

static double axial_diff_deg(double a_deg, double b_deg) {
  const double a = deg_to_rad(a_deg);
  const double b = deg_to_rad(b_deg);
  const double d = 0.5 * std::atan2(std::sin(2.0 * (a - b)), std::cos(2.0 * (a - b)));
  return std::abs(rad_to_deg(d));
}

static double weighted_axial_mean_deg(double a_deg, double b_deg, double wa, double wb) {
  const double a = deg_to_rad(a_deg);
  const double b = deg_to_rad(b_deg);
  const double x = wa * std::cos(2.0 * a) + wb * std::cos(2.0 * b);
  const double y = wa * std::sin(2.0 * a) + wb * std::sin(2.0 * b);
  if (x == 0.0 && y == 0.0)
    return normalize_axial_deg(a_deg);
  return normalize_axial_deg(rad_to_deg(0.5 * std::atan2(y, x)));
}

struct NonZeroStats {
  double mean  = 0.0;
  double stdev = 0.0;
  bool ok      = false;
};

static NonZeroStats compute_positive_stats(const std::vector<double>& img) {
  double sum  = 0.0;
  double sum2 = 0.0;
  size_t n    = 0;
  for (double v : img) {
    if (v > 0.0) {
      sum += v;
      sum2 += v * v;
      ++n;
    }
  }
  if (n < 10)
    return {};

  const double mean = sum / static_cast<double>(n);
  const double var  = std::max(0.0, sum2 / static_cast<double>(n) - mean * mean);
  NonZeroStats s;
  s.mean  = mean;
  s.stdev = std::sqrt(var);
  s.ok    = std::isfinite(s.mean) && std::isfinite(s.stdev);
  return s;
}

struct WeightedMoments {
  double x_mean = 0.0;
  double y_mean = 0.0;
  double cxx    = 0.0;
  double cxy    = 0.0;
  double cyy    = 0.0;
  double sum_w  = 0.0;
  bool ok       = false;
};

static WeightedMoments compute_weighted_cov(const std::vector<double>& img, int width, int height,
                                            double threshold) {
  WeightedMoments m;
  if (width <= 0 || height <= 0)
    return m;

  double sum_w = 0.0;
  double sum_x = 0.0;
  double sum_y = 0.0;
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      const double v =
          img[static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)];
      if (v > threshold) {
        const double w = v;
        sum_w += w;
        sum_x += w * static_cast<double>(x);
        sum_y += w * static_cast<double>(y);
      }
    }
  }

  if (sum_w <= 0.0)
    return m;

  const double x_mean = sum_x / sum_w;
  const double y_mean = sum_y / sum_w;

  double cxx = 0.0;
  double cxy = 0.0;
  double cyy = 0.0;
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      const double v =
          img[static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)];
      if (v > threshold) {
        const double w  = v;
        const double dx = static_cast<double>(x) - x_mean;
        const double dy = static_cast<double>(y) - y_mean;
        cxx += w * dx * dx;
        cxy += w * dx * dy;
        cyy += w * dy * dy;
      }
    }
  }

  m.x_mean = x_mean;
  m.y_mean = y_mean;
  m.cxx    = cxx / sum_w;
  m.cxy    = cxy / sum_w;
  m.cyy    = cyy / sum_w;
  m.sum_w  = sum_w;
  m.ok     = std::isfinite(m.cxx) && std::isfinite(m.cxy) && std::isfinite(m.cyy);
  return m;
}

static double angle_from_cov_pca(const WeightedMoments& m) {
  if (!m.ok)
    return 0.0;
  const double theta = 0.5 * std::atan2(2.0 * m.cxy, m.cxx - m.cyy);
  return normalize_axial_deg(rad_to_deg(theta));
}

static double angle_from_inertia(const std::vector<double>& img, int width, int height,
                                 const WeightedMoments& m, double threshold) {
  if (!m.ok)
    return 0.0;

  double mxx = 0.0;
  double myy = 0.0;
  double mxy = 0.0;
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      const double v =
          img[static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)];
      if (v > threshold) {
        const double w  = v;
        const double dx = static_cast<double>(x) - m.x_mean;
        const double dy = static_cast<double>(y) - m.y_mean;
        mxx += w * dx * dx;
        myy += w * dy * dy;
        mxy += w * dx * dy;
      }
    }
  }

  const double theta = 0.5 * std::atan2(2.0 * mxy, mxx - myy);
  return normalize_axial_deg(rad_to_deg(theta));
}

static double safe_var_from_threshold_sweep(double a0, double a_minus, double a_plus) {
  const double dm  = axial_diff_deg(a0, a_minus);
  const double dp  = axial_diff_deg(a0, a_plus);
  const double d   = std::max(dm, dp);
  const double var = d * d;
  if (!std::isfinite(var) || var <= 1e-12)
    return 1e-12;
  return var;
}

static double bilinear_sample(const std::vector<double>& img, int width, int height, double x,
                              double y) {
  if (width <= 0 || height <= 0)
    return 0.0;
  if (x < 0.0 || y < 0.0 || x > static_cast<double>(width - 1)
      || y > static_cast<double>(height - 1))
    return 0.0;

  const int x0 = static_cast<int>(std::floor(x));
  const int y0 = static_cast<int>(std::floor(y));
  const int x1 = std::min(x0 + 1, width - 1);
  const int y1 = std::min(y0 + 1, height - 1);

  const double fx = x - static_cast<double>(x0);
  const double fy = y - static_cast<double>(y0);

  auto at = [&](int xx, int yy) -> double {
    return img[static_cast<size_t>(yy) * static_cast<size_t>(width) + static_cast<size_t>(xx)];
  };

  const double v00 = at(x0, y0);
  const double v10 = at(x1, y0);
  const double v01 = at(x0, y1);
  const double v11 = at(x1, y1);

  const double v0 = (1.0 - fx) * v00 + fx * v10;
  const double v1 = (1.0 - fx) * v01 + fx * v11;
  return (1.0 - fy) * v0 + fy * v1;
}

static double median_of(std::vector<double> v) {
  if (v.empty())
    return 0.0;
  const size_t mid = v.size() / 2;
  std::nth_element(v.begin(), v.begin() + static_cast<long>(mid), v.end());
  const double m = v[mid];
  if (v.size() % 2 == 1)
    return m;
  std::nth_element(v.begin(), v.begin() + static_cast<long>(mid - 1), v.end());
  return 0.5 * (m + v[mid - 1]);
}

static double stdev_of(const std::vector<double>& v, double mean) {
  if (v.size() < 2)
    return 0.0;
  double sum = 0.0;
  for (double x : v) {
    const double d = x - mean;
    sum += d * d;
  }
  return std::sqrt(sum / static_cast<double>(v.size() - 1));
}

struct GaussFit {
  double A            = 0.0;
  double mu           = 0.0;
  double mu_err       = std::numeric_limits<double>::infinity();
  double sigma        = 1.0;
  double B            = 0.0;
  double sigma_err    = std::numeric_limits<double>::infinity();
  double chi2_ndof    = std::numeric_limits<double>::infinity();
  double rms_residual = std::numeric_limits<double>::infinity();
  double noise        = std::numeric_limits<double>::infinity();
  bool converged      = false;
};

static bool solve_4x4(std::array<std::array<double, 5>, 4>& a) {
  for (size_t i = 0; i < 4; ++i) {
    size_t piv  = i;
    double maxv = std::abs(a[i][i]);
    for (size_t r = i + 1; r < 4; ++r) {
      const double v = std::abs(a[r][i]);
      if (v > maxv) {
        maxv = v;
        piv  = r;
      }
    }
    if (maxv < 1e-18)
      return false;
    if (piv != i)
      std::swap(a[piv], a[i]);

    const double div = a[i][i];
    for (size_t c = i; c < 5; ++c)
      a[i][c] /= div;

    for (size_t r = 0; r < 4; ++r) {
      if (r == i)
        continue;
      const double f = a[r][i];
      for (size_t c = i; c < 5; ++c)
        a[r][c] -= f * a[i][c];
    }
  }
  return true;
}

static bool invert_4x4(const std::array<std::array<double, 4>, 4>& m,
                       std::array<std::array<double, 4>, 4>& inv) {
  std::array<std::array<double, 8>, 4> a{};
  for (size_t r = 0; r < 4; ++r) {
    for (size_t c = 0; c < 4; ++c)
      a[r][c] = m[r][c];
    for (size_t c = 0; c < 4; ++c)
      a[r][4 + c] = (r == c) ? 1.0 : 0.0;
  }

  for (size_t i = 0; i < 4; ++i) {
    size_t piv  = i;
    double maxv = std::abs(a[i][i]);
    for (size_t r = i + 1; r < 4; ++r) {
      const double v = std::abs(a[r][i]);
      if (v > maxv) {
        maxv = v;
        piv  = r;
      }
    }
    if (maxv < 1e-18)
      return false;
    if (piv != i)
      std::swap(a[piv], a[i]);

    const double div = a[i][i];
    for (size_t c = i; c < 8; ++c)
      a[i][c] /= div;

    for (size_t r = 0; r < 4; ++r) {
      if (r == i)
        continue;
      const double f = a[r][i];
      for (size_t c = i; c < 8; ++c)
        a[r][c] -= f * a[i][c];
    }
  }

  for (size_t r = 0; r < 4; ++r)
    for (size_t c = 0; c < 4; ++c)
      inv[r][c] = a[r][4 + c];
  return true;
}

static GaussFit fit_gaussian_1d(const std::vector<double>& x, const std::vector<double>& y,
                                double sigma_max) {
  GaussFit out;
  const int n = static_cast<int>(x.size());
  if (n < 7)
    return out;

  std::vector<double> y_edges;
  y_edges.reserve(static_cast<size_t>(n));
  const int edge_n = std::max(1, n / 5);
  for (int i = 0; i < edge_n; ++i) {
    y_edges.push_back(y[static_cast<size_t>(i)]);
    y_edges.push_back(y[static_cast<size_t>(n - 1 - i)]);
  }
  const double B0 = median_of(y_edges);

  double ymax = y[0];
  int imax    = 0;
  for (int i = 1; i < n; ++i) {
    if (y[static_cast<size_t>(i)] > ymax) {
      ymax = y[static_cast<size_t>(i)];
      imax = i;
    }
  }
  const double A0 = std::max(0.0, ymax - B0);

  double wsum = 0.0;
  double xsum = 0.0;
  for (int i = 0; i < n; ++i) {
    const double w = std::max(0.0, y[static_cast<size_t>(i)] - B0);
    wsum += w;
    xsum += w * x[static_cast<size_t>(i)];
  }
  const double mu0 = (wsum > 0.0) ? (xsum / wsum) : x[static_cast<size_t>(imax)];

  double var = 0.0;
  if (wsum > 0.0) {
    for (int i = 0; i < n; ++i) {
      const double w = std::max(0.0, y[static_cast<size_t>(i)] - B0);
      const double d = x[static_cast<size_t>(i)] - mu0;
      var += w * d * d;
    }
    var /= wsum;
  }
  double sigma0 = std::sqrt(std::max(1e-6, var));
  sigma0        = std::clamp(sigma0, 0.5, sigma_max);

  out.A     = A0;
  out.mu    = mu0;
  out.sigma = sigma0;
  out.B     = B0;

  auto eval = [&](double A, double mu, double sigma, double B, int i) -> double {
    const double dx = x[static_cast<size_t>(i)] - mu;
    const double e  = std::exp(-(dx * dx) / (2.0 * sigma * sigma));
    return A * e + B;
  };

  double lambda    = 1e-3;
  double best_chi2 = std::numeric_limits<double>::infinity();
  GaussFit best    = out;
  bool has_best    = false;

  for (int iter = 0; iter < 100; ++iter) {
    std::array<std::array<double, 4>, 4> jtj{};
    std::array<double, 4> jtr{};
    double chi2 = 0.0;

    for (int i = 0; i < n; ++i) {
      const double xi = x[static_cast<size_t>(i)];
      const double yi = y[static_cast<size_t>(i)];
      const double dx = xi - out.mu;
      const double s2 = out.sigma * out.sigma;
      const double e  = std::exp(-(dx * dx) / (2.0 * s2));
      const double fi = out.A * e + out.B;
      const double ri = yi - fi;
      chi2 += ri * ri;

      const double dA     = e;
      const double dmu    = out.A * e * (dx / s2);
      const double dsigma = out.A * e * ((dx * dx) / (out.sigma * s2));
      const double dB     = 1.0;

      const std::array<double, 4> J{dA, dmu, dsigma, dB};
      for (int r = 0; r < 4; ++r) {
        jtr[static_cast<size_t>(r)] += J[static_cast<size_t>(r)] * ri;
        for (int c = 0; c < 4; ++c)
          jtj[static_cast<size_t>(r)][static_cast<size_t>(c)] +=
              J[static_cast<size_t>(r)] * J[static_cast<size_t>(c)];
      }
    }

    if (!std::isfinite(chi2))
      break;

    if (!has_best || chi2 < best_chi2) {
      best_chi2 = chi2;
      best      = out;
      has_best  = true;
    }

    std::array<std::array<double, 5>, 4> a{};
    for (int r = 0; r < 4; ++r) {
      for (int c = 0; c < 4; ++c)
        a[static_cast<size_t>(r)][static_cast<size_t>(c)] =
            jtj[static_cast<size_t>(r)][static_cast<size_t>(c)];
      a[static_cast<size_t>(r)][static_cast<size_t>(r)] *= (1.0 + lambda);
      a[static_cast<size_t>(r)][4] = jtr[static_cast<size_t>(r)];
    }

    if (!solve_4x4(a))
      break;

    const double dA  = a[0][4];
    const double dmu = a[1][4];
    const double ds  = a[2][4];
    const double dB  = a[3][4];

    GaussFit trial = out;
    trial.A        = out.A + dA;
    trial.mu       = out.mu + dmu;
    trial.sigma    = out.sigma + ds;
    trial.B        = out.B + dB;

    if (!(trial.sigma > 1e-3) || !std::isfinite(trial.sigma))
      trial.sigma = std::max(1e-3, out.sigma);
    trial.sigma = std::min(trial.sigma, sigma_max);
    if (!std::isfinite(trial.A))
      break;

    double chi2_trial = 0.0;
    for (int i = 0; i < n; ++i) {
      const double ri =
          y[static_cast<size_t>(i)] - eval(trial.A, trial.mu, trial.sigma, trial.B, i);
      chi2_trial += ri * ri;
    }

    if (chi2_trial < chi2) {
      out = trial;
      lambda *= 0.3;
      lambda = std::max(1e-12, lambda);

      const double rel = std::max(
          {std::abs(dA) / (std::abs(out.A) + 1e-9), std::abs(dmu) / (std::abs(out.mu) + 1e-9),
           std::abs(ds) / (std::abs(out.sigma) + 1e-9), std::abs(dB) / (std::abs(out.B) + 1e-9)});
      if (rel < 1e-6) {
        out.converged = true;
        break;
      }
    } else {
      lambda *= 10.0;
      if (lambda > 1e12)
        break;
    }
  }

  if (has_best && !out.converged)
    out = best;

  double chi2 = 0.0;
  std::vector<double> residuals;
  residuals.reserve(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) {
    const double ri = y[static_cast<size_t>(i)] - eval(out.A, out.mu, out.sigma, out.B, i);
    residuals.push_back(ri);
    chi2 += ri * ri;
  }
  out.rms_residual = std::sqrt(chi2 / static_cast<double>(n));

  const int ndof  = n - 4;
  const double s2 = (ndof > 0) ? (chi2 / static_cast<double>(ndof)) : chi2;

  std::vector<double> y_tail;
  y_tail.reserve(static_cast<size_t>(n));
  const double absmax   = std::max(std::abs(x.front()), std::abs(x.back()));
  const double tail_thr = 0.7 * absmax;
  for (int i = 0; i < n; ++i)
    if (std::abs(x[static_cast<size_t>(i)]) >= tail_thr)
      y_tail.push_back(y[static_cast<size_t>(i)] - out.B);
  const double tail_mean = (y_tail.empty()) ? 0.0
                                            : (std::accumulate(y_tail.begin(), y_tail.end(), 0.0)
                                               / static_cast<double>(y_tail.size()));
  const double noise     = std::max(1e-9, stdev_of(y_tail, tail_mean));
  out.noise              = noise;
  double chi2n           = 0.0;
  for (double r : residuals)
    chi2n += (r * r) / (noise * noise);
  out.chi2_ndof = (ndof > 0) ? (chi2n / static_cast<double>(ndof)) : chi2n;

  std::array<std::array<double, 4>, 4> jtj{};
  for (int i = 0; i < n; ++i) {
    const double xi     = x[static_cast<size_t>(i)];
    const double dx     = xi - out.mu;
    const double s2fit  = out.sigma * out.sigma;
    const double e      = std::exp(-(dx * dx) / (2.0 * s2fit));
    const double dA     = e;
    const double dmu    = out.A * e * (dx / s2fit);
    const double dsigma = out.A * e * ((dx * dx) / (out.sigma * s2fit));
    const double dB     = 1.0;
    const std::array<double, 4> J{dA, dmu, dsigma, dB};
    for (int r = 0; r < 4; ++r)
      for (int c = 0; c < 4; ++c)
        jtj[static_cast<size_t>(r)][static_cast<size_t>(c)] +=
            J[static_cast<size_t>(r)] * J[static_cast<size_t>(c)];
  }

  std::array<std::array<double, 4>, 4> cov{};
  if (invert_4x4(jtj, cov)) {
    const double var_mu = cov[1][1] * s2;
    if (var_mu > 0.0 && std::isfinite(var_mu))
      out.mu_err = std::sqrt(var_mu);
    const double var_sigma = cov[2][2] * s2;
    if (var_sigma > 0.0 && std::isfinite(var_sigma))
      out.sigma_err = std::sqrt(var_sigma);
  }

  return out;
}

} // namespace

namespace riptide::exp2 {

double estimate_trace_angle(const std::vector<double>& img, int width, int height,
                            const TraceConfig& cfg) {
  return estimate_trace_angle(img, width, height, cfg, nullptr);
}

double estimate_trace_angle(const std::vector<double>& img, int width, int height,
                            const TraceConfig&, double* angle_err_deg_out) {
  if (static_cast<size_t>(std::max(0, width) * std::max(0, height)) != img.size())
    throw std::invalid_argument("estimate_trace_angle: dimensioni immagine incoerenti");

  const NonZeroStats s = compute_positive_stats(img);
  if (!s.ok) {
    if (angle_err_deg_out)
      *angle_err_deg_out = 0.0;
    return 0.0;
  }

  const double thr0      = s.mean + 2.0 * s.stdev;
  const double thr_minus = std::max(0.0, s.mean + 1.0 * s.stdev);
  const double thr_plus  = s.mean + 3.0 * s.stdev;

  auto angle_at = [&](double thr, bool warn) -> double {
    const WeightedMoments m = compute_weighted_cov(img, width, height, thr);
    if (!m.ok)
      return 0.0;
    const double a_pca = angle_from_cov_pca(m);
    const double a_I   = angle_from_inertia(img, width, height, m, thr);
    if (warn && axial_diff_deg(a_pca, a_I) > 5.0) {
      std::cerr << "[WARNING] exp2::estimate_trace_angle: PCA e inerzia differiscono di "
                << axial_diff_deg(a_pca, a_I) << " deg\n";
    }

    const double a_pca_m = angle_from_cov_pca(compute_weighted_cov(img, width, height, thr_minus));
    const double a_pca_p = angle_from_cov_pca(compute_weighted_cov(img, width, height, thr_plus));
    const double var_pca = safe_var_from_threshold_sweep(a_pca, a_pca_m, a_pca_p);

    const WeightedMoments m_minus = compute_weighted_cov(img, width, height, thr_minus);
    const WeightedMoments m_plus  = compute_weighted_cov(img, width, height, thr_plus);
    const double a_I_m            = angle_from_inertia(img, width, height, m_minus, thr_minus);
    const double a_I_p            = angle_from_inertia(img, width, height, m_plus, thr_plus);
    const double var_I            = safe_var_from_threshold_sweep(a_I, a_I_m, a_I_p);

    const double w_pca = 1.0 / var_pca;
    const double w_I   = 1.0 / var_I;
    return weighted_axial_mean_deg(a_pca, a_I, w_pca, w_I);
  };

  const double a0 = angle_at(thr0, true);
  const double am = angle_at(thr_minus, false);
  const double ap = angle_at(thr_plus, false);

  const double err = std::max(axial_diff_deg(a0, am), axial_diff_deg(a0, ap));
  if (angle_err_deg_out)
    *angle_err_deg_out = err;
  return normalize_axial_deg(a0);
}

TraceResult extract_trace_profile(const std::vector<double>& img, int width, int height,
                                  double angle_deg, const TraceConfig& cfg) {
  if (static_cast<size_t>(std::max(0, width) * std::max(0, height)) != img.size())
    throw std::invalid_argument("extract_trace_profile: dimensioni immagine incoerenti");

  TraceResult out;
  out.angle_deg = normalize_axial_deg(angle_deg);

  const NonZeroStats s = compute_positive_stats(img);
  const double thr     = s.ok ? (s.mean + 2.0 * s.stdev) : 0.0;
  WeightedMoments m    = compute_weighted_cov(img, width, height, thr);

  if (!m.ok) {
    m.x_mean = 0.5 * static_cast<double>(width - 1);
    m.y_mean = 0.5 * static_cast<double>(height - 1);
    m.ok     = true;
  }

  const double theta = deg_to_rad(out.angle_deg);
  const double ux    = std::cos(theta);
  const double uy    = std::sin(theta);
  const double nx    = -uy;
  const double ny    = ux;

  double tmin = std::numeric_limits<double>::infinity();
  double tmax = -std::numeric_limits<double>::infinity();
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      const double v =
          img[static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)];
      if (v > thr) {
        const double dx = static_cast<double>(x) - m.x_mean;
        const double dy = static_cast<double>(y) - m.y_mean;
        const double t  = dx * ux + dy * uy;
        tmin            = std::min(tmin, t);
        tmax            = std::max(tmax, t);
      }
    }
  }

  if (!std::isfinite(tmin) || !std::isfinite(tmax) || tmax <= tmin) {
    const std::array<std::pair<double, double>, 4> corners{{
        {0.0 - m.x_mean, 0.0 - m.y_mean},
        {static_cast<double>(width - 1) - m.x_mean, 0.0 - m.y_mean},
        {0.0 - m.x_mean, static_cast<double>(height - 1) - m.y_mean},
        {static_cast<double>(width - 1) - m.x_mean, static_cast<double>(height - 1) - m.y_mean},
    }};
    tmin = std::numeric_limits<double>::infinity();
    tmax = -std::numeric_limits<double>::infinity();
    for (const auto& c : corners) {
      const double t = c.first * ux + c.second * uy;
      tmin           = std::min(tmin, t);
      tmax           = std::max(tmax, t);
    }
  }

  const double t0   = tmin;
  const double step = static_cast<double>(std::max(1, cfg.slice_step));

  for (double t = tmin; t <= tmax; t += step) {
    SliceProfile sp;
    sp.t = t - t0;

    const double x0 = m.x_mean + ux * t;
    const double y0 = m.y_mean + uy * t;

    double s_lo = -static_cast<double>(cfg.gaussian_range);
    double s_hi = static_cast<double>(cfg.gaussian_range);

    auto tighten = [&](double ncomp, double coord, double lo_bound, double hi_bound) {
      if (std::abs(ncomp) < 1e-12)
        return;
      const double s1 = (lo_bound - coord) / ncomp;
      const double s2 = (hi_bound - coord) / ncomp;
      const double a  = std::min(s1, s2);
      const double b  = std::max(s1, s2);
      s_lo            = std::max(s_lo, a);
      s_hi            = std::min(s_hi, b);
    };

    tighten(nx, x0, 0.0, static_cast<double>(width - 1));
    tighten(ny, y0, 0.0, static_cast<double>(height - 1));

    const double max_abs = std::min({static_cast<double>(cfg.gaussian_range), s_hi, -s_lo});
    const int max_i      = static_cast<int>(std::floor(max_abs));
    if (max_i < 5) {
      sp.valid     = false;
      sp.near_edge = true;
      out.profile.push_back(sp);
      continue;
    }

    sp.near_edge = (max_i < cfg.gaussian_range);

    std::vector<double> xs;
    std::vector<double> ys;
    xs.reserve(static_cast<size_t>(2 * max_i + 1));
    ys.reserve(static_cast<size_t>(2 * max_i + 1));

    for (int si = -max_i; si <= max_i; ++si) {
      const double ss = static_cast<double>(si);
      xs.push_back(ss);
      double acc   = 0.0;
      int n_acc    = 0;
      const int hw = std::max(0, cfg.slice_width);
      for (int oi = -hw; oi <= hw; ++oi) {
        const double oo = static_cast<double>(oi);
        const double xx = x0 + ux * oo + nx * ss;
        const double yy = y0 + uy * oo + ny * ss;
        const double vv = bilinear_sample(img, width, height, xx, yy);
        acc += vv;
        ++n_acc;
      }
      ys.push_back((n_acc > 0) ? (acc / static_cast<double>(n_acc)) : 0.0);
    }

    GaussFit fit = fit_gaussian_1d(xs, ys, cfg.sigma_max);
    if (!fit.converged) {
      const double B0 = median_of(ys);
      double wsum     = 0.0;
      double xsum     = 0.0;
      for (size_t i = 0; i < xs.size(); ++i) {
        const double w = std::max(0.0, ys[i] - B0);
        wsum += w;
        xsum += w * xs[i];
      }
      const double mu0 = (wsum > 0.0) ? (xsum / wsum) : 0.0;
      double var       = 0.0;
      if (wsum > 0.0) {
        for (size_t i = 0; i < xs.size(); ++i) {
          const double w = std::max(0.0, ys[i] - B0);
          const double d = xs[i] - mu0;
          var += w * d * d;
        }
        var /= wsum;
      }
      fit.mu           = mu0;
      fit.mu_err       = std::numeric_limits<double>::infinity();
      fit.sigma        = std::clamp(std::sqrt(std::max(1e-6, var)), 0.5, cfg.sigma_max);
      fit.A            = *std::max_element(ys.begin(), ys.end()) - B0;
      fit.B            = B0;
      fit.sigma_err    = std::numeric_limits<double>::infinity();
      fit.chi2_ndof    = std::numeric_limits<double>::infinity();
      fit.rms_residual = std::numeric_limits<double>::infinity();
    }

    sp.amplitude = fit.A;
    sp.center    = fit.mu;
    {
      double ce = fit.mu_err;
      if (!std::isfinite(ce))
        ce = std::numeric_limits<double>::infinity();
      ce = ce * cfg.center_err_scale;
      if (!std::isfinite(ce) || ce <= 0.0)
        ce = cfg.center_err_floor;
      sp.center_err = std::max(cfg.center_err_floor, ce);
    }
    sp.sigma = fit.sigma;
    {
      double se = fit.sigma_err;
      if (!std::isfinite(se))
        se = std::numeric_limits<double>::infinity();
      se = se * cfg.sigma_err_scale;
      if (!std::isfinite(se) || se <= 0.0)
        se = cfg.sigma_err_floor;
      sp.sigma_err = std::max(cfg.sigma_err_floor, se);
    }
    sp.chi2_ndof = fit.chi2_ndof;

    const double snr = (fit.noise > 0.0 && std::isfinite(fit.noise)) ? (fit.A / fit.noise) : 0.0;
    sp.snr           = snr;
    sp.valid = std::isfinite(sp.sigma) && std::isfinite(sp.sigma_err) && (sp.sigma_err > 0.0)
            && (sp.chi2_ndof < 5.0) && (sp.sigma > 0.5) && (sp.sigma < cfg.sigma_max)
            && (snr > cfg.min_snr);

    out.profile.push_back(sp);
  }

  if (cfg.enable_trace_trim && !out.profile.empty()) {
    double snr_max = 0.0;
    for (const auto& sp : out.profile)
      if (std::isfinite(sp.snr))
        snr_max = std::max(snr_max, sp.snr);

    const double snr_thr = std::max(cfg.trace_trim_min_snr, cfg.trace_trim_frac * snr_max);

    int best_start = -1;
    int best_len   = 0;
    int cur_start  = -1;
    int cur_len    = 0;

    for (int i = 0; i < static_cast<int>(out.profile.size()); ++i) {
      const auto& sp  = out.profile[static_cast<size_t>(i)];
      const bool keep = std::isfinite(sp.snr) && (sp.snr >= snr_thr) && std::isfinite(sp.sigma)
                     && (sp.sigma > 0.5) && (sp.sigma < cfg.sigma_max);

      if (keep) {
        if (cur_start < 0) {
          cur_start = i;
          cur_len   = 1;
        } else {
          ++cur_len;
        }
        if (cur_len > best_len) {
          best_len   = cur_len;
          best_start = cur_start;
        }
      } else {
        cur_start = -1;
        cur_len   = 0;
      }
    }

    if (best_len >= cfg.trace_trim_min_slices && best_start >= 0) {
      const int pad = std::max(0, cfg.trace_trim_pad_slices);
      int start     = std::max(0, best_start - pad);
      int end = std::min(static_cast<int>(out.profile.size()) - 1, best_start + best_len - 1 + pad);

      for (int i = 0; i < static_cast<int>(out.profile.size()); ++i) {
        if (i < start || i > end) {
          out.profile[static_cast<size_t>(i)].in_trace = false;
          out.profile[static_cast<size_t>(i)].valid    = false;
        }
      }
    }
  }

  double wsum = 0.0;
  double ssum = 0.0;
  int n_valid = 0;
  for (const auto& sp : out.profile) {
    if (sp.valid && !sp.near_edge && std::isfinite(sp.sigma_err) && sp.sigma_err > 0.0) {
      const double w = 1.0 / (sp.sigma_err * sp.sigma_err);
      wsum += w;
      ssum += w * sp.sigma;
      ++n_valid;
    }
  }

  out.n_valid_slices = n_valid;
  if (wsum > 0.0 && n_valid >= 5) {
    out.sigma_mean     = ssum / wsum;
    out.sigma_mean_err = std::sqrt(1.0 / wsum);
    out.fwhm_mean      = 2.355 * out.sigma_mean;
  } else {
    out.sigma_mean     = std::numeric_limits<double>::quiet_NaN();
    out.sigma_mean_err = std::numeric_limits<double>::quiet_NaN();
    out.fwhm_mean      = std::numeric_limits<double>::quiet_NaN();
  }

  return out;
}

CentroidFitResult fit_centroid_line(const TraceResult& trace) {
  if (trace.profile.size() < 3)
    throw std::invalid_argument("fit_centroid_line: profilo troppo corto");

  std::vector<double> dt;
  dt.reserve(trace.profile.size());
  for (size_t i = 1; i < trace.profile.size(); ++i)
    dt.push_back(std::abs(trace.profile[i].t - trace.profile[i - 1].t));
  const double dt_med = median_of(dt);
  const double var_z  = (dt_med > 0.0) ? (0.25 * dt_med * dt_med) : 1.0;

  std::vector<riptide::TracePoint> pts;
  pts.reserve(trace.profile.size());

  for (const auto& sp : trace.profile) {
    riptide::TracePoint p{};
    p.t                = sp.t;
    p.r                = 0.0;
    p.x_src            = 0.0;
    p.y_src            = 0.0;
    p.z_src            = 0.0;
    p.mu_y             = sp.center;
    p.mu_z             = sp.t;
    const double var_y = (std::isfinite(sp.center_err) && sp.center_err > 0.0)
                           ? (sp.center_err * sp.center_err)
                           : 1.0;
    p.cov              = {var_y, 0.0, var_z};
    p.valid            = sp.valid && !sp.near_edge;
    p.n_hits           = 1000.0;
    p.n_hits_count     = 1000.0;
    pts.push_back(p);
  }

  auto count_valid = [&]() -> int {
    int n = 0;
    for (const auto& p : pts)
      if (p.valid)
        ++n;
    return n;
  };

  auto pull_of = [&](const riptide::TracePoint& p, const CentroidFitResult& fit) -> double {
    const double y      = p.mu_y;
    const double z      = p.mu_z;
    const double var_y  = p.cov.yy;
    const double var_zp = p.cov.zz;

    double u     = 0.0;
    double v     = 0.0;
    double var_u = 0.0;
    double var_v = 0.0;
    if (fit.axis == riptide::FitAxis::ZvsY) {
      u     = y;
      v     = z;
      var_u = var_y;
      var_v = var_zp;
    } else {
      u     = z;
      v     = y;
      var_u = var_zp;
      var_v = var_y;
    }

    const double denom   = std::sqrt(1.0 + fit.a * fit.a);
    const double resid   = (v - fit.a * u - fit.b) / denom;
    const double sigma_d = std::sqrt((fit.a * fit.a * var_u + var_v) / (1.0 + fit.a * fit.a));
    if (!(sigma_d > 0.0) || !std::isfinite(sigma_d))
      return 0.0;
    return resid / sigma_d;
  };

  CentroidFitResult fit = riptide::fit_trace(pts, 0.0);

  const double pull_max = 8.0;
  const int max_iter    = 2;
  const int min_points  = 30;

  for (int iter = 0; iter < max_iter; ++iter) {
    int removed = 0;
    for (auto& p : pts) {
      if (!p.valid)
        continue;
      const double pull = pull_of(p, fit);
      if (std::abs(pull) > pull_max) {
        p.valid = false;
        ++removed;
      }
    }
    if (removed == 0)
      break;
    if (count_valid() < min_points)
      break;
    fit = riptide::fit_trace(pts, 0.0);
  }

  return fit;
}

} // namespace riptide::exp2
