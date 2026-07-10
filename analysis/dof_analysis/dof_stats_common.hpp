#pragma once

#include "plot_style_common.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>
#include <vector>

// Funzioni statistiche condivise byte-identiche tra dof_map.cpp,
// dof_plot.cpp e magnification_map.cpp. Lo stile ROOT condiviso vive in
// plot_style_common.hpp; dof_plot.cpp ha una propria variante piu' corta e
// resta locale a quel file.

inline double weighted_std(const std::vector<double>& x, const std::vector<double>& w) {
  double sum_w = 0.0;
  double sum_x = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    sum_w += w[i];
    sum_x += w[i] * x[i];
  }
  if (sum_w <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  double mean = sum_x / sum_w;
  double var  = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    double d = x[i] - mean;
    var += w[i] * d * d;
  }
  var /= sum_w;
  return std::sqrt(std::max(0.0, var));
}

inline double weighted_mean(const std::vector<double>& x, const std::vector<double>& w) {
  double sum_w = 0.0;
  double sum_x = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    sum_w += w[i];
    sum_x += w[i] * x[i];
  }
  if (sum_w <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return sum_x / sum_w;
}

inline std::vector<size_t> select_core_indices_around_mean(const std::vector<double>& y,
                                                            const std::vector<double>& z,
                                                            const std::vector<double>& w,
                                                            double core_fraction) {
  if (y.empty() || z.size() != y.size() || w.size() != y.size()) {
    return {};
  }
  if (core_fraction >= 1.0) {
    std::vector<size_t> idx(y.size());
    std::iota(idx.begin(), idx.end(), 0);
    return idx;
  }
  if (core_fraction <= 0.0) {
    return {};
  }

  double sum_w = 0.0;
  double sum_y = 0.0;
  double sum_z = 0.0;
  for (size_t i = 0; i < y.size(); ++i) {
    sum_w += w[i];
    sum_y += w[i] * y[i];
    sum_z += w[i] * z[i];
  }
  if (sum_w <= 0.0) {
    return {};
  }
  double mu_y = sum_y / sum_w;
  double mu_z = sum_z / sum_w;

  double cov_yy = 0.0;
  double cov_zz = 0.0;
  double cov_yz = 0.0;
  for (size_t i = 0; i < y.size(); ++i) {
    double dy = y[i] - mu_y;
    double dz = z[i] - mu_z;
    cov_yy += w[i] * dy * dy;
    cov_zz += w[i] * dz * dz;
    cov_yz += w[i] * dy * dz;
  }
  cov_yy /= sum_w;
  cov_zz /= sum_w;
  cov_yz /= sum_w;

  double scale = 0.5 * (cov_yy + cov_zz);
  double eps   = 1e-9 + 1e-6 * std::max(0.0, scale);
  cov_yy += eps;
  cov_zz += eps;

  double det = cov_yy * cov_zz - cov_yz * cov_yz;
  if (!std::isfinite(det) || std::abs(det) < 1e-18) {
    cov_yz = 0.0;
    det    = cov_yy * cov_zz;
    if (!std::isfinite(det) || std::abs(det) < 1e-18) {
      return {};
    }
  }

  double inv_yy = cov_zz / det;
  double inv_zz = cov_yy / det;
  double inv_yz = -cov_yz / det;

  std::vector<size_t> idx(y.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
    double dya = y[a] - mu_y;
    double dza = z[a] - mu_z;
    double dyb = y[b] - mu_y;
    double dzb = z[b] - mu_z;
    double d2a = dya * dya * inv_yy + dza * dza * inv_zz + 2.0 * dya * dza * inv_yz;
    double d2b = dyb * dyb * inv_yy + dzb * dzb * inv_zz + 2.0 * dyb * dzb * inv_yz;
    return d2a < d2b;
  });

  std::vector<size_t> keep;
  keep.reserve(idx.size());
  double target = core_fraction * sum_w;
  double acc    = 0.0;
  for (size_t k : idx) {
    keep.push_back(k);
    acc += w[k];
    if (acc >= target) {
      break;
    }
  }
  std::sort(keep.begin(), keep.end());
  return keep;
}

inline void apply_selection(std::vector<double>& y0, std::vector<double>& z0,
                            std::vector<double>& dy, std::vector<double>& dz,
                            std::vector<double>& w, std::vector<double>& ysrc,
                            const std::vector<size_t>& keep) {
  auto pick = [&](std::vector<double>& v) {
    std::vector<double> out;
    out.reserve(keep.size());
    for (size_t i : keep) {
      out.push_back(v[i]);
    }
    v.swap(out);
  };
  pick(y0);
  pick(z0);
  pick(dy);
  pick(dz);
  pick(w);
  pick(ysrc);
}

inline double weighted_percentile(const std::vector<double>& x, const std::vector<double>& w,
                                  double p) {
  if (x.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  std::vector<size_t> idx(x.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return x[a] < x[b]; });

  double sum_w = 0.0;
  for (double wi : w) {
    sum_w += wi;
  }
  if (sum_w <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double target = p * sum_w;
  double acc    = 0.0;
  for (size_t j = 0; j < idx.size(); ++j) {
    acc += w[idx[j]];
    if (acc >= target) {
      return x[idx[j]];
    }
  }
  return x[idx.back()];
}

inline std::optional<double> quadratic_vertex_from_three_points(double x0, double y0, double x1,
                                                                 double y1, double x2, double y2) {
  double d0 = (x0 - x1) * (x0 - x2);
  double d1 = (x1 - x0) * (x1 - x2);
  double d2 = (x2 - x0) * (x2 - x1);
  if (std::abs(d0) < 1e-12 || std::abs(d1) < 1e-12 || std::abs(d2) < 1e-12) {
    return std::nullopt;
  }

  double a = (y0 / d0) + (y1 / d1) + (y2 / d2);
  double b = (-y0 * (x1 + x2) / d0) + (-y1 * (x0 + x2) / d1) + (-y2 * (x0 + x1) / d2);
  if (!std::isfinite(a) || !std::isfinite(b) || std::abs(a) < 1e-18 || a <= 0.0) {
    return std::nullopt;
  }
  return -b / (2.0 * a);
}
