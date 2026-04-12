/*
 * Copyright 2026 Giulio Mesini
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 */

#include "psf_interpolator.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

static int g_n_pass = 0;
static int g_n_fail = 0;

static bool check(bool condition, const std::string& msg, const std::string& test_name = "") {
  if (condition) {
    std::cout << "  [PASS] " << msg << "\n";
    ++g_n_pass;
  } else {
    std::cerr << "  [FAIL] " << msg;
    if (!test_name.empty())
      std::cerr << "  (in " << test_name << ")";
    std::cerr << "\n";
    ++g_n_fail;
  }
  return condition;
}

static bool near(double a, double b, double tol = 1e-6, const std::string& label = "",
                 const std::string& test = "") {
  bool ok = std::abs(a - b) <= tol;
  std::ostringstream msg;
  msg << label << ": got " << a << ", expected " << b << " (tol=" << tol
      << ", diff=" << std::abs(a - b) << ")";
  return check(ok, msg.str(), test);
}

static bool in_range(double v, double lo, double hi, const std::string& label = "",
                     const std::string& test = "") {
  bool ok = (v >= lo) && (v <= hi);
  std::ostringstream msg;
  msg << label << ": got " << v << ", expected in [" << lo << ", " << hi << "]";
  return check(ok, msg.str(), test);
}

static bool rel_close(double v, double expected, double rel_tol, const std::string& label = "",
                      const std::string& test = "") {
  double denom = std::max(1e-12, std::abs(expected));
  double rel   = std::abs(v - expected) / denom;
  std::ostringstream msg;
  msg << label << ": got " << v << ", expected " << expected << " (rel_tol=" << rel_tol
      << ", rel=" << rel << ")";
  return check(rel <= rel_tol, msg.str(), test);
}

static bool factor_close(double v, double expected, double max_factor,
                         const std::string& label = "", const std::string& test = "") {
  double lo = expected / max_factor;
  double hi = expected * max_factor;
  return in_range(v, lo, hi, label, test);
}

static riptide::TracePoint make_point(double t, double mu_y, double mu_z, double var_y,
                                      double var_z, double cov_yz = 0.0) {
  riptide::TracePoint p{};
  p.t            = t;
  p.mu_y         = mu_y;
  p.mu_z         = mu_z;
  p.cov          = {var_y, cov_yz, var_z};
  p.valid        = true;
  p.n_hits       = 1000;
  p.n_hits_count = 1000;
  return p;
}

static std::vector<riptide::TracePoint> make_perfect_trace(double a, double b, int N, double var_y,
                                                           double var_z, double cov_yz = 0.0) {
  std::vector<riptide::TracePoint> trace;
  trace.reserve(N);
  for (int i = 0; i < N; ++i) {
    double y = -5.0 + 10.0 * i / (N - 1);
    trace.push_back(make_point(static_cast<double>(i), y, a * y + b, var_y, var_z, cov_yz));
  }
  return trace;
}

static std::vector<riptide::TracePoint> make_trace_with_cov(double a, double b, int N, double cov,
                                                            int n_hits) {
  std::vector<riptide::TracePoint> trace;
  trace.reserve(N);
  for (int i = 0; i < N; ++i) {
    double y = -5.0 + 10.0 * i / (N - 1);
    riptide::TracePoint p{};
    p.t            = static_cast<double>(i);
    p.mu_y         = y;
    p.mu_z         = a * y + b;
    p.cov          = {cov, 0.0, cov};
    p.valid        = true;
    p.n_hits       = static_cast<double>(n_hits);
    p.n_hits_count = static_cast<double>(n_hits);
    trace.push_back(p);
  }
  return trace;
}

static std::vector<riptide::TracePoint>
make_trace_with_nhits(double a, double b, int N, double sigma2, int n_hits, bool scale_by_nhits) {
  double cov = scale_by_nhits ? sigma2 / static_cast<double>(n_hits) : sigma2;
  return make_trace_with_cov(a, b, N, cov, n_hits);
}

static std::vector<riptide::TracePoint> add_noise(std::vector<riptide::TracePoint> trace,
                                                  double sigma_noise, unsigned seed = 42) {
  std::mt19937 gen(seed);
  std::normal_distribution<double> gaus(0.0, sigma_noise);
  for (auto& p : trace) {
    p.mu_y += gaus(gen);
    p.mu_z += gaus(gen);
  }
  return trace;
}

static std::vector<riptide::TracePoint>
add_correlated_noise_along_normal(std::vector<riptide::TracePoint> trace, double a,
                                  double sigma_noise, double rho, unsigned seed = 42) {
  rho = std::clamp(rho, 0.0, 0.999999);
  std::mt19937 gen(seed);
  std::normal_distribution<double> gaus(0.0, sigma_noise);

  double norm = std::sqrt(1.0 + a * a);
  double ny   = -a / norm;
  double nz   = 1.0 / norm;

  double eps_prev = gaus(gen);
  for (size_t i = 0; i < trace.size(); ++i) {
    double xi  = gaus(gen);
    double eps = (i == 0) ? eps_prev : rho * eps_prev + std::sqrt(1.0 - rho * rho) * xi;
    eps_prev   = eps;
    trace[i].mu_y += ny * eps;
    trace[i].mu_z += nz * eps;
  }
  return trace;
}

struct MeanStats {
  double mean = 0.0;
  double rms  = 0.0;
};

static MeanStats mean_rms(const std::vector<double>& v) {
  MeanStats s{};
  if (v.empty())
    return s;
  double sum = std::accumulate(v.begin(), v.end(), 0.0);
  s.mean     = sum / static_cast<double>(v.size());
  double sq  = 0.0;
  for (double x : v)
    sq += (x - s.mean) * (x - s.mean);
  s.rms = std::sqrt(sq / static_cast<double>(v.size()));
  return s;
}

static double expected_chi2_ndof_ar1_linear_fit(int N, double rho) {
  if (N <= 2)
    return 0.0;

  rho = std::clamp(rho, 0.0, 0.999999);

  std::vector<double> y;
  y.reserve(static_cast<size_t>(N));
  for (int i = 0; i < N; ++i)
    y.push_back(-5.0 + 10.0 * static_cast<double>(i) / static_cast<double>(N - 1));

  double S1  = static_cast<double>(N);
  double Sy  = std::accumulate(y.begin(), y.end(), 0.0);
  double Syy = 0.0;
  for (double yi : y)
    Syy += yi * yi;

  double det = Syy * S1 - Sy * Sy;
  if (std::abs(det) < 1e-12)
    return 0.0;

  double inv00 = S1 / det;
  double inv01 = -Sy / det;
  double inv10 = -Sy / det;
  double inv11 = Syy / det;

  double M00 = 0.0, M01 = 0.0, M11 = 0.0;
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      double Rij = std::pow(rho, std::abs(i - j));
      M00 += y[static_cast<size_t>(i)] * Rij * y[static_cast<size_t>(j)];
      M01 += y[static_cast<size_t>(i)] * Rij;
      M11 += Rij;
    }
  }
  double M10 = M01;

  double trace_HR   = inv00 * M00 + inv01 * M10 + inv10 * M01 + inv11 * M11;
  double trace_R    = static_cast<double>(N);
  double trace_IH_R = trace_R - trace_HR;

  return trace_IH_R / static_cast<double>(N - 2);
}

static void test_TD1() {
  std::cout << "\n[TD1] Baseline: cov=errore sulla media, dati perfetti\n";
  const std::string T = "TD1";

  const int N         = 21;
  const double sigma2 = 0.01;
  const int n_hits    = 1000;

  auto trace = make_trace_with_nhits(1.0, 0.0, N, sigma2, n_hits, true);
  auto res   = riptide::fit_trace(trace);

  near(res.chi2, 0.0, 1e-9, "chi2", T);
  near(res.chi2_ndof, 0.0, 1e-9, "chi2/ndof", T);
  check(res.ndof == (N - 2), "ndof = N-2", T);
  check(res.converged, "converged", T);
}

static void test_TD_REG() {
  std::cout << "\n[TD_REG] Regressione: build_trace_3d scala cov con 1/n_hits_count\n";
  const std::string T = "TD_REG";

  const double sigma2 = 0.04;
  const int N_H       = 500;

  riptide::PSFDatabase db;
  riptide::LensConfig cfg{50.0, 120.0};

  std::vector<double> x_vals = {-10.0, 0.0, 10.0};
  for (double x : x_vals) {
    for (double r = 0.0; r <= 10.0; r += 1.0) {
      db[cfg].push_back({x, r, r, 0.0, sigma2, 0.0, sigma2, true, static_cast<double>(N_H),
                         static_cast<double>(N_H)});
    }
  }

  riptide::Point3D p1{-5.0, 5.0, 0.0};
  riptide::Point3D p2{+5.0, 5.0, 0.0};
  auto trace = riptide::build_trace_3d(p1, p2, cfg, db, 1.0);
  check(!trace.empty(), "traccia non vuota", T);

  double expected = sigma2 / static_cast<double>(N_H);
  for (const auto& pt : trace) {
    if (pt.n_hits_count <= 0.0)
      continue;
    check(std::abs(pt.cov.yy - expected) < 0.01 * expected,
          "cov.yy ≈ sigma2/N_H (got " + std::to_string(pt.cov.yy) + ", expected "
              + std::to_string(expected) + ")",
          T);
    check(std::abs(pt.cov.zz - expected) < 0.01 * expected,
          "cov.zz ≈ sigma2/N_H (got " + std::to_string(pt.cov.zz) + ", expected "
              + std::to_string(expected) + ")",
          T);
  }
}

static void test_TD2() {
  std::cout << "\n[TD2] Candidato 1: scaling covarianza vs n_hits (rumore i.i.d.)\n";
  const std::string T = "TD2";

  const int N                = 51;
  const double sigma2        = 0.01;
  const double sigma_noise   = std::sqrt(sigma2);
  const int n_hits           = 100;
  const int n_realizations   = 200;
  const unsigned seed_offset = 1000;

  std::vector<double> chi2_over, chi2_correct, chi2_under;
  chi2_over.reserve(n_realizations);
  chi2_correct.reserve(n_realizations);
  chi2_under.reserve(n_realizations);

  for (int k = 0; k < n_realizations; ++k) {
    unsigned seed = seed_offset + static_cast<unsigned>(k);

    auto tr_over =
        add_noise(make_trace_with_nhits(1.0, 0.0, N, sigma2, n_hits, true), sigma_noise, seed);
    auto tr_correct =
        add_noise(make_trace_with_nhits(1.0, 0.0, N, sigma2, n_hits, false), sigma_noise, seed);
    auto tr_under =
        add_noise(make_trace_with_cov(1.0, 0.0, N, sigma2 * n_hits, n_hits), sigma_noise, seed);

    chi2_over.push_back(riptide::fit_trace(tr_over).chi2_ndof);
    chi2_correct.push_back(riptide::fit_trace(tr_correct).chi2_ndof);
    chi2_under.push_back(riptide::fit_trace(tr_under).chi2_ndof);
  }

  MeanStats s_over    = mean_rms(chi2_over);
  MeanStats s_correct = mean_rms(chi2_correct);
  MeanStats s_under   = mean_rms(chi2_under);

  std::cout << "  mean chi2/ndof: over=" << s_over.mean << "  correct=" << s_correct.mean
            << "  under=" << s_under.mean << "\n";

  in_range(s_correct.mean, 0.5, 2.0, "chi2_correct in [0.5, 2.0]", T);
  check(s_over.mean > 10.0 * s_correct.mean, "chi2_over > 10 * chi2_correct", T);
  check(s_under.mean < 0.1 * s_correct.mean, "chi2_under < 0.1 * chi2_correct", T);
  check(s_over.mean > s_correct.mean && s_correct.mean > s_under.mean,
        "ordine: chi2_over >> chi2_correct >> chi2_under", T);
}

static void test_TD3() {
  std::cout << "\n[TD3] Candidato 2: residui correlati AR(1) con cov marginale corretta\n";
  const std::string T = "TD3";

  const int N                = 51;
  const double sigma2        = 0.01;
  const double sigma_noise   = std::sqrt(sigma2);
  const int n_realizations   = 500;
  const unsigned seed_offset = 20000;

  std::vector<double> rhos = {0.0, 0.5, 0.9};

  for (size_t ic = 0; ic < rhos.size(); ++ic) {
    double rho      = rhos[ic];
    double expected = expected_chi2_ndof_ar1_linear_fit(N, rho);

    std::vector<double> vals;
    vals.reserve(n_realizations);

    for (int k = 0; k < n_realizations; ++k) {
      unsigned seed = seed_offset + static_cast<unsigned>(1000 * ic + k);
      auto tr       = make_trace_with_cov(1.0, 0.0, N, sigma2, 1000);
      tr            = add_correlated_noise_along_normal(std::move(tr), 1.0, sigma_noise, rho, seed);
      vals.push_back(riptide::fit_trace(tr).chi2_ndof);
    }

    MeanStats s = mean_rms(vals);
    std::cout << "  rho=" << rho << "  mean chi2/ndof=" << s.mean << "  expected~" << expected
              << "\n";

    rel_close(s.mean, expected, 0.20, "AR(1) chi2/ndof (atteso con fit lineare)", T);
    if (rho == 0.0)
      in_range(s.mean, 0.5, 2.0, "baseline rho=0: chi2/ndof ~ 1", T);
  }
}

static void test_TD4() {
  std::cout << "\n[TD4] Combinato: cov distribuzione + rumore AR(1) (errore sulla media)\n";
  const std::string T = "TD4";

  const int N                = 51;
  const double sigma2        = 0.01;
  const int n_hits           = 100;
  const double rho           = 0.7;
  const int n_realizations   = 500;
  const unsigned seed_offset = 40000;

  const double sigma_noise_mean = std::sqrt(sigma2 / static_cast<double>(n_hits));
  const double expected = expected_chi2_ndof_ar1_linear_fit(N, rho) / static_cast<double>(n_hits);

  std::vector<double> vals;
  vals.reserve(n_realizations);
  for (int k = 0; k < n_realizations; ++k) {
    unsigned seed = seed_offset + static_cast<unsigned>(k);
    auto tr       = make_trace_with_cov(1.0, 0.0, N, sigma2, n_hits);
    tr = add_correlated_noise_along_normal(std::move(tr), 1.0, sigma_noise_mean, rho, seed);
    vals.push_back(riptide::fit_trace(tr).chi2_ndof);
  }

  MeanStats s = mean_rms(vals);
  std::cout << "  mean chi2/ndof=" << s.mean << "  expected~" << expected << "\n";
  factor_close(s.mean, expected, 3.0, "chi2/ndof entro fattore 3 dalla formula", T);
}

static void test_TD5() {
  std::cout << "\n[TD5] Monotonia su dt: fisica fissa, campionamento più fitto\n";
  const std::string T = "TD5";

  const double L      = 10.0;
  const double sigma2 = 0.01;
  const double sigma  = std::sqrt(sigma2);
  const int n_hits    = 1000;

  const double rho_dt1 = 0.8;
  const double l_corr  = -1.0 / std::log(rho_dt1);

  std::vector<int> Ns      = {11, 21, 41, 81};
  std::vector<double> mean = {};
  mean.reserve(Ns.size());

  const int n_realizations   = 300;
  const unsigned seed_offset = 60000;

  for (size_t i = 0; i < Ns.size(); ++i) {
    int N        = Ns[i];
    double dt    = L / static_cast<double>(N - 1);
    double rho   = std::exp(-dt / l_corr);
    double c_exp = expected_chi2_ndof_ar1_linear_fit(N, rho);

    std::vector<double> vals;
    vals.reserve(n_realizations);
    for (int k = 0; k < n_realizations; ++k) {
      unsigned seed = seed_offset + static_cast<unsigned>(1000 * i + k);
      auto tr       = make_trace_with_cov(1.0, 0.0, N, sigma2, n_hits);
      tr            = add_correlated_noise_along_normal(std::move(tr), 1.0, sigma, rho, seed);
      vals.push_back(riptide::fit_trace(tr).chi2_ndof);
    }

    MeanStats s = mean_rms(vals);
    mean.push_back(s.mean);
    std::cout << "  N=" << N << "  dt=" << dt << "  rho~" << rho << "  mean chi2/ndof=" << s.mean
              << "  expected~" << c_exp << "\n";
  }

  bool mono = true;
  for (size_t i = 0; i + 1 < mean.size(); ++i) {
    if (!(mean[i] > mean[i + 1] + 0.005)) {
      mono = false;
      break;
    }
  }
  check(mono, "chi2/ndof strettamente decrescente con N", T);
}

int main() {
  std::cout << "test_chi2_diagnostics — riptide chi2 diagnostic tests\n\n";

  test_TD_REG();
  test_TD1();
  test_TD2();
  test_TD3();
  test_TD4();
  test_TD5();

  std::cout << "\nTEST SUMMARY: " << g_n_pass << " PASS, " << g_n_fail << " FAIL\n";
  return (g_n_fail == 0) ? 0 : 1;
}
