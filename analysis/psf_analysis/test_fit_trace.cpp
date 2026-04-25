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

/*
 * test_fit_trace — Test unitari per fit_trace(), is_trace_valid() e compute_Q().
 *
 * Tutti i casi hanno soluzione analitica nota: non serve nessun file ROOT.
 * I dati vengono costruiti sinteticamente come vettori di TracePoint o
 * PSFDatabase, bypassando load_psf_database() e build_trace().
 *
 * I test usano una PSF analitica molto semplice:
 * mu_y(r) = r
 * mu_z(r) = 0
 * cov     = isotropa (var)
 *
 * In questo modo le tracce sul detector sono semplicemente
 * mu_y = sqrt(y0^2 + t^2)
 * mu_z = 0
 *
 * Questo permette di calcolare i risultati attesi di fit_trace.
 *
 * I test su compute_Q verificano il comportamento ad alto livello.
 */

#include "psf_interpolator.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Helpers
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

// Costruisce un TracePoint con valid=true.
// NOTA: valid=true è essenziale — fit_trace() scarta i punti con valid=false.
static riptide::TracePoint make_point(double t, double mu_y, double mu_z, double var_y,
                                      double var_z, double cov_yz = 0.0) {
  riptide::TracePoint p{};
  p.t            = t;
  p.mu_y         = mu_y;
  p.mu_z         = mu_z;
  p.cov          = {var_y, cov_yz, var_z};
  p.valid        = true;
  p.n_hits       = 1000; // Valido per default
  p.n_hits_count = 1000;
  return p;
}

// Costruisce N punti sulla retta z = a*y + b, tutti con valid=true.
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

static riptide::PSFDatabase make_synthetic_db(double x1 = 50.0, double x2 = 120.0, double sy = 0.1,
                                              double sz = 0.1) {
  riptide::PSFDatabase db;
  riptide::LensConfig cfg{x1, x2};
  // Griglia 1D (x=0, y in [0, 10]) per compatibilità con i vecchi test
  for (double y = 0; y <= 10.0; y += 1.0) {
    db[cfg].push_back({0.0, y, y, 0.0, sy * sy, 0.0, sz * sz, true, 1000, 1000});
  }
  return db;
}

// T1
static void test_T1() {
  std::cout << "\n[T1] Retta perfetta z = y  (a=1, b=0), cov isotropa\n";
  const std::string T = "T1";

  auto trace = make_perfect_trace(1.0, 0.0, 21, 0.01, 0.01);
  for (const auto& pt : trace)
    check(pt.valid, "punto valid=true", T);

  auto res = riptide::fit_trace(trace);
  near(res.a, 1.0, 1e-9, "a", T);
  near(res.b, 0.0, 1e-9, "b", T);
  near(res.chi2, 0.0, 1e-9, "chi2", T);
  near(res.chi2_ndof, 0.0, 1e-9, "chi2/ndof", T);
  check(res.converged, "converged", T);
  check(res.ndof == 19, "ndof = 19", T);
  check(res.n_points_used == 21, "n_points_used = 21", T);

  bool pulls_ok = true;
  for (size_t i = 0; i < res.pull.size(); ++i)
    if (std::abs(res.pull[i]) > 1e-9) {
      pulls_ok = false;
      break;
    }
  check(pulls_ok, "tutti i pull ≈ 0", T);
}

// T2
static void test_T2() {
  std::cout << "\n[T2] Retta perfetta z = 0.5*y + 3, cov isotropa\n";
  const std::string T = "T2";

  auto trace = make_perfect_trace(0.5, 3.0, 31, 0.04, 0.04);
  auto res   = riptide::fit_trace(trace);

  near(res.a, 0.5, 1e-9, "a", T);
  near(res.b, 3.0, 1e-9, "b", T);
  near(res.chi2, 0.0, 1e-9, "chi2", T);
  check(res.converged, "converged", T);
  check(res.ndof == 29, "ndof = 29", T);
  check(res.n_points_used == 31, "n_points_used = 31", T);
  check(res.sigma_a > 0.0, "sigma_a > 0", T);
  check(res.sigma_b > 0.0, "sigma_b > 0", T);
}

// T3
static void test_T3() {
  std::cout << "\n[T3] Retta perfetta z = 2*y - 1, cov anisotropa\n";
  const std::string T = "T3";

  auto trace = make_perfect_trace(2.0, -1.0, 21, 1.0, 1e-4);
  auto res   = riptide::fit_trace(trace);

  near(res.a, 2.0, 1e-7, "a", T);
  near(res.b, -1.0, 1e-7, "b", T);
  near(res.chi2, 0.0, 1e-9, "chi2", T);
  check(res.converged, "converged", T);
}

// T4
static void test_T4() {
  std::cout << "\n[T4] Outlier singolo Δz = 0.5 mm\n";
  const std::string T = "T4";

  const int N        = 21;
  const double sigma = 0.1, var = sigma * sigma, delta = 0.5;
  auto trace    = make_perfect_trace(0.0, 0.0, N, var, var);
  const int idx = N / 2;
  trace[idx].mu_z += delta;

  auto res = riptide::fit_trace(trace);

  double ep = delta / sigma;
  check(res.pull.size() == static_cast<size_t>(N), "size pull == N", T);
  check(res.n_points_used == N, "n_points_used == N", T);
  if (!res.pull.empty()) {
    double p = std::abs(res.pull[idx]);
    check(p > 0.8 * ep && p <= ep + 1.0,
          "pull outlier in [" + std::to_string(0.8 * ep) + ", " + std::to_string(ep + 1.0)
              + "], got " + std::to_string(p),
          T);
  }
  bool ok = true;
  for (int i = 0; i < N; ++i)
    if (i != idx && std::abs(res.pull[i]) > 2.0) {
      ok = false;
      break;
    }
  check(ok, "pull non-outlier < 2", T);
  check(res.chi2_ndof > 0.1 * delta * delta / var / (N - 2), "chi2/ndof dominato dall'outlier", T);
}

// T5
static void test_T5() {
  std::cout << "\n[T5] Eccezione su traccia con < 3 punti validi\n";
  const std::string T = "T5";

  {
    bool t = false;
    try {
      riptide::fit_trace({});
    } catch (const std::invalid_argument&) {
      t = true;
    }
    check(t, "traccia vuota → invalid_argument", T);
  }

  {
    auto tr = make_perfect_trace(1.0, 0.0, 2, 0.01, 0.01);
    bool t  = false;
    try {
      riptide::fit_trace(tr);
    } catch (const std::invalid_argument&) {
      t = true;
    }
    check(t, "2 punti validi → invalid_argument", T);
  }

  {
    auto tr = make_perfect_trace(1.0, 0.0, 3, 0.01, 0.01);
    for (auto& pt : tr)
      pt.valid = false;
    bool t = false;
    try {
      riptide::fit_trace(tr);
    } catch (const std::invalid_argument&) {
      t = true;
    }
    check(t, "3 punti tutti invalid → invalid_argument", T);
  }

  {
    auto tr = make_perfect_trace(1.0, 0.0, 3, 0.01, 0.01);
    bool t  = false;
    try {
      riptide::fit_trace(tr);
    } catch (...) {
      t = true;
    }
    check(!t, "3 punti validi → nessuna eccezione", T);
  }

  {
    auto tr = make_perfect_trace(1.0, 0.0, 5, 0.01, 0.01);
    for (int i = 2; i < 5; ++i)
      tr[i].valid = false;
    bool t = false;
    try {
      riptide::fit_trace(tr);
    } catch (const std::invalid_argument&) {
      t = true;
    }
    check(t, "5 punti con 2 valid → invalid_argument", T);
  }
}

// T6
static void test_T6() {
  std::cout << "\n[T6] Covarianza degenere (var=0)\n";
  const std::string T = "T6";

  auto trace = make_perfect_trace(0.3, 1.5, 11, 0.0, 0.0);
  bool threw = false;
  riptide::LineFitResult res{};
  try {
    res = riptide::fit_trace(trace);
  } catch (...) {
    threw = true;
  }

  check(!threw, "nessuna eccezione con var=0", T);
  check(res.converged, "converged", T);
  near(res.a, 0.3, 1e-7, "a", T);
  near(res.b, 1.5, 1e-7, "b", T);
}

// T7
static void test_T7() {
  std::cout << "\n[T7] Verifica formula sigma_d con cov non diagonale\n";
  const std::string T = "T7";

  const int N = 11;
  std::vector<riptide::TracePoint> trace;
  trace.reserve(N);
  for (int i = 0; i < N; ++i) {
    double y = -5.0 + 10.0 * i / (N - 1);
    trace.push_back(make_point(i, y, 0.0, 4.0, 0.09, 0.5));
  }
  const int idx   = N / 2;
  trace[idx].mu_z = 0.6;

  auto res   = riptide::fit_trace(trace);
  double esd = std::sqrt(0.09), ep = 0.6 / esd;

  check(std::abs(res.residual_sig[idx] - esd) < 0.05 * esd,
        "sigma_d ≈ sqrt(cov_zz)=" + std::to_string(esd) + ", got "
            + std::to_string(res.residual_sig[idx]),
        T);
  check(std::abs(res.pull[idx]) > 0.7 * ep,
        "pull > 70% di " + std::to_string(ep) + ", got " + std::to_string(std::abs(res.pull[idx])),
        T);
}

// T8: selezione automatica asse — traccia quasi verticale
// Con Δz >> Δy, il fit deve scegliere FitAxis::YvsZ e produrre chi2 ≈ 0.
static void test_T8() {
  std::cout << "\n[T8] Selezione automatica asse: traccia quasi verticale\n";
  const std::string T = "T8";

  const int N = 21;
  std::vector<riptide::TracePoint> trace;
  trace.reserve(N);
  for (int i = 0; i < N; ++i) {
    double z = -5.0 + 10.0 * i / (N - 1);
    double y = 0.02 * z;
    trace.push_back(make_point(static_cast<double>(i), y, z, 0.0001, 0.01));
  }

  riptide::LineFitResult res = riptide::fit_trace(trace);

  check(res.axis == riptide::FitAxis::YvsZ, "axis == YvsZ per traccia verticale", T);
  check(res.converged, "converged", T);
  near(res.chi2, 0.0, 1e-6, "chi2 ≈ 0 (punti collineari)", T);
  near(res.chi2_ndof, 0.0, 1e-6, "chi2/ndof ≈ 0", T);
  near(res.a, 0.02, 1e-6, "a ≈ 0.02 (pendenza y vs z)", T);
  near(res.b, 0.0, 1e-6, "b ≈ 0 (intercetta)", T);
}

// Suite TV
static void test_TV() {
  std::cout << "\n[TV1] is_trace_valid: tutti i punti validi\n";
  {
    const std::string T = "TV1";
    auto trace          = make_perfect_trace(1.0, 0.0, 10, 0.01, 0.01);
    check(riptide::is_trace_valid(trace, 0.75), "soglia 75% → valid", T);
    check(riptide::is_trace_valid(trace, 1.00), "soglia 100% → valid", T);
  }

  std::cout << "\n[TV2] is_trace_valid: 60% validi\n";
  {
    const std::string T = "TV2";
    auto trace          = make_perfect_trace(1.0, 0.0, 10, 0.01, 0.01);
    for (int i = 0; i < 4; ++i)
      trace[i].valid = false;
    // 6/10 = 60%, soglia richiesta 75%
    check(!riptide::is_trace_valid(trace, 0.75), "soglia 75% → invalid", T);
    // 6/10 = 60%, soglia richiesta 50%
    check(riptide::is_trace_valid(trace, 0.50), "soglia 50% → valid", T);
  }

  std::cout << "\n[TV3] is_trace_valid: traccia vuota\n";
  {
    const std::string T = "TV3";
    check(!riptide::is_trace_valid({}, 0.75), "vuota → invalid", T);
  }
}

// TQ1
static void test_TQ1() {
  std::cout << "\n[TQ1] compute_Q con PSF ideale (Q ≈ 0)\n";
  const std::string T = "TQ1";

  auto db = make_synthetic_db(50.0, 120.0, 0.1, 0.1);
  riptide::LensConfig cfg{50.0, 120.0};

  riptide::QConfig qcfg;
  qcfg.n_tracks = 10;
  qcfg.scint_x  = 10.0;
  qcfg.scint_y  = 10.0;
  qcfg.scint_z  = 0.0;
  qcfg.trace_dt = 0.5;

  riptide::QResult res;
  bool threw = false;
  try {
    res = riptide::compute_Q(cfg, db, qcfg, true);
  } catch (const std::exception& e) {
    std::cerr << "  ECCEZIONE: " << e.what() << "\n";
    threw = true;
  }

  check(!threw, "nessuna eccezione", T);
  if (threw)
    return;

  check(res.n_traces > 0, "almeno una traccia valida", T);
  check(res.n_failed == 0, "n_failed = 0", T);
  near(res.Q, 0.0, 1e-4, "Q ≈ 0", T);
}

// TQ2
static void test_TQ2() {
  std::cout << "\n[TQ2] compute_Q con scint personalizzato\n";
  const std::string T = "TQ2";

  auto db = make_synthetic_db(50.0, 120.0, 0.1, 0.1);
  riptide::LensConfig cfg{50.0, 120.0};

  riptide::QConfig qcfg;
  qcfg.scint_x  = 20.0;
  qcfg.scint_y  = 10.0;
  qcfg.scint_z  = 0.0;
  qcfg.n_tracks = 50;
  qcfg.trace_dt = 0.5;

  auto res = riptide::compute_Q(cfg, db, qcfg, true);
  check(res.n_traces > 0, "n_traces > 0", T);
  near(res.Q, 0.0, 1e-4, "Q ≈ 0", T);
}

// TQ3
static void test_TQ3() {
  std::cout << "\n[TQ3] compute_Q: config non presente → eccezione\n";
  const std::string T = "TQ3";

  auto db    = make_synthetic_db(50.0, 120.0, 0.1, 0.1);
  bool threw = false;
  try {
    riptide::compute_Q({99.0, 199.0}, db);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  check(threw, "std::invalid_argument per cfg non presente", T);
}

// TI2D — verifica interpolazione bilineare (x e y non-banali)
static void test_TI2D() {
  std::cout << "\n[TI2D] interpolate() bilineare su griglia 2D\n";
  const std::string T = "TI2D";

  // DB con due valori di x (0 e 10) e y in [0, 5]
  // mu_y(x, y) = x + y  → facile da verificare analiticamente
  riptide::PSFDatabase db;
  riptide::LensConfig cfg{50.0, 120.0};
  for (double xs : {0.0, 10.0}) {
    for (double ys = 0.0; ys <= 5.0; ys += 1.0) {
      db[cfg].push_back(
          {xs, ys, xs + ys, 0.0, 0.01, 0.0, 0.01, true, 1000.0, 1000.0});
    }
  }
  // I punti sono già ordinati per x poi y (x=0 prima, x=10 dopo).

  // Punto interno: x=5, y=2.5 → mu_y atteso = 5+2.5 = 7.5
  auto v = riptide::interpolate(5.0, 2.5, cfg, db);
  near(v.mu_y, 7.5, 1e-9, "mu_y bilineare (x=5, y=2.5)", T);
  check(v.on_detector, "on_detector=true per punto interno", T);

  // Punto al bordo x=0: deve restituire esattamente il valore della griglia
  auto v0 = riptide::interpolate(0.0, 3.0, cfg, db);
  near(v0.mu_y, 3.0, 1e-9, "mu_y bordo x=0, y=3", T);

  // Punto al bordo x=10
  auto v1 = riptide::interpolate(10.0, 2.0, cfg, db);
  near(v1.mu_y, 12.0, 1e-9, "mu_y bordo x=10, y=2", T);
}

// TQ4
static void test_TQ4() {
  std::cout << "\n[TQ4] Q monotona: PSF curva, sigma_z piccola → Q maggiore\n";
  const std::string T = "TQ4";

  auto make_curved = [](double cov_zz_val) {
    riptide::PSFDatabase db;
    riptide::LensConfig cfg{50.0, 120.0};
    for (double x = -30.0; x <= 30.0; x += 10.0) {
      for (double y = 0.0; y <= 10.0; y += 1.0) {
        db[cfg].push_back({x, y, y, 0.05 * y * y, 1e-6, 0.0, cov_zz_val, true, 1000, 1000});
      }
    }
    return db;
  };

  riptide::LensConfig cfg{50.0, 120.0};
  riptide::QConfig qcfg;
  qcfg.n_tracks = 50;
  qcfg.scint_x  = 60.0;
  qcfg.scint_y  = 20.0;
  qcfg.scint_z  = 20.0;
  qcfg.trace_dt = 1.0;

  auto rt = riptide::compute_Q(cfg, make_curved(1e-4), qcfg);
  auto rb = riptide::compute_Q(cfg, make_curved(1e-2), qcfg);

  check(rt.Q > rb.Q, "Q(sigma=0.01) > Q(sigma=0.1)", T);
}

// main
int main() {
  std::cout << "  test_fit_trace — riptide PSF unit tests\n";
  std::cout << "  (fit_trace, is_trace_valid, compute_Q)\n\n";

  test_T1();
  test_T2();
  test_T3();
  test_T4();
  test_T5();
  test_T6();
  test_T7();
  test_T8();
  test_TV();
  test_TI2D();
  test_TQ1();
  test_TQ2();
  test_TQ3();
  test_TQ4();

  std::cout << "\n--------------------------------------\n";
  std::cout << "TEST SUMMARY: " << g_n_pass << " PASS, " << g_n_fail << " FAIL\n";
  return (g_n_fail == 0) ? 0 : 1;
}
