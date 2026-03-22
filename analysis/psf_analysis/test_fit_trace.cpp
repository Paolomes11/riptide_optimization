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
 * INVARIANTI DI INIZIALIZZAZIONE
 *
 * TracePoint::valid = true
 *   fit_trace() filtra i punti con valid==false (fuori dal fotocatodo).
 *   La zero-inizializzazione di struct lascia valid=false → fit_trace()
 *   non trova punti utilizzabili e lancia std::invalid_argument.
 *   make_point() imposta valid=true esplicitamente.
 *
 * PSFPoint::on_detector = true
 *   build_trace() copia on_detector in TracePoint::valid.
 *   Se false → is_trace_valid() ritorna false → compute_Q() accumula
 *   n_invalid invece di n_traces → Q rimane 0 per motivo sbagliato.
 *   make_synthetic_db() imposta on_detector=true esplicitamente.
 *
 * NOTA SUL TEST TQ5 (temporal unfolding)
 *
 * L'unfolding z̃_i = μ_z_i + i·δz modifica la coordinata z di ogni punto.
 * Con μ_z = 0 per tutti i punti, z̃_i = i·δz è identico per qualsiasi DB.
 * Il chi2 del fit ODR su (μ_y_i, z̃_i) dipende SOLO da come μ_y_i si
 * distribuisce rispetto alla retta ottimale.
 *
 * Impossibilità di usare PSFDatabase + build_trace() per questo test:
 *   r(t) = sqrt(y0^2 + t^2) è simmetrica in t per qualsiasi y0>0, quindi
 *   μ_y_i = μ_y_{N-1-i} per qualsiasi DB con μ_y(r) monotona. Non è
 *   possibile costruire una traccia "lineare" (μ_y_i ∝ i) tramite
 *   build_trace() senza accesso a r(t) lineare in i, che richiederebbe y0=0
 *   con t ∈ [0, L] — ma build_trace usa t ∈ [-L/2, +L/2] → r = |t| → ripiegata.
 *
 * Soluzione: costruiamo TracePoint direttamente con μ_z già srotolato:
 *
 *   Traccia LINEARE:   μ_y_i = i·k, z̃_i = i·δz
 *     → punti (i·k, i·δz) su retta z̃ = (δz/k)·μ_y → chi2 ≈ 0
 *
 *   Traccia RIPIEGATA: μ_y_i = |i-N/2|·k, z̃_i = i·δz
 *     → punti (|i-N/2|·k, i·δz) NON collineari → chi2 >> 0
 *
 * Con unfolding OFF (μ_z = 0 per tutti): entrambe le tracce hanno z=0
 * → retta ottimale z=0 → chi2 = 0 per entrambe → delta_chi2 = 0.
 *
 * SUITE DI TEST
 *
 *   T1  Retta z=y, cov isotropa              → a=1, b=0, chi2=0, pull≈0
 *   T2  Retta z=0.5y+3, cov isotropa         → a=0.5, b=3, chi2=0
 *   T3  Retta z=2y-1, cov anisotropa          → parametri invariati
 *   T4  Outlier singolo Δz=0.5mm              → pull_outlier≈Δz/σ
 *   T5  < 3 punti validi                      → std::invalid_argument
 *   T6  Covarianza degenere (var=0)            → convergenza con floor
 *   T7  Cov non diagonale, outlier             → sigma_d=sqrt(cov_zz)
 *   TV1 is_trace_valid, tutti valid            → true
 *   TV2 is_trace_valid, 60% valid              → dipende dalla soglia
 *   TV3 is_trace_valid, traccia vuota          → false
 *   TQ1 PSF ideale, unfolding OFF              → Q≈0
 *   TQ2 y0_values esplicito, unfolding OFF     → n_traces=3, Q≈0
 *   TQ3 Config non presente                    → std::invalid_argument
 *   TQ4 PSF curva, due σ_z diverse             → Q(piccola)>Q(grande)
 *   TQ5 Unfolding ON: lineare vs ripiegata     → chi2(fold)>>chi2(lin),
 *                                                 delta_ON >> delta_OFF
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
  riptide::TracePoint pt{};
  pt.t     = t;
  pt.r     = std::abs(mu_y);
  pt.mu_y  = mu_y;
  pt.mu_z  = mu_z;
  pt.cov   = {var_y, cov_yz, var_z};
  pt.valid = true; // essenziale
  return pt;
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

// Database PSF sintetico con PSF ideale (mu_y=r, mu_z=0, on_detector=true).
// NOTA: on_detector=true è essenziale — build_trace() lo copia in TracePoint::valid.
static riptide::PSFDatabase make_synthetic_db(double x1 = 50.0, double x2 = 120.0,
                                              double r_max = 16.0, double dr = 0.1,
                                              double cov_val = 1e-6) {
  riptide::PSFDatabase db;
  riptide::LensConfig cfg{x1, x2};
  auto& pts = db[cfg];
  for (double r = 0.0; r <= r_max + dr * 1e-9; r += dr) {
    riptide::PSFPoint p;
    p.y_source    = r;
    p.mu_y        = r;
    p.mu_z        = 0.0;
    p.cov_yy      = cov_val;
    p.cov_yz      = 0.0;
    p.cov_zz      = cov_val;
    p.on_detector = true; // essenziale
    pts.push_back(p);
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
    check(riptide::is_trace_valid(trace, 0.50), "soglia 50% → valid", T);
    check(!riptide::is_trace_valid(trace, 0.75), "soglia 75% → invalid", T);
    check(!riptide::is_trace_valid(trace, 1.00), "soglia 100% → invalid", T);
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

  auto db = make_synthetic_db();
  riptide::LensConfig cfg{50.0, 120.0};

  riptide::QConfig qcfg;
  qcfg.y0_min                   = 0.0;
  qcfg.y0_max                   = 10.0;
  qcfg.dy0                      = 1.0;
  qcfg.trace_L                  = 10.0;
  qcfg.trace_dt                 = 0.5;
  qcfg.point_valid_fraction     = 0.50;
  qcfg.trace_valid_fraction     = 0.50;
  qcfg.apply_temporal_unfolding = false;

  riptide::QResult res;
  bool threw = false;
  try {
    res = riptide::compute_Q(cfg, db, qcfg);
  } catch (const std::exception& e) {
    std::cerr << "  ECCEZIONE: " << e.what() << "\n";
    threw = true;
  }

  check(!threw, "nessuna eccezione", T);
  if (threw)
    return;

  check(res.n_traces == 11, "n_traces = 11", T);
  check(res.n_failed == 0, "n_failed = 0", T);
  check(res.config_valid, "config_valid = true", T);
  near(res.Q, 0.0, 1e-4, "Q ≈ 0", T);
  check(res.y0_used.size() == 11, "y0_used.size() == 11", T);
  check(res.chi2_per_y0.size() == 11, "chi2_per_y0.size() == 11", T);

  double Q_sum = 0.0;
  for (double c : res.chi2_per_y0)
    Q_sum += c;
  near(res.Q, Q_sum, 1e-12, "Q == sum(chi2_per_y0)", T);
}

// TQ2
static void test_TQ2() {
  std::cout << "\n[TQ2] compute_Q con y0_values esplicito\n";
  const std::string T = "TQ2";

  auto db = make_synthetic_db();
  riptide::LensConfig cfg{50.0, 120.0};

  riptide::QConfig qcfg;
  qcfg.y0_values                = {2.0, 5.0, 8.0};
  qcfg.trace_L                  = 10.0;
  qcfg.trace_dt                 = 0.5;
  qcfg.point_valid_fraction     = 0.50;
  qcfg.trace_valid_fraction     = 0.50;
  qcfg.apply_temporal_unfolding = false;

  auto res = riptide::compute_Q(cfg, db, qcfg);
  check(res.n_traces == 3, "n_traces == 3", T);
  check(res.y0_used.size() == 3, "y0_used.size() == 3", T);
  near(res.Q, 0.0, 1e-4, "Q ≈ 0", T);
}

// TQ3
static void test_TQ3() {
  std::cout << "\n[TQ3] compute_Q: config non presente → eccezione\n";
  const std::string T = "TQ3";

  auto db    = make_synthetic_db(50.0, 120.0);
  bool threw = false;
  try {
    riptide::compute_Q({99.0, 199.0}, db);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  check(threw, "std::invalid_argument per cfg non presente", T);
}

// TQ4
static void test_TQ4() {
  std::cout << "\n[TQ4] Q monotona: PSF curva, sigma_z piccola → Q maggiore\n";
  const std::string T = "TQ4";

  auto make_curved = [](double cov_zz_val) {
    riptide::PSFDatabase db;
    riptide::LensConfig cfg{50.0, 120.0};
    auto& pts = db[cfg];
    for (double r = 0.0; r <= 16.0 + 1e-9; r += 0.1) {
      riptide::PSFPoint p;
      p.y_source    = r;
      p.mu_y        = r;
      p.mu_z        = 0.05 * r * r;
      p.cov_yy      = 1e-6;
      p.cov_yz      = 0.0;
      p.cov_zz      = cov_zz_val;
      p.on_detector = true;
      pts.push_back(p);
    }
    return db;
  };

  riptide::LensConfig cfg{50.0, 120.0};
  riptide::QConfig qcfg;
  qcfg.y0_min                   = 2.0;
  qcfg.y0_max                   = 8.0;
  qcfg.dy0                      = 2.0;
  qcfg.trace_L                  = 10.0;
  qcfg.trace_dt                 = 0.5;
  qcfg.point_valid_fraction     = 0.50;
  qcfg.trace_valid_fraction     = 0.50;
  qcfg.apply_temporal_unfolding = false;

  auto rt = riptide::compute_Q(cfg, make_curved(1e-4), qcfg);
  auto rl = riptide::compute_Q(cfg, make_curved(1.0), qcfg);

  check(rt.Q > rl.Q, "Q(σ_pic)=" + std::to_string(rt.Q) + " > Q(σ_gr)=" + std::to_string(rl.Q), T);
  check(rt.n_failed == 0, "n_failed==0 (tight)", T);
  check(rl.n_failed == 0, "n_failed==0 (loose)", T);
}

// TQ5: discriminatore del temporal unfolding
//
// Perché costruiamo TracePoint direttamente invece di usare PSFDatabase:
//   r(t) = sqrt(y0^2+t^2) è simmetrica per qualsiasi y0>0, quindi
//   mu_y_i = mu_y_{N-1-i} tramite build_trace(). Non è possibile ottenere
//   una traccia con mu_y_i ∝ i (lineare in i) passando da PSFDatabase,
//   perché r non è lineare in t.
//
//   Costruendo i TracePoint a mano con mu_z già srotolato (mu_z = i*dz),
//   simuliamo esattamente l'output che compute_Q() passerebbe a fit_trace()
//   dopo aver applicato apply_unfolding().
//
//   Traccia LINEARE:   mu_y_i = i*k,        mu_z_i = i*dz
//     → punti (i*k, i*dz) su retta z = (dz/k)*y → chi2 ≈ 0
//
//   Traccia RIPIEGATA: mu_y_i = |i-N/2|*k,  mu_z_i = i*dz
//     → punti non collineari → chi2 >> 0
//
//   Con unfolding OFF: mu_z = 0 per tutti → chi2 = 0 per entrambe.
//
static void test_TQ5() {
  std::cout << "\n[TQ5] Temporal unfolding: traccia ripiegata penalizzata vs lineare\n";
  const std::string T = "TQ5";

  const int N      = 21;
  const double L   = 10.0;
  const double dz  = L / (N - 1); // δz automatico = 0.5 mm/passo
  const double k   = 0.5;         // ampiezza mu_y [mm/passo]
  const double var = 0.01;        // sigma = 0.1 mm

  // Costruisce N TracePoint con mu_y e mu_z specificati, tutti valid=true.
  auto make_trace = [&](const std::vector<double>& mu_y_v, const std::vector<double>& mu_z_v) {
    std::vector<riptide::TracePoint> tr;
    tr.reserve(N);
    for (int i = 0; i < N; ++i) {
      riptide::TracePoint pt{};
      pt.t     = -L / 2.0 + i * (L / (N - 1));
      pt.r     = std::abs(mu_y_v[i]);
      pt.mu_y  = mu_y_v[i];
      pt.mu_z  = mu_z_v[i];
      pt.cov   = {var, 0.0, var};
      pt.valid = true;
      tr.push_back(pt);
    }
    return tr;
  };

  std::vector<double> z_unf(N), z_zero(N, 0.0);
  std::vector<double> mu_y_lin(N), mu_y_fold(N);
  for (int i = 0; i < N; ++i) {
    z_unf[i]     = i * dz;
    mu_y_lin[i]  = i * k;
    mu_y_fold[i] = std::abs(i - N / 2) * k;
  }

  // Con unfolding ON (mu_z già modificato con i*dz)
  auto res_lin_on  = riptide::fit_trace(make_trace(mu_y_lin, z_unf));
  auto res_fold_on = riptide::fit_trace(make_trace(mu_y_fold, z_unf));

  std::cout << "  [unf ON]  chi2_lin=" << res_lin_on.chi2 << "  chi2_fold=" << res_fold_on.chi2
            << "\n";

  check(res_fold_on.chi2 > res_lin_on.chi2 + 1.0,
        "unf ON: chi2(fold)=" + std::to_string(res_fold_on.chi2)
            + " >> chi2(lin)=" + std::to_string(res_lin_on.chi2),
        T);

  // Con unfolding OFF (mu_z = 0 per tutti → chi2 = 0 per entrambe)
  auto res_lin_off  = riptide::fit_trace(make_trace(mu_y_lin, z_zero));
  auto res_fold_off = riptide::fit_trace(make_trace(mu_y_fold, z_zero));

  std::cout << "  [unf OFF] chi2_lin=" << res_lin_off.chi2 << "  chi2_fold=" << res_fold_off.chi2
            << "\n";

  near(res_lin_off.chi2, 0.0, 1e-9, "unf OFF: chi2_lin ≈ 0", T);
  near(res_fold_off.chi2, 0.0, 1e-9, "unf OFF: chi2_fold ≈ 0", T);

  // La differenza con unfolding ON deve essere molto maggiore che con OFF
  double delta_on  = res_fold_on.chi2 - res_lin_on.chi2;
  double delta_off = res_fold_off.chi2 - res_lin_off.chi2;
  check(delta_on > delta_off + 1.0,
        "delta_chi2(unf ON)=" + std::to_string(delta_on)
            + " >> delta_chi2(unf OFF)=" + std::to_string(delta_off),
        T);
}

// main
int main() {
  std::cout << "  test_fit_trace — riptide PSF unit tests\n";
  std::cout << "  (fit_trace, is_trace_valid, compute_Q, temporal unfolding)\n\n";

  test_T1();
  test_T2();
  test_T3();
  test_T4();
  test_T5();
  test_T6();
  test_T7();
  test_TV();
  test_TQ1();
  test_TQ2();
  test_TQ3();
  test_TQ4();
  test_TQ5();

  std::cout << "\n  Risultato: " << g_n_pass << " PASS, " << g_n_fail << " FAIL\n";
  return (g_n_fail == 0) ? 0 : 1;
}