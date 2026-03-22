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
 * test_fit_trace — Test unitari per fit_trace() in psf_interpolator.
 *
 * Tutti i casi hanno soluzione analitica nota: non serve nessun file ROOT.
 * La traccia viene costruita sinteticamente direttamente come vettore di
 * TracePoint, bypassando load_psf_database() e build_trace().
 *
 * Casi testati:
 *   [T1]  Retta perfetta y=x,  cov isotropa uniforme
 *         → a=1, b=0, chi2/ndof≈0, pull tutti ≈ 0
 *
 *   [T2]  Retta z = 0.5*y + 3, cov isotropa uniforme
 *         → a=0.5, b=3, chi2/ndof≈0
 *
 *   [T3]  Retta perfetta z = 2*y - 1, cov anisotropa (sigma_y >> sigma_z)
 *         → parametri invariati rispetto al caso isotropo: la soluzione ODR
 *           deve coincidere con OLS quando i punti sono esattamente sulla retta
 *           (residui nulli indipendentemente dai pesi)
 *
 *   [T4]  Retta z = 0, punti esattamente su z=0 tranne uno spostato di delta
 *         → pull dell'outlier ≈ delta/sigma_d, pull degli altri ≈ 0
 *         → chi2/ndof ≈ (delta/sigma_d)^2 / (N-2)
 *
 *   [T5]  Traccia con meno di 3 punti → deve lanciare std::invalid_argument
 *
 *   [T6]  Covarianza degenere (sigma=0 su tutti i punti) → deve convergere
 *         grazie al floor interno sd2_floor, senza crash
 */

#include "psf_interpolator.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

//  Helpers
static int g_n_pass = 0;
static int g_n_fail = 0;

// Stampa il risultato di un singolo check e aggiorna i contatori globali.
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

// Confronta due double con tolleranza assoluta.
static bool near(double a, double b, double tol = 1e-6, const std::string& label = "",
                 const std::string& test = "") {
  bool ok = std::abs(a - b) <= tol;
  std::ostringstream msg;
  msg << label << ": got " << a << ", expected " << b << " (tol=" << tol
      << ", diff=" << std::abs(a - b) << ")";
  return check(ok, msg.str(), test);
}

// Costruisce una TracePoint con covarianza diagonale (cov_yz = 0).
static riptide::TracePoint make_point(double t, double mu_y, double mu_z, double var_y,
                                      double var_z, double cov_yz = 0.0) {
  riptide::TracePoint pt{};
  pt.t    = t;
  pt.r    = std::abs(mu_y); // non usato dal fit, ma riempito per completezza
  pt.mu_y = mu_y;
  pt.mu_z = mu_z;
  pt.cov  = {var_y, cov_yz, var_z};
  return pt;
}

// Costruisce N punti esattamente sulla retta z = a*y + b, con cov uniforme.
static std::vector<riptide::TracePoint> make_perfect_trace(double a, double b, int N, double var_y,
                                                           double var_z, double cov_yz = 0.0) {
  std::vector<riptide::TracePoint> trace;
  trace.reserve(N);
  for (int i = 0; i < N; ++i) {
    double y = -5.0 + 10.0 * i / (N - 1); // y in [-5, 5] mm
    double z = a * y + b;
    trace.push_back(make_point(static_cast<double>(i), y, z, var_y, var_z, cov_yz));
  }
  return trace;
}

//  Test T1: retta y=x (a=1, b=0), cov isotropa sigma=0.1 mm
static void test_T1() {
  std::cout << "\n[T1] Retta perfetta z = y  (a=1, b=0), cov isotropa\n";
  const std::string T = "T1";

  const double a_true = 1.0, b_true = 0.0;
  const double var = 0.01; // sigma = 0.1 mm
  auto trace       = make_perfect_trace(a_true, b_true, 21, var, var, 0.0);

  auto res = riptide::fit_trace(trace);

  near(res.a, a_true, 1e-9, "a", T);
  near(res.b, b_true, 1e-9, "b", T);
  near(res.chi2, 0.0, 1e-9, "chi2", T);
  near(res.chi2_ndof, 0.0, 1e-9, "chi2/ndof", T);
  check(res.converged, "converged", T);
  check(res.ndof == 19, "ndof = N-2 = 19", T);

  // Tutti i pull devono essere zero (punti esattamente sulla retta)
  bool pulls_ok = true;
  for (int i = 0; i < static_cast<int>(res.pull.size()); ++i) {
    if (std::abs(res.pull[i]) > 1e-9) {
      pulls_ok = false;
      std::cerr << "  pull[" << i << "] = " << res.pull[i] << " (atteso 0)\n";
    }
  }
  check(pulls_ok, "tutti i pull ≈ 0", T);
}

//  Test T2: retta z = 0.5*y + 3, cov isotropa
static void test_T2() {
  std::cout << "\n[T2] Retta perfetta z = 0.5*y + 3  (a=0.5, b=3), cov isotropa\n";
  const std::string T = "T2";

  const double a_true = 0.5, b_true = 3.0;
  const double var = 0.04;
  auto trace       = make_perfect_trace(a_true, b_true, 31, var, var, 0.0);

  auto res = riptide::fit_trace(trace);

  near(res.a, a_true, 1e-9, "a", T);
  near(res.b, b_true, 1e-9, "b", T);
  near(res.chi2, 0.0, 1e-9, "chi2", T);
  check(res.converged, "converged", T);
  check(res.ndof == 29, "ndof = 29", T);

  // sigma_a e sigma_b devono essere > 0 (funzione della covarianza)
  check(res.sigma_a > 0.0, "sigma_a > 0", T);
  check(res.sigma_b > 0.0, "sigma_b > 0", T);
}

//  Test T3: retta z = 2*y - 1, cov anisotropa (sigma_y >> sigma_z)
//  Quando tutti i punti giacciono ESATTAMENTE sulla retta, i residui
//  perpendicolari sono zero indipendentemente dai pesi.  Il fit deve
//  quindi restituire a=2, b=-1 anche con cov anisotropa.
static void test_T3() {
  std::cout << "\n[T3] Retta perfetta z = 2*y - 1, cov anisotropa (var_y=1.0, var_z=0.0001)\n";
  const std::string T = "T3";

  const double a_true = 2.0, b_true = -1.0;
  // var_y >> var_z: sigma_y = 1 mm, sigma_z = 0.01 mm
  auto trace = make_perfect_trace(a_true, b_true, 21, 1.0, 1e-4, 0.0);

  auto res = riptide::fit_trace(trace);

  // Con punti esatti sulla retta la soluzione e' unica indipendentemente dai pesi
  near(res.a, a_true, 1e-7, "a", T);
  near(res.b, b_true, 1e-7, "b", T);
  near(res.chi2, 0.0, 1e-9, "chi2", T);
  check(res.converged, "converged", T);
}

//  Test T4: outlier singolo
//  N punti su z=0, il punto centrale spostato di delta = 0.5 mm in z.
//  Il pull di quel punto deve essere ≈ delta / sigma_d.
//  Il chi2 deve essere ≈ (delta/sigma_d)^2  (gli altri punti hanno residuo ~0).
static void test_T4() {
  std::cout << "\n[T4] Outlier singolo: un punto spostato di delta = 0.5 mm\n";
  const std::string T = "T4";

  const int N        = 21;
  const double sigma = 0.1; // mm
  const double var   = sigma * sigma;
  const double delta = 0.5; // mm (spostamento in z del punto centrale)

  // a=0, b=0 → z=0; cov isotropa (n_hat = (0,1), sigma_d = sigma_z = sigma)
  auto trace = make_perfect_trace(0.0, 0.0, N, var, var, 0.0);

  // Sposta il punto centrale (indice N/2 = 10)
  const int idx_out = N / 2;
  trace[idx_out].mu_z += delta;

  auto res = riptide::fit_trace(trace);

  // Il pull dell'outlier deve essere vicino a delta/sigma
  // (il fit tirera' leggermente la retta verso l'outlier, quindi il pull
  //  effettivo e' leggermente < delta/sigma; la tolleranza e' del 10%)
  double expected_pull = delta / sigma;
  check(res.pull.size() == static_cast<size_t>(N), "size pull == N", T);
  if (!res.pull.empty()) {
    double p_out = std::abs(res.pull[idx_out]);
    check(p_out > 0.8 * expected_pull && p_out <= expected_pull + 1.0,
          "pull outlier nell'intervallo atteso [" + std::to_string(0.8 * expected_pull) + ", "
              + std::to_string(expected_pull + 1.0) + "], got " + std::to_string(p_out),
          T);
  }

  // Gli altri pull devono essere piccoli (< 2 in valore assoluto)
  bool others_ok = true;
  for (int i = 0; i < N; ++i) {
    if (i == idx_out)
      continue;
    if (std::abs(res.pull[i]) > 2.0) {
      others_ok = false;
      std::cerr << "  pull[" << i << "] = " << res.pull[i] << " (atteso < 2)\n";
    }
  }
  check(others_ok, "pull dei punti non-outlier < 2", T);

  // chi2/ndof deve essere dominato dall'outlier: circa delta^2/sigma^2 / (N-2)
  double chi2_expected_approx = delta * delta / var / (N - 2);
  check(res.chi2_ndof > 0.1 * chi2_expected_approx,
        "chi2/ndof > 10% del valore atteso " + std::to_string(chi2_expected_approx) + ", got "
            + std::to_string(res.chi2_ndof),
        T);
}

//  Test T5: eccezione su traccia troppo corta
static void test_T5() {
  std::cout << "\n[T5] Eccezione su traccia con < 3 punti\n";
  const std::string T = "T5";

  // Traccia vuota
  {
    bool threw = false;
    try {
      riptide::fit_trace({});
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    check(threw, "traccia vuota → std::invalid_argument", T);
  }

  // Traccia con 2 punti
  {
    auto trace2 = make_perfect_trace(1.0, 0.0, 2, 0.01, 0.01);
    bool threw  = false;
    try {
      riptide::fit_trace(trace2);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    check(threw, "2 punti → std::invalid_argument", T);
  }

  // Traccia con esattamente 3 punti: NON deve lanciare
  {
    auto trace3 = make_perfect_trace(1.0, 0.0, 3, 0.01, 0.01);
    bool threw  = false;
    try {
      riptide::fit_trace(trace3);
    } catch (...) {
      threw = true;
    }
    check(!threw, "3 punti → nessuna eccezione", T);
  }
}

//  Test T6: covarianza degenere (var = 0)
//  Con Sigma_i = 0 il floor interno sd2_floor = 1e-6 mm^2 deve evitare
//  divisioni per zero e permettere la convergenza.
static void test_T6() {
  std::cout << "\n[T6] Covarianza degenere (var=0): nessun crash, convergenza\n";
  const std::string T = "T6";

  const double a_true = 0.3, b_true = 1.5;
  auto trace = make_perfect_trace(a_true, b_true, 11, 0.0, 0.0, 0.0);

  bool threw = false;
  riptide::LineFitResult res{};
  try {
    res = riptide::fit_trace(trace);
  } catch (...) {
    threw = true;
  }
  check(!threw, "nessuna eccezione con var=0", T);
  check(res.converged, "converged anche con var=0", T);
  near(res.a, a_true, 1e-7, "a", T);
  near(res.b, b_true, 1e-7, "b", T);
}

//  Test T7: verifica formula sigma_d con cov non diagonale
//  Per una retta z = 0 (a=0, b=0), n_hat = (0, 1).
//  sigma_{d,i}^2 = 0^2 * cov_yy + 2*0*1*cov_yz + 1^2 * cov_zz = cov_zz
//  Quindi sigma_d = sqrt(cov_zz) indipendentemente da cov_yy e cov_yz.
//  Con punti su z=0 e uno spostato di delta:
//    pull_outlier ≈ delta / sqrt(cov_zz_outlier)
static void test_T7() {
  std::cout << "\n[T7] Verifica formula sigma_d^2 = n_hat^T Sigma n_hat con cov non diagonale\n";
  const std::string T = "T7";

  const int N         = 11;
  const double cov_yy = 4.0;  // sigma_y = 2 mm (irrilevante per n_hat=(0,1))
  const double cov_zz = 0.09; // sigma_z = 0.3 mm (questo e' sigma_d)
  const double cov_yz = 0.5;  // correlazione (irrilevante per n_hat=(0,1))
  const double delta  = 0.6;  // spostamento outlier in z [mm]

  // Tutti i punti su z=0, cov non diagonale identica
  std::vector<riptide::TracePoint> trace;
  trace.reserve(N);
  for (int i = 0; i < N; ++i) {
    double y = -5.0 + 10.0 * i / (N - 1);
    trace.push_back(make_point(i, y, 0.0, cov_yy, cov_zz, cov_yz));
  }
  // Outlier al centro
  const int idx   = N / 2;
  trace[idx].mu_z = delta;

  auto res = riptide::fit_trace(trace);

  // Con a ≈ 0 e n_hat ≈ (0,1): sigma_d_outlier ≈ sqrt(cov_zz) = 0.3 mm
  double expected_sigma_d = std::sqrt(cov_zz);
  double expected_pull    = delta / expected_sigma_d; // ≈ 2.0

  // Tolleranza piu' larga: il fit sposta leggermente la retta verso l'outlier
  check(std::abs(res.residual_sig[idx] - expected_sigma_d) < 0.05 * expected_sigma_d,
        "sigma_d outlier ≈ sqrt(cov_zz) = " + std::to_string(expected_sigma_d) + ", got "
            + std::to_string(res.residual_sig[idx]),
        T);

  check(std::abs(res.pull[idx]) > 0.7 * expected_pull,
        "pull outlier > 70% di " + std::to_string(expected_pull) + ", got "
            + std::to_string(std::abs(res.pull[idx])),
        T);
}

// Test TQ: compute_Q
// Non possiamo usare un PSFDatabase reale (nessun file ROOT disponibile nel
// test runner), ma possiamo costruire un PSFDatabase sintetico in memoria
// e verificare le proprietà analitiche di Q.
//
// Setup sintetico:
//   - Una sola configurazione lenti: x1=50, x2=120
//   - PSF "ideale": mu_y = y_source, mu_z = 0, cov = diag(eps, 0, eps)
//     (sistema che mappa perfettamente la posizione radiale sull'asse y)
//   - Traccia ideale a y0=5: r(t) = sqrt(25 + t^2)
//     → mu_y(t) = r(t), mu_z(t) = 0  → traccia curva su y, piatta su z
//     → fit z = a*y + b con tutti i z=0: b=0, a=0, chi2 ≈ 0
//   - Q = sum_y0 chi2 ≈ 0

// Costruisce un PSFDatabase sintetico con PSF mu_y = y_source, mu_z = 0.
static riptide::PSFDatabase make_synthetic_db(double x1 = 50.0, double x2 = 120.0,
                                              double r_max = 16.0, double dr = 0.1,
                                              double cov_val = 1e-6) {
  riptide::PSFDatabase db;
  riptide::LensConfig cfg{x1, x2};
  auto& pts = db[cfg];

  const double eps = dr * 1e-9;
  for (double r = 0.0; r <= r_max + eps; r += dr) {
    riptide::PSFPoint p;
    p.y_source = r;
    p.mu_y     = r;   // mappa lineare: mu_y = r
    p.mu_z     = 0.0; // nessuna deflessione in z
    p.cov_yy   = cov_val;
    p.cov_yz   = 0.0;
    p.cov_zz   = cov_val;
    pts.push_back(p);
  }
  return db;
}

static void test_TQ1() {
  std::cout << "\n[TQ1] compute_Q con PSF ideale (chi2 ≈ 0 per ogni traccia)\n";
  const std::string T = "TQ1";

  auto db = make_synthetic_db();
  riptide::LensConfig cfg{50.0, 120.0};

  riptide::QConfig qcfg;
  qcfg.y0_min   = 0.0;
  qcfg.y0_max   = 10.0;
  qcfg.dy0      = 1.0; // 11 tracce: y0 = 0, 1, ..., 10
  qcfg.trace_L  = 10.0;
  qcfg.trace_dt = 0.5; // 21 punti per traccia

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

  // Con PSF ideale (tutti i punti su una retta z=0), chi2 deve essere ≈ 0
  // La tolleranza è generosa perché cov_val=1e-6 porta chi2/ndof ≈ 0 ma Q = sum chi2
  near(res.Q, 0.0, 1e-4, "Q ≈ 0 (PSF ideale)", T);

  // Verifica struttura vettori per-traccia
  check(res.y0_used.size() == 11, "y0_used.size() == 11", T);
  check(res.chi2_per_y0.size() == 11, "chi2_per_y0.size() == 11", T);
  check(res.chi2_ndof_per_y0.size() == 11, "chi2_ndof.size() == 11", T);

  // Q deve essere la somma esatta dei chi2 per-traccia
  double Q_sum = 0.0;
  for (double c : res.chi2_per_y0)
    Q_sum += c;
  near(res.Q, Q_sum, 1e-12, "Q == sum(chi2_per_y0)", T);
}

static void test_TQ2() {
  std::cout << "\n[TQ2] compute_Q con y0_values esplicito\n";
  const std::string T = "TQ2";

  auto db = make_synthetic_db();
  riptide::LensConfig cfg{50.0, 120.0};

  riptide::QConfig qcfg;
  qcfg.y0_values = {2.0, 5.0, 8.0}; // solo 3 tracce
  qcfg.trace_L   = 10.0;
  qcfg.trace_dt  = 0.5;

  auto res = riptide::compute_Q(cfg, db, qcfg);

  check(res.n_traces == 3, "n_traces == 3", T);
  check(res.y0_used.size() == 3, "y0_used.size() == 3", T);
  near(res.Q, 0.0, 1e-4, "Q ≈ 0", T);
}

static void test_TQ3() {
  std::cout << "\n[TQ3] compute_Q con configurazione non presente nel db → eccezione\n";
  const std::string T = "TQ3";

  auto db = make_synthetic_db(50.0, 120.0);
  riptide::LensConfig cfg_wrong{99.0, 199.0}; // non esiste nel db

  bool threw = false;
  try {
    riptide::compute_Q(cfg_wrong, db);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  check(threw, "eccezione std::invalid_argument per cfg non presente", T);
}

static void test_TQ4() {
  std::cout << "\n[TQ4] Q monotona: PSF con dispersione crescente → Q cresce con cov\n";
  const std::string T = "TQ4";

  // PSF con cov_zz = sigma^2: mu_y = r, mu_z linearmente crescente con r
  // → le tracce non sono più piatte in z → chi2 > 0, proporzionale a 1/sigma^2
  // Costruiamo due DB: uno con sigma piccola (cov=1e-4) e uno con sigma grande (cov=1.0)
  // Con sigma piccola i pesi sono grandi → chi2 grande per lo stesso residuo
  // Con sigma grande i pesi sono piccoli → chi2 piccolo
  // Quindi Q_small_sigma > Q_large_sigma  [se i residui sono non nulli]

  // PSF con mu_z = 0.1 * y_source  (introduce curvatura non lineare)
  auto make_curved_db = [](double cov_zz_val) {
    riptide::PSFDatabase db;
    riptide::LensConfig cfg{50.0, 120.0};
    auto& pts = db[cfg];
    for (double r = 0.0; r <= 16.0 + 1e-9; r += 0.1) {
      riptide::PSFPoint p;
      p.y_source = r;
      p.mu_y     = r;
      p.mu_z     = 0.05 * r * r; // curvatura quadratica in z → non lineare
      p.cov_yy   = 1e-6;
      p.cov_yz   = 0.0;
      p.cov_zz   = cov_zz_val;
      pts.push_back(p);
    }
    return db;
  };

  riptide::LensConfig cfg{50.0, 120.0};
  riptide::QConfig qcfg;
  qcfg.y0_min   = 2.0;
  qcfg.y0_max   = 8.0;
  qcfg.dy0      = 2.0;
  qcfg.trace_L  = 10.0;
  qcfg.trace_dt = 0.5;

  auto db_tight = make_curved_db(1e-4);
  auto db_loose = make_curved_db(1.0);

  auto res_tight = riptide::compute_Q(cfg, db_tight, qcfg);
  auto res_loose = riptide::compute_Q(cfg, db_loose, qcfg);

  // Con PSF più "stretta" (cov_zz piccolo) i residui di curvatura pesano di più → Q maggiore
  check(res_tight.Q > res_loose.Q,
        "Q(sigma_piccola) > Q(sigma_grande) per traccia curva: " + std::to_string(res_tight.Q)
            + " > " + std::to_string(res_loose.Q),
        T);
  // Entrambe devono avere n_failed == 0
  check(res_tight.n_failed == 0, "n_failed == 0 (tight)", T);
  check(res_loose.n_failed == 0, "n_failed == 0 (loose)", T);
}

//  main
int main() {
  std::cout << "  test_fit_trace — riptide::fit_trace() unit tests\n";

  test_T1();
  test_T2();
  test_T3();
  test_T4();
  test_T5();
  test_T6();
  test_T7();
  test_TQ1();
  test_TQ2();
  test_TQ3();
  test_TQ4();

  std::cout << "  Risultato: " << g_n_pass << " PASS, " << g_n_fail << " FAIL\n";

  return (g_n_fail == 0) ? 0 : 1;
}
