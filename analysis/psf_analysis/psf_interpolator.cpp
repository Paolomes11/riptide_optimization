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

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace riptide {

//  Caricamento database ────────────────────────────────────────────────────────

PSFDatabase load_psf_database(const std::string& root_path) {
  TFile* f = TFile::Open(root_path.c_str(), "READ");
  if (!f || f->IsZombie())
    throw std::runtime_error("load_psf_database: impossibile aprire " + root_path);

  TTree* tree = (TTree*)f->Get("PSF");
  if (!tree)
    throw std::runtime_error("load_psf_database: TTree 'PSF' non trovato in " + root_path);

  Int_t config_id;
  Double_t x1, x2;
  Float_t y_source_f;
  Double_t mean_y, mean_z;
  Double_t cov_yy, cov_yz, cov_zz;
  Int_t n_hits_filtered;

  tree->SetBranchAddress("config_id", &config_id);
  tree->SetBranchAddress("x1", &x1);
  tree->SetBranchAddress("x2", &x2);
  tree->SetBranchAddress("y_source", &y_source_f);
  tree->SetBranchAddress("mean_y", &mean_y);
  tree->SetBranchAddress("mean_z", &mean_z);
  tree->SetBranchAddress("cov_yy", &cov_yy);
  tree->SetBranchAddress("cov_yz", &cov_yz);
  tree->SetBranchAddress("cov_zz", &cov_zz);
  tree->SetBranchAddress("n_hits_filtered", &n_hits_filtered);

  bool on_det_branch_exists = (tree->GetBranch("on_detector") != nullptr);
  bool on_det_val           = true;
  if (on_det_branch_exists) {
    tree->SetBranchAddress("on_detector", &on_det_val);
    std::cout << "load_psf_database: branch 'on_detector' trovato.\n";
  } else {
    std::cout << "load_psf_database: branch 'on_detector' assente, "
                 "uso fallback n_hits_filtered >= 10.\n";
  }

  PSFDatabase db;
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    bool on_det = on_det_branch_exists ? on_det_val : (n_hits_filtered >= 10);
    LensConfig cfg{x1, x2};
    db[cfg].push_back(
        {static_cast<double>(y_source_f), mean_y, mean_z, cov_yy, cov_yz, cov_zz, on_det});
  }

  for (auto& [cfg, pts] : db)
    std::sort(pts.begin(), pts.end(),
              [](const PSFPoint& a, const PSFPoint& b) { return a.y_source < b.y_source; });

  f->Close();
  std::cout << "PSF database caricato: " << db.size() << " configurazioni\n";
  return db;
}

//  Ricerca configurazione più vicina ───────────────────────────────────────────

LensConfig find_nearest_config(const LensConfig& cfg, const PSFDatabase& db) {
  if (db.empty())
    throw std::runtime_error("find_nearest_config: database vuoto");

  const LensConfig* best = nullptr;
  double best_dist       = std::numeric_limits<double>::max();

  for (const auto& [key, _] : db) {
    double d = std::hypot(key.x1 - cfg.x1, key.x2 - cfg.x2);
    if (d < best_dist) {
      best_dist = d;
      best      = &key;
    }
  }

  if (best_dist > 1e-4)
    std::cout << "[WARNING] Config (x1=" << cfg.x1 << ", x2=" << cfg.x2
              << ") non trovata. Più vicina: (x1=" << best->x1 << ", x2=" << best->x2
              << "), dist=" << best_dist << " mm\n";

  return *best;
}

//  Interpolazione ──────────────────────────────────────────────────────────────

PSFValue interpolate(double r, const LensConfig& cfg, const PSFDatabase& db) {
  auto it = db.find(cfg);
  if (it == db.end()) {
    std::ostringstream oss;
    oss << "interpolate: config (x1=" << cfg.x1 << ", x2=" << cfg.x2 << ") non trovata";
    throw std::out_of_range(oss.str());
  }

  const auto& pts = it->second;
  if (pts.empty())
    throw std::out_of_range("interpolate: nessun punto PSF per questa configurazione");

  auto upper = std::lower_bound(pts.begin(), pts.end(), r,
                                [](const PSFPoint& p, double v) { return p.y_source < v; });

  if (upper == pts.begin()) {
    const auto& p = pts.front();
    return {p.mu_y, p.mu_z, {p.cov_yy, p.cov_yz, p.cov_zz}, p.on_detector};
  }
  if (upper == pts.end()) {
    const auto& p = pts.back();
    return {p.mu_y, p.mu_z, {p.cov_yy, p.cov_yz, p.cov_zz}, p.on_detector};
  }

  const PSFPoint& hi = *upper;
  const PSFPoint& lo = *std::prev(upper);
  double dr          = hi.y_source - lo.y_source;
  double alpha       = (dr > 1e-12) ? (r - lo.y_source) / dr : 0.0;
  auto lerp          = [alpha](double a, double b) { return a + alpha * (b - a); };

  return {lerp(lo.mu_y, hi.mu_y),
          lerp(lo.mu_z, hi.mu_z),
          {lerp(lo.cov_yy, hi.cov_yy), lerp(lo.cov_yz, hi.cov_yz), lerp(lo.cov_zz, hi.cov_zz)},
          lo.on_detector && hi.on_detector};
}

//  Costruzione traccia ─────────────────────────────────────────────────────────

std::vector<TracePoint> build_trace(double y0, const LensConfig& cfg, const PSFDatabase& db,
                                    double L, double dt) {
  const int N = static_cast<int>(std::round(L / dt)) + 1;
  std::vector<TracePoint> trace;
  trace.reserve(N);

  for (int i = 0; i < N; ++i) {
    double t     = std::min(-L / 2.0 + i * dt, L / 2.0);
    double r     = std::hypot(y0, t);
    PSFValue psf = interpolate(r, cfg, db);
    trace.push_back({t, r, psf.mu_y, psf.mu_z, psf.cov, psf.on_detector});
  }
  return trace;
}

//  Validità traccia ────────────────────────────────────────────────────────────

bool is_trace_valid(const std::vector<TracePoint>& trace, double point_valid_fraction) {
  if (trace.empty())
    return false;
  int nv = 0;
  for (const auto& pt : trace)
    if (pt.valid)
      ++nv;
  return static_cast<double>(nv) / static_cast<double>(trace.size()) >= point_valid_fraction;
}

// ─── Fit lineare pesato ODR ───────────────────────────────────────────────────

static bool solve_wls(const std::vector<double>& y, const std::vector<double>& z,
                      const std::vector<double>& w, double& a_out, double& b_out, double& var_a,
                      double& var_b, double& cov_ab) {
  const int N = static_cast<int>(y.size());
  if (N < 2)
    return false;

  double S1 = 0, Sy = 0, Sz = 0, Syy = 0, Syz = 0;
  for (int i = 0; i < N; ++i) {
    S1 += w[i];
    Sy += w[i] * y[i];
    Sz += w[i] * z[i];
    Syy += w[i] * y[i] * y[i];
    Syz += w[i] * y[i] * z[i];
  }
  double det = Syy * S1 - Sy * Sy;
  if (std::abs(det) < 1e-30)
    return false;

  a_out  = (Syz * S1 - Sz * Sy) / det;
  b_out  = (Syy * Sz - Sy * Syz) / det;
  var_a  = S1 / det;
  var_b  = Syy / det;
  cov_ab = -Sy / det;
  return true;
}

LineFitResult fit_trace(const std::vector<TracePoint>& trace, int max_iter, double tol) {
  // Estrai solo i punti validi
  std::vector<double> vy, vz;
  std::vector<Cov2> vcov;
  vy.reserve(trace.size());
  vz.reserve(trace.size());
  vcov.reserve(trace.size());
  for (const auto& pt : trace) {
    if (pt.valid) {
      vy.push_back(pt.mu_y);
      vz.push_back(pt.mu_z);
      vcov.push_back(pt.cov);
    }
  }

  const int N = static_cast<int>(vy.size());
  if (N < 3)
    throw std::invalid_argument("fit_trace: punti validi insufficienti (trovati: "
                                + std::to_string(N)
                                + ", richiesti: 3). Verifica point_valid_fraction.");

  LineFitResult res{};
  res.n_iter        = 0;
  res.converged     = false;
  res.n_points_used = N;

  // Stima iniziale OLS
  {
    std::vector<double> w1(N, 1.0);
    double va, vb, cab;
    if (!solve_wls(vy, vz, w1, res.a, res.b, va, vb, cab)) {
      res.a = 0.0;
      res.b = vz[0];
    }
  }

  // Loop IRLS
  std::vector<double> ww(N);
  for (int iter = 0; iter < max_iter; ++iter) {
    res.n_iter    = iter + 1;
    double a_prev = res.a;
    double norm   = std::sqrt(1.0 + res.a * res.a);
    double ny     = -res.a / norm;
    double nz     = 1.0 / norm;

    for (int i = 0; i < N; ++i) {
      double sd2 = ny * ny * vcov[i].yy + 2.0 * ny * nz * vcov[i].yz + nz * nz * vcov[i].zz;
      ww[i]      = 1.0 / std::max(sd2, 1e-6);
    }

    double new_a, new_b;
    if (!solve_wls(vy, vz, ww, new_a, new_b, res.sigma_a, res.sigma_b, res.cov_ab)) {
      std::cerr << "[fit_trace] WLS singolare iter=" << iter + 1 << "\n";
      break;
    }
    res.a = new_a;
    res.b = new_b;
    if (std::abs(res.a - a_prev) < tol) {
      res.converged = true;
      break;
    }
  }

  // χ², residui, pull
  double normf = std::sqrt(1.0 + res.a * res.a);
  double nyf   = -res.a / normf;
  double nzf   = 1.0 / normf;
  res.chi2     = 0.0;
  res.residuals.resize(N);
  res.residual_sig.resize(N);
  res.pull.resize(N);

  for (int i = 0; i < N; ++i) {
    double d   = (res.a * vy[i] - vz[i] + res.b) / normf;
    double sd2 = nyf * nyf * vcov[i].yy + 2.0 * nyf * nzf * vcov[i].yz + nzf * nzf * vcov[i].zz;
    double sd  = std::sqrt(std::max(sd2, 1e-6));
    res.residuals[i]    = d;
    res.residual_sig[i] = sd;
    res.pull[i]         = d / sd;
    res.chi2 += d * d / std::max(sd2, 1e-6);
  }
  res.ndof      = N - 2;
  res.chi2_ndof = (res.ndof > 0) ? res.chi2 / res.ndof : 0.0;
  res.sigma_a   = std::sqrt(std::max(res.sigma_a, 0.0));
  res.sigma_b   = std::sqrt(std::max(res.sigma_b, 0.0));
  return res;
}

// ─── Temporal unfolding ───────────────────────────────────────────────────────
//
// Perché l'approccio i·δz funziona (e perché il codice precedente era rotto)
// ──────────────────────────────────────────────────────────────────────────────
// Se una traccia si ripiega: μ_y = [2,4,6,4,2], μ_z ≈ [0,0,0,0,0]
// Con unfolding z̃_i = μ_z_i + i·δz  i punti diventano (nel piano y,z̃):
//   (2,0) (4,δz) (6,2δz) (4,3δz) (2,4δz)
// I punti (2,0) e (2,4δz) hanno la stessa y ma z̃ diverse: nessuna retta
// z̃ = a·y + b può passare per entrambi → χ² elevato. ✓
//
// Per una traccia monotona non ripiegata i punti ruotano ma rimangono
// approssimativamente collineari — il fit ODR troverà una retta con χ² piccolo.
// Questo è il comportamento DESIDERATO: premiamo le tracce che si mappano
// in modo monotono e regolare sul detector.
//
// Il bug nel codice precedente era una dangling reference:
//   const auto& trace_for_fit = apply_unfolding(trace, qcfg);  // temporaneo distrutto!
// trace_for_fit puntava a memoria già liberata, quindi fit_trace riceveva dati
// casuali (o la traccia originale se il compilatore riutilizzava lo stack).
// La correzione: usare sempre una variabile std::vector con lifetime esplicita.
//
// Il parametro z_unfold_step controlla δz:
//   = 0 (default): δz = trace_L / (N-1)  → offset totale = trace_L mm
//   > 0: δz fisso in mm
// Disabilita con apply_temporal_unfolding = false.

static std::vector<TracePoint> apply_unfolding(const std::vector<TracePoint>& trace,
                                               const QConfig& qcfg) {
  std::vector<TracePoint> unfolded = trace; // COPIA esplicita, lifetime garantita
  const int N                      = static_cast<int>(unfolded.size());
  if (N < 2)
    return unfolded;

  double dz =
      (qcfg.z_unfold_step > 0.0) ? qcfg.z_unfold_step : qcfg.trace_L / static_cast<double>(N - 1);

  for (int i = 0; i < N; ++i)
    unfolded[i].mu_z += static_cast<double>(i) * dz;
  //  Cov invariata: l'offset i·dz è deterministico, non stocastico.
  //  Cov(z̃_i, z̃_j) = Cov(z_i + i·dz, z_j + j·dz) = Cov(z_i, z_j)

  return unfolded;
}

// ─── compute_Q ───────────────────────────────────────────────────────────────

QResult compute_Q(const LensConfig& cfg, const PSFDatabase& db, const QConfig& qcfg,
                  bool include_non_converged) {
  if (db.find(cfg) == db.end()) {
    std::ostringstream oss;
    oss << "compute_Q: config (x1=" << cfg.x1 << ", x2=" << cfg.x2 << ") non trovata.";
    throw std::invalid_argument(oss.str());
  }

  // Lista y0
  std::vector<double> y0_list;
  if (!qcfg.y0_values.empty()) {
    y0_list = qcfg.y0_values;
  } else {
    const double eps = qcfg.dy0 * 1e-9;
    for (double y0 = qcfg.y0_min; y0 <= qcfg.y0_max + eps; y0 += qcfg.dy0)
      y0_list.push_back(y0);
  }
  if (y0_list.empty())
    throw std::invalid_argument("compute_Q: lista di y0 vuota");

  const int total_y0 = static_cast<int>(y0_list.size());

  QResult res{};
  res.Q            = 0.0;
  res.n_traces     = 0;
  res.n_failed     = 0;
  res.n_invalid    = 0;
  res.config_valid = true;
  res.y0_used.reserve(total_y0);
  res.chi2_per_y0.reserve(total_y0);
  res.chi2_ndof_per_y0.reserve(total_y0);
  res.trace_valid_flags.reserve(total_y0);

  for (double y0 : y0_list) {
    // 1. Traccia fisica
    std::vector<TracePoint> trace;
    try {
      trace = build_trace(y0, cfg, db, qcfg.trace_L, qcfg.trace_dt);
    } catch (const std::exception& e) {
      res.warnings.push_back({QWarning::Kind::BuildTraceFailed, y0,
                              std::numeric_limits<double>::quiet_NaN(), e.what()});
      if (qcfg.verbose)
        std::cerr << "[compute_Q] build_trace failed y0=" << y0 << ": " << e.what() << "\n";
      ++res.n_failed;
      res.trace_valid_flags.push_back(false);
      continue;
    }

    // 2. Validità traccia (sulla traccia fisica)
    bool trace_ok = is_trace_valid(trace, qcfg.point_valid_fraction);
    res.trace_valid_flags.push_back(trace_ok);
    if (!trace_ok) {
      int nv = 0;
      for (const auto& pt : trace)
        if (pt.valid)
          ++nv;
      double frac = static_cast<double>(nv) / static_cast<double>(trace.size());
      res.warnings.push_back(
          {QWarning::Kind::TraceInvalid, y0, std::numeric_limits<double>::quiet_NaN(),
           "y0=" + std::to_string(y0) + " ha solo " + std::to_string(static_cast<int>(frac * 100))
               + "% punti on_detector"});
      if (qcfg.verbose)
        std::cerr << "[compute_Q] " << res.warnings.back().message << "\n";
      ++res.n_invalid;
      continue;
    }

    // 3. Temporal unfolding: z̃_i = μ_z_i + i·δz
    //    Operazione su una copia locale (trace fisica NON modificata).
    //    Il vettore `unfolded` ha lifetime garantita fino alla fine del blocco.
    std::vector<TracePoint> unfolded;         // lifetime esplicita — nessuna dangling ref
    const std::vector<TracePoint>* fit_input; // punta a trace o unfolded
    if (qcfg.apply_temporal_unfolding) {
      unfolded  = apply_unfolding(trace, qcfg);
      fit_input = &unfolded;
    } else {
      fit_input = &trace;
    }

    // 4. Fit ODR pesato
    LineFitResult fit;
    try {
      fit = fit_trace(*fit_input, qcfg.fit_max_iter, qcfg.fit_tol);
    } catch (const std::exception& e) {
      res.warnings.push_back(
          {QWarning::Kind::FitFailed, y0, std::numeric_limits<double>::quiet_NaN(), e.what()});
      if (qcfg.verbose)
        std::cerr << "[compute_Q] fit_trace threw y0=" << y0 << ": " << e.what() << "\n";
      ++res.n_failed;
      continue;
    }

    // 5. Convergenza
    if (!fit.converged && !include_non_converged) {
      res.warnings.push_back({QWarning::Kind::FitNotConverged, y0, fit.a,
                              "fit non convergente dopo " + std::to_string(fit.n_iter) + " iter"});
      if (qcfg.verbose)
        std::cerr << "[compute_Q] " << res.warnings.back().message << " y0=" << y0 << "\n";
      ++res.n_failed;
      continue;
    }

    res.Q += fit.chi2;
    res.y0_used.push_back(y0);
    res.chi2_per_y0.push_back(fit.chi2);
    res.chi2_ndof_per_y0.push_back(fit.chi2_ndof);
    ++res.n_traces;
  }

  // Validità configurazione
  int n_att = total_y0 - res.n_failed;
  if (n_att > 0) {
    double vf        = static_cast<double>(res.n_traces) / static_cast<double>(n_att);
    res.config_valid = (vf >= qcfg.trace_valid_fraction);
    if (!res.config_valid) {
      res.warnings.push_back(
          {QWarning::Kind::ConfigInvalid, -1.0, std::numeric_limits<double>::quiet_NaN(),
           "solo " + std::to_string(static_cast<int>(vf * 100)) + "% tracce valide"});
      if (qcfg.verbose)
        std::cerr << "[compute_Q] " << res.warnings.back().message << "\n";
    }
  } else {
    res.config_valid = false;
  }

  return res;
}

} // namespace riptide