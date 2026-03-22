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

//  Caricamento database

PSFDatabase load_psf_database(const std::string& root_path) {
  TFile* f = TFile::Open(root_path.c_str(), "READ");
  if (!f || f->IsZombie())
    throw std::runtime_error("load_psf_database: impossibile aprire " + root_path);

  TTree* tree = (TTree*)f->Get("PSF");
  if (!tree)
    throw std::runtime_error("load_psf_database: TTree 'PSF' non trovato in " + root_path);

  // Branch del TTree PSF — y_source è Float_t, il resto Double_t
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

  PSFDatabase db;
  const Long64_t N = tree->GetEntries();

  for (Long64_t i = 0; i < N; ++i) {
    tree->GetEntry(i);

    if (n_hits_filtered < 10)
      continue;

    LensConfig cfg{x1, x2};
    db[cfg].push_back({static_cast<double>(y_source_f), mean_y, mean_z, cov_yy, cov_yz, cov_zz});
  }

  // Ordina ogni vettore per y_source crescente (necessario per la ricerca binaria)
  for (auto& [cfg, pts] : db) {
    std::sort(pts.begin(), pts.end(),
              [](const PSFPoint& a, const PSFPoint& b) { return a.y_source < b.y_source; });
  }

  f->Close();

  std::cout << "PSF database caricato: " << db.size() << " configurazioni\n";
  return db;
}

//  Ricerca configurazione più vicina

LensConfig find_nearest_config(const LensConfig& cfg, const PSFDatabase& db) {
  if (db.empty())
    throw std::runtime_error("find_nearest_config: database vuoto");

  const LensConfig* best = nullptr;
  double best_dist       = std::numeric_limits<double>::max();

  for (const auto& [key, _] : db) {
    double dx1  = key.x1 - cfg.x1;
    double dx2  = key.x2 - cfg.x2;
    double dist = std::sqrt(dx1 * dx1 + dx2 * dx2);
    if (dist < best_dist) {
      best_dist = dist;
      best      = &key;
    }
  }

  if (best_dist > 1e-4) {
    std::cout << "[WARNING] Configurazione (x1=" << cfg.x1 << ", x2=" << cfg.x2
              << ") non trovata nel database.\n"
              << "          Configurazione più vicina: (x1=" << best->x1 << ", x2=" << best->x2
              << "), distanza=" << best_dist << " mm\n";
  }

  return *best;
}

//  Interpolazione

PSFValue interpolate(double r, const LensConfig& cfg, const PSFDatabase& db) {
  auto it = db.find(cfg);
  if (it == db.end()) {
    std::ostringstream oss;
    oss << "interpolate: configurazione (x1=" << cfg.x1 << ", x2=" << cfg.x2
        << ") non trovata nel database";
    throw std::out_of_range(oss.str());
  }

  const auto& pts = it->second;
  if (pts.empty())
    throw std::out_of_range("interpolate: nessun punto PSF per questa configurazione");

  // Cerca i due punti adiacenti con ricerca binaria su y_source
  // Lower bound: primo punto con y_source >= r
  auto upper = std::lower_bound(pts.begin(), pts.end(), r,
                                [](const PSFPoint& p, double val) { return p.y_source < val; });

  // Caso r <= primo punto: usa il primo punto senza estrapolazione
  if (upper == pts.begin()) {
    const auto& p = pts.front();
    return {p.mu_y, p.mu_z, {p.cov_yy, p.cov_yz, p.cov_zz}};
  }

  // Caso r >= ultimo punto: usa l'ultimo punto senza estrapolazione
  if (upper == pts.end()) {
    const auto& p = pts.back();
    return {p.mu_y, p.mu_z, {p.cov_yy, p.cov_yz, p.cov_zz}};
  }

  // Caso generale: interpolazione lineare tra lower e upper
  const PSFPoint& hi = *upper;
  const PSFPoint& lo = *std::prev(upper);

  double dr    = hi.y_source - lo.y_source;
  double alpha = (dr > 1e-12) ? (r - lo.y_source) / dr : 0.0;

  auto lerp = [alpha](double a, double b) { return a + alpha * (b - a); };

  return {lerp(lo.mu_y, hi.mu_y),
          lerp(lo.mu_z, hi.mu_z),
          {lerp(lo.cov_yy, hi.cov_yy), lerp(lo.cov_yz, hi.cov_yz), lerp(lo.cov_zz, hi.cov_zz)}};
}

//  Costruzione traccia

std::vector<TracePoint> build_trace(double y0, const LensConfig& cfg, const PSFDatabase& db,
                                    double L, double dt) {
  std::vector<TracePoint> trace;

  // Numero di passi: da -L/2 a +L/2 con step dt
  const int N = static_cast<int>(std::round(L / dt)) + 1;
  trace.reserve(N);

  for (int i = 0; i < N; ++i) {
    double t = -L / 2.0 + i * dt;
    // Clamp t all'intervallo [-L/2, L/2] per evitare errori floating point
    if (t > L / 2.0)
      t = L / 2.0;

    double r = std::sqrt(y0 * y0 + t * t);

    PSFValue psf = interpolate(r, cfg, db);

    trace.push_back({t, r, psf.mu_y, psf.mu_z, psf.cov});
  }

  return trace;
}

// Fit lineare pesato ODR della traccia media

static bool solve_wls(const std::vector<double>& y, const std::vector<double>& z,
                      const std::vector<double>& w, double& a_out, double& b_out, double& var_a,
                      double& var_b, double& cov_ab_out) {
  const int N = static_cast<int>(y.size());
  if (N < 2)
    return false;

  // Accumula le somme pesate per il sistema normale 2x2:
  //   [S_yy  S_y ] [a]   [S_yz]
  //   [S_y   S_1 ] [b] = [S_z ]
  double S1  = 0.0; // sum w_i
  double Sy  = 0.0; // sum w_i * y_i
  double Sz  = 0.0; // sum w_i * z_i
  double Syy = 0.0; // sum w_i * y_i^2
  double Syz = 0.0; // sum w_i * y_i * z_i

  for (int i = 0; i < N; ++i) {
    S1 += w[i];
    Sy += w[i] * y[i];
    Sz += w[i] * z[i];
    Syy += w[i] * y[i] * y[i];
    Syz += w[i] * y[i] * z[i];
  }

  // Determinante della matrice normale (deve essere > 0 per sistema non degenere)
  double det = Syy * S1 - Sy * Sy;
  if (std::abs(det) < 1e-30)
    return false; // sistema singolare (tutti y_i uguali o N < 2)

  // Soluzione di Cramer
  a_out = (Syz * S1 - Sz * Sy) / det;
  b_out = (Syy * Sz - Sy * Syz) / det;

  // Matrice di covarianza dei parametri: C = (X^T W X)^{-1}
  //   X^T W X = [[Syy, Sy], [Sy, S1]]
  var_a      = S1 / det;
  var_b      = Syy / det;
  cov_ab_out = -Sy / det;

  return true;
}

LineFitResult fit_trace(const std::vector<TracePoint>& trace, int max_iter, double tol) {
  const int N = static_cast<int>(trace.size());
  if (N < 3)
    throw std::invalid_argument(
        "fit_trace: la traccia deve contenere almeno 3 punti (trovati: " + std::to_string(N) + ")");

  // Estrai i vettori di coordinate dalla traccia
  std::vector<double> vy(N), vz(N);
  for (int i = 0; i < N; ++i) {
    vy[i] = trace[i].mu_y;
    vz[i] = trace[i].mu_z;
  }

  LineFitResult res{};
  res.n_iter    = 0;
  res.converged = false;

  // Step 1: stima iniziale con OLS non pesato
  // Risolvi z = a*y + b con pesi unitari per ottenere la direzione iniziale
  {
    std::vector<double> w_unit(N, 1.0);
    double var_a_tmp, var_b_tmp, cov_ab_tmp;
    bool ok = solve_wls(vy, vz, w_unit, res.a, res.b, var_a_tmp, var_b_tmp, cov_ab_tmp);
    if (!ok) {
      // Fallback: retta orizzontale
      res.a = 0.0;
      res.b = vz[0];
    }
  }

  // Step 2–3: loop IRLS (Iteratively Reweighted Least Squares)
  std::vector<double> weights(N);

  for (int iter = 0; iter < max_iter; ++iter) {
    res.n_iter = iter + 1;

    double a_prev  = res.a;
    double norm_sq = 1.0 + res.a * res.a;
    double norm    = std::sqrt(norm_sq);

    // Componenti del vettore normale alla retta z = a*y + b
    double n_y = -res.a / norm; //  n_y = -a / sqrt(1+a^2)
    double n_z = 1.0 / norm;    //  n_z =  1 / sqrt(1+a^2)

    // Calcola pesi w_i = 1 / sigma_{d,i}^2
    for (int i = 0; i < N; ++i) {
      const Cov2& cov = trace[i].cov;
      // sigma^2_{d,i} = n_hat^T Sigma_i n_hat
      double sd2 = n_y * n_y * cov.yy + 2.0 * n_y * n_z * cov.yz + n_z * n_z * cov.zz;

      // Protezione contro covarianza degenerata o nulla:
      // se sd2 e' troppo piccola usiamo un floor pari a (0.001 mm)^2
      // per evitare pesi infiniti che destabilizzerebbero il fit.
      const double sd2_floor = 1e-6; // (0.001 mm)^2
      weights[i]             = 1.0 / std::max(sd2, sd2_floor);
    }

    // Risolvi il WLS con i pesi aggiornati
    double new_a, new_b;
    bool ok = solve_wls(vy, vz, weights, new_a, new_b, res.sigma_a, res.sigma_b, res.cov_ab);
    if (!ok) {
      // Il sistema e' diventato singolare: mantieni la stima precedente
      std::cerr << "[fit_trace] WARNING: sistema WLS singolare all'iterazione " << iter + 1
                << ", uso stima precedente.\n";
      break;
    }

    res.a = new_a;
    res.b = new_b;

    // Controlla convergenza su variazione di a
    if (std::abs(res.a - a_prev) < tol) {
      res.converged = true;
      break;
    }
  }

  // Step 4: calcolo finale di chi^2, residui e pull
  double norm_final = std::sqrt(1.0 + res.a * res.a);
  double n_y_final  = -res.a / norm_final;
  double n_z_final  = 1.0 / norm_final;

  res.chi2 = 0.0;
  res.residuals.resize(N);
  res.residual_sig.resize(N);
  res.pull.resize(N);

  for (int i = 0; i < N; ++i) {
    // Distanza perpendicolare con segno (positivo se il punto e' "sopra" la retta)
    double d_i = (res.a * vy[i] - vz[i] + res.b) / norm_final;

    // sigma_{d,i} con il vettore normale finale
    const Cov2& cov = trace[i].cov;
    double sd2      = n_y_final * n_y_final * cov.yy + 2.0 * n_y_final * n_z_final * cov.yz
               + n_z_final * n_z_final * cov.zz;

    const double sd2_floor = 1e-6;
    double sd_i            = std::sqrt(std::max(sd2, sd2_floor));

    res.residuals[i]    = d_i;
    res.residual_sig[i] = sd_i;
    res.pull[i]         = d_i / sd_i;

    res.chi2 += (d_i * d_i) / std::max(sd2, sd2_floor);
  }

  // Gradi di liberta' e chi^2 ridotto
  res.ndof      = N - 2;
  res.chi2_ndof = (res.ndof > 0) ? res.chi2 / static_cast<double>(res.ndof) : 0.0;

  // Trasforma sigma_a, sigma_b da varianza a deviazione standard
  res.sigma_a = std::sqrt(std::max(res.sigma_a, 0.0));
  res.sigma_b = std::sqrt(std::max(res.sigma_b, 0.0));

  return res;
}

// compute_Q
QResult compute_Q(const LensConfig& cfg, const PSFDatabase& db, const QConfig& qcfg,
                  bool include_non_converged) {
  if (db.find(cfg) == db.end()) {
    std::ostringstream oss;
    oss << "compute_Q: configurazione (x1=" << cfg.x1 << ", x2=" << cfg.x2
        << ") non trovata nel database.";
    throw std::invalid_argument(oss.str());
  }

  // Costruisce la lista di y0
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

  QResult res{};
  res.Q        = 0.0;
  res.n_traces = 0;
  res.n_failed = 0;
  res.y0_used.reserve(y0_list.size());
  res.chi2_per_y0.reserve(y0_list.size());
  res.chi2_ndof_per_y0.reserve(y0_list.size());

  for (double y0 : y0_list) {
    // build_trace
    std::vector<TracePoint> trace;
    try {
      trace = build_trace(y0, cfg, db, qcfg.trace_L, qcfg.trace_dt);
    } catch (const std::exception& e) {
      QWarning w;
      w.kind    = QWarning::Kind::BuildTraceFailed;
      w.y0      = y0;
      w.a_final = std::numeric_limits<double>::quiet_NaN();
      w.message = e.what();
      res.warnings.push_back(std::move(w));
      if (qcfg.verbose)
        std::cerr << "[compute_Q] build_trace failed y0=" << y0 << ": " << e.what() << "\n";
      ++res.n_failed;
      continue;
    }

    // fit_trace
    LineFitResult fit;
    try {
      fit = fit_trace(trace, qcfg.fit_max_iter, qcfg.fit_tol);
    } catch (const std::exception& e) {
      QWarning w;
      w.kind    = QWarning::Kind::FitFailed;
      w.y0      = y0;
      w.a_final = std::numeric_limits<double>::quiet_NaN();
      w.message = e.what();
      res.warnings.push_back(std::move(w));
      if (qcfg.verbose)
        std::cerr << "[compute_Q] fit_trace threw y0=" << y0 << ": " << e.what() << "\n";
      ++res.n_failed;
      continue;
    }

    // Convergenza
    if (!fit.converged && !include_non_converged) {
      QWarning w;
      w.kind    = QWarning::Kind::FitNotConverged;
      w.y0      = y0;
      w.a_final = fit.a;
      w.message = "fit non convergente dopo " + std::to_string(fit.n_iter) + " iterazioni";
      res.warnings.push_back(std::move(w));
      if (qcfg.verbose)
        std::cerr << "[compute_Q] " << w.message << " y0=" << y0 << " a=" << fit.a << "\n";
      ++res.n_failed;
      continue;
    }

    res.Q += fit.chi2;
    res.y0_used.push_back(y0);
    res.chi2_per_y0.push_back(fit.chi2);
    res.chi2_ndof_per_y0.push_back(fit.chi2_ndof);
    ++res.n_traces;
  }

  return res;
}

} // namespace riptide
