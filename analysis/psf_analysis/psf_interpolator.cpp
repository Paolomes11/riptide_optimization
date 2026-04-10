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
  Double_t x_source_d, y_source_d;
  Double_t mean_y, mean_z;
  Double_t cov_yy, cov_yz, cov_zz;
  Double_t n_hits_filtered;
  Int_t n_hits_filtered_count = 0;
  Bool_t on_detector;

  tree->SetBranchAddress("config_id", &config_id);
  tree->SetBranchAddress("x1", &x1);
  tree->SetBranchAddress("x2", &x2);
  tree->SetBranchAddress("x_source", &x_source_d);
  tree->SetBranchAddress("y_source", &y_source_d);
  tree->SetBranchAddress("mean_y", &mean_y);
  tree->SetBranchAddress("mean_z", &mean_z);
  tree->SetBranchAddress("cov_yy", &cov_yy);
  tree->SetBranchAddress("cov_yz", &cov_yz);
  tree->SetBranchAddress("cov_zz", &cov_zz);
  tree->SetBranchAddress("n_hits_filtered", &n_hits_filtered);
  tree->SetBranchAddress("on_detector", &on_detector);

  bool has_count_branch = (tree->GetBranch("n_hits_filtered_count") != nullptr);
  if (has_count_branch) {
    tree->SetBranchAddress("n_hits_filtered_count", &n_hits_filtered_count);
  }

  PSFDatabase db;
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    LensConfig cfg{x1, x2};
    db[cfg].push_back(
        {x_source_d, y_source_d, mean_y, mean_z, cov_yy, cov_yz, cov_zz, (bool)on_detector,
         n_hits_filtered,
         has_count_branch ? static_cast<double>(n_hits_filtered_count) : n_hits_filtered});
  }

  // Sort per x crescente, poi y crescente per facilitare l'interpolazione 2D
  for (auto& [cfg, pts] : db) {
    std::sort(pts.begin(), pts.end(), [](const PSFPoint& a, const PSFPoint& b) {
      if (std::abs(a.x_source - b.x_source) > 1e-4)
        return a.x_source < b.x_source;
      return a.y_source < b.y_source;
    });
  }

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

//  Interpolazione 2D (bilineare su griglia x_src, y_src) ──────────────────────────

PSFValue interpolate(double x, double y, const LensConfig& cfg, const PSFDatabase& db) {
  auto it = db.find(cfg);
  if (it == db.end()) {
    std::ostringstream oss;
    oss << "interpolate: config (x1=" << cfg.x1 << ", x2=" << cfg.x2 << ") non trovata";
    throw std::out_of_range(oss.str());
  }

  const auto& pts = it->second;
  if (pts.empty())
    throw std::out_of_range("interpolate: nessun punto PSF per questa configurazione");

  // 1. Estrai valori x distinti presenti nel database
  std::vector<double> x_vals;
  for (const auto& p : pts) {
    if (x_vals.empty() || std::abs(p.x_source - x_vals.back()) > 1e-4)
      x_vals.push_back(p.x_source);
  }

  // 2. Trova l'intervallo in x
  auto it_x = std::lower_bound(x_vals.begin(), x_vals.end(), x);
  double x0, x1, alpha_x;
  if (it_x == x_vals.begin()) {
    x0 = x1 = x_vals.front();
    alpha_x = 0.0;
  } else if (it_x == x_vals.end()) {
    x0 = x1 = x_vals.back();
    alpha_x = 0.0;
  } else {
    x1      = *it_x;
    x0      = *std::prev(it_x);
    alpha_x = (x - x0) / (x1 - x0);
  }

  // Helper per interpolare in y per una x fissata
  auto get_y_interp = [&](double target_x) -> PSFValue {
    std::vector<const PSFPoint*> pts_at_x;
    for (const auto& p : pts) {
      if (std::abs(p.x_source - target_x) < 1e-4)
        pts_at_x.push_back(&p);
    }
    if (pts_at_x.empty())
      throw std::runtime_error("interpolate: nessun dato per x=" + std::to_string(target_x));

    auto it_y = std::lower_bound(pts_at_x.begin(), pts_at_x.end(), y,
                                 [](const PSFPoint* p, double v) { return p->y_source < v; });

    if (it_y == pts_at_x.begin()) {
      const auto& p = **it_y;
      return {p.mu_y,        p.mu_z,   {p.cov_yy, p.cov_yz, p.cov_zz},
              p.on_detector, p.n_hits, p.n_hits_count};
    }
    if (it_y == pts_at_x.end()) {
      const auto& p = **std::prev(it_y);
      return {p.mu_y,        p.mu_z,   {p.cov_yy, p.cov_yz, p.cov_zz},
              p.on_detector, p.n_hits, p.n_hits_count};
    }

    const auto& p1 = **it_y;
    const auto& p0 = **std::prev(it_y);
    double dy      = p1.y_source - p0.y_source;
    double alpha_y = (dy > 1e-12) ? (y - p0.y_source) / dy : 0.0;
    auto lerp_y    = [alpha_y](double a, double b) { return a + alpha_y * (b - a); };

    return {
        lerp_y(p0.mu_y, p1.mu_y),
        lerp_y(p0.mu_z, p1.mu_z),
        {lerp_y(p0.cov_yy, p1.cov_yy), lerp_y(p0.cov_yz, p1.cov_yz), lerp_y(p0.cov_zz, p1.cov_zz)},
        p0.on_detector && p1.on_detector,
        lerp_y(p0.n_hits, p1.n_hits),
        lerp_y(p0.n_hits_count, p1.n_hits_count)};
  };

  // 3. Esegui l'interpolazione bilineare
  PSFValue v0 = get_y_interp(x0);
  if (std::abs(x1 - x0) < 1e-4)
    return v0;

  PSFValue v1 = get_y_interp(x1);

  // Se siamo vicino all'asse (r < r_min), interpola linearmente tra (0,0) e il primo punto
  // NOTA: Assumiamo che se r=0, mu_y=0 e mu_z=0 per simmetria.
  // Troviamo il r_min per questa configurazione
  auto get_r_min = [&](double target_x) {
    for (const auto& p : pts) {
      if (std::abs(p.x_source - target_x) < 1e-4)
        return p.y_source;
    }
    return 0.0;
  };

  double r0_min = get_r_min(x0);
  double r1_min = get_r_min(x1);

  if (y < r0_min && r0_min > 1e-6) {
    double alpha = y / r0_min;
    v0.mu_y *= alpha;
    v0.mu_z = 0.0; // Forza simmetria perfetta lungo l'asse
    // Rendiamo la covarianza isotropa vicino all'asse per evitare artefatti di rotazione
    double sigma_avg = 0.5 * (v0.cov.yy + v0.cov.zz);
    v0.cov.yy        = sigma_avg;
    v0.cov.zz        = sigma_avg;
    v0.cov.yz        = 0.0;
    if (v0.n_hits_count_interp > 0.0)
      v0.on_detector = true;
  }
  if (y < r1_min && r1_min > 1e-6) {
    double alpha = y / r1_min;
    v1.mu_y *= alpha;
    v1.mu_z          = 0.0;
    double sigma_avg = 0.5 * (v1.cov.yy + v1.cov.zz);
    v1.cov.yy        = sigma_avg;
    v1.cov.zz        = sigma_avg;
    v1.cov.yz        = 0.0;
    if (v1.n_hits_count_interp > 0.0)
      v1.on_detector = true;
  }

  auto lerp_x = [alpha_x](double a, double b) { return a + alpha_x * (b - a); };

  return {
      lerp_x(v0.mu_y, v1.mu_y),
      lerp_x(v0.mu_z, v1.mu_z),
      {lerp_x(v0.cov.yy, v1.cov.yy), lerp_x(v0.cov.yz, v1.cov.yz), lerp_x(v0.cov.zz, v1.cov.zz)},
      v0.on_detector && v1.on_detector,
      lerp_x(v0.n_hits_interp, v1.n_hits_interp),
      lerp_x(v0.n_hits_count_interp, v1.n_hits_count_interp)};
}

//  Costruzione traccia 3D ──────────────────────────────────────────────────────

std::vector<TracePoint> build_trace_3d(const Point3D& p1, const Point3D& p2, const LensConfig& cfg,
                                       const PSFDatabase& db, double dt) {
  double dx   = p2.x - p1.x;
  double dy   = p2.y - p1.y;
  double dz   = p2.z - p1.z;
  double L    = std::sqrt(dx * dx + dy * dy + dz * dz);
  const int N = (L > 1e-9) ? static_cast<int>(std::ceil(L / dt)) + 1 : 1;

  std::vector<TracePoint> trace;
  trace.reserve(N);

  for (int i = 0; i < N; ++i) {
    double f = (N > 1) ? static_cast<double>(i) / (N - 1) : 0.0;
    double t = f * L;
    double x = p1.x + f * dx;
    double y = p1.y + f * dy;
    double z = p1.z + f * dz;

    // Simmetria circolare: r = dist dall'asse X, phi = angolo nel piano YZ
    double r   = std::hypot(y, z);
    double phi = std::atan2(z, y);

    // Interpola PSF(x, r) per una sorgente ruotata a phi=0
    PSFValue psf0 = interpolate(x, r, cfg, db);

    // Ruota il centro (mu_y, mu_z) di phi.
    // Per simmetria circolare perfetta, una sorgente a phi=0 produce mu_z=0.
    // Ignoriamo la componente mu_z residua del database per evitare "salti"
    // e asimmetrie numeriche quando la sorgente attraversa l'asse.
    double cos_p = std::cos(phi);
    double sin_p = std::sin(phi);

    double mu_y_rot = psf0.mu_y * cos_p;
    double mu_z_rot = psf0.mu_y * sin_p;

    // Ruota la matrice di covarianza C = R * C0 * R^T
    // C0 = [[yy, yz], [yz, zz]]
    // R = [[cos, -sin], [sin, cos]]
    double c0_yy = psf0.cov.yy;
    double c0_yz = psf0.cov.yz;
    double c0_zz = psf0.cov.zz;

    double c_yy = cos_p * cos_p * c0_yy - 2.0 * cos_p * sin_p * c0_yz + sin_p * sin_p * c0_zz;
    double c_zz = sin_p * sin_p * c0_yy + 2.0 * cos_p * sin_p * c0_yz + cos_p * cos_p * c0_zz;
    double c_yz =
        cos_p * sin_p * c0_yy + (cos_p * cos_p - sin_p * sin_p) * c0_yz - cos_p * sin_p * c0_zz;

    trace.push_back({t,
                     r,
                     x,
                     y,
                     z,
                     mu_y_rot,
                     mu_z_rot,
                     {c_yy, c_yz, c_zz},
                     psf0.on_detector,
                     psf0.n_hits_interp,
                     psf0.n_hits_count_interp});
  }
  return trace;
}

//  Wrapper di compatibilità ───────────────────────────────────────────────────

std::vector<TracePoint> build_trace(double y0, const LensConfig& cfg, const PSFDatabase& db,
                                    double L, double dt) {
  Point3D p1{-L / 2.0, y0, 0.0};
  Point3D p2{L / 2.0, y0, 0.0};
  return build_trace_3d(p1, p2, cfg, db, dt);
}

//  Validità traccia ────────────────────────────────────────────────────────────

bool is_trace_valid(const std::vector<TracePoint>& trace, double point_valid_fraction) {
  if (trace.empty())
    return false;
  int nv = 0;
  for (const auto& pt : trace) {
    if (pt.valid)
      ++nv;
  }
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

LineFitResult fit_trace(const std::vector<TracePoint>& trace, double min_hits_per_point,
                        int max_iter, double tol) {
  // Estrai solo i punti validi (con abbastanza fotoni)
  std::vector<double> vy, vz;
  std::vector<Cov2> vcov;
  vy.reserve(trace.size());
  vz.reserve(trace.size());
  vcov.reserve(trace.size());
  for (const auto& pt : trace) {
    if (pt.valid && pt.n_hits_count >= min_hits_per_point) {
      vy.push_back(pt.mu_y);
      vz.push_back(pt.mu_z);
      vcov.push_back(pt.cov);
    }
  }

  const int N = static_cast<int>(vy.size());
  if (N < 3)
    throw std::invalid_argument("fit_trace: punti validi insufficienti (trovati: "
                                + std::to_string(N) + ", richiesti: 3).");

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

// ─── Temporal unfolding ortogonale ───────────────────────────────────────────

static std::vector<TracePoint> apply_unfolding(const std::vector<TracePoint>& trace,
                                               const QConfig& qcfg) {
  std::vector<TracePoint> unfolded = trace;
  const int N                      = static_cast<int>(unfolded.size());
  if (N < 2)
    return unfolded;

  // 1. Identifica la direzione della traccia (vettore tra primo e ultimo punto valido)
  int i_first = -1, i_last = -1;
  for (int i = 0; i < N; ++i) {
    if (trace[i].valid) {
      if (i_first == -1)
        i_first = i;
      i_last = i;
    }
  }

  // Se non ci sono abbastanza punti validi, non possiamo calcolare una normale
  if (i_first == -1 || i_last == i_first)
    return unfolded;

  double dy_axis = trace[i_last].mu_y - trace[i_first].mu_y;
  double dz_axis = trace[i_last].mu_z - trace[i_first].mu_z;

  // Angolo della traccia sul detector
  double alpha = std::atan2(dz_axis, dy_axis);

  // Direzione ortogonale (normale): alpha + 90 gradi
  double theta = alpha + M_PI / 2.0;
  double nx    = std::cos(theta);
  double ny    = std::sin(theta);

  // 2. Calcola il passo di srotolamento (magnitudo dello spostamento)
  double L  = trace.back().t - trace.front().t;
  double dS = (qcfg.z_unfold_step > 0.0) ? qcfg.z_unfold_step : L / static_cast<double>(N - 1);

  // 3. Applica l'unfolding lungo la normale
  for (int i = 0; i < N; ++i) {
    double shift = static_cast<double>(i) * dS;
    unfolded[i].mu_y += shift * nx;
    unfolded[i].mu_z += shift * ny;
  }

  return unfolded;
}

// ─── compute_Q ───────────────────────────────────────────────────────────────

QResult compute_Q(const LensConfig& cfg, const PSFDatabase& db, const QConfig& qcfg,
                  bool include_non_converged) {
  if (db.find(cfg) == db.end()) {
    throw std::invalid_argument("compute_Q: config non trovata");
  }

  QResult res{};
  res.Q            = 0.0;
  res.n_traces     = 0;
  res.n_failed     = 0;
  res.n_invalid    = 0;
  res.config_valid = true;
  res.chi2_per_trace.reserve(qcfg.n_tracks);
  res.chi2_ndof_per_trace.reserve(qcfg.n_tracks);
  res.trace_valid_flags.reserve(qcfg.n_tracks);

  // Generazione tracce casuali nello scintillatore
  for (int i = 0; i < qcfg.n_tracks; ++i) {
    auto get_random_point = [&]() -> Point3D {
      // Scintillatore centrato in x=0, y=0, z=0
      // Superficie: 6 facce
      double x_hl = qcfg.scint_x / 2.0;
      double y_hl = qcfg.scint_y / 2.0;
      double z_hl = qcfg.scint_z / 2.0;

      int face = std::rand() % 6;
      double u = (static_cast<double>(std::rand()) / RAND_MAX) * 2.0 - 1.0;
      double v = (static_cast<double>(std::rand()) / RAND_MAX) * 2.0 - 1.0;

      if (face == 0)
        return {x_hl, u * y_hl, v * z_hl}; // +X
      if (face == 1)
        return {-x_hl, u * y_hl, v * z_hl}; // -X
      if (face == 2)
        return {u * x_hl, y_hl, v * z_hl}; // +Y
      if (face == 3)
        return {u * x_hl, -y_hl, v * z_hl}; // -Y
      if (face == 4)
        return {u * x_hl, v * y_hl, z_hl}; // +Z
      return {u * x_hl, v * y_hl, -z_hl};  // -Z
    };

    Point3D p1 = get_random_point();
    Point3D p2 = get_random_point();

    // Filtra tracce troppo corte (devono avere almeno 3 punti campionati)
    double dist =
        std::sqrt(std::pow(p2.x - p1.x, 2) + std::pow(p2.y - p1.y, 2) + std::pow(p2.z - p1.z, 2));
    if (dist < 2.5 * qcfg.trace_dt) {
      --i; // Riprova con un'altra coppia
      continue;
    }

    // 1. Traccia 3D
    std::vector<TracePoint> trace;
    try {
      trace = build_trace_3d(p1, p2, cfg, db, qcfg.trace_dt);
    } catch (...) {
      ++res.n_failed;
      res.trace_valid_flags.push_back(false);
      continue;
    }

    // 2. Validità traccia (75% dei punti devono essere on_detector)
    bool trace_ok = is_trace_valid(trace, 0.75);
    res.trace_valid_flags.push_back(trace_ok);
    if (!trace_ok) {
      ++res.n_invalid;
      continue;
    }

    // 3. Unfolding
    std::vector<TracePoint> unfolded;
    const std::vector<TracePoint>* fit_input;
    if (qcfg.apply_temporal_unfolding) {
      unfolded  = apply_unfolding(trace, qcfg);
      fit_input = &unfolded;
    } else {
      fit_input = &trace;
    }

    // 4. Fit
    LineFitResult fit;
    try {
      fit = fit_trace(*fit_input, qcfg.min_hits_per_point, qcfg.fit_max_iter, qcfg.fit_tol);
    } catch (...) {
      ++res.n_failed;
      continue;
    }

    if (!fit.converged && !include_non_converged) {
      ++res.n_failed;
      continue;
    }

    // Controllo validità Chi-squared ridotto
    if (!std::isfinite(fit.chi2_ndof) || fit.chi2_ndof < 0.0) {
      ++res.n_failed;
      continue;
    }

    res.Q += fit.chi2_ndof;
    res.chi2_per_trace.push_back(fit.chi2);
    res.chi2_ndof_per_trace.push_back(fit.chi2_ndof);
    ++res.n_traces;
  }

  // Normalizzazione finale di Q (media del Chi-squared ridotto)
  if (res.n_traces > 0) {
    res.Q /= static_cast<double>(res.n_traces);
    double vf        = static_cast<double>(res.n_traces) / static_cast<double>(qcfg.n_tracks);
    res.config_valid = (vf >= qcfg.trace_valid_fraction);
  } else {
    res.Q            = 0.0;
    res.config_valid = false;
  }

  return res;
}

CoverageResult compute_coverage(const LensConfig& cfg, const PSFDatabase& db, const QConfig& qcfg) {
  if (db.find(cfg) == db.end()) {
    throw std::invalid_argument("compute_coverage: config non trovata");
  }

  CoverageResult res{};
  res.coverage       = 0.0;
  res.n_y0_requested = qcfg.n_tracks;
  res.n_y0_evaluated = 0;
  res.config_valid   = false;

  double sum_frac = 0.0;

  for (int i = 0; i < qcfg.n_tracks; ++i) {
    auto get_random_point = [&]() -> Point3D {
      double x_hl = qcfg.scint_x / 2.0;
      double y_hl = qcfg.scint_y / 2.0;
      double z_hl = qcfg.scint_z / 2.0;
      int face    = std::rand() % 6;
      double u    = (static_cast<double>(std::rand()) / RAND_MAX) * 2.0 - 1.0;
      double v    = (static_cast<double>(std::rand()) / RAND_MAX) * 2.0 - 1.0;
      if (face == 0)
        return {x_hl, u * y_hl, v * z_hl};
      if (face == 1)
        return {-x_hl, u * y_hl, v * z_hl};
      if (face == 2)
        return {u * x_hl, y_hl, v * z_hl};
      if (face == 3)
        return {u * x_hl, -y_hl, v * z_hl};
      if (face == 4)
        return {u * x_hl, v * y_hl, z_hl};
      return {u * x_hl, v * y_hl, -z_hl};
    };

    Point3D p1 = get_random_point();
    Point3D p2 = get_random_point();

    std::vector<TracePoint> trace;
    try {
      trace = build_trace_3d(p1, p2, cfg, db, qcfg.trace_dt);
    } catch (...) {
      continue;
    }

    int nv = 0;
    for (const auto& pt : trace) {
      if (pt.n_hits_count >= qcfg.min_hits_per_point)
        ++nv;
    }
    sum_frac += static_cast<double>(nv) / static_cast<double>(trace.size());
    ++res.n_y0_evaluated;
  }

  if (res.n_y0_evaluated > 0) {
    res.coverage     = sum_frac / res.n_y0_evaluated;
    res.config_valid = true;
  }

  return res;
}

} // namespace riptide
