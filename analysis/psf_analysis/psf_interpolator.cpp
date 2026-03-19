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

} // namespace riptide
