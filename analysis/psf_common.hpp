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
 * psf_common.hpp — struct e algoritmo condivisi per la stima PSF
 *
 * Estratto verbatim da psf_extractor.cpp per essere riusato anche da
 * altri tool di analisi (es. psf_gaussianity) senza duplicare l'algoritmo
 * di filtro outlier Mahalanobis iterativo.
 */

#pragma once

#include <cmath>
#include <vector>

// Struttura dati per una singola hit
struct Hit {
  double y, z;
};

// Calcola media e matrice di covarianza con filtro outlier iterativo (3-sigma)
// Restituisce false se non ci sono hit sufficienti
struct PSFResult {
  double mean_y, mean_z;
  double cov_yy, cov_yz, cov_zz;
  double n_hits_raw;      // hit prima del filtro (pesate)
  double n_hits_filtered; // hit dopo il filtro (pesate)
  int n_hits_raw_count;
  int n_hits_filtered_count;
};

inline bool compute_psf(const std::vector<Hit>& raw_hits, PSFResult& result,
                        double sigma_cut = 3.0, int n_iter = 4, double weight = 1.0) {
  if (raw_hits.empty())
    return false;

  result.n_hits_raw       = static_cast<double>(raw_hits.size()) * weight;
  result.n_hits_raw_count = static_cast<int>(raw_hits.size());
  std::vector<Hit> hits   = raw_hits;

  double mean_y = 0, mean_z = 0;

  for (int iter = 0; iter < n_iter; ++iter) {
    size_t n = hits.size();
    if (n < 2)
      return false;

    // MEDIA
    mean_y = 0;
    mean_z = 0;
    for (auto& h : hits) {
      mean_y += h.y;
      mean_z += h.z;
    }
    mean_y /= n;
    mean_z /= n;

    // COVARIANZA
    double var_y = 0, var_z = 0, cov_yz = 0;
    for (auto& h : hits) {
      double dy = h.y - mean_y;
      double dz = h.z - mean_z;
      var_y += dy * dy;
      var_z += dz * dz;
      cov_yz += dy * dz;
    }

    var_y /= n;
    var_z /= n;
    cov_yz /= n;

    // INVERSA COVARIANZA (regolarizzazione Tikhonov: cov += ε·I)
    constexpr double kRegEps = 1e-4;
    double reg_yy = var_y + kRegEps;
    double reg_zz = var_z + kRegEps;
    double det    = reg_yy * reg_zz - cov_yz * cov_yz;

    double inv_yy = reg_zz / det;
    double inv_zz = reg_yy / det;
    double inv_yz = -cov_yz / det;

    // FILTRO ELLITTICO
    std::vector<Hit> filtered;
    filtered.reserve(hits.size());

    for (auto& h : hits) {
      double dy = h.y - mean_y;
      double dz = h.z - mean_z;

      double d2 = dy * dy * inv_yy + dz * dz * inv_zz + 2.0 * dy * dz * inv_yz;

      if (d2 <= sigma_cut * sigma_cut)
        filtered.push_back(h);
    }

    if (filtered.empty())
      break;

    hits.swap(filtered);
  }

  size_t n = hits.size();
  if (n < 3)
    return false;

  // MEDIA FINALE
  mean_y = 0;
  mean_z = 0;
  for (auto& h : hits) {
    mean_y += h.y;
    mean_z += h.z;
  }
  mean_y /= n;
  mean_z /= n;

  // COVARIANZA FINALE (stimatore campionario)
  double cov_yy = 0, cov_yz = 0, cov_zz = 0;
  for (auto& h : hits) {
    double dy = h.y - mean_y;
    double dz = h.z - mean_z;
    cov_yy += dy * dy;
    cov_yz += dy * dz;
    cov_zz += dz * dz;
  }

  double denom = static_cast<double>(n - 1);

  result.mean_y                = mean_y;
  result.mean_z                = mean_z;
  result.cov_yy                = cov_yy / denom;
  result.cov_yz                = cov_yz / denom;
  result.cov_zz                = cov_zz / denom;
  result.n_hits_filtered       = static_cast<double>(n) * weight;
  result.n_hits_filtered_count = static_cast<int>(n);

  return true;
}

// Soglia Bonferroni-corretta sulla coda esatta di chi^2_2 (P(D^2>x) = exp(-x/2)),
// per classificare n_tests ipotesi indipendenti (una per hit) controllando
// l'errore family-wise a alpha_fw. Valida quando media/covarianza sono stimate
// su una frazione h/n prossima a 1 (nessuna correzione finite-sample MCD
// necessaria, cfr. Hardin & Rocke 2005).
inline double bonferroni_chi2_2_cutoff(std::size_t n_tests, double alpha_fw = 0.05) {
  double alpha_per_test = alpha_fw / static_cast<double>(n_tests);
  return -2.0 * std::log(alpha_per_test);
}

// Distanza di Mahalanobis quadrata di una hit rispetto a media/covarianza
// stimate da compute_psf, con la stessa regolarizzazione Tikhonov (ε=1e-4)
// usata internamente al filtro outlier.
inline double mahalanobis_d2(const Hit& h, const PSFResult& r) {
  constexpr double kRegEps = 1e-4;
  double reg_yy = r.cov_yy + kRegEps;
  double reg_zz = r.cov_zz + kRegEps;
  double det    = reg_yy * reg_zz - r.cov_yz * r.cov_yz;

  double inv_yy = reg_zz / det;
  double inv_zz = reg_yy / det;
  double inv_yz = -r.cov_yz / det;

  double dy = h.y - r.mean_y;
  double dz = h.z - r.mean_z;

  return dy * dy * inv_yy + dz * dz * inv_zz + 2.0 * dy * dz * inv_yz;
}
