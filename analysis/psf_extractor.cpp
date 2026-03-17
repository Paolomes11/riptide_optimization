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
 * psf_extractor — estrae media e matrice di covarianza della PSF
 *
 * Per ogni combinazione (config_id, y_source) presente nel file lens.root,
 * calcola:
 *   - media y e z (con filtro outlier 2-sigma, 2 iterazioni)
 *   - matrice di covarianza 2x2: cov_yy, cov_yz, cov_zz
 *   - numero di hit dopo filtraggio
 *
 * Output: psf_data.root con un TTree "PSF" contenente una riga per run.
 *
 * Uso:
 *   ./psf_extractor [input.root] [output.root]
 *   default input:  output/lens_simulation/lens.root
 *   default output: output/psf/psf_data.root
 */

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// Struttura dati per una singola hit
struct Hit {
  double y, z;
};

// Calcola media e matrice di covarianza con filtro outlier iterativo (2-sigma)
// Restituisce false se non ci sono hit sufficienti
struct PSFResult {
  double mean_y, mean_z;
  double cov_yy, cov_yz, cov_zz;
  int n_hits_raw;      // hit prima del filtro
  int n_hits_filtered; // hit dopo il filtro
};

bool compute_psf(const std::vector<Hit>& raw_hits, PSFResult& result, double sigma_cut = 2.0,
                 int n_iter = 2) {
  if (raw_hits.empty())
    return false;

  result.n_hits_raw     = static_cast<int>(raw_hits.size());
  std::vector<Hit> hits = raw_hits;

  double mean_y = 0, mean_z = 0;
  double sig_y = 0, sig_z = 0;

  for (int iter = 0; iter < n_iter; ++iter) {
    // Calcola medie
    mean_y = 0;
    mean_z = 0;
    for (auto& h : hits) {
      mean_y += h.y;
      mean_z += h.z;
    }
    size_t n = hits.size();
    if (n == 0)
      return false;
    mean_y /= n;
    mean_z /= n;

    // Calcola deviazioni standard
    double var_y = 0, var_z = 0;
    for (auto& h : hits) {
      var_y += (h.y - mean_y) * (h.y - mean_y);
      var_z += (h.z - mean_z) * (h.z - mean_z);
    }
    sig_y = std::sqrt(var_y / n);
    sig_z = std::sqrt(var_z / n);

    // Filtra outlier
    std::vector<Hit> filtered;
    filtered.reserve(hits.size());
    for (auto& h : hits) {
      if (std::abs(h.y - mean_y) <= sigma_cut * sig_y
          && std::abs(h.z - mean_z) <= sigma_cut * sig_z)
        filtered.push_back(h);
    }
    hits = std::move(filtered);
  }

  size_t n = hits.size();
  if (n < 3)
    return false; // minimo per covarianza sensata

  // Ricalcola media finale sulle hit filtrate
  mean_y = 0;
  mean_z = 0;
  for (auto& h : hits) {
    mean_y += h.y;
    mean_z += h.z;
  }
  mean_y /= n;
  mean_z /= n;

  // Matrice di covarianza (estimatore campionario, divisione per n-1)
  double cov_yy = 0, cov_yz = 0, cov_zz = 0;
  for (auto& h : hits) {
    double dy = h.y - mean_y;
    double dz = h.z - mean_z;
    cov_yy += dy * dy;
    cov_yz += dy * dz;
    cov_zz += dz * dz;
  }
  double denom           = static_cast<double>(n - 1);
  result.mean_y          = mean_y;
  result.mean_z          = mean_z;
  result.cov_yy          = cov_yy / denom;
  result.cov_yz          = cov_yz / denom;
  result.cov_zz          = cov_zz / denom;
  result.n_hits_filtered = static_cast<int>(n);
  return true;
}

// Main
int main(int argc, char** argv) {
  std::string input_path  = "output/lens_simulation/lens.root";
  std::string output_path = "output/psf/psf_data.root";

  if (argc > 1)
    input_path = argv[1];
  if (argc > 2)
    output_path = argv[2];

  // Apri file di input
  TFile* fin = TFile::Open(input_path.c_str(), "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Errore: impossibile aprire " << input_path << "\n";
    return 1;
  }

  TTree* tConfig = dynamic_cast<TTree*>(fin->Get("Configurations"));
  TTree* tRuns   = dynamic_cast<TTree*>(fin->Get("Runs"));
  TTree* tHits   = dynamic_cast<TTree*>(fin->Get("Hits"));

  if (!tConfig || !tRuns || !tHits) {
    std::cerr << "Errore: TTree mancanti nel file di input\n";
    return 1;
  }

  // Leggi Configurations
  int cfg_id;
  double cfg_x1, cfg_x2;
  tConfig->SetBranchAddress("config_id", &cfg_id);
  tConfig->SetBranchAddress("x1", &cfg_x1);
  tConfig->SetBranchAddress("x2", &cfg_x2);

  // Mappa config_id -> (x1, x2)
  std::map<int, std::pair<double, double>> config_map;
  for (Long64_t i = 0; i < tConfig->GetEntries(); ++i) {
    tConfig->GetEntry(i);
    config_map[cfg_id] = {cfg_x1, cfg_x2};
  }
  std::cout << "Configurazioni lette: " << config_map.size() << "\n";

  // Leggi Runs e costruisci mappa run_id -> offset hits
  int run_id_r, run_cfg_id, n_hits_r;
  float y_source_r;
  tRuns->SetBranchAddress("run_id", &run_id_r);
  tRuns->SetBranchAddress("config_id", &run_cfg_id);
  tRuns->SetBranchAddress("x_source", &y_source_r);
  tRuns->SetBranchAddress("n_hits", &n_hits_r);

  struct RunInfo {
    int config_id;
    float y_source;
    Long64_t first_hit;
    int n_hits;
  };

  std::vector<RunInfo> runs;
  runs.reserve(static_cast<size_t>(tRuns->GetEntries()));
  Long64_t cumulative = 0;
  for (Long64_t i = 0; i < tRuns->GetEntries(); ++i) {
    tRuns->GetEntry(i);
    runs.push_back({run_cfg_id, y_source_r, cumulative, n_hits_r});
    cumulative += n_hits_r;
  }
  std::cout << "Run letti: " << runs.size() << "\n";
  std::cout << "Hit totali: " << cumulative << "\n";

  // Setup lettura hits
  float y_hit_r, z_hit_r;
  tHits->SetBranchAddress("y_hit", &y_hit_r);
  tHits->SetBranchAddress("z_hit", &z_hit_r);

  // Apri file di output
  TFile* fout = TFile::Open(output_path.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "Errore: impossibile creare " << output_path << "\n";
    return 1;
  }

  // Crea TTree di output
  TTree* tPSF = new TTree("PSF", "PSF mean and covariance per run");

  int out_config_id;
  double out_x1, out_x2;
  float out_y_source;
  double out_mean_y, out_mean_z;
  double out_cov_yy, out_cov_yz, out_cov_zz;
  int out_n_hits_raw, out_n_hits_filtered;

  tPSF->Branch("config_id", &out_config_id, "config_id/I");
  tPSF->Branch("x1", &out_x1, "x1/D");
  tPSF->Branch("x2", &out_x2, "x2/D");
  tPSF->Branch("y_source", &out_y_source, "y_source/F");
  tPSF->Branch("mean_y", &out_mean_y, "mean_y/D");
  tPSF->Branch("mean_z", &out_mean_z, "mean_z/D");
  tPSF->Branch("cov_yy", &out_cov_yy, "cov_yy/D");
  tPSF->Branch("cov_yz", &out_cov_yz, "cov_yz/D");
  tPSF->Branch("cov_zz", &out_cov_zz, "cov_zz/D");
  tPSF->Branch("n_hits_raw", &out_n_hits_raw, "n_hits_raw/I");
  tPSF->Branch("n_hits_filtered", &out_n_hits_filtered, "n_hits_filtered/I");

  // Loop su tutti i run
  int skipped = 0;
  std::vector<Hit> hits_buf;
  hits_buf.reserve(1000);

  for (size_t r = 0; r < runs.size(); ++r) {
    const auto& run = runs[r];

    // Trova x1, x2 per questa config
    auto it = config_map.find(run.config_id);
    if (it == config_map.end()) {
      ++skipped;
      continue;
    }

    // Carica le hit del run
    hits_buf.clear();
    for (Long64_t j = run.first_hit; j < run.first_hit + run.n_hits; ++j) {
      tHits->GetEntry(j);
      hits_buf.push_back({static_cast<double>(y_hit_r), static_cast<double>(z_hit_r)});
    }

    // Calcola PSF
    PSFResult res{};
    if (!compute_psf(hits_buf, res)) {
      ++skipped;
      continue;
    }

    // Riempi rami output
    out_config_id       = run.config_id;
    out_x1              = it->second.first;
    out_x2              = it->second.second;
    out_y_source        = run.y_source;
    out_mean_y          = res.mean_y;
    out_mean_z          = res.mean_z;
    out_cov_yy          = res.cov_yy;
    out_cov_yz          = res.cov_yz;
    out_cov_zz          = res.cov_zz;
    out_n_hits_raw      = res.n_hits_raw;
    out_n_hits_filtered = res.n_hits_filtered;
    tPSF->Fill();

    // Progresso ogni 5000 run
    if ((r + 1) % 5000 == 0)
      std::cout << "  Processati " << r + 1 << " / " << runs.size() << " run...\n";
  }

  fout->cd();
  tPSF->Write();
  fout->Close();
  fin->Close();

  std::cout << "\nRun processati: " << runs.size() - skipped << "\n";
  std::cout << "Run saltati (hit insufficienti): " << skipped << "\n";
  std::cout << "Output salvato in: " << output_path << "\n";

  return 0;
}
