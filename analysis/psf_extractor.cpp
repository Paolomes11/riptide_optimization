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
 *   - media y e z (con filtro outlier 3-sigma, 2 iterazioni)
 *   - matrice di covarianza 2x2: cov_yy, cov_yz, cov_zz
 *   - numero di hit dopo filtraggio
 *   - flag on_detector: true se n_hits_filtered >= min_hits (default 150)
 *
 * Il flag on_detector marca i run in cui i fotoni hanno effettivamente
 * raggiunto il fotocatodo in numero sufficiente per una stima statistica
 * affidabile. Run con pochi hit corrispondono a configurazioni in cui la
 * PSF cade parzialmente o totalmente fuori dal fotocatodo.
 *
 * Uso:
 *   ./psf_extractor [input.root] [output.root] [min_hits]
 *   default input:    output/lens_simulation/lens.root
 *   default output:   output/psf/psf_data.root
 *   default min_hits: 150
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

// Calcola media e matrice di covarianza con filtro outlier iterativo (3-sigma)
// Restituisce false se non ci sono hit sufficienti
struct PSFResult {
  double mean_y, mean_z;
  double cov_yy, cov_yz, cov_zz;
  int n_hits_raw;      // hit prima del filtro
  int n_hits_filtered; // hit dopo il filtro
};

bool compute_psf(const std::vector<Hit>& raw_hits, PSFResult& result, double sigma_cut = 3.0,
                 int n_iter = 4) {
  if (raw_hits.empty())
    return false;

  result.n_hits_raw     = static_cast<int>(raw_hits.size());
  std::vector<Hit> hits = raw_hits;

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

    // INVERSA COVARIANZA
    double det = var_y * var_z - cov_yz * cov_yz;
    if (det < 1e-12)
      det = 1e-12;

    double inv_yy = var_z / det;
    double inv_zz = var_y / det;
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
  int min_hits            = 150;

  if (argc > 1)
    input_path = argv[1];
  if (argc > 2)
    output_path = argv[2];
  if (argc > 3)
    min_hits = std::stoi(argv[3]);

  std::cout << "Soglia on_detector: n_hits_filtered >= " << min_hits << "\n";

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
  bool out_on_detector;

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
  tPSF->Branch("on_detector", &out_on_detector, "on_detector/O");

  // Loop su tutti i run
  int skipped     = 0;
  int flagged_off = 0;
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

    // Caso run con zero hit: salva la riga comunque con on_detector = false
    // e statistiche nulle, così il database è completo per y_source
    if (hits_buf.empty()) {
      out_config_id       = run.config_id;
      out_x1              = it->second.first;
      out_x2              = it->second.second;
      out_y_source        = run.y_source;
      out_mean_y          = 0.0;
      out_mean_z          = 0.0;
      out_cov_yy          = 0.0;
      out_cov_yz          = 0.0;
      out_cov_zz          = 0.0;
      out_n_hits_raw      = 0;
      out_n_hits_filtered = 0;
      out_on_detector     = false;
      tPSF->Fill();
      ++flagged_off;
      continue;
    }

    // Calcola PSF
    PSFResult res{};
    bool ok = compute_psf(hits_buf, res);

    if (!ok) {
      // Meno di 3 hit dopo il filtro: salva con on_detector = false
      out_config_id       = run.config_id;
      out_x1              = it->second.first;
      out_x2              = it->second.second;
      out_y_source        = run.y_source;
      out_mean_y          = 0.0;
      out_mean_z          = 0.0;
      out_cov_yy          = 0.0;
      out_cov_yz          = 0.0;
      out_cov_zz          = 0.0;
      out_n_hits_raw      = static_cast<int>(hits_buf.size());
      out_n_hits_filtered = 0;
      out_on_detector     = false;
      tPSF->Fill();
      ++flagged_off;
      continue;
    }

    // Imposta il flag: on_detector solo se ci sono abbastanza hit
    bool on_det = (res.n_hits_filtered >= min_hits);

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
    out_on_detector     = on_det;
    tPSF->Fill();

    if (!on_det)
      ++flagged_off;

    // Progresso ogni 5000 run
    if ((r + 1) % 5000 == 0)
      std::cout << "  Processati " << r + 1 << " / " << runs.size() << " run...\n";
  }

  fout->cd();
  tPSF->Write();
  fout->Close();
  fin->Close();

  std::cout << "\nRun processati:                        " << runs.size() - skipped << "\n";
  std::cout << "Run saltati (config non trovata):      " << skipped << "\n";
  std::cout << "Run con on_detector = false:           " << flagged_off << "\n";
  std::cout << "Output salvato in: " << output_path << "\n";

  return 0;
}
