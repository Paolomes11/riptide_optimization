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

#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

struct Hit {
  double y;
  double z;
};

int main(int argc, char** argv) {
  if (argc != 4) {
    std::cerr << "Uso: " << argv[0] << " x1_target x2_target y0_target\n";
    return 1;
  }

  double x1_target = std::stod(argv[1]);
  double x2_target = std::stod(argv[2]);
  double y0_target = std::stod(argv[3]);

  // apri file ROOT
  TFile* file = TFile::Open("output/lens_simulation/lens_20260319_170944.root");
  if (!file || file->IsZombie()) {
    std::cerr << "Errore nell'aprire il file ROOT\n";
    return 1;
  }

  TTree* tConfig = (TTree*)file->Get("Configurations");
  TTree* tRuns   = (TTree*)file->Get("Runs");
  TTree* tHits   = (TTree*)file->Get("Hits");
  if (!tConfig || !tRuns || !tHits) {
    std::cerr << "Impossibile leggere uno dei TTree\n";
    return 1;
  }

  // setup branch
  int config_id;
  double x1, x2;
  tConfig->SetBranchAddress("config_id", &config_id);
  tConfig->SetBranchAddress("x1", &x1);
  tConfig->SetBranchAddress("x2", &x2);

  int run_id, run_config, n_hits_run;
  float y_source;
  tRuns->SetBranchAddress("run_id", &run_id);
  tRuns->SetBranchAddress("config_id", &run_config);
  tRuns->SetBranchAddress("x_source", &y_source);
  tRuns->SetBranchAddress("n_hits", &n_hits_run);

  float y_hit, z_hit;
  tHits->SetBranchAddress("y_hit", &y_hit);
  tHits->SetBranchAddress("z_hit", &z_hit);

  // trova configurazione più vicina
  int target_config_id = -1;
  double best_dist     = std::numeric_limits<double>::max();
  double x1_sel = 0, x2_sel = 0;

  for (Long64_t i = 0; i < tConfig->GetEntries(); ++i) {
    tConfig->GetEntry(i);
    double dist =
        std::sqrt((x1 - x1_target) * (x1 - x1_target) + (x2 - x2_target) * (x2 - x2_target));
    if (dist < best_dist) {
      best_dist        = dist;
      target_config_id = config_id;
      x1_sel           = x1;
      x2_sel           = x2;
    }
  }

  if (target_config_id < 0) {
    std::cerr << "Nessuna configurazione trovata!\n";
    return 1;
  }
  std::cout << "Configurazione scelta: config_id=" << target_config_id << ", x1=" << x1_sel
            << ", x2=" << x2_sel << "\n";

  // trova run esatto
  struct RunInfo {
    Long64_t first_hit;
    int n_hits;
  };
  std::map<int, RunInfo> run_info_map;
  Long64_t cumulative_hits = 0;
  for (Long64_t i = 0; i < tRuns->GetEntries(); ++i) {
    tRuns->GetEntry(i);
    run_info_map[run_id] = {cumulative_hits, n_hits_run};
    cumulative_hits += n_hits_run;
  }

  int exact_run_id = -1;
  for (Long64_t i = 0; i < tRuns->GetEntries(); ++i) {
    tRuns->GetEntry(i);
    if (run_config == target_config_id && std::abs(y_source - y0_target) < 1e-3) {
      exact_run_id = run_id;
      break;
    }
  }

  if (exact_run_id < 0) {
    std::cerr << "Nessun run trovato con y0_target = " << y0_target << "\n";
    return 1;
  }
  std::cout << "Run selezionato: run_id=" << exact_run_id << ", y_source=" << y0_target << "\n";

  // raccogli tutte le hit
  std::vector<Hit> hits;
  // ri = run_info
  const auto& ri = run_info_map.at(exact_run_id);
  hits.reserve(ri.n_hits);
  for (Long64_t i = ri.first_hit; i < ri.first_hit + ri.n_hits; ++i) {
    tHits->GetEntry(i);
    hits.push_back({static_cast<double>(y_hit), static_cast<double>(z_hit)});
  }

  if (hits.empty()) {
    std::cerr << "Nessuna hit trovata!\n";
    return 1;
  }

  // calcolo robusto media e covarianza con taglio ellittico
  const double max_sigma         = 2.0; // soglia (raggio ellisse)
  std::vector<Hit> hits_filtered = hits;

  double y_mean = 0, z_mean = 0;
  double y_sigma = 0, z_sigma = 0;

  // più iterazioni → più stabile
  for (int iter = 0; iter < 4; ++iter) {
    size_t n = hits_filtered.size();
    if (n < 2)
      break;

    // media
    double y_sum = 0, z_sum = 0;
    for (auto& h : hits_filtered) {
      y_sum += h.y;
      z_sum += h.z;
    }
    y_mean = y_sum / n;
    z_mean = z_sum / n;

    // covarianza
    double y_var = 0, z_var = 0, yz_cov = 0;
    for (auto& h : hits_filtered) {
      double dy = h.y - y_mean;
      double dz = h.z - z_mean;
      y_var += dy * dy;
      z_var += dz * dz;
      yz_cov += dy * dz;
    }

    y_var /= n;
    z_var /= n;
    yz_cov /= n;

    y_sigma = std::sqrt(y_var);
    z_sigma = std::sqrt(z_var);

    // inversa della matrice di covarianza
    double det = y_var * z_var - yz_cov * yz_cov;
    if (det < 1e-12)
      det = 1e-12;

    double inv_yy = z_var / det;
    double inv_zz = y_var / det;
    double inv_yz = -yz_cov / det;

    // filtro ellittico (Mahalanobis)
    std::vector<Hit> temp;
    temp.reserve(hits_filtered.size());

    for (auto& h : hits_filtered) {
      double dy = h.y - y_mean;
      double dz = h.z - z_mean;

      double d2 = dy * dy * inv_yy + dz * dz * inv_zz + 2.0 * dy * dz * inv_yz;

      if (d2 <= max_sigma * max_sigma) {
        temp.push_back(h);
      }
    }

    // evita di svuotare tutto per instabilità numeriche
    if (temp.empty())
      break;

    hits_filtered.swap(temp);
  }

  std::cout << "Hit dopo filtraggio: " << hits_filtered.size() << ", media y=" << y_mean
            << ", sigma y=" << y_sigma << ", media z=" << z_mean << ", sigma z=" << z_sigma << "\n";

  // debug migliorato
  std::cout << "Prime 10 hit raccolte (y, z):\n";
  for (size_t i = 0; i < std::min<size_t>(10, hits.size()); ++i) {
    std::cout << "  (" << hits[i].y << ", " << hits[i].z << ")\n";
  }

  // canvas e istogramma
  double y_min = y_mean - max_sigma * y_sigma;
  double y_max = y_mean + max_sigma * y_sigma;
  double z_min = z_mean - max_sigma * z_sigma;
  double z_max = z_mean + max_sigma * z_sigma;

  TCanvas* c = new TCanvas("c", "Hits sul detector", 1200, 1000);
  gStyle->SetOptStat(0);
  c->SetGrid();

  TH2D* hist = new TH2D("hist",
                        ("Hits: x1=" + std::to_string(x1_sel) + ", x2=" + std::to_string(x2_sel)
                         + ", y0=" + std::to_string(y0_target))
                            .c_str(),
                        200, y_min, y_max, 200, z_min, z_max);

  for (auto& h : hits_filtered)
    hist->Fill(h.y, h.z);

  hist->GetXaxis()->SetTitle("Y [mm]");
  hist->GetYaxis()->SetTitle("Z [mm]");
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetTitleOffset(1.2);
  hist->Draw("COLZ");

  std::ostringstream filename;
  filename << "output/mean_covariance_maps/detector_hits_config_" << target_config_id << "_y0_"
           << y0_target << ".png";
  c->SaveAs(filename.str().c_str());
  std::cout << "Immagine salvata in: " << filename.str() << "\n";

  return 0;
}