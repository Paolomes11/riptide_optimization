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

#include <nlohmann/json.hpp>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TTree.h>
#include <omp.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

struct Config {
  double x1, x2;
};

int main() {
  using json = nlohmann::json;

  // --- Leggi configurazione ---
  std::ifstream f("config/config.json");
  if (!f.is_open()) {
    std::cerr << "Impossibile aprire config/config.json\n";
    return 1;
  }
  json config;
  f >> config;

  double x_min = config["x_min"];
  double x_max = config["x_max"];
  double dx    = config["dx"];
  double r1 = config["r1"], h1 = config["h1"];
  double r2 = config["r2"], h2 = config["h2"];

  // --- Apri file ROOT ---
  TFile file("output/events.root", "READ");
  if (!file.IsOpen()) {
    std::cerr << "Impossibile aprire output/events.root\n";
    return 1;
  }
  TTree* tree = (TTree*)file.Get("events");
  if (!tree) {
    std::cerr << "Tree 'events' non trovato nel file\n";
    return 1;
  }

  double x1_val, x2_val;
  int config_id;
  tree->SetBranchAddress("x1", &x1_val);
  tree->SetBranchAddress("x2", &x2_val);
  tree->SetBranchAddress("config_id", &config_id);

  // --- Ricostruisci lista configurazioni e limiti ---
  std::vector<Config> configs;
  double x1_min_all = 1e9, x1_max_all = -1e9;
  double x2_min_all = 1e9, x2_max_all = -1e9;

  for (double x1 = x_min - r1 + h1; x1 <= x_max - h2 - 1 - r1; x1 += dx) {
    if (x1 < x1_min_all)
      x1_min_all = x1;
    if (x1 > x1_max_all)
      x1_max_all = x1;

    double x2_min = x1 + r1 + r2 + 1;
    double x2_max = x_max + r2 - h2;
    if (x2_min < x2_min_all)
      x2_min_all = x2_min;
    if (x2_max > x2_max_all)
      x2_max_all = x2_max;

    for (double x2 = x2_min; x2 <= x2_max; x2 += dx) {
      configs.push_back({x1, x2});
    }
  }

  if (configs.empty()) {
    std::cerr << "Errore: lista configurazioni vuota!\n";
    return 1;
  }

  const double N_generated = 10000.0; // fotoni generati per config

  // --- Prepara vettore contatori ---
  std::vector<int> counts(configs.size(), 0);

  Long64_t nEntries = tree->GetEntries();
  std::cout << "Leggo " << nEntries << " eventi...\n";

// --- Loop sugli eventi, accumula contatori ---
#pragma omp parallel
  {
    std::vector<int> local_counts(configs.size(), 0); // buffer thread-local
#pragma omp for schedule(static)
    for (Long64_t i = 0; i < nEntries; ++i) {
      tree->GetEntry(i);
      if (config_id >= 0 && static_cast<size_t>(config_id) < configs.size()) {
        local_counts[config_id]++;
      }
    }
// Merge thread-local
#pragma omp critical
    for (size_t j = 0; j < counts.size(); ++j)
      counts[j] += local_counts[j];
  }

  // --- Crea istogramma 2D ---
  int bins_x = static_cast<int>(std::ceil((x1_max_all - x1_min_all) / dx));
  int bins_y = static_cast<int>(std::ceil((x2_max_all - x2_min_all) / dx));
  TH2D h_efficiency("h_efficiency", "Geometrical Efficiency;X1 [mm];X2 [mm]", bins_x, x1_min_all,
                    x1_max_all, bins_y, x2_min_all, x2_max_all);

  for (size_t cfg = 0; cfg < configs.size(); ++cfg) {
    double eff = counts[cfg] / N_generated;
    h_efficiency.Fill(configs[cfg].x1, configs[cfg].x2, eff);
  }

  // --- Disegna e salva ---
  TCanvas c("c", "Efficiency", 800, 600);
  gStyle->SetOptStat(0);
  h_efficiency.Draw("COLZ");
  c.SaveAs("output/efficiency2D.png");

  std::cout << "Mappa di efficienza salvata in output/efficiency2D.png\n";
  file.Close();
  return 0;
}
