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
#include <TPaletteAxis.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

struct Config {
  double x1;
  double x2;
};

int main() {
  using json = nlohmann::json;

  // ===============================
  // Lettura configurazione
  // ===============================
  std::ifstream f("config/config.json");
  if (!f.is_open()) {
    std::cerr << "Errore: impossibile aprire config/config.json\n";
    return 1;
  }

  json config;
  f >> config;

  double x_min = config["x_min"];
  double x_max = config["x_max"];
  double dx    = config["dx"];

  double r1 = config["r1"];
  double h1 = config["h1"];
  double r2 = config["r2"];
  double h2 = config["h2"];

  double lower_percentile = config.value("lower_percentile", 0.0);
  double upper_percentile = config.value("upper_percentile", 0.25);

  // ===============================
  // Apertura file ROOT
  // ===============================
  TFile file("output/events.root", "READ");
  if (!file.IsOpen()) {
    std::cerr << "Errore: impossibile aprire output/events.root\n";
    return 1;
  }

  TTree* tree = (TTree*)file.Get("events");
  if (!tree) {
    std::cerr << "Errore: TTree 'events' non trovato\n";
    return 1;
  }

  double x1_val;
  double x2_val;
  int config_id;

  tree->SetBranchAddress("x1", &x1_val);
  tree->SetBranchAddress("x2", &x2_val);
  tree->SetBranchAddress("config_id", &config_id);

  // ===============================
  // Ricostruzione configurazioni
  // ===============================
  std::vector<Config> configs;
  configs.reserve(10000);

  double x1_min_all = 1e9;
  double x1_max_all = -1e9;
  double x2_min_all = 1e9;
  double x2_max_all = -1e9;

  for (double x1 = x_min - r1 + h1; x1 <= x_max - h2 - 1 - r1; x1 += dx) {
    x1_min_all = std::min(x1_min_all, x1);
    x1_max_all = std::max(x1_max_all, x1);

    double x2_min = x1 + r1 + r2 + 1;
    double x2_max = x_max + r2 - h2;

    x2_min_all = std::min(x2_min_all, x2_min);
    x2_max_all = std::max(x2_max_all, x2_max);

    for (double x2 = x2_min; x2 <= x2_max; x2 += dx) {
      configs.push_back({x1, x2});
    }
  }

  if (configs.empty()) {
    std::cerr << "Errore: nessuna configurazione generata\n";
    return 1;
  }
  std::cout << "Numero configurazioni: " << configs.size() << std::endl;

  // ===============================
  // Conteggio eventi
  // ===============================
  const double N_generated = 10000.0;
  std::vector<int> counts(configs.size(), 0);
  Long64_t nEntries = tree->GetEntries();
  std::cout << "Eventi nel file: " << nEntries << std::endl;

  for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    if (config_id >= 0 && (size_t)config_id < configs.size()) {
      counts[config_id]++;
    }
  }

  // ===============================
  // Calcolo efficienze
  // ===============================
  std::vector<double> efficiencies(configs.size());
  for (size_t i = 0; i < configs.size(); ++i) {
    efficiencies[i] = counts[i] / N_generated;
  }

  // ===============================
  // Calcolo percentili
  // ===============================
  std::vector<double> sorted = efficiencies;
  std::sort(sorted.begin(), sorted.end());

  size_t N          = sorted.size();
  size_t low_index  = static_cast<size_t>(lower_percentile * N);
  size_t high_index = static_cast<size_t>((1.0 - upper_percentile) * N);
  if (high_index >= N)
    high_index = N - 1;

  double p_low  = sorted[low_index];
  double p_high = sorted[high_index];

  std::cout << "Percentile basso : " << p_low << std::endl;
  std::cout << "Percentile alto  : " << p_high << std::endl;

  // ===============================
  // Creazione istogramma
  // ===============================
  int bins_x = std::ceil((x1_max_all - x1_min_all) / dx);
  int bins_y = std::ceil((x2_max_all - x2_min_all) / dx);

  TH2D h_efficiency("h_efficiency", "Geometrical Efficiency;X1 [mm];X2 [mm]", bins_x, x1_min_all,
                    x1_max_all, bins_y, x2_min_all, x2_max_all);

  // ===============================
  // Riempimento SOLO dentro percentili
  // ===============================
  double actual_min = 1e9;
  double actual_max = -1e9;

  for (size_t i = 0; i < configs.size(); ++i) {
    double value = efficiencies[i];
    if (value < p_low || value > p_high)
      continue;

    h_efficiency.Fill(configs[i].x1, configs[i].x2, value);

    actual_min = std::min(actual_min, value);
    actual_max = std::max(actual_max, value);
  }

  // Aggiorna scala colori in base ai valori effettivi
  h_efficiency.SetMinimum(actual_min);
  h_efficiency.SetMaximum(actual_max);

  // ===============================
  // Plot con margini, color bar e griglia leggera
  // ===============================
  gStyle->SetOptStat(0);

  // colori della griglia molto leggeri
  gStyle->SetGridColor(kGray + 2); // grigio chiaro
  gStyle->SetGridStyle(2);         // linee tratteggiate
  gStyle->SetGridWidth(1);

  // ===============================
  // Margini e griglia
  // ===============================
  TCanvas canvas("canvas", "Efficiency map", 900, 700);
  canvas.SetLeftMargin(0.12);
  canvas.SetRightMargin(0.18); // margine più moderato
  canvas.SetTopMargin(0.08);
  canvas.SetBottomMargin(0.12);

  canvas.SetGridx();
  canvas.SetGridy();

  // disegna l'istogramma con COLZ
  h_efficiency.Draw("COLZ");

  // forza ROOT a creare la palette
  canvas.Update();

  // ottieni la palette generata da COLZ
  TPaletteAxis* palette = (TPaletteAxis*)h_efficiency.GetListOfFunctions()->FindObject("palette");
  if (palette) {
    palette->SetX1NDC(0.88); // bordo sinistro palette vicino al grafico
    palette->SetX2NDC(0.92); // bordo destro palette
    palette->SetTitle("Efficienza");
    palette->SetTitleOffset(1.0); // titolo appena spostato, non troppo lontano
    palette->SetLabelSize(0.03);
    palette->SetTitleSize(0.035);
  }

  // aggiorna canvas e salva
  canvas.Update();
  canvas.SaveAs("output/efficiency2D.png");

  std::cout << "Mappa salvata in output/efficiency2D.png\n";

  file.Close();
  return 0;
}