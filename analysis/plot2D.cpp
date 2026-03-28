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
#include <map>
#include <vector>

struct Config {
  double x1;
  double x2;
};

struct CliConfig {
  std::string input_file  = "output/events.root";
  std::string config_file = "config/config.json";
  std::string output_png  = "output/efficiency2D.png";
  double lower_percentile = -1.0; // -1 significa "usa dal config"
  double upper_percentile = -1.0;
};

static CliConfig parse_args(int argc, char** argv) {
  CliConfig cfg;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    auto next       = [&]() -> std::string {
      if (i + 1 >= argc) {
        std::cerr << "Argomento mancante dopo " << arg << "\n";
        std::exit(1);
      }
      return argv[++i];
    };

    if (arg == "--input" || arg == "-i")
      cfg.input_file = next();
    else if (arg == "--config" || arg == "-c")
      cfg.config_file = next();
    else if (arg == "--output" || arg == "-o")
      cfg.output_png = next();
    else if (arg == "--low")
      cfg.lower_percentile = std::stod(next());
    else if (arg == "--high")
      cfg.upper_percentile = std::stod(next());
    else if (arg == "--help" || arg == "-h") {
      std::cout
          << "Uso: plot2D [opzioni]\n"
          << "Opzioni:\n"
          << "  -i, --input <file>    File ROOT di input (default: output/events.root)\n"
          << "  -c, --config <file>   File JSON di configurazione (default: config/config.json)\n"
          << "  -o, --output <file>   Nome del file PNG in uscita (default: "
             "output/efficiency2D.png)\n"
          << "  --low <val>           Percentile inferiore (sovrascrive config)\n"
          << "  --high <val>          Percentile superiore (sovrascrive config)\n";
      std::exit(0);
    } else {
      std::cerr << "Opzione sconosciuta: " << arg << "\n";
      std::exit(1);
    }
  }
  return cfg;
}

int main(int argc, char** argv) {
  using json = nlohmann::json;

  CliConfig cli = parse_args(argc, argv);

  // ===============================
  // Lettura configurazione
  // ===============================
  std::ifstream f(cli.config_file);
  if (!f.is_open()) {
    std::cerr << "Errore: impossibile aprire " << cli.config_file << "\n";
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

  double lower_percentile =
      (cli.lower_percentile >= 0) ? cli.lower_percentile : config.value("lower_percentile", 0.0);
  double upper_percentile =
      (cli.upper_percentile >= 0) ? cli.upper_percentile : config.value("upper_percentile", 0.0);

  // ===============================
  // Apertura file ROOT
  // ===============================
  TFile file(cli.input_file.c_str(), "READ");
  if (!file.IsOpen()) {
    std::cerr << "Errore: impossibile aprire " << cli.input_file << "\n";
    return 1;
  }

  TTree* tree_config = (TTree*)file.Get("Configurations");
  TTree* tree_eff    = (TTree*)file.Get("Efficiency");
  if (!tree_config || !tree_eff) {
    std::cerr << "Errore: TTree 'Configurations' or 'Efficiency' non trovato\n";
    return 1;
  }

  double x1_val;
  double x2_val;
  int config_id;
  int n_photons;
  int n_hits;

  tree_config->SetBranchAddress("x1", &x1_val);
  tree_config->SetBranchAddress("x2", &x2_val);
  tree_config->SetBranchAddress("config_id", &config_id);

  tree_eff->SetBranchAddress("config_id", &config_id);
  tree_eff->SetBranchAddress("n_photons", &n_photons);
  tree_eff->SetBranchAddress("n_hits", &n_hits);

  // Mappa per associare config_id a (x1, x2)
  std::map<int, Config> config_map;
  for (int i = 0; i < tree_config->GetEntries(); ++i) {
    tree_config->GetEntry(i);
    config_map[config_id] = {x1_val, x2_val};
  }

  // ===============================
  // Ricostruzione configurazioni
  // ===============================
  double x1_min_all = 1e9;
  double x1_max_all = -1e9;
  double x2_min_all = 1e9;
  double x2_max_all = -1e9;

  for (auto const& [id, cfg] : config_map) {
    x1_min_all = std::min(x1_min_all, cfg.x1);
    x1_max_all = std::max(x1_max_all, cfg.x1);
    x2_min_all = std::min(x2_min_all, cfg.x2);
    x2_max_all = std::max(x2_max_all, cfg.x2);
  }

  // Creazione istogramma
  int n_bins_x1 = std::round((x1_max_all - x1_min_all) / dx) + 1;
  int n_bins_x2 = std::round((x2_max_all - x2_min_all) / dx) + 1;

  TH2D h_efficiency("h_efficiency", "Geometrical Efficiency;X1 [mm];X2 [mm]", n_bins_x1,
                    x1_min_all - dx / 2.0, x1_max_all + dx / 2.0, n_bins_x2, x2_min_all - dx / 2.0,
                    x2_max_all + dx / 2.0);

  // Riempimento istogramma
  std::vector<double> effs;
  for (int i = 0; i < tree_eff->GetEntries(); ++i) {
    tree_eff->GetEntry(i);
    effs.push_back((double)n_hits / n_photons);
  }

  std::vector<double> sorted = effs;
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

  double actual_min = 1e9;
  double actual_max = -1e9;
  for (int i = 0; i < tree_eff->GetEntries(); ++i) {
    tree_eff->GetEntry(i);
    double eff = (double)n_hits / n_photons;
    if (eff < p_low || eff > p_high)
      continue;

    if (config_map.count(config_id)) {
      h_efficiency.Fill(config_map[config_id].x1, config_map[config_id].x2, eff);
      actual_min = std::min(actual_min, eff);
      actual_max = std::max(actual_max, eff);
    }
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
  canvas.SaveAs(cli.output_png.c_str());

  std::cout << "Mappa salvata in " << cli.output_png << "\n";

  file.Close();
  return 0;
}