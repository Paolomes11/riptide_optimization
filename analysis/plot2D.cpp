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
  std::string lens1_id;
  std::string lens2_id;
};

struct CliConfig {
  std::string input_file  = "output/events.root";
  std::string config_file = "config/config.json";
  std::string output_png  = "output/efficiency2D.png";
  std::string lens1_id    = ""; // Se vuoto, usa la prima trovata
  std::string lens2_id    = "";
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
    else if (arg == "--lens1")
      cfg.lens1_id = next();
    else if (arg == "--lens2")
      cfg.lens2_id = next();
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
          << "  --lens1 <id>          ID della prima lente da plottare (default: prima trovata)\n"
          << "  --lens2 <id>          ID della seconda lente da plottare (default: prima trovata)\n"
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

  double x_min = config.value("x_min", 0.0);
  double x_max = config.value("x_max", 200.0);
  double dx    = config.value("dx", 1.0);

  double r1 = config.value("r1", 0.0);
  double h1 = config.value("h1", 0.0);
  double r2 = config.value("r2", 0.0);
  double h2 = config.value("h2", 0.0);

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
  bool is_opt_mode   = true;

  if (!tree_eff) {
    tree_eff    = (TTree*)file.Get("Runs");
    is_opt_mode = false;
  }

  if (!tree_config || !tree_eff) {
    std::cerr << "Errore: TTree 'Configurations' or 'Efficiency/Runs' non trovato\n";
    return 1;
  }

  double x1_val;
  double x2_val;
  int config_id;
  int n_photons_val;
  int n_hits_val;
  char l1_id_buf[256];
  char l2_id_buf[256];

  tree_config->SetBranchAddress("x1", &x1_val);
  tree_config->SetBranchAddress("x2", &x2_val);
  tree_config->SetBranchAddress("config_id", &config_id);
  if (tree_config->GetBranch("lens75_id")) {
    tree_config->SetBranchAddress("lens75_id", l1_id_buf);
    tree_config->SetBranchAddress("lens60_id", l2_id_buf);
  }

  if (is_opt_mode) {
    tree_eff->SetBranchAddress("config_id", &config_id);
    tree_eff->SetBranchAddress("n_photons", &n_photons_val);
    tree_eff->SetBranchAddress("n_hits", &n_hits_val);
  } else {
    // Lens simulation mode
    tree_eff->SetBranchAddress("config_id", &config_id);
    tree_eff->SetBranchAddress("n_hits", &n_hits_val);
    n_photons_val = config.value("n_photons", 10000);
  }

  // Mappa per associare config_id a Config
  std::map<int, Config> config_map;
  std::vector<std::pair<std::string, std::string>> available_pairs;

  for (int i = 0; i < tree_config->GetEntries(); ++i) {
    tree_config->GetEntry(i);
    std::string s1 = tree_config->GetBranch("lens75_id") ? l1_id_buf : "default";
    std::string s2 = tree_config->GetBranch("lens60_id") ? l2_id_buf : "default";
    config_map[config_id] = {x1_val, x2_val, s1, s2};

    bool found = false;
    for (auto const& p : available_pairs) {
      if (p.first == s1 && p.second == s2) {
        found = true;
        break;
      }
    }
    if (!found)
      available_pairs.push_back({s1, s2});
  }

  // Selezione lente
  std::string sel_l1 = cli.lens1_id;
  std::string sel_l2 = cli.lens2_id;

  if (sel_l1.empty() || sel_l2.empty()) {
    if (!available_pairs.empty()) {
      sel_l1 = available_pairs[0].first;
      sel_l2 = available_pairs[0].second;
    }
  }

  std::cout << "Plotting pair: " << sel_l1 << " & " << sel_l2 << std::endl;

  // Filtra config_map per la coppia selezionata
  std::map<int, Config> filtered_config;
  for (auto const& [id, cfg] : config_map) {
    if (cfg.lens1_id == sel_l1 && cfg.lens2_id == sel_l2) {
      filtered_config[id] = cfg;
    }
  }

  if (filtered_config.empty()) {
    std::cerr << "Nessuna configurazione trovata per " << sel_l1 << " & " << sel_l2 << std::endl;
    return 1;
  }

  // ===============================
  // Ricostruzione configurazioni
  // ===============================
  double x1_min_all = 1e9;
  double x1_max_all = -1e9;
  double x2_min_all = 1e9;
  double x2_max_all = -1e9;

  for (auto const& [id, cfg] : filtered_config) {
    x1_min_all = std::min(x1_min_all, cfg.x1);
    x1_max_all = std::max(x1_max_all, cfg.x1);
    x2_min_all = std::min(x2_min_all, cfg.x2);
    x2_max_all = std::max(x2_max_all, cfg.x2);
  }

  // Creazione istogramma
  int n_bins_x1 = std::round((x1_max_all - x1_min_all) / dx) + 1;
  int n_bins_x2 = std::round((x2_max_all - x2_min_all) / dx) + 1;

  TH2D h_efficiency("h_efficiency",
                    Form("Efficiency: %s + %s;X1 [mm];X2 [mm]", sel_l1.c_str(), sel_l2.c_str()),
                    n_bins_x1, x1_min_all - dx / 2.0, x1_max_all + dx / 2.0, n_bins_x2,
                    x2_min_all - dx / 2.0, x2_max_all + dx / 2.0);

  // Riempimento istogramma
  std::map<int, double> efficiency_sum;
  std::map<int, int> n_entries;

  for (int i = 0; i < tree_eff->GetEntries(); ++i) {
    tree_eff->GetEntry(i);
    if (filtered_config.count(config_id)) {
      double eff = (double)n_hits_val / n_photons_val;
      efficiency_sum[config_id] += eff;
      n_entries[config_id]++;
    }
  }

  std::vector<double> effs;
  for (auto const& [id, eff] : efficiency_sum) {
    double final_eff = eff;
    if (!is_opt_mode)
      final_eff /= n_entries[id]; // Media sulle posizioni sorgente
    effs.push_back(final_eff);
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

  for (auto const& [id, eff_sum] : efficiency_sum) {
    double eff = eff_sum;
    if (!is_opt_mode)
      eff /= n_entries[id];

    if (eff < p_low || eff > p_high)
      continue;

    if (filtered_config.count(id)) {
      h_efficiency.Fill(filtered_config[id].x1, filtered_config[id].x2, eff);
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