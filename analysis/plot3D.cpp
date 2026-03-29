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
#include <TH3D.h>
#include <TPolyLine3D.h>
#include <TStyle.h>
#include <TTree.h>
#include <TGraph2D.h>
#include <TAxis.h>
#include <TText.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
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
  std::string output_png  = "output/efficiency3D.png";
  std::string lens1_id    = ""; // Se vuoto, mostra tutte in un grid
  bool show_2d            = false;
  double lower_percentile = -1.0;
  double upper_percentile = -1.0;
};

static CliConfig parse_args(int argc, char** argv) {
  CliConfig cfg;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    auto next       = [&]() -> std::string {
      if (i + 1 >= argc) {
        std::cerr << "Missing argument after " << arg << "\n";
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
    else if (arg == "--2d")
      cfg.show_2d = true;
    else if (arg == "--low")
      cfg.lower_percentile = std::stod(next());
    else if (arg == "--high")
      cfg.upper_percentile = std::stod(next());
    else if (arg == "--help" || arg == "-h") {
      std::cout << "Usage: plot3D [options]\n"
                << "Options:\n"
                << "  -i, --input <file>    Input ROOT file\n"
                << "  -c, --config <file>   Config JSON file\n"
                << "  -o, --output <file>   Output PNG file\n"
                << "  --lens1 <id>          ID of the first lens to plot specifically\n"
                << "  --2d                  Show a grid of 2D plots for each lens2 (only with "
                   "--lens1)\n"
                << "  --low <val>           Lower percentile\n"
                << "  --high <val>          Upper percentile\n";
      std::exit(0);
    }
  }
  return cfg;
}

int main(int argc, char** argv) {
  using json = nlohmann::json;
  CliConfig cli = parse_args(argc, argv);

  std::ifstream f(cli.config_file);
  if (!f.is_open()) return 1;
  json config;
  f >> config;

  double dx = config.value("dx", 1.0);
  double lower_p = (cli.lower_percentile >= 0) ? cli.lower_percentile : config.value("lower_percentile", 0.0);
  double upper_p = (cli.upper_percentile >= 0) ? cli.upper_percentile : config.value("upper_percentile", 0.0);

  TFile file(cli.input_file.c_str(), "READ");
  if (!file.IsOpen()) return 1;

  TTree* tree_config = (TTree*)file.Get("Configurations");
  TTree* tree_eff = (TTree*)file.Get("Efficiency");
  bool is_opt = true;
  if (!tree_eff) { tree_eff = (TTree*)file.Get("Runs"); is_opt = false; }
  if (!tree_config || !tree_eff) return 1;

  double x1_val, x2_val;
  int config_id, n_photons_val, n_hits_val;
  char l1_id_buf[256], l2_id_buf[256];

  tree_config->SetBranchAddress("x1", &x1_val);
  tree_config->SetBranchAddress("x2", &x2_val);
  tree_config->SetBranchAddress("config_id", &config_id);
  if (tree_config->GetBranch("lens75_id")) {
    tree_config->SetBranchAddress("lens75_id", l1_id_buf);
    tree_config->SetBranchAddress("lens60_id", l2_id_buf);
  }

  std::map<int, Config> config_map;
  std::set<std::string> lens1_set, lens2_set;
  for (int i = 0; i < tree_config->GetEntries(); ++i) {
    tree_config->GetEntry(i);
    std::string s1 = tree_config->GetBranch("lens75_id") ? l1_id_buf : "default";
    std::string s2 = tree_config->GetBranch("lens60_id") ? l2_id_buf : "default";
    config_map[config_id] = {x1_val, x2_val, s1, s2};
    lens1_set.insert(s1);
    lens2_set.insert(s2);
  }

  std::vector<std::string> lens1_list;
  if (!cli.lens1_id.empty()) {
    if (lens1_set.count(cli.lens1_id)) {
      lens1_list.push_back(cli.lens1_id);
    } else {
      std::cerr << "Lens 1 ID " << cli.lens1_id << " not found in file.\n";
      return 1;
    }
  } else {
    lens1_list.assign(lens1_set.begin(), lens1_set.end());
  }

  std::vector<std::string> lens2_list(lens2_set.begin(), lens2_set.end());
  std::map<std::string, int> lens2_idx;
  for (int i = 0; i < lens2_list.size(); ++i) lens2_idx[lens2_list[i]] = i;

  if (is_opt) {
    tree_eff->SetBranchAddress("config_id", &config_id);
    tree_eff->SetBranchAddress("n_photons", &n_photons_val);
    tree_eff->SetBranchAddress("n_hits", &n_hits_val);
  } else {
    tree_eff->SetBranchAddress("config_id", &config_id);
    tree_eff->SetBranchAddress("n_hits", &n_hits_val);
    n_photons_val = config.value("n_photons", 10000);
  }

  std::map<int, double> efficiency_sum;
  std::map<int, int> n_entries;
  std::vector<double> all_effs;

  for (int i = 0; i < tree_eff->GetEntries(); ++i) {
    tree_eff->GetEntry(i);
    double eff = (double)n_hits_val / n_photons_val;
    efficiency_sum[config_id] += eff;
    n_entries[config_id]++;
  }

  for (auto const& [id, sum] : efficiency_sum) {
    double eff = is_opt ? sum : sum / n_entries[id];
    all_effs.push_back(eff);
  }

  std::sort(all_effs.begin(), all_effs.end());
  double p_low = all_effs[static_cast<size_t>(lower_p * all_effs.size())];
  double p_high = all_effs[std::min(all_effs.size() - 1, static_cast<size_t>((1.0 - upper_p) * all_effs.size()))];

  double x1_min = 1e9, x1_max = -1e9, x2_min = 1e9, x2_max = -1e9;
  for (auto const& [id, cfg] : config_map) {
    x1_min = std::min(x1_min, cfg.x1); x1_max = std::max(x1_max, cfg.x1);
    x2_min = std::min(x2_min, cfg.x2); x2_max = std::max(x2_max, cfg.x2);
  }

  int n_x1 = std::round((x1_max - x1_min) / dx) + 1;
  int n_x2 = std::round((x2_max - x2_min) / dx) + 1;
  int n_l2 = lens2_list.size();

  if (cli.show_2d && !cli.lens1_id.empty()) {
    std::string l1 = cli.lens1_id;
    TCanvas c2d("c2d", Form("2D Efficiency for Lens 1: %s", l1.c_str()), 1200, 900);
    int n_rows_2d = std::sqrt(n_l2);
    int n_cols_2d = (n_l2 + n_rows_2d - 1) / n_rows_2d;
    c2d.Divide(n_cols_2d, n_rows_2d);

    std::vector<TH2D*> h2s;
    for (int j = 0; j < n_l2; ++j) {
      c2d.cd(j + 1);
      std::string l2 = lens2_list[j];
      TH2D* h2       = new TH2D(Form("h2_%d", j), Form("%s;X1 [mm];X2 [mm]", l2.c_str()), n_x1,
                                x1_min - dx / 2, x1_max + dx / 2, n_x2, x2_min - dx / 2, x2_max + dx / 2);

      for (auto const& [id, sum] : efficiency_sum) {
        if (config_map[id].lens1_id == l1 && config_map[id].lens2_id == l2) {
          double eff = is_opt ? sum : sum / n_entries[id];
          if (eff >= p_low && eff <= p_high) {
            h2->Fill(config_map[id].x1, config_map[id].x2, eff);
          }
        }
      }
      h2->SetMinimum(p_low);
      h2->SetMaximum(p_high);
      h2->SetStats(0);
      h2->Draw("COLZ");
      h2s.push_back(h2);
    }
    c2d.SaveAs(cli.output_png.c_str());
    file.Close();
    return 0;
  }

  TCanvas canvas("c", "3D Efficiency", 1200, 800);
  if (lens1_list.size() > 1) {
    int n_rows = std::sqrt(lens1_list.size());
    int n_cols = (lens1_list.size() + n_rows - 1) / n_rows;
    canvas.Divide(n_cols, n_rows);
  }

  std::vector<TH3D*> hists;
  std::vector<TPolyLine3D*> lines;

  for (int i = 0; i < lens1_list.size(); ++i) {
    if (lens1_list.size() > 1)
      canvas.cd(i + 1);
    else
      canvas.cd();

    // Imposta l'angolo di visuale per invertire la percezione di X1 (Y axis)
    // Phi default è 30. Portandolo a 150-210 ruotiamo il cubo.
    gPad->SetTheta(30);
    gPad->SetPhi(150);

    std::string l1 = lens1_list[i];
    TH3D* h =
        new TH3D(Form("h_%d", i), Form("Lens 1: %s;Lens 2;X1 [mm];X2 [mm]", l1.c_str()), n_l2, -0.5,
                 n_l2 - 0.5, n_x1, x1_min - dx / 2, x1_max + dx / 2, n_x2, x2_min - dx / 2,
                 x2_max + dx / 2);

    // Imposta le label per le lenti 2
    for (int j = 0; j < n_l2; ++j) {
      h->GetXaxis()->SetBinLabel(j + 1, lens2_list[j].c_str());
    }
    h->GetXaxis()->SetLabelSize(0.05);

    std::map<int, double> sum_x1, sum_x2, sum_w;

    for (auto const& [id, sum] : efficiency_sum) {
      if (config_map[id].lens1_id != l1) continue;
      double eff = is_opt ? sum : sum / n_entries[id];
      if (eff < p_low || eff > p_high) continue;
      
      int l2_i = lens2_idx[config_map[id].lens2_id];
      h->Fill(l2_i, config_map[id].x1, config_map[id].x2, eff);
      
      sum_x1[l2_i] += config_map[id].x1 * eff;
      sum_x2[l2_i] += config_map[id].x2 * eff;
      sum_w[l2_i] += eff;
    }

    h->SetMinimum(p_low);
    h->SetMaximum(p_high);
    h->Draw("BOX2Z");
    gPad->Update(); // Forza la creazione del sistema di coordinate 3D

    std::vector<double> lx, ly, lz;
    for (int j = 0; j < n_l2; ++j) {
      if (sum_w[j] > 0) {
        lx.push_back(j);
        ly.push_back(sum_x1[j] / sum_w[j]);
        lz.push_back(sum_x2[j] / sum_w[j]);
      }
    }

    if (!lx.empty()) {
      TPolyLine3D* line = new TPolyLine3D(lx.size());
      for (size_t k = 0; k < lx.size(); ++k) {
        line->SetPoint(k, lx[k], ly[k], lz[k]);
      }
      line->SetLineColor(kRed + 1);
      line->SetLineWidth(4);
      line->Draw("SAME");
      lines.push_back(line);
    }

    hists.push_back(h);
  }

  canvas.SaveAs(cli.output_png.c_str());
  file.Close();
  return 0;
}
