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
#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TPolyLine3D.h>
#include <TStyle.h>
#include <TTree.h>

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
  std::string lens1_id    = "";
  double lower_percentile = -1.0;
  double upper_percentile = -1.0;
};

static CliConfig parse_args(int argc, char** argv) {
  CliConfig cfg;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    auto next       = [&]() { return argv[++i]; };
    if (arg == "--input" || arg == "-i")
      cfg.input_file = next();
    else if (arg == "--config" || arg == "-c")
      cfg.config_file = next();
    else if (arg == "--output" || arg == "-o")
      cfg.output_png = next();
    else if (arg == "--lens1")
      cfg.lens1_id = next();
    else if (arg == "--low")
      cfg.lower_percentile = std::stod(next());
    else if (arg == "--high")
      cfg.upper_percentile = std::stod(next());
  }
  return cfg;
}

int main(int argc, char** argv) {
  using json    = nlohmann::json;
  CliConfig cli = parse_args(argc, argv);

  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
  gStyle->SetCanvasPreferGL(kFALSE);

  std::ifstream f(cli.config_file);
  if (!f.is_open())
    return 1;
  json config;
  f >> config;

  double dx = config.value("dx", 1.0);
  double lower_p =
      (cli.lower_percentile >= 0) ? cli.lower_percentile : config.value("lower_percentile", 0.0);
  double upper_p =
      (cli.upper_percentile >= 0) ? cli.upper_percentile : config.value("upper_percentile", 0.0);

  TFile file(cli.input_file.c_str(), "READ");
  TTree* tree_config = (TTree*)file.Get("Configurations");
  TTree* tree_eff    = (TTree*)file.Get("Efficiency");
  if (!tree_eff)
    tree_eff = (TTree*)file.Get("Runs");

  double x1_val, x2_val;
  int config_id, n_photons_val;
  double n_hits_val;
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
    std::string s1        = tree_config->GetBranch("lens75_id") ? l1_id_buf : "default";
    std::string s2        = tree_config->GetBranch("lens60_id") ? l2_id_buf : "default";
    config_map[config_id] = {x1_val, x2_val, s1, s2};
    lens1_set.insert(s1);
    lens2_set.insert(s2);
  }

  std::vector<std::string> lens1_list;
  if (!cli.lens1_id.empty())
    lens1_list.push_back(cli.lens1_id);
  else
    lens1_list.assign(lens1_set.begin(), lens1_set.end());

  std::vector<std::string> lens2_list(lens2_set.begin(), lens2_set.end());
  std::sort(lens2_list.rbegin(), lens2_list.rend());
  std::map<std::string, int> lens2_idx;
  for (int i = 0; i < lens2_list.size(); ++i)
    lens2_idx[lens2_list[i]] = i;

  tree_eff->SetBranchAddress("config_id", &config_id);
  tree_eff->SetBranchAddress("n_photons", &n_photons_val);
  tree_eff->SetBranchAddress("n_hits", &n_hits_val);

  std::map<int, double> efficiency_sum;
  std::map<int, int> n_entries;
  std::vector<double> all_effs;
  for (int i = 0; i < tree_eff->GetEntries(); ++i) {
    tree_eff->GetEntry(i);
    double eff = (double)n_hits_val / n_photons_val;
    efficiency_sum[config_id] += eff;
    n_entries[config_id]++;
  }
  for (auto const& [id, sum] : efficiency_sum)
    all_effs.push_back(sum / n_entries[id]);
  std::sort(all_effs.begin(), all_effs.end());
  double p_low  = all_effs[static_cast<size_t>(lower_p * all_effs.size())];
  double p_high = all_effs[std::min(all_effs.size() - 1,
                                    static_cast<size_t>((1.0 - upper_p) * all_effs.size()))];

  double x1_min = 1e9, x1_max = -1e9, x2_min = 1e9, x2_max = -1e9;
  for (auto const& [id, cfg] : config_map) {
    x1_min = std::min(x1_min, cfg.x1);
    x1_max = std::max(x1_max, cfg.x1);
    x2_min = std::min(x2_min, cfg.x2);
    x2_max = std::max(x2_max, cfg.x2);
  }

  int n_x1 = std::round((x1_max - x1_min) / dx) + 1;
  int n_x2 = std::round((x2_max - x2_min) / dx) + 1;
  int n_l2 = lens2_list.size();

  TCanvas canvas("c", "3D Efficiency", 1300, 800);
  if (lens1_list.size() > 1) {
    int n_rows = std::sqrt(lens1_list.size());
    int n_cols = (lens1_list.size() + n_rows - 1) / n_rows;
    canvas.Divide(n_cols, n_rows);
  }

  std::vector<TH3D*> hists;
  std::vector<TPolyLine3D*> lines;
  std::vector<TPolyLine3D*> grid_cage; // Memoria per la griglia

  for (int i = 0; i < lens1_list.size(); ++i) {
    if (lens1_list.size() > 1)
      canvas.cd(i + 1);
    else
      canvas.cd();

    gPad->SetTheta(30);
    gPad->SetPhi(150);
    gPad->SetRightMargin(0.15);

    std::string l1 = lens1_list[i];
    int n_x_bins   = 2 * n_l2 - 1;
    TH3D* h = new TH3D(Form("h_%d", i), Form("Lens 1: %s;Lens 2;X1 [mm];X2 [mm]", l1.c_str()),
                       n_x_bins, -0.5, n_x_bins - 0.5, n_x1, x1_min - dx / 2, x1_max + dx / 2, n_x2,
                       x2_min - dx / 2, x2_max + dx / 2);

    // --- ASSE X1 (Y) ---
    for (int b = 1; b <= n_x1; ++b)
      h->GetYaxis()->SetBinLabel(b, "");
    for (double v = 40; v <= 160; v += 20) {
      double x1_p = x1_max - (v - x1_min);
      int bin     = h->GetYaxis()->FindBin(x1_p);
      if (bin >= 1 && bin <= n_x1)
        h->GetYaxis()->SetBinLabel(bin, Form("%.0f", v));
    }
    // --- ASSE LENS 2 (X) ---
    h->GetXaxis()->CenterLabels();
    for (int j = 0; j < n_l2; ++j)
      h->GetXaxis()->SetBinLabel(2 * j + 1, lens2_list[j].c_str());

    h->GetYaxis()->SetLabelSize(0.04);
    h->GetXaxis()->SetLabelSize(0.04);

    std::map<int, double> max_eff_val;
    std::map<int, double> x1_at_max, x2_at_max;

    for (auto const& [id, sum] : efficiency_sum) {
      if (config_map[id].lens1_id != l1)
        continue;
      double eff = sum / n_entries[id];
      if (eff < p_low || eff > p_high)
        continue;
      int l2_plot_idx = 2 * lens2_idx[config_map[id].lens2_id];
      double x1_p     = x1_max - (config_map[id].x1 - x1_min);
      h->Fill(l2_plot_idx, x1_p, config_map[id].x2, eff);
      if (eff > max_eff_val[l2_plot_idx]) {
        max_eff_val[l2_plot_idx] = eff;
        x1_at_max[l2_plot_idx]   = x1_p;
        x2_at_max[l2_plot_idx]   = config_map[id].x2;
      }
    }

    h->SetMinimum(p_low);
    h->SetMaximum(p_high);
    h->SetFillColorAlpha(kAzure, 0.35);
    h->Draw("BOX2 Z");
    h->GetZaxis()->SetTitle("Efficienza [a.d.]");
    h->GetZaxis()->CenterTitle(kTRUE);
    h->GetZaxis()->SetTitleOffset(1.6);
    gPad->Update();

    // --- AGGIUNTA GRIGLIA DI SEPARAZIONE (GRID CAGE) ---
    double z_min = x2_min - dx / 2, z_max = x2_max + dx / 2;
    double x_min_p = -0.5, x_max_p = n_x_bins - 0.5;
    double y_min_p = x1_min - dx / 2, y_max_p = x1_max + dx / 2;

    // Disegna separatori tra le lenti
    for (int j = 0; j < n_l2; ++j) {
      double left  = 2 * j - 0.5;
      double right = 2 * j + 0.5;

      // 1. Disegna un riquadro sul pavimento per ogni lente
      TPolyLine3D* box = new TPolyLine3D(5);
      box->SetPoint(0, left, y_min_p, z_min);
      box->SetPoint(1, right, y_min_p, z_min);
      box->SetPoint(2, right, y_max_p, z_min);
      box->SetPoint(3, left, y_max_p, z_min);
      box->SetPoint(4, left, y_min_p, z_min);
      box->SetLineColor(kBlack);
      box->SetLineWidth(1);
      box->SetLineStyle(1);
      box->Draw("SAME");
      grid_cage.push_back(box);

      // 2. Disegna linee verticali agli angoli del settore lente (muro posteriore)
      TPolyLine3D* v1 = new TPolyLine3D(2);
      v1->SetPoint(0, left, y_max_p, z_min);
      v1->SetPoint(1, left, y_max_p, z_max);
      v1->SetLineColor(kGray + 1);
      v1->SetLineStyle(2);
      v1->Draw("SAME");
      grid_cage.push_back(v1);

      TPolyLine3D* v2 = new TPolyLine3D(2);
      v2->SetPoint(0, right, y_max_p, z_min);
      v2->SetPoint(1, right, y_max_p, z_max);
      v2->SetLineColor(kGray + 1);
      v2->SetLineStyle(2);
      v2->Draw("SAME");
      grid_cage.push_back(v2);
    }

    // --- DISEGNO LINEA ROSSA ---
    std::vector<double> lx, ly, lz;
    for (int j = 0; j < n_l2; ++j) {
      int idx = 2 * j;
      if (max_eff_val.count(idx)) {
        lx.push_back(h->GetXaxis()->GetBinCenter(idx + 1));
        ly.push_back(x1_at_max[idx]);
        lz.push_back(x2_at_max[idx]);
      }
    }
    if (!lx.empty()) {
      TPolyLine3D* line = new TPolyLine3D(lx.size());
      for (size_t k = 0; k < lx.size(); ++k)
        line->SetPoint(k, lx[k], ly[k], lz[k]);
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