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

// beam_scan_plot.cpp
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <vector>

struct DetectorHit {
  double y, z;
};

std::string format_double(double val, int precision = 1) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(precision) << val;
  return oss.str();
}

int main() {
  std::filesystem::path input_file = "output/lens_simulation/lens.root";
  std::filesystem::path output_dir = input_file.parent_path();

  TFile f(input_file.c_str());
  if (f.IsZombie()) {
    std::cerr << "Cannot open file " << input_file << std::endl;
    return 1;
  }

  TTree* hits    = (TTree*)f.Get("Hits");
  TTree* runs    = (TTree*)f.Get("Runs");
  TTree* configs = (TTree*)f.Get("Configurations"); // contiene x1 e x2

  if (!hits || !runs || !configs) {
    std::cerr << "Hits, Runs, or Configurations trees not found!" << std::endl;
    return 1;
  }

  // Lettura configurazioni
  int cfg_id_cfg;
  double cfg_x1, cfg_x2;
  configs->SetBranchAddress("config_id", &cfg_id_cfg);
  configs->SetBranchAddress("x1", &cfg_x1);
  configs->SetBranchAddress("x2", &cfg_x2);

  std::map<int, std::pair<double, double>> config_values;
  Long64_t n_cfg = configs->GetEntries();
  for (Long64_t i = 0; i < n_cfg; i++) {
    configs->GetEntry(i);
    config_values[cfg_id_cfg] = {cfg_x1, cfg_x2};
  }

  // Lettura runs
  int run_id, config_id;
  double x_source;
  runs->SetBranchAddress("run_id", &run_id);
  runs->SetBranchAddress("config_id", &config_id);
  runs->SetBranchAddress("x_source", &x_source);

  std::map<int, std::vector<std::pair<int, double>>> config_to_runs;
  Long64_t n_runs = runs->GetEntries();
  for (Long64_t i = 0; i < n_runs; i++) {
    runs->GetEntry(i);
    config_to_runs[config_id].push_back({run_id, x_source});
  }

  // Lettura hits
  int hit_run_id;
  double y_hit, z_hit;
  hits->SetBranchAddress("run_id", &hit_run_id);
  hits->SetBranchAddress("y_hit", &y_hit);
  hits->SetBranchAddress("z_hit", &z_hit);

  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.05);

  // Selezione casuale di 10 configurazioni
  std::vector<int> all_configs;
  for (auto& kv : config_to_runs)
    all_configs.push_back(kv.first);
  std::shuffle(all_configs.begin(), all_configs.end(), std::mt19937{std::random_device{}()});
  if (all_configs.size() > 10)
    all_configs.resize(10);

  TCanvas* c     = new TCanvas("c", "Beam positions", 1200, 600);
  int file_index = 1;

  for (auto cfg_id : all_configs) {
    auto& runs_vec = config_to_runs[cfg_id];

    c->Clear();
    c->cd();
    // Dividiamo il canvas in due grafici, uno per y e uno per z
    c->Divide(2, 1);

    // Colori disponibili
    std::vector<int> colors = {kRed,      kBlue,       kGreen + 2, kMagenta,
                               kCyan + 1, kOrange + 7, kViolet,    kTeal + 2};

    // --- Grafico Y ---
    c->cd(1);
    TGraph* g_y_all = new TGraph();
    int point_index = 0;

    for (size_t i = 0; i < runs_vec.size(); i++) {
      int run_i = runs_vec[i].first;
      double x0 = runs_vec[i].second; // posizione iniziale del fascio

      std::vector<double> y_hits_vec;
      Long64_t n_hits = hits->GetEntries();
      for (Long64_t j = 0; j < n_hits; j++) {
        hits->GetEntry(j);
        if (hit_run_id == run_i)
          y_hits_vec.push_back(y_hit);
      }

      int color = colors[i % colors.size()]; // colori ciclici
      for (auto y_val : y_hits_vec) {
        g_y_all->SetPoint(point_index++, x0, y_val);
      }
    }
    g_y_all->SetMarkerStyle(20);
    g_y_all->SetMarkerSize(1);
    g_y_all->SetMarkerColor(kRed);
    g_y_all->SetTitle(("Posizione Y - x1=" + format_double(config_values[cfg_id].first)
                       + " mm, x2=" + format_double(config_values[cfg_id].second) + " mm")
                          .c_str());
    g_y_all->GetXaxis()->SetTitle("Posizione iniziale fascio [mm]");
    g_y_all->GetYaxis()->SetTitle("Posizione y fotoni [mm]");
    g_y_all->Draw("AP");

    // --- Grafico Z ---
    c->cd(2);
    TGraph* g_z_all = new TGraph();
    point_index     = 0;

    for (size_t i = 0; i < runs_vec.size(); i++) {
      int run_i = runs_vec[i].first;
      double x0 = runs_vec[i].second;

      std::vector<double> z_hits_vec;
      Long64_t n_hits = hits->GetEntries();
      for (Long64_t j = 0; j < n_hits; j++) {
        hits->GetEntry(j);
        if (hit_run_id == run_i)
          z_hits_vec.push_back(z_hit);
      }

      int color = colors[i % colors.size()];
      for (auto z_val : z_hits_vec) {
        g_z_all->SetPoint(point_index++, x0, z_val);
      }
    }
    g_z_all->SetMarkerStyle(20);
    g_z_all->SetMarkerSize(1);
    g_z_all->SetMarkerColor(kBlue);
    g_z_all->SetTitle(("Posizione Z - x1=" + format_double(config_values[cfg_id].first)
                       + " mm, x2=" + format_double(config_values[cfg_id].second) + " mm")
                          .c_str());
    g_z_all->GetXaxis()->SetTitle("Posizione iniziale fascio [mm]");
    g_z_all->GetYaxis()->SetTitle("Posizione z fotoni [mm]");
    g_z_all->Draw("AP");

    std::string filename = (output_dir / (std::to_string(file_index) + ".png")).string();
    c->SaveAs(filename.c_str());
    file_index++;

    delete g_y_all;
    delete g_z_all;
  }

  delete c;
  return 0;
}