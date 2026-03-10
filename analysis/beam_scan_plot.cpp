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

#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TTree.h>

#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

struct HitWithRun {
  int run_id;
  double y, z;
};

std::string format_double(double val, int precision = 1) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(precision) << val;
  return oss.str();
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " x1 x2 [input_file=root_file]" << std::endl;
    return 1;
  }

  double x1_target                 = std::stod(argv[1]);
  double x2_target                 = std::stod(argv[2]);
  std::filesystem::path input_file = (argc > 3) ? argv[3] : "output/lens_simulation/lens.root";
  std::filesystem::path output_dir = input_file.parent_path();

  TFile f(input_file.c_str());
  if (f.IsZombie()) {
    std::cerr << "Cannot open file " << input_file << std::endl;
    return 1;
  }

  TTree* hits    = (TTree*)f.Get("Hits");
  TTree* runs    = (TTree*)f.Get("Runs");
  TTree* configs = (TTree*)f.Get("Configurations");
  if (!hits || !runs || !configs) {
    std::cerr << "Hits, Runs, or Configurations trees not found!" << std::endl;
    return 1;
  }

  // --- Configurazioni ---
  int cfg_id_cfg;
  double cfg_x1, cfg_x2;
  configs->SetBranchAddress("config_id", &cfg_id_cfg);
  configs->SetBranchAddress("x1", &cfg_x1);
  configs->SetBranchAddress("x2", &cfg_x2);

  int selected_cfg_id = -1;
  for (Long64_t i = 0; i < configs->GetEntries(); i++) {
    configs->GetEntry(i);
    if (std::abs(cfg_x1 - x1_target) < 1e-6 && std::abs(cfg_x2 - x2_target) < 1e-6) {
      selected_cfg_id = cfg_id_cfg;
      break;
    }
  }

  if (selected_cfg_id == -1) {
    std::cerr << "Configuration with x1=" << x1_target << " and x2=" << x2_target << " not found!"
              << std::endl;
    return 1;
  }
  std::cout << "Selected configuration id: " << selected_cfg_id << std::endl;

  // --- Runs per la configurazione ---
  int run_id, config_id;
  double x_source;
  runs->SetBranchAddress("run_id", &run_id);
  runs->SetBranchAddress("config_id", &config_id);
  runs->SetBranchAddress("x_source", &x_source);

  std::map<int, double> runs_map; // run_id -> x_source
  for (Long64_t i = 0; i < runs->GetEntries(); i++) {
    runs->GetEntry(i);
    if (config_id == selected_cfg_id)
      runs_map[run_id] = x_source;
  }

  if (runs_map.empty()) {
    std::cerr << "No runs found for this configuration!" << std::endl;
    return 1;
  }

  // --- Hits ---
  int hit_run_id;
  double y_hit, z_hit;
  hits->SetBranchAddress("run_id", &hit_run_id);
  hits->SetBranchAddress("y_hit", &y_hit);
  hits->SetBranchAddress("z_hit", &z_hit);

  // Aggregazione per x_source
  std::map<double, std::vector<double>> y_hits_map;
  std::map<double, std::vector<double>> z_hits_map;

  for (Long64_t j = 0; j < hits->GetEntries(); j++) {
    hits->GetEntry(j);
    auto it = runs_map.find(hit_run_id);
    if (it != runs_map.end()) {
      double xs = it->second;
      y_hits_map[xs].push_back(y_hit);
      z_hits_map[xs].push_back(z_hit);
    }
  }

  // --- Calcolo medie e deviazioni ---
  std::vector<double> xs_vec, y_mean_vec, y_err_vec, z_mean_vec, z_err_vec;
  for (auto& [xs, yvec] : y_hits_map) {
    double y_sum = 0;
    for (double v : yvec)
      y_sum += v;
    double y_mean = y_sum / yvec.size();
    double y_std  = 0;
    for (double v : yvec)
      y_std += (v - y_mean) * (v - y_mean);
    y_std = std::sqrt(y_std / yvec.size());

    xs_vec.push_back(xs);
    y_mean_vec.push_back(y_mean);
    y_err_vec.push_back(y_std);

    auto& zvec   = z_hits_map[xs];
    double z_sum = 0;
    for (double v : zvec)
      z_sum += v;
    double z_mean = z_sum / zvec.size();
    double z_std  = 0;
    for (double v : zvec)
      z_std += (v - z_mean) * (v - z_mean);
    z_std = std::sqrt(z_std / zvec.size());

    z_mean_vec.push_back(z_mean);
    z_err_vec.push_back(z_std);
  }

  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetPadGridX(true);
  gStyle->SetPadGridY(true);

  TCanvas* c = new TCanvas("c", "Beam positions", 1400, 700);
  c->Divide(2, 1, 0.01, 0.01);

  // --- Grafico Y ---
  c->cd(1);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.08);
  gPad->SetBottomMargin(0.12);

  TGraphErrors* g_y =
      new TGraphErrors(xs_vec.size(), xs_vec.data(), y_mean_vec.data(), nullptr, y_err_vec.data());
  g_y->SetMarkerStyle(20);
  g_y->SetMarkerSize(1);
  g_y->SetMarkerColor(kRed);
  g_y->SetTitle(("Y positions - x1=" + format_double(x1_target)
                 + " mm, x2=" + format_double(x2_target) + " mm")
                    .c_str());
  g_y->GetXaxis()->SetTitle("Posizione iniziale fascio [mm]");
  g_y->GetYaxis()->SetTitle("Posizione y fotoni [mm]");
  g_y->Draw("AP");

  // --- Grafico Z ---
  c->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.08);
  gPad->SetBottomMargin(0.12);

  TGraphErrors* g_z =
      new TGraphErrors(xs_vec.size(), xs_vec.data(), z_mean_vec.data(), nullptr, z_err_vec.data());
  g_z->SetMarkerStyle(20);
  g_z->SetMarkerSize(1);
  g_z->SetMarkerColor(kBlue);
  g_z->SetTitle(("Z positions - x1=" + format_double(x1_target)
                 + " mm, x2=" + format_double(x2_target) + " mm")
                    .c_str());
  g_z->GetXaxis()->SetTitle("Posizione iniziale fascio [mm]");
  g_z->GetYaxis()->SetTitle("Posizione z fotoni [mm]");
  g_z->Draw("AP");

  std::string filename =
      (output_dir
       / ("beam_x1_" + format_double(x1_target, 2) + "_x2_" + format_double(x2_target, 2) + ".png"))
          .string();
  c->SaveAs(filename.c_str());

  delete g_y;
  delete g_z;
  delete c;

  std::cout << "Plot saved to " << filename << std::endl;

  return 0;
}