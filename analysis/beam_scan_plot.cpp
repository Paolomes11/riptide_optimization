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
#include <TF2.h>
#include <TFile.h>
#include <TGraph2DErrors.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
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

struct Point3D {
  double x_src, y_src;
  double mean_y, err_y;
  double mean_z, err_z;
};

std::string format_double(double val, int precision = 1) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(precision) << val;
  return oss.str();
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " x1 x2 [input_file=root_file] [x0_1,x0_2,...]"
              << std::endl;
    return 1;
  }

  double x1_target                 = std::stod(argv[1]);
  double x2_target                 = std::stod(argv[2]);
  std::filesystem::path input_file = (argc > 3) ? argv[3] : "output/lens_simulation/lens.root";
  std::filesystem::path output_dir = input_file.parent_path();

  std::vector<double> x0_filters;
  if (argc > 4) {
    std::string x0_str = argv[4];
    std::stringstream ss(x0_str);
    std::string item;
    while (std::getline(ss, item, ',')) {
      x0_filters.push_back(std::stod(item));
    }
  }

  TFile f(input_file.c_str());
  if (f.IsZombie()) {
    std::cerr << "Cannot open file " << input_file << std::endl;
    return 1;
  }

  TTree* runs    = (TTree*)f.Get("Runs");
  TTree* configs = (TTree*)f.Get("Configurations");
  if (!runs || !configs) {
    std::cerr << "Required trees not found in " << input_file << "!" << std::endl;
    return 1;
  }

  // Configurazioni
  int cfg_id_cfg;
  double cfg_x1, cfg_x2;
  configs->SetBranchAddress("config_id", &cfg_id_cfg);
  configs->SetBranchAddress("x1", &cfg_x1);
  configs->SetBranchAddress("x2", &cfg_x2);

  int selected_cfg_id = -1;
  double best_dist    = 1e18;
  double found_x1 = 0, found_x2 = 0;

  for (Long64_t i = 0; i < configs->GetEntries(); i++) {
    configs->GetEntry(i);
    double dist = std::hypot(cfg_x1 - x1_target, cfg_x2 - x2_target);
    if (dist < best_dist) {
      best_dist       = dist;
      selected_cfg_id = cfg_id_cfg;
      found_x1        = cfg_x1;
      found_x2        = cfg_x2;
    }
  }

  if (selected_cfg_id == -1) {
    std::cerr << "No configurations found in " << input_file << "!" << std::endl;
    return 1;
  }

  if (best_dist > 1e-4) {
    std::cout << "[WARNING] Config (x1=" << x1_target << ", x2=" << x2_target
              << ") non trovata. Uso la più vicina: (x1=" << found_x1 << ", x2=" << found_x2
              << "), dist=" << best_dist << " mm" << std::endl;
  } else {
    std::cout << "Selected configuration id: " << selected_cfg_id << " (x1=" << found_x1
              << ", x2=" << found_x2 << ")" << std::endl;
  }

  // Runs per la configurazione
  int run_id, config_id, n_hits_run;
  float x_source, y_source;
  std::vector<float>* y_hits_ptr = nullptr;
  std::vector<float>* z_hits_ptr = nullptr;

  runs->SetBranchAddress("run_id", &run_id);
  runs->SetBranchAddress("config_id", &config_id);
  runs->SetBranchAddress("x_source", &x_source);
  runs->SetBranchAddress("y_source", &y_source);
  runs->SetBranchAddress("n_hits", &n_hits_run);
  runs->SetBranchAddress("y_hits", &y_hits_ptr);
  runs->SetBranchAddress("z_hits", &z_hits_ptr);

  // Aggregazione per (x_source, y_source)
  std::map<std::pair<double, double>, std::vector<double>> y_hits_map;
  std::map<std::pair<double, double>, std::vector<double>> z_hits_map;

  for (Long64_t i = 0; i < runs->GetEntries(); i++) {
    runs->GetEntry(i);
    if (config_id != selected_cfg_id)
      continue;

    double x_src_rounded = std::round(static_cast<double>(x_source) * 10.0) / 10.0;
    double y_src_rounded = std::round(static_cast<double>(y_source) * 10.0) / 10.0;
    auto key             = std::make_pair(x_src_rounded, y_src_rounded);
    auto& yvec           = y_hits_map[key];
    auto& zvec           = z_hits_map[key];

    if (y_hits_ptr && z_hits_ptr) {
      for (size_t j = 0; j < y_hits_ptr->size(); ++j) {
        yvec.push_back(static_cast<double>((*y_hits_ptr)[j]));
        zvec.push_back(static_cast<double>((*z_hits_ptr)[j]));
      }
    }
  }

  if (y_hits_map.empty()) {
    std::cerr << "No runs found for this configuration!" << std::endl;
    return 1;
  }

  // Calcolo medie e deviazioni robuste (filtro 3σ)
  std::vector<Point3D> points;
  const double max_sigma    = 3.0;
  const size_t min_hits_req = 10; // Soglia minima di hit nel detector

  for (auto& [src_pos, yvec_orig] : y_hits_map) {
    auto& zvec_orig = z_hits_map[src_pos];

    if (yvec_orig.size() < min_hits_req) {
      std::cout << "[INFO] Saltato punto (x=" << src_pos.first << ", y=" << src_pos.second
                << ") per pochi hit: " << yvec_orig.size() << std::endl;
      continue;
    }

    std::vector<double> yvec = yvec_orig;
    std::vector<double> zvec = zvec_orig;
    double y_mean = 0, z_mean = 0;
    double y_sigma = 0, z_sigma = 0;

    // due iterazioni per stabilizzare media e sigma
    for (int iter = 0; iter < 2; ++iter) {
      double y_sum = 0, z_sum = 0;
      for (double v : yvec)
        y_sum += v;
      for (double v : zvec)
        z_sum += v;

      size_t n = yvec.size();
      if (n < min_hits_req)
        break;
      y_mean = y_sum / n;
      z_mean = z_sum / n;

      double y_var = 0, z_var = 0;
      for (double v : yvec)
        y_var += (v - y_mean) * (v - y_mean);
      for (double v : zvec)
        z_var += (v - z_mean) * (v - z_mean);
      y_sigma = std::sqrt(y_var / n);
      z_sigma = std::sqrt(z_var / n);

      std::vector<double> y_temp, z_temp;
      for (size_t i = 0; i < yvec.size(); ++i) {
        if (std::abs(yvec[i] - y_mean) <= max_sigma * y_sigma
            && std::abs(zvec[i] - z_mean) <= max_sigma * z_sigma) {
          y_temp.push_back(yvec[i]);
          z_temp.push_back(zvec[i]);
        }
      }
      yvec = y_temp;
      zvec = z_temp;
    }

    if (yvec.size() >= min_hits_req) {
      points.push_back({src_pos.first, src_pos.second, y_mean, y_sigma, z_mean, z_sigma});
    }
  }

  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.04);
  gStyle->SetPadGridX(true);
  gStyle->SetPadGridY(true);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetFrameFillColor(kWhite);

  TCanvas* c = new TCanvas("c", "Beam positions 3D", 1600, 800);
  c->Divide(2, 1, 0.01, 0.01);

  // Funzioni di fit (piano: z = a + b*x + c*y)
  TF2* fit_y = new TF2("fit_y", "[0] + [1]*x + [2]*y", -50, 50, -50, 50);
  fit_y->SetParameters(0, 0, 0);
  fit_y->SetParNames("Intercept", "Slope_y0", "Slope_x0");
  fit_y->SetNpx(50);
  fit_y->SetNpy(50);

  TF2* fit_z = new TF2("fit_z", "[0] + [1]*x + [2]*y", -50, 50, -50, 50);
  fit_z->SetParameters(0, 0, 0);
  fit_z->SetParNames("Intercept", "Slope_y0", "Slope_x0");
  fit_z->SetNpx(50);
  fit_z->SetNpy(50);

  // Grafico Y
  c->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.1);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.12);

  TGraph2DErrors* g_y = new TGraph2DErrors(points.size());
  for (int i = 0; i < points.size(); ++i) {
    // X = y_source, Y = x_source, Z = mean_y
    g_y->SetPoint(i, points[i].y_src, points[i].x_src, points[i].mean_y);
    g_y->SetPointError(i, 0, 0, points[i].err_y);
  }

  // Eseguo il fit
  g_y->Fit(fit_y, "QW");

  g_y->SetMarkerStyle(20);
  g_y->SetMarkerSize(0.8);
  g_y->SetMarkerColor(kRed);
  g_y->SetLineColor(kRed);
  g_y->SetTitle(("Centroid Y in detector (x1=" + format_double(found_x1)
                 + ", x2=" + format_double(found_x2) + ")")
                    .c_str());
  g_y->GetXaxis()->SetTitle("y_{0} sorgente [mm]");
  g_y->GetYaxis()->SetTitle("x_{0} sorgente [mm]");
  g_y->GetZaxis()->SetTitle("y detector [mm]");
  g_y->GetXaxis()->SetTitleOffset(1.5);
  g_y->GetYaxis()->SetTitleOffset(1.5);
  g_y->GetZaxis()->SetTitleOffset(1.5);
  g_y->Draw("P ERR");

  // Disegno il piano di fit con trasparenza
  fit_y->SetRange(g_y->GetXaxis()->GetXmin(), g_y->GetYaxis()->GetXmin(),
                  g_y->GetXaxis()->GetXmax(), g_y->GetYaxis()->GetXmax());
  fit_y->SetLineColorAlpha(kOrange + 1, 0.1); // 0.4 di trasparenza
  fit_y->SetFillColorAlpha(kOrange + 1, 0.05);
  //fit_y->Draw("SURF4 SAME");

  // Grafico Z
  c->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.1);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.12);

  TGraph2DErrors* g_z = new TGraph2DErrors(points.size());
  for (int i = 0; i < points.size(); ++i) {
    g_z->SetPoint(i, points[i].y_src, points[i].x_src, points[i].mean_z);
    g_z->SetPointError(i, 0, 0, points[i].err_z);
  }

  // Eseguo il fit
  g_z->Fit(fit_z, "Q");

  g_z->SetMarkerStyle(20);
  g_z->SetMarkerSize(0.8);
  g_z->SetMarkerColor(kBlue);
  g_z->SetLineColor(kBlue);
  g_z->SetTitle(("Centroid Z in detector (x1=" + format_double(found_x1)
                 + ", x2=" + format_double(found_x2) + ")")
                    .c_str());
  g_z->GetXaxis()->SetTitle("y_{0} sorgente [mm]");
  g_z->GetYaxis()->SetTitle("x_{0} sorgente [mm]");
  g_z->GetZaxis()->SetTitle("z detector [mm]");
  g_z->GetXaxis()->SetTitleOffset(1.5);
  g_z->GetYaxis()->SetTitleOffset(1.5);
  g_z->GetZaxis()->SetTitleOffset(1.5);
  g_z->Draw("P ERR");

  // Disegno il piano di fit con trasparenza
  fit_z->SetRange(g_z->GetXaxis()->GetXmin(), g_z->GetYaxis()->GetXmin(),
                  g_z->GetXaxis()->GetXmax(), g_z->GetYaxis()->GetXmax());
  fit_z->SetLineColorAlpha(kCyan + 1, 0.1);
  fit_z->SetFillColorAlpha(kCyan + 1, 0.05);
  // fit_z->Draw("SURF4 SAME");

  // Stampa parametri sul terminale
  auto print_fit = [](const char* name, TF2* f) {
    std::cout << "Fit Plane " << name << ": z = " << f->GetParameter(0) << " + ("
              << f->GetParameter(1) << ")*y0 + (" << f->GetParameter(2) << ")*x0" << std::endl;
  };
  print_fit("Y", fit_y);
  print_fit("Z", fit_z);

  std::string filename = (output_dir
                          / ("beam_scan_3D_x1_" + format_double(found_x1, 2) + "_x2_"
                             + format_double(found_x2, 2) + ".png"))
                             .string();
  c->SaveAs(filename.c_str());

  // Pulizia (ROOT handles deletion of objects attached to pads/canvas usually,
  // but explicitly deleting if they are pointers is safer if not added to a list)
  // delete g_y; // TGraph2D is often managed by the pad
  // delete g_z;
  delete c;

  std::cout << "Plot saved to " << filename << std::endl;

  // --- Parte 2: Plot 2D se richiesti x0_filters ---
  if (!x0_filters.empty()) {
    TCanvas* c2d = new TCanvas("c2d", "Beam positions 2D slices", 1400, 700);
    c2d->Divide(2, 1, 0.01, 0.01);

    TMultiGraph* mg_y = new TMultiGraph();
    mg_y->SetTitle(("Y Detector vs Y0 Source (x1=" + format_double(found_x1)
                    + ", x2=" + format_double(found_x2) + ");y_{0} sorgente [mm];y detector [mm]")
                       .c_str());

    TMultiGraph* mg_z = new TMultiGraph();
    mg_z->SetTitle(("Z Detector vs Y0 Source (x1=" + format_double(found_x1)
                    + ", x2=" + format_double(found_x2) + ");y_{0} sorgente [mm];z detector [mm]")
                       .c_str());

    TLegend* leg = new TLegend(0.1, 0.7, 0.3, 0.9);

    int colors[]  = {kRed,        kBlue,     kGreen + 2, kMagenta,
                     kOrange + 7, kCyan + 2, kAzure + 1, kGray + 2};
    int color_idx = 0;

    for (double x0_target : x0_filters) {
      std::vector<double> y0_src, y_det, y_err, z_det, z_err;

      for (const auto& p : points) {
        if (std::abs(p.x_src - x0_target) < 0.1) {
          y0_src.push_back(p.y_src);
          y_det.push_back(p.mean_y);
          y_err.push_back(p.err_y);
          z_det.push_back(p.mean_z);
          z_err.push_back(p.err_z);
        }
      }

      if (!y0_src.empty()) {
        int color = colors[color_idx % 8];
        color_idx++;

        TGraphErrors* g2y =
            new TGraphErrors(y0_src.size(), y0_src.data(), y_det.data(), nullptr, y_err.data());
        g2y->SetMarkerStyle(20);
        g2y->SetMarkerColor(color);
        g2y->SetLineColor(color);
        mg_y->Add(g2y, "PL");

        TGraphErrors* g2z =
            new TGraphErrors(y0_src.size(), y0_src.data(), z_det.data(), nullptr, z_err.data());
        g2z->SetMarkerStyle(20);
        g2z->SetMarkerColor(color);
        g2z->SetLineColor(color);
        mg_z->Add(g2z, "PL");

        leg->AddEntry(g2y, ("x0 = " + format_double(x0_target)).c_str(), "PL");
      }
    }

    c2d->cd(1);
    gPad->SetLeftMargin(0.15);
    mg_y->Draw("A");
    leg->Draw();

    c2d->cd(2);
    gPad->SetLeftMargin(0.15);
    mg_z->Draw("A");
    leg->Draw();

    std::string filename_2d = (output_dir
                               / ("beam_scan_2D_x1_" + format_double(found_x1, 2) + "_x2_"
                                  + format_double(found_x2, 2) + ".png"))
                                  .string();
    c2d->SaveAs(filename_2d.c_str());
    delete c2d;
    std::cout << "2D Slices saved to " << filename_2d << std::endl;
  }

  return 0;
}