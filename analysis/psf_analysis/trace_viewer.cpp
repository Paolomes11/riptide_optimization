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

/*
 * trace_viewer -- visualizza la traccia media sul detector con ellissi di covarianza
 *
 * Uso:
 *   trace_viewer --x1 <val> --x2 <val> --y0 <val> [--psf <path>] [--output <path>]
 *                [--dt <val>] [--L <val>] [--sigma <val>]
 */

#include "psf_interpolator.hpp"

#include <TCanvas.h>
#include <TEllipse.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TStyle.h>

#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

//  Parsing CLI 

struct Config {
  double x1            = -1.0;
  double x2            = -1.0;
  double y0            = -1.0;
  double dt            = 0.1;
  double L             = 10.0;
  double sigma_scale   = 1.0;
  std::string psf_path = "output/psf/psf_data.root";
  std::string output   = "";
};

Config parse_args(int argc, char** argv) {
  Config cfg;
  for (int i = 1; i < argc - 1; ++i) {
    std::string key = argv[i];
    std::string val = argv[i + 1];
    if (key == "--x1") {
      cfg.x1 = std::stod(val);
      ++i;
    } else if (key == "--x2") {
      cfg.x2 = std::stod(val);
      ++i;
    } else if (key == "--y0") {
      cfg.y0 = std::stod(val);
      ++i;
    } else if (key == "--dt") {
      cfg.dt = std::stod(val);
      ++i;
    } else if (key == "--L") {
      cfg.L = std::stod(val);
      ++i;
    } else if (key == "--sigma") {
      cfg.sigma_scale = std::stod(val);
      ++i;
    } else if (key == "--psf") {
      cfg.psf_path = val;
      ++i;
    } else if (key == "--output") {
      cfg.output = val;
      ++i;
    }
  }
  return cfg;
}

//  Calcolo ellisse di covarianza 

struct Ellipse {
  double cy, cz;
  double a, b;
  double theta_deg;
};

Ellipse cov_to_ellipse(double cy, double cz, const riptide::Cov2& cov, double n_sigma = 1.0) {
  double trace   = cov.yy + cov.zz;
  double det     = cov.yy * cov.zz - cov.yz * cov.yz;
  double disc    = std::sqrt(std::max(0.0, trace * trace / 4.0 - det));
  double lambda1 = trace / 2.0 + disc;
  double lambda2 = trace / 2.0 - disc;

  double a = n_sigma * std::sqrt(std::max(0.0, lambda1));
  double b = n_sigma * std::sqrt(std::max(0.0, lambda2));

  double theta_rad = 0.0;
  if (std::abs(cov.yz) > 1e-12 || std::abs(cov.yy - lambda1) > 1e-12)
    theta_rad = std::atan2(lambda1 - cov.yy, cov.yz);

  return {cy, cz, a, b, theta_rad * 180.0 / M_PI};
}

//  main 

int main(int argc, char** argv) {
  Config cfg = parse_args(argc, argv);

  if (cfg.x1 < 0 || cfg.x2 < 0 || cfg.y0 < 0) {
    std::cerr << "Uso: trace_viewer --x1 <val> --x2 <val> --y0 <val>\n"
              << "     [--psf <path>]     default: output/psf/psf_data.root\n"
              << "     [--output <path>]  default: output/psf_analysis/trace_*.png\n"
              << "     [--dt <val>]       step traccia in mm, default: 0.1\n"
              << "     [--L <val>]        lunghezza traccia in mm, default: 10.0\n"
              << "     [--sigma <val>]    scala ellissi, default: 1.0\n";
    return 1;
  }

  // Carica database PSF
  riptide::PSFDatabase db;
  try {
    db = riptide::load_psf_database(cfg.psf_path);
  } catch (const std::exception& e) {
    std::cerr << "Errore: " << e.what() << "\n";
    return 1;
  }

  // Trova configurazione piu' vicina e costruisce la traccia
  riptide::LensConfig lens{cfg.x1, cfg.x2};
  riptide::LensConfig actual = riptide::find_nearest_config(lens, db);

  std::vector<riptide::TracePoint> trace;
  try {
    trace = riptide::build_trace(cfg.y0, actual, db, cfg.L, cfg.dt);
  } catch (const std::exception& e) {
    std::cerr << "Errore nella costruzione della traccia: " << e.what() << "\n";
    return 1;
  }

  if (trace.empty()) {
    std::cerr << "Errore: traccia vuota\n";
    return 1;
  }

  std::cout << "Traccia costruita: " << trace.size() << " punti\n";
  std::cout << "  Richiesta: x1=" << cfg.x1 << " mm, x2=" << cfg.x2 << " mm\n";
  std::cout << "  Usata:     x1=" << actual.x1 << " mm, x2=" << actual.x2 << " mm, y0=" << cfg.y0
            << " mm\n";

  //  Calcola range degli assi 
  double y_min = 1e9, y_max = -1e9;
  double z_min = 1e9, z_max = -1e9;
  double max_sigma_y = 0, max_sigma_z = 0;

  for (const auto& pt : trace) {
    y_min       = std::min(y_min, pt.mu_y);
    y_max       = std::max(y_max, pt.mu_y);
    z_min       = std::min(z_min, pt.mu_z);
    z_max       = std::max(z_max, pt.mu_z);
    max_sigma_y = std::max(max_sigma_y, std::sqrt(pt.cov.yy));
    max_sigma_z = std::max(max_sigma_z, std::sqrt(pt.cov.zz));
  }

  // Usa la sigma maggiore tra y e z come unita' di margine per entrambi gli assi.
  // Cosi' se sigma_z e' quasi zero il range z rimane comunque aperto quanto sigma_y,
  // e le ellissi sono visibili con le proporzioni fisiche corrette.
  // double max_sigma = std::max({max_sigma_y, max_sigma_z, 1e-4});
  double margin_y    = std::max(2.0 * max_sigma_y * cfg.sigma_scale, 0.001);
  double margin_z    = std::max(2.0 * max_sigma_z * cfg.sigma_scale, 0.001);

  double py_min = y_min - margin_y;
  double py_max = y_max + margin_y;

  // Range z centrato e simmetrico attorno al valore medio della traccia
  double z_center = (z_min + z_max) / 2.0;
  // double z_half   = std::max(std::abs(z_max - z_center), margin_z);
  double pz_min   = z_center - margin_z;
  double pz_max   = z_center + margin_z;

  //  Canvas 
  gStyle->SetOptStat(0);
  gStyle->SetPadGridX(true);
  gStyle->SetPadGridY(true);
  gStyle->SetGridColor(17);

  TCanvas* canvas = new TCanvas("canvas", "Trace Viewer", 1000, 800);
  canvas->SetLeftMargin(0.13);
  canvas->SetRightMargin(0.05);
  canvas->SetTopMargin(0.12);
  canvas->SetBottomMargin(0.12);

  TH2D* frame = new TH2D("frame", "", 100, py_min, py_max, 100, pz_min, pz_max);
  frame->GetXaxis()->SetTitle("y_{det} [mm]");
  frame->GetYaxis()->SetTitle("z_{det} [mm]");
  frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetTitleOffset(1.5);
  frame->GetXaxis()->SetTitleSize(0.048);
  frame->GetYaxis()->SetTitleSize(0.048);
  frame->Draw();

  // Ellissi di covarianza
  std::vector<TEllipse*> ellipses;
  for (const auto& pt : trace) {
    auto el = cov_to_ellipse(pt.mu_y, pt.mu_z, pt.cov, cfg.sigma_scale);
    if (el.a < 1e-6 || el.b < 1e-6)
      continue;
    TEllipse* e = new TEllipse(el.cy, el.cz, el.a, el.b, 0, 360, el.theta_deg);
    e->SetFillColorAlpha(kBlue - 9, 0.30);
    e->SetLineColor(kBlue - 7);
    e->SetLineWidth(1);
    e->Draw("same");
    ellipses.push_back(e);
  }

  // Centroidi
  std::vector<double> gy, gz;
  gy.reserve(trace.size());
  gz.reserve(trace.size());
  for (const auto& pt : trace) {
    gy.push_back(pt.mu_y);
    gz.push_back(pt.mu_z);
  }
  TGraph* g = new TGraph(static_cast<int>(gy.size()), gy.data(), gz.data());
  g->SetLineColor(kRed);
  g->SetLineWidth(2);
  g->SetMarkerColor(kRed);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.5);
  g->Draw("LP same");

  // Titolo
  TLatex lat;
  lat.SetNDC();
  lat.SetTextAlign(22);
  lat.SetTextSize(0.038);
  std::ostringstream title;
  title << std::fixed << std::setprecision(1) << "Traccia sul detector  "
        << "x_{1}=" << actual.x1 << " mm, "
        << "x_{2}=" << actual.x2 << " mm, "
        << "y_{0}=" << cfg.y0 << " mm";
  lat.DrawLatex(0.54, 0.95, title.str().c_str());

  // Legenda
  lat.SetTextSize(0.030);
  lat.SetTextAlign(12);
  lat.SetTextColor(kRed);
  lat.DrawLatex(0.15, 0.90, "Centroidi #mu_{y}, #mu_{z}");
  lat.SetTextColor(kBlue - 7);
  std::ostringstream leg;
  leg << std::fixed << std::setprecision(1) << "Ellissi " << cfg.sigma_scale << "#sigma";
  lat.DrawLatex(0.50, 0.90, leg.str().c_str());

  // Info sigma
  lat.SetTextColor(kBlack);
  lat.SetTextSize(0.026);
  lat.SetTextAlign(12);
  std::ostringstream sinfo;
  sinfo << std::fixed << std::setprecision(4) << "#sigma_{y,max}=" << max_sigma_y << " mm    "
        << "#sigma_{z,max}=" << max_sigma_z << " mm";
  lat.DrawLatex(0.15, 0.04, sinfo.str().c_str());

  canvas->Update();

  // Salva
  if (!cfg.output.empty()) {
    std::filesystem::create_directories(std::filesystem::path(cfg.output).parent_path());
    canvas->SaveAs(cfg.output.c_str());
    std::cout << "Plot salvato in: " << cfg.output << "\n";
  } else {
    std::ostringstream auto_name;
    auto_name << "output/psf_analysis/trace"
              << "_x1_" << std::fixed << std::setprecision(1) << actual.x1 << "_x2_" << actual.x2
              << "_y0_" << cfg.y0 << ".png";
    std::filesystem::create_directories("output/psf_analysis");
    canvas->SaveAs(auto_name.str().c_str());
    std::cout << "Plot salvato in: " << auto_name.str() << "\n";
  }

  delete g;
  for (auto* e : ellipses)
    delete e;
  delete frame;
  delete canvas;

  return 0;
}
