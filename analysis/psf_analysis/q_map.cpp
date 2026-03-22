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
 * q_map — Mappa 2D della funzione di qualità Q(x1, x2)
 *
 * Per ogni configurazione (x1, x2) presente in psf_data.root, calcola
 * Q = sum_{y0} chi^2(y0, x1, x2) tramite compute_Q() e produce una
 * mappa TH2D.
 *
 * Uso:
 *   q_map [--psf <path>] [--config <path>] [--output <path>]
 *         [--y0-min <val>] [--y0-max <val>] [--dy0 <val>]
 *         [--L <val>] [--dt <val>]
 *         [--log] [--norm] [--tsv <path>]
 */

#include "psf_interpolator.hpp"

#include <nlohmann/json.hpp>

#include <TAxis.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TGaxis.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMarker.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// Parsing CLI

struct CliConfig {
  std::string psf_path    = "output/psf/psf_data.root";
  std::string config_path = "config/config.json";
  std::string output_path = "output/psf_analysis/q_map.png";
  std::string tsv_path    = "";

  double y0_min = -1.0;
  double y0_max = -1.0;
  double dy0    = -1.0;
  double L      = 10.0;
  double dt     = 0.1;

  bool log_scale = false;
  bool normalize = false;
};

static CliConfig parse_args(int argc, char** argv) {
  CliConfig cfg;
  for (int i = 1; i < argc; ++i) {
    std::string key = argv[i];
    auto next       = [&]() -> std::string {
      if (i + 1 >= argc) {
        std::cerr << "Argomento mancante dopo " << key << "\n";
        std::exit(1);
      }
      return argv[++i];
    };
    if (key == "--psf")
      cfg.psf_path = next();
    else if (key == "--config")
      cfg.config_path = next();
    else if (key == "--output")
      cfg.output_path = next();
    else if (key == "--tsv")
      cfg.tsv_path = next();
    else if (key == "--y0-min")
      cfg.y0_min = std::stod(next());
    else if (key == "--y0-max")
      cfg.y0_max = std::stod(next());
    else if (key == "--dy0")
      cfg.dy0 = std::stod(next());
    else if (key == "--L")
      cfg.L = std::stod(next());
    else if (key == "--dt")
      cfg.dt = std::stod(next());
    else if (key == "--log")
      cfg.log_scale = true;
    else if (key == "--norm")
      cfg.normalize = true;
    else {
      std::cerr << "Opzione sconosciuta: " << key << "\n";
      std::exit(1);
    }
  }
  return cfg;
}

// Identico a trace_viewer per uniformità visiva.
static void apply_style() {
  gStyle->Reset();
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleFont(42, "");
  gStyle->SetStatFont(42);

  gStyle->SetTextSize(0.040);
  gStyle->SetLabelSize(0.038, "XYZ");
  gStyle->SetTitleSize(0.044, "XYZ");
  gStyle->SetTitleSize(0.046, "");

  gStyle->SetTitleOffset(1.55, "Y");
  gStyle->SetTitleOffset(1.20, "X");
  gStyle->SetTitleOffset(1.40, "Z");

  gStyle->SetTickLength(0.018, "X");
  gStyle->SetTickLength(0.018, "Y");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetGridColor(kGray + 1);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0); // Titolo gestito manualmente con TLatex

  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);

  gStyle->SetNdivisions(505, "X");
  gStyle->SetNdivisions(505, "Y");

  // kBird è la palette publication-standard per mappe 2D
  gStyle->SetPalette(kBird);
  gStyle->SetNumberContours(255);
}

// Helpers
static std::string fmt(double v, int n = 1) {
  std::ostringstream o;
  o << std::fixed << std::setprecision(n) << v;
  return o.str();
}

// Restituisce il colore ROOT dalla palette corrente per t ∈ [0, 1].
static Int_t palette_color(double t) {
  int n   = gStyle->GetNumberOfColors();
  int idx = static_cast<int>(std::clamp(t, 0.0, 1.0) * (n - 1));
  return gStyle->GetColorPalette(idx);
}

// main
int main(int argc, char** argv) {
  CliConfig cli = parse_args(argc, argv);

  // 1. Legge config.json
  using json = nlohmann::json;
  std::ifstream f_cfg(cli.config_path);
  if (!f_cfg.is_open()) {
    std::cerr << "Errore: impossibile aprire " << cli.config_path << "\n";
    return 1;
  }
  json jcfg;
  f_cfg >> jcfg;

  const double dx   = jcfg["dx"];
  const double r1   = jcfg["r1"];
  const double h1   = jcfg["h1"];
  const double r2   = jcfg["r2"];
  const double h2   = jcfg["h2"];
  const double xmin = jcfg["x_min"];
  const double xmax = jcfg["x_max"];

  const double x1_lo = xmin - r1 + h1;
  const double x1_hi = xmax - h2 - 3.0 - r1;
  const double x2_lo = x1_lo + r1 + r2 + 3.0;
  const double x2_hi = xmax + r2 - h2;

  // 2. Parametri di campionamento y0
  riptide::QConfig qcfg;
  qcfg.y0_min       = (cli.y0_min >= 0.0) ? cli.y0_min : 0.0;
  qcfg.y0_max       = (cli.y0_max >= 0.0) ? cli.y0_max : 10.0;
  qcfg.dy0          = (cli.dy0 > 0.0) ? cli.dy0 : 1.0;
  qcfg.trace_L      = cli.L;
  qcfg.trace_dt     = cli.dt;
  qcfg.fit_max_iter = 20;
  qcfg.fit_tol      = 1e-8;

  int n_y0 = static_cast<int>(std::round((qcfg.y0_max - qcfg.y0_min) / qcfg.dy0)) + 1;
  std::cout << "Campionamento y0: [" << qcfg.y0_min << ", " << qcfg.y0_max << "] passo " << qcfg.dy0
            << " mm  →  " << n_y0 << " tracce per configurazione\n";

  // 3. Carica PSF database
  riptide::PSFDatabase db;
  try {
    db = riptide::load_psf_database(cli.psf_path);
  } catch (const std::exception& e) {
    std::cerr << "Errore: " << e.what() << "\n";
    return 1;
  }
  std::cout << "Configurazioni nel database: " << db.size() << "\n";

  // 4. Calcola Q per ogni configurazione
  struct QEntry {
    double x1, x2, Q;
    int n_traces, n_failed;
  };
  std::vector<QEntry> entries;
  entries.reserve(db.size());

  int config_idx = 0;
  int total_cfgs = static_cast<int>(db.size());
  int print_step = std::max(1, total_cfgs / 20);

  for (const auto& [cfg, _] : db) {
    riptide::QResult res = riptide::compute_Q(cfg, db, qcfg);

    double Q_val = res.Q;
    if (cli.normalize && res.n_traces > 0)
      Q_val /= static_cast<double>(res.n_traces);

    entries.push_back({cfg.x1, cfg.x2, Q_val, res.n_traces, res.n_failed});

    if ((config_idx + 1) % print_step == 0 || config_idx + 1 == total_cfgs) {
      std::cout << "  [" << std::setw(3) << (config_idx + 1) << "/" << total_cfgs << "]"
                << "  x1=" << fmt(cfg.x1) << "  x2=" << fmt(cfg.x2) << "  Q=" << std::scientific
                << std::setprecision(3) << Q_val << std::defaultfloat << "\n";
    }
    ++config_idx;
  }

  // 5. Minimo di Q
  auto it_best = std::min_element(entries.begin(), entries.end(),
                                  [](const QEntry& a, const QEntry& b) { return a.Q < b.Q; });
  std::cout << "\n★  Configurazione ottimale:\n"
            << "   x1 = " << fmt(it_best->x1, 2) << " mm\n"
            << "   x2 = " << fmt(it_best->x2, 2) << " mm\n"
            << "   Q  = " << it_best->Q << "\n";

  // 6. Esportazione TSV opzionale
  if (!cli.tsv_path.empty()) {
    std::filesystem::create_directories(std::filesystem::path(cli.tsv_path).parent_path());
    std::ofstream tsv(cli.tsv_path);
    tsv << "x1\tx2\tQ\tn_traces\tn_failed\n";
    for (const auto& e : entries)
      tsv << e.x1 << "\t" << e.x2 << "\t" << e.Q << "\t" << e.n_traces << "\t" << e.n_failed
          << "\n";
    std::cout << "TSV salvato in: " << cli.tsv_path << "\n";
  }

  // Sezione grafica
  apply_style();

  // 7. Costruisce il TH2D
  int bins_x = std::max(1, static_cast<int>(std::ceil((x1_hi - x1_lo) / dx)) + 1);
  int bins_y = std::max(1, static_cast<int>(std::ceil((x2_hi - x2_lo) / dx)) + 1);
  double hx  = dx / 2.0;

  TH2D h_Q("h_Q", "", bins_x, x1_lo - hx, x1_hi + hx, bins_y, x2_lo - hx, x2_hi + hx);

  for (const auto& e : entries) {
    double val = e.Q;
    if (val <= 0.0)
      val = std::numeric_limits<double>::min();
    h_Q.Fill(e.x1, e.x2, val);
  }

  h_Q.GetXaxis()->SetTitle("x_{1}  [mm]  (lente 75 mm)");
  h_Q.GetYaxis()->SetTitle("x_{2}  [mm]  (lente 60 mm)");
  h_Q.GetXaxis()->CenterTitle(true);
  h_Q.GetYaxis()->CenterTitle(true);
  h_Q.GetXaxis()->SetTitleOffset(1.20);
  h_Q.GetYaxis()->SetTitleOffset(1.55);
  h_Q.GetXaxis()->SetNdivisions(506);
  h_Q.GetYaxis()->SetNdivisions(505);

  // Range colori: percentili 5–95
  {
    std::vector<double> q_vals;
    q_vals.reserve(entries.size());
    for (const auto& e : entries)
      if (e.Q > 0.0)
        q_vals.push_back(e.Q);
    std::sort(q_vals.begin(), q_vals.end());
    if (!q_vals.empty()) {
      size_t N    = q_vals.size();
      double z_lo = q_vals[static_cast<size_t>(0.05 * N)];
      double z_hi = q_vals[static_cast<size_t>(0.95 * N)];
      h_Q.SetMinimum(z_lo);
      h_Q.SetMaximum(z_hi);
    }
  }

  // 8. Canvas e layout — identico a trace_viewer
  TCanvas* canvas = new TCanvas("canvas", "Q map", 1200, 1000);
  canvas->SetLeftMargin(0.0);
  canvas->SetRightMargin(0.0);
  canvas->SetTopMargin(0.0);
  canvas->SetBottomMargin(0.0);

  TPad* pad_plot = new TPad("pad_plot", "", 0.00, 0.12, 0.88, 1.00);
  TPad* pad_cb   = new TPad("pad_cb", "", 0.88, 0.12, 0.96, 1.00);
  TPad* pad_info = new TPad("pad_info", "", 0.00, 0.00, 1.00, 0.12);

  // pad_plot: margini identici a trace_viewer
  pad_plot->SetLeftMargin(0.13);
  pad_plot->SetRightMargin(0.015);
  pad_plot->SetTopMargin(0.11);
  pad_plot->SetBottomMargin(0.13);
  pad_plot->SetGridx();
  pad_plot->SetGridy();
  pad_plot->SetFrameLineWidth(2);
  if (cli.log_scale)
    pad_plot->SetLogz();

  // pad_cb: color bar manuale — NON usa COLZ, nessuna palette automatica
  pad_cb->SetLeftMargin(0.25);
  pad_cb->SetRightMargin(0.30);
  pad_cb->SetTopMargin(0.11);    // allineato con pad_plot
  pad_cb->SetBottomMargin(0.13); // allineato con pad_plot

  pad_info->SetLeftMargin(0.02);
  pad_info->SetRightMargin(0.02);
  pad_info->SetTopMargin(0.05);
  pad_info->SetBottomMargin(0.05);

  canvas->cd();
  pad_plot->Draw();
  pad_cb->Draw();
  pad_info->Draw();

  // 9. Plot principale: disegna con "COL" (non "COLZ") per evitare
  pad_plot->cd();
  h_Q.Draw("COL");

  // Marker stella sul minimo
  {
    TMarker* mk = new TMarker(it_best->x1, it_best->x2, 29);
    mk->SetMarkerColor(kRed);
    mk->SetMarkerSize(2.8);
    mk->Draw("same");

    TLatex lbl;
    lbl.SetTextFont(42);
    lbl.SetTextSize(0.032);
    lbl.SetTextColor(kRed + 1);
    lbl.SetTextAlign(12);
    std::string lbl_str = " #bf{min} (" + fmt(it_best->x1, 1) + ", " + fmt(it_best->x2, 1) + ")";
    lbl.DrawLatex(it_best->x1, it_best->x2, lbl_str.c_str());
  }

  // Titolo manuale nel margine superiore del pad_plot (NDC).
  {
    TLatex title_lat;
    title_lat.SetNDC();
    title_lat.SetTextFont(42);
    title_lat.SetTextSize(0.048);
    title_lat.SetTextAlign(22);
    title_lat.SetTextColor(kBlack);
    const std::string z_label = cli.normalize ? "Q / n_{tracce}" : "Q(x_{1}, x_{2})";
    // Titolo centrato tra i margini sinistro e destro del pad
    title_lat.DrawLatex(0.535, 0.953,
                        ("Mappa della funzione di qualit#grave{a}  " + z_label).c_str());
  }

  pad_plot->RedrawAxis();

  // 10. Color bar manuale nel pad_cb — tecnica identica a trace_viewer
  pad_cb->cd();
  pad_cb->Range(0.0, 0.0, 1.0, 1.0);

  const int NB = 255;
  double cb_x0 = pad_cb->GetLeftMargin();
  double cb_x1 = 1.0 - pad_cb->GetRightMargin();
  // Allinea i limiti verticali della barra con l'area utile del pad_plot
  // espressi nelle coordinate NDC del pad_cb stesso.
  double cb_y0 = pad_cb->GetBottomMargin();
  double cb_y1 = 1.0 - pad_cb->GetTopMargin();

  for (int i = 0; i < NB; ++i) {
    double f0  = static_cast<double>(i) / NB;
    double f1  = static_cast<double>(i + 1) / NB;
    double yb0 = cb_y0 + f0 * (cb_y1 - cb_y0);
    double yb1 = cb_y0 + f1 * (cb_y1 - cb_y0);

    TBox* box = new TBox(cb_x0, yb0, cb_x1, yb1);
    box->SetFillColor(gStyle->GetColorPalette(i));
    box->SetLineWidth(0);
    box->Draw();
  }

  // Asse della color bar con range fisico di Q
  double cb_zmin = h_Q.GetMinimum();
  double cb_zmax = h_Q.GetMaximum();
  TGaxis* cb_axis =
      new TGaxis(cb_x1, cb_y0, cb_x1, cb_y1, cb_zmin, cb_zmax, 505, cli.log_scale ? "+LG" : "+L");
  cb_axis->SetLabelFont(42);
  cb_axis->SetLabelSize(0.18);
  cb_axis->SetTickSize(0.35);
  cb_axis->SetLabelOffset(0.03);

  const std::string z_label_cb = cli.normalize ? "Q / n_{tracce}" : "Q(x_{1}, x_{2})";
  cb_axis->SetTitle(z_label_cb.c_str());
  cb_axis->SetTitleFont(42);
  cb_axis->SetTitleSize(0.20);
  cb_axis->SetTitleOffset(0.55);
  cb_axis->Draw();

  // 11. Info panel — tre colonne, identico a trace_viewer
  pad_info->cd();
  pad_info->SetFillColor(TColor::GetColor(245, 245, 248));

  TLine* sep_line = new TLine(0.0, 0.97, 1.0, 0.97);
  sep_line->SetNDC(true);
  sep_line->SetLineColor(kGray + 1);
  sep_line->SetLineWidth(1);
  sep_line->Draw();

  TLatex info;
  info.SetNDC(true);
  info.SetTextFont(42);

  // Coordinate NDC del pad_info (altezza 12% del canvas)
  const double col1 = 0.03, col2 = 0.54;
  const double hdr = 0.82, row1 = 0.52, row2 = 0.18;

  // Header colonne
  info.SetTextSize(0.13);
  info.SetTextColor(kGray + 2);
  info.SetTextAlign(12);
  info.DrawLatex(col1, hdr, "CAMPIONAMENTO y_{0}");
  info.DrawLatex(col2, hdr, "OTTIMO");

  // Dati
  info.SetTextSize(0.20);
  info.SetTextColor(kBlack);

  // Colonna 1: campionamento y0 e traccia
  {
    info.DrawLatex(col1, row1,
                   ("y_{0} #in [" + fmt(qcfg.y0_min, 1) + ", " + fmt(qcfg.y0_max, 1)
                    + "] mm   #Deltay_{0} = " + fmt(qcfg.dy0, 1) + " mm")
                       .c_str());
    info.DrawLatex(col1, row2,
                   ("L = " + fmt(qcfg.trace_L, 1) + " mm   #Deltat = " + fmt(qcfg.trace_dt, 2)
                    + " mm   (" + std::to_string(n_y0) + " tracce/config)")
                       .c_str());
  }

  // Colonna 2: minimo di Q
  {
    info.DrawLatex(col2, row1,
                   ("#bf{x_{1}^{*}} = " + fmt(it_best->x1, 1)
                    + " mm,   #bf{x_{2}^{*}} = " + fmt(it_best->x2, 1) + " mm")
                       .c_str());
    std::ostringstream q_str;
    q_str << std::scientific << std::setprecision(3) << it_best->Q;
    info.DrawLatex(col2, row2, ("#bf{Q_{min}} = " + q_str.str()).c_str());
  }

  // 12. Salva
  canvas->Update();

  std::filesystem::create_directories(std::filesystem::path(cli.output_path).parent_path());
  canvas->Print(cli.output_path.c_str());
  std::cout << "\nMappa salvata in: " << cli.output_path << "\n";

  delete canvas;
  return 0;
}