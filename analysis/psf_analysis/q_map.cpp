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
 * TEMPORAL UNFOLDING (attivo per default)
 *
 * Per default il calcolo di Q utilizza il temporal unfolding: a ogni punto i
 * della traccia viene aggiunto un offset artificiale in z:
 *
 *     z̃_i = μ_{z,i} + i · δz_unfold
 *
 * con δz_unfold = trace_L / (N-1) per default (offset totale = trace_L mm).
 * Questo rende visibile qualsiasi non-linearità (curvatura, ripiegamento)
 * che sarebbe altrimenti mascherata dalla simmetria z≈0.
 *
 * Opzioni:
 *   --unfold-dz <val>   passo di srotolamento fisso [mm/passo] (0 = automatico)
 *   --no-unfold         disabilita il temporal unfolding (usa z fisico)
 *
 * Gestione celle invalide:
 *   Celle grigie = configurazioni in cui meno del 75% delle tracce è on_detector.
 *
 * Uso:
 *   q_map [--psf <path>] [--config <path>] [--output <path>]
 *         [--y0-min <val>] [--y0-max <val>] [--dy0 <val>]
 *         [--L <val>] [--dt <val>]
 *         [--unfold-dz <val>]
 *         [--no-unfold]
 *         [--point-frac <val>]   (soglia punti on_detector per traccia, default 0.75)
 *         [--trace-frac <val>]   (soglia tracce valide per config, default 0.75)
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
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
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

  double point_frac = 0.75;
  double trace_frac = 0.75;

  // Temporal unfolding
  double unfold_dz = 0.0;   // 0 = automatico
  bool no_unfold   = false; // se true, disabilita l'unfolding

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
    else if (key == "--point-frac")
      cfg.point_frac = std::stod(next());
    else if (key == "--trace-frac")
      cfg.trace_frac = std::stod(next());
    else if (key == "--unfold-dz")
      cfg.unfold_dz = std::stod(next());
    else if (key == "--no-unfold")
      cfg.no_unfold = true;
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

// Stile ROOT
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
  gStyle->SetOptTitle(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetNdivisions(505, "Y");
  gStyle->SetPalette(kBird);
  gStyle->SetNumberContours(255);
}

static std::string fmt(double v, int n = 1) {
  std::ostringstream o;
  o << std::fixed << std::setprecision(n) << v;
  return o.str();
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

  // 2. Parametri QConfig
  riptide::QConfig qcfg;
  qcfg.y0_min               = (cli.y0_min >= 0.0) ? cli.y0_min : 0.0;
  qcfg.y0_max               = (cli.y0_max >= 0.0) ? cli.y0_max : 10.0;
  qcfg.dy0                  = (cli.dy0 > 0.0) ? cli.dy0 : 1.0;
  qcfg.trace_L              = cli.L;
  qcfg.trace_dt             = cli.dt;
  qcfg.fit_max_iter         = 20;
  qcfg.fit_tol              = 1e-8;
  qcfg.point_valid_fraction = cli.point_frac;
  qcfg.trace_valid_fraction = cli.trace_frac;

  // Temporal unfolding
  qcfg.apply_temporal_unfolding = !cli.no_unfold;
  qcfg.z_unfold_step            = cli.unfold_dz;

  int n_y0 = static_cast<int>(std::round((qcfg.y0_max - qcfg.y0_min) / qcfg.dy0)) + 1;

  // Calcola il dz effettivo per la stampa diagnostica
  int n_trace_pts     = static_cast<int>(std::round(qcfg.trace_L / qcfg.trace_dt)) + 1;
  double dz_effective = (qcfg.z_unfold_step > 0.0)
                          ? qcfg.z_unfold_step
                          : qcfg.trace_L / static_cast<double>(n_trace_pts - 1);

  std::cout << "Campionamento y0: [" << qcfg.y0_min << ", " << qcfg.y0_max << "] passo " << qcfg.dy0
            << " mm  →  " << n_y0 << " tracce per configurazione\n";
  std::cout << "Soglie validità: punti " << cli.point_frac * 100 << "% | tracce "
            << cli.trace_frac * 100 << "%\n";
  if (qcfg.apply_temporal_unfolding) {
    std::cout << "Temporal unfolding: ATTIVO  δz = " << fmt(dz_effective, 4)
              << " mm/passo  (offset totale ≈ " << fmt(dz_effective * (n_trace_pts - 1), 2)
              << " mm)\n";
  } else {
    std::cout << "Temporal unfolding: DISATTIVATO\n";
  }

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
    double x1, x2;
    double Q;
    int n_traces;
    int n_failed;
    int n_invalid;
    bool config_valid;
  };
  std::vector<QEntry> entries;
  entries.reserve(db.size());

  int config_idx    = 0;
  int total_cfgs    = static_cast<int>(db.size());
  int print_step    = std::max(1, total_cfgs / 20);
  int n_cfg_invalid = 0;

  for (const auto& [cfg, _] : db) {
    riptide::QResult res;
    try {
      res = riptide::compute_Q(cfg, db, qcfg);
    } catch (const std::exception& e) {
      std::cerr << "  [WARN] compute_Q failed x1=" << cfg.x1 << " x2=" << cfg.x2 << ": " << e.what()
                << "\n";
      entries.push_back({cfg.x1, cfg.x2, 0.0, 0, 0, 0, false});
      ++n_cfg_invalid;
      ++config_idx;
      continue;
    }

    double Q_val = res.Q;
    if (cli.normalize && res.n_traces > 0)
      Q_val /= static_cast<double>(res.n_traces);

    if (!res.config_valid)
      ++n_cfg_invalid;

    entries.push_back(
        {cfg.x1, cfg.x2, Q_val, res.n_traces, res.n_failed, res.n_invalid, res.config_valid});

    if ((config_idx + 1) % print_step == 0 || config_idx + 1 == total_cfgs) {
      std::cout << "  [" << std::setw(3) << config_idx + 1 << "/" << total_cfgs << "]"
                << "  x1=" << fmt(cfg.x1) << "  x2=" << fmt(cfg.x2);
      if (res.config_valid)
        std::cout << "  Q=" << std::scientific << std::setprecision(3) << Q_val;
      else
        std::cout << "  [INVALIDA: " << res.n_traces << "/" << n_y0 << " tracce valide]";
      std::cout << std::defaultfloat << "\n";
    }
    ++config_idx;
  }

  std::cout << "\nConfigurazioni valide:   " << total_cfgs - n_cfg_invalid << "\n";
  std::cout << "Configurazioni invalide: " << n_cfg_invalid << "\n";

  // 5. Minimo di Q
  auto it_best =
      std::min_element(entries.begin(), entries.end(), [](const QEntry& a, const QEntry& b) {
        if (!a.config_valid)
          return false;
        if (!b.config_valid)
          return true;
        return a.Q < b.Q;
      });

  if (!it_best->config_valid) {
    std::cerr << "ERRORE: nessuna configurazione valida trovata!\n";
    return 1;
  }
  std::cout << "\n★  Configurazione ottimale:\n"
            << "   x1 = " << fmt(it_best->x1, 2) << " mm\n"
            << "   x2 = " << fmt(it_best->x2, 2) << " mm\n"
            << "   Q  = " << it_best->Q << "\n";

  // 6. Esportazione TSV opzionale
  if (!cli.tsv_path.empty()) {
    std::filesystem::create_directories(std::filesystem::path(cli.tsv_path).parent_path());
    std::ofstream tsv(cli.tsv_path);
    tsv << "x1\tx2\tQ\tn_traces\tn_failed\tn_invalid\tconfig_valid\n";
    for (const auto& e : entries)
      tsv << e.x1 << "\t" << e.x2 << "\t" << (e.config_valid ? std::to_string(e.Q) : "NaN") << "\t"
          << e.n_traces << "\t" << e.n_failed << "\t" << e.n_invalid << "\t"
          << (e.config_valid ? "1" : "0") << "\n";
    std::cout << "TSV salvato in: " << cli.tsv_path << "\n";
  }

  // Sezione grafica
  apply_style();

  // 7. TH2D valide + TH2D invalide
  int bins_x = std::max(1, static_cast<int>(std::round((x1_hi - x1_lo) / dx)) + 1);
  int bins_y = std::max(1, static_cast<int>(std::round((x2_hi - x2_lo) / dx)) + 1);
  double hx  = dx / 2.0;

  double ax_x1_lo = x1_lo - hx, ax_x1_hi = x1_hi + hx;
  double ax_x2_lo = x2_lo - hx, ax_x2_hi = x2_hi + hx;

  TH2D h_Q("h_Q", "", bins_x, ax_x1_lo, ax_x1_hi, bins_y, ax_x2_lo, ax_x2_hi);
  TH2D h_inv("h_inv", "", bins_x, ax_x1_lo, ax_x1_hi, bins_y, ax_x2_lo, ax_x2_hi);
  h_Q.Reset();
  h_inv.Reset();

  for (const auto& e : entries) {
    int bx = h_Q.GetXaxis()->FindFixBin(e.x1);
    int by = h_Q.GetYaxis()->FindFixBin(e.x2);
    if (e.config_valid)
      h_Q.SetBinContent(bx, by, e.Q);
    else
      h_inv.SetBinContent(bx, by, 1.0);
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
      if (e.config_valid && e.Q > 0.0)
        q_vals.push_back(e.Q);
    std::sort(q_vals.begin(), q_vals.end());
    if (!q_vals.empty()) {
      size_t N = q_vals.size();
      h_Q.SetMinimum(q_vals[static_cast<size_t>(0.05 * N)]);
      h_Q.SetMaximum(q_vals[static_cast<size_t>(0.95 * N)]);
    }
  }

  Int_t invalid_color = TColor::GetColor(80, 80, 80);
  h_inv.SetFillColor(invalid_color);
  h_inv.SetLineColor(invalid_color);

  // 8. Canvas e layout
  TCanvas* canvas = new TCanvas("canvas", "Q map", 1200, 1000);
  canvas->SetLeftMargin(0.0);
  canvas->SetRightMargin(0.0);
  canvas->SetTopMargin(0.0);
  canvas->SetBottomMargin(0.0);

  TPad* pad_plot = new TPad("pad_plot", "", 0.00, 0.12, 0.88, 1.00);
  TPad* pad_cb   = new TPad("pad_cb", "", 0.88, 0.12, 0.96, 1.00);
  TPad* pad_info = new TPad("pad_info", "", 0.00, 0.00, 1.00, 0.12);

  pad_plot->SetLeftMargin(0.13);
  pad_plot->SetRightMargin(0.015);
  pad_plot->SetTopMargin(0.11);
  pad_plot->SetBottomMargin(0.13);
  pad_plot->SetGridx();
  pad_plot->SetGridy();
  pad_plot->SetFrameLineWidth(2);
  if (cli.log_scale)
    pad_plot->SetLogz();

  pad_cb->SetLeftMargin(0.25);
  pad_cb->SetRightMargin(0.30);
  pad_cb->SetTopMargin(0.11);
  pad_cb->SetBottomMargin(0.13);

  pad_info->SetLeftMargin(0.02);
  pad_info->SetRightMargin(0.02);
  pad_info->SetTopMargin(0.05);
  pad_info->SetBottomMargin(0.05);

  canvas->cd();
  pad_plot->Draw();
  pad_cb->Draw();
  pad_info->Draw();

  // 9. Plot principale
  pad_plot->cd();
  h_Q.Draw("COL");
  h_inv.Draw("BOX same");

  // Marker sul minimo
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
    lbl.DrawLatex(it_best->x1, it_best->x2,
                  (" #bf{min} (" + fmt(it_best->x1, 1) + ", " + fmt(it_best->x2, 1) + ")").c_str());
  }

  // Etichetta unfolding (angolo in alto a sinistra del pad)
  {
    TLatex unf_lbl;
    unf_lbl.SetNDC(true);
    unf_lbl.SetTextFont(42);
    unf_lbl.SetTextSize(0.028);
    unf_lbl.SetTextColor(kGray + 2);
    unf_lbl.SetTextAlign(12);
    if (qcfg.apply_temporal_unfolding) {
      unf_lbl.DrawLatex(0.145, 0.915,
                        ("#delta z_{unfold} = " + fmt(dz_effective, 6) + " mm/passo").c_str());
    } else {
      unf_lbl.DrawLatex(0.145, 0.915, "temporal unfolding: OFF");
    }
  }

  // Titolo
  {
    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextSize(0.046);
    title.SetTextAlign(22);
    title.SetTextColor(kBlack);
    const std::string z_lbl = cli.normalize ? "Q / n_{tracce}" : "Q(x_{1}, x_{2})";
    title.DrawLatex(0.535, 0.953, ("Mappa della funzione di qualit#grave{a}  " + z_lbl).c_str());
  }

  pad_plot->RedrawAxis();

  // 10. Color bar manuale
  pad_cb->cd();
  pad_cb->Range(0.0, 0.0, 1.0, 1.0);

  const int NB = 255;
  double cb_x0 = pad_cb->GetLeftMargin();
  double cb_x1 = 1.0 - pad_cb->GetRightMargin();
  double cb_y0 = pad_cb->GetBottomMargin();
  double cb_y1 = 1.0 - pad_cb->GetTopMargin();

  for (int i = 0; i < NB; ++i) {
    double f0  = static_cast<double>(i) / NB;
    double f1  = static_cast<double>(i + 1) / NB;
    double yb0 = cb_y0 + f0 * (cb_y1 - cb_y0);
    double yb1 = cb_y0 + f1 * (cb_y1 - cb_y0);
    TBox* box  = new TBox(cb_x0, yb0, cb_x1, yb1);
    box->SetFillColor(gStyle->GetColorPalette(i));
    box->SetLineWidth(0);
    box->Draw();
  }

  // Rettangolino grigio per N/A
  {
    double inv_y0 = std::max(0.0, cb_y0 - 0.09);
    double inv_y1 = cb_y0 - 0.01;
    if (inv_y1 > inv_y0) {
      TBox* inv_box = new TBox(cb_x0, inv_y0, cb_x1, inv_y1);
      inv_box->SetFillColor(invalid_color);
      inv_box->SetLineColor(invalid_color);
      inv_box->Draw();
      TLatex inv_lbl;
      inv_lbl.SetNDC();
      inv_lbl.SetTextFont(42);
      inv_lbl.SetTextSize(0.12);
      inv_lbl.SetTextAlign(12);
      inv_lbl.SetTextColor(kWhite);
      inv_lbl.DrawLatex(cb_x0 + 0.04, (inv_y0 + inv_y1) / 2.0, "N/A");
    }
  }

  TGaxis* cb_axis = new TGaxis(cb_x1, cb_y0, cb_x1, cb_y1, h_Q.GetMinimum(), h_Q.GetMaximum(), 505,
                               cli.log_scale ? "+LG" : "+L");
  cb_axis->SetLabelFont(42);
  cb_axis->SetLabelSize(0.18);
  cb_axis->SetTickSize(0.35);
  cb_axis->SetLabelOffset(0.03);
  const std::string z_lbl_cb = cli.normalize ? "Q / n_{tracce}" : "Q(x_{1}, x_{2})";
  cb_axis->SetTitle(z_lbl_cb.c_str());
  cb_axis->SetTitleFont(42);
  cb_axis->SetTitleSize(0.20);
  cb_axis->SetTitleOffset(0.55);
  cb_axis->Draw();

  // 11. Info panel
  pad_info->cd();
  pad_info->SetFillColor(TColor::GetColor(245, 245, 248));

  TLine* sep = new TLine(0.0, 0.97, 1.0, 0.97);
  sep->SetNDC(true);
  sep->SetLineColor(kGray + 1);
  sep->SetLineWidth(1);
  sep->Draw();

  TLatex info;
  info.SetNDC(true);
  info.SetTextFont(42);

  const double col1 = 0.03, col2 = 0.52;
  const double hdr = 0.82, row1 = 0.52, row2 = 0.18;

  info.SetTextSize(0.13);
  info.SetTextColor(kGray + 2);
  info.SetTextAlign(12);
  info.DrawLatex(col1, hdr, "CAMPIONAMENTO y_{0}");
  info.DrawLatex(col2, hdr,
                 ("OTTIMO  (valide: " + std::to_string(total_cfgs - n_cfg_invalid) + "/"
                  + std::to_string(total_cfgs) + " config)")
                     .c_str());

  info.SetTextSize(0.20);
  info.SetTextColor(kBlack);

  info.DrawLatex(col1, row1,
                 ("y_{0} #in [" + fmt(qcfg.y0_min, 1) + ", " + fmt(qcfg.y0_max, 1)
                  + "] mm   #Deltay_{0} = " + fmt(qcfg.dy0, 1) + " mm")
                     .c_str());
  info.DrawLatex(col1, row2,
                 ("L = " + fmt(qcfg.trace_L, 1) + " mm   #Deltat = " + fmt(qcfg.trace_dt, 2)
                  + " mm   soglie: " + fmt(cli.point_frac * 100, 0) + "% pts | "
                  + fmt(cli.trace_frac * 100, 0) + "% tr.")
                     .c_str());

  info.DrawLatex(col2, row1,
                 ("#bf{x_{1}^{*}} = " + fmt(it_best->x1, 1)
                  + " mm,   #bf{x_{2}^{*}} = " + fmt(it_best->x2, 1) + " mm")
                     .c_str());
  {
    std::ostringstream qs;
    qs << std::scientific << std::setprecision(3) << it_best->Q;
    info.DrawLatex(col2, row2, ("#bf{Q_{min}} = " + qs.str()).c_str());
  }

  // 12. Salva
  canvas->Update();
  std::filesystem::create_directories(std::filesystem::path(cli.output_path).parent_path());
  canvas->Print(cli.output_path.c_str());
  std::cout << "\nMappa salvata in: " << cli.output_path << "\n";

  delete canvas;
  return 0;
}