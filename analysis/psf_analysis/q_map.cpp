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
 *         oppure mappa di copertura geometrica (--coverage)
 *
 * MODALITÀ Q (default)
 *
 * Calcola Q(x1,x2) = Σ_{y0} chi²(y0, x1, x2) e ne produce la mappa 2D.
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
 * MODALITÀ COVERAGE (--coverage)
 *
 * Produce una mappa della copertura geometrica:
 *
 *   coverage(x1,x2) = (1/N_y0) · Σ_{y0} [ n_valid(y0) / N_pts_per_trace ]
 *
 * dove n_valid(y0) è il numero di punti della traccia che cadono on_detector.
 * La mappa mostra in percentuale quanta luce arriva effettivamente sul
 * fotocatodo per ogni configurazione, mediata su tutte le posizioni sorgente
 * campionate.
 *
 * Non esegue nessun fit ODR: è una misura puramente geometrica.
 * Utile come diagnostica complementare a Q: configurazioni ad alta copertura
 * ma alto Q raccolgono molti fotoni ma li focalizzano male.
 *
 * In modalità coverage:
 *   - --unfold-dz, --no-unfold, --point-frac, --trace-frac vengono ignorati
 *   - --log e --norm si applicano normalmente
 *   - il minimo di Q è sostituito dal massimo di coverage
 *
 * GESTIONE CELLE INVALIDE
 *
 * Celle grigie = configurazioni in cui build_trace ha fallito per tutti i y0.
 *
 * Uso:
 *   q_map [--psf <path>] [--config <path>] [--output <path>]
 *         [--y0-min <val>] [--y0-max <val>] [--dy0 <val>]
 *         [--L <val>] [--dt <val>]
 *         [--unfold-dz <val>]        (solo modalità Q)
 *         [--no-unfold]              (solo modalità Q)
 *         [--point-frac <val>]       (soglia punti on_detector, default 0.75)
 *         [--trace-frac <val>]       (soglia tracce valide, default 0.75)
 *         [--coverage]               ← nuova modalità: mappa di copertura
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
  std::string output_path = ""; // vuoto = default automatico per modalità
  std::string tsv_path    = "";

  double scint_x = 60.0;
  double scint_y = 20.0;
  double scint_z = 20.0;
  int n_tracks   = 100;

  double dt = 0.1;

  double min_hits   = 10.0;
  double trace_frac = 0.75;

  // Temporal unfolding (solo modalità Q)
  double unfold_dz = 0.0;
  bool no_unfold   = false;

  // Modalità
  bool coverage_mode  = false; // se true: mappa copertura invece di Q
  bool dist_to_target = false;
  double dist_target  = -1.0;

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
    else if (key == "--scint-x")
      cfg.scint_x = std::stod(next());
    else if (key == "--scint-y")
      cfg.scint_y = std::stod(next());
    else if (key == "--scint-z")
      cfg.scint_z = std::stod(next());
    else if (key == "--n-tracks")
      cfg.n_tracks = std::stoi(next());
    else if (key == "--dt")
      cfg.dt = std::stod(next());
    else if (key == "--min-hits")
      cfg.min_hits = std::stod(next());
    else if (key == "--trace-frac")
      cfg.trace_frac = std::stod(next());
    else if (key == "--unfold-dz")
      cfg.unfold_dz = std::stod(next());
    else if (key == "--no-unfold")
      cfg.no_unfold = true;
    else if (key == "--coverage")
      cfg.coverage_mode = true;
    else if (key == "--dist-to-target") {
      cfg.dist_to_target = true;
      if (i + 1 < argc) {
        std::string maybe_val = argv[i + 1];
        if (!maybe_val.empty() && maybe_val.rfind("--", 0) != 0) {
          cfg.dist_target = std::stod(next());
        }
      }
    } else if (key == "--target") {
      cfg.dist_target    = std::stod(next());
      cfg.dist_to_target = true;
    } else if (key == "--log")
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
  gStyle->SetNumberContours(255);
}

static std::string fmt(double v, int n = 1) {
  std::ostringstream o;
  o << std::fixed << std::setprecision(n) << v;
  return o.str();
}

// Disegno canvas comune
// Crea canvas + 3 pad (plot | colorbar | info) e li ritorna via puntatori.
// I pad vengono disegnati sul canvas ma NON inizializzati ulteriormente.
struct PadLayout {
  TCanvas* canvas;
  TPad* pad_plot;
  TPad* pad_cb;
  TPad* pad_info;
};

static PadLayout make_canvas(bool log_z) {
  PadLayout pl;
  pl.canvas = new TCanvas("canvas", "q_map", 1200, 1000);
  pl.canvas->SetLeftMargin(0.0);
  pl.canvas->SetRightMargin(0.0);
  pl.canvas->SetTopMargin(0.0);
  pl.canvas->SetBottomMargin(0.0);

  pl.pad_plot = new TPad("pad_plot", "", 0.00, 0.12, 0.88, 1.00);
  pl.pad_cb   = new TPad("pad_cb", "", 0.88, 0.12, 0.96, 1.00);
  pl.pad_info = new TPad("pad_info", "", 0.00, 0.00, 1.00, 0.12);

  pl.pad_plot->SetLeftMargin(0.10);
  pl.pad_plot->SetRightMargin(0.015);
  pl.pad_plot->SetTopMargin(0.11);
  pl.pad_plot->SetBottomMargin(0.13);
  pl.pad_plot->SetGridx();
  pl.pad_plot->SetGridy();
  pl.pad_plot->SetFrameLineWidth(2);
  if (log_z)
    pl.pad_plot->SetLogz();

  pl.pad_cb->SetLeftMargin(0.25);
  pl.pad_cb->SetRightMargin(0.30);
  pl.pad_cb->SetTopMargin(0.11);
  pl.pad_cb->SetBottomMargin(0.13);

  pl.pad_info->SetLeftMargin(0.02);
  pl.pad_info->SetRightMargin(0.02);
  pl.pad_info->SetTopMargin(0.05);
  pl.pad_info->SetBottomMargin(0.05);

  pl.canvas->cd();
  pl.pad_plot->Draw();
  pl.pad_cb->Draw();
  pl.pad_info->Draw();

  return pl;
}

// Disegna la colorbar manuale nel pad_cb
static void draw_colorbar(TPad* pad_cb, double vmin, double vmax, bool log_scale,
                          const std::string& title, Int_t invalid_color,
                          bool show_invalid_box = true) {
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

  if (show_invalid_box) {
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

  TGaxis* cb_axis =
      new TGaxis(cb_x1, cb_y0, cb_x1, cb_y1, vmin, vmax, 505, log_scale ? "+LG" : "+L");
  cb_axis->SetLabelFont(42);
  cb_axis->SetLabelSize(0.18);
  cb_axis->SetTickSize(0.35);
  cb_axis->SetLabelOffset(0.03);
  cb_axis->SetTitle(title.c_str());
  cb_axis->SetTitleFont(42);
  cb_axis->SetTitleSize(0.20);
  cb_axis->CenterTitle(kTRUE);
  cb_axis->SetTitleOffset(1.8);
  cb_axis->Draw();
}

// Disegno info panel comune
static void draw_info_panel_q(TPad* pad_info, const riptide::QConfig& qcfg, int total_cfgs,
                              int n_invalid, double best_x1, double best_x2, double best_metric,
                              double best_Q_raw, double best_rho, double best_target,
                              bool dist_to_target) {
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
  info.DrawLatex(col1, hdr, "TRACCE CASUALI 3D NELLO SCINTILLATORE");
  info.DrawLatex(col2, hdr,
                 ("OTTIMO  (valide: " + std::to_string(total_cfgs - n_invalid) + "/"
                  + std::to_string(total_cfgs) + " config)")
                     .c_str());

  info.SetTextSize(0.20);
  info.SetTextColor(kBlack);

  info.DrawLatex(col1, row1,
                 ("Scint: " + fmt(qcfg.scint_x, 0) + "x" + fmt(qcfg.scint_y, 0) + "x"
                  + fmt(qcfg.scint_z, 0) + " mm^3   #tracce = " + std::to_string(qcfg.n_tracks))
                     .c_str());
  info.DrawLatex(col1, row2,
                 ("#Deltat = " + fmt(qcfg.trace_dt, 2)
                  + " mm   min_hits = " + fmt(qcfg.min_hits_per_point, 0)
                  + "   soglia tr: " + fmt(qcfg.trace_valid_fraction * 100, 0) + "%")
                     .c_str());

  info.DrawLatex(
      col2, row1,
      ("#bf{x_{1}^{*}} = " + fmt(best_x1, 1) + " mm,   #bf{x_{2}^{*}} = " + fmt(best_x2, 1) + " mm")
          .c_str());
  {
    std::ostringstream ss_metric;
    ss_metric << std::scientific << std::setprecision(3) << best_metric;
    std::ostringstream ss_qraw;
    ss_qraw << std::scientific << std::setprecision(3) << best_Q_raw;
    std::ostringstream ss_tgt;
    ss_tgt << std::scientific << std::setprecision(3) << best_target;
    std::ostringstream ss_rho;
    ss_rho << std::fixed << std::setprecision(3) << best_rho;

    std::string line;
    if (dist_to_target) {
      line = "#bf{|Q-target|_{min}} = " + ss_metric.str() + "   Q = " + ss_qraw.str()
           + "   Q_{target} = " + ss_tgt.str() + "   #hat{#rho} = " + ss_rho.str();
    } else {
      line = "#bf{Q_{min}} = " + ss_qraw.str() + "   Q_{target} = " + ss_tgt.str()
           + "   #hat{#rho} = " + ss_rho.str();
    }
    info.DrawLatex(col2, row2, line.c_str());
  }
}

static void draw_info_panel_coverage(TPad* pad_info, const riptide::QConfig& qcfg, int total_cfgs,
                                     int n_invalid, double best_x1, double best_x2,
                                     double best_cov) {
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
  info.DrawLatex(col1, hdr, "TRACCE CASUALI 3D  [MAPPA DI COPERTURA]");
  info.DrawLatex(col2, hdr,
                 ("MASSIMO COPERTURA  (valide: " + std::to_string(total_cfgs - n_invalid) + "/"
                  + std::to_string(total_cfgs) + " config)")
                     .c_str());

  info.SetTextSize(0.20);
  info.SetTextColor(kBlack);

  info.DrawLatex(col1, row1,
                 ("Scint: " + fmt(qcfg.scint_x, 0) + "x" + fmt(qcfg.scint_y, 0) + "x"
                  + fmt(qcfg.scint_z, 0) + " mm^3   #tracce = " + std::to_string(qcfg.n_tracks))
                     .c_str());
  info.DrawLatex(col1, row2, ("#Deltat = " + fmt(qcfg.trace_dt, 2) + " mm").c_str());

  info.DrawLatex(
      col2, row1,
      ("#bf{x_{1}^{*}} = " + fmt(best_x1, 1) + " mm,   #bf{x_{2}^{*}} = " + fmt(best_x2, 1) + " mm")
          .c_str());
  info.DrawLatex(col2, row2, ("#bf{copertura_{max}} = " + fmt(best_cov * 100.0, 1) + " %").c_str());
}

// main
int main(int argc, char** argv) {
  CliConfig cli = parse_args(argc, argv);

  // Default output path dipende dalla modalità
  if (cli.output_path.empty()) {
    cli.output_path = cli.coverage_mode ? "output/psf_analysis/coverage_map.png"
                                        : "output/psf_analysis/q_map.png";
  }

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
  const double xmin = jcfg["x_min"];
  const double xmax = jcfg["x_max"];

  // 2. QConfig (condivisa tra le due modalità)
  riptide::QConfig qcfg;
  qcfg.scint_x              = cli.scint_x;
  qcfg.scint_y              = cli.scint_y;
  qcfg.scint_z              = cli.scint_z;
  qcfg.n_tracks             = cli.n_tracks;
  qcfg.trace_dt             = cli.dt;
  qcfg.fit_max_iter         = 20;
  qcfg.fit_tol              = 1e-8;
  qcfg.min_hits_per_point   = cli.min_hits;
  qcfg.trace_valid_fraction = cli.trace_frac;

  // Diagnostica console
  if (cli.coverage_mode) {
    std::cout << "Modalità: MAPPA DI COPERTURA (--coverage)\n";
    std::cout << "Scintillatore: " << qcfg.scint_x << "x" << qcfg.scint_y << "x" << qcfg.scint_z
              << " mm^3\n";
    std::cout << "Tracce: " << qcfg.n_tracks << ", dt=" << qcfg.trace_dt << " mm\n";
  } else {
    std::cout << "Modalità: MAPPA Q\n";
    std::cout << "Scintillatore: " << qcfg.scint_x << "x" << qcfg.scint_y << "x" << qcfg.scint_z
              << " mm^3\n";
    std::cout << "Tracce: " << qcfg.n_tracks << ", dt=" << qcfg.trace_dt << " mm\n";
    std::cout << "Soglie: min_hits=" << qcfg.min_hits_per_point
              << ", trace_frac=" << qcfg.trace_valid_fraction << "\n";
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

  if (db.empty()) {
    std::cerr << "Errore: database vuoto\n";
    return 1;
  }

  // Calcola i limiti x1, x2 direttamente dal database
  double x1_lo = std::numeric_limits<double>::max();
  double x1_hi = -std::numeric_limits<double>::max();
  double x2_lo = std::numeric_limits<double>::max();
  double x2_hi = -std::numeric_limits<double>::max();

  for (const auto& [cfg, _] : db) {
    x1_lo = std::min(x1_lo, cfg.x1);
    x1_hi = std::max(x1_hi, cfg.x1);
    x2_lo = std::min(x2_lo, cfg.x2);
    x2_hi = std::max(x2_hi, cfg.x2);
  }

  // Dimensioni griglia istogramma
  int bins_x      = std::max(1, static_cast<int>(std::round((x1_hi - x1_lo) / dx)) + 1);
  int bins_y      = std::max(1, static_cast<int>(std::round((x2_hi - x2_lo) / dx)) + 1);
  double hx       = dx / 2.0;
  double ax_x1_lo = x1_lo - hx, ax_x1_hi = x1_hi + hx;
  double ax_x2_lo = x2_lo - hx, ax_x2_hi = x2_hi + hx;

  Int_t invalid_color = TColor::GetColor(80, 80, 80);
  int total_cfgs      = static_cast<int>(db.size());

  // MODALITÀ COVERAGE
  if (cli.coverage_mode) {
    gStyle->SetPalette(kRainBow);
    gStyle->SetPalette(kViridis);
    gStyle->SetNumberContours(255);

    struct CovEntry {
      double x1, x2;
      double coverage; // ∈ [0,1] oppure -1 se invalida
      int n_evaluated;
      bool config_valid;
    };
    std::vector<CovEntry> entries;
    entries.reserve(db.size());

    int n_cfg_invalid = 0;
    int config_idx    = 0;
    int print_step    = std::max(1, total_cfgs / 20);

    for (const auto& [cfg, _] : db) {
      riptide::CoverageResult res;
      try {
        res = riptide::compute_coverage(cfg, db, qcfg);
      } catch (const std::exception& e) {
        std::cerr << "  [WARN] compute_coverage failed x1=" << cfg.x1 << " x2=" << cfg.x2 << ": "
                  << e.what() << "\n";
        entries.push_back({cfg.x1, cfg.x2, -1.0, 0, false});
        ++n_cfg_invalid;
        ++config_idx;
        continue;
      }

      if (!res.config_valid)
        ++n_cfg_invalid;
      double cov_pct = res.config_valid ? res.coverage : -1.0;
      entries.push_back({cfg.x1, cfg.x2, cov_pct, res.n_y0_evaluated, res.config_valid});

      if ((config_idx + 1) % print_step == 0 || config_idx + 1 == total_cfgs) {
        std::cout << "  [" << std::setw(3) << config_idx + 1 << "/" << total_cfgs << "]"
                  << "  x1=" << fmt(cfg.x1) << "  x2=" << fmt(cfg.x2);
        if (res.config_valid)
          std::cout << "  cov=" << fmt(res.coverage * 100.0, 1) << "%";
        else
          std::cout << "  [INVALIDA]";
        std::cout << "\n";
      }
      ++config_idx;
    }

    std::cout << "\nConfigurazioni valide:   " << total_cfgs - n_cfg_invalid << "\n";
    std::cout << "Configurazioni invalide: " << n_cfg_invalid << "\n";

    // Trova massimo copertura
    auto it_best =
        std::max_element(entries.begin(), entries.end(), [](const CovEntry& a, const CovEntry& b) {
          if (!a.config_valid)
            return true;
          if (!b.config_valid)
            return false;
          return a.coverage < b.coverage;
        });

    if (!it_best->config_valid) {
      std::cerr << "ERRORE: nessuna configurazione valida trovata!\n";
      return 1;
    }
    std::cout << "\n★  Configurazione con copertura massima:\n"
              << "   x1       = " << fmt(it_best->x1, 2) << " mm\n"
              << "   x2       = " << fmt(it_best->x2, 2) << " mm\n"
              << "   copertura = " << fmt(it_best->coverage * 100.0, 2) << " %\n";

    // Esportazione TSV
    if (!cli.tsv_path.empty()) {
      std::filesystem::create_directories(std::filesystem::path(cli.tsv_path).parent_path());
      std::ofstream tsv(cli.tsv_path);
      tsv << "x1\tx2\tcoverage_pct\tn_y0_evaluated\tconfig_valid\n";
      for (const auto& e : entries)
        tsv << e.x1 << "\t" << e.x2 << "\t" << (e.config_valid ? fmt(e.coverage * 100.0, 4) : "NaN")
            << "\t" << e.n_evaluated << "\t" << (e.config_valid ? "1" : "0") << "\n";
      std::cout << "TSV salvato in: " << cli.tsv_path << "\n";
    }

    // Plot coverage
    apply_style();
    gStyle->SetPalette(kViridis);
    gStyle->SetNumberContours(255);

    TH2D h_cov("h_cov", "", bins_x, ax_x1_lo, ax_x1_hi, bins_y, ax_x2_lo, ax_x2_hi);
    TH2D h_inv("h_inv", "", bins_x, ax_x1_lo, ax_x1_hi, bins_y, ax_x2_lo, ax_x2_hi);
    h_cov.Reset();
    h_inv.Reset();

    for (const auto& e : entries) {
      int bx = h_cov.GetXaxis()->FindFixBin(e.x1);
      int by = h_cov.GetYaxis()->FindFixBin(e.x2);
      if (e.config_valid)
        h_cov.SetBinContent(bx, by, e.coverage * 100.0); // in %
      else
        h_inv.SetBinContent(bx, by, 1.0);
    }

    h_cov.GetXaxis()->SetTitle("x_{1}  [mm]  (lente 75 mm)");
    h_cov.GetYaxis()->SetTitle("x_{2}  [mm]  (lente 60 mm)");
    h_cov.GetXaxis()->CenterTitle(true);
    h_cov.GetYaxis()->CenterTitle(true);
    h_cov.GetXaxis()->SetTitleOffset(1.20);
    h_cov.GetYaxis()->SetTitleOffset(1.55);
    h_cov.GetXaxis()->SetNdivisions(506);
    h_cov.GetYaxis()->SetNdivisions(505);

    // Range colori [0, 100] — percentuale intera
    h_cov.SetMinimum(0.0);
    h_cov.SetMaximum(100.0);

    h_inv.SetFillColor(invalid_color);
    h_inv.SetLineColor(invalid_color);

    auto pl = make_canvas(cli.log_scale);

    pl.pad_plot->cd();
    h_cov.Draw("COL");
    for (int iy = 1; iy <= h_inv.GetNbinsY(); ++iy) {
      for (int ix = 1; ix <= h_inv.GetNbinsX(); ++ix) {
        if (h_inv.GetBinContent(ix, iy) > 0.5) {
          TBox* b = new TBox(h_inv.GetXaxis()->GetBinLowEdge(ix),
                             h_inv.GetYaxis()->GetBinLowEdge(iy),
                             h_inv.GetXaxis()->GetBinUpEdge(ix),
                             h_inv.GetYaxis()->GetBinUpEdge(iy));
          b->SetFillColor(invalid_color);
          b->SetLineColor(invalid_color);
          b->Draw("same");
        }
      }
    }

    // Marker stella sul massimo
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
      lbl.DrawLatex(
          it_best->x1, it_best->x2,
          (" #bf{max} (" + fmt(it_best->x1, 1) + ", " + fmt(it_best->x2, 1) + ")").c_str());
    }

    // Titolo
    {
      TLatex title;
      title.SetNDC();
      title.SetTextFont(42);
      title.SetTextSize(0.046);
      title.SetTextAlign(22);
      title.SetTextColor(kBlack);
      title.DrawLatex(0.535, 0.953, "Copertura geometrica  coverage(x_{1}, x_{2})  [%]");
    }

    pl.pad_plot->RedrawAxis();

    draw_colorbar(pl.pad_cb, 0.0, 100.0, cli.log_scale, "copertura [%]", invalid_color,
                  /*show_invalid_box=*/true);

    draw_info_panel_coverage(pl.pad_info, qcfg, total_cfgs, n_cfg_invalid, it_best->x1, it_best->x2,
                             it_best->coverage);

    pl.canvas->Update();
    std::filesystem::create_directories(std::filesystem::path(cli.output_path).parent_path());
    pl.canvas->Print(cli.output_path.c_str());
    std::cout << "\nMappa salvata in: " << cli.output_path << "\n";
    delete pl.canvas;
    return 0;
  }

  // MODALITÀ Q (default)
  gStyle->SetPalette(kBird);
  gStyle->SetNumberContours(255);

  struct QEntry {
    double x1, x2;
    double metric;
    double Q_raw;
    double Q_target;
    double rho_est;
    int n_traces;
    int n_failed;
    int n_invalid;
    bool config_valid;
  };
  std::vector<QEntry> entries;
  entries.reserve(db.size());

  int config_idx    = 0;
  int print_step    = std::max(1, total_cfgs / 20);
  int n_cfg_invalid = 0;

  for (const auto& [cfg, _] : db) {
    riptide::QResult res;
    try {
      res = riptide::compute_Q(cfg, db, qcfg);
    } catch (const std::exception& e) {
      std::cerr << "  [WARN] compute_Q failed x1=" << cfg.x1 << " x2=" << cfg.x2 << ": " << e.what()
                << "\n";
      entries.push_back({cfg.x1, cfg.x2, 0.0, 0.0, 1.0, 0.0, 0, 0, 0, false});
      ++n_cfg_invalid;
      ++config_idx;
      continue;
    }

    double target = (cli.dist_target > 0.0) ? cli.dist_target : res.Q_target;
    double metric = cli.dist_to_target ? std::abs(res.Q - target) : res.Q;

    if (!res.config_valid)
      ++n_cfg_invalid;

    entries.push_back({cfg.x1, cfg.x2, metric, res.Q, target, res.rho_estimate, res.n_traces,
                       res.n_failed, res.n_invalid, res.config_valid});

    if ((config_idx + 1) % print_step == 0 || config_idx + 1 == total_cfgs) {
      std::cout << "  [" << std::setw(3) << config_idx + 1 << "/" << total_cfgs << "]"
                << "  x1=" << fmt(cfg.x1) << "  x2=" << fmt(cfg.x2);
      if (res.config_valid) {
        if (cli.dist_to_target) {
          std::cout << "  |Q-target|=" << std::scientific << std::setprecision(3) << metric
                    << "  Q=" << std::scientific << std::setprecision(3) << res.Q
                    << "  rho=" << std::fixed << std::setprecision(3) << res.rho_estimate;
        } else {
          std::cout << "  Q=" << std::scientific << std::setprecision(3) << metric
                    << "  rho=" << std::fixed << std::setprecision(3) << res.rho_estimate;
        }
      } else {
        std::cout << "  [INVALIDA]";
      }
      std::cout << std::defaultfloat << "\n";
    }
    ++config_idx;
  }

  std::cout << "\nConfigurazioni valide:   " << total_cfgs - n_cfg_invalid << "\n";
  std::cout << "Configurazioni invalide: " << n_cfg_invalid << "\n";

  // Minimo di Q
  auto it_best =
      std::min_element(entries.begin(), entries.end(), [](const QEntry& a, const QEntry& b) {
        if (!a.config_valid)
          return false;
        if (!b.config_valid)
          return true;
        return a.metric < b.metric;
      });

  if (!it_best->config_valid) {
    std::cerr << "ERRORE: nessuna configurazione valida trovata!\n";
    return 1;
  }
  std::cout << "\n★  Configurazione ottimale:\n"
            << "   x1 = " << fmt(it_best->x1, 2) << " mm\n"
            << "   x2 = " << fmt(it_best->x2, 2) << " mm\n";
  if (cli.dist_to_target) {
    std::cout << "   |Q-target| = " << it_best->metric << "\n"
              << "   Q_raw      = " << it_best->Q_raw << "\n"
              << "   Q_target   = " << it_best->Q_target << "\n";
  } else {
    std::cout << "   Q          = " << it_best->Q_raw << "\n"
              << "   Q_target   = " << it_best->Q_target << "\n";
  }
  std::cout << "   rho_hat    = " << it_best->rho_est << "\n";

  // Esportazione TSV
  if (!cli.tsv_path.empty()) {
    std::filesystem::create_directories(std::filesystem::path(cli.tsv_path).parent_path());
    std::ofstream tsv(cli.tsv_path);
    tsv << "x1\tx2\tmetric\tQ_raw\tQ_target\trho_hat\tn_traces\tn_failed\tn_invalid\tconfig_"
           "valid\n";
    for (const auto& e : entries)
      tsv << e.x1 << "\t" << e.x2 << "\t" << (e.config_valid ? std::to_string(e.metric) : "NaN")
          << "\t" << (e.config_valid ? std::to_string(e.Q_raw) : "NaN") << "\t"
          << (e.config_valid ? std::to_string(e.Q_target) : "NaN") << "\t"
          << (e.config_valid ? std::to_string(e.rho_est) : "NaN") << "\t" << e.n_traces << "\t"
          << e.n_failed << "\t" << e.n_invalid << "\t" << (e.config_valid ? "1" : "0") << "\n";
    std::cout << "TSV salvato in: " << cli.tsv_path << "\n";
  }

  // Plot Q
  apply_style();
  gStyle->SetPalette(kBird);
  gStyle->SetNumberContours(255);

  TH2D h_Q("h_Q", "", bins_x, ax_x1_lo, ax_x1_hi, bins_y, ax_x2_lo, ax_x2_hi);
  TH2D h_inv("h_inv", "", bins_x, ax_x1_lo, ax_x1_hi, bins_y, ax_x2_lo, ax_x2_hi);
  h_Q.Reset();
  h_inv.Reset();

  for (const auto& e : entries) {
    int bx = h_Q.GetXaxis()->FindFixBin(e.x1);
    int by = h_Q.GetYaxis()->FindFixBin(e.x2);
    if (e.config_valid)
      h_Q.SetBinContent(bx, by, e.metric);
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
      if (e.config_valid && std::isfinite(e.metric))
        q_vals.push_back(e.metric);
    std::sort(q_vals.begin(), q_vals.end());
    if (!q_vals.empty()) {
      size_t N = q_vals.size();
      h_Q.SetMinimum(q_vals[static_cast<size_t>(0.00 * N)]);
      h_Q.SetMaximum(q_vals[static_cast<size_t>(0.95 * N)]);
    }
  }

  h_inv.SetFillColor(invalid_color);
  h_inv.SetLineColor(invalid_color);

  auto pl = make_canvas(cli.log_scale);

  pl.pad_plot->cd();
  h_Q.Draw("COL");
  for (int iy = 1; iy <= h_inv.GetNbinsY(); ++iy) {
    for (int ix = 1; ix <= h_inv.GetNbinsX(); ++ix) {
      if (h_inv.GetBinContent(ix, iy) > 0.5) {
        TBox* b = new TBox(h_inv.GetXaxis()->GetBinLowEdge(ix),
                           h_inv.GetYaxis()->GetBinLowEdge(iy),
                           h_inv.GetXaxis()->GetBinUpEdge(ix),
                           h_inv.GetYaxis()->GetBinUpEdge(iy));
        b->SetFillColor(invalid_color);
        b->SetLineColor(invalid_color);
        b->Draw("same");
      }
    }
  }

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
    lbl.DrawLatex(it_best->x1, it_best->x2,
                  (" #bf{min} (" + fmt(it_best->x1, 1) + ", " + fmt(it_best->x2, 1) + ")").c_str());
  }

  pl.pad_plot->RedrawAxis();

  std::string z_lbl_cb;
  if (cli.dist_to_target)
    z_lbl_cb = "|#LT#chi^{2}_{red}#GT - Q_{target}|";
  else
    z_lbl_cb = "#LT#chi^{2}_{red}#GT";
  draw_colorbar(pl.pad_cb, h_Q.GetMinimum(), h_Q.GetMaximum(), cli.log_scale, z_lbl_cb,
                invalid_color, /*show_invalid_box=*/true);

  draw_info_panel_q(pl.pad_info, qcfg, total_cfgs, n_cfg_invalid, it_best->x1, it_best->x2,
                    it_best->metric, it_best->Q_raw, it_best->rho_est, it_best->Q_target,
                    cli.dist_to_target);

  pl.canvas->Update();
  std::filesystem::create_directories(std::filesystem::path(cli.output_path).parent_path());
  pl.canvas->Print(cli.output_path.c_str());
  std::cout << "\nMappa salvata in: " << cli.output_path << "\n";
  delete pl.canvas;
  return 0;
}
