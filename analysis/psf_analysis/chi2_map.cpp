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
 * chi2_map — Mappa 2D del chi-quadro ottenuto dal fit di un piano
 *            alle posizioni medie sul detector (mu_y, mu_z) in funzione
 *            della posizione sorgente (x_src, y_src).
 */

#include "psf_interpolator.hpp"

#include <nlohmann/json.hpp>

#include <TAxis.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TF2.h>
#include <TGaxis.h>
#include <TGraph2DErrors.h>
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
  std::string output_path = "output/psf_analysis/chi2_map.png";
  std::string tsv_path    = "";

  double min_hits  = 10.0;
  bool log_scale   = false;
  bool use_reduced = true;
  bool dist_to_n   = false;
  double dist_n    = 1.0;

  // Parametri percentili per la scala colori
  double perc_low  = 0.0;
  double perc_high = 95.0;
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
    else if (key == "--min-hits")
      cfg.min_hits = std::stod(next());
    else if (key == "--log")
      cfg.log_scale = true;
    else if (key == "--dist-to-n") {
      cfg.dist_to_n = true;
      cfg.dist_n    = std::stod(next());
    } else if (key == "--dist-to-one") {
      cfg.dist_to_n = true;
      cfg.dist_n    = 1.0;
    } else if (key == "--no-reduced")
      cfg.use_reduced = false;
    else if (key == "--p-low")
      cfg.perc_low = std::stod(next());
    else if (key == "--p-high")
      cfg.perc_high = std::stod(next());
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

struct PadLayout {
  TCanvas* canvas;
  TPad* pad_plot;
  TPad* pad_cb;
  TPad* pad_info;
};

static PadLayout make_canvas(bool log_z) {
  PadLayout pl;
  pl.canvas = new TCanvas("canvas", "chi2_map", 1200, 1000);
  pl.canvas->SetLeftMargin(0.0);
  pl.canvas->SetRightMargin(0.0);
  pl.canvas->SetTopMargin(0.0);
  pl.canvas->SetBottomMargin(0.0);

  pl.pad_plot = new TPad("pad_plot", "", 0.00, 0.12, 0.88, 1.00);
  pl.pad_cb   = new TPad("pad_cb", "", 0.88, 0.12, 0.96, 1.00);
  pl.pad_info = new TPad("pad_info", "", 0.00, 0.00, 1.00, 0.12);

  pl.pad_plot->SetLeftMargin(0.13);
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

static void draw_colorbar(TPad* pad_cb, double vmin, double vmax, bool log_scale,
                          const std::string& title, Int_t invalid_color) {
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

  double inv_y0 = std::max(0.0, cb_y0 - 0.09);
  double inv_y1 = cb_y0 - 0.01;
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

  TGaxis* cb_axis =
      new TGaxis(cb_x1, cb_y0, cb_x1, cb_y1, vmin, vmax, 505, log_scale ? "+LG" : "+L");
  cb_axis->SetLabelFont(42);
  cb_axis->SetLabelSize(0.18);
  cb_axis->SetTickSize(0.35);
  cb_axis->SetLabelOffset(0.03);
  cb_axis->SetTitle(title.c_str());
  cb_axis->SetTitleFont(42);
  cb_axis->SetTitleSize(0.20);
  cb_axis->SetTitleOffset(0.55);
  cb_axis->Draw();
}

static void draw_info_panel(TPad* pad_info, int total_cfgs, int n_invalid, double best_x1,
                            double best_x2, double best_metric, double best_chi2_red,
                            double min_hits, bool use_reduced, bool dist_to_n, double dist_n) {
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
  info.DrawLatex(col1, hdr, "FIT PIANO: mu_{y,z} = a + b*y_{0} + c*x_{0}");
  info.DrawLatex(col2, hdr,
                 ("OTTIMO  (valide: " + std::to_string(total_cfgs - n_invalid) + "/"
                  + std::to_string(total_cfgs) + " config)")
                     .c_str());

  info.SetTextSize(0.20);
  info.SetTextColor(kBlack);
  info.DrawLatex(col1, row1, ("Min hits PSF: " + fmt(min_hits, 0)).c_str());
  if (dist_to_n) {
    info.DrawLatex(col1, row2,
                   ("Metrica: |#chi^{2}/ndof - " + fmt(dist_n, 3) + "|  (Y+Z)").c_str());
  } else {
    info.DrawLatex(col1, row2,
                   (use_reduced ? "Metrica: #chi^{2}/ndof (Y+Z)" : "Metrica: #chi^{2} (Y+Z)"));
  }
  info.DrawLatex(
      col2, row1,
      ("#bf{x_{1}^{*}} = " + fmt(best_x1, 1) + " mm,   #bf{x_{2}^{*}} = " + fmt(best_x2, 1) + " mm")
          .c_str());

  std::ostringstream ms;
  ms << std::fixed << std::setprecision(4) << best_metric;
  if (dist_to_n) {
    std::ostringstream cr;
    cr << std::fixed << std::setprecision(3) << best_chi2_red;
    info.DrawLatex(col2, row2,
                   ("#bf{min |#chi^{2}/ndof - " + fmt(dist_n, 3) + "|} = " + ms.str()
                    + "   (#chi^{2}/ndof = " + cr.str() + ")")
                       .c_str());
  } else {
    info.DrawLatex(col2, row2, ("#bf{min} = " + ms.str()).c_str());
  }
}

int main(int argc, char** argv) {
  CliConfig cli = parse_args(argc, argv);

  using json = nlohmann::json;
  std::ifstream f_cfg(cli.config_path);
  if (!f_cfg.is_open()) {
    std::cerr << "Errore: impossibile aprire " << cli.config_path << "\n";
    return 1;
  }
  json jcfg;
  f_cfg >> jcfg;
  const double dx = jcfg["dx"];

  riptide::PSFDatabase db;
  try {
    db = riptide::load_psf_database(cli.psf_path);
  } catch (const std::exception& e) {
    std::cerr << "Errore: " << e.what() << "\n";
    return 1;
  }
  if (db.empty())
    return 1;

  double x1_lo = 1e18, x1_hi = -1e18, x2_lo = 1e18, x2_hi = -1e18;
  for (const auto& [cfg, _] : db) {
    x1_lo = std::min(x1_lo, cfg.x1);
    x1_hi = std::max(x1_hi, cfg.x1);
    x2_lo = std::min(x2_lo, cfg.x2);
    x2_hi = std::max(x2_hi, cfg.x2);
  }
  int bins_x      = std::max(1, static_cast<int>(std::round((x1_hi - x1_lo) / dx)) + 1);
  int bins_y      = std::max(1, static_cast<int>(std::round((x2_hi - x2_lo) / dx)) + 1);
  double hx       = dx / 2.0;
  double ax_x1_lo = x1_lo - hx, ax_x1_hi = x1_hi + hx;
  double ax_x2_lo = x2_lo - hx, ax_x2_hi = x2_hi + hx;

  TF2* fit_func = new TF2("fit_plane", "[0] + [1]*x + [2]*y", -100, 100, -100, 100);

  struct Chi2Entry {
    double x1, x2;
    double metric;
    double chi2_red;
    bool valid;
  };
  std::vector<Chi2Entry> entries;
  entries.reserve(db.size());

  int config_idx = 0;
  int total_cfgs = db.size();
  int print_step = std::max(1, total_cfgs / 20);

  for (const auto& [cfg, points] : db) {
    std::vector<const riptide::PSFPoint*> valid_points;
    for (const auto& p : points) {
      if (p.on_detector && p.n_hits_count >= cli.min_hits)
        valid_points.push_back(&p);
    }

    if (valid_points.size() < 10) {
      entries.push_back({cfg.x1, cfg.x2, 0.0, 0.0, false});
    } else {
      bool old_reg = TH1::AddDirectoryStatus();
      TH1::AddDirectory(false);
      TGraph2DErrors g_y(valid_points.size()), g_z(valid_points.size());

      for (size_t i = 0; i < valid_points.size(); ++i) {
        const auto* p = valid_points[i];
        // Converti covarianza della distribuzione in errore sulla media.
        double n_count = std::max(1.0, p->n_hits_count);
        double err_y   = std::sqrt(std::max(0.0, p->cov_yy) / n_count);
        double err_z   = std::sqrt(std::max(0.0, p->cov_zz) / n_count);
        g_y.SetPoint(i, p->y_source, p->x_source, p->mu_y);
        g_y.SetPointError(i, 0, 0, err_y);
        g_z.SetPoint(i, p->y_source, p->x_source, p->mu_z);
        g_z.SetPointError(i, 0, 0, err_z);
      }

      fit_func->SetParameters(0, 0, 0);
      g_y.Fit(fit_func, "QW");
      double c2y     = fit_func->GetChisquare();
      double ndf_y   = static_cast<double>(std::max(1, fit_func->GetNDF()));
      double c2y_red = c2y / ndf_y;

      fit_func->SetParameters(0, 0, 0);
      g_z.Fit(fit_func, "QW");
      double c2z     = fit_func->GetChisquare();
      double ndf_z   = static_cast<double>(std::max(1, fit_func->GetNDF()));
      double c2z_red = c2z / ndf_z;

      double c2_tot     = c2y + c2z;
      double ndf_tot    = std::max(1.0, ndf_y + ndf_z);
      double c2_tot_red = c2_tot / ndf_tot;

      double metric = 0.0;
      if (cli.dist_to_n) {
        metric = std::abs(c2_tot_red - cli.dist_n);
      } else {
        metric = cli.use_reduced ? c2_tot_red : c2_tot;
      }

      if (!std::isfinite(metric) || !std::isfinite(c2_tot_red)) {
        entries.push_back({cfg.x1, cfg.x2, 0.0, 0.0, false});
      } else {
        entries.push_back({cfg.x1, cfg.x2, metric, c2_tot_red, true});
      }
      TH1::AddDirectory(old_reg);
    }

    if ((config_idx + 1) % print_step == 0 || config_idx + 1 == total_cfgs) {
      std::string metric_label = "chi2";
      if (cli.dist_to_n) {
        metric_label = "|chi2/ndof-n|";
      } else if (cli.use_reduced) {
        metric_label = "chi2/ndof";
      }
      std::cout << "  [" << std::setw(3) << config_idx + 1 << "/" << total_cfgs
                << "] x1=" << fmt(cfg.x1) << " " << metric_label << "="
                << (entries.back().valid ? fmt(entries.back().metric, 4) : "N/A") << "\n";
    }
    config_idx++;
  }

  auto it_best =
      std::min_element(entries.begin(), entries.end(), [](const Chi2Entry& a, const Chi2Entry& b) {
        if (!a.valid)
          return false;
        if (!b.valid)
          return true;
        if (a.metric != b.metric)
          return a.metric < b.metric;
        return a.chi2_red < b.chi2_red;
      });

  // --- CALCOLO LIMITI PERCENTILI ---
  std::vector<double> v;
  for (const auto& e : entries)
    if (e.valid)
      v.push_back(e.metric);
  std::sort(v.begin(), v.end());

  double z_min = 0.0, z_max = 1.0;
  if (!v.empty()) {
    size_t i_lo = static_cast<size_t>((cli.perc_low / 100.0) * (v.size() - 1));
    size_t i_hi = static_cast<size_t>((cli.perc_high / 100.0) * (v.size() - 1));
    z_min       = v[std::clamp(i_lo, (size_t)0, v.size() - 1)];
    z_max       = v[std::clamp(i_hi, (size_t)0, v.size() - 1)];
    if (z_min == z_max)
      z_max += 1.0;
  }
  if (cli.log_scale && z_min <= 0.0) {
    double min_pos = std::numeric_limits<double>::infinity();
    for (double vv : v) {
      if (vv > 0.0)
        min_pos = std::min(min_pos, vv);
    }
    if (std::isfinite(min_pos))
      z_min = min_pos;
    else
      z_min = 1e-6;
  }

  apply_style();
  gStyle->SetPalette(kBird);

  TH2D h_chi2("h_chi2", "", bins_x, ax_x1_lo, ax_x1_hi, bins_y, ax_x2_lo, ax_x2_hi);
  TH2D h_inv("h_inv", "", bins_x, ax_x1_lo, ax_x1_hi, bins_y, ax_x2_lo, ax_x2_hi);

  for (const auto& e : entries) {
    int bx = h_chi2.GetXaxis()->FindFixBin(e.x1);
    int by = h_chi2.GetYaxis()->FindFixBin(e.x2);
    if (e.valid)
      h_chi2.SetBinContent(bx, by, e.metric);
    else
      h_inv.SetBinContent(bx, by, 1.0);
  }

  h_chi2.SetMinimum(z_min);
  h_chi2.SetMaximum(z_max);
  h_chi2.GetXaxis()->SetTitle("x_{1} [mm]");
  h_chi2.GetYaxis()->SetTitle("x_{2} [mm]");

  auto pl = make_canvas(cli.log_scale);
  pl.pad_plot->cd();
  h_chi2.Draw("COL");
  h_inv.SetFillColor(TColor::GetColor(80, 80, 80));
  h_inv.Draw("BOX same");

  TMarker* mk = new TMarker(it_best->x1, it_best->x2, 29);
  mk->SetMarkerColor(kRed);
  mk->SetMarkerSize(2.8);
  mk->Draw("same");

  TLatex title;
  title.SetNDC();
  title.SetTextFont(42);
  title.SetTextSize(0.046);
  title.SetTextAlign(22);
  if (cli.dist_to_n) {
    title.DrawLatex(
        0.535, 0.953,
        ("Distanza da #chi^{2}/ndof = " + fmt(cli.dist_n, 3) + "  (fit piano #mu_{y,z})").c_str());
  } else {
    title.DrawLatex(0.535, 0.953, "Linearit#grave{a} della risposta (fit piano #mu_{y,z})");
  }

  std::string cb_title = "#chi^{2}";
  if (cli.dist_to_n) {
    cb_title = "|#chi^{2}/ndof - n|";
  } else if (cli.use_reduced) {
    cb_title = "#chi^{2}/ndof";
  }
  draw_colorbar(pl.pad_cb, z_min, z_max, cli.log_scale, cb_title, TColor::GetColor(80, 80, 80));

  draw_info_panel(
      pl.pad_info, total_cfgs,
      std::count_if(entries.begin(), entries.end(), [](const Chi2Entry& e) { return !e.valid; }),
      it_best->x1, it_best->x2, it_best->metric, it_best->chi2_red, cli.min_hits, cli.use_reduced,
      cli.dist_to_n, cli.dist_n);

  pl.canvas->Update();
  std::filesystem::create_directories(std::filesystem::path(cli.output_path).parent_path());
  pl.canvas->Print(cli.output_path.c_str());

  return 0;
}
