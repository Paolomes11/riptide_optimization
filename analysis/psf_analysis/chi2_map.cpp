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
  std::string output_path = "output/psf_analysis/chi2_map.png";
  std::string tsv_path    = "";

  double min_hits  = 10.0;
  bool log_scale   = false;
  bool use_reduced = true;
  bool dist_to_n   = false;
  double dist_n    = 1.0;

  bool corr_map        = false;
  bool adaptive_target = false;

  bool apply_corr = true;
  double grid_dx  = -1.0;
  double grid_dy  = -1.0;

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
    else if (key == "--corr-map")
      cfg.corr_map = true;
    else if (key == "--adaptive-target")
      cfg.adaptive_target = true;
    else if (key == "--no-corr")
      cfg.apply_corr = false;
    else if (key == "--grid-dx")
      cfg.grid_dx = std::stod(next());
    else if (key == "--grid-dy")
      cfg.grid_dy = std::stod(next());
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

struct PlaneFitResult {
  double a = 0.0, b = 0.0, c = 0.0;
  double chi2      = 0.0;
  int ndof         = 0;
  double chi2_ndof = 0.0;
  bool valid       = false;

  std::vector<double> residuals;
  std::vector<double> residual_sig;
};

static bool solve_plane_wls(const std::vector<double>& x_src, const std::vector<double>& y_src,
                            const std::vector<double>& mu, const std::vector<double>& sigma,
                            PlaneFitResult& res) {
  const int N = static_cast<int>(x_src.size());
  if (N < 4 || static_cast<int>(y_src.size()) != N || static_cast<int>(mu.size()) != N
      || static_cast<int>(sigma.size()) != N) {
    res.valid = false;
    return false;
  }

  double M[3][3] = {};
  double rhs[3]  = {};

  for (int i = 0; i < N; ++i) {
    const double sig2 = std::max(sigma[i] * sigma[i], 1e-30);
    const double w    = 1.0 / sig2;
    const double yi   = y_src[i];
    const double xi   = x_src[i];
    const double mi   = mu[i];

    const double row[3] = {1.0, yi, xi};
    for (int r = 0; r < 3; ++r) {
      rhs[r] += w * row[r] * mi;
      for (int c = 0; c < 3; ++c)
        M[r][c] += w * row[r] * row[c];
    }
  }

  const double det = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
                   - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                   + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);

  if (!std::isfinite(det) || std::abs(det) < 1e-30) {
    res.valid = false;
    return false;
  }

  double inv[3][3];
  inv[0][0] = (M[1][1] * M[2][2] - M[1][2] * M[2][1]) / det;
  inv[0][1] = -(M[0][1] * M[2][2] - M[0][2] * M[2][1]) / det;
  inv[0][2] = (M[0][1] * M[1][2] - M[0][2] * M[1][1]) / det;
  inv[1][0] = -(M[1][0] * M[2][2] - M[1][2] * M[2][0]) / det;
  inv[1][1] = (M[0][0] * M[2][2] - M[0][2] * M[2][0]) / det;
  inv[1][2] = -(M[0][0] * M[1][2] - M[0][2] * M[1][0]) / det;
  inv[2][0] = (M[1][0] * M[2][1] - M[1][1] * M[2][0]) / det;
  inv[2][1] = -(M[0][0] * M[2][1] - M[0][1] * M[2][0]) / det;
  inv[2][2] = (M[0][0] * M[1][1] - M[0][1] * M[1][0]) / det;

  double theta[3] = {};
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c)
      theta[r] += inv[r][c] * rhs[c];

  if (!std::isfinite(theta[0]) || !std::isfinite(theta[1]) || !std::isfinite(theta[2])) {
    res.valid = false;
    return false;
  }

  res.a = theta[0];
  res.b = theta[1];
  res.c = theta[2];

  res.chi2 = 0.0;
  res.residuals.resize(static_cast<size_t>(N));
  res.residual_sig.resize(static_cast<size_t>(N));
  for (int i = 0; i < N; ++i) {
    const double mu_hat = res.a + res.b * y_src[i] + res.c * x_src[i];
    const double r_i    = mu[i] - mu_hat;
    const double sig2   = std::max(sigma[i] * sigma[i], 1e-30);
    res.chi2 += r_i * r_i / sig2;
    res.residuals[static_cast<size_t>(i)]    = r_i;
    res.residual_sig[static_cast<size_t>(i)] = sigma[i];
  }

  res.ndof      = N - 3;
  res.chi2_ndof = (res.ndof > 0) ? res.chi2 / static_cast<double>(res.ndof) : 0.0;
  res.valid     = std::isfinite(res.chi2_ndof);
  return res.valid;
}

static double estimate_min_step(std::vector<double> v) {
  v.erase(std::remove_if(v.begin(), v.end(), [](double x) { return !std::isfinite(x); }), v.end());
  if (v.size() < 2)
    return 0.0;
  std::sort(v.begin(), v.end());
  v.erase(std::unique(v.begin(), v.end()), v.end());
  double step = std::numeric_limits<double>::infinity();
  for (size_t i = 1; i < v.size(); ++i) {
    const double d = v[i] - v[i - 1];
    if (d > 0.0 && d < step)
      step = d;
  }
  return std::isfinite(step) ? step : 0.0;
}

static double estimate_residual_correlation(const std::vector<double>& x_src,
                                            const std::vector<double>& y_src,
                                            const std::vector<double>& residuals, double dx_thresh,
                                            double dy_thresh, int& n_pairs_out) {
  n_pairs_out = 0;
  if (dx_thresh <= 0.0 || dy_thresh <= 0.0)
    return 0.0;

  const int N = static_cast<int>(residuals.size());
  if (N < 2)
    return 0.0;

  double mean_r = 0.0;
  for (int i = 0; i < N; ++i)
    mean_r += residuals[static_cast<size_t>(i)];
  mean_r /= static_cast<double>(N);

  double num = 0.0, den = 0.0;

  for (int i = 0; i < N; ++i) {
    for (int j = i + 1; j < N; ++j) {
      const double dx = std::abs(x_src[static_cast<size_t>(i)] - x_src[static_cast<size_t>(j)]);
      const double dy = std::abs(y_src[static_cast<size_t>(i)] - y_src[static_cast<size_t>(j)]);
      if (dx <= dx_thresh && dy <= dy_thresh) {
        const double ri = residuals[static_cast<size_t>(i)] - mean_r;
        const double rj = residuals[static_cast<size_t>(j)] - mean_r;
        num += ri * rj;
        den += 0.5 * (ri * ri + rj * rj);
        ++n_pairs_out;
      }
    }
  }

  if (n_pairs_out == 0 || den < 1e-30)
    return 0.0;
  return std::clamp(num / den, 0.0, 0.999);
}

static double compute_inflation_factor(double rho, int N, double n_neighbors_avg) {
  if (N <= 3)
    return 1.0;
  const double frac = n_neighbors_avg / static_cast<double>(N);
  return 1.0 + 2.0 * rho * frac * static_cast<double>(N - 3) / static_cast<double>(N - 1);
}

static int count_distinct(std::vector<double> v, double eps = 1e-3) {
  v.erase(std::remove_if(v.begin(), v.end(), [](double x) { return !std::isfinite(x); }), v.end());
  if (v.empty())
    return 0;
  std::sort(v.begin(), v.end());
  int n = 1;
  for (size_t i = 1; i < v.size(); ++i) {
    if (v[i] - v[i - 1] > eps)
      ++n;
  }
  return n;
}

static double compute_chi2_target(double rho, int N_x, int N_y) {
  if (rho <= 0.0 || (N_x <= 2 && N_y <= 2))
    return 1.0;

  rho = std::clamp(rho, 0.0, 0.999999);

  const bool use_x = N_x > 2;
  const bool use_y = N_y > 2;

  if (use_x && !use_y)
    return riptide::expected_chi2_ndof_ar1(N_x, rho);
  if (!use_x && use_y)
    return riptide::expected_chi2_ndof_ar1(N_y, rho);

  const double tx = riptide::expected_chi2_ndof_ar1(N_x, rho);
  const double ty = riptide::expected_chi2_ndof_ar1(N_y, rho);
  const double wx = static_cast<double>(N_x) / static_cast<double>(N_x + N_y);
  const double wy = 1.0 - wx;
  const double t  = wx * tx + wy * ty;
  return std::isfinite(t) ? t : 1.0;
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
                            double best_chi2_red_raw, double best_rho, double best_chi2_target,
                            double min_hits, bool use_reduced, bool dist_to_n, double dist_n,
                            bool corr_mode, bool adaptive_target) {
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
  if (corr_mode) {
    info.DrawLatex(col1, row2, "Metrica: #rho residui vicini (Y+Z)");
  } else if (adaptive_target) {
    info.DrawLatex(col1, row2, "Metrica: |#chi^{2}/ndof - #chi^{2}_{target}(#hat{#rho})|  (Y+Z)");
  } else if (dist_to_n) {
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
  if (corr_mode) {
    std::ostringstream cr;
    cr << std::fixed << std::setprecision(3) << best_chi2_red;
    info.DrawLatex(
        col2, row2,
        ("#bf{max #rho} = " + ms.str() + "   (#chi^{2}/ndof = " + cr.str() + ")").c_str());
  } else if (adaptive_target) {
    std::ostringstream cr;
    std::ostringstream ct;
    std::ostringstream rr;
    cr << std::fixed << std::setprecision(3) << best_chi2_red_raw;
    ct << std::fixed << std::setprecision(3) << best_chi2_target;
    rr << std::fixed << std::setprecision(3) << best_rho;
    info.DrawLatex(col2, row2,
                   ("#bf{min} = " + ms.str() + "   (#chi^{2}/ndof = " + cr.str()
                    + "   target = " + ct.str() + "   #hat{#rho} = " + rr.str() + ")")
                       .c_str());
  } else if (dist_to_n) {
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
  if (cli.corr_map && cli.output_path == "output/psf_analysis/chi2_map.png")
    cli.output_path = "output/psf_analysis/corr_map.png";
  if (cli.adaptive_target && !cli.corr_map && cli.output_path == "output/psf_analysis/chi2_map.png")
    cli.output_path = "output/psf_analysis/chi2_map_adaptive.png";

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

  struct Chi2Entry {
    double x1, x2;
    double metric;
    double chi2_red;
    double chi2_red_raw;
    double rho;
    double chi2_target;
    double infl;
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
      entries.push_back({cfg.x1, cfg.x2, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, false});
    } else {
      std::vector<double> xs, ys, mu_y_vec, mu_z_vec, sig_y_vec, sig_z_vec;
      xs.reserve(valid_points.size());
      ys.reserve(valid_points.size());
      mu_y_vec.reserve(valid_points.size());
      mu_z_vec.reserve(valid_points.size());
      sig_y_vec.reserve(valid_points.size());
      sig_z_vec.reserve(valid_points.size());

      for (const auto* p : valid_points) {
        const double n_count = std::max(1.0, p->n_hits_count);
        const double err_y   = std::sqrt(std::max(0.0, p->cov_yy) / n_count);
        const double err_z   = std::sqrt(std::max(0.0, p->cov_zz) / n_count);

        xs.push_back(p->x_source);
        ys.push_back(p->y_source);
        mu_y_vec.push_back(p->mu_y);
        mu_z_vec.push_back(p->mu_z);
        sig_y_vec.push_back(std::max(err_y, 1e-6));
        sig_z_vec.push_back(std::max(err_z, 1e-6));
      }

      PlaneFitResult fit_y, fit_z;
      const bool ok_y = solve_plane_wls(xs, ys, mu_y_vec, sig_y_vec, fit_y);
      const bool ok_z = solve_plane_wls(xs, ys, mu_z_vec, sig_z_vec, fit_z);

      if (!ok_y || !ok_z) {
        entries.push_back({cfg.x1, cfg.x2, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, false});
      } else {
        const double c2_tot   = fit_y.chi2 + fit_z.chi2;
        const int ndof_tot_i  = fit_y.ndof + fit_z.ndof;
        const double ndof_tot = static_cast<double>(ndof_tot_i);
        const double c2_tot_red_raw =
            (ndof_tot_i > 0) ? c2_tot / ndof_tot : std::numeric_limits<double>::quiet_NaN();

        double c2_tot_red  = c2_tot_red_raw;
        double rho         = 0.0;
        double chi2_target = 1.0;
        double infl        = 1.0;

        double dx_grid = cli.grid_dx;
        double dy_grid = cli.grid_dy;
        if (dx_grid <= 0.0)
          dx_grid = estimate_min_step(xs);
        if (dy_grid <= 0.0)
          dy_grid = estimate_min_step(ys);

        bool corr_ok = false;
        if (dx_grid > 0.0 && dy_grid > 0.0 && fit_y.residuals.size() == xs.size()
            && fit_z.residuals.size() == xs.size()) {
          std::vector<double> rstd_y(xs.size()), rstd_z(xs.size());
          for (size_t i = 0; i < xs.size(); ++i) {
            const double sy = std::max(fit_y.residual_sig[i], 1e-30);
            const double sz = std::max(fit_z.residual_sig[i], 1e-30);
            rstd_y[i]       = fit_y.residuals[i] / sy;
            rstd_z[i]       = fit_z.residuals[i] / sz;
          }

          int n_pairs_y = 0;
          int n_pairs_z = 0;
          const double rho_y =
              estimate_residual_correlation(xs, ys, rstd_y, dx_grid, dy_grid, n_pairs_y);
          const double rho_z =
              estimate_residual_correlation(xs, ys, rstd_z, dx_grid, dy_grid, n_pairs_z);

          rho = 0.5 * (rho_y + rho_z);

          const double n_pairs_avg =
              0.5 * (static_cast<double>(n_pairs_y) + static_cast<double>(n_pairs_z));
          const double neighborsAvg =
              (xs.empty()) ? 0.0 : (2.0 * n_pairs_avg / static_cast<double>(xs.size()));

          infl    = compute_inflation_factor(rho, static_cast<int>(xs.size()), neighborsAvg);
          corr_ok = std::isfinite(rho) && std::isfinite(infl) && infl > 0.0
                 && (n_pairs_y + n_pairs_z) > 0;
        }

        const int N_x = count_distinct(xs);
        const int N_y = count_distinct(ys);
        chi2_target   = compute_chi2_target(rho, N_x, N_y);

        const bool apply_corr_effective = cli.apply_corr && !cli.adaptive_target;
        if (apply_corr_effective && std::isfinite(infl) && infl > 1e-6)
          c2_tot_red = c2_tot_red_raw / infl;

        double metric = 0.0;
        if (cli.corr_map) {
          metric = rho;
        } else if (cli.dist_to_n) {
          metric = std::abs(c2_tot_red - cli.dist_n);
        } else if (cli.adaptive_target) {
          metric = std::abs(c2_tot_red_raw - chi2_target);
        } else {
          metric = cli.use_reduced ? c2_tot_red : c2_tot;
        }

        const bool valid =
            cli.corr_map ? (corr_ok && std::isfinite(metric)) : std::isfinite(metric);
        if (!valid || !std::isfinite(c2_tot_red)) {
          entries.push_back({cfg.x1, cfg.x2, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, false});
        } else {
          entries.push_back(
              {cfg.x1, cfg.x2, metric, c2_tot_red, c2_tot_red_raw, rho, chi2_target, infl, true});
        }
      }
    }

    if ((config_idx + 1) % print_step == 0 || config_idx + 1 == total_cfgs) {
      std::string metric_label = "chi2";
      if (cli.corr_map) {
        metric_label = "rho";
      } else if (cli.dist_to_n) {
        metric_label = "|chi2/ndof-n|";
      } else if (cli.adaptive_target) {
        metric_label = "|chi2/ndof-target|";
      } else if (cli.use_reduced) {
        metric_label = "chi2/ndof";
      }
      std::cout << "  [" << std::setw(3) << config_idx + 1 << "/" << total_cfgs
                << "] x1=" << fmt(cfg.x1) << " " << metric_label << "="
                << (entries.back().valid ? fmt(entries.back().metric, 4) : "N/A") << "\n";
    }
    config_idx++;
  }

  auto it_best = cli.corr_map ? std::max_element(entries.begin(), entries.end(),
                                                 [](const Chi2Entry& a, const Chi2Entry& b) {
                                                   if (!a.valid)
                                                     return true;
                                                   if (!b.valid)
                                                     return false;
                                                   if (a.metric != b.metric)
                                                     return a.metric < b.metric;
                                                   return a.chi2_red < b.chi2_red;
                                                 })
                              : std::min_element(entries.begin(), entries.end(),
                                                 [](const Chi2Entry& a, const Chi2Entry& b) {
                                                   if (!a.valid)
                                                     return false;
                                                   if (!b.valid)
                                                     return true;
                                                   if (a.metric != b.metric)
                                                     return a.metric < b.metric;
                                                   return a.chi2_red < b.chi2_red;
                                                 });

  if (!cli.tsv_path.empty()) {
    std::filesystem::create_directories(std::filesystem::path(cli.tsv_path).parent_path());
    std::ofstream out(cli.tsv_path);
    if (!out.is_open()) {
      std::cerr << "Errore: impossibile aprire " << cli.tsv_path << "\n";
      return 1;
    }

    out << "x1\tx2\tmetric\tchi2_red\tchi2_red_raw\trho\tchi2_target\tinfl\tvalid\n";
    for (const auto& e : entries) {
      out << std::fixed << std::setprecision(6) << e.x1 << "\t" << e.x2 << "\t" << e.metric << "\t"
          << e.chi2_red << "\t" << e.chi2_red_raw << "\t" << e.rho << "\t" << e.chi2_target << "\t"
          << e.infl << "\t" << (e.valid ? 1 : 0) << "\n";
    }
  }

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
  if (cli.corr_map) {
    title.DrawLatex(0.535, 0.953, "Correlazione residui (fit piano #mu_{y,z})");
  } else if (cli.adaptive_target) {
    title.DrawLatex(0.535, 0.953, "Distanza da target adattivo (fit piano #mu_{y,z})");
  } else if (cli.dist_to_n) {
    title.DrawLatex(
        0.535, 0.953,
        ("Distanza da #chi^{2}/ndof = " + fmt(cli.dist_n, 3) + "  (fit piano #mu_{y,z})").c_str());
  } else {
    title.DrawLatex(0.535, 0.953, "Linearit#grave{a} della risposta (fit piano #mu_{y,z})");
  }

  std::string cb_title = "#chi^{2}";
  if (cli.corr_map) {
    cb_title = "#rho";
  } else if (cli.adaptive_target) {
    cb_title = "|#chi^{2}/ndof - #chi^{2}_{target}|";
  } else if (cli.dist_to_n) {
    cb_title = "|#chi^{2}/ndof - n|";
  } else if (cli.use_reduced) {
    cb_title = "#chi^{2}/ndof";
  }
  draw_colorbar(pl.pad_cb, z_min, z_max, cli.log_scale, cb_title, TColor::GetColor(80, 80, 80));

  draw_info_panel(
      pl.pad_info, total_cfgs,
      std::count_if(entries.begin(), entries.end(), [](const Chi2Entry& e) { return !e.valid; }),
      it_best->x1, it_best->x2, it_best->metric, it_best->chi2_red, it_best->chi2_red_raw,
      it_best->rho, it_best->chi2_target, cli.min_hits, cli.use_reduced, cli.dist_to_n, cli.dist_n,
      cli.corr_map, cli.adaptive_target);

  pl.canvas->Update();
  std::filesystem::create_directories(std::filesystem::path(cli.output_path).parent_path());
  pl.canvas->Print(cli.output_path.c_str());

  return 0;
}
