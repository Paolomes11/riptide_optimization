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

#include <TAttFill.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TEllipse.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMarker.h>
#include <TPaveText.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

//  Parsing CLI

struct Config {
  double x1  = -1.0;
  double x2  = -1.0;
  double y0  = -1.0;
  double p1x = 1e9, p1y = 1e9, p1z = 1e9;
  double p2x = 1e9, p2y = 1e9, p2z = 1e9;
  double dt            = 0.1;
  double L             = 10.0;
  double sigma_scale   = 1.0;
  std::string psf_path = "output/psf/psf_data.root";
  std::string output   = "";
  bool fit             = false;
};

Config parse_args(int argc, char** argv) {
  Config cfg;
  for (int i = 1; i < argc; ++i) {
    std::string key = argv[i];
    auto next       = [&]() -> std::string {
      if (i + 1 >= argc)
        return "";
      return argv[++i];
    };
    if (key == "--x1")
      cfg.x1 = std::stod(next());
    else if (key == "--x2")
      cfg.x2 = std::stod(next());
    else if (key == "--y0")
      cfg.y0 = std::stod(next());
    else if (key == "--p1x")
      cfg.p1x = std::stod(next());
    else if (key == "--p1y")
      cfg.p1y = std::stod(next());
    else if (key == "--p1z")
      cfg.p1z = std::stod(next());
    else if (key == "--p2x")
      cfg.p2x = std::stod(next());
    else if (key == "--p2y")
      cfg.p2y = std::stod(next());
    else if (key == "--p2z")
      cfg.p2z = std::stod(next());
    else if (key == "--dt")
      cfg.dt = std::stod(next());
    else if (key == "--L")
      cfg.L = std::stod(next());
    else if (key == "--sigma")
      cfg.sigma_scale = std::stod(next());
    else if (key == "--psf")
      cfg.psf_path = next();
    else if (key == "--output")
      cfg.output = next();
    else if (key == "--fit")
      cfg.fit = true;
  }
  return cfg;
}

//  Calcolo ellisse di covarianza

struct Ellipse {
  double cy, cz, a, b, theta_deg;
};

Ellipse cov_to_ellipse(double cy, double cz, const riptide::Cov2& cov, double n_sigma = 1.0) {
  double trace = cov.yy + cov.zz;
  double det   = cov.yy * cov.zz - cov.yz * cov.yz;
  double disc  = std::sqrt(std::max(0.0, trace * trace / 4.0 - det));
  double l1    = trace / 2.0 + disc;
  double l2    = trace / 2.0 - disc;
  double a     = n_sigma * std::sqrt(std::max(0.0, l1));
  double b     = n_sigma * std::sqrt(std::max(0.0, l2));
  double theta = 0.0;
  if (std::abs(cov.yz) > 1e-12 || std::abs(cov.yy - l1) > 1e-12)
    theta = std::atan2(l1 - cov.yy, cov.yz);
  return {cy, cz, a, b, theta * 180.0 / M_PI};
}

//  Palette sequenziale personalizzata

Int_t create_trace_palette() {
  gStyle->SetPalette(kRainBow);
  gStyle->SetNumberContours(255);
  return 0;
}

// Restituisce il colore interpolato nella palette custom per t in [0,1]
Int_t palette_color(double t) {
  int n   = gStyle->GetNumberOfColors();
  int idx = static_cast<int>(std::clamp(t, 0.0, 1.0) * (n - 1));
  return gStyle->GetColorPalette(idx);
}

//  Stile globale publication-quality

void apply_style() {
  gStyle->Reset();
  // Font: Helvetica (42) per tutto
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleFont(42, ""); // titolo pad
  gStyle->SetStatFont(42);

  // Dimensioni
  gStyle->SetTextSize(0.040);
  gStyle->SetLabelSize(0.038, "XYZ");
  gStyle->SetTitleSize(0.044, "XYZ");
  gStyle->SetTitleSize(0.046, "");

  // Offset assi
  gStyle->SetTitleOffset(1.55, "Y");
  gStyle->SetTitleOffset(1.20, "X");
  gStyle->SetTitleOffset(1.15, "Z");

  // Tick
  gStyle->SetTickLength(0.020, "X");
  gStyle->SetTickLength(0.020, "Y");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  // Griglia leggera
  gStyle->SetGridColor(kGray);
  gStyle->SetGridStyle(3); // tratteggio fine
  gStyle->SetGridWidth(1);

  // Nessuna stat box, nessun titolo frame
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // Canvas / Pad
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);

  // Divisioni asse
  gStyle->SetNdivisions(505, "X");
  gStyle->SetNdivisions(505, "Y");
}

//  Helpers per testo NDC

// Disegna una TLatex NDC con stile fisso
void draw_label(TLatex& lat, double x, double y, const std::string& text, double size = 0.033,
                int align = 12, int color = kBlack) {
  lat.SetNDC();
  lat.SetTextSize(size);
  lat.SetTextAlign(align);
  lat.SetTextColor(color);
  lat.DrawLatex(x, y, text.c_str());
}

// Formatta un double con n cifre decimali
std::string fmt(double v, int n = 1) {
  std::ostringstream o;
  o << std::fixed << std::setprecision(n) << v;
  return o.str();
}

//  main

int main(int argc, char** argv) {
  Config cfg = parse_args(argc, argv);

  bool use_3d = (cfg.p1x < 1e8 && cfg.p1y < 1e8 && cfg.p1z < 1e8 && cfg.p2x < 1e8 && cfg.p2y < 1e8
                 && cfg.p2z < 1e8);

  if (cfg.x1 < 0 || cfg.x2 < 0 || (cfg.y0 < 0 && !use_3d)) {
    std::cerr << "Uso (2D): trace_viewer --x1 <val> --x2 <val> --y0 <val> [--L <val>]\n"
              << "Uso (3D): trace_viewer --x1 <val> --x2 <val> --p1x <v> --p1y <v> --p1z <v> --p2x "
                 "<v> --p2y <v> --p2z <v>\n"
              << "Opzioni:\n"
              << "     [--psf <path>]     default: output/psf/psf_data.root\n"
              << "     [--output <path>]  default: output/psf_analysis/trace_*.png\n"
              << "     [--dt <val>]       step traccia [mm], default: 0.1\n"
              << "     [--sigma <val>]    scala ellissi, default: 1.0\n";
    return 1;
  }

  // Carica PSF database
  riptide::PSFDatabase db;
  try {
    db = riptide::load_psf_database(cfg.psf_path);
  } catch (const std::exception& e) {
    std::cerr << "Errore: " << e.what() << "\n";
    return 1;
  }

  // Costruisce la traccia
  riptide::LensConfig lens{cfg.x1, cfg.x2};
  riptide::LensConfig actual = riptide::find_nearest_config(lens, db);

  riptide::Point3D p1, p2;
  if (use_3d) {
    p1 = {cfg.p1x, cfg.p1y, cfg.p1z};
    p2 = {cfg.p2x, cfg.p2y, cfg.p2z};
  } else {
    p1 = {-cfg.L / 2.0, cfg.y0, 0.0};
    p2 = {cfg.L / 2.0, cfg.y0, 0.0};
  }

  std::vector<riptide::TracePoint> trace;
  try {
    trace = riptide::build_trace_3d(p1, p2, actual, db, cfg.dt);
  } catch (const std::exception& e) {
    std::cerr << "Errore nella costruzione della traccia: " << e.what() << "\n";
    return 1;
  }

  if (trace.empty()) {
    std::cerr << "Errore: traccia vuota\n";
    return 1;
  }

  double total_L = std::hypot(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);

  // Eventuale fit
  riptide::LineFitResult fit_res;
  bool has_fit = false;
  if (cfg.fit) {
    try {
      fit_res = riptide::fit_trace(trace);
      has_fit = true;
    } catch (const std::exception& e) {
      std::cerr << "Fit fallito: " << e.what() << "\n";
    }
  }

  const int N = static_cast<int>(trace.size());

  std::vector<double> gy, gz;
  gy.reserve(N);
  gz.reserve(N);
  for (const auto& pt : trace) {
    gy.push_back(pt.mu_y);
    gz.push_back(pt.mu_z);
  }

  double y_min       = std::numeric_limits<double>::infinity();
  double y_max       = -std::numeric_limits<double>::infinity();
  double z_min       = std::numeric_limits<double>::infinity();
  double z_max       = -std::numeric_limits<double>::infinity();
  double max_sigma_y = 0.0, max_sigma_z = 0.0;

  for (const auto& pt : trace) {
    max_sigma_y = std::max(max_sigma_y, std::sqrt(std::max(0.0, pt.cov.yy)));
    max_sigma_z = std::max(max_sigma_z, std::sqrt(std::max(0.0, pt.cov.zz)));

    auto el = cov_to_ellipse(pt.mu_y, pt.mu_z, pt.cov, cfg.sigma_scale);
    if (el.a > 0.0 && el.b > 0.0) {
      double theta = el.theta_deg * M_PI / 180.0;
      double c     = std::cos(theta);
      double s     = std::sin(theta);
      double dy    = std::sqrt(el.a * el.a * c * c + el.b * el.b * s * s);
      double dz    = std::sqrt(el.a * el.a * s * s + el.b * el.b * c * c);
      y_min        = std::min(y_min, pt.mu_y - dy);
      y_max        = std::max(y_max, pt.mu_y + dy);
      z_min        = std::min(z_min, pt.mu_z - dz);
      z_max        = std::max(z_max, pt.mu_z + dz);
    } else {
      y_min = std::min(y_min, pt.mu_y);
      y_max = std::max(y_max, pt.mu_y);
      z_min = std::min(z_min, pt.mu_z);
      z_max = std::max(z_max, pt.mu_z);
    }
  }

  if (has_fit) {
    double y_start = y_min;
    double y_end   = y_max;
    double z_start = fit_res.a * y_start + fit_res.b;
    double z_end   = fit_res.a * y_end + fit_res.b;
    y_min          = std::min({y_min, y_start, y_end});
    y_max          = std::max({y_max, y_start, y_end});
    z_min          = std::min({z_min, z_start, z_end});
    z_max          = std::max({z_max, z_start, z_end});
  }

  if (!std::isfinite(y_min) || !std::isfinite(y_max) || !std::isfinite(z_min)
      || !std::isfinite(z_max)) {
    std::cerr << "Errore: bounds non finiti\n";
    return 1;
  }

  double y_range  = std::max(0.0, y_max - y_min);
  double z_range  = std::max(0.0, z_max - z_min);
  double margin_y = std::max({0.05 * y_range, 3.0 * max_sigma_y * cfg.sigma_scale, 0.005});
  double margin_z = std::max({0.05 * z_range, 3.0 * max_sigma_z * cfg.sigma_scale, 0.005});

  double py_min = y_min - margin_y;
  double py_max = y_max + margin_y;
  double pz_min = z_min - margin_z;
  double pz_max = z_max + margin_z;

  //  Applica stile e crea palette

  apply_style();

  //  Canvas principale: proporzione 4:3, risoluzione alta per pubblicazione

  TCanvas* canvas = new TCanvas("canvas", "Trace Viewer", 1200, 900);
  canvas->SetLeftMargin(0.0);
  canvas->SetRightMargin(0.0);
  canvas->SetTopMargin(0.0);
  canvas->SetBottomMargin(0.0);

  // Layout a due pannelli: plot (sinistra, ~82%) + color bar verticale (destra, ~6%)
  // + pannello info (fondo, ~12%)
  //
  //  ┌────────────────────────────────┬──┐
  //  │                                │CB│  ← plot + colorbar (88% altezza)
  //  ├────────────────────────────────┴──┤
  //  │           info panel              │  ← 12% altezza
  //  └───────────────────────────────────┘

  TPad* pad_plot = new TPad("pad_plot", "", 0.00, 0.12, 0.88, 1.00);
  TPad* pad_cb   = new TPad("pad_cb", "", 0.88, 0.12, 0.96, 1.00);
  TPad* pad_info = new TPad("pad_info", "", 0.00, 0.00, 1.00, 0.12);

  // Bordi e margini del pad principale
  pad_plot->SetLeftMargin(0.10);
  pad_plot->SetRightMargin(0.015);
  pad_plot->SetTopMargin(0.11);
  pad_plot->SetBottomMargin(0.13);
  pad_plot->SetGridx();
  pad_plot->SetGridy();
  pad_plot->SetFrameLineWidth(2);

  pad_cb->SetLeftMargin(0.25);
  pad_cb->SetRightMargin(0.30);
  pad_cb->SetTopMargin(0.06);
  pad_cb->SetBottomMargin(0.13);

  pad_info->SetLeftMargin(0.02);
  pad_info->SetRightMargin(0.02);
  pad_info->SetTopMargin(0.05);
  pad_info->SetBottomMargin(0.05);

  canvas->cd();
  pad_plot->Draw();
  pad_cb->Draw();
  pad_info->Draw();

  //  Pad principale: frame + ellissi + centroidi

  pad_plot->cd();

  TH2D* frame = new TH2D("frame", "", 200, py_min, py_max, 200, pz_min, pz_max);
  frame->GetXaxis()->SetTitle("y_{det}  [mm]");
  frame->GetYaxis()->SetTitle("z_{det}  [mm]");
  frame->GetXaxis()->CenterTitle(true);
  frame->GetYaxis()->CenterTitle(true);
  frame->GetYaxis()->SetTitleOffset(1.55);
  frame->GetXaxis()->SetTitleOffset(1.20);
  frame->GetXaxis()->SetNdivisions(506);
  frame->GetYaxis()->SetNdivisions(505);
  frame->Draw("AXIS");

  // Indice del punto centrale (t più vicino a 0, ovvero il verde della palette)
  int i_center = N / 2;

  // Ordine di rendering: dal centro verso i bordi.
  // Il verde (centro) viene disegnato per primo e finisce sotto;
  // blu (i=0) e rosso (i=N-1) vengono disegnati per ultimi e rimangono visibili.
  // L'indice i/(N-1) mappa direttamente blu→rosso indipendentemente da t.
  std::vector<int> draw_order;
  draw_order.reserve(N);
  for (int d = 0; i_center - d >= 0 || i_center + d < N; ++d) {
    if (i_center - d >= 0)
      draw_order.push_back(i_center - d);
    if (d > 0 && i_center + d < N)
      draw_order.push_back(i_center + d);
  }

  std::vector<TEllipse*> ellipses;
  ellipses.reserve(N);

  for (int i : draw_order) {
    const auto& pt = trace[i];
    double t_norm  = (N > 1) ? static_cast<double>(i) / (N - 1) : 0.5;
    auto el        = cov_to_ellipse(pt.mu_y, pt.mu_z, pt.cov, cfg.sigma_scale);
    if (el.a < 1e-6 || el.b < 1e-6)
      continue;
    Int_t col   = palette_color(t_norm);
    TEllipse* e = new TEllipse(el.cy, el.cz, el.a, el.b, 0, 360, el.theta_deg);
    if (pt.valid) {
      e->SetFillColorAlpha(col, 0.50);
      e->SetLineColorAlpha(col, 0.75);
      e->SetLineStyle(1);
    } else {
      e->SetFillColorAlpha(col, 0.15); // Molto trasparente se "off detector"
      e->SetLineColorAlpha(col, 0.40);
      e->SetLineStyle(3); // Tratteggiato
    }
    e->SetLineWidth(1);
    e->Draw("same");
    ellipses.push_back(e);
  }

  // Linea del centroide
  TGraph* g_line = new TGraph(N, gy.data(), gz.data());
  g_line->SetLineColor(kBlack);
  g_line->SetLineWidth(2);
  g_line->Draw("L same");

  // Disegna fit se presente
  if (has_fit) {
    TLine* l_fit = nullptr;
    if (fit_res.axis == riptide::FitAxis::ZvsY) {
      double y_start = py_min;
      double y_end   = py_max;
      double z_start = fit_res.a * y_start + fit_res.b;
      double z_end   = fit_res.a * y_end + fit_res.b;
      if (std::abs(fit_res.a) > 1e-15) {
        if (z_start < pz_min) { y_start = (pz_min - fit_res.b) / fit_res.a; z_start = pz_min; }
        if (z_start > pz_max) { y_start = (pz_max - fit_res.b) / fit_res.a; z_start = pz_max; }
        if (z_end   < pz_min) { y_end   = (pz_min - fit_res.b) / fit_res.a; z_end   = pz_min; }
        if (z_end   > pz_max) { y_end   = (pz_max - fit_res.b) / fit_res.a; z_end   = pz_max; }
      }
      l_fit = new TLine(y_start, z_start, y_end, z_end);
    } else {
      double z_start = pz_min;
      double z_end   = pz_max;
      double y_start = fit_res.a * z_start + fit_res.b;
      double y_end   = fit_res.a * z_end + fit_res.b;
      if (std::abs(fit_res.a) > 1e-15) {
        if (y_start < py_min) { z_start = (py_min - fit_res.b) / fit_res.a; y_start = py_min; }
        if (y_start > py_max) { z_start = (py_max - fit_res.b) / fit_res.a; y_start = py_max; }
        if (y_end   < py_min) { z_end   = (py_min - fit_res.b) / fit_res.a; y_end   = py_min; }
        if (y_end   > py_max) { z_end   = (py_max - fit_res.b) / fit_res.a; y_end   = py_max; }
      }
      l_fit = new TLine(y_start, z_start, y_end, z_end);
    }
    l_fit->SetLineColor(kRed + 1);
    l_fit->SetLineStyle(2);
    l_fit->SetLineWidth(3);
    l_fit->Draw("same");
  }

  // Marker colorati: stessa logica centro → bordi
  for (int i : draw_order) {
    const auto& pt = trace[i];
    double t_norm  = (N > 1) ? static_cast<double>(i) / (N - 1) : 0.5;
    Int_t col      = palette_color(t_norm);
    TMarker* mk    = new TMarker(pt.mu_y, pt.mu_z, pt.valid ? 20 : 24);
    mk->SetMarkerColor(col);
    mk->SetMarkerSize(pt.valid ? 0.65 : 0.45);
    mk->Draw("same");
  }

  // Marker stella a inizio (blu) e fine (rosso) traccia
  {
    TMarker* mk_start = new TMarker(gy.front(), gz.front(), 29);
    mk_start->SetMarkerColor(palette_color(0.0));
    mk_start->SetMarkerSize(2.2);
    mk_start->Draw("same");

    TMarker* mk_end = new TMarker(gy.back(), gz.back(), 29);
    mk_end->SetMarkerColor(palette_color(1.0));
    mk_end->SetMarkerSize(2.2);
    mk_end->Draw("same");
  }

  // Titolo nel pad
  TLatex title_lat;
  title_lat.SetNDC();
  title_lat.SetTextFont(42);
  title_lat.SetTextSize(0.045);
  title_lat.SetTextAlign(22);
  title_lat.SetTextColor(kBlack);
  std::string title_str = "Traccia sul Detector";
  // Y = 1 - TopMargin/2  →  1 - 0.055 = 0.945, ben dentro il margine
  title_lat.DrawLatex(0.535, 0.945, title_str.c_str());

  // Etichette inizio/fine traccia
  TLatex lbl;
  lbl.SetNDC(false); // coordinate world
  lbl.SetTextFont(42);
  lbl.SetTextSize(0.030);
  lbl.SetTextColor(kBlack);

  pad_plot->RedrawAxis();

  //  Pad color bar: asse verticale con la palette (t da -L/2 a +L/2)

  pad_cb->cd();

  // Imposta range world del pad a [0,1]×[0,1] — TBox usa sempre coordinate world.
  // I margini del pad definiscono l'area utile della barra in unità world.
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

    TBox* box = new TBox(cb_x0, yb0, cb_x1, yb1); // world coords [0,1]
    box->SetFillColor(gStyle->GetColorPalette(i));
    box->SetLineWidth(0);
    box->Draw();
  }

  // Asse della color bar: disegnato in world coords del pad [0,1]×[0,1].
  // Il range fisico (t) è mappato su [cb_y0, cb_y1].
  double t_start = 0.0;
  double t_end   = total_L;

  TGaxis* cb_axis = new TGaxis(cb_x1, cb_y0,   // (xstart, ystart) in world coords
                               cb_x1, cb_y1,   // (xend,   yend)   in world coords
                               t_start, t_end, // wmin, wmax
                               505, "+L"       // ndiv, "+L" = label a destra dell'asse
  );
  cb_axis->SetLabelFont(42);
  cb_axis->SetLabelSize(0.18);
  cb_axis->SetTickSize(0.35);
  cb_axis->SetLabelOffset(0.03);
  cb_axis->SetTitle("t  [mm]");
  cb_axis->SetTitleFont(42);
  cb_axis->SetTitleSize(0.20);
  cb_axis->CenterTitle(kTRUE);
  cb_axis->SetTitleOffset(1.8);
  cb_axis->Draw();

  //  Pad info: riepilogo parametri fisici in due colonne

  pad_info->cd();

  // Sfondo grigio chiarissimo per distinguere l'area info
  pad_info->SetFillColor(TColor::GetColor(245, 245, 248));
  pad_info->SetFrameLineWidth(0);

  // Linea separatrice superiore
  TLine* sep_line = new TLine(0.0, 0.98, 1.0, 0.98);
  sep_line->SetNDC(true);
  sep_line->SetLineColor(kGray + 1);
  sep_line->SetLineWidth(1);
  sep_line->Draw();

  TLatex info;
  info.SetNDC(true);
  info.SetTextFont(42);

  // Tre colonne, due righe. Coordinate NDC del pad_info.
  double col1_x = 0.03, col2_x = 0.34, col3_x = 0.66;
  double hdr_y  = 0.80; // header colonna
  double row1_y = 0.52; // prima riga dati
  double row2_y = 0.18; // seconda riga dati

  // Header colonne (grigio, corsivo, piccolo)
  info.SetTextSize(0.13);
  info.SetTextColor(kGray + 2);
  info.SetTextAlign(12);
  info.DrawLatex(col1_x, hdr_y, "CONFIGURAZIONE LENTI");
  info.DrawLatex(col2_x, hdr_y, "PARAMETRI SORGENTE");
  info.DrawLatex(col3_x, hdr_y, "STATISTICHE PSF");

  // Dati (nero, dimensione normale)
  info.SetTextSize(0.18);
  info.SetTextColor(kBlack);

  // Colonna 1
  info.DrawLatex(col1_x, row1_y, ("x_{1} (lente 75 mm) = " + fmt(actual.x1) + " mm").c_str());
  info.DrawLatex(col1_x, row2_y, ("x_{2} (lente 60 mm) = " + fmt(actual.x2) + " mm").c_str());

  // Colonna 2
  if (use_3d) {
    info.DrawLatex(col2_x, row1_y,
                   ("P1: (" + fmt(p1.x) + "," + fmt(p1.y) + "," + fmt(p1.z) + ") mm").c_str());
    info.DrawLatex(col2_x, row2_y,
                   ("P2: (" + fmt(p2.x) + "," + fmt(p2.y) + "," + fmt(p2.z) + ") mm").c_str());
  } else {
    info.DrawLatex(col2_x, row1_y,
                   ("y_{0} = " + fmt(cfg.y0) + " mm,   L = " + fmt(cfg.L)
                    + " mm,   #Deltat = " + fmt(cfg.dt, 2) + " mm")
                       .c_str());
    info.DrawLatex(
        col2_x, row2_y,
        (std::to_string(N) + " punti traccia,   ellissi " + fmt(cfg.sigma_scale, 1) + "#sigma")
            .c_str());
  }

  // Colonna 3
  if (has_fit) {
    info.DrawLatex(col3_x, row1_y,
                   ("#chi^{2}/ndof = " + fmt(fit_res.chi2_ndof, 3) + "   (N = "
                    + std::to_string(fit_res.n_points_used) + "/" + std::to_string(N) + ")")
                       .c_str());
    std::string fit_str = (fit_res.axis == riptide::FitAxis::ZvsY)
                            ? "fit: z = " + fmt(fit_res.a, 3) + "y + " + fmt(fit_res.b, 3)
                            : "fit: y = " + fmt(fit_res.a, 3) + "z + " + fmt(fit_res.b, 3);
    info.DrawLatex(col3_x, row2_y, fit_str.c_str());
  } else {
    info.DrawLatex(col3_x, row1_y, ("#sigma_{y,max} = " + fmt(max_sigma_y, 4) + " mm").c_str());
    info.DrawLatex(col3_x, row2_y, ("#sigma_{z,max} = " + fmt(max_sigma_z, 4) + " mm").c_str());
  }

  //  Output

  canvas->Update();

  std::string out_path;
  if (!cfg.output.empty()) {
    std::filesystem::create_directories(std::filesystem::path(cfg.output).parent_path());
    out_path = cfg.output;
  } else {
    std::filesystem::create_directories("output/psf_analysis");
    std::ostringstream auto_name;
    auto_name << "output/psf_analysis/trace"
              << "_x1_" << std::fixed << std::setprecision(1) << actual.x1 << "_x2_" << actual.x2;
    if (use_3d) {
      auto_name << "_3d_" << std::hex << (static_cast<int>(p1.y + p2.z) & 0xFFF);
    } else {
      auto_name << "_y0_" << cfg.y0;
    }
    auto_name << ".png";
    out_path = auto_name.str();
  }

  // Salva in PNG ad alta risoluzione (300 DPI equiv. per canvas 1200×900)
  canvas->Print(out_path.c_str());
  std::cout << "Plot salvato in: " << out_path << "\n";

  // Cleanup
  for (auto* e : ellipses)
    delete e;
  delete g_line;
  delete frame;
  delete canvas;

  return 0;
}
