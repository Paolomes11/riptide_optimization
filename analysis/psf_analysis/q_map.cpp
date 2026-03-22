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
 * mappa TH2D analoga a efficiency2D.png di plot2D.
 *
 * Funzionalità aggiuntive rispetto a plot2D:
 *   - Marker sulla configurazione ottimale (minimo di Q)
 *   - Scala logaritmica opzionale (--log)
 *   - Normalizzazione Q/n_traces (chi^2 medio per traccia, --norm)
 *   - Esportazione TSV con i valori di Q per uso esterno (--tsv)
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
#include <TPaletteAxis.h>
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

  double y0_min = -1.0; // -1 = usa config.json
  double y0_max = -1.0;
  double dy0    = -1.0;
  double L      = 10.0;
  double dt     = 0.1;

  bool log_scale = false;
  bool normalize = false; // Q / n_traces
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
      std::cerr << "Uso: q_map [--psf path] [--config path] [--output path]\n"
                << "          [--tsv path] [--y0-min val] [--y0-max val] [--dy0 val]\n"
                << "          [--L val] [--dt val] [--log] [--norm]\n";
      std::exit(1);
    }
  }
  return cfg;
}

// Stile ROOT publication-quality
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

  // Palette: kBird è publication-standard e percettivamente uniforme
  gStyle->SetPalette(kBird);
  gStyle->SetNumberContours(255);
}

// Formattatore double
static std::string fmt(double v, int n = 1) {
  std::ostringstream o;
  o << std::fixed << std::setprecision(n) << v;
  return o.str();
}

// main
int main(int argc, char** argv) {
  CliConfig cli = parse_args(argc, argv);

  // 1. Legge config.json per i parametri di griglia
  using json = nlohmann::json;
  std::ifstream f_cfg(cli.config_path);
  if (!f_cfg.is_open()) {
    std::cerr << "Errore: impossibile aprire " << cli.config_path << "\n";
    return 1;
  }
  json jcfg;
  f_cfg >> jcfg;

  // Parametri griglia lenti (usati solo per ricostruire i bin dell'istogramma)
  const double dx   = jcfg["dx"];
  const double r1   = jcfg["r1"];
  const double h1   = jcfg["h1"];
  const double r2   = jcfg["r2"];
  const double h2   = jcfg["h2"];
  const double xmin = jcfg["x_min"];
  const double xmax = jcfg["x_max"];

  // Limiti x1 e x2 identici a quelli usati in lens_scan / optimizer
  const double x1_lo = xmin - r1 + h1;
  const double x1_hi = xmax - h2 - 3.0 - r1;
  const double x2_lo = x1_lo + r1 + r2 + 3.0; // x2_min del primo x1
  const double x2_hi = xmax + r2 - h2;

  // 2. Parametri di campionamento y0
  riptide::QConfig qcfg;
  qcfg.y0_min   = (cli.y0_min >= 0.0) ? cli.y0_min : 0.0;
  qcfg.y0_max   = (cli.y0_max >= 0.0) ? cli.y0_max : 10.0;
  qcfg.dy0      = (cli.dy0 > 0.0) ? cli.dy0 : 1.0;
  qcfg.trace_L  = cli.L;
  qcfg.trace_dt = cli.dt;

  std::cout << "Campionamento y0: [" << qcfg.y0_min << ", " << qcfg.y0_max << "] passo " << qcfg.dy0
            << " mm  →  "
            << static_cast<int>(std::round((qcfg.y0_max - qcfg.y0_min) / qcfg.dy0)) + 1
            << " tracce per configurazione\n";
  std::cout << "Traccia: L=" << qcfg.trace_L << " mm, dt=" << qcfg.trace_dt << " mm\n";

  // 3. Carica il database PSF
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

  // Statistiche globali sulle esclusioni (per il report finale)
  struct ExclusionSummary {
    double x1, x2;
    int n_failed;
    int n_total;
    // y0 esclusi raggruppati per tipo
    std::vector<double> y0_not_converged;
    std::vector<double> y0_build_failed;
    std::vector<double> y0_fit_failed;
  };
  std::vector<ExclusionSummary> exclusions; // solo cfg con n_failed > 0

  int config_idx = 0;
  int total_cfgs = static_cast<int>(db.size());
  int print_step = std::max(1, total_cfgs / 20);

  // QConfig con verbose=false: gli warning vengono restituiti nel QResult
  qcfg.fit_max_iter = 20;
  qcfg.fit_tol      = 1e-8;
  // (verbose resta false — non vogliamo spam su stderr)

  for (const auto& [cfg, _] : db) {
    riptide::QResult res = riptide::compute_Q(cfg, db, qcfg);

    double Q_val = res.Q;
    if (cli.normalize && res.n_traces > 0)
      Q_val /= static_cast<double>(res.n_traces);

    entries.push_back({cfg.x1, cfg.x2, Q_val, res.n_traces, res.n_failed});

    // Raccoglie info sulle esclusioni per il report finale
    if (res.n_failed > 0) {
      ExclusionSummary es;
      es.x1       = cfg.x1;
      es.x2       = cfg.x2;
      es.n_failed = res.n_failed;
      es.n_total  = res.n_traces + res.n_failed;
      for (const auto& w : res.warnings) {
        switch (w.kind) {
        case riptide::QWarning::Kind::FitNotConverged:
          es.y0_not_converged.push_back(w.y0);
          break;
        case riptide::QWarning::Kind::BuildTraceFailed:
          es.y0_build_failed.push_back(w.y0);
          break;
        case riptide::QWarning::Kind::FitFailed:
          es.y0_fit_failed.push_back(w.y0);
          break;
        }
      }
      exclusions.push_back(std::move(es));
    }

    if ((config_idx + 1) % print_step == 0 || config_idx + 1 == total_cfgs) {
      // Riga di progresso: mostra anche le esclusioni di questa cfg
      std::string excl_str = (res.n_failed > 0)
                               ? "  [" + std::to_string(res.n_failed) + "/"
                                     + std::to_string(res.n_traces + res.n_failed) + " y0 esclusi]"
                               : "";
      std::cout << "  [" << std::setw(3) << (config_idx + 1) << "/" << total_cfgs << "]"
                << "  x1=" << fmt(cfg.x1) << "  x2=" << fmt(cfg.x2) << "  Q=" << std::scientific
                << std::setprecision(3) << Q_val << std::defaultfloat << excl_str << "\n";
    }
    ++config_idx;
  }

  // 4b. Report finale sulle esclusioni
  if (!exclusions.empty()) {
    std::cout << "\n─── Riepilogo esclusioni (" << exclusions.size()
              << " configurazioni con y0 esclusi) ───\n";

    // Conta tipi di esclusione globali
    int tot_not_conv = 0, tot_build = 0, tot_fit = 0;
    for (const auto& es : exclusions) {
      tot_not_conv += static_cast<int>(es.y0_not_converged.size());
      tot_build += static_cast<int>(es.y0_build_failed.size());
      tot_fit += static_cast<int>(es.y0_fit_failed.size());
    }

    std::cout << "  Tipo esclusione          | Occorrenze\n";
    std::cout << "  ─────────────────────────┼───────────\n";
    if (tot_not_conv > 0)
      std::cout << "  Fit non convergente      | " << tot_not_conv << "\n";
    if (tot_build > 0)
      std::cout << "  build_trace fallita      | " << tot_build << "\n";
    if (tot_fit > 0)
      std::cout << "  fit_trace eccezione      | " << tot_fit << "\n";

    // Mostra le cfg più problematiche (quelle con più y0 esclusi)
    std::sort(exclusions.begin(), exclusions.end(),
              [](const ExclusionSummary& a, const ExclusionSummary& b) {
                return a.n_failed > b.n_failed;
              });

    const int max_show = std::min(10, static_cast<int>(exclusions.size()));
    std::cout << "\n  Top-" << max_show << " configurazioni con più esclusioni:\n";
    std::cout << "  " << std::left << std::setw(8) << "x1 [mm]" << std::setw(10) << "x2 [mm]"
              << std::setw(14) << "esclusi/totali"
              << "y0 esclusi\n";
    std::cout << "  " << std::string(60, '-') << "\n";

    for (int i = 0; i < max_show; ++i) {
      const auto& es = exclusions[i];
      std::cout << "  " << std::left << std::setw(8) << fmt(es.x1, 1) << std::setw(10)
                << fmt(es.x2, 1) << std::setw(14)
                << (std::to_string(es.n_failed) + "/" + std::to_string(es.n_total));

      // Elenca i y0 esclusi in forma compatta (max 8 per riga)
      auto print_y0_list = [](const std::vector<double>& v, const std::string& tag) {
        if (v.empty())
          return;
        std::cout << tag << ": ";
        int shown = 0;
        for (double y : v) {
          if (shown++ >= 8) {
            std::cout << "...";
            break;
          }
          std::cout << y;
          if (shown < static_cast<int>(v.size()) && shown < 8)
            std::cout << ", ";
        }
        std::cout << "  ";
      };
      print_y0_list(es.y0_not_converged, "no-conv");
      print_y0_list(es.y0_build_failed, "build");
      print_y0_list(es.y0_fit_failed, "fit-exc");
      std::cout << "\n";
    }
    std::cout << "─────────────────────────────────────────────────────────────\n\n";
  } else {
    std::cout << "\n✓  Nessuna esclusione: tutte le tracce hanno converguto.\n\n";
  }

  // 5. Trova il minimo di Q
  auto it_best = std::min_element(entries.begin(), entries.end(),
                                  [](const QEntry& a, const QEntry& b) { return a.Q < b.Q; });
  std::cout << "\n★  Configurazione ottimale:\n"
            << "   x1 = " << fmt(it_best->x1, 2) << " mm\n"
            << "   x2 = " << fmt(it_best->x2, 2) << " mm\n"
            << "   Q  = " << it_best->Q << "\n";

  // 6. Eventuale esportazione TSV
  if (!cli.tsv_path.empty()) {
    std::filesystem::create_directories(std::filesystem::path(cli.tsv_path).parent_path());
    std::ofstream tsv(cli.tsv_path);
    tsv << "x1\tx2\tQ\tn_traces\tn_failed\n";
    for (const auto& e : entries)
      tsv << e.x1 << "\t" << e.x2 << "\t" << e.Q << "\t" << e.n_traces << "\t" << e.n_failed
          << "\n";
    std::cout << "TSV salvato in: " << cli.tsv_path << "\n";
  }

  // 7. Costruisce l'istogramma 2D
  apply_style();

  // Dimensioni della griglia (stessa logica di plot2D)
  int bins_x = std::max(1, static_cast<int>(std::ceil((x1_hi - x1_lo) / dx)) + 1);
  int bins_y = std::max(1, static_cast<int>(std::ceil((x2_hi - x2_lo) / dx)) + 1);

  // Margini dei bin: ±dx/2 per centrare i punti nelle celle
  double hx = dx / 2.0;
  TH2D h_Q("h_Q", "", bins_x, x1_lo - hx, x1_hi + hx, bins_y, x2_lo - hx, x2_hi + hx);

  // Scala logaritmica: ROOT gestisce il log sull'asse Z tramite SetLogz()
  // ma richiede che il TH2D contenga valori > 0 e che il range sia impostato.
  for (const auto& e : entries) {
    double val = e.Q;
    if (cli.log_scale) {
      // Protezione contro Q == 0 (non dovrebbe accadere con dati reali)
      val = (val > 0.0) ? val : std::numeric_limits<double>::min();
    }
    h_Q.Fill(e.x1, e.x2, val);
  }

  // Imposta le etichette degli assi
  h_Q.GetXaxis()->SetTitle("x_{1}  [mm]  (lente 75 mm)");
  h_Q.GetYaxis()->SetTitle("x_{2}  [mm]  (lente 60 mm)");
  h_Q.GetXaxis()->CenterTitle(true);
  h_Q.GetYaxis()->CenterTitle(true);
  h_Q.GetXaxis()->SetTitleOffset(1.20);
  h_Q.GetYaxis()->SetTitleOffset(1.55);

  // Range colori: usa i percentili 5–95 per evitare che outlier estremi
  // dominino la scala cromatica, rendendo invisibili le differenze centrali.
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

  // 8. Canvas e layout
  TCanvas* canvas = new TCanvas("canvas", "Q map", 1200, 1000);
  canvas->SetLeftMargin(0.0);
  canvas->SetRightMargin(0.0);
  canvas->SetTopMargin(0.0);
  canvas->SetBottomMargin(0.0);

  TPad* pad_plot = new TPad("pad_plot", "", 0.00, 0.13, 0.88, 1.00);
  TPad* pad_cb   = new TPad("pad_cb", "", 0.88, 0.13, 0.97, 1.00);
  TPad* pad_info = new TPad("pad_info", "", 0.00, 0.00, 1.00, 0.13);

  pad_plot->SetLeftMargin(0.14);
  pad_plot->SetRightMargin(0.015);
  pad_plot->SetTopMargin(0.10);
  pad_plot->SetBottomMargin(0.13);
  pad_plot->SetGridx();
  pad_plot->SetGridy();
  pad_plot->SetFrameLineWidth(2);
  if (cli.log_scale)
    pad_plot->SetLogz();

  pad_cb->SetLeftMargin(0.25);
  pad_cb->SetRightMargin(0.28);
  pad_cb->SetTopMargin(0.10);
  pad_cb->SetBottomMargin(0.13);

  pad_info->SetLeftMargin(0.02);
  pad_info->SetRightMargin(0.02);
  pad_info->SetTopMargin(0.08);
  pad_info->SetBottomMargin(0.08);

  canvas->cd();
  pad_plot->Draw();
  pad_cb->Draw();
  pad_info->Draw();

  // 9. Disegna la mappa
  pad_plot->cd();
  h_Q.Draw("COLZ");
  canvas->Update(); // forza ROOT a generare la palette

  // Riposiziona la palette automatica generata da COLZ
  TPaletteAxis* pal = (TPaletteAxis*)h_Q.GetListOfFunctions()->FindObject("palette");
  if (pal) {
    pal->SetX1NDC(0.86);
    pal->SetX2NDC(0.92);
    pal->SetY1NDC(pad_plot->GetBottomMargin());
    pal->SetY2NDC(1.0 - pad_plot->GetTopMargin());
    pal->GetAxis()->SetLabelFont(42);
    pal->GetAxis()->SetLabelSize(0.035);
    pal->GetAxis()->SetTitleOffset(1.40);
  }

  // Marker stella sul minimo di Q
  {
    TMarker* mk = new TMarker(it_best->x1, it_best->x2, 29); // stella
    mk->SetMarkerColor(kRed);
    mk->SetMarkerSize(2.8);
    mk->Draw("same");

    // Etichetta testo vicino al minimo
    TLatex lbl;
    lbl.SetTextFont(42);
    lbl.SetTextSize(0.032);
    lbl.SetTextColor(kRed + 1);
    lbl.SetTextAlign(12);
    std::string lbl_str = " #bf{min} (" + fmt(it_best->x1, 1) + ", " + fmt(it_best->x2, 1) + ")";
    lbl.DrawLatex(it_best->x1, it_best->x2, lbl_str.c_str());
  }

  // Titolo nel margine superiore del pad_plot
  TLatex title_lat;
  title_lat.SetNDC();
  title_lat.SetTextFont(42);
  title_lat.SetTextSize(0.048);
  title_lat.SetTextAlign(22);
  title_lat.SetTextColor(kBlack);
  std::string z_label = cli.normalize ? "Q / n_{tracce}" : "Q(x_{1}, x_{2})";
  title_lat.DrawLatex(0.535, 0.955,
                      ("Mappa della funzione di qualit#grave{a}  " + z_label).c_str());

  pad_plot->RedrawAxis();

  // 10. Color bar manuale nel pad_cb
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

  // Asse della color bar
  double cb_zmin = h_Q.GetMinimum();
  double cb_zmax = h_Q.GetMaximum();
  TGaxis* cb_axis =
      new TGaxis(cb_x1, cb_y0, cb_x1, cb_y1, cb_zmin, cb_zmax, 505, cli.log_scale ? "+LG" : "+L");
  cb_axis->SetLabelFont(42);
  cb_axis->SetLabelSize(0.16);
  cb_axis->SetTickSize(0.35);
  cb_axis->SetLabelOffset(0.03);
  cb_axis->SetTitle(z_label.c_str());
  cb_axis->SetTitleFont(42);
  cb_axis->SetTitleSize(0.18);
  cb_axis->SetTitleOffset(0.60);
  cb_axis->Draw();

  // 11. Info pad
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

  // Layout tre colonne
  const double col1 = 0.03, col2 = 0.34, col3 = 0.67;
  const double hdr = 0.80, row1 = 0.50, row2 = 0.15;

  info.SetTextSize(0.13);
  info.SetTextColor(kGray + 2);
  info.SetTextAlign(12);
  info.DrawLatex(col1, hdr, "PARAMETRI PSF");
  info.DrawLatex(col2, hdr, "CAMPIONAMENTO y_{0}");
  info.DrawLatex(col3, hdr, "OTTIMO");

  info.SetTextSize(0.20);
  info.SetTextColor(kBlack);

  // Colonna 1: file sorgente
  {
    std::string psf_base = std::filesystem::path(cli.psf_path).filename().string();
    info.DrawLatex(col1, row1, ("Sorgente PSF: " + psf_base).c_str());
    info.DrawLatex(
        col1, row2,
        (std::to_string(total_cfgs) + " configurazioni   dx = " + fmt(dx, 1) + " mm").c_str());
  }

  // Colonna 2: campionamento
  {
    int n_y0 = static_cast<int>(std::round((qcfg.y0_max - qcfg.y0_min) / qcfg.dy0)) + 1;
    info.DrawLatex(col2, row1,
                   ("y_{0} #in [" + fmt(qcfg.y0_min, 1) + ", " + fmt(qcfg.y0_max, 1)
                    + "] mm   #Deltay_{0} = " + fmt(qcfg.dy0, 1) + " mm")
                       .c_str());
    info.DrawLatex(col2, row2,
                   ("L = " + fmt(qcfg.trace_L, 1) + " mm   #Deltat = " + fmt(qcfg.trace_dt, 2)
                    + " mm   (" + std::to_string(n_y0) + " tracce/config)")
                       .c_str());
  }

  // Colonna 3: minimo
  {
    info.DrawLatex(col3, row1,
                   ("#bf{x_{1}^{*}} = " + fmt(it_best->x1, 1)
                    + " mm,   #bf{x_{2}^{*}} = " + fmt(it_best->x2, 1) + " mm")
                       .c_str());
    std::ostringstream q_str;
    q_str << std::scientific << std::setprecision(3) << it_best->Q;
    info.DrawLatex(col3, row2, ("#bf{Q_{min}} = " + q_str.str()).c_str());
  }

  // 12. Salva
  canvas->Update();

  std::filesystem::create_directories(std::filesystem::path(cli.output_path).parent_path());
  canvas->Print(cli.output_path.c_str());
  std::cout << "\nMappa salvata in: " << cli.output_path << "\n";

  delete canvas;
  return 0;
}
