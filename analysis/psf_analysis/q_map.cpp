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

#include "plot_style_common.hpp"
#include "psf_interpolator.hpp"

#include <nlohmann/json.hpp>

#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TStyle.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

#include <omp.h>

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

  int n_jobs = 0; // 0 = tutti i core disponibili
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
    else if (key == "--jobs")
      cfg.n_jobs = std::stoi(next());
    else {
      std::cerr << "Opzione sconosciuta: " << key << "\n";
      std::exit(1);
    }
  }
  return cfg;
}

static std::string fmt(double v, int n = 1) {
  std::ostringstream o;
  o << std::fixed << std::setprecision(n) << v;
  return o.str();
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

  // Numero di thread OpenMP
  if (cli.n_jobs > 0)
    omp_set_num_threads(cli.n_jobs);
  int n_omp_threads = (cli.n_jobs > 0) ? cli.n_jobs : omp_get_max_threads();

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
  std::cout << "Thread OpenMP: " << n_omp_threads << "\n";

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

  std::string l1_label = "L1", l2_label = "L2";
  {
    const auto& key = db.begin()->first;
    if (!key.l1_id.empty())
      l1_label = "Lente " + key.l1_id;
    if (!key.l2_id.empty())
      l2_label = "Lente " + key.l2_id;
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
    struct CovEntry {
      double x1, x2;
      double coverage; // ∈ [0,1] oppure -1 se invalida
      int n_evaluated;
      bool config_valid;
    };
    std::vector<riptide::LensConfig> cov_configs;
    cov_configs.reserve(db.size());
    for (const auto& [k, _] : db) cov_configs.push_back(k);

    std::vector<CovEntry> entries(total_cfgs);

    int n_cfg_invalid = 0;
    int print_step    = std::max(1, total_cfgs / 20);

    std::atomic<int> cov_n_invalid{0};
    std::atomic<int> cov_processed{0};
    std::mutex cov_log_mtx;

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < total_cfgs; ++i) {
      const auto& cfg = cov_configs[i];
      riptide::QConfig qcfg_local = qcfg;
      qcfg_local.seed = qcfg.seed ^ (static_cast<unsigned>(i) * 2654435761u);
      riptide::CoverageResult res;
      bool failed = false;
      try {
        res = riptide::compute_coverage(cfg, db, qcfg_local);
      } catch (const std::exception& e) {
        std::lock_guard<std::mutex> lk(cov_log_mtx);
        std::cerr << "  [WARN] compute_coverage failed x1=" << cfg.x1 << " x2=" << cfg.x2 << ": "
                  << e.what() << "\n";
        entries[i] = {cfg.x1, cfg.x2, -1.0, 0, false};
        ++cov_n_invalid;
        failed = true;
      }
      if (failed) { ++cov_processed; continue; }

      if (!res.config_valid) ++cov_n_invalid;
      double cov_pct = res.config_valid ? res.coverage : -1.0;
      entries[i] = {cfg.x1, cfg.x2, cov_pct, res.n_y0_evaluated, res.config_valid};

      int done = ++cov_processed;
      if (done % print_step == 0 || done == total_cfgs) {
        std::lock_guard<std::mutex> lk(cov_log_mtx);
        std::cout << "  [" << std::setw(3) << done << "/" << total_cfgs << "]"
                  << "  x1=" << fmt(cfg.x1) << "  x2=" << fmt(cfg.x2);
        if (res.config_valid)
          std::cout << "  cov=" << fmt(res.coverage * 100.0, 1) << "%";
        else
          std::cout << "  [INVALIDA]";
        std::cout << "\n";
      }
    }
    n_cfg_invalid = cov_n_invalid.load();

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
    set_root_style();
    set_viridis_palette(false);

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

    h_cov.GetXaxis()->SetTitle(("x_{1}  [mm]  (" + l1_label + ")").c_str());
    h_cov.GetYaxis()->SetTitle(("x_{2}  [mm]  (" + l2_label + ")").c_str());
    h_cov.GetXaxis()->CenterTitle(true);
    h_cov.GetYaxis()->CenterTitle(true);
    h_cov.GetXaxis()->SetTitleOffset(1.20);
    h_cov.GetYaxis()->SetTitleOffset(1.50);
    h_cov.GetXaxis()->SetNdivisions(506);
    h_cov.GetYaxis()->SetNdivisions(505);
    h_cov.GetZaxis()->SetTitle("copertura [%]");
    apply_zaxis_style(&h_cov);

    // Range colori [0, 100] — percentuale intera
    h_cov.SetMinimum(0.0);
    h_cov.SetMaximum(100.0);

    h_inv.SetFillColor(invalid_color);
    h_inv.SetLineColor(invalid_color);

    TCanvas* c = make_map_canvas("coverage_map", cli.log_scale);
    h_cov.Draw("COLZ");
    h_inv.Draw("BOX same");
    c->Update();
    draw_na_legend(&h_cov, invalid_color);

    draw_map_title("Copertura geometrica  coverage(x_{1}, x_{2})  [%]");

    std::filesystem::create_directories(std::filesystem::path(cli.output_path).parent_path());
    c->Print(cli.output_path.c_str());
    std::cout << "\nMappa salvata in: " << cli.output_path << "\n";
    return 0;
  }

  // MODALITÀ Q (default)
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
  std::vector<riptide::LensConfig> q_configs;
  q_configs.reserve(db.size());
  for (const auto& [k, _] : db) q_configs.push_back(k);

  std::vector<QEntry> entries(total_cfgs);

  int print_step    = std::max(1, total_cfgs / 20);
  int n_cfg_invalid = 0;

  std::atomic<int> q_n_invalid{0};
  std::atomic<int> q_processed{0};
  std::mutex q_log_mtx;

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < total_cfgs; ++i) {
    const auto& cfg = q_configs[i];

    // Seed deterministico e distinto per configurazione (costante di Knuth)
    riptide::QConfig qcfg_local = qcfg;
    qcfg_local.seed = qcfg.seed ^ (static_cast<unsigned>(i) * 2654435761u);

    riptide::QResult res;
    bool failed = false;
    try {
      res = riptide::compute_Q(cfg, db, qcfg_local);
    } catch (const std::exception& e) {
      std::lock_guard<std::mutex> lk(q_log_mtx);
      std::cerr << "  [WARN] compute_Q failed x1=" << cfg.x1 << " x2=" << cfg.x2 << ": " << e.what()
                << "\n";
      entries[i] = {cfg.x1, cfg.x2, 0.0, 0.0, 1.0, 0.0, 0, 0, 0, false};
      ++q_n_invalid;
      failed = true;
    }
    if (failed) { ++q_processed; continue; }

    double target = (cli.dist_target > 0.0) ? cli.dist_target : res.Q_target;
    double metric = cli.dist_to_target ? std::abs(res.Q - target) : res.Q;

    if (!res.config_valid) ++q_n_invalid;
    entries[i] = {cfg.x1, cfg.x2, metric, res.Q, target, res.rho_estimate, res.n_traces,
                  res.n_failed, res.n_invalid, res.config_valid};

    int done = ++q_processed;
    if (done % print_step == 0 || done == total_cfgs) {
      std::lock_guard<std::mutex> lk(q_log_mtx);
      std::cout << "  [" << std::setw(3) << done << "/" << total_cfgs << "]"
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
  }
  n_cfg_invalid = q_n_invalid.load();

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
  set_root_style();
  set_viridis_palette(true);

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

  h_Q.GetXaxis()->SetTitle(("x_{1}  [mm]  (" + l1_label + ")").c_str());
  h_Q.GetYaxis()->SetTitle(("x_{2}  [mm]  (" + l2_label + ")").c_str());
  h_Q.GetXaxis()->CenterTitle(true);
  h_Q.GetYaxis()->CenterTitle(true);
  h_Q.GetXaxis()->SetTitleOffset(1.20);
  h_Q.GetYaxis()->SetTitleOffset(1.50);
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

  std::string z_lbl_cb;
  if (cli.dist_to_target)
    z_lbl_cb = "|#LT#chi_{red}^{2}#GT - Q_{target}|";
  else
    z_lbl_cb = "#LT#chi_{red}^{2}#GT";
  h_Q.GetZaxis()->SetTitle(z_lbl_cb.c_str());
  apply_zaxis_style(&h_Q);

  TCanvas* c = make_map_canvas("q_map", cli.log_scale);
  h_Q.Draw("COLZ");
  h_inv.Draw("BOX same");
  c->Update();
  draw_na_legend(&h_Q, invalid_color);

  if (cli.dist_to_target)
    draw_map_title("Risoluzione tracce 3D  -  distanza da target");
  else
    draw_map_title("Risoluzione tracce 3D  #LT#chi^{2}_{red}#GT(x_{1}, x_{2})");

  std::filesystem::create_directories(std::filesystem::path(cli.output_path).parent_path());
  c->Print(cli.output_path.c_str());
  std::cout << "\nMappa salvata in: " << cli.output_path << "\n";
  return 0;
}
