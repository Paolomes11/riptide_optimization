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
 * exp1_main — Pipeline di analisi per immagini FITS a 16 bit
 *
 * Uso:
 *   exp1_main [--data-dir <path>] [--ssd <mount>]
 *             [--good <subdir>] [--bad <subdir>] [--bg <subdir>]
 *             [--output <path>] [--frames N]
 *             [--roi x0 y0 x1 y1]
 *             [--sigma N] [--iter N]
 *             [--method sigma|mean|median]
 *             [--no-root] [--no-png]
 *
 * Struttura attesa dati:
 *   <data-dir>/
 *     good/         ← 10 frame .fit configurazione "buona"
 *     bad/          ← 10 frame .fit configurazione "cattiva"
 *     background/   ← 10 frame .fit fondo
 *
 * Output (in --output, default: output/exp1/):
 *   stacked_means.png    — mappe 2D delle 3 serie stackate
 *   sigma_maps.png       — mappe 2D delle incertezze
 *   diff_maps.png        — differenza segnale - fondo (good e bad)
 *   profiles.png         — profili integrati X e Y + differenza good-bad
 *   summary.png          — pannello riassuntivo con ΔI e significatività
 */

#include "fit_reader.hpp"
#include "frame_stats.hpp"
#include "image_analysis.hpp"

#include <TFile.h>
#include <TTree.h>

#include <nlohmann/json.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>

// Parsing CLI minimale

struct CliOptions {
  std::filesystem::path data_dir = "analysis/exp1/data";
  std::string good_dir           = "good";
  std::string bad_dir            = "bad";
  std::string bg_dir             = "background";
  std::filesystem::path output   = "output/exp1";
  std::string ssd_mount          = "";
  size_t max_frames              = 0;       // 0 = tutti
  std::string method             = "sigma"; // sigma | mean | median
  double clip_sigma              = 3.0;
  int clip_iter                  = 3;
  double z_min                   = 0.005; // 0.5%
  double z_max                   = 0.995; // 99.5%
  bool save_png                  = true;
  bool save_root                 = true;
  exp1::ROI roi; // default: intero frame

  std::string lens_path   = "";
  std::string config_path = "config/config.json";
  double x1_good          = 0.0;
  double x2_good          = 0.0;
  double x1_bad           = 0.0;
  double x2_bad           = 0.0;
};

static void print_usage(const char* prog) {
  std::cout << "Uso: " << prog << " [opzioni]\n\n"
            << "Opzioni:\n"
            << "  --data-dir <path>      Cartella radice dati (default: analysis/exp1/data)\n"
            << "  --ssd <mount>          Mount point SSD (sovrascrive --data-dir)\n"
            << "                           Atteso: <mount>/exp1/{good,bad,background}/\n"
            << "  --good <subdir>        Sottocartella good (default: good)\n"
            << "  --bad  <subdir>        Sottocartella bad  (default: bad)\n"
            << "  --bg   <subdir>        Sottocartella bg   (default: background)\n"
            << "  --output <path>        Cartella output    (default: output/exp1)\n"
            << "  --frames N             Massimo frame per serie (0 = tutti, default: 0)\n"
            << "  --method [sigma|mean|median]  Metodo stacking (default: sigma)\n"
            << "  --sigma <val>          Soglia sigma-clipping (default: 3.0)\n"
            << "  --iter  <N>            Iterazioni sigma-clipping (default: 3)\n"
            << "  --z-min <perc>         Percentile inferiore scala Z (0.0-1.0, default: 0.005)\n"
            << "  --z-max <perc>         Percentile superiore scala Z (0.0-1.0, default: 0.995)\n"
            << "  --roi x0 y0 x1 y1     ROI per integrazione (default: intero frame)\n"
            << "                           (x1=-1, y1=-1 = fino al bordo)\n"
            << "  --no-root              Non salvare file .root\n"
            << "  --no-png               Non salvare file .png\n"
            << "  --lens <path>          lens.root per confronto simulazione (opzionale)\n"
            << "  --sim-config <path>    config.json per n_photons [config/config.json]\n"
            << "  --x1-good F            x1 configurazione good [mm]\n"
            << "  --x2-good F            x2 configurazione good [mm]\n"
            << "  --x1-bad  F            x1 configurazione bad  [mm]\n"
            << "  --x2-bad  F            x2 configurazione bad  [mm]\n"
            << "  --help                 Mostra questo messaggio\n";
}

static CliOptions parse_args(int argc, char** argv) {
  CliOptions opt;
  bool roi_set = false;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    auto next       = [&]() -> std::string {
      if (i + 1 >= argc) {
        std::cerr << "Argomento mancante dopo " << arg << "\n";
        std::exit(1);
      }
      return argv[++i];
    };
    auto next_int = [&]() { return std::stoi(next()); };
    auto next_dbl = [&]() { return std::stod(next()); };

    if (arg == "--help" || arg == "-h") {
      print_usage(argv[0]);
      std::exit(0);
    } else if (arg == "--data-dir")
      opt.data_dir = next();
    else if (arg == "--ssd") {
      std::string mount = next();
      opt.data_dir      = std::filesystem::path(mount) / "exp1";
    } else if (arg == "--good")
      opt.good_dir = next();
    else if (arg == "--bad")
      opt.bad_dir = next();
    else if (arg == "--bg")
      opt.bg_dir = next();
    else if (arg == "--output")
      opt.output = next();
    else if (arg == "--frames")
      opt.max_frames = static_cast<size_t>(next_int());
    else if (arg == "--method")
      opt.method = next();
    else if (arg == "--sigma")
      opt.clip_sigma = next_dbl();
    else if (arg == "--iter")
      opt.clip_iter = next_int();
    else if (arg == "--z-min")
      opt.z_min = next_dbl();
    else if (arg == "--z-max")
      opt.z_max = next_dbl();
    else if (arg == "--roi") {
      opt.roi.x0 = next_int();
      opt.roi.y0 = next_int();
      opt.roi.x1 = next_int();
      opt.roi.y1 = next_int();
      roi_set    = true;
    } else if (arg == "--no-root")
      opt.save_root = false;
    else if (arg == "--no-png")
      opt.save_png = false;
    else if (arg == "--lens")
      opt.lens_path = next();
    else if (arg == "--sim-config")
      opt.config_path = next();
    else if (arg == "--x1-good")
      opt.x1_good = next_dbl();
    else if (arg == "--x2-good")
      opt.x2_good = next_dbl();
    else if (arg == "--x1-bad")
      opt.x1_bad = next_dbl();
    else if (arg == "--x2-bad")
      opt.x2_bad = next_dbl();
    else {
      std::cerr << "Opzione sconosciuta: " << arg << "\n";
      print_usage(argv[0]);
      std::exit(1);
    }
  }

  (void)roi_set; // ROI con tutti -1 = intero frame
  return opt;
}

// Funzione di stacking secondo il metodo scelto

static exp1::StackedImage do_stack(const std::vector<exp1::FitsFrame>& frames,
                                   const CliOptions& opt) {
  exp1::StackConfig cfg;
  cfg.n_sigma = opt.clip_sigma;
  cfg.n_iter  = opt.clip_iter;

  if (opt.method == "mean") {
    std::cout << "[STACK] Metodo: mean (senza sigma-clipping)\n";
    return exp1::mean_stack(frames);
  } else if (opt.method == "median") {
    std::cout << "[STACK] Metodo: median + MAD\n";
    return exp1::median_stack(frames);
  } else {
    std::cout << "[STACK] Metodo: sigma-clipping (" << opt.clip_sigma << "σ, " << opt.clip_iter
              << " iterazioni)\n";
    return exp1::sigma_clip_stack(frames, cfg);
  }
}

// ─── Fattore correttivo geometrico sorgente rettangolare → sferica ──────────
//
// La sorgente GPS è un rettangolo 2·src_xh × 2·src_zh mm a Y = src_y mm sopra
// l'asse ottico (X), che emette in un emisfero verso le lenti (a X = x1_i).
// La sorgente sperimentale è un bulbo semisferico di pochi mm centrato in
// origine (0,0,0), assimilato a una sorgente puntiforme.
//
// L'angolo solido parassiale da un punto (xs, src_y, zs) verso l'apertura
// della prima lente in x = x1 è proporzionale a (x1-xs)/d³, con d = distanza.
// Mediando su z ∈ [-src_zh, src_zh]:
//
//   Ω̄_z(xs; x1) ∝ (x1-xs) / [(a²+src_y²)·√(a²+src_y²+src_zh²)]   a = x1-xs
//
// Per src_y = 5 mm, src_zh = 5 mm i termini costanti diventano 25 e 50.
// L'antiderivata di [2·src_zh·a / ((a²+src_y²)√(a²+src_y²+src_zh²))] rispetto ad a
// (sostituzione t = √(a²+src_y²+src_zh²)) è:
//
//   L(a) = ln|(√(a²+50) − 5) / (√(a²+50) + 5)|   [src_y=src_zh=5 → 50=25+25]
//
// Mediando anche su xs ∈ [-src_xh, src_xh]:
//
//   <Ω_rect>(x1) ∝ (1/(4·src_xh·src_zh)) · [L(x1+src_xh) − L(x1−src_xh)]
//
// Ω_sphere(x1) ∝ 1/x1²   (sorgente puntiforme in origine)
//
// Fattore correttivo per config i:  f_i = 4·src_xh·src_zh / (x1_i² · ΔL_i)
// Correzione al rapporto:           C = f_good / f_bad

// L(a) = ln|(√(a²+src_y²+src_zh²)−src_zh) / (√(a²+src_y²+src_zh²)+src_zh)|
// hardcoded per src_y = src_zh = 5 mm → costanti 25+25 = 50
static double geom_L(double a) {
  double t = std::sqrt(a * a + 50.0);
  return std::log((t - 5.0) / (t + 5.0));
}

static double geom_correction(double x1, double src_xh, double src_zh) {
  double delta = geom_L(x1 + src_xh) - geom_L(x1 - src_xh);
  if (delta <= 0.0) return 1.0;
  return (4.0 * src_xh * src_zh) / (x1 * x1 * delta);
}

// main

int main(int argc, char** argv) {
  CliOptions opt = parse_args(argc, argv);

  std::cout << "═══════════════════════════════════════════════════════════\n";
  std::cout << "  exp1 — Analisi immagini FITS 16 bit\n";
  std::cout << "═══════════════════════════════════════════════════════════\n";
  std::cout << "  Dati: " << opt.data_dir.string() << "\n";
  std::cout << "  Output: " << opt.output.string() << "\n";

  try {
    // 1. Lettura dati (carichiamo tutto in RAM per l'analisi frame-to-frame)

    std::cout << "\n--- Caricamento serie FITS ---\n";
    auto good_frames = exp1::read_fits_stack(opt.data_dir / opt.good_dir, opt.max_frames);
    auto bad_frames  = exp1::read_fits_stack(opt.data_dir / opt.bad_dir, opt.max_frames);
    auto bg_frames   = exp1::read_fits_stack(opt.data_dir / opt.bg_dir, opt.max_frames);

    // 2. Analisi frame-to-frame (nuovo metodo con incertezze sperimentali)

    std::cout << "\n--- Analisi frame-to-frame ---\n";
    std::optional<exp1::FrameByFrameResult> fbf;
    try {
      fbf =
          exp1::analyze_frame_by_frame(good_frames, bad_frames, bg_frames, opt.roi, opt.clip_sigma);
    } catch (const std::exception& e) {
      std::cerr << "  [WARNING] Analisi frame-to-frame fallita: " << e.what() << "\n";
    }

    // 3. Stacking (metodo standard) e pulizia RAM

    std::cout << "\n--- Stacking serie GOOD ---\n";
    auto good_stack = do_stack(good_frames, opt);
    good_frames.clear();
    good_frames.shrink_to_fit();

    std::cout << "\n--- Stacking serie BAD ---\n";
    auto bad_stack = do_stack(bad_frames, opt);
    bad_frames.clear();
    bad_frames.shrink_to_fit();

    std::cout << "\n--- Stacking serie BACKGROUND ---\n";
    auto bg_stack = do_stack(bg_frames, opt);
    bg_frames.clear();
    bg_frames.shrink_to_fit();

    // 4. Sottrazione fondo

    std::cout << "\n Sottrazione fondo\n";
    auto good_diff = exp1::subtract_background(good_stack, bg_stack);
    auto bad_diff  = exp1::subtract_background(bad_stack, bg_stack);

    // 5. Integrazione

    std::cout << "\n Integrazione\n";
    auto good_int = exp1::integrate(good_diff, opt.roi);
    auto bad_int  = exp1::integrate(bad_diff, opt.roi);

    std::cout << "  ROI usata: [" << good_int.roi_used.x0 << "," << good_int.roi_used.y0 << "] → ["
              << good_int.roi_used.x1 << "," << good_int.roi_used.y1 << "]  (" << good_int.n_pixels
              << " pixel)\n";

    // 6. Confronto

    auto comparison = exp1::compare(good_int, bad_int);

    // 7. Confronto simulazione (opzionale) — solo events.root da optimization_main

    std::optional<exp1::SimEff1> maybe_sim;
    if (!opt.lens_path.empty() && opt.x1_good != 0.0 && opt.x1_bad != 0.0) {

      // Parametri sorgente per la correzione geometrica
      double src_xh = 30.0; // source_thickness/2 [mm]
      double src_zh = 5.0;  // source_width/2     [mm]
      {
        std::ifstream fcfg(opt.config_path);
        if (fcfg) {
          nlohmann::json j;
          fcfg >> j;
          src_xh = j.value("source_thickness", 60.0) / 2.0;
          src_zh = j.value("source_width",     10.0) / 2.0;
        }
      }

      TFile ev_file(opt.lens_path.c_str(), "READ");
      if (ev_file.IsZombie()) {
        std::cerr << "[WARNING] Impossibile aprire il file ROOT: " << opt.lens_path << "\n";
      } else {
        TTree* tc  = static_cast<TTree*>(ev_file.Get("Configurations"));
        TTree* tef = static_cast<TTree*>(ev_file.Get("Efficiency"));

        if (!tc || !tef) {
          std::cerr << "[WARNING] Tree 'Configurations' o 'Efficiency' non trovato in "
                    << opt.lens_path << "\n";
        } else {
          // Leggi Configurations → mappa config_id → (x1, x2)
          struct CfgPos { double x1, x2; };
          std::map<int, CfgPos> cfg_map;
          {
            double x1_val, x2_val;
            int cfg_id_val;
            tc->SetBranchAddress("config_id", &cfg_id_val);
            tc->SetBranchAddress("x1", &x1_val);
            tc->SetBranchAddress("x2", &x2_val);
            for (Long64_t i = 0; i < tc->GetEntries(); ++i) {
              tc->GetEntry(i);
              cfg_map[cfg_id_val] = {x1_val, x2_val};
            }
          }

          // Leggi Efficiency: una riga per configurazione
          std::map<int, double> eff_map;
          std::map<int, double> total_ph;
          {
            int eff_cfg_id, n_photons_branch;
            double n_hits_val;
            tef->SetBranchAddress("config_id",  &eff_cfg_id);
            tef->SetBranchAddress("n_photons",  &n_photons_branch);
            tef->SetBranchAddress("n_hits",     &n_hits_val);
            for (Long64_t i = 0; i < tef->GetEntries(); ++i) {
              tef->GetEntry(i);
              if (n_photons_branch > 0) {
                eff_map[eff_cfg_id]  = n_hits_val / n_photons_branch;
                total_ph[eff_cfg_id] = static_cast<double>(n_photons_branch);
              }
            }
          }

          int id_good = -1, id_bad = -1;
          for (const auto& [id, pos] : cfg_map) {
            if (std::abs(pos.x1 - opt.x1_good) < 1e-3 && std::abs(pos.x2 - opt.x2_good) < 1e-3)
              id_good = id;
            if (std::abs(pos.x1 - opt.x1_bad) < 1e-3 && std::abs(pos.x2 - opt.x2_bad) < 1e-3)
              id_bad = id;
          }

          if (id_good >= 0 && id_bad >= 0 && eff_map.count(id_good) && eff_map.count(id_bad)) {
            exp1::SimEff1 s;
            s.eta_good = eff_map[id_good];
            s.eta_bad  = eff_map[id_bad];
            if (s.eta_bad > 0.0) {
              s.ratio_sim = s.eta_good / s.eta_bad;

              double n_g  = total_ph[id_good];
              double n_b  = total_ph[id_bad];
              double se_g = (n_g > 0.0) ? std::sqrt(s.eta_good * (1.0 - s.eta_good) / n_g) : 0.0;
              double se_b = (n_b > 0.0) ? std::sqrt(s.eta_bad  * (1.0 - s.eta_bad)  / n_b) : 0.0;
              s.sigma_ratio_sim = s.ratio_sim
                                  * std::sqrt((se_g / s.eta_good) * (se_g / s.eta_good)
                                              + (se_b / s.eta_bad) * (se_b / s.eta_bad));

              // Fattore correttivo geometrico sorgente rettangolare → sferica a origine
              double f_good = geom_correction(opt.x1_good, src_xh, src_zh);
              double f_bad  = geom_correction(opt.x1_bad,  src_xh, src_zh);
              s.c_geom               = f_good / f_bad;
              s.ratio_sim_corr       = s.ratio_sim       * s.c_geom;
              s.sigma_ratio_sim_corr = s.sigma_ratio_sim * s.c_geom;

              maybe_sim = s;
            }
          } else {
            std::cerr << "[WARNING] Configurazioni (x1=" << opt.x1_good << ",x2=" << opt.x2_good
                      << ") o (x1=" << opt.x1_bad << ",x2=" << opt.x2_bad
                      << ") non trovate nel file ROOT\n";
          }
        }
        ev_file.Close();
      }
    }

    // 8. Output ROOT

    std::cout << "\n Produzione output\n";
    exp1::AnalysisConfig acfg;
    acfg.output_dir       = opt.output;
    acfg.roi              = opt.roi;
    acfg.save_png         = opt.save_png;
    acfg.save_root        = opt.save_root;
    acfg.z_min_percentile = opt.z_min;
    acfg.z_max_percentile = opt.z_max;

    exp1::produce_output(good_diff, bad_diff, good_stack, bad_stack, bg_stack, comparison, fbf,
                         maybe_sim, acfg);

  } catch (const std::exception& e) {
    std::cerr << "\n[ERRORE] " << e.what() << "\n";
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "\n[ERRORE] Eccezione sconosciuta\n";
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}