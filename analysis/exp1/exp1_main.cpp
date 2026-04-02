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

#include <filesystem>
#include <iostream>
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

// main

int main(int argc, char** argv) {
  CliOptions opt = parse_args(argc, argv);

  std::cout << "═══════════════════════════════════════════════════════════\n";
  std::cout << "  exp1 — Analisi immagini FITS 16 bit\n";
  std::cout << "═══════════════════════════════════════════════════════════\n";
  std::cout << "  Dati: " << opt.data_dir.string() << "\n";
  std::cout << "  Output: " << opt.output.string() << "\n";

  try {
    // 1. Elaborazione serie per serie (stacking + clear) per risparmiare RAM

    std::cout << "\n--- Elaborazione serie GOOD ---\n";
    auto good_frames = exp1::read_fits_stack(opt.data_dir / opt.good_dir, opt.max_frames);
    auto good_stack  = do_stack(good_frames, opt);
    good_frames.clear();
    good_frames.shrink_to_fit();

    std::cout << "\n--- Elaborazione serie BAD ---\n";
    auto bad_frames = exp1::read_fits_stack(opt.data_dir / opt.bad_dir, opt.max_frames);
    auto bad_stack  = do_stack(bad_frames, opt);
    bad_frames.clear();
    bad_frames.shrink_to_fit();

    std::cout << "\n--- Elaborazione serie BACKGROUND ---\n";
    auto bg_frames = exp1::read_fits_stack(opt.data_dir / opt.bg_dir, opt.max_frames);
    auto bg_stack  = do_stack(bg_frames, opt);
    bg_frames.clear();
    bg_frames.shrink_to_fit();

    // 2. Sottrazione fondo

    std::cout << "\n Sottrazione fondo\n";
    auto good_diff = exp1::subtract_background(good_stack, bg_stack);
    auto bad_diff  = exp1::subtract_background(bad_stack, bg_stack);

    // 3. Integrazione

    std::cout << "\n Integrazione\n";
    auto good_int = exp1::integrate(good_diff, opt.roi);
    auto bad_int  = exp1::integrate(bad_diff, opt.roi);

    std::cout << "  ROI usata: [" << good_int.roi_used.x0 << "," << good_int.roi_used.y0 << "] → ["
              << good_int.roi_used.x1 << "," << good_int.roi_used.y1 << "]  (" << good_int.n_pixels
              << " pixel)\n";

    // 4. Confronto

    auto comparison = exp1::compare(good_int, bad_int);

    // 5. Output ROOT

    std::cout << "\n Produzione output\n";
    exp1::AnalysisConfig acfg;
    acfg.output_dir       = opt.output;
    acfg.roi              = opt.roi;
    acfg.save_png         = opt.save_png;
    acfg.save_root        = opt.save_root;
    acfg.z_min_percentile = opt.z_min;
    acfg.z_max_percentile = opt.z_max;

    exp1::produce_output(good_diff, bad_diff, good_stack, bad_stack, bg_stack, comparison, acfg);

  } catch (const std::exception& e) {
    std::cerr << "\n[ERRORE] " << e.what() << "\n";
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "\n[ERRORE] Eccezione sconosciuta\n";
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}