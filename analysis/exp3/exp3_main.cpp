/* Copyright 2026 Giulio Mesini — Apache License 2.0 */

#include "exp3_analysis.hpp"
#include "exp3b_analysis.hpp"
#include "exp3_config.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------
struct CliOptions {
  std::string mode;           // calib | measure | exp3b | plots
  fs::path config_path  = "config/exp3/exp3_config.json";
  fs::path data_dir     = "data/exp3";
  fs::path data_dir_b   = "data/exp3b";
  fs::path output_dir   = "output/exp3";
  fs::path output_dir_b = "output/exp3b";
  std::string lens_config = "all";   // good | bad | all
  bool no_png  = false;
  bool verbose = false;
  bool help    = false;
};

static void print_usage(const char* prog) {
  std::cout
      << "Uso: " << prog << " --mode <modalita> [opzioni]\n"
      << "\nModalita:\n"
      << "  calib     Calibrazione via linea di riferimento\n"
      << "  measure   Misura Q sulle 3 tracce parallele (usa calib gia' eseguita)\n"
      << "  exp3b     Misura M(r) tramite colonna di dot\n"
      << "  plots     Produce tutti i grafici ROOT\n"
      << "\nOpzioni:\n"
      << "  --config <path>         Config JSON [default: config/exp3/exp3_config.json]\n"
      << "  --data-dir <path>       Root dati Exp3 [default: data/exp3/]\n"
      << "  --data-dir-b <path>     Root dati Exp3b [default: data/exp3b/]\n"
      << "  --output <path>         Root output Exp3 [default: output/exp3/]\n"
      << "  --output-b <path>       Root output Exp3b [default: output/exp3b/]\n"
      << "  --lens-config <s>       good|bad|all [default: all]\n"
      << "  --no-png                Disabilita salvataggio PNG\n"
      << "  --verbose               Output dettagliato\n"
      << "  --help                  Mostra questo messaggio\n"
      << "\nStruttura dati attesa:\n"
      << "  data/exp3/calib/d{dist}/signal|background\n"
      << "  data/exp3/{good|bad}/d{dist}/signal|background\n"
      << "  data/exp3b/{good|bad}/d{dist}/signal|background\n"
      << "\nOutput prodotti (calib):\n"
      << "  output/exp3/calib/scale_d{dist}mm.json\n"
      << "Output prodotti (measure):\n"
      << "  output/exp3/{good|bad}/q_exp_map.tsv\n"
      << "Output prodotti (exp3b):\n"
      << "  output/exp3b/{good|bad}/M_summary.tsv\n"
      << "Output prodotti (plots):\n"
      << "  output/exp3/q_vs_r_d{dist}.png\n"
      << "  output/exp3/q_vs_dax_r{idx}.png\n"
      << "  output/exp3b/M_profile_d{dist}_{config}.png\n"
      << "  output/exp3b/M_vs_dax.png\n"
      << "  output/exp3b/M_nonlinearity_d{dist}_{config}.png\n";
}

static CliOptions parse_args(int argc, char** argv) {
  CliOptions opt;
  if (argc < 2) { print_usage(argv[0]); std::exit(1); }

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) {
        std::cerr << "Argomento mancante dopo " << arg << "\n";
        std::exit(1);
      }
      return argv[++i];
    };

    if      (arg == "--help" || arg == "-h") { opt.help = true; }
    else if (arg == "--mode")       { opt.mode        = next(); }
    else if (arg == "--config")     { opt.config_path = next(); }
    else if (arg == "--data-dir")   { opt.data_dir    = next(); }
    else if (arg == "--data-dir-b") { opt.data_dir_b  = next(); }
    else if (arg == "--output")     { opt.output_dir  = next(); }
    else if (arg == "--output-b")   { opt.output_dir_b = next(); }
    else if (arg == "--lens-config") { opt.lens_config = next(); }
    else if (arg == "--no-png")     { opt.no_png  = true; }
    else if (arg == "--verbose")    { opt.verbose = true; }
    else {
      std::cerr << "Opzione sconosciuta: " << arg << "\n";
      print_usage(argv[0]);
      std::exit(1);
    }
  }
  return opt;
}

// ---------------------------------------------------------------------------
// Helper: formatta "d150" da distanza nominale
// ---------------------------------------------------------------------------
static std::string dist_str(double d) {
  std::ostringstream s;
  s << "d" << std::fixed << std::setprecision(0) << d;
  return s.str();
}

static std::vector<std::string> active_configs(const std::string& lens_config) {
  if (lens_config == "good") return {"good"};
  if (lens_config == "bad")  return {"bad"};
  return {"good", "bad"};
}

// ---------------------------------------------------------------------------
// Fase calib
// ---------------------------------------------------------------------------
static void run_calib(const CliOptions& opt, const riptide::exp3::Exp3Config& cfg) {
  namespace exp3 = riptide::exp3;
  fs::path calib_data = opt.data_dir / "calib";
  fs::path calib_out  = opt.output_dir / "calib";
  fs::create_directories(calib_out);

  size_t n = cfg.axial_distances_nominal_mm.size();
  for (size_t di = 0; di < n; ++di) {
    double d_nom = cfg.axial_distances_nominal_mm[di];
    fs::path sig_dir = calib_data / dist_str(d_nom) / "signal";
    fs::path bg_dir  = calib_data / dist_str(d_nom) / "background";

    if (!fs::exists(sig_dir)) {
      std::cerr << "[exp3 calib] WARNING: non trovato: " << sig_dir << "\n";
      continue;
    }
    try {
      auto calib = exp3::calibrate_line(sig_dir, bg_dir, cfg);
      calib.d_ax_mm = cfg.axial_distances_measured_mm[di];

      std::ostringstream fname;
      fname << "scale_d" << std::fixed << std::setprecision(0) << d_nom << "mm.json";
      exp3::save_line_calib(calib, calib_out / fname.str());

      if (opt.verbose)
        std::cout << "[exp3 calib] d=" << d_nom
                  << " mm  scale=" << std::fixed << std::setprecision(6)
                  << calib.scale_mm_per_sens_px
                  << " mm/px  L_sens=" << std::fixed << std::setprecision(1)
                  << calib.L_px_sensor << " px\n";
    } catch (const std::exception& e) {
      std::cerr << "[exp3 calib] ERROR d=" << d_nom << ": " << e.what() << "\n";
    }
  }
}

// ---------------------------------------------------------------------------
// Fase measure (Exp3)
// ---------------------------------------------------------------------------
static void run_measure(const CliOptions& opt, const riptide::exp3::Exp3Config& cfg) {
  namespace exp3 = riptide::exp3;
  fs::path calib_out = opt.output_dir / "calib";

  for (const auto& lc : active_configs(opt.lens_config)) {
    std::vector<exp3::QResult> all_results;
    size_t n = cfg.axial_distances_nominal_mm.size();

    for (size_t di = 0; di < n; ++di) {
      double d_nom  = cfg.axial_distances_nominal_mm[di];
      double d_meas = cfg.axial_distances_measured_mm[di];

      // Carica calibrazione
      std::ostringstream fname;
      fname << "scale_d" << std::fixed << std::setprecision(0) << d_nom << "mm.json";
      fs::path calib_path = calib_out / fname.str();

      exp3::LineCalib calib;
      try {
        calib = exp3::load_line_calib(calib_path);
      } catch (const std::exception&) {
        // Fallback: scala display usata direttamente da detect_lines_auto
        if (opt.verbose)
          std::cerr << "[exp3 measure] INFO: calib assente per d=" << d_nom
                    << ", uso scala display come fallback\n";
      }

      fs::path sig_dir = opt.data_dir / lc / dist_str(d_nom) / "signal";
      fs::path bg_dir  = opt.data_dir / lc / dist_str(d_nom) / "background";

      if (!fs::exists(sig_dir)) {
        if (opt.verbose)
          std::cerr << "[exp3 measure] WARNING: non trovato: " << sig_dir << "\n";
        continue;
      }

      try {
        auto results = exp3::run_measurement_parallel_lines(sig_dir, bg_dir, cfg, calib, d_meas);
        if (opt.verbose)
          for (const auto& r : results)
            std::cout << "[exp3 measure] " << lc << " d=" << d_nom
                      << " r_idx=" << r.r_idx
                      << " chi2=" << std::fixed << std::setprecision(3) << r.chi2_ndof
                      << (r.warning ? " [!]" : "") << "\n";
        for (auto& r : all_results) (void)r;
        all_results.insert(all_results.end(), results.begin(), results.end());
      } catch (const std::exception& e) {
        std::cerr << "[exp3 measure] ERROR " << lc << " d=" << d_nom
                  << ": " << e.what() << "\n";
      }
    }

    if (!all_results.empty()) {
      fs::path tsv_path = opt.output_dir / lc / "q_exp_map.tsv";
      try {
        exp3::write_q_results_tsv(lc, all_results, tsv_path);
        std::cout << "[exp3 measure] Salvato: " << tsv_path << "\n";
      } catch (const std::exception& e) {
        std::cerr << "[exp3 measure] WARNING TSV: " << e.what() << "\n";
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Fase exp3b
// ---------------------------------------------------------------------------
static void run_exp3b(const CliOptions& opt, const riptide::exp3::Exp3Config& cfg) {
  namespace exp3b = riptide::exp3b;
  namespace exp3  = riptide::exp3;
  fs::path calib_out = opt.output_dir / "calib";

  for (const auto& lc : active_configs(opt.lens_config)) {
    std::vector<exp3b::MResult> all_results;
    size_t n = cfg.axial_distances_nominal_mm.size();

    for (size_t di = 0; di < n; ++di) {
      double d_nom  = cfg.axial_distances_nominal_mm[di];
      double d_meas = cfg.axial_distances_measured_mm[di];

      std::ostringstream fname;
      fname << "scale_d" << std::fixed << std::setprecision(0) << d_nom << "mm.json";
      exp3::LineCalib calib;
      try {
        calib = exp3::load_line_calib(calib_out / fname.str());
      } catch (const std::exception& e) {
        std::cerr << "[exp3b] WARNING: calib non trovata per d=" << d_nom
                  << ": " << e.what() << "\n";
        continue;
      }

      fs::path sig_dir = opt.data_dir_b / lc / dist_str(d_nom) / "signal";
      fs::path bg_dir  = opt.data_dir_b / lc / dist_str(d_nom) / "background";

      if (!fs::exists(sig_dir)) {
        if (opt.verbose)
          std::cerr << "[exp3b] WARNING: non trovato: " << sig_dir << "\n";
        continue;
      }

      try {
        auto r = exp3b::run_exp3b(sig_dir, bg_dir, cfg, calib, d_meas);
        if (opt.verbose)
          std::cout << "[exp3b] " << lc << " d=" << d_nom
                    << " M=" << std::fixed << std::setprecision(4) << r.M_global
                    << " n_dots=" << r.n_dots << "\n";
        all_results.push_back(r);
      } catch (const std::exception& e) {
        std::cerr << "[exp3b] ERROR " << lc << " d=" << d_nom
                  << ": " << e.what() << "\n";
      }
    }

    if (!all_results.empty()) {
      fs::path tsv_path = opt.output_dir_b / lc / "M_summary.tsv";
      try {
        exp3b::write_M_summary_tsv(lc, all_results, tsv_path);
        std::cout << "[exp3b] Salvato: " << tsv_path << "\n";
      } catch (const std::exception& e) {
        std::cerr << "[exp3b] WARNING TSV: " << e.what() << "\n";
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Fase plots
// ---------------------------------------------------------------------------
static void run_plots(const CliOptions& opt, const riptide::exp3::Exp3Config& cfg) {
  namespace exp3  = riptide::exp3;
  namespace exp3b = riptide::exp3b;

  if (opt.no_png) return;

  // Carica i QResult da TSV
  auto load_q = [&](const std::string& lc) -> std::vector<exp3::QResult> {
    std::vector<exp3::QResult> v;
    fs::path p = opt.output_dir / lc / "q_exp_map.tsv";
    std::ifstream ifs(p);
    if (!ifs) return v;
    std::string line;
    std::getline(ifs, line); // header
    while (std::getline(ifs, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      std::string config; int warn;
      exp3::QResult r;
      if (!(ss >> config >> r.d_ax_mm >> r.r_mm >> r.r_idx
                >> r.chi2_ndof >> r.n_valid_slices >> warn)) continue;
      r.warning = (warn != 0);
      v.push_back(r);
    }
    return v;
  };

  auto good_q = load_q("good");
  auto bad_q  = load_q("bad");

  // confronto Q_exp(d_ax) vs Q_sim per good e bad
  try {
    exp3::produce_q_comparison(good_q, bad_q, cfg.q_comparison,
                               opt.output_dir / "q_comparison.png");
    std::cout << "[plots] Salvato: " << opt.output_dir / "q_comparison.png" << "\n";
  } catch (const std::exception& e) {
    std::cerr << "[plots] WARNING q_comparison: " << e.what() << "\n";
  }

  // q_vs_dax per ogni r_idx
  int n_r = static_cast<int>(cfg.radial_offsets_px.size());
  for (int ri = 0; ri < n_r; ++ri) {
    std::ostringstream fname;
    fname << "q_vs_dax_r" << ri << ".png";
    try {
      exp3::produce_q_vs_dax(good_q, bad_q, ri, opt.output_dir / fname.str());
    } catch (const std::exception& e) {
      std::cerr << "[plots] WARNING q_vs_dax r=" << ri << ": " << e.what() << "\n";
    }
  }

  // Carica MResult da TSV e produce i plot Exp3b
  auto load_m = [&](const std::string& lc) -> std::vector<exp3b::MResult> {
    std::vector<exp3b::MResult> v;
    fs::path p = opt.output_dir_b / lc / "M_summary.tsv";
    std::ifstream ifs(p);
    if (!ifs) return v;
    std::string line;
    std::getline(ifs, line);
    while (std::getline(ifs, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      std::string config;
      exp3b::MResult r;
      if (!(ss >> config >> r.d_ax_mm >> r.M_global >> r.M_rms_residual
                >> r.chi2_ndof >> r.n_dots)) continue;
      v.push_back(r);
    }
    return v;
  };

  auto good_m = load_m("good");
  auto bad_m  = load_m("bad");

  // M_vs_dax
  try {
    exp3b::produce_M_vs_dax(good_m, bad_m, opt.output_dir_b / "M_vs_dax.png");
  } catch (const std::exception& e) {
    std::cerr << "[plots] WARNING M_vs_dax: " << e.what() << "\n";
  }

  // M_profile e M_nonlinearity per ogni distanza + config
  // (i dati di r_mm e M_local non sono nel TSV di summary: questi plot
  //  richiedono di ri-eseguire run_exp3b con i FITS, oppure un TSV più ricco.
  //  Per ora vengono saltati se MResult ha n_dots < 2.)
  for (const auto& lc : active_configs(opt.lens_config)) {
    exp3::LineCalib dummy_calib;
    for (size_t di = 0; di < cfg.axial_distances_nominal_mm.size(); ++di) {
      double d_nom  = cfg.axial_distances_nominal_mm[di];
      double d_meas = cfg.axial_distances_measured_mm[di];

      std::ostringstream cfname;
      cfname << "scale_d" << std::fixed << std::setprecision(0) << d_nom << "mm.json";
      try { dummy_calib = exp3::load_line_calib(opt.output_dir / "calib" / cfname.str()); }
      catch (...) { continue; }

      fs::path sig_dir = opt.data_dir_b / lc / dist_str(d_nom) / "signal";
      fs::path bg_dir  = opt.data_dir_b / lc / dist_str(d_nom) / "background";
      if (!fs::exists(sig_dir)) continue;

      try {
        auto mr = exp3b::run_exp3b(sig_dir, bg_dir, cfg, dummy_calib, d_meas);
        if (mr.n_dots < 2) continue;

        std::ostringstream pf, nf;
        pf << "M_profile_d" << std::fixed << std::setprecision(0) << d_nom
           << "_" << lc << ".png";
        nf << "M_nonlinearity_d" << std::fixed << std::setprecision(0) << d_nom
           << "_" << lc << ".png";

        exp3b::produce_M_profile(mr, lc, opt.output_dir_b / pf.str());
        exp3b::produce_M_nonlinearity(mr, lc, opt.output_dir_b / nf.str());
      } catch (const std::exception& e) {
        std::cerr << "[plots] WARNING M_profile " << lc << " d=" << d_nom
                  << ": " << e.what() << "\n";
      }
    }
  }
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char** argv) {
  auto opt = parse_args(argc, argv);

  if (opt.help) { print_usage(argv[0]); return 0; }

  if (opt.mode.empty()) {
    std::cerr << "[exp3] Errore: specificare --mode calib|measure|exp3b|plots\n\n";
    print_usage(argv[0]);
    return 1;
  }

  riptide::exp3::Exp3Config cfg;
  try {
    cfg = riptide::exp3::load_exp3_config(opt.config_path);
  } catch (const std::exception& e) {
    std::cerr << "[exp3] Errore caricamento config: " << e.what() << "\n";
    return 1;
  }

  fs::create_directories(opt.output_dir);
  fs::create_directories(opt.output_dir_b);

  if (opt.verbose) {
    std::cout << "[exp3] Mode:   " << opt.mode << "\n"
              << "[exp3] Config: " << opt.config_path << "\n"
              << "[exp3] Data:   " << opt.data_dir << "\n"
              << "[exp3] Output: " << opt.output_dir << "\n"
              << "[exp3] Distanze assiali: "
              << cfg.axial_distances_nominal_mm.size() << "\n"
              << "[exp3] Offset radiali:   "
              << cfg.radial_offsets_px.size() << "\n";
  }

  try {
    if      (opt.mode == "calib")   { run_calib(opt, cfg); }
    else if (opt.mode == "measure") { run_measure(opt, cfg); }
    else if (opt.mode == "exp3b")   { run_exp3b(opt, cfg); }
    else if (opt.mode == "plots")   { run_plots(opt, cfg); }
    else {
      std::cerr << "[exp3] Modalita' sconosciuta: " << opt.mode << "\n";
      print_usage(argv[0]);
      return 1;
    }
  } catch (const std::exception& e) {
    std::cerr << "[exp3] Errore fatale: " << e.what() << "\n";
    return 1;
  }

  std::cout << "[exp3] Completato.\n";
  return 0;
}
