/* Copyright 2026 Giulio Mesini — Apache License 2.0 */

#include "exp3_analysis.hpp"
#include "exp3_config.hpp"
#include "homography.hpp"

#include "fits_io.hpp"
#include "stacking.hpp"

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

// ---------------------------------------------------------------------------
// Opzioni CLI
// ---------------------------------------------------------------------------
struct CliOptions {
  // Modalità (obbligatorio uno tra calibrate/analyze/full)
  bool mode_calibrate = false;
  bool mode_analyze   = false;
  bool mode_full      = false;

  fs::path config_path  = "config/exp3/exp3_config.json";
  fs::path data_dir     = "data/exp3";
  fs::path output_dir   = "output/exp3";
  std::string lens_config = "all";   // "good" | "bad" | "all"

  bool no_png    = false;
  bool no_root   = true;   // ROOT disabilitato di default (come exp2)
  bool verbose   = false;
  bool help      = false;
};

static void print_usage(const char* prog) {
  std::cout
      << "Uso: " << prog << " <modalita> [opzioni]\n"
      << "\nModalita (obbligatorio scegliere una):\n"
      << "  --calibrate       Esegue solo la calibrazione geometrica\n"
      << "  --analyze         Esegue l'analisi (calibrazione gia' eseguita)\n"
      << "  --full            Calibrazione + analisi in sequenza\n"
      << "\nOpzioni comuni:\n"
      << "  --config <path>         Config JSON [default: config/exp3/exp3_config.json]\n"
      << "  --data-dir <path>       Root cartella dati [default: data/exp3/]\n"
      << "  --output <path>         Root cartella output [default: output/exp3/]\n"
      << "  --lens-config <s>       good|bad|all [default: all]\n"
      << "  --no-png                Disabilita salvataggio PNG\n"
      << "  --no-root               Disabilita salvataggio ROOT (default: gia' disabilitato)\n"
      << "  --verbose               Output dettagliato per ogni acquisizione\n"
      << "  --help                  Mostra questo messaggio\n"
      << "\nStruttura dati attesa:\n"
      << "  data/exp3/calib/d{dist}/     FITS griglia calibrazione\n"
      << "  data/exp3/calib/background/  FITS background calibrazione\n"
      << "  data/exp3/{good|bad}/d{dist}/theta{theta}/  FITS segnale\n"
      << "  data/exp3/{good|bad}/background/            FITS background\n"
      << "\nOutput prodotti:\n"
      << "  output/exp3/calib/homography_d{dist}mm.json\n"
      << "  output/exp3/calib/calib_report.png\n"
      << "  output/exp3/{good|bad}/trace_d{dist}_theta{angle}.png\n"
      << "  output/exp3/{good|bad}/Q_profile.png\n"
      << "  output/exp3/summary.png\n"
      << "  output/exp3/results.tsv\n"
      << "  output/exp3/summary.tsv\n";
}

static CliOptions parse_args(int argc, char** argv) {
  CliOptions opt;

  if (argc < 2) {
    print_usage(argv[0]);
    std::exit(1);
  }

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];

    auto next = [&]() -> std::string {
      if (i + 1 >= argc) {
        std::cerr << "Argomento mancante dopo " << arg << "\n";
        std::exit(1);
      }
      return argv[++i];
    };

    if (arg == "--help" || arg == "-h") {
      opt.help = true;
    } else if (arg == "--calibrate") {
      opt.mode_calibrate = true;
    } else if (arg == "--analyze") {
      opt.mode_analyze = true;
    } else if (arg == "--full") {
      opt.mode_full = true;
    } else if (arg == "--config") {
      opt.config_path = next();
    } else if (arg == "--data-dir") {
      opt.data_dir = next();
    } else if (arg == "--output") {
      opt.output_dir = next();
    } else if (arg == "--lens-config") {
      opt.lens_config = next();
    } else if (arg == "--no-png") {
      opt.no_png = true;
    } else if (arg == "--no-root") {
      opt.no_root = true;
    } else if (arg == "--verbose") {
      opt.verbose = true;
    } else {
      std::cerr << "Opzione sconosciuta: " << arg << "\n";
      print_usage(argv[0]);
      std::exit(1);
    }
  }

  return opt;
}

// ---------------------------------------------------------------------------
// Fase di calibrazione
// ---------------------------------------------------------------------------
static std::vector<riptide::exp3::Homography>
run_calibration(const CliOptions& opt,
                const riptide::exp3::Exp3Config& cfg,
                std::vector<std::vector<riptide::exp3::CalibPoint>>* pts_out = nullptr) {
  namespace exp3 = riptide::exp3;

  fs::path calib_data = opt.data_dir / "calib";
  fs::path calib_out  = opt.output_dir / "calib";
  fs::create_directories(calib_out);

  std::vector<exp3::Homography> homographies;
  std::vector<double> axial_dists;
  std::vector<std::vector<exp3::CalibPoint>> all_pts;

  size_t n_dist = cfg.axial_distances_nominal_mm.size();
  for (size_t di = 0; di < n_dist; ++di) {
    double d_nom = cfg.axial_distances_nominal_mm[di];

    std::ostringstream d_str;
    d_str << "d" << std::fixed << std::setprecision(0) << d_nom;

    fs::path sig_dir = calib_data / d_str.str() / "signal";
    fs::path bg_dir  = calib_data / d_str.str() / "background";

    if (!fs::exists(sig_dir)) {
      std::cerr << "[exp3] WARNING: directory calibrazione non trovata: "
                << sig_dir << "\n";
      continue;
    }

    try {
      auto H = exp3::calibrate_distance(sig_dir, bg_dir, cfg);

      // Salva omografia
      std::ostringstream hfile;
      hfile << "homography_d" << std::fixed << std::setprecision(0) << d_nom << "mm.json";
      exp3::save_homography(H, calib_out / hfile.str());

      if (opt.verbose)
        std::cout << "[exp3] d=" << d_nom << " mm: n_pts=" << H.n_points
                  << " RMS=" << std::fixed << std::setprecision(3) << H.rms_residual
                  << " px\n";

      homographies.push_back(H);
      axial_dists.push_back(d_nom);

      // Recupera i punti per il report (ricalcola solo la detection)
      std::vector<exp3::CalibPoint> pts;
      try {
        auto sig_frames = riptide::fits::read_fits_stack(sig_dir);
        auto bg_frames  = fs::exists(bg_dir) ? riptide::fits::read_fits_stack(bg_dir)
                                             : std::vector<riptide::fits::FitsFrame>{};
        riptide::stack::StackConfig sc;
        sc.n_sigma = cfg.stack_n_sigma; sc.n_iter = cfg.stack_n_iter;
        sc.method  = riptide::stack::StackMethod::SigmaClip;
        auto sig_st = riptide::stack::sigma_clip_stack(sig_frames, sc);
        riptide::stack::StackedImage bg_st;
        if (!bg_frames.empty())
          bg_st = riptide::stack::mean_stack(bg_frames);
        else {
          bg_st.width = sig_st.width; bg_st.height = sig_st.height;
          bg_st.mean.assign(static_cast<size_t>(sig_st.width * sig_st.height), 0.0);
        }
        std::vector<double> diff(sig_st.npixels());
        for (size_t i = 0; i < diff.size(); ++i)
          diff[i] = sig_st.mean[i] - bg_st.mean[i];
        pts = exp3::detect_calibration_dots(diff, sig_st.width, sig_st.height,
                                             cfg.calibration_grid_step_px,
                                             cfg.display.width_px, cfg.display.height_px);
      } catch (...) {}
      all_pts.push_back(pts);

    } catch (const std::exception& e) {
      std::cerr << "[exp3] ERROR calibrazione d=" << d_nom << ": " << e.what() << "\n";
      continue;
    }
  }

  // Produce calib_report.png
  if (!homographies.empty() && !opt.no_png) {
    try {
      exp3::produce_calibration_report(homographies, axial_dists, all_pts,
                                        calib_out / "calib_report.png");
      std::cout << "[exp3] Salvato: " << (calib_out / "calib_report.png") << "\n";
    } catch (const std::exception& e) {
      std::cerr << "[exp3] WARNING calib_report: " << e.what() << "\n";
    }
  }

  if (pts_out) *pts_out = all_pts;
  return homographies;
}

// ---------------------------------------------------------------------------
// Fase di analisi
// ---------------------------------------------------------------------------
static void run_analysis(const CliOptions& opt,
                          const riptide::exp3::Exp3Config& cfg) {
  namespace exp3 = riptide::exp3;

  fs::path calib_out = opt.output_dir / "calib";

  exp3::OutputConfig out_cfg;
  out_cfg.output_dir  = opt.output_dir;
  out_cfg.save_png    = !opt.no_png;
  out_cfg.save_root   = !opt.no_root;
  out_cfg.verbose     = opt.verbose;

  // Filtra configurazioni lenti
  std::vector<exp3::LensConfig> lens_to_run;
  for (const auto& lc : cfg.lens_configs) {
    if (opt.lens_config == "all" || opt.lens_config == lc.label)
      lens_to_run.push_back(lc);
  }

  if (lens_to_run.empty()) {
    std::cerr << "[exp3] WARNING: nessuna configurazione lenti da analizzare.\n";
    return;
  }

  std::vector<exp3::ConfigResult> all_results;
  for (const auto& lc : lens_to_run) {
    std::cout << "[exp3] Analisi configurazione: " << lc.label << "\n";
    try {
      auto result = exp3::analyze_config(opt.data_dir, calib_out, lc, cfg, out_cfg);
      std::cout << "[exp3]   Q_exp_global = " << std::fixed << std::setprecision(4)
                << result.Q_exp_global
                << "  Q_sim = " << result.Q_sim
                << "  R = " << result.R << "\n";
      all_results.push_back(result);
    } catch (const std::exception& e) {
      std::cerr << "[exp3] ERROR analisi " << lc.label << ": " << e.what() << "\n";
    }
  }

  if (all_results.empty()) return;

  // Summary PNG
  if (!opt.no_png && all_results.size() > 1) {
    try {
      exp3::produce_results_summary(all_results, opt.output_dir / "summary.png");
      std::cout << "[exp3] Salvato: " << (opt.output_dir / "summary.png") << "\n";
    } catch (const std::exception& e) {
      std::cerr << "[exp3] WARNING summary: " << e.what() << "\n";
    }
  }

  // TSV
  try {
    exp3::write_results_tsv(all_results, opt.output_dir / "results.tsv");
    std::cout << "[exp3] Salvato: " << (opt.output_dir / "results.tsv") << "\n";
  } catch (const std::exception& e) {
    std::cerr << "[exp3] WARNING results.tsv: " << e.what() << "\n";
  }

  try {
    exp3::write_summary_tsv(all_results, opt.output_dir / "summary.tsv");
    std::cout << "[exp3] Salvato: " << (opt.output_dir / "summary.tsv") << "\n";
  } catch (const std::exception& e) {
    std::cerr << "[exp3] WARNING summary.tsv: " << e.what() << "\n";
  }
}

// ---------------------------------------------------------------------------
int main(int argc, char** argv) {
  auto opt = parse_args(argc, argv);

  if (opt.help) {
    print_usage(argv[0]);
    return 0;
  }

  // Verifica che almeno una modalità sia selezionata
  if (!opt.mode_calibrate && !opt.mode_analyze && !opt.mode_full) {
    std::cerr << "[exp3] Errore: specificare --calibrate, --analyze o --full\n\n";
    print_usage(argv[0]);
    return 1;
  }

  // Carica configurazione
  riptide::exp3::Exp3Config cfg;
  try {
    cfg = riptide::exp3::load_exp3_config(opt.config_path);
  } catch (const std::exception& e) {
    std::cerr << "[exp3] Errore caricamento config: " << e.what() << "\n";
    return 1;
  }

  fs::create_directories(opt.output_dir);

  if (opt.verbose) {
    std::cout << "[exp3] Config: " << opt.config_path << "\n"
              << "[exp3] Data:   " << opt.data_dir << "\n"
              << "[exp3] Output: " << opt.output_dir << "\n"
              << "[exp3] Distanze assiali: " << cfg.axial_distances_nominal_mm.size() << "\n"
              << "[exp3] Orientazioni:     " << cfg.orientations_deg.size() << "\n"
              << "[exp3] Lenti:            " << cfg.lens_configs.size() << "\n";
  }

  try {
    if (opt.mode_calibrate || opt.mode_full) {
      std::cout << "[exp3] === Calibrazione ===\n";
      run_calibration(opt, cfg);
    }

    if (opt.mode_analyze || opt.mode_full) {
      std::cout << "[exp3] === Analisi ===\n";
      run_analysis(opt, cfg);
    }
  } catch (const std::exception& e) {
    std::cerr << "[exp3] Errore fatale: " << e.what() << "\n";
    return 1;
  }

  std::cout << "[exp3] Completato.\n";
  return 0;
}
