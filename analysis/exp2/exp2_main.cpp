/* Copyright 2026 Giulio Mesini — Apache License 2.0 */

#include "exp2_analysis.hpp"

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

struct CliOptions {
  std::filesystem::path data_dir = "analysis/exp2/data";
  std::string good_focus_dir     = "good_focus";
  std::string good_nofocus_dir   = "good_nofocus";
  std::string bad_focus_dir      = "bad_focus";
  std::string bad_nofocus_dir    = "bad_nofocus";
  std::string bg_dir             = "background";
  std::string bg_alt_dir         = "background1";

  std::filesystem::path output = "output/exp2";

  size_t max_frames  = 0;
  std::string method = "sigma";
  double clip_sigma  = 3.0;
  int clip_iter      = 3;

  int slice_width           = 5;
  int slice_step            = 3;
  double snr_min            = 5.0;
  double center_err_floor   = 0.2;
  double center_err_scale   = 1.0;
  double sigma_err_floor    = 0.2;
  double sigma_err_scale    = 1.0;
  bool enable_trace_trim    = true;
  int trace_trim_pad_slices = 10;
  int trace_trim_min_slices = 50;
  double min_aspect_ratio          = 1.3;
  double trace_trim_max_center_err = 5.0;

  double x1_good       = 0.0;
  double x2_good       = 0.0;
  double xdet_good     = 0.0;
  double xdet_good_opt = 0.0;

  double x1_bad       = 0.0;
  double x2_bad       = 0.0;
  double xdet_bad     = 0.0;
  double xdet_bad_opt = 0.0;

  double z_min = 0.005;
  double z_max = 0.995;

  bool save_png  = true;
  bool save_root = true;
};

static void print_usage(const char* prog) {
  std::cout << "Uso: " << prog << " [--data-dir <path>] [--ssd <mount>]\n"
            << "          [--good-focus <subdir>] [--good-nofocus <subdir>]\n"
            << "          [--bad-focus <subdir>]  [--bad-nofocus <subdir>]\n"
            << "          [--bg <subdir>] [--bg-alt <subdir>] [--output <path>]\n"
            << "          [--frames N] [--method sigma|mean|median]\n"
            << "          [--sigma N] [--iter N]\n"
            << "          [--slice-width N] [--slice-step N] [--snr-min F]\n"
            << "          [--center-err-floor F] [--center-err-scale F]\n"
            << "          [--sigma-err-floor F]  [--sigma-err-scale F]\n"
            << "          [--no-trim]\n"
            << "          [--trim-pad-slices N] [--trim-min-slices N]\n"
            << "          [--min-aspect-ratio F] [--trim-max-center-err F]\n"
            << "          [--x1-good F] [--x2-good F] [--xdet-good F]\n"
            << "          [--x1-bad F]  [--x2-bad F]  [--xdet-bad F]\n"
            << "          [--xdet-good-opt F] [--xdet-bad-opt F]\n"
            << "          [--z-min F] [--z-max F]\n"
            << "          [--no-root] [--no-png] [--help]\n";
}

static CliOptions parse_args(int argc, char** argv) {
  CliOptions opt;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];

    auto next = [&]() -> std::string {
      if (i + 1 >= argc) {
        std::cerr << "Argomento mancante dopo " << arg << "\n";
        std::exit(1);
      }
      return argv[++i];
    };
    auto next_int = [&]() { return std::stoi(next()); };
    auto next_szt = [&]() { return static_cast<size_t>(std::stoul(next())); };
    auto next_dbl = [&]() { return std::stod(next()); };

    if (arg == "--help" || arg == "-h") {
      print_usage(argv[0]);
      std::exit(0);
    } else if (arg == "--data-dir")
      opt.data_dir = next();
    else if (arg == "--ssd")
      opt.data_dir = std::filesystem::path(next()) / "exp2";
    else if (arg == "--good-focus")
      opt.good_focus_dir = next();
    else if (arg == "--good-nofocus")
      opt.good_nofocus_dir = next();
    else if (arg == "--bad-focus")
      opt.bad_focus_dir = next();
    else if (arg == "--bad-nofocus")
      opt.bad_nofocus_dir = next();
    else if (arg == "--bg")
      opt.bg_dir = next();
    else if (arg == "--bg-alt")
      opt.bg_alt_dir = next();
    else if (arg == "--output")
      opt.output = next();
    else if (arg == "--frames")
      opt.max_frames = next_szt();
    else if (arg == "--method")
      opt.method = next();
    else if (arg == "--sigma")
      opt.clip_sigma = next_dbl();
    else if (arg == "--iter")
      opt.clip_iter = next_int();
    else if (arg == "--slice-width")
      opt.slice_width = next_int();
    else if (arg == "--slice-step")
      opt.slice_step = next_int();
    else if (arg == "--snr-min")
      opt.snr_min = next_dbl();
    else if (arg == "--center-err-floor")
      opt.center_err_floor = next_dbl();
    else if (arg == "--center-err-scale")
      opt.center_err_scale = next_dbl();
    else if (arg == "--sigma-err-floor")
      opt.sigma_err_floor = next_dbl();
    else if (arg == "--sigma-err-scale")
      opt.sigma_err_scale = next_dbl();
    else if (arg == "--no-trim")
      opt.enable_trace_trim = false;
    else if (arg == "--trim-pad-slices")
      opt.trace_trim_pad_slices = next_int();
    else if (arg == "--trim-min-slices")
      opt.trace_trim_min_slices = next_int();
    else if (arg == "--min-aspect-ratio")
      opt.min_aspect_ratio = next_dbl();
    else if (arg == "--trim-max-center-err")
      opt.trace_trim_max_center_err = next_dbl();
    else if (arg == "--x1-good")
      opt.x1_good = next_dbl();
    else if (arg == "--x2-good")
      opt.x2_good = next_dbl();
    else if (arg == "--xdet-good")
      opt.xdet_good = next_dbl();
    else if (arg == "--x1-bad")
      opt.x1_bad = next_dbl();
    else if (arg == "--x2-bad")
      opt.x2_bad = next_dbl();
    else if (arg == "--xdet-bad")
      opt.xdet_bad = next_dbl();
    else if (arg == "--xdet-good-opt")
      opt.xdet_good_opt = next_dbl();
    else if (arg == "--xdet-bad-opt")
      opt.xdet_bad_opt = next_dbl();
    else if (arg == "--z-min")
      opt.z_min = next_dbl();
    else if (arg == "--z-max")
      opt.z_max = next_dbl();
    else if (arg == "--no-root")
      opt.save_root = false;
    else if (arg == "--no-png")
      opt.save_png = false;
    else {
      std::cerr << "Opzione sconosciuta: " << arg << "\n";
      print_usage(argv[0]);
      std::exit(1);
    }
  }

  return opt;
}

static riptide::stack::StackMethod parse_method(const std::string& s) {
  if (s == "mean")
    return riptide::stack::StackMethod::Mean;
  if (s == "median")
    return riptide::stack::StackMethod::Median;
  return riptide::stack::StackMethod::SigmaClip;
}

static std::string method_str(riptide::stack::StackMethod m) {
  switch (m) {
  case riptide::stack::StackMethod::SigmaClip:
    return "sigma";
  case riptide::stack::StackMethod::Mean:
    return "mean";
  case riptide::stack::StackMethod::Median:
    return "median";
  default:
    return "sigma";
  }
}

int main(int argc, char** argv) {
  const CliOptions opt = parse_args(argc, argv);

  std::cout << "═══════════════════════════════════════════════════════════\n";
  std::cout << "  exp2 — Analisi tracce laser su immagini FITS 16 bit\n";
  std::cout << "═══════════════════════════════════════════════════════════\n";
  std::cout << "  Dati:   " << opt.data_dir.string() << "\n";
  std::cout << "  Output: " << opt.output.string() << "\n";

  try {
    riptide::stack::StackConfig stack_cfg;
    stack_cfg.method     = parse_method(opt.method);
    stack_cfg.max_frames = opt.max_frames;
    stack_cfg.n_sigma    = opt.clip_sigma;
    stack_cfg.n_iter     = opt.clip_iter;

    riptide::exp2::TraceConfig trace_cfg;
    trace_cfg.slice_width           = opt.slice_width;
    trace_cfg.slice_step            = opt.slice_step;
    trace_cfg.min_snr               = opt.snr_min;
    trace_cfg.center_err_floor      = opt.center_err_floor;
    trace_cfg.center_err_scale      = opt.center_err_scale;
    trace_cfg.sigma_err_floor       = opt.sigma_err_floor;
    trace_cfg.sigma_err_scale       = opt.sigma_err_scale;
    trace_cfg.enable_trace_trim     = opt.enable_trace_trim;
    trace_cfg.trace_trim_pad_slices = opt.trace_trim_pad_slices;
    trace_cfg.trace_trim_min_slices = opt.trace_trim_min_slices;
    trace_cfg.min_aspect_ratio           = opt.min_aspect_ratio;
    trace_cfg.trace_trim_max_center_err  = opt.trace_trim_max_center_err;

    riptide::exp2::OutputConfig out_cfg;
    out_cfg.output_dir       = opt.output;
    out_cfg.save_png         = opt.save_png;
    out_cfg.save_root        = opt.save_root;
    out_cfg.z_min_percentile = opt.z_min;
    out_cfg.z_max_percentile = opt.z_max;

    const auto good_focus_dir   = opt.data_dir / opt.good_focus_dir;
    const auto good_nofocus_dir = opt.data_dir / opt.good_nofocus_dir;
    const auto bad_focus_dir    = opt.data_dir / opt.bad_focus_dir;
    const auto bad_nofocus_dir  = opt.data_dir / opt.bad_nofocus_dir;
    const auto bg_dir           = opt.data_dir / opt.bg_dir;
    const auto bg_alt_dir       = opt.data_dir / opt.bg_alt_dir;

    riptide::exp2::OpticsParams opt_good;
    opt_good.x1            = opt.x1_good;
    opt_good.x2            = opt.x2_good;
    opt_good.x_det         = opt.xdet_good;
    opt_good.x_det_optimal = opt.xdet_good_opt;

    riptide::exp2::OpticsParams opt_bad;
    opt_bad.x1            = opt.x1_bad;
    opt_bad.x2            = opt.x2_bad;
    opt_bad.x_det         = opt.xdet_bad;
    opt_bad.x_det_optimal = opt.xdet_bad_opt;

    if (opt_good.x_det == 0.0 && opt_good.x_det_optimal > 0.0) {
      std::cerr << "[WARNING] xdet_good non fornito: assumo xdet_good = xdet_good_opt ("
                << opt_good.x_det_optimal << " mm)\n";
      opt_good.x_det = opt_good.x_det_optimal;
    }
    if (opt_bad.x_det == 0.0 && opt_bad.x_det_optimal > 0.0) {
      std::cerr << "[WARNING] xdet_bad non fornito: assumo xdet_bad = xdet_bad_opt ("
                << opt_bad.x_det_optimal << " mm)\n";
      opt_bad.x_det = opt_bad.x_det_optimal;
    }

    riptide::exp2::OpticsParams opt_good_nf = opt_good;
    opt_good_nf.x_det                       = 0.0;
    opt_good_nf.x_det_optimal               = 0.0;

    riptide::exp2::OpticsParams opt_bad_nf = opt_bad;
    opt_bad_nf.x_det                       = 0.0;
    opt_bad_nf.x_det_optimal               = 0.0;

    std::vector<riptide::exp2::ConfigResult> results;
    results.reserve(4);

    std::cout << "\n--- Configurazione stacking: " << method_str(stack_cfg.method) << " ---\n";

    auto analyze_with_bg_fallback =
        [&](const std::filesystem::path& signal_dir, riptide::exp2::ConfigLabel label,
            const riptide::exp2::OpticsParams& optics) -> riptide::exp2::ConfigResult {
      try {
        return riptide::exp2::analyze_config(signal_dir, bg_dir, label, optics, stack_cfg,
                                             trace_cfg);
      } catch (const std::exception& e) {
        const std::string msg = e.what();
        const bool dim_mismatch =
            (msg.find("dimensioni diverse tra signal e background") != std::string::npos);
        const bool alt_is_different = (bg_alt_dir != bg_dir);
        const bool alt_exists =
            std::filesystem::exists(bg_alt_dir) && std::filesystem::is_directory(bg_alt_dir);
        if (dim_mismatch && alt_is_different && alt_exists) {
          std::cerr << "[WARNING] " << riptide::exp2::config_label_str(label)
                    << ": mismatch dimensioni con bg '" << bg_dir.string()
                    << "', provo fallback con '" << bg_alt_dir.string() << "'\n";
          return riptide::exp2::analyze_config(signal_dir, bg_alt_dir, label, optics, stack_cfg,
                                               trace_cfg);
        }
        throw;
      }
    };

    results.push_back(
        analyze_with_bg_fallback(good_focus_dir, riptide::exp2::ConfigLabel::GoodFocus, opt_good));
    results.push_back(analyze_with_bg_fallback(
        good_nofocus_dir, riptide::exp2::ConfigLabel::GoodNoFocus, opt_good_nf));
    results.push_back(
        analyze_with_bg_fallback(bad_focus_dir, riptide::exp2::ConfigLabel::BadFocus, opt_bad));
    results.push_back(analyze_with_bg_fallback(bad_nofocus_dir,
                                               riptide::exp2::ConfigLabel::BadNoFocus, opt_bad_nf));

    for (const auto& r : results)
      riptide::exp2::produce_config_output(r, out_cfg);
    riptide::exp2::produce_summary(results, out_cfg);

    std::cout << "\n--- Summary (stdout) ---\n";
    int ref_width = 0;
    for (const auto& r : results)
      ref_width = std::max(ref_width, r.signal_stack.width);

    std::cout << std::left << std::setw(14) << "config" << std::right << std::setw(10) << "W"
              << std::setw(10) << "H" << std::setw(10) << "xscale" << std::setw(8) << "N"
              << std::setw(18) << "sigma_px" << std::setw(18) << "sigma_refpx"
              << std::setw(12) << "sig_minor" << std::setw(12) << "sig_major"
              << std::setw(10) << "aspect" << std::setw(14) << "chi2/ndof"
              << std::setw(14) << "metric_ref" << std::setw(14) << "dx_det" << "\n";

    for (const auto& r : results) {
      const double scale =
          (ref_width > 0 && r.signal_stack.width > 0)
              ? (static_cast<double>(ref_width) / static_cast<double>(r.signal_stack.width))
              : 1.0;
      const double sigma_ref     = r.trace.sigma_mean * scale;
      const double sigma_ref_err = r.trace.sigma_mean_err * scale;
      const double fwhm_ref      = r.trace.fwhm_mean * scale;
      const double metric_ref    = r.trace.trace_detected
                                       ? (r.centroid_fit.chi2_ndof + fwhm_ref)
                                       : std::numeric_limits<double>::quiet_NaN();

      std::cout << std::left << std::setw(14) << riptide::exp2::config_label_str(r.label)
                << std::right << std::setw(10) << r.signal_stack.width << std::setw(10)
                << r.signal_stack.height << std::setw(10) << std::fixed << std::setprecision(2)
                << scale << std::setw(8) << r.trace.n_valid_slices << std::setw(8) << " "
                << std::setw(7) << std::fixed << std::setprecision(2) << r.trace.sigma_mean << "±"
                << std::setw(4) << std::fixed << std::setprecision(2) << r.trace.sigma_mean_err
                << std::setw(8) << " " << std::setw(7) << std::fixed << std::setprecision(2)
                << sigma_ref << "±" << std::setw(4) << std::fixed << std::setprecision(2)
                << sigma_ref_err
                << std::setw(12) << std::fixed << std::setprecision(2) << r.trace.sigma_minor
                << std::setw(12) << std::fixed << std::setprecision(2) << r.trace.sigma_major
                << std::setw(10) << std::fixed << std::setprecision(2) << r.trace.aspect_ratio
                << std::setw(14) << std::fixed << std::setprecision(3)
                << (r.trace.trace_detected ? r.centroid_fit.chi2_ndof
                                           : std::numeric_limits<double>::quiet_NaN())
                << std::setw(14) << std::fixed << std::setprecision(3) << metric_ref;

      if (std::isfinite(r.dof_delta_mm))
        std::cout << std::setw(14) << std::fixed << std::setprecision(2) << r.dof_delta_mm << "\n";
      else
        std::cout << std::setw(14) << "NaN"
                  << "\n";
    }

    std::cout << "\nOutput:\n";
    std::cout << "  " << (out_cfg.output_dir / "summary.png").string() << "\n";
    if (out_cfg.save_root)
      std::cout << "  " << (out_cfg.output_dir / "exp2_analysis.root").string() << "\n";

    return 0;
  } catch (const std::exception& e) {
    std::cerr << "\n[ERROR] " << e.what() << "\n";
    return 1;
  }
}
