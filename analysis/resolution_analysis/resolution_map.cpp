#include <nlohmann/json.hpp>

#include <TBox.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TTree.h>

#include <omp.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

struct CliConfig {
  std::string input_file  = "output/psf_dof_simulation/psf_dof.root";
  std::string config_file = "config/config.json";
  std::string output_dir  = "output/resolution_analysis";
  std::optional<double> k;
  std::optional<double> scan_min;
  std::optional<double> scan_max;
  std::optional<double> scan_step;
  double lower_percentile = -1.0;
  double upper_percentile = -1.0;
  std::string tsv_out;
  std::string dump_csv;
  int max_entries  = -1;
  int entry_stride = 1;
  int entry_offset = 0;
  int n_jobs = 0; // 0 = tutti i core disponibili
};

static CliConfig parse_args(int argc, char** argv) {
  CliConfig cfg;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    auto next       = [&]() -> std::string {
      if (i + 1 >= argc) {
        std::cerr << "Argomento mancante dopo " << arg << "\n";
        std::exit(1);
      }
      return argv[++i];
    };

    if (arg == "--input" || arg == "-i")
      cfg.input_file = next();
    else if (arg == "--config" || arg == "-c")
      cfg.config_file = next();
    else if (arg == "--output" || arg == "-o")
      cfg.output_dir = next();
    else if (arg == "--k")
      cfg.k = std::stod(next());
    else if (arg == "--scan-min")
      cfg.scan_min = std::stod(next());
    else if (arg == "--scan-max")
      cfg.scan_max = std::stod(next());
    else if (arg == "--scan-step")
      cfg.scan_step = std::stod(next());
    else if (arg == "--low")
      cfg.lower_percentile = std::stod(next());
    else if (arg == "--high")
      cfg.upper_percentile = std::stod(next());
    else if (arg == "--tsv")
      cfg.tsv_out = next();
    else if (arg == "--dump-csv")
      cfg.dump_csv = next();
    else if (arg == "--max-entries")
      cfg.max_entries = std::stoi(next());
    else if (arg == "--entry-stride")
      cfg.entry_stride = std::stoi(next());
    else if (arg == "--entry-offset")
      cfg.entry_offset = std::stoi(next());
    else if (arg == "--jobs")
      cfg.n_jobs = std::stoi(next());
    else if (arg == "--help" || arg == "-h") {
      std::cout << "Uso:\n"
                << "  resolution_map --input psf_dof.root --config config/config.json\n"
                << "                 [--output dir/] [--k 1.414]\n"
                << "                 [--scan-min 180] [--scan-max 400] [--scan-step 0.5]\n"
                << "                 [--low 0.0] [--high 0.0]\n"
                << "                 [--tsv output/resolution_map.tsv]\n"
                << "                 [--dump-csv diag.csv] [--max-entries N]\n";
      std::exit(0);
    } else {
      std::cerr << "Opzione sconosciuta: " << arg << "\n";
      std::exit(1);
    }
  }
  return cfg;
}

static void set_root_style() {
  gStyle->Reset();
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleFont(42, "");
  gStyle->SetStatFont(42);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetNumberContours(255);
}

struct ConfigInfo {
  int config_id    = -1;
  double x1        = 0.0;
  double x2        = 0.0;
  double x_virtual = 0.0;
};

struct RunRow {
  int config_id   = -1;
  int run_id      = -1;
  double x_src    = 0.0;
  double y_src    = 0.0;
  double n_hits   = 0.0;
  double mu_y     = 0.0;
  double sigma_z  = 0.0;
  double sigma_dz = 0.0;
  double cov_z_dz = 0.0;
  double mu_dy    = 0.0;
  double sigma_y  = 0.0;
  double sigma_dy = 0.0;
  double cov_y_dy = 0.0;
};

struct Aggregated {
  double sum_dof      = 0.0;
  int n_dof           = 0;
  double sum_delta_y  = 0.0;
  int n_delta_y       = 0;
  double sum_focus    = 0.0;
  int    n_focus      = 0;
  double sum_EE80     = 0.0;
  int    n_EE80       = 0;
  double dof_mean     = 0.0;
  double delta_y_mean = 0.0;
  double focus_mean   = std::numeric_limits<double>::quiet_NaN();
  double EE80_mean    = std::numeric_limits<double>::quiet_NaN();
};

static std::optional<std::pair<double, double>>
percentile_range(std::vector<double> values, double low_frac, double high_frac) {
  if (values.empty()) {
    return std::nullopt;
  }
  if (!(low_frac >= 0.0 && low_frac < 1.0) || !(high_frac >= 0.0 && high_frac < 1.0)) {
    return std::nullopt;
  }
  if (low_frac + high_frac >= 1.0) {
    return std::nullopt;
  }

  std::sort(values.begin(), values.end());
  size_t N    = values.size();
  size_t i_lo = static_cast<size_t>(low_frac * static_cast<double>(N));
  if (i_lo >= N) {
    i_lo = N - 1;
  }
  size_t i_hi = static_cast<size_t>((1.0 - high_frac) * static_cast<double>(N));
  if (i_hi >= N) {
    i_hi = N - 1;
  }
  if (i_hi < i_lo) {
    return std::nullopt;
  }
  return std::make_pair(values[i_lo], values[i_hi]);
}

static std::vector<double> build_scan(double scan_min, double scan_max, double scan_step) {
  std::vector<double> x;
  for (double v = scan_min; v <= scan_max + 1e-12; v += scan_step) {
    x.push_back(v);
  }
  return x;
}

static std::pair<double, double> compute_sigma_z_min_and_focus(const RunRow& r,
                                                               const ConfigInfo& cfg,
                                                               const std::vector<double>& x_scan) {
  double best_sigma = std::numeric_limits<double>::infinity();
  double best_x     = std::numeric_limits<double>::quiet_NaN();

  double s0_2  = r.sigma_z * r.sigma_z;
  double sdz_2 = r.sigma_dz * r.sigma_dz;
  double cov   = r.cov_z_dz;
  double x0    = cfg.x_virtual;

  for (double x_det : x_scan) {
    double dx = x_det - x0;
    double v  = s0_2 + 2.0 * cov * dx + sdz_2 * dx * dx;
    if (!std::isfinite(v)) {
      continue;
    }
    double sigma = std::sqrt(std::max(0.0, v));
    if (sigma < best_sigma) {
      best_sigma = sigma;
      best_x     = x_det;
    }
  }

  return {best_sigma, best_x};
}

static double compute_dof_from_curve(const RunRow& r, const ConfigInfo& cfg,
                                     const std::vector<double>& x_scan, double k) {
  if (x_scan.size() < 2) {
    return 0.0;
  }

  auto [sigma_min, x_focus] = compute_sigma_z_min_and_focus(r, cfg, x_scan);
  if (!std::isfinite(sigma_min) || !std::isfinite(x_focus)) {
    return 0.0;
  }

  double thr   = k * sigma_min;
  double s0_2  = r.sigma_z * r.sigma_z;
  double sdz_2 = r.sigma_dz * r.sigma_dz;
  double cov   = r.cov_z_dz;
  double x0    = cfg.x_virtual;

  std::vector<double> sigma(x_scan.size(), 0.0);
  for (size_t i = 0; i < x_scan.size(); ++i) {
    double dx = x_scan[i] - x0;
    double v  = s0_2 + 2.0 * cov * dx + sdz_2 * dx * dx;
    sigma[i]  = std::sqrt(std::max(0.0, v));
  }

  size_t i_min =
      std::distance(sigma.begin(), std::min_element(sigma.begin(), sigma.end(),
                                                    [](double a, double b) { return a < b; }));
  int i_lo = static_cast<int>(i_min);
  while (i_lo - 1 >= 0 && sigma[static_cast<size_t>(i_lo - 1)] < thr) {
    --i_lo;
  }
  int i_hi = static_cast<int>(i_min);
  while (i_hi + 1 < static_cast<int>(sigma.size()) && sigma[static_cast<size_t>(i_hi + 1)] < thr) {
    ++i_hi;
  }

  return x_scan[static_cast<size_t>(i_hi)] - x_scan[static_cast<size_t>(i_lo)];
}

static std::optional<double> compute_delta_y_min_at_focus(const RunRow& r, const ConfigInfo& cfg,
                                                          const std::vector<double>& x_scan,
                                                          double k) {
  auto [sigma_min, x_focus] = compute_sigma_z_min_and_focus(r, cfg, x_scan);
  if (!std::isfinite(sigma_min) || !std::isfinite(x_focus)) {
    return std::nullopt;
  }

  double dx = x_focus - cfg.x_virtual;
  if (!(std::abs(r.y_src) > 1e-12)) {
    return std::nullopt;
  }

  double y_bar_focus = r.mu_y + r.mu_dy * dx;
  if (!std::isfinite(y_bar_focus) || std::abs(y_bar_focus) < 1e-9) {
    return std::nullopt;
  }

  double Mbar = y_bar_focus / r.y_src;
  if (!std::isfinite(Mbar) || std::abs(Mbar) < 1e-18) {
    return std::nullopt;
  }

  return k * sigma_min / std::abs(Mbar);
}

static std::optional<double> compute_EE80_at_focus(const RunRow& r, const ConfigInfo& cfg,
                                                    const std::vector<double>& x_scan) {
  auto [sigma_z_min, x_focus] = compute_sigma_z_min_and_focus(r, cfg, x_scan);
  if (!std::isfinite(sigma_z_min) || !std::isfinite(x_focus)) return std::nullopt;

  double dx = x_focus - cfg.x_virtual;
  double v_y =
      r.sigma_y * r.sigma_y + 2.0 * r.cov_y_dy * dx + r.sigma_dy * r.sigma_dy * dx * dx;
  double sigma_y_focus = std::sqrt(std::max(0.0, v_y));

  // approssimazione 2D gaussiana isotropa: R80 = sqrt(-2*ln(0.2))*sigma ≈ 1.7941*sigma
  double sigma_rms = std::sqrt((sigma_y_focus * sigma_y_focus + sigma_z_min * sigma_z_min) / 2.0);
  constexpr double k_EE80 = 2.0 * 1.7941;
  return k_EE80 * sigma_rms;
}

int main(int argc, char** argv) {
  using json = nlohmann::json;

  CliConfig cli = parse_args(argc, argv);
  set_root_style();
  std::filesystem::create_directories(cli.output_dir);

  std::ifstream jf(cli.config_file);
  if (!jf.is_open()) {
    std::cerr << "Errore: impossibile aprire " << cli.config_file << "\n";
    return 1;
  }
  json config;
  jf >> config;

  double scan_min  = cli.scan_min.value_or(config.value("dof_x_scan_min", 180.0));
  double scan_max  = cli.scan_max.value_or(config.value("dof_x_scan_max", 400.0));
  double scan_step = cli.scan_step.value_or(config.value("dof_x_scan_step", 0.5));
  double k         = cli.k.value_or(config.value("dof_k_threshold", std::sqrt(2.0)));
  double dx_grid   = config.value("dx", 1.0);
  double lower_percentile =
      (cli.lower_percentile >= 0.0) ? cli.lower_percentile : config.value("lower_percentile", 0.0);
  double upper_percentile =
      (cli.upper_percentile >= 0.0) ? cli.upper_percentile : config.value("upper_percentile", 0.0);
  double lens_det_gap = config.value("lens_det_gap", 0.0);

  if (!(scan_step > 0.0)) {
    std::cerr << "Errore: scan_step deve essere > 0\n";
    return 1;
  }
  if (!(scan_max > scan_min)) {
    std::cerr << "Errore: scan_max deve essere > scan_min\n";
    return 1;
  }
  if (!(k > 0.0)) {
    std::cerr << "Errore: k deve essere > 0\n";
    return 1;
  }
  if (!(lower_percentile >= 0.0 && lower_percentile < 1.0)
      || !(upper_percentile >= 0.0 && upper_percentile < 1.0)
      || (lower_percentile + upper_percentile >= 1.0)) {
    std::cerr << "Errore: percentili non validi (richiesto 0<=low<1, 0<=high<1, low+high<1)\n";
    return 1;
  }

  TFile file(cli.input_file.c_str(), "READ");
  if (!file.IsOpen()) {
    std::cerr << "Errore: impossibile aprire " << cli.input_file << "\n";
    return 1;
  }

  TTree* tree_cfg  = (TTree*)file.Get("PsfDofConfigs");
  TTree* tree_runs = (TTree*)file.Get("PsfDofRuns");
  if (!tree_cfg || !tree_runs) {
    std::cerr << "Errore: TTree 'PsfDofConfigs' o 'PsfDofRuns' non trovato\n";
    return 1;
  }

  int config_id_cfg    = 0;
  double x1_cfg        = 0.0;
  double x2_cfg        = 0.0;
  double x_virtual_cfg = 0.0;
  tree_cfg->SetBranchAddress("config_id", &config_id_cfg);
  tree_cfg->SetBranchAddress("x1", &x1_cfg);
  tree_cfg->SetBranchAddress("x2", &x2_cfg);
  tree_cfg->SetBranchAddress("x_virtual", &x_virtual_cfg);

  std::unordered_map<int, ConfigInfo> config_map;
  config_map.reserve(static_cast<size_t>(tree_cfg->GetEntries()));
  double x1_min_data = std::numeric_limits<double>::infinity();
  double x1_max_data = -std::numeric_limits<double>::infinity();
  double x2_min_data = std::numeric_limits<double>::infinity();
  double x2_max_data = -std::numeric_limits<double>::infinity();

  for (int i = 0; i < tree_cfg->GetEntries(); ++i) {
    tree_cfg->GetEntry(i);
    ConfigInfo info{config_id_cfg, x1_cfg, x2_cfg, x_virtual_cfg};
    config_map[config_id_cfg] = info;
    x1_min_data               = std::min(x1_min_data, x1_cfg);
    x1_max_data               = std::max(x1_max_data, x1_cfg);
    x2_min_data               = std::min(x2_min_data, x2_cfg);
    x2_max_data               = std::max(x2_max_data, x2_cfg);
  }

  std::unordered_map<int, std::vector<double>> x_scan_per_cfg;
  x_scan_per_cfg.reserve(config_map.size());
  for (const auto& [id, c] : config_map) {
    double cfg_scan_min = std::max(scan_min, c.x2 + lens_det_gap);
    auto cfg_scan = build_scan(cfg_scan_min, scan_max, scan_step);
    if (!cfg_scan.empty()) {
      x_scan_per_cfg[id] = std::move(cfg_scan);
    }
  }

  RunRow rr;
  tree_runs->SetBranchAddress("config_id", &rr.config_id);
  tree_runs->SetBranchAddress("run_id", &rr.run_id);
  tree_runs->SetBranchAddress("x_source", &rr.x_src);
  tree_runs->SetBranchAddress("y_source", &rr.y_src);
  tree_runs->SetBranchAddress("n_hits", &rr.n_hits);
  tree_runs->SetBranchAddress("mu_y", &rr.mu_y);
  tree_runs->SetBranchAddress("sigma_z", &rr.sigma_z);
  tree_runs->SetBranchAddress("sigma_dz", &rr.sigma_dz);
  tree_runs->SetBranchAddress("cov_z_dz", &rr.cov_z_dz);
  tree_runs->SetBranchAddress("mu_dy", &rr.mu_dy);
  tree_runs->SetBranchAddress("sigma_y",  &rr.sigma_y);
  tree_runs->SetBranchAddress("sigma_dy", &rr.sigma_dy);
  tree_runs->SetBranchAddress("cov_y_dy", &rr.cov_y_dy);

  if (cli.n_jobs > 0)
    omp_set_num_threads(cli.n_jobs);

  int total_entries = tree_runs->GetEntries();
  int stride        = std::max(1, cli.entry_stride);
  int offset        = std::max(0, std::min(cli.entry_offset, stride - 1));
  int i_start       = offset;
  int n_to_process  = (cli.max_entries > 0)
                          ? std::min(i_start + cli.max_entries * stride, total_entries)
                          : total_entries;

  // Fase 1: pre-caricamento sequenziale (GetEntry non è thread-safe)
  struct PreloadedRow {
    RunRow rr;
    const ConfigInfo* cfg;
    const std::vector<double>* x_scan;
  };
  std::vector<PreloadedRow> rows;
  rows.reserve(static_cast<size_t>((n_to_process - i_start + stride - 1) / stride));

  for (int i = i_start; i < n_to_process; i += stride) {
    tree_runs->GetEntry(i);
    if (!(rr.n_hits > 0.0))
      continue;
    auto it = config_map.find(rr.config_id);
    if (it == config_map.end())
      continue;
    auto scan_it = x_scan_per_cfg.find(rr.config_id);
    if (scan_it == x_scan_per_cfg.end())
      continue;
    rows.push_back({rr, &it->second, &scan_it->second});
  }

  // Fase 2: calcolo parallelo per riga
  struct RowResult {
    int config_id;
    double dof;
    bool dof_valid;
    double focus;
    bool focus_valid;
    std::optional<double> delta;
    std::optional<double> ee80;
    // campi per --dump-csv
    double x_focus_analytical;
    bool focus_in_scan;
    double sigma_min_comp;
    double mag_M;
  };
  std::vector<RowResult> row_results(rows.size());

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < static_cast<int>(rows.size()); ++i) {
    const auto& row = rows[i];
    const RunRow& r = row.rr;
    const ConfigInfo& c = *row.cfg;
    const std::vector<double>& xs = *row.x_scan;
    RowResult& res = row_results[i];

    res.config_id = r.config_id;

    double dof_val  = compute_dof_from_curve(r, c, xs, k);
    res.dof         = dof_val;
    res.dof_valid   = std::isfinite(dof_val) && dof_val > 0.0;

    res.focus_valid = false;
    res.focus       = std::numeric_limits<double>::quiet_NaN();
    if (r.sigma_dz > 1e-12) {
      double xf = c.x_virtual - r.cov_z_dz / (r.sigma_dz * r.sigma_dz);
      if (std::isfinite(xf)) {
        res.focus       = xf;
        res.focus_valid = true;
      }
    }

    res.delta = compute_delta_y_min_at_focus(r, c, xs, k);
    res.ee80  = compute_EE80_at_focus(r, c, xs);

    res.x_focus_analytical = res.focus;
    res.focus_in_scan      = res.focus_valid && res.focus >= scan_min && res.focus <= scan_max;

    auto [sm, xf_comp] = compute_sigma_z_min_and_focus(r, c, xs);
    res.sigma_min_comp  = sm;

    res.mag_M = std::numeric_limits<double>::quiet_NaN();
    if (std::isfinite(xf_comp) && std::abs(r.y_src) > 1e-12) {
      double dx_f      = xf_comp - c.x_virtual;
      double y_bar_foc = r.mu_y + r.mu_dy * dx_f;
      if (std::isfinite(y_bar_foc))
        res.mag_M = y_bar_foc / r.y_src;
    }
  }

  // Fase 3: aggregazione sequenziale e dump CSV
  std::unordered_map<int, Aggregated> agg;
  agg.reserve(config_map.size());

  std::ofstream dump_file;
  bool do_dump = !cli.dump_csv.empty();
  if (do_dump) {
    std::filesystem::path dp(cli.dump_csv);
    if (dp.has_parent_path())
      std::filesystem::create_directories(dp.parent_path());
    dump_file.open(cli.dump_csv);
    dump_file << "config_id\tx1\tx2\tx_virtual\t"
              << "run_id\ty_src\tn_hits\t"
              << "sigma_z\tsigma_dz\tcov_z_dz\t"
              << "x_focus_analytical\tfocus_in_scan\t"
              << "sigma_min_computed\tdof\tdelta_y\tmag_M\n";
  }

  for (int i = 0; i < static_cast<int>(rows.size()); ++i) {
    const RunRow& r     = rows[i].rr;
    const ConfigInfo& c = *rows[i].cfg;
    const RowResult& res = row_results[i];

    if (res.dof_valid) {
      agg[res.config_id].sum_dof += res.dof;
      agg[res.config_id].n_dof  += 1;
    }
    if (res.focus_valid) {
      agg[res.config_id].sum_focus += res.focus;
      agg[res.config_id].n_focus  += 1;
    }
    if (res.delta.has_value() && std::isfinite(*res.delta) && *res.delta > 0.0) {
      agg[res.config_id].sum_delta_y += *res.delta;
      agg[res.config_id].n_delta_y  += 1;
    }
    if (res.ee80.has_value() && std::isfinite(*res.ee80) && *res.ee80 > 0.0) {
      agg[res.config_id].sum_EE80 += *res.ee80;
      agg[res.config_id].n_EE80  += 1;
    }

    if (do_dump) {
      dump_file << c.config_id << "\t" << c.x1 << "\t" << c.x2 << "\t" << c.x_virtual << "\t"
                << r.run_id << "\t" << r.y_src << "\t" << r.n_hits << "\t"
                << r.sigma_z << "\t" << r.sigma_dz << "\t" << r.cov_z_dz << "\t"
                << res.x_focus_analytical << "\t" << (res.focus_in_scan ? 1 : 0) << "\t"
                << res.sigma_min_comp << "\t"
                << (std::isfinite(res.dof) ? res.dof : std::numeric_limits<double>::quiet_NaN()) << "\t"
                << (res.delta.has_value() ? *res.delta : std::numeric_limits<double>::quiet_NaN()) << "\t"
                << res.mag_M << "\n";
    }
  }

  for (auto& [id, a] : agg) {
    a.dof_mean     = (a.n_dof > 0) ? (a.sum_dof / static_cast<double>(a.n_dof)) : 0.0;
    a.delta_y_mean = (a.n_delta_y > 0) ? (a.sum_delta_y / static_cast<double>(a.n_delta_y)) : 0.0;
    a.focus_mean   = (a.n_focus > 0) ? (a.sum_focus / static_cast<double>(a.n_focus))
                                      : std::numeric_limits<double>::quiet_NaN();
    a.EE80_mean    = (a.n_EE80 > 0) ? (a.sum_EE80 / static_cast<double>(a.n_EE80))
                                     : std::numeric_limits<double>::quiet_NaN();
  }

  int n_bins_x1 = std::round((x1_max_data - x1_min_data) / dx_grid) + 1;
  int n_bins_x2 = std::round((x2_max_data - x2_min_data) / dx_grid) + 1;

  TH2D h_dof_mean("h_dof_mean", ";x1 [mm];x2 [mm];DoF_{mean} [mm]", n_bins_x1,
                  x1_min_data - dx_grid / 2.0, x1_max_data + dx_grid / 2.0, n_bins_x2,
                  x2_min_data - dx_grid / 2.0, x2_max_data + dx_grid / 2.0);
  TH2D h_delta_y_mean("h_delta_y_mean", ";x1 [mm];x2 [mm];#delta y_{min,mean} [mm]", n_bins_x1,
                      x1_min_data - dx_grid / 2.0, x1_max_data + dx_grid / 2.0, n_bins_x2,
                      x2_min_data - dx_grid / 2.0, x2_max_data + dx_grid / 2.0);
  TH2D h_EE80_mean("h_EE80_mean", ";x1 [mm];x2 [mm];EE80_{mean} [mm]", n_bins_x1,
                   x1_min_data - dx_grid / 2.0, x1_max_data + dx_grid / 2.0, n_bins_x2,
                   x2_min_data - dx_grid / 2.0, x2_max_data + dx_grid / 2.0);

  TH2D h_mask_cfg("h_mask_cfg", "", n_bins_x1, x1_min_data - dx_grid / 2.0,
                  x1_max_data + dx_grid / 2.0, n_bins_x2, x2_min_data - dx_grid / 2.0,
                  x2_max_data + dx_grid / 2.0);
  TH2D h_mask_dof("h_mask_dof", "", n_bins_x1, x1_min_data - dx_grid / 2.0,
                  x1_max_data + dx_grid / 2.0, n_bins_x2, x2_min_data - dx_grid / 2.0,
                  x2_max_data + dx_grid / 2.0);
  TH2D h_mask_delta("h_mask_delta", "", n_bins_x1, x1_min_data - dx_grid / 2.0,
                    x1_max_data + dx_grid / 2.0, n_bins_x2, x2_min_data - dx_grid / 2.0,
                    x2_max_data + dx_grid / 2.0);
  TH2D h_mask_EE80("h_mask_EE80", "", n_bins_x1, x1_min_data - dx_grid / 2.0,
                   x1_max_data + dx_grid / 2.0, n_bins_x2, x2_min_data - dx_grid / 2.0,
                   x2_max_data + dx_grid / 2.0);
  TH2D h_focus_warn("h_focus_warn", "", n_bins_x1, x1_min_data - dx_grid / 2.0,
                    x1_max_data + dx_grid / 2.0, n_bins_x2, x2_min_data - dx_grid / 2.0,
                    x2_max_data + dx_grid / 2.0);
  for (int bx = 1; bx <= n_bins_x1; ++bx) {
    for (int by = 1; by <= n_bins_x2; ++by) {
      h_mask_cfg.SetBinContent(bx, by, 1.0);
      h_mask_dof.SetBinContent(bx, by, 1.0);
      h_mask_delta.SetBinContent(bx, by, 1.0);
      h_mask_EE80.SetBinContent(bx, by, 1.0);
    }
  }

  double dof_min = std::numeric_limits<double>::infinity();
  double dof_max  = -std::numeric_limits<double>::infinity();
  double dy_min   = std::numeric_limits<double>::infinity();
  double dy_max   = -std::numeric_limits<double>::infinity();
  double ee80_min = std::numeric_limits<double>::infinity();
  double ee80_max = -std::numeric_limits<double>::infinity();

  for (const auto& [id, c] : config_map) {
    int bx = h_dof_mean.GetXaxis()->FindBin(c.x1);
    int by = h_dof_mean.GetYaxis()->FindBin(c.x2);
    h_mask_cfg.SetBinContent(bx, by, 0.0);

    auto it = agg.find(id);
    if (it == agg.end()) {
      continue;
    }

    if (std::isfinite(it->second.focus_mean) && it->second.focus_mean < c.x2 + lens_det_gap) {
      h_focus_warn.SetBinContent(bx, by, 1.0);
    }

    if (it->second.n_dof > 0) {
      h_dof_mean.SetBinContent(bx, by, it->second.dof_mean);
      h_mask_dof.SetBinContent(bx, by, 0.0);
      dof_min = std::min(dof_min, it->second.dof_mean);
      dof_max = std::max(dof_max, it->second.dof_mean);
    }
    if (it->second.n_delta_y > 0) {
      h_delta_y_mean.SetBinContent(bx, by, it->second.delta_y_mean);
      h_mask_delta.SetBinContent(bx, by, 0.0);
      dy_min = std::min(dy_min, it->second.delta_y_mean);
      dy_max = std::max(dy_max, it->second.delta_y_mean);
    }
    if (it->second.n_EE80 > 0 && std::isfinite(it->second.EE80_mean)) {
      h_EE80_mean.SetBinContent(bx, by, it->second.EE80_mean);
      h_mask_EE80.SetBinContent(bx, by, 0.0);
      ee80_min = std::min(ee80_min, it->second.EE80_mean);
      ee80_max = std::max(ee80_max, it->second.EE80_mean);
    }
  }

  if (std::isfinite(dof_min) && std::isfinite(dof_max)) {
    h_dof_mean.SetMinimum(dof_min);
    h_dof_mean.SetMaximum(dof_max);
  }
  if (std::isfinite(dy_min) && std::isfinite(dy_max)) {
    h_delta_y_mean.SetMinimum(dy_min);
    h_delta_y_mean.SetMaximum(dy_max);
  }
  if (std::isfinite(ee80_min) && std::isfinite(ee80_max)) {
    h_EE80_mean.SetMinimum(ee80_min);
    h_EE80_mean.SetMaximum(ee80_max);
  }

  {
    std::vector<double> dof_vals;
    dof_vals.reserve(static_cast<size_t>(n_bins_x1 * n_bins_x2));
    std::vector<double> dy_vals;
    dy_vals.reserve(static_cast<size_t>(n_bins_x1 * n_bins_x2));

    for (int by = 1; by <= n_bins_x2; ++by) {
      for (int bx = 1; bx <= n_bins_x1; ++bx) {
        if (h_mask_dof.GetBinContent(bx, by) < 0.5) {
          double v = h_dof_mean.GetBinContent(bx, by);
          if (std::isfinite(v)) {
            dof_vals.push_back(v);
          }
        }
        if (h_mask_delta.GetBinContent(bx, by) < 0.5) {
          double v = h_delta_y_mean.GetBinContent(bx, by);
          if (std::isfinite(v)) {
            dy_vals.push_back(v);
          }
        }
      }
    }

    auto r_dof = percentile_range(std::move(dof_vals), lower_percentile, upper_percentile);
    if (r_dof.has_value() && r_dof->second > r_dof->first) {
      h_dof_mean.SetMinimum(r_dof->first);
      h_dof_mean.SetMaximum(r_dof->second);
    }

    auto r_dy = percentile_range(std::move(dy_vals), lower_percentile, upper_percentile);
    if (r_dy.has_value() && r_dy->second > r_dy->first) {
      h_delta_y_mean.SetMinimum(r_dy->first);
      h_delta_y_mean.SetMaximum(r_dy->second);
    }
  }

  auto build_warn_boxes = [&]() {
    std::vector<std::unique_ptr<TBox>> boxes;
    boxes.reserve(static_cast<size_t>(n_bins_x1 * n_bins_x2));
    for (int bx = 1; bx <= n_bins_x1; ++bx) {
      double x_lo = h_focus_warn.GetXaxis()->GetBinLowEdge(bx);
      double x_hi = h_focus_warn.GetXaxis()->GetBinUpEdge(bx);
      for (int by = 1; by <= n_bins_x2; ++by) {
        if (h_focus_warn.GetBinContent(bx, by) < 0.5) {
          continue;
        }
        double y_lo = h_focus_warn.GetYaxis()->GetBinLowEdge(by);
        double y_hi = h_focus_warn.GetYaxis()->GetBinUpEdge(by);
        auto box    = std::make_unique<TBox>(x_lo, y_lo, x_hi, y_hi);
        box->SetFillColorAlpha(kOrange + 1, 0.35);
        box->SetLineColor(0);
        box->Draw("same");
        boxes.push_back(std::move(box));
      }
    }
    return boxes;
  };

  auto build_mask_boxes = [&](TH2D& h_mask) {
    std::vector<std::unique_ptr<TBox>> boxes;
    boxes.reserve(static_cast<size_t>(n_bins_x1 * n_bins_x2));
    for (int bx = 1; bx <= n_bins_x1; ++bx) {
      double x_lo = h_mask.GetXaxis()->GetBinLowEdge(bx);
      double x_hi = h_mask.GetXaxis()->GetBinUpEdge(bx);
      for (int by = 1; by <= n_bins_x2; ++by) {
        if (h_mask.GetBinContent(bx, by) < 0.5) {
          continue;
        }
        double y_lo = h_mask.GetYaxis()->GetBinLowEdge(by);
        double y_hi = h_mask.GetYaxis()->GetBinUpEdge(by);
        auto box    = std::make_unique<TBox>(x_lo, y_lo, x_hi, y_hi);
        box->SetFillColor(kBlack);
        box->SetLineColor(kBlack);
        box->Draw();
        boxes.push_back(std::move(box));
      }
    }
    return boxes;
  };

  auto save_map = [&](TH2D& h, TH2D& h_mask, const std::string& name, int palette) {
    gStyle->SetPalette(palette);
    TCanvas c(("c_" + name).c_str(), name.c_str(), 1100, 900);
    c.SetLeftMargin(0.16);
    c.SetBottomMargin(0.14);
    c.SetRightMargin(0.16);
    c.SetTopMargin(0.08);
    h.Draw("COLZ");
    auto warn_boxes = build_warn_boxes();
    auto boxes      = build_mask_boxes(h_mask);
    std::string out = (std::filesystem::path(cli.output_dir) / (name + ".png")).string();
    c.SaveAs(out.c_str());
  };

  save_map(h_dof_mean, h_mask_cfg, "resolution_dof_mean_map", kViridis);
  save_map(h_delta_y_mean, h_mask_delta, "resolution_delta_y_min_mean_map", kViridis);
  save_map(h_EE80_mean, h_mask_EE80, "resolution_EE80_mean_map", kViridis);

  if (!cli.tsv_out.empty()) {
    std::filesystem::path tsv_path = cli.tsv_out;
    if (tsv_path.has_parent_path()) {
      std::filesystem::create_directories(tsv_path.parent_path());
    }
    std::ofstream out(tsv_path);
    out << "x1\tx2\tdof_mean\tdelta_y_min_mean\tEE80_mean\tconfig_id\tn_runs_dof\tn_runs_delta_y\tn_runs_EE80\n";
    for (const auto& [id, c] : config_map) {
      auto it = agg.find(id);
      if (it == agg.end()) {
        continue;
      }
      out << c.x1 << "\t" << c.x2 << "\t" << it->second.dof_mean << "\t" << it->second.delta_y_mean
          << "\t" << it->second.EE80_mean << "\t" << id << "\t" << it->second.n_dof << "\t"
          << it->second.n_delta_y << "\t" << it->second.n_EE80 << "\n";
    }
  }

  return 0;
}
