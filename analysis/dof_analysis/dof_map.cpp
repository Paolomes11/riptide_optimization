#include <nlohmann/json.hpp>

#include <TBox.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TColor.h>
#include <TFile.h>
#include <TH2D.h>
#include <TMarker.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

struct CliConfig {
  std::string input_file  = "output/dof_simulation/focal.root";
  std::string config_file = "config/config.json";
  std::string output_dir  = "output/dof_analysis";
  std::optional<double> scan_min;
  std::optional<double> scan_max;
  std::optional<double> scan_step;
  std::optional<double> k_threshold;
  std::optional<double> core_fraction;
  std::optional<double> m_target;
  std::string tsv_out;
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
    else if (arg == "--scan-min")
      cfg.scan_min = std::stod(next());
    else if (arg == "--scan-max")
      cfg.scan_max = std::stod(next());
    else if (arg == "--scan-step")
      cfg.scan_step = std::stod(next());
    else if (arg == "--k")
      cfg.k_threshold = std::stod(next());
    else if (arg == "--core-fraction")
      cfg.core_fraction = std::stod(next());
    else if (arg == "--m-target")
      cfg.m_target = std::stod(next());
    else if (arg == "--tsv")
      cfg.tsv_out = next();
    else if (arg == "--help" || arg == "-h") {
      std::cout << "Uso:\n"
                << "  dof_map --input focal.root --config config/config.json [--output dir/]\n"
                << "          [--scan-min 100] [--scan-max 350] [--scan-step 0.5] [--k 1.414]\n"
                << "          [--core-fraction 1.0] [--m-target 0.1333]\n"
                << "          [--tsv output/dof_map.tsv]\n";
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
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetNumberContours(255);
  gStyle->SetGridColor(kGray + 1);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
}

static double weighted_std(const std::vector<double>& x, const std::vector<double>& w) {
  double sum_w = 0.0;
  double sum_x = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    sum_w += w[i];
    sum_x += w[i] * x[i];
  }
  if (sum_w <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  double mean = sum_x / sum_w;
  double var  = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    double d = x[i] - mean;
    var += w[i] * d * d;
  }
  var /= sum_w;
  return std::sqrt(std::max(0.0, var));
}

static double weighted_mean(const std::vector<double>& x, const std::vector<double>& w) {
  double sum_w = 0.0;
  double sum_x = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    sum_w += w[i];
    sum_x += w[i] * x[i];
  }
  if (sum_w <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return sum_x / sum_w;
}

static std::vector<size_t> select_core_indices_around_mean(const std::vector<double>& y,
                                                           const std::vector<double>& z,
                                                           const std::vector<double>& w,
                                                           double core_fraction) {
  if (y.empty() || z.size() != y.size() || w.size() != y.size()) {
    return {};
  }
  if (core_fraction >= 1.0) {
    std::vector<size_t> idx(y.size());
    std::iota(idx.begin(), idx.end(), 0);
    return idx;
  }
  if (core_fraction <= 0.0) {
    return {};
  }

  double sum_w = 0.0;
  double sum_y = 0.0;
  double sum_z = 0.0;
  for (size_t i = 0; i < y.size(); ++i) {
    sum_w += w[i];
    sum_y += w[i] * y[i];
    sum_z += w[i] * z[i];
  }
  if (sum_w <= 0.0) {
    return {};
  }
  double mu_y = sum_y / sum_w;
  double mu_z = sum_z / sum_w;

  double cov_yy = 0.0;
  double cov_zz = 0.0;
  double cov_yz = 0.0;
  for (size_t i = 0; i < y.size(); ++i) {
    double dy = y[i] - mu_y;
    double dz = z[i] - mu_z;
    cov_yy += w[i] * dy * dy;
    cov_zz += w[i] * dz * dz;
    cov_yz += w[i] * dy * dz;
  }
  cov_yy /= sum_w;
  cov_zz /= sum_w;
  cov_yz /= sum_w;

  double scale = 0.5 * (cov_yy + cov_zz);
  double eps   = 1e-9 + 1e-6 * std::max(0.0, scale);
  cov_yy += eps;
  cov_zz += eps;

  double det = cov_yy * cov_zz - cov_yz * cov_yz;
  if (!std::isfinite(det) || std::abs(det) < 1e-18) {
    cov_yz = 0.0;
    det    = cov_yy * cov_zz;
    if (!std::isfinite(det) || std::abs(det) < 1e-18) {
      return {};
    }
  }

  double inv_yy = cov_zz / det;
  double inv_zz = cov_yy / det;
  double inv_yz = -cov_yz / det;

  std::vector<size_t> idx(y.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
    double dya = y[a] - mu_y;
    double dza = z[a] - mu_z;
    double dyb = y[b] - mu_y;
    double dzb = z[b] - mu_z;
    double d2a = dya * dya * inv_yy + dza * dza * inv_zz + 2.0 * dya * dza * inv_yz;
    double d2b = dyb * dyb * inv_yy + dzb * dzb * inv_zz + 2.0 * dyb * dzb * inv_yz;
    return d2a < d2b;
  });

  std::vector<size_t> keep;
  keep.reserve(idx.size());
  double target = core_fraction * sum_w;
  double acc    = 0.0;
  for (size_t k : idx) {
    keep.push_back(k);
    acc += w[k];
    if (acc >= target) {
      break;
    }
  }
  std::sort(keep.begin(), keep.end());
  return keep;
}

static void apply_selection(std::vector<double>& y0, std::vector<double>& z0,
                            std::vector<double>& dy, std::vector<double>& dz,
                            std::vector<double>& w, std::vector<double>& ysrc,
                            const std::vector<size_t>& keep) {
  auto pick = [&](std::vector<double>& v) {
    std::vector<double> out;
    out.reserve(keep.size());
    for (size_t i : keep) {
      out.push_back(v[i]);
    }
    v.swap(out);
  };
  pick(y0);
  pick(z0);
  pick(dy);
  pick(dz);
  pick(w);
  pick(ysrc);
}

static double weighted_percentile(const std::vector<double>& x, const std::vector<double>& w,
                                  double p) {
  if (x.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  std::vector<size_t> idx(x.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return x[a] < x[b]; });

  double sum_w = 0.0;
  for (double wi : w) {
    sum_w += wi;
  }
  if (sum_w <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double target = p * sum_w;
  double acc    = 0.0;
  for (size_t j = 0; j < idx.size(); ++j) {
    acc += w[idx[j]];
    if (acc >= target) {
      return x[idx[j]];
    }
  }
  return x[idx.back()];
}

static std::optional<double> quadratic_vertex_from_three_points(double x0, double y0, double x1,
                                                                double y1, double x2, double y2) {
  double d0 = (x0 - x1) * (x0 - x2);
  double d1 = (x1 - x0) * (x1 - x2);
  double d2 = (x2 - x0) * (x2 - x1);
  if (std::abs(d0) < 1e-12 || std::abs(d1) < 1e-12 || std::abs(d2) < 1e-12) {
    return std::nullopt;
  }

  double a = (y0 / d0) + (y1 / d1) + (y2 / d2);
  double b = (-y0 * (x1 + x2) / d0) + (-y1 * (x0 + x2) / d1) + (-y2 * (x0 + x1) / d2);
  if (!std::isfinite(a) || !std::isfinite(b) || std::abs(a) < 1e-18) {
    return std::nullopt;
  }
  return -b / (2.0 * a);
}

struct ConfigInfo {
  int config_id    = -1;
  double x1        = 0.0;
  double x2        = 0.0;
  double x_virtual = 0.0;
};

struct ResultRow {
  int config_id           = -1;
  double x1               = 0.0;
  double x2               = 0.0;
  double x_focus          = 0.0;
  double x_focus_scan     = 0.0;
  double dof              = 0.0;
  double M                = 0.0;
  double M_abs_err        = 0.0;
  double stripe_width     = 0.0;
  int within_photocathode = 0;
  double n_rays           = 0.0;
  double n_rays_core      = 0.0;
  bool focus_before_lens2 = false;
};

static ResultRow analyze_config(const ConfigInfo& cfg, double n_rays, const std::vector<double>& y0,
                                const std::vector<double>& z0, const std::vector<double>& dy,
                                const std::vector<double>& dz, const std::vector<double>& w,
                                const std::vector<double>& ysrc, const std::vector<double>& x_scan,
                                double k_threshold, double core_fraction, double m_target,
                                double sigma_y_src_theory) {
  ResultRow out;
  out.config_id = cfg.config_id;
  out.x1        = cfg.x1;
  out.x2        = cfg.x2;
  out.n_rays    = n_rays;

  if (x_scan.size() < 2) {
    return out;
  }

  std::vector<double> tmp(z0.size(), 0.0);
  std::vector<double> sigma_z(x_scan.size(), 0.0);
  for (size_t ix = 0; ix < x_scan.size(); ++ix) {
    double dx_det = x_scan[ix] - cfg.x_virtual;
    for (size_t i = 0; i < z0.size(); ++i) {
      tmp[i] = z0[i] + dz[i] * dx_det;
    }
    sigma_z[ix] = weighted_std(tmp, w);
  }

  size_t i_min =
      std::distance(sigma_z.begin(), std::min_element(sigma_z.begin(), sigma_z.end(),
                                                      [](double a, double b) { return a < b; }));
  out.x_focus_scan   = x_scan[i_min];
  out.x_focus        = out.x_focus_scan;
  double sigma_z_min = sigma_z[i_min];

  if (i_min == 0 && x_scan.size() >= 3) {
    double x0  = x_scan[0];
    double x1  = x_scan[1];
    double x2  = x_scan[2];
    double y0q = sigma_z[0] * sigma_z[0];
    double y1q = sigma_z[1] * sigma_z[1];
    double y2q = sigma_z[2] * sigma_z[2];
    auto xv    = quadratic_vertex_from_three_points(x0, y0q, x1, y1q, x2, y2q);
    if (xv.has_value() && std::isfinite(*xv)) {
      out.x_focus = *xv;
    }
  }
  out.focus_before_lens2 = (out.x_focus < cfg.x2);

  double thr = k_threshold * sigma_z_min;
  int i_lo   = static_cast<int>(i_min);
  while (i_lo - 1 >= 0 && sigma_z[static_cast<size_t>(i_lo - 1)] < thr) {
    --i_lo;
  }
  int i_hi = static_cast<int>(i_min);
  while (i_hi + 1 < static_cast<int>(sigma_z.size())
         && sigma_z[static_cast<size_t>(i_hi + 1)] < thr) {
    ++i_hi;
  }
  double x_lo = x_scan[static_cast<size_t>(i_lo)];
  double x_hi = x_scan[static_cast<size_t>(i_hi)];
  out.dof     = x_hi - x_lo;

  std::vector<double> y0_f   = y0;
  std::vector<double> z0_f   = z0;
  std::vector<double> dy_f   = dy;
  std::vector<double> dz_f   = dz;
  std::vector<double> w_f    = w;
  std::vector<double> ysrc_f = ysrc;

  double n_w_before = 0.0;
  for (double wi : w_f) {
    n_w_before += wi;
  }

  std::vector<double> y_focus(y0_f.size(), 0.0);
  std::vector<double> z_focus(z0_f.size(), 0.0);
  double dx_focus = out.x_focus - cfg.x_virtual;
  for (size_t i = 0; i < y0_f.size(); ++i) {
    y_focus[i] = y0_f[i] + dy_f[i] * dx_focus;
    z_focus[i] = z0_f[i] + dz_f[i] * dx_focus;
  }

  if (core_fraction < 1.0) {
    auto keep = select_core_indices_around_mean(y_focus, z_focus, w_f, core_fraction);
    apply_selection(y0_f, z0_f, dy_f, dz_f, w_f, ysrc_f, keep);

    if (y0_f.size() >= 3) {
      std::vector<double> tmp2(z0_f.size(), 0.0);
      std::vector<double> sigma_z2(x_scan.size(), 0.0);
      for (size_t ix = 0; ix < x_scan.size(); ++ix) {
        double dx_det2 = x_scan[ix] - cfg.x_virtual;
        for (size_t i = 0; i < z0_f.size(); ++i) {
          tmp2[i] = z0_f[i] + dz_f[i] * dx_det2;
        }
        sigma_z2[ix] = weighted_std(tmp2, w_f);
      }

      size_t i_min2    = std::distance(sigma_z2.begin(),
                                       std::min_element(sigma_z2.begin(), sigma_z2.end(),
                                                        [](double a, double b) { return a < b; }));
      out.x_focus_scan = x_scan[i_min2];
      out.x_focus      = out.x_focus_scan;
      sigma_z_min      = sigma_z2[i_min2];

      if (i_min2 == 0 && x_scan.size() >= 3) {
        double x0  = x_scan[0];
        double x1  = x_scan[1];
        double x2  = x_scan[2];
        double y0q = sigma_z2[0] * sigma_z2[0];
        double y1q = sigma_z2[1] * sigma_z2[1];
        double y2q = sigma_z2[2] * sigma_z2[2];
        auto xv    = quadratic_vertex_from_three_points(x0, y0q, x1, y1q, x2, y2q);
        if (xv.has_value() && std::isfinite(*xv)) {
          out.x_focus = *xv;
        }
      }
      out.focus_before_lens2 = (out.x_focus < cfg.x2);

      double thr2 = k_threshold * sigma_z_min;
      int i_lo2   = static_cast<int>(i_min2);
      while (i_lo2 - 1 >= 0 && sigma_z2[static_cast<size_t>(i_lo2 - 1)] < thr2) {
        --i_lo2;
      }
      int i_hi2 = static_cast<int>(i_min2);
      while (i_hi2 + 1 < static_cast<int>(sigma_z2.size())
             && sigma_z2[static_cast<size_t>(i_hi2 + 1)] < thr2) {
        ++i_hi2;
      }
      double x_lo2 = x_scan[static_cast<size_t>(i_lo2)];
      double x_hi2 = x_scan[static_cast<size_t>(i_hi2)];
      out.dof      = x_hi2 - x_lo2;
    }
  }

  y_focus.assign(y0_f.size(), 0.0);
  z_focus.assign(z0_f.size(), 0.0);
  dx_focus = out.x_focus - cfg.x_virtual;
  for (size_t i = 0; i < y0_f.size(); ++i) {
    y_focus[i] = y0_f[i] + dy_f[i] * dx_focus;
    z_focus[i] = z0_f[i] + dz_f[i] * dx_focus;
  }

  double sigma_y_det = weighted_std(y_focus, w_f);
  out.M              = (sigma_y_src_theory > 0.0) ? (sigma_y_det / sigma_y_src_theory) : 0.0;
  out.M_abs_err      = std::abs(out.M - m_target);

  double p10              = weighted_percentile(y_focus, w_f, 0.10);
  double p90              = weighted_percentile(y_focus, w_f, 0.90);
  out.stripe_width        = p90 - p10;
  out.within_photocathode = (out.stripe_width <= 16.0) ? 1 : 0;

  double n_w_after = 0.0;
  for (double wi : w_f) {
    n_w_after += wi;
  }
  out.n_rays_core = (core_fraction < 1.0) ? n_w_after : n_w_before;

  return out;
}

int main(int argc, char** argv) {
  using json = nlohmann::json;

  CliConfig cli = parse_args(argc, argv);
  set_root_style();
  std::filesystem::create_directories(cli.output_dir);

  std::ifstream f(cli.config_file);
  if (!f.is_open()) {
    std::cerr << "Errore: impossibile aprire " << cli.config_file << "\n";
    return 1;
  }
  json config;
  f >> config;

  double x_min = config.value("x_min", 0.0);
  double x_max = config.value("x_max", 200.0);
  double dx    = config.value("dx", 1.0);

  double scan_min      = cli.scan_min.value_or(config.value("dof_x_scan_min", 100.0));
  double scan_max      = cli.scan_max.value_or(config.value("dof_x_scan_max", 350.0));
  double scan_step     = cli.scan_step.value_or(config.value("dof_x_scan_step", 0.5));
  double k_threshold   = cli.k_threshold.value_or(config.value("dof_k_threshold", std::sqrt(2.0)));
  double core_fraction = cli.core_fraction.value_or(config.value("dof_core_fraction", 1.0));
  double m_target      = cli.m_target.value_or(config.value("m_target", 1.0 / 7.5));
  double source_halfy       = config.value("dof_source_halfy", 5.0);
  double sigma_y_src_theory = source_halfy / std::sqrt(3.0);
  if (!(core_fraction > 0.0 && core_fraction <= 1.0)) {
    std::cerr << "Errore: core_fraction deve essere in (0, 1]\n";
    return 1;
  }
  if (!(scan_step > 0.0)) {
    std::cerr << "Errore: scan_step deve essere > 0\n";
    return 1;
  }
  if (!(scan_max > scan_min)) {
    std::cerr << "Errore: scan_max deve essere > scan_min\n";
    return 1;
  }

  TFile file(cli.input_file.c_str(), "READ");
  if (!file.IsOpen()) {
    std::cerr << "Errore: impossibile aprire " << cli.input_file << "\n";
    return 1;
  }

  TTree* tree_cfg  = (TTree*)file.Get("FocalConfigurations");
  TTree* tree_rays = (TTree*)file.Get("FocalRays");
  if (!tree_cfg || !tree_rays) {
    std::cerr << "Errore: TTree 'FocalConfigurations' o 'FocalRays' non trovato\n";
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
  for (int i = 0; i < tree_cfg->GetEntries(); ++i) {
    tree_cfg->GetEntry(i);
    config_map[config_id_cfg] = {config_id_cfg, x1_cfg, x2_cfg, x_virtual_cfg};
  }

  std::vector<double> x_scan;
  double x_virtual_min = std::numeric_limits<double>::infinity();
  double x_virtual_max = -std::numeric_limits<double>::infinity();
  for (const auto& [id, c] : config_map) {
    x_virtual_min = std::min(x_virtual_min, c.x_virtual);
    x_virtual_max = std::max(x_virtual_max, c.x_virtual);
  }
  if (std::isfinite(x_virtual_max)) {
    scan_min = std::max(scan_min, x_virtual_max);
  }

  for (double x = scan_min; x <= scan_max + 1e-12; x += scan_step) {
    x_scan.push_back(x);
  }
  if (x_scan.empty()) {
    std::cerr << "Errore: scan range vuoto (controlla scan_min/scan_max/scan_step)\n";
    return 1;
  }

  int config_id_rays              = 0;
  double n_rays                   = 0.0;
  std::vector<float>* y_hits_f    = nullptr;
  std::vector<float>* z_hits_f    = nullptr;
  std::vector<float>* dy_hits_f   = nullptr;
  std::vector<float>* dz_hits_f   = nullptr;
  std::vector<float>* w_hits_f    = nullptr;
  std::vector<float>* ysrc_hits_f = nullptr;

  tree_rays->SetBranchAddress("config_id", &config_id_rays);
  tree_rays->SetBranchAddress("n_rays", &n_rays);
  tree_rays->SetBranchAddress("y_hits", &y_hits_f);
  tree_rays->SetBranchAddress("z_hits", &z_hits_f);
  tree_rays->SetBranchAddress("dy_hits", &dy_hits_f);
  tree_rays->SetBranchAddress("dz_hits", &dz_hits_f);
  tree_rays->SetBranchAddress("weight_hits", &w_hits_f);
  tree_rays->SetBranchAddress("y_source_hits", &ysrc_hits_f);

  std::vector<ResultRow> results;
  results.reserve(static_cast<size_t>(tree_rays->GetEntries()));

  std::vector<double> y0, z0, dy, dz, w, ysrc;
  for (int i = 0; i < tree_rays->GetEntries(); ++i) {
    tree_rays->GetEntry(i);

    auto it = config_map.find(config_id_rays);
    if (it == config_map.end()) {
      continue;
    }

    if (!y_hits_f || !z_hits_f || !dy_hits_f || !dz_hits_f || !w_hits_f || !ysrc_hits_f) {
      std::cerr << "Errore: uno o piu' branch vettoriali non sono disponibili in FocalRays\n";
      return 1;
    }

    size_t n = y_hits_f ? y_hits_f->size() : 0;
    y0.resize(n);
    z0.resize(n);
    dy.resize(n);
    dz.resize(n);
    w.resize(n);
    ysrc.resize(n);
    for (size_t k = 0; k < n; ++k) {
      y0[k]   = (*y_hits_f)[k];
      z0[k]   = (*z_hits_f)[k];
      dy[k]   = (*dy_hits_f)[k];
      dz[k]   = (*dz_hits_f)[k];
      w[k]    = (*w_hits_f)[k];
      ysrc[k] = (*ysrc_hits_f)[k];
    }

    results.push_back(analyze_config(it->second, n_rays, y0, z0, dy, dz, w, ysrc, x_scan,
                                     k_threshold, core_fraction, m_target, sigma_y_src_theory));
  }

  if (!cli.tsv_out.empty()) {
    std::filesystem::path tsv_path = cli.tsv_out;
    if (tsv_path.has_parent_path()) {
      std::filesystem::create_directories(tsv_path.parent_path());
    }
    std::ofstream out(tsv_path);
    out << "x1\tx2\tx_focus\tx_focus_scan\tdof\tM\tm_target\tM_abs_err\tstripe_width\t"
           "within_photocathode\tn_rays\tn_rays_core\tcore_fraction\tconfig_id\tfocus_before_lens2\n";
    for (const auto& r : results) {
      out << r.x1 << "\t" << r.x2 << "\t" << r.x_focus << "\t" << r.x_focus_scan << "\t" << r.dof
          << "\t" << r.M << "\t" << m_target << "\t" << r.M_abs_err << "\t" << r.stripe_width
          << "\t" << r.within_photocathode << "\t" << r.n_rays << "\t" << r.n_rays_core << "\t"
          << core_fraction << "\t" << r.config_id << "\t" << (r.focus_before_lens2 ? 1 : 0)
          << "\n";
    }
  }

  double x1_min_data = std::numeric_limits<double>::infinity();
  double x1_max_data = -std::numeric_limits<double>::infinity();
  double x2_min_data = std::numeric_limits<double>::infinity();
  double x2_max_data = -std::numeric_limits<double>::infinity();
  for (const auto& [id, c] : config_map) {
    x1_min_data = std::min(x1_min_data, c.x1);
    x1_max_data = std::max(x1_max_data, c.x1);
    x2_min_data = std::min(x2_min_data, c.x2);
    x2_max_data = std::max(x2_max_data, c.x2);
  }

  int n_bins_x1 = std::round((x1_max_data - x1_min_data) / dx) + 1;
  int n_bins_x2 = std::round((x2_max_data - x2_min_data) / dx) + 1;

  TH2D h_focus("h_focus", ";x1 [mm];x2 [mm];x* [mm]", n_bins_x1, x1_min_data - dx / 2.0,
               x1_max_data + dx / 2.0, n_bins_x2, x2_min_data - dx / 2.0, x2_max_data + dx / 2.0);
  TH2D h_dof("h_dof", ";x1 [mm];x2 [mm];DoF [mm]", n_bins_x1, x1_min_data - dx / 2.0,
             x1_max_data + dx / 2.0, n_bins_x2, x2_min_data - dx / 2.0, x2_max_data + dx / 2.0);
  TH2D h_M("h_M", ";x1 [mm];x2 [mm];M", n_bins_x1, x1_min_data - dx / 2.0, x1_max_data + dx / 2.0,
           n_bins_x2, x2_min_data - dx / 2.0, x2_max_data + dx / 2.0);
  TH2D h_M_abs_err("h_M_abs_err", ";x1 [mm];x2 [mm];|M - M_{target}|", n_bins_x1,
                   x1_min_data - dx / 2.0, x1_max_data + dx / 2.0, n_bins_x2,
                   x2_min_data - dx / 2.0, x2_max_data + dx / 2.0);
  TH2D h_stripe("h_stripe", ";x1 [mm];x2 [mm];stripe width [mm]", n_bins_x1, x1_min_data - dx / 2.0,
                x1_max_data + dx / 2.0, n_bins_x2, x2_min_data - dx / 2.0, x2_max_data + dx / 2.0);
  TH2D h_margin("h_margin", ";x1 [mm];x2 [mm];16 - stripe width [mm]", n_bins_x1,
                x1_min_data - dx / 2.0, x1_max_data + dx / 2.0, n_bins_x2, x2_min_data - dx / 2.0,
                x2_max_data + dx / 2.0);

  TH2D h_mask("h_mask", "", n_bins_x1, x1_min_data - dx / 2.0, x1_max_data + dx / 2.0, n_bins_x2,
              x2_min_data - dx / 2.0, x2_max_data + dx / 2.0);
  for (int bx = 1; bx <= n_bins_x1; ++bx) {
    for (int by = 1; by <= n_bins_x2; ++by) {
      h_mask.SetBinContent(bx, by, 1.0);
    }
  }

  TH2D h_warn("h_warn", "", n_bins_x1, x1_min_data - dx / 2.0, x1_max_data + dx / 2.0, n_bins_x2,
              x2_min_data - dx / 2.0, x2_max_data + dx / 2.0);

  const double x_photocathode      = 180.0;
  const ResultRow* best_focus_near = nullptr;
  const ResultRow* best_dof_max    = nullptr;
  const ResultRow* best_M_neutral  = nullptr;
  const ResultRow* best_M_target   = nullptr;

  double focus_min  = std::numeric_limits<double>::infinity();
  double focus_max  = -std::numeric_limits<double>::infinity();
  double dof_min    = std::numeric_limits<double>::infinity();
  double dof_max    = -std::numeric_limits<double>::infinity();
  double M_min      = std::numeric_limits<double>::infinity();
  double M_max      = -std::numeric_limits<double>::infinity();
  double M_err_min  = std::numeric_limits<double>::infinity();
  double M_err_max  = -std::numeric_limits<double>::infinity();
  double stripe_min = std::numeric_limits<double>::infinity();
  double stripe_max = -std::numeric_limits<double>::infinity();
  double margin_min = std::numeric_limits<double>::infinity();
  double margin_max = -std::numeric_limits<double>::infinity();

  for (const auto& r : results) {
    int bin_x = h_focus.GetXaxis()->FindBin(r.x1);
    int bin_y = h_focus.GetYaxis()->FindBin(r.x2);

    h_mask.SetBinContent(bin_x, bin_y, 0.0);
    if (r.focus_before_lens2) {
      h_warn.SetBinContent(bin_x, bin_y, 1.0);
    }

    h_focus.SetBinContent(bin_x, bin_y, r.x_focus);
    h_dof.SetBinContent(bin_x, bin_y, r.dof);
    h_M.SetBinContent(bin_x, bin_y, r.M);
    h_M_abs_err.SetBinContent(bin_x, bin_y, r.M_abs_err);
    h_stripe.SetBinContent(bin_x, bin_y, r.stripe_width);
    double margin = 16.0 - r.stripe_width;
    h_margin.SetBinContent(bin_x, bin_y, margin);

    focus_min  = std::min(focus_min, r.x_focus);
    focus_max  = std::max(focus_max, r.x_focus);
    dof_min    = std::min(dof_min, r.dof);
    dof_max    = std::max(dof_max, r.dof);
    M_min      = std::min(M_min, r.M);
    M_max      = std::max(M_max, r.M);
    M_err_min  = std::min(M_err_min, r.M_abs_err);
    M_err_max  = std::max(M_err_max, r.M_abs_err);
    stripe_min = std::min(stripe_min, r.stripe_width);
    stripe_max = std::max(stripe_max, r.stripe_width);
    margin_min = std::min(margin_min, margin);
    margin_max = std::max(margin_max, margin);

    if (!best_focus_near
        || std::abs(r.x_focus - x_photocathode)
               < std::abs(best_focus_near->x_focus - x_photocathode)) {
      best_focus_near = &r;
    }
    if (!best_dof_max || r.dof > best_dof_max->dof) {
      best_dof_max = &r;
    }
    if (!best_M_neutral || std::abs(r.M - 1.0) < std::abs(best_M_neutral->M - 1.0)) {
      best_M_neutral = &r;
    }
    if (!best_M_target || r.M_abs_err < best_M_target->M_abs_err) {
      best_M_target = &r;
    }
  }

  if (std::isfinite(focus_min) && std::isfinite(focus_max)) {
    h_focus.SetMinimum(focus_min);
    h_focus.SetMaximum(focus_max);
  }
  if (std::isfinite(dof_min) && std::isfinite(dof_max)) {
    h_dof.SetMinimum(dof_min);
    h_dof.SetMaximum(dof_max);
  }
  if (std::isfinite(M_min) && std::isfinite(M_max)) {
    h_M.SetMinimum(M_min);
    h_M.SetMaximum(M_max);
  }
  if (std::isfinite(M_err_min) && std::isfinite(M_err_max)) {
    h_M_abs_err.SetMinimum(M_err_min);
    h_M_abs_err.SetMaximum(M_err_max);
  }
  if (std::isfinite(stripe_min) && std::isfinite(stripe_max)) {
    h_stripe.SetMinimum(stripe_min);
    h_stripe.SetMaximum(stripe_max);
  }
  if (std::isfinite(margin_min) && std::isfinite(margin_max)) {
    double max_abs = std::max(std::abs(margin_min), std::abs(margin_max));
    if (max_abs <= 0.0) {
      max_abs = 1.0;
    }
    h_margin.SetMinimum(-max_abs);
    h_margin.SetMaximum(+max_abs);
  }

  auto build_warn_boxes = [&]() {
    std::vector<std::unique_ptr<TBox>> boxes;
    boxes.reserve(static_cast<size_t>(n_bins_x1 * n_bins_x2));
    for (int bx = 1; bx <= n_bins_x1; ++bx) {
      double x_lo = h_warn.GetXaxis()->GetBinLowEdge(bx);
      double x_hi = h_warn.GetXaxis()->GetBinUpEdge(bx);
      for (int by = 1; by <= n_bins_x2; ++by) {
        if (h_warn.GetBinContent(bx, by) < 0.5) {
          continue;
        }
        double y_lo = h_warn.GetYaxis()->GetBinLowEdge(by);
        double y_hi = h_warn.GetYaxis()->GetBinUpEdge(by);
        auto box    = std::make_unique<TBox>(x_lo, y_lo, x_hi, y_hi);
        box->SetFillColorAlpha(kOrange + 1, 0.35);
        box->SetLineColor(0);
        box->Draw("same");
        boxes.push_back(std::move(box));
      }
    }
    return boxes;
  };

  auto build_mask_boxes = [&]() {
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
        box->SetFillColor(kWhite);
        box->SetLineColor(kWhite);
        box->Draw();
        boxes.push_back(std::move(box));
      }
    }
    return boxes;
  };

  static const std::unordered_map<std::string, std::string> kMapTitles = {
      {"dof_focus_map",        "Posizione di fuoco ottimale"},
      {"dof_dof_map",          "Profondit#grave{a} di campo (DoF)"},
      {"dof_M_map",            "Magnificazione M"},
      {"dof_stripe_map",       "Larghezza striscia"},
      {"dof_M_error_map",      "Errore magnificazione |M#minusM_{tgt}|"},
      {"dof_within_photocathode", "Margine fotocatodo"},
  };

  auto save_map = [&](TH2D& h, const std::string& name, int palette, const ResultRow* star) {
    gStyle->SetPalette(palette);
    TCanvas c(("c_" + name).c_str(), name.c_str(), 1100, 900);
    c.SetLeftMargin(0.10);
    c.SetBottomMargin(0.14);
    c.SetRightMargin(0.16);
    c.SetTopMargin(0.10);
    c.SetGridx();
    c.SetGridy();
    h.GetZaxis()->CenterTitle(kTRUE);
    h.GetZaxis()->SetTitleOffset(1.6);
    h.Draw("COLZ");
    c.Update();
    if (star) {
      TMarker m(star->x1, star->x2, 29);
      m.SetMarkerSize(2.0);
      m.SetMarkerColor(kBlack);
      m.Draw("same");
    }
    auto boxes      = build_mask_boxes();
    auto warn_boxes = build_warn_boxes();
    auto it         = kMapTitles.find(name);
    if (it != kMapTitles.end()) {
      TLatex tit;
      tit.SetNDC();
      tit.SetTextFont(42);
      tit.SetTextSize(0.042);
      tit.SetTextAlign(22);
      tit.DrawLatex(0.50, 0.955, it->second.c_str());
    }
    std::string out = (std::filesystem::path(cli.output_dir) / (name + ".png")).string();
    c.SaveAs(out.c_str());
  };

  save_map(h_focus, "dof_focus_map", kRainBow, best_focus_near);
  save_map(h_dof, "dof_dof_map", kBird, best_dof_max);
  save_map(h_M, "dof_M_map", kViridis, best_M_neutral);
  save_map(h_stripe, "dof_stripe_map", kViridis, nullptr);
  save_map(h_M_abs_err, "dof_M_error_map", kViridis, best_M_target);

  {
    const int nRGBs     = 3;
    double stops[nRGBs] = {0.0, 0.5, 1.0};
    double red[nRGBs]   = {1.0, 1.0, 0.0};
    double green[nRGBs] = {0.0, 1.0, 1.0};
    double blue[nRGBs]  = {0.0, 1.0, 0.0};
    TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, 255);
    gStyle->SetNumberContours(255);

    TCanvas c("c_within", "within", 1100, 900);
    c.SetLeftMargin(0.10);
    c.SetBottomMargin(0.14);
    c.SetRightMargin(0.16);
    c.SetTopMargin(0.10);
    c.SetGridx();
    c.SetGridy();
    h_margin.GetZaxis()->CenterTitle(kTRUE);
    h_margin.GetZaxis()->SetTitleOffset(1.6);
    h_margin.Draw("COLZ");
    c.Update();
    auto boxes      = build_mask_boxes();
    auto warn_boxes = build_warn_boxes();
    TLatex tit_w;
    tit_w.SetNDC();
    tit_w.SetTextFont(42);
    tit_w.SetTextSize(0.042);
    tit_w.SetTextAlign(22);
    tit_w.DrawLatex(0.50, 0.955, "Margine fotocatodo");
    std::string out =
        (std::filesystem::path(cli.output_dir) / "dof_within_photocathode.png").string();
    c.SaveAs(out.c_str());
  }

  return 0;
}
