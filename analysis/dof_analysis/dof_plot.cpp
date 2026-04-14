#include <nlohmann/json.hpp>

#include <TBox.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <string>
#include <utility>
#include <vector>

struct CliConfig {
  std::string input_file = "output/dof_simulation/focal.root";
  std::string output_dir = "output/dof_analysis";
  double x1              = std::numeric_limits<double>::quiet_NaN();
  double x2              = std::numeric_limits<double>::quiet_NaN();
  std::optional<double> scan_min;
  std::optional<double> scan_max;
  std::optional<double> scan_step;
  std::optional<double> k_threshold;
  double core_fraction = 1.0;
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
    else if (arg == "--output" || arg == "-o")
      cfg.output_dir = next();
    else if (arg == "--x1")
      cfg.x1 = std::stod(next());
    else if (arg == "--x2")
      cfg.x2 = std::stod(next());
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
    else if (arg == "--help" || arg == "-h") {
      std::cout << "Uso:\n"
                << "  dof_plot --input focal.root --x1 50.0 --x2 120.0 [--output dir/]\n"
                << "           [--scan-min 100] [--scan-max 350] [--scan-step 0.5] [--k 1.414]\n"
                << "           [--core-fraction 1.0]\n";
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

static std::vector<size_t> select_core_indices_around_mean(const std::vector<double>& x,
                                                           const std::vector<double>& w,
                                                           double core_fraction) {
  if (x.empty()) {
    return {};
  }
  if (core_fraction >= 1.0) {
    std::vector<size_t> idx(x.size());
    std::iota(idx.begin(), idx.end(), 0);
    return idx;
  }
  if (core_fraction <= 0.0) {
    return {};
  }

  double mu = weighted_mean(x, w);
  if (!std::isfinite(mu)) {
    return {};
  }

  double sum_w = 0.0;
  for (double wi : w) {
    sum_w += wi;
  }
  if (sum_w <= 0.0) {
    return {};
  }

  std::vector<size_t> idx(x.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(),
            [&](size_t a, size_t b) { return std::abs(x[a] - mu) < std::abs(x[b] - mu); });

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

struct FocusEllipse {
  double mean_y    = 0.0;
  double mean_z    = 0.0;
  double r1        = 0.0;
  double r2        = 0.0;
  double angle_deg = 0.0;
};

static FocusEllipse weighted_covariance_ellipse(const std::vector<double>& y,
                                                const std::vector<double>& z,
                                                const std::vector<double>& w) {
  FocusEllipse out;
  if (y.empty()) {
    return out;
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
    return out;
  }

  out.mean_y = sum_y / sum_w;
  out.mean_z = sum_z / sum_w;

  double vyy = 0.0;
  double vzz = 0.0;
  double vyz = 0.0;
  for (size_t i = 0; i < y.size(); ++i) {
    double dy = y[i] - out.mean_y;
    double dz = z[i] - out.mean_z;
    vyy += w[i] * dy * dy;
    vzz += w[i] * dz * dz;
    vyz += w[i] * dy * dz;
  }
  vyy /= sum_w;
  vzz /= sum_w;
  vyz /= sum_w;

  double tr   = vyy + vzz;
  double det  = vyy * vzz - vyz * vyz;
  double disc = tr * tr - 4.0 * det;
  disc        = std::max(0.0, disc);

  double l1 = 0.5 * (tr + std::sqrt(disc));
  double l2 = 0.5 * (tr - std::sqrt(disc));

  out.r1 = std::sqrt(std::max(0.0, l1));
  out.r2 = std::sqrt(std::max(0.0, l2));

  out.angle_deg = 0.5 * std::atan2(2.0 * vyz, vyy - vzz) * 180.0 / M_PI;
  return out;
}

int main(int argc, char** argv) {
  using json = nlohmann::json;

  CliConfig cli = parse_args(argc, argv);
  if (!std::isfinite(cli.x1) || !std::isfinite(cli.x2)) {
    std::cerr << "Errore: --x1 e --x2 sono obbligatori\n";
    return 1;
  }
  if (!(cli.core_fraction > 0.0 && cli.core_fraction <= 1.0)) {
    std::cerr << "Errore: --core-fraction deve essere in (0, 1]\n";
    return 1;
  }

  set_root_style();
  std::filesystem::create_directories(cli.output_dir);

  json cfg_json;
  {
    std::ifstream f_cfg("config/config.json");
    if (f_cfg.is_open()) {
      f_cfg >> cfg_json;
    }
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
  char lens75_id_buf[256];
  char lens60_id_buf[256];

  tree_cfg->SetBranchAddress("config_id", &config_id_cfg);
  tree_cfg->SetBranchAddress("x1", &x1_cfg);
  tree_cfg->SetBranchAddress("x2", &x2_cfg);
  tree_cfg->SetBranchAddress("x_virtual", &x_virtual_cfg);
  tree_cfg->SetBranchAddress("lens75_id", lens75_id_buf);
  tree_cfg->SetBranchAddress("lens60_id", lens60_id_buf);

  int best_config_id    = -1;
  double best_x1        = 0.0;
  double best_x2        = 0.0;
  double best_x_virtual = 0.0;
  std::string best_l75;
  std::string best_l60;
  double best_d2 = std::numeric_limits<double>::infinity();

  for (int i = 0; i < tree_cfg->GetEntries(); ++i) {
    tree_cfg->GetEntry(i);
    double d2 = (x1_cfg - cli.x1) * (x1_cfg - cli.x1) + (x2_cfg - cli.x2) * (x2_cfg - cli.x2);
    if (d2 < best_d2) {
      best_d2        = d2;
      best_config_id = config_id_cfg;
      best_x1        = x1_cfg;
      best_x2        = x2_cfg;
      best_x_virtual = x_virtual_cfg;
      best_l75       = lens75_id_buf;
      best_l60       = lens60_id_buf;
    }
  }

  if (best_config_id < 0) {
    std::cerr << "Errore: nessuna configurazione trovata nel file\n";
    return 1;
  }

  double scan_step    = cli.scan_step.value_or(cfg_json.value("dof_x_scan_step", 0.5));
  double scan_max     = cli.scan_max.value_or(cfg_json.value("dof_x_scan_max", 400.0));
  double scan_min_cfg = cfg_json.value("dof_x_scan_min", best_x_virtual);
  double scan_min     = cli.scan_min.value_or(scan_min_cfg);
  double k_threshold  = cli.k_threshold.value_or(cfg_json.value("dof_k_threshold", std::sqrt(2.0)));

  if (scan_min < best_x_virtual) {
    scan_min = best_x_virtual;
  }
  if (!(scan_step > 0.0) || !(scan_max > scan_min)) {
    std::cerr << "Errore: scan range non valido\n";
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

  bool found_rays = false;
  std::vector<double> y0, z0, dy, dz, w, ysrc;
  for (int i = 0; i < tree_rays->GetEntries(); ++i) {
    tree_rays->GetEntry(i);
    if (config_id_rays != best_config_id) {
      continue;
    }
    found_rays = true;
    size_t n   = y_hits_f ? y_hits_f->size() : 0;
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
    break;
  }

  if (!found_rays) {
    std::cerr << "Errore: rays non trovati per config_id=" << best_config_id << "\n";
    return 1;
  }

  std::vector<double> x_scan;
  for (double x = scan_min; x <= scan_max + 1e-12; x += scan_step) {
    x_scan.push_back(x);
  }

  std::vector<double> sigma_z(x_scan.size(), 0.0);
  std::vector<double> z_tmp(y0.size(), 0.0);
  for (size_t ix = 0; ix < x_scan.size(); ++ix) {
    double x_det  = x_scan[ix];
    double dx_det = x_det - best_x_virtual;
    for (size_t i = 0; i < z0.size(); ++i) {
      z_tmp[i] = z0[i] + dz[i] * dx_det;
    }
    sigma_z[ix] = weighted_std(z_tmp, w);
  }

  size_t i_min =
      std::distance(sigma_z.begin(), std::min_element(sigma_z.begin(), sigma_z.end(),
                                                      [](double a, double b) { return a < b; }));
  double x_focus     = x_scan[i_min];
  double sigma_z_min = sigma_z[i_min];

  std::vector<double> y_focus(y0.size(), 0.0);
  std::vector<double> z_focus(z0.size(), 0.0);
  double dx_focus = x_focus - best_x_virtual;
  for (size_t i = 0; i < y0.size(); ++i) {
    y_focus[i] = y0[i] + dy[i] * dx_focus;
    z_focus[i] = z0[i] + dz[i] * dx_focus;
  }

  double n_rays_before = 0.0;
  for (double wi : w) {
    n_rays_before += wi;
  }

  if (cli.core_fraction < 1.0) {
    auto keep = select_core_indices_around_mean(y_focus, w, cli.core_fraction);
    apply_selection(y0, z0, dy, dz, w, ysrc, keep);

    x_scan.clear();
    for (double x = scan_min; x <= scan_max + 1e-12; x += scan_step) {
      x_scan.push_back(x);
    }

    sigma_z.assign(x_scan.size(), 0.0);
    z_tmp.assign(y0.size(), 0.0);
    for (size_t ix = 0; ix < x_scan.size(); ++ix) {
      double x_det  = x_scan[ix];
      double dx_det = x_det - best_x_virtual;
      for (size_t i = 0; i < z0.size(); ++i) {
        z_tmp[i] = z0[i] + dz[i] * dx_det;
      }
      sigma_z[ix] = weighted_std(z_tmp, w);
    }

    i_min =
        std::distance(sigma_z.begin(), std::min_element(sigma_z.begin(), sigma_z.end(),
                                                        [](double a, double b) { return a < b; }));
    x_focus     = x_scan[i_min];
    sigma_z_min = sigma_z[i_min];

    y_focus.assign(y0.size(), 0.0);
    z_focus.assign(z0.size(), 0.0);
    dx_focus = x_focus - best_x_virtual;
    for (size_t i = 0; i < y0.size(); ++i) {
      y_focus[i] = y0[i] + dy[i] * dx_focus;
      z_focus[i] = z0[i] + dz[i] * dx_focus;
    }
  }

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
  double dof  = x_hi - x_lo;

  double sigma_y_src = weighted_std(ysrc, w);
  double sigma_y_det = weighted_std(y_focus, w);
  double M           = (sigma_y_src > 0.0) ? (sigma_y_det / sigma_y_src) : 0.0;

  double p10               = weighted_percentile(y_focus, w, 0.10);
  double p90               = weighted_percentile(y_focus, w, 0.90);
  double stripe_width      = p90 - p10;
  bool within_photocathode = (stripe_width <= 16.0);

  std::cout << "Focal configuration (nearest):\n";
  std::cout << "  config_id = " << best_config_id << "\n";
  std::cout << "  x1 = " << best_x1 << " mm, x2 = " << best_x2 << " mm\n";
  std::cout << "  x_virtual = " << best_x_virtual << " mm\n";
  std::cout << "  lenses = " << best_l75 << " & " << best_l60 << "\n";
  std::cout << "  n_rays (weighted) = " << n_rays << "\n\n";

  std::cout << "Results:\n";
  std::cout << "  x* = " << x_focus << " mm\n";
  std::cout << "  sigma_z(x*) = " << sigma_z_min << " mm\n";
  std::cout << "  DoF(k=" << k_threshold << ") = " << dof << " mm  (" << x_lo << " .. " << x_hi
            << ")\n";
  std::cout << "  sigma_y_src = " << sigma_y_src << " mm\n";
  std::cout << "  sigma_y_det(x*) = " << sigma_y_det << " mm\n";
  std::cout << "  M(x*) = " << M << "\n";
  std::cout << "  stripe_width(P90-P10) = " << stripe_width << " mm\n";
  std::cout << "  within_photocathode(16mm) = " << (within_photocathode ? 1 : 0) << "\n";
  if (cli.core_fraction < 1.0) {
    double n_rays_after = 0.0;
    for (double wi : w) {
      n_rays_after += wi;
    }
    std::cout << "  core_fraction = " << cli.core_fraction << "\n";
    std::cout << "  n_rays_core (weighted) = " << n_rays_after << " (from " << n_rays_before
              << ")\n";
  }

  double sigma_z_max = *std::max_element(sigma_z.begin(), sigma_z.end());

  {
    TCanvas c("c_sigma", "sigma_z", 1100, 750);
    c.SetLeftMargin(0.16);
    c.SetBottomMargin(0.14);
    c.SetRightMargin(0.10);
    c.SetTopMargin(0.08);

    TGraph g(static_cast<int>(x_scan.size()));
    for (int i = 0; i < static_cast<int>(x_scan.size()); ++i) {
      g.SetPoint(i, x_scan[static_cast<size_t>(i)], sigma_z[static_cast<size_t>(i)]);
    }
    g.SetTitle(";x_{det} [mm];#sigma_{z}(x_{det}) [mm]");
    g.SetLineWidth(2);
    g.Draw("AL");

    TBox box(x_lo, 0.0, x_hi, sigma_z_max);
    box.SetFillColorAlpha(kGray + 1, 0.20);
    box.SetLineColor(0);
    box.Draw("same");

    TLine line(x_focus, 0.0, x_focus, sigma_z_max);
    line.SetLineColor(kRed + 1);
    line.SetLineWidth(2);
    line.Draw("same");

    TPaveText pt(0.18, 0.78, 0.62, 0.92, "NDC");
    pt.SetFillColor(0);
    pt.SetTextFont(42);
    pt.SetTextAlign(12);
    pt.AddText(Form("x* = %.2f mm", x_focus));
    pt.AddText(Form("#sigma_{z}(x*) = %.3f mm", sigma_z_min));
    pt.AddText(Form("DoF = %.1f mm", dof));
    pt.Draw("same");

    std::string out = (std::filesystem::path(cli.output_dir)
                       / ("sigma_z_config_" + std::to_string(best_config_id) + ".png"))
                          .string();
    c.SaveAs(out.c_str());
  }

  {
    double y_abs_max = 0.0;
    double z_abs_max = 0.0;
    for (size_t i = 0; i < y_focus.size(); ++i) {
      y_abs_max = std::max(y_abs_max, std::abs(y_focus[i]));
      z_abs_max = std::max(z_abs_max, std::abs(z_focus[i]));
    }
    double y_range = std::max(12.0, 1.2 * y_abs_max);
    double z_range = std::max(12.0, 1.2 * z_abs_max);

    int nbins = 200;
    TH2D h("h_focus", ";y [mm];z [mm]", nbins, -y_range, y_range, nbins, -z_range, z_range);
    for (size_t i = 0; i < y_focus.size(); ++i) {
      h.Fill(y_focus[i], z_focus[i], w[i]);
    }

    TCanvas c("c_focus", "focus", 900, 850);
    c.SetLeftMargin(0.16);
    c.SetBottomMargin(0.14);
    c.SetRightMargin(0.16);
    c.SetTopMargin(0.08);
    h.Draw("COLZ");

    FocusEllipse el = weighted_covariance_ellipse(y_focus, z_focus, w);
    TEllipse e(el.mean_y, el.mean_z, el.r1, el.r2, 0, 360, el.angle_deg);
    e.SetLineColor(kRed + 1);
    e.SetLineWidth(2);
    e.SetFillStyle(0);
    e.Draw("same");

    TBox photocathode(-8.0, -8.0, 8.0, 8.0);
    photocathode.SetFillStyle(0);
    photocathode.SetLineColor(kGreen + 2);
    photocathode.SetLineWidth(3);
    photocathode.Draw("same");

    TPaveText pt(0.18, 0.80, 0.62, 0.92, "NDC");
    pt.SetFillColor(0);
    pt.SetTextFont(42);
    pt.SetTextAlign(12);
    pt.AddText(Form("M = %.3f", M));
    pt.AddText(Form("stripe = %.2f mm", stripe_width));
    pt.Draw("same");

    std::string out = (std::filesystem::path(cli.output_dir)
                       / ("focus_yz_config_" + std::to_string(best_config_id) + ".png"))
                          .string();
    c.SaveAs(out.c_str());
  }

  {
    double min_v = std::numeric_limits<double>::infinity();
    double max_v = -std::numeric_limits<double>::infinity();
    for (double v : ysrc) {
      min_v = std::min(min_v, v);
      max_v = std::max(max_v, v);
    }
    for (double v : y_focus) {
      min_v = std::min(min_v, v);
      max_v = std::max(max_v, v);
    }
    double pad = 0.10 * (max_v - min_v);
    if (pad <= 0.0) {
      pad = 1.0;
    }

    TH1D h_src("h_src", ";y [mm];norm.", 120, min_v - pad, max_v + pad);
    TH1D h_det("h_det", ";y [mm];norm.", 120, min_v - pad, max_v + pad);
    for (size_t i = 0; i < y_focus.size(); ++i) {
      h_src.Fill(ysrc[i], w[i]);
      h_det.Fill(y_focus[i], w[i]);
    }
    if (h_src.Integral() > 0.0) {
      h_src.Scale(1.0 / h_src.Integral());
    }
    if (h_det.Integral() > 0.0) {
      h_det.Scale(1.0 / h_det.Integral());
    }

    h_src.SetLineColor(kBlue + 1);
    h_src.SetLineWidth(2);
    h_det.SetLineColor(kRed + 1);
    h_det.SetLineWidth(2);

    TCanvas c("c_ycmp", "y_compare", 1100, 750);
    c.SetLeftMargin(0.16);
    c.SetBottomMargin(0.14);
    c.SetRightMargin(0.10);
    c.SetTopMargin(0.08);
    h_src.SetMaximum(1.25 * std::max(h_src.GetMaximum(), h_det.GetMaximum()));
    h_src.Draw("HIST");
    h_det.Draw("HIST SAME");

    TLegend leg(0.62, 0.78, 0.88, 0.90);
    leg.SetTextFont(42);
    leg.SetBorderSize(0);
    leg.SetFillColor(0);
    leg.AddEntry(&h_src, "y_{src}", "l");
    leg.AddEntry(&h_det, "y_{det}(x*)", "l");
    leg.Draw("same");

    std::string out = (std::filesystem::path(cli.output_dir)
                       / ("y_compare_config_" + std::to_string(best_config_id) + ".png"))
                          .string();
    c.SaveAs(out.c_str());
  }

  {
    double sigma_y_src_local = weighted_std(ysrc, w);
    double M_focus_local     = (sigma_y_src_local > 0.0) ? (sigma_y_det / sigma_y_src_local) : 0.0;

    std::vector<double> sigma_y_over_M(x_scan.size(), 0.0);
    std::vector<double> y_tmp(y0.size(), 0.0);
    for (size_t ix = 0; ix < x_scan.size(); ++ix) {
      double x_det  = x_scan[ix];
      double dx_det = x_det - best_x_virtual;
      for (size_t i = 0; i < y0.size(); ++i) {
        y_tmp[i] = y0[i] + dy[i] * dx_det;
      }
      double sy          = weighted_std(y_tmp, w);
      sigma_y_over_M[ix] = (M_focus_local > 0.0) ? (sy / M_focus_local) : 0.0;
    }

    double y_right_max = *std::max_element(sigma_y_over_M.begin(), sigma_y_over_M.end());
    if (y_right_max <= 0.0) {
      y_right_max = 1.0;
    }
    double scale = sigma_z_max / y_right_max;

    TCanvas c("c_combo", "sigma_combo", 1100, 750);
    c.SetLeftMargin(0.16);
    c.SetBottomMargin(0.14);
    c.SetRightMargin(0.16);
    c.SetTopMargin(0.08);

    TGraph g1(static_cast<int>(x_scan.size()));
    TGraph g2(static_cast<int>(x_scan.size()));
    for (int i = 0; i < static_cast<int>(x_scan.size()); ++i) {
      double x_det = x_scan[static_cast<size_t>(i)];
      g1.SetPoint(i, x_det, sigma_z[static_cast<size_t>(i)]);
      g2.SetPoint(i, x_det, sigma_y_over_M[static_cast<size_t>(i)] * scale);
    }

    g1.SetTitle(";x_{det} [mm];#sigma_{z}(x_{det}) [mm]");
    g1.SetLineColor(kBlue + 1);
    g1.SetLineWidth(2);
    g1.Draw("AL");

    g2.SetLineColor(kRed + 1);
    g2.SetLineStyle(7);
    g2.SetLineWidth(2);
    g2.Draw("L SAME");

    TLine line(x_focus, 0.0, x_focus, sigma_z_max);
    line.SetLineColor(kBlack);
    line.SetLineWidth(2);
    line.Draw("same");

    TGaxis axis(scan_max, 0.0, scan_max, sigma_z_max, 0.0, y_right_max, 510, "+L");
    axis.SetLineColor(kRed + 1);
    axis.SetLabelColor(kRed + 1);
    axis.SetTitleColor(kRed + 1);
    axis.SetTitle("#sigma_{y}(x_{det})/M(x*) [mm]");
    axis.SetTitleFont(42);
    axis.SetLabelFont(42);
    axis.Draw();

    TLegend leg(0.62, 0.78, 0.88, 0.90);
    leg.SetTextFont(42);
    leg.SetBorderSize(0);
    leg.SetFillColor(0);
    leg.AddEntry(&g1, "#sigma_{z}(x_{det})", "l");
    leg.AddEntry(&g2, "#sigma_{y}(x_{det})/M (scaled)", "l");
    leg.Draw("same");

    std::string out = (std::filesystem::path(cli.output_dir)
                       / ("sigma_combo_config_" + std::to_string(best_config_id) + ".png"))
                          .string();
    c.SaveAs(out.c_str());
  }

  return 0;
}
