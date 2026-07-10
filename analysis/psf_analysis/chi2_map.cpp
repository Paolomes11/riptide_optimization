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
 * chi2_map — Mappa 2D del chi-quadro ottenuto dal fit di un piano
 *            alle posizioni medie sul detector (mu_y, mu_z) in funzione
 *            della posizione sorgente (x_src, y_src).
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

#include <omp.h>

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

// Parsing CLI
struct CliConfig {
  std::string psf_path    = "output/psf/psf_data.root";
  std::string config_path = "config/config.json";
  std::string output_path = "output/psf_analysis/chi2_map.png";
  std::string tsv_path    = "";

  double min_hits  = 10.0;
  bool log_scale   = false;
  bool use_reduced = true;
  bool dist_to_n   = false;
  double dist_n    = 1.0;

  bool corr_map        = false;
  bool adaptive_target = false;

  bool apply_corr = true;
  double grid_dx  = -1.0;
  double grid_dy  = -1.0;

  // Parametri percentili per la scala colori
  double perc_low  = 0.0;
  double perc_high = 95.0;

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
    else if (key == "--min-hits")
      cfg.min_hits = std::stod(next());
    else if (key == "--log")
      cfg.log_scale = true;
    else if (key == "--dist-to-n") {
      cfg.dist_to_n = true;
      cfg.dist_n    = std::stod(next());
    } else if (key == "--dist-to-one") {
      cfg.dist_to_n = true;
      cfg.dist_n    = 1.0;
    } else if (key == "--no-reduced")
      cfg.use_reduced = false;
    else if (key == "--corr-map")
      cfg.corr_map = true;
    else if (key == "--adaptive-target")
      cfg.adaptive_target = true;
    else if (key == "--no-corr")
      cfg.apply_corr = false;
    else if (key == "--grid-dx")
      cfg.grid_dx = std::stod(next());
    else if (key == "--grid-dy")
      cfg.grid_dy = std::stod(next());
    else if (key == "--p-low")
      cfg.perc_low = std::stod(next());
    else if (key == "--p-high")
      cfg.perc_high = std::stod(next());
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

struct PlaneFitResult {
  double a = 0.0, b = 0.0, c = 0.0;
  double chi2      = 0.0;
  int ndof         = 0;
  double chi2_ndof = 0.0;
  bool valid       = false;

  std::vector<double> residuals;
  std::vector<double> residual_sig;
};

static bool solve_plane_wls(const std::vector<double>& x_src, const std::vector<double>& y_src,
                            const std::vector<double>& mu, const std::vector<double>& sigma,
                            PlaneFitResult& res) {
  const int N = static_cast<int>(x_src.size());
  if (N < 4 || static_cast<int>(y_src.size()) != N || static_cast<int>(mu.size()) != N
      || static_cast<int>(sigma.size()) != N) {
    res.valid = false;
    return false;
  }

  double M[3][3] = {};
  double rhs[3]  = {};

  for (int i = 0; i < N; ++i) {
    const double sig2 = std::max(sigma[i] * sigma[i], 1e-30);
    const double w    = 1.0 / sig2;
    const double yi   = y_src[i];
    const double xi   = x_src[i];
    const double mi   = mu[i];

    const double row[3] = {1.0, yi, xi};
    for (int r = 0; r < 3; ++r) {
      rhs[r] += w * row[r] * mi;
      for (int c = 0; c < 3; ++c)
        M[r][c] += w * row[r] * row[c];
    }
  }

  const double det = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
                   - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                   + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);

  const double scale =
      std::abs(M[0][0]) * std::abs(M[1][1]) * std::abs(M[2][2]);
  if (!std::isfinite(det) || std::abs(det) < 1e-12 * std::max(scale, 1.0)) {
    res.valid = false;
    return false;
  }

  double inv[3][3];
  inv[0][0] = (M[1][1] * M[2][2] - M[1][2] * M[2][1]) / det;
  inv[0][1] = -(M[0][1] * M[2][2] - M[0][2] * M[2][1]) / det;
  inv[0][2] = (M[0][1] * M[1][2] - M[0][2] * M[1][1]) / det;
  inv[1][0] = -(M[1][0] * M[2][2] - M[1][2] * M[2][0]) / det;
  inv[1][1] = (M[0][0] * M[2][2] - M[0][2] * M[2][0]) / det;
  inv[1][2] = -(M[0][0] * M[1][2] - M[0][2] * M[1][0]) / det;
  inv[2][0] = (M[1][0] * M[2][1] - M[1][1] * M[2][0]) / det;
  inv[2][1] = -(M[0][0] * M[2][1] - M[0][1] * M[2][0]) / det;
  inv[2][2] = (M[0][0] * M[1][1] - M[0][1] * M[1][0]) / det;

  double theta[3] = {};
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c)
      theta[r] += inv[r][c] * rhs[c];

  if (!std::isfinite(theta[0]) || !std::isfinite(theta[1]) || !std::isfinite(theta[2])) {
    res.valid = false;
    return false;
  }

  res.a = theta[0];
  res.b = theta[1];
  res.c = theta[2];

  res.chi2 = 0.0;
  res.residuals.resize(static_cast<size_t>(N));
  res.residual_sig.resize(static_cast<size_t>(N));
  for (int i = 0; i < N; ++i) {
    const double mu_hat = res.a + res.b * y_src[i] + res.c * x_src[i];
    const double r_i    = mu[i] - mu_hat;
    const double sig2   = std::max(sigma[i] * sigma[i], 1e-30);
    res.chi2 += r_i * r_i / sig2;
    res.residuals[static_cast<size_t>(i)]    = r_i;
    res.residual_sig[static_cast<size_t>(i)] = sigma[i];
  }

  res.ndof      = N - 3;
  res.chi2_ndof = (res.ndof > 0) ? res.chi2 / static_cast<double>(res.ndof) : 0.0;
  res.valid     = std::isfinite(res.chi2_ndof);
  return res.valid;
}

static double estimate_min_step(std::vector<double> v) {
  v.erase(std::remove_if(v.begin(), v.end(), [](double x) { return !std::isfinite(x); }), v.end());
  if (v.size() < 2)
    return 0.0;
  std::sort(v.begin(), v.end());
  v.erase(std::unique(v.begin(), v.end()), v.end());
  double step = std::numeric_limits<double>::infinity();
  for (size_t i = 1; i < v.size(); ++i) {
    const double d = v[i] - v[i - 1];
    if (d > 0.0 && d < step)
      step = d;
  }
  return std::isfinite(step) ? step : 0.0;
}

static double estimate_residual_correlation(const std::vector<double>& x_src,
                                            const std::vector<double>& y_src,
                                            const std::vector<double>& residuals, double dx_thresh,
                                            double dy_thresh, int& n_pairs_out) {
  n_pairs_out = 0;
  if (dx_thresh <= 0.0 || dy_thresh <= 0.0)
    return 0.0;

  const int N = static_cast<int>(residuals.size());
  if (N < 2)
    return 0.0;

  double mean_r = 0.0;
  for (int i = 0; i < N; ++i)
    mean_r += residuals[static_cast<size_t>(i)];
  mean_r /= static_cast<double>(N);

  double num = 0.0, den = 0.0;

  for (int i = 0; i < N; ++i) {
    for (int j = i + 1; j < N; ++j) {
      const double dx = std::abs(x_src[static_cast<size_t>(i)] - x_src[static_cast<size_t>(j)]);
      const double dy = std::abs(y_src[static_cast<size_t>(i)] - y_src[static_cast<size_t>(j)]);
      if (dx <= dx_thresh && dy <= dy_thresh) {
        const double ri = residuals[static_cast<size_t>(i)] - mean_r;
        const double rj = residuals[static_cast<size_t>(j)] - mean_r;
        num += ri * rj;
        den += 0.5 * (ri * ri + rj * rj);
        ++n_pairs_out;
      }
    }
  }

  if (n_pairs_out == 0 || den < 1e-30)
    return 0.0;
  return std::clamp(num / den, 0.0, 0.999);
}

static double compute_inflation_factor(double rho, int N, double n_neighbors_avg) {
  if (N <= 3)
    return 1.0;
  const double frac = n_neighbors_avg / static_cast<double>(N);
  return 1.0 + 2.0 * rho * frac * static_cast<double>(N - 3) / static_cast<double>(N - 1);
}

static int count_distinct(std::vector<double> v, double eps = 1e-3) {
  v.erase(std::remove_if(v.begin(), v.end(), [](double x) { return !std::isfinite(x); }), v.end());
  if (v.empty())
    return 0;
  std::sort(v.begin(), v.end());
  int n = 1;
  for (size_t i = 1; i < v.size(); ++i) {
    if (v[i] - v[i - 1] > eps)
      ++n;
  }
  return n;
}

static double compute_chi2_target(double rho, int N_x, int N_y) {
  if (rho <= 0.0 || (N_x <= 2 && N_y <= 2))
    return 1.0;

  rho = std::clamp(rho, 0.0, 0.999999);

  const bool use_x = N_x > 2;
  const bool use_y = N_y > 2;

  if (use_x && !use_y)
    return riptide::expected_chi2_ndof_ar1(N_x, rho);
  if (!use_x && use_y)
    return riptide::expected_chi2_ndof_ar1(N_y, rho);

  const double tx = riptide::expected_chi2_ndof_ar1(N_x, rho);
  const double ty = riptide::expected_chi2_ndof_ar1(N_y, rho);
  const double wx = static_cast<double>(N_x) / static_cast<double>(N_x + N_y);
  const double wy = 1.0 - wx;
  const double t  = wx * tx + wy * ty;
  return std::isfinite(t) ? t : 1.0;
}

int main(int argc, char** argv) {
  CliConfig cli = parse_args(argc, argv);
  if (cli.corr_map && cli.output_path == "output/psf_analysis/chi2_map.png")
    cli.output_path = "output/psf_analysis/corr_map.png";
  if (cli.adaptive_target && !cli.corr_map && cli.output_path == "output/psf_analysis/chi2_map.png")
    cli.output_path = "output/psf_analysis/chi2_map_adaptive.png";

  using json = nlohmann::json;
  std::ifstream f_cfg(cli.config_path);
  if (!f_cfg.is_open()) {
    std::cerr << "Errore: impossibile aprire " << cli.config_path << "\n";
    return 1;
  }
  json jcfg;
  f_cfg >> jcfg;
  const double dx = jcfg["dx"];

  riptide::PSFDatabase db;
  try {
    db = riptide::load_psf_database(cli.psf_path);
  } catch (const std::exception& e) {
    std::cerr << "Errore: " << e.what() << "\n";
    return 1;
  }
  if (db.empty())
    return 1;

  double x1_lo = 1e18, x1_hi = -1e18, x2_lo = 1e18, x2_hi = -1e18;
  for (const auto& [cfg, _] : db) {
    x1_lo = std::min(x1_lo, cfg.x1);
    x1_hi = std::max(x1_hi, cfg.x1);
    x2_lo = std::min(x2_lo, cfg.x2);
    x2_hi = std::max(x2_hi, cfg.x2);
  }
  int bins_x      = std::max(1, static_cast<int>(std::round((x1_hi - x1_lo) / dx)) + 1);
  int bins_y      = std::max(1, static_cast<int>(std::round((x2_hi - x2_lo) / dx)) + 1);
  double hx       = dx / 2.0;
  double ax_x1_lo = x1_lo - hx, ax_x1_hi = x1_hi + hx;
  double ax_x2_lo = x2_lo - hx, ax_x2_hi = x2_hi + hx;

  struct Chi2Entry {
    double x1, x2;
    double metric;
    double chi2_red;
    double chi2_red_raw;
    double rho;
    double chi2_target;
    double infl;
    bool valid;
  };
  if (cli.n_jobs > 0)
    omp_set_num_threads(cli.n_jobs);

  std::vector<std::pair<riptide::LensConfig, const std::vector<riptide::PSFPoint>*>> db_vec;
  db_vec.reserve(db.size());
  for (const auto& [k, v] : db)
    db_vec.push_back({k, &v});

  int total_cfgs = static_cast<int>(db_vec.size());
  int print_step = std::max(1, total_cfgs / 20);

  std::vector<Chi2Entry> entries(total_cfgs, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, false});
  std::atomic<int> processed{0};
  std::mutex log_mtx;

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < total_cfgs; ++i) {
    const auto& cfg    = db_vec[i].first;
    const auto& points = *db_vec[i].second;

    std::vector<const riptide::PSFPoint*> valid_points;
    for (const auto& p : points) {
      if (p.on_detector && p.n_hits_count >= cli.min_hits)
        valid_points.push_back(&p);
    }

    if (valid_points.size() < 10) {
      entries[i] = {cfg.x1, cfg.x2, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, false};
    } else {
      std::vector<double> xs, ys, mu_y_vec, mu_z_vec, sig_y_vec, sig_z_vec;
      xs.reserve(valid_points.size());
      ys.reserve(valid_points.size());
      mu_y_vec.reserve(valid_points.size());
      mu_z_vec.reserve(valid_points.size());
      sig_y_vec.reserve(valid_points.size());
      sig_z_vec.reserve(valid_points.size());

      for (const auto* p : valid_points) {
        const double n_count = std::max(1.0, p->n_hits_count);
        const double err_y   = std::sqrt(std::max(0.0, p->cov_yy) / n_count);
        const double err_z   = std::sqrt(std::max(0.0, p->cov_zz) / n_count);

        xs.push_back(p->x_source);
        ys.push_back(p->y_source);
        mu_y_vec.push_back(p->mu_y);
        mu_z_vec.push_back(p->mu_z);
        sig_y_vec.push_back(std::max(err_y, 1e-6));
        sig_z_vec.push_back(std::max(err_z, 1e-6));
      }

      PlaneFitResult fit_y, fit_z;
      const bool ok_y = solve_plane_wls(xs, ys, mu_y_vec, sig_y_vec, fit_y);
      const bool ok_z = solve_plane_wls(xs, ys, mu_z_vec, sig_z_vec, fit_z);

      if (!ok_y || !ok_z) {
        entries[i] = {cfg.x1, cfg.x2, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, false};
      } else {
        const double c2_tot   = fit_y.chi2 + fit_z.chi2;
        const int ndof_tot_i  = fit_y.ndof + fit_z.ndof;
        const double ndof_tot = static_cast<double>(ndof_tot_i);
        const double c2_tot_red_raw =
            (ndof_tot_i > 0) ? c2_tot / ndof_tot : std::numeric_limits<double>::quiet_NaN();

        double c2_tot_red  = c2_tot_red_raw;
        double rho         = 0.0;
        double chi2_target = 1.0;
        double infl        = 1.0;

        double dx_grid = cli.grid_dx;
        double dy_grid = cli.grid_dy;
        if (dx_grid <= 0.0)
          dx_grid = estimate_min_step(xs);
        if (dy_grid <= 0.0)
          dy_grid = estimate_min_step(ys);

        bool corr_ok = false;
        if (dx_grid > 0.0 && dy_grid > 0.0 && fit_y.residuals.size() == xs.size()
            && fit_z.residuals.size() == xs.size()) {
          std::vector<double> rstd_y(xs.size()), rstd_z(xs.size());
          for (size_t j = 0; j < xs.size(); ++j) {
            const double sy = std::max(fit_y.residual_sig[j], 1e-30);
            const double sz = std::max(fit_z.residual_sig[j], 1e-30);
            rstd_y[j]       = fit_y.residuals[j] / sy;
            rstd_z[j]       = fit_z.residuals[j] / sz;
          }

          int n_pairs_y = 0;
          int n_pairs_z = 0;
          const double rho_y =
              estimate_residual_correlation(xs, ys, rstd_y, dx_grid, dy_grid, n_pairs_y);
          const double rho_z =
              estimate_residual_correlation(xs, ys, rstd_z, dx_grid, dy_grid, n_pairs_z);

          rho = 0.5 * (rho_y + rho_z);

          const double n_pairs_avg =
              0.5 * (static_cast<double>(n_pairs_y) + static_cast<double>(n_pairs_z));
          const double neighborsAvg =
              (xs.empty()) ? 0.0 : (2.0 * n_pairs_avg / static_cast<double>(xs.size()));

          infl    = compute_inflation_factor(rho, static_cast<int>(xs.size()), neighborsAvg);
          corr_ok = std::isfinite(rho) && std::isfinite(infl) && infl > 0.0
                 && (n_pairs_y + n_pairs_z) > 0;
        }

        const int N_x = count_distinct(xs);
        const int N_y = count_distinct(ys);
        chi2_target   = compute_chi2_target(rho, N_x, N_y);

        const bool apply_corr_effective = cli.apply_corr && !cli.adaptive_target;
        if (apply_corr_effective && std::isfinite(infl) && infl > 1e-6)
          c2_tot_red = c2_tot_red_raw / infl;

        double metric = 0.0;
        if (cli.corr_map) {
          metric = rho;
        } else if (cli.dist_to_n) {
          metric = std::abs(c2_tot_red - cli.dist_n);
        } else if (cli.adaptive_target) {
          metric = std::abs(c2_tot_red_raw - chi2_target);
        } else {
          metric = cli.use_reduced ? c2_tot_red : c2_tot;
        }

        const bool valid =
            cli.corr_map ? (corr_ok && std::isfinite(metric)) : std::isfinite(metric);
        if (!valid || !std::isfinite(c2_tot_red)) {
          entries[i] = {cfg.x1, cfg.x2, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, false};
        } else {
          entries[i] =
              {cfg.x1, cfg.x2, metric, c2_tot_red, c2_tot_red_raw, rho, chi2_target, infl, true};
        }
      }
    }

    int done = ++processed;
    if (done % print_step == 0 || done == total_cfgs) {
      std::string metric_label = "chi2";
      if (cli.corr_map) {
        metric_label = "rho";
      } else if (cli.dist_to_n) {
        metric_label = "|chi2/ndof-n|";
      } else if (cli.adaptive_target) {
        metric_label = "|chi2/ndof-target|";
      } else if (cli.use_reduced) {
        metric_label = "chi2/ndof";
      }
      std::lock_guard<std::mutex> lk(log_mtx);
      std::cout << "  [" << std::setw(3) << done << "/" << total_cfgs
                << "] x1=" << fmt(cfg.x1) << " " << metric_label << "="
                << (entries[i].valid ? fmt(entries[i].metric, 4) : "N/A") << "\n";
    }
  }

  if (!cli.tsv_path.empty()) {
    std::filesystem::create_directories(std::filesystem::path(cli.tsv_path).parent_path());
    std::ofstream out(cli.tsv_path);
    if (!out.is_open()) {
      std::cerr << "Errore: impossibile aprire " << cli.tsv_path << "\n";
      return 1;
    }

    out << "x1\tx2\tmetric\tchi2_red\tchi2_red_raw\trho\tchi2_target\tinfl\tvalid\n";
    for (const auto& e : entries) {
      out << std::fixed << std::setprecision(6) << e.x1 << "\t" << e.x2 << "\t" << e.metric << "\t"
          << e.chi2_red << "\t" << e.chi2_red_raw << "\t" << e.rho << "\t" << e.chi2_target << "\t"
          << e.infl << "\t" << (e.valid ? 1 : 0) << "\n";
    }
  }

  // --- CALCOLO LIMITI PERCENTILI ---
  std::vector<double> v;
  for (const auto& e : entries)
    if (e.valid)
      v.push_back(e.metric);
  std::sort(v.begin(), v.end());

  double z_min = 0.0, z_max = 1.0;
  if (!v.empty()) {
    size_t i_lo = static_cast<size_t>((cli.perc_low / 100.0) * (v.size() - 1));
    size_t i_hi = static_cast<size_t>((cli.perc_high / 100.0) * (v.size() - 1));
    z_min       = v[std::clamp(i_lo, (size_t)0, v.size() - 1)];
    z_max       = v[std::clamp(i_hi, (size_t)0, v.size() - 1)];
    if (z_min == z_max)
      z_max += 1.0;
  }
  if (cli.log_scale && z_min <= 0.0) {
    double min_pos = std::numeric_limits<double>::infinity();
    for (double vv : v) {
      if (vv > 0.0)
        min_pos = std::min(min_pos, vv);
    }
    if (std::isfinite(min_pos))
      z_min = min_pos;
    else
      z_min = 1e-6;
  }

  set_root_style();
  set_viridis_palette(!cli.corr_map);

  TH2D h_chi2("h_chi2", "", bins_x, ax_x1_lo, ax_x1_hi, bins_y, ax_x2_lo, ax_x2_hi);
  TH2D h_inv("h_inv", "", bins_x, ax_x1_lo, ax_x1_hi, bins_y, ax_x2_lo, ax_x2_hi);

  for (const auto& e : entries) {
    int bx = h_chi2.GetXaxis()->FindFixBin(e.x1);
    int by = h_chi2.GetYaxis()->FindFixBin(e.x2);
    if (e.valid)
      h_chi2.SetBinContent(bx, by, e.metric);
    else
      h_inv.SetBinContent(bx, by, 1.0);
  }

  h_chi2.SetMinimum(z_min);
  h_chi2.SetMaximum(z_max);
  h_chi2.GetXaxis()->SetTitle("x_{1} [mm]");
  h_chi2.GetYaxis()->SetTitle("x_{2} [mm]");

  std::string cb_title = "#chi^{2}";
  if (cli.corr_map) {
    cb_title = "#rho  [a.d.]";
  } else if (cli.adaptive_target) {
    cb_title = "|#chi^{2}/ndof - #chi_{target}^{2}|";
  } else if (cli.dist_to_n) {
    cb_title = "|#chi^{2}/ndof - n|";
  } else if (cli.use_reduced) {
    cb_title = "#chi^{2}/ndof";
  }
  h_chi2.GetZaxis()->SetTitle(cb_title.c_str());
  apply_zaxis_style(&h_chi2);

  const Int_t invalid_color = TColor::GetColor(80, 80, 80);
  TCanvas* c                = make_map_canvas("chi2_map", cli.log_scale);
  h_chi2.Draw("COLZ");
  h_inv.SetFillColor(invalid_color);
  h_inv.Draw("BOX same");
  c->Update();
  draw_na_legend(&h_chi2, invalid_color);

  if (cli.corr_map) {
    draw_map_title("Correlazione residui (fit piano #mu_{y,z})");
  } else if (cli.adaptive_target) {
    draw_map_title("Linearit#grave{a} della risposta  -  distanza da target");
  } else if (cli.dist_to_n) {
    draw_map_title("Distanza da #chi^{2}/ndof = " + fmt(cli.dist_n, 3)
                    + "  (fit piano #mu_{y,z})");
  } else {
    draw_map_title("Linearit#grave{a} della risposta (fit piano #mu_{y,z})");
  }

  std::filesystem::create_directories(std::filesystem::path(cli.output_path).parent_path());
  c->Print(cli.output_path.c_str());

  return 0;
}
