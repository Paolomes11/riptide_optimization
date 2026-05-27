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

#include "exp3b_analysis.hpp"

#include "fits_io.hpp"
#include "stacking.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace riptide::exp3b {

// ---------------------------------------------------------------------------
// Helper: make_diff e stack_frames (ridondanti rispetto ad exp3, ma isolati)
// ---------------------------------------------------------------------------
static std::vector<double> make_diff(const stack::StackedImage& sig,
                                     const stack::StackedImage& bg) {
    std::vector<double> d(sig.npixels());
    for (size_t i = 0; i < d.size(); ++i)
        d[i] = sig.mean[i] - bg.mean[i];
    return d;
}

static void stack_frames(const std::vector<fits::FitsFrame>& sig_frames,
                          const std::vector<fits::FitsFrame>& bg_frames,
                          const exp3::StackingConfig& sc,
                          stack::StackedImage& sig_out,
                          stack::StackedImage& bg_out) {
    stack::StackConfig cfg;
    cfg.n_sigma    = sc.n_sigma;
    cfg.n_iter     = sc.n_iter;
    cfg.min_frames = sc.min_frames;
    cfg.method = (static_cast<int>(sig_frames.size()) >= sc.min_frames)
                 ? stack::StackMethod::SigmaClip
                 : stack::StackMethod::Mean;

    sig_out = stack::sigma_clip_stack(sig_frames, cfg);
    if (!bg_frames.empty()) {
        bg_out = stack::mean_stack(bg_frames);
    } else {
        bg_out.width  = sig_out.width;
        bg_out.height = sig_out.height;
        bg_out.mean.assign(static_cast<size_t>(sig_out.width * sig_out.height), 0.0);
        bg_out.sigma.assign(static_cast<size_t>(sig_out.width * sig_out.height), 0.0);
        bg_out.count.assign(static_cast<size_t>(sig_out.width * sig_out.height), 1);
    }
}

// ---------------------------------------------------------------------------
// Fit gaussiano 2D (centroide pesato) su una finestra dell'immagine.
// Restituisce true se il SNR è sufficiente, e scrive (cx, cy, sigma_px).
// ---------------------------------------------------------------------------
static bool fit_dot_centroid_2d(const std::vector<double>& img, int W, int H,
                                 int cx0, int cy0, int half,
                                 double& cx_out, double& cy_out,
                                 double& sigma_out) {
    double sum_w = 0.0, sum_wx = 0.0, sum_wy = 0.0;
    double sum_w2 = 0.0;
    double peak = 0.0;

    for (int dy = -half; dy <= half; ++dy) {
        for (int dx = -half; dx <= half; ++dx) {
            int x = cx0 + dx, y = cy0 + dy;
            if (x < 0 || x >= W || y < 0 || y >= H) continue;
            double v = img[static_cast<size_t>(y * W + x)];
            if (v < 0.0) continue;
            peak = std::max(peak, v);
            sum_w  += v;
            sum_wx += v * static_cast<double>(x);
            sum_wy += v * static_cast<double>(y);
            sum_w2 += v * v;
        }
    }

    if (peak < 1e-9 || sum_w < 1e-12) return false;

    // SNR: confronta picco con std dei pixel di bordo
    double bg_est = 0.0;
    double bg2    = 0.0;
    int n_border  = 0;
    for (int dy = -half; dy <= half; ++dy) {
        for (int dx = -half; dx <= half; ++dx) {
            if (std::abs(dx) != half && std::abs(dy) != half) continue;
            int x = cx0 + dx, y = cy0 + dy;
            if (x < 0 || x >= W || y < 0 || y >= H) continue;
            double v = img[static_cast<size_t>(y * W + x)];
            bg_est += v; bg2 += v * v; ++n_border;
        }
    }
    if (n_border > 0) {
        bg_est /= static_cast<double>(n_border);
        double var   = bg2 / static_cast<double>(n_border) - bg_est * bg_est;
        double noise = (var > 0.0) ? std::sqrt(var) : 1.0;
        if ((peak - bg_est) < 3.0 * noise) return false;
    }

    cx_out = sum_wx / sum_w;
    cy_out = sum_wy / sum_w;

    // Approssima sigma come RMS della distribuzione dei pesi
    double var_x = 0.0;
    for (int dy = -half; dy <= half; ++dy) {
        for (int dx = -half; dx <= half; ++dx) {
            int x = cx0 + dx, y = cy0 + dy;
            if (x < 0 || x >= W || y < 0 || y >= H) continue;
            double v = img[static_cast<size_t>(y * W + x)];
            if (v < 0.0) continue;
            double dd = static_cast<double>(x) - cx_out;
            var_x += v * dd * dd;
        }
    }
    sigma_out = (sum_w > 1e-12) ? std::sqrt(var_x / sum_w) : 1.0;
    if (sigma_out < 0.1) sigma_out = 0.5;

    return true;
}

// ---------------------------------------------------------------------------
// run_exp3b
// ---------------------------------------------------------------------------
MResult run_exp3b(const fs::path& signal_dir,
                  const fs::path& bg_dir,
                  const Exp3Config& cfg,
                  const LineCalib& calib,
                  double d_ax_mm)
{
    auto sig_frames = fits::read_fits_stack(signal_dir);
    if (sig_frames.empty())
        throw std::runtime_error("run_exp3b: nessun frame in " + signal_dir.string());

    std::vector<fits::FitsFrame> bg_frames;
    if (fs::exists(bg_dir))
        bg_frames = fits::read_fits_stack(bg_dir);

    stack::StackedImage sig_stacked, bg_stacked;
    stack_frames(sig_frames, bg_frames, cfg.stacking, sig_stacked, bg_stacked);

    auto diff = make_diff(sig_stacked, bg_stacked);
    int W = sig_stacked.width;
    int H = sig_stacked.height;

    // Posizione X del centro asse ottico sul sensore
    int cx_sens = static_cast<int>(cfg.optical_axis_center_px[0] *
                                    static_cast<double>(W) / cfg.display.width_px);
    cx_sens = std::max(10, std::min(W - 11, cx_sens));

    // Numero di dot attesi
    int n_dots_expected = cfg.display.height_px / cfg.dot_column_step_px;
    int dot_half        = 10; // finestra di ricerca ±10 px

    std::vector<DotCentroid> dots;
    dots.reserve(static_cast<size_t>(n_dots_expected));

    for (int k = 0; k < n_dots_expected; ++k) {
        // Posizione teorica nel display
        double y_disp_px = static_cast<double>(k) * static_cast<double>(cfg.dot_column_step_px);
        double y_disp_mm = y_disp_px * cfg.display.mm_per_px_y;

        // Converti in coordinate sensore approssimative
        int cy_seed = static_cast<int>(y_disp_px * static_cast<double>(H) /
                                        cfg.display.height_px);
        cy_seed = std::max(dot_half, std::min(H - dot_half - 1, cy_seed));

        double cx_found, cy_found, sigma_found;
        if (!fit_dot_centroid_2d(diff, W, H, cx_sens, cy_seed, dot_half,
                                  cx_found, cy_found, sigma_found))
            continue;

        DotCentroid dc;
        dc.y_disp_mm = y_disp_mm;
        dc.y_sens_px = cy_found;
        dc.y_sens_mm = cy_found * calib.scale_mm_per_sens_px;
        dc.sigma_px  = sigma_found;
        dots.push_back(dc);
    }

    MResult result;
    result.d_ax_mm = d_ax_mm;
    result.n_dots  = static_cast<int>(dots.size());

    if (result.n_dots < 3) {
        result.M_global = std::numeric_limits<double>::quiet_NaN();
        return result;
    }

    // Fit lineare: y_sens_mm = M * y_disp_mm + q
    double sum_x = 0.0, sum_y = 0.0, sum_xx = 0.0, sum_xy = 0.0;
    int n = result.n_dots;
    for (const auto& dc : dots) {
        sum_x  += dc.y_disp_mm;
        sum_y  += dc.y_sens_mm;
        sum_xx += dc.y_disp_mm * dc.y_disp_mm;
        sum_xy += dc.y_disp_mm * dc.y_sens_mm;
    }
    double denom = static_cast<double>(n) * sum_xx - sum_x * sum_x;
    if (std::abs(denom) < 1e-12) {
        result.M_global = std::numeric_limits<double>::quiet_NaN();
        return result;
    }
    double M = (static_cast<double>(n) * sum_xy - sum_x * sum_y) / denom;
    double q = (sum_y - M * sum_x) / static_cast<double>(n);

    result.M_global = M;
    result.q_offset = q;

    // Residui e chi2
    double chi2 = 0.0, rms2 = 0.0;
    for (const auto& dc : dots) {
        double res    = dc.y_sens_mm - M * dc.y_disp_mm - q;
        double sigma  = dc.sigma_px * calib.scale_mm_per_sens_px;
        if (sigma < 1e-9) sigma = 0.5 * calib.scale_mm_per_sens_px;
        chi2 += (res * res) / (sigma * sigma);
        rms2 += res * res;
    }
    result.chi2_ndof       = (n > 2) ? chi2 / static_cast<double>(n - 2) : 0.0;
    result.M_rms_residual  = std::sqrt(rms2 / static_cast<double>(n));

    // M_local e r_mm
    result.r_mm.reserve(static_cast<size_t>(n));
    result.M_local.reserve(static_cast<size_t>(n - 1));
    for (const auto& dc : dots)
        result.r_mm.push_back(dc.y_sens_mm);
    for (int i = 0; i < n - 1; ++i) {
        double dy_sens = dots[static_cast<size_t>(i + 1)].y_sens_mm -
                         dots[static_cast<size_t>(i)].y_sens_mm;
        double dy_disp = dots[static_cast<size_t>(i + 1)].y_disp_mm -
                         dots[static_cast<size_t>(i)].y_disp_mm;
        result.M_local.push_back(std::abs(dy_disp) > 1e-9 ? dy_sens / dy_disp : 0.0);
    }

    return result;
}

// ---------------------------------------------------------------------------
// write_M_summary_tsv
// ---------------------------------------------------------------------------
void write_M_summary_tsv(const std::string& config,
                          const std::vector<MResult>& results,
                          const fs::path& path) {
    fs::create_directories(path.parent_path());
    std::ofstream ofs(path);
    if (!ofs)
        throw std::runtime_error("write_M_summary_tsv: impossibile aprire " + path.string());

    ofs << std::fixed << std::setprecision(6);
    ofs << "config\td_ax_mm\tM_global\tM_rms_residual\tchi2_ndof\tn_dots\n";
    for (const auto& r : results) {
        ofs << config          << '\t'
            << r.d_ax_mm       << '\t'
            << r.M_global      << '\t'
            << r.M_rms_residual << '\t'
            << r.chi2_ndof     << '\t'
            << r.n_dots        << '\n';
    }
}

} // namespace riptide::exp3b
