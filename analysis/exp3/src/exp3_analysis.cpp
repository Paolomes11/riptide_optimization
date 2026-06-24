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

#include "exp3_analysis.hpp"

// exp2 riuso diretto
#include "exp2_analysis.hpp"
#include "trace_extractor.hpp"

// exp_common
#include "fits_io.hpp"
#include "stacking.hpp"

// nlohmann/json per save/load LineCalib
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace riptide::exp3 {

// ---------------------------------------------------------------------------
// Helper interno: sottrae stacked background da stacked signal
// ---------------------------------------------------------------------------
static std::vector<double> make_diff(const stack::StackedImage& sig,
                                     const stack::StackedImage& bg) {
    std::vector<double> diff(sig.npixels());
    for (size_t i = 0; i < diff.size(); ++i)
        diff[i] = sig.mean[i] - bg.mean[i];
    return diff;
}

// ---------------------------------------------------------------------------
// Helper interno: stack segnale + background con configurazione standard
// ---------------------------------------------------------------------------
static void stack_frames(const std::vector<fits::FitsFrame>& sig_frames,
                          const std::vector<fits::FitsFrame>& bg_frames,
                          const StackingConfig& sc,
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
// Helper: centroide 1D pesato sul vettore P nel range [x0, x1]
// ---------------------------------------------------------------------------
static double centroid_1d(const std::vector<double>& P, int x0, int x1) {
    double sw = 0.0, swx = 0.0;
    x0 = std::max(0, x0);
    x1 = std::min(static_cast<int>(P.size()) - 1, x1);
    for (int x = x0; x <= x1; ++x) {
        double v = std::max(0.0, P[static_cast<size_t>(x)]);
        sw  += v;
        swx += v * static_cast<double>(x);
    }
    return (sw > 1e-12) ? swx / sw : static_cast<double>((x0 + x1) / 2);
}

// ---------------------------------------------------------------------------
// Calibrazione via linea di riferimento
// ---------------------------------------------------------------------------
LineCalib calibrate_line(const fs::path& signal_dir,
                         const fs::path& bg_dir,
                         const Exp3Config& cfg) {
    auto sig_frames = fits::read_fits_stack(signal_dir);
    if (sig_frames.empty())
        throw std::runtime_error("calibrate_line: nessun frame in " + signal_dir.string());

    std::vector<fits::FitsFrame> bg_frames;
    if (fs::exists(bg_dir))
        bg_frames = fits::read_fits_stack(bg_dir);

    stack::StackedImage sig_stacked, bg_stacked;
    stack_frames(sig_frames, bg_frames, cfg.stacking, sig_stacked, bg_stacked);

    auto diff = make_diff(sig_stacked, bg_stacked);
    int W = sig_stacked.width;
    int H = sig_stacked.height;

    // Posizione approssimativa del centro nel sensore
    int cy = static_cast<int>(cfg.optical_axis_center_px[1] *
                               static_cast<double>(H) / cfg.display.height_px);
    cy = std::max(0, std::min(H - 1, cy));

    // Profilo colonna P(x) integrato su una banda verticale attorno al centro
    int y_band = 30;
    int y_lo = std::max(0, cy - y_band);
    int y_hi = std::min(H, cy + y_band);

    std::vector<double> P(static_cast<size_t>(W), 0.0);
    for (int y = y_lo; y < y_hi; ++y)
        for (int x = 0; x < W; ++x)
            P[static_cast<size_t>(x)] += diff[static_cast<size_t>(y) * static_cast<size_t>(W) +
                                              static_cast<size_t>(x)];

    double p_peak = *std::max_element(P.begin(), P.end());
    if (p_peak < 1e-9)
        throw std::runtime_error("calibrate_line: nessun segnale rilevato in " +
                                 signal_dir.string());

    double threshold = p_peak * 0.3;

    int x_left_idx  = -1;
    int x_right_idx = -1;
    for (int x = 0; x < W; ++x)
        if (P[static_cast<size_t>(x)] > threshold) { x_left_idx = x; break; }
    for (int x = W - 1; x >= 0; --x)
        if (P[static_cast<size_t>(x)] > threshold) { x_right_idx = x; break; }

    if (x_left_idx < 0 || x_right_idx < 0 || x_right_idx <= x_left_idx)
        throw std::runtime_error("calibrate_line: impossibile rilevare gli estremi della linea");

    double x_left  = centroid_1d(P, x_left_idx  - 10, x_left_idx  + 10);
    double x_right = centroid_1d(P, x_right_idx - 10, x_right_idx + 10);

    double L_sens_px = x_right - x_left;
    if (L_sens_px < 1.0)
        throw std::runtime_error("calibrate_line: lunghezza linea sensore troppo piccola");

    double L_disp_px = static_cast<double>(cfg.calib_line_length_px);

    LineCalib calib;
    calib.L_px_display         = L_disp_px;
    calib.L_px_sensor          = L_sens_px;
    calib.scale_mm_per_sens_px = (L_disp_px * cfg.display.mm_per_px_x) / L_sens_px;
    calib.d_ax_mm              = 0.0;

    return calib;
}

// ---------------------------------------------------------------------------
// Serializzazione LineCalib
// ---------------------------------------------------------------------------
void save_line_calib(const LineCalib& calib, const fs::path& path) {
    fs::create_directories(path.parent_path());
    nlohmann::json j;
    j["scale_mm_per_sens_px"] = calib.scale_mm_per_sens_px;
    j["L_px_display"]         = calib.L_px_display;
    j["L_px_sensor"]          = calib.L_px_sensor;
    j["d_ax_mm"]              = calib.d_ax_mm;
    std::ofstream ofs(path);
    if (!ofs)
        throw std::runtime_error("save_line_calib: impossibile aprire " + path.string());
    ofs << j.dump(4) << "\n";
}

LineCalib load_line_calib(const fs::path& path) {
    std::ifstream ifs(path);
    if (!ifs)
        throw std::runtime_error("load_line_calib: file non trovato: " + path.string());
    nlohmann::json j;
    try { ifs >> j; } catch (const nlohmann::json::parse_error& e) {
        throw std::runtime_error(std::string("load_line_calib: JSON malformato: ") + e.what());
    }
    LineCalib calib;
    calib.scale_mm_per_sens_px = j.value("scale_mm_per_sens_px", 0.0);
    calib.L_px_display         = j.value("L_px_display",         0.0);
    calib.L_px_sensor          = j.value("L_px_sensor",          0.0);
    calib.d_ax_mm              = j.value("d_ax_mm",              0.0);
    return calib;
}

// ---------------------------------------------------------------------------
// Segmentazione delle tre tracce parallele
// ---------------------------------------------------------------------------
std::vector<LineROI> segment_parallel_lines(const std::vector<double>& diff,
                                             int W, int H,
                                             const Exp3Config& cfg,
                                             const LineCalib& /*calib*/) {
    std::vector<LineROI> rois;
    rois.reserve(cfg.radial_offsets_px.size());

    for (int i = 0; i < static_cast<int>(cfg.radial_offsets_px.size()); ++i) {
        int offset_disp = cfg.radial_offsets_px[i];

        // Converti display → sensore approssimativamente
        int seed_y = static_cast<int>(
            (cfg.optical_axis_center_px[1] + static_cast<double>(offset_disp)) *
            static_cast<double>(H) / cfg.display.height_px);
        seed_y = std::max(20, std::min(H - 21, seed_y));

        // Cerca la riga con massima intensità integrata in [seed_y-50, seed_y+50]
        int search_half = 50;
        int y_best  = seed_y;
        double best = -1e300;

        for (int y = std::max(0, seed_y - search_half);
                 y < std::min(H, seed_y + search_half); ++y) {
            double row_sum = 0.0;
            for (int x = 0; x < W; ++x)
                row_sum += diff[static_cast<size_t>(y) * static_cast<size_t>(W) +
                                static_cast<size_t>(x)];
            if (row_sum > best) { best = row_sum; y_best = y; }
        }

        LineROI roi;
        roi.y_center   = y_best;
        roi.half_width = 20;
        roi.line_idx   = i;
        // offset radiale in mm fisici (usando la scala display)
        roi.r_mm = static_cast<double>(offset_disp) * cfg.display.mm_per_px_y;
        rois.push_back(roi);
    }

    return rois;
}

// ---------------------------------------------------------------------------
// Helper: smoothing gaussiano 1D (kernel troncato a ±3σ, boundary REFLECT)
// Boundary REFLECT: idx < 0 → -idx; idx >= n → 2*(n-1) - idx.
// Molto meglio di CLAMP per stimare il baseline con kernel larghi (σ ~ centinaia di px).
// ---------------------------------------------------------------------------
static std::vector<double> gaussian_smooth_1d(const std::vector<double>& S, double sigma) {
    const int n    = static_cast<int>(S.size());
    const int half = std::min(static_cast<int>(3.0 * sigma + 0.5), n / 2);
    std::vector<double> kernel(static_cast<size_t>(2 * half + 1));
    double ksum = 0.0;
    for (int k = -half; k <= half; ++k) {
        kernel[static_cast<size_t>(k + half)] = std::exp(-0.5 * k * k / (sigma * sigma));
        ksum += kernel[static_cast<size_t>(k + half)];
    }
    for (auto& v : kernel) v /= ksum;

    std::vector<double> out(static_cast<size_t>(n), 0.0);
    for (int y = 0; y < n; ++y) {
        for (int k = -half; k <= half; ++k) {
            int idx = y + k;
            if (idx < 0)   idx = -idx;           // rifletti a sinistra
            if (idx >= n)  idx = 2 * (n - 1) - idx; // rifletti a destra
            idx = std::clamp(idx, 0, n - 1);     // sicurezza per doppia riflessione
            out[static_cast<size_t>(y)] += kernel[static_cast<size_t>(k + half)] * S[static_cast<size_t>(idx)];
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// Auto-rilevamento tracce tramite profilo di riga (somma firmata + multi-scala)
// ---------------------------------------------------------------------------
std::vector<LineROI> detect_lines_auto(const std::vector<double>& diff,
                                        int W, int H,
                                        const Exp3Config& cfg) {
    // 1. Somma firmata per riga — nessun clipping (né per-pixel né per-riga).
    //    Se il diff è negativo (bg ≥ signal), clippare a 0 azzererebbe tutto e il
    //    confronto con la soglia assoluta farebbe break immediato senza trovare nulla.
    //    La mediana gestisce il baseline nel peak-finding a valle.
    std::vector<double> S(static_cast<size_t>(H), 0.0);
    for (int y = 0; y < H; ++y) {
        double row_sum = 0.0;
        for (int x = 0; x < W; ++x)
            row_sum += diff[static_cast<size_t>(y) * static_cast<size_t>(W)
                           + static_cast<size_t>(x)];
        S[static_cast<size_t>(y)] = row_sum;
    }

    // 2. Detrending: rimuovi il baseline lento (gradienti flat-field, vignettatura, profilo fascio).
    //    σ_bg = H/8 ≈ 700 px cattura variazioni su scala ~ 100s di px senza toccare i picchi
    //    di linea (FWHM ≤ 200 px).  La mediana nel peak-finding a valle centra il profilo
    //    residuo indipendentemente dal segno del diff.
    const double sigma_bg = static_cast<double>(H) / 8.0;
    const auto   S_bg     = gaussian_smooth_1d(S, sigma_bg);
    std::vector<double> S_det(static_cast<size_t>(H));
    for (int y = 0; y < H; ++y)
        S_det[static_cast<size_t>(y)] = S[static_cast<size_t>(y)] - S_bg[static_cast<size_t>(y)];

    // 3. Multi-scala: σ = 5.5 px (fuoco netto), 25 px (lieve defocus), 75 px (forte defocus)
    const std::vector<double> sigmas    = {5.5, 25.0, 75.0};
    const int    sep       = cfg.trace_extraction.line_min_separation_px;
    const int    max_lines = cfg.trace_extraction.max_lines;
    const double min_snr   = cfg.trace_extraction.min_snr;

    std::vector<int> peaks;
    peaks.reserve(static_cast<size_t>(max_lines));

    for (double sigma : sigmas) {
        auto Ss = gaussian_smooth_1d(S_det, sigma);

        // Soglia MAD robusta sul profilo smoothed a questa scala.
        // min_snr è volutamente basso (≤ 1.5): il MAD risulta inflazionato di 3–10× rispetto
        // al rumore di shot teorico (RFPN + disallineamento esposizioni) perciò anche con
        // min_snr=1.4 si ottiene un vero SNR di fotoni >> 3σ per tutti i picchi reali.
        // Il filtro secondario (aspect ratio blob, n_valid_slices) protegge dai falsi positivi.
        std::vector<double> sorted_Ss = Ss;
        std::nth_element(sorted_Ss.begin(),
                         sorted_Ss.begin() + static_cast<ptrdiff_t>(sorted_Ss.size() / 2),
                         sorted_Ss.end());
        const double med = sorted_Ss[sorted_Ss.size() / 2];
        std::vector<double> abs_dev(static_cast<size_t>(H));
        for (int i = 0; i < H; ++i)
            abs_dev[static_cast<size_t>(i)] = std::abs(Ss[static_cast<size_t>(i)] - med);
        std::nth_element(abs_dev.begin(),
                         abs_dev.begin() + static_cast<ptrdiff_t>(abs_dev.size() / 2),
                         abs_dev.end());
        const double threshold = min_snr * abs_dev[abs_dev.size() / 2] * 1.4826;

        std::vector<double> work = Ss;
        for (auto& v : work) v -= med;
        for (int n = 0; n < max_lines; ++n) {
            auto it = std::max_element(work.begin(), work.end());
            if (*it < threshold) break;
            const int y_peak = static_cast<int>(std::distance(work.begin(), it));

            // Scarta se già trovato (a questa scala o da scala precedente)
            bool dup = false;
            for (int p : peaks)
                if (std::abs(p - y_peak) < sep) { dup = true; break; }

            if (!dup)
                peaks.push_back(y_peak);

            const int y0 = std::max(0, y_peak - sep);
            const int y1 = std::min(H - 1, y_peak + sep);
            for (int k = y0; k <= y1; ++k)
                work[static_cast<size_t>(k)] = 0.0;
        }
    }

    // 3. Ordina per Y crescente; tronca a max_lines
    std::sort(peaks.begin(), peaks.end());
    if (static_cast<int>(peaks.size()) > max_lines)
        peaks.resize(static_cast<size_t>(max_lines));

    // 4. Costruisce LineROI
    std::vector<LineROI> rois;
    rois.reserve(peaks.size());

    const double H_disp = static_cast<double>(cfg.display.height_px);
    const double axis_y = cfg.optical_axis_center_px[1];

    for (int i = 0; i < static_cast<int>(peaks.size()); ++i) {
        const int    y_best = peaks[static_cast<size_t>(i)];
        const double y_disp = static_cast<double>(y_best) / static_cast<double>(H) * H_disp;
        const double r_disp = y_disp - axis_y;

        LineROI roi;
        roi.y_center   = y_best;
        roi.half_width = cfg.trace_extraction.half_width_px;
        roi.line_idx   = i;
        roi.r_mm       = r_disp * cfg.display.mm_per_px_y;
        rois.push_back(roi);
    }

    return rois;
}

// ---------------------------------------------------------------------------
// Pipeline di misura per le tre tracce parallele
// ---------------------------------------------------------------------------
std::vector<QResult> run_measurement_parallel_lines(
    const fs::path& signal_dir,
    const fs::path& bg_dir,
    const Exp3Config& cfg,
    const LineCalib& calib,
    double d_ax_mm)
{
    auto sig_frames = fits::read_fits_stack(signal_dir);
    if (sig_frames.empty())
        throw std::runtime_error("run_measurement_parallel_lines: nessun frame in " +
                                 signal_dir.string());

    std::vector<fits::FitsFrame> bg_frames;
    if (fs::exists(bg_dir))
        bg_frames = fits::read_fits_stack(bg_dir);

    stack::StackedImage sig_stacked, bg_stacked;
    stack_frames(sig_frames, bg_frames, cfg.stacking, sig_stacked, bg_stacked);

    auto diff = make_diff(sig_stacked, bg_stacked);
    int W = sig_stacked.width;
    int H = sig_stacked.height;

    auto rois = cfg.radial_offsets_px.empty()
                ? detect_lines_auto(diff, W, H, cfg)
                : segment_parallel_lines(diff, W, H, cfg, calib);

    std::vector<QResult> results;
    results.reserve(rois.size());

    exp2::TraceConfig tc;
    tc.min_snr = cfg.trace_extraction.min_snr;

    for (const auto& roi : rois) {
        QResult qr;
        qr.d_ax_mm = d_ax_mm;
        qr.r_mm    = roi.r_mm;
        qr.r_idx   = roi.line_idx;
        qr.warning = false;

        try {
            // Ritaglia il diff alla ROI
            int y_lo   = std::max(0, roi.y_center - roi.half_width);
            int y_hi   = std::min(H, roi.y_center + roi.half_width + 1);
            int crop_H = y_hi - y_lo;

            std::vector<double> crop(static_cast<size_t>(W * crop_H));
            for (int y = y_lo; y < y_hi; ++y)
                for (int x = 0; x < W; ++x)
                    crop[static_cast<size_t>((y - y_lo) * W + x)] =
                        diff[static_cast<size_t>(y * W + x)];

            // Stima angolo + estrazione profilo con raffinamento iterativo
            double angle = exp2::estimate_trace_angle(crop, W, crop_H, tc, nullptr);
            auto trace   = exp2::extract_trace_profile(crop, W, crop_H, angle, tc);

            for (int iter = 0; iter < 2 && trace.n_valid_slices >= 5; ++iter) {
                auto fit_tmp = exp2::fit_centroid_line(trace);
                if (!fit_tmp.converged) break;
                double delta = std::atan(fit_tmp.a) * 180.0 / M_PI;
                if (std::abs(delta) < 0.05 || std::abs(delta) > 10.0) break;
                angle += delta;
                trace = exp2::extract_trace_profile(crop, W, crop_H, angle, tc);
            }

            auto fit = exp2::fit_centroid_line(trace);
            qr.chi2_ndof      = fit.chi2_ndof;
            qr.n_valid_slices = trace.n_valid_slices;

            if (!trace.trace_detected ||
                trace.n_valid_slices < cfg.trace_extraction.min_valid_slices)
                qr.warning = true;
        } catch (const std::exception& e) {
            // Linea non caratterizzabile (troppo defocused, blob circolare, ecc.)
            qr.chi2_ndof      = 0.0;
            qr.n_valid_slices = 0;
            qr.warning        = true;
            std::cerr << "[exp3 measure] WARNING r_idx=" << roi.line_idx
                      << " y=" << roi.y_center << ": " << e.what() << "\n";
        }

        results.push_back(qr);
    }

    return results;
}

// ---------------------------------------------------------------------------
// Carica Q_sim da TSV — nearest-neighbor entro 10 mm, config_valid=1
// Colonne TSV: x1 x2 metric Q_raw Q_target rho_hat n_traces n_failed n_invalid config_valid
// Restituisce la colonna "metric" (chi2/ndof della sim) del punto più vicino.
// ---------------------------------------------------------------------------
double load_Q_sim(const fs::path& tsv_path, double x1_mm, double x2_mm) {
    std::ifstream ifs(tsv_path);
    if (!ifs)
        return std::numeric_limits<double>::quiet_NaN();

    std::string line;
    std::getline(ifs, line);  // header

    const double max_dist_sq = 10.0 * 10.0;
    double best_dist_sq = max_dist_sq + 1.0;
    double best_q       = std::numeric_limits<double>::quiet_NaN();

    while (std::getline(ifs, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);

        double x1_f, x2_f, metric, q_raw, q_target, rho_hat;
        int    n_traces, n_failed, n_invalid, config_valid;

        if (!(ss >> x1_f >> x2_f >> metric >> q_raw >> q_target >> rho_hat
                 >> n_traces >> n_failed >> n_invalid >> config_valid))
            continue;
        if (config_valid != 1 || !std::isfinite(metric))
            continue;

        double dx = x1_f - x1_mm, dy = x2_f - x2_mm;
        double dist_sq = dx * dx + dy * dy;
        if (dist_sq < best_dist_sq) {
            best_dist_sq = dist_sq;
            best_q       = metric;
        }
    }
    return best_q;
}

// ---------------------------------------------------------------------------
// Scrittura TSV
// ---------------------------------------------------------------------------
void write_q_results_tsv(const std::string& config,
                          const std::vector<QResult>& results,
                          const fs::path& path) {
    fs::create_directories(path.parent_path());
    std::ofstream ofs(path);
    if (!ofs)
        throw std::runtime_error("write_q_results_tsv: impossibile aprire " + path.string());

    ofs << std::fixed << std::setprecision(6);
    ofs << "config\td_ax_mm\tr_mm\tr_idx\tchi2_ndof\tn_valid_slices\twarning\n";

    for (const auto& r : results) {
        ofs << config        << '\t'
            << r.d_ax_mm     << '\t'
            << r.r_mm        << '\t'
            << r.r_idx       << '\t'
            << r.chi2_ndof   << '\t'
            << r.n_valid_slices << '\t'
            << (r.warning ? 1 : 0) << '\n';
    }
}

} // namespace riptide::exp3
