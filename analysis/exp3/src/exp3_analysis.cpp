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
#include "homography.hpp"

// exp2 riuso diretto
#include "exp2_analysis.hpp"
#include "trace_extractor.hpp"

// exp_common
#include "fits_io.hpp"
#include "stacking.hpp"

// ROOT
#include <TAxis.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace riptide::exp3 {

// ---------------------------------------------------------------------------
// Stile ROOT condiviso (mirror di exp2)
// ---------------------------------------------------------------------------
static void apply_riptide_style() {
    gStyle->Reset();
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetTitleFont(42, "");
    gStyle->SetStatFont(42);
    gStyle->SetTextSize(0.040f);
    gStyle->SetLabelSize(0.038f, "XYZ");
    gStyle->SetTitleSize(0.044f, "XYZ");
    gStyle->SetTitleSize(0.046f, "");
    gStyle->SetTitleOffset(1.55f, "Y");
    gStyle->SetTitleOffset(1.20f, "X");
    gStyle->SetTickLength(0.018f, "X");
    gStyle->SetTickLength(0.018f, "Y");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetPadColor(0);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetNdivisions(505, "X");
    gStyle->SetNdivisions(505, "Y");
    gStyle->SetNumberContours(255);
    gStyle->SetGridColor(kGray + 1);
    gStyle->SetGridStyle(3);
    gStyle->SetGridWidth(1);
    TGaxis::SetMaxDigits(4);
}

// ---------------------------------------------------------------------------
// Utilità: percentile su vettore
// ---------------------------------------------------------------------------
static double percentile(std::vector<double> v, double p) {
    if (v.empty())
        return 0.0;
    std::sort(v.begin(), v.end());
    double idx = p * static_cast<double>(v.size() - 1);
    size_t lo  = static_cast<size_t>(idx);
    size_t hi  = std::min(lo + 1, v.size() - 1);
    double frac = idx - static_cast<double>(lo);
    return v[lo] * (1.0 - frac) + v[hi] * frac;
}

// ---------------------------------------------------------------------------
// Costruisce diff double da due StackedImage
// ---------------------------------------------------------------------------
static std::vector<double> make_diff(const stack::StackedImage& sig,
                                     const stack::StackedImage& bg) {
    std::vector<double> diff(sig.npixels());
    for (size_t i = 0; i < diff.size(); ++i)
        diff[i] = sig.mean[i] - bg.mean[i];
    return diff;
}

// ---------------------------------------------------------------------------
// Struttura interna con tutti i dati di una singola analisi traccia
// ---------------------------------------------------------------------------
struct TraceAnalysisData {
    SingleTraceResult result;
    std::vector<double> diff;
    int width  = 0;
    int height = 0;
    exp2::TraceResult  trace;
    exp2::CentroidFitResult fit;
};

// ---------------------------------------------------------------------------
// Funzione interna che esegue l'analisi completa e restituisce anche i dati grezzi
// ---------------------------------------------------------------------------
static TraceAnalysisData analyze_single_trace_full(const fs::path& signal_dir,
                                                    const fs::path& bg_dir,
                                                    double theta_nominal_deg,
                                                    double d_ax_nominal_mm,
                                                    double d_ax_measured_mm,
                                                    const Homography& /*H*/,
                                                    const Exp3Config& cfg) {
    TraceAnalysisData out;
    out.result.d_ax_nominal_mm  = d_ax_nominal_mm;
    out.result.d_ax_measured_mm = d_ax_measured_mm;
    out.result.theta_nominal_deg = theta_nominal_deg;
    out.result.valid = false;

    // Carica frame
    std::vector<fits::FitsFrame> sig_frames, bg_frames;
    try {
        sig_frames = fits::read_fits_stack(signal_dir);
    } catch (const std::exception& e) {
        out.result.warning = std::string("Signal dir error: ") + e.what();
        return out;
    }
    try {
        bg_frames = fits::read_fits_stack(bg_dir);
    } catch (const std::exception& e) {
        out.result.warning = std::string("Background dir error: ") + e.what();
        return out;
    }

    if (sig_frames.empty()) {
        out.result.warning = "Nessun frame segnale trovato";
        return out;
    }
    if (bg_frames.empty()) {
        out.result.warning = "Nessun frame background trovato";
        return out;
    }

    out.width  = sig_frames[0].width();
    out.height = sig_frames[0].height();

    // Stacking
    stack::StackConfig sc;
    sc.n_sigma    = cfg.stack_n_sigma;
    sc.n_iter     = cfg.stack_n_iter;
    sc.min_frames = cfg.stack_min_frames;
    bool use_clip = (static_cast<int>(sig_frames.size()) >= cfg.stack_min_frames);
    sc.method = use_clip ? stack::StackMethod::SigmaClip : stack::StackMethod::Mean;
    if (!use_clip)
        out.result.warning = "Meno di " + std::to_string(cfg.stack_min_frames) +
                             " frame segnale: uso media semplice";

    stack::StackedImage sig_stacked = stack::sigma_clip_stack(sig_frames, sc);
    stack::StackedImage bg_stacked  = stack::mean_stack(bg_frames);

    out.diff = make_diff(sig_stacked, bg_stacked);

    // Configurazione estrattore
    exp2::TraceConfig tc;
    tc.min_snr = cfg.min_snr;

    // Stima angolo
    double angle_err = 0.0;
    double angle = exp2::estimate_trace_angle(out.diff, out.width, out.height, tc, &angle_err);

    // Estrazione profilo con raffinamento iterativo (come in exp2)
    out.trace = exp2::extract_trace_profile(out.diff, out.width, out.height, angle, tc);
    for (int iter = 0; iter < 2 && out.trace.n_valid_slices >= tc.trace_trim_min_slices; ++iter) {
        auto fit_tmp = exp2::fit_centroid_line(out.trace);
        if (!fit_tmp.converged) break;
        double delta = std::atan(fit_tmp.a) * 180.0 / M_PI;
        if (std::abs(delta) < 0.05 || std::abs(delta) > 10.0) break;
        angle += delta;
        out.trace = exp2::extract_trace_profile(out.diff, out.width, out.height, angle, tc);
    }

    out.fit = exp2::fit_centroid_line(out.trace);
    out.result.theta_measured_deg = angle;
    out.result.chi2_ndof          = out.fit.chi2_ndof;
    out.result.n_valid_slices     = out.trace.n_valid_slices;

    // Validazione
    double angle_diff = std::abs(angle - theta_nominal_deg);
    // Normalizza la differenza angolare in [0, 90] (simmetria assiale)
    while (angle_diff > 90.0) angle_diff -= 90.0;
    angle_diff = std::min(angle_diff, 90.0 - angle_diff);

    bool enough_slices = (out.trace.n_valid_slices >= cfg.min_valid_slices);
    bool angle_ok      = (angle_diff < cfg.angle_tolerance_deg);
    out.result.valid   = enough_slices && angle_ok && out.trace.trace_detected;

    if (!out.trace.trace_detected)
        out.result.warning = "Blob troppo circolare, traccia non rilevata";
    else if (!enough_slices)
        out.result.warning = "Slice valide insufficienti (" +
                             std::to_string(out.trace.n_valid_slices) + ")";
    else if (!angle_ok)
        out.result.warning = "Angolo misurato fuori tolleranza (" +
                             std::to_string(angle_diff) + " deg)";

    return out;
}

// ---------------------------------------------------------------------------
// API pubblica
// ---------------------------------------------------------------------------
SingleTraceResult analyze_single_trace(const fs::path& signal_dir,
                                       const fs::path& bg_dir,
                                       double theta_nominal_deg,
                                       double d_ax_nominal_mm,
                                       double d_ax_measured_mm,
                                       const Homography& H,
                                       const Exp3Config& cfg) {
    return analyze_single_trace_full(signal_dir, bg_dir, theta_nominal_deg,
                                     d_ax_nominal_mm, d_ax_measured_mm, H, cfg).result;
}

// ---------------------------------------------------------------------------
// Fit gaussiano 2D per il rilevamento dei dot di calibrazione.
// Restituisce (cx, cy) del centroide e false se il fit non converge.
// ---------------------------------------------------------------------------
static bool fit_dot_centroid(const std::vector<double>& img, int W, int H,
                              int cx0, int cy0, int half,
                              double& cx_out, double& cy_out) {
    // Raccoglie pixel della finestra
    double sum_w = 0.0, sum_wx = 0.0, sum_wy = 0.0;
    double peak = 0.0;
    for (int dy = -half; dy <= half; ++dy) {
        for (int dx = -half; dx <= half; ++dx) {
            int x = cx0 + dx, y = cy0 + dy;
            if (x < 0 || x >= W || y < 0 || y >= H) continue;
            double v = img[static_cast<size_t>(y) * static_cast<size_t>(W) + static_cast<size_t>(x)];
            if (v < 0.0) continue;
            peak = std::max(peak, v);
            sum_w  += v;
            sum_wx += v * static_cast<double>(x);
            sum_wy += v * static_cast<double>(y);
        }
    }

    double bg_est = 0.0;
    if (peak < 1e-6 || sum_w < 1e-12) return false;

    // SNR minimo: il picco deve essere almeno 3× la std locale di background
    // (stima grezza: std dei pixel ai bordi della finestra)
    double sum2 = 0.0;
    int n_border = 0;
    for (int dy = -half; dy <= half; ++dy) {
        for (int dx = -half; dx <= half; ++dx) {
            if (std::abs(dx) != half && std::abs(dy) != half) continue;
            int x = cx0 + dx, y = cy0 + dy;
            if (x < 0 || x >= W || y < 0 || y >= H) continue;
            double v = img[static_cast<size_t>(y) * static_cast<size_t>(W) + static_cast<size_t>(x)];
            bg_est += v;
            sum2   += v * v;
            ++n_border;
        }
    }
    if (n_border > 0) {
        bg_est /= static_cast<double>(n_border);
        double var = sum2 / static_cast<double>(n_border) - bg_est * bg_est;
        double noise = (var > 0.0) ? std::sqrt(var) : 1.0;
        if ((peak - bg_est) < 3.0 * noise) return false;
    }

    cx_out = sum_wx / sum_w;
    cy_out = sum_wy / sum_w;
    return true;
}

// ---------------------------------------------------------------------------
std::vector<CalibPoint> detect_calibration_dots(const std::vector<double>& diff_image,
                                                int width, int height,
                                                int grid_step_px,
                                                int display_width_px,
                                                int display_height_px) {
    std::vector<CalibPoint> pts;
    int half = grid_step_px / 2;

    // Itera sui nodi teorici della griglia display
    for (int gy = grid_step_px; gy < display_height_px - grid_step_px / 2; gy += grid_step_px) {
        for (int gx = grid_step_px; gx < display_width_px - grid_step_px / 2; gx += grid_step_px) {
            // Posizione teorica nel sensore: assume allineamento approssimativo
            // (la prima stima è la stessa coordinata scalata al sensore)
            double scale_x = static_cast<double>(width)  / static_cast<double>(display_width_px);
            double scale_y = static_cast<double>(height) / static_cast<double>(display_height_px);
            int cx0 = static_cast<int>(static_cast<double>(gx) * scale_x);
            int cy0 = static_cast<int>(static_cast<double>(gy) * scale_y);

            double cx_found, cy_found;
            if (fit_dot_centroid(diff_image, width, height, cx0, cy0, half,
                                  cx_found, cy_found)) {
                CalibPoint p;
                p.disp_x = static_cast<double>(gx);
                p.disp_y = static_cast<double>(gy);
                p.sens_x = cx_found;
                p.sens_y = cy_found;
                pts.push_back(p);
            }
        }
    }

    return pts;
}

// ---------------------------------------------------------------------------
Homography calibrate_distance(const fs::path& signal_dir,
                              const fs::path& bg_dir,
                              const Exp3Config& cfg) {
    auto sig_frames = fits::read_fits_stack(signal_dir);
    auto bg_frames  = fits::read_fits_stack(bg_dir);

    if (sig_frames.empty())
        throw std::runtime_error("calibrate_distance: nessun frame in " + signal_dir.string());

    int W = sig_frames[0].width();
    int H = sig_frames[0].height();

    stack::StackConfig sc;
    sc.n_sigma    = cfg.stack_n_sigma;
    sc.n_iter     = cfg.stack_n_iter;
    sc.min_frames = cfg.stack_min_frames;
    sc.method = bg_frames.empty() ? stack::StackMethod::Mean : stack::StackMethod::SigmaClip;

    auto sig_stacked = stack::sigma_clip_stack(sig_frames, sc);
    stack::StackedImage bg_stacked;
    if (!bg_frames.empty())
        bg_stacked = stack::mean_stack(bg_frames);
    else {
        bg_stacked.width  = W;
        bg_stacked.height = H;
        bg_stacked.mean.assign(static_cast<size_t>(W * H), 0.0);
        bg_stacked.sigma.assign(static_cast<size_t>(W * H), 0.0);
        bg_stacked.count.assign(static_cast<size_t>(W * H), 1);
    }

    auto diff = make_diff(sig_stacked, bg_stacked);

    auto pts = detect_calibration_dots(diff, W, H,
                                       cfg.calibration_grid_step_px,
                                       cfg.display.width_px,
                                       cfg.display.height_px);

    if (static_cast<int>(pts.size()) < 4)
        throw std::runtime_error("calibrate_distance: troppi pochi dot trovati (" +
                                 std::to_string(pts.size()) + ")");

    return compute_homography(pts);
}

// ---------------------------------------------------------------------------
double load_Q_sim(const fs::path& tsv_path, double x1_mm, double x2_mm) {
    std::ifstream ifs(tsv_path);
    if (!ifs)
        return std::numeric_limits<double>::quiet_NaN();

    std::string line;
    std::getline(ifs, line);  // header

    while (std::getline(ifs, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        double x1_f, x2_f, q;
        if (!(ss >> x1_f >> x2_f >> q)) continue;
        if (std::abs(x1_f - x1_mm) < 0.5 && std::abs(x2_f - x2_mm) < 0.5)
            return q;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

// ---------------------------------------------------------------------------
// Output ROOT: singola traccia
// ---------------------------------------------------------------------------
static void produce_trace_output(const TraceAnalysisData& data,
                                  const fs::path& output_dir,
                                  double d_ax_mm, double theta_deg,
                                  const OutputConfig& out_cfg) {
    if (!out_cfg.save_png && !out_cfg.save_root) return;
    if (data.diff.empty()) return;

    apply_riptide_style();

    int W = data.width, H = data.height;
    std::vector<double> flat = data.diff;

    double z_lo = percentile(flat, out_cfg.z_min_percentile);
    double z_hi = percentile(flat, out_cfg.z_max_percentile);
    if (z_hi <= z_lo) z_hi = z_lo + 1.0;

    // Costruisce titolo
    std::ostringstream title_oss;
    title_oss << "d_{ax}=" << std::fixed << std::setprecision(0) << d_ax_mm
              << " mm  #theta=" << std::fixed << std::setprecision(0) << theta_deg << "#circ";

    auto* c = new TCanvas("trace_canvas", title_oss.str().c_str(), 1600, 600);
    c->Divide(2, 1);

    // Pad 1: immagine diff con traccia
    c->cd(1);
    gPad->SetLeftMargin(0.16f);
    gPad->SetBottomMargin(0.14f);

    auto* h2 = new TH2D("hdiff", title_oss.str().c_str(),
                         W, 0, W, H, 0, H);
    for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
            h2->SetBinContent(x + 1, y + 1,
                flat[static_cast<size_t>(y) * static_cast<size_t>(W) + static_cast<size_t>(x)]);

    gStyle->SetPalette(kViridis);
    h2->GetZaxis()->SetRangeUser(z_lo, z_hi);
    h2->GetXaxis()->SetTitle("x [pixel]");
    h2->GetYaxis()->SetTitle("y [pixel]");
    h2->Draw("COLZ");

    // Pad 2: centroidi + fit ODR
    c->cd(2);
    gPad->SetLeftMargin(0.16f);
    gPad->SetBottomMargin(0.14f);
    gPad->SetGrid();

    // Raccoglie centroidi validi
    std::vector<double> t_vals, c_vals, c_errs;
    for (const auto& sl : data.trace.profile) {
        if (!sl.valid || !sl.in_trace) continue;
        t_vals.push_back(sl.t);
        c_vals.push_back(sl.center);
        c_errs.push_back(sl.center_err > 0.0 ? sl.center_err : 0.5);
    }

    if (!t_vals.empty()) {
        auto* gr = new TGraphErrors(static_cast<int>(t_vals.size()),
                                    t_vals.data(), c_vals.data(),
                                    nullptr, c_errs.data());
        gr->SetTitle((title_oss.str() + ";t [pixel];centroide [pixel]").c_str());
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(0.5f);
        gr->SetMarkerColor(kBlue + 1);
        gr->SetLineColor(kBlue + 1);
        gr->Draw("AP");

        // Linea ODR
        if (data.fit.converged && !t_vals.empty()) {
            double t_min = *std::min_element(t_vals.begin(), t_vals.end());
            double t_max = *std::max_element(t_vals.begin(), t_vals.end());
            auto* line = new TLine(t_min, data.fit.a * t_min + data.fit.b,
                                    t_max, data.fit.a * t_max + data.fit.b);
            line->SetLineColor(kRed);
            line->SetLineWidth(2);
            line->Draw();
        }
    }

    // Pannello info
    auto* info = new TPaveText(0.62, 0.70, 0.98, 0.98, "NDC");
    info->SetFillColor(0);
    info->SetBorderSize(1);
    info->SetTextFont(42);
    info->SetTextSize(0.035f);
    info->AddText(("#chi^{2}/ndof = " + [&](){
        std::ostringstream s;
        s << std::fixed << std::setprecision(2) << data.result.chi2_ndof;
        return s.str();
    }()).c_str());
    info->AddText(("n_{slices} = " + std::to_string(data.result.n_valid_slices)).c_str());
    {
        std::ostringstream s;
        s << std::fixed << std::setprecision(1) << data.result.theta_measured_deg;
        info->AddText(("#theta_{meas} = " + s.str() + "#circ").c_str());
    }
    info->Draw();

    // Salva
    std::ostringstream fname;
    fname << "trace_d" << std::fixed << std::setprecision(0) << d_ax_mm
          << "_theta" << std::fixed << std::setprecision(0) << theta_deg;

    fs::create_directories(output_dir);
    if (out_cfg.save_png)
        c->SaveAs((output_dir / (fname.str() + ".png")).c_str());

    delete c;
}

// ---------------------------------------------------------------------------
// Output ROOT: rapporto calibrazione
// ---------------------------------------------------------------------------
static void produce_calib_report(
    const std::vector<Homography>& homographies,
    const std::vector<double>& axial_dists,
    const std::vector<std::vector<CalibPoint>>& calib_pts_per_dist,
    const fs::path& output_path) {

    apply_riptide_style();

    int n = static_cast<int>(homographies.size());
    if (n == 0) return;

    int cols = std::min(n, 3);
    int rows = (n + cols - 1) / cols;

    auto* c = new TCanvas("calib_report", "Calibration Report",
                           static_cast<int>(600 * cols),
                           static_cast<int>(500 * rows));
    c->Divide(cols, rows);

    for (int i = 0; i < n; ++i) {
        c->cd(i + 1);
        gPad->SetLeftMargin(0.16f);
        gPad->SetBottomMargin(0.14f);
        gPad->SetGrid();

        const auto& H   = homographies[static_cast<size_t>(i)];
        double d_ax     = (i < static_cast<int>(axial_dists.size())) ? axial_dists[static_cast<size_t>(i)] : 0.0;
        const auto& pts = calib_pts_per_dist[static_cast<size_t>(i)];

        if (pts.empty()) continue;

        // Residui
        std::vector<double> xs, ys, res;
        for (const auto& p : pts) {
            double w  = H.H[6] * p.disp_x + H.H[7] * p.disp_y + H.H[8];
            double px = (H.H[0] * p.disp_x + H.H[1] * p.disp_y + H.H[2]) / w;
            double py = (H.H[3] * p.disp_x + H.H[4] * p.disp_y + H.H[5]) / w;
            xs.push_back(p.sens_x);
            ys.push_back(p.sens_y);
            double r = std::sqrt((px - p.sens_x) * (px - p.sens_x) +
                                  (py - p.sens_y) * (py - p.sens_y));
            res.push_back(r);
        }

        auto* gr = new TGraph(static_cast<int>(xs.size()), xs.data(), ys.data());
        std::ostringstream t;
        t << "d_{ax}=" << std::fixed << std::setprecision(0) << d_ax << " mm;x_{sens} [px];y_{sens} [px]";
        gr->SetTitle(t.str().c_str());
        gr->SetMarkerStyle(20);
        gr->SetMarkerColor(kBlue + 1);
        gr->Draw("AP");

        double rms = H.rms_residual;
        auto* pave = new TPaveText(0.55, 0.70, 0.98, 0.98, "NDC");
        pave->SetFillColor(0);
        pave->SetBorderSize(1);
        pave->SetTextFont(42);
        pave->SetTextSize(0.04f);
        {
            std::ostringstream s;
            s << "n_{pts}=" << H.n_points;
            pave->AddText(s.str().c_str());
        }
        {
            std::ostringstream s;
            s << "RMS=" << std::fixed << std::setprecision(2) << rms << " px";
            pave->AddText(s.str().c_str());
        }
        pave->Draw();
    }

    fs::create_directories(output_path.parent_path());
    c->SaveAs(output_path.c_str());
    delete c;
}

// ---------------------------------------------------------------------------
// Output ROOT: profilo Q(d_ax)
// ---------------------------------------------------------------------------
static void produce_Q_profile(const std::vector<AxialResult>& axial,
                               double Q_sim,
                               const std::string& config_label,
                               const fs::path& output_path) {
    apply_riptide_style();

    std::vector<double> xv, yv, yev;
    for (const auto& ar : axial) {
        if (ar.n_valid_orientations == 0) continue;
        xv.push_back(ar.d_ax_measured_mm);
        yv.push_back(ar.Q_exp);
        yev.push_back(ar.Q_exp_sigma);
    }

    if (xv.empty()) return;

    auto* c = new TCanvas("Q_profile", ("Q profile — " + config_label).c_str(), 800, 600);
    c->SetLeftMargin(0.16f);
    c->SetBottomMargin(0.14f);
    c->SetGrid();

    auto* gr = new TGraphErrors(static_cast<int>(xv.size()),
                                 xv.data(), yv.data(), nullptr, yev.data());
    std::string title = config_label + ";d_{ax} [mm];Q = #chi^{2}/ndof";
    gr->SetTitle(title.c_str());
    gr->SetMarkerStyle(21);
    gr->SetMarkerColor(kBlue + 1);
    gr->SetLineColor(kBlue + 1);
    gr->SetLineWidth(2);
    gr->Draw("APE");

    double x_min = *std::min_element(xv.begin(), xv.end()) - 5.0;
    double x_max = *std::max_element(xv.begin(), xv.end()) + 5.0;

    auto* leg = new TLegend(0.60, 0.75, 0.95, 0.95);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(gr, "Q_{exp} #pm #sigma_{#theta}", "PE");

    if (std::isfinite(Q_sim)) {
        auto* lsim = new TLine(x_min, Q_sim, x_max, Q_sim);
        lsim->SetLineColor(kRed);
        lsim->SetLineWidth(2);
        lsim->SetLineStyle(2);
        lsim->Draw();
        leg->AddEntry(lsim, "Q_{sim}", "L");
    }
    leg->Draw();

    fs::create_directories(output_path.parent_path());
    c->SaveAs(output_path.c_str());
    delete c;
}

// ---------------------------------------------------------------------------
// Output ROOT: pannello summary (good vs bad)
// ---------------------------------------------------------------------------
static void produce_summary_panel(const std::vector<ConfigResult>& results,
                                   const fs::path& output_path) {
    apply_riptide_style();

    auto* c = new TCanvas("summary", "Summary", 1200, 600);
    c->Divide(2, 1);

    // Raccoglie valori
    std::vector<std::string> labels;
    std::vector<double> Q_exps, Q_sims, Rs;
    for (const auto& r : results) {
        labels.push_back(r.config_label);
        Q_exps.push_back(r.Q_exp_global);
        Q_sims.push_back(r.Q_sim);
        Rs.push_back(r.R);
    }

    int n = static_cast<int>(results.size());

    // Pad 1: Q_exp_global e Q_sim
    c->cd(1);
    gPad->SetLeftMargin(0.16f);
    gPad->SetBottomMargin(0.14f);
    gPad->SetGrid();

    auto* gr_exp = new TGraph(n);
    auto* gr_sim = new TGraph(n);
    for (int i = 0; i < n; ++i) {
        gr_exp->SetPoint(i, static_cast<double>(i) + 0.8, Q_exps[static_cast<size_t>(i)]);
        gr_sim->SetPoint(i, static_cast<double>(i) + 1.2, Q_sims[static_cast<size_t>(i)]);
    }
    gr_exp->SetTitle(";configurazione;Q = #chi^{2}/ndof");
    gr_exp->SetMarkerStyle(21); gr_exp->SetMarkerColor(kBlue + 1); gr_exp->SetMarkerSize(1.5f);
    gr_sim->SetMarkerStyle(22); gr_sim->SetMarkerColor(kRed);       gr_sim->SetMarkerSize(1.5f);

    auto* mg1 = new TGraph();
    mg1->SetTitle(";configurazione;Q = #chi^{2}/ndof");
    gr_exp->GetXaxis()->SetLimits(0.0, static_cast<double>(n) + 1.0);
    gr_exp->Draw("AP");
    gr_sim->Draw("P SAME");

    auto* leg1 = new TLegend(0.55, 0.75, 0.95, 0.95);
    leg1->SetBorderSize(0); leg1->SetFillStyle(0);
    leg1->AddEntry(gr_exp, "Q_{exp}", "P");
    leg1->AddEntry(gr_sim, "Q_{sim}", "P");
    leg1->Draw();

    // Pad 2: R = Q_exp / Q_sim
    c->cd(2);
    gPad->SetLeftMargin(0.16f);
    gPad->SetBottomMargin(0.14f);
    gPad->SetGrid();

    auto* gr_R = new TGraph(n);
    for (int i = 0; i < n; ++i)
        gr_R->SetPoint(i, static_cast<double>(i) + 1.0, Rs[static_cast<size_t>(i)]);
    gr_R->SetTitle(";configurazione;R = Q_{exp}/Q_{sim}");
    gr_R->SetMarkerStyle(21); gr_R->SetMarkerColor(kViolet + 1); gr_R->SetMarkerSize(1.5f);
    gr_R->GetXaxis()->SetLimits(0.0, static_cast<double>(n) + 1.0);
    gr_R->Draw("AP");

    // Linea R = 1
    auto* l1 = new TLine(0.0, 1.0, static_cast<double>(n) + 1.0, 1.0);
    l1->SetLineColor(kGray + 2); l1->SetLineWidth(2); l1->SetLineStyle(2);
    l1->Draw();

    delete mg1;

    fs::create_directories(output_path.parent_path());
    c->SaveAs(output_path.c_str());
    delete c;
}

// ---------------------------------------------------------------------------
// Pipeline completa per una configurazione lenti
// ---------------------------------------------------------------------------
ConfigResult analyze_config(const fs::path& data_root,
                             const fs::path& calib_root,
                             const LensConfig& lens,
                             const Exp3Config& cfg,
                             const OutputConfig& out_cfg) {
    ConfigResult result;
    result.config_label = lens.label;
    result.lens         = lens;

    fs::path config_data = data_root / lens.label;
    fs::path trace_out   = out_cfg.output_dir / lens.label;

    size_t n_dist = cfg.axial_distances_nominal_mm.size();

    for (size_t di = 0; di < n_dist; ++di) {
        double d_nom = cfg.axial_distances_nominal_mm[di];
        double d_meas = cfg.axial_distances_measured_mm[di];

        // Formato directory: d150, d165, ...
        std::ostringstream d_str;
        d_str << "d" << std::fixed << std::setprecision(0) << d_nom;

        // Carica omografia calibrata
        std::ostringstream hfile;
        hfile << "homography_d" << std::fixed << std::setprecision(0) << d_nom << "mm.json";
        fs::path hpath = calib_root / hfile.str();

        Homography H;
        bool h_loaded = false;
        try {
            H = load_homography(hpath);
            h_loaded = true;
        } catch (const std::exception& e) {
            if (out_cfg.verbose)
                fprintf(stderr, "[exp3] WARNING: omografia non trovata per d=%.0f: %s\n",
                        d_nom, e.what());
        }

        AxialResult ar;
        ar.d_ax_nominal_mm  = d_nom;
        ar.d_ax_measured_mm = d_meas;

        for (double theta : cfg.orientations_deg) {
            std::ostringstream t_str;
            t_str << "theta" << std::setfill('0') << std::setw(3)
                  << std::fixed << std::setprecision(0) << static_cast<int>(theta);

            fs::path sig_dir = config_data / d_str.str() / t_str.str() / "signal";
            fs::path bg_dir  = config_data / d_str.str() / t_str.str() / "background";

            if (!fs::exists(sig_dir)) {
                if (out_cfg.verbose)
                    fprintf(stderr, "[exp3] WARNING: directory non trovata: %s\n",
                            sig_dir.c_str());
                SingleTraceResult r;
                r.d_ax_nominal_mm  = d_nom;
                r.d_ax_measured_mm = d_meas;
                r.theta_nominal_deg = theta;
                r.valid   = false;
                r.warning = "Directory segnale non trovata";
                ar.traces.push_back(r);
                continue;
            }

            auto data = analyze_single_trace_full(sig_dir, bg_dir, theta,
                                                   d_nom, d_meas, H, cfg);

            if (!data.result.warning.empty() && out_cfg.verbose)
                fprintf(stderr, "[exp3] WARNING (d=%.0f theta=%.0f %s): %s\n",
                        d_nom, theta, lens.label.c_str(), data.result.warning.c_str());

            // Produce PNG inline
            if (h_loaded && out_cfg.save_png && data.result.valid)
                produce_trace_output(data, trace_out, d_nom, theta, out_cfg);

            ar.traces.push_back(data.result);
        }

        // Media chi2 sulle orientazioni valide
        std::vector<double> chi2_vals;
        for (const auto& tr : ar.traces)
            if (tr.valid) chi2_vals.push_back(tr.chi2_ndof);

        ar.n_valid_orientations = static_cast<int>(chi2_vals.size());
        if (!chi2_vals.empty()) {
            double sum = 0.0, sum2 = 0.0;
            for (double v : chi2_vals) { sum += v; sum2 += v * v; }
            ar.Q_exp = sum / static_cast<double>(chi2_vals.size());
            double var = sum2 / static_cast<double>(chi2_vals.size()) - ar.Q_exp * ar.Q_exp;
            ar.Q_exp_sigma = (var > 0.0) ? std::sqrt(var) : 0.0;
        }

        result.axial.push_back(ar);
    }

    // Media pesata su d_ax (peso = n_valid_orientations)
    double sum_w = 0.0, sum_wq = 0.0;
    for (const auto& ar : result.axial) {
        if (ar.n_valid_orientations > 0) {
            double w = static_cast<double>(ar.n_valid_orientations);
            sum_w  += w;
            sum_wq += w * ar.Q_exp;
        }
    }
    result.Q_exp_global = (sum_w > 0.0) ? sum_wq / sum_w : 0.0;

    // Confronto simulazione
    result.Q_sim = load_Q_sim(cfg.q_map_tsv, lens.x1_mm, lens.x2_mm);
    result.R     = (std::isfinite(result.Q_sim) && result.Q_sim > 1e-12)
                   ? result.Q_exp_global / result.Q_sim
                   : std::numeric_limits<double>::quiet_NaN();

    // Produce Q_profile.png
    if (out_cfg.save_png)
        produce_Q_profile(result.axial, result.Q_sim, lens.label,
                           out_cfg.output_dir / lens.label / "Q_profile.png");

    return result;
}

// ---------------------------------------------------------------------------
// Funzione esposta per produrre il calib_report: wrappa produce_calib_report
// ---------------------------------------------------------------------------
void produce_calibration_report(const std::vector<Homography>& homographies,
                                 const std::vector<double>& axial_dists,
                                 const std::vector<std::vector<CalibPoint>>& pts_per_dist,
                                 const fs::path& output_path) {
    produce_calib_report(homographies, axial_dists, pts_per_dist, output_path);
}

void produce_results_summary(const std::vector<ConfigResult>& results,
                              const fs::path& output_path) {
    produce_summary_panel(results, output_path);
}

// ---------------------------------------------------------------------------
// Scrittura TSV
// ---------------------------------------------------------------------------
void write_results_tsv(const std::vector<ConfigResult>& results,
                        const fs::path& path) {
    fs::create_directories(path.parent_path());
    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("write_results_tsv: impossibile aprire " + path.string());

    ofs << std::fixed << std::setprecision(6);
    ofs << "label\tx1_mm\tx2_mm\td_ax_nominal\td_ax_measured\t"
           "theta_nominal\ttheta_measured\tchi2_ndof\tn_valid_slices\tvalid\twarning\n";

    for (const auto& cr : results) {
        for (const auto& ar : cr.axial) {
            for (const auto& tr : ar.traces) {
                ofs << cr.config_label << '\t'
                    << cr.lens.x1_mm << '\t' << cr.lens.x2_mm << '\t'
                    << tr.d_ax_nominal_mm << '\t' << tr.d_ax_measured_mm << '\t'
                    << tr.theta_nominal_deg << '\t' << tr.theta_measured_deg << '\t'
                    << tr.chi2_ndof << '\t' << tr.n_valid_slices << '\t'
                    << (tr.valid ? 1 : 0) << '\t' << tr.warning << '\n';
            }
        }
    }
}

void write_summary_tsv(const std::vector<ConfigResult>& results,
                        const fs::path& path) {
    fs::create_directories(path.parent_path());
    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("write_summary_tsv: impossibile aprire " + path.string());

    ofs << std::fixed << std::setprecision(6);
    ofs << "label\tx1_mm\tx2_mm\tQ_exp_global\tQ_sim\tR\n";
    for (const auto& cr : results) {
        ofs << cr.config_label << '\t'
            << cr.lens.x1_mm << '\t' << cr.lens.x2_mm << '\t'
            << cr.Q_exp_global << '\t' << cr.Q_sim << '\t' << cr.R << '\n';
    }
}

} // namespace riptide::exp3
