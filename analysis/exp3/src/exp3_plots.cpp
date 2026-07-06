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
#include "exp3b_analysis.hpp"
#include "fits_io.hpp"

// ROOT
#include <TAxis.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1F.h>
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
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <vector>

namespace fs = std::filesystem;

namespace {

// ---------------------------------------------------------------------------
// Stile RIPTIDE condiviso
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
    TGaxis::SetMaxDigits(4);
}

static void set_pad_margins() {
    gPad->SetLeftMargin(0.16f);
    gPad->SetBottomMargin(0.14f);
    gPad->SetGrid();
}

static std::string fmt(double v, int prec = 1) {
    std::ostringstream s;
    s << std::fixed << std::setprecision(prec) << v;
    return s.str();
}

// Percentile robusto su TH2D (ignora pixel ~0) per una scala colore che non
// venga dominata da hot pixel/outlier. Adattato da exp1/image_analysis.cpp.
static double get_th2d_percentile(TH2D* h, double percentile) {
    if (!h || percentile <= 0.0 || percentile >= 1.0) return 0.0;
    std::vector<double> vals;
    vals.reserve(static_cast<size_t>(h->GetNbinsX() * h->GetNbinsY()));
    for (int iy = 1; iy <= h->GetNbinsY(); ++iy)
        for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
            double v = h->GetBinContent(ix, iy);
            if (v > 1e-9) vals.push_back(v);
        }
    if (vals.empty()) return 0.0;
    size_t idx = static_cast<size_t>(static_cast<double>(vals.size()) * percentile);
    std::nth_element(vals.begin(), vals.begin() + static_cast<long>(idx), vals.end());
    return vals[idx];
}

// Rebinna un TH2D se troppo grande per il rendering PNG (safety memoria).
// Adattato da exp1/image_analysis.cpp.
static TH2D* rebin_for_display(TH2D* h) {
    if (!h) return nullptr;
    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();
    int rx = (nx > 2000) ? nx / 1024 : 1;
    int ry = (ny > 2000) ? ny / 1024 : 1;
    if (rx <= 1 && ry <= 1) return static_cast<TH2D*>(h->Clone());
    auto* h_rebinned = static_cast<TH2D*>(
        h->Rebin2D(rx, ry, (std::string(h->GetName()) + "_disp").c_str()));
    h_rebinned->Scale(1.0 / (static_cast<double>(rx) * ry));
    return h_rebinned;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// Exp3: Q_exp(r) per d_ax fisso
// ---------------------------------------------------------------------------
namespace riptide::exp3 {

void produce_q_vs_r(const std::vector<QResult>& good,
                    const std::vector<QResult>& bad,
                    double d_ax_mm,
                    const fs::path& output_path) {
    apply_riptide_style();

    // Filtra per d_ax
    auto pick = [&](const std::vector<QResult>& v) {
        std::vector<double> xs, ys;
        for (const auto& r : v)
            if (std::abs(r.d_ax_mm - d_ax_mm) < 0.6) {
                xs.push_back(r.r_mm);
                ys.push_back(r.chi2_ndof);
            }
        return std::make_pair(xs, ys);
    };

    auto [gx, gy] = pick(good);
    auto [bx, by] = pick(bad);

    if (gx.empty() && bx.empty()) return;

    auto* c = new TCanvas("qvsr", ("Q(r) d_{ax}=" + fmt(d_ax_mm) + " mm").c_str(), 800, 600);
    c->SetLeftMargin(0.16f);
    c->SetBottomMargin(0.14f);
    set_pad_margins();

    auto* leg = new TLegend(0.60, 0.75, 0.95, 0.95);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);

    bool first = true;
    auto draw_graph = [&](std::vector<double>& xs, std::vector<double>& ys,
                           int color, int style, const char* label) {
        if (xs.empty()) return;
        auto* gr = new TGraph(static_cast<int>(xs.size()), xs.data(), ys.data());
        gr->SetTitle(";r [mm];Q = #chi^{2}/ndof");
        gr->SetMarkerStyle(style);
        gr->SetMarkerColor(color);
        gr->SetLineColor(color);
        gr->SetLineWidth(2);
        gr->SetMarkerSize(1.2f);
        gr->Draw(first ? "AP" : "P SAME");
        leg->AddEntry(gr, label, "P");
        first = false;
    };

    draw_graph(gx, gy, kBlue + 1, 21, "good");
    draw_graph(bx, by, kRed,      25, "bad");
    leg->Draw();

    auto* pave = new TPaveText(0.16, 0.88, 0.58, 0.98, "NDC");
    pave->SetFillColor(0); pave->SetBorderSize(0); pave->SetTextFont(42);
    pave->SetTextSize(0.035f);
    pave->AddText(("d_{ax} = " + fmt(d_ax_mm) + " mm").c_str());
    pave->Draw();

    fs::create_directories(output_path.parent_path());
    c->SaveAs(output_path.c_str());
    delete c;
}

// ---------------------------------------------------------------------------
// Exp3: Q_exp(d_ax) per r_idx fisso
// ---------------------------------------------------------------------------
void produce_q_vs_dax(const std::vector<QResult>& good,
                      const std::vector<QResult>& bad,
                      int r_idx,
                      const fs::path& output_path) {
    apply_riptide_style();

    auto pick = [&](const std::vector<QResult>& v) {
        std::vector<double> xs, ys;
        for (const auto& r : v)
            if (r.r_idx == r_idx) {
                xs.push_back(r.d_ax_mm);
                ys.push_back(r.chi2_ndof);
            }
        return std::make_pair(xs, ys);
    };

    auto [gx, gy] = pick(good);
    auto [bx, by] = pick(bad);

    if (gx.empty() && bx.empty()) return;

    std::string r_mm_label = (!good.empty() ? fmt(good.front().r_mm) :
                               (!bad.empty() ? fmt(bad.front().r_mm) : "?"));
    for (const auto& r : good) if (r.r_idx == r_idx) { r_mm_label = fmt(r.r_mm); break; }

    auto* c = new TCanvas("qvsdax",
                           ("Q(d_{ax}) r=" + r_mm_label + " mm").c_str(), 800, 600);
    c->SetLeftMargin(0.16f);
    c->SetBottomMargin(0.14f);
    set_pad_margins();

    auto* leg = new TLegend(0.60, 0.75, 0.95, 0.95);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42);

    bool first = true;
    auto draw_graph = [&](std::vector<double>& xs, std::vector<double>& ys,
                           int color, int style, const char* label) {
        if (xs.empty()) return;
        auto* gr = new TGraph(static_cast<int>(xs.size()), xs.data(), ys.data());
        gr->SetTitle(";d_{ax} [mm];Q = #chi^{2}/ndof");
        gr->SetMarkerStyle(style);
        gr->SetMarkerColor(color);
        gr->SetLineColor(color);
        gr->SetLineWidth(2);
        gr->SetMarkerSize(1.2f);
        gr->Draw(first ? "APL" : "PL SAME");
        leg->AddEntry(gr, label, "PL");
        first = false;
    };

    draw_graph(gx, gy, kBlue + 1, 21, "good");
    draw_graph(bx, by, kRed,      25, "bad");
    leg->Draw();

    auto* pave = new TPaveText(0.16, 0.88, 0.60, 0.98, "NDC");
    pave->SetFillColor(0); pave->SetBorderSize(0); pave->SetTextFont(42);
    pave->SetTextSize(0.035f);
    pave->AddText(("r = " + r_mm_label + " mm  (r_{idx}=" +
                   std::to_string(r_idx) + ")").c_str());
    pave->Draw();

    fs::create_directories(output_path.parent_path());
    c->SaveAs(output_path.c_str());
    delete c;
}

// ---------------------------------------------------------------------------
// Exp3: confronto Q_exp(d_ax) vs Q_sim per good e bad  (pannello singolo)
// ---------------------------------------------------------------------------
void produce_q_comparison(const std::vector<QResult>& good,
                          const std::vector<QResult>& bad,
                          const QComparisonConfig& cmp,
                          const fs::path& output_path) {
    apply_riptide_style();

    // Filtra misure valide. Nota: alcune distanze reali (es. d_ax=150mm per
    // "good", d_ax=210mm per "bad") non compaiono nel grafico aggregato
    // perche' quelle specifiche acquisizioni FITS hanno warning=true
    // (n_valid_slices sotto soglia min_valid_slices, o traccia non
    // caratterizzabile/troppo defocused) - filtro di qualita' dato
    // legittimo, non un bug di plotting.
    auto filter_valid = [](const std::vector<QResult>& v) {
        std::vector<QResult> out;
        for (const auto& r : v)
            if (!r.warning && r.n_valid_slices > 0)
                out.push_back(r);
        return out;
    };
    auto gv = filter_valid(good);
    auto bv = filter_valid(bad);
    if (gv.empty() && bv.empty()) return;

    // Aggrega per d_ax: media, min, max di chi2/ndof
    struct DAXStats { double d_ax, mean, lo, hi; };
    auto group_stats = [](const std::vector<QResult>& v) {
        std::map<double, std::vector<double>> by_dax;
        for (const auto& r : v)
            by_dax[r.d_ax_mm].push_back(r.chi2_ndof);
        std::vector<DAXStats> out;
        for (const auto& [d, vals] : by_dax) {
            double s = 0.0;
            for (double q : vals) s += q;
            double mean = s / static_cast<double>(vals.size());
            double lo = *std::min_element(vals.begin(), vals.end());
            double hi = *std::max_element(vals.begin(), vals.end());
            out.push_back({d, mean, lo, hi});
        }
        return out;
    };
    auto gs = group_stats(gv);
    auto bs = group_stats(bv);

    // Q_sim: preferisce i valori fissi da config (banco exp3 a fuoco fisso
    // 180mm), non ricavabili da q_map_tsv perche' quest'ultimo viene da run
    // a fuoco mobile e non e' confrontabile con l'esperimento. Fallback sul
    // lookup nearest-neighbor in q_map_tsv se non impostati.
    double q_sim_good = std::isfinite(cmp.q_sim_good)
        ? cmp.q_sim_good
        : load_Q_sim(cmp.q_map_tsv, cmp.good_x1_mm, cmp.good_x2_mm);
    double q_sim_bad = std::isfinite(cmp.q_sim_bad)
        ? cmp.q_sim_bad
        : load_Q_sim(cmp.q_map_tsv, cmp.bad_x1_mm, cmp.bad_x2_mm);

    const double x_lo = 143.0, x_hi = 217.0;
    // -------------------------------------------------------------------
    // Canvas a pannello singolo
    // -------------------------------------------------------------------
    auto* c = new TCanvas("qcomp", "Q confronto sim vs exp", 820, 620);
    c->SetFillColor(0);
    c->SetLeftMargin(0.16f);
    c->SetRightMargin(0.06f);
    c->SetTopMargin(0.12f);
    c->SetBottomMargin(0.14f);
    c->SetGrid();

    double y_max = 1.5;
    for (const auto& s : bs) y_max = std::max(y_max, s.hi);
    for (const auto& s : gs) y_max = std::max(y_max, s.hi);
    if (std::isfinite(q_sim_bad))  y_max = std::max(y_max, q_sim_bad);
    if (std::isfinite(q_sim_good)) y_max = std::max(y_max, q_sim_good);
    // Margine ampio: lascia una fascia libera in alto per la legenda,
    // cosi' non interseca mai le curve dati.
    y_max *= 1.55;

    auto* frame = static_cast<TH1F*>(
        c->DrawFrame(x_lo, 0.0, x_hi, y_max, ";d_{ax} [mm];Q = #chi^{2}/ndof"));
    frame->GetYaxis()->SetTitleOffset(1.30f);
    frame->GetXaxis()->SetTitleOffset(1.10f);

    auto* leg = new TLegend(0.18, 0.76, 0.92, 0.885);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42);
    leg->SetTextSize(0.032f);
    leg->SetNColumns(2);

    auto draw_series = [&](const std::vector<DAXStats>& sv, int color, int marker,
                            const char* label) {
        if (sv.empty()) return;
        int n = static_cast<int>(sv.size());
        std::vector<double> xs(n), ys(n), exl(n, 0.0), exh(n, 0.0), eyl(n), eyh(n);
        for (int i = 0; i < n; ++i) {
            xs[i]  = sv[i].d_ax;
            ys[i]  = sv[i].mean;
            eyl[i] = sv[i].mean - sv[i].lo;
            eyh[i] = sv[i].hi   - sv[i].mean;
        }
        auto* gr = new TGraphAsymmErrors(n,
            xs.data(), ys.data(),
            exl.data(), exh.data(),
            eyl.data(), eyh.data());
        gr->SetMarkerStyle(marker);
        gr->SetMarkerColor(color);
        gr->SetLineColor(color);
        gr->SetLineWidth(2);
        gr->SetMarkerSize(1.2f);
        gr->SetFillColorAlpha(color, 0.15f);
        gr->Draw("PL E2 SAME");
        leg->AddEntry(gr, label, "PLE");
    };

    draw_series(gs, kBlue + 1, 21,
        Form("exp good (x_{1}=%g, x_{2}=%g mm)", cmp.good_x1_mm, cmp.good_x2_mm));
    draw_series(bs, kRed, 25,
        Form("exp bad  (x_{1}=%g, x_{2}=%g mm)", cmp.bad_x1_mm,  cmp.bad_x2_mm));

    if (std::isfinite(q_sim_good)) {
        auto* l = new TLine(x_lo, q_sim_good, x_hi, q_sim_good);
        l->SetLineColor(kBlue + 1); l->SetLineStyle(7); l->SetLineWidth(2); l->Draw();
        leg->AddEntry(l, Form("sim good  Q_{sim}=%.2f", q_sim_good), "L");
    }
    if (std::isfinite(q_sim_bad)) {
        auto* l = new TLine(x_lo, q_sim_bad, x_hi, q_sim_bad);
        l->SetLineColor(kRed); l->SetLineStyle(7); l->SetLineWidth(2); l->Draw();
        leg->AddEntry(l, Form("sim bad   Q_{sim}=%.2f", q_sim_bad), "L");
    }
    leg->Draw();

    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextSize(0.045f);
    title.SetTextAlign(22);
    title.DrawLatex(0.5, 0.955, "LA4464 / LA4464R  --  Q simulato vs sperimentale");

    fs::create_directories(output_path.parent_path());
    c->SaveAs(output_path.c_str());
    delete c;
}

// ---------------------------------------------------------------------------
// Conversione FITS grezzi -> PNG di ispezione visiva
// ---------------------------------------------------------------------------

// Titolo leggibile a partire dal path relativo, atteso nella forma
// <config>/d<d_ax>/<signal|background>/<timestamp>/<file>. Se la struttura
// non corrisponde (dataset diverso da exp3), ripiega sul nome file grezzo.
static std::string build_fits_preview_title(const fs::path& rel) {
    std::vector<std::string> parts;
    for (const auto& p : rel) parts.push_back(p.string());
    if (parts.size() < 4) return rel.filename().string();

    const std::string& config   = parts[0];
    const std::string& d_dir    = parts[1];
    const std::string& category = parts[2];
    const std::string& filename = parts.back();

    if (d_dir.size() < 2 || d_dir[0] != 'd') return rel.filename().string();
    double d_ax_mm = 0.0;
    try {
        d_ax_mm = std::stod(d_dir.substr(1));
    } catch (...) {
        return rel.filename().string();
    }

    std::string exposure_s, frame_no;
    std::stringstream ss(filename);
    std::string token;
    while (std::getline(ss, token, '_')) {
        if (exposure_s.empty() && token.size() > 3 &&
            token.substr(token.size() - 3) == "sec")
            exposure_s = token.substr(0, token.size() - 3);
        if (token.rfind("frame", 0) == 0) {
            std::size_t dot = token.find('.');
            std::string digits = token.substr(5, dot == std::string::npos ? std::string::npos
                                                                            : dot - 5);
            try { frame_no = std::to_string(std::stoi(digits)); } catch (...) {}
        }
    }

    std::ostringstream out;
    out << config << "  --  d_{ax} = " << std::fixed << std::setprecision(0) << d_ax_mm
        << " mm  --  " << category;
    if (!exposure_s.empty()) out << "  --  esp. " << exposure_s << " s";
    if (!frame_no.empty())   out << "  --  frame " << frame_no;
    return out.str();
}

int produce_fits_previews(const fs::path& data_dir, const fs::path& output_dir) {
    if (!fs::exists(data_dir))
        throw std::runtime_error("produce_fits_previews: data_dir non trovata: " + data_dir.string());

    apply_riptide_style();
    gStyle->SetPalette(kViridis);

    const fs::path preview_root = output_dir / "fits_preview";
    int n_written = 0;

    for (const auto& entry : fs::recursive_directory_iterator(data_dir)) {
        if (!entry.is_regular_file()) continue;
        std::string ext = entry.path().extension().string();
        if (ext != ".fit" && ext != ".fits" && ext != ".fts") continue;

        riptide::fits::FitsFrame frame;
        try {
            frame = riptide::fits::read_fits(entry.path());
        } catch (const std::exception& e) {
            std::cerr << "[fits_preview] WARNING: impossibile leggere " << entry.path()
                      << ": " << e.what() << "\n";
            continue;
        }

        auto* h = new TH2D("h_fits", entry.path().stem().string().c_str(),
                           frame.width(), 0.0, static_cast<double>(frame.width()),
                           frame.height(), 0.0, static_cast<double>(frame.height()));
        for (int y = 0; y < frame.height(); ++y)
            for (int x = 0; x < frame.width(); ++x)
                h->SetBinContent(x + 1, y + 1, frame.pixel(x, y));
        h->GetXaxis()->SetTitle("Pixel X [px]");
        h->GetYaxis()->SetTitle("Pixel Y [px]");
        h->GetXaxis()->CenterTitle(true);
        h->GetYaxis()->CenterTitle(true);
        h->GetXaxis()->SetTitleOffset(1.20f);
        h->GetYaxis()->SetTitleOffset(1.55f);

        TH2D* h_disp = rebin_for_display(h);
        double vmin = get_th2d_percentile(h_disp, 0.01);
        double vmax = get_th2d_percentile(h_disp, 0.99);
        if (vmax <= vmin) { vmin = 0.0; vmax = h_disp->GetMaximum(); }
        h_disp->SetTitle("");
        h_disp->SetMinimum(vmin);
        h_disp->SetMaximum(vmax);
        h_disp->GetZaxis()->SetTitle("ADU");
        h_disp->GetZaxis()->CenterTitle(true);
        h_disp->GetZaxis()->SetTitleOffset(1.35f);
        h_disp->GetZaxis()->SetTitleSize(0.045f);
        h_disp->GetZaxis()->SetLabelSize(0.032f);

        auto* c = new TCanvas("c_fits_preview", "", 900, 700);
        c->SetLeftMargin(0.16f);
        c->SetRightMargin(0.20f);
        c->SetTopMargin(0.13f);
        c->SetBottomMargin(0.14f);
        h_disp->Draw("COLZ");

        fs::path rel = fs::relative(entry.path(), data_dir);
        TLatex t;
        t.SetNDC();
        t.SetTextFont(42);
        t.SetTextSize(0.032f);
        t.SetTextAlign(22);
        t.DrawLatex(0.53, 0.955, build_fits_preview_title(rel).c_str());

        fs::path png_path = preview_root / rel;
        png_path.replace_extension(".png");
        fs::create_directories(png_path.parent_path());
        c->SaveAs(png_path.c_str());

        delete c;
        delete h_disp;
        delete h;
        ++n_written;
    }

    return n_written;
}

} // namespace riptide::exp3

// ---------------------------------------------------------------------------
// Exp3b plots
// ---------------------------------------------------------------------------
namespace riptide::exp3b {

void produce_M_profile(const MResult& result,
                       const std::string& config_label,
                       const fs::path& output_path) {
    apply_riptide_style();
    if (result.n_dots < 2) return;

    int n = result.n_dots;
    // Ricostruisce y_disp_mm dai dati: non disponibile direttamente in MResult,
    // usiamo r_mm come y_sens_mm e ricalcoliamo i residui dal fit
    std::vector<double> y_s = result.r_mm; // y_sens_mm
    // y_disp_mm non è salvato in MResult; lo ricostruiamo invertendo il fit
    // y_disp = (y_sens - q_offset) / M_global
    double M = result.M_global;
    double q = result.q_offset;

    std::vector<double> y_d(static_cast<size_t>(n)), res(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        y_d[static_cast<size_t>(i)] = (std::abs(M) > 1e-9)
            ? (y_s[static_cast<size_t>(i)] - q) / M
            : static_cast<double>(i);
        res[static_cast<size_t>(i)] = y_s[static_cast<size_t>(i)] -
                                       M * y_d[static_cast<size_t>(i)] - q;
    }

    auto* c = new TCanvas("Mprof",
                           ("M profile  " + config_label + "  d_{ax}=" +
                            fmt(result.d_ax_mm) + " mm").c_str(), 800, 900);
    c->Divide(1, 2);

    // Pad 1: y_sens vs y_disp + fit
    c->cd(1);
    gPad->SetLeftMargin(0.16f);
    gPad->SetBottomMargin(0.14f);
    gPad->SetGrid();

    auto* gr = new TGraph(n, y_d.data(), y_s.data());
    gr->SetTitle(";y_{disp} [mm];y_{sens} [mm]");
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue + 1);
    gr->SetMarkerSize(0.8f);
    gr->Draw("AP");

    double xlo = *std::min_element(y_d.begin(), y_d.end());
    double xhi = *std::max_element(y_d.begin(), y_d.end());
    auto* lfit = new TLine(xlo, M * xlo + q, xhi, M * xhi + q);
    lfit->SetLineColor(kRed); lfit->SetLineWidth(2); lfit->Draw();

    auto* pave = new TPaveText(0.55, 0.15, 0.95, 0.45, "NDC");
    pave->SetFillColor(0); pave->SetBorderSize(1); pave->SetTextFont(42);
    pave->SetTextSize(0.04f);
    pave->AddText(("M = " + fmt(M, 4)).c_str());
    pave->AddText(("#chi^{2}/ndof = " + fmt(result.chi2_ndof, 3)).c_str());
    pave->AddText(("n_{dots} = " + std::to_string(n)).c_str());
    pave->Draw();

    // Pad 2: residui
    c->cd(2);
    gPad->SetLeftMargin(0.16f);
    gPad->SetBottomMargin(0.14f);
    gPad->SetGrid();

    auto* gr_res = new TGraph(n, y_d.data(), res.data());
    gr_res->SetTitle(";y_{disp} [mm];residuo [mm]");
    gr_res->SetMarkerStyle(20);
    gr_res->SetMarkerColor(kBlue + 1);
    gr_res->SetMarkerSize(0.8f);
    gr_res->Draw("AP");

    auto* lzero = new TLine(xlo, 0.0, xhi, 0.0);
    lzero->SetLineColor(kGray + 2); lzero->SetLineWidth(2);
    lzero->SetLineStyle(2); lzero->Draw();

    fs::create_directories(output_path.parent_path());
    c->SaveAs(output_path.c_str());
    delete c;
}

void produce_M_vs_dax(const std::vector<MResult>& good,
                      const std::vector<MResult>& bad,
                      const fs::path& output_path) {
    apply_riptide_style();

    auto pick = [](const std::vector<MResult>& v) {
        std::vector<double> xs, ys;
        for (const auto& r : v)
            if (std::isfinite(r.M_global)) {
                xs.push_back(r.d_ax_mm);
                ys.push_back(r.M_global);
            }
        return std::make_pair(xs, ys);
    };

    auto [gx, gy] = pick(good);
    auto [bx, by] = pick(bad);
    if (gx.empty() && bx.empty()) return;

    auto* c = new TCanvas("Mvsdax", "M_global vs d_{ax}", 800, 600);
    c->SetLeftMargin(0.16f);
    c->SetBottomMargin(0.14f);
    set_pad_margins();

    auto* leg = new TLegend(0.60, 0.75, 0.95, 0.95);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42);

    bool first = true;
    auto draw_g = [&](std::vector<double>& xs, std::vector<double>& ys,
                       int color, int style, const char* label) {
        if (xs.empty()) return;
        auto* gr = new TGraph(static_cast<int>(xs.size()), xs.data(), ys.data());
        gr->SetTitle(";d_{ax} [mm];M_{global}");
        gr->SetMarkerStyle(style);
        gr->SetMarkerColor(color);
        gr->SetLineColor(color);
        gr->SetLineWidth(2);
        gr->SetMarkerSize(1.2f);
        gr->Draw(first ? "APL" : "PL SAME");
        leg->AddEntry(gr, label, "PL");
        first = false;
    };

    draw_g(gx, gy, kBlue + 1, 21, "good");
    draw_g(bx, by, kRed,      25, "bad");
    leg->Draw();

    fs::create_directories(output_path.parent_path());
    c->SaveAs(output_path.c_str());
    delete c;
}

void produce_M_nonlinearity(const MResult& result,
                             const std::string& config_label,
                             const fs::path& output_path) {
    apply_riptide_style();
    int n_local = static_cast<int>(result.M_local.size());
    if (n_local < 2) return;

    // Ascissa: r_mm dei punti intermedi
    std::vector<double> x_mid(static_cast<size_t>(n_local));
    const auto& r = result.r_mm;
    for (int i = 0; i < n_local; ++i)
        x_mid[static_cast<size_t>(i)] = (r[static_cast<size_t>(i)] +
                                          r[static_cast<size_t>(i + 1)]) * 0.5;

    std::vector<double> ml = result.M_local;

    // Media e deviazione standard di M_local
    double mean_ml = 0.0;
    for (double v : ml) mean_ml += v;
    mean_ml /= static_cast<double>(n_local);
    double var_ml = 0.0;
    for (double v : ml) var_ml += (v - mean_ml) * (v - mean_ml);
    double std_ml = (n_local > 1) ? std::sqrt(var_ml / static_cast<double>(n_local - 1)) : 0.0;

    auto* c = new TCanvas("Mnl",
                           ("M_{local}(r)  " + config_label +
                            "  d_{ax}=" + fmt(result.d_ax_mm) + " mm").c_str(),
                           800, 600);
    c->SetLeftMargin(0.16f);
    c->SetBottomMargin(0.14f);
    set_pad_margins();

    auto* gr = new TGraph(n_local, x_mid.data(), ml.data());
    gr->SetTitle(";r [mm];M_{local}");
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue + 1);
    gr->SetLineColor(kBlue + 1);
    gr->SetLineWidth(2);
    gr->SetMarkerSize(1.0f);
    gr->Draw("APL");

    double x_lo = *std::min_element(x_mid.begin(), x_mid.end());
    double x_hi = *std::max_element(x_mid.begin(), x_mid.end());

    // Banda ±1σ
    auto* l_mean  = new TLine(x_lo, mean_ml, x_hi, mean_ml);
    l_mean->SetLineColor(kGray + 2); l_mean->SetLineWidth(2);
    l_mean->SetLineStyle(2); l_mean->Draw();

    auto* l_up = new TLine(x_lo, mean_ml + std_ml, x_hi, mean_ml + std_ml);
    l_up->SetLineColor(kGray + 1); l_up->SetLineStyle(3); l_up->SetLineWidth(1);
    l_up->Draw();

    auto* l_dn = new TLine(x_lo, mean_ml - std_ml, x_hi, mean_ml - std_ml);
    l_dn->SetLineColor(kGray + 1); l_dn->SetLineStyle(3); l_dn->SetLineWidth(1);
    l_dn->Draw();

    auto* pave = new TPaveText(0.55, 0.75, 0.95, 0.95, "NDC");
    pave->SetFillColor(0); pave->SetBorderSize(1); pave->SetTextFont(42);
    pave->SetTextSize(0.04f);
    pave->AddText(("#bar{M} = " + fmt(mean_ml, 4)).c_str());
    pave->AddText(("#sigma_{M} = " + fmt(std_ml, 4)).c_str());
    pave->Draw();

    fs::create_directories(output_path.parent_path());
    c->SaveAs(output_path.c_str());
    delete c;
}

} // namespace riptide::exp3b
