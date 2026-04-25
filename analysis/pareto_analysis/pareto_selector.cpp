/*
 * pareto_selector — Selezione ottimale multi-criterio tramite analisi del
 *                   fronte di Pareto e funzione di merito scalare pesata.
 *
 * Aggrega: events.root (η), q_map.tsv (Q), chi2_map.tsv (χ²),
 *          dof_map.tsv (DoF, M, M_abs_err, x_focus)
 */

#include "pareto_core.hpp"

#include <TROOT.h>
#include <TAxis.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace riptide::pareto;

// ── CLI ──────────────────────────────────────────────────────────────────────

struct CliConfig {
    std::string events_path  = "output/optimization/events.root";
    std::string qmap_path    = "output/psf_analysis/q_map.tsv";
    std::string chi2map_path = "output/psf_analysis/chi2_map.tsv";
    std::string dofmap_path  = "output/dof_analysis/dof_map.tsv";
    std::string output_path  = "output/pareto_analysis/pareto_plot.png";
    std::string tsv_path     = "output/pareto_analysis/pareto_results.tsv";
    std::string lens75_id    = "";
    std::string lens60_id    = "";
    FilterConfig fc;
    WeightConfig wc;
};

static void print_help() {
    std::cout <<
        "Uso: pareto_selector [opzioni]\n\n"
        "  --events   PATH   File ROOT efficienza  [output/optimization/events.root]\n"
        "  --qmap     PATH   TSV qualità Q         [output/psf_analysis/q_map.tsv]\n"
        "  --chi2map  PATH   TSV linearità chi2    [output/psf_analysis/chi2_map.tsv]\n"
        "  --dofmap   PATH   TSV DoF               [output/dof_analysis/dof_map.tsv]\n"
        "  --output   PATH   PNG output            [output/pareto_analysis/pareto_plot.png]\n"
        "  --tsv      PATH   TSV risultati         [output/pareto_analysis/pareto_results.tsv]\n"
        "  --eta-frac  F     Soglia eta/eta_max    [0.75]\n"
        "  --x-det     F     Posizione detector mm [180.0]\n"
        "  --focus-tol F     Tolleranza fuoco mm   [15.0]\n"
        "  --dof-min   F     DoF minima mm (0=off) [0.0]\n"
        "  --w-eta     F     Peso eta in Mtot      [0.35]\n"
        "  --w-Q       F     Peso Q in Mtot        [0.40]\n"
        "  --w-dof     F     Peso DoF in Mtot      [0.15]\n"
        "  --w-M       F     Peso M in Mtot        [0.10]\n"
        "  --lens75-id STR   Filtra lente L1       [\"\"]\n"
        "  --lens60-id STR   Filtra lente L2       [\"\"]\n"
        "  --help            Stampa uso ed esce\n";
}

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
        if      (key == "--events")    cfg.events_path  = next();
        else if (key == "--qmap")      cfg.qmap_path    = next();
        else if (key == "--chi2map")   cfg.chi2map_path = next();
        else if (key == "--dofmap")    cfg.dofmap_path  = next();
        else if (key == "--output")    cfg.output_path  = next();
        else if (key == "--tsv")       cfg.tsv_path     = next();
        else if (key == "--eta-frac")  cfg.fc.eta_frac  = std::stod(next());
        else if (key == "--x-det")     cfg.fc.x_det     = std::stod(next());
        else if (key == "--focus-tol") cfg.fc.focus_tol = std::stod(next());
        else if (key == "--dof-min")   cfg.fc.dof_min   = std::stod(next());
        else if (key == "--w-eta")     cfg.wc.w_eta     = std::stod(next());
        else if (key == "--w-Q")       cfg.wc.w_Q       = std::stod(next());
        else if (key == "--w-dof")     cfg.wc.w_dof     = std::stod(next());
        else if (key == "--w-M")       cfg.wc.w_M       = std::stod(next());
        else if (key == "--lens75-id") cfg.lens75_id    = next();
        else if (key == "--lens60-id") cfg.lens60_id    = next();
        else if (key == "--help")      { print_help(); std::exit(0); }
        else {
            std::cerr << "Opzione sconosciuta: " << key << "\n";
            std::exit(1);
        }
    }
    return cfg;
}

// ── TSV loader ───────────────────────────────────────────────────────────────

// Ritorna vector<map<string,string>>, una mappa per ogni riga dati (escluso header).
static std::vector<std::map<std::string, std::string>> load_tsv(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "[ERRORE] Impossibile aprire TSV: " << path << "\n";
        std::exit(1);
    }
    std::vector<std::map<std::string, std::string>> rows;
    std::string line;
    std::vector<std::string> headers;
    bool first = true;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::vector<std::string> fields;
        std::istringstream ss(line);
        std::string tok;
        while (std::getline(ss, tok, '\t'))
            fields.push_back(tok);
        if (first) {
            headers = fields;
            first   = false;
            continue;
        }
        std::map<std::string, std::string> row;
        for (size_t i = 0; i < headers.size() && i < fields.size(); ++i)
            row[headers[i]] = fields[i];
        rows.push_back(std::move(row));
    }
    return rows;
}

// ── ROOT style (da q_map.cpp) ────────────────────────────────────────────────

static void apply_style() {
    gStyle->Reset();
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetTitleFont(42, "");
    gStyle->SetStatFont(42);
    gStyle->SetTextSize(0.040);
    gStyle->SetLabelSize(0.038, "XYZ");
    gStyle->SetTitleSize(0.044, "XYZ");
    gStyle->SetTitleSize(0.046, "");
    gStyle->SetTitleOffset(1.55, "Y");
    gStyle->SetTitleOffset(1.20, "X");
    gStyle->SetTickLength(0.018, "X");
    gStyle->SetTickLength(0.018, "Y");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetGridColor(kGray + 1);
    gStyle->SetGridStyle(3);
    gStyle->SetGridWidth(1);
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
}

static void draw_colorbar(TPad* pad_cb, double vmin, double vmax, const std::string& title) {
    pad_cb->cd();
    pad_cb->Range(0.0, 0.0, 1.0, 1.0);

    const int NB   = 255;
    double cb_x0   = pad_cb->GetLeftMargin();
    double cb_x1   = 1.0 - pad_cb->GetRightMargin();
    double cb_y0   = pad_cb->GetBottomMargin();
    double cb_y1   = 1.0 - pad_cb->GetTopMargin();

    for (int i = 0; i < NB; ++i) {
        double f0  = static_cast<double>(i) / NB;
        double f1  = static_cast<double>(i + 1) / NB;
        double yb0 = cb_y0 + f0 * (cb_y1 - cb_y0);
        double yb1 = cb_y0 + f1 * (cb_y1 - cb_y0);
        TBox* box  = new TBox(cb_x0, yb0, cb_x1, yb1);
        box->SetFillColor(gStyle->GetColorPalette(i));
        box->SetLineWidth(0);
        box->Draw();
    }

    TGaxis* ax = new TGaxis(cb_x1, cb_y0, cb_x1, cb_y1, vmin, vmax, 505, "+L");
    ax->SetLabelFont(42);
    ax->SetLabelSize(0.18);
    ax->SetTickSize(0.35);
    ax->SetLabelOffset(0.03);
    ax->SetTitle(title.c_str());
    ax->SetTitleFont(42);
    ax->SetTitleSize(0.20);
    ax->SetTitleOffset(0.55);
    ax->Draw();
}

// ── main ─────────────────────────────────────────────────────────────────────

int main(int argc, char** argv) {
    // Batch mode prima di qualsiasi ROOT call
    gROOT->SetBatch(true);

    auto cli = parse_args(argc, argv);

    // ── 1. Carica events.root ────────────────────────────────────────────────
    TFile root_file(cli.events_path.c_str(), "READ");
    if (!root_file.IsOpen()) {
        std::cerr << "[ERRORE] Impossibile aprire: " << cli.events_path << "\n";
        return 1;
    }

    TTree* tree_config = (TTree*)root_file.Get("Configurations");
    TTree* tree_eff    = (TTree*)root_file.Get("Efficiency");
    if (!tree_config || !tree_eff) {
        std::cerr << "[ERRORE] TTree 'Configurations' o 'Efficiency' non trovato in "
                  << cli.events_path << "\n";
        return 1;
    }

    double x1_val = 0.0, x2_val = 0.0;
    int  cfg_id        = 0;
    int  n_photons_val = 0;
    int  n_hits_val    = 0;
    char l1_buf[256]   = {};
    char l2_buf[256]   = {};

    tree_config->SetBranchAddress("x1",        &x1_val);
    tree_config->SetBranchAddress("x2",        &x2_val);
    tree_config->SetBranchAddress("config_id", &cfg_id);
    const bool has_lens = (tree_config->GetBranch("lens75_id") != nullptr);
    if (has_lens) {
        tree_config->SetBranchAddress("lens75_id", l1_buf);
        tree_config->SetBranchAddress("lens60_id", l2_buf);
    }

    tree_eff->SetBranchAddress("config_id", &cfg_id);
    tree_eff->SetBranchAddress("n_photons",  &n_photons_val);
    tree_eff->SetBranchAddress("n_hits",     &n_hits_val);

    struct RootCfg { double x1, x2; std::string l1, l2; };
    std::map<int, RootCfg> cfg_map;

    for (Long64_t i = 0; i < tree_config->GetEntries(); ++i) {
        tree_config->GetEntry(i);
        std::string s1 = has_lens ? l1_buf : "default";
        std::string s2 = has_lens ? l2_buf : "default";
        if (!cli.lens75_id.empty() && s1 != cli.lens75_id) continue;
        if (!cli.lens60_id.empty() && s2 != cli.lens60_id) continue;
        cfg_map[cfg_id] = {x1_val, x2_val, s1, s2};
    }

    std::map<int, double> eta_sum;
    std::map<int, int>    eta_cnt;
    for (Long64_t i = 0; i < tree_eff->GetEntries(); ++i) {
        tree_eff->GetEntry(i);
        if (!cfg_map.count(cfg_id) || n_photons_val <= 0) continue;
        eta_sum[cfg_id] += static_cast<double>(n_hits_val) / n_photons_val;
        eta_cnt[cfg_id]++;
    }

    const int n_events_root = (int)cfg_map.size();

    // ── 2. Carica TSV ────────────────────────────────────────────────────────
    auto qmap_rows    = load_tsv(cli.qmap_path);
    auto chi2map_rows = load_tsv(cli.chi2map_path);
    auto dofmap_rows  = load_tsv(cli.dofmap_path);

    // ── 3. Join inner su (x1, x2) ────────────────────────────────────────────
    int n_found_q    = 0;
    int n_found_chi2 = 0;
    int n_found_dof  = 0;
    int n_joined     = 0;

    std::vector<ConfigData> configs;
    configs.reserve(cfg_map.size());

    for (const auto& [id, rc] : cfg_map) {
        double eta = 0.0;
        if (eta_cnt.count(id) && eta_cnt.at(id) > 0)
            eta = eta_sum.at(id) / eta_cnt.at(id);

        // Cerca in q_map.tsv
        const std::map<std::string, std::string>* qrow = nullptr;
        for (const auto& r : qmap_rows) {
            try {
                if (coords_match(rc.x1, rc.x2, std::stod(r.at("x1")), std::stod(r.at("x2")))) {
                    qrow = &r;
                    break;
                }
            } catch (...) {}
        }

        // Cerca in chi2_map.tsv
        const std::map<std::string, std::string>* chi2row = nullptr;
        for (const auto& r : chi2map_rows) {
            try {
                if (coords_match(rc.x1, rc.x2, std::stod(r.at("x1")), std::stod(r.at("x2")))) {
                    chi2row = &r;
                    break;
                }
            } catch (...) {}
        }

        // Cerca in dof_map.tsv
        const std::map<std::string, std::string>* dofrow = nullptr;
        for (const auto& r : dofmap_rows) {
            try {
                if (coords_match(rc.x1, rc.x2, std::stod(r.at("x1")), std::stod(r.at("x2")))) {
                    dofrow = &r;
                    break;
                }
            } catch (...) {}
        }

        // Contatori indipendenti per il riepilogo join
        if (qrow)    ++n_found_q;
        if (chi2row) ++n_found_chi2;
        if (dofrow)  ++n_found_dof;

        if (!qrow) {
            std::cerr << "[WARN] (x1=" << rc.x1 << ", x2=" << rc.x2
                      << "): non trovata in " << cli.qmap_path << "\n";
            continue;
        }
        if (!chi2row) {
            std::cerr << "[WARN] (x1=" << rc.x1 << ", x2=" << rc.x2
                      << "): non trovata in " << cli.chi2map_path << "\n";
            continue;
        }
        if (!dofrow) {
            std::cerr << "[WARN] (x1=" << rc.x1 << ", x2=" << rc.x2
                      << "): non trovata in " << cli.dofmap_path << "\n";
            continue;
        }

        // Validità: salta configurazioni marcate come non valide
        try { if (qrow->at("config_valid") == "0") continue; } catch (...) {}
        try { if (chi2row->at("valid") == "0")      continue; } catch (...) {}

        ConfigData cd;
        cd.x1  = rc.x1;
        cd.x2  = rc.x2;
        cd.eta = eta;

        try { cd.Q         = std::stod(qrow->at("metric")); }       catch (...) { continue; }
        try { cd.chi2      = std::stod(chi2row->at("metric")); }      catch (...) { continue; }
        try { cd.DoF       = std::stod(dofrow->at("dof")); }         catch (...) { continue; }
        try { cd.M         = std::stod(dofrow->at("M")); }           catch (...) { continue; }
        try { cd.M_abs_err = std::stod(dofrow->at("M_abs_err")); }   catch (...) { continue; }
        try { cd.x_focus   = std::stod(dofrow->at("x_focus")); }     catch (...) { continue; }

        // Salta valori non finiti (es. "NaN" nel TSV per config non valide)
        if (!std::isfinite(cd.Q) || !std::isfinite(cd.chi2) ||
            !std::isfinite(cd.DoF) || !std::isfinite(cd.M_abs_err))
            continue;

        configs.push_back(cd);
        ++n_joined;
    }

    std::cout << "[JOIN] Configurazioni in events.root:   " << n_events_root  << "\n"
              << "[JOIN] Trovate in q_map.tsv:            " << n_found_q      << "\n"
              << "[JOIN] Trovate in chi2_map.tsv:         " << n_found_chi2   << "\n"
              << "[JOIN] Trovate in dof_map.tsv:          " << n_found_dof    << "\n"
              << "[JOIN] Configurazioni complete (join):  " << n_joined       << "\n\n";

    if (configs.empty()) {
        std::cerr << "[ERRORE] Join vuoto: nessuna configurazione completa trovata.\n";
        return 1;
    }

    // ── 4. Pipeline filtri ───────────────────────────────────────────────────
    std::cout << "[FILTER] Totale prima dei filtri:       " << configs.size() << "\n";

    configs = apply_eta_filter(configs, cli.fc);
    std::cout << "[FILTER] Dopo filtro eta (>=" << cli.fc.eta_frac << "*eta_max):   "
              << configs.size() << "\n";

    configs = apply_focus_filter(configs, cli.fc);
    std::cout << "[FILTER] Dopo filtro fuoco (|x*-" << cli.fc.x_det << "|<="
              << cli.fc.focus_tol << "): " << configs.size() << "\n";

    if (cli.fc.dof_min > 0.0) {
        configs = apply_dof_filter(configs, cli.fc);
        std::cout << "[FILTER] Dopo filtro DoF (>=" << cli.fc.dof_min << " mm):     "
                  << configs.size() << "\n";
    }

    std::cout << "\n";

    if (configs.empty()) {
        std::cerr << "[WARN] Nessuna configurazione sopravvissuta ai filtri. Nessun output prodotto.\n";
        return 0;
    }

    // ── 5. Calcolo Mtot e fronte di Pareto ──────────────────────────────────
    compute_mtot(configs, cli.wc);
    compute_pareto_front(configs);

    // Fronte ordinato per rank
    std::vector<ConfigData*> front;
    for (auto& c : configs)
        if (c.on_pareto) front.push_back(&c);
    std::sort(front.begin(), front.end(),
              [](const ConfigData* a, const ConfigData* b) {
                  return a->pareto_rank < b->pareto_rank;
              });

    std::cout << "[PARETO] Configurazioni sul fronte:     " << front.size() << "\n";
    if (!front.empty()) {
        const ConfigData& best = *front[0];
        std::cout << "[PARETO] Configurazione raccomandata: x1=" << best.x1
                  << " x2=" << best.x2 << "\n"
                  << "         eta=" << std::fixed << std::setprecision(3) << best.eta
                  << "  Q=" << best.Q
                  << "  chi2=" << best.chi2
                  << "  DoF=" << best.DoF
                  << "  M=" << best.M
                  << "  Mtot=" << best.Mtot << "\n\n";
    }

    // Normalizzatori per assi e TSV
    double eta_max = 0.0, Q_max = 0.0, DoF_max = 0.0;
    for (const auto& c : configs) {
        eta_max = std::max(eta_max, c.eta);
        Q_max   = std::max(Q_max,   c.Q);
        DoF_max = std::max(DoF_max,  c.DoF);
    }
    if (eta_max <= 0.0) eta_max = 1.0;
    if (Q_max   <= 0.0) Q_max   = 1.0;
    if (DoF_max <= 0.0) DoF_max = 1.0;

    // ── 6. TSV output ─────────────────────────────────────────────────────────
    {
        std::filesystem::create_directories(
            std::filesystem::path(cli.tsv_path).parent_path());
        std::ofstream tsv(cli.tsv_path);
        if (!tsv.is_open()) {
            std::cerr << "[ERRORE] Impossibile aprire TSV output: " << cli.tsv_path << "\n";
        } else {
            tsv << "x1\tx2\teta\teta_norm\tQ\tchi2\tDoF\tM\tM_abs_err\t"
                   "x_focus\ton_pareto\tMtot\tpareto_rank\n";
            for (const auto& c : configs) {
                tsv << std::fixed << std::setprecision(6)
                    << c.x1       << "\t" << c.x2         << "\t"
                    << c.eta      << "\t" << c.eta / eta_max << "\t"
                    << c.Q        << "\t" << c.chi2         << "\t"
                    << c.DoF      << "\t" << c.M            << "\t"
                    << c.M_abs_err << "\t" << c.x_focus     << "\t"
                    << (c.on_pareto ? 1 : 0) << "\t"
                    << c.Mtot     << "\t" << c.pareto_rank  << "\n";
            }
        }
    }

    // ── 7. Grafico ROOT ───────────────────────────────────────────────────────
    apply_style();

    TCanvas* canvas = new TCanvas("pareto_canvas", "Pareto Selector", 1200, 900);
    canvas->SetLeftMargin(0.0);
    canvas->SetRightMargin(0.0);
    canvas->SetTopMargin(0.0);
    canvas->SetBottomMargin(0.0);

    TPad* pad_top = new TPad("pad_top", "", 0.00, 0.30, 0.88, 1.0);
    TPad* pad_cb  = new TPad("pad_cb",  "", 0.88, 0.30, 0.96, 1.0);
    TPad* pad_bot = new TPad("pad_bot", "", 0.00, 0.00, 1.00, 0.30);

    pad_top->SetLeftMargin(0.16);
    pad_top->SetRightMargin(0.05);
    pad_top->SetTopMargin(0.08);
    pad_top->SetBottomMargin(0.14);
    pad_top->SetGridx();
    pad_top->SetGridy();

    pad_cb->SetLeftMargin(0.25);
    pad_cb->SetRightMargin(0.30);
    pad_cb->SetTopMargin(0.08);
    pad_cb->SetBottomMargin(0.14);

    pad_bot->SetLeftMargin(0.02);
    pad_bot->SetRightMargin(0.02);
    pad_bot->SetTopMargin(0.05);
    pad_bot->SetBottomMargin(0.05);

    canvas->cd();
    pad_top->Draw();
    pad_cb->Draw();
    pad_bot->Draw();

    // ── Pad superiore: scatter plot ───────────────────────────────────────────
    pad_top->cd();

    // Calcola range assi
    double xmin_p = 1e9, xmax_p = -1e9, ymin_p = 1e9, ymax_p = -1e9;
    for (const auto& c : configs) {
        if (c.Q <= 0.0) continue;
        double xv = c.eta / eta_max;
        double yv = Q_max  / c.Q;
        xmin_p = std::min(xmin_p, xv);
        xmax_p = std::max(xmax_p, xv);
        ymin_p = std::min(ymin_p, yv);
        ymax_p = std::max(ymax_p, yv);
    }
    if (xmin_p > xmax_p) { xmin_p = 0.0; xmax_p = 1.0; }
    if (ymin_p > ymax_p) { ymin_p = 0.5; ymax_p = 1.5; }
    double xpad = 0.08 * (xmax_p - xmin_p);
    double ypad = 0.10 * (ymax_p - ymin_p);
    double xlo = std::max(0.0, xmin_p - xpad);
    double xhi = std::min(1.05, xmax_p + xpad);
    double ylo = std::max(0.0, ymin_p - ypad);
    double yhi = ymax_p + ypad;

    // Frame vuoto per assi
    TH2D* frame = new TH2D("frame", "", 100, xlo, xhi, 100, ylo, yhi);
    frame->GetXaxis()->SetTitle("#eta / #eta_{max}");
    frame->GetYaxis()->SetTitle("Q_{max}/Q");
    frame->GetXaxis()->SetTitleOffset(1.2);
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->Draw("AXIS");

    // Punti di background (non sul fronte) — marker size proporzionale a DoF
    std::vector<TGraph*> g_bg_list;
    for (const auto& c : configs) {
        if (c.on_pareto || c.Q <= 0.0) continue;
        double sz = 0.6 + 1.4 * (c.DoF / DoF_max);
        sz = std::max(0.6, std::min(2.0, sz));
        TGraph* gp = new TGraph(1);
        gp->SetPoint(0, c.eta / eta_max, Q_max / c.Q);
        gp->SetMarkerStyle(kOpenCircle);
        gp->SetMarkerColor(kGray + 1);
        gp->SetLineColor(kGray + 1);
        gp->SetMarkerSize(sz);
        gp->Draw("P SAME");
        g_bg_list.push_back(gp);
    }

    // Punti sul fronte (viridis per |M-1|) — tutti tranne il best
    gStyle->SetPalette(kViridis);
    const int ncolors = gStyle->GetNumberOfColors();

    double M_diff_max = 0.0;
    for (const auto* p : front)
        M_diff_max = std::max(M_diff_max, p->M_abs_err);
    if (M_diff_max <= 0.0) M_diff_max = 1.0;

    std::vector<TGraph*> g_front_list;
    for (const auto* p : front) {
        if (p->pareto_rank == 1) continue;  // best disegnato dopo
        if (p->Q <= 0.0) continue;
        double frac      = p->M_abs_err / M_diff_max;
        int    cidx      = static_cast<int>(frac * (ncolors - 1));
        Int_t  root_col  = gStyle->GetColorPalette(cidx);

        TGraph* gp = new TGraph(1);
        gp->SetPoint(0, p->eta / eta_max, Q_max / p->Q);
        gp->SetMarkerStyle(kFullCircle);
        gp->SetMarkerColor(root_col);
        gp->SetLineColor(kRed);
        gp->SetLineWidth(2);
        gp->SetMarkerSize(1.5);
        gp->Draw("P SAME");
        g_front_list.push_back(gp);
    }

    // Punto best — stella rossa
    TGraph* g_best = nullptr;
    if (!front.empty() && front[0]->Q > 0.0) {
        const ConfigData& best = *front[0];
        g_best = new TGraph(1);
        g_best->SetPoint(0, best.eta / eta_max, Q_max / best.Q);
        g_best->SetMarkerStyle(kFullStar);
        g_best->SetMarkerColor(kRed);
        g_best->SetMarkerSize(3.0);
        g_best->Draw("P SAME");

        // Label vicino al best
        double bx_ndc = pad_top->GetLeftMargin() +
                        (best.eta / eta_max - xlo) / (xhi - xlo) *
                        (1.0 - pad_top->GetLeftMargin() - pad_top->GetRightMargin());
        double by_ndc = pad_top->GetBottomMargin() +
                        (Q_max / best.Q - ylo) / (yhi - ylo) *
                        (1.0 - pad_top->GetBottomMargin() - pad_top->GetTopMargin());
        double pw = 0.20, ph = 0.10;
        double px1 = bx_ndc + 0.02;
        double py1 = by_ndc + 0.02;
        if (px1 + pw > 0.95) px1 = bx_ndc - pw - 0.02;
        if (py1 + ph > 0.95) py1 = by_ndc - ph - 0.02;
        py1 = std::max(0.0, py1);
        px1 = std::max(0.0, px1);

        TPaveText* pt = new TPaveText(px1, py1, px1 + pw, py1 + ph, "NDC");
        pt->SetFillColor(0);
        pt->SetBorderSize(1);
        pt->SetTextFont(42);
        pt->SetTextSize(0.030);
        pt->AddText("BEST");
        pt->AddText(Form("x1=%.1f  x2=%.1f", best.x1, best.x2));
        pt->AddText(Form("Mtot=%.3f", best.Mtot));
        pt->Draw();
    }

    // Legenda
    TLegend* leg = new TLegend(0.60, 0.68, 0.93, 0.91);
    leg->SetTextFont(42);
    leg->SetTextSize(0.030);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    if (!g_bg_list.empty())
        leg->AddEntry(g_bg_list[0], "Filtrati (non fronte)", "p");
    if (!g_front_list.empty())
        leg->AddEntry(g_front_list[0], "Fronte di Pareto", "p");
    if (g_best)
        leg->AddEntry(g_best, "Raccomandato", "p");
    leg->Draw();

    // ── Color bar |M−1| ───────────────────────────────────────────────────────
    draw_colorbar(pad_cb, 0.0, M_diff_max, "|M#minus{}M_{tgt}|");
    pad_top->cd();

    // ── Pad inferiore: tabella top-5 ──────────────────────────────────────────
    pad_bot->cd();

    TPaveText* table = new TPaveText(0.01, 0.01, 0.99, 0.99, "NDC");
    table->SetFillColor(0);
    table->SetBorderSize(0);
    table->SetTextFont(42);
    table->SetTextSize(0.065);
    table->SetTextAlign(12);
    table->AddText("Top-5 configurazioni sul fronte di Pareto (per M_{tot}):");

    int shown = 0;
    for (const auto* p : front) {
        if (shown >= 5) break;
        table->AddText(Form(
            "#%d  x1=%5.1f  x2=%6.1f  |  #eta=%.3f  Q=%.3f  #chi^{2}=%.2f"
            "  DoF=%.1f  M=%.3f  #DeltaM=%.4f  |  Mtot=%.3f",
            p->pareto_rank, p->x1, p->x2,
            p->eta, p->Q, p->chi2, p->DoF, p->M, p->M_abs_err, p->Mtot));
        ++shown;
    }
    table->Draw();

    // ── Salva ────────────────────────────────────────────────────────────────
    std::filesystem::create_directories(
        std::filesystem::path(cli.output_path).parent_path());
    canvas->SaveAs(cli.output_path.c_str());
    delete canvas;

    // ── 8. Stdout finale ─────────────────────────────────────────────────────
    const std::string sep(48, '=');
    std::cout << "\n" << sep << "\n"
              << " RIPTIDE - Pareto Selector - Risultati finali\n"
              << sep << "\n"
              << " Configurazioni analizzate:      " << configs.size() << "\n"
              << " Sul fronte di Pareto:           " << front.size()   << "\n";

    if (!front.empty()) {
        const ConfigData& b = *front[0];
        std::cout << " Configurazione raccomandata (#1):\n"
                  << "   x1 = " << b.x1 << " mm   x2 = " << b.x2 << " mm\n"
                  << "   eta= " << std::fixed << std::setprecision(3) << b.eta
                  << "     Q  = " << b.Q
                  << "    chi2= " << b.chi2 << "\n"
                  << "   DoF= " << b.DoF << " mm   M  = " << b.M
                  << "    Mtot= " << b.Mtot << "\n";
        for (int r = 1; r < (int)front.size() && r < 3; ++r) {
            const ConfigData& fp = *front[r];
            std::cout << " #" << (r + 1) << ": x1=" << fp.x1 << " x2=" << fp.x2
                      << " | Mtot=" << std::fixed << std::setprecision(3) << fp.Mtot << "\n";
        }
    }
    std::cout << sep << "\n"
              << " Output: " << cli.output_path << "\n"
              << "         " << cli.tsv_path    << "\n"
              << sep << "\n";

    return 0;
}
