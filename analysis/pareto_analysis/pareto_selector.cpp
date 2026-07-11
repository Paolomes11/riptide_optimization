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
#include <TColor.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMarker.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace riptide::pareto;

// ── CLI ──────────────────────────────────────────────────────────────────────

struct CliConfig {
    std::string events_path      = "output/optimization/events.root";
    std::string qmap_path        = "output/psf_analysis/q_map.tsv";
    std::string chi2map_path     = "output/psf_analysis/chi2_map.tsv";
    std::string dofmap_path      = "output/dof_analysis/dof_map.tsv";
    std::string resolution_path  = "";
    std::string output_path      = "output/pareto_analysis/pareto_plot.png";
    std::string tsv_path         = "output/pareto_analysis/pareto_results.tsv";
    std::string l1_id        = "";
    std::string l2_id        = "";
    FilterConfig fc;
    WeightConfig wc;

    // Modalità weight-sweep (ternario + mappa (x1,x2))
    bool        weight_sweep     = false;
    double      weight_step      = 0.05;
    std::string weight_sweep_tsv = "output/pareto_analysis/weight_sweep_results.tsv";
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
        "  --mobile-focus    Disattiva filtro x_det/focus_tol (detector statico);\n"
        "                    vincolo focale delegato a --dof-min           [off]\n"
        "  --resolution PATH TSV risoluzione (opt) [\"\"]\n"
        "  --ee80-max  F     EE80 max mm (0=off)   [0.0]\n"
        "  --w-eta     F     Peso eta in Mtot      [0.35]\n"
        "  --w-Q       F     Peso Q in Mtot        [0.40]\n"
        "  --w-dof     F     Peso DoF in Mtot      [0.15]\n"
        "  --w-M       F     Peso M in Mtot        [0.10]\n"
        "  --l1-id STR   Filtra lente L1       [\"\"]\n"
        "  --l2-id STR   Filtra lente L2       [\"\"]\n"
        "\n"
        "  --weight-sweep         Modalità sweep pesi (ternario + mappa x1,x2) invece\n"
        "                         del plot singolo-peso [off]\n"
        "  --weight-step F        Passo griglia baricentrica pesi        [0.05]\n"
        "  --weight-sweep-tsv PATH TSV aggregato sweep\n"
        "                         [output/pareto_analysis/weight_sweep_results.tsv]\n"
        "\n"
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
        else if (key == "--dof-min")    cfg.fc.dof_min      = std::stod(next());
        else if (key == "--mobile-focus") cfg.fc.mobile_focus = true;
        else if (key == "--resolution") cfg.resolution_path = next();
        else if (key == "--ee80-max")   cfg.fc.ee80_max     = std::stod(next());
        else if (key == "--w-eta")     cfg.wc.w_eta     = std::stod(next());
        else if (key == "--w-Q")       cfg.wc.w_Q       = std::stod(next());
        else if (key == "--w-dof")     cfg.wc.w_dof     = std::stod(next());
        else if (key == "--w-M")       cfg.wc.w_M       = std::stod(next());
        else if (key == "--l1-id") cfg.l1_id    = next();
        else if (key == "--l2-id") cfg.l2_id    = next();
        else if (key == "--weight-sweep")     cfg.weight_sweep     = true;
        else if (key == "--weight-step")      cfg.weight_step      = std::stod(next());
        else if (key == "--weight-sweep-tsv") cfg.weight_sweep_tsv = next();
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

// ── TSV header validation ─────────────────────────────────────────────────────

static void validate_tsv_headers(
    const std::vector<std::map<std::string, std::string>>& rows,
    const std::vector<std::string>& required,
    const std::string& path)
{
    if (rows.empty()) return;
    for (const auto& col : required) {
        if (!rows[0].count(col)) {
            std::cerr << "[ERRORE] Colonna '" << col << "' mancante in " << path << "\n";
            std::exit(1);
        }
    }
}

// ── Indice per join O(1) su (x1,x2) ──────────────────────────────────────────

using CoordIndex = std::unordered_map<int64_t, const std::map<std::string, std::string>*>;

static CoordIndex build_coord_index(const std::vector<std::map<std::string, std::string>>& rows) {
    CoordIndex idx;
    idx.reserve(rows.size());
    for (const auto& r : rows)
        idx[coord_key(std::stod(r.at("x1")), std::stod(r.at("x2")))] = &r;
    return idx;
}

static const std::map<std::string, std::string>* find_row(const CoordIndex& idx, double x1, double x2) {
    auto it = idx.find(coord_key(x1, x2));
    return it != idx.end() ? it->second : nullptr;
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
    gStyle->SetTitleOffset(1.10, "Y");
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
    ax->Draw();

    TLatex lbl;
    lbl.SetNDC();
    lbl.SetTextFont(42);
    lbl.SetTextSize(0.15);
    lbl.SetTextColor(kWhite);
    lbl.SetTextAlign(22);
    lbl.SetTextAngle(90.0);
    double cx = (cb_x0 + cb_x1) / 2.0;
    double cy = (cb_y0 + cb_y1) / 2.0;
    lbl.DrawLatex(cx, cy, title.c_str());
}

// ── Weight sweep: plotting ───────────────────────────────────────────────────

static Int_t categorical_color(int category_id, int n_cat) {
    const Float_t h = (n_cat > 0) ? 360.0f * (category_id % std::max(1, n_cat)) / std::max(1, n_cat) : 0.0f;
    Float_t r = 0.0f, g = 0.0f, b = 0.0f;
    TColor::HSV2RGB(h, 0.65f, 0.95f, r, g, b);
    return TColor::GetColor(r, g, b);
}

// Plot ternario: sviluppo (net) del tetraedro dei pesi (w_eta, w_Q, w_dof, w_M).
// La faccia centrale (w_M=0) e' il triangolo A/B/C di ternary_to_xy; le altre 3
// facce (w_eta=0, w_Q=0, w_dof=0) sono ribaltate verso l'esterno sui lati BC/AC/AB
// condivisi con la centrale — sviluppo standard del tetraedro per diagrammi di
// composizione quaternari (v. README/Bibliografia). Ciascuna faccia e' tassellata
// esplicitamente sul reticolo baricentrico (TGraph+Draw("F"), coppie "up"/"down"),
// senza buchi/sovrapposizioni per costruzione (dimostrato: i 3 apici dei lembi sono
// a due a due a distanza esatta 2, v. T14).
static void draw_ternary_plot(const std::vector<SweepPoint>& sweep,
                               int n_cat, double weight_step,
                               const std::string& out_path) {
    const int n = static_cast<int>(std::lround(1.0 / weight_step));

    // Vertici della faccia centrale (w_M=0): A=w_eta, B=w_Q, C=w_dof.
    const auto A = ternary_to_xy(1.0, 0.0, 0.0);
    const auto B = ternary_to_xy(0.0, 1.0, 0.0);
    const auto C = ternary_to_xy(0.0, 0.0, 1.0);

    // Apici dei 3 lembi ribaltati: ciascuno e' il vertice w_M=1 sulla faccia
    // opposta al vertice riflesso (dimostrato pairwise a distanza 2, v. T14).
    const auto A_flap = reflect_across_line(A, B, C);  // faccia w_eta=0
    const auto B_flap = reflect_across_line(B, A, C);  // faccia w_Q=0
    const auto C_flap = reflect_across_line(C, A, B);  // faccia w_dof=0

    const FaceAffine identity;                                    // faccia centrale
    const FaceAffine face_eta = solve_face_affine(B, C, A_flap);   // locale (w_Q, w_dof, w_M)
    const FaceAffine face_Q   = solve_face_affine(A, C, B_flap);   // locale (w_eta, w_dof, w_M)
    const FaceAffine face_dof = solve_face_affine(A, B, C_flap);   // locale (w_eta, w_Q, w_M)

    const double title_h = 0.05;
    const double info_h  = 0.14;

    TCanvas* canvas = new TCanvas("ternary_canvas", "Pareto Weight Sweep - Ternario", 1300, 1300);
    canvas->SetLeftMargin(0.0);
    canvas->SetRightMargin(0.0);
    canvas->SetTopMargin(0.0);
    canvas->SetBottomMargin(0.0);

    TPad* pad_net  = new TPad("pad_net",  "", 0.0, info_h, 1.0, 1.0 - title_h);
    TPad* pad_info = new TPad("pad_info", "", 0.0, 0.0,    1.0, info_h);
    canvas->cd();
    pad_net->Draw();
    pad_info->Draw();

    pad_net->cd();
    pad_net->SetLeftMargin(0.02);
    pad_net->SetRightMargin(0.02);
    pad_net->SetTopMargin(0.02);
    pad_net->SetBottomMargin(0.02);

    // Il net occupa un rombo [-0.5,1.5]x[-h,h] con h=sqrt(3)/2; margine 0.15.
    const double h_tri = std::sqrt(3.0) / 2.0;
    TH2D* frame = new TH2D("frame_net", "", 2, -0.65, 1.65, 2, -h_tri - 0.15, h_tri + 0.15);
    frame->GetXaxis()->SetLabelSize(0);
    frame->GetYaxis()->SetLabelSize(0);
    frame->GetXaxis()->SetTickLength(0);
    frame->GetYaxis()->SetTickLength(0);
    frame->Draw("AXIS");

    // Tassella una faccia: (get_i,get_j) estraggono gli indici di reticolo delle
    // prime 2 coordinate baricentriche locali (la terza e' n-i-j per costruzione);
    // on_face seleziona i punti dello sweep con il peso escluso da questa faccia
    // ~0; affine mappa il triangolo canonico locale sulla posizione globale.
    auto tessellate_face = [&](const std::function<int(const SweepPoint&)>& get_i,
                                const std::function<int(const SweepPoint&)>& get_j,
                                const std::function<bool(const SweepPoint&)>& on_face,
                                const FaceAffine& affine) {
        std::map<std::pair<int, int>, int> cat_at;
        for (const auto& sp : sweep) {
            if (!on_face(sp)) continue;
            cat_at[{get_i(sp), get_j(sp)}] = sp.category_id;
        }
        auto lattice_xy = [&](int i, int j) {
            const double a = (n > 0) ? static_cast<double>(i) / n : 0.0;
            const double b = (n > 0) ? static_cast<double>(j) / n : 0.0;
            const double c = 1.0 - a - b;
            const auto local = ternary_to_xy(a, b, c);
            return affine.apply(local.first, local.second);
        };
        auto draw_tri = [&](std::pair<double, double> p0, std::pair<double, double> p1,
                             std::pair<double, double> p2, Int_t col) {
            double xs[4] = {p0.first, p1.first, p2.first, p0.first};
            double ys[4] = {p0.second, p1.second, p2.second, p0.second};
            TGraph* tri = new TGraph(4, xs, ys);
            tri->SetFillColor(col);
            tri->SetLineColor(col);
            tri->Draw("F SAME");
        };
        for (int i = 0; n > 0 && i <= n; ++i) {
            for (int j = 0; i + j <= n; ++j) {
                if (i + j <= n - 1 && cat_at.count({i, j}))
                    draw_tri(lattice_xy(i, j), lattice_xy(i + 1, j), lattice_xy(i, j + 1),
                              categorical_color(cat_at.at({i, j}), n_cat));
                if (i + j <= n - 2 && cat_at.count({i + 1, j + 1}))
                    draw_tri(lattice_xy(i + 1, j), lattice_xy(i, j + 1), lattice_xy(i + 1, j + 1),
                              categorical_color(cat_at.at({i + 1, j + 1}), n_cat));
            }
        }
    };

    auto ridx = [n](double w) { return static_cast<int>(std::lround(w * n)); };

    tessellate_face(
        [ridx](const SweepPoint& sp) { return ridx(sp.wc.w_eta); },
        [ridx](const SweepPoint& sp) { return ridx(sp.wc.w_Q); },
        [ridx](const SweepPoint& sp) { return ridx(sp.wc.w_M) == 0; },
        identity);
    tessellate_face(
        [ridx](const SweepPoint& sp) { return ridx(sp.wc.w_Q); },
        [ridx](const SweepPoint& sp) { return ridx(sp.wc.w_dof); },
        [ridx](const SweepPoint& sp) { return ridx(sp.wc.w_eta) == 0; },
        face_eta);
    tessellate_face(
        [ridx](const SweepPoint& sp) { return ridx(sp.wc.w_eta); },
        [ridx](const SweepPoint& sp) { return ridx(sp.wc.w_dof); },
        [ridx](const SweepPoint& sp) { return ridx(sp.wc.w_Q) == 0; },
        face_Q);
    tessellate_face(
        [ridx](const SweepPoint& sp) { return ridx(sp.wc.w_eta); },
        [ridx](const SweepPoint& sp) { return ridx(sp.wc.w_Q); },
        [ridx](const SweepPoint& sp) { return ridx(sp.wc.w_dof) == 0; },
        face_dof);

    // Contorni delle 4 facce, per rendere visibile lo sviluppo geometrico.
    auto draw_outline = [](std::pair<double, double> p0, std::pair<double, double> p1,
                            std::pair<double, double> p2) {
        double xs[4] = {p0.first, p1.first, p2.first, p0.first};
        double ys[4] = {p0.second, p1.second, p2.second, p0.second};
        TGraph* g = new TGraph(4, xs, ys);
        g->SetLineColor(kBlack);
        g->SetLineWidth(2);
        g->Draw("L SAME");
    };
    draw_outline(A, B, C);
    draw_outline(B, C, A_flap);
    draw_outline(A, C, B_flap);
    draw_outline(A, B, C_flap);

    // Etichette vertici — TLatex puro (nessun box di sfondo), spinte radialmente
    // verso l'esterno dal centroide del net per non tagliare mai i contorni.
    const double net_cx = 0.5, net_cy = std::sqrt(3.0) / 6.0;
    auto label_vertex = [&](std::pair<double, double> p, const char* txt) {
        double dx = p.first - net_cx, dy = p.second - net_cy;
        const double len = std::hypot(dx, dy);
        if (len < 1e-9) { dx = 0.0; dy = 1.0; }
        else { dx /= len; dy /= len; }
        TLatex* lbl = new TLatex(p.first + 0.12 * dx, p.second + 0.12 * dy, txt);
        lbl->SetTextFont(62);
        lbl->SetTextAlign(22);
        lbl->SetTextSize(0.032);
        lbl->Draw();
    };
    label_vertex(A, "w_{#eta}");
    label_vertex(B, "w_{Q}");
    label_vertex(C, "w_{DoF}");
    label_vertex(A_flap, "w_{M}");
    label_vertex(B_flap, "w_{M}");
    label_vertex(C_flap, "w_{M}");

    // Titolo canvas
    canvas->cd();
    TLatex title;
    title.SetNDC();
    title.SetTextFont(62);
    title.SetTextAlign(22);
    title.SetTextSize(0.026);
    title.DrawLatex(0.5, 1.0 - title_h / 2.0,
        "Sweep pesi Pareto - sviluppo (net) del tetraedro dei pesi (w_{#eta}, w_{Q}, w_{DoF}, w_{M})");

    // Legenda: spiegazione fissa del significato dei colori, non elenco per categoria
    // (con n_cat~15-20 l'elenco sarebbe illeggibile; i due plot ternario+mappa si
    // spiegano a vicenda in tesi).
    pad_info->cd();
    TLatex info;
    info.SetNDC();
    info.SetTextFont(42);
    info.SetTextAlign(12);
    info.SetTextSize(0.22);
    info.DrawLatex(0.02, 0.65,
        "Colore = configurazione (x_{1},x_{2}) vincente");
    info.DrawLatex(0.02, 0.30,
        "per quella combinazione di pesi (v. pareto_ternary_map).");

    std::filesystem::create_directories(std::filesystem::path(out_path).parent_path());
    canvas->SaveAs(out_path.c_str());
    delete canvas;
}

// Mappa fisica (x1,x2): punti del fronte colorati per categoria vincente (stessa
// palette del ternario), non vincenti in grigio chiaro — dominio nativo di
// q_map/chi2_map/dof_map. Sovrappone la maschera del dominio invalido (stesso schema
// TBox grigio di q_map.cpp) per renderla direttamente confrontabile con quei plot, e
// una legenda esplicita (colore->categoria, marker->stato).
static void draw_secondary_map(const std::vector<ConfigData>& front,
                                const std::vector<SweepPoint>& sweep,
                                const std::vector<std::map<std::string, std::string>>& qmap_rows,
                                int n_cat,
                                const std::string& out_path) {
    std::unordered_map<int64_t, int> winner_cat;
    for (const auto& sp : sweep)
        winner_cat[coord_key(sp.x1_winner, sp.x2_winner)] = sp.category_id;

    apply_style();

    const double legend_h = 0.16;

    TCanvas* canvas = new TCanvas("sweep_map_canvas", "Pareto Weight Sweep - Mappa (x1,x2)", 1100, 1050);
    TPad* pad_plot = new TPad("pad_plot", "", 0.0, legend_h, 1.0, 1.0);
    TPad* pad_leg  = new TPad("pad_leg",  "", 0.0, 0.0,       1.0, legend_h);
    canvas->cd();
    pad_plot->Draw();
    pad_leg->Draw();

    pad_plot->cd();
    pad_plot->SetLeftMargin(0.16);
    pad_plot->SetRightMargin(0.04);
    pad_plot->SetTopMargin(0.10);
    pad_plot->SetBottomMargin(0.14);

    // Dominio (x1,x2) e passo griglia dedotti da TUTTE le righe di q_map (non solo
    // quelle sopravvissute al join), altrimenti il dominio invalido sarebbe
    // indistinguibile da "nessun dato".
    std::vector<double> x1_vals, x2_vals;
    double x1min = 1e9, x1max = -1e9, x2min = 1e9, x2max = -1e9;
    for (const auto& row : qmap_rows) {
        try {
            const double x1 = std::stod(row.at("x1"));
            const double x2 = std::stod(row.at("x2"));
            x1_vals.push_back(x1);
            x2_vals.push_back(x2);
            x1min = std::min(x1min, x1); x1max = std::max(x1max, x1);
            x2min = std::min(x2min, x2); x2max = std::max(x2max, x2);
        } catch (...) {}
    }
    if (x1min >= x1max) { x1min -= 1.0; x1max += 1.0; }
    if (x2min >= x2max) { x2min -= 1.0; x2max += 1.0; }

    auto infer_step = [](std::vector<double> v) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end(),
                             [](double a, double b) { return std::abs(a - b) < 1e-9; }),
                 v.end());
        double step = 1.0;
        for (size_t i = 1; i < v.size(); ++i) step = std::min(step, v[i] - v[i - 1]);
        return (v.size() >= 2) ? step : 1.0;
    };
    const double dx = infer_step(x1_vals);
    const double dy = infer_step(x2_vals);

    const int bins_x = std::max(1, static_cast<int>(std::round((x1max - x1min) / dx)) + 1);
    const int bins_y = std::max(1, static_cast<int>(std::round((x2max - x2min) / dy)) + 1);
    const double ax_x1_lo = x1min - dx / 2.0, ax_x1_hi = x1max + dx / 2.0;
    const double ax_x2_lo = x2min - dy / 2.0, ax_x2_hi = x2max + dy / 2.0;

    TH2D* frame = new TH2D("frame_map", "", bins_x, ax_x1_lo, ax_x1_hi, bins_y, ax_x2_lo, ax_x2_hi);
    frame->GetXaxis()->SetTitle("x1  [mm]");
    frame->GetYaxis()->SetTitle("x2  [mm]");
    frame->GetXaxis()->SetTitleOffset(1.2);
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->Draw("AXIS");

    // Riempie l'INTERO dominio simulato (non solo le celle invalide) cosi' il
    // "triangolo" delle configurazioni valide e' visibile come area piena,
    // comparabile a q_map.png/dof_map.png: valido = grigio chiaro, invalido = nero.
    const Int_t valid_color   = TColor::GetColor(224, 224, 224);
    const Int_t invalid_color = kBlack;
    TH2D h_has("h_has_map", "", bins_x, ax_x1_lo, ax_x1_hi, bins_y, ax_x2_lo, ax_x2_hi);
    TH2D h_inv("h_inv_map", "", bins_x, ax_x1_lo, ax_x1_hi, bins_y, ax_x2_lo, ax_x2_hi);
    for (const auto& row : qmap_rows) {
        try {
            const double x1 = std::stod(row.at("x1"));
            const double x2 = std::stod(row.at("x2"));
            const int bx = h_has.GetXaxis()->FindFixBin(x1);
            const int by = h_has.GetYaxis()->FindFixBin(x2);
            h_has.SetBinContent(bx, by, 1.0);
            if (row.at("config_valid") == "0")
                h_inv.SetBinContent(bx, by, 1.0);
        } catch (...) {}
    }
    for (int iy = 1; iy <= h_has.GetNbinsY(); ++iy) {
        for (int ix = 1; ix <= h_has.GetNbinsX(); ++ix) {
            if (h_has.GetBinContent(ix, iy) <= 0.5) continue;  // nessuna riga q_map per questa cella
            const bool invalid = h_inv.GetBinContent(ix, iy) > 0.5;
            TBox* b = new TBox(h_has.GetXaxis()->GetBinLowEdge(ix), h_has.GetYaxis()->GetBinLowEdge(iy),
                                h_has.GetXaxis()->GetBinUpEdge(ix),  h_has.GetYaxis()->GetBinUpEdge(iy));
            b->SetFillColor(invalid ? invalid_color : valid_color);
            b->SetLineColor(invalid ? invalid_color : valid_color);
            b->Draw("same");
        }
    }

    for (const auto& c : front) {
        auto it = winner_cat.find(coord_key(c.x1, c.x2));
        TGraph* gp = new TGraph(1);
        gp->SetPoint(0, c.x1, c.x2);
        if (it != winner_cat.end()) {
            gp->SetMarkerStyle(kFullCircle);
            gp->SetMarkerColor(categorical_color(it->second, n_cat));
            gp->SetMarkerSize(1.6);
        } else {
            gp->SetMarkerStyle(kOpenCircle);
            gp->SetMarkerColor(kBlack);
            gp->SetMarkerSize(0.8);
        }
        gp->Draw("P SAME");
    }

    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextSize(0.035);
    title.SetTextAlign(22);
    title.DrawLatex(0.5, 0.96, "Configurazioni vincenti nello sweep dei pesi");

    // Legenda: spiegazione fissa del significato di colori/marker/sfondo, non
    // elenco per categoria (con n_cat~15-20 l'elenco sarebbe illeggibile; questo
    // plot e pareto_ternary.png si spiegano a vicenda in tesi tramite lo stesso
    // mapping categoria->colore). Stile TLegend coerente con quella di
    // draw_pareto_frontier (SetBorderSize(0)+SetFillStyle(0)); swatch grigio
    // neutro per il marker pieno (non un colore specifico: il colore reale
    // varia per categoria, v. categorical_color).
    pad_leg->cd();

    TGraph* proxy_winner = new TGraph(1);
    proxy_winner->SetPoint(0, 0, 0);
    proxy_winner->SetMarkerStyle(kFullCircle);
    proxy_winner->SetMarkerColor(kGray + 2);
    proxy_winner->SetMarkerSize(1.6);

    TGraph* proxy_front = new TGraph(1);
    proxy_front->SetPoint(0, 0, 0);
    proxy_front->SetMarkerStyle(kOpenCircle);
    proxy_front->SetMarkerColor(kBlack);
    proxy_front->SetMarkerSize(1.4);

    TBox* proxy_valid = new TBox(0, 0, 0, 0);
    proxy_valid->SetFillColor(valid_color);
    proxy_valid->SetLineColor(kBlack);

    TBox* proxy_invalid = new TBox(0, 0, 0, 0);
    proxy_invalid->SetFillColor(invalid_color);
    proxy_invalid->SetLineColor(kBlack);

    TLegend* leg = new TLegend(0.02, 0.03, 0.99, 0.97);
    leg->SetTextFont(42);
    leg->SetTextAlign(12);
    leg->SetTextSize(0.11);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetMargin(0.05);
    leg->AddEntry(proxy_winner, "Configurazione (x_{1},x_{2}) vincente in almeno una combinazione di pesi (colore = categoria, v. pareto_ternary)", "p");
    leg->AddEntry(proxy_front, "Sul fronte di Pareto, mai vincente nello sweep dei pesi", "p");
    leg->AddEntry(proxy_valid, "Dominio simulato valido", "f");
    leg->AddEntry(proxy_invalid, "Dominio non valido (vincoli geometrici della lente)", "f");
    leg->Draw();

    std::filesystem::create_directories(std::filesystem::path(out_path).parent_path());
    canvas->SaveAs(out_path.c_str());
    delete canvas;
}

// Orchestrazione della modalità --weight-sweep: fronte+normalizzatori calcolati
// una sola volta (indipendenti dai pesi), poi sweep economico sulla griglia.
static void run_weight_sweep_mode(std::vector<ConfigData>& configs,
                                   const std::vector<std::map<std::string, std::string>>& qmap_rows,
                                   const CliConfig& cli) {
    mark_pareto_front(configs);
    const MtotNormalizers normalizers = compute_mtot_normalizers(configs);

    std::vector<ConfigData> front_copy;
    for (const auto& c : configs)
        if (c.on_pareto) front_copy.push_back(c);

    std::cout << "[SWEEP] Configurazioni sul fronte: " << front_copy.size() << "\n";

    std::filesystem::create_directories(std::filesystem::path(cli.weight_sweep_tsv).parent_path());
    if (front_copy.empty()) {
        std::cerr << "[ERRORE] Fronte di Pareto vuoto: impossibile eseguire lo sweep.\n";
        std::ofstream tsv(cli.weight_sweep_tsv);
        if (tsv.is_open())
            tsv << "w_eta\tw_Q\tw_dof\tw_M\tx1_winner\tx2_winner\tMtot_winner\tcategory_id\n";
        return;
    }

    const auto grid  = generate_barycentric_grid(cli.weight_step);
    auto sweep        = run_weight_sweep(front_copy, normalizers, grid);
    assign_category_ids(sweep);

    int n_cat = 0;
    for (const auto& sp : sweep)
        n_cat = std::max(n_cat, sp.category_id + 1);

    std::cout << "[SWEEP] Punti griglia pesi (step=" << cli.weight_step << "): " << sweep.size()
              << "\n[SWEEP] Configurazioni vincenti distinte: " << n_cat << "\n\n";

    {
        std::ofstream tsv(cli.weight_sweep_tsv);
        if (!tsv.is_open()) {
            std::cerr << "[ERRORE] Impossibile aprire TSV sweep: " << cli.weight_sweep_tsv << "\n";
        } else {
            tsv << "w_eta\tw_Q\tw_dof\tw_M\tx1_winner\tx2_winner\tMtot_winner\tcategory_id\n";
            for (const auto& sp : sweep) {
                tsv << std::fixed << std::setprecision(6)
                    << sp.wc.w_eta << "\t" << sp.wc.w_Q << "\t" << sp.wc.w_dof << "\t" << sp.wc.w_M << "\t"
                    << sp.x1_winner << "\t" << sp.x2_winner << "\t" << sp.Mtot_winner << "\t"
                    << sp.category_id << "\n";
            }
        }
    }

    apply_style();
    draw_ternary_plot(sweep, n_cat, cli.weight_step, cli.output_path);

    const std::filesystem::path out_path(cli.output_path);
    const std::string secondary_path =
        (out_path.parent_path() / (out_path.stem().string() + "_map" + out_path.extension().string())).string();
    draw_secondary_map(front_copy, sweep, qmap_rows, n_cat, secondary_path);

    std::cout << "[SWEEP] Output: " << cli.weight_sweep_tsv << "\n"
              << "         " << cli.output_path << "\n"
              << "         " << secondary_path  << "\n";
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
    double n_hits_val  = 0.0;
    char l1_buf[4096]  = {};
    char l2_buf[4096]  = {};

    tree_config->SetBranchAddress("x1",        &x1_val);
    tree_config->SetBranchAddress("x2",        &x2_val);
    tree_config->SetBranchAddress("config_id", &cfg_id);
    const bool has_lens = (tree_config->GetBranch("l1_id") != nullptr);
    if (has_lens) {
        tree_config->SetBranchAddress("l1_id", l1_buf);
        tree_config->SetBranchAddress("l2_id", l2_buf);
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
        if (!cli.l1_id.empty() && s1 != cli.l1_id) continue;
        if (!cli.l2_id.empty() && s2 != cli.l2_id) continue;
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

    validate_tsv_headers(qmap_rows,    {"x1", "x2", "metric"},                          cli.qmap_path);
    validate_tsv_headers(chi2map_rows, {"x1", "x2", "metric"},                          cli.chi2map_path);
    validate_tsv_headers(dofmap_rows,  {"x1", "x2", "dof", "M", "M_abs_err", "x_focus"}, cli.dofmap_path);

    // ── 3. Join inner su (x1, x2) ────────────────────────────────────────────
    int n_found_q    = 0;
    int n_found_chi2 = 0;
    int n_found_dof  = 0;
    int n_joined     = 0;

    const CoordIndex qmap_idx    = build_coord_index(qmap_rows);
    const CoordIndex chi2map_idx = build_coord_index(chi2map_rows);
    const CoordIndex dofmap_idx  = build_coord_index(dofmap_rows);

    std::vector<ConfigData> configs;
    configs.reserve(cfg_map.size());

    for (const auto& [id, rc] : cfg_map) {
        double eta = 0.0;
        if (eta_cnt.count(id) && eta_cnt.at(id) > 0)
            eta = eta_sum.at(id) / eta_cnt.at(id);

        const std::map<std::string, std::string>* qrow    = find_row(qmap_idx,    rc.x1, rc.x2);
        const std::map<std::string, std::string>* chi2row = find_row(chi2map_idx, rc.x1, rc.x2);
        const std::map<std::string, std::string>* dofrow  = find_row(dofmap_idx,  rc.x1, rc.x2);

        // Contatori indipendenti per il riepilogo join
        if (qrow)    ++n_found_q;
        if (chi2row) ++n_found_chi2;
        if (dofrow)  ++n_found_dof;

        if (!qrow || !chi2row || !dofrow) continue;

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

        // x_lo/x_hi: colonne opzionali (TSV vecchi non le hanno) — fallback simmetrico ±DoF/2
        {
            bool has_lo_hi = false;
            double x_lo_raw = 0.0, x_hi_raw = 0.0;
            try {
                x_lo_raw = std::stod(dofrow->at("x_lo"));
                x_hi_raw = std::stod(dofrow->at("x_hi"));
                has_lo_hi = true;
            } catch (...) {}
            resolve_dof_bounds(cd, has_lo_hi, x_lo_raw, x_hi_raw);
        }

        // Salta valori non finiti (es. "NaN" nel TSV per config non valide)
        if (!std::isfinite(cd.Q) || !std::isfinite(cd.chi2) ||
            !std::isfinite(cd.DoF) || !std::isfinite(cd.M_abs_err))
            continue;

        configs.push_back(cd);
        ++n_joined;
    }

    std::cout << "[JOIN] Configurazioni (dopo filtro lente): " << n_events_root << "\n"
              << "[JOIN] Trovate in q_map.tsv:            " << n_found_q      << "\n"
              << "[JOIN] Trovate in chi2_map.tsv:         " << n_found_chi2   << "\n"
              << "[JOIN] Trovate in dof_map.tsv:          " << n_found_dof    << "\n"
              << "[JOIN] Configurazioni complete (join):  " << n_joined       << "\n\n";

    if (n_events_root - n_found_q > 0)
        std::cerr << "[WARN] " << (n_events_root - n_found_q)
                  << " configurazioni da events.root assenti in " << cli.qmap_path
                  << " (griglia PSF incompleta — ri-eseguire lens scan con lens_gap_margin aggiornato)\n\n";

    if (configs.empty()) {
        std::cerr << "[ERRORE] Join vuoto: nessuna configurazione completa trovata.\n";
        return 1;
    }

    // ── 3b. Join opzionale con TSV di resolution (EE80) ─────────────────────
    if (!cli.resolution_path.empty()) {
        auto res_rows = load_tsv(cli.resolution_path);
        validate_tsv_headers(res_rows, {"x1", "x2", "EE80_mean"}, cli.resolution_path);
        const CoordIndex res_idx = build_coord_index(res_rows);
        int n_ee80_found = 0;
        for (auto& cd : configs) {
            const auto* r = find_row(res_idx, cd.x1, cd.x2);
            if (r) {
                try { cd.EE80 = std::stod(r->at("EE80_mean")); ++n_ee80_found; } catch (...) {}
            }
        }
        std::cout << "[JOIN] Trovate in resolution.tsv (EE80): " << n_ee80_found
                  << "/" << configs.size() << "\n\n";
    }

    // ── 4. Pipeline filtri ───────────────────────────────────────────────────
    std::cout << "[FILTER] Totale prima dei filtri:       " << configs.size() << "\n";

    configs = apply_eta_filter(configs, cli.fc);
    std::cout << "[FILTER] Dopo filtro eta (>=" << cli.fc.eta_frac << "*eta_max):   "
              << configs.size() << "\n";

    configs = apply_focus_filter(configs, cli.fc);
    if (cli.fc.mobile_focus) {
        std::cout << "[FILTER] Filtro fuoco fisso disattivato (--mobile-focus): "
                     "vincolo focale delegato a DoF>=dof_min\n";
    } else {
        std::cout << "[FILTER] Dopo filtro fuoco (|x*-" << cli.fc.x_det << "|<="
                  << cli.fc.focus_tol << "): " << configs.size() << "\n";
    }

    if (cli.fc.dof_min > 0.0) {
        configs = apply_dof_filter(configs, cli.fc);
        std::cout << "[FILTER] Dopo filtro DoF (>=" << cli.fc.dof_min << " mm):     "
                  << configs.size() << "\n";
    }

    if (cli.fc.ee80_max > 0.0) {
        configs = apply_ee80_filter(configs, cli.fc);
        std::cout << "[FILTER] Dopo filtro EE80 (<=" << cli.fc.ee80_max << " mm):   "
                  << configs.size() << "\n";
    }

    std::cout << "\n";

    if (configs.empty()) {
        std::cerr << "[WARN] Nessuna configurazione sopravvissuta ai filtri.\n";
        if (cli.weight_sweep) {
            std::filesystem::create_directories(
                std::filesystem::path(cli.weight_sweep_tsv).parent_path());
            std::ofstream tsv(cli.weight_sweep_tsv);
            if (tsv.is_open())
                tsv << "w_eta\tw_Q\tw_dof\tw_M\tx1_winner\tx2_winner\tMtot_winner\tcategory_id\n";
            return 0;
        }
        std::filesystem::create_directories(
            std::filesystem::path(cli.tsv_path).parent_path());
        std::ofstream tsv(cli.tsv_path);
        if (tsv.is_open())
            tsv << "x1\tx2\teta\teta_norm\tQ\tchi2\tDoF\tM\tM_abs_err\t"
                   "x_focus\tEE80\ton_pareto\tMtot\tpareto_rank\n";
        return 0;
    }

    if (cli.weight_sweep) {
        run_weight_sweep_mode(configs, qmap_rows, cli);
        return 0;
    }

    // ── 5. Calcolo Mtot e fronte di Pareto ──────────────────────────────────
    {
        double w_sum = cli.wc.w_eta + cli.wc.w_Q + cli.wc.w_dof + cli.wc.w_M;
        if (std::abs(w_sum - 1.0) > 0.01)
            std::cerr << "[WARN] Pesi non sommano a 1.0 (somma=" << std::fixed
                      << std::setprecision(3) << w_sum
                      << "); normalizzati internamente.\n\n";
    }
    compute_pareto_front(configs, cli.wc);

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
    double Q_min   = std::numeric_limits<double>::max();
    for (const auto& c : configs) {
        eta_max = std::max(eta_max, c.eta);
        Q_max   = std::max(Q_max,   c.Q);
        DoF_max = std::max(DoF_max,  c.DoF);
        if (c.Q > 0.0) Q_min = std::min(Q_min, c.Q);
    }
    if (eta_max <= 0.0) eta_max = 1.0;
    if (Q_max   <= 0.0) Q_max   = 1.0;
    if (DoF_max <= 0.0) DoF_max = 1.0;
    if (Q_min >= 1e30 || Q_min <= 0.0) Q_min = 1.0;

    // ── 6. TSV output ─────────────────────────────────────────────────────────
    {
        std::filesystem::create_directories(
            std::filesystem::path(cli.tsv_path).parent_path());
        std::ofstream tsv(cli.tsv_path);
        if (!tsv.is_open()) {
            std::cerr << "[ERRORE] Impossibile aprire TSV output: " << cli.tsv_path << "\n";
        } else {
            tsv << "x1\tx2\teta\teta_norm\tQ\tchi2\tDoF\tM\tM_abs_err\t"
                   "x_focus\tEE80\ton_pareto\tMtot\tpareto_rank\n";
            for (const auto& c : configs) {
                tsv << std::fixed << std::setprecision(6)
                    << c.x1       << "\t" << c.x2         << "\t"
                    << c.eta      << "\t" << c.eta / eta_max << "\t"
                    << c.Q        << "\t" << c.chi2         << "\t"
                    << c.DoF      << "\t" << c.M            << "\t"
                    << c.M_abs_err << "\t" << c.x_focus     << "\t"
                    << c.EE80     << "\t"
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

    TPad* pad_top = new TPad("pad_top", "", 0.00, 0.00, 0.84, 1.0);
    TPad* pad_cb  = new TPad("pad_cb",  "", 0.84, 0.00, 1.00, 1.0);

    pad_top->SetLeftMargin(0.12);
    pad_top->SetRightMargin(0.02);
    pad_top->SetTopMargin(0.12);
    pad_top->SetBottomMargin(0.14);
    pad_top->SetGridx();
    pad_top->SetGridy();

    pad_cb->SetLeftMargin(0.10);
    pad_cb->SetRightMargin(0.38);
    pad_cb->SetTopMargin(0.08);
    pad_cb->SetBottomMargin(0.14);

    canvas->cd();
    pad_top->Draw();
    pad_cb->Draw();

    // ── Pad superiore: scatter plot ───────────────────────────────────────────
    pad_top->cd();
    {
      TLatex main_t;
      main_t.SetNDC();
      main_t.SetTextFont(42);
      main_t.SetTextSize(0.040);
      main_t.SetTextAlign(22);
      main_t.DrawLatex(0.50, 0.96, "Analisi del Fronte di Pareto");
    }

    // Calcola range assi
    double xmin_p = 1e9, xmax_p = -1e9, ymin_p = 1e9, ymax_p = -1e9;
    for (const auto& c : configs) {
        if (c.Q <= 0.0) continue;
        double xv = c.eta / eta_max;
        double yv = Q_min  / c.Q;
        xmin_p = std::min(xmin_p, xv);
        xmax_p = std::max(xmax_p, xv);
        ymin_p = std::min(ymin_p, yv);
        ymax_p = std::max(ymax_p, yv);
    }
    if (xmin_p >= xmax_p) { xmin_p = 0.0; xmax_p = 1.0; }
    if (ymin_p >= ymax_p) { ymin_p = 0.5; ymax_p = 1.5; }
    double xpad = 0.08 * (xmax_p - xmin_p);
    double ypad = 0.10 * (ymax_p - ymin_p);
    double xlo = std::max(0.0, xmin_p - xpad);
    double xhi = std::min(1.05, xmax_p + xpad);
    double ylo = std::max(0.0, ymin_p - ypad);
    double yhi = ymax_p + ypad;

    // Frame vuoto per assi
    TH2D* frame = new TH2D("frame", "", 100, xlo, xhi, 100, ylo, yhi);
    frame->GetXaxis()->SetTitle("#eta/#eta_{max}  [a.d.]");
    frame->GetYaxis()->SetTitle("#DeltaQ_{min}/#DeltaQ  [a.d.]");
    frame->GetXaxis()->SetTitleOffset(1.2);
    frame->GetYaxis()->SetTitleOffset(1.10);
    frame->Draw("AXIS");

    // Punti di background (non sul fronte) — marker size proporzionale a DoF
    std::vector<TGraph*> g_bg_list;
    for (const auto& c : configs) {
        if (c.on_pareto || c.Q <= 0.0) continue;
        double sz = 0.6 + 1.4 * (c.DoF / DoF_max);
        sz = std::max(0.6, std::min(2.0, sz));
        TGraph* gp = new TGraph(1);
        gp->SetPoint(0, c.eta / eta_max, Q_min / c.Q);
        gp->SetMarkerStyle(kOpenCircle);
        gp->SetMarkerColor(kGray + 1);
        gp->SetLineColor(kGray + 1);
        gp->SetMarkerSize(sz);
        gp->Draw("P SAME");
        g_bg_list.push_back(gp);
    }

    // Punti sul fronte (viridis per M assoluto) — tutti tranne il best
    gStyle->SetPalette(kViridis);
    const int ncolors = gStyle->GetNumberOfColors();

    double M_min = std::numeric_limits<double>::max();
    double M_max = -std::numeric_limits<double>::max();
    for (const auto* p : front) {
        M_min = std::min(M_min, p->M);
        M_max = std::max(M_max, p->M);
    }
    if (M_max <= M_min) { M_min -= 0.5; M_max += 0.5; }

    auto m_color = [&](double m) {
        double frac = (m - M_min) / (M_max - M_min);
        int    cidx = static_cast<int>(std::clamp(frac, 0.0, 1.0) * (ncolors - 1));
        return gStyle->GetColorPalette(cidx);
    };

    std::vector<TGraph*> g_front_list;
    for (const auto* p : front) {
        if (p->pareto_rank == 1) continue;  // best disegnato dopo
        if (p->Q <= 0.0) continue;
        Int_t root_col = m_color(p->M);

        TGraph* gp = new TGraph(1);
        gp->SetPoint(0, p->eta / eta_max, Q_min / p->Q);
        gp->SetMarkerStyle(kFullCircle);
        gp->SetMarkerColor(root_col);
        gp->SetLineColor(kRed);
        gp->SetLineWidth(2);
        gp->SetMarkerSize(1.5);
        gp->Draw("P SAME");
        g_front_list.push_back(gp);
    }

    // Punto best — stella con interno colorato per M e bordo rosso spesso e ben
    // visibile: stella rossa piena leggermente più grande disegnata sotto, e
    // stella colorata per M più piccola sopra (il bordo rosso e' l'anello
    // visibile tra le due dimensioni — piu' robusto di kOpenStar, il cui
    // contorno ROOT disegna con una linea troppo sottile per essere notata).
    TGraph* g_best = nullptr;
    if (!front.empty() && front[0]->Q > 0.0) {
        const ConfigData& best = *front[0];
        g_best = new TGraph(1);
        g_best->SetPoint(0, best.eta / eta_max, Q_min / best.Q);
        g_best->SetMarkerStyle(kFullStar);
        g_best->SetMarkerColor(kRed);
        g_best->SetMarkerSize(4.6);
        g_best->Draw("P SAME");

        TGraph* g_best_fill = new TGraph(1);
        g_best_fill->SetPoint(0, best.eta / eta_max, Q_min / best.Q);
        g_best_fill->SetMarkerStyle(kFullStar);
        g_best_fill->SetMarkerColor(m_color(best.M));
        g_best_fill->SetMarkerSize(3.2);
        g_best_fill->Draw("P SAME");
    }

    // Legenda — angolo in alto a destra per non intersecare i dati (concentrati
    // tipicamente in basso a sinistra/al centro nel piano eta_norm/Q_norm)
    TLegend* leg = new TLegend(0.62, 0.55, 0.97, 0.72);
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

    // ── Color bar M ──────────────────────────────────────────────────────────
    draw_colorbar(pad_cb, M_min, M_max, "M  [a.d.]");
    pad_top->cd();

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
