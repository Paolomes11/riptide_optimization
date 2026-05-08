/* Copyright 2026 Giulio Mesini — Apache License 2.0 */

#include "exp3_analysis.hpp"
#include "homography.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace riptide::exp3;

// ── Helpers ──────────────────────────────────────────────────────────────────

static int g_n_pass = 0;
static int g_n_fail = 0;

static bool check(bool cond, const std::string& msg, const std::string& T = "") {
    if (cond) {
        std::cout << "  [PASS] " << msg << "\n";
        ++g_n_pass;
    } else {
        std::cerr << "  [FAIL] " << msg;
        if (!T.empty()) std::cerr << "  (in " << T << ")";
        std::cerr << "\n";
        ++g_n_fail;
    }
    return cond;
}

static bool near(double a, double b, double tol, const std::string& label,
                 const std::string& T = "") {
    bool ok = std::abs(a - b) <= tol;
    std::ostringstream msg;
    msg << label << ": got " << a << ", expected " << b
        << " (tol=" << tol << ", diff=" << std::abs(a - b) << ")";
    return check(ok, msg.str(), T);
}

// ── Dati sintetici ────────────────────────────────────────────────────────────

// Griglia di punti display su un rettangolo rows×cols a partire da (x0,y0) con passo step.
static std::vector<std::pair<double, double>> grid_pts(double x0, double y0,
                                                        double step, int rows, int cols) {
    std::vector<std::pair<double, double>> pts;
    pts.reserve(static_cast<size_t>(rows * cols));
    for (int r = 0; r < rows; ++r)
        for (int c = 0; c < cols; ++c)
            pts.push_back({x0 + c * step, y0 + r * step});
    return pts;
}

// Punti CalibPoint da trasformazione affine: [xs,ys] = A·[xd,yd] + [tx,ty]
static std::vector<CalibPoint> make_affine_pts(
    double a00, double a01, double tx,
    double a10, double a11, double ty,
    const std::vector<std::pair<double, double>>& disp) {
    std::vector<CalibPoint> pts;
    pts.reserve(disp.size());
    for (const auto& [xd, yd] : disp) {
        CalibPoint p;
        p.disp_x = xd;
        p.disp_y = yd;
        p.sens_x = a00 * xd + a01 * yd + tx;
        p.sens_y = a10 * xd + a11 * yd + ty;
        pts.push_back(p);
    }
    return pts;
}

// Punti CalibPoint da omografia proiettiva H (row-major 3×3).
static std::vector<CalibPoint> make_proj_pts(
    const std::array<double, 9>& H,
    const std::vector<std::pair<double, double>>& disp) {
    std::vector<CalibPoint> pts;
    pts.reserve(disp.size());
    for (const auto& [xd, yd] : disp) {
        double w  = H[6] * xd + H[7] * yd + H[8];
        CalibPoint p;
        p.disp_x = xd;
        p.disp_y = yd;
        p.sens_x = (H[0] * xd + H[1] * yd + H[2]) / w;
        p.sens_y = (H[3] * xd + H[4] * yd + H[5]) / w;
        pts.push_back(p);
    }
    return pts;
}

// Immagine sintetica con gaussiane 2D (row-major).
static std::vector<double> make_gauss_image(
    int W, int H,
    const std::vector<std::pair<double, double>>& centers,
    double amplitude, double sigma) {
    std::vector<double> img(static_cast<size_t>(W * H), 0.0);
    double inv2s2 = 1.0 / (2.0 * sigma * sigma);
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            double val = 0.0;
            for (const auto& [cx, cy] : centers) {
                double dx = x - cx, dy = y - cy;
                val += amplitude * std::exp(-(dx * dx + dy * dy) * inv2s2);
            }
            img[static_cast<size_t>(y) * static_cast<size_t>(W) + static_cast<size_t>(x)] = val;
        }
    }
    return img;
}

// Scrive un TSV in /tmp/ e restituisce il path.
static std::filesystem::path write_tmp_tsv(const std::string& content,
                                            const std::string& fname) {
    std::filesystem::path p = fname;
    std::ofstream ofs(p);
    ofs << content;
    return p;
}

// Rileva NaN tramite bit IEEE 754 — robusto con -ffast-math che rompe std::isnan.
static bool is_nan_bits(double v) {
    uint64_t bits;
    __builtin_memcpy(&bits, &v, sizeof(bits));
    return (bits & 0x7FF0000000000000ULL) == 0x7FF0000000000000ULL
        && (bits & 0x000FFFFFFFFFFFFFULL) != 0ULL;
}

// PRNG lineare congruente per rumore riproducibile (no <random> per semplicità).
static double lcg_gauss(unsigned& s) {
    s = s * 1664525u + 1013904223u;
    double u1 = (s & 0xFFFFu) / 65535.0 + 1e-12;
    s = s * 1664525u + 1013904223u;
    double u2 = (s & 0xFFFFu) / 65535.0 + 1e-12;
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
}

// ── Nodi standard per i test detect_calibration_dots ─────────────────────────
// display 400×300, grid_step=100 → gx∈{100,200,300}, gy∈{100,200} → 6 nodi
static const int DISP_W    = 400;
static const int DISP_H    = 300;
static const int GRID_STEP = 100;
static const std::vector<std::pair<double, double>> STD_NODES = {
    {100.0, 100.0}, {200.0, 100.0}, {300.0, 100.0},
    {100.0, 200.0}, {200.0, 200.0}, {300.0, 200.0}};

// ── Cluster H: compute_homography / apply_homography ─────────────────────────

static void test_H1() {
    std::cout << "\n[H1] Identità: 4 punti disp==sens\n";
    const std::string T = "H1";
    auto pts = make_affine_pts(1, 0, 0, 0, 1, 0,
                               grid_pts(100, 100, 100, 2, 2));
    auto hom = compute_homography(pts);
    near(hom.H[0], 1.0, 1e-9, "H[0]", T);
    near(hom.H[4], 1.0, 1e-9, "H[4]", T);
    near(hom.H[8], 1.0, 1e-9, "H[8]", T);
    near(hom.H[1], 0.0, 1e-9, "H[1]", T);
    near(hom.H[2], 0.0, 1e-9, "H[2]", T);
    near(hom.H[3], 0.0, 1e-9, "H[3]", T);
    near(hom.H[5], 0.0, 1e-9, "H[5]", T);
    near(hom.rms_residual, 0.0, 1e-9, "rms", T);
    check(hom.n_points == 4, "n_points=4", T);
}

static void test_H2() {
    std::cout << "\n[H2] Traslazione pura tx=50, ty=-30\n";
    const std::string T = "H2";
    auto pts = make_affine_pts(1, 0, 50, 0, 1, -30,
                               grid_pts(100, 100, 100, 3, 3));
    auto hom = compute_homography(pts);
    near(hom.H[0], 1.0,   1e-8, "H[0]", T);
    near(hom.H[4], 1.0,   1e-8, "H[4]", T);
    near(hom.H[2], 50.0,  1e-8, "H[2]", T);
    near(hom.H[5], -30.0, 1e-8, "H[5]", T);
    near(hom.H[1], 0.0,   1e-8, "H[1]", T);
    near(hom.H[3], 0.0,   1e-8, "H[3]", T);
    near(hom.rms_residual, 0.0, 1e-8, "rms", T);
}

static void test_H3() {
    std::cout << "\n[H3] Scaling isotropico s=2.5\n";
    const std::string T = "H3";
    auto pts = make_affine_pts(2.5, 0, 0, 0, 2.5, 0,
                               grid_pts(40, 40, 40, 3, 3));
    auto hom = compute_homography(pts);
    near(hom.H[0], 2.5, 1e-8, "H[0]", T);
    near(hom.H[4], 2.5, 1e-8, "H[4]", T);
    near(hom.H[1], 0.0, 1e-8, "H[1]", T);
    near(hom.H[3], 0.0, 1e-8, "H[3]", T);
    near(hom.H[2], 0.0, 1e-8, "H[2]", T);
    near(hom.H[5], 0.0, 1e-8, "H[5]", T);
    near(hom.rms_residual, 0.0, 1e-8, "rms", T);
}

static void test_H4() {
    std::cout << "\n[H4] Rotazione 30°\n";
    const std::string T = "H4";
    const double c30 = std::sqrt(3.0) / 2.0;
    const double s30 = 0.5;
    auto pts = make_affine_pts(c30, -s30, 0, s30, c30, 0,
                               grid_pts(100, 100, 50, 3, 3));
    auto hom = compute_homography(pts);
    near(hom.H[0],  c30, 1e-7, "H[0]=cos30", T);
    near(hom.H[1], -s30, 1e-7, "H[1]=-sin30", T);
    near(hom.H[3],  s30, 1e-7, "H[3]=sin30", T);
    near(hom.H[4],  c30, 1e-7, "H[4]=cos30", T);
    near(hom.rms_residual, 0.0, 1e-7, "rms", T);
}

static void test_H5() {
    std::cout << "\n[H5] Prospettiva + apply_homography\n";
    const std::string T = "H5";
    std::array<double, 9> H_ref = {1, 0, 0,
                                    0, 1, 0,
                                    1e-4, 0, 1};
    auto disp = grid_pts(100, 100, 50, 3, 3);
    auto pts  = make_proj_pts(H_ref, disp);
    auto hom  = compute_homography(pts);

    // Verifica su un punto non nel training
    double xd = 175.0, yd = 225.0;
    double w_ref  = H_ref[6] * xd + H_ref[7] * yd + H_ref[8];
    double xs_ref = (H_ref[0] * xd + H_ref[1] * yd + H_ref[2]) / w_ref;
    double ys_ref = (H_ref[3] * xd + H_ref[4] * yd + H_ref[5]) / w_ref;
    double pix_mm = 0.01;
    auto phys = apply_homography(hom, xd, yd, pix_mm);
    near(phys.x_mm, xs_ref * pix_mm, 1e-6, "x_mm", T);
    near(phys.y_mm, ys_ref * pix_mm, 1e-6, "y_mm", T);
    near(hom.rms_residual, 0.0, 1e-7, "rms training", T);
}

static void test_H6() {
    std::cout << "\n[H6] Eccezione con < 4 punti\n";
    const std::string T = "H6";
    std::vector<CalibPoint> three(3);
    three[0] = {0, 0, 0, 0};
    three[1] = {1, 0, 1, 0};
    three[2] = {0, 1, 0, 1};
    bool threw = false;
    try { compute_homography(three); }
    catch (const std::invalid_argument&) { threw = true; }
    check(threw, "3 punti → invalid_argument", T);

    std::vector<CalibPoint> zero;
    threw = false;
    try { compute_homography(zero); }
    catch (const std::invalid_argument&) { threw = true; }
    check(threw, "0 punti → invalid_argument", T);
}

static void test_H7() {
    std::cout << "\n[H7] 25 punti con rumore sigma=0.1 px\n";
    const std::string T = "H7";
    const double sigma = 0.1;
    auto disp = grid_pts(50, 50, 50, 5, 5);
    auto clean = make_affine_pts(1.8, 0.1, 20, -0.05, 1.6, -10, disp);
    unsigned seed = 42u;
    for (auto& p : clean) {
        p.sens_x += sigma * lcg_gauss(seed);
        p.sens_y += sigma * lcg_gauss(seed);
    }
    auto hom = compute_homography(clean);
    check(hom.rms_residual < 3.0 * sigma, "rms < 3*sigma", T);
    check(hom.n_points == 25, "n_points=25", T);
}

// ── Cluster S: save/load JSON ─────────────────────────────────────────────────

static void test_S1() {
    std::cout << "\n[S1] Roundtrip save→load\n";
    const std::string T = "S1";
    Homography H_orig;
    H_orig.H = {1.5, 0.1, -20.3, -0.05, 1.4, 15.7, 1e-5, 2e-5, 1.0};
    H_orig.rms_residual = 0.123456789;
    H_orig.n_points = 17;

    std::filesystem::path tmp = "/tmp/riptide_test_H.json";
    save_homography(H_orig, tmp);
    auto H_loaded = load_homography(tmp);

    for (int i = 0; i < 9; ++i) {
        std::ostringstream lbl;
        lbl << "H[" << i << "]";
        near(H_loaded.H[static_cast<size_t>(i)],
             H_orig.H[static_cast<size_t>(i)], 1e-12, lbl.str(), T);
    }
    near(H_loaded.rms_residual, H_orig.rms_residual, 1e-12, "rms", T);
    check(H_loaded.n_points == 17, "n_points=17", T);
}

static void test_S2() {
    std::cout << "\n[S2] Creazione directory nidificate\n";
    const std::string T = "S2";
    Homography H;
    H.H.fill(0.0);
    H.H[0] = H.H[4] = H.H[8] = 1.0;
    std::filesystem::path nested = "/tmp/riptide_test_dir/sub/H.json";
    save_homography(H, nested);
    check(std::filesystem::exists(nested), "file creato in path nidificata", T);
}

static void test_S3() {
    std::cout << "\n[S3] Load file inesistente → eccezione\n";
    const std::string T = "S3";
    bool threw = false;
    try { load_homography("/tmp/riptide_NONEXISTENT_XXXXX.json"); }
    catch (const std::runtime_error&) { threw = true; }
    check(threw, "file mancante → runtime_error", T);
}

static void test_S4() {
    std::cout << "\n[S4] Pipeline compute→save→load→apply\n";
    const std::string T = "S4";
    // sens = 2*disp + (10, 5)
    auto pts = make_affine_pts(2, 0, 10, 0, 2, 5,
                               grid_pts(50, 50, 50, 3, 3));
    auto H_comp = compute_homography(pts);
    std::filesystem::path tmp = "/tmp/riptide_test_H_s4.json";
    save_homography(H_comp, tmp);
    auto H_load = load_homography(tmp);

    // Punto verifica: (150, 200) → sens (310, 405) → mm con pixsize=0.01
    double pix_mm = 0.01;
    auto phys = apply_homography(H_load, 150.0, 200.0, pix_mm);
    near(phys.x_mm, 310.0 * pix_mm, 1e-5, "x_mm", T);
    near(phys.y_mm, 405.0 * pix_mm, 1e-5, "y_mm", T);
}

// ── Cluster D: detect_calibration_dots ───────────────────────────────────────

static void test_D1() {
    std::cout << "\n[D1] Griglia 3×2 gaussiane perfette, scala 1:1\n";
    const std::string T = "D1";
    auto img = make_gauss_image(DISP_W, DISP_H, STD_NODES, 100.0, 3.0);
    auto pts = detect_calibration_dots(img, DISP_W, DISP_H,
                                       GRID_STEP, DISP_W, DISP_H);
    check(pts.size() == 6, "trovati 6 dot", T);
    // Ordina per (disp_y, disp_x) per confronto deterministico
    std::sort(pts.begin(), pts.end(), [](const CalibPoint& a, const CalibPoint& b) {
        return a.disp_y < b.disp_y || (a.disp_y == b.disp_y && a.disp_x < b.disp_x);
    });
    for (size_t i = 0; i < pts.size() && i < STD_NODES.size(); ++i) {
        std::ostringstream lx, ly;
        lx << "dot[" << i << "].sens_x";
        ly << "dot[" << i << "].sens_y";
        near(pts[i].sens_x, STD_NODES[i].first,  1e-6, lx.str(), T);
        near(pts[i].sens_y, STD_NODES[i].second, 1e-6, ly.str(), T);
    }
}

static void test_D2() {
    std::cout << "\n[D2] Immagine zero: nessun dot\n";
    const std::string T = "D2";
    std::vector<double> img(static_cast<size_t>(DISP_W * DISP_H), 0.0);
    auto pts = detect_calibration_dots(img, DISP_W, DISP_H,
                                       GRID_STEP, DISP_W, DISP_H);
    check(pts.empty(), "0 dot trovati", T);
}

static void test_D3() {
    std::cout << "\n[D3] Gaussiane spostate (+3.7, -2.1)\n";
    const std::string T = "D3";
    const double ox = 3.7, oy = -2.1;
    std::vector<std::pair<double, double>> shifted;
    for (const auto& [cx, cy] : STD_NODES)
        shifted.push_back({cx + ox, cy + oy});
    auto img = make_gauss_image(DISP_W, DISP_H, shifted, 100.0, 3.0);
    auto pts = detect_calibration_dots(img, DISP_W, DISP_H,
                                       GRID_STEP, DISP_W, DISP_H);
    check(pts.size() == 6, "trovati 6 dot", T);
    std::sort(pts.begin(), pts.end(), [](const CalibPoint& a, const CalibPoint& b) {
        return a.disp_y < b.disp_y || (a.disp_y == b.disp_y && a.disp_x < b.disp_x);
    });
    for (size_t i = 0; i < pts.size() && i < STD_NODES.size(); ++i) {
        std::ostringstream lx, ly;
        lx << "dot[" << i << "].sens_x≈" << STD_NODES[i].first + ox;
        ly << "dot[" << i << "].sens_y≈" << STD_NODES[i].second + oy;
        near(pts[i].sens_x, STD_NODES[i].first  + ox, 1e-4, lx.str(), T);
        near(pts[i].sens_y, STD_NODES[i].second + oy, 1e-4, ly.str(), T);
    }
}

static void test_D4() {
    std::cout << "\n[D4] Picco sotto soglia 1e-6: nessun dot\n";
    const std::string T = "D4";
    // amplitude=1e-7 → peak < 1e-6 → rifiuto immediato nella fit_dot_centroid
    auto img = make_gauss_image(DISP_W, DISP_H, STD_NODES, 1e-7, 3.0);
    auto pts = detect_calibration_dots(img, DISP_W, DISP_H,
                                       GRID_STEP, DISP_W, DISP_H);
    check(pts.empty(), "0 dot trovati (peak < 1e-6)", T);
}

static void test_D5() {
    std::cout << "\n[D5] Scala non unitaria: sensor 300×225, display 400×300\n";
    const std::string T = "D5";
    // scale=0.75: gaussiane a 75px di distanza nel sensore, finestra half=50 → nessuna sovrapposizione
    const int SW = 300, SH = 225;
    std::vector<std::pair<double, double>> sens_centers;
    for (const auto& [gx, gy] : STD_NODES)
        sens_centers.push_back({gx * 0.75, gy * 0.75});
    auto img = make_gauss_image(SW, SH, sens_centers, 100.0, 3.0);
    auto pts = detect_calibration_dots(img, SW, SH,
                                       GRID_STEP, DISP_W, DISP_H);
    check(pts.size() == 6, "trovati 6 dot", T);
    std::sort(pts.begin(), pts.end(), [](const CalibPoint& a, const CalibPoint& b) {
        return a.disp_y < b.disp_y || (a.disp_y == b.disp_y && a.disp_x < b.disp_x);
    });
    for (size_t i = 0; i < pts.size() && i < STD_NODES.size(); ++i) {
        std::ostringstream lx, ly;
        lx << "dot[" << i << "].sens_x≈" << STD_NODES[i].first * 0.75;
        ly << "dot[" << i << "].sens_y≈" << STD_NODES[i].second * 0.75;
        near(pts[i].sens_x, STD_NODES[i].first  * 0.75, 1e-4, lx.str(), T);
        near(pts[i].sens_y, STD_NODES[i].second * 0.75, 1e-4, ly.str(), T);
    }
}

static void test_D6() {
    std::cout << "\n[D6] Pipeline detect → compute_homography\n";
    const std::string T = "D6";
    auto img = make_gauss_image(DISP_W, DISP_H, STD_NODES, 100.0, 3.0);
    auto pts = detect_calibration_dots(img, DISP_W, DISP_H,
                                       GRID_STEP, DISP_W, DISP_H);
    if (!check(pts.size() >= 4, "almeno 4 dot per DLT", T)) return;
    auto hom = compute_homography(pts);
    check(hom.rms_residual < 0.01, "rms pipeline < 0.01 px", T);
}

// ── Cluster Q: load_Q_sim ─────────────────────────────────────────────────────

static const std::filesystem::path Q_TSV = "/tmp/riptide_test_qsim.tsv";

static void setup_Q_tsv() {
    write_tmp_tsv(
        "x1\tx2\tq\n"
        "100.0\t200.0\t1.234\n"
        "110.0\t210.0\t5.678\n"
        "50.0\t60.0\t9.999\n",
        Q_TSV.string());
}

static void test_Q1() {
    std::cout << "\n[Q1] Match esatto\n";
    const std::string T = "Q1";
    double q = load_Q_sim(Q_TSV, 100.0, 200.0);
    near(q, 1.234, 1e-12, "Q_sim", T);
}

static void test_Q2() {
    std::cout << "\n[Q2] Match con tolleranza (Δ=0.3 < 0.5)\n";
    const std::string T = "Q2";
    double q = load_Q_sim(Q_TSV, 100.3, 200.4);
    near(q, 1.234, 1e-12, "Q_sim con tolleranza", T);
}

static void test_Q3() {
    std::cout << "\n[Q3] Nessun match → NaN\n";
    const std::string T = "Q3";
    double q = load_Q_sim(Q_TSV, 200.0, 300.0);
    check(is_nan_bits(q), "nessun match → NaN", T);
}

static void test_Q4() {
    std::cout << "\n[Q4] File non esiste → NaN\n";
    const std::string T = "Q4";
    double q = load_Q_sim("/tmp/riptide_NONEXISTENT_qsim.tsv", 1.0, 2.0);
    check(is_nan_bits(q), "file mancante → NaN", T);
}

static void test_Q5() {
    std::cout << "\n[Q5] Righe commentate e vuote saltate\n";
    const std::string T = "Q5";
    auto p = write_tmp_tsv(
        "x1\tx2\tq\n"
        "# commento\n"
        "\n"
        "50.0\t60.0\t9.999\n",
        "/tmp/riptide_test_qsim_comments.tsv");
    double q = load_Q_sim(p, 50.0, 60.0);
    near(q, 9.999, 1e-12, "Q_sim dopo commenti", T);
}

static void test_Q6() {
    std::cout << "\n[Q6] Confine tolleranza 0.5\n";
    const std::string T = "Q6";
    // Δ=0.499 → match
    double q_in = load_Q_sim(Q_TSV, 100.499, 200.0);
    check(!is_nan_bits(q_in), "Δ=0.499 → match trovato", T);
    // Δ=0.501 → no match
    double q_out = load_Q_sim(Q_TSV, 100.501, 200.0);
    check(is_nan_bits(q_out), "Δ=0.501 → nessun match", T);
}

// ── main ──────────────────────────────────────────────────────────────────────

int main() {
    std::cout << "test_exp3_homography — riptide exp3 unit tests\n";

    test_H1();
    test_H2();
    test_H3();
    test_H4();
    test_H5();
    test_H6();
    test_H7();

    test_S1();
    test_S2();
    test_S3();
    test_S4();

    test_D1();
    test_D2();
    test_D3();
    test_D4();
    test_D5();
    test_D6();

    setup_Q_tsv();
    test_Q1();
    test_Q2();
    test_Q3();
    test_Q4();
    test_Q5();
    test_Q6();

    std::cout << "\nTEST SUMMARY: " << g_n_pass << " PASS, " << g_n_fail << " FAIL\n";
    return (g_n_fail == 0) ? 0 : 1;
}
