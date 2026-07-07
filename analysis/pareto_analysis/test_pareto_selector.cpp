/*
 * test_pareto_selector — Test unitari per pareto_core.hpp
 *
 * T1: fronte di Pareto su caso 4-punti noto
 * T2: apply_eta_filter
 * T3: apply_focus_filter
 * T4: compute_mtot + pareto_rank ordering
 * T5: coords_match con tolleranza
 * T6: apply_dof_filter
 * T7: compute_pareto_front casi degeneri (singolo punto, punti identici)
 * T8: apply_mtot con normalizzatori precalcolati == compute_mtot (regressione)
 * T9: generate_barycentric_grid(0.05): 1771 punti, ognuno somma a 1.0
 * T10: sweep pesi sul caso a 3 punti di T4 == vincitore rank-1 di T4
 * T11: ternary_to_xy sanity sui vertici
 * T12: apply_focus_filter con mobile_focus=true è no-op
 * T13: resolve_dof_bounds — valori TSV vs fallback ±DoF/2 attorno a x_focus
 * T14: reflect_across_line — apici dei 3 lembi ribaltati del net del tetraedro
 *      a due a due a distanza esatta 2.0 (tassellazione senza buchi/sovrapposizioni)
 * T15: solve_face_affine — vertici canonici locali mappati esattamente sui
 *      vertici globali noti
 *
 * Ritorna 0 se tutti i test passano, 1 altrimenti.
 */

#include "pareto_core.hpp"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using namespace riptide::pareto;

static int g_n_pass = 0;
static int g_n_fail = 0;

static bool check(bool condition, const std::string& msg) {
    if (condition) {
        std::cout << "  [PASS] " << msg << "\n";
        ++g_n_pass;
    } else {
        std::cerr << "  [FAIL] " << msg << "\n";
        ++g_n_fail;
    }
    return condition;
}

static bool near(double a, double b, double tol, const std::string& label) {
    return check(std::abs(a - b) < tol, label + " (" + std::to_string(a) + " vs " + std::to_string(b) + ")");
}

// ── helpers per costruire ConfigData minimali ────────────────────────────────

static ConfigData make_config(double x1, double x2, double eta, double Q,
                               double x_focus = 180.0, double DoF = 50.0,
                               double M = 1.0, double M_abs_err = 0.0) {
    ConfigData c;
    c.x1 = x1; c.x2 = x2;
    c.eta = eta; c.Q = Q;
    c.x_focus = x_focus; c.DoF = DoF;
    c.M = M; c.M_abs_err = M_abs_err;
    c.chi2 = 1.0;
    return c;
}

// ── T1: fronte di Pareto su caso 4-punti noto ────────────────────────────────

static void test_T1() {
    std::cout << "\n[T1] Fronte di Pareto su caso 4-punti noto\n";

    // A: alta eta, alta Q → sul fronte (alta eta unica)
    // B: bassa eta, bassa Q → sul fronte (bassa Q unica)
    // C: dominata da A (stessa eta, Q peggiore)
    // D: dominata da B (eta peggiore, Q peggiore)
    std::vector<ConfigData> cfgs = {
        make_config(1.0, 1.0, 0.9, 1.0),  // A
        make_config(2.0, 2.0, 0.5, 0.5),  // B
        make_config(3.0, 3.0, 0.9, 1.5),  // C — A domina C
        make_config(4.0, 4.0, 0.4, 0.6),  // D — B domina D
    };

    WeightConfig wc;
    compute_pareto_front(cfgs, wc);

    check(cfgs[0].on_pareto,  "A è sul fronte");
    check(cfgs[1].on_pareto,  "B è sul fronte");
    check(!cfgs[2].on_pareto, "C non è sul fronte (dominata da A)");
    check(!cfgs[3].on_pareto, "D non è sul fronte (dominata da B)");
}

// ── T2: apply_eta_filter ─────────────────────────────────────────────────────

static void test_T2() {
    std::cout << "\n[T2] apply_eta_filter con eta_frac=0.75\n";

    // 10 config con eta = 0.1, 0.2, ..., 1.0
    // eta_max = 1.0, soglia = 0.75 → passano eta ∈ {0.8, 0.9, 1.0}
    std::vector<ConfigData> cfgs;
    for (int i = 1; i <= 10; ++i)
        cfgs.push_back(make_config((double)i, 0.0, i * 0.1, 1.0, 180.0, 50.0));

    FilterConfig fc;
    fc.eta_frac  = 0.75;
    fc.x_det     = 180.0;
    fc.focus_tol = 1000.0;  // disabilitato
    fc.dof_min   = 0.0;

    auto result = apply_eta_filter(cfgs, fc);

    check(result.size() == 3, "N sopravvissuti = 3 (eta ∈ {0.8, 0.9, 1.0})");
    if (result.size() == 3) {
        near(result[0].eta, 0.8, 1e-9, "result[0].eta = 0.8");
        near(result[1].eta, 0.9, 1e-9, "result[1].eta = 0.9");
        near(result[2].eta, 1.0, 1e-9, "result[2].eta = 1.0");
    }
}

// ── T3: apply_focus_filter ───────────────────────────────────────────────────

static void test_T3() {
    std::cout << "\n[T3] apply_focus_filter con x_det=180, focus_tol=15\n";

    // 5 config con x_focus ∈ {160, 170, 180, 190, 200}
    // |x_focus - 180| ≤ 15 → passano {170, 180, 190}; droppate {160, 200}
    std::vector<double> xf = {160.0, 170.0, 180.0, 190.0, 200.0};
    std::vector<ConfigData> cfgs;
    for (double f : xf)
        cfgs.push_back(make_config(0.0, 0.0, 0.9, 1.0, f, 50.0));

    FilterConfig fc;
    fc.eta_frac  = 0.0;  // disabilitato
    fc.x_det     = 180.0;
    fc.focus_tol = 15.0;
    fc.dof_min   = 0.0;

    auto result = apply_focus_filter(cfgs, fc);

    check(result.size() == 3, "N sopravvissuti = 3 (x_focus ∈ {170,180,190})");
    if (result.size() == 3) {
        near(result[0].x_focus, 170.0, 1e-9, "result[0].x_focus = 170");
        near(result[1].x_focus, 180.0, 1e-9, "result[1].x_focus = 180");
        near(result[2].x_focus, 190.0, 1e-9, "result[2].x_focus = 190");
    }
}

// ── T4: compute_mtot + pareto_rank ordering ──────────────────────────────────

static void test_T4() {
    std::cout << "\n[T4] compute_mtot + pareto_rank: P2 > P3 > P1\n";

    // Valori già nell'intervallo [0,1] con massimo = 1.0 per ogni dimensione
    // → normalizzazione interna non altera i valori attesi
    //
    // Formula: Mtot = nw_eta*(η/η_max) + nw_Q*(1−Q/Q_max) + nw_dof*(DoF/DoF_max) + nw_M*(1−M/M_max)
    // Tutti i termini ∈ [0,1]; Mtot ∈ [0,1].  [vecchia formula con sottrazione: RIMOSSA]
    //
    // max values across {P1,P2,P3}: eta_max=1.0, Q_max=1.0, DoF_max=1.0, M_max=1.0
    //
    // P1: eta=1.0, Q=1.0, DoF=0.5, M_abs_err=0.0
    //   Mtot = 0.35*1.0 + 0.40*(1−1.0) + 0.15*0.5 + 0.10*(1−0.0) = 0.525
    // P2: eta=0.8, Q=0.6, DoF=1.0, M_abs_err=0.5
    //   Mtot = 0.35*0.8 + 0.40*(1−0.6) + 0.15*1.0 + 0.10*(1−0.5) = 0.640
    // P3: eta=0.6, Q=0.3, DoF=0.8, M_abs_err=1.0
    //   Mtot = 0.35*0.6 + 0.40*(1−0.3) + 0.15*0.8 + 0.10*(1−1.0) = 0.610
    //
    // Ordine atteso: P2 (0.640) > P3 (0.610) > P1 (0.525)
    // Tutti e 3 sono sul fronte (nessuno domina l'altro in (eta, Q, ΔM))

    std::vector<ConfigData> cfgs;
    cfgs.push_back(make_config(1.0, 1.0, 1.0, 1.0, 180.0, 0.5, 1.0, 0.0));  // P1
    cfgs.push_back(make_config(2.0, 2.0, 0.8, 0.6, 180.0, 1.0, 1.0, 0.5));  // P2
    cfgs.push_back(make_config(3.0, 3.0, 0.6, 0.3, 180.0, 0.8, 1.0, 1.0));  // P3

    WeightConfig wc;  // pesi di default: w_eta=0.35, w_Q=0.40, w_dof=0.15, w_M=0.10

    compute_mtot(cfgs, wc);

    // Verifica Mtot calcolati (tolleranza 1e-9)
    near(cfgs[0].Mtot, 0.525, 1e-9, "Mtot P1 = 0.525");
    near(cfgs[1].Mtot, 0.640, 1e-9, "Mtot P2 = 0.640");
    near(cfgs[2].Mtot, 0.610, 1e-9, "Mtot P3 = 0.610");

    compute_pareto_front(cfgs, wc);

    // Tutti e 3 sul fronte
    check(cfgs[0].on_pareto, "P1 sul fronte");
    check(cfgs[1].on_pareto, "P2 sul fronte");
    check(cfgs[2].on_pareto, "P3 sul fronte");

    // Ranking: P2=1 (0.640), P3=2 (0.610), P1=3 (0.525)
    check(cfgs[1].pareto_rank == 1, "P2 ha rank 1 (Mtot massimo = 0.640)");
    check(cfgs[2].pareto_rank == 2, "P3 ha rank 2 (Mtot = 0.610)");
    check(cfgs[0].pareto_rank == 3, "P1 ha rank 3 (Mtot minimo = 0.525)");
}

// ── T5: coords_match con tolleranza ──────────────────────────────────────────

static void test_T5() {
    std::cout << "\n[T5] coords_match con tolleranza 1e-3 mm\n";

    std::vector<std::pair<double,double>> nominals = {
        {100.0, 200.0}, {110.0, 210.0}, {120.0, 220.0}
    };

    // Sub-test A: delta = 0.0005 < 1e-3 → tutti matchano
    {
        double delta = 0.0005;
        int matches = 0;
        for (const auto& [x1, x2] : nominals)
            if (coords_match(x1, x2, x1 + delta, x2 + delta))
                ++matches;
        check(matches == 3, "Sub-test A (delta=0.0005): tutti e 3 matchano");
    }

    // Sub-test B: delta = 0.002 > 1e-3 → nessuno matcha
    {
        double delta = 0.002;
        int matches = 0;
        for (const auto& [x1, x2] : nominals)
            if (coords_match(x1, x2, x1 + delta, x2 + delta))
                ++matches;
        check(matches == 0, "Sub-test B (delta=0.002): nessun match");
    }
}

// ── T6: apply_dof_filter ─────────────────────────────────────────────────────

static void test_T6() {
    std::cout << "\n[T6] apply_dof_filter\n";

    // 5 config con DoF ∈ {0, 10, 20, 30, 40}, dof_min=15 → passano {20,30,40}
    std::vector<ConfigData> cfgs;
    for (int d : {0, 10, 20, 30, 40})
        cfgs.push_back(make_config(0.0, 0.0, 0.9, 1.0, 180.0, (double)d));

    FilterConfig fc;
    fc.dof_min = 15.0;

    auto result = apply_dof_filter(cfgs, fc);

    check(result.size() == 3, "N sopravvissuti = 3 (DoF ∈ {20,30,40})");
    if (result.size() == 3) {
        near(result[0].DoF, 20.0, 1e-9, "result[0].DoF = 20");
        near(result[1].DoF, 30.0, 1e-9, "result[1].DoF = 30");
        near(result[2].DoF, 40.0, 1e-9, "result[2].DoF = 40");
    }

    // Path disabilitato: dof_min=0 → tutti passano
    fc.dof_min = 0.0;
    auto all = apply_dof_filter(cfgs, fc);
    check(all.size() == 5, "dof_min=0: tutti e 5 passano");
}

// ── T7: compute_pareto_front casi degeneri ────────────────────────────────────

static void test_T7() {
    std::cout << "\n[T7] compute_pareto_front casi degeneri\n";

    WeightConfig wc;

    // Sub-test A: singolo punto → sul fronte con rank=1
    {
        std::vector<ConfigData> cfgs = { make_config(1.0, 1.0, 0.8, 0.5) };
        compute_pareto_front(cfgs, wc);
        check(cfgs[0].on_pareto,       "A (singolo punto): on_pareto=true");
        check(cfgs[0].pareto_rank == 1, "A (singolo punto): pareto_rank=1");
    }

    // Sub-test B: due punti identici (stessi eta, Q) → nessuno domina l'altro
    {
        std::vector<ConfigData> cfgs = {
            make_config(1.0, 1.0, 0.7, 0.6),
            make_config(2.0, 2.0, 0.7, 0.6),
        };
        compute_pareto_front(cfgs, wc);
        check(cfgs[0].on_pareto, "B (identici): cfgs[0] sul fronte");
        check(cfgs[1].on_pareto, "B (identici): cfgs[1] sul fronte");
    }
}

// ── T8: apply_mtot con normalizzatori precalcolati == compute_mtot ──────────

static void test_T8() {
    std::cout << "\n[T8] apply_mtot(normalizers) == compute_mtot (regressione bit-per-bit)\n";

    std::vector<ConfigData> cfgs_a;
    cfgs_a.push_back(make_config(1.0, 1.0, 1.0, 1.0, 180.0, 0.5, 1.0, 0.0));
    cfgs_a.push_back(make_config(2.0, 2.0, 0.8, 0.6, 180.0, 1.0, 1.0, 0.5));
    cfgs_a.push_back(make_config(3.0, 3.0, 0.6, 0.3, 180.0, 0.8, 1.0, 1.0));
    std::vector<ConfigData> cfgs_b = cfgs_a;

    WeightConfig wc;
    compute_mtot(cfgs_a, wc);

    const MtotNormalizers norm = compute_mtot_normalizers(cfgs_b);
    apply_mtot(cfgs_b, norm, wc);

    for (size_t i = 0; i < cfgs_a.size(); ++i)
        near(cfgs_a[i].Mtot, cfgs_b[i].Mtot, 1e-12,
             "Mtot identico per punto " + std::to_string(i));
}

// ── T9: generate_barycentric_grid ────────────────────────────────────────────

static void test_T9() {
    std::cout << "\n[T9] generate_barycentric_grid(0.05): 1771 punti, somme=1.0\n";

    auto grid = generate_barycentric_grid(0.05);
    check(grid.size() == 1771, "grid.size() == 1771 (C(23,3))");

    bool all_sum_to_one = true;
    for (const auto& wc : grid) {
        const double s = wc.w_eta + wc.w_Q + wc.w_dof + wc.w_M;
        if (std::abs(s - 1.0) >= 1e-9) { all_sum_to_one = false; break; }
    }
    check(all_sum_to_one, "ogni vettore di pesi somma esattamente a 1.0");
}

// ── T10: sweep pesi == vincitore rank-1 di T4 ────────────────────────────────

static void test_T10() {
    std::cout << "\n[T10] run_weight_sweep sul caso T4: vincitore == P2 (rank 1)\n";

    std::vector<ConfigData> cfgs;
    cfgs.push_back(make_config(1.0, 1.0, 1.0, 1.0, 180.0, 0.5, 1.0, 0.0));  // P1
    cfgs.push_back(make_config(2.0, 2.0, 0.8, 0.6, 180.0, 1.0, 1.0, 0.5));  // P2
    cfgs.push_back(make_config(3.0, 3.0, 0.6, 0.3, 180.0, 0.8, 1.0, 1.0));  // P3

    const MtotNormalizers norm = compute_mtot_normalizers(cfgs);

    WeightConfig wc;  // default: w_eta=0.35, w_Q=0.40, w_dof=0.15, w_M=0.10
    auto grid = generate_barycentric_grid(0.05);

    // Il punto di default è esattamente sulla griglia con step=0.05
    // (0.35,0.40,0.15,0.10)*20 = (7,8,3,2), somma=20=n.
    const bool grid_has_default = std::any_of(grid.begin(), grid.end(), [&wc](const WeightConfig& g) {
        return std::abs(g.w_eta - wc.w_eta) < 1e-9 && std::abs(g.w_Q - wc.w_Q) < 1e-9 &&
               std::abs(g.w_dof - wc.w_dof) < 1e-9 && std::abs(g.w_M - wc.w_M) < 1e-9;
    });
    check(grid_has_default, "il peso di default è presente nella griglia step=0.05");

    auto sweep = run_weight_sweep(cfgs, norm, grid);

    const SweepPoint* found = nullptr;
    for (const auto& sp : sweep) {
        if (std::abs(sp.wc.w_eta - wc.w_eta) < 1e-9 &&
            std::abs(sp.wc.w_Q   - wc.w_Q)   < 1e-9 &&
            std::abs(sp.wc.w_dof - wc.w_dof) < 1e-9 &&
            std::abs(sp.wc.w_M   - wc.w_M)   < 1e-9) {
            found = &sp;
            break;
        }
    }
    if (check(found != nullptr, "punto di sweep per il peso di default trovato")) {
        near(found->x1_winner, 2.0, 1e-9, "x1_winner == P2.x1 (2.0)");
        near(found->x2_winner, 2.0, 1e-9, "x2_winner == P2.x2 (2.0)");
        near(found->Mtot_winner, 0.640, 1e-9, "Mtot_winner == 0.640 (P2, rank 1 in T4)");
    }
}

// ── T11: ternary_to_xy sanity sui vertici ────────────────────────────────────

static void test_T11() {
    std::cout << "\n[T11] ternary_to_xy sanity sui vertici\n";

    auto [xa, ya] = ternary_to_xy(1.0, 0.0, 0.0);
    near(xa, 0.0, 1e-9, "vertice A=(1,0,0) -> x=0");
    near(ya, 0.0, 1e-9, "vertice A=(1,0,0) -> y=0");

    auto [xb, yb] = ternary_to_xy(0.0, 1.0, 0.0);
    near(xb, 1.0, 1e-9, "vertice B=(0,1,0) -> x=1");
    near(yb, 0.0, 1e-9, "vertice B=(0,1,0) -> y=0");

    auto [xc, yc] = ternary_to_xy(0.0, 0.0, 1.0);
    near(xc, 0.5, 1e-9, "vertice C=(0,0,1) -> x=0.5");
    near(yc, std::sqrt(3.0) / 2.0, 1e-9, "vertice C=(0,0,1) -> y=sqrt(3)/2");
}

// ── T12: apply_focus_filter con mobile_focus=true è no-op ───────────────────

static void test_T12() {
    std::cout << "\n[T12] apply_focus_filter con mobile_focus=true è no-op\n";

    // Stesso set di T3 (x_focus fuori tolleranza per 160/200), ma con mobile_focus
    // attivo il filtro non deve scartare nulla, indipendentemente da x_det/focus_tol.
    std::vector<double> xf = {160.0, 170.0, 180.0, 190.0, 200.0};
    std::vector<ConfigData> cfgs;
    for (double f : xf)
        cfgs.push_back(make_config(0.0, 0.0, 0.9, 1.0, f, 50.0));

    FilterConfig fc;
    fc.x_det        = 180.0;
    fc.focus_tol    = 15.0;
    fc.mobile_focus = true;

    auto result = apply_focus_filter(cfgs, fc);
    check(result.size() == cfgs.size(), "N sopravvissuti == N input (nessun filtro applicato)");
}

// ── T13: resolve_dof_bounds ──────────────────────────────────────────────────

static void test_T13() {
    std::cout << "\n[T13] resolve_dof_bounds: valori TSV vs fallback ±DoF/2\n";

    // Colonne x_lo/x_hi presenti nel TSV: usate direttamente (asimmetriche)
    ConfigData cd1 = make_config(0.0, 0.0, 0.9, 1.0, /*x_focus*/180.0, /*DoF*/50.0);
    resolve_dof_bounds(cd1, /*has_x_lo_hi*/true, /*x_lo_raw*/160.0, /*x_hi_raw*/205.0);
    near(cd1.x_lo, 160.0, 1e-9, "x_lo == valore TSV (160.0), non il fallback simmetrico");
    near(cd1.x_hi, 205.0, 1e-9, "x_hi == valore TSV (205.0), non il fallback simmetrico");

    // Colonne assenti (TSV vecchio): fallback simmetrico attorno a x_focus
    ConfigData cd2 = make_config(0.0, 0.0, 0.9, 1.0, /*x_focus*/180.0, /*DoF*/50.0);
    resolve_dof_bounds(cd2, /*has_x_lo_hi*/false, 0.0, 0.0);
    near(cd2.x_lo, 155.0, 1e-9, "x_lo == x_focus - DoF/2 (fallback, 180-25)");
    near(cd2.x_hi, 205.0, 1e-9, "x_hi == x_focus + DoF/2 (fallback, 180+25)");
}

// ── T14: reflect_across_line — apici del net a distanza esatta 2.0 ──────────

static double dist(std::pair<double, double> P, std::pair<double, double> Q) {
    return std::sqrt((P.first - Q.first) * (P.first - Q.first) +
                      (P.second - Q.second) * (P.second - Q.second));
}

static void test_T14() {
    std::cout << "\n[T14] reflect_across_line: apici dei lembi a distanza 2.0\n";

    const auto A = ternary_to_xy(1.0, 0.0, 0.0);
    const auto B = ternary_to_xy(0.0, 1.0, 0.0);
    const auto C = ternary_to_xy(0.0, 0.0, 1.0);

    // Lembo w_eta=0 (condivide BC): ribalta A attorno a BC.
    const auto apex_eta = reflect_across_line(A, B, C);
    // Lembo w_Q=0 (condivide AC): ribalta B attorno ad AC.
    const auto apex_Q = reflect_across_line(B, A, C);
    // Lembo w_dof=0 (condivide AB): ribalta C attorno ad AB.
    const auto apex_dof = reflect_across_line(C, A, B);

    near(dist(apex_eta, apex_Q), 2.0, 1e-9, "dist(apex_eta, apex_Q) == 2.0");
    near(dist(apex_eta, apex_dof), 2.0, 1e-9, "dist(apex_eta, apex_dof) == 2.0");
    near(dist(apex_Q, apex_dof), 2.0, 1e-9, "dist(apex_Q, apex_dof) == 2.0");
}

// ── T15: solve_face_affine — vertici canonici mappati sui vertici globali ────

static void test_T15() {
    std::cout << "\n[T15] solve_face_affine: L0/L1/L2 -> G0/G1/G2 esatti\n";

    const std::pair<double, double> G0{1.5, -2.7};
    const std::pair<double, double> G1{3.2, 0.4};
    const std::pair<double, double> G2{-0.8, 4.1};

    const FaceAffine fa = solve_face_affine(G0, G1, G2);
    const auto p0 = fa.apply(0.0, 0.0);
    const auto p1 = fa.apply(1.0, 0.0);
    const auto p2 = fa.apply(0.5, std::sqrt(3.0) / 2.0);

    near(p0.first, G0.first, 1e-9, "L0.x -> G0.x");
    near(p0.second, G0.second, 1e-9, "L0.y -> G0.y");
    near(p1.first, G1.first, 1e-9, "L1.x -> G1.x");
    near(p1.second, G1.second, 1e-9, "L1.y -> G1.y");
    near(p2.first, G2.first, 1e-9, "L2.x -> G2.x");
    near(p2.second, G2.second, 1e-9, "L2.y -> G2.y");
}

// ── main ─────────────────────────────────────────────────────────────────────

int main() {
    std::cout << "test_pareto_selector — riptide pareto core tests\n";

    test_T1();
    test_T2();
    test_T3();
    test_T4();
    test_T5();
    test_T6();
    test_T7();
    test_T8();
    test_T9();
    test_T10();
    test_T11();
    test_T12();
    test_T13();
    test_T14();
    test_T15();

    std::cout << "\nTEST SUMMARY: " << g_n_pass << " PASS, " << g_n_fail << " FAIL\n";
    return (g_n_fail == 0) ? 0 : 1;
}
