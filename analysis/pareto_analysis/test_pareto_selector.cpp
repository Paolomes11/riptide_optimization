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
    // P1: eta=1.0, Q=1.0, DoF=0.5, M_abs_err=0.0
    //   Mtot = 0.35*1.0 - 0.40*1.0 + 0.15*0.5 - 0.10*0.0 = 0.025
    // P2: eta=0.8, Q=0.6, DoF=1.0, M_abs_err=0.5
    //   Mtot = 0.35*0.8 - 0.40*0.6 + 0.15*1.0 - 0.10*0.5 = 0.140
    // P3: eta=0.6, Q=0.3, DoF=0.8, M_abs_err=1.0
    //   Mtot = 0.35*0.6 - 0.40*0.3 + 0.15*0.8 - 0.10*1.0 = 0.110
    //
    // Ordine atteso: P2 (0.140) > P3 (0.110) > P1 (0.025)
    // Tutti e 3 sono sul fronte (nessuno domina l'altro in (eta, Q))

    std::vector<ConfigData> cfgs;
    cfgs.push_back(make_config(1.0, 1.0, 1.0, 1.0, 180.0, 0.5, 1.0, 0.0));  // P1
    cfgs.push_back(make_config(2.0, 2.0, 0.8, 0.6, 180.0, 1.0, 1.0, 0.5));  // P2
    cfgs.push_back(make_config(3.0, 3.0, 0.6, 0.3, 180.0, 0.8, 1.0, 1.0));  // P3

    WeightConfig wc;  // pesi di default: w_eta=0.35, w_Q=0.40, w_dof=0.15, w_M=0.10

    compute_mtot(cfgs, wc);

    // Verifica Mtot calcolati (tolleranza 1e-9)
    near(cfgs[0].Mtot, 0.025, 1e-9, "Mtot P1 = 0.025");
    near(cfgs[1].Mtot, 0.140, 1e-9, "Mtot P2 = 0.140");
    near(cfgs[2].Mtot, 0.110, 1e-9, "Mtot P3 = 0.110");

    compute_pareto_front(cfgs, wc);

    // Tutti e 3 sul fronte
    check(cfgs[0].on_pareto, "P1 sul fronte");
    check(cfgs[1].on_pareto, "P2 sul fronte");
    check(cfgs[2].on_pareto, "P3 sul fronte");

    // Ranking: P2=1, P3=2, P1=3
    check(cfgs[1].pareto_rank == 1, "P2 ha rank 1 (Mtot massimo)");
    check(cfgs[2].pareto_rank == 2, "P3 ha rank 2");
    check(cfgs[0].pareto_rank == 3, "P1 ha rank 3 (Mtot minimo)");
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

    std::cout << "\nTEST SUMMARY: " << g_n_pass << " PASS, " << g_n_fail << " FAIL\n";
    return (g_n_fail == 0) ? 0 : 1;
}
