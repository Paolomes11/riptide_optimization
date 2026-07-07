#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace riptide::pareto {

struct ConfigData {
    double x1, x2;
    double eta;        // efficienza geometrica [0,1]
    double Q;          // qualità ottica (da minimizzare; colonna "metric" del TSV)
    double chi2;       // chi²/ndof (colonna "chi2_red" del TSV)
    double DoF;        // profondità di campo [mm]
    double M;          // magnificazione
    double M_abs_err;  // |M - m_target| (colonna "M_abs_err" del TSV)
    double x_focus;    // piano di fuoco [mm]
    double x_lo = 0.0; // limite inferiore intervallo DoF [mm] (x_focus - DoF/2 se non disponibile nel TSV)
    double x_hi = 0.0; // limite superiore intervallo DoF [mm] (x_focus + DoF/2 se non disponibile nel TSV)
    double EE80 = std::numeric_limits<double>::quiet_NaN(); // EE80 medio [mm] (da resolution TSV)
    bool   on_pareto   = false;
    double Mtot        = 0.0;
    int    pareto_rank = 0;  // 1=best Mtot sul fronte, 0=non sul fronte
};

struct FilterConfig {
    double eta_frac  = 0.75;   // soglia η relativa a η_max
    double x_det     = 180.0;  // posizione detector [mm]
    double focus_tol = 15.0;   // tolleranza fuoco [mm]
    double dof_min   = 0.0;    // DoF minima [mm] (0 = disabilitato)
    double ee80_max  = 0.0;    // EE80 massimo [mm] (0 = disabilitato)
    bool   mobile_focus = false; // se true, il filtro x_det/focus_tol (pensato per detector statico) è no-op
};

struct WeightConfig {
    double w_eta = 0.35;
    double w_Q   = 0.40;
    double w_dof = 0.15;
    double w_M   = 0.10;
};

// Matching coordinate con tolleranza
inline bool coords_match(double x1a, double x2a, double x1b, double x2b,
                         double tol = 1e-3) {
    return std::abs(x1a - x1b) < tol && std::abs(x2a - x2b) < tol;
}

// Chiave intera per il join su (x1,x2) via unordered_map (O(1) invece di scansione lineare).
// Sicura perché il passo di griglia (dx=1.0mm, vedi config/config.json) è molto più largo
// della risoluzione di arrotondamento (1e-4mm) usata qui.
inline int64_t coord_key(double x1, double x2) {
    const int64_t k1 = llround(x1 * 1e4);
    const int64_t k2 = llround(x2 * 1e4);
    const uint64_t u1 = static_cast<uint64_t>(k1) & 0xffffffffULL;
    const uint64_t u2 = static_cast<uint64_t>(k2) & 0xffffffffULL;
    return static_cast<int64_t>((u1 << 32) | u2);
}

// A domina B nel senso di Pareto (η↑, Q↓, ΔM↓):
// A.eta >= B.eta AND A.Q <= B.Q AND A.M_abs_err <= B.M_abs_err
// con almeno una disuguaglianza stretta
inline bool pareto_dominates(const ConfigData& A, const ConfigData& B) {
    return (A.eta >= B.eta && A.Q <= B.Q && A.M_abs_err <= B.M_abs_err)
        && (A.eta > B.eta  || A.Q < B.Q  || A.M_abs_err < B.M_abs_err);
}

inline std::vector<ConfigData> apply_eta_filter(const std::vector<ConfigData>& configs,
                                                 const FilterConfig& fc) {
    if (configs.empty()) return {};
    double eta_max = 0.0;
    for (const auto& c : configs)
        eta_max = std::max(eta_max, c.eta);
    double threshold = fc.eta_frac * eta_max;
    std::vector<ConfigData> result;
    for (const auto& c : configs)
        if (c.eta >= threshold)
            result.push_back(c);
    return result;
}

inline std::vector<ConfigData> apply_focus_filter(const std::vector<ConfigData>& configs,
                                                    const FilterConfig& fc) {
    // In fuoco mobile non esiste un detector a posizione fissa da confrontare con x_focus:
    // il vincolo di qualità focale è delegato interamente a apply_dof_filter (dof_min).
    if (fc.mobile_focus) return configs;
    std::vector<ConfigData> result;
    for (const auto& c : configs)
        if (std::abs(c.x_focus - fc.x_det) <= fc.focus_tol)
            result.push_back(c);
    return result;
}

// Popola cd.x_lo/x_hi: usa i valori forniti (colonne opzionali del TSV dof_map) se
// presenti, altrimenti fallback simmetrico attorno a x_focus (TSV generati prima
// dell'introduzione delle colonne x_lo/x_hi). Richiede cd.x_focus/cd.DoF già impostati.
inline void resolve_dof_bounds(ConfigData& cd, bool has_x_lo_hi, double x_lo_raw, double x_hi_raw) {
    if (has_x_lo_hi) {
        cd.x_lo = x_lo_raw;
        cd.x_hi = x_hi_raw;
    } else {
        cd.x_lo = cd.x_focus - cd.DoF / 2.0;
        cd.x_hi = cd.x_focus + cd.DoF / 2.0;
    }
}

// Applica DoF >= dof_min; se dof_min <= 0 ritorna tutti
inline std::vector<ConfigData> apply_dof_filter(const std::vector<ConfigData>& configs,
                                                  const FilterConfig& fc) {
    if (fc.dof_min <= 0.0) return configs;
    std::vector<ConfigData> result;
    for (const auto& c : configs)
        if (c.DoF >= fc.dof_min)
            result.push_back(c);
    return result;
}

// Applica EE80 <= ee80_max; se ee80_max <= 0 o EE80 non disponibile ritorna tutti
inline std::vector<ConfigData> apply_ee80_filter(const std::vector<ConfigData>& configs,
                                                   const FilterConfig& fc) {
    if (fc.ee80_max <= 0.0) return configs;
    std::vector<ConfigData> result;
    for (const auto& c : configs)
        if (std::isfinite(c.EE80) && c.EE80 <= fc.ee80_max)
            result.push_back(c);
    return result;
}

// Normalizzatori indipendenti dai pesi (max di eta/Q/DoF/M_abs_err sul set).
// Vanno calcolati sull'intero set filtrato (non solo sul fronte) per restare
// equivalenti a compute_mtot: cambiare il set su cui si calcolano i massimi
// cambierebbe la normalizzazione e quindi i valori di Mtot.
struct MtotNormalizers {
    double eta_max   = 1.0;
    double Q_max     = 1.0;
    double DoF_max   = 1.0;
    double M_err_max = 1.0;
};

inline MtotNormalizers compute_mtot_normalizers(const std::vector<ConfigData>& configs) {
    MtotNormalizers n{0.0, 0.0, 0.0, 0.0};
    for (const auto& c : configs) {
        n.eta_max   = std::max(n.eta_max,   c.eta);
        n.Q_max     = std::max(n.Q_max,     c.Q);
        n.DoF_max   = std::max(n.DoF_max,   c.DoF);
        n.M_err_max = std::max(n.M_err_max, c.M_abs_err);
    }
    if (n.eta_max   <= 0.0) n.eta_max   = 1.0;
    if (n.Q_max     <= 0.0) n.Q_max     = 1.0;
    if (n.DoF_max   <= 0.0) n.DoF_max   = 1.0;
    if (n.M_err_max <= 0.0) n.M_err_max = 1.0;
    return n;
}

// Calcola Mtot per tutti i config (modifica in-place) dati normalizzatori precalcolati.
// Formula: Mtot = w_eta*(η/η_max) + w_Q*(1−Q/Q_max) + w_dof*(DoF/DoF_max) + w_M*(1−M_abs_err/M_abs_err_max)
// Tutti i termini ∈ [0,1]; Mtot ∈ [0,1] dove 1 è il punto ideale.
inline void apply_mtot(std::vector<ConfigData>& configs, const MtotNormalizers& n,
                       const WeightConfig& wc) {
    double w_sum = wc.w_eta + wc.w_Q + wc.w_dof + wc.w_M;
    if (w_sum <= 0.0) w_sum = 1.0;
    const double nw_eta = wc.w_eta / w_sum;
    const double nw_Q   = wc.w_Q   / w_sum;
    const double nw_dof = wc.w_dof / w_sum;
    const double nw_M   = wc.w_M   / w_sum;

    for (auto& c : configs) {
        double eta_n = c.eta       / n.eta_max;
        double Q_n   = c.Q         / n.Q_max;
        double dof_n = c.DoF       / n.DoF_max;
        double M_n   = c.M_abs_err / n.M_err_max;
        c.Mtot = nw_eta * eta_n + nw_Q * (1.0 - Q_n) + nw_dof * dof_n + nw_M * (1.0 - M_n);
    }
}

// Calcola Mtot con normalizzazione interna al vettore (equivalente a
// apply_mtot(configs, compute_mtot_normalizers(configs), wc)).
inline void compute_mtot(std::vector<ConfigData>& configs, const WeightConfig& wc) {
    if (configs.empty()) return;
    apply_mtot(configs, compute_mtot_normalizers(configs), wc);
}

// Marca on_pareto=true per i punti non dominati. Dipende solo da eta/Q/M_abs_err,
// non dai pesi: può essere calcolato una sola volta e riusato per qualunque sweep di pesi.
inline void mark_pareto_front(std::vector<ConfigData>& configs) {
    const size_t n = configs.size();
    for (size_t i = 0; i < n; ++i) {
        configs[i].on_pareto = true;
        for (size_t j = 0; j < n; ++j) {
            if (i != j && pareto_dominates(configs[j], configs[i])) {
                configs[i].on_pareto = false;
                break;
            }
        }
    }
}

// Assegna pareto_rank in ordine decrescente di Mtot (1=best) ai punti on_pareto.
// Richiede che Mtot sia già stato calcolato (compute_mtot/apply_mtot) e il fronte
// già marcato (mark_pareto_front).
inline void assign_pareto_ranks(std::vector<ConfigData>& configs) {
    std::vector<size_t> front_idx;
    front_idx.reserve(configs.size());
    for (size_t i = 0; i < configs.size(); ++i)
        if (configs[i].on_pareto) front_idx.push_back(i);

    std::sort(front_idx.begin(), front_idx.end(),
              [&](size_t a, size_t b) { return configs[a].Mtot > configs[b].Mtot; });

    for (size_t r = 0; r < front_idx.size(); ++r)
        configs[front_idx[r]].pareto_rank = static_cast<int>(r + 1);
}

// Marca on_pareto=true per i punti non dominati; poi assegna pareto_rank
// in ordine decrescente di Mtot (1=best) ai punti sul fronte.
inline void compute_pareto_front(std::vector<ConfigData>& configs, const WeightConfig& wc) {
    compute_mtot(configs, wc);
    mark_pareto_front(configs);
    assign_pareto_ranks(configs);
}

// ── Weight sweep (spazio dei pesi, sul fronte già calcolato) ────────────────

struct SweepPoint {
    WeightConfig wc;
    double x1_winner   = 0.0;
    double x2_winner   = 0.0;
    double Mtot_winner = 0.0;
    int    category_id = -1;
};

// Griglia baricentrica sui 4 pesi (w_eta,w_Q,w_dof,w_M), ogni vettore somma
// esattamente a 1 (aritmetica intera: n = round(1/step)).
// Per step=0.05 (n=20): C(23,3) = 1771 punti.
inline std::vector<WeightConfig> generate_barycentric_grid(double step) {
    std::vector<WeightConfig> grid;
    const int n = static_cast<int>(std::lround(1.0 / step));
    if (n <= 0) return grid;

    for (int i = 0; i <= n; ++i) {
        for (int j = 0; i + j <= n; ++j) {
            for (int k = 0; i + j + k <= n; ++k) {
                const int l = n - i - j - k;
                WeightConfig wc;
                wc.w_eta = static_cast<double>(i) / n;
                wc.w_Q   = static_cast<double>(j) / n;
                wc.w_dof = static_cast<double>(k) / n;
                wc.w_M   = static_cast<double>(l) / n;
                grid.push_back(wc);
            }
        }
    }
    return grid;
}

// Per ogni vettore di pesi della griglia, ricalcola Mtot sul fronte (già filtrato
// e marcato) e trova il vincitore. Economico: O(grid.size() * front.size()).
inline std::vector<SweepPoint> run_weight_sweep(std::vector<ConfigData> front,
                                                 const MtotNormalizers& normalizers,
                                                 const std::vector<WeightConfig>& grid) {
    std::vector<SweepPoint> results;
    results.reserve(grid.size());
    for (const auto& wc : grid) {
        apply_mtot(front, normalizers, wc);
        auto best = std::max_element(front.begin(), front.end(),
            [](const ConfigData& a, const ConfigData& b) { return a.Mtot < b.Mtot; });

        SweepPoint sp;
        sp.wc = wc;
        if (best != front.end()) {
            sp.x1_winner   = best->x1;
            sp.x2_winner   = best->x2;
            sp.Mtot_winner = best->Mtot;
        }
        results.push_back(sp);
    }
    return results;
}

// Deduplica i punti (x1,x2) vincenti e assegna ID sequenziali stabili per la colorazione.
inline void assign_category_ids(std::vector<SweepPoint>& sweep) {
    std::unordered_map<int64_t, int> id_of_key;
    for (auto& sp : sweep) {
        const int64_t key = coord_key(sp.x1_winner, sp.x2_winner);
        auto [it, inserted] = id_of_key.try_emplace(key, static_cast<int>(id_of_key.size()));
        sp.category_id = it->second;
    }
}

// Trasformazione baricentrica (a=w_eta, b=w_Q, c=w_dof; a+b+c=1) -> cartesiana,
// per il plot ternario. Vertici: A=(0,0), B=(1,0), C=(0.5, sqrt(3)/2).
inline std::pair<double, double> ternary_to_xy(double a, double b, double c) {
    (void)a;  // ridondante (a = 1-b-c) ma incluso per chiarezza dei chiamanti
    return {b + 0.5 * c, (std::sqrt(3.0) / 2.0) * c};
}

// ── Sviluppo (net) del tetraedro dei pesi ────────────────────────────────────
//
// Il tetraedro dei pesi (w_eta,w_Q,w_dof,w_M) ha 4 facce triangolari. La faccia
// w_M=0 è il triangolo centrale A=(0,0)/B=(1,0)/C=(0.5,sqrt(3)/2) di ternary_to_xy.
// Le altre 3 facce (w_eta=0, w_Q=0, w_dof=0) condividono ciascuna un lato con la
// faccia centrale (rispettivamente BC, AC, AB) e si "aprono" ribaltando il
// triangolo centrale attorno a quel lato condiviso — tecnica standard di
// sviluppo del tetraedro regolare per diagrammi di composizione quaternari.

// Riflette il punto P rispetto alla retta passante per Q1,Q2 (proiezione + specchio).
inline std::pair<double, double> reflect_across_line(std::pair<double, double> P,
                                                       std::pair<double, double> Q1,
                                                       std::pair<double, double> Q2) {
    const double dx = Q2.first  - Q1.first;
    const double dy = Q2.second - Q1.second;
    const double len2 = dx * dx + dy * dy;
    if (len2 <= 0.0) return P;  // Q1==Q2, retta degenere
    const double vx = P.first  - Q1.first;
    const double vy = P.second - Q1.second;
    const double t = (vx * dx + vy * dy) / len2;
    const double fx = Q1.first  + t * dx;
    const double fy = Q1.second + t * dy;
    return {2.0 * fx - P.first, 2.0 * fy - P.second};
}

// Trasformazione affine (x,y) -> (a*x+b*y+e, c*x+d*y+f) che manda il triangolo
// canonico locale L0=(0,0), L1=(1,0), L2=(0.5,sqrt(3)/2) (stessi vertici usati da
// ternary_to_xy) sui 3 vertici globali noti G0,G1,G2, in ordine corrispondente.
// Forma chiusa (non iterativa): risolta direttamente dalle 3 corrispondenze punto-punto.
struct FaceAffine {
    double a = 1.0, b = 0.0, c = 0.0, d = 1.0, e = 0.0, f = 0.0;

    std::pair<double, double> apply(double x, double y) const {
        return {a * x + b * y + e, c * x + d * y + f};
    }
};

inline FaceAffine solve_face_affine(std::pair<double, double> G0,
                                     std::pair<double, double> G1,
                                     std::pair<double, double> G2) {
    FaceAffine fa;
    fa.e = G0.first;
    fa.f = G0.second;
    fa.a = G1.first  - fa.e;
    fa.c = G1.second - fa.f;
    const double h = std::sqrt(3.0) / 2.0;
    fa.b = (G2.first  - fa.e - 0.5 * fa.a) / h;
    fa.d = (G2.second - fa.f - 0.5 * fa.c) / h;
    return fa;
}

}  // namespace riptide::pareto
