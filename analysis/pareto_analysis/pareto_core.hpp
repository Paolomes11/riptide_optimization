#pragma once

#include <algorithm>
#include <cmath>
#include <string>
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
    bool   on_pareto   = false;
    double Mtot        = 0.0;
    int    pareto_rank = 0;  // 1=best Mtot sul fronte, 0=non sul fronte
};

struct FilterConfig {
    double eta_frac  = 0.75;   // soglia η relativa a η_max
    double x_det     = 180.0;  // posizione detector [mm]
    double focus_tol = 15.0;   // tolleranza fuoco [mm]
    double dof_min   = 0.0;    // DoF minima [mm] (0 = disabilitato)
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

// A domina B nel senso di Pareto (η↑, Q↓):
// A.eta >= B.eta AND A.Q <= B.Q con almeno una disuguaglianza stretta
inline bool pareto_dominates(const ConfigData& A, const ConfigData& B) {
    return (A.eta >= B.eta && A.Q <= B.Q) && (A.eta > B.eta || A.Q < B.Q);
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
    std::vector<ConfigData> result;
    for (const auto& c : configs)
        if (std::abs(c.x_focus - fc.x_det) <= fc.focus_tol)
            result.push_back(c);
    return result;
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

// Calcola Mtot per tutti i config (modifica in-place).
// Formula: Mtot = w_eta*(η/η_max) - w_Q*(Q/Q_max) + w_dof*(DoF/DoF_max) - w_M*(M_abs_err/M_abs_err_max)
// Normalizzazione interna al vettore; se max==0 usa 1.0 come guardia.
inline void compute_mtot(std::vector<ConfigData>& configs, const WeightConfig& wc) {
    if (configs.empty()) return;

    double eta_max     = 0.0;
    double Q_max       = 0.0;
    double DoF_max     = 0.0;
    double M_err_max   = 0.0;
    for (const auto& c : configs) {
        eta_max   = std::max(eta_max,   c.eta);
        Q_max     = std::max(Q_max,     c.Q);
        DoF_max   = std::max(DoF_max,   c.DoF);
        M_err_max = std::max(M_err_max, c.M_abs_err);
    }
    if (eta_max   <= 0.0) eta_max   = 1.0;
    if (Q_max     <= 0.0) Q_max     = 1.0;
    if (DoF_max   <= 0.0) DoF_max   = 1.0;
    if (M_err_max <= 0.0) M_err_max = 1.0;

    for (auto& c : configs) {
        double eta_n = c.eta       / eta_max;
        double Q_n   = c.Q         / Q_max;
        double dof_n = c.DoF       / DoF_max;
        double M_n   = c.M_abs_err / M_err_max;
        c.Mtot = wc.w_eta * eta_n - wc.w_Q * Q_n + wc.w_dof * dof_n - wc.w_M * M_n;
    }
}

// Marca on_pareto=true per i punti non dominati; poi assegna pareto_rank
// in ordine decrescente di Mtot (1=best) ai punti sul fronte.
// Assume che compute_mtot sia già stato chiamato.
inline void compute_pareto_front(std::vector<ConfigData>& configs) {
    const size_t n = configs.size();

    // Determina il fronte: un punto è sul fronte se nessun altro lo domina
    for (size_t i = 0; i < n; ++i) {
        configs[i].on_pareto = true;
        for (size_t j = 0; j < n; ++j) {
            if (i != j && pareto_dominates(configs[j], configs[i])) {
                configs[i].on_pareto = false;
                break;
            }
        }
    }

    // Raccoglie indici del fronte e li ordina per Mtot decrescente
    std::vector<size_t> front_idx;
    front_idx.reserve(n);
    for (size_t i = 0; i < n; ++i)
        if (configs[i].on_pareto) front_idx.push_back(i);

    std::sort(front_idx.begin(), front_idx.end(),
              [&](size_t a, size_t b) { return configs[a].Mtot > configs[b].Mtot; });

    for (size_t r = 0; r < front_idx.size(); ++r)
        configs[front_idx[r]].pareto_rank = static_cast<int>(r + 1);
}

}  // namespace riptide::pareto
