#pragma once
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

struct PairStats {
  std::string lens1_id, lens2_id;
  double mean     = 0.0;
  double variance = 0.0;
  double score    = 0.0;
  size_t n_points = 0;
};

// CfgT must expose: int focus_valid, std::string lens1_id, std::string lens2_id.
// always_divide=true  → eff = sum / n_entries[id]  (lens-sim mode o plot3D)
// always_divide=false → eff = sum                   (opt-mode, single entry per config)
template <typename CfgT>
std::vector<PairStats> compute_pair_ranking(const std::map<int, double>& eff_sum,
                                             const std::map<int, int>&    n_entries,
                                             const std::map<int, CfgT>&   config_map,
                                             bool                          always_divide) {
  std::map<std::pair<std::string, std::string>, std::vector<double>> pair_effs;
  for (auto const& [id, sum] : eff_sum) {
    auto cfg_it = config_map.find(id);
    if (cfg_it == config_map.end()) continue;
    auto const& cfg = cfg_it->second;
    if (cfg.focus_valid == 0) continue;
    double eff;
    if (always_divide) {
      auto ne_it = n_entries.find(id);
      if (ne_it == n_entries.end() || ne_it->second == 0) continue;
      eff = sum / static_cast<double>(ne_it->second);
    } else {
      eff = sum;
    }
    pair_effs[{cfg.lens1_id, cfg.lens2_id}].push_back(eff);
  }

  std::vector<PairStats> stats;
  stats.reserve(pair_effs.size());
  for (auto const& [pair, ev] : pair_effs) {
    PairStats ps;
    ps.lens1_id = pair.first;
    ps.lens2_id = pair.second;
    ps.n_points = ev.size();
    double s    = 0.0;
    for (double e : ev) s += e;
    ps.mean     = s / static_cast<double>(ev.size());
    double vs   = 0.0;
    for (double e : ev) vs += (e - ps.mean) * (e - ps.mean);
    ps.variance = vs / static_cast<double>(ev.size());
    // Lower Confidence Bound: penalizza le coppie con std_eta alta rispetto
    // alla media (Markowitz mean-variance / bandit LCB, k=0.5 calibrato sulla
    // distribuzione empirica del CV). Vedi README §Riferimenti metodologici.
    constexpr double kRiskAversion = 0.5;
    ps.score = ps.mean - kRiskAversion * std::sqrt(ps.variance);
    stats.push_back(ps);
  }
  std::sort(stats.begin(), stats.end(),
            [](const PairStats& a, const PairStats& b) { return a.score > b.score; });
  return stats;
}

inline void write_ranking_csv(const std::vector<PairStats>& stats,
                               const std::string&             csv_path) {
  std::ofstream out(csv_path);
  out << "rank,lens1_id,lens2_id,score,mean_eta,std_eta,variance,n_grid_points\n";
  out << std::fixed << std::setprecision(6);
  for (size_t r = 0; r < stats.size(); ++r) {
    auto const& ps = stats[r];
    out << (r + 1) << "," << ps.lens1_id << "," << ps.lens2_id << "," << ps.score << ","
        << ps.mean << "," << std::sqrt(ps.variance) << "," << ps.variance << ","
        << ps.n_points << "\n";
  }
}

inline void print_ranking_summary(const std::vector<PairStats>& stats, size_t top = 5) {
  top = std::min(stats.size(), top);
  std::cout << "\nRanking coppie di lenti (top " << top << "):\n";
  std::cout << std::left << std::setw(6) << "rank" << std::setw(22) << "lens1"
            << std::setw(22) << "lens2" << std::right << std::setw(12) << "score"
            << std::setw(12) << "mean_eta" << std::setw(12) << "std_eta" << std::setw(8)
            << "N\n";
  for (size_t r = 0; r < top; ++r) {
    auto const& ps = stats[r];
    std::cout << std::left << std::setw(6) << (r + 1) << std::setw(22) << ps.lens1_id
              << std::setw(22) << ps.lens2_id << std::right << std::fixed
              << std::setprecision(6) << std::setw(12) << ps.score << std::setw(12)
              << ps.mean << std::setw(12) << std::sqrt(ps.variance) << std::setw(8)
              << ps.n_points << "\n";
  }
}
