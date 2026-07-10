/*
 * psf_dof_extractor — legge PsfDofConfigs + PsfDofRuns da psf_dof.root,
 * calcola DoF / EE80 per config_id e scrive:
 *   - <output_data.root>  : TTree "PsfDofSummary" compatto
 *   - <output.tsv>        : stesso formato di resolution_map.tsv
 *
 * CLI: ./psf_dof_extractor <input.root> <output_data.root> <output.tsv> <config.json>
 *
 * Formule identiche a resolution_map.cpp.
 */

#include <nlohmann/json.hpp>

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// ─── Strutture dati ──────────────────────────────────────────────────────────

struct ConfigInfo {
  int    config_id = -1;
  double x1        = 0.0;
  double x2        = 0.0;
  double x_virtual = 0.0;
};

struct RunRow {
  int    config_id = -1;
  int    run_id    = -1;
  double x_src     = 0.0;
  double y_src     = 0.0;
  double n_hits    = 0.0;
  double mu_y      = 0.0;
  double sigma_z   = 0.0;
  double sigma_dz  = 0.0;
  double cov_z_dz  = 0.0;
  double mu_dy     = 0.0;
  double sigma_y   = 0.0;
  double sigma_dy  = 0.0;
  double cov_y_dy  = 0.0;
};

struct Aggregated {
  double sum_dof     = 0.0; int n_dof     = 0;
  double sum_EE80    = 0.0; int n_EE80    = 0;
};

// ─── Formule (da resolution_map.cpp) ─────────────────────────────────────────

static std::vector<double> build_scan(double scan_min, double scan_max, double scan_step) {
  std::vector<double> x;
  for (double v = scan_min; v <= scan_max + 1e-12; v += scan_step)
    x.push_back(v);
  return x;
}

static std::pair<double, double> sigma_z_min_and_focus(const RunRow& r, const ConfigInfo& c,
                                                        const std::vector<double>& x_scan) {
  double best_sigma = std::numeric_limits<double>::infinity();
  double best_x     = std::numeric_limits<double>::quiet_NaN();
  double s0_2 = r.sigma_z * r.sigma_z;
  double sdz2 = r.sigma_dz * r.sigma_dz;
  double cov  = r.cov_z_dz;
  double x0   = c.x_virtual;
  for (double x_det : x_scan) {
    double dx    = x_det - x0;
    double v     = s0_2 + 2.0 * cov * dx + sdz2 * dx * dx;
    double sigma = std::sqrt(std::max(0.0, v));
    if (!std::isfinite(sigma)) continue;
    if (sigma < best_sigma) { best_sigma = sigma; best_x = x_det; }
  }
  return {best_sigma, best_x};
}

static double compute_dof(const RunRow& r, const ConfigInfo& c,
                           const std::vector<double>& x_scan, double k) {
  if (x_scan.size() < 2) return 0.0;
  auto [sigma_min, x_focus] = sigma_z_min_and_focus(r, c, x_scan);
  if (!std::isfinite(sigma_min) || !std::isfinite(x_focus)) return 0.0;

  double thr  = k * sigma_min;
  double s0_2 = r.sigma_z * r.sigma_z;
  double sdz2 = r.sigma_dz * r.sigma_dz;
  double cov  = r.cov_z_dz;
  double x0   = c.x_virtual;

  std::vector<double> sigma(x_scan.size());
  for (size_t i = 0; i < x_scan.size(); ++i) {
    double dx = x_scan[i] - x0;
    sigma[i]  = std::sqrt(std::max(0.0, s0_2 + 2.0*cov*dx + sdz2*dx*dx));
  }

  size_t i_min = std::distance(sigma.begin(),
                               std::min_element(sigma.begin(), sigma.end()));
  int i_lo = static_cast<int>(i_min);
  int i_hi = static_cast<int>(i_min);
  while (i_lo - 1 >= 0 && sigma[static_cast<size_t>(i_lo-1)] < thr) --i_lo;
  while (i_hi + 1 < static_cast<int>(sigma.size()) && sigma[static_cast<size_t>(i_hi+1)] < thr) ++i_hi;
  return x_scan[static_cast<size_t>(i_hi)] - x_scan[static_cast<size_t>(i_lo)];
}

static std::optional<double> compute_EE80(const RunRow& r, const ConfigInfo& c,
                                           const std::vector<double>& x_scan) {
  auto [sigma_z_min, x_focus] = sigma_z_min_and_focus(r, c, x_scan);
  if (!std::isfinite(sigma_z_min) || !std::isfinite(x_focus)) return std::nullopt;
  double dx           = x_focus - c.x_virtual;
  double v_y          = r.sigma_y*r.sigma_y + 2.0*r.cov_y_dy*dx + r.sigma_dy*r.sigma_dy*dx*dx;
  double sigma_y_foc  = std::sqrt(std::max(0.0, v_y));
  double sigma_rms    = std::sqrt((sigma_y_foc*sigma_y_foc + sigma_z_min*sigma_z_min) / 2.0);
  constexpr double k_EE80 = 2.0 * 1.7941;
  return k_EE80 * sigma_rms;
}

// ─── main ────────────────────────────────────────────────────────────────────

int main(int argc, char** argv) {
  if (argc < 5) {
    std::cerr << "Uso: psf_dof_extractor <input.root> <output_data.root>"
              << " <output.tsv> <config.json>\n";
    return 1;
  }
  const std::string in_path  = argv[1];
  const std::string out_root = argv[2];
  const std::string out_tsv  = argv[3];
  const std::string cfg_path = argv[4];

  // Leggi config.json
  std::ifstream jf(cfg_path);
  if (!jf.is_open()) {
    std::cerr << "Errore: impossibile aprire " << cfg_path << "\n"; return 1;
  }
  nlohmann::json config;
  jf >> config;

  double scan_min     = config.value("dof_x_scan_min",  180.0);
  double scan_max     = config.value("dof_x_scan_max",  400.0);
  double scan_step    = config.value("dof_x_scan_step",   0.5);
  double k            = config.value("dof_k_threshold", std::sqrt(2.0));
  double lens_det_gap = config.value("lens_det_gap",      0.0);

  if (!(scan_step > 0.0) || !(scan_max > scan_min) || !(k > 0.0)) {
    std::cerr << "Errore: parametri scan/k non validi nella config\n"; return 1;
  }

  // Apri file di input
  TFile* fin = TFile::Open(in_path.c_str(), "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Errore: impossibile aprire " << in_path << "\n"; return 1;
  }
  TTree* tree_cfg  = dynamic_cast<TTree*>(fin->Get("PsfDofConfigs"));
  TTree* tree_runs = dynamic_cast<TTree*>(fin->Get("PsfDofRuns"));
  if (!tree_cfg || !tree_runs) {
    std::cerr << "Errore: TTree PsfDofConfigs o PsfDofRuns non trovato\n"; return 1;
  }

  // Leggi configurazioni
  int    cid = 0; double cx1 = 0, cx2 = 0, cx_virt = 0;
  tree_cfg->SetBranchAddress("config_id", &cid);
  tree_cfg->SetBranchAddress("x1",        &cx1);
  tree_cfg->SetBranchAddress("x2",        &cx2);
  tree_cfg->SetBranchAddress("x_virtual", &cx_virt);

  std::unordered_map<int, ConfigInfo> config_map;
  config_map.reserve(static_cast<size_t>(tree_cfg->GetEntries()));
  for (Long64_t i = 0; i < tree_cfg->GetEntries(); ++i) {
    tree_cfg->GetEntry(i);
    config_map[cid] = {cid, cx1, cx2, cx_virt};
  }
  std::cout << "Configurazioni: " << config_map.size() << "\n";

  // Costruisci x_scan per config (accounting per lens_det_gap)
  std::unordered_map<int, std::vector<double>> x_scan_map;
  x_scan_map.reserve(config_map.size());
  for (const auto& [id, c] : config_map) {
    double cfg_scan_min = std::max(scan_min, c.x2 + lens_det_gap);
    auto xs = build_scan(cfg_scan_min, scan_max, scan_step);
    if (!xs.empty()) x_scan_map[id] = std::move(xs);
  }

  // Leggi run e calcola metriche
  RunRow rr{};
  tree_runs->SetBranchAddress("config_id", &rr.config_id);
  tree_runs->SetBranchAddress("run_id",    &rr.run_id);
  tree_runs->SetBranchAddress("x_source",  &rr.x_src);
  tree_runs->SetBranchAddress("y_source",  &rr.y_src);
  tree_runs->SetBranchAddress("n_hits",    &rr.n_hits);
  tree_runs->SetBranchAddress("mu_y",      &rr.mu_y);
  tree_runs->SetBranchAddress("sigma_z",   &rr.sigma_z);
  tree_runs->SetBranchAddress("sigma_dz",  &rr.sigma_dz);
  tree_runs->SetBranchAddress("cov_z_dz",  &rr.cov_z_dz);
  tree_runs->SetBranchAddress("mu_dy",     &rr.mu_dy);
  tree_runs->SetBranchAddress("sigma_y",   &rr.sigma_y);
  tree_runs->SetBranchAddress("sigma_dy",  &rr.sigma_dy);
  tree_runs->SetBranchAddress("cov_y_dy",  &rr.cov_y_dy);

  std::unordered_map<int, Aggregated> agg;
  agg.reserve(config_map.size());

  Long64_t n_runs = tree_runs->GetEntries();
  for (Long64_t i = 0; i < n_runs; ++i) {
    tree_runs->GetEntry(i);
    if (!(rr.n_hits > 0.0)) continue;
    auto it_c = config_map.find(rr.config_id);
    if (it_c == config_map.end()) continue;
    auto it_s = x_scan_map.find(rr.config_id);
    if (it_s == x_scan_map.end()) continue;

    const ConfigInfo&         c  = it_c->second;
    const std::vector<double>& xs = it_s->second;

    double dof_val = compute_dof(rr, c, xs, k);
    if (std::isfinite(dof_val) && dof_val > 0.0) {
      agg[rr.config_id].sum_dof += dof_val;
      agg[rr.config_id].n_dof  += 1;
    }
    auto ee80 = compute_EE80(rr, c, xs);
    if (ee80.has_value() && std::isfinite(*ee80) && *ee80 > 0.0) {
      agg[rr.config_id].sum_EE80 += *ee80;
      agg[rr.config_id].n_EE80  += 1;
    }

    if ((i + 1) % 50000 == 0)
      std::cout << "  Processati " << i+1 << " / " << n_runs << " run...\n";
  }
  fin->Close();

  // Scrivi ROOT di output
  TFile* fout = TFile::Open(out_root.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "Errore: impossibile creare " << out_root << "\n"; return 1;
  }
  TTree* tSumm = new TTree("PsfDofSummary", "Aggregated DoF/EE80 per config_id");
  int    o_config_id;
  double o_x1, o_x2, o_dof_mean, o_EE80_mean;
  int    o_n_dof, o_n_EE80;
  tSumm->Branch("config_id",        &o_config_id,    "config_id/I");
  tSumm->Branch("x1",               &o_x1,            "x1/D");
  tSumm->Branch("x2",               &o_x2,            "x2/D");
  tSumm->Branch("dof_mean",         &o_dof_mean,      "dof_mean/D");
  tSumm->Branch("EE80_mean",        &o_EE80_mean,     "EE80_mean/D");
  tSumm->Branch("n_runs_dof",       &o_n_dof,         "n_runs_dof/I");
  tSumm->Branch("n_runs_EE80",      &o_n_EE80,        "n_runs_EE80/I");

  // Scrivi TSV di output
  std::ofstream tsv(out_tsv);
  if (!tsv.is_open()) {
    std::cerr << "Errore: impossibile creare " << out_tsv << "\n"; return 1;
  }
  tsv << "x1\tx2\tdof_mean\tEE80_mean\t"
         "config_id\tn_runs_dof\tn_runs_EE80\n";

  constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

  for (const auto& [id, c] : config_map) {
    auto it = agg.find(id);

    double dof_mean  = (it != agg.end() && it->second.n_dof > 0)
                       ? it->second.sum_dof / it->second.n_dof : 0.0;
    double ee80_mean = (it != agg.end() && it->second.n_EE80 > 0)
                       ? it->second.sum_EE80 / it->second.n_EE80 : kNaN;
    int n_dof_out   = (it != agg.end()) ? it->second.n_dof     : 0;
    int n_EE80_out  = (it != agg.end()) ? it->second.n_EE80    : 0;

    o_config_id    = id;
    o_x1           = c.x1;
    o_x2           = c.x2;
    o_dof_mean     = dof_mean;
    o_EE80_mean    = ee80_mean;
    o_n_dof        = n_dof_out;
    o_n_EE80       = n_EE80_out;
    tSumm->Fill();

    auto fmt = [](double v) -> std::string {
      if (!std::isfinite(v)) return "nan";
      std::ostringstream oss;
      oss << v;
      return oss.str();
    };
    tsv << c.x1 << "\t" << c.x2 << "\t"
        << fmt(dof_mean) << "\t"
        << fmt(ee80_mean) << "\t"
        << id << "\t"
        << n_dof_out << "\t" << n_EE80_out << "\n";
  }

  fout->cd();
  tSumm->Write();
  fout->Close();

  std::cout << "Configurazioni aggregate: " << agg.size()
            << " / " << config_map.size() << "\n";
  std::cout << "ROOT: " << out_root << "\n";
  std::cout << "TSV:  " << out_tsv  << "\n";
  return 0;
}
