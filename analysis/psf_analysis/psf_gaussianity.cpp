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

/*
 * psf_gaussianity — diagnostica di gaussianita' della PSF sul detector.
 *
 * Per ogni run in lens.root (Configurations+Runs), stima media/covarianza
 * locale con compute_psf (filtro Mahalanobis iterativo, psf_common.hpp),
 * poi calcola la distanza di Mahalanobis quadrata D^2 di ogni hit grezzo
 * rispetto a quella stima. Sotto H0 (gaussiana bivariata) D^2 ~ chi^2_2
 * indipendentemente dalla covarianza specifica del run, quindi i valori
 * D^2 di tutti i run possono essere raggruppati in un unico campione
 * (pooling) per costruire un Q-Q plot contro chi^2_2 e stimare la curtosi
 * in eccesso multivariata di Mardia: gamma2 = mean(D^4)/8 - 1.
 *
 * Uso:
 *   ./psf_gaussianity [--input lens.root] [--output out.png]
 *                     [--sigma-cut 3.0] [--n-iter 4] [--min-hits 30]
 *                     [--contamination-alpha 0.05] [--hz-subsample-size 2000]
 *
 * Il test di Mardia calcolato su tutti gli hit grezzi (pool completo) e' noto
 * essere estremamente sensibile a qualunque contaminazione, per quanto piccola
 * (luce spuria non fisica, es. riflessione totale interna). Per distinguere
 * "nucleo non gaussiano" da "nucleo gaussiano + coda di contaminazione" si
 * classifica ogni hit con un cutoff Bonferroni-corretto sulla coda esatta di
 * chi^2_2 (bonferroni_chi2_2_cutoff, psf_common.hpp): D^2 > cutoff => flaggato
 * come contaminazione a livello di confidenza family-wise 1-contamination_alpha.
 * Mardia gamma2 e' poi ricalcolato solo sul nucleo (D^2 <= cutoff). In aggiunta,
 * il test di Henze-Zirkler (piu' potente ma O(n^2)) viene applicato come
 * diagnostica secondaria su un sottocampione casuale del run rappresentativo.
 */

#include "psf_common.hpp"

#include <spdlog/spdlog.h>

#include <TCanvas.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

// Densita' chi^2 con 2 gradi di liberta' (forma chiusa: f(x) = 0.5*exp(-x/2))
static double chi2_2_pdf(double x) {
  return 0.5 * std::exp(-0.5 * x);
}

// Curtosi in eccesso multivariata di Mardia (p=2) su un pool di D^2: b2 = mean(D^4),
// gamma2 = b2/8 - 1, se(gamma2) = sqrt(64/N)/8. Atteso 0 per gaussiana bivariata.
struct MardiaStats {
  double b2, gamma2, se_gamma2;
};

static MardiaStats mardia_stats(const std::vector<double>& d2) {
  double sum_d4 = 0.0;
  for (double d2v : d2)
    sum_d4 += d2v * d2v;
  double n         = static_cast<double>(d2.size());
  MardiaStats stats;
  stats.b2        = sum_d4 / n;
  stats.gamma2    = stats.b2 / 8.0 - 1.0;
  stats.se_gamma2 = std::sqrt(64.0 / n) / 8.0;
  return stats;
}

// Statistica di Henze-Zirkler (1990) per normalita' bivariata (p=2), con smoothing
// beta = (1/sqrt(2))*[n(2p+1)/4]^(1/(p+4)). Media/covarianza stimate internamente
// sul campione fornito (plug-in, non robusto: il campione va gia' filtrato a monte).
// Restituisce {HZ, p-value} col null asintotico lognormale.
struct HZResult {
  double hz, p_value;
  bool ok;
};

static HZResult henze_zirkler(const std::vector<Hit>& hits) {
  HZResult res{0.0, 1.0, false};
  int n = static_cast<int>(hits.size());
  if (n < 10)
    return res;

  double mean_y = 0.0, mean_z = 0.0;
  for (const auto& h : hits) {
    mean_y += h.y;
    mean_z += h.z;
  }
  mean_y /= n;
  mean_z /= n;

  double var_y = 0.0, var_z = 0.0, cov_yz = 0.0;
  for (const auto& h : hits) {
    double dy = h.y - mean_y, dz = h.z - mean_z;
    var_y += dy * dy;
    var_z += dz * dz;
    cov_yz += dy * dz;
  }
  var_y /= (n - 1);
  var_z /= (n - 1);
  cov_yz /= (n - 1);

  constexpr double kRegEps = 1e-4;
  double reg_yy = var_y + kRegEps, reg_zz = var_z + kRegEps;
  double det    = reg_yy * reg_zz - cov_yz * cov_yz;
  if (det <= 0.0)
    return res;
  double inv_yy = reg_zz / det, inv_zz = reg_yy / det, inv_yz = -cov_yz / det;

  auto d2_to_mean = [&](const Hit& h) {
    double dy = h.y - mean_y, dz = h.z - mean_z;
    return dy * dy * inv_yy + dz * dz * inv_zz + 2.0 * dy * dz * inv_yz;
  };
  auto d2_pair = [&](const Hit& a, const Hit& b) {
    double dy = a.y - b.y, dz = a.z - b.z;
    return dy * dy * inv_yy + dz * dz * inv_zz + 2.0 * dy * dz * inv_yz;
  };

  constexpr double p = 2.0;
  double beta = (1.0 / std::sqrt(2.0)) *
                std::pow(static_cast<double>(n) * (2.0 * p + 1.0) / 4.0, 1.0 / (p + 4.0));
  double b2 = beta * beta;

  std::vector<double> di(n);
  for (int i = 0; i < n; ++i)
    di[i] = d2_to_mean(hits[i]);

  double term_pairs = 0.0;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      term_pairs += std::exp(-0.5 * b2 * d2_pair(hits[i], hits[j]));
  term_pairs /= (static_cast<double>(n) * n);

  double term_mean = 0.0;
  for (int i = 0; i < n; ++i)
    term_mean += std::exp(-b2 / (2.0 * (1.0 + b2)) * di[i]);
  term_mean *= 2.0 * std::pow(1.0 + b2, -p / 2.0) / n;

  double term_const = std::pow(1.0 + 2.0 * b2, -p / 2.0);

  double hz = n * (term_pairs - term_mean + term_const);

  // Media/varianza asintotica del null (Henze & Zirkler 1990), p=2.
  double w_b  = (1.0 + b2) * (1.0 + 3.0 * b2);
  double mu   = 1.0 - std::pow(1.0 + 2.0 * b2, -p / 2.0) *
                        (1.0 + p * b2 / (1.0 + 2.0 * b2) +
                         p * (p + 2.0) * b2 * b2 / (2.0 * (1.0 + 2.0 * b2) * (1.0 + 2.0 * b2)));
  double sigma2 =
      2.0 * std::pow(1.0 + 4.0 * b2, -p / 2.0) +
      2.0 * std::pow(1.0 + 2.0 * b2, -p) *
          (1.0 + 2.0 * p * b2 * b2 / ((1.0 + 2.0 * b2) * (1.0 + 2.0 * b2)) +
           3.0 * p * (p + 2.0) * b2 * b2 * b2 * b2 /
               (4.0 * (1.0 + 2.0 * b2) * (1.0 + 2.0 * b2) * (1.0 + 2.0 * b2) * (1.0 + 2.0 * b2))) -
      4.0 * std::pow(w_b, -p / 2.0) *
          (1.0 + 3.0 * p * b2 * b2 / (2.0 * w_b) +
           p * (p + 2.0) * b2 * b2 * b2 * b2 / (2.0 * w_b * w_b));

  if (sigma2 <= 0.0 || mu <= 0.0)
    return res;

  double sigma2_ln = std::log(1.0 + sigma2 / (mu * mu));
  double mu_ln     = std::log(mu) - 0.5 * sigma2_ln;
  double z         = (std::log(hz) - mu_ln) / std::sqrt(sigma2_ln);
  double p_value   = 0.5 * std::erfc(z / std::sqrt(2.0));

  res.hz      = hz;
  res.p_value = p_value;
  res.ok      = true;
  return res;
}

struct CliConfig {
  std::string input_path  = "output/lens_simulation/lens.root";
  std::string output_path = "output/psf_analysis/psf_gaussianity.png";
  double sigma_cut         = 3.0;
  int n_iter               = 4;
  int min_hits             = 30;
  double contamination_alpha = 0.05;
  int hz_subsample_size      = 2000;
};

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
    if (key == "--input")
      cfg.input_path = next();
    else if (key == "--output")
      cfg.output_path = next();
    else if (key == "--sigma-cut")
      cfg.sigma_cut = std::stod(next());
    else if (key == "--n-iter")
      cfg.n_iter = std::stoi(next());
    else if (key == "--min-hits")
      cfg.min_hits = std::stoi(next());
    else if (key == "--contamination-alpha")
      cfg.contamination_alpha = std::stod(next());
    else if (key == "--hz-subsample-size")
      cfg.hz_subsample_size = std::stoi(next());
    else {
      std::cerr << "Opzione sconosciuta: " << key << "\n";
      std::exit(1);
    }
  }
  return cfg;
}

static void apply_style() {
  gStyle->Reset();
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleFont(42, "");
  gStyle->SetStatFont(42);
  gStyle->SetTitleOffset(1.5, "Y");
  gStyle->SetTitleOffset(1.2, "X");
  gStyle->SetTitleOffset(1.4, "Z");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetLineWidth(2);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}

// Autovalori/autovettore principale di una matrice 2x2 simmetrica
// [[var_y, cov_yz],[cov_yz, var_z]] -> semiassi (sqrt(lambda1,2)) e angolo
// (in gradi) dell'asse maggiore rispetto all'asse y.
struct EllipseParams {
  double a1, a2;   // semiassi (1 sigma)
  double angle_deg;
};

static EllipseParams ellipse_from_cov(double var_y, double var_z, double cov_yz) {
  double tr    = var_y + var_z;
  double det   = var_y * var_z - cov_yz * cov_yz;
  double disc  = std::sqrt(std::max(0.0, (tr * tr) / 4.0 - det));
  double lam1  = tr / 2.0 + disc;
  double lam2  = std::max(0.0, tr / 2.0 - disc);
  double theta = 0.5 * std::atan2(2.0 * cov_yz, var_y - var_z);
  return {std::sqrt(std::max(0.0, lam1)), std::sqrt(lam2), theta * 180.0 / TMath::Pi()};
}

int main(int argc, char** argv) {
  CliConfig cli = parse_args(argc, argv);

  TFile* fin = TFile::Open(cli.input_path.c_str(), "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Errore: impossibile aprire " << cli.input_path << "\n";
    return 1;
  }

  TTree* tConfig = dynamic_cast<TTree*>(fin->Get("Configurations"));
  TTree* tRuns   = dynamic_cast<TTree*>(fin->Get("Runs"));
  if (!tConfig || !tRuns) {
    std::cerr << "Errore: TTree Configurations o Runs mancanti in " << cli.input_path << "\n";
    return 1;
  }

  // Preferisce il piano virtuale non troncato (PsfDofRuns/y_virtual_hits,z_virtual_hits,
  // Tecnica E) se disponibile; altrimenti fallback sugli hit fisici del fotocatodo
  // (Runs/y_hits,z_hits), troncati dall'apertura del sensore.
  TTree* tVRuns    = dynamic_cast<TTree*>(fin->Get("PsfDofRuns"));
  bool use_virtual = tVRuns && tVRuns->GetBranch("y_virtual_hits") && tVRuns->GetBranch("z_virtual_hits");
  TTree* tHits     = use_virtual ? tVRuns : tRuns;
  if (use_virtual) {
    spdlog::info("Uso hit non troncati dal piano virtuale (PsfDofRuns/y_virtual_hits,z_virtual_hits)");
  } else {
    spdlog::warn(
        "PsfDofRuns/y_virtual_hits assenti in {} - fallback su Runs/y_hits (fisici, troncati "
        "dal fotocatodo). Rigenera lens.root con lens_save_virtual_hits=true per l'analisi "
        "non troncata.",
        cli.input_path);
  }

  int cfg_id;
  double cfg_x1, cfg_x2;
  char cfg_l1_id[256] = {}, cfg_l2_id[256] = {};
  tConfig->SetBranchAddress("config_id", &cfg_id);
  tConfig->SetBranchAddress("x1", &cfg_x1);
  tConfig->SetBranchAddress("x2", &cfg_x2);
  bool has_cfg_ids = (tConfig->GetBranch("l1_id") && tConfig->GetBranch("l2_id"));
  if (has_cfg_ids) {
    tConfig->SetBranchAddress("l1_id", cfg_l1_id);
    tConfig->SetBranchAddress("l2_id", cfg_l2_id);
  }

  std::map<int, std::tuple<double, double, std::string, std::string>> config_map;
  for (Long64_t i = 0; i < tConfig->GetEntries(); ++i) {
    tConfig->GetEntry(i);
    config_map[cfg_id] = {cfg_x1, cfg_x2, has_cfg_ids ? cfg_l1_id : "", has_cfg_ids ? cfg_l2_id : ""};
  }

  int run_id_r, run_cfg_id;
  double n_hits_r, x_source_r, y_source_r;
  std::vector<float>* y_hits_ptr = nullptr;
  std::vector<float>* z_hits_ptr = nullptr;
  tHits->SetBranchAddress("run_id", &run_id_r);
  tHits->SetBranchAddress("config_id", &run_cfg_id);
  tHits->SetBranchAddress("x_source", &x_source_r);
  tHits->SetBranchAddress("y_source", &y_source_r);
  tHits->SetBranchAddress("n_hits", &n_hits_r);
  tHits->SetBranchAddress(use_virtual ? "y_virtual_hits" : "y_hits", &y_hits_ptr);
  tHits->SetBranchAddress(use_virtual ? "z_virtual_hits" : "z_hits", &z_hits_ptr);

  std::vector<double> d2_pool;
  int n_runs_used = 0, n_runs_skipped = 0;

  // Run "rappresentativo": quello con piu' hit grezzi, usato per lo scatter
  // 2D e l'istogramma marginale (mostra la PSF meglio popolata).
  std::vector<Hit> best_hits;
  PSFResult best_res{};
  double best_x1 = 0.0, best_x2 = 0.0;
  std::string best_l1, best_l2;

  Long64_t n_runs = tHits->GetEntries();
  for (Long64_t r = 0; r < n_runs; ++r) {
    tHits->GetEntry(r);
    auto it = config_map.find(run_cfg_id);
    if (it == config_map.end())
      continue;

    std::vector<Hit> hits_buf;
    if (y_hits_ptr && z_hits_ptr) {
      hits_buf.reserve(y_hits_ptr->size());
      for (size_t j = 0; j < y_hits_ptr->size(); ++j)
        hits_buf.push_back(
            {static_cast<double>((*y_hits_ptr)[j]), static_cast<double>((*z_hits_ptr)[j])});
    }

    if (static_cast<int>(hits_buf.size()) < cli.min_hits) {
      ++n_runs_skipped;
      continue;
    }

    PSFResult res{};
    bool ok = compute_psf(hits_buf, res, cli.sigma_cut, cli.n_iter, 1.0);
    if (!ok) {
      ++n_runs_skipped;
      continue;
    }

    for (const auto& h : hits_buf)
      d2_pool.push_back(mahalanobis_d2(h, res));

    ++n_runs_used;

    if (hits_buf.size() > best_hits.size()) {
      best_hits = hits_buf;
      best_res  = res;
      best_x1   = std::get<0>(it->second);
      best_x2   = std::get<1>(it->second);
      best_l1   = std::get<2>(it->second);
      best_l2   = std::get<3>(it->second);
    }
  }

  fin->Close();

  if (d2_pool.empty() || best_hits.empty()) {
    std::cerr << "Errore: nessun run valido (min_hits=" << cli.min_hits << ")\n";
    return 1;
  }

  std::cout << "Run usati: " << n_runs_used << ", run scartati: " << n_runs_skipped << "\n";
  std::cout << "Hit pooled per Q-Q plot: " << d2_pool.size() << "\n";

  // --- Q-Q plot: D^2 pooled (ordinato) vs quantili teorici chi^2_2 ---
  std::sort(d2_pool.begin(), d2_pool.end());
  int N = static_cast<int>(d2_pool.size());

  std::vector<double> q_theory(N), q_sample(N);
  for (int i = 0; i < N; ++i) {
    double p     = (i + 0.5) / static_cast<double>(N);
    q_theory[i]  = TMath::ChisquareQuantile(p, 2);
    q_sample[i]  = d2_pool[i];
  }

  // --- Curtosi multivariata di Mardia (p=2) sul campione completo (include
  // eventuale contaminazione: riferimento per confronto, comportamento invariato) ---
  MardiaStats mardia_full = mardia_stats(d2_pool);

  std::cout << "Mardia b2 (campione completo) = " << mardia_full.b2
            << ", gamma2 = " << mardia_full.gamma2 << " +/- " << mardia_full.se_gamma2 << "\n";

  // --- Split nucleo/contaminazione: cutoff Bonferroni-corretto su chi^2_2 ---
  double d2_cutoff = bonferroni_chi2_2_cutoff(d2_pool.size(), cli.contamination_alpha);
  auto it_cutoff    = std::upper_bound(d2_pool.begin(), d2_pool.end(), d2_cutoff);
  std::vector<double> d2_core(d2_pool.begin(), it_cutoff);
  int n_core           = static_cast<int>(d2_core.size());
  int n_contaminated   = N - n_core;
  double contam_frac   = static_cast<double>(n_contaminated) / static_cast<double>(N);
  MardiaStats mardia_core = mardia_stats(d2_core);

  std::cout << "Cutoff Bonferroni D^2 (alpha_fw=" << cli.contamination_alpha << ", N=" << N
            << ") = " << d2_cutoff << "\n";
  std::cout << "Contaminazione: " << n_contaminated << "/" << N << " (" << (contam_frac * 100.0)
            << "%)\n";
  std::cout << "Mardia b2 (nucleo) = " << mardia_core.b2 << ", gamma2 = " << mardia_core.gamma2
            << " +/- " << mardia_core.se_gamma2 << "\n";

  // --- Henze-Zirkler (diagnostica secondaria) su sottocampione del run
  // rappresentativo, filtrato allo stesso cutoff nucleo/contaminazione ---
  std::vector<Hit> best_core;
  best_core.reserve(best_hits.size());
  for (const auto& h : best_hits)
    if (mahalanobis_d2(h, best_res) <= d2_cutoff)
      best_core.push_back(h);

  std::vector<Hit> hz_sample;
  int hz_n_req = std::min(cli.hz_subsample_size, static_cast<int>(best_core.size()));
  hz_sample.reserve(hz_n_req);
  std::mt19937 hz_rng(42);
  std::sample(best_core.begin(), best_core.end(), std::back_inserter(hz_sample), hz_n_req, hz_rng);
  HZResult hz = henze_zirkler(hz_sample);

  if (hz.ok)
    std::cout << "Henze-Zirkler (sottocampione n=" << hz_sample.size() << ") = " << hz.hz
              << ", p-value = " << hz.p_value << "\n";
  else
    std::cout << "Henze-Zirkler: sottocampione insufficiente o degenere, saltato\n";

  // --- Canvas ---
  apply_style();
  TCanvas* c = new TCanvas("c", "psf_gaussianity", 1500, 750);
  c->SetLeftMargin(0.0);
  c->SetRightMargin(0.0);
  c->SetTopMargin(0.0);
  c->SetBottomMargin(0.0);

  TPad* pad_scatter = new TPad("pad_scatter", "", 0.00, 0.24, 0.42, 1.00);
  TPad* pad_marg    = new TPad("pad_marg", "", 0.42, 0.24, 0.70, 1.00);
  TPad* pad_qq      = new TPad("pad_qq", "", 0.70, 0.24, 1.00, 1.00);
  TPad* pad_info    = new TPad("pad_info", "", 0.00, 0.00, 1.00, 0.24);

  pad_scatter->SetLeftMargin(0.16);
  pad_scatter->SetRightMargin(0.03);
  pad_scatter->SetTopMargin(0.10);
  pad_scatter->SetBottomMargin(0.14);

  // Margine sinistro maggiore di pad_scatter: pad_marg e pad_qq sono piu' stretti
  // (larghezza 0.28 e 0.30 del canvas contro 0.42 di pad_scatter), quindi a parita' di
  // frazione avrebbero meno spazio assoluto in pixel per titolo+numeri dell'asse Y,
  // causando sovrapposizione (specialmente con notazione scientifica sui conteggi).
  pad_marg->SetLeftMargin(0.24);
  pad_marg->SetRightMargin(0.03);
  pad_marg->SetTopMargin(0.10);
  pad_marg->SetBottomMargin(0.14);

  pad_qq->SetLeftMargin(0.22);
  pad_qq->SetRightMargin(0.04);
  pad_qq->SetTopMargin(0.10);
  pad_qq->SetBottomMargin(0.14);

  pad_info->SetLeftMargin(0.02);
  pad_info->SetRightMargin(0.02);

  c->cd();
  pad_scatter->Draw();
  pad_marg->Draw();
  pad_qq->Draw();
  pad_info->Draw();

  // --- Pannello 1: scatter 2D con ellissi 1/2/3 sigma ---
  pad_scatter->cd();
  std::vector<double> ys(best_hits.size()), zs(best_hits.size());
  for (size_t i = 0; i < best_hits.size(); ++i) {
    ys[i] = best_hits[i].y;
    zs[i] = best_hits[i].z;
  }

  EllipseParams ep = ellipse_from_cov(best_res.cov_yy, best_res.cov_zz, best_res.cov_yz);

  // Range esplicito = unione tra estensione del nucleo (best_core, D^2 <= cutoff Bonferroni)
  // e bounding box dell'ellisse 3 sigma (ruotata di ep.angle_deg), + margine 10%. Si usa il
  // nucleo e non l'intero best_hits per non far esplodere lo zoom per via di pochi hit di
  // contaminazione (riflessioni spurie) molto lontani dal centro: quei punti restano comunque
  // disegnati nello scatter, solo clippati dal pad se fuori range.
  double theta_rad = ep.angle_deg * TMath::Pi() / 180.0;
  double bbox_y_3s = std::sqrt(std::pow(3.0 * ep.a1 * std::cos(theta_rad), 2) +
                                std::pow(3.0 * ep.a2 * std::sin(theta_rad), 2));
  double bbox_z_3s = std::sqrt(std::pow(3.0 * ep.a1 * std::sin(theta_rad), 2) +
                                std::pow(3.0 * ep.a2 * std::cos(theta_rad), 2));

  std::vector<double> core_ys(best_core.size()), core_zs(best_core.size());
  for (size_t i = 0; i < best_core.size(); ++i) {
    core_ys[i] = best_core[i].y;
    core_zs[i] = best_core[i].z;
  }

  double y_lo = std::min(*std::min_element(core_ys.begin(), core_ys.end()), best_res.mean_y - bbox_y_3s);
  double y_hi = std::max(*std::max_element(core_ys.begin(), core_ys.end()), best_res.mean_y + bbox_y_3s);
  double z_lo = std::min(*std::min_element(core_zs.begin(), core_zs.end()), best_res.mean_z - bbox_z_3s);
  double z_hi = std::max(*std::max_element(core_zs.begin(), core_zs.end()), best_res.mean_z + bbox_z_3s);

  double y_margin1 = 0.10 * std::max(1.0, y_hi - y_lo);
  double z_margin1 = 0.10 * std::max(1.0, z_hi - z_lo);
  y_lo -= y_margin1;
  y_hi += y_margin1;
  z_lo -= z_margin1;
  z_hi += z_margin1;

  TH2F* h_frame1 = new TH2F("h_frame1", "", 100, y_lo, y_hi, 100, z_lo, z_hi);
  h_frame1->GetXaxis()->SetTitle("y_{hit} [mm]");
  h_frame1->GetYaxis()->SetTitle("z_{hit} [mm]");
  h_frame1->GetXaxis()->SetTitleFont(42);
  h_frame1->GetYaxis()->SetTitleFont(42);
  h_frame1->Draw("AXIS");

  TGraph* g_scatter = new TGraph(static_cast<int>(ys.size()), ys.data(), zs.data());
  g_scatter->SetMarkerStyle(6);
  g_scatter->SetMarkerColor(kAzure + 2);
  g_scatter->Draw("P same");

  int colors[3] = {kGreen + 2, kOrange + 1, kRed + 1};
  for (int n = 1; n <= 3; ++n) {
    TEllipse* el = new TEllipse(best_res.mean_y, best_res.mean_z, ep.a1 * n, ep.a2 * n, 0, 360,
                                ep.angle_deg);
    el->SetFillStyle(0);
    el->SetLineColor(colors[n - 1]);
    el->SetLineWidth(2);
    el->Draw("same");
  }

  TLatex tl1;
  tl1.SetNDC();
  tl1.SetTextFont(42);
  tl1.SetTextSize(0.045);
  tl1.SetTextAlign(22);
  tl1.DrawLatex(0.55, 0.955, "Distribuzione hit sul detector");

  // --- Pannello 2: istogramma marginale (asse y) con fit gaussiano ---
  pad_marg->cd();
  // Range dal nucleo (best_core), non dall'intero best_hits, per lo stesso motivo del
  // pannello 1 (zoom sulla gaussiana invece che sulla contaminazione). L'istogramma
  // riempie comunque tutti gli hit di best_hits.
  double y_min = *std::min_element(core_ys.begin(), core_ys.end());
  double y_max = *std::max_element(core_ys.begin(), core_ys.end());
  double pad_w = 0.15 * std::max(1.0, y_max - y_min);
  TH1D* h_marg = new TH1D("h_marg", "", 60, y_min - pad_w, y_max + pad_w);
  for (double y : ys)
    h_marg->Fill(y);
  h_marg->SetLineColor(kAzure + 2);
  h_marg->SetLineWidth(2);
  h_marg->GetXaxis()->SetTitle("y_{hit} [mm]");
  h_marg->GetYaxis()->SetTitle("Conteggi");
  h_marg->GetXaxis()->SetTitleFont(42);
  h_marg->GetYaxis()->SetTitleFont(42);
  h_marg->Draw("HIST");

  TF1* f_gaus = new TF1("f_gaus", "gaus", y_min - pad_w, y_max + pad_w);
  f_gaus->SetLineColor(kRed + 1);
  f_gaus->SetLineWidth(2);
  h_marg->Fit(f_gaus, "QR");
  f_gaus->Draw("same");

  TLatex tl2;
  tl2.SetNDC();
  tl2.SetTextFont(42);
  tl2.SetTextSize(0.045);
  tl2.SetTextAlign(22);
  tl2.DrawLatex(0.55, 0.955, "Profilo marginale + fit gaussiano");

  // --- Pannello 3: Q-Q plot D^2 vs chi^2_2 ---
  pad_qq->cd();

  // Banda di confidenza 95% puntuale (asintotica, delta method sulle statistiche
  // d'ordine): se(q_theory_(i)) = sqrt(p_i(1-p_i)/n) / f_chi2_2(q_theory_(i)).
  std::vector<double> band_lo(N), band_hi(N);
  for (int i = 0; i < N; ++i) {
    double p    = (i + 0.5) / static_cast<double>(N);
    double dens = chi2_2_pdf(q_theory[i]);
    double se   = std::sqrt(p * (1.0 - p) / static_cast<double>(N)) / std::max(dens, 1e-12);
    band_lo[i]  = q_theory[i] - 1.96 * se;
    band_hi[i]  = q_theory[i] + 1.96 * se;
  }

  // Range esplicito: unione tra quantili teorici/campionari e banda, limitata al nucleo
  // (primi n_core punti, D^2 <= cutoff Bonferroni) + margine 10%. Il pool completo include
  // la contaminazione (D^2 fino a migliaia), che schiaccerebbe il nucleo vicino a zero se
  // usato per definire il range - i punti disegnati restano tutti (g_qq, g_band, line_id),
  // quelli oltre il nuovo range vengono clippati dal pad.
  double band_hi_max = *std::max_element(band_hi.begin(), band_hi.begin() + n_core);
  double band_lo_min = *std::min_element(band_lo.begin(), band_lo.begin() + n_core);
  double qq_lo = std::min(0.0, band_lo_min);
  double qq_hi = std::max({q_theory[n_core - 1], q_sample[n_core - 1], band_hi_max});
  double qq_margin = 0.10 * std::max(1.0, qq_hi - qq_lo);
  qq_lo -= qq_margin;
  qq_hi += qq_margin;

  TH2F* h_frame3 = new TH2F("h_frame3", "", 100, qq_lo, qq_hi, 100, qq_lo, qq_hi);
  h_frame3->GetXaxis()->SetTitle("Quantile teorico #chi^{2}_{2}");
  h_frame3->GetYaxis()->SetTitle("D^{2} campionario (ordinato)");
  h_frame3->GetXaxis()->SetTitleFont(42);
  h_frame3->GetYaxis()->SetTitleFont(42);
  h_frame3->Draw("AXIS");

  // Banda ombreggiata 95%: poligono chiuso, bordo superiore poi inferiore invertito.
  std::vector<double> band_x(2 * N), band_y(2 * N);
  for (int i = 0; i < N; ++i) {
    band_x[i]              = q_theory[i];
    band_y[i]               = band_hi[i];
    band_x[2 * N - 1 - i]  = q_theory[i];
    band_y[2 * N - 1 - i]  = band_lo[i];
  }
  TGraph* g_band = new TGraph(2 * N, band_x.data(), band_y.data());
  g_band->SetFillColorAlpha(kAzure + 2, 0.25);
  g_band->SetLineColorAlpha(kAzure + 2, 0.0);
  g_band->Draw("F same");

  TGraph* g_qq = new TGraph(N, q_theory.data(), q_sample.data());
  g_qq->SetMarkerStyle(7);
  g_qq->SetMarkerColor(kAzure + 2);
  g_qq->Draw("P same");

  TLine* line_id = new TLine(0, 0, qq_hi, qq_hi);
  line_id->SetLineColor(kRed + 1);
  line_id->SetLineStyle(2);
  line_id->SetLineWidth(2);
  line_id->Draw("same");

  // Cutoff Bonferroni nucleo/contaminazione (orizzontale, sull'asse D^2 campionario)
  TLine* line_cutoff = new TLine(qq_lo, d2_cutoff, qq_hi, d2_cutoff);
  line_cutoff->SetLineColor(kMagenta + 2);
  line_cutoff->SetLineStyle(7);
  line_cutoff->SetLineWidth(2);
  line_cutoff->Draw("same");

  TLatex tl3;
  tl3.SetNDC();
  tl3.SetTextFont(42);
  tl3.SetTextSize(0.045);
  tl3.SetTextAlign(22);
  tl3.DrawLatex(0.55, 0.955, "Q-Q plot vs #chi^{2}_{2}");

  // --- Pannello metadati ---
  pad_info->cd();
  TPaveText* pt = new TPaveText(0.02, 0.05, 0.98, 0.95, "NDC");
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->SetTextFont(42);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.15);

  std::ostringstream line;
  line << "Coppia: " << best_l1 << " + " << best_l2 << "   x_{1}=" << best_x1
       << " mm, x_{2}=" << best_x2 << " mm    |    N_{run}=" << n_runs_used
       << ", N_{hit} pooled=" << N;
  pt->AddText(line.str().c_str());

  std::ostringstream line_contam;
  line_contam << "Cutoff Bonferroni D^{2} > " << std::fixed << std::setprecision(2) << d2_cutoff
              << " (#alpha_{fw}=" << cli.contamination_alpha << ")   |   Contaminazione: "
              << n_contaminated << "/" << N << " (" << std::setprecision(3) << contam_frac * 100.0
              << "%)";
  pt->AddText(line_contam.str().c_str());

  std::ostringstream line3;
  line3 << "Mardia #gamma_{2} nucleo (N_{core}=" << n_core << ") = " << std::fixed
        << std::setprecision(4) << mardia_core.gamma2 << " #pm " << mardia_core.se_gamma2
        << "   (atteso 0 per gaussiana bivariata)";
  pt->AddText(line3.str().c_str());
  pt->Draw();

  c->Update();
  std::filesystem::create_directories(std::filesystem::path(cli.output_path).parent_path());
  c->Print(cli.output_path.c_str());

  std::cout << "Output salvato in: " << cli.output_path << "\n";
  return 0;
}
