#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <vector>

double moving_average(const std::vector<double>& data, int idx, int window) {
  int half  = window / 2;
  int start = std::max(0, idx - half);
  int end   = std::min(static_cast<int>(data.size()) - 1, idx + half);

  double sum = 0;
  for (int i = start; i <= end; i++)
    sum += data[i];
  return sum / (end - start + 1);
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "Usage: ./geant4_heatmap lens.root\n";
    return 1;
  }

  gStyle->SetOptStat(0);

  TFile f(argv[1]);
  if (f.IsZombie()) {
    std::cout << "Cannot open file\n";
    return 1;
  }

  TTree* tConfig = (TTree*)f.Get("Configurations");
  TTree* tRuns   = (TTree*)f.Get("Runs");
  TTree* tHits   = (TTree*)f.Get("Hits");

  if (!tConfig || !tRuns || !tHits) {
    std::cout << "Missing trees!\n";
    return 1;
  }

  // =============================
  // LETTURA CONFIGURAZIONI
  // =============================
  int config_id;
  double x1, x2;
  tConfig->SetBranchAddress("config_id", &config_id);
  tConfig->SetBranchAddress("x1", &x1);
  tConfig->SetBranchAddress("x2", &x2);

  std::map<int, std::pair<double, double>> config_positions;
  std::set<double> x1_unique, x2_unique;

  for (Long64_t i = 0; i < tConfig->GetEntries(); i++) {
    tConfig->GetEntry(i);
    config_positions[config_id] = {x1, x2};
    x1_unique.insert(x1);
    x2_unique.insert(x2);
  }

  std::vector<double> x1_vals(x1_unique.begin(), x1_unique.end());
  std::vector<double> x2_vals(x2_unique.begin(), x2_unique.end());

  int nbins_x1 = x1_vals.size();
  int nbins_x2 = x2_vals.size();

  // =============================
  // RUNS
  // =============================
  int run_id, run_config_id;
  double x_source;
  tRuns->SetBranchAddress("run_id", &run_id);
  tRuns->SetBranchAddress("config_id", &run_config_id);
  tRuns->SetBranchAddress("x_source", &x_source);

  std::map<int, std::vector<std::pair<int, double>>> config_runs;
  for (Long64_t i = 0; i < tRuns->GetEntries(); i++) {
    tRuns->GetEntry(i);
    config_runs[run_config_id].push_back({run_id, x_source});
  }

  // =============================
  // HITS
  // =============================
  int hit_run_id;
  double y_hit;
  tHits->SetBranchAddress("run_id", &hit_run_id);
  tHits->SetBranchAddress("y_hit", &y_hit);

  std::map<int, std::vector<double>> run_hits;
  for (Long64_t i = 0; i < tHits->GetEntries(); i++) {
    tHits->GetEntry(i);
    run_hits[hit_run_id].push_back(y_hit);
  }

  // =============================
  // CREAZIONE HEATMAP
  // =============================
  TH2D* h = new TH2D("heatmap", "RMSE from Geant4;X1 (mm);X2 (mm)", nbins_x1, x1_vals.front() - 0.5,
                     x1_vals.back() + 0.5, nbins_x2, x2_vals.front() - 0.5, x2_vals.back() + 0.5);

  auto find_bin = [](const std::vector<double>& vals, double v) {
    auto it = std::lower_bound(vals.begin(), vals.end(), v);
    return std::distance(vals.begin(), it) + 1;
  };

  // =============================
  // LOOP CONFIGURAZIONI CON MEDIA MOBILE
  // =============================
  int window = 3; // finestra della media mobile (regolabile)
  for (auto& cfg : config_runs) {
    int cfg_id = cfg.first;
    auto& runs = cfg.second;

    std::vector<double> source_positions;
    std::vector<double> mean_detector;

    for (auto& r : runs) {
      int rid    = r.first;
      auto& hits = run_hits[rid];
      if (hits.empty())
        continue;

      double sum = 0;
      for (double v : hits)
        sum += v;
      double mean = sum / hits.size();

      source_positions.push_back(r.second);
      mean_detector.push_back(mean);
    }

    if (source_positions.size() < 2)
      continue;

    // Applica media mobile
    std::vector<double> smoothed;
    for (size_t i = 0; i < mean_detector.size(); i++) {
      smoothed.push_back(moving_average(mean_detector, i, window));
    }

    // Fit lineare su dati smoothed
    TGraph g(smoothed.size());
    for (size_t i = 0; i < smoothed.size(); i++)
      g.SetPoint(i, source_positions[i], smoothed[i]);
    TF1 fit("fit", "pol1");
    g.Fit(&fit, "Q");

    double sumsq = 0;
    for (size_t i = 0; i < smoothed.size(); i++) {
      double fitted = fit.Eval(source_positions[i]);
      sumsq += pow(smoothed[i] - fitted, 2);
    }

    double rmse = sqrt(sumsq / smoothed.size());

    auto [pos_x1, pos_x2] = config_positions[cfg_id];
    int bin_x             = find_bin(x1_vals, pos_x1);
    int bin_y             = find_bin(x2_vals, pos_x2);
    h->SetBinContent(bin_x, bin_y, rmse);
  }

  // =============================
  // PLOT
  // =============================
  TCanvas* c = new TCanvas("c", "heatmap", 1000, 700);
  h->SetContour(255);
  h->Draw("COLZ");
  c->SaveAs("output/lens_simulation/geant4_heatmap.png");

  std::cout << "Heatmap created from Geant4 data (smoothed).\n";
}