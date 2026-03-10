#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

struct Hit {
  double y;
  double z;
};

int main(int argc, char** argv) {
  if (argc != 4) {
    std::cerr << "Uso: " << argv[0] << " x1_target x2_target y0_target\n";
    return 1;
  }

  double x1_target = std::stod(argv[1]);
  double x2_target = std::stod(argv[2]);
  double y0_target = std::stod(argv[3]);

  // ---------------------- apri file ROOT ----------------------
  TFile* file = TFile::Open("output/lens_simulation/lens.root");
  if (!file || file->IsZombie()) {
    std::cerr << "Errore nell'aprire il file ROOT\n";
    return 1;
  }

  TTree* tConfig = (TTree*)file->Get("Configurations");
  TTree* tRuns   = (TTree*)file->Get("Runs");
  TTree* tHits   = (TTree*)file->Get("Hits");
  if (!tConfig || !tRuns || !tHits) {
    std::cerr << "Impossibile leggere uno dei TTree\n";
    return 1;
  }

  // ---------------------- setup branch ----------------------
  int config_id;
  double x1, x2;
  tConfig->SetBranchAddress("config_id", &config_id);
  tConfig->SetBranchAddress("x1", &x1);
  tConfig->SetBranchAddress("x2", &x2);

  int run_id, run_config;
  double y_source;
  tRuns->SetBranchAddress("run_id", &run_id);
  tRuns->SetBranchAddress("config_id", &run_config);
  tRuns->SetBranchAddress("x_source", &y_source);

  int hit_run;
  double y_hit, z_hit;
  tHits->SetBranchAddress("run_id", &hit_run);
  tHits->SetBranchAddress("y_hit", &y_hit);
  tHits->SetBranchAddress("z_hit", &z_hit);

  // ---------------------- trova configurazione più vicina ----------------------
  int target_config_id = -1;
  double best_dist     = std::numeric_limits<double>::max();
  double x1_sel = 0, x2_sel = 0;
  for (Long64_t i = 0; i < tConfig->GetEntries(); ++i) {
    tConfig->GetEntry(i);
    double dist =
        std::sqrt((x1 - x1_target) * (x1 - x1_target) + (x2 - x2_target) * (x2 - x2_target));
    if (dist < best_dist) {
      best_dist        = dist;
      target_config_id = config_id;
      x1_sel           = x1;
      x2_sel           = x2;
    }
  }
  if (target_config_id < 0) {
    std::cerr << "Nessuna configurazione trovata!\n";
    return 1;
  }
  std::cout << "Configurazione scelta: config_id=" << target_config_id << ", x1=" << x1_sel
            << ", x2=" << x2_sel << "\n";

  // ---------------------- trova run esatto ----------------------
  int exact_run_id = -1;
  for (Long64_t i = 0; i < tRuns->GetEntries(); ++i) {
    tRuns->GetEntry(i);
    if (run_config == target_config_id && std::abs(y_source - y0_target) < 1e-6) {
      exact_run_id = run_id;
      break;
    }
  }
  if (exact_run_id < 0) {
    std::cerr << "Nessun run trovato con y0_target = " << y0_target << "\n";
    return 1;
  }
  std::cout << "Run selezionato: run_id=" << exact_run_id << ", y_source=" << y0_target << "\n";

  // ---------------------- raccogli tutte le hit ----------------------
  std::vector<Hit> hits;
  for (Long64_t i = 0; i < tHits->GetEntries(); ++i) {
    tHits->GetEntry(i);
    if (hit_run == exact_run_id)
      hits.push_back({y_hit, z_hit});
  }
  if (hits.empty()) {
    std::cerr << "Nessuna hit trovata!\n";
    return 1;
  }

  // ---------------------- calcola media e deviazione standard ----------------------
  double y_sum = 0, z_sum = 0;
  for (auto& h : hits) {
    y_sum += h.y;
    z_sum += h.z;
  }
  double y_mean = y_sum / hits.size();
  double z_mean = z_sum / hits.size();

  double y_var = 0, z_var = 0;
  for (auto& h : hits) {
    y_var += (h.y - y_mean) * (h.y - y_mean);
    z_var += (h.z - z_mean) * (h.z - z_mean);
  }
  double y_sigma = std::sqrt(y_var / hits.size());
  double z_sigma = std::sqrt(z_var / hits.size());

  // ---------------------- filtra outlier (>3σ) ----------------------
  std::vector<Hit> hits_zoom;
  for (auto& h : hits) {
    if (std::abs(h.y - y_mean) <= 3 * y_sigma && std::abs(h.z - z_mean) <= 3 * z_sigma)
      hits_zoom.push_back(h);
  }

  // ridefinisci range per zoom
  double y_min = y_mean - 0.5 * y_sigma;
  double y_max = y_mean + 0.5 * y_sigma;
  double z_min = z_mean - 3 * z_sigma;
  double z_max = z_mean + 3 * z_sigma;

  // ---------------------- canvas e istogramma ----------------------
  TCanvas* c = new TCanvas("c", "Hits sul detector", 1200, 1000); // canvas più grande
  gStyle->SetOptStat(0);
  c->SetGrid();

  TH2D* hist = new TH2D("hist",
                        ("Hits: x1=" + std::to_string(x1_sel) + ", x2=" + std::to_string(x2_sel)
                         + ", y0=" + std::to_string(y0_target))
                            .c_str(),
                        200, y_min, y_max, 200, z_min, z_max);

  for (auto& h : hits_zoom)
    hist->Fill(h.y, h.z);

  hist->GetXaxis()->SetTitle("Y [mm]");
  hist->GetYaxis()->SetTitle("Z [mm]");
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetTitleOffset(1.2); // evita taglio
  hist->Draw("COLZ");

  std::ostringstream filename;
  filename << "output/mean_covariance_maps/detector_hits_config_" << target_config_id << "_y0_"
           << y0_target << ".png";
  c->SaveAs(filename.str().c_str());
  std::cout << "Immagine salvata in: " << filename.str() << "\n";

  return 0;
}