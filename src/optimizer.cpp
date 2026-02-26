#include "optimizer.hpp"
#include "action_initialization.hpp"
#include "detector_construction.hpp"
#include "efficiency_collector.hpp"

#include <G4RunManager.hh>
#include <G4UImanager.hh>

#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>
#include <fstream>

#include <TFile.h>
#include <TH2D.h>

#include <CLHEP/Random/Random.h>
#include <chrono>

namespace riptide {

void run_optimization(G4RunManager* run_manager, const std::filesystem::path& macro_file,
                      EfficiencyCollector* collector) {
  using json = nlohmann::json;

  // Legge i parametri del file di configurazione
  std::ifstream f("config/config.json");
  if (!f.is_open()) {
    throw std::runtime_error("Impossibile aprire config/config.json");
  }
  json config;
  f >> config;

  // Trova detector
  auto det = const_cast<DetectorConstruction*>(
      dynamic_cast<const DetectorConstruction*>(run_manager->GetUserDetectorConstruction()));
  if (!det) {
    throw std::runtime_error("DetectorConstruction not found!");
  }

  // Parametri di ottimizzazione
  double x_min  = config["x_min"];
  double x_max  = config["x_max"];
  double dx     = config["dx"];
  int n_photons = config["n_photons"];

  double r1 = config["r1"], h1 = config["h1"];
  double r2 = config["r2"], h2 = config["h2"];

  auto UImanager                     = G4UImanager::GetUIpointer();
  std::filesystem::path macro_to_run = macro_file.empty() ? "macros/run1.mac" : macro_file;

  double best_eff = -1.0;
  double best_x1 = 0.0, best_x2 = 0.0;

  // CREA TH2D
  double x1_min_fill = x_min - r1 + h1;
  double x1_max_fill = x_max - h2 - 1 - r1;
  double x2_min_fill = x1_min_fill + r1 + r2 + 1;
  double x2_max_fill = x_max + r2 - h2;
  int nx             = static_cast<int>((x1_max_fill - x1_min_fill) / dx) + 1;
  int ny             = static_cast<int>((x2_max_fill - x2_min_fill) / dx) + 1;
  auto h2_eff        = new TH2D("h2_eff", "Geometrical efficiency;x1 [mm];x2 [mm]", nx, x1_min_fill,
                                x1_max_fill, ny, x2_min_fill, x2_max_fill);

  for (double x1 = x_min - r1 + h1; x1 <= x_max - h2 - 1 - r1; x1 += dx) {
    double x2_min = x1 + r1 + r2 + 1;
    double x2_max = x_max + r2 - h2;

    for (double x2 = x2_min; x2 <= x2_max; x2 += dx) {
      det->SetLensPositions(x1, x2);
      run_manager->GeometryHasBeenModified();
      collector->reset();

      // CAMBIA SEED QUI
      long seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
      CLHEP::HepRandom::setTheSeed(seed);

      UImanager->ApplyCommand("/control/execute " + macro_to_run.string());

      double eff = collector->computeEfficiency(n_photons);
      if (eff > best_eff) {
        best_eff = eff;
        best_x1  = x1;
        best_x2  = x2;
      }

      // Riempie la TH2D
      h2_eff->Fill(x1, x2, eff);

      spdlog::info("x1={} mm, x2={} mm -> eff={}", x1, x2, eff);
    }
  }

  spdlog::info("Best configuration: x1={} mm, x2={} mm, eff={}", best_x1, best_x2, best_eff);

  // SALVA TH2D
  TFile out("output/optimization.root", "RECREATE");
  h2_eff->Write();
  out.Close();

  delete h2_eff;
}

} // namespace riptide