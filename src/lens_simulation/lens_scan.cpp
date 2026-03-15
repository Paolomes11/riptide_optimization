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

#include "lens_scan.hpp"
#include "action_initialization.hpp"
#include "detector_construction.hpp"
#include "event_action.hpp"

#include <G4AnalysisManager.hh>
#include <G4RunManager.hh>
#include <G4UImanager.hh>

#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace riptide {

void lens_scan(G4RunManager* run_manager, const std::filesystem::path& macro_file,
               const std::string& root_output_file, const std::filesystem::path& config_file) {
  using json = nlohmann::json;

  // Crea la directory di output se non esiste
  std::filesystem::path output_path(root_output_file);
  if (output_path.has_parent_path()) {
    std::filesystem::create_directories(output_path.parent_path());
    spdlog::info("Output directory: {}", output_path.parent_path().string());
  }

  // Legge i parametri del file di configurazione
  std::ifstream f(config_file);
  if (!f.is_open())
    throw std::runtime_error("Impossibile aprire config: " + config_file.string());
  json config;
  f >> config;

  auto det = const_cast<DetectorConstruction*>(
      dynamic_cast<const DetectorConstruction*>(run_manager->GetUserDetectorConstruction()));
  if (!det)
    throw std::runtime_error("DetectorConstruction not found!");

  auto UImanager = G4UImanager::GetUIpointer();
  std::filesystem::path macro_to_run =
      macro_file.empty() ? "macros/lens_simulation.mac" : macro_file;

  // Apertura file ROOT e creazione delle Ntuple
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->OpenFile(root_output_file);

  // Compressione LZ4 livello 4 — encoding: algoritmo*100 + livello (404 = LZ4 lv4)
  analysisManager->SetCompressionLevel(404);

  // Ntuple 0: configurazioni lenti
  analysisManager->CreateNtuple("Configurations", "Lens configurations");
  analysisManager->CreateNtupleIColumn("config_id");
  analysisManager->CreateNtupleDColumn("x1");
  analysisManager->CreateNtupleDColumn("x2");
  analysisManager->FinishNtuple(0);

  // Ntuple 1: posizione sorgente / run
  analysisManager->CreateNtuple("Runs", "Source positions per configuration");
  analysisManager->CreateNtupleIColumn("run_id");
  analysisManager->CreateNtupleIColumn("config_id");
  analysisManager->CreateNtupleFColumn("x_source");
  analysisManager->CreateNtupleIColumn("n_hits");
  analysisManager->FinishNtuple(1);

  // Ntuple 2: hit/fotoni
  analysisManager->CreateNtuple("Hits", "Detected photons");
  analysisManager->CreateNtupleFColumn("y_hit");
  analysisManager->CreateNtupleFColumn("z_hit");
  analysisManager->FinishNtuple(2);

  // Parametri di scansione lenti e sorgente
  // Tiene conto di 30mm di scintillatore e 3mm di margine
  double x_min = config["x_min"];
  // Tiene conto di 3mm di margine
  double x_max = config["x_max"];
  double dx    = config["dx"];

  double r1 = config["r1"], h1 = config["h1"];
  double r2 = config["r2"], h2 = config["h2"];

  // Parametri sorgente GPS
  double source_y_min = 0.0;
  double source_y_max = 10.0;
  double source_dy    = 0.1;

  int config_counter = 0;
  int run_counter    = 0;

  for (double x1 = x_min - r1 + h1; x1 <= x_max - h2 - 1 - r1; x1 += dx) {
    double x2_min = x1 + r1 + r2 + 1;
    double x2_max = x_max + r2 - h2;

    for (double x2 = x2_min; x2 <= x2_max; x2 += dx) {
      det->SetLensPositions(x1, x2);
      run_manager->GeometryHasBeenModified();

      // Salva la configurazione delle lenti
      analysisManager->FillNtupleIColumn(0, 0, config_counter);
      analysisManager->FillNtupleDColumn(0, 1, x1);
      analysisManager->FillNtupleDColumn(0, 2, x2);
      analysisManager->AddNtupleRow(0);

      // Loop sulla posizione sorgente
      for (double y_source = source_y_min; y_source <= source_y_max; y_source += source_dy) {
        // Imposta posizione del centro del GPS
        UImanager->ApplyCommand("/gps/pos/centre 0 " + std::to_string(y_source) + " 0 mm");

        // Identificatore del run corrente
        int run_id = run_counter++;

        auto* eventAction = dynamic_cast<EventAction*>(
            const_cast<G4UserEventAction*>(run_manager->GetUserEventAction()));
        if (eventAction) {
          eventAction->runID = run_id; // assegna il run corrente
        }

        // Cambia seed
        long seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        CLHEP::HepRandom::setTheSeed(seed);

        // Esegue la macro
        UImanager->ApplyCommand("/control/execute " + macro_to_run.string());

        // legge n_hits scritte da EventAction durante questo run
        int n_hits = eventAction ? eventAction->GetLastRunHitCount() : 0;

        // Salva run nella Ntuple Runs
        analysisManager->FillNtupleIColumn(1, 0, run_id);
        analysisManager->FillNtupleIColumn(1, 1, config_counter);
        analysisManager->FillNtupleFColumn(1, 2, static_cast<float>(y_source));
        analysisManager->FillNtupleIColumn(1, 3, n_hits);
        analysisManager->AddNtupleRow(1);

        spdlog::info("Run done: config_id={}, y_source={} mm, run_id={}", config_counter, y_source,
                     run_id);
      }

      config_counter++;
    }
  }

  // Scrive e chiude file
  analysisManager->Write();
  analysisManager->CloseFile();

  spdlog::info("Optimization completed");
}

} // namespace riptide