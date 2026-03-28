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

  // Legge i parametri dal file di configurazione (path parametrico)
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

  // Ottieni dimensioni reali dal detector (possono essere Thorlabs o GDML default)
  double h1      = det->GetLens75Thickness();
  double offset1 = det->GetLens75CenterOffset();
  double h2      = det->GetLens60Thickness();
  double offset2 = det->GetLens60CenterOffset();

  // Parametri sorgente GPS
  double source_y_min = 0.0;
  double source_y_max = 10.0 * sqrt(2.0);
  double source_dy    = 0.1;

  // Offset per config_id e run_id: permette il merge corretto di chunk paralleli
  int config_id_offset = config.value("config_id_offset", 0);
  int run_id_offset    = config.value("run_id_offset", 0);

  int config_counter = config_id_offset;
  int run_counter    = run_id_offset;

  // Costruisce la lista di coppie (x1, x2) da simulare.
  // Se il config contiene "pairs" (generato da run.sh per chunk bilanciati),
  // usa quella lista. Altrimenti usa il loop geometrico completo.
  struct LensPair {
    double x1, x2;
  };
  std::vector<LensPair> pairs;

  if (config.contains("pairs")) {
    for (const auto& p : config["pairs"]) {
      pairs.push_back({p[0].get<double>(), p[1].get<double>()});
    }
    spdlog::info("Loaded {} pairs from config", pairs.size());
  } else {
    // Loop geometrico adattivo: evita intersezioni basandosi sulle dimensioni reali
    // delle lenti (GDML default o Thorlabs) e mantiene l'ordine L1 < L2.
    const double margin = 3.0; // 3mm di spazio minimo tra le lenti per lens_scan

    double x1_start = x_min;
    double x1_end   = x_max;

    if (config.contains("x1_start"))
      x1_start = config["x1_start"];
    if (config.contains("x1_end"))
      x1_end = config["x1_end"];

    for (double x1 = x1_start; x1 <= x1_end + 1e-9; x1 += dx) {
      // Calcola x2_min per evitare collisioni:
      // back_l1 + margin < front_l2
      // (x1 + offset1 + h1/2) + margin < (x2 + offset2 - h2/2)
      double x2_min_collision = x1 + (offset1 - offset2) + (h1 + h2) / 2.0 + margin;
      double x2_start         = std::max(x1 + dx, x2_min_collision);

      for (double x2 = x2_start; x2 <= x_max + 1e-9; x2 += dx) {
        pairs.push_back({x1, x2});
      }
    }
    spdlog::info("Generated {} pairs from adaptive geometry loop", pairs.size());
  }

  for (const auto& pair : pairs) {
    double x1 = pair.x1;
    double x2 = pair.x2;

    det->SetLensPositions(x1, x2);
    run_manager->GeometryHasBeenModified();

    // Salva la configurazione delle lenti
    analysisManager->FillNtupleIColumn(0, 0, config_counter);
    analysisManager->FillNtupleDColumn(0, 1, x1);
    analysisManager->FillNtupleDColumn(0, 2, x2);
    analysisManager->AddNtupleRow(0);

    // Loop sulla posizione sorgente
    for (double y_source = source_y_min; y_source <= source_y_max + 1e-9; y_source += source_dy) {
      UImanager->ApplyCommand("/gps/pos/centre 0 " + std::to_string(y_source) + " 0 mm");

      int run_id = run_counter++;

      auto* eventAction = dynamic_cast<EventAction*>(
          const_cast<G4UserEventAction*>(run_manager->GetUserEventAction()));
      if (eventAction) {
        eventAction->runID = run_id;
      }

      long seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
      CLHEP::HepRandom::setTheSeed(seed);

      UImanager->ApplyCommand("/control/execute " + macro_to_run.string());

      int n_hits = eventAction ? eventAction->GetLastRunHitCount() : 0;

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

  // Scrive e chiude file
  analysisManager->Write();
  analysisManager->CloseFile();

  spdlog::info("Total runs completed: {}", run_counter - run_id_offset);
  spdlog::info("Optimization completed");
}

} // namespace riptide