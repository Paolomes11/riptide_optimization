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

#include "optimizer.hpp"
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

void run_optimization(G4RunManager* run_manager, const std::filesystem::path& macro_file,
                      const std::string& root_output_file) {
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

  // Apertura file di output e creazione Ntuple
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->OpenFile(root_output_file);

  analysisManager->CreateNtuple("events", "Eventi fotoni");
  analysisManager->CreateNtupleDColumn("x1");
  analysisManager->CreateNtupleDColumn("x2");
  analysisManager->CreateNtupleIColumn("config_id");
  analysisManager->FinishNtuple();

  // Parametri di ottimizzazione
  double x_min = config["x_min"];
  double x_max = config["x_max"];
  double dx    = config["dx"];

  double r1 = config["r1"], h1 = config["h1"];
  double r2 = config["r2"], h2 = config["h2"];

  auto UImanager                     = G4UImanager::GetUIpointer();
  std::filesystem::path macro_to_run = macro_file.empty() ? "macros/run1.mac" : macro_file;

  int config_counter = 0;
  for (double x1 = x_min - r1 + h1; x1 <= x_max - h2 - 1 - r1; x1 += dx) {
    double x2_min = x1 + r1 + r2 + 1;
    double x2_max = x_max + r2 - h2;

    for (double x2 = x2_min; x2 <= x2_max; x2 += dx) {
      det->SetLensPositions(x1, x2);
      run_manager->GeometryHasBeenModified();

      // Imposta l'identificatore della configurazione
      auto* eventAction = dynamic_cast<EventAction*>(
          const_cast<G4UserEventAction*>(run_manager->GetUserEventAction()));
      if (eventAction) {
        eventAction->SetConfigId(config_counter);
      }

      // CAMBIA SEED QUI
      long seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
      CLHEP::HepRandom::setTheSeed(seed);

      UImanager->ApplyCommand("/control/execute " + macro_to_run.string());

      spdlog::info("Simulation done for x1={} mm, x2={} mm, config_id={}", x1, x2, config_counter);
      config_counter++;
    }
  }

  // Scrive e chiude il file di output
  analysisManager->Write();
  analysisManager->CloseFile();

  spdlog::info("Optimization completed");
}

} // namespace riptide
