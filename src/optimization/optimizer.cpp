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
                      const std::string& root_output_file,
                      const std::filesystem::path& config_file) {
  using json = nlohmann::json;

  // Crea la directory di output se non esiste
  std::filesystem::path output_path(root_output_file);
  if (output_path.has_parent_path()) {
    std::filesystem::create_directories(output_path.parent_path());
    spdlog::info("Output directory: {}", output_path.parent_path().string());
  }

  // Legge i parametri dal file di configurazione (path parametrico)
  std::ifstream f(config_file);
  if (!f.is_open()) {
    throw std::runtime_error("Impossibile aprire config: " + config_file.string());
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
  if (!analysisManager->OpenFile(root_output_file)) {
    throw std::runtime_error("Impossibile aprire il file ROOT di output: " + root_output_file);
  }

  // Compressione LZ4 livello 4 — encoding: algoritmo*100 + livello (404 = LZ4 lv4)
  analysisManager->SetCompressionLevel(404);

  // Ntuple 0: configurazioni lenti
  analysisManager->CreateNtuple("Configurations", "Lens configurations");
  analysisManager->CreateNtupleIColumn("config_id");
  analysisManager->CreateNtupleDColumn("x1");
  analysisManager->CreateNtupleDColumn("x2");
  analysisManager->FinishNtuple(0);

  // Ntuple 1: Efficienza geometrica per configurazione
  analysisManager->CreateNtuple("Efficiency", "Geometrical efficiency");
  analysisManager->CreateNtupleIColumn("config_id");
  analysisManager->CreateNtupleIColumn("n_photons");
  analysisManager->CreateNtupleIColumn("n_hits");
  analysisManager->FinishNtuple(1);

  // Parametri di ottimizzazione
  double x_min = config["x_min"];
  double x_max = config["x_max"];
  double dx    = config["dx"];

  double r1 = config["r1"], h1 = config["h1"];
  double r2 = config["r2"], h2 = config["h2"];

  // Parametri sorgente GPS
  double source_width  = config.value("source_width", 60.0);
  double source_height = config.value("source_height", 10.0);
  int n_photons        = config.value("n_photons", 10000);

  auto UImanager                     = G4UImanager::GetUIpointer();
  std::filesystem::path macro_to_run = macro_file.empty() ? "macros/optimization.mac" : macro_file;

  // Configurazione sorgente GPS (Rettangolo con emissione isotropa)
  UImanager->ApplyCommand("/gps/particle opticalphoton");
  UImanager->ApplyCommand("/gps/energy 2.5 eV");
  UImanager->ApplyCommand("/gps/pos/type Plane");
  UImanager->ApplyCommand("/gps/pos/shape Rectangle");
  UImanager->ApplyCommand("/gps/pos/centre 0 " + std::to_string(source_height / 2.0) + " 0 mm");
  UImanager->ApplyCommand("/gps/pos/halfx " + std::to_string(source_height / 2.0) + " mm");
  UImanager->ApplyCommand("/gps/pos/halfy " + std::to_string(source_width / 2.0) + " mm");
  UImanager->ApplyCommand("/gps/pos/rot1 0 1 0"); // local X = global Y
  UImanager->ApplyCommand("/gps/pos/rot2 0 0 1"); // local Y = global Z

  UImanager->ApplyCommand("/gps/ang/type iso");
  UImanager->ApplyCommand("/gps/direction 1 0 0");
  UImanager->ApplyCommand("/gps/ang/mintheta 0 deg");
  UImanager->ApplyCommand("/gps/ang/maxtheta 90 deg");

  int config_counter = 0;
  // Loop geometrico
  for (double x1 = x_min - r1 + h1; x1 <= x_max - h2 - 1 - r1; x1 += dx) {
    double x2_min = x1 + r1 + r2 + 1;
    double x2_max = x_max + r2 - h2;

    for (double x2 = x2_min; x2 <= x2_max; x2 += dx) {
      det->SetLensPositions(x1, x2);
      run_manager->GeometryHasBeenModified();

      // Salva la configurazione delle lenti (Ntuple 0)
      analysisManager->FillNtupleIColumn(0, 0, config_counter);
      analysisManager->FillNtupleDColumn(0, 1, x1);
      analysisManager->FillNtupleDColumn(0, 2, x2);
      analysisManager->AddNtupleRow(0);

      // Imposta l'identificatore della configurazione e resetta il contatore hit
      auto* eventAction = dynamic_cast<EventAction*>(
          const_cast<G4UserEventAction*>(run_manager->GetUserEventAction()));
      if (eventAction) {
        eventAction->SetConfigId(config_counter);
        eventAction->ResetLastRunHitCount();
      }

      // CAMBIA SEED QUI
      long seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
      CLHEP::HepRandom::setTheSeed(seed);

      // Esegue la macro con n_photons
      if (macro_file.empty()) {
        UImanager->ApplyCommand("/run/beamOn " + std::to_string(n_photons));
      } else {
        UImanager->ApplyCommand("/control/execute " + macro_to_run.string());
      }

      // Salva l'efficienza (Ntuple 1)
      int n_hits = eventAction ? eventAction->GetLastRunHitCount() : 0;
      analysisManager->FillNtupleIColumn(1, 0, config_counter);
      analysisManager->FillNtupleIColumn(1, 1, n_photons);
      analysisManager->FillNtupleIColumn(1, 2, n_hits);
      analysisManager->AddNtupleRow(1);

      spdlog::info("Config done: x1={} mm, x2={} mm, config_id={}, hits={}/{}", x1, x2,
                   config_counter, n_hits, n_photons);
      config_counter++;
    }
  }

  // Scrive e chiude il file di output
  analysisManager->Write();
  analysisManager->CloseFile();

  spdlog::info("Optimization completed");
}

} // namespace riptide
