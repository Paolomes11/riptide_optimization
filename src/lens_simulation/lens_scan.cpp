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
#include "importance_sampling.hpp"
#include "primary_generator_action.hpp"

#include <G4AnalysisManager.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4UImanager.hh>

#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace riptide {

void lens_scan(G4RunManager* run_manager, const std::filesystem::path& macro_file,
               const std::string& root_output_file, const std::filesystem::path& config_file,
               const std::string& lens75_id_arg, const std::string& lens60_id_arg) {
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
  try {
    f >> config;
  } catch (const json::parse_error& e) {
    throw std::runtime_error("Config JSON malformato: " + std::string(e.what()));
  }

  auto det = const_cast<DetectorConstruction*>(
      dynamic_cast<const DetectorConstruction*>(run_manager->GetUserDetectorConstruction()));
  if (!det)
    throw std::runtime_error("DetectorConstruction not found!");

  auto UImanager = G4UImanager::GetUIpointer();
  std::filesystem::path macro_to_run =
      macro_file.empty() ? "macros/lens_simulation.mac" : macro_file;
  if (!std::filesystem::exists(macro_to_run))
    throw std::runtime_error("Macro non trovata: " + macro_to_run.string());

  // Ottieni i modelli di lenti correnti dal detector o dagli argomenti
  std::string lens75_id = !lens75_id_arg.empty() ? lens75_id_arg : det->GetLens75Id();
  std::string lens60_id = !lens60_id_arg.empty() ? lens60_id_arg : det->GetLens60Id();

  // Se sono state passate lenti diverse da quelle attuali nel detector, aggiornale
  if (lens75_id != det->GetLens75Id() || lens60_id != det->GetLens60Id()) {
    det->SetLenses(lens75_id, lens60_id);
  }

  // Apertura file ROOT e creazione delle Ntuple
  auto analysisManager = G4AnalysisManager::Instance();
  if (!analysisManager->OpenFile(root_output_file))
    throw std::runtime_error("lens_scan: impossibile aprire il file ROOT di output: " +
                             root_output_file);

  // Compressione LZ4 livello 4 — encoding: algoritmo*100 + livello (404 = LZ4 lv4)
  analysisManager->SetCompressionLevel(404);

  // Ntuple 0: configurazioni lenti
  analysisManager->CreateNtuple("Configurations", "Lens configurations");
  analysisManager->CreateNtupleIColumn("config_id");
  analysisManager->CreateNtupleDColumn("x1");
  analysisManager->CreateNtupleDColumn("x2");
  analysisManager->CreateNtupleSColumn("lens75_id");
  analysisManager->CreateNtupleSColumn("lens60_id");
  analysisManager->FinishNtuple(0);

  // Vettori bound per le hit (necessari per CreateNtupleFColumn con std::vector)
  std::vector<float> y_hits_vec;
  std::vector<float> z_hits_vec;

  // Ntuple 1: posizione sorgente / run + Hits (vettori)
  analysisManager->CreateNtuple("Runs", "PSF data");
  analysisManager->CreateNtupleIColumn("config_id");
  analysisManager->CreateNtupleIColumn("run_id");
  analysisManager->CreateNtupleDColumn("x_source");
  analysisManager->CreateNtupleDColumn("y_source");
  analysisManager->CreateNtupleDColumn("n_hits");
  analysisManager->CreateNtupleFColumn("y_hits", y_hits_vec);
  analysisManager->CreateNtupleFColumn("z_hits", z_hits_vec);
  analysisManager->FinishNtuple(1);

  // Parametri di scansione lenti e sorgente
  double x_min = config["x_min"];
  double x_max = config["x_max"];
  double dx    = config["dx"];

  // Parametri sorgente GPS (mappa 3D PSF)
  double source_y_min = config.value("source_y_min", 0.0);
  double source_y_max = config.value("source_y_max", 10.0 * sqrt(2.0));
  double source_dy    = config.value("source_dy", 1.0);

  double source_x_min = config.value("source_x_min", -30.0);
  double source_x_max = config.value("source_x_max", 30.0);
  double source_dx    = config.value("source_dx", 5.0);

  int config_id_offset = config.value("config_id_offset", 0);
  int run_id_offset    = config.value("run_id_offset", 0);

  int config_counter = config_id_offset;
  int run_counter    = run_id_offset;

  spdlog::info("Simulating lens pair: {} & {}", lens75_id, lens60_id);

  // Ricalcola h dopo il cambio lenti
  double h1 = det->GetLens75Thickness();
  double h2 = det->GetLens60Thickness();

  std::vector<std::pair<double, double>> pairs;
  if (config.contains("pairs")) {
    for (const auto& p : config["pairs"]) {
      pairs.push_back({p[0].get<double>(), p[1].get<double>()});
    }
  } else {
    const double margin = config.value("lens_gap_margin", 1.0);
    double x1_start     = config.value("x1_start", x_min);
    double x1_end       = config.value("x1_end", x_max);

    for (double x1 = x1_start; x1 <= x1_end + 1e-9; x1 += dx) {
      // Ora che SetLensPositions lavora con i centri geometrici,
      // la collisione è data semplicemente dalla somma dei semi-spessori
      double x2_min_collision = x1 + (h1 + h2) / 2.0 + margin;
      double x2_start         = std::max(x1 + dx, x2_min_collision);
      for (double x2 = x2_start; x2 <= x_max + 1e-9; x2 += dx) {
        pairs.push_back({x1, x2});
      }
    }
  }

  for (const auto& pair : pairs) {
    double x1 = pair.first;
    double x2 = pair.second;

    det->SetLensPositions(x1, x2);
    run_manager->GeometryHasBeenModified();

    analysisManager->FillNtupleIColumn(0, 0, config_counter);
    analysisManager->FillNtupleDColumn(0, 1, x1);
    analysisManager->FillNtupleDColumn(0, 2, x2);
    analysisManager->FillNtupleSColumn(0, 3, lens75_id);
    analysisManager->FillNtupleSColumn(0, 4, lens60_id);
    analysisManager->AddNtupleRow(0);

    for (double x_source = source_x_min; x_source <= source_x_max + 1e-9; x_source += source_dx) {
      for (double y_source = source_y_min; y_source <= source_y_max + 1e-9; y_source += source_dy) {
        UImanager->ApplyCommand("/gps/pos/centre " + std::to_string(x_source) + " "
                                + std::to_string(y_source) + " 0 mm");

        // Calcola il cono statico per questo run (sorgente puntiforme)
        if (config.value("use_importance_sampling", false)) {
          G4ThreeVector pos(x_source, y_source, 0);
          G4ThreeVector axis;
          double maxTheta;
          auto params = det->GetLens75Params();
          axis        = (G4ThreeVector(params.x, 0, 0) - pos).unit();
          ImportanceSamplingHelper::CalculateCone(pos, params, axis, maxTheta);

          // Angolo solido emisferico originale (piatto)
          double originalOmega = 2.0 * M_PI;
          double cosMaxTheta   = std::cos(maxTheta);
          double newOmega      = 2.0 * M_PI * (1.0 - cosMaxTheta);
          double weight        = newOmega / originalOmega;

          if (run_counter == 0) {
            spdlog::info("Configured Importance Sampling: theta_max = {:.2f} deg, weight = {:.4f}",
                         maxTheta * 180.0 / M_PI, weight);
          }

          auto* primaryGen =
              const_cast<PrimaryGeneratorAction*>(static_cast<const PrimaryGeneratorAction*>(
                  run_manager->GetUserPrimaryGeneratorAction()));
          if (primaryGen) {
            primaryGen->SetStaticCone(axis, maxTheta);
          }
        }

        int run_id        = run_counter++;
        auto* eventAction = dynamic_cast<EventAction*>(
            const_cast<G4UserEventAction*>(run_manager->GetUserEventAction()));
        if (eventAction) {
          eventAction->runID = run_id;
          eventAction->ResetLastRunHitCount();
          eventAction->ClearEventHits();
        }

        long seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        CLHEP::HepRandom::setTheSeed(seed);

        UImanager->ApplyCommand("/control/execute " + macro_to_run.string());

        if (eventAction) {
          const auto& hits = eventAction->GetEventHits();
          y_hits_vec.clear();
          z_hits_vec.clear();
          y_hits_vec.reserve(hits.size());
          z_hits_vec.reserve(hits.size());
          for (const auto& h : hits) {
            y_hits_vec.push_back(h.y);
            z_hits_vec.push_back(h.z);
          }

          double n_hits = eventAction ? eventAction->GetLastRunHitCount() : 0.0;
          analysisManager->FillNtupleIColumn(1, 0, config_counter);
          analysisManager->FillNtupleIColumn(1, 1, run_id);
          analysisManager->FillNtupleDColumn(1, 2, x_source);
          analysisManager->FillNtupleDColumn(1, 3, y_source);
          analysisManager->FillNtupleDColumn(1, 4, n_hits);
          analysisManager->AddNtupleRow(1);
        }

        spdlog::info("Run done: config_id={}, x_src={}, y_src={}, n_hits={}", config_counter,
                     x_source, y_source, eventAction ? eventAction->GetLastRunHitCount() : 0);
      }
    }
    config_counter++;
  }

  analysisManager->Write();
  analysisManager->CloseFile();
  spdlog::info("Simulation completed. Total configs: {}, Total runs: {}", config_counter,
               run_counter);
}

} // namespace riptide
