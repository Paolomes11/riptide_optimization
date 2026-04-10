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

void run_optimization(G4RunManager* run_manager, const std::filesystem::path& macro_file,
                      const std::string& root_output_file, const std::filesystem::path& config_file,
                      bool all_lenses, const std::string& lens75_id, const std::string& lens60_id) {
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

  // Carica i modelli di lenti se richiesto
  struct LensModel {
    std::string id75;
    std::string id60;
  };
  std::vector<LensModel> models;

  if (all_lenses) {
    LensCutter cutter("lens_cutter/lens_data/thorlabs_biconvex.tsv");
    const auto& lenses = cutter.get_lenses();
    for (const auto& l1 : lenses) {
      for (const auto& l2 : lenses) {
        models.push_back({l1.id, l2.id});
      }
    }
    spdlog::warn("Optimization of ALL lens combinations enabled ({} combinations)", models.size());
  } else if (!lens75_id.empty() && !lens60_id.empty()) {
    models.push_back({lens75_id, lens60_id});
  } else {
    // Usa i modelli correnti dal detector se disponibili, altrimenti default
    std::string current_id75 = det->GetLens75Id();
    std::string current_id60 = det->GetLens60Id();
    models.push_back({current_id75, current_id60});
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
  analysisManager->CreateNtupleSColumn("lens75_id");
  analysisManager->CreateNtupleSColumn("lens60_id");
  analysisManager->FinishNtuple(0);

  // Ntuple 1: Efficienza geometrica per configurazione
  analysisManager->CreateNtuple("Efficiency", "Geometrical efficiency");
  analysisManager->CreateNtupleIColumn("config_id");
  analysisManager->CreateNtupleIColumn("n_photons");
  analysisManager->CreateNtupleDColumn("n_hits");
  analysisManager->FinishNtuple(1);

  // Parametri di ottimizzazione
  double x_min = config["x_min"];
  double x_max = config["x_max"];
  double dx    = config["dx"];

  // Parametri sorgente GPS
  double source_width     = config.value("source_width", 60.0);
  double source_height    = config.value("source_height", 10.0);
  double source_thickness = config.value("source_thickness", 30.0);
  int n_photons           = config.value("n_photons", 10000);

  auto UImanager                     = G4UImanager::GetUIpointer();
  std::filesystem::path macro_to_run = macro_file.empty() ? "macros/optimization.mac" : macro_file;

  // Configurazione sorgente GPS
  UImanager->ApplyCommand("/gps/particle opticalphoton");
  UImanager->ApplyCommand("/gps/energy 2.5 eV");
  UImanager->ApplyCommand("/gps/pos/type Plane");
  UImanager->ApplyCommand("/gps/pos/shape Rectangle");
  UImanager->ApplyCommand("/gps/pos/centre 0 " + std::to_string(source_height) + " 0 mm");
  UImanager->ApplyCommand("/gps/pos/halfx " + std::to_string(source_thickness / 2.0) + " mm");
  UImanager->ApplyCommand("/gps/pos/halfy " + std::to_string(source_width / 2.0) + " mm");
  UImanager->ApplyCommand("/gps/pos/rot1 1 0 0");
  UImanager->ApplyCommand("/gps/pos/rot2 0 0 1");

  UImanager->ApplyCommand("/gps/ang/type iso");
  UImanager->ApplyCommand("/gps/ang/rot1 0 -1 0");
  UImanager->ApplyCommand("/gps/ang/rot2 0 0 1");
  UImanager->ApplyCommand("/gps/ang/mintheta 0 deg");
  UImanager->ApplyCommand("/gps/ang/maxtheta 90 deg");

  int config_id_offset = config.value("config_id_offset", 0);
  int config_counter   = config_id_offset;

  for (const auto& model : models) {
    spdlog::info("Optimizing lens pair: {} & {}", model.id75, model.id60);
    det->SetLenses(model.id75, model.id60);

    double h1 = det->GetLens75Thickness();
    double h2 = det->GetLens60Thickness();

    std::vector<std::pair<double, double>> pairs;
    if (config.contains("pairs")) {
      for (const auto& p : config["pairs"]) {
        pairs.push_back({p[0].get<double>(), p[1].get<double>()});
      }
    } else {
      const double margin = 1.0;
      for (double x1 = x_min; x1 <= x_max; x1 += dx) {
        // Ora che SetLensPositions lavora con i centri geometrici,
        // la collisione è data semplicemente dalla somma dei semi-spessori
        double x2_min_collision = x1 + (h1 + h2) / 2.0 + margin;
        double x2_start         = std::max(x1 + dx, x2_min_collision);
        for (double x2 = x2_start; x2 <= x_max; x2 += dx) {
          pairs.push_back({x1, x2});
        }
      }
    }

    for (const auto& pair : pairs) {
      double x1 = pair.first;
      double x2 = pair.second;

      det->SetLensPositions(x1, x2);
      run_manager->GeometryHasBeenModified();

      // Calcola il cono globale per la sorgente estesa
      if (config.value("use_importance_sampling", false)) {
        double sw                          = source_width;
        double st                          = source_thickness;
        double sh                          = source_height;
        std::vector<G4ThreeVector> corners = {
            G4ThreeVector(-st / 2, sh, -sw / 2), G4ThreeVector(st / 2, sh, -sw / 2),
            G4ThreeVector(-st / 2, sh, sw / 2), G4ThreeVector(st / 2, sh, sw / 2)};
        G4ThreeVector axis;
        double maxTheta;
        ImportanceSamplingHelper::CalculateGlobalCone(corners, det->GetLens75Params(), axis,
                                                      maxTheta);

        // Angolo solido emisferico originale (piatto)
        double originalOmega = 2.0 * M_PI;
        double cosMaxTheta   = std::cos(maxTheta);
        double newOmega      = 2.0 * M_PI * (1.0 - cosMaxTheta);
        double weight        = newOmega / originalOmega;

        spdlog::info("Configured Importance Sampling: theta_max = {:.2f} deg, weight = {:.4f}",
                     maxTheta * 180.0 / M_PI, weight);

        auto* primaryGen =
            const_cast<PrimaryGeneratorAction*>(static_cast<const PrimaryGeneratorAction*>(
                run_manager->GetUserPrimaryGeneratorAction()));
        if (primaryGen) {
          primaryGen->SetStaticCone(axis, maxTheta);
        }
      }

      analysisManager->FillNtupleIColumn(0, 0, config_counter);
      analysisManager->FillNtupleDColumn(0, 1, x1);
      analysisManager->FillNtupleDColumn(0, 2, x2);
      analysisManager->FillNtupleSColumn(0, 3, model.id75);
      analysisManager->FillNtupleSColumn(0, 4, model.id60);
      analysisManager->AddNtupleRow(0);

      auto* eventAction = dynamic_cast<EventAction*>(
          const_cast<G4UserEventAction*>(run_manager->GetUserEventAction()));
      if (eventAction) {
        eventAction->SetConfigId(config_counter);
        eventAction->ResetLastRunHitCount();
      }

      long seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
      CLHEP::HepRandom::setTheSeed(seed);

      if (macro_file.empty()) {
        UImanager->ApplyCommand("/run/beamOn " + std::to_string(n_photons));
      } else {
        UImanager->ApplyCommand("/control/execute " + macro_to_run.string());
      }

      double n_hits = eventAction ? eventAction->GetLastRunHitCount() : 0.0;
      analysisManager->FillNtupleIColumn(1, 0, config_counter);
      analysisManager->FillNtupleIColumn(1, 1, n_photons);
      analysisManager->FillNtupleDColumn(1, 2, n_hits);
      analysisManager->AddNtupleRow(1);

      spdlog::info("Run done: config_id={}, hits={}/{}", config_counter, n_hits, n_photons);
      config_counter++;
    }
  }

  analysisManager->Write();
  analysisManager->CloseFile();
  spdlog::info("Optimization completed. Total configs: {}", config_counter);
}

} // namespace riptide
