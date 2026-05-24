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
#include "focus_map.hpp"
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
#include <set>
#include <sstream>

namespace riptide {

void run_optimization(G4RunManager* run_manager, const std::filesystem::path& macro_file,
                      const std::string& root_output_file, const std::filesystem::path& config_file,
                      bool all_lenses, const std::string& l1_id, const std::string& l2_id,
                      const std::string& focus_tsv, const std::string& lens_subset) {
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
  try {
    f >> config;
  } catch (const json::parse_error& e) {
    throw std::runtime_error("Config JSON malformato: " + std::string(e.what()));
  }

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
    std::set<std::string> subset_ids;
    if (!lens_subset.empty()) {
      std::istringstream ss(lens_subset);
      std::string tok;
      while (std::getline(ss, tok, ','))
        if (!tok.empty()) subset_ids.insert(tok);
    }
    for (const auto& l1 : lenses) {
      if (!subset_ids.empty() && !subset_ids.count(l1.id)) continue;
      for (const auto& l2 : lenses) {
        if (!subset_ids.empty() && !subset_ids.count(l2.id)) continue;
        models.push_back({l1.id, l2.id});
      }
    }
    spdlog::warn("Optimization of ALL lens combinations enabled ({} combinations)", models.size());
  } else if (!l1_id.empty() && !l2_id.empty()) {
    models.push_back({l1_id, l2_id});
  } else {
    // Usa i modelli correnti dal detector se disponibili, altrimenti default
    std::string current_id75 = det->GetL1Id();
    std::string current_id60 = det->GetL2Id();
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
  analysisManager->CreateNtupleSColumn("l1_id");
  analysisManager->CreateNtupleSColumn("l2_id");
  analysisManager->CreateNtupleDColumn("x_det");      // posizione detector per questa config
  analysisManager->CreateNtupleIColumn("focus_valid"); // 1=valida, 0=fuoco invalido
  analysisManager->FinishNtuple(0);

  // Ntuple 1: Efficienza geometrica per configurazione
  analysisManager->CreateNtuple("Efficiency", "Geometrical efficiency");
  analysisManager->CreateNtupleIColumn("config_id");
  analysisManager->CreateNtupleIColumn("n_photons");
  analysisManager->CreateNtupleDColumn("n_hits");
  analysisManager->FinishNtuple(1);

  // Parametri di ottimizzazione
  double x_min = config["x_min"];
  double dx    = config["dx"];
  double x_max;
  double global_x_det   = 0.0;
  double lens_det_gap   = config.value("lens_det_gap", 0.0);
  bool   use_focus_map  = !focus_tsv.empty();

  if (use_focus_map) {
    // focus-tsv mode: usa l'intera griglia; validità verificata per-config
    x_max = config.contains("x_max") ? config["x_max"].get<double>()
                                      : config["x_det"].get<double>() - lens_det_gap;
    global_x_det = config.value("x_det", 0.0);
    det->SetDetectorPosition(global_x_det);
    spdlog::info("focus-tsv mode: x_max={:.1f}, TSV={}", x_max, focus_tsv);
  } else if (config.contains("x_det")) {
    double x_det = config["x_det"].get<double>();
    x_max        = x_det - lens_det_gap;
    global_x_det = x_det;
    spdlog::info("x_max = x_det({:.1f}) - lens_det_gap({:.1f}) = {:.1f} mm", x_det, lens_det_gap,
                 x_max);
    det->SetDetectorPosition(x_det);
  } else {
    x_max        = config["x_max"];
    global_x_det = config.value("x_det", 0.0);
  }

  // Carica mappa fuochi se richiesto
  FocusMap focus_map;
  if (use_focus_map)
    focus_map = load_focus_map(focus_tsv);

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

    double h1 = det->GetL1Thickness();
    double h2 = det->GetL2Thickness();

    std::vector<std::pair<double, double>> pairs;
    if (config.contains("pairs")) {
      for (const auto& p : config["pairs"]) {
        pairs.push_back({p[0].get<double>(), p[1].get<double>()});
      }
    } else if (use_focus_map) {
      // Coppie estratte direttamente dal TSV per garantire allineamento
      double x1_start = config.value("x1_start", x_min);
      double x1_end   = config.value("x1_end", x_max);
      for (auto& p : get_pairs_from_focus_map(focus_map)) {
        if (p.first >= x1_start - 1e-6 && p.first <= x1_end + 1e-6 &&
            p.second <= x_max + 1e-6)
          pairs.push_back(p);
      }
      spdlog::info("focus-tsv: {} coppie caricate dal TSV", pairs.size());
    } else {
      const double margin = config.value("lens_gap_margin", 1.0);
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

      // Determina posizione detector e validità per questa config
      double x_det_config = global_x_det;
      int    focus_valid  = 1;
      if (use_focus_map) {
        auto fopt = lookup_focus(focus_map, x1, x2);
        if (!fopt) {
          focus_valid  = 0;
          x_det_config = 0.0;
          spdlog::warn("focus-tsv: nessuna corrispondenza per ({:.1f},{:.1f}) — config invalida",
                       x1, x2);
        } else {
          x_det_config = *fopt;
          if (std::isnan(x_det_config) || x_det_config <= 0.0) {
            focus_valid = 0;
            spdlog::warn("focus-tsv: fuoco {:.1f} invalido (NaN/non positivo)", x_det_config);
          } else {
            if (x_det_config < x2 + lens_det_gap) {
              spdlog::warn(
                  "focus-tsv: fuoco {:.1f} prima del limite {:.1f}, clamped a {:.1f}",
                  x_det_config, x2 + lens_det_gap, x2 + lens_det_gap);
              x_det_config = x2 + lens_det_gap;
            }
            det->SetDetectorPosition(x_det_config);
          }
        }
      }

      det->SetLensPositions(x1, x2);
      run_manager->GeometryHasBeenModified();

      analysisManager->FillNtupleIColumn(0, 0, config_counter);
      analysisManager->FillNtupleDColumn(0, 1, x1);
      analysisManager->FillNtupleDColumn(0, 2, x2);
      analysisManager->FillNtupleSColumn(0, 3, model.id75);
      analysisManager->FillNtupleSColumn(0, 4, model.id60);
      analysisManager->FillNtupleDColumn(0, 5, x_det_config);
      analysisManager->FillNtupleIColumn(0, 6, focus_valid);
      analysisManager->AddNtupleRow(0);

      if (focus_valid) {
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
          ImportanceSamplingHelper::CalculateGlobalCone(corners, det->GetL1Params(), axis,
                                                        maxTheta);

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
      } else {
        // Config invalida: riga sentinel in Efficiency (n_photons=0, n_hits=0)
        analysisManager->FillNtupleIColumn(1, 0, config_counter);
        analysisManager->FillNtupleIColumn(1, 1, 0);
        analysisManager->FillNtupleDColumn(1, 2, 0.0);
        analysisManager->AddNtupleRow(1);

        spdlog::info("Skipped config_id={} (focus_valid=0, x_det={:.1f})", config_counter,
                     x_det_config);
      }

      config_counter++;
    }
  }

  analysisManager->Write();
  analysisManager->CloseFile();
  spdlog::info("Optimization completed. Total configs: {}", config_counter);
}

} // namespace riptide
