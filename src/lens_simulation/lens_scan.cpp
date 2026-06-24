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
#include "focus_map.hpp"
#include "importance_sampling.hpp"
#include "primary_generator_action.hpp"
#include "spot_grid.hpp"

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
               const std::string& l1_id_arg, const std::string& l2_id_arg,
               const std::string& focus_tsv) {
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
  std::string l1_id = !l1_id_arg.empty() ? l1_id_arg : det->GetL1Id();
  std::string l2_id = !l2_id_arg.empty() ? l2_id_arg : det->GetL2Id();

  // Se sono state passate lenti diverse da quelle attuali nel detector, aggiornale
  if (l1_id != det->GetL1Id() || l2_id != det->GetL2Id()) {
    det->SetLenses(l1_id, l2_id);
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
  analysisManager->CreateNtupleSColumn("l1_id");
  analysisManager->CreateNtupleSColumn("l2_id");
  analysisManager->CreateNtupleDColumn("x_det");      // posizione detector per questa config
  analysisManager->CreateNtupleIColumn("focus_valid"); // 1=valida, 0=fuoco invalido
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
  double dx    = config["dx"];
  double x_max;
  double global_x_det  = 0.0;
  double lens_det_gap  = config.value("lens_det_gap", 0.0);
  bool   use_focus_map = !focus_tsv.empty();

  if (use_focus_map) {
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

  // Parametri sorgente GPS (mappa 3D PSF)
  double source_y_min = config.value("source_y_min", 0.0);
  double source_y_max = config.value("source_y_max", 10.0 * sqrt(2.0));
  double source_dy    = config.value("source_dy", 1.0);

  double source_x_min = config.value("source_x_min", -30.0);
  double source_x_max = config.value("source_x_max", 30.0);
  double source_dx    = config.value("source_dx", 5.0);

  int lens_n_photons   = config.value("lens_n_photons", 1000);
  int config_id_offset = config.value("config_id_offset", 0);
  int run_id_offset    = config.value("run_id_offset", 0);

  int config_counter = config_id_offset;
  int run_counter    = run_id_offset;

  spdlog::info("Simulating lens pair: {} & {}", l1_id, l2_id);

  // Ricalcola h dopo il cambio lenti
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

  auto t_scan_start    = std::chrono::steady_clock::now();
  int  total_pairs     = static_cast<int>(pairs.size());
  int  completed_pairs = 0;
  int  next_progress_pct = 10;
  int  n_src_y = static_cast<int>((source_y_max - source_y_min) / source_dy) + 1;

  auto log_progress = [&](const char* tag) {
    int pct = total_pairs > 0 ? (completed_pairs * 100 / total_pairs) : 0;
    if (pct < next_progress_pct)
      return;
    double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t_scan_start).count();
    double eta = elapsed * (total_pairs - completed_pairs) / completed_pairs;
    spdlog::info("[PROGRESS] {} {}% ({}/{} coppie) — elapsed {:.0f}s ETA {:.0f}s",
                 tag, pct, completed_pairs, total_pairs, elapsed, eta);
    spdlog::default_logger()->flush();
    next_progress_pct = (pct / 10 + 1) * 10;
  };

  bool use_is  = config.value("use_importance_sampling", false);

  auto* eventAction = dynamic_cast<EventAction*>(
      const_cast<G4UserEventAction*>(run_manager->GetUserEventAction()));
  if (!eventAction)
    throw std::runtime_error("EventAction non trovato in lens_scan");

  auto* primaryGen = const_cast<PrimaryGeneratorAction*>(
      static_cast<const PrimaryGeneratorAction*>(run_manager->GetUserPrimaryGeneratorAction()));
  if (!primaryGen)
    throw std::runtime_error("PrimaryGeneratorAction non trovato in lens_scan");

  auto spots  = build_spot_grid(source_x_min, source_x_max, source_dx,
                                source_y_min, source_y_max, source_dy);
  int  n_spots = static_cast<int>(spots.size());
  spdlog::info("[lens] Inizio scan: {} coppie, {} spot/coppia, 1 BeamOn/coppia",
               total_pairs, n_spots);

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
    analysisManager->FillNtupleSColumn(0, 3, l1_id);
    analysisManager->FillNtupleSColumn(0, 4, l2_id);
    analysisManager->FillNtupleDColumn(0, 5, x_det_config);
    analysisManager->FillNtupleIColumn(0, 6, focus_valid);
    analysisManager->AddNtupleRow(0);

    if (focus_valid) {
      // Aggiorna coni IS per questa coppia (dipendono da x1)
      if (use_is) {
        auto params = det->GetL1Params();
        for (auto& sc : spots) {
          G4ThreeVector axis = (G4ThreeVector(params.x, 0, 0) - sc.pos).unit();
          double maxTheta;
          ImportanceSamplingHelper::CalculateCone(sc.pos, params, axis, maxTheta);
          sc.is_axis   = axis;
          sc.is_theta  = maxTheta;
          sc.is_weight = 1.0 - std::cos(maxTheta);
          sc.use_is    = true;
        }
      }

      // Configura EventAction per accumulo per-spot
      eventAction->InitSpotMode(n_spots, source_x_min, source_dx,
                                source_y_min, source_dy, n_src_y);

      // Configura cycling PGA e lancia 1 BeamOn per tutti gli spot
      primaryGen->ConfigureSpotCycling(spots);
      primaryGen->ResetCycling();

      long seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
      CLHEP::HepRandom::setTheSeed(seed);

      // Macro configura GPS (particle type, energy, angular distribution)
      // La posizione viene sovrascritta dal cycling PGA per ogni evento
      UImanager->ApplyCommand("/control/execute " + macro_to_run.string());
      UImanager->ApplyCommand("/run/beamOn "
                              + std::to_string(static_cast<long long>(n_spots) * lens_n_photons));

      primaryGen->DisableCycling();

      // Scrivi n_spots righe nel Runs ntuple (struttura ROOT invariata)
      const auto& all_hits    = eventAction->GetSpotHits();
      const auto& all_weights = eventAction->GetSpotWeights();

      for (int spot_idx = 0; spot_idx < n_spots; ++spot_idx) {
        double x_source      = spots[spot_idx].pos.x();
        double y_source      = spots[spot_idx].pos.y();
        const auto& hits     = all_hits[spot_idx];
        y_hits_vec.clear();
        z_hits_vec.clear();
        y_hits_vec.reserve(hits.size());
        z_hits_vec.reserve(hits.size());
        for (const auto& h : hits) {
          y_hits_vec.push_back(h.y);
          z_hits_vec.push_back(h.z);
        }
        int run_id = run_counter++;
        analysisManager->FillNtupleIColumn(1, 0, config_counter);
        analysisManager->FillNtupleIColumn(1, 1, run_id);
        analysisManager->FillNtupleDColumn(1, 2, x_source);
        analysisManager->FillNtupleDColumn(1, 3, y_source);
        analysisManager->FillNtupleDColumn(1, 4, all_weights[spot_idx]);
        analysisManager->AddNtupleRow(1);
      }

      eventAction->ResetSpotAccumulators();
      spdlog::info("config_id={} completata: {} spot, BeamOn={}", config_counter, n_spots,
                   static_cast<long long>(n_spots) * lens_n_photons);
    } else {
      spdlog::info("Skipped config_id={} (focus_valid=0, x_det={:.1f})", config_counter,
                   x_det_config);
    }

    completed_pairs++;
    log_progress("lens");
    config_counter++;
  }

  analysisManager->Write();
  analysisManager->CloseFile();
  spdlog::info("Simulation completed. Total configs: {}, Total runs: {}", config_counter,
               run_counter);
}

} // namespace riptide
