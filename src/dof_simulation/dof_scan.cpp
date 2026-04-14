#include "dof_scan.hpp"

#include "dof_event_action.hpp"
#include "dof_stepping_action.hpp"
#include "importance_sampling.hpp"
#include "primary_generator_action.hpp"

#include "optimization/detector_construction.hpp"

#include <G4AnalysisManager.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4UImanager.hh>

#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>

#include <CLHEP/Random/Random.h>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace riptide {

void run_dof_scan(G4RunManager* run_manager, const std::filesystem::path& macro_file,
                  const std::string& root_output_file, const std::filesystem::path& config_file,
                  bool all_lenses, const std::string& lens75_id, const std::string& lens60_id) {
  using json = nlohmann::json;

  std::filesystem::path output_path(root_output_file);
  if (output_path.has_parent_path()) {
    std::filesystem::create_directories(output_path.parent_path());
    spdlog::info("Output directory: {}", output_path.parent_path().string());
  }

  std::ifstream f(config_file);
  if (!f.is_open()) {
    throw std::runtime_error("Impossibile aprire config: " + config_file.string());
  }
  json config;
  f >> config;

  auto det = const_cast<DetectorConstruction*>(
      dynamic_cast<const DetectorConstruction*>(run_manager->GetUserDetectorConstruction()));
  if (!det) {
    throw std::runtime_error("DetectorConstruction not found!");
  }

  auto* eventAction = dynamic_cast<DofEventAction*>(
      const_cast<G4UserEventAction*>(run_manager->GetUserEventAction()));
  if (!eventAction) {
    throw std::runtime_error("DofEventAction not found!");
  }

  auto* steppingAction = dynamic_cast<DofSteppingAction*>(
      const_cast<G4UserSteppingAction*>(run_manager->GetUserSteppingAction()));
  if (!steppingAction) {
    throw std::runtime_error("DofSteppingAction not found!");
  }

  auto* primaryGen = dynamic_cast<PrimaryGeneratorAction*>(
      const_cast<G4VUserPrimaryGeneratorAction*>(run_manager->GetUserPrimaryGeneratorAction()));
  if (!primaryGen) {
    throw std::runtime_error("PrimaryGeneratorAction not found!");
  }

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
    spdlog::warn("DoF scan of ALL lens combinations enabled ({} combinations)", models.size());
  } else if (!lens75_id.empty() && !lens60_id.empty()) {
    models.push_back({lens75_id, lens60_id});
  } else {
    models.push_back({det->GetLens75Id(), det->GetLens60Id()});
  }

  auto analysisManager = G4AnalysisManager::Instance();
  if (!analysisManager->OpenFile(root_output_file)) {
    throw std::runtime_error("Impossibile aprire il file ROOT di output: " + root_output_file);
  }
  analysisManager->SetCompressionLevel(404);

  analysisManager->CreateNtuple("FocalConfigurations", "DoF focal configurations");
  analysisManager->CreateNtupleIColumn("config_id");
  analysisManager->CreateNtupleDColumn("x1");
  analysisManager->CreateNtupleDColumn("x2");
  analysisManager->CreateNtupleDColumn("x_virtual");
  analysisManager->CreateNtupleSColumn("lens75_id");
  analysisManager->CreateNtupleSColumn("lens60_id");
  analysisManager->FinishNtuple(0);

  analysisManager->CreateNtuple("FocalRays", "Rays at virtual plane");
  analysisManager->CreateNtupleIColumn("config_id");
  analysisManager->CreateNtupleDColumn("n_rays");
  analysisManager->CreateNtupleFColumn("y_hits", eventAction->YHits());
  analysisManager->CreateNtupleFColumn("z_hits", eventAction->ZHits());
  analysisManager->CreateNtupleFColumn("dy_hits", eventAction->DyHits());
  analysisManager->CreateNtupleFColumn("dz_hits", eventAction->DzHits());
  analysisManager->CreateNtupleFColumn("weight_hits", eventAction->WeightHits());
  analysisManager->CreateNtupleFColumn("y_source_hits", eventAction->YSourceHits());
  analysisManager->FinishNtuple(1);

  double x_min = config["x_min"];
  double x_max = config["x_max"];
  double dx    = config["dx"];

  int n_photons = config.value("dof_n_photons", 10000);

  bool has_virtual_mm     = config.contains("dof_x_virtual_mm");
  double x_virtual_mm     = config.value("dof_x_virtual_mm", 180.0);
  double x_virtual_offset = config.value("dof_x_virtual_offset", 30.0);

  double source_halfy    = config.value("dof_source_halfy", 5.0);
  double source_halfz    = config.value("dof_source_halfz", 30.0);
  double source_y_centre = config.value("dof_source_y_centre", 5.0);

  auto UImanager = G4UImanager::GetUIpointer();

  UImanager->ApplyCommand("/gps/particle opticalphoton");
  UImanager->ApplyCommand("/gps/energy 2.5 eV");
  UImanager->ApplyCommand("/gps/pos/type Plane");
  UImanager->ApplyCommand("/gps/pos/shape Rectangle");
  UImanager->ApplyCommand("/gps/pos/centre 0 " + std::to_string(source_y_centre) + " 0 mm");
  UImanager->ApplyCommand("/gps/pos/halfx " + std::to_string(source_halfz) + " mm");
  UImanager->ApplyCommand("/gps/pos/halfy " + std::to_string(source_halfy) + " mm");
  UImanager->ApplyCommand("/gps/pos/rot1 0 0 1");
  UImanager->ApplyCommand("/gps/pos/rot2 0 1 0");

  UImanager->ApplyCommand("/gps/ang/type iso");
  UImanager->ApplyCommand("/gps/ang/rot1 0 -1 0");
  UImanager->ApplyCommand("/gps/ang/rot2 0 0 1");
  UImanager->ApplyCommand("/gps/ang/mintheta 0 deg");
  UImanager->ApplyCommand("/gps/ang/maxtheta 90 deg");

  std::vector<G4ThreeVector> source_points;
  source_points.reserve(4);
  source_points.push_back(G4ThreeVector(
      0.0 * CLHEP::mm, (source_y_centre - source_halfy) * CLHEP::mm, (-source_halfz) * CLHEP::mm));
  source_points.push_back(G4ThreeVector(
      0.0 * CLHEP::mm, (source_y_centre - source_halfy) * CLHEP::mm, (+source_halfz) * CLHEP::mm));
  source_points.push_back(G4ThreeVector(
      0.0 * CLHEP::mm, (source_y_centre + source_halfy) * CLHEP::mm, (-source_halfz) * CLHEP::mm));
  source_points.push_back(G4ThreeVector(
      0.0 * CLHEP::mm, (source_y_centre + source_halfy) * CLHEP::mm, (+source_halfz) * CLHEP::mm));

  int config_id_offset = config.value("config_id_offset", 0);
  int config_counter   = config_id_offset;

  for (const auto& model : models) {
    spdlog::info("DoF scan lens pair: {} & {}", model.id75, model.id60);
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

      double x_virtual = has_virtual_mm ? x_virtual_mm : (x2 + x_virtual_offset);
      steppingAction->SetVirtualPlane(x_virtual);

      if (config.value("use_importance_sampling", false)) {
        G4ThreeVector axis;
        double maxTheta = 0.0;
        ImportanceSamplingHelper::CalculateGlobalCone(source_points, det->GetLens75Params(), axis,
                                                      maxTheta);
        primaryGen->SetStaticCone(axis, maxTheta);
      }

      analysisManager->FillNtupleIColumn(0, 0, config_counter);
      analysisManager->FillNtupleDColumn(0, 1, x1);
      analysisManager->FillNtupleDColumn(0, 2, x2);
      analysisManager->FillNtupleDColumn(0, 3, x_virtual);
      analysisManager->FillNtupleSColumn(0, 4, model.id75);
      analysisManager->FillNtupleSColumn(0, 5, model.id60);
      analysisManager->AddNtupleRow(0);

      eventAction->SetConfigId(config_counter);
      eventAction->ClearRays();

      long seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
      CLHEP::HepRandom::setTheSeed(seed);

      if (macro_file.empty()) {
        UImanager->ApplyCommand("/run/beamOn " + std::to_string(n_photons));
      } else {
        UImanager->ApplyCommand("/control/execute " + macro_file.string());
      }

      analysisManager->FillNtupleIColumn(1, 0, config_counter);
      analysisManager->FillNtupleDColumn(1, 1, eventAction->GetWeightedRayCount());
      analysisManager->AddNtupleRow(1);

      spdlog::info("Run done: config_id={}, rays_w={:.1f}, photons={}", config_counter,
                   eventAction->GetWeightedRayCount(), n_photons);

      eventAction->ClearRays();
      config_counter++;
    }
  }

  analysisManager->Write();
  analysisManager->CloseFile();
  spdlog::info("DoF scan completed. Total configs: {}", config_counter);
}

} // namespace riptide
