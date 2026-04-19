#include "psf_dof_scan.hpp"

#include "importance_sampling.hpp"
#include "primary_generator_action.hpp"
#include "psf_dof_event_action.hpp"
#include "psf_dof_stepping_action.hpp"

#include "optimization/detector_construction.hpp"

#include <G4AnalysisManager.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4UImanager.hh>

#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>

#include <CLHEP/Random/Random.h>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace riptide {

struct PsfDofStats {
  double n_hits    = 0.0;
  double mu_y      = 0.0;
  double mu_z      = 0.0;
  double sigma_y   = 0.0;
  double sigma_z   = 0.0;
  double sigma_yz  = 0.0;
  double mu_dy     = 0.0;
  double mu_dz     = 0.0;
  double sigma_dy  = 0.0;
  double sigma_dz  = 0.0;
  double cov_y_dy  = 0.0;
  double cov_z_dz  = 0.0;
  bool has_moments = false;
};

static PsfDofStats compute_stats(const PsfDofMoments& m) {
  PsfDofStats out;
  if (!(m.sum_w > 0.0)) {
    return out;
  }
  out.has_moments = true;
  out.n_hits      = m.sum_w;
  out.mu_y        = m.sum_y / m.sum_w;
  out.mu_z        = m.sum_z / m.sum_w;
  out.mu_dy       = m.sum_dy / m.sum_w;
  out.mu_dz       = m.sum_dz / m.sum_w;

  double e_yy    = m.sum_yy / m.sum_w;
  double e_zz    = m.sum_zz / m.sum_w;
  double e_yz    = m.sum_yz / m.sum_w;
  double e_dy_dy = m.sum_dy_dy / m.sum_w;
  double e_dz_dz = m.sum_dz_dz / m.sum_w;
  double e_y_dy  = m.sum_y_dy / m.sum_w;
  double e_z_dz  = m.sum_z_dz / m.sum_w;

  double var_y  = e_yy - out.mu_y * out.mu_y;
  double var_z  = e_zz - out.mu_z * out.mu_z;
  double cov_yz = e_yz - out.mu_y * out.mu_z;

  double var_dy = e_dy_dy - out.mu_dy * out.mu_dy;
  double var_dz = e_dz_dz - out.mu_dz * out.mu_dz;

  out.sigma_y  = std::sqrt(std::max(0.0, var_y));
  out.sigma_z  = std::sqrt(std::max(0.0, var_z));
  out.sigma_yz = cov_yz;

  out.sigma_dy = std::sqrt(std::max(0.0, var_dy));
  out.sigma_dz = std::sqrt(std::max(0.0, var_dz));

  out.cov_y_dy = e_y_dy - out.mu_y * out.mu_dy;
  out.cov_z_dz = e_z_dz - out.mu_z * out.mu_dz;

  return out;
}

void run_psf_dof_scan(G4RunManager* run_manager, const std::filesystem::path& macro_file,
                      const std::string& root_output_file, const std::filesystem::path& config_file,
                      const std::string& lens75_id_arg, const std::string& lens60_id_arg) {
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

  auto* eventAction = dynamic_cast<PsfDofEventAction*>(
      const_cast<G4UserEventAction*>(run_manager->GetUserEventAction()));
  if (!eventAction) {
    throw std::runtime_error("PsfDofEventAction not found!");
  }

  auto* steppingAction = dynamic_cast<PsfDofSteppingAction*>(
      const_cast<G4UserSteppingAction*>(run_manager->GetUserSteppingAction()));
  if (!steppingAction) {
    throw std::runtime_error("PsfDofSteppingAction not found!");
  }

  auto* primaryGen = dynamic_cast<PrimaryGeneratorAction*>(
      const_cast<G4VUserPrimaryGeneratorAction*>(run_manager->GetUserPrimaryGeneratorAction()));
  if (!primaryGen) {
    throw std::runtime_error("PrimaryGeneratorAction not found!");
  }

  std::string lens75_id = !lens75_id_arg.empty() ? lens75_id_arg : det->GetLens75Id();
  std::string lens60_id = !lens60_id_arg.empty() ? lens60_id_arg : det->GetLens60Id();
  det->SetLenses(lens75_id, lens60_id);

  auto analysisManager = G4AnalysisManager::Instance();
  if (!analysisManager->OpenFile(root_output_file)) {
    throw std::runtime_error("Impossibile aprire il file ROOT di output: " + root_output_file);
  }
  analysisManager->SetCompressionLevel(404);

  bool save_hits     = config.value("psf_dof_save_hits", false);
  int n_photons      = config.value("psf_dof_n_photons", 10000);
  double x_virtual   = config.value("psf_dof_x_virtual_mm", 180.0);
  int config_id_off  = config.value("config_id_offset", 0);
  int run_id_off     = config.value("run_id_offset", 0);
  int config_counter = config_id_off;
  int run_counter    = run_id_off;

  eventAction->SetSaveHits(save_hits);

  analysisManager->CreateNtuple("PsfDofConfigs", "PSF+DoF configurations");
  analysisManager->CreateNtupleIColumn("config_id");
  analysisManager->CreateNtupleDColumn("x1");
  analysisManager->CreateNtupleDColumn("x2");
  analysisManager->CreateNtupleDColumn("x_virtual");
  analysisManager->CreateNtupleSColumn("lens75_id");
  analysisManager->CreateNtupleSColumn("lens60_id");
  analysisManager->FinishNtuple(0);

  analysisManager->CreateNtuple("PsfDofRuns", "PSF+DoF runs");
  analysisManager->CreateNtupleIColumn("config_id");
  analysisManager->CreateNtupleIColumn("run_id");
  analysisManager->CreateNtupleDColumn("x_source");
  analysisManager->CreateNtupleDColumn("y_source");
  analysisManager->CreateNtupleDColumn("n_hits");
  analysisManager->CreateNtupleDColumn("mu_y");
  analysisManager->CreateNtupleDColumn("mu_z");
  analysisManager->CreateNtupleDColumn("sigma_y");
  analysisManager->CreateNtupleDColumn("sigma_z");
  analysisManager->CreateNtupleDColumn("sigma_yz");
  analysisManager->CreateNtupleDColumn("mu_dy");
  analysisManager->CreateNtupleDColumn("mu_dz");
  analysisManager->CreateNtupleDColumn("sigma_dy");
  analysisManager->CreateNtupleDColumn("sigma_dz");
  analysisManager->CreateNtupleDColumn("cov_y_dy");
  analysisManager->CreateNtupleDColumn("cov_z_dz");
  if (save_hits) {
    analysisManager->CreateNtupleFColumn("y_hits", eventAction->YHits());
    analysisManager->CreateNtupleFColumn("z_hits", eventAction->ZHits());
    analysisManager->CreateNtupleFColumn("dy_hits", eventAction->DyHits());
    analysisManager->CreateNtupleFColumn("dz_hits", eventAction->DzHits());
    analysisManager->CreateNtupleFColumn("w_hits", eventAction->WeightHits());
  }
  analysisManager->FinishNtuple(1);

  double x_min = config["x_min"];
  double x_max = config["x_max"];
  double dx    = config["dx"];

  double source_y_min = config.value("source_y_min", 0.0);
  double source_y_max = config.value("source_y_max", 10.0 * std::sqrt(2.0));
  double source_dy    = config.value("source_dy", 1.0);

  double source_x_min = config.value("source_x_min", -30.0);
  double source_x_max = config.value("source_x_max", 30.0);
  double source_dx    = config.value("source_dx", 5.0);

  auto UImanager = G4UImanager::GetUIpointer();

  UImanager->ApplyCommand("/gps/particle opticalphoton");
  UImanager->ApplyCommand("/gps/energy 2.5 eV");
  UImanager->ApplyCommand("/gps/pos/type Point");
  UImanager->ApplyCommand("/gps/ang/type iso");
  UImanager->ApplyCommand("/gps/ang/rot1 0 -1 0");
  UImanager->ApplyCommand("/gps/ang/rot2 0 0 1");
  UImanager->ApplyCommand("/gps/ang/mintheta 0 deg");
  UImanager->ApplyCommand("/gps/ang/maxtheta 90 deg");

  steppingAction->SetVirtualPlane(x_virtual);
  {
    auto lens1      = det->GetLens75Params();
    double x1_front = lens1.x - 0.5 * lens1.tc - 1e-3;
    double r1_lens  = 0.5 * det->GetLens75Diameter();
    double x2_front = det->GetLens60X() - 0.5 * det->GetLens60Thickness() - 1e-3;
    double r2_lens  = 0.5 * det->GetLens60Diameter();
    steppingAction->SetLensAperturePlanes(x1_front, r1_lens, x2_front, r2_lens);
  }

  double h1 = det->GetLens75Thickness();
  double h2 = det->GetLens60Thickness();

  std::vector<std::pair<double, double>> pairs;
  if (config.contains("pairs")) {
    for (const auto& p : config["pairs"]) {
      pairs.push_back({p[0].get<double>(), p[1].get<double>()});
    }
  } else {
    const double margin = 3.0;
    double x1_start     = config.value("x1_start", x_min);
    double x1_end       = config.value("x1_end", x_max);

    for (double x1 = x1_start; x1 <= x1_end + 1e-9; x1 += dx) {
      double x2_min_collision = x1 + (h1 + h2) / 2.0 + margin;
      double x2_start         = std::max(x1 + dx, x2_min_collision);
      for (double x2 = x2_start; x2 <= x_max + 1e-9; x2 += dx) {
        pairs.push_back({x1, x2});
      }
    }
  }

  spdlog::info("PSF+DoF scan lens pair: {} & {}", lens75_id, lens60_id);

  for (const auto& pair : pairs) {
    double x1 = pair.first;
    double x2 = pair.second;

    det->SetLensPositions(x1, x2);
    run_manager->GeometryHasBeenModified();

    analysisManager->FillNtupleIColumn(0, 0, config_counter);
    analysisManager->FillNtupleDColumn(0, 1, x1);
    analysisManager->FillNtupleDColumn(0, 2, x2);
    analysisManager->FillNtupleDColumn(0, 3, x_virtual);
    analysisManager->FillNtupleSColumn(0, 4, lens75_id);
    analysisManager->FillNtupleSColumn(0, 5, lens60_id);
    analysisManager->AddNtupleRow(0);

    for (double x_source = source_x_min; x_source <= source_x_max + 1e-9; x_source += source_dx) {
      for (double y_source = source_y_min; y_source <= source_y_max + 1e-9; y_source += source_dy) {
        UImanager->ApplyCommand("/gps/pos/centre " + std::to_string(x_source) + " "
                                + std::to_string(y_source) + " 0 mm");

        if (config.value("use_importance_sampling", false)) {
          G4ThreeVector pos(x_source, y_source, 0);
          auto params = det->GetLens75Params();
          G4ThreeVector axis;
          double maxTheta = 0.0;
          axis            = (G4ThreeVector(params.x, 0, 0) - pos).unit();
          ImportanceSamplingHelper::CalculateCone(pos, params, axis, maxTheta);
          primaryGen->SetStaticCone(axis, maxTheta);
        }

        int run_id = run_counter++;
        eventAction->ClearRays();

        long seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        CLHEP::HepRandom::setTheSeed(seed);

        if (macro_file.empty()) {
          UImanager->ApplyCommand("/run/beamOn " + std::to_string(n_photons));
        } else {
          UImanager->ApplyCommand("/control/execute " + macro_file.string());
        }

        auto stats = compute_stats(eventAction->Moments());

        analysisManager->FillNtupleIColumn(1, 0, config_counter);
        analysisManager->FillNtupleIColumn(1, 1, run_id);
        analysisManager->FillNtupleDColumn(1, 2, x_source);
        analysisManager->FillNtupleDColumn(1, 3, y_source);
        analysisManager->FillNtupleDColumn(1, 4, stats.n_hits);
        analysisManager->FillNtupleDColumn(1, 5, stats.mu_y);
        analysisManager->FillNtupleDColumn(1, 6, stats.mu_z);
        analysisManager->FillNtupleDColumn(1, 7, stats.sigma_y);
        analysisManager->FillNtupleDColumn(1, 8, stats.sigma_z);
        analysisManager->FillNtupleDColumn(1, 9, stats.sigma_yz);
        analysisManager->FillNtupleDColumn(1, 10, stats.mu_dy);
        analysisManager->FillNtupleDColumn(1, 11, stats.mu_dz);
        analysisManager->FillNtupleDColumn(1, 12, stats.sigma_dy);
        analysisManager->FillNtupleDColumn(1, 13, stats.sigma_dz);
        analysisManager->FillNtupleDColumn(1, 14, stats.cov_y_dy);
        analysisManager->FillNtupleDColumn(1, 15, stats.cov_z_dz);
        analysisManager->AddNtupleRow(1);

        spdlog::info("Run done: config_id={}, run_id={}, x_src={}, y_src={}, hits_w={:.1f}",
                     config_counter, run_id, x_source, y_source, stats.n_hits);
      }
    }

    config_counter++;
  }

  analysisManager->Write();
  analysisManager->CloseFile();
  spdlog::info("PSF+DoF scan completed. Total configs: {}, Total runs: {}", config_counter,
               run_counter);
}

} // namespace riptide
