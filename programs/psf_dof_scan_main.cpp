#include "physics_list.hpp"
#include "psf_dof_action_initialization.hpp"
#include "psf_dof_detector_construction.hpp"
#include "psf_dof_scan.hpp"

#include <G4RunManager.hh>
#include <G4UIExecutive.hh>
#include <G4UImanager.hh>
#include <G4VisExecutive.hh>

#include <lyra/lyra.hpp>
#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>

#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>

int main(int argc, char** argv) {
  std::filesystem::path geometry_path = "geometry/dof_geometry.gdml";
  std::filesystem::path macro_file;
  std::filesystem::path macro_vis   = "macros/vis.mac";
  std::filesystem::path macro_run   = "macros/lens_simulation.mac";
  std::filesystem::path config_file = "config/config.json";
  std::string lens75_id;
  std::string lens60_id;
  bool visualize = false;
  bool batch     = false;
  bool psf_dof   = false;
  bool show_help = false;
  bool use_ssd   = false;

  std::string root_output_file = "output/psf_dof_simulation/psf_dof.root";
  std::string ssd_mount        = "/mnt/external_ssd";

  auto cli = lyra::cli() | lyra::help(show_help)
           | lyra::opt(geometry_path, "geometry")["-g"]["--geometry"]("Path to GDML geometry file")
                 .required()
           | lyra::opt(macro_file, "macro")["-m"]["--macro"]("Path to macro file (default: none)")
           | lyra::opt(config_file, "config")["--config"](
                 "Path to config JSON file (default: config/config.json)")
           | lyra::opt(lens75_id, "id")["--lens75-id"]("Thorlabs ID for lens 1 (75mm)")
           | lyra::opt(lens60_id, "id")["--lens60-id"]("Thorlabs ID for lens 2 (60mm)")
           | lyra::opt(root_output_file, "output")["--output"]("Path to ROOT output file")
           | lyra::opt(ssd_mount, "ssd-mount")["--ssd-mount"](
                 "Mount point of external SSD (default: /mnt/external_ssd)")
           | lyra::opt(use_ssd)["--ssd"]("Write output to external SSD (uses --ssd-mount)")
           | lyra::opt(visualize)["-v"]["--visualize"]("Enable visualization mode")
           | lyra::opt(batch)["-b"]["--batch"]("Enable batch mode (no visualization)")
           | lyra::opt(psf_dof)["-p"]["--psf-dof"]("Enable PSF+DoF scan mode");

  auto result = cli.parse({argc, argv});
  if (!result) {
    spdlog::error("CLI parsing error: {}", result.message());
    std::cerr << cli << '\n';
    return EXIT_FAILURE;
  }

  if (show_help) {
    std::cout << cli << '\n';
    return EXIT_SUCCESS;
  }

  if (visualize && batch) {
    spdlog::error("Cannot use --visualize and --batch together.");
    return EXIT_FAILURE;
  }

  if (use_ssd) {
    std::string timestamp = []() {
      auto t  = std::time(nullptr);
      auto tm = *std::localtime(&t);
      std::ostringstream oss;
      oss << std::put_time(&tm, "%Y%m%d_%H%M%S");
      return oss.str();
    }();
    root_output_file = ssd_mount + "/riptide/runs/run_" + timestamp + "/psf_dof.root";
    spdlog::info("SSD mode: output -> {}", root_output_file);
  }

  spdlog::info("Output file : {}", root_output_file);
  spdlog::info("Config file : {}", config_file.string());

  try {
    G4RunManager run_manager{};

    run_manager.GeometryHasBeenModified(true);
    if (!lens75_id.empty() && !lens60_id.empty()) {
      run_manager.SetUserInitialization(
          new riptide::PsfDofDetectorConstruction(geometry_path.string(), lens75_id, lens60_id));
    } else {
      run_manager.SetUserInitialization(
          new riptide::PsfDofDetectorConstruction(geometry_path.string()));
    }

    bool use_importance_sampling = false;
    std::ifstream f(config_file);
    if (f.is_open()) {
      auto config             = nlohmann::json::parse(f);
      use_importance_sampling = config.value("use_importance_sampling", false);
      if (use_importance_sampling) {
        spdlog::info("Geometric importance sampling enabled");
      }
    }

    run_manager.SetUserInitialization(new riptide::PhysicsList());
    run_manager.SetUserInitialization(
        new riptide::PsfDofActionInitialization(use_importance_sampling));
    run_manager.Initialize();

    auto UImanager              = G4UImanager::GetUIpointer();
    G4UIExecutive* ui           = nullptr;
    G4VisExecutive* vis_manager = nullptr;

    if (psf_dof) {
      spdlog::info("Running PSF+DoF scan");
      riptide::run_psf_dof_scan(&run_manager, macro_file, root_output_file, config_file, lens75_id,
                                lens60_id);
      return EXIT_SUCCESS;
    }

    if (visualize) {
      ui          = new G4UIExecutive(argc, argv);
      vis_manager = new G4VisExecutive();
      vis_manager->Initialize();

      UImanager->ApplyCommand("/control/execute " + macro_vis.string());
      if (!macro_file.empty()) {
        UImanager->ApplyCommand("/control/execute " + macro_file.string());
      } else {
        UImanager->ApplyCommand("/control/execute " + macro_run.string());
      }
      ui->SessionStart();
      delete ui;
      delete vis_manager;
      return EXIT_SUCCESS;
    }

    if (batch || !visualize) {
      if (!macro_file.empty()) {
        UImanager->ApplyCommand("/control/execute " + macro_file.string());
      } else {
        UImanager->ApplyCommand("/control/execute " + macro_run.string());
      }
    }
  } catch (const std::exception& e) {
    spdlog::error("Error: {}", e.what());
    return EXIT_FAILURE;
  } catch (...) {
    spdlog::error("Unknown error occurred");
    return EXIT_FAILURE;
  }
}
