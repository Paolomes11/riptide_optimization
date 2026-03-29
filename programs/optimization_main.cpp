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

/*
    Un'altra versione del main con una migliore gestione dei flag da riga di comando usando lyra e
   spdlog.
*/

#include "action_initialization.hpp"
#include "detector_construction.hpp"
#include "optimizer.hpp"
#include "physics_list.hpp"

#include <G4RunManager.hh>
#include <G4UIExecutive.hh>
#include <G4UImanager.hh>
#include <G4VisExecutive.hh>

#include <lyra/lyra.hpp>
#include <spdlog/spdlog.h>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>

int main(int argc, char** argv) {
  std::filesystem::path geometry_path = "geometry/main.gdml";
  std::filesystem::path macro_file;
  std::filesystem::path macro_vis   = "macros/vis.mac";
  std::filesystem::path macro_run   = "macros/optimization.mac";
  std::filesystem::path config_file = "config/config.json";
  std::string lens75_id;
  std::string lens60_id;
  bool visualize                    = false;
  bool batch                        = false;
  bool optimize                     = false;
  bool all_lenses                   = false;
  bool show_help                    = false;
  bool use_ssd                      = false;

  // Path di output: default locale, sovrascrivibile da CLI o --ssd
  std::string root_output_file = "output/events.root";
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
           | lyra::opt(optimize)["-o"]["--optimize"]("Enable geometrical efficiency optimization")
           | lyra::opt(all_lenses)["--all-lenses"](
                 "Simulate all combinations of Thorlabs lenses (may use a lot of disk space!)");

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

  // Risolvi il path di output: --ssd sovrascrive --output con un path datato sull'SSD
  if (use_ssd) {
    std::string timestamp = []() {
      auto t  = std::time(nullptr);
      auto tm = *std::localtime(&t);
      std::ostringstream oss;
      oss << std::put_time(&tm, "%Y%m%d_%H%M%S");
      return oss.str();
    }();
    root_output_file = ssd_mount + "/riptide/runs/run_" + timestamp + "/events.root";
    spdlog::info("SSD mode: output -> {}", root_output_file);
  }

  spdlog::info("Output file : {}", root_output_file);
  spdlog::info("Config file : {}", config_file.string());

  // Try catch errors
  try {
    // Crea il run manager
    G4RunManager run_manager{};

    run_manager.GeometryHasBeenModified(true);
    if (!lens75_id.empty() && !lens60_id.empty()) {
      run_manager.SetUserInitialization(
          new riptide::DetectorConstruction(geometry_path.string(), lens75_id, lens60_id, 83.9,
                                            153.4));
    } else {
      run_manager.SetUserInitialization(
          new riptide::DetectorConstruction(geometry_path.string(), 83.9, 153.4));
    }
    run_manager.SetUserInitialization(new riptide::PhysicsList());
    run_manager.SetUserInitialization(new riptide::ActionInitialization());
    run_manager.Initialize();

    auto UImanager              = G4UImanager::GetUIpointer();
    G4UIExecutive* ui           = nullptr;
    G4VisExecutive* vis_manager = nullptr;

    if (optimize) {
      if (all_lenses && !use_ssd) {
        spdlog::warn("Option --all-lenses used without --ssd. Output file might be very large!");
      }
      spdlog::info("Running optimization");
      riptide::run_optimization(&run_manager, macro_file, root_output_file, config_file, all_lenses);
      return EXIT_SUCCESS;
    }

    // Visualizzazione
    if (visualize) {
      spdlog::info("Visualization mode enabled");
      ui          = new G4UIExecutive(argc, argv);
      vis_manager = new G4VisExecutive();
      vis_manager->Initialize();

      // Prima esegui la macro grafica standard
      UImanager->ApplyCommand("/control/execute " + macro_vis.string());

      // Poi esegui la macro di simulazione (o quella specificata)
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

    // Batch mode
    if (batch || !visualize) {
      spdlog::info("Batch mode enabled");
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