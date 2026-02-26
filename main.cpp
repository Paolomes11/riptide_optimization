/*
    Un'altra versione del main con una migliore gestione dei flag da riga di comando usando lyra e
   spdlog.
*/

#include "action_initialization.hpp"
#include "detector_construction.hpp"
#include "efficiency_collector.hpp"
#include "optimizer.hpp"
#include "physics_list.hpp"

#include <G4RunManager.hh>
#include <G4UIExecutive.hh>
#include <G4UImanager.hh>
#include <G4VisExecutive.hh>

#include <lyra/lyra.hpp>
#include <spdlog/spdlog.h>

int main(int argc, char** argv) {
  std::filesystem::path geometry_path = "geometry/main.gdml";
  std::filesystem::path macro_file;
  std::filesystem::path macro_vis = "macros/vis.mac";  // macro grafica standard
  std::filesystem::path macro_run = "macros/run1.mac"; // macro simulazione standard
  bool visualize                  = false;
  bool batch                      = false;
  bool optimize                   = false;
  bool show_help                  = false;

  auto cli = lyra::cli() | lyra::help(show_help)
           | lyra::opt(geometry_path, "geometry")["-g"]["--geometry"]("Path to GDML geometry file")
                 .required()
           | lyra::opt(macro_file, "macro")["-m"]["--macro"]("Path to macro file (default: none)")
           | lyra::opt(visualize)["-v"]["--visualize"]("Enable visualization mode")
           | lyra::opt(batch)["-b"]["--batch"]("Enable batch mode (no visualization)")
           | lyra::opt(optimize)["-o"]["--optimize"]("Enable geometrical efficiency optimization");

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

  // Try catch errors
  try {
    // Crea il run manager
    G4RunManager run_manager{};
    auto collector = std::make_unique<riptide::EfficiencyCollector>();

    // Imposta il collector come istanza globale per il sensitive detector
    riptide::EfficiencyCollector::SetInstance(collector.get());

    run_manager.SetUserInitialization(new riptide::DetectorConstruction(geometry_path.string()));
    run_manager.SetUserInitialization(new riptide::PhysicsList());
    G4UIExecutive* ui           = nullptr;
    G4VisExecutive* vis_manager = nullptr;
    auto UImanager              = G4UImanager::GetUIpointer();
    run_manager.SetUserInitialization(new riptide::ActionInitialization(collector.get()));
    run_manager.Initialize();

    if (optimize) {
      spdlog::info("Running optimization");
      // Funzione separata, passa run manager e macro scelta
      riptide::run_optimization(&run_manager, macro_file, collector.get());
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