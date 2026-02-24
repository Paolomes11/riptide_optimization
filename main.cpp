#include "physics_list.hpp"
#include "detector_construction.hpp"
#include "action_initialization.hpp"

#include <G4RunManager.hh>
#include <G4UImanager.hh>
#include <G4VisExecutive.hh>
#include <G4UIExecutive.hh>

int main(int argc, char **argv)
{
    try
    {
        // Posizioni delle lenti di default
        double lens75_x = 83.9;
        double lens60_x = 153.4;
        G4String macro_file = "run1.mac";

        // Legge i parametri da riga di comando
        if (argc >= 3)
        {
            lens75_x = std::stod(argv[1]);
            lens60_x = std::stod(argv[2]);
        }
        if (argc >= 4)
        {
            macro_file = argv[3];
        }

        std::cout << "Using lens positions: "
                  << "lens75_x = " << lens75_x << " mm, "
                  << "lens60_x = " << lens60_x << " mm\n";
        std::cout << "Macro file: " << macro_file << "\n";

        // Create the run manager
        G4RunManager run_manager{};

        // Set mandatory initialization classes
        run_manager.SetUserInitialization(new riptide::DetectorConstruction("geometry/main.gdml", lens75_x, lens60_x));
        run_manager.SetUserInitialization(new riptide::PhysicsList());
        run_manager.SetUserInitialization(new riptide::ActionInitialization());

        // Initialize G4 kernel
        run_manager.Initialize();

        G4UIExecutive *ui = nullptr;
        if (argc == 1)
        {
            ui = new G4UIExecutive(argc, argv);
        }

        // Initialize visualization with the default graphics system
        auto visManager = new G4VisExecutive(argc, argv);
        visManager->Initialize();

        // Get the pointer to the User Interface manager
        auto UImanager = G4UImanager::GetUIpointer();

        // Process macro or start UI session
        if (ui)
        {
            UImanager->ApplyCommand("/control/execute macro/run1.mac");
            ui->SessionStart();
        }
        else
        {
            UImanager->ApplyCommand("/control/execute macro/" + macro_file);
        }

        // Clean up
        delete visManager;
        delete ui;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    catch (...)
    {
        std::cerr << "Unknown error occurred \n";
        return 1;
    }
}