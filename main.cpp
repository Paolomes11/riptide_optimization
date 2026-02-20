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
        // Create the run manager
        G4RunManager run_manager{};

        // Set mandatory initialization classes
        run_manager.SetUserInitialization(new riptide::DetectorConstruction("geometry/main.gdml"));
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
        // auto visManager = new G4VisExecutive(argc, argv);
        // visManager->Initialize();

        // Get the pointer to the User Interface manager
        auto UImanager = G4UImanager::GetUIpointer();

        // Process macro or start UI session
        if (!ui)
        {
            // batch mode
            if (argc > 1)
            {
                G4String macroFile = argv[1];
                UImanager->ApplyCommand("/control/execute " + macroFile);
            }
        }
        else
        {
            // interactive mode
            UImanager->ApplyCommand("/control/execute run1.mac");
            ui->SessionStart();
            delete ui;
        }

        // Clean up
        // delete visManager;
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