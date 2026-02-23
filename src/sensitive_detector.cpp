#include "sensitive_detector.hpp"

#include <G4AnalysisManager.hh>
#include <G4EventManager.hh>
#include <G4OpticalPhoton.hh>
#include <G4SystemOfUnits.hh>
#include <G4VProcess.hh>

#include <fmt/format.h>
#include <tuple>
#include <iostream>

namespace riptide
{

    G4bool SensitivePhotocathode::ProcessHits(G4Step *step, G4TouchableHistory *history)
    {
        std::ignore = history;
        auto *track = step->GetTrack();

        // Considera solo i fotoni ottici
        if (track->GetDefinition() != G4OpticalPhoton::Definition())
        {
            return false;
        }

        // Considera solo i fotoni che attraversano la superficie del fotocatodo
        auto *pre = step->GetPreStepPoint();
        if (pre->GetStepStatus() != fGeomBoundary)
        {
            return false;
        }

        // Log degli ingressi dei fotoni riflessi dal fotocatodo nel fotocatodo
        // static std::map<int, int> entries;
        // auto *post = step->GetPostStepPoint();
        // if (post->GetStepStatus() == fGeomBoundary)
        // {
        //     int event_id =
        //         G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

        //     entries[event_id]++;

        //     if (entries[event_id] > 1)
        //     {
        //         std::cout << "Fotone rientrato nel detector! Event "
        //                   << event_id
        //                   << " ingressi = "
        //                   << entries[event_id]
        //                   << std::endl;
        //     }
        // }

        auto event_id = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
        auto position = step->GetPostStepPoint()->GetPosition();

        G4AnalysisManager *am = G4AnalysisManager::Instance();
        am->FillNtupleIColumn(0, event_id);
        am->FillNtupleDColumn(1, position.x());
        am->FillNtupleDColumn(2, position.y());
        am->FillNtupleDColumn(3, position.z());
        am->AddNtupleRow(0);

        // auto energy = track->GetKineticEnergy();
        // auto direction = track->GetMomentumDirection();
        // auto volume = track->GetVolume()->GetName();
        // auto process = post_step->GetProcessDefinedStep()->GetProcessName();

        return true;
    }

} // namespace riptide
