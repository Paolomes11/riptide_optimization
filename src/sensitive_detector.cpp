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
        auto *post_step = step->GetPostStepPoint();

        if (track->GetDefinition() != G4OpticalPhoton::Definition())
        {
            return false;
        }

        auto event_id = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
        auto position = post_step->GetPosition();

        std::cout << event_id << " " << position.x() << " " << position.y() << " " << position.z() << "\n";

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
