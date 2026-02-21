#include "sensitive_detector.hpp"

// #include <G4AnalysisManager.hh>
// #include <G4EventManager.hh>
#include <G4OpticalPhoton.hh>
#include <G4SystemOfUnits.hh>
#include <G4SDManager.hh>
// #include <G4VProcess.hh>

// #include <fmt/format.h>
// #include <tuple>
#include <iostream>

namespace riptide
{

    SensitivePhotocathode::SensitivePhotocathode(G4String const &name, int detector_id)
        : G4VSensitiveDetector{name}, m_detector_id{detector_id}
    {
        collectionName.insert(name + "_hits");
    }

    void SensitivePhotocathode::Initialize(G4HCofThisEvent *hce)
    {
        m_hits_collection = new PhotocathodeHitsCollection(SensitiveDetectorName, collectionName[0]);
        if (m_hc_id < 0)
        {
            m_hc_id = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
        }
        hce->AddHitsCollection(m_hc_id, m_hits_collection);
    }

    G4bool SensitivePhotocathode::ProcessHits(G4Step *step, G4TouchableHistory *history)
    {
        std::ignore = history;
        auto *track = step->GetTrack();
        auto *post_step = step->GetPostStepPoint();

        if (track->GetDefinition() != G4OpticalPhoton::Definition())
        {
            return false;
        }

        // Controlla che m_hits_collection non sia nullptr
        if (!m_hits_collection)
        {
            G4ExceptionDescription msg;
            msg << "Hits collection is null in ProcessHits";
            G4Exception("SensitivePhotocathode::ProcessHits", "PHOTOCATHODE001", JustWarning, msg);
            return false;
        }

        auto *hit = new PhotocathodeHit();
        hit->set_global_position(post_step->GetPosition());
        hit->set_time(post_step->GetGlobalTime());
        hit->set_energy(track->GetKineticEnergy());
        hit->set_detector_id(m_detector_id);

        m_hits_collection->insert(hit);

        // assorbe il fotone
        track->SetTrackStatus(fStopAndKill);

        return true;
    }

} // namespace riptide
