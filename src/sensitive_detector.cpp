#include "sensitive_detector.hpp"
#include "efficiency_collector.hpp"

#include <G4AnalysisManager.hh>
#include <G4EventManager.hh>
#include <G4OpticalPhoton.hh>
#include <G4SystemOfUnits.hh>
#include <G4VProcess.hh>

#include <fmt/format.h>
#include <iostream>
#include <tuple>

namespace riptide {

G4bool SensitivePhotocathode::ProcessHits(G4Step* step, G4TouchableHistory* history) {
  std::ignore = history;
  auto* track = step->GetTrack();

  // Considera solo i fotoni ottici
  if (track->GetDefinition() != G4OpticalPhoton::Definition()) {
    std::cout << "SensitiveDetector: Not an optical photon, skipping" << std::endl;
    return false;
  }

  // Considera solo i fotoni che attraversano la superficie del fotocatodo
  auto* pre = step->GetPreStepPoint();
  if (pre->GetStepStatus() != fGeomBoundary) {
    std::cout << "SensitiveDetector: Not at geometry boundary, skipping" << std::endl;
    return false;
  }

  // Registra l'evento come evento con hit
  auto event_id = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

  // Registra l'evento con hit usando il collector globale
  auto collector = riptide::EfficiencyCollector::GetInstance();
  if (collector) {
    collector->recordEventWithHit(event_id);
  } else {
  }

  // // Salva i dati del fotone colpito nell'ntuple
  // auto position = step->GetPostStepPoint()->GetPosition();

  // G4AnalysisManager* am = G4AnalysisManager::Instance();
  // am->FillNtupleIColumn(0, event_id);
  // am->FillNtupleDColumn(1, position.x());
  // am->FillNtupleDColumn(2, position.y());
  // am->FillNtupleDColumn(3, position.z());
  // am->AddNtupleRow(0);

  return true;
}

} // namespace riptide
