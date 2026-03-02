#include "sensitive_detector.hpp"
#include "detector_construction.hpp"
#include "event_action.hpp"

#include <G4OpticalPhoton.hh>
#include <G4RunManager.hh>

#include <iostream>

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

  // Recupera EventAction tramite il puntatore statico
  auto* eventAction = EventAction::GetEventAction();
  if (!eventAction) {
    std::cerr << "SensitivePhotocathode: EventAction not found!" << std::endl;
    return false;
  }

  // Recupera DetectorConstruction per leggere le posizioni correnti delle lenti
  auto* det = dynamic_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  if (!det) {
    std::cerr << "SensitivePhotocathode: DetectorConstruction not found!" << std::endl;
    return false;
  }

  double current_x1 = det->GetLens75X();
  double current_x2 = det->GetLens60X();

  // Registra il fotone come "hit"
  eventAction->AddPhotonHit(current_x1, current_x2);

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
