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
    return false;
  }

  // Considera solo i fotoni che attraversano la superficie del fotocatodo
  auto* pre = step->GetPreStepPoint();
  if (pre->GetStepStatus() != fGeomBoundary) {
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

  // Registra il fotone come "hit"
  eventAction->AddPhotonHit();

  track->SetTrackStatus(fStopAndKill);

  return true;
}

} // namespace riptide
