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


namespace riptide {

G4bool SensitivePhotocathode::ProcessHits(G4Step* step, G4TouchableHistory*) {
  auto* pre = step->GetPreStepPoint();

  if (pre->GetStepStatus() != fGeomBoundary)
    return false;

  auto* eventAction = EventAction::GetEventAction();
  if (!eventAction)
    return false;

  double weight = step->GetTrack()->GetWeight();
  eventAction->AddPhotonHit(weight);
  step->GetTrack()->SetTrackStatus(fStopAndKill);
  return true;
}

} // namespace riptide
