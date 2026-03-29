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

#include "stepping_action.hpp"

#include <G4OpticalPhoton.hh>
#include <G4Step.hh>
#include <G4Track.hh>

namespace riptide {

SteppingAction::SteppingAction() : G4UserSteppingAction() {}

void SteppingAction::UserSteppingAction(const G4Step* step) {
  auto* track = step->GetTrack();

  // Considera solo i fotoni ottici
  if (track->GetDefinition() != G4OpticalPhoton::Definition()) {
    return;
  }

  // Uccide i fotoni che vanno nella direzione sbagliata (X negativa)
  // Assumendo che la sorgente sia a X=0 e l'apparato in +X
  auto pos      = track->GetPosition();
  auto momentum = track->GetMomentum();

  // Se il fotone ha superato la sorgente (X < -1 mm) e va verso -X, lo uccidiamo
  if (pos.x() < -1.0 && momentum.x() < 0) {
    track->SetTrackStatus(fStopAndKill);
    return;
  }

  // Se il fotone è troppo lontano dall'asse ottico (Y o Z > 150 mm), lo uccidiamo
  // I diametri delle lenti sono max 77mm, quindi 150mm è un raggio molto conservativo.
  if (std::abs(pos.y()) > 150.0 || std::abs(pos.z()) > 150.0) {
    track->SetTrackStatus(fStopAndKill);
    return;
  }
}

} // namespace riptide
