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

#include "event_action.hpp"
#include <G4AnalysisManager.hh>
#include <G4Event.hh>
#include <iostream>

namespace riptide {

void EventAction::AddPhotonHit(double y, double z) {
  PhotonHit photon;
  photon.y = y;
  photon.z = z;
  eventHits.push_back(photon);
}

void EventAction::BeginOfEventAction(const G4Event* /*event*/) {
  eventHits.clear();
  s_currentEventAction = this; // imposta il puntatore statico
}

// Alla fine dell'evento, scrivi tutti i fotoni nel TTree
void EventAction::EndOfEventAction(const G4Event* /*event*/) {
  // Ottieni il manager di analisi ROOT
  auto* analysisManager = G4AnalysisManager::Instance();

  // Scrive tutte le hit accumulate nell'evento corrente
  for (const auto& photon : eventHits) {
    analysisManager->FillNtupleFColumn(2, 0, photon.y);
    analysisManager->FillNtupleFColumn(2, 1, photon.z);
    analysisManager->AddNtupleRow(2);
  }

  m_lastRunHitCount += static_cast<int>(eventHits.size());

  // Pulisce il vettore per il prossimo evento
  eventHits.clear();
}

} // namespace riptide
