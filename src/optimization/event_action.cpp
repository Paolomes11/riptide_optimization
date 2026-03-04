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

// Funzione che viene chiamata dal SensitiveDetector per registrare un hit
void EventAction::AddPhotonHit(double lens_x1, double lens_x2) {
  PhotonHit photon;
  photon.x1        = lens_x1;
  photon.x2        = lens_x2;
  photon.config_id = config_id;
  eventHits.push_back(photon);
}

// Funzione per impostare l'identificatore della configurazione
void EventAction::SetConfigId(int config_id) {
  this->config_id = config_id;
}

void EventAction::BeginOfEventAction(const G4Event* /*event*/) {
  eventHits.clear();
  s_currentEventAction = this; // imposta il puntatore statico
}

// Alla fine dell'evento, scrivi tutti i fotoni nel TTree
void EventAction::EndOfEventAction(const G4Event* event) {
  (void)event; // non serve l'evento qui

  auto analysisManager = G4AnalysisManager::Instance();

  for (auto& photon : eventHits) {
    analysisManager->FillNtupleDColumn(0, photon.x1);
    analysisManager->FillNtupleDColumn(1, photon.x2);
    analysisManager->FillNtupleIColumn(2, photon.config_id);

    analysisManager->AddNtupleRow();
  }

  eventHits.clear(); // prepara per il prossimo evento
}

} // namespace riptide
