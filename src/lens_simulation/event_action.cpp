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
void EventAction::AddPhotonHit(double y, double z, double weight) {
  eventHits.push_back({static_cast<float>(y), static_cast<float>(z)});
  m_lastRunHitCount += weight;
}

void EventAction::BeginOfEventAction(const G4Event* /*event*/) {
  s_currentEventAction = this; // imposta il puntatore statico
}

// Alla fine dell'evento, aggiorna il conteggio totale dei fotoni
void EventAction::EndOfEventAction(const G4Event* /*event*/) {
  // m_lastRunHitCount è aggiornato in AddPhotonHit e tiene conto dei pesi.
  // Non puliamo eventHits qui, lo facciamo manualmente in lens_scan/optimizer
  // all'inizio di ogni run per permettere l'accumulo vettoriale.
}

} // namespace riptide
