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
#include "G4Event.hh"
#include <iostream>

namespace riptide {

EventAction::EventAction(const std::string& output_file_name) {
  m_file = new TFile(output_file_name.c_str(), "RECREATE");
  m_tree = new TTree("events", "Eventi fotoni");

  // Definisce i rami dell'albero
  m_tree->Branch("x1", &x1, "x1/D");
  m_tree->Branch("x2", &x2, "x2/D");
  m_tree->Branch("config_id", &config_id, "config_id/I");
}

EventAction::~EventAction() {
  m_file->cd();
  m_tree->Write();
  m_file->Close();
  delete m_file;
}

// Funzione che viene chiamata dal SensitiveDetector per registrare un hit
void EventAction::AddPhotonHit(double lens_x1, double lens_x2) {
  PhotonHit photon;
  photon.x1       = lens_x1;
  photon.x2       = lens_x2;
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

  for (auto& photon : eventHits) {
    x1       = photon.x1;
    x2       = photon.x2;
    m_tree->Fill();
  }

  eventHits.clear(); // prepara per il prossimo evento
}

} // namespace riptide
