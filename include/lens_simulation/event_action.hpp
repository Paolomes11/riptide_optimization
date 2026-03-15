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

#ifndef RIPTIDE_EVENT_ACTION_HPP
#define RIPTIDE_EVENT_ACTION_HPP

#include <G4UserEventAction.hh>

#include <vector>

namespace riptide {

// Struttura per un hit/fotone
struct PhotonHit {
  // Posizione del fotone
  float y;
  float z;
};

class EventAction : public G4UserEventAction {
  // Vector per memorizzare tutti i fotoni di un evento
  std::vector<PhotonHit> eventHits;

  // Per evitare run_id
  int m_lastRunHitCount = 0;

  // puntatore statico
  static inline EventAction* s_currentEventAction = nullptr;

 public:
  EventAction()          = default;
  virtual ~EventAction() = default;

  // Funzione che viene chiamata dal SensitiveDetector per registrare un hit
  void AddPhotonHit(double y, double z);

  // Getter per ottenere tutti gli hit di questo evento
  const std::vector<PhotonHit>& GetEventHits() const {
    return eventHits;
  }

  // Pulisce gli hit dopo averli salvati
  void ClearEventHits() {
    eventHits.clear();
  }

  // n_hits getter e clearer
  int GetLastRunHitCount() const {
    return m_lastRunHitCount;
  }

  void ResetLastRunHitCount() {
    m_lastRunHitCount = 0;
  }

  virtual void EndOfEventAction(const G4Event* event) override;
  virtual void BeginOfEventAction(const G4Event* event) override;

  // Metodo statico per ottenere l'EventAction corrente
  static EventAction* GetEventAction() {
    return s_currentEventAction;
  }

  // Identificatore del run corrente, da impostare in lens_scan.cpp
  int runID = -1;
};

} // namespace riptide

#endif