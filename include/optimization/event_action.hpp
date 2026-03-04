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
  double x1;
  double x2;
  int config_id; // Identificatore della configurazione delle lenti
};

class EventAction : public G4UserEventAction {
  // Variabili da salvare per ogni fotone
  int config_id; // Identificatore della configurazione delle lenti

  // Vector per memorizzare tutti i fotoni di un evento
  std::vector<PhotonHit> eventHits;

  // puntatore statico
  static inline EventAction* s_currentEventAction = nullptr;

 public:
  EventAction() = default;
  virtual ~EventAction() = default;

  // Funzione che viene chiamata dal SensitiveDetector per registrare un hit
  void AddPhotonHit(double lens_x1, double lens_x2);

  // Funzione per impostare l'identificatore della configurazione
  void SetConfigId(int config_id);

  virtual void EndOfEventAction(const G4Event* event) override;
  virtual void BeginOfEventAction(const G4Event* event) override;

  // Metodo statico per ottenere l'EventAction corrente
  static EventAction* GetEventAction() {
    return s_currentEventAction;
  }
};

} // namespace riptide

#endif