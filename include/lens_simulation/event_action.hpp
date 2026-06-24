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

#include "spot_grid.hpp"

#include <G4UserEventAction.hh>

#include <cmath>
#include <vector>

namespace riptide {

// Struttura per un hit/fotone
struct PhotonHit {
  float y;
  float z;
};

class EventAction : public G4UserEventAction {
  // Buffer legacy (un run = uno spot)
  std::vector<PhotonHit> eventHits;
  double m_lastRunHitCount = 0;

  // Buffer spot-mode (un run = tutti gli spot)
  bool m_spotMode = false;
  std::vector<std::vector<PhotonHit>> m_spotHits;
  std::vector<double> m_spotWeights;

  // Parametri griglia sorgente (per calcolo spot_id in ProcessHits)
  double m_xSrcMin = 0.0;
  double m_dxSrc   = 1.0;
  double m_ySrcMin = 0.0;
  double m_dySrc   = 1.0;
  int m_nySrc      = 1;

  static inline EventAction* s_currentEventAction = nullptr;

 public:
  EventAction()          = default;
  virtual ~EventAction() = default;

  // Legacy: hit senza spot_id (usato da optimizer e da SD in legacy mode)
  void AddPhotonHit(double y, double z, double weight = 1.0);

  // Spot-mode: hit con spot_id pre-calcolato
  void AddPhotonHit(int spot_id, double y, double z, double weight);

  bool IsSpotMode() const {
    return m_spotMode;
  }

  // Calcola spot_id da coordinate sorgente (in mm)
  int CalcSpotId(double x_src_mm, double y_src_mm) const {
    return calc_spot_id(x_src_mm, y_src_mm, m_xSrcMin, m_dxSrc, m_ySrcMin, m_dySrc, m_nySrc);
  }

  // Abilita spot-mode e inizializza buffer
  void InitSpotMode(int n_spots, double x_src_min, double dx_src, double y_src_min, double dy_src,
                    int ny_src);

  // Resetta buffer spot (chiamare tra una pair e l'altra)
  void ResetSpotAccumulators();

  // Getters spot-mode
  const std::vector<std::vector<PhotonHit>>& GetSpotHits() const {
    return m_spotHits;
  }
  const std::vector<double>& GetSpotWeights() const {
    return m_spotWeights;
  }

  // Getter legacy
  const std::vector<PhotonHit>& GetEventHits() const {
    return eventHits;
  }
  void ClearEventHits() {
    eventHits.clear();
  }
  double GetLastRunHitCount() const {
    return m_lastRunHitCount;
  }
  void ResetLastRunHitCount() {
    m_lastRunHitCount = 0;
  }

  virtual void EndOfEventAction(const G4Event* event) override;
  virtual void BeginOfEventAction(const G4Event* event) override;

  static EventAction* GetEventAction() {
    return s_currentEventAction;
  }

  int runID = -1;
};

} // namespace riptide

#endif