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

// Accumulatore momenti del secondo ordine per il piano virtuale (Tecnica E)
struct VirtualPlaneMoments {
  double sum_w     = 0.0;
  double sum_y     = 0.0;
  double sum_z     = 0.0;
  double sum_dy    = 0.0;
  double sum_dz    = 0.0;
  double sum_yy    = 0.0;
  double sum_zz    = 0.0;
  double sum_yz    = 0.0;
  double sum_dy_dy = 0.0;
  double sum_dz_dz = 0.0;
  double sum_y_dy  = 0.0;
  double sum_z_dz  = 0.0;
};

class EventAction : public G4UserEventAction {
  double m_lastRunHitCount = 0;

  // Buffer spot-mode (un run = tutti gli spot)
  bool m_spotMode = false;
  std::vector<std::vector<PhotonHit>> m_spotHits;
  std::vector<double> m_spotWeights;

  // Accumulatori piano virtuale (Tecnica E)
  std::vector<VirtualPlaneMoments> m_spotVirtualMoments;
  std::vector<int> m_spotKilledVLens1;
  std::vector<int> m_spotKilledVLens2;
  std::vector<int> m_spotKilledVBack;

  // Hit grezzi al piano virtuale (opt-in, default disattivo — vedi SetSaveVirtualHits)
  bool m_saveVirtualHits = false;
  std::vector<std::vector<float>> m_spotVirtualHitsY;
  std::vector<std::vector<float>> m_spotVirtualHitsZ;

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

  // Spot-mode: hit con spot_id pre-calcolato
  void AddPhotonHit(int spot_id, double y, double z, double weight);

  // Piano virtuale (Tecnica E): accumulo momenti ponderati
  void AddVirtualHit(int spot_id, double y, double z, double dy, double dz, double w);
  void AddKilledVLens1(int spot_id);
  void AddKilledVLens2(int spot_id);
  void AddKilledVBack(int spot_id);

  // Abilita/disabilita il salvataggio degli hit grezzi al piano virtuale (default: off)
  void SetSaveVirtualHits(bool enable) {
    m_saveVirtualHits = enable;
  }

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

  // Getters piano virtuale
  const std::vector<VirtualPlaneMoments>& GetSpotVirtualMoments() const {
    return m_spotVirtualMoments;
  }
  const std::vector<int>& GetSpotKilledVLens1() const {
    return m_spotKilledVLens1;
  }
  const std::vector<int>& GetSpotKilledVLens2() const {
    return m_spotKilledVLens2;
  }
  const std::vector<int>& GetSpotKilledVBack() const {
    return m_spotKilledVBack;
  }
  const std::vector<std::vector<float>>& GetSpotVirtualHitsY() const {
    return m_spotVirtualHitsY;
  }
  const std::vector<std::vector<float>>& GetSpotVirtualHitsZ() const {
    return m_spotVirtualHitsZ;
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
};

} // namespace riptide

#endif