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

void EventAction::AddPhotonHit(double y, double z, double weight) {
  eventHits.push_back({static_cast<float>(y), static_cast<float>(z)});
  m_lastRunHitCount += weight;
}

void EventAction::AddPhotonHit(int spot_id, double y, double z, double weight) {
  if (spot_id >= 0 && spot_id < static_cast<int>(m_spotHits.size())) {
    m_spotHits[spot_id].push_back({static_cast<float>(y), static_cast<float>(z)});
    m_spotWeights[spot_id] += weight;
  }
}

void EventAction::InitSpotMode(int n_spots, double x_src_min, double dx_src,
                               double y_src_min, double dy_src, int ny_src) {
  m_spotMode = true;
  m_xSrcMin  = x_src_min;
  m_dxSrc    = dx_src;
  m_ySrcMin  = y_src_min;
  m_dySrc    = dy_src;
  m_nySrc    = ny_src;
  m_spotHits.assign(n_spots, {});
  m_spotWeights.assign(n_spots, 0.0);
  m_spotVirtualMoments.assign(n_spots, VirtualPlaneMoments{});
  m_spotKilledVLens1.assign(n_spots, 0);
  m_spotKilledVLens2.assign(n_spots, 0);
  m_spotKilledVBack.assign(n_spots, 0);
  m_spotVirtualHitsY.assign(n_spots, {});
  m_spotVirtualHitsZ.assign(n_spots, {});
}

void EventAction::ResetSpotAccumulators() {
  m_spotMode = false;
  m_spotHits.clear();
  m_spotWeights.clear();
  m_spotVirtualMoments.clear();
  m_spotKilledVLens1.clear();
  m_spotKilledVLens2.clear();
  m_spotKilledVBack.clear();
  m_spotVirtualHitsY.clear();
  m_spotVirtualHitsZ.clear();
}

void EventAction::AddVirtualHit(int spot_id, double y, double z,
                                 double dy, double dz, double w) {
  if (spot_id < 0 || spot_id >= static_cast<int>(m_spotVirtualMoments.size())) return;
  auto& m    = m_spotVirtualMoments[spot_id];
  m.sum_w   += w;
  m.sum_y   += w * y;
  m.sum_z   += w * z;
  m.sum_dy  += w * dy;
  m.sum_dz  += w * dz;
  m.sum_yy  += w * y * y;
  m.sum_zz  += w * z * z;
  m.sum_yz  += w * y * z;
  m.sum_dy_dy += w * dy * dy;
  m.sum_dz_dz += w * dz * dz;
  m.sum_y_dy  += w * y * dy;
  m.sum_z_dz  += w * z * dz;

  if (m_saveVirtualHits) {
    m_spotVirtualHitsY[spot_id].push_back(static_cast<float>(y));
    m_spotVirtualHitsZ[spot_id].push_back(static_cast<float>(z));
  }
}

void EventAction::AddKilledVLens1(int spot_id) {
  if (spot_id >= 0 && spot_id < static_cast<int>(m_spotKilledVLens1.size()))
    ++m_spotKilledVLens1[spot_id];
}

void EventAction::AddKilledVLens2(int spot_id) {
  if (spot_id >= 0 && spot_id < static_cast<int>(m_spotKilledVLens2.size()))
    ++m_spotKilledVLens2[spot_id];
}

void EventAction::AddKilledVBack(int spot_id) {
  if (spot_id >= 0 && spot_id < static_cast<int>(m_spotKilledVBack.size()))
    ++m_spotKilledVBack[spot_id];
}

void EventAction::BeginOfEventAction(const G4Event* /*event*/) {
  s_currentEventAction = this;
}

void EventAction::EndOfEventAction(const G4Event* /*event*/) {
  // Gli hit si accumulano nel run (sia legacy che spot-mode).
  // La pulizia è esplicita: ClearEventHits() o ResetSpotAccumulators().
}

} // namespace riptide
