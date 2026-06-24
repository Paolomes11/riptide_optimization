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
}

void EventAction::ResetSpotAccumulators() {
  m_spotMode = false;
  m_spotHits.clear();
  m_spotWeights.clear();
}

void EventAction::BeginOfEventAction(const G4Event* /*event*/) {
  s_currentEventAction = this;
}

void EventAction::EndOfEventAction(const G4Event* /*event*/) {
  // Gli hit si accumulano nel run (sia legacy che spot-mode).
  // La pulizia è esplicita: ClearEventHits() o ResetSpotAccumulators().
}

} // namespace riptide
