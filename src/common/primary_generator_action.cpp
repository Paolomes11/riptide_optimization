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

#include "primary_generator_action.hpp"

#include <G4Event.hh>
#include <G4GeneralParticleSource.hh>
#include <G4PrimaryParticle.hh>
#include <G4PrimaryVertex.hh>
#include <G4SPSPosDistribution.hh>
#include <G4SingleParticleSource.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>

namespace riptide {

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : m_gps{new G4GeneralParticleSource()}
    , m_useImportanceSampling{false}
    , m_paramsProvider{nullptr}
    , m_hasStaticCone{false}
    , m_staticAxis{0, 0, 0}
    , m_staticMaxTheta{0.0} {
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  delete m_gps;
}

void PrimaryGeneratorAction::ConfigureSpotCycling(const std::vector<SpotConfig>& spots) {
  m_spotConfigs    = spots;
  m_currentSpotIdx = 0;
}

void PrimaryGeneratorAction::ResetCycling() {
  m_currentSpotIdx = 0;
}

void PrimaryGeneratorAction::DisableCycling() {
  m_spotConfigs.clear();
  m_currentSpotIdx = 0;
}

static void apply_cone_sampling(G4Event* event, const G4ThreeVector& axis, double maxTheta) {
  G4PrimaryVertex* vertex = event->GetPrimaryVertex();
  if (!vertex)
    return;
  G4PrimaryParticle* particle = vertex->GetPrimary();
  if (!particle)
    return;

  double cosMaxTheta = std::cos(maxTheta);
  double rand1       = G4UniformRand();
  double rand2       = G4UniformRand();

  double cosTheta = 1.0 - rand1 * (1.0 - cosMaxTheta);
  double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
  double phi      = 2.0 * M_PI * rand2;

  G4ThreeVector localDir(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
  G4ThreeVector newDir = localDir;
  newDir.rotateUz(axis);
  particle->SetMomentumDirection(newDir);

  double newOmega = 2.0 * M_PI * (1.0 - cosMaxTheta);
  particle->SetWeight(newOmega / (2.0 * M_PI));
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  if (!m_spotConfigs.empty()) {
    const auto& spot = m_spotConfigs[m_currentSpotIdx];

    // Imposta posizione GPS per questo spot direttamente via C++ API
    m_gps->GetCurrentSource()->GetPosDist()->SetCentreCoords(spot.pos);

    m_gps->GeneratePrimaryVertex(event);

    if (spot.use_is) {
      apply_cone_sampling(event, spot.is_axis, spot.is_theta);
    }

    m_currentSpotIdx = (m_currentSpotIdx + 1) % static_cast<int>(m_spotConfigs.size());
    return;
  }

  // Modalità legacy
  m_gps->GeneratePrimaryVertex(event);

  if (m_useImportanceSampling) {
    G4PrimaryVertex* vertex = event->GetPrimaryVertex();
    if (!vertex)
      return;

    G4PrimaryParticle* particle = vertex->GetPrimary();
    if (!particle)
      return;

    G4ThreeVector axis;
    double maxTheta;

    if (m_hasStaticCone) {
      axis     = m_staticAxis;
      maxTheta = m_staticMaxTheta;
    } else if (m_paramsProvider) {
      auto params       = m_paramsProvider();
      G4ThreeVector pos = vertex->GetPosition();
      axis              = (params.x > pos.x() ? G4ThreeVector(1, 0, 0) : G4ThreeVector(-1, 0, 0));
      ImportanceSamplingHelper::CalculateCone(pos, params, axis, maxTheta);
    } else {
      return;
    }

    apply_cone_sampling(event, axis, maxTheta);
  }
}

} // namespace riptide