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

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  // Genera il primo vertice dell'interazione particella
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

    // Genera una nuova direzione uniforme nell'angolo solido del cono
    double cosMaxTheta = std::cos(maxTheta);
    double rand1       = G4UniformRand();
    double rand2       = G4UniformRand();

    double cosTheta = 1.0 - rand1 * (1.0 - cosMaxTheta);
    double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
    double phi      = 2.0 * M_PI * rand2;

    G4ThreeVector localDir(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);

    // Ruota la direzione locale per allinearla con l'asse del cono
    G4ThreeVector newDir = localDir;
    newDir.rotateUz(axis);

    particle->SetMomentumDirection(newDir);

    // Peso: Rapporto tra l'angolo solido del cono e l'angolo solido originale (emisfero = 2pi)
    double originalOmega = 2.0 * M_PI; // maxtheta era 90 gradi nelle macro originali
    double newOmega      = 2.0 * M_PI * (1.0 - cosMaxTheta);
    double weight        = newOmega / originalOmega;

    particle->SetWeight(weight);
  }
}

} // namespace riptide