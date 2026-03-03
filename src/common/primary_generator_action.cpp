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

#include <G4GeneralParticleSource.hh>

namespace riptide {

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : m_gps{new G4GeneralParticleSource()} {
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  delete m_gps;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  // Genera il primo vertice dell'interazione particella
  m_gps->GeneratePrimaryVertex(event);
}

} // namespace riptide