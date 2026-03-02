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

#ifndef RIPTIDE_PRIMARY_GENERATOR_ACTION_HPP
#define RIPTIDE_PRIMARY_GENERATOR_ACTION_HPP

#include <G4VUserPrimaryGeneratorAction.hh>

// Forward declarations for compilation speedup
class G4Event;
class G4GeneralParticleSource;

namespace riptide {

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
  // m_prefix is a member variable
  G4GeneralParticleSource* m_gps{nullptr};

 public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event* event) override;
};

} // namespace riptide

#endif // RIPTIDE_PRIMARY_GENERATOR_ACTION_HPP
