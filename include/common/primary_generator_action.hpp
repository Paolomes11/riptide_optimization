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

#include "importance_sampling.hpp"
#include <G4VUserPrimaryGeneratorAction.hh>
#include <functional>

// Forward declarations for compilation speedup
class G4Event;
class G4GeneralParticleSource;

namespace riptide {

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event* event) override;

  void SetImportanceSampling(bool enable,
                             std::function<ImportanceSamplingHelper::LensParams()> provider) {
    m_useImportanceSampling = enable;
    m_paramsProvider        = provider;
  }

  void SetStaticCone(const G4ThreeVector& axis, double maxTheta) {
    m_staticAxis     = axis;
    m_staticMaxTheta = maxTheta;
    m_hasStaticCone  = true;
  }

 private:
  G4GeneralParticleSource* m_gps{nullptr};

  bool m_useImportanceSampling = false;
  std::function<ImportanceSamplingHelper::LensParams()> m_paramsProvider;

  bool m_hasStaticCone = false;
  G4ThreeVector m_staticAxis;
  double m_staticMaxTheta = 0.0;
};

} // namespace riptide

#endif // RIPTIDE_PRIMARY_GENERATOR_ACTION_HPP
