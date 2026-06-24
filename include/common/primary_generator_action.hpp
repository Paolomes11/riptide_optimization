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
#include <G4ThreeVector.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <functional>
#include <vector>

// Forward declarations for compilation speedup
class G4Event;
class G4GeneralParticleSource;

namespace riptide {

struct SpotConfig {
  G4ThreeVector pos;       // posizione sorgente in Geant4 units (mm)
  G4ThreeVector is_axis;   // asse cono IS
  double is_theta  = 0.0;  // angolo massimo cono IS
  double is_weight = 1.0;  // peso IS per colonna ntuple (1-cos(theta_max))
  bool use_is      = false;
};

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
    m_staticAxis            = axis;
    m_staticMaxTheta        = maxTheta;
    m_hasStaticCone         = true;
    m_useImportanceSampling = true;
  }

  // Cycling: genera fotoni uno per spot ciclicamente su tutti gli spot.
  // Il BeamOn deve usare n_spots * n_photons_per_spot come conteggio.
  void ConfigureSpotCycling(const std::vector<SpotConfig>& spots);
  void ResetCycling();     // azzera m_currentSpotIdx (chiamare prima di BeamOn)
  void DisableCycling();   // torna alla modalità legacy

 private:
  G4GeneralParticleSource* m_gps{nullptr};

  bool m_useImportanceSampling = false;
  std::function<ImportanceSamplingHelper::LensParams()> m_paramsProvider;

  bool m_hasStaticCone = false;
  G4ThreeVector m_staticAxis;
  double m_staticMaxTheta = 0.0;

  // Cycling spot
  std::vector<SpotConfig> m_spotConfigs;
  int m_currentSpotIdx = 0;
};

} // namespace riptide

#endif // RIPTIDE_PRIMARY_GENERATOR_ACTION_HPP
