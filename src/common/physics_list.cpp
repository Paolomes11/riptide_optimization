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

#include "physics_list.hpp"

#include <G4EmStandardPhysics.hh>
#include <G4OpticalParameters.hh>
#include <G4OpticalPhysics.hh>
#include <G4SystemOfUnits.hh>

namespace riptide {

PhysicsList::PhysicsList() {
  G4VModularPhysicsList::defaultCutValue = 1.0 * mm;

  RegisterPhysics(new G4EmStandardPhysics{0});

  RegisterPhysics(new G4OpticalPhysics{0});

  // Disabilita processi ottici non necessari (fotoni generati via GPS)
  auto* optParams = G4OpticalParameters::Instance();
  optParams->SetProcessActivation("Scintillation", false);
  optParams->SetProcessActivation("Cerenkov", false);
  optParams->SetProcessActivation("WLS", false);
  optParams->SetProcessActivation("WLS2", false);
}

} // namespace riptide