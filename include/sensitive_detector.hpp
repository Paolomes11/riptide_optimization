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

#ifndef RIPTIDE_SENSITIVE_PHOTOCATHODE_HPP
#define RIPTIDE_SENSITIVE_PHOTOCATHODE_HPP

#include <G4VSensitiveDetector.hh>

namespace riptide {

class SensitivePhotocathode : public G4VSensitiveDetector {
 public:
  using G4VSensitiveDetector::G4VSensitiveDetector;
  bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
};

} // namespace riptide

#endif