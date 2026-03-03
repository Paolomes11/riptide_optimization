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

#ifndef RIPTIDE_DETECTOR_CONSTRUCTION_HPP
#define RIPTIDE_DETECTOR_CONSTRUCTION_HPP

#include <G4GDMLParser.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VUserDetectorConstruction.hh>

#include <filesystem>

namespace riptide {

class DetectorConstruction : public G4VUserDetectorConstruction {
  G4GDMLParser m_parser{};
  std::filesystem::path m_geometry_path;

  // Posizioni delle lenti (per ottimizzazione)
  double m_lens75_x;
  double m_lens60_x;

  // Volumi fisici
  G4VPhysicalVolume* m_world       = nullptr;
  G4VPhysicalVolume* m_lens75_phys = nullptr;
  G4VPhysicalVolume* m_lens60_phys = nullptr;

 public:
  DetectorConstruction(std::filesystem::path geometry_path, double lens75_x = 83.9,
                       double lens60_x = 153.4);
  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;

  // Setter per ottimizzazione
  void SetLensPositions(double lens75_x, double lens60_x);

  // Getter per ottimizzazione
  double GetLens75X() const {
    return m_lens75_x;
  }
  double GetLens60X() const {
    return m_lens60_x;
  }
};

} // namespace riptide

#endif // RIPTIDE_DETECTOR_CONSTRUCTION_HPP