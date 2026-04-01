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

#include "lens_cutter.hpp"
#include <filesystem>
#include <optional>

namespace riptide {

class DetectorConstruction : public G4VUserDetectorConstruction {
  G4GDMLParser m_parser{};
  std::filesystem::path m_geometry_path;

  // Posizioni delle lenti (per ottimizzazione)
  double m_lens75_x;
  double m_lens60_x;

  // Dimensioni e offset (aggiornate in Construct)
  double m_lens75_thickness     = 12.5;
  double m_lens60_thickness     = 16.3;
  double m_lens75_center_offset = -32.35; // (zcut1+zcut2)/2 = (-38.6-26.1)/2
  double m_lens60_center_offset = 22.75;  // (zcut1+zcut2)/2 = (14.6+30.9)/2
  double m_lens75_rotation_deg  = 90.0;   // 90 = rot_90y GDML default
  double m_lens60_rotation_deg  = 90.0;

  // IDs per lenti Thorlabs (opzionali)
  std::optional<std::string> m_lens75_id;
  std::optional<std::string> m_lens60_id;

  // Volumi fisici
  G4VPhysicalVolume* m_world       = nullptr;
  G4VPhysicalVolume* m_lens75_phys = nullptr;
  G4VPhysicalVolume* m_lens60_phys = nullptr;

 public:
  DetectorConstruction(std::filesystem::path geometry_path, double lens75_x = 83.9,
                       double lens60_x = 153.4);

  // Nuovo costruttore con supporto LensCutter
  DetectorConstruction(std::filesystem::path geometry_path, std::string lens75_id,
                       std::string lens60_id, double lens75_x = 83.9, double lens60_x = 153.4);
  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;

  // Setter per ottimizzazione
  void SetLensPositions(double lens75_x, double lens60_x);

  // Nuovo: Cambia i modelli delle lenti al volo
  void SetLenses(const std::string& lens75_id, const std::string& lens60_id);

  // Getter per dimensioni e posizioni
  double GetLens75Thickness() const;
  double GetLens60Thickness() const;
  // Getter per debug/verifica rotazione
  double GetLens75RotationDeg() const {
    return m_lens75_rotation_deg;
  }
  double GetLens60RotationDeg() const {
    return m_lens60_rotation_deg;
  }

  // Ritorna l'offset del centro del solido rispetto alla posizione x impostata
  // Per i solidi GDML originali, non sono centrati in 0.
  // Per i solidi Thorlabs, sono centrati in 0.
  double GetLens75CenterOffset() const;
  double GetLens60CenterOffset() const;

  double GetLens75X() const {
    return m_lens75_x;
  }
  double GetLens60X() const {
    return m_lens60_x;
  }

  std::string GetLens75Id() const {
    return m_lens75_id.value_or("LB4553");
  }
  std::string GetLens60Id() const {
    return m_lens60_id.value_or("LB4592");
  }
};

} // namespace riptide

#endif // RIPTIDE_DETECTOR_CONSTRUCTION_HPP