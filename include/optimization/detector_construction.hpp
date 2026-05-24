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

#include "importance_sampling.hpp"
#include "lens_cutter.hpp"
#include <filesystem>
#include <optional>

namespace riptide {

class DetectorConstruction : public G4VUserDetectorConstruction {
  G4GDMLParser m_parser{};
  std::filesystem::path m_geometry_path;

  // Posizioni delle lenti (per ottimizzazione)
  double m_l1_x;
  double m_l2_x;

  // Dimensioni e offset (aggiornate in Construct)
  double m_l1_thickness      = 12.5;
  double m_l2_thickness      = 16.3;
  double m_l1_edge_thickness = 1.0;
  double m_l1_center_offset  = -32.35; // (zcut1+zcut2)/2 = (-38.6-26.1)/2
  double m_l2_center_offset  = 22.75;  // (zcut1+zcut2)/2 = (14.6+30.9)/2
  double m_l1_rotation_deg   = 90.0;   // 90 = rot_90y GDML default
  double m_l2_rotation_deg   = 90.0;
  double m_l1_diameter       = 50.8; // Default 2 inch
  double m_l2_diameter       = 50.8;
  double m_l1_R1             = 75.0; // Default
  double m_l1_R2             = 75.0;
  bool m_l1_is_biconvex      = true;

  // IDs per lenti Thorlabs (opzionali)
  std::optional<std::string> m_l1_id;
  std::optional<std::string> m_l2_id;

  // Volumi fisici
  G4VPhysicalVolume* m_world            = nullptr;
  G4VPhysicalVolume* m_l1_phys      = nullptr;
  G4VPhysicalVolume* m_l2_phys      = nullptr;
  G4VPhysicalVolume* m_photocathode_phys = nullptr;

 public:
  DetectorConstruction(std::filesystem::path geometry_path, double l1_x = 83.9,
                       double l2_x = 153.4);

  // Nuovo costruttore con supporto LensCutter
  DetectorConstruction(std::filesystem::path geometry_path, std::string l1_id,
                       std::string l2_id, double l1_x = 83.9, double l2_x = 153.4);
  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;

  // Setter per ottimizzazione
  void SetLensPositions(double l1_x, double l2_x);
  void SetDetectorPosition(double x_det);

  // Nuovo: Cambia i modelli delle lenti al volo
  void SetLenses(const std::string& l1_id, const std::string& l2_id);

  // Getter per dimensioni e posizioni
  double GetL1Thickness() const;
  double GetL2Thickness() const;
  // Getter per debug/verifica rotazione
  double GetL1RotationDeg() const {
    return m_l1_rotation_deg;
  }
  double GetL2RotationDeg() const {
    return m_l2_rotation_deg;
  }

  // Ritorna l'offset del centro del solido rispetto alla posizione x impostata
  // Per i solidi GDML originali, non sono centrati in 0.
  // Per i solidi Thorlabs, sono centrati in 0.
  double GetL1CenterOffset() const;
  double GetL2CenterOffset() const;
  double GetL1X() const {
    return m_l1_x;
  }
  double GetL1Diameter() const {
    return m_l1_diameter;
  }
  double GetL2Diameter() const {
    return m_l2_diameter;
  }

  // Restituisce i parametri completi per l'importance sampling
  ImportanceSamplingHelper::LensParams GetL1Params() const;
  double GetL2X() const {
    return m_l2_x;
  }

  std::string GetL1Id() const {
    return m_l1_id.value_or("LB4553");
  }
  std::string GetL2Id() const {
    return m_l2_id.value_or("LB4592");
  }
};

} // namespace riptide

#endif // RIPTIDE_DETECTOR_CONSTRUCTION_HPP
