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

#include "detector_construction.hpp"
#include "sensitive_detector.hpp"

#include <G4LogicalVolumeStore.hh>
#include <G4PhysicalVolumeStore.hh>

namespace riptide {

DetectorConstruction::DetectorConstruction(std::filesystem::path geometry_path, double lens75_x,
                                           double lens60_x)
    : m_geometry_path{std::move(geometry_path)}
    , m_lens75_x{lens75_x}
    , m_lens60_x{lens60_x} {
  if (!std::filesystem::exists(m_geometry_path)
      || !std::filesystem::is_regular_file(m_geometry_path)) {
    throw std::runtime_error("Geometry path is invalid");
  }
}

G4VPhysicalVolume* DetectorConstruction::Construct() {
  // Legge GDML solo la prima volta
  if (!m_world) {
    m_parser.Read(m_geometry_path.string());
    m_world = m_parser.GetWorldVolume();

    // Trova i volumi fisici delle lenti
    auto pv_store = G4PhysicalVolumeStore::GetInstance();

    for (auto pv : *pv_store) {
      auto name = pv->GetName();

      if (name == "lens75_x_phys") {
        m_lens75_phys = pv;
      } else if (name == "lens60_x_phys") {
        m_lens60_phys = pv;
      }
    }

    if (!m_lens75_phys || !m_lens60_phys) {
      throw std::runtime_error("Lens physical volumes not found in GDML");
    }
  }

  // Imposta posizione iniziale
  SetLensPositions(m_lens75_x, m_lens60_x);

  return m_world;
}

void DetectorConstruction::ConstructSDandField() {
  auto volume_store = G4LogicalVolumeStore::GetInstance();

  for (auto lv : *volume_store) {
    auto const& aux_list = m_parser.GetVolumeAuxiliaryInformation(lv);
    auto const it        = std::find_if(aux_list.begin(), aux_list.end(),
                                        [](const auto& aux) { return aux.type == "sd_name"; });
    if (it != aux_list.end()) {
      if (it->value.empty()) {
        throw std::runtime_error("Empty sd_name value for volume");
      }

      lv->SetSensitiveDetector(new SensitivePhotocathode{it->value});
    }
  }
}

void DetectorConstruction::SetLensPositions(double lens75_x, double lens60_x) {
  m_lens75_x = lens75_x;
  m_lens60_x = lens60_x;

  if (m_lens75_phys)
    m_lens75_phys->SetTranslation(G4ThreeVector(lens75_x, 0., 0.));

  if (m_lens60_phys)
    m_lens60_phys->SetTranslation(G4ThreeVector(lens60_x, 0., 0.));
}

} // namespace riptide