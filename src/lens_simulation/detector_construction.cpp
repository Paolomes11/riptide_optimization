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
#include <G4SystemOfUnits.hh>
#include <spdlog/spdlog.h>

namespace riptide {

DetectorConstruction::DetectorConstruction(std::filesystem::path geometry_path, double lens75_x,
                                           double lens60_x)
    : m_geometry_path{std::move(geometry_path)}
    , m_lens75_x{lens75_x}
    , m_lens60_x{lens60_x}
    , m_lens75_id{std::nullopt}
    , m_lens60_id{std::nullopt} {
  if (!std::filesystem::exists(m_geometry_path)
      || !std::filesystem::is_regular_file(m_geometry_path)) {
    throw std::runtime_error("Geometry path is invalid");
  }
}

DetectorConstruction::DetectorConstruction(std::filesystem::path geometry_path,
                                           std::string lens75_id, std::string lens60_id,
                                           double lens75_x, double lens60_x)
    : m_geometry_path{std::move(geometry_path)}
    , m_lens75_x{lens75_x}
    , m_lens60_x{lens60_x}
    , m_lens75_id{std::move(lens75_id)}
    , m_lens60_id{std::move(lens60_id)} {
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

    // Supporto per LensCutter: sostituisci i solidi se gli ID sono forniti
    if (m_lens75_id || m_lens60_id) {
      LensCutter cutter("lens_cutter/lens_data/thorlabs_biconvex.tsv");
      cutter.load_catalog("lens_cutter/lens_data/thorlabs_planoconvex.tsv");

      auto replace_solid = [&](G4VPhysicalVolume* pv, const std::string& id, bool is_lens75) {
        const auto* lens = cutter.get_lens_by_id(id);
        if (lens) {
          auto lv        = pv->GetLogicalVolume();
          auto new_solid = lens->to_g4_solid("_custom");
          lv->SetSolid(new_solid);

          // Convenzione Thorlabs: rotation=0 -> lato curvo verso la sorgente (-X)
          // La nostra rotazione base 90y orienta Z local verso +X global.
          // Quindi per avere il lato curvo (Z>0) verso la sorgente (-X),
          // la rotazione di default deve essere 180 rispetto a quella base.

          auto* rot = new G4RotationMatrix();
          rot->rotateY(90. * CLHEP::deg); // asse Z -> X

          // Se rotation_deg = 0 (Thorlabs), vogliamo lato curvo verso sorgente (-X)
          // Questo significa ruotare di altri 180 gradi.
          double final_rot = 180.0 + lens->rotation_deg;
          rot->rotateY(final_rot * CLHEP::deg);
          pv->SetRotation(rot);

          // Calcola l'offset del centro geometrico nel sistema globale X
          // Dopo rotY(90) e rotY(180+rot), il local Z maps a:
          // rot=0   -> local Z maps to -X
          // rot=180 -> local Z maps to +X
          double sign   = (std::abs(lens->rotation_deg) < 1e-6) ? 1.0 : -1.0;
          double offset = sign * lens->get_center_offset();

          if (is_lens75) {
            m_lens75_thickness     = lens->center_thickness;
            m_lens75_center_offset = offset;
          } else {
            m_lens60_thickness     = lens->center_thickness;
            m_lens60_center_offset = offset;
          }

          spdlog::info("Replaced solid for {} with Thorlabs lens {} (offset: {:.3f} mm)",
                       pv->GetName(), id, offset);
        } else {
          spdlog::warn("Thorlabs lens ID {} not found, keeping GDML solid", id);
        }
      };

      if (m_lens75_id)
        replace_solid(m_lens75_phys, *m_lens75_id, true);
      if (m_lens60_id)
        replace_solid(m_lens60_phys, *m_lens60_id, false);
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

  // Controllo collisioni (opzionale, logga un warning)
  double dist               = std::abs(m_lens60_x - m_lens75_x);
  double sum_half_thickness = (m_lens75_thickness + m_lens60_thickness) / 2.0;
  if (dist < sum_half_thickness) {
    spdlog::warn("Lenses are overlapping! Distance: {:.3f} mm, sum of half-thicknesses: {:.3f} mm",
                 dist, sum_half_thickness);
  }

  if (m_lens75_phys) {
    m_lens75_phys->SetTranslation(G4ThreeVector(m_lens75_x - m_lens75_center_offset, 0, 0));
  }
  if (m_lens60_phys) {
    m_lens60_phys->SetTranslation(G4ThreeVector(m_lens60_x - m_lens60_center_offset, 0, 0));
  }
}

void DetectorConstruction::SetLenses(const std::string& lens75_id, const std::string& lens60_id) {
  m_lens75_id = lens75_id;
  m_lens60_id = lens60_id;

  LensCutter cutter("lens_cutter/lens_data/thorlabs_biconvex.tsv");
  cutter.load_catalog("lens_cutter/lens_data/thorlabs_planoconvex.tsv");

  auto replace_solid = [&](G4VPhysicalVolume* pv, const std::string& id, bool is_lens75) {
    if (!pv) {
      spdlog::warn("Physical volume is null");
      return;
    }
    const auto* lens = cutter.get_lens_by_id(id);
    if (lens) {
      auto lv        = pv->GetLogicalVolume();
      auto new_solid = lens->to_g4_solid(is_lens75 ? "_75" : "_60");
      lv->SetSolid(new_solid);

      // Convenzione Thorlabs: rotation=0 -> lato curvo verso la sorgente (-X)
      // La nostra rotazione base 90y orienta Z local verso +X global.
      // Quindi per avere il lato curvo (Z>0) verso la sorgente (-X),
      // la rotazione di default deve essere 180 rispetto a quella base.

      auto* rot = new G4RotationMatrix();
      rot->rotateY(90. * CLHEP::deg); // asse Z -> X

      // Se rotation_deg = 0 (Thorlabs), vogliamo lato curvo verso sorgente (-X)
      // Questo significa ruotare di altri 180 gradi.
      double final_rot = 180.0 + lens->rotation_deg;
      rot->rotateY(final_rot * CLHEP::deg);
      pv->SetRotation(rot);

      // Calcola l'offset del centro geometrico nel sistema globale X
      // Dopo rotY(90) e rotY(180+rot), il local Z maps a:
      // rot=0   -> local Z maps to -X
      // rot=180 -> local Z maps to +X
      double sign   = (std::abs(lens->rotation_deg) < 1e-6) ? 1.0 : -1.0;
      double offset = sign * lens->get_center_offset();

      if (is_lens75) {
        m_lens75_thickness      = lens->center_thickness;
        m_lens75_edge_thickness = lens->edge_thickness;
        m_lens75_center_offset  = offset;
        m_lens75_diameter       = lens->diameter;
        m_lens75_R1             = lens->radius_of_curvature;
        m_lens75_R2 = (lens->type == LensType::Biconvex) ? lens->radius_of_curvature : 1e12; // Flat
        m_lens75_is_biconvex = (lens->type == LensType::Biconvex);
      } else {
        m_lens60_thickness     = lens->center_thickness;
        m_lens60_center_offset = offset;
      }
      spdlog::info("Updated solid for {} with Thorlabs lens {} (offset: {:.3f} mm)", pv->GetName(),
                   id, offset);
    } else {
      spdlog::warn("Thorlabs lens ID {} not found", id);
    }
  };

  if (m_lens75_phys)
    replace_solid(m_lens75_phys, lens75_id, true);
  if (m_lens60_phys)
    replace_solid(m_lens60_phys, lens60_id, false);
}

ImportanceSamplingHelper::LensParams DetectorConstruction::GetLens75Params() const {
  return {m_lens75_x,         m_lens75_diameter,       m_lens75_R1,         m_lens75_R2,
          m_lens75_thickness, m_lens75_edge_thickness, m_lens75_is_biconvex};
}

double DetectorConstruction::GetLens75Thickness() const {
  return m_lens75_thickness;
}

double DetectorConstruction::GetLens60Thickness() const {
  return m_lens60_thickness;
}

double DetectorConstruction::GetLens75CenterOffset() const {
  return m_lens75_center_offset;
}

double DetectorConstruction::GetLens60CenterOffset() const {
  return m_lens60_center_offset;
}

} // namespace riptide