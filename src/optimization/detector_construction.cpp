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
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <spdlog/spdlog.h>

namespace riptide {

DetectorConstruction::DetectorConstruction(std::filesystem::path geometry_path, double l1_x,
                                           double l2_x)
    : m_geometry_path{std::move(geometry_path)}
    , m_l1_x{l1_x}
    , m_l2_x{l2_x}
    , m_l1_id{std::nullopt}
    , m_l2_id{std::nullopt} {
  if (!std::filesystem::exists(m_geometry_path)
      || !std::filesystem::is_regular_file(m_geometry_path)) {
    throw std::runtime_error("Geometry path is invalid");
  }
}

DetectorConstruction::DetectorConstruction(std::filesystem::path geometry_path,
                                           std::string l1_id, std::string l2_id,
                                           double l1_x, double l2_x)
    : m_geometry_path{std::move(geometry_path)}
    , m_l1_x{l1_x}
    , m_l2_x{l2_x}
    , m_l1_id{std::move(l1_id)}
    , m_l2_id{std::move(l2_id)} {
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

      if (name == "l1_x_phys") {
        m_l1_phys = pv;
      } else if (name == "l2_x_phys") {
        m_l2_phys = pv;
      } else if (name == "photocathode_x_phys") {
        m_photocathode_phys = pv;
      }
    }

    if (!m_l1_phys || !m_l2_phys) {
      throw std::runtime_error("Lens physical volumes not found in GDML");
    }
    if (!m_photocathode_phys) {
      spdlog::warn("Photocathode physical volume not found in GDML — SD disabled");
    }

    // Supporto per LensCutter: sostituisci i solidi se gli ID sono forniti
    if (m_l1_id || m_l2_id) {
      LensCutter cutter("lens_cutter/lens_data/thorlabs_biconvex.tsv");
      cutter.load_catalog("lens_cutter/lens_data/thorlabs_planoconvex.tsv");

      auto replace_solid = [&](G4VPhysicalVolume* pv, const std::string& id, bool is_l1) {
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

          if (is_l1) {
            m_l1_thickness     = lens->center_thickness;
            m_l1_center_offset = offset;
            m_l1_rotation_deg  = 90.0 + lens->rotation_deg;
            m_l1_diameter      = lens->diameter;
          } else {
            m_l2_thickness     = lens->center_thickness;
            m_l2_center_offset = offset;
            m_l2_rotation_deg  = 90.0 + lens->rotation_deg;
            m_l2_diameter      = lens->diameter;
          }

          spdlog::info("Replaced solid for {} with Thorlabs lens {} (offset: {:.3f} mm)",
                       pv->GetName(), id, offset);
        } else {
          spdlog::warn("Thorlabs lens ID {} not found, keeping GDML solid", id);
        }
      };

      if (m_l1_id)
        replace_solid(m_l1_phys, *m_l1_id, true);
      if (m_l2_id)
        replace_solid(m_l2_phys, *m_l2_id, false);
    }
  }

  // Imposta posizione iniziale
  SetLensPositions(m_l1_x, m_l2_x);

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

void DetectorConstruction::SetLensPositions(double l1_x, double l2_x) {
  m_l1_x = l1_x;
  m_l2_x = l2_x;

  // Controllo collisioni (opzionale, logga un warning)
  double dist               = std::abs(m_l2_x - m_l1_x);
  double sum_half_thickness = (m_l1_thickness + m_l2_thickness) / 2.0;
  if (dist < sum_half_thickness) {
    spdlog::warn("Lenses are overlapping! Distance: {:.3f} mm, sum of half-thicknesses: {:.3f} mm",
                 dist, sum_half_thickness);
  }

  if (m_l1_phys)
    m_l1_phys->SetTranslation(G4ThreeVector(m_l1_x - m_l1_center_offset, 0, 0));
  if (m_l2_phys)
    m_l2_phys->SetTranslation(G4ThreeVector(m_l2_x - m_l2_center_offset, 0, 0));

  // Forza il ricalcolo della geometria per evitare che fotoni
  // vengano tracciati con la geometria della configurazione precedente
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void DetectorConstruction::SetDetectorPosition(double x_det) {
  if (!m_photocathode_phys) {
    spdlog::warn("SetDetectorPosition: m_photocathode_phys is null, skipping");
    return;
  }
  m_photocathode_phys->SetTranslation(G4ThreeVector(x_det, 0, 0));
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  spdlog::info("Detector (photocathode) position set to x = {:.2f} mm", x_det);
}

void DetectorConstruction::SetLenses(const std::string& l1_id, const std::string& l2_id) {
  m_l1_id = l1_id;
  m_l2_id = l2_id;

  LensCutter cutter("lens_cutter/lens_data/thorlabs_biconvex.tsv");
  cutter.load_catalog("lens_cutter/lens_data/thorlabs_planoconvex.tsv");

  auto replace_solid = [&](G4VPhysicalVolume* pv, const std::string& id, bool is_l1) {
    const auto* lens = cutter.get_lens_by_id(id);
    if (lens) {
      auto lv        = pv->GetLogicalVolume();
      auto new_solid = lens->to_g4_solid(is_l1 ? "_75" : "_60");
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

      if (is_l1) {
        m_l1_thickness      = lens->center_thickness;
        m_l1_edge_thickness = lens->edge_thickness;
        m_l1_center_offset  = offset;
        m_l1_rotation_deg   = 90.0 + lens->rotation_deg;
        m_l1_diameter       = lens->diameter;
        m_l1_R1             = lens->radius_of_curvature;
        m_l1_R2 = (lens->type == LensType::Biconvex) ? lens->radius_of_curvature : 1e12; // Flat
        m_l1_is_biconvex = (lens->type == LensType::Biconvex);
      } else {
        m_l2_thickness     = lens->center_thickness;
        m_l2_center_offset = offset;
        m_l2_rotation_deg  = 90.0 + lens->rotation_deg;
        m_l2_diameter      = lens->diameter;
      }
      spdlog::info("Updated solid for {} with Thorlabs lens {} (offset: {:.3f} mm)", pv->GetName(),
                   id, offset);
    } else {
      spdlog::warn("Thorlabs lens ID {} not found", id);
    }
  };

  if (m_l1_phys)
    replace_solid(m_l1_phys, l1_id, true);
  if (m_l2_phys)
    replace_solid(m_l2_phys, l2_id, false);
}

ImportanceSamplingHelper::LensParams DetectorConstruction::GetL1Params() const {
  return {m_l1_x,         m_l1_diameter,       m_l1_R1,         m_l1_R2,
          m_l1_thickness, m_l1_edge_thickness, m_l1_is_biconvex};
}

double DetectorConstruction::GetL1Thickness() const {
  return m_l1_thickness;
}

double DetectorConstruction::GetL2Thickness() const {
  return m_l2_thickness;
}

double DetectorConstruction::GetL1CenterOffset() const {
  return m_l1_center_offset;
}

double DetectorConstruction::GetL2CenterOffset() const {
  return m_l2_center_offset;
}

} // namespace riptide
