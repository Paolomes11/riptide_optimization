#include "lens_cutter.hpp"

#include <G4Ellipsoid.hh>
#include <G4LogicalVolume.hh>
#include <G4NistManager.hh>
#include <G4RotationMatrix.hh>
#include <G4SystemOfUnits.hh>
#include <G4Tubs.hh>
#include <G4UnionSolid.hh>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace riptide {

LensCutter::LensCutter(const std::filesystem::path& data_path) {
  load_data(data_path);
}

void LensCutter::load_catalog(const std::filesystem::path& data_path) {
  load_data(data_path);
}

void LensCutter::load_data(const std::filesystem::path& data_path) {
  std::ifstream file(data_path);
  if (!file.is_open())
    throw std::runtime_error("Could not open lens data file: " + data_path.string());

  std::string header_line;
  if (!std::getline(file, header_line))
    return;

  // Rileva se il TSV contiene la colonna Rotation_deg (plano-convesse)
  bool has_rotation      = (header_line.find("Rotation_deg") != std::string::npos);
  LensType detected_type = has_rotation ? LensType::PlanoConvex : LensType::Biconvex;

  std::string line;
  while (std::getline(file, line)) {
    if (line.empty())
      continue;
    std::stringstream ss(line);
    Lens lens;
    lens.type = detected_type;

    if (!std::getline(ss, lens.id, '\t'))
      continue;

    auto read_val = [&](double& target) {
      std::string val;
      if (std::getline(ss, val, '\t'))
        target = std::stod(val);
    };

    read_val(lens.diameter);
    read_val(lens.focal_length);
    read_val(lens.radius_of_curvature);
    read_val(lens.center_thickness);
    read_val(lens.edge_thickness);
    read_val(lens.back_focal_length);

    if (has_rotation) {
      read_val(lens.rotation_deg);
    }

    m_lenses.push_back(lens);
    // Dopo m_lenses.push_back(lens):
    if (lens.type == LensType::PlanoConvex) {
      double s = lens.center_thickness - lens.edge_thickness;
      if (s <= 0.0 || s >= lens.radius_of_curvature) {
        throw std::runtime_error("Lente " + lens.id
                                 + ": geometria non valida (s=" + std::to_string(s)
                                 + ", R=" + std::to_string(lens.radius_of_curvature) + ")");
      }
    }
  }
}

const Lens* LensCutter::get_lens_by_id(const std::string& id) const {
  for (const auto& lens : m_lenses) {
    if (lens.id == id)
      return &lens;
  }
  return nullptr;
}

std::string Lens::to_gdml_solid(const std::string& suffix) const {
  if (type == LensType::PlanoConvex) {
    double r     = diameter / 2.0;
    double R     = radius_of_curvature;
    double tc    = center_thickness;
    double te    = edge_thickness;
    double s     = tc - te;
    double shift = te / 2.0 - (R - s);
    if (s <= 0.0 || s >= radius_of_curvature)
      throw std::runtime_error(id + ": parametri geometrici non validi per plano-convessa");

    std::stringstream ss;
    ss << std::fixed << std::setprecision(3);
    ss << "  <tube name=\"" << id << "_pc_tube" << suffix << "\" rmin=\"0\" rmax=\"" << r
       << "\" z=\"" << te << "\" deltaphi=\"2*pi\" startphi=\"0\" aunit=\"rad\" lunit=\"mm\" />\n";

    ss << "  <ellipsoid name=\"" << id << "_pc_cap" << suffix << "\" ax=\"" << R << "\" by=\"" << R
       << "\" cz=\"" << R << "\" zcut1=\"" << (R - s) << "\" zcut2=\"" << R
       << "\" lunit=\"mm\" />\n";

    ss << "  <union name=\"" << id << "_sol" << suffix << "\">\n"
       << "    <first ref=\"" << id << "_pc_tube" << suffix << "\" />\n"
       << "    <second ref=\"" << id << "_pc_cap" << suffix << "\" />\n"
       << "    <position name=\"" << id << "_pos_pc" << suffix << "\" z=\"" << shift
       << "\" lunit=\"mm\" />\n"
       << "  </union>\n";
    return ss.str();
  }

  double r  = diameter / 2.0;
  double R  = radius_of_curvature;
  double tc = center_thickness;
  double te = edge_thickness;
  double s  = (tc - te) / 2.0;

  std::stringstream ss;
  ss << std::fixed << std::setprecision(3);

  // 1. Central Tube
  ss << "  <tube name=\"" << id << "_central_tube" << suffix << "\" rmin=\"0\" rmax=\"" << r
     << "\" z=\"" << te << "\" deltaphi=\"2*pi\" startphi=\"0\" aunit=\"rad\" lunit=\"mm\" />\n";

  // 2. Spherical segment (one side)
  // axes ax, by, cz are all R for a spherical lens
  ss << "  <ellipsoid name=\"" << id << "_spherical_segment" << suffix << "\" ax=\"" << R
     << "\" by=\"" << R << "\" cz=\"" << R << "\" zcut1=\"" << R - s << "\" zcut2=\"" << R
     << "\" lunit=\"mm\" />\n";

  // 3. Union (side 1 + tube)
  // The ellipsoid's center must be shifted so its cut face (at R-s) matches the tube's face (at
  // te/2) Shift = te/2 - (R - s)
  double shift = te / 2.0 - (R - s);
  ss << "  <union name=\"" << id << "_union_1" << suffix << "\">\n";
  ss << "    <first ref=\"" << id << "_central_tube" << suffix << "\" />\n";
  ss << "    <second ref=\"" << id << "_spherical_segment" << suffix << "\" />\n";
  ss << "    <position name=\"" << id << "_pos_1" << suffix << "\" z=\"" << shift
     << "\" lunit=\"mm\" />\n";
  ss << "  </union>\n";

  // 4. Final Union (side 2 + union_1)
  // The other side is identical but flipped and shifted in opposite direction
  // Shift = -(te/2 - (R - s)) = (R - s) - te/2
  double shift2 = -shift;
  ss << "  <union name=\"" << id << "_sol" << suffix << "\">\n";
  ss << "    <first ref=\"" << id << "_union_1" << suffix << "\" />\n";
  ss << "    <second ref=\"" << id << "_spherical_segment" << suffix << "\" />\n";
  ss << "    <position name=\"" << id << "_pos_2" << suffix << "\" z=\"" << shift2
     << "\" lunit=\"mm\" />\n";
  ss << "    <rotation name=\"" << id << "_rot_2" << suffix << "\" x=\"180\" unit=\"deg\" />\n";
  ss << "  </union>\n";

  return ss.str();
}

// ── Geometria plano-convessa ─────────────────────────────────────────────────
//
//  Il solido è composto da:
//    1. Un cilindro piatto di spessore `edge_thickness` e raggio r  (corpo base)
//    2. Una calotta sferica di raggio R che sporge dal lato curvo
//       per un'altezza sagitta s = center_thickness - edge_thickness
//
//  L'asse del solido è Z. Il lato piano è a z = -te/2, il lato curvo a z = +te/2 + s.
//  La rotazione rotation_deg (attorno a Y) è applicata dal chiamante al volume fisico.

static G4VSolid* make_planoconvex_solid(const riptide::Lens& lens, const std::string& suffix) {
  double r  = lens.diameter / 2.0 * mm;
  double R  = lens.radius_of_curvature * mm;
  double tc = lens.center_thickness * mm;
  double te = lens.edge_thickness * mm;
  double s  = tc - te; // sagitta

  // Validazione geometrica
  if (s <= 0.0)
    throw std::runtime_error(lens.id + ": center_thickness <= edge_thickness");
  if (s >= R)
    throw std::runtime_error(lens.id + ": sagitta s=" + std::to_string(s / mm)
                             + " >= R=" + std::to_string(R / mm) + " mm — geometria impossibile");

  // Verifica che la calotta non sporga oltre il raggio della lente:
  // la corda alla base della calotta è 2*sqrt(R^2 - (R-s)^2) = 2*sqrt(2Rs - s^2)
  // deve essere <= diametro
  double chord = 2.0 * std::sqrt(2.0 * R * s - s * s);
  if (chord > 2.0 * r)
    throw std::runtime_error(lens.id + ": calotta (chord=" + std::to_string(chord / mm)
                             + " mm) sporge oltre il diametro (" + std::to_string(2.0 * r / mm)
                             + " mm)");

  // 1. Cilindro base
  auto tube = new G4Tubs(lens.id + "_pc_tube" + suffix, 0., r, te / 2.0, 0., 360. * deg);

  // 2. Calotta: l'ellissoide è centrato nell'origine della sfera.
  //    zcut1 = R - s  (faccia piatta della calotta, coincide con la faccia del cilindro)
  //    zcut2 = R      (apice della calotta)
  //    Shift per allineare zcut1 con la faccia superiore del cilindro (z = +te/2):
  //      z_centro_sfera = te/2 - (R - s)   [può essere negativo per lenti sottili]
  double shift = te / 2.0 - (R - s);
  auto cap     = new G4Ellipsoid(lens.id + "_pc_cap" + suffix, R, R, R, (R - s), R);

  return new G4UnionSolid(lens.id + "_sol" + suffix, tube, cap, nullptr,
                          G4ThreeVector(0., 0., shift));
}

std::string Lens::to_gdml_structure(const std::string& material, const std::string& suffix) const {
  std::stringstream ss;
  ss << "  <volume name=\"" << id << "_log" << suffix << "\">\n";
  ss << "    <materialref ref=\"" << material << "\" />\n";
  ss << "    <solidref ref=\"" << id << "_sol" << suffix << "\" />\n";
  ss << "  </volume>\n";
  return ss.str();
}

G4VSolid* Lens::to_g4_solid(const std::string& suffix) const {
  if (type == LensType::PlanoConvex)
    return make_planoconvex_solid(*this, suffix);

  // ── Biconvessa (codice esistente, invariato) ─────────────────────────────
  double r  = diameter / 2.0 * mm;
  double R  = radius_of_curvature * mm;
  double tc = center_thickness * mm;
  double te = edge_thickness * mm;
  double s  = (tc - te) / 2.0;

  auto tube      = new G4Tubs(id + "_central_tube" + suffix, 0., r, te / 2.0, 0., 360. * deg);
  auto ellipsoid = new G4Ellipsoid(id + "_spherical_segment" + suffix, R, R, R, R - s, R);

  double shift = te / 2.0 - (R - s);
  auto union1  = new G4UnionSolid(id + "_union_1" + suffix, tube, ellipsoid, nullptr,
                                  G4ThreeVector(0, 0, shift));

  double shift2 = -shift;
  G4RotationMatrix rot;
  rot.rotateX(180. * deg);
  return new G4UnionSolid(id + "_sol" + suffix, union1, ellipsoid,
                          G4Transform3D(rot, G4ThreeVector(0, 0, shift2)));
}

double Lens::get_center_offset() const {
  if (type == LensType::Biconvex)
    return 0.0;
  double s = center_thickness - edge_thickness;
  return s / 2.0;
}

G4LogicalVolume* Lens::to_g4_logical(const std::string& material_name,
                                     const std::string& suffix) const {
  auto nist     = G4NistManager::Instance();
  auto material = nist->FindOrBuildMaterial(material_name);
  if (!material) {
    throw std::runtime_error("Could not find material: " + material_name);
  }

  auto solid = to_g4_solid(suffix);
  return new G4LogicalVolume(solid, material, id + "_log" + suffix);
}

} // namespace riptide
