#include "lens_cutter.hpp"

#include <G4Tubs.hh>
#include <G4Ellipsoid.hh>
#include <G4UnionSolid.hh>
#include <G4NistManager.hh>
#include <G4LogicalVolume.hh>
#include <G4SystemOfUnits.hh>
#include <G4RotationMatrix.hh>
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

void LensCutter::load_data(const std::filesystem::path& data_path) {
    std::ifstream file(data_path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open lens data file: " + data_path.string());
    }

    std::string line;
    // Skip header
    if (std::getline(file, line)) {
        // First line is header
    }

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        Lens lens;
        if (std::getline(ss, lens.id, '\t')) {
            std::string val;
            auto read_val = [&](double& target) {
                if (std::getline(ss, val, '\t')) {
                    target = std::stod(val);
                }
            };
            read_val(lens.diameter);
            read_val(lens.focal_length);
            read_val(lens.radius_of_curvature);
            read_val(lens.center_thickness);
            read_val(lens.edge_thickness);
            read_val(lens.back_focal_length);
            m_lenses.push_back(lens);
        }
    }
}

const Lens* LensCutter::get_lens_by_id(const std::string& id) const {
    for (const auto& lens : m_lenses) {
        if (lens.id == id) return &lens;
    }
    return nullptr;
}

std::string Lens::to_gdml_solid(const std::string& suffix) const {
    double r = diameter / 2.0;
    double R = radius_of_curvature;
    double tc = center_thickness;
    double te = edge_thickness;
    double s = (tc - te) / 2.0;

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
    // The ellipsoid's center must be shifted so its cut face (at R-s) matches the tube's face (at te/2)
    // Shift = te/2 - (R - s)
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

std::string Lens::to_gdml_structure(const std::string& material, const std::string& suffix) const {
    std::stringstream ss;
    ss << "  <volume name=\"" << id << "_log" << suffix << "\">\n";
    ss << "    <materialref ref=\"" << material << "\" />\n";
    ss << "    <solidref ref=\"" << id << "_sol" << suffix << "\" />\n";
    ss << "  </volume>\n";
    return ss.str();
}

G4VSolid* Lens::to_g4_solid(const std::string& suffix) const {
    double r  = diameter / 2.0 * mm;
    double R  = radius_of_curvature * mm;
    double tc = center_thickness * mm;
    double te = edge_thickness * mm;
    double s  = (tc - te) / 2.0;

    // 1. Central Tube
    auto tube = new G4Tubs(id + "_central_tube" + suffix, 0., r, te / 2.0, 0., 360. * deg);

    // 2. Spherical segment
    auto ellipsoid = new G4Ellipsoid(id + "_spherical_segment" + suffix, R, R, R, R - s, R);

    // 3. Union 1
    double shift = te / 2.0 - (R - s);
    auto union1  = new G4UnionSolid(id + "_union_1" + suffix, tube, ellipsoid, nullptr,
                                   G4ThreeVector(0, 0, shift));

    // 4. Final Union
    double shift2 = -shift;
    G4RotationMatrix rot;
    rot.rotateX(180. * deg);
    auto final_sol = new G4UnionSolid(id + "_sol" + suffix, union1, ellipsoid,
                                     G4Transform3D(rot, G4ThreeVector(0, 0, shift2)));

    return final_sol;
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
