#ifndef LENS_CUTTER_HPP
#define LENS_CUTTER_HPP

#include <G4VSolid.hh>
#include <G4LogicalVolume.hh>
#include <filesystem>
#include <string>
#include <vector>

namespace riptide {

struct Lens {
    std::string id;
    double diameter;
    double focal_length;
    double radius_of_curvature;
    double center_thickness;
    double edge_thickness;
    double back_focal_length;

    // Generates a GDML snippet for this lens (solids and structure)
    std::string to_gdml_solid(const std::string& suffix = "") const;
    std::string to_gdml_structure(const std::string& material = "G4_SILICON_DIOXIDE",
                                  const std::string& suffix = "") const;

    // Generates G4 solids and logical volumes directly
    G4VSolid* to_g4_solid(const std::string& suffix = "") const;
    G4LogicalVolume* to_g4_logical(const std::string& material_name = "G4_SILICON_DIOXIDE",
                                   const std::string& suffix = "") const;
};

class LensCutter {
public:
    explicit LensCutter(const std::filesystem::path& data_path);

    const std::vector<Lens>& get_lenses() const { return m_lenses; }
    const Lens* get_lens_by_id(const std::string& id) const;

private:
    std::vector<Lens> m_lenses;
    void load_data(const std::filesystem::path& data_path);
};

} // namespace riptide

#endif // LENS_CUTTER_HPP
