#ifndef LENS_CUTTER_HPP
#define LENS_CUTTER_HPP

#include <G4LogicalVolume.hh>
#include <G4VSolid.hh>
#include <filesystem>
#include <string>
#include <vector>

namespace riptide {

// Tipo geometrico della lente
enum class LensType {
  Biconvex,    // simmetrica — nessuna rotazione necessaria
  PlanoConvex, // asimmetrica — rotation_deg specifica il lato esposto
};

struct Lens {
  std::string id;
  LensType type = LensType::Biconvex;
  double diameter;
  double focal_length;
  double radius_of_curvature;
  double center_thickness;
  double edge_thickness;
  double back_focal_length;

  // Rotazione attorno all'asse Y [deg] da applicare al volume fisico.
  // Per lenti biconvesse = 0 (simmetria), per plano-convesse indica quale
  // faccia è rivolta verso la sorgente:
  //   0   → lato curvo verso sorgente (default Thorlabs per ridurre aberrazioni)
  //   180 → lato piano verso sorgente
  double rotation_deg = 0.0;

  // Generates a GDML snippet for this lens (solids and structure)
  std::string to_gdml_solid(const std::string& suffix = "") const;
  std::string to_gdml_structure(const std::string& material = "G4_SILICON_DIOXIDE",
                                const std::string& suffix   = "") const;

  // Generates G4 solids and logical volumes directly
  G4VSolid* to_g4_solid(const std::string& suffix = "") const;
  G4LogicalVolume* to_g4_logical(const std::string& material_name = "G4_SILICON_DIOXIDE",
                                 const std::string& suffix        = "") const;

  // Stringa leggibile del tipo
  std::string type_str() const {
    return (type == LensType::Biconvex) ? "biconvex" : "planoconvex";
  }

  // Offset del centro geometrico rispetto all'origine del solido [mm].
  // Per biconvesse è 0. Per plano-convesse è s/2.
  double get_center_offset() const;
};

class LensCutter {
 public:
  explicit LensCutter(const std::filesystem::path& data_path);

  // Carica un secondo file TSV e aggiunge le lenti al catalogo esistente.
  // Può essere chiamato più volte per cataloghi multipli.
  void load_catalog(const std::filesystem::path& data_path);

  const std::vector<Lens>& get_lenses() const {
    return m_lenses;
  }
  const Lens* get_lens_by_id(const std::string& id) const;

 private:
  std::vector<Lens> m_lenses;
  void load_data(const std::filesystem::path& data_path);
};

} // namespace riptide

#endif // LENS_CUTTER_HPP
