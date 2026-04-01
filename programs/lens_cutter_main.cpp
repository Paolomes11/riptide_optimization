#include "lens_cutter.hpp"

#include <lyra/lyra.hpp>
#include <spdlog/spdlog.h>
#include <iostream>

int main(int argc, char** argv) {
  std::filesystem::path data_path = "lens_cutter/lens_data/thorlabs_biconvex.tsv";
  std::string lens_id;
  bool list_lenses = false;
  bool show_help   = false;
  std::string extra_catalog;

  auto cli =
      lyra::cli() | lyra::help(show_help)
      | lyra::opt(data_path, "data")["-d"]["--data"]("Path to lens data TSV")
      | lyra::opt(lens_id, "id")["-i"]["--id"]("Lens ID to generate GDML for")
      | lyra::opt(list_lenses)["-l"]["--list"]("List available lenses")
      | lyra::opt(extra_catalog, "catalog")["--catalog"]("Path to additional lens catalog TSV");

  auto result = cli.parse({argc, argv});
  if (!result) {
    spdlog::error("CLI parsing error: {}", result.message());
    std::cerr << cli << '\n';
    return EXIT_FAILURE;
  }

  if (show_help) {
    std::cout << cli << '\n';
    return EXIT_SUCCESS;
  }

  try {
    riptide::LensCutter cutter(data_path);

    if (!extra_catalog.empty()) {
      cutter.load_catalog(extra_catalog);
    }

    if (list_lenses) {
      std::cout << "Available lenses in " << data_path << ":\n";
      std::cout << std::left << std::setw(12) << "ID" << std::setw(14) << "Type" << std::setw(8)
                << "Dia" << std::setw(8) << "f" << std::setw(8) << "R" << std::setw(8) << "tc"
                << std::setw(8) << "te" << std::setw(10) << "Rot[deg]"
                << "\n";
      std::cout << std::string(76, '-') << "\n";
      for (const auto& lens : cutter.get_lenses()) {
        std::cout << std::left << std::setw(12) << lens.id << std::setw(14) << lens.type_str()
                  << std::setw(8) << lens.diameter << std::setw(8) << lens.focal_length
                  << std::setw(8) << lens.radius_of_curvature << std::setw(8)
                  << lens.center_thickness << std::setw(8) << lens.edge_thickness << std::setw(10)
                  << lens.rotation_deg << "\n";
      }
      return EXIT_SUCCESS;
    }

    if (!lens_id.empty()) {
      const auto* lens = cutter.get_lens_by_id(lens_id);
      if (!lens) {
        spdlog::error("Lens with ID {} not found", lens_id);
        return EXIT_FAILURE;
      }

      std::cout << "<!-- GDML Solids for " << lens_id << " -->\n";
      std::cout << lens->to_gdml_solid();
      std::cout << "\n<!-- GDML Structure for " << lens_id << " -->\n";
      std::cout << lens->to_gdml_structure();
      return EXIT_SUCCESS;
    }

    std::cout << "Use --list to see available lenses or --id <ID> to generate GDML.\n";

  } catch (const std::exception& e) {
    spdlog::error("Error: {}", e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
