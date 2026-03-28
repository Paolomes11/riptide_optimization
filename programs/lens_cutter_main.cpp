#include "lens_cutter.hpp"

#include <iostream>
#include <lyra/lyra.hpp>
#include <spdlog/spdlog.h>

int main(int argc, char** argv) {
    std::filesystem::path data_path = "lens_cutter/lens_data/thorlabs_biconvex.tsv";
    std::string lens_id;
    bool list_lenses = false;
    bool show_help   = false;

    auto cli = lyra::cli() | lyra::help(show_help)
             | lyra::opt(data_path, "data")["-d"]["--data"]("Path to lens data TSV")
             | lyra::opt(lens_id, "id")["-i"]["--id"]("Lens ID to generate GDML for")
             | lyra::opt(list_lenses)["-l"]["--list"]("List available lenses");

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

        if (list_lenses) {
            std::cout << "Available lenses in " << data_path << ":\n";
            std::cout << std::left << std::setw(10) << "ID" << std::setw(10) << "Dia" << std::setw(10)
                      << "f" << std::setw(10) << "R" << std::setw(10) << "tc" << std::setw(10) << "te"
                      << "\n";
            std::cout << std::string(60, '-') << "\n";
            for (const auto& lens : cutter.get_lenses()) {
                std::cout << std::left << std::setw(10) << lens.id << std::setw(10) << lens.diameter
                          << std::setw(10) << lens.focal_length << std::setw(10)
                          << lens.radius_of_curvature << std::setw(10) << lens.center_thickness
                          << std::setw(10) << lens.edge_thickness << "\n";
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
