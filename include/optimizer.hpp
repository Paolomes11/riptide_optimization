#ifndef RIPTIDE_OPTIMIZER_HPP
#define RIPTIDE_OPTIMIZER_HPP

#include <filesystem>

class G4RunManager;

namespace riptide {
void run_optimization(G4RunManager* run_manager, const std::filesystem::path& macro_file);
}

#endif // RIPTIDE_OPTIMIZER_HPP