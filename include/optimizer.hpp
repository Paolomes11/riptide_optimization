#ifndef RIPTIDE_OPTIMIZER_HPP
#define RIPTIDE_OPTIMIZER_HPP

#include <filesystem>

class G4RunManager;
class EfficiencyCollector;

namespace riptide {
void run_optimization(G4RunManager* run_manager, const std::filesystem::path& macro_file,
                      EfficiencyCollector* collector);
}

#endif // RIPTIDE_OPTIMIZER_HPP