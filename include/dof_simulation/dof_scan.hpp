#ifndef RIPTIDE_DOF_SCAN_HPP
#define RIPTIDE_DOF_SCAN_HPP

#include <filesystem>
#include <string>

class G4RunManager;

namespace riptide {

void run_dof_scan(G4RunManager* run_manager, const std::filesystem::path& macro_file,
                  const std::string& root_output_file,
                  const std::filesystem::path& config_file = "config/config.json",
                  bool all_lenses = false, const std::string& lens75_id = "",
                  const std::string& lens60_id = "");

} // namespace riptide

#endif // RIPTIDE_DOF_SCAN_HPP
