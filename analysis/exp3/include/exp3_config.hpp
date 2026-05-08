/* Copyright 2026 Giulio Mesini — Apache License 2.0 */

#ifndef RIPTIDE_EXP3_CONFIG_HPP
#define RIPTIDE_EXP3_CONFIG_HPP

#include <filesystem>
#include <string>
#include <vector>

namespace riptide::exp3 {

/// Parametri fisici del display (smartphone).
struct DisplayConfig {
  double mm_per_px_x  = 0.0;
  double mm_per_px_y  = 0.0;
  int    width_px     = 1920;
  int    height_px    = 1080;
  double wavelength_nm = 525.0;
};

/// Configurazione ottica (posizione delle due lenti).
struct LensConfig {
  double x1_mm = 0.0;
  double x2_mm = 0.0;
  std::string label;  // "good" o "bad"
};

/// Configurazione completa dell'esperimento 3.
struct Exp3Config {
  DisplayConfig display;

  /// Distanze assiali nominali [mm] (valore impostato sul banco).
  std::vector<double> axial_distances_nominal_mm;
  /// Distanze assiali misurate col righello [mm].
  std::vector<double> axial_distances_measured_mm;
  /// Orientazioni nominali delle tracce [gradi].
  std::vector<double> orientations_deg;

  /// Passo della griglia di calibrazione [pixel display].
  int calibration_grid_step_px = 50;

  // -- Stacking --
  double stack_n_sigma  = 3.0;
  int    stack_n_iter   = 3;
  int    stack_min_frames = 2;

  // -- Estrazione traccia --
  double min_snr            = 5.0;
  int    min_valid_slices   = 30;
  double angle_tolerance_deg = 10.0;

  // -- Confronto simulazione --
  std::string q_map_tsv;

  /// Configurazioni lenti da analizzare.
  std::vector<LensConfig> lens_configs;
};

/**
 * Carica la configurazione da un file JSON.
 *
 * @param path  Path al file JSON
 * @return      Exp3Config popolato
 * @throws      std::runtime_error se il file non esiste o il JSON è malformato
 */
Exp3Config load_exp3_config(const std::filesystem::path& path);

} // namespace riptide::exp3

#endif // RIPTIDE_EXP3_CONFIG_HPP
