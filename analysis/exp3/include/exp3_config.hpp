/* Copyright 2026 Giulio Mesini — Apache License 2.0 */

#ifndef RIPTIDE_EXP3_CONFIG_HPP
#define RIPTIDE_EXP3_CONFIG_HPP

#include <array>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>

namespace riptide::exp3 {

/// Parametri fisici del display (smartphone).
struct DisplayConfig {
  double mm_per_px_x   = 0.0;
  double mm_per_px_y   = 0.0;
  int    width_px      = 1080;
  int    height_px     = 2340;
  double wavelength_nm = 525.0;
};

/// Parametri di stacking sigma-clip.
struct StackingConfig {
  double n_sigma    = 3.0;
  int    n_iter     = 3;
  int    min_frames = 2;
};

/// Soglie per l'estrazione della traccia.
struct TraceExtractionConfig {
  double min_snr               = 5.0;
  int    min_valid_slices      = 20;
  int    max_lines             = 3;
  int    half_width_px         = 40;
  int    line_min_separation_px = 100;
};

/// Confronto con la simulazione (buona vs cattiva configurazione lenti).
struct QComparisonConfig {
  std::string q_map_tsv;
  double good_x1_mm = 0.0;
  double good_x2_mm = 0.0;
  double bad_x1_mm  = 0.0;
  double bad_x2_mm  = 0.0;

  /// Q_sim di riferimento a fuoco fisso (banco exp3), non ricavabili da
  /// q_map_tsv perche' quest'ultimo proviene da run a fuoco mobile. Se
  /// finiti, hanno precedenza sul lookup nearest-neighbor in q_map_tsv.
  double q_sim_good = std::numeric_limits<double>::quiet_NaN();
  double q_sim_bad  = std::numeric_limits<double>::quiet_NaN();
};

/// Configurazione completa dell'esperimento 3.
struct Exp3Config {
  DisplayConfig display;

  /// Distanze assiali nominali [mm] (valore impostato sul banco).
  std::vector<double> axial_distances_nominal_mm;
  /// Distanze assiali misurate col righello [mm].
  std::vector<double> axial_distances_measured_mm;

  /// Offset radiali delle tre tracce parallele [pixel display].
  std::vector<int> radial_offsets_px;
  /// Passo della colonna di dot per Exp3b [pixel display].
  int dot_column_step_px = 50;
  /// Lunghezza della linea di calibrazione [pixel display].
  int calib_line_length_px = 800;
  /// Centro dell'asse ottico nel display [pixel display].
  std::array<double, 2> optical_axis_center_px = {540.0, 1170.0};

  StackingConfig        stacking;
  TraceExtractionConfig trace_extraction;
  QComparisonConfig     q_comparison;
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
