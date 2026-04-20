/* Copyright 2026 Giulio Mesini — Apache License 2.0 */

#ifndef RIPTIDE_EXP2_ANALYSIS_HPP
#define RIPTIDE_EXP2_ANALYSIS_HPP

#include "stacking.hpp"
#include "trace_extractor.hpp"

#include <filesystem>
#include <limits>
#include <string>
#include <vector>

namespace riptide::exp2 {

/// Etichette delle 4 configurazioni sperimentali (2 qualità × 2 posizioni focali).
enum class ConfigLabel { GoodFocus, GoodNoFocus, BadFocus, BadNoFocus };

/**
 * Restituisce una stringa stabile per ConfigLabel (usata in output filenames e summary).
 *
 * @param label Etichetta configurazione
 * @return      Stringa in snake_case
 */
std::string config_label_str(ConfigLabel label);

/// Parametri ottici associati alla configurazione (per confronto DoF).
struct OpticsParams {
  double x1            = 0.0;
  double x2            = 0.0;
  double x_det         = 0.0;
  double x_det_optimal = 0.0;
};

/// Risultato completo della pipeline per una configurazione.
struct ConfigResult {
  ConfigLabel label;
  OpticsParams optics;
  riptide::stack::StackedImage signal_stack;
  std::vector<double> diff;
  TraceResult trace;
  CentroidFitResult centroid_fit;
  double dof_delta_mm = std::numeric_limits<double>::quiet_NaN();
};

/**
 * Pipeline sequenziale di analisi di una singola configurazione exp2.
 *
 * @param signal_dir     Directory contenente i frame con laser
 * @param background_dir Directory contenente i frame di background
 * @param label          Etichetta configurazione
 * @param optics         Parametri ottici associati
 * @param stack_cfg      Config sigma-clipping (o metodo mean/median)
 * @param trace_cfg      Config estrazione traccia
 * @return               ConfigResult completo
 */
ConfigResult analyze_config(const std::filesystem::path& signal_dir,
                            const std::filesystem::path& background_dir, ConfigLabel label,
                            const OpticsParams& optics,
                            const riptide::stack::StackConfig& stack_cfg = {},
                            const TraceConfig& trace_cfg                 = {});

/// Configurazione output ROOT/PNG.
struct OutputConfig {
  std::filesystem::path output_dir = "output/exp2";
  bool save_png                    = true;
  bool save_root                   = true;
  double z_min_percentile          = 0.005;
  double z_max_percentile          = 0.995;
};

/**
 * Produce output ROOT/PNG separato per una configurazione.
 *
 * @param result Risultato analisi configurazione
 * @param cfg    Configurazione output
 */
void produce_config_output(const ConfigResult& result, const OutputConfig& cfg);

/**
 * Produce pannello summary (numeric comparison) su tutte le configurazioni.
 *
 * @param results Risultati delle 4 configurazioni
 * @param cfg     Configurazione output
 */
void produce_summary(const std::vector<ConfigResult>& results, const OutputConfig& cfg);

} // namespace riptide::exp2

#endif // RIPTIDE_EXP2_ANALYSIS_HPP
