/* Copyright 2026 Giulio Mesini — Apache License 2.0 */

#ifndef RIPTIDE_EXP3B_ANALYSIS_HPP
#define RIPTIDE_EXP3B_ANALYSIS_HPP

#include "exp3_analysis.hpp"   // LineCalib, Exp3Config

#include <filesystem>
#include <string>
#include <vector>

namespace riptide::exp3b {

using riptide::exp3::Exp3Config;
using riptide::exp3::LineCalib;

/// Centroide di un singolo dot nella colonna.
struct DotCentroid {
  double y_disp_mm = 0.0; ///< posizione display [mm] dal bordo superiore
  double y_sens_px = 0.0; ///< centroide sensore [px]
  double y_sens_mm = 0.0; ///< centroide sensore [mm]
  double sigma_px  = 0.5; ///< incertezza fit gaussiano [px]
};

/// Risultato della misura di M(r) per una distanza assiale.
struct MResult {
  double d_ax_mm      = 0.0;
  double M_global     = 0.0; ///< pendenza fit lineare globale
  double q_offset     = 0.0; ///< intercetta
  double chi2_ndof    = 0.0;
  double M_rms_residual = 0.0;
  int    n_dots       = 0;
  std::vector<double> r_mm;    ///< posizioni y_sens_mm dei dot
  std::vector<double> M_local; ///< M incrementale tra dot adiacenti
};

/**
 * Pipeline Exp3b: misura di M(r) tramite colonna di dot.
 *
 * 1. Stack sigma-clip + sottrazione background
 * 2. Fit gaussiano 2D su ogni dot (seed da posizione teorica)
 * 3. Fit lineare globale y_sens vs y_disp → M_global
 * 4. M_local[i] = (y_sens[i+1]-y_sens[i]) / (y_disp[i+1]-y_disp[i])
 *
 * @param signal_dir  Cartella frame segnale FITS (colonna di dot)
 * @param bg_dir      Cartella frame background FITS
 * @param cfg         Configurazione
 * @param calib       Calibrazione (scale_mm_per_sens_px)
 * @param d_ax_mm     Distanza assiale misurata [mm]
 * @return            MResult con M_global, M_local, chi2_ndof
 */
MResult run_exp3b(const std::filesystem::path& signal_dir,
                  const std::filesystem::path& bg_dir,
                  const Exp3Config& cfg,
                  const LineCalib& calib,
                  double d_ax_mm);

/**
 * Scrive M_summary.tsv.
 *
 * Colonne: config  d_ax_mm  M_global  M_rms_residual  chi2_ndof  n_dots
 */
void write_M_summary_tsv(const std::string& config,
                          const std::vector<MResult>& results,
                          const std::filesystem::path& path);

// ---------------------------------------------------------------------------
// Plot (implementate in exp3_plots.cpp)
// ---------------------------------------------------------------------------

/**
 * Produce M_profile_d{dist}.png:
 * y_sens vs y_disp con fit lineare + pannello residui.
 */
void produce_M_profile(const MResult& result,
                       const std::string& config_label,
                       const std::filesystem::path& output_path);

/**
 * Produce M_vs_dax.png:
 * M_global(d_ax) good vs bad con M_sim sovrapposto (linea tratteggiata).
 */
void produce_M_vs_dax(const std::vector<MResult>& good,
                      const std::vector<MResult>& bad,
                      const std::filesystem::path& output_path);

/**
 * Produce M_nonlinearity_d{dist}.png:
 * M_local(r) con banda ±1σ.
 */
void produce_M_nonlinearity(const MResult& result,
                             const std::string& config_label,
                             const std::filesystem::path& output_path);

} // namespace riptide::exp3b

#endif // RIPTIDE_EXP3B_ANALYSIS_HPP
