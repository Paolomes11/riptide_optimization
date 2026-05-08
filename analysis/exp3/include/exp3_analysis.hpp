/* Copyright 2026 Giulio Mesini — Apache License 2.0 */

#ifndef RIPTIDE_EXP3_ANALYSIS_HPP
#define RIPTIDE_EXP3_ANALYSIS_HPP

#include "exp3_config.hpp"
#include "homography.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace riptide::exp3 {

/// Risultato del fit ODR per una singola acquisizione (d_ax, θ).
struct SingleTraceResult {
  double d_ax_nominal_mm  = 0.0;
  double d_ax_measured_mm = 0.0;
  double theta_nominal_deg  = 0.0;
  double theta_measured_deg = 0.0;
  double chi2_ndof          = 0.0;
  int    n_valid_slices     = 0;
  bool   valid              = false;
  std::string warning;
};

/// Risultato per una distanza assiale (media su orientazioni).
struct AxialResult {
  double d_ax_nominal_mm  = 0.0;
  double d_ax_measured_mm = 0.0;
  double Q_exp            = 0.0;
  double Q_exp_sigma      = 0.0;
  int    n_valid_orientations = 0;
  std::vector<SingleTraceResult> traces;
};

/// Risultato globale per una configurazione lenti.
struct ConfigResult {
  std::string config_label;
  LensConfig  lens;
  std::vector<AxialResult> axial;
  double Q_exp_global = 0.0;
  double Q_sim        = 0.0;
  double R            = 0.0;
};

/// Configurazione output ROOT/PNG.
struct OutputConfig {
  std::filesystem::path output_dir  = "output/exp3";
  bool save_png                     = true;
  bool save_root                    = false;
  double z_min_percentile           = 0.005;
  double z_max_percentile           = 0.995;
  bool verbose                      = false;
};

/**
 * Analisi di una singola acquisizione (d_ax, θ):
 * stack → sottrazione bg → stima angolo → estrazione profilo → fit ODR.
 *
 * Riusa riptide::exp2::estimate_trace_angle, extract_trace_profile, fit_centroid_line.
 *
 * @param signal_dir        Cartella con i frame .fit del segnale
 * @param bg_dir            Cartella con i frame di background
 * @param theta_nominal_deg Angolo nominale della traccia [gradi]
 * @param d_ax_nominal_mm   Distanza assiale nominale [mm]
 * @param d_ax_measured_mm  Distanza assiale misurata [mm]
 * @param H                 Omografia calibrata per questa distanza
 * @param cfg               Configurazione
 * @return                  SingleTraceResult con chi2_ndof e flag di validità
 */
SingleTraceResult analyze_single_trace(const std::filesystem::path& signal_dir,
                                       const std::filesystem::path& bg_dir,
                                       double theta_nominal_deg,
                                       double d_ax_nominal_mm,
                                       double d_ax_measured_mm,
                                       const Homography& H,
                                       const Exp3Config& cfg);

/**
 * Calibrazione geometrica per una distanza assiale:
 * legge i frame della griglia, rileva i dot, calcola H.
 *
 * @param signal_dir  Cartella con i frame FITS della griglia (segnale)
 * @param bg_dir      Cartella con i frame di background
 * @param cfg         Configurazione
 * @return            Omografia calibrata
 */
Homography calibrate_distance(const std::filesystem::path& signal_dir,
                              const std::filesystem::path& bg_dir,
                              const Exp3Config& cfg);

/**
 * Rilevamento dot della griglia di calibrazione.
 *
 * Per ogni nodo della griglia teorica (passo grid_step_px), esegue un fit
 * gaussiano 2D in una finestra centrata sulla posizione attesa e, se il fit
 * converge con SNR sufficiente, aggiunge il punto di calibrazione.
 *
 * @param diff_image         Immagine differenza segnale-background (row-major)
 * @param width              Larghezza immagine [pixel]
 * @param height             Altezza immagine [pixel]
 * @param grid_step_px       Passo griglia nel display [pixel display]
 * @param display_width_px   Larghezza display [pixel]
 * @param display_height_px  Altezza display [pixel]
 * @return                   Lista di CalibPoint (attesi→trovati)
 */
std::vector<CalibPoint> detect_calibration_dots(const std::vector<double>& diff_image,
                                                int width, int height,
                                                int grid_step_px,
                                                int display_width_px,
                                                int display_height_px);

/**
 * Pipeline completa per una configurazione lenti.
 *
 * Itera su tutte le distanze assiali e orientazioni, producendo output ROOT/PNG
 * opzionali inline. Assume che le omografie siano già state calcolate e salvate
 * in out_cfg.output_dir/calib/homography_d{dist}mm.json.
 *
 * Struttura dati attesa:
 *   data_root/{label}/d{dist}/theta{theta}/  ← frame segnale
 *   data_root/{label}/background/            ← frame background
 *
 * @param data_root   Root cartella dati (es. data/exp3/)
 * @param calib_root  Root cartella calibrazione (es. data/exp3/calib/ o output/exp3/calib/)
 * @param lens        Configurazione lenti da analizzare
 * @param cfg         Configurazione esperimento
 * @param out_cfg     Configurazione output
 * @return            ConfigResult con Q_exp_global e Q_sim
 */
ConfigResult analyze_config(const std::filesystem::path& data_root,
                            const std::filesystem::path& calib_root,
                            const LensConfig& lens,
                            const Exp3Config& cfg,
                            const OutputConfig& out_cfg = {});

/**
 * Carica Q_sim da un file TSV prodotto da q_map per una coppia (x1, x2).
 *
 * Usa la prima riga con |x1_file - x1| < 0.5 mm e |x2_file - x2| < 0.5 mm.
 * Restituisce NaN se nessuna riga corrisponde.
 *
 * @param tsv_path  Path al file TSV
 * @param x1_mm    Posizione prima lente [mm]
 * @param x2_mm    Posizione seconda lente [mm]
 * @return          Q_sim (o NaN se non trovato)
 */
double load_Q_sim(const std::filesystem::path& tsv_path, double x1_mm, double x2_mm);

/**
 * Produce calib_report.png: scatter dei punti di calibrazione e RMS residui per distanza.
 *
 * @param homographies    Vettore di omografie (una per distanza assiale)
 * @param axial_dists     Distanze assiali corrispondenti [mm]
 * @param pts_per_dist    Punti di calibrazione per ogni distanza
 * @param output_path     Path del file PNG di output
 */
void produce_calibration_report(const std::vector<Homography>& homographies,
                                const std::vector<double>& axial_dists,
                                const std::vector<std::vector<CalibPoint>>& pts_per_dist,
                                const std::filesystem::path& output_path);

/**
 * Produce summary.png: confronto Q_exp_global e Q_sim per tutte le configurazioni.
 *
 * @param results      Vettore di risultati (una per configurazione)
 * @param output_path  Path del file PNG di output
 */
void produce_results_summary(const std::vector<ConfigResult>& results,
                             const std::filesystem::path& output_path);

/**
 * Scrive results.tsv: tabella completa per ogni singola traccia.
 *
 * @param results  Vettore di risultati
 * @param path     Path del file TSV di output
 */
void write_results_tsv(const std::vector<ConfigResult>& results,
                       const std::filesystem::path& path);

/**
 * Scrive summary.tsv: Q_exp_global, Q_sim, R per ogni configurazione.
 *
 * @param results  Vettore di risultati
 * @param path     Path del file TSV di output
 */
void write_summary_tsv(const std::vector<ConfigResult>& results,
                       const std::filesystem::path& path);

} // namespace riptide::exp3

#endif // RIPTIDE_EXP3_ANALYSIS_HPP
