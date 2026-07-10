/* Copyright 2026 Giulio Mesini — Apache License 2.0 */

#ifndef RIPTIDE_EXP3_ANALYSIS_HPP
#define RIPTIDE_EXP3_ANALYSIS_HPP

#include "exp3_config.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace riptide::exp3 {

/// Scala sensor↔display derivata dalla linea di calibrazione.
struct LineCalib {
  double scale_mm_per_sens_px = 0.0; ///< mm fisici per pixel sensore
  double L_px_display         = 0.0; ///< lunghezza nominale linea display [px]
  double L_px_sensor          = 0.0; ///< lunghezza misurata sul sensore [px]
  double d_ax_mm              = 0.0; ///< distanza assiale di questa calibrazione
};

/// Regione di interesse (ROI) per una singola traccia orizzontale.
struct LineROI {
  int    y_center   = 0;
  int    half_width = 20;
  int    line_idx   = 0;  ///< 0=center, 1=+Δ, 2=+2Δ
  double r_mm       = 0.0; ///< offset radiale in mm fisici
};

/// Risultato della misura di Q per una singola traccia (d_ax, r).
struct QResult {
  double d_ax_mm        = 0.0;
  double r_mm           = 0.0;
  int    r_idx          = 0;
  double chi2_ndof      = 0.0;
  int    n_valid_slices = 0;
  bool   warning        = false;
};

/// Configurazione output ROOT/PNG.
struct OutputConfig {
  std::filesystem::path output_dir = "output/exp3";
  bool save_png  = true;
  bool save_root = false;
  bool verbose   = false;
};

/**
 * Calibrazione via linea di riferimento a lunghezza nota.
 *
 * Legge i frame FITS, esegue stack+diff, trova gli estremi della linea
 * orizzontale e calcola la scala mm/px del sensore.
 *
 * @param signal_dir  Cartella frame segnale FITS
 * @param bg_dir      Cartella frame background FITS
 * @param cfg         Configurazione
 * @return            LineCalib con scale_mm_per_sens_px e L_px_sensor
 * @throws            std::runtime_error se il segnale è assente o la linea non rilevata
 */
LineCalib calibrate_line(const std::filesystem::path& signal_dir,
                         const std::filesystem::path& bg_dir,
                         const Exp3Config& cfg);

/**
 * Serializza LineCalib in JSON.
 *
 * @param calib  Risultato di calibrate_line
 * @param path   File di destinazione
 */
void save_line_calib(const LineCalib& calib, const std::filesystem::path& path);

/**
 * Deserializza LineCalib da JSON.
 *
 * @param path  File JSON scritto da save_line_calib
 * @return      LineCalib
 * @throws      std::runtime_error se il file non esiste o è malformato
 */
LineCalib load_line_calib(const std::filesystem::path& path);

/**
 * Segmenta le tracce orizzontali parallele nell'immagine diff usando
 * radial_offsets_px come seed di ricerca (modalità legacy).
 */
std::vector<LineROI> segment_parallel_lines(const std::vector<double>& diff,
                                             int W, int H,
                                             const Exp3Config& cfg,
                                             const LineCalib& calib);

/**
 * Auto-rileva fino a cfg.trace_extraction.max_lines tracce orizzontali
 * nell'immagine diff trovando i picchi del profilo di riga.
 *
 * Non richiede radial_offsets_px né calibrazione. Usata quando
 * radial_offsets_px è vuoto nel config.
 *
 * @param diff  Immagine differenza segnale-background (row-major)
 * @param W     Larghezza immagine sensore [pixel]
 * @param H     Altezza immagine sensore [pixel]
 * @param cfg   Configurazione
 * @return      Vettore di LineROI ordinato per y_center crescente
 */
std::vector<LineROI> detect_lines_auto(const std::vector<double>& diff,
                                        int W, int H,
                                        const Exp3Config& cfg);

/**
 * Pipeline completa per la misura di Q sulle tre tracce parallele.
 *
 * Per ogni LineROI: ritaglia l'immagine diff, stima l'angolo, estrae il
 * profilo centroide, fit ODR → chi2/ndof.
 *
 * @param signal_dir  Cartella frame segnale FITS
 * @param bg_dir      Cartella frame background FITS
 * @param cfg         Configurazione
 * @param calib       Calibrazione per questa distanza assiale
 * @param d_ax_mm     Distanza assiale misurata [mm]
 * @return            Vettore di QResult (uno per traccia)
 */
std::vector<QResult> run_measurement_parallel_lines(
    const std::filesystem::path& signal_dir,
    const std::filesystem::path& bg_dir,
    const Exp3Config& cfg,
    const LineCalib& calib,
    double d_ax_mm);

/**
 * Carica Q_sim da TSV per una coppia (x1, x2).
 *
 * Restituisce NaN se la coppia non è trovata o il file non esiste.
 */
double load_Q_sim(const std::filesystem::path& tsv_path, double x1_mm, double x2_mm);

/**
 * Scrive q_exp_map.tsv: una riga per QResult.
 *
 * Colonne: config  d_ax_mm  r_mm  r_idx  chi2_ndof  n_valid_slices  warning
 *
 * @param config   Label della configurazione lenti ("good" o "bad")
 * @param results  Vettore di QResult
 * @param path     File di destinazione
 */
void write_q_results_tsv(const std::string& config,
                          const std::vector<QResult>& results,
                          const std::filesystem::path& path);

// ---------------------------------------------------------------------------
// Plot (implementate in exp3_plots.cpp)
// ---------------------------------------------------------------------------

/**
 * Produce q_vs_dax_r{r}.png: Q_exp(d_ax) per r_idx fisso, good vs bad.
 */
void produce_q_vs_dax(const std::vector<QResult>& good,
                      const std::vector<QResult>& bad,
                      int r_idx,
                      const std::filesystem::path& output_path);

/**
 * Produce q_comparison.png: Q_exp(d_ax) per good e bad con Q_sim come riferimento.
 *
 * Le misure valide (warning=false, n_valid_slices>0) vengono aggregate per d_ax
 * (media ± range min-max). Q_sim è cercato nel TSV via nearest-neighbor (≤10 mm).
 */
void produce_q_comparison(const std::vector<QResult>& good,
                          const std::vector<QResult>& bad,
                          const QComparisonConfig& cmp,
                          const std::filesystem::path& output_path);

/**
 * Converte ricorsivamente ogni file .fit/.fits/.fts sotto data_dir in un
 * PNG di ispezione visiva (scala colore percentile, palette kViridis),
 * mantenendo la struttura di sottocartelle sotto output_dir/fits_preview/.
 *
 * @param data_dir    Cartella radice con i frame FITS grezzi (scansione ricorsiva)
 * @param output_dir  Cartella di output exp3 (i PNG finiscono in output_dir/fits_preview/...)
 * @return            Numero di PNG generati
 */
int produce_fits_previews(const std::filesystem::path& data_dir,
                          const std::filesystem::path& output_dir);

} // namespace riptide::exp3

#endif // RIPTIDE_EXP3_ANALYSIS_HPP
