/* Copyright 2026 Giulio Mesini — Apache License 2.0 */

#ifndef RIPTIDE_EXP2_TRACE_EXTRACTOR_HPP
#define RIPTIDE_EXP2_TRACE_EXTRACTOR_HPP

#include "psf_interpolator.hpp"

#include <vector>

namespace riptide::exp2 {

/// Parametri del fit gaussiano trasversale in un singolo slice della traccia.
struct SliceProfile {
  double t          = 0.0;
  double amplitude  = 0.0;
  double snr        = 0.0;
  double center     = 0.0;
  double center_err = 0.0;
  double sigma      = 0.0;
  double sigma_err  = 0.0;
  double chi2_ndof  = 0.0;
  bool valid        = false;
  bool near_edge    = false;
  bool in_trace     = true;
};

/// Risultato dell'estrazione della traccia da un'immagine stacked.
struct TraceResult {
  double angle_deg     = 0.0;
  double angle_err_deg = 0.0;
  std::vector<SliceProfile> profile;
  double sigma_mean     = 0.0;
  double sigma_mean_err = 0.0;
  double fwhm_mean      = 0.0;
  int n_valid_slices    = 0;
};

/// Risultato del fit ODR lineare sulla posizione del centroide.
using CentroidFitResult = riptide::LineFitResult;

/// Configurazione dell'estrattore.
struct TraceConfig {
  int slice_width           = 5;
  int slice_step            = 3;
  double min_snr            = 5.0;
  int gaussian_range        = 20;
  double sigma_max          = 30.0;
  double center_err_floor   = 0.2;
  double center_err_scale   = 1.0;
  double sigma_err_floor    = 0.2;
  double sigma_err_scale    = 1.0;
  bool enable_trace_trim    = true;
  double trace_trim_frac    = 0.25;
  double trace_trim_min_snr = 3.0;
  int trace_trim_pad_slices = 10;
  int trace_trim_min_slices = 50;
};

/**
 * Stima automatica della direzione della traccia tramite PCA e momento di inerzia
 * sull'immagine dopo sottrazione background.
 *
 * @param img    Immagine differenza (signal - background), layout row-major
 * @param width  Larghezza immagine [pixel]
 * @param height Altezza immagine [pixel]
 * @param cfg    Configurazione
 * @return       Angolo in gradi [-90, 90], con 0 = orizzontale
 */
double estimate_trace_angle(const std::vector<double>& img, int width, int height,
                            const TraceConfig& cfg = {});

/**
 * Variante che restituisce anche l'incertezza sull'angolo stimata tramite sweep della soglia.
 *
 * @param img              Immagine differenza (signal - background), layout row-major
 * @param width            Larghezza immagine [pixel]
 * @param height           Altezza immagine [pixel]
 * @param cfg              Configurazione
 * @param angle_err_deg_out Puntatore (opzionale) dove scrivere l'errore in gradi
 * @return                 Angolo in gradi [-90, 90], con 0 = orizzontale
 */
double estimate_trace_angle(const std::vector<double>& img, int width, int height,
                            const TraceConfig& cfg, double* angle_err_deg_out);

/**
 * Estrae il profilo trasversale della traccia slice per slice.
 *
 * @param img       Immagine differenza (signal - background), layout row-major
 * @param width     Larghezza immagine [pixel]
 * @param height    Altezza immagine [pixel]
 * @param angle_deg Angolo traccia in gradi [-90, 90]
 * @param cfg       Configurazione
 * @return          TraceResult con profilo e statistiche aggregate
 */
TraceResult extract_trace_profile(const std::vector<double>& img, int width, int height,
                                  double angle_deg, const TraceConfig& cfg = {});

/**
 * Esegue il fit ODR lineare sui centroidi delle slice valide.
 *
 * @param trace TraceResult con il profilo estratto
 * @return      Risultato fit ODR riptide::fit_trace()
 */
CentroidFitResult fit_centroid_line(const TraceResult& trace);

} // namespace riptide::exp2

#endif // RIPTIDE_EXP2_TRACE_EXTRACTOR_HPP
