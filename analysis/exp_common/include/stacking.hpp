/* Copyright 2026 Giulio Mesini — Apache License 2.0 */

#ifndef RIPTIDE_EXP_COMMON_STACKING_HPP
#define RIPTIDE_EXP_COMMON_STACKING_HPP

/**
 * stacking — Stacking statistico di uno stack di frame FITS
 *
 * Per ogni pixel (x, y), raccoglie i valori dai N frame e calcola:
 *   - Media e deviazione standard naive
 *   - Sigma-clipping iterativo
 *   - Mediana + MAD
 */

#include "fits_io.hpp"

#include <cstddef>
#include <vector>

namespace riptide::stack {

/// Immagine stacked: mean e sigma per pixel, più count di frame contribuenti.
struct StackedImage {
  int width   = 0;
  int height  = 0;
  int nframes = 0;
  std::vector<double> mean;
  std::vector<double> sigma;
  std::vector<int> count;

  /// Accesso sicuro alla media.
  double pixel_mean(int x, int y) const {
    if (x < 0 || x >= width || y < 0 || y >= height)
      return 0.0;
    return mean[static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)];
  }

  /// Accesso sicuro alla sigma.
  double pixel_sigma(int x, int y) const {
    if (x < 0 || x >= width || y < 0 || y >= height)
      return 0.0;
    return sigma[static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)];
  }

  /// Numero totale di pixel.
  size_t npixels() const {
    return static_cast<size_t>(width) * static_cast<size_t>(height);
  }
};

/// Metodo stacking selezionabile da CLI.
enum class StackMethod { SigmaClip, Mean, Median };

/// Parametri stacking.
struct StackConfig {
  StackMethod method = StackMethod::SigmaClip;
  size_t max_frames  = 0;
  double n_sigma     = 3.0;
  int n_iter         = 3;
  int min_frames     = 3;
};

/**
 * Esegue lo stacking statistico sigma-clipping di una sequenza di frame.
 *
 * Tutti i frame devono avere le stesse dimensioni (width × height).
 *
 * @param frames  Stack di frame (deve averne almeno 1)
 * @param cfg     Parametri sigma-clipping
 * @return        StackedImage con mean, sigma, count per pixel
 * @throws        std::invalid_argument se frames è vuoto o dimensioni incoerenti
 */
StackedImage sigma_clip_stack(const std::vector<riptide::fits::FitsFrame>& frames,
                              const StackConfig& cfg = {});

/**
 * Media semplice (senza sigma-clipping).
 *
 * @param frames  Stack di frame (deve averne almeno 1)
 * @return        StackedImage con mean e sigma naive
 */
StackedImage mean_stack(const std::vector<riptide::fits::FitsFrame>& frames);

/**
 * Mediana per-pixel + MAD come stima robusta di sigma.
 *
 * @param frames  Stack di frame (deve averne almeno 1)
 * @return        StackedImage con mean=mediana e sigma≈1.4826*MAD
 */
StackedImage median_stack(const std::vector<riptide::fits::FitsFrame>& frames);

} // namespace riptide::stack

#endif // RIPTIDE_EXP_COMMON_STACKING_HPP
