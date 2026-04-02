/*
 * Copyright 2026 Giulio Mesini
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 */

#ifndef EXP1_FRAME_STATS_HPP
#define EXP1_FRAME_STATS_HPP

/**
 * frame_stats — Stacking statistico di uno stack di frame FITS
 *
 * Per ogni pixel (x, y), raccoglie i valori dai 10 frame e calcola:
 *   1. Media e deviazione standard naive sull'intero stack
 *   2. Sigma-clipping iterativo: esclude i valori oltre n_sigma * σ
 *   3. Ricalcola media e σ sui pixel sopravvissuti al clipping
 *
 * Il risultato è una StackedImage con:
 *   - mean[i]:   media clippata per pixel i
 *   - sigma[i]:  deviazione standard clippata (= incertezza sul segnale)
 *   - count[i]:  numero di frame contribuenti dopo il clipping
 *
 * Il sigma-clipping è lo standard in astronomia/imaging scientifico per
 * eliminare raggi cosmici, hot pixel e artefatti impulsivi.
 */

#include "fit_reader.hpp"

#include <cstdint>
#include <vector>

namespace exp1 {

// Immagine stacked

struct StackedImage {
  int width   = 0;
  int height  = 0;
  int nframes = 0; // numero di frame nello stack

  // Layout row-major: indice i = y * width + x
  std::vector<double> mean;  // media clippata [ADU, 0..65535]
  std::vector<double> sigma; // deviazione standard clippata [ADU]
  std::vector<int> count;    // frame contribuenti dopo sigma-clipping

  // Accesso sicuro
  double pixel_mean(int x, int y) const {
    if (x < 0 || x >= width || y < 0 || y >= height)
      return 0.0;
    return mean[static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)];
  }
  double pixel_sigma(int x, int y) const {
    if (x < 0 || x >= width || y < 0 || y >= height)
      return 0.0;
    return sigma[static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)];
  }
  size_t npixels() const {
    return static_cast<size_t>(width) * static_cast<size_t>(height);
  }
};

// Parametri stacking

struct StackConfig {
  double n_sigma = 3.0; // soglia sigma-clipping
  int n_iter     = 3;   // numero iterazioni sigma-clipping
  int min_frames = 3;   // frame minimi per pixel valido (altrimenti sigma = 0)
};

// API pubblica

/**
 * Esegue lo stacking statistico di una sequenza di frame.
 *
 * Tutti i frame devono avere le stesse dimensioni (width × height).
 * L'algoritmo per pixel:
 *   1. Raccoglie i valori dei N frame
 *   2. Calcola μ e σ
 *   3. Rimuove i valori con |v - μ| > n_sigma * σ
 *   4. Ripete n_iter volte
 *   5. Calcola μ finale e σ finale sui valori sopravvissuti
 *
 * @param frames  Stack di frame (deve averne almeno 2)
 * @param cfg     Parametri sigma-clipping
 * @return        StackedImage con mean, sigma, count per pixel
 * @throws        std::invalid_argument se frames è vuoto o dimensioni incoerenti
 */
StackedImage sigma_clip_stack(const std::vector<FitsFrame>& frames, const StackConfig& cfg = {});

/**
 * Media semplice (senza sigma-clipping) — usata per confronto/debug.
 */
StackedImage mean_stack(const std::vector<FitsFrame>& frames);

/**
 * Mediana per-pixel — robusta agli outlier, ma più lenta.
 * Utile come sanity check rispetto al sigma-clipping.
 */
StackedImage median_stack(const std::vector<FitsFrame>& frames);

} // namespace exp1

#endif // EXP1_FRAME_STATS_HPP