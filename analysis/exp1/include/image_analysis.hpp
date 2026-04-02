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

#ifndef EXP1_IMAGE_ANALYSIS_HPP
#define EXP1_IMAGE_ANALYSIS_HPP

/**
 * image_analysis — Sottrazione fondo, integrazione e plotting ROOT
 *
 * Pipeline:
 *
 *   1. SOTTRAZIONE FONDO (statistica)
 *      diff = stacked_signal - stacked_background
 *      σ_diff = sqrt(σ_signal² + σ_background²)   [propagazione quadratica]
 *
 *   2. INTEGRAZIONE
 *      integral = Σ_{pixel ∈ ROI} diff[pixel]
 *      σ_integral = sqrt(Σ_{pixel ∈ ROI} σ_diff[pixel]²)
 *
 *   3. CONFRONTO good vs bad
 *      ΔI = integral_good - integral_bad
 *      σ_ΔI = sqrt(σ_good² + σ_bad²)
 *      significance = ΔI / σ_ΔI  [in unità σ]
 *
 *   4. OUTPUT ROOT
 *      - TH2D mappa 2D del segnale netto (colori)
 *      - TH2D mappa σ (incertezza)
 *      - TH1D profilo integrato lungo X e lungo Y
 *      - Pannello riassuntivo con ΔI e significatività
 */

#include "frame_stats.hpp"

#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>

#include <filesystem>
#include <optional>
#include <string>

namespace exp1 {

// ROI (Region of Interest)

struct ROI {
  int x0 = 0, y0 = 0;   // angolo in alto a sinistra (inclusivo)
  int x1 = -1, y1 = -1; // angolo in basso a destra (inclusivo, -1 = tutto)

  // Ritorna true se il pixel è nella ROI (dopo aver risolto i -1)
  bool contains(int x, int y, int width, int height) const;

  // Risolve i valori -1 con le dimensioni reali dell'immagine
  ROI resolve(int width, int height) const;
};

// Immagine differenza (signal - background)

struct DiffImage {
  int width = 0, height = 0;

  std::vector<double> diff;  // segnale netto [ADU]
  std::vector<double> sigma; // incertezza propagata [ADU]

  double pixel_diff(int x, int y) const;
  double pixel_sigma(int x, int y) const;
  size_t npixels() const {
    return static_cast<size_t>(width) * static_cast<size_t>(height);
  }
};

// Risultato integrazione

struct IntegralResult {
  double integral;       // Σ diff [ADU·pixel]
  double sigma_integral; // sqrt(Σ σ²) [ADU·pixel]
  int n_pixels;          // pixel della ROI contribuenti
  ROI roi_used;          // ROI effettivamente usata
};

// Risultato confronto good vs bad

struct ComparisonResult {
  IntegralResult good;
  IntegralResult bad;

  double delta;       // integral_good - integral_bad [ADU·pixel]
  double sigma_delta; // sqrt(σ_good² + σ_bad²)
  double ratio;       // integral_good / integral_bad
  double sigma_ratio; // incertezza propagata sul rapporto
};

// Configurazione analisi

struct AnalysisConfig {
  ROI roi; // ROI per l'integrazione (default: intero frame)
  bool save_png                    = true;
  bool save_root                   = true;
  std::filesystem::path output_dir = "output/exp1";

  // Scaling percentili per la visualizzazione (Z-axis)
  double z_min_percentile = 0.005; // 0.5% inferiore (ignora pixel morti/neri)
  double z_max_percentile = 0.995; // 99.5% superiore (ignora hot pixel)
};

// API pubblica

/**
 * Sottrae il fondo da un'immagine stacked con propagazione dell'incertezza.
 *
 * @param signal      Immagine segnale (good o bad) già stacked
 * @param background  Immagine fondo già stacked
 * @return            DiffImage con diff = signal.mean - background.mean
 *                    e sigma = sqrt(signal.sigma² + background.sigma²)
 * @throws            std::invalid_argument se le dimensioni non coincidono
 */
DiffImage subtract_background(const StackedImage& signal, const StackedImage& background);

/**
 * Integra il segnale netto nella ROI specificata.
 */
IntegralResult integrate(const DiffImage& diff, const ROI& roi = {});

/**
 * Calcola ΔI = integral_good - integral_bad con significatività statistica.
 */
ComparisonResult compare(const IntegralResult& good, const IntegralResult& bad);

/**
 * Pipeline completa: produce tutte le mappe ROOT e il pannello di confronto.
 *
 * @param good_diff   DiffImage (good - background)
 * @param bad_diff    DiffImage (bad  - background)
 * @param good_stack  StackedImage good (per istogrammi diagnostici)
 * @param bad_stack   StackedImage bad
 * @param bg_stack    StackedImage background
 * @param comparison  Risultato del confronto
 * @param cfg         Configurazione output
 */
void produce_output(const DiffImage& good_diff, const DiffImage& bad_diff,
                    const StackedImage& good_stack, const StackedImage& bad_stack,
                    const StackedImage& bg_stack, const ComparisonResult& comparison,
                    const AnalysisConfig& cfg);

// Helpers per ROOT

/**
 * Converte uno StackedImage in TH2D ROOT (mean o sigma).
 * Il chiamante è responsabile della vita dell'oggetto.
 */
TH2D* stacked_to_th2d(const StackedImage& img, const std::string& name, const std::string& title,
                      bool use_sigma = false);

/**
 * Converte una DiffImage in TH2D ROOT.
 */
TH2D* diff_to_th2d(const DiffImage& img, const std::string& name, const std::string& title,
                   bool use_sigma = false);

/**
 * Proiezione orizzontale (integrazione su Y) di una DiffImage → TH1D.
 * Utile per visualizzare il profilo del segnale lungo X.
 */
TH1D* project_x(const DiffImage& img, const std::string& name, const std::string& title,
                const ROI& roi = {});

/**
 * Proiezione verticale (integrazione su X) di una DiffImage → TH1D.
 */
TH1D* project_y(const DiffImage& img, const std::string& name, const std::string& title,
                const ROI& roi = {});

/**
 * Applica lo stile ROOT publication-quality standard del progetto RIPTIDE.
 */
void apply_riptide_style();

} // namespace exp1

#endif // EXP1_IMAGE_ANALYSIS_HPP