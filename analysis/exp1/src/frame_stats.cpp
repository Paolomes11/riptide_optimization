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

#include "frame_stats.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

namespace exp1 {

// Helper: media e sigma su un vettore di double

static void compute_mean_sigma(const std::vector<double>& vals, double& mean, double& sigma) {
  if (vals.empty()) {
    mean = sigma = 0.0;
    return;
  }
  double sum = 0.0;
  for (double v : vals)
    sum += v;
  mean = sum / static_cast<double>(vals.size());

  if (vals.size() < 2) {
    sigma = 0.0;
    return;
  }
  double var = 0.0;
  for (double v : vals)
    var += (v - mean) * (v - mean);
  sigma = std::sqrt(var / static_cast<double>(vals.size() - 1));
}

// sigma_clip_stack

StackedImage sigma_clip_stack(const std::vector<FitsFrame>& frames, const StackConfig& cfg) {
  if (frames.empty())
    throw std::invalid_argument("sigma_clip_stack: stack vuoto");

  const int W = frames[0].width();
  const int H = frames[0].height();
  const int N = static_cast<int>(frames.size());

  for (int f = 1; f < N; ++f) {
    if (frames[static_cast<size_t>(f)].width() != W || frames[static_cast<size_t>(f)].height() != H)
      throw std::invalid_argument("sigma_clip_stack: dimensioni incoerenti al frame "
                                  + std::to_string(f));
  }

  StackedImage result;
  result.width   = W;
  result.height  = H;
  result.nframes = N;
  result.mean.resize(static_cast<size_t>(W) * static_cast<size_t>(H), 0.0);
  result.sigma.resize(static_cast<size_t>(W) * static_cast<size_t>(H), 0.0);
  result.count.resize(static_cast<size_t>(W) * static_cast<size_t>(H), 0);

  // Scratch buffer per i valori di un pixel su tutti i frame
  std::vector<double> pixel_vals(static_cast<size_t>(N));

  const size_t total = static_cast<size_t>(W) * static_cast<size_t>(H);
  for (size_t i = 0; i < total; ++i) {
    int x = static_cast<int>(i % static_cast<size_t>(W));
    int y = static_cast<int>(i / static_cast<size_t>(W));

    // Raccoglie i valori del pixel i dai N frame
    for (int f = 0; f < N; ++f)
      pixel_vals[static_cast<size_t>(f)] = frames[static_cast<size_t>(f)].pixel(x, y);

    // Sigma-clipping iterativo
    std::vector<double> working(pixel_vals);

    for (int iter = 0; iter < cfg.n_iter && working.size() > 1; ++iter) {
      double mu, sg;
      compute_mean_sigma(working, mu, sg);

      if (sg < 1e-10) // tutti i valori identici: nessun outlier
        break;

      std::vector<double> clipped;
      clipped.reserve(working.size());
      for (double v : working) {
        if (std::abs(v - mu) <= cfg.n_sigma * sg)
          clipped.push_back(v);
      }
      if (clipped.empty())
        break; // non rimuovere tutto
      working.swap(clipped);
    }

    double mu, sg;
    compute_mean_sigma(working, mu, sg);

    result.mean[i]  = mu;
    result.sigma[i] = sg;
    result.count[i] = static_cast<int>(working.size());

    // Se troppo pochi frame sopravvivono, poni sigma = 0 (valore non affidabile)
    if (result.count[i] < cfg.min_frames)
      result.sigma[i] = 0.0;
  }

  int total_px   = W * H;
  int clipped_px = 0;
  for (int c : result.count)
    if (c < N)
      ++clipped_px;

  std::cout << "[STACK] Sigma-clipping completato: " << clipped_px << "/" << total_px
            << " pixel con almeno un frame rimosso\n";

  return result;
}

// mean_stack

StackedImage mean_stack(const std::vector<FitsFrame>& frames) {
  if (frames.empty())
    throw std::invalid_argument("mean_stack: stack vuoto");

  const int W = frames[0].width();
  const int H = frames[0].height();
  const int N = static_cast<int>(frames.size());

  StackedImage result;
  result.width       = W;
  result.height      = H;
  result.nframes     = N;
  const size_t total = static_cast<size_t>(W) * static_cast<size_t>(H);
  result.mean.resize(total, 0.0);
  result.sigma.resize(total, 0.0);
  result.count.resize(total, N);

  for (size_t i = 0; i < total; ++i) {
    int x = static_cast<int>(i % static_cast<size_t>(W));
    int y = static_cast<int>(i / static_cast<size_t>(W));

    std::vector<double> vals(static_cast<size_t>(N));
    for (int f = 0; f < N; ++f)
      vals[static_cast<size_t>(f)] = frames[static_cast<size_t>(f)].pixel(x, y);

    double mu, sg;
    compute_mean_sigma(vals, mu, sg);
    result.mean[i]  = mu;
    result.sigma[i] = sg;
  }

  return result;
}

// median_stack

StackedImage median_stack(const std::vector<FitsFrame>& frames) {
  if (frames.empty())
    throw std::invalid_argument("median_stack: stack vuoto");

  const int W = frames[0].width();
  const int H = frames[0].height();
  const int N = static_cast<int>(frames.size());

  StackedImage result;
  result.width       = W;
  result.height      = H;
  result.nframes     = N;
  const size_t total = static_cast<size_t>(W) * static_cast<size_t>(H);
  result.mean.resize(total, 0.0); // useremo mean per la mediana
  result.sigma.resize(total, 0.0);
  result.count.resize(total, N);

  std::vector<double> vals(static_cast<size_t>(N));
  for (size_t i = 0; i < total; ++i) {
    int x = static_cast<int>(i % static_cast<size_t>(W));
    int y = static_cast<int>(i / static_cast<size_t>(W));

    for (int f = 0; f < N; ++f)
      vals[static_cast<size_t>(f)] = frames[static_cast<size_t>(f)].pixel(x, y);

    std::sort(vals.begin(), vals.end());
    double median =
        (N % 2 == 0)
            ? (vals[static_cast<size_t>(N / 2 - 1)] + vals[static_cast<size_t>(N / 2)]) / 2.0
            : vals[static_cast<size_t>(N / 2)];

    // MAD (Median Absolute Deviation) come stima robusta di sigma
    std::vector<double> absd(static_cast<size_t>(N));
    for (int f = 0; f < N; ++f)
      absd[static_cast<size_t>(f)] = std::abs(vals[static_cast<size_t>(f)] - median);
    std::sort(absd.begin(), absd.end());
    double mad = (N % 2 == 0)
                   ? (absd[static_cast<size_t>(N / 2 - 1)] + absd[static_cast<size_t>(N / 2)]) / 2.0
                   : absd[static_cast<size_t>(N / 2)];
    // Conversione MAD → σ per distribuzione normale: σ ≈ 1.4826 * MAD
    double sigma = 1.4826 * mad;

    result.mean[i]  = median;
    result.sigma[i] = sigma;
  }

  return result;
}

} // namespace exp1