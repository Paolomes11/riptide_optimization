/* Copyright 2026 Giulio Mesini — Apache License 2.0 */

#include "stacking.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

namespace riptide::stack {

static void compute_mean_sigma(const std::vector<double>& vals, double& mean, double& sigma) {
  if (vals.empty()) {
    mean  = 0.0;
    sigma = 0.0;
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

StackedImage sigma_clip_stack(const std::vector<riptide::fits::FitsFrame>& frames,
                              const StackConfig& cfg) {
  if (frames.empty())
    throw std::invalid_argument("sigma_clip_stack: stack vuoto");

  const int W = frames[0].width();
  const int H = frames[0].height();
  const int N = static_cast<int>(frames.size());

  for (int f = 1; f < N; ++f) {
    const auto& fr = frames[static_cast<size_t>(f)];
    if (fr.width() != W || fr.height() != H)
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

  std::vector<double> pixel_vals(static_cast<size_t>(N));
  const size_t total = static_cast<size_t>(W) * static_cast<size_t>(H);

  for (size_t i = 0; i < total; ++i) {
    const int x = static_cast<int>(i % static_cast<size_t>(W));
    const int y = static_cast<int>(i / static_cast<size_t>(W));

    for (int f = 0; f < N; ++f)
      pixel_vals[static_cast<size_t>(f)] = frames[static_cast<size_t>(f)].pixel(x, y);

    std::vector<double> working(pixel_vals);
    for (int iter = 0; iter < cfg.n_iter && working.size() > 1; ++iter) {
      double mu = 0.0;
      double sg = 0.0;
      compute_mean_sigma(working, mu, sg);

      if (sg < 1e-10)
        break;

      std::vector<double> clipped;
      clipped.reserve(working.size());
      for (double v : working) {
        if (std::abs(v - mu) <= cfg.n_sigma * sg)
          clipped.push_back(v);
      }
      if (clipped.empty())
        break;
      working.swap(clipped);
    }

    double mu = 0.0;
    double sg = 0.0;
    compute_mean_sigma(working, mu, sg);

    result.mean[i]  = mu;
    result.sigma[i] = sg;
    result.count[i] = static_cast<int>(working.size());

    if (result.count[i] < cfg.min_frames)
      result.sigma[i] = 0.0;
  }

  const int total_px = W * H;
  int clipped_px     = 0;
  for (int c : result.count)
    if (c < N)
      ++clipped_px;

  std::cout << "[STACK] Sigma-clipping completato: " << clipped_px << "/" << total_px
            << " pixel con almeno un frame rimosso\n";

  return result;
}

StackedImage mean_stack(const std::vector<riptide::fits::FitsFrame>& frames) {
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
    const int x = static_cast<int>(i % static_cast<size_t>(W));
    const int y = static_cast<int>(i / static_cast<size_t>(W));

    std::vector<double> vals(static_cast<size_t>(N));
    for (int f = 0; f < N; ++f)
      vals[static_cast<size_t>(f)] = frames[static_cast<size_t>(f)].pixel(x, y);

    double mu = 0.0;
    double sg = 0.0;
    compute_mean_sigma(vals, mu, sg);
    result.mean[i]  = mu;
    result.sigma[i] = sg;
  }

  return result;
}

StackedImage median_stack(const std::vector<riptide::fits::FitsFrame>& frames) {
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
  result.mean.resize(total, 0.0);
  result.sigma.resize(total, 0.0);
  result.count.resize(total, N);

  std::vector<double> vals(static_cast<size_t>(N));
  for (size_t i = 0; i < total; ++i) {
    const int x = static_cast<int>(i % static_cast<size_t>(W));
    const int y = static_cast<int>(i / static_cast<size_t>(W));

    for (int f = 0; f < N; ++f)
      vals[static_cast<size_t>(f)] = frames[static_cast<size_t>(f)].pixel(x, y);

    std::sort(vals.begin(), vals.end());
    const double median =
        (N % 2 == 0)
            ? (vals[static_cast<size_t>(N / 2 - 1)] + vals[static_cast<size_t>(N / 2)]) / 2.0
            : vals[static_cast<size_t>(N / 2)];

    std::vector<double> absd(static_cast<size_t>(N));
    for (int f = 0; f < N; ++f)
      absd[static_cast<size_t>(f)] = std::abs(vals[static_cast<size_t>(f)] - median);
    std::sort(absd.begin(), absd.end());
    const double mad =
        (N % 2 == 0)
            ? (absd[static_cast<size_t>(N / 2 - 1)] + absd[static_cast<size_t>(N / 2)]) / 2.0
            : absd[static_cast<size_t>(N / 2)];

    result.mean[i]  = median;
    result.sigma[i] = 1.4826 * mad;
  }

  return result;
}

} // namespace riptide::stack
