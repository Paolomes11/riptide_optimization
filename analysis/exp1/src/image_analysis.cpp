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

#include "image_analysis.hpp"

#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace exp1 {

// ROI

ROI ROI::resolve(int width, int height) const {
  ROI r = *this;
  if (r.x1 < 0)
    r.x1 = width - 1;
  if (r.y1 < 0)
    r.y1 = height - 1;
  r.x0 = std::max(0, r.x0);
  r.y0 = std::max(0, r.y0);
  r.x1 = std::min(width - 1, r.x1);
  r.y1 = std::min(height - 1, r.y1);
  return r;
}

bool ROI::contains(int x, int y, int width, int height) const {
  ROI r = resolve(width, height);
  return x >= r.x0 && x <= r.x1 && y >= r.y0 && y <= r.y1;
}

// DiffImage accessors

double DiffImage::pixel_diff(int x, int y) const {
  if (x < 0 || x >= width || y < 0 || y >= height)
    return 0.0;
  return diff[static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)];
}

double DiffImage::pixel_sigma(int x, int y) const {
  if (x < 0 || x >= width || y < 0 || y >= height)
    return 0.0;
  return sigma[static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)];
}

// subtract_background

DiffImage subtract_background(const StackedImage& signal, const StackedImage& background) {
  if (signal.width != background.width || signal.height != background.height)
    throw std::invalid_argument("subtract_background: dimensioni diverse tra segnale e fondo ("
                                + std::to_string(signal.width) + "×" + std::to_string(signal.height)
                                + " vs " + std::to_string(background.width) + "×"
                                + std::to_string(background.height) + ")");

  DiffImage out;
  out.width      = signal.width;
  out.height     = signal.height;
  const size_t N = out.npixels();
  out.diff.resize(N);
  out.sigma.resize(N);

  for (size_t i = 0; i < N; ++i) {
    // diff = segnale - fondo
    out.diff[i] = signal.mean[i] - background.mean[i];

    // Propagazione incertezza quadratica: σ_diff = sqrt(σ_s² + σ_bg²)
    double ss    = signal.sigma[i];
    double sb    = background.sigma[i];
    out.sigma[i] = std::sqrt(ss * ss + sb * sb);
  }

  return out;
}

// analyze_frame_by_frame

FrameByFrameResult analyze_frame_by_frame(const std::vector<FitsFrame>& good_frames,
                                          const std::vector<FitsFrame>& bad_frames,
                                          const std::vector<FitsFrame>& bg_frames, const ROI& roi,
                                          double clip_sigma) {
  size_t N = good_frames.size();
  if (bad_frames.size() != N || bg_frames.size() != N) {
    throw std::invalid_argument("analyze_frame_by_frame: numero di frame incoerente ("
                                + std::to_string(N) + " vs " + std::to_string(bad_frames.size())
                                + " vs " + std::to_string(bg_frames.size()) + ")");
  }

  std::vector<double> ratios;
  ratios.reserve(N);

  for (size_t k = 0; k < N; ++k) {
    ROI r_good = roi.resolve(good_frames[k].width(), good_frames[k].height());
    ROI r_bad  = roi.resolve(bad_frames[k].width(), bad_frames[k].height());
    ROI r_bg   = roi.resolve(bg_frames[k].width(), bg_frames[k].height());

    auto integrate_frame = [](const FitsFrame& f, const ROI& r) {
      double sum = 0;
      for (int y = r.y0; y <= r.y1; ++y) {
        for (int x = r.x0; x <= r.x1; ++x) {
          sum += f.pixel(x, y);
        }
      }
      return sum;
    };

    double sum_good = integrate_frame(good_frames[k], r_good);
    double sum_bad  = integrate_frame(bad_frames[k], r_bad);
    double sum_bg   = integrate_frame(bg_frames[k], r_bg);

    double Ig_k = sum_good - sum_bg;
    double Ib_k = sum_bad - sum_bg;

    if (Ib_k != 0) {
      ratios.push_back(Ig_k / Ib_k);
    }
  }

  if (ratios.empty()) {
    throw std::runtime_error("analyze_frame_by_frame: nessun rapporto calcolato (I_bad = 0?)");
  }

  // Outlier rejection (clipping a n_sigma)
  if (clip_sigma > 0 && ratios.size() > 2) {
    double sum_r = std::accumulate(ratios.begin(), ratios.end(), 0.0);
    double mean  = sum_r / ratios.size();
    double sq_sum =
        std::inner_product(ratios.begin(), ratios.end(), ratios.begin(), 0.0, std::plus<>(),
                           [mean](double a, double b) { return (a - mean) * (b - mean); });
    double stdev = std::sqrt(sq_sum / (ratios.size() - 1));

    std::vector<double> filtered;
    for (double r : ratios) {
      if (std::abs(r - mean) <= clip_sigma * stdev) {
        filtered.push_back(r);
      }
    }
    ratios = std::move(filtered);
  }

  size_t n_used = ratios.size();
  if (n_used < 2) {
    throw std::runtime_error("analyze_frame_by_frame: troppi pochi frame validi ("
                             + std::to_string(n_used) + ")");
  }

  double mean_R = std::accumulate(ratios.begin(), ratios.end(), 0.0) / n_used;
  double sq_sum =
      std::inner_product(ratios.begin(), ratios.end(), ratios.begin(), 0.0, std::plus<>(),
                         [mean_R](double a, double b) { return (a - mean_R) * (b - mean_R); });

  double var_R        = sq_sum / (n_used - 1);
  double sigma_R      = std::sqrt(var_R);
  double sigma_mean_R = sigma_R / std::sqrt(n_used);
  double significance = (mean_R - 1.0) / sigma_mean_R;

  FrameByFrameResult res;
  res.mean_ratio       = mean_R;
  res.sigma_ratio_std  = sigma_R;
  res.sigma_ratio_mean = sigma_mean_R;
  res.significance     = significance;
  res.n_frames_used    = static_cast<int>(n_used);
  res.ratios           = ratios;

  return res;
}

// integrate

IntegralResult integrate(const DiffImage& diff, const ROI& roi) {
  IntegralResult res{};
  res.roi_used = roi.resolve(diff.width, diff.height);

  double sum_val = 0.0;
  double sum_var = 0.0;
  int n_pixels   = 0;

  for (int y = res.roi_used.y0; y <= res.roi_used.y1; ++y) {
    for (int x = res.roi_used.x0; x <= res.roi_used.x1; ++x) {
      double v = diff.pixel_diff(x, y);
      double s = diff.pixel_sigma(x, y);
      sum_val += v;
      sum_var += s * s;
      ++n_pixels;
    }
  }

  res.integral       = sum_val;
  res.sigma_integral = std::sqrt(sum_var);
  res.n_pixels       = n_pixels;

  return res;
}

// compare

ComparisonResult compare(const IntegralResult& good, const IntegralResult& bad) {
  ComparisonResult res;
  res.good        = good;
  res.bad         = bad;
  res.delta       = good.integral - bad.integral;
  res.sigma_delta = std::sqrt(good.sigma_integral * good.sigma_integral
                              + bad.sigma_integral * bad.sigma_integral);

  // Calcolo rapporto e propagazione incertezza
  if (std::abs(bad.integral) > 1e-9) {
    res.ratio        = good.integral / bad.integral;
    double rel_good  = (std::abs(good.integral) > 1e-9) ? good.sigma_integral / good.integral : 0.0;
    double rel_bad   = bad.sigma_integral / bad.integral;
    res.sigma_ratio  = std::abs(res.ratio) * std::sqrt(rel_good * rel_good + rel_bad * rel_bad);
    res.significance = (res.sigma_ratio > 1e-12) ? (res.ratio - 1.0) / res.sigma_ratio : 0.0;
  } else {
    res.ratio        = 0.0;
    res.sigma_ratio  = 0.0;
    res.significance = 0.0;
  }

  return res;
}

// apply_riptide_style

void apply_riptide_style() {
  gStyle->Reset();
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleFont(42, "");
  gStyle->SetStatFont(42);
  gStyle->SetTextSize(0.040);
  gStyle->SetLabelSize(0.038, "XYZ");
  gStyle->SetTitleSize(0.044, "XYZ");
  gStyle->SetTitleSize(0.046, "");
  gStyle->SetTitleOffset(1.55, "Y");
  gStyle->SetTitleOffset(1.20, "X");
  gStyle->SetTickLength(0.018, "X");
  gStyle->SetTickLength(0.018, "Y");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetNdivisions(505, "Y");
  gStyle->SetNumberContours(255);
}

// Helpers TH2D

TH2D* stacked_to_th2d(const StackedImage& img, const std::string& name, const std::string& title,
                      bool use_sigma) {
  auto* h = new TH2D(name.c_str(), title.c_str(), img.width, 0.0, static_cast<double>(img.width),
                     img.height, 0.0, static_cast<double>(img.height));

  for (int y = 0; y < img.height; ++y) {
    for (int x = 0; x < img.width; ++x) {
      size_t i   = static_cast<size_t>(y) * static_cast<size_t>(img.width) + static_cast<size_t>(x);
      double val = use_sigma ? img.sigma[i] : img.mean[i];
      h->SetBinContent(x + 1, y + 1, val);
    }
  }

  h->GetXaxis()->SetTitle("Pixel X");
  h->GetYaxis()->SetTitle("Pixel Y");
  h->GetXaxis()->CenterTitle(true);
  h->GetYaxis()->CenterTitle(true);
  h->GetXaxis()->SetTitleOffset(1.20);
  h->GetYaxis()->SetTitleOffset(1.55);
  return h;
}

TH2D* diff_to_th2d(const DiffImage& img, const std::string& name, const std::string& title,
                   bool use_sigma) {
  auto* h = new TH2D(name.c_str(), title.c_str(), img.width, 0.0, static_cast<double>(img.width),
                     img.height, 0.0, static_cast<double>(img.height));

  for (int y = 0; y < img.height; ++y) {
    for (int x = 0; x < img.width; ++x) {
      double val = use_sigma ? img.pixel_sigma(x, y) : img.pixel_diff(x, y);
      h->SetBinContent(x + 1, y + 1, val);
    }
  }

  h->GetXaxis()->SetTitle("Pixel X");
  h->GetYaxis()->SetTitle("Pixel Y");
  h->GetXaxis()->CenterTitle(true);
  h->GetYaxis()->CenterTitle(true);
  h->GetXaxis()->SetTitleOffset(1.20);
  h->GetYaxis()->SetTitleOffset(1.55);
  return h;
}

TH1D* project_x(const DiffImage& img, const std::string& name, const std::string& title,
                const ROI& roi) {
  ROI r   = roi.resolve(img.width, img.height);
  auto* h = new TH1D(name.c_str(), title.c_str(), img.width, 0.0, static_cast<double>(img.width));

  for (int x = r.x0; x <= r.x1; ++x) {
    double sum = 0.0;
    for (int y = r.y0; y <= r.y1; ++y)
      sum += img.pixel_diff(x, y);
    h->SetBinContent(x + 1, sum);
  }

  h->GetXaxis()->SetTitle("Pixel X");
  h->GetYaxis()->SetTitle("Intensit#grave{a} integrata [ADU]");
  h->GetXaxis()->CenterTitle(true);
  h->GetYaxis()->CenterTitle(true);
  h->GetXaxis()->SetTitleOffset(1.20);
  h->GetYaxis()->SetTitleOffset(1.55);
  return h;
}

TH1D* project_y(const DiffImage& img, const std::string& name, const std::string& title,
                const ROI& roi) {
  ROI r   = roi.resolve(img.width, img.height);
  auto* h = new TH1D(name.c_str(), title.c_str(), img.height, 0.0, static_cast<double>(img.height));

  for (int y = r.y0; y <= r.y1; ++y) {
    double sum = 0.0;
    for (int x = r.x0; x <= r.x1; ++x)
      sum += img.pixel_diff(x, y);
    h->SetBinContent(y + 1, sum);
  }

  h->GetXaxis()->SetTitle("Pixel Y");
  h->GetYaxis()->SetTitle("Intensit#grave{a} integrata [ADU]");
  h->GetXaxis()->CenterTitle(true);
  h->GetYaxis()->CenterTitle(true);
  h->GetXaxis()->SetTitleOffset(1.20);
  h->GetYaxis()->SetTitleOffset(1.55);
  return h;
}

// Helper: format double

static std::string fmtd(double v, int prec = 2) {
  std::ostringstream o;
  o << std::fixed << std::setprecision(prec) << v;
  return o.str();
}

static std::string fmt_sci(double v, int prec = 3) {
  std::ostringstream o;
  o << std::scientific << std::setprecision(prec) << v;
  return o.str();
}

/**
 * Helper: Calcola un percentile di un istogramma TH2D per lo scaling della Z.
 * Utile per ignorare hot pixel o outlier che appiattiscono la scala colori.
 */
static double get_th2d_percentile(TH2D* h, double percentile) {
  if (!h || percentile <= 0.0 || percentile >= 1.0)
    return 0.0;
  std::vector<double> vals;
  vals.reserve(static_cast<size_t>(h->GetNbinsX() * h->GetNbinsY()));
  for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
    for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
      double v = h->GetBinContent(ix, iy);
      if (v > 1e-9) // ignora zeri/pixel neri
        vals.push_back(v);
    }
  }
  if (vals.empty())
    return 0.0;
  size_t idx = static_cast<size_t>(static_cast<double>(vals.size()) * percentile);
  std::nth_element(vals.begin(), vals.begin() + static_cast<long>(idx), vals.end());
  return vals[idx];
}

/**
 * Helper: Rebinna un TH2D se è troppo grande per la visualizzazione in PNG.
 * Questo riduce drasticamente l'occupazione di memoria durante il rendering del canvas.
 */
static TH2D* rebin_for_display(TH2D* h) {
  if (!h)
    return nullptr;
  int nx = h->GetNbinsX();
  int ny = h->GetNbinsY();
  int rx = 1;
  int ry = 1;

  if (nx > 2000)
    rx = nx / 1024;
  if (ny > 2000)
    ry = ny / 1024;

  if (rx <= 1 && ry <= 1)
    return static_cast<TH2D*>(h->Clone());

  TH2D* h_rebinned =
      static_cast<TH2D*>(h->Rebin2D(rx, ry, (std::string(h->GetName()) + "_disp").c_str()));
  return h_rebinned;
}

// produce_output

void produce_output(const DiffImage& good_diff, const DiffImage& bad_diff,
                    const StackedImage& good_stack, const StackedImage& bad_stack,
                    const StackedImage& bg_stack, const ComparisonResult& comparison,
                    const std::optional<FrameByFrameResult>& fbf, const AnalysisConfig& cfg) {
  std::filesystem::create_directories(cfg.output_dir);

  apply_riptide_style();
  gStyle->SetPalette(kViridis);

  // Apertura TFile per salvataggio istogrammi a piena risoluzione
  // Salviamo in un unico file invece di usare Canvas::Print(".root") che fallisce sopra 1GB
  std::unique_ptr<TFile> root_file;
  if (cfg.save_root) {
    root_file = std::make_unique<TFile>((cfg.output_dir / "exp1_analysis.root").string().c_str(),
                                        "RECREATE");
    if (root_file->IsZombie())
      std::cerr << "[WARNING] Impossibile creare il file ROOT di output\n";
  }

  // Pagina 1: mappe 2D delle immagini stacked

  {
    TCanvas c("c_stacked", "Stacked Images", 1800, 600);
    c.Divide(3, 1, 0.005, 0.005);

    auto* hg  = stacked_to_th2d(good_stack, "hg_mean", "Good (media)", false);
    auto* hb  = stacked_to_th2d(bad_stack, "hb_mean", "Bad (media)", false);
    auto* hbg = stacked_to_th2d(bg_stack, "hbg_mean", "Background (media)", false);

    if (root_file) {
      hg->Write();
      hb->Write();
      hbg->Write();
    }

    // Stessa scala colori per confronto diretto: usiamo i percentili per ignorare outlier
    double vmin = std::min({get_th2d_percentile(hg, cfg.z_min_percentile),
                            get_th2d_percentile(hb, cfg.z_min_percentile),
                            get_th2d_percentile(hbg, cfg.z_min_percentile)});
    double vmax = std::max({get_th2d_percentile(hg, cfg.z_max_percentile),
                            get_th2d_percentile(hb, cfg.z_max_percentile),
                            get_th2d_percentile(hbg, cfg.z_max_percentile)});
    if (vmax <= vmin) {
      vmin = 0.0;
      vmax = std::max({hg->GetMaximum(), hb->GetMaximum(), hbg->GetMaximum()});
    }

    auto draw_pad = [&](int idx, TH2D* h, const std::string& lbl) {
      auto* pad = static_cast<TPad*>(c.GetPad(idx));
      pad->SetLeftMargin(0.16);
      pad->SetRightMargin(0.14);
      pad->SetTopMargin(0.11);
      pad->SetBottomMargin(0.14);
      pad->cd();

      // Rebinning per il canvas (risparmia memoria durante il Print PNG)
      TH2D* h_disp = rebin_for_display(h);
      h_disp->SetMinimum(vmin);
      h_disp->SetMaximum(vmax);
      h_disp->Draw("COLZ");

      TLatex t;
      t.SetNDC();
      t.SetTextFont(42);
      t.SetTextSize(0.055);
      t.SetTextAlign(22);
      t.DrawLatex(0.53, 0.952, lbl.c_str());
      pad->Update();
      return h_disp; // per cancellarlo dopo
    };

    auto* hgd  = draw_pad(1, hg, "Good (media #sigma-clip)");
    auto* hbd  = draw_pad(2, hb, "Bad  (media #sigma-clip)");
    auto* hbgd = draw_pad(3, hbg, "Background (media #sigma-clip)");

    c.Update();
    if (cfg.save_png)
      c.Print((cfg.output_dir / "stacked_means.png").string().c_str());

    delete hg;
    delete hb;
    delete hbg;
    delete hgd;
    delete hbd;
    delete hbgd;
  }

  // Pagina 2: mappe sigma

  {
    TCanvas c("c_sigma", "Sigma Maps", 1800, 600);
    c.Divide(3, 1, 0.005, 0.005);

    auto* hg  = stacked_to_th2d(good_stack, "hg_sig", "Good #sigma", true);
    auto* hb  = stacked_to_th2d(bad_stack, "hb_sig", "Bad #sigma", true);
    auto* hbg = stacked_to_th2d(bg_stack, "hbg_sig", "Bg #sigma", true);

    if (root_file) {
      hg->Write();
      hb->Write();
      hbg->Write();
    }

    // Scaling robusto per le mappe sigma
    double smin = std::min({get_th2d_percentile(hg, cfg.z_min_percentile),
                            get_th2d_percentile(hb, cfg.z_min_percentile),
                            get_th2d_percentile(hbg, cfg.z_min_percentile)});
    double smax = std::max({get_th2d_percentile(hg, cfg.z_max_percentile),
                            get_th2d_percentile(hb, cfg.z_max_percentile),
                            get_th2d_percentile(hbg, cfg.z_max_percentile)});
    if (smax <= smin) {
      smin = 0.0;
      smax = std::max({hg->GetMaximum(), hb->GetMaximum(), hbg->GetMaximum()});
    }

    auto draw_pad = [&](int idx, TH2D* h, const std::string& lbl) {
      auto* pad = static_cast<TPad*>(c.GetPad(idx));
      pad->SetLeftMargin(0.16);
      pad->SetRightMargin(0.14);
      pad->SetTopMargin(0.11);
      pad->SetBottomMargin(0.14);
      pad->cd();

      TH2D* h_disp = rebin_for_display(h);
      h_disp->SetMinimum(smin);
      h_disp->SetMaximum(smax);
      h_disp->Draw("COLZ");

      TLatex t;
      t.SetNDC();
      t.SetTextFont(42);
      t.SetTextSize(0.055);
      t.SetTextAlign(22);
      t.DrawLatex(0.53, 0.952, lbl.c_str());
      pad->Update();
      return h_disp;
    };

    auto* hgd  = draw_pad(1, hg, "#sigma_{good}  [ADU]");
    auto* hbd  = draw_pad(2, hb, "#sigma_{bad}   [ADU]");
    auto* hbgd = draw_pad(3, hbg, "#sigma_{bg}    [ADU]");

    c.Update();
    if (cfg.save_png)
      c.Print((cfg.output_dir / "sigma_maps.png").string().c_str());

    delete hg;
    delete hb;
    delete hbg;
    delete hgd;
    delete hbd;
    delete hbgd;
  }

  // Pagina 3: mappe differenza (signal - background)

  {
    TCanvas c("c_diff", "Difference Maps", 1400, 600);
    c.Divide(2, 1, 0.005, 0.005);

    auto* hgd = diff_to_th2d(good_diff, "hgd", "Good - Background", false);
    auto* hbd = diff_to_th2d(bad_diff, "hbd", "Bad  - Background", false);

    if (root_file) {
      hgd->Write();
      hbd->Write();
    }

    // Scala simmetrica attorno a 0 per enfatizzare il segnale netto
    double absmax = std::max({std::abs(hgd->GetMinimum()), hgd->GetMaximum(),
                              std::abs(hbd->GetMinimum()), hbd->GetMaximum()});
    gStyle->SetPalette(kCool); // palette divergente per valori negativi

    auto draw_pad = [&](int idx, TH2D* h, const std::string& lbl) {
      auto* pad = static_cast<TPad*>(c.GetPad(idx));
      pad->SetLeftMargin(0.16);
      pad->SetRightMargin(0.14);
      pad->SetTopMargin(0.11);
      pad->SetBottomMargin(0.14);
      pad->cd();

      TH2D* h_disp = rebin_for_display(h);
      h_disp->SetMinimum(-absmax);
      h_disp->SetMaximum(absmax);
      h_disp->GetZaxis()->SetTitle("[ADU]");
      h_disp->Draw("COLZ");

      TLatex t;
      t.SetNDC();
      t.SetTextFont(42);
      t.SetTextSize(0.055);
      t.SetTextAlign(22);
      t.DrawLatex(0.53, 0.952, lbl.c_str());
      pad->Update();
      return h_disp;
    };

    auto* hgdd = draw_pad(1, hgd, "Good #minus Background  [ADU]");
    auto* hbdd = draw_pad(2, hbd, "Bad  #minus Background  [ADU]");

    c.Update();
    if (cfg.save_png)
      c.Print((cfg.output_dir / "diff_maps.png").string().c_str());

    delete hgd;
    delete hbd;
    delete hgdd;
    delete hbdd;
  }

  // Pagina 4: profili integrati (proiezione X e Y)

  {
    TCanvas c("c_profiles", "Integrated Profiles", 1200, 800);
    c.Divide(1, 2, 0.005, 0.005);

    ROI roi = cfg.roi;

    auto* hgx = project_x(good_diff, "hgx", "Good: profilo X", roi);
    auto* hbx = project_x(bad_diff, "hbx", "Bad:  profilo X", roi);
    auto* hgy = project_y(good_diff, "hgy", "Good: profilo Y", roi);
    auto* hby = project_y(bad_diff, "hby", "Bad:  profilo Y", roi);

    if (root_file) {
      hgx->Write();
      hbx->Write();
      hgy->Write();
      hby->Write();
    }

    hgx->SetLineColor(kBlue + 1);
    hbx->SetLineColor(kRed + 1);
    hgy->SetLineColor(kBlue + 1);
    hby->SetLineColor(kRed + 1);
    for (auto* h : {hgx, hbx, hgy, hby})
      h->SetLineWidth(2);

    // Proiezione X: Good vs Bad sovrapposti
    auto setup_pad = [](TPad* p) {
      p->SetLeftMargin(0.16);
      p->SetRightMargin(0.05);
      p->SetTopMargin(0.10);
      p->SetBottomMargin(0.14);
      p->SetGridx();
      p->SetGridy();
    };

    setup_pad(static_cast<TPad*>(c.GetPad(1)));
    c.cd(1);
    double xmax = std::max(hgx->GetMaximum(), hbx->GetMaximum());
    double xmin = std::min(hgx->GetMinimum(), hbx->GetMinimum());
    hgx->SetMinimum(xmin - 0.05 * std::abs(xmin));
    hgx->SetMaximum(xmax + 0.05 * std::abs(xmax));
    hgx->Draw("HIST");
    hbx->Draw("HIST SAME");
    TLatex t1;
    t1.SetNDC();
    t1.SetTextFont(42);
    t1.SetTextSize(0.050);
    t1.SetTextAlign(22);
    t1.DrawLatex(0.53, 0.952, "Profilo integrato X (Good vs Bad)");

    setup_pad(static_cast<TPad*>(c.GetPad(2)));
    c.cd(2);
    double ymax = std::max(hgy->GetMaximum(), hby->GetMaximum());
    double ymin = std::min(hgy->GetMinimum(), hby->GetMinimum());
    hgy->SetMinimum(ymin - 0.05 * std::abs(ymin));
    hgy->SetMaximum(ymax + 0.05 * std::abs(ymax));
    hgy->Draw("HIST");
    hby->Draw("HIST SAME");
    TLatex t2;
    t2.SetNDC();
    t2.SetTextFont(42);
    t2.SetTextSize(0.050);
    t2.SetTextAlign(22);
    t2.DrawLatex(0.53, 0.952, "Profilo integrato Y (Good vs Bad)");

    c.Update();
    if (cfg.save_png)
      c.Print((cfg.output_dir / "profiles.png").string().c_str());

    delete hgx;
    delete hbx;
    delete hgy;
    delete hby;
  }

  // Pagina 5: pannello riassuntivo

  {
    TCanvas c("c_summary", "Summary", 900, 700);
    c.SetLeftMargin(0.02);
    c.SetRightMargin(0.02);
    c.SetTopMargin(0.02);
    c.SetBottomMargin(0.02);
    c.cd();

    // Sfondo grigio chiaro
    c.SetFillColor(TColor::GetColor(248, 248, 252));

    // Titolo
    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextSize(0.055);
    title.SetTextAlign(22);
    title.SetTextColor(TColor::GetColor(26, 58, 107)); // darkblue RIPTIDE
    title.DrawLatex(0.50, 0.92, "Exp1 — Analisi Statistica Immagini FITS");

    // Linea separatrice
    TLine sep(0.05, 0.88, 0.95, 0.88);
    sep.SetLineColor(TColor::GetColor(26, 58, 107));
    sep.SetLineWidth(2);
    sep.Draw();

    // Corpo testo
    TPaveText pt(0.05, 0.10, 0.95, 0.86, "NDC");
    pt.SetFillColor(0);
    pt.SetBorderSize(0);
    pt.SetTextFont(42);
    pt.SetTextSize(0.040);
    pt.SetTextAlign(12);

    auto add = [&](const std::string& s) { pt.AddText(s.c_str()); };

    add(" ");
    add("Stacking: #sigma-clipping (3#sigma, 3 iterazioni) su 10 frame per serie");
    add(" ");
    add("  Integrale netto Good  (segnale - fondo):");
    add(("    I_{good} = " + fmt_sci(comparison.good.integral, 4) + " ADU #cdot px   (#pm "
         + fmt_sci(comparison.good.sigma_integral, 3) + ")")
            .c_str());
    add(" ");
    add("  Integrale netto Bad   (segnale - fondo):");
    add(("    I_{bad}  = " + fmt_sci(comparison.bad.integral, 4) + " ADU #cdot px   (#pm "
         + fmt_sci(comparison.bad.sigma_integral, 3) + ")")
            .c_str());
    add(" ");
    add("  #bf{Risultati:}");
    add(("    #DeltaI = I_{good} - I_{bad} = " + fmt_sci(comparison.delta, 4) + " ADU #cdot px")
            .c_str());
    add(("    #sigma_{#DeltaI} = " + fmt_sci(comparison.sigma_delta, 3) + " ADU #cdot px").c_str());
    add(" ");
    add(("    #bf{Rapporto:} I_{good} / I_{bad} = " + fmtd(comparison.ratio, 3) + " #pm "
         + fmtd(comparison.sigma_ratio, 3))
            .c_str());
    add(("    #bf{Significativit#grave{a}:} (R-1)/#sigma_{R} = " + fmtd(comparison.significance, 2)
         + " #sigma")
            .c_str());

    if (fbf) {
      add(" ");
      add("  #bf{Metodo 2: Frame-to-Frame Fluctuations}");
      add(("    Rapporto medio R = " + fmtd(fbf->mean_ratio, 3) + " #pm "
           + fmtd(fbf->sigma_ratio_mean, 3))
              .c_str());
      add(("    Incertezza esp. #sigma_{R} = " + fmtd(fbf->sigma_ratio_std, 3) + " ("
           + std::to_string(fbf->n_frames_used) + " frame)")
              .c_str());
      add(("    Significativit#grave{a} S = " + fmtd(fbf->significance, 2) + " #sigma").c_str());
    }

    add(" ");
    add(("    N pixel ROI: " + std::to_string(comparison.good.n_pixels)).c_str());

    pt.Draw();

    c.Update();
    if (cfg.save_png)
      c.Print((cfg.output_dir / "summary.png").string().c_str());

    if (root_file) {
      c.Write("c_summary"); // Salviamo il canvas del riassunto nel file ROOT
    }
  }

  if (root_file) {
    root_file->Close();
  }

  std::cout << "\n[RESULT] ─────────────────────────────────────────────────────────\n";
  std::cout << "  I_good  = " << fmt_sci(comparison.good.integral, 4) << " ADU·px  (±"
            << fmt_sci(comparison.good.sigma_integral, 3) << ")\n";
  std::cout << "  I_bad   = " << fmt_sci(comparison.bad.integral, 4) << " ADU·px  (±"
            << fmt_sci(comparison.bad.sigma_integral, 3) << ")\n";
  std::cout << "  ΔI      = " << fmt_sci(comparison.delta, 4) << " ADU·px\n";
  std::cout << "  σ_ΔI    = " << fmt_sci(comparison.sigma_delta, 3) << " ADU·px\n";
  std::cout << "  Rapporto = " << fmtd(comparison.ratio, 3) << " ± "
            << fmtd(comparison.sigma_ratio, 3) << "\n";
  std::cout << "  Signif.  = " << fmtd(comparison.significance, 2) << " σ\n";

  if (fbf) {
    std::cout << "  [Metodo 2: Frame-to-Frame Fluctuations]\n";
    std::cout << "  R̄        = " << fmtd(fbf->mean_ratio, 3) << " ± "
              << fmtd(fbf->sigma_ratio_mean, 3) << "\n";
    std::cout << "  σ_R      = " << fmtd(fbf->sigma_ratio_std, 3) << " (dev. std sperimentale)\n";
    std::cout << "  N        = " << fbf->n_frames_used << " frame utilizzati\n";
    std::cout << "  Signif.  = " << fmtd(fbf->significance, 2) << " σ\n";
    std::cout << "  Valori R_k: ";
    for (size_t i = 0; i < fbf->ratios.size(); ++i) {
      std::cout << fmtd(fbf->ratios[i], 3) << (i == fbf->ratios.size() - 1 ? "" : ", ");
    }
    std::cout << "\n";
  }

  std::cout << "─────────────────────────────────────────────────────────────────\n";
  std::cout << "Output salvato in: " << cfg.output_dir.string() << "/\n";
}

} // namespace exp1