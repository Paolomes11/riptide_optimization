/* Copyright 2026 Giulio Mesini — Apache License 2.0 */

#include "exp2_analysis.hpp"

#include "fits_io.hpp"

#include <TAxis.h>
#include <TBox.h>
#include <TGaxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TEllipse.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>

namespace {

static std::string fmtd(double v, int prec = 2) {
  std::ostringstream o;
  o << std::fixed << std::setprecision(prec) << v;
  return o.str();
}

static void apply_riptide_style() {
  gStyle->Reset();
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleFont(42, "");
  gStyle->SetStatFont(42);
  gStyle->SetTextSize(0.040f);
  gStyle->SetLabelSize(0.038f, "XYZ");
  gStyle->SetTitleSize(0.044f, "XYZ");
  gStyle->SetTitleSize(0.046f, "");
  gStyle->SetTitleOffset(1.55f, "Y");
  gStyle->SetTitleOffset(1.20f, "X");
  gStyle->SetTickLength(0.018f, "X");
  gStyle->SetTickLength(0.018f, "Y");
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
  gStyle->SetGridColor(kGray + 1);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
  TGaxis::SetMaxDigits(4);
}

static double vector_percentile(std::vector<double> vals, double p, bool ignore_zeros) {
  if (vals.empty() || p <= 0.0 || p >= 1.0)
    return 0.0;
  if (ignore_zeros) {
    vals.erase(std::remove_if(vals.begin(), vals.end(), [](double v) { return v <= 1e-12; }),
               vals.end());
  }
  if (vals.empty())
    return 0.0;
  const size_t idx = static_cast<size_t>(static_cast<double>(vals.size()) * p);
  const size_t i   = std::min(idx, vals.size() - 1);
  std::nth_element(vals.begin(), vals.begin() + static_cast<long>(i), vals.end());
  return vals[i];
}

static double get_th2d_percentile(TH2D* h, double percentile) {
  if (!h || percentile <= 0.0 || percentile >= 1.0)
    return 0.0;
  std::vector<double> vals;
  vals.reserve(static_cast<size_t>(h->GetNbinsX() * h->GetNbinsY()));
  for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
    for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
      const double v = h->GetBinContent(ix, iy);
      if (v > 1e-12)
        vals.push_back(v);
    }
  }
  if (vals.empty())
    return 0.0;
  const size_t idx = static_cast<size_t>(static_cast<double>(vals.size()) * percentile);
  const size_t i   = std::min(idx, vals.size() - 1);
  std::nth_element(vals.begin(), vals.begin() + static_cast<long>(i), vals.end());
  return vals[i];
}

static TH2D* rebin_for_display(TH2D* h) {
  if (!h)
    return nullptr;
  const int nx = h->GetNbinsX();
  const int ny = h->GetNbinsY();
  int rx       = 1;
  int ry       = 1;
  if (nx > 2000)
    rx = nx / 1024;
  if (ny > 2000)
    ry = ny / 1024;
  if (rx <= 1 && ry <= 1)
    return static_cast<TH2D*>(h->Clone());
  return static_cast<TH2D*>(h->Rebin2D(rx, ry, (std::string(h->GetName()) + "_disp").c_str()));
}

static TH2D* stacked_mean_to_th2d(const riptide::stack::StackedImage& img, const std::string& name,
                                  const std::string& title) {
  const int w = img.width;
  const int h = img.height;
  auto* hist  = new TH2D(name.c_str(), title.c_str(), w, 0, w, h, 0, h);
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      const double v =
          img.mean[static_cast<size_t>(y) * static_cast<size_t>(w) + static_cast<size_t>(x)];
      hist->SetBinContent(x + 1, y + 1, v);
    }
  }
  hist->GetXaxis()->SetTitle("x [px]");
  hist->GetYaxis()->SetTitle("y [px]");
  hist->GetXaxis()->CenterTitle(true);
  hist->GetYaxis()->CenterTitle(true);
  return hist;
}

static TH2D* diff_to_th2d(const std::vector<double>& diff, int width, int height,
                          const std::string& name, const std::string& title) {
  auto* hist = new TH2D(name.c_str(), title.c_str(), width, 0, width, height, 0, height);
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      const double v =
          diff[static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)];
      hist->SetBinContent(x + 1, y + 1, v);
    }
  }
  hist->GetXaxis()->SetTitle("x [px]");
  hist->GetYaxis()->SetTitle("y [px]");
  hist->GetXaxis()->CenterTitle(true);
  hist->GetYaxis()->CenterTitle(true);
  return hist;
}

static void set_pad_margins(TPad* pad) {
  pad->SetLeftMargin(0.16f);
  pad->SetRightMargin(0.14f);
  pad->SetTopMargin(0.11f);
  pad->SetBottomMargin(0.14f);
}

static void set_diverging_palette() {
  static bool done = false;
  if (done)
    return;
  done             = true;
  Double_t stops[] = {0.0, 0.5, 1.0};
  Double_t r[]     = {0.15, 1.0, 0.75};
  Double_t g[]     = {0.30, 1.0, 0.15};
  Double_t b[]     = {0.75, 1.0, 0.15};
  TColor::CreateGradientColorTable(3, stops, r, g, b, 255);
  gStyle->SetNumberContours(255);
}

static void draw_trace_line(int width, int height, double angle_deg) {
  const double cx    = 0.5 * static_cast<double>(width - 1);
  const double cy    = 0.5 * static_cast<double>(height - 1);
  const double theta = angle_deg * (M_PI / 180.0);
  const double ux    = std::cos(theta);
  const double uy    = std::sin(theta);
  const double L     = 0.8 * static_cast<double>(std::min(width, height));
  const double dx    = 0.5 * L * ux;
  const double dy    = 0.5 * L * uy;

  TLine line(cx - dx, cy - dy, cx + dx, cy + dy);
  line.SetLineColor(kRed);
  line.SetLineWidth(2);
  line.Draw("SAME");
}

static riptide::stack::StackedImage
stack_frames(const std::vector<riptide::fits::FitsFrame>& frames,
             const riptide::stack::StackConfig& cfg) {
  if (cfg.method == riptide::stack::StackMethod::Mean)
    return riptide::stack::mean_stack(frames);
  if (cfg.method == riptide::stack::StackMethod::Median)
    return riptide::stack::median_stack(frames);
  return riptide::stack::sigma_clip_stack(frames, cfg);
}

static double compute_pull(const riptide::exp2::SliceProfile& sp,
                           const riptide::exp2::CentroidFitResult& fit, double var_z) {
  const double y = sp.center;
  const double z = sp.t;
  const double var_y =
      (std::isfinite(sp.center_err) && sp.center_err > 0.0) ? (sp.center_err * sp.center_err) : 1.0;

  double u     = 0.0;
  double v     = 0.0;
  double var_u = 0.0;
  double var_v = 0.0;

  if (fit.axis == riptide::FitAxis::ZvsY) {
    u     = y;
    v     = z;
    var_u = var_y;
    var_v = var_z;
  } else {
    u     = z;
    v     = y;
    var_u = var_z;
    var_v = var_y;
  }

  const double denom   = std::sqrt(1.0 + fit.a * fit.a);
  const double resid   = (v - fit.a * u - fit.b) / denom;
  const double sigma_d = std::sqrt((fit.a * fit.a * var_u + var_v) / (1.0 + fit.a * fit.a));
  if (!(sigma_d > 0.0) || !std::isfinite(sigma_d))
    return 0.0;
  return resid / sigma_d;
}

static void produce_blob_image(const riptide::exp2::ConfigResult& result,
                               const riptide::exp2::OutputConfig& cfg) {
  const std::string lbl = riptide::exp2::config_label_str(result.label);
  const int w           = result.signal_stack.width;
  const int h           = result.signal_stack.height;
  if (w <= 0 || h <= 0)
    return;

  // Centro di massa del blob (solo valori positivi)
  double wsum = 0.0, cx = 0.0, cy = 0.0;
  for (int iy = 0; iy < h; ++iy) {
    for (int ix = 0; ix < w; ++ix) {
      const double v =
          std::max(0.0, result.diff[static_cast<size_t>(iy) * static_cast<size_t>(w)
                                    + static_cast<size_t>(ix)]);
      wsum += v;
      cx += ix * v;
      cy += iy * v;
    }
  }
  if (!(wsum > 0.0))
    return;
  cx /= wsum;
  cy /= wsum;

  // Regione di crop: ±4*sigma_major (minimo 40 px)
  const double margin = std::max(40.0, 4.0 * result.trace.sigma_major);
  const int x0        = std::max(0, static_cast<int>(cx - margin));
  const int x1        = std::min(w - 1, static_cast<int>(cx + margin));
  const int y0        = std::max(0, static_cast<int>(cy - margin));
  const int y1        = std::min(h - 1, static_cast<int>(cy + margin));
  if (x1 <= x0 || y1 <= y0)
    return;

  gStyle->SetPalette(kViridis);
  TCanvas c(("c_blob_" + lbl).c_str(), ("Blob " + lbl).c_str(), 900, 800);
  set_pad_margins(static_cast<TPad*>(&c));

  auto* hc = new TH2D(("h_blob_" + lbl).c_str(), "",
                      x1 - x0 + 1, static_cast<double>(x0), static_cast<double>(x1 + 1),
                      y1 - y0 + 1, static_cast<double>(y0), static_cast<double>(y1 + 1));
  for (int iy = y0; iy <= y1; ++iy) {
    for (int ix = x0; ix <= x1; ++ix) {
      hc->SetBinContent(ix - x0 + 1, iy - y0 + 1,
                        result.diff[static_cast<size_t>(iy) * static_cast<size_t>(w)
                                    + static_cast<size_t>(ix)]);
    }
  }
  hc->GetXaxis()->SetTitle("x [px]");
  hc->GetYaxis()->SetTitle("y [px]");
  hc->GetXaxis()->CenterTitle(true);
  hc->GetYaxis()->CenterTitle(true);
  hc->GetZaxis()->SetTitle("ADU");
  hc->GetZaxis()->CenterTitle(kTRUE);
  hc->GetZaxis()->SetTitleOffset(1.6);

  const double vmax = get_th2d_percentile(hc, cfg.z_max_percentile);
  hc->SetMinimum(0.0);
  if (vmax > 0.0)
    hc->SetMaximum(vmax);
  hc->Draw("COLZ");

  // Ellisse 2σ (semiassi sigma_major, sigma_minor)
  if (result.trace.sigma_minor > 0.0 && result.trace.sigma_major > 0.0) {
    TEllipse ell(cx, cy,
                 2.0 * result.trace.sigma_major, 2.0 * result.trace.sigma_minor,
                 0.0, 360.0,
                 result.trace.angle_deg);
    ell.SetLineColor(kRed);
    ell.SetLineWidth(2);
    ell.SetFillStyle(0);
    ell.Draw("SAME");
  }

  // Linea della traccia attraverso il centroide
  {
    const double theta = result.trace.angle_deg * (M_PI / 180.0);
    const double ux    = std::cos(theta);
    const double uy    = std::sin(theta);
    const double L     = 0.8 * std::min(x1 - x0, y1 - y0);
    const double dx    = 0.5 * L * ux;
    const double dy    = 0.5 * L * uy;
    TLine tline(cx - dx, cy - dy, cx + dx, cy + dy);
    tline.SetLineColor(kOrange + 1);
    tline.SetLineWidth(2);
    tline.SetLineStyle(2);
    tline.Draw("SAME");
  }

  TLatex main_t;
  main_t.SetNDC();
  main_t.SetTextFont(42);
  main_t.SetTextSize(0.042f);
  main_t.SetTextAlign(22);
  main_t.DrawLatex(0.50, 0.972f, ("Profilo 2D del blob  [" + lbl + "]").c_str());

  TLatex lbl_txt;
  lbl_txt.SetNDC();
  lbl_txt.SetTextFont(42);
  lbl_txt.SetTextSize(0.034f);
  lbl_txt.SetTextAlign(11);
  lbl_txt.DrawLatex(0.18, 0.90,
                    ("#sigma_{minor}=" + fmtd(result.trace.sigma_minor, 1)
                     + "  #sigma_{major}=" + fmtd(result.trace.sigma_major, 1)
                     + " px  #theta=" + fmtd(result.trace.angle_deg, 1) + "#circ")
                        .c_str());

  if (cfg.save_png)
    c.Print((cfg.output_dir / ("blob_" + lbl + ".png")).string().c_str());
}

} // namespace

namespace riptide::exp2 {

std::string config_label_str(ConfigLabel label) {
  switch (label) {
  case ConfigLabel::GoodFocus:
    return "good_focus";
  case ConfigLabel::GoodNoFocus:
    return "good_nofocus";
  case ConfigLabel::BadFocus:
    return "bad_focus";
  case ConfigLabel::BadNoFocus:
    return "bad_nofocus";
  default:
    return "unknown";
  }
}

ConfigResult analyze_config(const std::filesystem::path& signal_dir,
                            const std::filesystem::path& background_dir, ConfigLabel label,
                            const OpticsParams& optics,
                            const riptide::stack::StackConfig& stack_cfg,
                            const TraceConfig& trace_cfg) {
  auto signal_frames = riptide::fits::read_fits_stack(signal_dir, stack_cfg.max_frames);
  auto bg_frames     = riptide::fits::read_fits_stack(background_dir, stack_cfg.max_frames);

  auto signal_stack = stack_frames(signal_frames, stack_cfg);
  auto bg_stack     = stack_frames(bg_frames, stack_cfg);

  if (signal_stack.width != bg_stack.width || signal_stack.height != bg_stack.height)
    throw std::runtime_error("exp2::analyze_config: dimensioni diverse tra signal e background");

  const int width  = signal_stack.width;
  const int height = signal_stack.height;
  const size_t n   = static_cast<size_t>(width) * static_cast<size_t>(height);

  std::vector<double> diff(n, 0.0);
  for (size_t i = 0; i < n; ++i)
    diff[i] = signal_stack.mean[i] - bg_stack.mean[i];

  double angle_err    = 0.0;
  double current_angle = estimate_trace_angle(diff, width, height, trace_cfg, &angle_err);

  TraceResult trace   = extract_trace_profile(diff, width, height, current_angle, trace_cfg);
  trace.angle_deg     = current_angle;
  trace.angle_err_deg = angle_err;

  CentroidFitResult centroid_fit{};
  if (trace.trace_detected) {
    centroid_fit = fit_centroid_line(trace);

    // P2: raffinamento iterativo dell'angolo (fino a 2 passaggi).
    // Il pendio ODR rivela l'errore angolare residuo: delta ≈ atan(slope).
    // Ref: tecnica standard nei beam profiler (ISO 11146) e tracker di tracce.
    for (int refine = 0; refine < 2; ++refine) {
      if (!centroid_fit.converged || centroid_fit.n_points_used < 20)
        break;

      double delta_deg = 0.0;
      if (centroid_fit.axis == riptide::FitAxis::YvsZ) {
        // center = a·t + b  →  a = tan(δ)  →  δ = atan(a)
        delta_deg = (180.0 / M_PI) * std::atan(centroid_fit.a);
      } else {
        // t = a·center + b  →  1/a = tan(δ)  →  δ = atan(1/a)
        if (std::abs(centroid_fit.a) > 1e-6)
          delta_deg = (180.0 / M_PI) * std::atan(1.0 / centroid_fit.a);
      }

      // Ignora correzioni trascurabili (<0.05°) o anomale (>10°).
      if (std::abs(delta_deg) < 0.05 || std::abs(delta_deg) > 10.0)
        break;

      current_angle += delta_deg;
      while (current_angle > 90.0)
        current_angle -= 180.0;
      while (current_angle <= -90.0)
        current_angle += 180.0;

      TraceResult refined = extract_trace_profile(diff, width, height, current_angle, trace_cfg);
      refined.angle_deg     = current_angle;
      refined.angle_err_deg = angle_err;
      if (!refined.trace_detected)
        break;

      trace        = std::move(refined);
      centroid_fit = fit_centroid_line(trace);
    }
  }

  const double dof_delta = (optics.x_det_optimal > 0.0) ? (optics.x_det - optics.x_det_optimal)
                                                        : std::numeric_limits<double>::quiet_NaN();

  ConfigResult res;
  res.label        = label;
  res.optics       = optics;
  res.signal_stack = std::move(signal_stack);
  res.diff         = std::move(diff);
  res.trace        = std::move(trace);
  res.centroid_fit = std::move(centroid_fit);
  res.dof_delta_mm = dof_delta;
  return res;
}

void produce_config_output(const ConfigResult& result, const OutputConfig& cfg) {
  std::filesystem::create_directories(cfg.output_dir);
  apply_riptide_style();

  std::unique_ptr<TFile> root_file;
  if (cfg.save_root) {
    root_file =
        std::make_unique<TFile>((cfg.output_dir / "exp2_analysis.root").string().c_str(), "UPDATE");
    if (root_file->IsZombie()) {
      std::cerr << "[WARNING] Impossibile creare il file ROOT di output\n";
      root_file.reset();
    }
  }

  const std::string lbl = config_label_str(result.label);
  const int width       = result.signal_stack.width;
  const int height      = result.signal_stack.height;

  {
    gStyle->SetPalette(kViridis);
    TCanvas c(("c_stacked_" + lbl).c_str(), ("Stacked " + lbl).c_str(), 900, 800);
    set_pad_margins(static_cast<TPad*>(&c));

    auto* h = stacked_mean_to_th2d(result.signal_stack, "h_stacked_" + lbl, "Stacked mean");
    if (root_file) {
      root_file->cd();
      h->Write();
    }

    const double vmax = get_th2d_percentile(h, cfg.z_max_percentile);

    TH2D* h_disp = rebin_for_display(h);
    if (vmax > 0.0) {
      h_disp->SetMinimum(0.0);
      h_disp->SetMaximum(vmax);
    }

    h_disp->Draw("COLZ");
    draw_trace_line(width, height, result.trace.angle_deg);

    if (cfg.save_png)
      c.Print((cfg.output_dir / ("stacked_" + lbl + ".png")).string().c_str());
    if (root_file) {
      root_file->cd();
      c.Write(("stacked_canvas_" + lbl).c_str());
    }
  }

  {
    set_diverging_palette();
    TCanvas c(("c_diff_" + lbl).c_str(), ("Diff " + lbl).c_str(), 900, 800);
    set_pad_margins(static_cast<TPad*>(&c));

    auto* h = diff_to_th2d(result.diff, width, height, "h_diff_" + lbl, "Signal - Background");
    if (root_file) {
      root_file->cd();
      h->Write();
    }

    std::vector<double> tmp = result.diff;
    const double p5         = vector_percentile(tmp, 0.05, false);
    const double p95        = vector_percentile(tmp, 0.95, false);
    const double v          = std::max(std::abs(p5), std::abs(p95));

    TH2D* h_disp = rebin_for_display(h);
    if (v > 0.0 && std::isfinite(v)) {
      h_disp->SetMinimum(-v);
      h_disp->SetMaximum(v);
    }
    h_disp->Draw("COLZ");
    draw_trace_line(width, height, result.trace.angle_deg);

    if (cfg.save_png)
      c.Print((cfg.output_dir / ("diff_" + lbl + ".png")).string().c_str());
    if (root_file) {
      root_file->cd();
      c.Write(("diff_canvas_" + lbl).c_str());
    }
  }

  {
    TCanvas c(("c_profile_" + lbl).c_str(), ("Profile " + lbl).c_str(), 900, 650);
    set_pad_margins(static_cast<TPad*>(&c));

    std::vector<double> t;
    std::vector<double> s;
    std::vector<double> se;
    for (const auto& sp : result.trace.profile) {
      if (sp.valid && std::isfinite(sp.sigma)) {
        t.push_back(sp.t);
        s.push_back(sp.sigma);
        se.push_back(std::isfinite(sp.sigma_err) ? sp.sigma_err : 0.0);
      }
    }

    TGraphErrors* g = nullptr;
    if (t.empty()) {
      TLatex nodata;
      nodata.SetNDC();
      nodata.SetTextFont(42);
      nodata.SetTextSize(0.050f);
      nodata.SetTextAlign(22);
      nodata.DrawLatex(0.5, 0.5, "Nessuna slice valida");
    } else {
      g = new TGraphErrors(static_cast<int>(t.size()));
      for (int i = 0; i < static_cast<int>(t.size()); ++i) {
        g->SetPoint(i, t[static_cast<size_t>(i)], s[static_cast<size_t>(i)]);
        g->SetPointError(i, 0.0, se[static_cast<size_t>(i)]);
      }
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.9f);
      g->SetLineWidth(2);
      g->GetXaxis()->SetTitle("t [px]");
      g->GetYaxis()->SetTitle("#sigma(t) [px]");
      g->Draw("AP");

      const double tmin = *std::min_element(t.begin(), t.end());
      const double tmax = *std::max_element(t.begin(), t.end());

      TLine mean_line(tmin, result.trace.sigma_mean, tmax, result.trace.sigma_mean);
      mean_line.SetLineColor(kRed + 1);
      mean_line.SetLineStyle(2);
      mean_line.SetLineWidth(2);
      mean_line.Draw("SAME");

      TBox box(tmin, result.trace.sigma_mean - result.trace.sigma_mean_err, tmax,
               result.trace.sigma_mean + result.trace.sigma_mean_err);
      box.SetFillColorAlpha(kRed + 1, 0.15f);
      box.SetLineColor(kRed + 1);
      box.SetLineWidth(1);
      box.Draw("SAME");

      TLatex txt;
      txt.SetNDC();
      txt.SetTextFont(42);
      txt.SetTextSize(0.040f);
      txt.SetTextAlign(13);
      txt.DrawLatex(0.18, 0.93,
                    ("#sigma_{mean} = " + fmtd(result.trace.sigma_mean, 2) + " #pm "
                     + fmtd(result.trace.sigma_mean_err, 2)
                     + " px    FWHM = " + fmtd(result.trace.fwhm_mean, 2)
                     + " px    N = " + std::to_string(result.trace.n_valid_slices) + " slice")
                        .c_str());
    }

    if (cfg.save_png)
      c.Print((cfg.output_dir / ("profile_" + lbl + ".png")).string().c_str());
    if (root_file) {
      root_file->cd();
      if (g)
        g->Write(("sigma_profile_" + lbl).c_str());
      c.Write(("profile_canvas_" + lbl).c_str());
    }
  }

  {
    TCanvas c(("c_centroid_" + lbl).c_str(), ("Centroid " + lbl).c_str(), 900, 750);

    TLatex c_title;
    c_title.SetNDC();
    c_title.SetTextFont(42);
    c_title.SetTextSize(0.038f);
    c_title.SetTextAlign(22);
    c_title.DrawLatex(0.50f, 0.973f, ("Centroide della traccia  [" + lbl + "]").c_str());

    auto* pad_top = new TPad(("pad_top_" + lbl).c_str(), "", 0.0, 0.30, 1.0, 1.0);
    auto* pad_bot = new TPad(("pad_bot_" + lbl).c_str(), "", 0.0, 0.00, 1.0, 0.30);
    pad_top->SetBottomMargin(0.02f);
    pad_bot->SetTopMargin(0.02f);
    pad_bot->SetBottomMargin(0.35f);
    pad_top->SetLeftMargin(0.12f);
    pad_top->SetRightMargin(0.06f);
    pad_top->SetTopMargin(0.06f);
    pad_bot->SetLeftMargin(0.12f);
    pad_bot->SetRightMargin(0.06f);
    pad_top->SetGridx();
    pad_top->SetGridy();

    pad_top->Draw();
    pad_bot->Draw();

    std::vector<double> dt;
    dt.reserve(result.trace.profile.size());
    for (size_t i = 1; i < result.trace.profile.size(); ++i)
      dt.push_back(std::abs(result.trace.profile[i].t - result.trace.profile[i - 1].t));
    const double dt_med = dt.empty() ? 1.0 : vector_percentile(dt, 0.5, false);
    const double var_z  = 0.25 * dt_med * dt_med;

    std::vector<double> t;
    std::vector<double> y;
    std::vector<double> ye;
    std::vector<double> pull;
    for (const auto& sp : result.trace.profile) {
      if (sp.valid && !sp.near_edge && std::isfinite(sp.center)) {
        t.push_back(sp.t);
        y.push_back(sp.center);
        ye.push_back((std::isfinite(sp.center_err) && sp.center_err > 0.0) ? sp.center_err : 0.0);
        pull.push_back(compute_pull(sp, result.centroid_fit, var_z));
      }
    }

    pad_top->cd();
    auto* g = new TGraphErrors(static_cast<int>(t.size()));
    for (int i = 0; i < static_cast<int>(t.size()); ++i) {
      g->SetPoint(i, t[static_cast<size_t>(i)], y[static_cast<size_t>(i)]);
      g->SetPointError(i, 0.0, ye[static_cast<size_t>(i)]);
    }
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.9f);
    g->SetLineWidth(2);
    g->GetXaxis()->SetLabelSize(0.0);
    g->GetXaxis()->SetTitleSize(0.0);
    g->GetYaxis()->SetTitle("centroide [px]");
    g->GetYaxis()->SetTitleSize(0.06f);
    g->GetYaxis()->SetTitleOffset(0.9f);
    g->Draw("AP");

    const double tmin = t.empty() ? 0.0 : *std::min_element(t.begin(), t.end());
    const double tmax = t.empty() ? 1.0 : *std::max_element(t.begin(), t.end());

    double y1 = 0.0;
    double y2 = 0.0;
    if (result.centroid_fit.axis == riptide::FitAxis::YvsZ) {
      y1 = result.centroid_fit.a * tmin + result.centroid_fit.b;
      y2 = result.centroid_fit.a * tmax + result.centroid_fit.b;
    } else {
      if (std::abs(result.centroid_fit.a) > 1e-12) {
        y1 = (tmin - result.centroid_fit.b) / result.centroid_fit.a;
        y2 = (tmax - result.centroid_fit.b) / result.centroid_fit.a;
      } else {
        y1 = 0.0;
        y2 = 0.0;
      }
    }
    TLine fit_line(tmin, y1, tmax, y2);
    fit_line.SetLineColor(kRed);
    fit_line.SetLineWidth(2);
    fit_line.Draw("SAME");

    TLatex txt;
    txt.SetNDC();
    txt.SetTextFont(42);
    txt.SetTextSize(0.045f);
    txt.SetTextAlign(13);
    txt.DrawLatex(0.18, 0.82,
                  ("#chi^{2}/ndof = " + fmtd(result.centroid_fit.chi2_ndof, 3)
                   + "  (ndof = " + std::to_string(result.centroid_fit.ndof) + ")")
                      .c_str());

    pad_bot->cd();
    auto* gp = new TGraph(static_cast<int>(t.size()));
    for (int i = 0; i < static_cast<int>(t.size()); ++i)
      gp->SetPoint(i, t[static_cast<size_t>(i)], pull[static_cast<size_t>(i)]);
    gp->SetMarkerStyle(20);
    gp->SetMarkerSize(0.7f);
    gp->SetLineWidth(2);
    gp->GetXaxis()->SetTitle("t  [px]");
    gp->GetXaxis()->SetTitleSize(0.10f);
    gp->GetXaxis()->SetTitleOffset(1.1f);
    gp->GetYaxis()->SetTitle("pull  [a.d.]");
    gp->GetYaxis()->SetTitleSize(0.10f);
    gp->GetYaxis()->SetTitleOffset(0.55f);
    gp->GetYaxis()->SetNdivisions(505);
    gp->Draw("AP");

    TLine l0(tmin, 0.0, tmax, 0.0);
    l0.SetLineColor(kBlack);
    l0.SetLineWidth(2);
    l0.Draw("SAME");

    for (double yy : {1.0, -1.0, 2.0, -2.0}) {
      TLine l(tmin, yy, tmax, yy);
      l.SetLineColor(kGray + 2);
      l.SetLineStyle(2);
      l.SetLineWidth(2);
      l.Draw("SAME");
    }

    if (cfg.save_png)
      c.Print((cfg.output_dir / ("centroid_" + lbl + ".png")).string().c_str());
    if (root_file) {
      root_file->cd();
      g->Write(("centroid_graph_" + lbl).c_str());
      gp->Write(("centroid_pull_" + lbl).c_str());
      c.Write(("centroid_canvas_" + lbl).c_str());
    }
  }

  produce_blob_image(result, cfg);
}

void produce_summary(const std::vector<ConfigResult>& results, const OutputConfig& cfg) {
  std::filesystem::create_directories(cfg.output_dir);
  apply_riptide_style();

  std::unique_ptr<TFile> root_file;
  if (cfg.save_root) {
    root_file =
        std::make_unique<TFile>((cfg.output_dir / "exp2_analysis.root").string().c_str(), "UPDATE");
    if (root_file->IsZombie())
      root_file.reset();
  }

  if (results.size() < 4)
    throw std::invalid_argument("produce_summary: risultati insufficienti");

  auto find = [&](ConfigLabel l) -> const ConfigResult* {
    for (const auto& r : results)
      if (r.label == l)
        return &r;
    return nullptr;
  };

  const auto* good_focus   = find(ConfigLabel::GoodFocus);
  const auto* good_nofocus = find(ConfigLabel::GoodNoFocus);
  const auto* bad_focus    = find(ConfigLabel::BadFocus);
  const auto* bad_nofocus  = find(ConfigLabel::BadNoFocus);

  int ref_width = 0;
  for (const auto& r : results)
    ref_width = std::max(ref_width, r.signal_stack.width);

  // sigma_minor è l'asse minore della distribuzione 2D dell'immagine:
  // invariante per rotazione, non influenzato dal trimming, confrontabile tra blob e streak.
  auto sigma_scaled = [&](const ConfigResult& r) -> std::pair<double, double> {
    if (!(r.trace.sigma_minor > 0.0))
      return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    if (ref_width <= 0 || r.signal_stack.width <= 0)
      return {r.trace.sigma_minor, 0.0};
    const double scale = static_cast<double>(ref_width) / static_cast<double>(r.signal_stack.width);
    return {r.trace.sigma_minor * scale, 0.0};
  };

  double sigma_min  = std::numeric_limits<double>::infinity();
  double sigma_max  = -std::numeric_limits<double>::infinity();
  ConfigLabel best  = results.front().label;
  ConfigLabel worst = results.front().label;

  for (const auto& r : results) {
    const auto s = sigma_scaled(r);
    if (std::isfinite(s.first) && s.first < sigma_min) {
      sigma_min = s.first;
      best      = r.label;
    }
    if (std::isfinite(s.first) && s.first > sigma_max) {
      sigma_max = s.first;
      worst     = r.label;
    }
  }

  TCanvas c("c_summary", "exp2 summary", 900, 700);
  auto* pad = static_cast<TPad*>(&c);
  set_pad_margins(pad);

  TPaveText box(0.06, 0.10, 0.94, 0.94, "NDC");
  box.SetFillColor(0);
  box.SetTextFont(42);
  box.SetTextSize(0.032f);
  box.SetTextAlign(12);
  box.SetBorderSize(0);

  box.AddText("exp2 — Summary (parametri | risultati)");
  box.AddText(" ");

  for (const auto& r : results) {
    const bool is_best  = (r.label == best);
    const bool is_worst = (r.label == worst);
    const int col       = is_best ? (kGreen + 2) : (is_worst ? (kRed + 1) : kBlack);

    std::ostringstream l1;
    l1 << "[" << config_label_str(r.label) << "]"
       << "    x1=" << fmtd(r.optics.x1, 1) << " mm"
       << "  x2=" << fmtd(r.optics.x2, 1) << " mm"
       << "  x_det=" << fmtd(r.optics.x_det, 1) << " mm";

    const std::string blob_tag = r.trace.trace_detected ? "" : "  [blob, no trace]";

    std::ostringstream l2;
    l2 << "            "
       << "#sigma_{minor}=" << fmtd(r.trace.sigma_minor, 2) << " px"
       << "   #sigma_{major}=" << fmtd(r.trace.sigma_major, 2) << " px"
       << "   aspect=" << fmtd(r.trace.aspect_ratio, 2) << blob_tag;

    std::ostringstream l3;
    if (r.trace.trace_detected)
      l3 << "            "
         << "#sigma_mean=" << fmtd(r.trace.sigma_mean, 2) << " #pm "
         << fmtd(r.trace.sigma_mean_err, 2) << " px   FWHM=" << fmtd(r.trace.fwhm_mean, 2)
         << " px   N=" << r.trace.n_valid_slices;
    else
      l3 << "            "
         << "#sigma_mean=N/A (traccia non rilevata)   aspect_ratio="
         << fmtd(r.trace.aspect_ratio, 2);

    std::ostringstream l4;
    l4 << "            ";
    if (r.trace.trace_detected)
      l4 << "#chi^{2}/ndof(centroide)=" << fmtd(r.centroid_fit.chi2_ndof, 2);
    else
      l4 << "#chi^{2}/ndof(centroide)=N/A";
    l4 << "   x_det_opt="
       << ((r.optics.x_det_optimal > 0.0) ? fmtd(r.optics.x_det_optimal, 1) : std::string("NaN"))
       << " mm"
       << "   #Deltax_det=" << (std::isfinite(r.dof_delta_mm) ? fmtd(r.dof_delta_mm, 2) : "NaN")
       << " mm";

    auto* t1 = box.AddText(l1.str().c_str());
    t1->SetTextColor(static_cast<Color_t>(col));
    auto* t2 = box.AddText(l2.str().c_str());
    t2->SetTextColor(static_cast<Color_t>(col));
    auto* t3 = box.AddText(l3.str().c_str());
    t3->SetTextColor(static_cast<Color_t>(col));
    auto* t4 = box.AddText(l4.str().c_str());
    t4->SetTextColor(static_cast<Color_t>(col));
    box.AddText(" ");
  }

  auto ratio_line = [&](const std::string& name, const ConfigResult* num, const ConfigResult* den) {
    if (!num || !den)
      return;
    const auto sn = sigma_scaled(*num);
    const auto sd = sigma_scaled(*den);
    if (!(sn.first > 0.0) || !(sd.first > 0.0))
      return;
    const double r = sn.first / sd.first;
    const double sr =
        r * std::sqrt(std::pow(sn.second / sn.first, 2.0) + std::pow(sd.second / sd.first, 2.0));
    const double sig = (sr > 0.0) ? ((r - 1.0) / sr) : 0.0;
    std::ostringstream o;
    o << name << " = " << fmtd(r, 3) << " ± " << fmtd(sr, 3) << "   (" << fmtd(sig, 2)
      << " #sigma)";
    box.AddText(o.str().c_str());
  };

  box.AddText("Rapporti chiave (#sigma_{minor}):");
  ratio_line("#sigma_{minor}(bad,focus) / #sigma_{minor}(good,focus)", bad_focus, good_focus);
  ratio_line("#sigma_{minor}(good,focus) / #sigma_{minor}(good,nofocus)", good_focus, good_nofocus);
  ratio_line("#sigma_{minor}(bad,focus) / #sigma_{minor}(bad,nofocus)", bad_focus, bad_nofocus);

  box.Draw();

  if (cfg.save_png)
    c.Print((cfg.output_dir / "summary.png").string().c_str());
  if (root_file) {
    root_file->cd();
    c.Write("summary_canvas");
  }
}

} // namespace riptide::exp2
