#pragma once

#include <TCanvas.h>
#include <TColor.h>
#include <TH2.h>
#include <TLatex.h>
#include <TPave.h>
#include <TPaletteAxis.h>
#include <TStyle.h>

#include <string>

// Stile ROOT condiviso e palette per tutti i grafici mappa (metriche di
// scelta: qualita', DoF, efficienza, risoluzione). Fonde le versioni un
// tempo duplicate in dof_stats_common.hpp/chi2_map.cpp/q_map.cpp/resolution_map.cpp.

inline void set_root_style() {
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
  gStyle->SetTitleOffset(1.20, "X");
  gStyle->SetTitleOffset(1.50, "Y");
  gStyle->SetTitleOffset(1.40, "Z");
  gStyle->SetTickLength(0.018, "X");
  gStyle->SetTickLength(0.018, "Y");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetNdivisions(505, "Y");
  gStyle->SetNumberContours(255);
  gStyle->SetGridColor(kGray + 1);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
}

// Palette direzionale condivisa per le mappe di scelta (qualita', DoF,
// efficienza, risoluzione). invert_for_low_is_good=true quando il valore
// migliore e' il piu' basso (es. chi2, errore assoluto, diametro spot).
inline void set_viridis_palette(bool invert_for_low_is_good) {
  gStyle->SetPalette(kViridis);
  gStyle->SetNumberContours(255);
  if (invert_for_low_is_good) {
    TColor::InvertPalette();
  }
}

// Palette divergente CVD-safe (ColorBrewer RdBu/PuOr) per le mappe con
// segno (errore firmato, margine, correlazione), al posto del
// rosso-bianco-verde non colorblind-friendly.
inline void set_diverging_palette_blue_white_orange() {
  const int nRGBs     = 3;
  double stops[nRGBs] = {0.0, 0.5, 1.0};
  double red[nRGBs]   = {33.0 / 255.0, 1.0, 230.0 / 255.0};
  double green[nRGBs] = {102.0 / 255.0, 1.0, 97.0 / 255.0};
  double blue[nRGBs]  = {172.0 / 255.0, 1.0, 1.0 / 255.0};
  TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, 255);
  gStyle->SetNumberContours(255);
}

// Canvas standard per le mappe (x1,x2) — layout unico condiviso da
// dof_map.cpp/magnification_map.cpp/resolution_map.cpp/chi2_map.cpp/q_map.cpp/plot2D.cpp.
inline TCanvas* make_map_canvas(const std::string& name, bool log_z = false) {
  TCanvas* c = new TCanvas(("c_" + name).c_str(), name.c_str(), 1100, 900);
  c->SetLeftMargin(0.14);
  c->SetBottomMargin(0.14);
  c->SetRightMargin(0.20);
  c->SetTopMargin(0.10);
  c->SetGridx();
  c->SetGridy();
  if (log_z) {
    c->SetLogz();
  }
  return c;
}

inline void apply_zaxis_style(TH2* h) {
  h->GetZaxis()->CenterTitle(kTRUE);
  h->GetZaxis()->SetTitleOffset(1.6);
}

// Banner titolo condiviso — da chiamare dopo Draw("COLZ").
inline void draw_map_title(const std::string& text) {
  TLatex tit;
  tit.SetNDC();
  tit.SetTextFont(42);
  tit.SetTextSize(0.042);
  tit.SetTextAlign(22);
  tit.DrawLatex(0.50, 0.955, text.c_str());
}

// Legenda "N/A" dedicata per le configurazioni invalide, agganciata alle
// coordinate NDC della TPaletteAxis auto-generata da Draw("COLZ"). Va
// chiamata dopo canvas->Update() (necessario perche' ROOT crei la
// TPaletteAxis prima che le sue coordinate NDC siano disponibili).
inline void draw_na_legend(TH2* h, Int_t invalid_color, const std::string& label = "N/A") {
  auto* pal =
      dynamic_cast<TPaletteAxis*>(h->GetListOfFunctions()->FindObject("palette"));
  if (!pal) {
    return;
  }
  double x1     = pal->GetX1NDC();
  double x2     = pal->GetX2NDC();
  double y1     = pal->GetY1NDC();
  double box_h  = 0.045;
  double gap    = 0.025;
  double box_y2 = y1 - gap;
  double box_y1 = box_y2 - box_h;
  if (box_y1 < 0.02) {
    box_y1 = 0.02;
    box_y2 = box_y1 + box_h;
  }

  auto* box = new TPave(x1, box_y1, x2, box_y2, 0, "NDC");
  box->SetFillColor(invalid_color);
  box->SetLineColor(invalid_color);
  box->Draw();

  TLatex lbl;
  lbl.SetNDC();
  lbl.SetTextFont(42);
  lbl.SetTextSize(0.028);
  lbl.SetTextAlign(22);
  lbl.SetTextColor(kWhite);
  lbl.DrawLatex((x1 + x2) / 2.0, (box_y1 + box_y2) / 2.0 - 0.01, label.c_str());
}
