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

/*
 * exp1_root_viewer — Ridisegna le immagini PNG a partire dal file .root di exp1
 *
 * Scopo: Permette di regolare la scala dei colori (percentili) senza dover
 * rieseguire l'intero stacking dei frame FITS (molto più veloce).
 */

#include <TCanvas.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

// Helper percentile (copiato da image_analysis.cpp)
static double get_th2d_percentile(TH2D* h, double percentile) {
  if (!h || percentile <= 0.0 || percentile >= 1.0)
    return 0.0;
  std::vector<double> vals;
  vals.reserve(static_cast<size_t>(h->GetNbinsX() * h->GetNbinsY()));
  for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
    for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
      double v = h->GetBinContent(ix, iy);
      if (v > 1e-9)
        vals.push_back(v);
    }
  }
  if (vals.empty())
    return 0.0;
  size_t idx = static_cast<size_t>(static_cast<double>(vals.size()) * percentile);
  std::nth_element(vals.begin(), vals.begin() + static_cast<long>(idx), vals.end());
  return vals[idx];
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "Uso: " << argv[0] << " <file.root> [z_min_perc] [z_max_perc]\n";
    std::cout << "Default: z_min=0.005, z_max=0.995\n";
    return 1;
  }

  std::string filename = argv[1];
  double z_min_perc    = (argc > 2) ? std::stod(argv[2]) : 0.005;
  double z_max_perc    = (argc > 3) ? std::stod(argv[3]) : 0.995;

  TFile* f = TFile::Open(filename.c_str(), "READ");
  if (!f || f->IsZombie())
    return 1;

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(kViridis);
  gStyle->SetGridColor(kGray + 1);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
  TGaxis::SetMaxDigits(4);

  // 1. Stacked Means
  TH2D* hg  = (TH2D*)f->Get("hg_mean");
  TH2D* hb  = (TH2D*)f->Get("hb_mean");
  TH2D* hbg = (TH2D*)f->Get("hbg_mean");

  hg->SetTitle("");
  hb->SetTitle("");
  hbg->SetTitle("");

  if (hg && hb && hbg) {
    TCanvas c("c_stacked", "Stacked", 1800, 600);
    c.Divide(3, 1);
    double vmin =
        std::min({get_th2d_percentile(hg, z_min_perc), get_th2d_percentile(hb, z_min_perc),
                  get_th2d_percentile(hbg, z_min_perc)});
    double vmax =
        std::max({get_th2d_percentile(hg, z_max_perc), get_th2d_percentile(hb, z_max_perc),
                  get_th2d_percentile(hbg, z_max_perc)});

    auto draw = [&](int i, TH2D* h, const char* t) {
      c.cd(i);
      gPad->SetLeftMargin(0.14);
      gPad->SetRightMargin(0.16);
      gPad->SetTopMargin(0.10);
      gPad->SetGridx();
      gPad->SetGridy();
      h->SetMinimum(vmin);
      h->SetMaximum(vmax);
      h->GetZaxis()->SetTitle("ADU");
      h->GetZaxis()->CenterTitle(kTRUE);
      h->GetZaxis()->SetTitleOffset(1.6);
      h->GetZaxis()->SetMaxDigits(2);
      h->Draw("COLZ");
      TLatex l;
      l.SetNDC();
      l.SetTextFont(42);
      l.SetTextSize(0.055);
      l.SetTextAlign(22);
      l.DrawLatex(0.50, 0.935, t);
    };
    draw(1, hg, "Good");
    draw(2, hb, "Bad");
    draw(3, hbg, "Background");
    c.Print("output/exp1/restacked_means.png");
  }

  // 2. Diff Maps
  TH2D* hgd = (TH2D*)f->Get("hgd");
  TH2D* hbd = (TH2D*)f->Get("hbd");
  if (hgd && hbd) {
    TCanvas c("c_diff", "Difference", 1200, 600);
    c.Divide(2, 1);
    gStyle->SetPalette(kViridis);
    double vmax =
        std::max(get_th2d_percentile(hgd, z_max_perc), get_th2d_percentile(hbd, z_max_perc));
    auto draw = [&](int i, TH2D* h, const char* t) {
      c.cd(i);
      gPad->SetLeftMargin(0.14);
      gPad->SetRightMargin(0.16);
      gPad->SetTopMargin(0.10);
      gPad->SetGridx();
      gPad->SetGridy();
      h->SetMinimum(0.0);
      h->SetMaximum(vmax);
      h->GetZaxis()->SetTitle("ADU  (#times10^{4})");
      h->GetZaxis()->CenterTitle(kTRUE);
      h->GetZaxis()->SetTitleOffset(1.6);
      h->Draw("COLZ");
      TLatex l;
      l.SetNDC();
      l.SetTextFont(42);
      l.SetTextSize(0.055);
      l.SetTextAlign(22);
      l.DrawLatex(0.50, 0.935, t);
    };
    draw(1, hgd, "Good #minus Background");
    draw(2, hbd, "Bad #minus Background");
    c.Print("output/exp1/diff_maps.png");
  }

  f->Close();
  return 0;
}
