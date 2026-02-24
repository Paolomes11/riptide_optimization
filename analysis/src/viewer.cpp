#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TTree.h>

#include <iostream>

int main() {
  // Disabilita il box delle statistiche
  gStyle->SetOptStat(0);

  // Percorso del file ROOT
  const char* filename = "output/output.root";

  // Apri file
  TFile* file = TFile::Open(filename, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Errore: impossibile aprire " << filename << std::endl;
    return 1;
  }

  // Ottieni il tree
  TTree* tree = nullptr;
  file->GetObject("Hits", tree);

  if (!tree) {
    std::cerr << "Errore: TTree 'Hits' non trovato." << std::endl;
    return 1;
  }

  // Canvas
  TCanvas* c1 = new TCanvas("c1", "Hits in detector", 800, 700);

  // Abilita griglia
  c1->SetGrid();
  gStyle->SetGridStyle(3); // linee tratteggiate
  gStyle->SetGridColor(kGray);
  gStyle->SetGridWidth(1);

  // Istogramma 2D (puoi regolare i range se vuoi)
  TH2D* h2 = new TH2D("h2", "Hit map; z (mm); y (mm)", 200, -8, 8, 200, -8, 8);

  // Riempimento istogramma (equivalente a Draw)
  tree->Draw("m_y:m_z>>h2", "", "colz");

  // Disegna
  h2->Draw("colz");

  c1->Update();
  c1->SaveAs("output/hits_map.png");

  std::cout << "Plot salvato in output/hits_map.png" << std::endl;

  file->Close();
  return 0;
}