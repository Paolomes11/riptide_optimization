#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <set>

double compute_geometric_efficiency(const std::string& root_file, int n_photons_shot) {
  TFile* file = TFile::Open(root_file.c_str(), "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Errore: impossibile aprire " << root_file << std::endl;
    return 0.0;
  }

  TTree* tree = nullptr;
  file->GetObject("Hits", tree);
  if (!tree) {
    std::cerr << "Errore: TTree 'Hits' non trovato\n";
    file->Close();
    return 0.0;
  }

  int eventID;
  tree->SetBranchAddress("m_event", &eventID);

  std::set<int> events_with_hits;

  Long64_t nentries = tree->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {
    tree->GetEntry(i);
    events_with_hits.insert(eventID); // contiamo ogni evento una sola volta
  }

  file->Close();

  double efficiency = static_cast<double>(events_with_hits.size()) / n_photons_shot;

  std::cout << "Eventi con almeno un hit: " << events_with_hits.size() << " su " << n_photons_shot
            << " → efficienza geometrica = " << efficiency << std::endl;

  return efficiency;
}