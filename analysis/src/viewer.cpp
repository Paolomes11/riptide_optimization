#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>

#include <iostream>
#include <vector>

int main(int argc, char **argv)
{
    std::string filename = "output/riptide_000.root";
    if (argc > 1)
        filename = argv[1];

    std::cout << "Apro file: " << filename << std::endl;

    TFile file(filename.c_str());
    if (file.IsZombie())
    {
        std::cerr << "Errore apertura file\n";
        return 1;
    }

    file.ls(); // DEBUG: mostra contenuto file

    auto *tree = static_cast<TTree *>(file.Get("photon_hits"));
    if (!tree)
    {
        std::cerr << "Ntuple 'photon_hits' non trovato\n";
        return 1;
    }

    std::cout << "Tree trovato\n";
    tree->Print(); // DEBUG: struttura ntuple

    std::vector<double> *hit_x = nullptr;
    std::vector<double> *hit_y = nullptr;

    tree->SetBranchAddress("hit_x", &hit_x);
    tree->SetBranchAddress("hit_y", &hit_y);

    TH2D hist("hits", "Photon hits;X;Y", 200, -100, 100, 200, -100, 100);

    Long64_t nentries = tree->GetEntries();
    std::cout << "Numero eventi nel tree: " << nentries << std::endl;

    size_t total_hits = 0;

    for (Long64_t i = 0; i < nentries; i++)
    {
        tree->GetEntry(i);

        if (!hit_x || !hit_y)
        {
            std::cout << "Evento " << i << " vettori null\n";
            continue;
        }

        if (i < 5) // DEBUG: stampa primi eventi
        {
            std::cout << "Evento " << i
                      << " hits: " << hit_x->size() << std::endl;
        }

        for (size_t j = 0; j < hit_x->size(); j++)
        {
            hist.Fill(hit_x->at(j), hit_y->at(j));
            total_hits++;
        }
    }

    std::cout << "Hit totali riempiti nell'istogramma: "
              << total_hits << std::endl;

    std::cout << "Entries histogram: " << hist.GetEntries() << std::endl;

    TCanvas c;
    hist.Draw("COLZ");
    c.SaveAs("photon_hits.png");

    std::cout << "Immagine salvata: photon_hits.png\n";
}