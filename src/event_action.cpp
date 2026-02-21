#include "event_action.hpp"
#include "photocathode_hit.hpp"

#include <G4HCofThisEvent.hh>
#include <G4Event.hh>
#include <G4SystemOfUnits.hh>
#include <G4AnalysisManager.hh>

namespace riptide
{
    void EventAction::BeginOfEventAction(G4Event const *)
    {
        // Pulisco i vettori per il nuovo evento
        m_hit_detector.clear();
        m_hit_x.clear();
        m_hit_y.clear();
        m_hit_z.clear();
        m_hit_time.clear();
        m_hit_energy.clear();
    }

    void EventAction::EndOfEventAction(G4Event const *event)
    {
        // Prende le hits dalle collezioni di hits
        // hce = hits collection of this event
        auto *hce = event->GetHCofThisEvent();
        if (hce)
        {
            int n_collections = hce->GetNumberOfCollections();
            for (int i = 0; i < n_collections; ++i)
            {
                auto *hc = dynamic_cast<PhotocathodeHitsCollection *>(hce->GetHC(i));
                if (!hc)
                    continue;

                auto n_hits = hc->GetSize();
                for (std::size_t j = 0; j < n_hits; ++j)
                {
                    auto *hit = static_cast<PhotocathodeHit *>(hc->GetHit(j));
                    m_hit_detector.push_back(hit->detector_id());
                    m_hit_x.push_back(hit->global_position().x() / mm);
                    m_hit_y.push_back(hit->global_position().y() / mm);
                    m_hit_z.push_back(hit->global_position().z() / mm);
                    m_hit_time.push_back(hit->time() / ns);
                    m_hit_energy.push_back(hit->energy() / eV);
                }
            }
        }

        // Debug: mostra cosa verrà scritto nel file ROOT
        std::cout << "=== DEBUG: Dati che verranno scritti nel ROOT ===\n";
        std::cout << "Evento ID: " << event->GetEventID() << "\n";
        std::cout << "Numero di hits: " << m_hit_detector.size() << "\n";
        
        if (!m_hit_detector.empty()) {
            std::cout << "Dettagli hits:\n";
            for (size_t i = 0; i < m_hit_detector.size(); ++i) {
                std::cout << "  Hit " << i << ": detector=" << m_hit_detector[i] 
                          << ", pos=(" << m_hit_x[i] << ", " << m_hit_y[i] << ", " << m_hit_z[i] << ") mm"
                          << ", tempo=" << m_hit_time[i] << " ns"
                          << ", energia=" << m_hit_energy[i] << " eV\n";
            }
        } else {
            std::cout << "  Nessun hit rilevato in questo evento\n";
        }
        std::cout << "===============================================\n\n";

        // Scrivi i dati dei fotoni rilevati nell'ntuple
        auto *am = G4AnalysisManager::Instance();
        am->FillNtupleIColumn(0, 0, event->GetEventID());
        am->FillNtupleIColumn(0, 1, static_cast<int>(m_hit_detector.size()));
        am->AddNtupleRow();
    }
}