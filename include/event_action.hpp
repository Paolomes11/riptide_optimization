#ifndef RIPTTIDE_EVENT_ACTION_HPP
#define RIPTTIDE_EVENT_ACTION_HPP

#include <G4UserEventAction.hh>

#include <vector>

namespace riptide
{

    class EventAction : public G4UserEventAction
    {
        // Detected photon hits
        std::vector<int> m_hit_detector;
        std::vector<double> m_hit_x;
        std::vector<double> m_hit_y;
        std::vector<double> m_hit_z;
        std::vector<double> m_hit_time;
        std::vector<double> m_hit_energy;

    public:
        void BeginOfEventAction(const G4Event *event) override;
        void EndOfEventAction(const G4Event *event) override;

        // Getters per accedere ai dati dei fotoni rilevati
        std::vector<int> &hit_detector() { return m_hit_detector; }
        std::vector<double> &hit_x() { return m_hit_x; }
        std::vector<double> &hit_y() { return m_hit_y; }
        std::vector<double> &hit_z() { return m_hit_z; }
        std::vector<double> &hit_time() { return m_hit_time; }
        std::vector<double> &hit_energy() { return m_hit_energy; }
    };

}

#endif // RIPTTIDE_EVENT_ACTION_HPP