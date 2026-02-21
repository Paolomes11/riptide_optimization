#ifndef RIPTIDE_RUN_ACTION_HPP
#define RIPTIDE_RUN_ACTION_HPP

#include <G4UserRunAction.hh>

#include <vector>

namespace riptide
{
    // Forward declaration
    class EventAction;

    class RunAction : public G4UserRunAction
    {
        EventAction *m_event_action{nullptr};

        // Dummy vectors
        std::vector<int> m_dummy_hit_detector;
        std::vector<double> m_dummy_hit_x;
        std::vector<double> m_dummy_hit_y;
        std::vector<double> m_dummy_hit_z;
        std::vector<double> m_dummy_hit_time;
        std::vector<double> m_dummy_hit_energy;

        static std::string s_filename;

    public:
        explicit RunAction(EventAction *ea = nullptr);
        void BeginOfRunAction(G4Run const *run) override;
        void EndOfRunAction(G4Run const *run) override;
    };

} // namespace riptide

#endif // RIPTIDE_RUN_ACTION_HPP