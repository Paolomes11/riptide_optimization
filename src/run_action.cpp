#include "run_action.hpp"
#include "event_action.hpp"

#include <G4AnalysisManager.hh>
#include <G4Run.hh>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <filesystem>

namespace
{

    std::string generate_unique_filename()
    {
        for (int i = 0; i < 999; ++i)
        {
            std::string filename = fmt::format("riptide_{:03}.root", i);
            if (!std::filesystem::exists(filename))
            {
                return filename;
            }
        }
        throw std::runtime_error("Maximum number of output files reached (999)");
    }

} // namespace

namespace riptide
{

    std::string RunAction::s_filename;

    RunAction::RunAction(EventAction *ea)
        : m_event_action{ea}
    {
    }

    void RunAction::BeginOfRunAction(G4Run const *)
    {
        auto *am = G4AnalysisManager::Instance();
        am->SetDefaultFileType("root");
        am->SetVerboseLevel(1);
        am->SetNtupleMerging(true);

        if (m_event_action == nullptr)
        {
            s_filename = generate_unique_filename();
        }

        auto stem = std::filesystem::path(s_filename).stem().string();
        am->OpenFile(stem);

        am->CreateNtuple("photon_hits", "Optical photon hits on detector");
        am->CreateNtupleIColumn("event_id");
        am->CreateNtupleIColumn("n_hits");

        auto &det = m_event_action ? m_event_action->hit_detector() : m_dummy_hit_detector;
        auto &hx = m_event_action ? m_event_action->hit_x() : m_dummy_hit_x;
        auto &hy = m_event_action ? m_event_action->hit_y() : m_dummy_hit_y;
        auto &hz = m_event_action ? m_event_action->hit_z() : m_dummy_hit_z;
        auto &ht = m_event_action ? m_event_action->hit_time() : m_dummy_hit_time;
        auto &he = m_event_action ? m_event_action->hit_energy() : m_dummy_hit_energy;

        am->CreateNtupleIColumn("hit_detector", det);
        am->CreateNtupleDColumn("hit_x", hx);
        am->CreateNtupleDColumn("hit_y", hy);
        am->CreateNtupleDColumn("hit_z", hz);
        am->CreateNtupleDColumn("hit_time", ht);
        am->CreateNtupleDColumn("hit_energy", he);

        am->FinishNtuple();

        spdlog::info("ROOT output file opened: {}", s_filename);
    }

    void RunAction::EndOfRunAction(const G4Run *)
    {
        auto *am = G4AnalysisManager::Instance();
        am->Write();
        am->CloseFile();

        spdlog::info("ROOT output file closed: {}", s_filename);
    }

} // namespace riptide