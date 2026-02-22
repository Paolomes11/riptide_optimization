#include "run_action.hpp"

#include <G4AnalysisManager.hh>

RunAction::RunAction() : G4UserRunAction()
{
}

RunAction::~RunAction()
{
}

void RunAction::BeginOfRunAction(const G4Run *)
{
    G4AnalysisManager *am = G4AnalysisManager::Instance();
    am->OpenFile("output.root");

    std::cout << "Creating ntuple..." << std::endl;

    am->CreateNtuple("Hits", "Hits");
    am->CreateNtupleIColumn("m_event");
    am->CreateNtupleDColumn("m_x");
    am->CreateNtupleDColumn("m_y");
    am->CreateNtupleDColumn("m_z");
    am->FinishNtuple(0);
}

void RunAction::EndOfRunAction(const G4Run *)
{
    G4AnalysisManager *am = G4AnalysisManager::Instance();
    std::cout << "Writing and closing file..." << std::endl;
    am->Write();
    am->CloseFile();
}