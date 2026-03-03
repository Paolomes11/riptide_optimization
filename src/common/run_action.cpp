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

#include "run_action.hpp"

#include <G4AnalysisManager.hh>

RunAction::RunAction()
    : G4UserRunAction() {
}

RunAction::~RunAction() {
}

void RunAction::BeginOfRunAction(const G4Run*) {
  // G4AnalysisManager* am = G4AnalysisManager::Instance();
  // am->OpenFile("output/output.root");

  // std::cout << "Creating ntuple..." << std::endl;

  // am->CreateNtuple("Hits", "Hits");
  // am->CreateNtupleIColumn("m_event");
  // am->CreateNtupleDColumn("m_x");
  // am->CreateNtupleDColumn("m_y");
  // am->CreateNtupleDColumn("m_z");
  // am->FinishNtuple(0);
}

void RunAction::EndOfRunAction(const G4Run*) {
  // G4AnalysisManager* am = G4AnalysisManager::Instance();
  // std::cout << "Writing and closing file..." << std::endl;
  // am->Write();
  // am->CloseFile();
}
