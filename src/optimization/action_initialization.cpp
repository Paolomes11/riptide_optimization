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

#include "action_initialization.hpp"
#include "event_action.hpp"
#include "primary_generator_action.hpp"
#include "run_action.hpp"
#include <iostream>

namespace riptide {

ActionInitialization::ActionInitialization(const std::string& output_file)
    : m_output_file(output_file) {
  std::cout << "ActionInitialization: Constructor called, output file: " << output_file
            << std::endl;
}

void ActionInitialization::BuildForMaster() const {
  std::cout << "ActionInitialization: BuildForMaster called" << std::endl;
}

void ActionInitialization::Build() const {
  // Azioni standard
  SetUserAction(new PrimaryGeneratorAction());
  SetUserAction(new RunAction());

  // EventAction con file ROOT
  SetUserAction(new EventAction(m_output_file));
}

} // namespace riptide
