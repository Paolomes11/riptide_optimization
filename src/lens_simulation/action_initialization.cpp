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
#include "stepping_action.hpp"

namespace riptide {

void ActionInitialization::BuildForMaster() const {
  SetUserAction(new RunAction());
}

void ActionInitialization::Build() const {
  // Azioni standard
  SetUserAction(new PrimaryGeneratorAction());
  SetUserAction(new RunAction());

  // SteppingAction per ottimizzazione uccisione fotoni
  SetUserAction(new SteppingAction());

  // EventAction con file ROOT
  SetUserAction(new EventAction());
}

} // namespace riptide
