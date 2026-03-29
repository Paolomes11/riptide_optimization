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

#ifndef RIPTIDE_STEPPING_ACTION_HPP
#define RIPTIDE_STEPPING_ACTION_HPP

#include <G4UserSteppingAction.hh>

namespace riptide {

class SteppingAction : public G4UserSteppingAction {
 public:
  SteppingAction();
  ~SteppingAction() override = default;

  void UserSteppingAction(const G4Step* step) override;
};

} // namespace riptide

#endif
