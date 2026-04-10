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

#ifndef RIPTIDE_ACTION_INITIALIZATION_HPP
#define RIPTIDE_ACTION_INITIALIZATION_HPP

#include <G4VUserActionInitialization.hh>
#include <string>

namespace riptide {

class ActionInitialization : public G4VUserActionInitialization {
  bool m_useImportanceSampling = false;

 public:
  ActionInitialization(bool useImportanceSampling = false)
      : m_useImportanceSampling(useImportanceSampling) {
  }
  virtual ~ActionInitialization() = default;

  void BuildForMaster() const override;
  void Build() const override;
};

} // namespace riptide

#endif // RIPTIDE_ACTION_INITIALIZATION_HPP