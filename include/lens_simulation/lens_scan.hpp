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

#ifndef RIPTIDE_LENS_SCAN_HPP
#define RIPTIDE_LENS_SCAN_HPP

#include <filesystem>
#include <string>

class G4RunManager;

namespace riptide {
void lens_scan(G4RunManager* run_manager, const std::filesystem::path& macro_file,
                      const std::string& root_output_file);
}

#endif // RIPTIDE_LENS_SCAN_HPP