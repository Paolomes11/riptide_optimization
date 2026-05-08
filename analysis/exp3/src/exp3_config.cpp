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

#include "exp3_config.hpp"

#include <nlohmann/json.hpp>

#include <fstream>
#include <stdexcept>

namespace riptide::exp3 {

Exp3Config load_exp3_config(const std::filesystem::path& path) {
    std::ifstream ifs(path);
    if (!ifs)
        throw std::runtime_error("load_exp3_config: file non trovato: " + path.string());

    nlohmann::json j;
    try {
        ifs >> j;
    } catch (const nlohmann::json::parse_error& e) {
        throw std::runtime_error(std::string("load_exp3_config: JSON malformato: ") + e.what());
    }

    Exp3Config cfg;

    // Display
    if (j.contains("display")) {
        const auto& d = j["display"];
        cfg.display.mm_per_px_x   = d.value("mm_per_px_x",   0.0);
        cfg.display.mm_per_px_y   = d.value("mm_per_px_y",   0.0);
        cfg.display.width_px      = d.value("width_px",      1920);
        cfg.display.height_px     = d.value("height_px",     1080);
        cfg.display.wavelength_nm = d.value("wavelength_nm", 525.0);
    }

    // Distanze assiali
    if (j.contains("axial_distances_nominal_mm"))
        cfg.axial_distances_nominal_mm =
            j["axial_distances_nominal_mm"].get<std::vector<double>>();
    if (j.contains("axial_distances_measured_mm"))
        cfg.axial_distances_measured_mm =
            j["axial_distances_measured_mm"].get<std::vector<double>>();

    // Se measured non specificato, copia da nominal
    if (cfg.axial_distances_measured_mm.empty())
        cfg.axial_distances_measured_mm = cfg.axial_distances_nominal_mm;

    // Assicura stesso numero di elementi
    if (cfg.axial_distances_measured_mm.size() != cfg.axial_distances_nominal_mm.size())
        throw std::runtime_error(
            "load_exp3_config: axial_distances_nominal_mm e axial_distances_measured_mm "
            "devono avere lo stesso numero di elementi");

    // Orientazioni
    if (j.contains("orientations_deg"))
        cfg.orientations_deg = j["orientations_deg"].get<std::vector<double>>();

    // Calibrazione
    cfg.calibration_grid_step_px = j.value("calibration_grid_step_px", 50);

    // Stacking
    cfg.stack_n_sigma    = j.value("stack_n_sigma",    3.0);
    cfg.stack_n_iter     = j.value("stack_n_iter",     3);
    cfg.stack_min_frames = j.value("stack_min_frames", 2);

    // Estrazione traccia
    cfg.min_snr             = j.value("min_snr",             5.0);
    cfg.min_valid_slices    = j.value("min_valid_slices",    30);
    cfg.angle_tolerance_deg = j.value("angle_tolerance_deg", 10.0);

    // Q map simulazione
    cfg.q_map_tsv = j.value("q_map_tsv", std::string("output/psf_analysis/q_tsv11.tsv"));

    // Configurazioni lenti
    if (j.contains("lens_configs")) {
        for (const auto& lj : j["lens_configs"]) {
            LensConfig lc;
            lc.x1_mm = lj.value("x1_mm", 0.0);
            lc.x2_mm = lj.value("x2_mm", 0.0);
            lc.label  = lj.value("label", std::string(""));
            cfg.lens_configs.push_back(lc);
        }
    }

    return cfg;
}

} // namespace riptide::exp3
