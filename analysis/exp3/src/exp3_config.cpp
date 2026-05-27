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
        cfg.display.width_px      = d.value("width_px",      1080);
        cfg.display.height_px     = d.value("height_px",     2340);
        cfg.display.wavelength_nm = d.value("wavelength_nm", 525.0);
    }

    // Distanze assiali
    if (j.contains("axial_distances_nominal_mm"))
        cfg.axial_distances_nominal_mm =
            j["axial_distances_nominal_mm"].get<std::vector<double>>();
    if (j.contains("axial_distances_measured_mm"))
        cfg.axial_distances_measured_mm =
            j["axial_distances_measured_mm"].get<std::vector<double>>();

    if (cfg.axial_distances_measured_mm.empty())
        cfg.axial_distances_measured_mm = cfg.axial_distances_nominal_mm;

    if (cfg.axial_distances_measured_mm.size() != cfg.axial_distances_nominal_mm.size())
        throw std::runtime_error(
            "load_exp3_config: axial_distances_nominal_mm e axial_distances_measured_mm "
            "devono avere lo stesso numero di elementi");

    // Nuovi campi geometrici
    if (j.contains("radial_offsets_px"))
        cfg.radial_offsets_px = j["radial_offsets_px"].get<std::vector<int>>();
    else
        cfg.radial_offsets_px = {0, 390, 780};

    cfg.dot_column_step_px   = j.value("dot_column_step_px",   50);
    cfg.calib_line_length_px = j.value("calib_line_length_px", 800);

    if (j.contains("optical_axis_center_px") &&
        j["optical_axis_center_px"].is_array() &&
        j["optical_axis_center_px"].size() == 2) {
        cfg.optical_axis_center_px[0] = j["optical_axis_center_px"][0].get<double>();
        cfg.optical_axis_center_px[1] = j["optical_axis_center_px"][1].get<double>();
    }

    // Stacking (oggetto annidato o campi flat per retrocompatibilità)
    if (j.contains("stacking")) {
        const auto& s = j["stacking"];
        cfg.stacking.n_sigma    = s.value("n_sigma",    3.0);
        cfg.stacking.n_iter     = s.value("n_iter",     3);
        cfg.stacking.min_frames = s.value("min_frames", 2);
    } else {
        cfg.stacking.n_sigma    = j.value("stack_n_sigma",    3.0);
        cfg.stacking.n_iter     = j.value("stack_n_iter",     3);
        cfg.stacking.min_frames = j.value("stack_min_frames", 2);
    }

    // Estrazione traccia
    if (j.contains("trace_extraction")) {
        const auto& t = j["trace_extraction"];
        cfg.trace_extraction.min_snr          = t.value("min_snr",          5.0);
        cfg.trace_extraction.min_valid_slices = t.value("min_valid_slices", 20);
    } else {
        cfg.trace_extraction.min_snr          = j.value("min_snr",          5.0);
        cfg.trace_extraction.min_valid_slices = j.value("min_valid_slices", 20);
    }

    // Q comparison
    if (j.contains("q_comparison")) {
        const auto& q = j["q_comparison"];
        cfg.q_comparison.q_map_tsv  = q.value("q_map_tsv",  std::string("output/psf_analysis/q_map.tsv"));
        cfg.q_comparison.good_x1_mm = q.value("good_x1_mm", 0.0);
        cfg.q_comparison.good_x2_mm = q.value("good_x2_mm", 0.0);
        cfg.q_comparison.bad_x1_mm  = q.value("bad_x1_mm",  0.0);
        cfg.q_comparison.bad_x2_mm  = q.value("bad_x2_mm",  0.0);
    } else {
        // Retrocompatibilità: q_map_tsv al livello radice
        cfg.q_comparison.q_map_tsv = j.value("q_map_tsv",
                                              std::string("output/psf_analysis/q_tsv11.tsv"));
    }

    return cfg;
}

} // namespace riptide::exp3
