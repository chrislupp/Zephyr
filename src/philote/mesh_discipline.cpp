// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include "vlm/philote/mesh_discipline.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace vlm {

// =============================================================================
// Precomputation
// =============================================================================

void MeshDiscipline::PrecomputeChordFracs() {
    chord_fracs_.resize(nx_ + 1);
    for (int64_t i = 0; i <= nx_; ++i) {
        if (cosine_chordwise_spacing_) {
            chord_fracs_[i] =
                0.5 * (1.0 - std::cos(M_PI * static_cast<double>(i) /
                                      static_cast<double>(nx_)));
        } else {
            chord_fracs_[i] = static_cast<double>(i) / static_cast<double>(nx_);
        }
    }
}

// =============================================================================
// Philote lifecycle
// =============================================================================

void MeshDiscipline::Initialize() {
    ExplicitDiscipline::Initialize();

    AddOption("nx", "int");
    AddOption("ny", "int");
    AddOption("cosine_chordwise_spacing", "string");
    AddOption("cosine_spanwise_spacing", "string");
}

void MeshDiscipline::SetOptions(
    const google::protobuf::Struct& options_struct) {
    const auto& fields = options_struct.fields();

    auto it = fields.find("nx");
    if (it != fields.end())
        nx_ = static_cast<int64_t>(it->second.number_value());

    it = fields.find("ny");
    if (it != fields.end())
        ny_ = static_cast<int64_t>(it->second.number_value());

    it = fields.find("cosine_chordwise_spacing");
    if (it != fields.end()) {
        std::string val = it->second.string_value();
        std::transform(val.begin(), val.end(), val.begin(), ::tolower);
        cosine_chordwise_spacing_ = (val == "true" || val == "1");
    }

    it = fields.find("cosine_spanwise_spacing");
    if (it != fields.end()) {
        std::string val = it->second.string_value();
        std::transform(val.begin(), val.end(), val.begin(), ::tolower);
        cosine_spanwise_spacing_ = (val == "true" || val == "1");
    }

    PrecomputeChordFracs();

    ExplicitDiscipline::SetOptions(options_struct);
}

void MeshDiscipline::Setup() {
    int64_t n_stations = ny_ + 1;

    // Per-spanwise-station inputs
    AddInput("chord", {n_stations}, "m");
    AddInput("twist", {n_stations}, "rad");
    AddInput("xle", {n_stations}, "m");
    AddInput("yle", {n_stations}, "m");
    AddInput("zle", {n_stations}, "m");

    // Flat mesh output
    int64_t num_values = (nx_ + 1) * n_stations * 3;
    AddOutput("mesh", {num_values}, "m");
}

void MeshDiscipline::SetupPartials() {
    DeclarePartials("mesh", "chord");
    DeclarePartials("mesh", "twist");
    DeclarePartials("mesh", "xle");
    DeclarePartials("mesh", "yle");
    DeclarePartials("mesh", "zle");
}

// =============================================================================
// Compute
// =============================================================================

void MeshDiscipline::Compute(const philote::Variables& inputs,
                             philote::Variables& outputs) {
    const auto& chord = inputs.at("chord");
    const auto& twist = inputs.at("twist");
    const auto& xle = inputs.at("xle");
    const auto& yle = inputs.at("yle");
    const auto& zle = inputs.at("zle");
    auto& mesh = outputs.at("mesh");

    for (int64_t j = 0; j <= ny_; ++j) {
        double cos_tw = std::cos(twist(j));
        double sin_tw = std::sin(twist(j));

        for (int64_t i = 0; i <= nx_; ++i) {
            double cf = chord_fracs_[i];
            int64_t idx = (i * (ny_ + 1) + j) * 3;

            mesh(idx + 0) = xle(j) + cf * chord(j) * cos_tw;
            mesh(idx + 1) = yle(j);
            mesh(idx + 2) = zle(j) - cf * chord(j) * sin_tw;
        }
    }
}

// =============================================================================
// ComputePartials
// =============================================================================

void MeshDiscipline::ComputePartials(const philote::Variables& inputs,
                                     philote::Partials& partials) {
    const auto& chord = inputs.at("chord");
    const auto& twist = inputs.at("twist");

    auto& dm_dc = partials[{"mesh", "chord"}];
    auto& dm_dt = partials[{"mesh", "twist"}];
    auto& dm_dxle = partials[{"mesh", "xle"}];
    auto& dm_dyle = partials[{"mesh", "yle"}];
    auto& dm_dzle = partials[{"mesh", "zle"}];

    int64_t n_st = ny_ + 1;
    int64_t n_mesh = (nx_ + 1) * n_st * 3;

    // Zero all partials
    for (int64_t k = 0; k < n_mesh * n_st; ++k) {
        dm_dc(k) = 0.0;
        dm_dt(k) = 0.0;
        dm_dxle(k) = 0.0;
        dm_dyle(k) = 0.0;
        dm_dzle(k) = 0.0;
    }

    for (int64_t j = 0; j <= ny_; ++j) {
        double cos_tw = std::cos(twist(j));
        double sin_tw = std::sin(twist(j));
        double c_j = chord(j);

        for (int64_t i = 0; i <= nx_; ++i) {
            double cf = chord_fracs_[i];
            int64_t out_base = (i * n_st + j) * 3;

            // x component (out_base + 0)
            int64_t row_x = (out_base + 0) * n_st;
            dm_dc(row_x + j) = cf * cos_tw;
            dm_dt(row_x + j) = -cf * c_j * sin_tw;
            dm_dxle(row_x + j) = 1.0;

            // y component (out_base + 1)
            int64_t row_y = (out_base + 1) * n_st;
            dm_dyle(row_y + j) = 1.0;

            // z component (out_base + 2)
            int64_t row_z = (out_base + 2) * n_st;
            dm_dc(row_z + j) = -cf * sin_tw;
            dm_dt(row_z + j) = -cf * c_j * cos_tw;
            dm_dzle(row_z + j) = 1.0;
        }
    }
}

}  // namespace vlm
