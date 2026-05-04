/*
    Project Pegasus
    Version 0.2.0

    Copyright 2024-2026, Christopher A. Lupp

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "pegasus/philote/mesh_discipline.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace pegasus {

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
    AddOption("y_duplicated", "string");
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

    it = fields.find("y_duplicated");
    if (it != fields.end()) {
        std::string val = it->second.string_value();
        std::transform(val.begin(), val.end(), val.begin(), ::tolower);
        y_duplicated_ = (val == "true" || val == "1");
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
    AddInput("zle", {n_stations}, "m");

    // Span parameterization
    AddInput("halfspan", {1}, "m");
    AddInput("eta", {n_stations}, "");

    // Flat mesh output
    int64_t num_values = (nx_ + 1) * n_stations * 3;
    AddOutput("mesh", {num_values}, "m");

    // Reference quantities
    AddOutput("reference_area", {1}, "m**2");
    AddOutput("reference_span", {1}, "m");
    AddOutput("reference_chord", {1}, "m");
}

void MeshDiscipline::SetupPartials() {
    DeclarePartials("mesh", "chord");
    DeclarePartials("mesh", "twist");
    DeclarePartials("mesh", "xle");
    DeclarePartials("mesh", "zle");
    DeclarePartials("mesh", "halfspan");
    DeclarePartials("mesh", "eta");

    // Reference area partials
    DeclarePartials("reference_area", "chord");
    DeclarePartials("reference_area", "zle");
    DeclarePartials("reference_area", "halfspan");
    DeclarePartials("reference_area", "eta");

    // Reference span partials
    DeclarePartials("reference_span", "halfspan");
    DeclarePartials("reference_span", "eta");

    // Reference chord (MAC) partials
    DeclarePartials("reference_chord", "chord");
    DeclarePartials("reference_chord", "zle");
    DeclarePartials("reference_chord", "halfspan");
    DeclarePartials("reference_chord", "eta");
}

// =============================================================================
// Compute
// =============================================================================

void MeshDiscipline::Compute(const philote::Variables& inputs,
                             philote::Variables& outputs) {
    const auto& chord = inputs.at("chord");
    const auto& twist = inputs.at("twist");
    const auto& xle = inputs.at("xle");
    const auto& zle = inputs.at("zle");
    const auto& eta = inputs.at("eta");
    double halfspan = inputs.at("halfspan")(0);
    auto& mesh = outputs.at("mesh");

    for (int64_t j = 0; j <= ny_; ++j) {
        double cos_tw = std::cos(twist(j));
        double sin_tw = std::sin(twist(j));
        double y_j = eta(j) * halfspan;

        for (int64_t i = 0; i <= nx_; ++i) {
            double cf = chord_fracs_[i];
            int64_t idx = (i * (ny_ + 1) + j) * 3;

            mesh(idx + 0) = xle(j) + cf * chord(j) * cos_tw;
            mesh(idx + 1) = y_j;
            mesh(idx + 2) = zle(j) - cf * chord(j) * sin_tw;
        }
    }

    // Reference area (planform area via trapezoidal rule)
    double sym = y_duplicated_ ? 2.0 : 1.0;
    double s_ref_half = 0.0;
    double mac_num = 0.0;

    for (int64_t j = 0; j < ny_; ++j) {
        double dy = (eta(j + 1) - eta(j)) * halfspan;
        double dz = zle(j + 1) - zle(j);
        double ds = std::sqrt(dy * dy + dz * dz);
        double avg_c = 0.5 * (chord(j) + chord(j + 1));
        s_ref_half += avg_c * ds;
        mac_num +=
            0.5 * (chord(j) * chord(j) + chord(j + 1) * chord(j + 1)) * ds;
    }

    outputs.at("reference_area")(0) = sym * s_ref_half;
    outputs.at("reference_span")(0) = sym * (eta(ny_) - eta(0)) * halfspan;
    outputs.at("reference_chord")(0) =
        (s_ref_half > 0.0) ? (mac_num / s_ref_half) : 0.0;
}

// =============================================================================
// ComputePartials
// =============================================================================

void MeshDiscipline::ComputePartials(const philote::Variables& inputs,
                                     philote::Partials& partials) {
    const auto& chord = inputs.at("chord");
    const auto& twist = inputs.at("twist");
    const auto& eta = inputs.at("eta");
    const auto& zle = inputs.at("zle");
    double halfspan = inputs.at("halfspan")(0);

    auto& dm_dc = partials[{"mesh", "chord"}];
    auto& dm_dt = partials[{"mesh", "twist"}];
    auto& dm_dxle = partials[{"mesh", "xle"}];
    auto& dm_dzle = partials[{"mesh", "zle"}];
    auto& dm_deta = partials[{"mesh", "eta"}];
    auto& dm_dhs = partials[{"mesh", "halfspan"}];

    int64_t n_st = ny_ + 1;
    int64_t n_mesh = (nx_ + 1) * n_st * 3;

    // Zero all mesh partials (n_mesh x n_st arrays)
    for (int64_t k = 0; k < n_mesh * n_st; ++k) {
        dm_dc(k) = 0.0;
        dm_dt(k) = 0.0;
        dm_dxle(k) = 0.0;
        dm_dzle(k) = 0.0;
        dm_deta(k) = 0.0;
    }
    // halfspan partial is n_mesh x 1
    for (int64_t k = 0; k < n_mesh; ++k) {
        dm_dhs(k) = 0.0;
    }

    for (int64_t j = 0; j <= ny_; ++j) {
        double cos_tw = std::cos(twist(j));
        double sin_tw = std::sin(twist(j));
        double c_j = chord(j);
        double eta_j = eta(j);

        for (int64_t i = 0; i <= nx_; ++i) {
            double cf = chord_fracs_[i];
            int64_t out_base = (i * n_st + j) * 3;

            // x component (out_base + 0)
            int64_t row_x = (out_base + 0) * n_st;
            dm_dc(row_x + j) = cf * cos_tw;
            dm_dt(row_x + j) = -cf * c_j * sin_tw;
            dm_dxle(row_x + j) = 1.0;

            // y component (out_base + 1): mesh_y = eta(j) * halfspan
            int64_t row_y = (out_base + 1) * n_st;
            dm_deta(row_y + j) = halfspan;
            dm_dhs(out_base + 1) = eta_j;

            // z component (out_base + 2)
            int64_t row_z = (out_base + 2) * n_st;
            dm_dc(row_z + j) = -cf * sin_tw;
            dm_dt(row_z + j) = -cf * c_j * cos_tw;
            dm_dzle(row_z + j) = 1.0;
        }
    }

    // =========================================================================
    // Reference quantity partials
    // =========================================================================
    double sym = y_duplicated_ ? 2.0 : 1.0;

    // Precompute per-segment values
    std::vector<double> ds(ny_), avg_c(ny_), deta_seg(ny_), dz_seg(ny_);
    double s_ref_half = 0.0;
    double mac_num = 0.0;

    for (int64_t j = 0; j < ny_; ++j) {
        deta_seg[j] = eta(j + 1) - eta(j);
        dz_seg[j] = zle(j + 1) - zle(j);
        double dy = deta_seg[j] * halfspan;
        ds[j] = std::sqrt(dy * dy + dz_seg[j] * dz_seg[j]);
        avg_c[j] = 0.5 * (chord(j) + chord(j + 1));
        s_ref_half += avg_c[j] * ds[j];
        mac_num +=
            0.5 * (chord(j) * chord(j) + chord(j + 1) * chord(j + 1)) * ds[j];
    }

    // --- reference_area partials ---
    auto& dSref_dc = partials[{"reference_area", "chord"}];
    auto& dSref_dz = partials[{"reference_area", "zle"}];
    auto& dSref_dhs = partials[{"reference_area", "halfspan"}];
    auto& dSref_deta = partials[{"reference_area", "eta"}];

    for (int64_t j = 0; j < n_st; ++j) {
        dSref_dc(j) = 0.0;
        dSref_dz(j) = 0.0;
        dSref_deta(j) = 0.0;
    }
    dSref_dhs(0) = 0.0;

    for (int64_t j = 0; j < ny_; ++j) {
        // dS_ref/dchord[j] and dS_ref/dchord[j+1]
        dSref_dc(j) += sym * 0.5 * ds[j];
        dSref_dc(j + 1) += sym * 0.5 * ds[j];

        if (ds[j] > 0.0) {
            // d(ds)/d(halfspan) = deta^2 * halfspan / ds
            double dds_dhs = deta_seg[j] * deta_seg[j] * halfspan / ds[j];
            dSref_dhs(0) += sym * avg_c[j] * dds_dhs;

            // d(ds)/d(eta[j]) = -deta * halfspan^2 / ds
            // d(ds)/d(eta[j+1]) = +deta * halfspan^2 / ds
            double eta_factor = deta_seg[j] * halfspan * halfspan / ds[j];
            dSref_deta(j) += sym * avg_c[j] * (-eta_factor);
            dSref_deta(j + 1) += sym * avg_c[j] * eta_factor;

            // d(ds)/d(zle)
            double dds_dz_j = -dz_seg[j] / ds[j];
            double dds_dz_j1 = dz_seg[j] / ds[j];
            dSref_dz(j) += sym * avg_c[j] * dds_dz_j;
            dSref_dz(j + 1) += sym * avg_c[j] * dds_dz_j1;
        }
    }

    // --- reference_span partials ---
    // ref_span = sym * (eta[ny] - eta[0]) * halfspan
    auto& dspan_dhs = partials[{"reference_span", "halfspan"}];
    auto& dspan_deta = partials[{"reference_span", "eta"}];

    dspan_dhs(0) = sym * (eta(ny_) - eta(0));
    for (int64_t j = 0; j < n_st; ++j) {
        dspan_deta(j) = 0.0;
    }
    dspan_deta(0) = -sym * halfspan;
    dspan_deta(ny_) += sym * halfspan;

    // --- reference_chord (MAC) partials via quotient rule ---
    // MAC = mac_num / s_ref_half
    auto& dMAC_dc = partials[{"reference_chord", "chord"}];
    auto& dMAC_dz = partials[{"reference_chord", "zle"}];
    auto& dMAC_dhs = partials[{"reference_chord", "halfspan"}];
    auto& dMAC_deta = partials[{"reference_chord", "eta"}];

    for (int64_t j = 0; j < n_st; ++j) {
        dMAC_dc(j) = 0.0;
        dMAC_dz(j) = 0.0;
        dMAC_deta(j) = 0.0;
    }
    dMAC_dhs(0) = 0.0;

    if (s_ref_half > 0.0) {
        double denom2 = s_ref_half * s_ref_half;

        for (int64_t j = 0; j < ny_; ++j) {
            // d_mac_num/dchord[j] and [j+1]
            double dnum_dc_j = chord(j) * ds[j];
            double dnum_dc_j1 = chord(j + 1) * ds[j];
            // d_s_ref_half/dchord[j] and [j+1]
            double dden_dc_j = 0.5 * ds[j];
            double dden_dc_j1 = 0.5 * ds[j];

            dMAC_dc(j) +=
                (dnum_dc_j * s_ref_half - mac_num * dden_dc_j) / denom2;
            dMAC_dc(j + 1) +=
                (dnum_dc_j1 * s_ref_half - mac_num * dden_dc_j1) / denom2;

            if (ds[j] > 0.0) {
                double num_coeff =
                    0.5 * (chord(j) * chord(j) + chord(j + 1) * chord(j + 1));
                double den_coeff = avg_c[j];

                // Partials through ds w.r.t. halfspan
                double dds_dhs = deta_seg[j] * deta_seg[j] * halfspan / ds[j];
                dMAC_dhs(0) += (num_coeff * dds_dhs * s_ref_half -
                                mac_num * den_coeff * dds_dhs) /
                               denom2;

                // Partials through ds w.r.t. eta
                double eta_factor = deta_seg[j] * halfspan * halfspan / ds[j];

                double dds_deta_j = -eta_factor;
                double dds_deta_j1 = eta_factor;

                dMAC_deta(j) += (num_coeff * dds_deta_j * s_ref_half -
                                 mac_num * den_coeff * dds_deta_j) /
                                denom2;
                dMAC_deta(j + 1) += (num_coeff * dds_deta_j1 * s_ref_half -
                                     mac_num * den_coeff * dds_deta_j1) /
                                    denom2;

                // Partials through ds w.r.t. zle
                double dds_dz_j = -dz_seg[j] / ds[j];
                double dds_dz_j1 = dz_seg[j] / ds[j];

                dMAC_dz(j) += (num_coeff * dds_dz_j * s_ref_half -
                               mac_num * den_coeff * dds_dz_j) /
                              denom2;
                dMAC_dz(j + 1) += (num_coeff * dds_dz_j1 * s_ref_half -
                                   mac_num * den_coeff * dds_dz_j1) /
                                  denom2;
            }
        }
    }
}

}  // namespace pegasus
