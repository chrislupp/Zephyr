/*
    Project Pegasus
    Version 0.1.0

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

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

namespace pegasus {
namespace {

// Helper: create a protobuf Struct with the given options
google::protobuf::Struct MakeOptions(int nx, int ny, bool y_duplicated = true) {
    google::protobuf::Struct opts;
    auto& fields = *opts.mutable_fields();
    fields["nx"].set_number_value(static_cast<double>(nx));
    fields["ny"].set_number_value(static_cast<double>(ny));
    fields["cosine_chordwise_spacing"].set_string_value("false");
    fields["cosine_spanwise_spacing"].set_string_value("false");
    fields["y_duplicated"].set_string_value(y_duplicated ? "true" : "false");
    return opts;
}

// =============================================================================
// Reference quantity values for a rectangular wing
// =============================================================================

TEST(MeshDisciplineRefTest, RectangularWingArea) {
    // Rectangular wing: chord=1, half-span=5, y-duplicated
    // Expected: S_ref = 2 * chord * half_span = 10.0
    const int nx = 2;
    const int ny = 4;
    const int n_st = ny + 1;
    const double chord_val = 1.0;
    const double half_span = 5.0;

    MeshDiscipline disc;
    disc.Initialize();
    disc.SetOptions(MakeOptions(nx, ny, true));
    disc.Setup();
    disc.SetupPartials();

    philote::Variables inputs;
    inputs["chord"].resize(n_st);
    inputs["twist"].resize(n_st);
    inputs["xle"].resize(n_st);
    inputs["yle"].resize(n_st);
    inputs["zle"].resize(n_st);

    for (int j = 0; j < n_st; ++j) {
        inputs["chord"](j) = chord_val;
        inputs["twist"](j) = 0.0;
        inputs["xle"](j) = 0.0;
        inputs["yle"](j) = half_span * static_cast<double>(j) / ny;
        inputs["zle"](j) = 0.0;
    }

    philote::Variables outputs;
    int64_t n_mesh = (nx + 1) * n_st * 3;
    outputs["mesh"].resize(n_mesh);
    outputs["reference_area"].resize(1);
    outputs["reference_span"].resize(1);
    outputs["reference_chord"].resize(1);

    disc.Compute(inputs, outputs);

    EXPECT_NEAR(outputs["reference_area"](0), 10.0, 1e-12);
    EXPECT_NEAR(outputs["reference_span"](0), 10.0, 1e-12);
    // MAC of a rectangular wing = chord
    EXPECT_NEAR(outputs["reference_chord"](0), 1.0, 1e-12);
}

TEST(MeshDisciplineRefTest, TaperedWingArea) {
    // Linearly tapered wing: root chord=2, tip chord=1, half-span=4
    // S_ref_half = 0.5*(2+1)*4 = 6, S_ref = 12
    // MAC = integral(c^2 dy) / integral(c dy)
    //     = sum 0.5*(c_j^2 + c_{j+1}^2)*ds / S_ref_half
    const int nx = 1;
    const int ny = 4;
    const int n_st = ny + 1;
    const double half_span = 4.0;

    MeshDiscipline disc;
    disc.Initialize();
    disc.SetOptions(MakeOptions(nx, ny, true));
    disc.Setup();
    disc.SetupPartials();

    philote::Variables inputs;
    inputs["chord"].resize(n_st);
    inputs["twist"].resize(n_st);
    inputs["xle"].resize(n_st);
    inputs["yle"].resize(n_st);
    inputs["zle"].resize(n_st);

    for (int j = 0; j < n_st; ++j) {
        double eta = static_cast<double>(j) / ny;
        inputs["chord"](j) = 2.0 - 1.0 * eta;  // linear taper
        inputs["twist"](j) = 0.0;
        inputs["xle"](j) = 0.0;
        inputs["yle"](j) = half_span * eta;
        inputs["zle"](j) = 0.0;
    }

    philote::Variables outputs;
    outputs["mesh"].resize((nx + 1) * n_st * 3);
    outputs["reference_area"].resize(1);
    outputs["reference_span"].resize(1);
    outputs["reference_chord"].resize(1);

    disc.Compute(inputs, outputs);

    EXPECT_NEAR(outputs["reference_area"](0), 12.0, 1e-12);
    EXPECT_NEAR(outputs["reference_span"](0), 8.0, 1e-12);

    // Compute expected MAC numerically
    double mac_num = 0.0, mac_den = 0.0;
    for (int j = 0; j < ny; ++j) {
        double c0 = inputs["chord"](j);
        double c1 = inputs["chord"](j + 1);
        double ds = half_span / ny;
        mac_num += 0.5 * (c0 * c0 + c1 * c1) * ds;
        mac_den += 0.5 * (c0 + c1) * ds;
    }
    double expected_mac = mac_num / mac_den;
    EXPECT_NEAR(outputs["reference_chord"](0), expected_mac, 1e-12);
}

// =============================================================================
// Finite-difference verification of partials
// =============================================================================

TEST(MeshDisciplineRefTest, PartialsFiniteDifference) {
    const int nx = 1;
    const int ny = 3;
    const int n_st = ny + 1;
    const double half_span = 5.0;

    MeshDiscipline disc;
    disc.Initialize();
    disc.SetOptions(MakeOptions(nx, ny, true));
    disc.Setup();
    disc.SetupPartials();

    // Non-uniform chord for interesting partials
    philote::Variables inputs;
    inputs["chord"].resize(n_st);
    inputs["twist"].resize(n_st);
    inputs["xle"].resize(n_st);
    inputs["yle"].resize(n_st);
    inputs["zle"].resize(n_st);

    for (int j = 0; j < n_st; ++j) {
        double eta = static_cast<double>(j) / ny;
        inputs["chord"](j) = 2.0 - 0.8 * eta;
        inputs["twist"](j) = 0.0;
        inputs["xle"](j) = 0.0;
        inputs["yle"](j) = half_span * eta;
        inputs["zle"](j) = 0.1 * eta * eta;  // slight dihedral
    }

    // Compute analytical partials
    philote::Partials partials;
    partials[{"reference_area", "chord"}].resize(n_st);
    partials[{"reference_area", "yle"}].resize(n_st);
    partials[{"reference_area", "zle"}].resize(n_st);
    partials[{"reference_span", "yle"}].resize(n_st);
    partials[{"reference_chord", "chord"}].resize(n_st);
    partials[{"reference_chord", "yle"}].resize(n_st);
    partials[{"reference_chord", "zle"}].resize(n_st);

    // Also need mesh partials (they get zeroed inside)
    int64_t n_mesh = (nx + 1) * n_st * 3;
    partials[{"mesh", "chord"}].resize(n_mesh * n_st);
    partials[{"mesh", "twist"}].resize(n_mesh * n_st);
    partials[{"mesh", "xle"}].resize(n_mesh * n_st);
    partials[{"mesh", "yle"}].resize(n_mesh * n_st);
    partials[{"mesh", "zle"}].resize(n_mesh * n_st);

    disc.ComputePartials(inputs, partials);

    // Helper lambda to compute outputs for perturbed inputs
    auto compute_refs = [&](philote::Variables& inp) {
        philote::Variables out;
        out["mesh"].resize(n_mesh);
        out["reference_area"].resize(1);
        out["reference_span"].resize(1);
        out["reference_chord"].resize(1);
        disc.Compute(inp, out);
        return out;
    };

    const double eps = 1e-7;

    // Check dS_ref/dchord, dMAC/dchord
    for (int j = 0; j < n_st; ++j) {
        auto inp_p = inputs;
        inp_p["chord"](j) += eps;
        auto out_p = compute_refs(inp_p);

        auto inp_m = inputs;
        inp_m["chord"](j) -= eps;
        auto out_m = compute_refs(inp_m);

        double fd_area =
            (out_p["reference_area"](0) - out_m["reference_area"](0)) /
            (2.0 * eps);
        double fd_mac =
            (out_p["reference_chord"](0) - out_m["reference_chord"](0)) /
            (2.0 * eps);

        EXPECT_NEAR(partials[{"reference_area", "chord"}](j), fd_area, 1e-5)
            << "dS_ref/dchord[" << j << "]";
        EXPECT_NEAR(partials[{"reference_chord", "chord"}](j), fd_mac, 1e-5)
            << "dMAC/dchord[" << j << "]";
    }

    // Check dS_ref/dyle, dspan/dyle, dMAC/dyle
    for (int j = 0; j < n_st; ++j) {
        auto inp_p = inputs;
        inp_p["yle"](j) += eps;
        auto out_p = compute_refs(inp_p);

        auto inp_m = inputs;
        inp_m["yle"](j) -= eps;
        auto out_m = compute_refs(inp_m);

        double fd_area =
            (out_p["reference_area"](0) - out_m["reference_area"](0)) /
            (2.0 * eps);
        double fd_span =
            (out_p["reference_span"](0) - out_m["reference_span"](0)) /
            (2.0 * eps);
        double fd_mac =
            (out_p["reference_chord"](0) - out_m["reference_chord"](0)) /
            (2.0 * eps);

        EXPECT_NEAR(partials[{"reference_area", "yle"}](j), fd_area, 1e-5)
            << "dS_ref/dyle[" << j << "]";
        EXPECT_NEAR(partials[{"reference_span", "yle"}](j), fd_span, 1e-5)
            << "dspan/dyle[" << j << "]";
        EXPECT_NEAR(partials[{"reference_chord", "yle"}](j), fd_mac, 1e-5)
            << "dMAC/dyle[" << j << "]";
    }

    // Check dS_ref/dzle, dMAC/dzle
    for (int j = 0; j < n_st; ++j) {
        auto inp_p = inputs;
        inp_p["zle"](j) += eps;
        auto out_p = compute_refs(inp_p);

        auto inp_m = inputs;
        inp_m["zle"](j) -= eps;
        auto out_m = compute_refs(inp_m);

        double fd_area =
            (out_p["reference_area"](0) - out_m["reference_area"](0)) /
            (2.0 * eps);
        double fd_mac =
            (out_p["reference_chord"](0) - out_m["reference_chord"](0)) /
            (2.0 * eps);

        EXPECT_NEAR(partials[{"reference_area", "zle"}](j), fd_area, 1e-5)
            << "dS_ref/dzle[" << j << "]";
        EXPECT_NEAR(partials[{"reference_chord", "zle"}](j), fd_mac, 1e-5)
            << "dMAC/dzle[" << j << "]";
    }
}

}  // namespace
}  // namespace pegasus
