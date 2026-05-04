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
    fields["y_duplicated"].set_string_value(y_duplicated ? "true" : "false");
    return opts;
}

// Helper: set up a MeshDiscipline through the public base-class interface
void SetupDiscipline(MeshDiscipline& disc, int nx, int ny,
                     bool y_duplicated = true) {
    philote::Discipline& dbase = disc;
    dbase.Initialize();
    dbase.SetOptions(MakeOptions(nx, ny, y_duplicated));
    dbase.Setup();
    dbase.SetupPartials();
}

// Helper: build inputs for a given chord/eta/zle distribution
philote::Variables MakeInputs(int ny, const std::vector<double>& chord_vals,
                              double halfspan,
                              const std::vector<double>& eta_vals,
                              const std::vector<double>& zle_vals) {
    int n_st = ny + 1;
    philote::Variables inputs;
    inputs["chord"] =
        philote::Variable(philote::kInput, {static_cast<size_t>(n_st)});
    inputs["twist"] =
        philote::Variable(philote::kInput, {static_cast<size_t>(n_st)});
    inputs["xle"] =
        philote::Variable(philote::kInput, {static_cast<size_t>(n_st)});
    inputs["zle"] =
        philote::Variable(philote::kInput, {static_cast<size_t>(n_st)});
    inputs["halfspan"] = philote::Variable(philote::kInput, {1});
    inputs["eta"] =
        philote::Variable(philote::kInput, {static_cast<size_t>(n_st)});

    for (int j = 0; j < n_st; ++j) {
        inputs["chord"](j) = chord_vals[j];
        inputs["twist"](j) = 0.0;
        inputs["xle"](j) = 0.0;
        inputs["eta"](j) = eta_vals[j];
        inputs["zle"](j) = zle_vals[j];
    }
    inputs["halfspan"](0) = halfspan;
    return inputs;
}

// Helper: build outputs matching the discipline's declared sizes
philote::Variables MakeOutputs(int nx, int ny) {
    int n_st = ny + 1;
    size_t n_mesh = static_cast<size_t>((nx + 1) * n_st * 3);
    philote::Variables outputs;
    outputs["mesh"] = philote::Variable(philote::kOutput, {n_mesh});
    outputs["reference_area"] = philote::Variable(philote::kOutput, {1});
    outputs["reference_span"] = philote::Variable(philote::kOutput, {1});
    outputs["reference_chord"] = philote::Variable(philote::kOutput, {1});
    return outputs;
}

// Helper: build partials matching the discipline's declared structure
philote::Partials MakePartials(int nx, int ny) {
    int n_st = ny + 1;
    size_t n_mesh = static_cast<size_t>((nx + 1) * n_st * 3);
    size_t sn = static_cast<size_t>(n_st);
    philote::Partials partials;

    // Mesh partials (dense: n_mesh rows x n_st cols)
    partials[{"mesh", "chord"}] =
        philote::Variable(philote::kOutput, {n_mesh * sn});
    partials[{"mesh", "twist"}] =
        philote::Variable(philote::kOutput, {n_mesh * sn});
    partials[{"mesh", "xle"}] =
        philote::Variable(philote::kOutput, {n_mesh * sn});
    partials[{"mesh", "zle"}] =
        philote::Variable(philote::kOutput, {n_mesh * sn});
    partials[{"mesh", "eta"}] =
        philote::Variable(philote::kOutput, {n_mesh * sn});
    // halfspan is scalar: n_mesh x 1
    partials[{"mesh", "halfspan"}] =
        philote::Variable(philote::kOutput, {n_mesh});

    // Reference quantity partials (1 row x n_st cols each)
    partials[{"reference_area", "chord"}] =
        philote::Variable(philote::kOutput, {sn});
    partials[{"reference_area", "zle"}] =
        philote::Variable(philote::kOutput, {sn});
    partials[{"reference_area", "halfspan"}] =
        philote::Variable(philote::kOutput, {1});
    partials[{"reference_area", "eta"}] =
        philote::Variable(philote::kOutput, {sn});

    partials[{"reference_span", "halfspan"}] =
        philote::Variable(philote::kOutput, {1});
    partials[{"reference_span", "eta"}] =
        philote::Variable(philote::kOutput, {sn});

    partials[{"reference_chord", "chord"}] =
        philote::Variable(philote::kOutput, {sn});
    partials[{"reference_chord", "zle"}] =
        philote::Variable(philote::kOutput, {sn});
    partials[{"reference_chord", "halfspan"}] =
        philote::Variable(philote::kOutput, {1});
    partials[{"reference_chord", "eta"}] =
        philote::Variable(philote::kOutput, {sn});

    return partials;
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
    SetupDiscipline(disc, nx, ny, true);

    std::vector<double> chords(n_st, chord_val);
    std::vector<double> etas(n_st), zles(n_st, 0.0);
    for (int j = 0; j < n_st; ++j)
        etas[j] = static_cast<double>(j) / ny;

    auto inputs = MakeInputs(ny, chords, half_span, etas, zles);
    auto outputs = MakeOutputs(nx, ny);

    philote::ExplicitDiscipline& base = disc;
    base.Compute(inputs, outputs);

    EXPECT_NEAR(outputs["reference_area"](0), 10.0, 1e-12);
    EXPECT_NEAR(outputs["reference_span"](0), 10.0, 1e-12);
    // MAC of a rectangular wing = chord
    EXPECT_NEAR(outputs["reference_chord"](0), 1.0, 1e-12);
}

TEST(MeshDisciplineRefTest, TaperedWingArea) {
    // Linearly tapered wing: root chord=2, tip chord=1, half-span=4
    // S_ref_half = 0.5*(2+1)*4 = 6, S_ref = 12
    const int nx = 1;
    const int ny = 4;
    const int n_st = ny + 1;
    const double half_span = 4.0;

    MeshDiscipline disc;
    SetupDiscipline(disc, nx, ny, true);

    std::vector<double> chords(n_st), etas(n_st), zles(n_st, 0.0);
    for (int j = 0; j < n_st; ++j) {
        double eta = static_cast<double>(j) / ny;
        chords[j] = 2.0 - 1.0 * eta;
        etas[j] = eta;
    }

    auto inputs = MakeInputs(ny, chords, half_span, etas, zles);
    auto outputs = MakeOutputs(nx, ny);

    philote::ExplicitDiscipline& base = disc;
    base.Compute(inputs, outputs);

    EXPECT_NEAR(outputs["reference_area"](0), 12.0, 1e-12);
    EXPECT_NEAR(outputs["reference_span"](0), 8.0, 1e-12);

    // Compute expected MAC numerically
    double mac_num = 0.0, mac_den = 0.0;
    for (int j = 0; j < ny; ++j) {
        double c0 = chords[j], c1 = chords[j + 1];
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
    SetupDiscipline(disc, nx, ny, true);

    // Non-uniform chord for interesting partials
    std::vector<double> chords(n_st), etas(n_st), zles(n_st);
    for (int j = 0; j < n_st; ++j) {
        double e = static_cast<double>(j) / ny;
        chords[j] = 2.0 - 0.8 * e;
        etas[j] = e;
        zles[j] = 0.1 * e * e;  // slight dihedral
    }

    auto inputs = MakeInputs(ny, chords, half_span, etas, zles);
    auto partials = MakePartials(nx, ny);

    philote::ExplicitDiscipline& base = disc;
    base.ComputePartials(inputs, partials);

    // Helper lambda to compute reference outputs
    auto compute_refs = [&](philote::Variables& inp) {
        auto out = MakeOutputs(nx, ny);
        base.Compute(inp, out);
        return out;
    };

    const double eps = 1e-7;

    // Extract partial references
    auto& dArea_dc = partials[{"reference_area", "chord"}];
    auto& dArea_dz = partials[{"reference_area", "zle"}];
    auto& dArea_dhs = partials[{"reference_area", "halfspan"}];
    auto& dArea_deta = partials[{"reference_area", "eta"}];
    auto& dSpan_dhs = partials[{"reference_span", "halfspan"}];
    auto& dSpan_deta = partials[{"reference_span", "eta"}];
    auto& dMAC_dc = partials[{"reference_chord", "chord"}];
    auto& dMAC_dz = partials[{"reference_chord", "zle"}];
    auto& dMAC_dhs = partials[{"reference_chord", "halfspan"}];
    auto& dMAC_deta = partials[{"reference_chord", "eta"}];

    // Check dS_ref/dchord, dMAC/dchord
    for (int j = 0; j < n_st; ++j) {
        auto inp_p = MakeInputs(ny, chords, half_span, etas, zles);
        inp_p["chord"](j) += eps;
        auto out_p = compute_refs(inp_p);

        auto inp_m = MakeInputs(ny, chords, half_span, etas, zles);
        inp_m["chord"](j) -= eps;
        auto out_m = compute_refs(inp_m);

        double fd_area =
            (out_p["reference_area"](0) - out_m["reference_area"](0)) /
            (2.0 * eps);
        double fd_mac =
            (out_p["reference_chord"](0) - out_m["reference_chord"](0)) /
            (2.0 * eps);

        EXPECT_NEAR(dArea_dc(j), fd_area, 1e-5) << "dS_ref/dchord[" << j << "]";
        EXPECT_NEAR(dMAC_dc(j), fd_mac, 1e-5) << "dMAC/dchord[" << j << "]";
    }

    // Check partials w.r.t. halfspan
    {
        auto inp_p = MakeInputs(ny, chords, half_span + eps, etas, zles);
        auto out_p = compute_refs(inp_p);
        auto inp_m = MakeInputs(ny, chords, half_span - eps, etas, zles);
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

        EXPECT_NEAR(dArea_dhs(0), fd_area, 1e-5) << "dS_ref/d(halfspan)";
        EXPECT_NEAR(dSpan_dhs(0), fd_span, 1e-5) << "dspan/d(halfspan)";
        EXPECT_NEAR(dMAC_dhs(0), fd_mac, 1e-5) << "dMAC/d(halfspan)";
    }

    // Check partials w.r.t. eta
    for (int j = 0; j < n_st; ++j) {
        auto inp_p = MakeInputs(ny, chords, half_span, etas, zles);
        inp_p["eta"](j) += eps;
        auto out_p = compute_refs(inp_p);

        auto inp_m = MakeInputs(ny, chords, half_span, etas, zles);
        inp_m["eta"](j) -= eps;
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

        EXPECT_NEAR(dArea_deta(j), fd_area, 1e-5) << "dS_ref/deta[" << j << "]";
        EXPECT_NEAR(dSpan_deta(j), fd_span, 1e-5) << "dspan/deta[" << j << "]";
        EXPECT_NEAR(dMAC_deta(j), fd_mac, 1e-5) << "dMAC/deta[" << j << "]";
    }

    // Check dS_ref/dzle, dMAC/dzle
    for (int j = 0; j < n_st; ++j) {
        auto inp_p = MakeInputs(ny, chords, half_span, etas, zles);
        inp_p["zle"](j) += eps;
        auto out_p = compute_refs(inp_p);

        auto inp_m = MakeInputs(ny, chords, half_span, etas, zles);
        inp_m["zle"](j) -= eps;
        auto out_m = compute_refs(inp_m);

        double fd_area =
            (out_p["reference_area"](0) - out_m["reference_area"](0)) /
            (2.0 * eps);
        double fd_mac =
            (out_p["reference_chord"](0) - out_m["reference_chord"](0)) /
            (2.0 * eps);

        EXPECT_NEAR(dArea_dz(j), fd_area, 1e-5) << "dS_ref/dzle[" << j << "]";
        EXPECT_NEAR(dMAC_dz(j), fd_mac, 1e-5) << "dMAC/dzle[" << j << "]";
    }

    // Check mesh partials w.r.t. halfspan (spot-check y components)
    {
        auto& dm_dhs = partials[{"mesh", "halfspan"}];
        auto inp_p = MakeInputs(ny, chords, half_span + eps, etas, zles);
        auto out_p = MakeOutputs(nx, ny);
        base.Compute(inp_p, out_p);

        auto inp_m = MakeInputs(ny, chords, half_span - eps, etas, zles);
        auto out_m = MakeOutputs(nx, ny);
        base.Compute(inp_m, out_m);

        int n_mesh = (nx + 1) * n_st * 3;
        for (int k = 0; k < n_mesh; ++k) {
            double fd = (out_p["mesh"](k) - out_m["mesh"](k)) / (2.0 * eps);
            EXPECT_NEAR(dm_dhs(k), fd, 1e-5)
                << "d(mesh[" << k << "])/d(halfspan)";
        }
    }

    // Check mesh partials w.r.t. eta
    {
        auto& dm_deta = partials[{"mesh", "eta"}];
        int n_mesh = (nx + 1) * n_st * 3;

        for (int j = 0; j < n_st; ++j) {
            auto inp_p = MakeInputs(ny, chords, half_span, etas, zles);
            inp_p["eta"](j) += eps;
            auto out_p = MakeOutputs(nx, ny);
            base.Compute(inp_p, out_p);

            auto inp_m = MakeInputs(ny, chords, half_span, etas, zles);
            inp_m["eta"](j) -= eps;
            auto out_m = MakeOutputs(nx, ny);
            base.Compute(inp_m, out_m);

            for (int k = 0; k < n_mesh; ++k) {
                double fd = (out_p["mesh"](k) - out_m["mesh"](k)) / (2.0 * eps);
                double analytical = dm_deta(k * n_st + j);
                EXPECT_NEAR(analytical, fd, 1e-5)
                    << "d(mesh[" << k << "])/d(eta[" << j << "])";
            }
        }
    }
}

}  // namespace
}  // namespace pegasus
