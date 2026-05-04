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

#include "pegasus/gradients/mesh_sensitivity.h"

#include "pegasus/core/types.h"
#include "pegasus/geometry/aircraft.h"
#include "pegasus/geometry/section.h"
#include "pegasus/geometry/surface.h"
#include "pegasus/solver/vortex_lattice.h"
#include "pegasus/solver/vortex_lattice_solver.h"

#include <cmath>

#include <gtest/gtest.h>

namespace pegasus {
namespace gradients {
namespace {

class MeshSensitivityTest : public ::testing::Test {
protected:
    void SetUp() override {
        Surface wing("Wing");
        wing.SetChordwisePanelCount(4);
        wing.SetSpanwisePanelCount(8);
        wing.SetYDuplicate(0.0);

        Section root(Vec3(0.0, 0.0, 0.0), 1.0);
        Section tip(Vec3(0.0, 5.0, 0.0), 1.0);
        wing.AddSection(root);
        wing.AddSection(tip);

        aircraft_.SetTitle("Test Wing");
        aircraft_.AddSurface(wing);
        aircraft_.SetReferenceArea(10.0);
        aircraft_.SetReferenceChord(1.0);
        aircraft_.SetReferenceSpan(10.0);

        solver_.SetAircraft(aircraft_);
    }

    Aircraft aircraft_;
    VortexLatticeSolver solver_;
};

TEST_F(MeshSensitivityTest, ExtractMeshArrays) {
    const VortexLattice& mesh = solver_.GetMesh();
    MeshArrays arrays = ExtractMeshArrays(mesh);

    EXPECT_EQ(arrays.num_vortices, mesh.GetNumVortices());
    EXPECT_EQ(static_cast<Index>(arrays.rv1.size()), 3 * arrays.num_vortices);
    EXPECT_EQ(static_cast<Index>(arrays.rv2.size()), 3 * arrays.num_vortices);
    EXPECT_EQ(static_cast<Index>(arrays.rc.size()), 3 * arrays.num_vortices);
    EXPECT_EQ(static_cast<Index>(arrays.normal.size()),
              3 * arrays.num_vortices);
    EXPECT_EQ(static_cast<Index>(arrays.chord.size()), arrays.num_vortices);
    EXPECT_EQ(static_cast<Index>(arrays.surface_idx.size()),
              arrays.num_vortices);
    EXPECT_EQ(static_cast<Index>(arrays.component_idx.size()),
              arrays.num_vortices);

    // Verify positions are finite
    for (Index i = 0; i < 3 * arrays.num_vortices; ++i) {
        EXPECT_TRUE(std::isfinite(arrays.rv1[i]));
        EXPECT_TRUE(std::isfinite(arrays.rv2[i]));
        EXPECT_TRUE(std::isfinite(arrays.rc[i]));
        EXPECT_TRUE(std::isfinite(arrays.normal[i]));
    }
}

TEST_F(MeshSensitivityTest, BuildSectionStripMappingSingleSurface) {
    const VortexLattice& mesh = solver_.GetMesh();

    // Surface 0 is the original half, surface 1 is the y-duplicate
    SectionStripMapping mapping = BuildSectionStripMapping(aircraft_, mesh, 0);

    EXPECT_EQ(mapping.surface_index, 0);
    EXPECT_EQ(mapping.num_sections, 2);
    EXPECT_GT(mapping.num_strips, 0);
    EXPECT_EQ(static_cast<Index>(mapping.strip_at_section.size()),
              mapping.num_sections);
    EXPECT_EQ(mapping.strip_at_section[0], 0);
}

TEST_F(MeshSensitivityTest, BuildSectionStripMappingYDuplicate) {
    const VortexLattice& mesh = solver_.GetMesh();

    // Surface 1 is the y-duplicate
    SectionStripMapping mapping = BuildSectionStripMapping(aircraft_, mesh, 1);

    EXPECT_EQ(mapping.surface_index, 1);
    EXPECT_EQ(mapping.num_sections, 2);
    EXPECT_GT(mapping.num_strips, 0);
}

TEST_F(MeshSensitivityTest, BuildDefaultChordwiseFractions) {
    std::vector<double> xvr, xcp;
    BuildDefaultChordwiseFractions(4, xvr, xcp);

    EXPECT_EQ(xvr.size(), 4u);
    EXPECT_EQ(xcp.size(), 4u);

    // Verify 1/4 - 3/4 chord rule within each panel
    for (size_t i = 0; i < xvr.size(); ++i) {
        EXPECT_LT(xvr[i], xcp[i]);
        EXPECT_GE(xvr[i], 0.0);
        EXPECT_LE(xcp[i], 1.0);
    }

    // First vortex should be near the leading edge
    EXPECT_LT(xvr[0], 0.2);
}

TEST_F(MeshSensitivityTest, StripInterpolationFraction) {
    // Single strip: should be 0.5
    double fc = ComputeStripInterpolationFraction(0, 1);
    EXPECT_DOUBLE_EQ(fc, 0.5);

    // Multiple strips: first should be near 0, last near 1
    double fc_first = ComputeStripInterpolationFraction(0, 8);
    double fc_last = ComputeStripInterpolationFraction(7, 8);
    EXPECT_LT(fc_first, 0.15);
    EXPECT_GT(fc_last, 0.85);
}

TEST_F(MeshSensitivityTest, StripEdgeFractions) {
    double f1, f2;
    ComputeStripEdgeFractions(0, 4, &f1, &f2);

    EXPECT_GE(f1, 0.0);
    EXPECT_LE(f2, 1.0);
    EXPECT_LT(f1, f2);

    // Last strip should end at 1.0
    ComputeStripEdgeFractions(3, 4, &f1, &f2);
    EXPECT_NEAR(f2, 1.0, 1e-10);
}

TEST_F(MeshSensitivityTest, BuildSpanwiseFractions) {
    std::vector<double> f1, f2, fc;
    BuildSpanwiseFractions(4, SpacingType::kCosine, f1, f2, fc);

    EXPECT_EQ(f1.size(), 4u);
    EXPECT_EQ(f2.size(), 4u);
    EXPECT_EQ(fc.size(), 4u);

    // First edge at 0, last edge at 1
    EXPECT_NEAR(f1[0], 0.0, 1e-10);
    EXPECT_NEAR(f2[3], 1.0, 1e-10);

    // Centers between edges
    for (size_t i = 0; i < 4; ++i) {
        EXPECT_GE(fc[i], f1[i]);
        EXPECT_LE(fc[i], f2[i]);
    }
}

TEST_F(MeshSensitivityTest, SectionStripMappingGetStripRange) {
    SectionStripMapping mapping;
    mapping.first_strip = 10;
    mapping.num_strips = 8;
    mapping.num_sections = 3;
    mapping.strip_at_section = {0, 4, 8};

    auto [first, count] = mapping.GetStripRange(0);
    EXPECT_EQ(first, 10);
    EXPECT_EQ(count, 4);

    auto [first2, count2] = mapping.GetStripRange(1);
    EXPECT_EQ(first2, 14);
    EXPECT_EQ(count2, 4);

    // Out of range
    auto [first3, count3] = mapping.GetStripRange(2);
    EXPECT_EQ(count3, 0);
}

TEST_F(MeshSensitivityTest, MeshGradientsResizeAndZero) {
    MeshGradients grads;
    grads.Resize(10);

    EXPECT_EQ(grads.d_rv1.size(), 30u);
    EXPECT_EQ(grads.d_rv2.size(), 30u);
    EXPECT_EQ(grads.d_rc.size(), 30u);
    EXPECT_EQ(grads.d_normal.size(), 30u);

    // Should be zero-initialized
    for (size_t i = 0; i < 30; ++i) {
        EXPECT_DOUBLE_EQ(grads.d_rv1[i], 0.0);
    }

    grads.d_rv1[0] = 1.0;
    grads.Zero();
    EXPECT_DOUBLE_EQ(grads.d_rv1[0], 0.0);
}

TEST_F(MeshSensitivityTest, PropagateLeadingEdgeGradient) {
    double d_rle[3] = {1.0, 2.0, 3.0};
    double d_le_left[3] = {0.0, 0.0, 0.0};
    double d_le_right[3] = {0.0, 0.0, 0.0};

    PropagateLeadingEdgeGradient(d_rle, 0.3, d_le_left, d_le_right);

    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(d_le_left[i], 0.7 * d_rle[i], 1e-12);
        EXPECT_NEAR(d_le_right[i], 0.3 * d_rle[i], 1e-12);
    }
}

TEST_F(MeshSensitivityTest, PropagateChordGradient) {
    double d_rv1[3] = {1.0, 0.0, 0.0};
    double d_rv2[3] = {0.5, 0.0, 0.0};
    double d_rc[3] = {0.8, 0.0, 0.0};
    double d_chord_left = 0.0;
    double d_chord_right = 0.0;

    PropagateChordGradient(d_rv1, d_rv2, d_rc, 0.25, 0.75, 0.0, 1.0, 0.5,
                           &d_chord_left, &d_chord_right);

    // Verify accumulation is non-zero
    EXPECT_NE(d_chord_left, 0.0);
    EXPECT_NE(d_chord_right, 0.0);
}

TEST_F(MeshSensitivityTest, ComputeSectionGradientFromMesh) {
    const VortexLattice& mesh = solver_.GetMesh();
    SectionStripMapping mapping = BuildSectionStripMapping(aircraft_, mesh, 0);

    std::vector<double> xvr, xcp;
    Index num_chordwise = mesh.GetSurfaceInfo(0).num_chordwise;
    BuildDefaultChordwiseFractions(num_chordwise, xvr, xcp);

    // Create mesh gradients with a simple perturbation
    MeshGradients grads;
    grads.Resize(mesh.GetNumVortices());
    // Set a small gradient on the first vortex normal
    if (mesh.GetNumVortices() > 0) {
        grads.d_normal[2] = 1.0;
    }

    SectionGradients result =
        ComputeSectionGradientFromMesh(grads, mapping, 0, mesh, xvr, xcp);

    // Results should be finite
    EXPECT_TRUE(std::isfinite(result.d_le_x));
    EXPECT_TRUE(std::isfinite(result.d_le_y));
    EXPECT_TRUE(std::isfinite(result.d_le_z));
    EXPECT_TRUE(std::isfinite(result.d_chord));
    EXPECT_TRUE(std::isfinite(result.d_twist));
}

TEST_F(MeshSensitivityTest, PropagateTwistGradient) {
    double d_normal[3] = {0.0, 0.0, 1.0};
    double bound_dir[3] = {0.0, 1.0, 0.0};

    double d_twist =
        PropagateTwistGradient(d_normal, 0.05, 0.0, 1.0, bound_dir, 0.0);
    EXPECT_TRUE(std::isfinite(d_twist));
}

// ============================================================================
// Multi-surface configuration
// ============================================================================

TEST(MeshSensitivityMultiSurfaceTest, TwoSurfaces) {
    Surface wing("Wing");
    wing.SetChordwisePanelCount(4);
    wing.SetSpanwisePanelCount(6);
    wing.SetYDuplicate(0.0);

    Section wing_root(Vec3(0.0, 0.0, 0.0), 1.0);
    Section wing_tip(Vec3(0.0, 5.0, 0.0), 1.0);
    wing.AddSection(wing_root);
    wing.AddSection(wing_tip);

    Surface tail("Tail");
    tail.SetChordwisePanelCount(4);
    tail.SetSpanwisePanelCount(4);
    tail.SetYDuplicate(0.0);

    Section tail_root(Vec3(4.0, 0.0, 0.0), 0.5);
    Section tail_tip(Vec3(4.0, 1.5, 0.0), 0.3);
    tail.AddSection(tail_root);
    tail.AddSection(tail_tip);

    Aircraft aircraft("Wing-Tail");
    aircraft.AddSurface(wing);
    aircraft.AddSurface(tail);
    aircraft.SetReferenceArea(10.0);
    aircraft.SetReferenceChord(1.0);
    aircraft.SetReferenceSpan(10.0);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    const VortexLattice& mesh = solver.GetMesh();
    MeshArrays arrays = ExtractMeshArrays(mesh);

    // Verify surfaces are indexed correctly
    bool found_surface_0 = false, found_surface_1 = false;
    for (Index i = 0; i < arrays.num_vortices; ++i) {
        if (arrays.surface_idx[i] == 0)
            found_surface_0 = true;
        if (arrays.surface_idx[i] == 1)
            found_surface_1 = true;
    }
    EXPECT_TRUE(found_surface_0);
    EXPECT_TRUE(found_surface_1);

    // Build mappings for both wing surfaces (original + y-dup)
    SectionStripMapping mapping0 = BuildSectionStripMapping(aircraft, mesh, 0);
    EXPECT_GT(mapping0.num_strips, 0);

    SectionStripMapping mapping1 = BuildSectionStripMapping(aircraft, mesh, 1);
    EXPECT_GT(mapping1.num_strips, 0);
}

}  // namespace
}  // namespace gradients
}  // namespace pegasus
