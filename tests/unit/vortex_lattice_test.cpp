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

#include "pegasus/solver/vortex_lattice.h"

#include "pegasus/geometry/aircraft.h"
#include "pegasus/geometry/section.h"
#include "pegasus/geometry/surface.h"

#include <cmath>

#include <gtest/gtest.h>

namespace pegasus {
namespace {

constexpr double kTestTolerance = 1.0e-6;

// ============================================================================
// VortexLattice Tests
// ============================================================================

class VortexLatticeTest : public ::testing::Test {
protected:
    // Helper to create a simple rectangular wing
    Aircraft CreateRectangularWing(Scalar span, Scalar chord, int num_chordwise,
                                   int num_spanwise) {
        Aircraft aircraft("Rectangular Wing");

        Surface wing("Wing");
        wing.SetChordwisePanelCount(num_chordwise);
        wing.SetSpanwisePanelCount(num_spanwise);
        wing.SetChordwiseSpacing(SpacingType::kCosine);
        wing.SetSpanwiseSpacing(SpacingType::kCosine);

        // Root section at y=0
        wing.AddSection(Section(Vec3(0, 0, 0), chord, 0.0));
        // Tip section at y=span/2
        wing.AddSection(Section(Vec3(0, span / 2, 0), chord, 0.0));

        // Y-duplicate for full wing
        wing.SetYDuplicate(0.0);

        aircraft.AddSurface(wing);
        aircraft.SetReferenceArea(span * chord);
        aircraft.SetReferenceChord(chord);
        aircraft.SetReferenceSpan(span);

        return aircraft;
    }

    // Helper to create a tapered wing
    Aircraft CreateTaperedWing(Scalar span, Scalar root_chord,
                               Scalar taper_ratio, int num_chordwise,
                               int num_spanwise) {
        Aircraft aircraft("Tapered Wing");

        Scalar tip_chord = root_chord * taper_ratio;

        Surface wing("Wing");
        wing.SetChordwisePanelCount(num_chordwise);
        wing.SetSpanwisePanelCount(num_spanwise);
        wing.SetChordwiseSpacing(SpacingType::kCosine);
        wing.SetSpanwiseSpacing(SpacingType::kCosine);

        // Root section at y=0
        wing.AddSection(Section(Vec3(0, 0, 0), root_chord, 0.0));
        // Tip section at y=span/2
        wing.AddSection(Section(Vec3(0, span / 2, 0), tip_chord, 0.0));

        // Y-duplicate for full wing
        wing.SetYDuplicate(0.0);

        aircraft.AddSurface(wing);

        return aircraft;
    }
};

TEST_F(VortexLatticeTest, EmptyAircraft) {
    Aircraft aircraft;
    VortexLattice mesh;

    mesh.GenerateFromAircraft(aircraft);

    EXPECT_EQ(mesh.GetNumVortices(), 0);
    EXPECT_EQ(mesh.GetNumStrips(), 0);
    EXPECT_EQ(mesh.GetNumSurfaces(), 0);
}

TEST_F(VortexLatticeTest, RectangularWingBasicCounts) {
    int num_chordwise = 6;
    int num_spanwise = 8;
    Aircraft aircraft =
        CreateRectangularWing(4.0, 1.0, num_chordwise, num_spanwise);

    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    // Y-duplicate creates two surfaces from one
    EXPECT_EQ(mesh.GetNumSurfaces(), 2);

    // Should have strips and vortices generated
    EXPECT_GT(mesh.GetNumStrips(), 0);
    EXPECT_GT(mesh.GetNumVortices(), 0);

    // Each strip should have num_chordwise vortices
    for (Index j = 0; j < mesh.GetNumStrips(); ++j) {
        EXPECT_EQ(mesh.GetStrip(j).num_vortices, num_chordwise);
    }
}

TEST_F(VortexLatticeTest, RectangularWingGeometry) {
    Scalar span = 4.0;
    Scalar chord = 1.0;
    int num_chordwise = 4;
    int num_spanwise = 4;
    Aircraft aircraft =
        CreateRectangularWing(span, chord, num_chordwise, num_spanwise);

    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    // Check that all strips have consistent chord
    for (Index j = 0; j < mesh.GetNumStrips(); ++j) {
        const Strip& strip = mesh.GetStrip(j);
        EXPECT_NEAR(strip.chord_center, chord, kTestTolerance);
    }

    // Check that vortices are within wing bounds
    for (Index i = 0; i < mesh.GetNumVortices(); ++i) {
        const HorseshoeVortex& vor = mesh.GetVortex(i);

        // X should be between 0 and chord
        EXPECT_GE(vor.control_point.x(), 0.0);
        EXPECT_LE(vor.control_point.x(), chord);

        // Y should be within wing span
        EXPECT_GE(vor.control_point.y(), -span / 2);
        EXPECT_LE(vor.control_point.y(), span / 2);
    }
}

TEST_F(VortexLatticeTest, VortexQuarterChord) {
    // In standard VLM, bound vortex is at 1/4 chord of each panel
    Scalar span = 2.0;
    Scalar chord = 1.0;
    int num_chordwise = 1;  // Single panel
    int num_spanwise = 1;   // Single strip

    Aircraft aircraft("Test");
    Surface wing("Wing");
    wing.SetChordwisePanelCount(num_chordwise);
    wing.SetSpanwisePanelCount(num_spanwise);
    wing.AddSection(Section(Vec3(0, 0, 0), chord, 0.0));
    wing.AddSection(Section(Vec3(0, span, 0), chord, 0.0));
    aircraft.AddSurface(wing);

    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    ASSERT_EQ(mesh.GetNumVortices(), 1);

    const HorseshoeVortex& vor = mesh.GetVortex(0);

    // Bound vortex should be at approximately 1/4 chord
    EXPECT_NEAR(vor.bound_vortex_center.x(), 0.25 * chord, 0.1);

    // Control point should be at approximately 3/4 chord
    EXPECT_NEAR(vor.control_point.x(), 0.75 * chord, 0.1);
}

TEST_F(VortexLatticeTest, ControlPointThreeQuarterChord) {
    // In standard VLM, control point is at 3/4 chord of each panel
    Scalar span = 2.0;
    Scalar chord = 1.0;
    int num_chordwise = 4;
    int num_spanwise = 2;

    Aircraft aircraft("Test");
    Surface wing("Wing");
    wing.SetChordwisePanelCount(num_chordwise);
    wing.SetSpanwisePanelCount(num_spanwise);
    wing.SetChordwiseSpacing(
        SpacingType::kUniform);  // Uniform for easy verification
    wing.AddSection(Section(Vec3(0, 0, 0), chord, 0.0));
    wing.AddSection(Section(Vec3(0, span, 0), chord, 0.0));
    aircraft.AddSurface(wing);

    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    // With uniform spacing and 4 panels, panel width = 0.25
    // Control points should be at panel edges + 0.75 * panel_width
    // Panel 0: [0, 0.25], CP at 0.1875
    // Panel 1: [0.25, 0.5], CP at 0.4375
    // etc.

    // Just check that control points are aft of bound vortices
    for (Index i = 0; i < mesh.GetNumVortices(); ++i) {
        const HorseshoeVortex& vor = mesh.GetVortex(i);
        EXPECT_GT(vor.control_point.x(), vor.bound_vortex_center.x());
    }
}

TEST_F(VortexLatticeTest, NormalVectorsFlatWing) {
    // For a flat wing (no twist, no camber), normals should point up (+z)
    Scalar span = 2.0;
    Scalar chord = 1.0;

    Aircraft aircraft("Test");
    Surface wing("Wing");
    wing.SetChordwisePanelCount(4);
    wing.SetSpanwisePanelCount(4);
    wing.AddSection(Section(Vec3(0, 0, 0), chord, 0.0));
    wing.AddSection(Section(Vec3(0, span, 0), chord, 0.0));
    aircraft.AddSurface(wing);

    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    for (Index i = 0; i < mesh.GetNumVortices(); ++i) {
        const HorseshoeVortex& vor = mesh.GetVortex(i);

        // Normal should be approximately (0, 0, 1) for flat wing
        // Allow some tolerance for numerical effects
        EXPECT_NEAR(vor.normal_at_control_point.z(), 1.0, 0.1);
    }
}

TEST_F(VortexLatticeTest, TwistedWingNormals) {
    // Wing with twist should have tilted normals
    Scalar span = 2.0;
    Scalar chord = 1.0;
    Scalar twist = 5.0;  // degrees

    Aircraft aircraft("Test");
    Surface wing("Wing");
    wing.SetChordwisePanelCount(4);
    wing.SetSpanwisePanelCount(4);

    // Root with positive twist (nose up)
    wing.AddSection(Section(Vec3(0, 0, 0), chord, twist));
    // Tip with zero twist
    wing.AddSection(Section(Vec3(0, span, 0), chord, 0.0));

    aircraft.AddSurface(wing);

    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    // Root strips should have normals tilted forward (positive x-component)
    // Tip strips should have normals pointing up

    // Check that normals vary spanwise
    const Strip& root_strip = mesh.GetStrip(0);
    const Strip& tip_strip = mesh.GetStrip(mesh.GetNumStrips() - 1);

    // Root should have larger twist angle than tip
    EXPECT_GT(std::abs(root_strip.twist_angle),
              std::abs(tip_strip.twist_angle));
}

TEST_F(VortexLatticeTest, YDuplicateSymmetry) {
    // Y-duplicated surface should have symmetric geometry
    Scalar span = 4.0;
    Scalar chord = 1.0;
    int num_chordwise = 4;
    int num_spanwise = 4;

    Aircraft aircraft =
        CreateRectangularWing(span, chord, num_chordwise, num_spanwise);

    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    // Should have two surfaces
    ASSERT_EQ(mesh.GetNumSurfaces(), 2);

    const auto& info1 = mesh.GetSurfaceInfo(0);
    const auto& info2 = mesh.GetSurfaceInfo(1);

    // Both surfaces should have same number of elements
    EXPECT_EQ(info1.num_vortices, info2.num_vortices);
    EXPECT_EQ(info1.num_strips, info2.num_strips);

    // Check that y-coordinates are mirrored
    for (Index j = 0; j < info1.num_strips; ++j) {
        const Strip& strip1 = mesh.GetStrip(info1.first_strip_index + j);
        const Strip& strip2 = mesh.GetStrip(info2.first_strip_index + j);

        // Y should be negated
        EXPECT_NEAR(strip1.leading_edge_center.y(),
                    -strip2.leading_edge_center.y(), kTestTolerance);

        // X and Z should be same
        EXPECT_NEAR(strip1.leading_edge_center.x(),
                    strip2.leading_edge_center.x(), kTestTolerance);
        EXPECT_NEAR(strip1.leading_edge_center.z(),
                    strip2.leading_edge_center.z(), kTestTolerance);
    }
}

TEST_F(VortexLatticeTest, TaperedWingChordDistribution) {
    Scalar span = 4.0;
    Scalar root_chord = 1.0;
    Scalar taper_ratio = 0.5;
    int num_chordwise = 4;
    int num_spanwise = 8;

    Aircraft aircraft = CreateTaperedWing(span, root_chord, taper_ratio,
                                          num_chordwise, num_spanwise);

    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    // Get the first surface (half-wing)
    const auto& info = mesh.GetSurfaceInfo(0);

    // Chord should decrease from root to tip
    Scalar prev_chord = root_chord;
    for (Index j = 0; j < info.num_strips; ++j) {
        const Strip& strip = mesh.GetStrip(info.first_strip_index + j);

        // Current chord should be less than or equal to previous
        // (with some tolerance for interpolation)
        EXPECT_LE(strip.chord_center, prev_chord + 0.01);
        prev_chord = strip.chord_center;
    }

    // Last strip should have chord close to tip chord
    const Strip& tip_strip =
        mesh.GetStrip(info.first_strip_index + info.num_strips - 1);
    Scalar tip_chord = root_chord * taper_ratio;
    EXPECT_NEAR(tip_strip.chord_center, tip_chord, 0.1);
}

TEST_F(VortexLatticeTest, SurfaceInfoConsistency) {
    int num_chordwise = 6;
    int num_spanwise = 10;
    Aircraft aircraft =
        CreateRectangularWing(4.0, 1.0, num_chordwise, num_spanwise);

    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    // Verify surface info is consistent
    for (Index s = 0; s < mesh.GetNumSurfaces(); ++s) {
        const auto& info = mesh.GetSurfaceInfo(s);

        EXPECT_EQ(info.num_chordwise, num_chordwise);
        EXPECT_EQ(info.num_spanwise, num_spanwise);
        EXPECT_EQ(info.num_vortices, info.num_chordwise * info.num_strips);

        // Verify indices are within bounds
        EXPECT_LT(info.first_vortex_index + info.num_vortices,
                  mesh.GetNumVortices() + 1);
        EXPECT_LT(info.first_strip_index + info.num_strips,
                  mesh.GetNumStrips() + 1);
    }
}

TEST_F(VortexLatticeTest, StripVortexIndicesConsistent) {
    int num_chordwise = 4;
    int num_spanwise = 5;
    Aircraft aircraft =
        CreateRectangularWing(4.0, 1.0, num_chordwise, num_spanwise);

    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    // Each strip should have num_chordwise vortices
    for (Index j = 0; j < mesh.GetNumStrips(); ++j) {
        const Strip& strip = mesh.GetStrip(j);
        EXPECT_EQ(strip.num_vortices, num_chordwise);

        // Verify vortex indices are valid
        EXPECT_LT(strip.first_vortex_index + strip.num_vortices,
                  mesh.GetNumVortices() + 1);

        // Verify vortices reference back to this strip
        for (Index k = 0; k < strip.num_vortices; ++k) {
            const HorseshoeVortex& vor =
                mesh.GetVortex(strip.first_vortex_index + k);
            EXPECT_EQ(vor.strip_index, j);
        }
    }
}

TEST_F(VortexLatticeTest, ClearAndRegenerate) {
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 4, 4);

    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    Index num_vor_1 = mesh.GetNumVortices();
    Index num_strip_1 = mesh.GetNumStrips();

    // Clear and regenerate
    mesh.Clear();
    EXPECT_EQ(mesh.GetNumVortices(), 0);
    EXPECT_EQ(mesh.GetNumStrips(), 0);

    mesh.GenerateFromAircraft(aircraft);
    EXPECT_EQ(mesh.GetNumVortices(), num_vor_1);
    EXPECT_EQ(mesh.GetNumStrips(), num_strip_1);
}

TEST_F(VortexLatticeTest, CoreRadiusFractions) {
    VortexLattice mesh;
    EXPECT_NEAR(mesh.GetVortexCoreRadiusChordFraction(), 0.75, 1e-6);
    EXPECT_NEAR(mesh.GetVortexCoreRadiusWidthFraction(), 0.25, 1e-6);

    mesh.SetVortexCoreRadiusChordFraction(0.5);
    mesh.SetVortexCoreRadiusWidthFraction(0.1);
    EXPECT_NEAR(mesh.GetVortexCoreRadiusChordFraction(), 0.5, 1e-6);
    EXPECT_NEAR(mesh.GetVortexCoreRadiusWidthFraction(), 0.1, 1e-6);
}

TEST_F(VortexLatticeTest, GetMeshInputNullByDefault) {
    VortexLattice mesh;
    EXPECT_EQ(mesh.GetMeshInput(), nullptr);
}

TEST_F(VortexLatticeTest, MutableVortexAccess) {
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 2, 2);
    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    // Get mutable vortex and modify it
    HorseshoeVortex& vor = mesh.GetVortex(0);
    Vec3 orig_cp = vor.control_point;
    vor.control_point = Vec3(99, 99, 99);
    EXPECT_NEAR(mesh.GetVortex(0).control_point.x(), 99.0, 1e-6);
    vor.control_point = orig_cp;  // restore
}

TEST_F(VortexLatticeTest, MutableStripAccess) {
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 2, 2);
    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    Strip& strip = mesh.GetStrip(0);
    EXPECT_GE(strip.num_vortices, 1);
}

TEST_F(VortexLatticeTest, GetVorticesMutable) {
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 2, 2);
    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    auto& vortices = mesh.GetVortices();
    EXPECT_FALSE(vortices.empty());
}

TEST_F(VortexLatticeTest, GetStripsMutable) {
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 2, 2);
    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    auto& strips = mesh.GetStrips();
    EXPECT_FALSE(strips.empty());
}

TEST_F(VortexLatticeTest, ControlSurfaceNamesEmpty) {
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 2, 2);
    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    EXPECT_TRUE(mesh.GetControlSurfaceNames().empty());
}

TEST_F(VortexLatticeTest, ControlSurfaceInfluencesSize) {
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 2, 2);
    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    auto& influences = mesh.GetControlSurfaceInfluences();
    EXPECT_EQ(static_cast<Index>(influences.size()), mesh.GetNumVortices());
}

TEST_F(VortexLatticeTest, RecalculateNormals) {
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 2, 2);
    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    // Normals should already be computed, save one
    Vec3 normal0 = mesh.GetVortex(0).normal_at_control_point;

    // Recalculate and verify they're the same
    mesh.RecalculateNormals();
    EXPECT_NEAR((mesh.GetVortex(0).normal_at_control_point - normal0).norm(),
                0.0, 1e-10);
}

}  // namespace
}  // namespace pegasus
