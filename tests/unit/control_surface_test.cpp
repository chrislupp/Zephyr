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

#include "pegasus/aero/force_calculator.h"
#include "pegasus/aero/stability_derivatives.h"
#include "pegasus/geometry/aircraft.h"
#include "pegasus/geometry/section.h"
#include "pegasus/geometry/surface.h"
#include "pegasus/io/avl_reader.h"
#include "pegasus/solver/vortex_lattice.h"
#include "pegasus/solver/vortex_lattice_solver.h"

#include <cmath>
#include <string>

#include <gtest/gtest.h>

namespace pegasus {
namespace {

constexpr Scalar kTestTolerance = 1e-10;

// ============================================================================
// Control Surface Definition Tests
// ============================================================================

class ControlSurfaceDefinitionTest : public ::testing::Test {};

TEST_F(ControlSurfaceDefinitionTest, NameDefaultsToEmpty) {
    ControlSurfaceDefinition control;

    // Name is std::string, so it default-constructs to empty
    EXPECT_TRUE(control.name.empty());
}

TEST_F(ControlSurfaceDefinitionTest, SetValues) {
    ControlSurfaceDefinition control;
    control.name = "aileron";
    control.hinge_x_fraction = 0.75;
    control.hinge_vector = Vec3(0, 1, 0);
    control.sign = 1.0;
    control.gain = 1.0;

    EXPECT_EQ(control.name, "aileron");
    EXPECT_NEAR(control.hinge_x_fraction, 0.75, kTestTolerance);
    EXPECT_EQ(control.hinge_vector, Vec3(0, 1, 0));
    EXPECT_NEAR(control.sign, 1.0, kTestTolerance);
    EXPECT_NEAR(control.gain, 1.0, kTestTolerance);
}

// ============================================================================
// Section Control Surface Tests
// ============================================================================

class SectionControlSurfaceTest : public ::testing::Test {};

TEST_F(SectionControlSurfaceTest, NoControlSurfacesByDefault) {
    Section section(Vec3(0, 0, 0), 1.0, 0.0);

    EXPECT_EQ(section.GetNumControlSurfaces(), 0);
    EXPECT_TRUE(section.GetControlSurfaces().empty());
}

TEST_F(SectionControlSurfaceTest, AddSingleControlSurface) {
    Section section(Vec3(0, 0, 0), 1.0, 0.0);

    ControlSurfaceDefinition aileron;
    aileron.name = "aileron";
    aileron.hinge_x_fraction = 0.75;
    aileron.hinge_vector = Vec3(0, 1, 0);
    aileron.gain = 1.0;
    aileron.sign = 1.0;

    section.AddControlSurface(aileron);

    EXPECT_EQ(section.GetNumControlSurfaces(), 1);
    const auto& controls = section.GetControlSurfaces();
    EXPECT_EQ(controls[0].name, "aileron");
    EXPECT_NEAR(controls[0].hinge_x_fraction, 0.75, kTestTolerance);
}

TEST_F(SectionControlSurfaceTest, AddMultipleControlSurfaces) {
    Section section(Vec3(0, 0, 0), 1.0, 0.0);

    ControlSurfaceDefinition flap;
    flap.name = "flap";
    flap.hinge_x_fraction = 0.70;
    flap.hinge_vector = Vec3(0, 1, 0);
    flap.gain = 1.0;
    flap.sign = 1.0;

    ControlSurfaceDefinition aileron;
    aileron.name = "aileron";
    aileron.hinge_x_fraction = 0.80;
    aileron.hinge_vector = Vec3(0, 1, 0);
    aileron.gain = 1.0;
    aileron.sign = -1.0;  // Antisymmetric

    section.AddControlSurface(flap);
    section.AddControlSurface(aileron);

    EXPECT_EQ(section.GetNumControlSurfaces(), 2);
    EXPECT_EQ(section.GetControlSurfaces()[0].name, "flap");
    EXPECT_EQ(section.GetControlSurfaces()[1].name, "aileron");
}

// ============================================================================
// Aircraft Control Surface Collection Tests
// ============================================================================

class AircraftControlSurfaceTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a wing with control surfaces
        Surface wing("Wing");
        wing.SetChordwisePanelCount(4);
        wing.SetSpanwisePanelCount(8);

        // Root section with elevator
        Section root(Vec3(0, 0, 0), 1.0, 0.0);
        ControlSurfaceDefinition elevator;
        elevator.name = "elevator";
        elevator.hinge_x_fraction = 0.75;
        elevator.hinge_vector = Vec3(0, 1, 0);
        elevator.gain = 1.0;
        elevator.sign = 1.0;
        root.AddControlSurface(elevator);

        // Tip section with elevator (same control surface spans both sections)
        Section tip(Vec3(0, 4.0, 0), 1.0, 0.0);
        tip.AddControlSurface(elevator);

        wing.AddSection(root);
        wing.AddSection(tip);

        aircraft_.AddSurface(wing);
    }

    Aircraft aircraft_{"Test Aircraft"};
};

TEST_F(AircraftControlSurfaceTest, GetControlSurfaceNames) {
    auto names = aircraft_.GetControlSurfaceNames();

    ASSERT_EQ(names.size(), 1);
    EXPECT_EQ(names[0], "elevator");
}

TEST_F(AircraftControlSurfaceTest, GetNumControlSurfaces) {
    EXPECT_EQ(aircraft_.GetNumControlSurfaces(), 1);
}

TEST_F(AircraftControlSurfaceTest, MultipleControlSurfaces) {
    // Add another surface with different control surfaces
    Surface htail("HTail");
    htail.SetChordwisePanelCount(4);
    htail.SetSpanwisePanelCount(4);

    Section htail_root(Vec3(3.0, 0, 0), 0.5, 0.0);
    Section htail_tip(Vec3(3.0, 1.0, 0), 0.5, 0.0);

    ControlSurfaceDefinition rudder;
    rudder.name = "rudder";
    rudder.hinge_x_fraction = 0.6;
    rudder.hinge_vector = Vec3(0, 0, 1);
    rudder.gain = 1.0;
    rudder.sign = 1.0;

    htail_root.AddControlSurface(rudder);
    htail_tip.AddControlSurface(rudder);

    htail.AddSection(htail_root);
    htail.AddSection(htail_tip);

    aircraft_.AddSurface(htail);

    auto names = aircraft_.GetControlSurfaceNames();
    EXPECT_EQ(names.size(), 2);

    // Names should be sorted
    bool has_elevator = false;
    bool has_rudder = false;
    for (const auto& name : names) {
        if (name == "elevator")
            has_elevator = true;
        if (name == "rudder")
            has_rudder = true;
    }
    EXPECT_TRUE(has_elevator);
    EXPECT_TRUE(has_rudder);
}

// ============================================================================
// Mesh Control Surface Influence Tests
// ============================================================================

class MeshControlSurfaceTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a simple wing with a trailing edge control surface
        Surface wing("Wing");
        wing.SetChordwisePanelCount(4);  // 4 chordwise panels
        wing.SetSpanwisePanelCount(4);   // 4 spanwise panels

        // Control surface from 75% chord to trailing edge
        ControlSurfaceDefinition aileron;
        aileron.name = "aileron";
        aileron.hinge_x_fraction = 0.75;
        aileron.hinge_vector = Vec3(0, 1, 0);
        aileron.gain = 1.0;
        aileron.sign = 1.0;

        Section root(Vec3(0, 0, 0), 1.0, 0.0);
        root.AddControlSurface(aileron);

        Section tip(Vec3(0, 2.0, 0), 1.0, 0.0);
        tip.AddControlSurface(aileron);

        wing.AddSection(root);
        wing.AddSection(tip);

        aircraft_.AddSurface(wing);
        mesh_.GenerateFromAircraft(aircraft_);
    }

    Aircraft aircraft_{"Test Aircraft"};
    VortexLattice mesh_;
};

TEST_F(MeshControlSurfaceTest, ControlSurfaceNamesInMesh) {
    const auto& names = mesh_.GetControlSurfaceNames();

    ASSERT_EQ(names.size(), 1);
    EXPECT_EQ(names[0], "aileron");
}

TEST_F(MeshControlSurfaceTest, ControlInfluencesHaveCorrectSize) {
    const auto& influences = mesh_.GetControlSurfaceInfluences();

    // Should have one influence vector per vortex element
    EXPECT_EQ(influences.size(), mesh_.GetNumVortices());
}

TEST_F(MeshControlSurfaceTest, TrailingEdgePanelsHaveInfluence) {
    // With 4 chordwise panels and hinge at 75%, the last panel (75-100% chord)
    // should have significant control surface influence
    const auto& influences = mesh_.GetControlSurfaceInfluences();

    int panels_with_influence = 0;
    for (const auto& vortex_influences : influences) {
        for (const auto& inf : vortex_influences) {
            if (inf.name == "aileron" && std::abs(inf.gain) > 0.01) {
                ++panels_with_influence;
            }
        }
    }

    // At least some panels should have aileron influence
    EXPECT_GT(panels_with_influence, 0);
}

TEST_F(MeshControlSurfaceTest, LeadingEdgePanelsNoInfluence) {
    // Panels ahead of the hinge line should have zero or near-zero influence
    const auto& influences = mesh_.GetControlSurfaceInfluences();

    // Check the first chordwise panel of each strip (should be near LE)
    for (Index i = 0; i < mesh_.GetNumStrips(); ++i) {
        const Strip& strip = mesh_.GetStrip(i);
        Index first_vortex = strip.first_vortex_index;

        // First vortex in strip is at leading edge
        const auto& le_influences = influences[first_vortex];
        for (const auto& inf : le_influences) {
            if (inf.name == "aileron") {
                // LE panel should have zero or very small influence
                EXPECT_LT(std::abs(inf.gain), 0.5)
                    << "Leading edge panel should not have significant control "
                       "influence";
            }
        }
    }
}

TEST_F(MeshControlSurfaceTest, InfluenceHingeAxisCorrect) {
    const auto& influences = mesh_.GetControlSurfaceInfluences();

    for (const auto& vortex_influences : influences) {
        for (const auto& inf : vortex_influences) {
            if (inf.name == "aileron" && std::abs(inf.gain) > 0.01) {
                // Hinge axis should be approximately spanwise (0, 1, 0)
                EXPECT_NEAR(inf.hinge_axis.x(), 0.0, 0.1);
                EXPECT_GT(std::abs(inf.hinge_axis.y()), 0.9);
                EXPECT_NEAR(inf.hinge_axis.z(), 0.0, 0.1);
            }
        }
    }
}

// ============================================================================
// Y-Duplicate Control Surface Tests
// ============================================================================

class YDuplicateControlSurfaceTest : public ::testing::Test {
protected:
    void SetUp() override {
        Surface wing("Wing");
        wing.SetChordwisePanelCount(4);
        wing.SetSpanwisePanelCount(4);
        wing.SetYDuplicate(0.0);  // Mirror about Y=0

        // Aileron with SgnDup = -1 (antisymmetric deflection)
        ControlSurfaceDefinition aileron;
        aileron.name = "aileron";
        aileron.hinge_x_fraction = 0.75;
        aileron.hinge_vector = Vec3(0, 1, 0);
        aileron.gain = 1.0;
        aileron.sign = -1.0;  // Antisymmetric

        Section root(Vec3(0, 0, 0), 1.0, 0.0);
        root.AddControlSurface(aileron);

        Section tip(Vec3(0, 2.0, 0), 1.0, 0.0);
        tip.AddControlSurface(aileron);

        wing.AddSection(root);
        wing.AddSection(tip);

        aircraft_.AddSurface(wing);
        mesh_.GenerateFromAircraft(aircraft_);
    }

    Aircraft aircraft_{"Test Aircraft"};
    VortexLattice mesh_;
};

TEST_F(YDuplicateControlSurfaceTest, MirroredSurfaceCreated) {
    // With y-duplicate, should have 2 surfaces in mesh
    EXPECT_EQ(mesh_.GetNumSurfaces(), 2);
}

TEST_F(YDuplicateControlSurfaceTest, MirroredInfluencesExist) {
    const auto& influences = mesh_.GetControlSurfaceInfluences();

    // Should have influences for all vortices (both original and mirrored)
    EXPECT_EQ(influences.size(), mesh_.GetNumVortices());
}

TEST_F(YDuplicateControlSurfaceTest, MirroredControlSurfaceInfluences) {
    // For control surfaces with sign = -1 (antisymmetric), the gain computation
    // is: new_gain = -orig_gain * reflection_sign where reflection_sign = sign
    // from ControlSurfaceDefinition
    //
    // So if orig_gain = 1.0 and sign = -1: new_gain = -1.0 * -1.0 = 1.0
    // The antisymmetric behavior comes from the multiplication, not sign flip
    // in storage
    const auto& influences = mesh_.GetControlSurfaceInfluences();

    const SurfaceMeshInfo& orig = mesh_.GetSurfaceInfo(0);
    const SurfaceMeshInfo& mirror = mesh_.GetSurfaceInfo(1);

    // Verify that mirrored influences exist and have matching control names
    bool found_matching_influence = false;
    for (Index i = 0; i < orig.num_vortices && i < mirror.num_vortices; ++i) {
        Index orig_idx = orig.first_vortex_index + i;
        Index mirror_idx = mirror.first_vortex_index + i;

        const auto& orig_inf = influences[orig_idx];
        const auto& mirror_inf = influences[mirror_idx];

        for (size_t j = 0; j < orig_inf.size() && j < mirror_inf.size(); ++j) {
            if (orig_inf[j].name == "aileron" &&
                mirror_inf[j].name == "aileron") {
                found_matching_influence = true;
                // Verify both have non-zero gains on trailing edge panels
                if (std::abs(orig_inf[j].gain) > 0.01) {
                    EXPECT_GT(std::abs(mirror_inf[j].gain), 0.01)
                        << "Mirrored surface should also have control "
                           "influence";
                }
            }
        }
    }
    EXPECT_TRUE(found_matching_influence)
        << "Should find matching aileron influences on both surfaces";
}

// ============================================================================
// AVL Reader Control Surface Parsing Tests
// ============================================================================

class AvlControlSurfaceParsingTest : public ::testing::Test {};

TEST_F(AvlControlSurfaceParsingTest, ParseBasicControlSurface) {
    std::string avl_content = R"(
Test Aircraft
0.0
0 0 0
10.0 1.0 10.0
0.25 0.0 0.0

SURFACE
Wing
4 1.0

SECTION
0.0 0.0 0.0 1.0 0.0
CONTROL
elevator 1.0 0.75 0 1 0 1

SECTION
0.0 2.0 0.0 1.0 0.0
CONTROL
elevator 1.0 0.75 0 1 0 1
)";

    AvlReader reader;
    Aircraft aircraft = reader.ParseString(avl_content);

    EXPECT_EQ(aircraft.GetNumControlSurfaces(), 1);

    auto names = aircraft.GetControlSurfaceNames();
    ASSERT_EQ(names.size(), 1);
    EXPECT_EQ(names[0], "elevator");

    // Check the control surface was added to sections
    const Surface& wing = aircraft.GetSurface(0);
    EXPECT_EQ(wing.GetSection(0).GetNumControlSurfaces(), 1);
    EXPECT_EQ(wing.GetSection(0).GetControlSurfaces()[0].name, "elevator");
}

TEST_F(AvlControlSurfaceParsingTest, ParseMultipleControlSurfaces) {
    std::string avl_content = R"(
Test Aircraft
0.0
0 0 0
10.0 1.0 10.0
0.25 0.0 0.0

SURFACE
Wing
4 1.0

SECTION
0.0 0.0 0.0 1.0 0.0
CONTROL
flap 1.0 0.70 0 1 0 1
CONTROL
aileron 1.0 0.80 0 1 0 -1

SECTION
0.0 2.0 0.0 1.0 0.0
CONTROL
flap 1.0 0.70 0 1 0 1
CONTROL
aileron 1.0 0.80 0 1 0 -1
)";

    AvlReader reader;
    Aircraft aircraft = reader.ParseString(avl_content);

    EXPECT_EQ(aircraft.GetNumControlSurfaces(), 2);

    const Surface& wing = aircraft.GetSurface(0);
    EXPECT_EQ(wing.GetSection(0).GetNumControlSurfaces(), 2);
}

TEST_F(AvlControlSurfaceParsingTest, ParseControlSurfaceHingeVector) {
    std::string avl_content = R"(
Test Aircraft
0.0
0 0 0
10.0 1.0 10.0
0.25 0.0 0.0

SURFACE
Fin
4 1.0

SECTION
0.0 0.0 0.0 0.5 0.0
CONTROL
rudder 1.0 0.60 0 0 1 1

SECTION
0.0 0.0 1.0 0.5 0.0
CONTROL
rudder 1.0 0.60 0 0 1 1
)";

    AvlReader reader;
    Aircraft aircraft = reader.ParseString(avl_content);

    const Surface& fin = aircraft.GetSurface(0);
    const auto& controls = fin.GetSection(0).GetControlSurfaces();

    ASSERT_EQ(controls.size(), 1);
    EXPECT_EQ(controls[0].name, "rudder");

    // Hinge vector should be vertical for rudder
    EXPECT_NEAR(controls[0].hinge_vector.x(), 0.0, kTestTolerance);
    EXPECT_NEAR(controls[0].hinge_vector.y(), 0.0, kTestTolerance);
    EXPECT_NEAR(controls[0].hinge_vector.z(), 1.0, kTestTolerance);
}

TEST_F(AvlControlSurfaceParsingTest, ParseControlSurfaceGain) {
    std::string avl_content = R"(
Test Aircraft
0.0
0 0 0
10.0 1.0 10.0
0.25 0.0 0.0

SURFACE
Wing
4 1.0

SECTION
0.0 0.0 0.0 1.0 0.0
CONTROL
spoiler -0.5 0.50 0 1 0 1

SECTION
0.0 2.0 0.0 1.0 0.0
CONTROL
spoiler -0.5 0.50 0 1 0 1
)";

    AvlReader reader;
    Aircraft aircraft = reader.ParseString(avl_content);

    const Surface& wing = aircraft.GetSurface(0);
    const auto& controls = wing.GetSection(0).GetControlSurfaces();

    ASSERT_EQ(controls.size(), 1);
    EXPECT_EQ(controls[0].name, "spoiler");
    EXPECT_NEAR(controls[0].gain, -0.5, kTestTolerance);
}

// ============================================================================
// Control Surface Effectiveness Tests (Integration)
// ============================================================================

class ControlSurfaceEffectivenessTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a rectangular wing with elevator
        Surface wing("Wing");
        wing.SetChordwisePanelCount(8);
        wing.SetSpanwisePanelCount(12);
        wing.SetYDuplicate(0.0);  // Full span wing

        ControlSurfaceDefinition elevator;
        elevator.name = "elevator";
        elevator.hinge_x_fraction = 0.75;
        elevator.hinge_vector = Vec3(0, 1, 0);
        elevator.gain = 1.0;
        elevator.sign = 1.0;  // Symmetric deflection

        Section root(Vec3(0, 0, 0), 1.0, 0.0);
        root.AddControlSurface(elevator);

        Section tip(Vec3(0, 5.0, 0), 1.0, 0.0);
        tip.AddControlSurface(elevator);

        wing.AddSection(root);
        wing.AddSection(tip);

        aircraft_.AddSurface(wing);
        aircraft_.SetReferenceArea(10.0);  // Full span area
        aircraft_.SetReferenceChord(1.0);
        aircraft_.SetReferenceSpan(10.0);

        solver_.SetAircraft(aircraft_);
    }

    Aircraft aircraft_{"Control Test"};
    VortexLatticeSolver solver_;
};

TEST_F(ControlSurfaceEffectivenessTest, ControlSurfacesPropagateToMesh) {
    const VortexLattice& mesh = solver_.GetMesh();

    EXPECT_EQ(mesh.GetControlSurfaceNames().size(), 1);
    EXPECT_EQ(mesh.GetControlSurfaceNames()[0], "elevator");
}

TEST_F(ControlSurfaceEffectivenessTest, InfluenceMatrixBuiltWithControls) {
    // Solve to ensure AIC is built
    FlightCondition condition;
    condition.alpha = 5.0 * kDegreesToRadians;

    SolutionState solution = solver_.Solve(condition);

    EXPECT_TRUE(solution.is_valid);
    EXPECT_GT(solution.circulation.size(), 0);
}

// ============================================================================
// Control Derivatives Struct Tests
// ============================================================================

class ControlDerivativesTest : public ::testing::Test {};

TEST_F(ControlDerivativesTest, DefaultValues) {
    ControlDerivatives derivs;

    EXPECT_TRUE(derivs.name.empty());
    EXPECT_NEAR(derivs.deflection, 0.0, kTestTolerance);
    EXPECT_NEAR(derivs.gain, 1.0, kTestTolerance);

    // All derivatives should default to zero
    EXPECT_NEAR(derivs.CL_d, 0.0, kTestTolerance);
    EXPECT_NEAR(derivs.CY_d, 0.0, kTestTolerance);
    EXPECT_NEAR(derivs.CD_d, 0.0, kTestTolerance);
    EXPECT_NEAR(derivs.Cl_d, 0.0, kTestTolerance);
    EXPECT_NEAR(derivs.Cm_d, 0.0, kTestTolerance);
    EXPECT_NEAR(derivs.Cn_d, 0.0, kTestTolerance);
    EXPECT_NEAR(derivs.Ch_d, 0.0, kTestTolerance);
}

TEST_F(ControlDerivativesTest, SetValues) {
    ControlDerivatives derivs;
    derivs.name = "elevator";
    derivs.CL_d = 0.5;
    derivs.Cm_d = -1.2;

    EXPECT_EQ(derivs.name, "elevator");
    EXPECT_NEAR(derivs.CL_d, 0.5, kTestTolerance);
    EXPECT_NEAR(derivs.Cm_d, -1.2, kTestTolerance);
}

// ============================================================================
// FlightCondition Control Deflection Tests
// ============================================================================

class FlightConditionControlTest : public ::testing::Test {};

TEST_F(FlightConditionControlTest, DefaultNoDeflections) {
    FlightCondition condition;

    EXPECT_TRUE(condition.control_deflections.empty());
}

TEST_F(FlightConditionControlTest, SetControlDeflection) {
    FlightCondition condition;
    condition.control_deflections.push_back(
        {"elevator", 5.0 * kDegreesToRadians});
    condition.control_deflections.push_back(
        {"aileron", 10.0 * kDegreesToRadians});

    EXPECT_EQ(condition.control_deflections.size(), 2);
    EXPECT_EQ(condition.control_deflections[0].first, "elevator");
    EXPECT_NEAR(condition.control_deflections[0].second,
                5.0 * kDegreesToRadians, kTestTolerance);
}

// ============================================================================
// Stability Derivatives Control Storage Tests
// ============================================================================

class StabilityDerivativesControlTest : public ::testing::Test {};

TEST_F(StabilityDerivativesControlTest, ControlDerivativesVector) {
    StabilityDerivatives derivs;

    EXPECT_TRUE(derivs.control_derivatives.empty());

    ControlDerivatives elevator_derivs;
    elevator_derivs.name = "elevator";
    elevator_derivs.CL_d = 0.4;
    elevator_derivs.Cm_d = -1.0;

    derivs.control_derivatives.push_back(elevator_derivs);

    EXPECT_EQ(derivs.control_derivatives.size(), 1);
    EXPECT_EQ(derivs.control_derivatives[0].name, "elevator");
}

}  // namespace
}  // namespace pegasus
