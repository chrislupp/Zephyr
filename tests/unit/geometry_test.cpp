/*
    Project Pegasus
    Version 0.0.0

    Copyright 2024, Christopher A. Lupp

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

#include "pegasus/geometry/aircraft.h"
#include "pegasus/geometry/airfoil.h"
#include "pegasus/geometry/section.h"
#include "pegasus/geometry/surface.h"
#include "pegasus/io/avl_reader.h"

#include <cmath>

#include <gtest/gtest.h>

namespace pegasus {
namespace {

constexpr double kTestTolerance = 1.0e-6;

// ============================================================================
// Airfoil Tests
// ============================================================================

class AirfoilTest : public ::testing::Test {};

TEST_F(AirfoilTest, FlatPlateCreation) {
    Airfoil airfoil = Airfoil::CreateFlatPlate();

    EXPECT_TRUE(airfoil.IsValid());
    EXPECT_EQ(airfoil.GetName(), "Flat Plate");
    EXPECT_NEAR(airfoil.GetMaxCamber(), 0.0, kTestTolerance);
    EXPECT_NEAR(airfoil.GetMaxThickness(), 0.0, kTestTolerance);

    // Flat plate should have zero camber everywhere
    EXPECT_NEAR(airfoil.EvaluateCamber(0.0), 0.0, kTestTolerance);
    EXPECT_NEAR(airfoil.EvaluateCamber(0.5), 0.0, kTestTolerance);
    EXPECT_NEAR(airfoil.EvaluateCamber(1.0), 0.0, kTestTolerance);

    // Slope should be zero
    EXPECT_NEAR(airfoil.EvaluateCamberSlope(0.5), 0.0, kTestTolerance);
}

TEST_F(AirfoilTest, Naca0012Creation) {
    Airfoil airfoil = Airfoil::CreateFromNaca4Digit(12);

    EXPECT_TRUE(airfoil.IsValid());
    EXPECT_EQ(airfoil.GetName(), "NACA 12");
    EXPECT_NEAR(airfoil.GetMaxCamber(), 0.0, kTestTolerance);
    EXPECT_NEAR(airfoil.GetMaxThickness(), 0.12, kTestTolerance);

    // Symmetric airfoil - zero camber
    EXPECT_NEAR(airfoil.EvaluateCamber(0.3), 0.0, kTestTolerance);
    EXPECT_NEAR(airfoil.EvaluateCamber(0.5), 0.0, kTestTolerance);
}

TEST_F(AirfoilTest, Naca2412Creation) {
    Airfoil airfoil = Airfoil::CreateFromNaca4Digit(2412);

    EXPECT_TRUE(airfoil.IsValid());
    EXPECT_EQ(airfoil.GetName(), "NACA 2412");
    EXPECT_NEAR(airfoil.GetMaxCamber(), 0.02, kTestTolerance);
    EXPECT_NEAR(airfoil.GetMaxCamberPosition(), 0.4, kTestTolerance);
    EXPECT_NEAR(airfoil.GetMaxThickness(), 0.12, kTestTolerance);

    // Cambered airfoil should have positive camber
    Scalar camber_at_max = airfoil.EvaluateCamber(0.4);
    EXPECT_GT(camber_at_max, 0.0);
    EXPECT_NEAR(camber_at_max, 0.02, 0.001);  // Max camber at 40% chord
}

TEST_F(AirfoilTest, Naca4DigitFromString) {
    Airfoil airfoil = Airfoil::CreateFromNaca4Digit("4412");

    EXPECT_EQ(airfoil.GetName(), "NACA 4412");
    EXPECT_NEAR(airfoil.GetMaxCamber(), 0.04, kTestTolerance);
    EXPECT_NEAR(airfoil.GetMaxCamberPosition(), 0.4, kTestTolerance);
    EXPECT_NEAR(airfoil.GetMaxThickness(), 0.12, kTestTolerance);
}

TEST_F(AirfoilTest, CamberSlopeAtLeadingEdge) {
    // For NACA 2412, slope at LE should be positive (camber increasing)
    Airfoil airfoil = Airfoil::CreateFromNaca4Digit(2412);

    Scalar slope_le = airfoil.EvaluateCamberSlope(0.0);
    EXPECT_GT(slope_le, 0.0);  // Positive slope at LE

    Scalar slope_te = airfoil.EvaluateCamberSlope(1.0);
    EXPECT_LT(slope_te, 0.0);  // Negative slope at TE (camber decreasing)
}

// ============================================================================
// Section Tests
// ============================================================================

class SectionTest : public ::testing::Test {};

TEST_F(SectionTest, DefaultConstruction) {
    Section section;

    EXPECT_EQ(section.GetChord(), 1.0);
    EXPECT_EQ(section.GetTwistAngle(), 0.0);
    EXPECT_NEAR(section.GetClAlphaFactor(), 1.0, kTestTolerance);
    EXPECT_FALSE(section.HasAirfoil());
}

TEST_F(SectionTest, ParameterizedConstruction) {
    Vec3 le_pos(1.0, 2.0, 0.5);
    Scalar chord = 0.3;
    Scalar twist = 5.0;

    Section section(le_pos, chord, twist);

    EXPECT_EQ(section.GetLeadingEdgePosition(), le_pos);
    EXPECT_EQ(section.GetChord(), chord);
    EXPECT_EQ(section.GetTwistAngle(), twist);
    EXPECT_NEAR(section.GetTwistAngleRadians(), twist * kDegreesToRadians,
                kTestTolerance);
}

TEST_F(SectionTest, AirfoilAssignment) {
    Section section;
    auto airfoil =
        std::make_shared<Airfoil>(Airfoil::CreateFromNaca4Digit(2412));
    section.SetAirfoil(airfoil);

    EXPECT_TRUE(section.HasAirfoil());
    EXPECT_EQ(section.GetAirfoil().GetName(), "NACA 2412");
}

TEST_F(SectionTest, ControlSurface) {
    Section section;

    ControlSurfaceDefinition control;
    control.name = "aileron";
    control.hinge_x_fraction = 0.75;
    control.hinge_vector = Vec3(0, 1, 0);
    control.gain = 1.0;
    control.sign = 1.0;

    section.AddControlSurface(control);

    EXPECT_EQ(section.GetNumControlSurfaces(), 1);
    EXPECT_EQ(section.GetControlSurfaces()[0].name, "aileron");
}

// ============================================================================
// Surface Tests
// ============================================================================

class SurfaceTest : public ::testing::Test {};

TEST_F(SurfaceTest, BasicConstruction) {
    Surface surface("Main Wing");

    EXPECT_EQ(surface.GetName(), "Main Wing");
    EXPECT_EQ(surface.GetNumSections(), 0);
    EXPECT_TRUE(surface.HasWake());
    EXPECT_TRUE(surface.ContributesToLoads());
}

TEST_F(SurfaceTest, AddSections) {
    Surface surface("Wing");

    surface.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    surface.AddSection(Section(Vec3(0, 2.0, 0), 0.8, -2.0));

    EXPECT_EQ(surface.GetNumSections(), 2);
}

TEST_F(SurfaceTest, PlanformAreaRectangular) {
    Surface surface("Wing");

    // Create a rectangular wing: chord=1.0, span=4.0
    surface.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    surface.AddSection(Section(Vec3(0, 4.0, 0), 1.0, 0.0));

    Scalar area = surface.ComputePlanformArea();
    EXPECT_NEAR(area, 4.0, kTestTolerance);  // chord * span = 1.0 * 4.0
}

TEST_F(SurfaceTest, PlanformAreaTapered) {
    Surface surface("Wing");

    // Tapered wing: root chord=2.0, tip chord=1.0, span=4.0
    surface.AddSection(Section(Vec3(0, 0, 0), 2.0, 0.0));
    surface.AddSection(Section(Vec3(0, 4.0, 0), 1.0, 0.0));

    Scalar area = surface.ComputePlanformArea();
    // Trapezoidal area = 0.5 * (c_root + c_tip) * span = 0.5 * 3.0 * 4.0 = 6.0
    EXPECT_NEAR(area, 6.0, kTestTolerance);
}

TEST_F(SurfaceTest, YDuplicateDoublesArea) {
    Surface surface("Wing");

    surface.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    surface.AddSection(Section(Vec3(0, 2.0, 0), 1.0, 0.0));

    Scalar area_half = surface.ComputePlanformArea();
    EXPECT_NEAR(area_half, 2.0, kTestTolerance);

    surface.SetYDuplicate(0.0);
    Scalar area_full = surface.ComputePlanformArea();
    EXPECT_NEAR(area_full, 4.0, kTestTolerance);  // Doubled
}

// ============================================================================
// Aircraft Tests
// ============================================================================

class AircraftTest : public ::testing::Test {};

TEST_F(AircraftTest, DefaultValues) {
    Aircraft aircraft;

    EXPECT_EQ(aircraft.GetNumSurfaces(), 0);
    EXPECT_EQ(aircraft.GetNumBodies(), 0);
    EXPECT_NEAR(aircraft.GetMachNumber(), 0.0, kTestTolerance);
}

TEST_F(AircraftTest, AddSurface) {
    Aircraft aircraft("Test Aircraft");

    Surface wing("Wing");
    wing.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    wing.AddSection(Section(Vec3(0, 4.0, 0), 1.0, 0.0));

    aircraft.AddSurface(wing);

    EXPECT_EQ(aircraft.GetNumSurfaces(), 1);
    EXPECT_EQ(aircraft.GetTitle(), "Test Aircraft");
}

// ============================================================================
// AVL Reader Tests
// ============================================================================

class AvlReaderTest : public ::testing::Test {};

TEST_F(AvlReaderTest, ParseSimpleGeometry) {
    // Simple AVL content
    std::string avl_content = R"(
Test Wing
0.0
0 0 0
10.0 1.0 10.0
0.5 0.0 0.0
0.0

SURFACE
Wing
12 1.0 20 1.0

SECTION
0.0 0.0 0.0 1.0 0.0

SECTION
0.0 5.0 0.0 1.0 0.0
)";

    AvlReader reader;
    Aircraft aircraft = reader.ParseString(avl_content);

    EXPECT_EQ(aircraft.GetTitle(), "Test Wing");
    EXPECT_NEAR(aircraft.GetMachNumber(), 0.0, kTestTolerance);
    EXPECT_NEAR(aircraft.GetReferenceArea(), 10.0, kTestTolerance);
    EXPECT_NEAR(aircraft.GetReferenceChord(), 1.0, kTestTolerance);
    EXPECT_NEAR(aircraft.GetReferenceSpan(), 10.0, kTestTolerance);

    ASSERT_EQ(aircraft.GetNumSurfaces(), 1);
    const Surface& wing = aircraft.GetSurface(0);
    EXPECT_EQ(wing.GetName(), "Wing");
    EXPECT_EQ(wing.GetChordwisePanelCount(), 12);
    EXPECT_EQ(wing.GetSpanwisePanelCount(), 20);

    ASSERT_EQ(wing.GetNumSections(), 2);
}

TEST_F(AvlReaderTest, ParseWithNacaAirfoil) {
    std::string avl_content = R"(
NACA Test
0.1
0 0 0
1.0 0.25 4.0
0.0 0.0 0.0
0.0

SURFACE
Wing
8 1.0

SECTION
0.0 0.0 0.0 0.5 0.0
NACA
2412

SECTION
0.0 2.0 0.0 0.5 0.0
NACA
2412
)";

    AvlReader reader;
    Aircraft aircraft = reader.ParseString(avl_content);

    const Surface& wing = aircraft.GetSurface(0);
    ASSERT_EQ(wing.GetNumSections(), 2);

    const Section& root = wing.GetSection(0);
    EXPECT_TRUE(root.HasAirfoil());
    EXPECT_EQ(root.GetAirfoil().GetName(), "NACA 2412");
}

TEST_F(AvlReaderTest, ParseWithYDuplicate) {
    std::string avl_content = R"(
Symmetric Wing
0.0
0 0 0
2.0 0.5 4.0
0.0 0.0 0.0

SURFACE
Wing
6 1.0
YDUPLICATE
0.0

SECTION
0.0 0.0 0.0 0.5 0.0

SECTION
0.0 2.0 0.0 0.5 0.0
)";

    AvlReader reader;
    Aircraft aircraft = reader.ParseString(avl_content);

    const Surface& wing = aircraft.GetSurface(0);
    EXPECT_TRUE(wing.IsYDuplicated());
    EXPECT_NEAR(wing.GetYDuplicate().value(), 0.0, kTestTolerance);
}

TEST_F(AvlReaderTest, ParseWithControlSurface) {
    std::string avl_content = R"(
Control Test
0.0
0 0 0
1.0 0.25 2.0
0.0 0.0 0.0

SURFACE
Wing
6 1.0

SECTION
0.0 0.0 0.0 0.5 0.0
CONTROL
aileron 1.0 0.75 0 1 0 1.0

SECTION
0.0 1.0 0.0 0.5 0.0
CONTROL
aileron 1.0 0.75 0 1 0 -1.0
)";

    AvlReader reader;
    Aircraft aircraft = reader.ParseString(avl_content);

    auto control_names = aircraft.GetControlSurfaceNames();
    ASSERT_EQ(control_names.size(), 1);
    EXPECT_EQ(control_names[0], "aileron");

    const Surface& wing = aircraft.GetSurface(0);
    EXPECT_EQ(wing.GetSection(0).GetNumControlSurfaces(), 1);
    EXPECT_EQ(wing.GetSection(0).GetControlSurfaces()[0].name, "aileron");
}

}  // namespace
}  // namespace pegasus
