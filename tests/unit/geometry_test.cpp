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
// Airfoil Additional Tests
// ============================================================================

TEST(AirfoilExtraTest, InvalidNacaString) {
    EXPECT_THROW((void)Airfoil::CreateFromNaca4Digit("xyz"),
                 std::invalid_argument);
}

TEST(AirfoilExtraTest, NacaCodeOutOfRange) {
    EXPECT_THROW((void)Airfoil::CreateFromNaca4Digit(-1),
                 std::invalid_argument);
    EXPECT_THROW((void)Airfoil::CreateFromNaca4Digit(10000),
                 std::invalid_argument);
}

TEST(AirfoilExtraTest, NacaZeroThickness) {
    Airfoil airfoil = Airfoil::CreateFromNaca4Digit(0);
    EXPECT_TRUE(airfoil.IsValid());
    EXPECT_NEAR(airfoil.GetMaxThickness(), 0.0, kTestTolerance);
    EXPECT_NEAR(airfoil.GetMaxCamber(), 0.0, kTestTolerance);
}

TEST(AirfoilExtraTest, EvaluateThicknessNaca0012) {
    Airfoil airfoil = Airfoil::CreateFromNaca4Digit(12);
    Scalar t_mid = airfoil.EvaluateThickness(0.3);
    EXPECT_GT(t_mid, 0.0);
    EXPECT_NEAR(airfoil.EvaluateThickness(0.0), 0.0, 0.01);
}

TEST(AirfoilExtraTest, EvaluateCamberClamping) {
    Airfoil airfoil = Airfoil::CreateFromNaca4Digit(2412);
    // Clamping outside [0,1]
    Scalar below = airfoil.EvaluateCamber(-0.5);
    Scalar at_zero = airfoil.EvaluateCamber(0.0);
    EXPECT_NEAR(below, at_zero, kTestTolerance);

    Scalar above = airfoil.EvaluateCamber(1.5);
    Scalar at_one = airfoil.EvaluateCamber(1.0);
    EXPECT_NEAR(above, at_one, kTestTolerance);
}

TEST(AirfoilExtraTest, InvalidAirfoilReturnsZero) {
    Airfoil airfoil;  // default-constructed, not valid
    EXPECT_FALSE(airfoil.IsValid());
    EXPECT_NEAR(airfoil.EvaluateCamber(0.5), 0.0, kTestTolerance);
    EXPECT_NEAR(airfoil.EvaluateCamberSlope(0.5), 0.0, kTestTolerance);
    EXPECT_NEAR(airfoil.EvaluateThickness(0.5), 0.0, kTestTolerance);
}

TEST(AirfoilExtraTest, SetName) {
    Airfoil airfoil = Airfoil::CreateFlatPlate();
    airfoil.SetName("Custom");
    EXPECT_EQ(airfoil.GetName(), "Custom");
}

TEST(AirfoilExtraTest, CreateFromCoordinatesCamberLine) {
    // Simple triangular camber line
    std::vector<Scalar> x = {0.0, 0.5, 1.0};
    std::vector<Scalar> y = {0.0, 0.05, 0.0};

    Airfoil airfoil = Airfoil::CreateFromCoordinates(x, y, true);
    EXPECT_TRUE(airfoil.IsValid());
    EXPECT_GT(airfoil.GetMaxCamber(), 0.0);
    EXPECT_NEAR(airfoil.GetMaxCamberPosition(), 0.5, kTestTolerance);
}

TEST(AirfoilExtraTest, CreateFromCoordinatesMismatch) {
    std::vector<Scalar> x = {0.0, 0.5, 1.0};
    std::vector<Scalar> y = {0.0, 0.05};
    EXPECT_THROW((void)Airfoil::CreateFromCoordinates(x, y, true),
                 std::invalid_argument);
}

TEST(AirfoilExtraTest, CreateFromCoordinatesTooFew) {
    std::vector<Scalar> x = {0.5};
    std::vector<Scalar> y = {0.0};
    EXPECT_THROW((void)Airfoil::CreateFromCoordinates(x, y, true),
                 std::invalid_argument);
}

TEST(AirfoilExtraTest, CreateFromCoordinatesFullAirfoil) {
    // Simple symmetric diamond airfoil (Selig order: TE upper -> LE -> TE
    // lower)
    std::vector<Scalar> x = {1.0, 0.5, 0.0, 0.5, 1.0};
    std::vector<Scalar> y = {0.0, 0.05, 0.0, -0.05, 0.0};

    Airfoil airfoil = Airfoil::CreateFromCoordinates(x, y, false);
    EXPECT_TRUE(airfoil.IsValid());
    EXPECT_GT(airfoil.GetMaxThickness(), 0.0);
    EXPECT_NEAR(airfoil.GetMaxCamber(), 0.0, 0.01);
}

TEST(AirfoilExtraTest, CreateFromCoordinateFileNotFound) {
    EXPECT_THROW((void)Airfoil::CreateFromCoordinateFile("nonexistent.dat"),
                 std::runtime_error);
}

TEST(AirfoilExtraTest, CamberXY) {
    Airfoil airfoil = Airfoil::CreateFromNaca4Digit(2412);
    const auto& cx = airfoil.GetCamberX();
    const auto& cy = airfoil.GetCamberY();
    EXPECT_EQ(cx.size(), cy.size());
    EXPECT_GT(cx.size(), 2u);
    EXPECT_NEAR(cx.front(), 0.0, kTestTolerance);
    EXPECT_NEAR(cx.back(), 1.0, kTestTolerance);
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
// Surface Additional Tests
// ============================================================================

TEST(SurfaceExtraTest, MeanAerodynamicChordRectangular) {
    Surface surface("Wing");
    surface.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    surface.AddSection(Section(Vec3(0, 4.0, 0), 1.0, 0.0));

    Scalar mac = surface.ComputeMeanAerodynamicChord();
    EXPECT_NEAR(mac, 1.0, kTestTolerance);
}

TEST(SurfaceExtraTest, MeanAerodynamicChordTapered) {
    Surface surface("Wing");
    surface.AddSection(Section(Vec3(0, 0, 0), 2.0, 0.0));
    surface.AddSection(Section(Vec3(0, 4.0, 0), 1.0, 0.0));

    Scalar mac = surface.ComputeMeanAerodynamicChord();
    // MAC = integral(c^2 dy) / integral(c dy)
    // For linear taper: MAC > avg_chord
    EXPECT_GT(mac, 1.5);
    EXPECT_LT(mac, 2.0);
}

TEST(SurfaceExtraTest, MeanAerodynamicChordSingleSection) {
    Surface surface("Wing");
    surface.AddSection(Section(Vec3(0, 0, 0), 1.5, 0.0));

    Scalar mac = surface.ComputeMeanAerodynamicChord();
    EXPECT_NEAR(mac, 1.5, kTestTolerance);
}

TEST(SurfaceExtraTest, MeanAerodynamicChordEmpty) {
    Surface surface("Wing");
    EXPECT_NEAR(surface.ComputeMeanAerodynamicChord(), 0.0, kTestTolerance);
}

TEST(SurfaceExtraTest, ComputeSpan) {
    Surface surface("Wing");
    surface.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    surface.AddSection(Section(Vec3(0, 5.0, 0), 1.0, 0.0));

    EXPECT_NEAR(surface.ComputeSpan(), 5.0, kTestTolerance);
}

TEST(SurfaceExtraTest, ComputeSpanYDuplicate) {
    Surface surface("Wing");
    surface.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    surface.AddSection(Section(Vec3(0, 5.0, 0), 1.0, 0.0));
    surface.SetYDuplicate(0.0);

    EXPECT_NEAR(surface.ComputeSpan(), 10.0, kTestTolerance);
}

TEST(SurfaceExtraTest, ComputeSpanEmpty) {
    Surface surface("Wing");
    EXPECT_NEAR(surface.ComputeSpan(), 0.0, kTestTolerance);
}

TEST(SurfaceExtraTest, PlanformAreaSingleSection) {
    Surface surface("Wing");
    surface.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));

    EXPECT_NEAR(surface.ComputePlanformArea(), 0.0, kTestTolerance);
}

TEST(SurfaceExtraTest, HasViscousDragDefault) {
    Surface surface("Wing");
    EXPECT_FALSE(surface.HasViscousDrag());
}

TEST(SurfaceExtraTest, HasViscousDragNonZero) {
    Surface surface("Wing");
    std::array<Scalar, 6> polar = {0.01, 0.0, 0.1, 0.0, 0.0, 0.0};
    surface.SetClCdPolar(polar);
    EXPECT_TRUE(surface.HasViscousDrag());

    auto retrieved = surface.GetClCdPolar();
    EXPECT_NEAR(retrieved[0], 0.01, kTestTolerance);
}

TEST(SurfaceExtraTest, ChordwiseSpacingNumeric) {
    Surface surface("Wing");

    surface.SetChordwiseSpacing(0.0);
    EXPECT_EQ(surface.GetChordwiseSpacing(), SpacingType::kUniform);

    surface.SetChordwiseSpacing(1.0);
    EXPECT_EQ(surface.GetChordwiseSpacing(), SpacingType::kCosine);

    surface.SetChordwiseSpacing(-1.0);
    EXPECT_EQ(surface.GetChordwiseSpacing(), SpacingType::kSine);
}

TEST(SurfaceExtraTest, SpanwiseSpacingNumeric) {
    Surface surface("Wing");

    surface.SetSpanwiseSpacing(0.0);
    EXPECT_EQ(surface.GetSpanwiseSpacing(), SpacingType::kUniform);

    surface.SetSpanwiseSpacing(1.0);
    EXPECT_EQ(surface.GetSpanwiseSpacing(), SpacingType::kCosine);

    surface.SetSpanwiseSpacing(-1.0);
    EXPECT_EQ(surface.GetSpanwiseSpacing(), SpacingType::kSine);
}

TEST(SurfaceExtraTest, ClearYDuplicate) {
    Surface surface("Wing");
    surface.SetYDuplicate(0.0);
    EXPECT_TRUE(surface.IsYDuplicated());

    surface.ClearYDuplicate();
    EXPECT_FALSE(surface.IsYDuplicated());
}

TEST(SurfaceExtraTest, AerodynamicOptions) {
    Surface surface("Wing");

    surface.SetHasWake(false);
    EXPECT_FALSE(surface.HasWake());

    surface.SetContributesToLoads(false);
    EXPECT_FALSE(surface.ContributesToLoads());

    surface.SetAffectedByAlphaBeta(false);
    EXPECT_FALSE(surface.IsAffectedByAlphaBeta());
}

TEST(SurfaceExtraTest, CoordinateTransforms) {
    Surface surface("Wing");

    surface.SetAngleOffset(5.0);
    EXPECT_NEAR(surface.GetAngleOffset(), 5.0, kTestTolerance);

    Vec3 trans(1.0, 2.0, 3.0);
    surface.SetTranslation(trans);
    EXPECT_EQ(surface.GetTranslation(), trans);

    Vec3 scale(2.0, 2.0, 2.0);
    surface.SetScale(scale);
    EXPECT_EQ(surface.GetScale(), scale);
}

TEST(SurfaceExtraTest, ComponentIndex) {
    Surface surface("Wing");
    surface.SetComponentIndex(3);
    EXPECT_EQ(surface.GetComponentIndex(), 3);
}

TEST(SurfaceExtraTest, SectionAccessors) {
    Surface surface("Wing");
    surface.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    surface.AddSection(Section(Vec3(0, 5.0, 0), 0.5, 0.0));

    const auto& sections = surface.GetSections();
    EXPECT_EQ(sections.size(), 2u);

    EXPECT_NEAR(surface.GetSection(0).GetChord(), 1.0, kTestTolerance);
    EXPECT_NEAR(surface.GetSection(1).GetChord(), 0.5, kTestTolerance);

    surface.GetSectionMutable(0).SetChord(0.8);
    EXPECT_NEAR(surface.GetSection(0).GetChord(), 0.8, kTestTolerance);

    auto& mut_sections = surface.GetSectionsMutable();
    EXPECT_EQ(mut_sections.size(), 2u);
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

TEST_F(AircraftTest, FindSurfaceByName) {
    Aircraft aircraft;
    Surface wing("Wing");
    wing.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    Surface tail("Tail");
    tail.AddSection(Section(Vec3(3, 0, 0), 0.5, 0.0));

    aircraft.AddSurface(wing);
    aircraft.AddSurface(tail);

    const Surface* found = aircraft.FindSurfaceByName("Tail");
    ASSERT_NE(found, nullptr);
    EXPECT_EQ(found->GetName(), "Tail");

    EXPECT_EQ(aircraft.FindSurfaceByName("Missing"), nullptr);
}

TEST_F(AircraftTest, ReferenceQuantities) {
    Aircraft aircraft;

    aircraft.SetReferenceArea(10.0);
    aircraft.SetReferenceChord(1.5);
    aircraft.SetReferenceSpan(8.0);
    Vec3 mrp(0.25, 0.0, 0.0);
    aircraft.SetMomentReferencePoint(mrp);

    EXPECT_NEAR(aircraft.GetReferenceArea(), 10.0, kTestTolerance);
    EXPECT_NEAR(aircraft.GetReferenceChord(), 1.5, kTestTolerance);
    EXPECT_NEAR(aircraft.GetReferenceSpan(), 8.0, kTestTolerance);
    EXPECT_EQ(aircraft.GetMomentReferencePoint(), mrp);
}

TEST_F(AircraftTest, FlowConditions) {
    Aircraft aircraft;

    aircraft.SetMachNumber(0.3);
    EXPECT_NEAR(aircraft.GetMachNumber(), 0.3, kTestTolerance);

    aircraft.SetProfileDragCoefficient(0.01);
    EXPECT_NEAR(aircraft.GetProfileDragCoefficient(), 0.01, kTestTolerance);
}

TEST_F(AircraftTest, Symmetry) {
    Aircraft aircraft;

    aircraft.SetYSymmetry(SymmetryType::kSymmetric);
    EXPECT_EQ(aircraft.GetYSymmetry(), SymmetryType::kSymmetric);

    aircraft.SetZSymmetry(SymmetryType::kAntisymmetric);
    EXPECT_EQ(aircraft.GetZSymmetry(), SymmetryType::kAntisymmetric);

    aircraft.SetZSymmetryPlane(1.0);
    EXPECT_NEAR(aircraft.GetZSymmetryPlane(), 1.0, kTestTolerance);
}

TEST_F(AircraftTest, Bodies) {
    Aircraft aircraft;
    Body body("Fuselage");
    aircraft.AddBody(body);

    EXPECT_EQ(aircraft.GetNumBodies(), 1);
    EXPECT_EQ(aircraft.GetBodies().size(), 1u);
    EXPECT_EQ(aircraft.GetBodies()[0].GetName(), "Fuselage");
}

TEST_F(AircraftTest, SurfaceAccessors) {
    Aircraft aircraft;
    Surface wing("Wing");
    wing.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    aircraft.AddSurface(wing);

    const auto& surfaces = aircraft.GetSurfaces();
    EXPECT_EQ(surfaces.size(), 1u);

    EXPECT_EQ(aircraft.GetSurface(0).GetName(), "Wing");

    auto& mut = aircraft.GetSurfacesMutable();
    mut[0].SetName("Modified Wing");
    EXPECT_EQ(aircraft.GetSurface(0).GetName(), "Modified Wing");
}

TEST_F(AircraftTest, MoveSurface) {
    Aircraft aircraft;
    Surface wing("Wing");
    wing.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    aircraft.AddSurface(std::move(wing));

    EXPECT_EQ(aircraft.GetNumSurfaces(), 1);
    EXPECT_EQ(aircraft.GetSurface(0).GetName(), "Wing");
}

TEST_F(AircraftTest, ControlSurfaceNames) {
    Aircraft aircraft;
    Surface wing("Wing");

    Section root(Vec3(0, 0, 0), 1.0, 0.0);
    ControlSurfaceDefinition aileron;
    aileron.name = "aileron";
    aileron.hinge_x_fraction = 0.75;
    aileron.hinge_vector = Vec3(0, 1, 0);
    aileron.gain = 1.0;
    aileron.sign = 1.0;
    root.AddControlSurface(aileron);
    wing.AddSection(root);

    Section tip(Vec3(0, 4.0, 0), 0.8, 0.0);
    tip.AddControlSurface(aileron);
    wing.AddSection(tip);

    aircraft.AddSurface(wing);

    auto names = aircraft.GetControlSurfaceNames();
    EXPECT_EQ(names.size(), 1u);
    EXPECT_EQ(names[0], "aileron");
    EXPECT_EQ(aircraft.GetNumControlSurfaces(), 1u);
}

TEST_F(AircraftTest, DesignVariableNames) {
    Aircraft aircraft;
    Surface wing("Wing");

    Section root(Vec3(0, 0, 0), 1.0, 0.0);
    DesignVariableDefinition dvar;
    dvar.name = "twist";
    dvar.weight = 1.0;
    root.AddDesignVariable(dvar);
    wing.AddSection(root);

    aircraft.AddSurface(wing);

    auto names = aircraft.GetDesignVariableNames();
    EXPECT_EQ(names.size(), 1u);
    EXPECT_EQ(names[0], "twist");
    EXPECT_EQ(aircraft.GetNumDesignVariables(), 1u);
}

TEST_F(AircraftTest, ComputeTotalVortexCount) {
    Aircraft aircraft;
    Surface wing("Wing");
    wing.SetChordwisePanelCount(4);
    wing.SetSpanwisePanelCount(6);
    wing.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    wing.AddSection(Section(Vec3(0, 5.0, 0), 1.0, 0.0));
    aircraft.AddSurface(wing);

    // 4 * 6 = 24 vortices
    EXPECT_EQ(aircraft.ComputeTotalVortexCount(), 24);
}

TEST_F(AircraftTest, ComputeTotalVortexCountYDuplicate) {
    Aircraft aircraft;
    Surface wing("Wing");
    wing.SetChordwisePanelCount(4);
    wing.SetSpanwisePanelCount(6);
    wing.SetYDuplicate(0.0);
    wing.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    wing.AddSection(Section(Vec3(0, 5.0, 0), 1.0, 0.0));
    aircraft.AddSurface(wing);

    // 4 * 6 * 2 = 48
    EXPECT_EQ(aircraft.ComputeTotalVortexCount(), 48);
}

TEST_F(AircraftTest, ComputeTotalStripCount) {
    Aircraft aircraft;
    Surface wing("Wing");
    wing.SetSpanwisePanelCount(8);
    wing.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    wing.AddSection(Section(Vec3(0, 5.0, 0), 1.0, 0.0));
    aircraft.AddSurface(wing);

    EXPECT_EQ(aircraft.ComputeTotalStripCount(), 8);
}

TEST_F(AircraftTest, ComputeTotalStripCountYDuplicate) {
    Aircraft aircraft;
    Surface wing("Wing");
    wing.SetSpanwisePanelCount(8);
    wing.SetYDuplicate(0.0);
    wing.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    wing.AddSection(Section(Vec3(0, 5.0, 0), 1.0, 0.0));
    aircraft.AddSurface(wing);

    EXPECT_EQ(aircraft.ComputeTotalStripCount(), 16);
}

TEST_F(AircraftTest, ComputeTotalPlanformArea) {
    Aircraft aircraft;
    Surface wing("Wing");
    wing.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    wing.AddSection(Section(Vec3(0, 5.0, 0), 1.0, 0.0));

    Surface tail("Tail");
    tail.AddSection(Section(Vec3(4, 0, 0), 0.5, 0.0));
    tail.AddSection(Section(Vec3(4, 1.5, 0), 0.3, 0.0));

    aircraft.AddSurface(wing);
    aircraft.AddSurface(tail);

    // Wing: 5.0, Tail: 0.5*(0.5+0.3)*1.5 = 0.6
    Scalar area = aircraft.ComputeTotalPlanformArea();
    EXPECT_NEAR(area, 5.6, kTestTolerance);
}

TEST_F(AircraftTest, ComputeTotalSpan) {
    Aircraft aircraft;
    Surface wing("Wing");
    wing.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    wing.AddSection(Section(Vec3(0, 5.0, 0), 1.0, 0.0));

    Surface tail("Tail");
    tail.AddSection(Section(Vec3(4, 0, 0), 0.5, 0.0));
    tail.AddSection(Section(Vec3(4, 1.5, 0), 0.3, 0.0));

    aircraft.AddSurface(wing);
    aircraft.AddSurface(tail);

    EXPECT_NEAR(aircraft.ComputeTotalSpan(), 5.0, kTestTolerance);
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

// ============================================================================
// AVL Reader Additional Tests
// ============================================================================

TEST(AvlReaderExtraTest, ReadFileNotFound) {
    AvlReader reader;
    EXPECT_THROW((void)reader.ReadFile("nonexistent.avl"), std::runtime_error);
}

TEST(AvlReaderExtraTest, ParseEmptyContent) {
    AvlReader reader;
    EXPECT_THROW((void)reader.ParseString(""), std::runtime_error);
}

TEST(AvlReaderExtraTest, ParseWithComments) {
    std::string content = R"(
# This is a comment
! Another comment
Test Wing
0.0
0 0 0
10.0 1.0 10.0
0.0 0.0 0.0
0.0

SURFACE
Wing
6 1.0

SECTION
0.0 0.0 0.0 1.0 0.0

SECTION
0.0 5.0 0.0 1.0 0.0
)";
    AvlReader reader;
    Aircraft aircraft = reader.ParseString(content);
    EXPECT_EQ(aircraft.GetTitle(), "Test Wing");
}

TEST(AvlReaderExtraTest, ParseWithProfileDrag) {
    std::string content = R"(
Wing with drag
0.0
0 0 0
10.0 1.0 10.0
0.0 0.0 0.0
0.02

SURFACE
Wing
6 1.0

SECTION
0.0 0.0 0.0 1.0 0.0

SECTION
0.0 5.0 0.0 1.0 0.0
)";
    AvlReader reader;
    Aircraft aircraft = reader.ParseString(content);
    EXPECT_NEAR(aircraft.GetProfileDragCoefficient(), 0.02, 1e-6);
}

TEST(AvlReaderExtraTest, ParseWithSymmetry) {
    std::string content = R"(
Symmetric Wing
0.0
1 -1 0.5
10.0 1.0 10.0
0.0 0.0 0.0
0.0

SURFACE
Wing
6 1.0

SECTION
0.0 0.0 0.0 1.0 0.0

SECTION
0.0 5.0 0.0 1.0 0.0
)";
    AvlReader reader;
    Aircraft aircraft = reader.ParseString(content);
    EXPECT_EQ(aircraft.GetYSymmetry(), SymmetryType::kSymmetric);
    EXPECT_EQ(aircraft.GetZSymmetry(), SymmetryType::kAntisymmetric);
    EXPECT_NEAR(aircraft.GetZSymmetryPlane(), 0.5, 1e-6);
}

TEST(AvlReaderExtraTest, ParseSurfaceModifiers) {
    std::string content = R"(
Modifiers Test
0.0
0 0 0
10.0 1.0 10.0
0.0 0.0 0.0
0.0

SURFACE
Wing
6 1.0 8 1.0
SCALE
2.0 2.0 2.0
TRANSLATE
1.0 0.0 0.0
ANGLE
5.0
COMPONENT
2
NOWAKE
NOLOAD
NOALBE
CDCL
0.01 0.0 0.005 0.0 0.0 0.0

SECTION
0.0 0.0 0.0 1.0 0.0

SECTION
0.0 5.0 0.0 1.0 0.0
)";
    AvlReader reader;
    Aircraft aircraft = reader.ParseString(content);

    const Surface& wing = aircraft.GetSurface(0);
    EXPECT_EQ(wing.GetScale(), Vec3(2.0, 2.0, 2.0));
    EXPECT_EQ(wing.GetTranslation(), Vec3(1.0, 0.0, 0.0));
    EXPECT_NEAR(wing.GetAngleOffset(), 5.0, 1e-6);
    EXPECT_EQ(wing.GetComponentIndex(), 2);
    EXPECT_FALSE(wing.HasWake());
    EXPECT_FALSE(wing.ContributesToLoads());
    EXPECT_FALSE(wing.IsAffectedByAlphaBeta());
    EXPECT_TRUE(wing.HasViscousDrag());
    EXPECT_EQ(wing.GetSpanwisePanelCount(), 8);
}

TEST(AvlReaderExtraTest, ParseSectionSpanwise) {
    std::string content = R"(
Section Spanwise
0.0
0 0 0
10.0 1.0 10.0
0.0 0.0 0.0
0.0

SURFACE
Wing
6 1.0

SECTION
0.0 0.0 0.0 1.0 3.0 4 -1.0
CLAF
1.5

SECTION
0.0 5.0 0.0 0.8 -2.0 4 1.0
)";
    AvlReader reader;
    Aircraft aircraft = reader.ParseString(content);

    const Surface& wing = aircraft.GetSurface(0);
    EXPECT_EQ(wing.GetNumSections(), 2);
    EXPECT_NEAR(wing.GetSection(0).GetTwistAngle(), 3.0, 1e-6);
    EXPECT_NEAR(wing.GetSection(0).GetClAlphaFactor(), 1.5, 1e-6);
    EXPECT_NEAR(wing.GetSection(1).GetChord(), 0.8, 1e-6);
}

TEST(AvlReaderExtraTest, ParseWithDesignVariable) {
    std::string content = R"(
Design Test
0.0
0 0 0
10.0 1.0 10.0
0.0 0.0 0.0
0.0

SURFACE
Wing
6 1.0

SECTION
0.0 0.0 0.0 1.0 0.0
DESIGN
twist 1.0

SECTION
0.0 5.0 0.0 1.0 0.0
DESIGN
twist 0.5
)";
    AvlReader reader;
    Aircraft aircraft = reader.ParseString(content);

    auto dv_names = aircraft.GetDesignVariableNames();
    EXPECT_EQ(dv_names.size(), 1u);
    EXPECT_EQ(dv_names[0], "twist");
}

TEST(AvlReaderExtraTest, ParseWithBody) {
    std::string content = R"(
Body Test
0.0
0 0 0
10.0 1.0 10.0
0.0 0.0 0.0
0.0

BODY
Fuselage
20 1.0
SCALE
1.0 1.0 1.0
TRANSLATE
0.0 0.0 0.0
YDUPLICATE
BSEC
0.0 0.0 0.0 0.2 0.2
BSEC
2.0 0.0 0.0 0.5 0.5
BSEC
5.0 0.0 0.0 0.1 0.1

SURFACE
Wing
6 1.0

SECTION
0.0 0.0 0.0 1.0 0.0

SECTION
0.0 5.0 0.0 1.0 0.0
)";
    AvlReader reader;
    Aircraft aircraft = reader.ParseString(content);

    EXPECT_EQ(aircraft.GetNumBodies(), 1);
    EXPECT_EQ(aircraft.GetBodies()[0].GetName(), "Fuselage");
    EXPECT_TRUE(aircraft.GetBodies()[0].IsYDuplicated());
    EXPECT_EQ(aircraft.GetNumSurfaces(), 1);
}

TEST(AvlReaderExtraTest, ParseBodyAtEnd) {
    std::string content = R"(
Body at End
0.0
0 0 0
10.0 1.0 10.0
0.0 0.0 0.0
0.0

BODY
Fuselage
20 1.0
BSEC
0.0 0.0 0.0 0.2 0.2
)";
    AvlReader reader;
    Aircraft aircraft = reader.ParseString(content);

    EXPECT_EQ(aircraft.GetNumBodies(), 1);
    EXPECT_EQ(aircraft.GetBodies()[0].GetName(), "Fuselage");
}

TEST(AvlReaderExtraTest, ParseVerbose) {
    std::string content = R"(
Verbose Test
0.0
0 0 0
10.0 1.0 10.0
0.0 0.0 0.0
0.0

SURFACE
Wing
6 1.0

SECTION
0.0 0.0 0.0 1.0 0.0

SECTION
0.0 5.0 0.0 1.0 0.0
)";
    AvlReader reader;
    reader.SetVerbose(true);
    Aircraft aircraft = reader.ParseString(content);
    EXPECT_EQ(aircraft.GetTitle(), "Verbose Test");
}

TEST(AvlReaderExtraTest, ParseMomentReference) {
    std::string content = R"(
MomRef Test
0.5
0 0 0
10.0 1.0 10.0
0.25 0.0 0.1
0.0

SURFACE
Wing
6 1.0

SECTION
0.0 0.0 0.0 1.0 0.0

SECTION
0.0 5.0 0.0 1.0 0.0
)";
    AvlReader reader;
    Aircraft aircraft = reader.ParseString(content);

    EXPECT_NEAR(aircraft.GetMachNumber(), 0.5, 1e-6);
    Vec3 mrp = aircraft.GetMomentReferencePoint();
    EXPECT_NEAR(mrp.x(), 0.25, 1e-6);
    EXPECT_NEAR(mrp.z(), 0.1, 1e-6);
}

TEST(AvlReaderExtraTest, ParseMultipleSurfaces) {
    std::string content = R"(
Multi Surface
0.0
0 0 0
10.0 1.0 10.0
0.0 0.0 0.0
0.0

SURFACE
Wing
6 1.0

SECTION
0.0 0.0 0.0 1.0 0.0

SECTION
0.0 5.0 0.0 1.0 0.0

SURFACE
Tail
4 1.0

SECTION
4.0 0.0 0.0 0.5 0.0

SECTION
4.0 1.5 0.0 0.3 0.0
)";
    AvlReader reader;
    Aircraft aircraft = reader.ParseString(content);

    EXPECT_EQ(aircraft.GetNumSurfaces(), 2);
    EXPECT_EQ(aircraft.GetSurface(0).GetName(), "Wing");
    EXPECT_EQ(aircraft.GetSurface(1).GetName(), "Tail");
}

}  // namespace
}  // namespace pegasus
