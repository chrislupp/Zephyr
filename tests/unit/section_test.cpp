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

#include "pegasus/geometry/section.h"

#include <memory>

#include <gtest/gtest.h>

namespace pegasus {
namespace {

constexpr double kTol = 1.0e-10;

// ============================================================================
// Section Construction
// ============================================================================

TEST(SectionUnitTest, DefaultConstructor) {
    Section s;
    EXPECT_NEAR(s.GetChord(), 1.0, kTol);
    EXPECT_NEAR(s.GetTwistAngle(), 0.0, kTol);
    EXPECT_FALSE(s.HasAirfoil());
}

TEST(SectionUnitTest, ParameterizedConstructor) {
    Vec3 pos(1.0, 2.0, 3.0);
    Section s(pos, 0.5, 5.0);
    EXPECT_NEAR(s.GetLeadingEdgePosition().x(), 1.0, kTol);
    EXPECT_NEAR(s.GetLeadingEdgePosition().y(), 2.0, kTol);
    EXPECT_NEAR(s.GetLeadingEdgePosition().z(), 3.0, kTol);
    EXPECT_NEAR(s.GetChord(), 0.5, kTol);
    EXPECT_NEAR(s.GetTwistAngle(), 5.0, kTol);
    EXPECT_FALSE(s.HasAirfoil());
}

TEST(SectionUnitTest, ConstructorWithAirfoil) {
    auto airfoil =
        std::make_shared<Airfoil>(Airfoil::CreateFromNaca4Digit("2412"));
    Section s(Vec3::Zero(), 1.0, 0.0, airfoil);
    EXPECT_TRUE(s.HasAirfoil());
}

// ============================================================================
// Airfoil Access
// ============================================================================

TEST(SectionUnitTest, GetAirfoilWithoutCustomReturnsFlat) {
    Section s(Vec3::Zero(), 1.0, 0.0);
    // Should return flat plate when no airfoil set
    const Airfoil& a = s.GetAirfoil();
    EXPECT_NEAR(a.EvaluateThickness(0.5), 0.0, kTol);
}

TEST(SectionUnitTest, GetAirfoilMutableCreatesFlat) {
    Section s(Vec3::Zero(), 1.0, 0.0);
    EXPECT_FALSE(s.HasAirfoil());
    // Mutable access should create a flat plate
    Airfoil& a = s.GetAirfoilMutable();
    EXPECT_TRUE(s.HasAirfoil());
    EXPECT_NEAR(a.EvaluateThickness(0.5), 0.0, kTol);
}

TEST(SectionUnitTest, GetAirfoilWithCustomReturnsCustom) {
    auto airfoil =
        std::make_shared<Airfoil>(Airfoil::CreateFromNaca4Digit("2412"));
    Section s(Vec3::Zero(), 1.0, 0.0, airfoil);
    const Airfoil& a = s.GetAirfoil();
    // NACA 2412 has non-zero thickness at 0.5
    EXPECT_GT(a.EvaluateThickness(0.5), 0.0);
}

TEST(SectionUnitTest, SetAirfoil) {
    Section s;
    EXPECT_FALSE(s.HasAirfoil());
    auto airfoil =
        std::make_shared<Airfoil>(Airfoil::CreateFromNaca4Digit("0012"));
    s.SetAirfoil(airfoil);
    EXPECT_TRUE(s.HasAirfoil());
}

// ============================================================================
// Setters and Getters
// ============================================================================

TEST(SectionUnitTest, SetLeadingEdgePosition) {
    Section s;
    Vec3 pos(5.0, 10.0, -1.0);
    s.SetLeadingEdgePosition(pos);
    EXPECT_NEAR(s.GetLeadingEdgePosition().x(), 5.0, kTol);
    EXPECT_NEAR(s.GetLeadingEdgePosition().y(), 10.0, kTol);
    EXPECT_NEAR(s.GetLeadingEdgePosition().z(), -1.0, kTol);
}

TEST(SectionUnitTest, SetChord) {
    Section s;
    s.SetChord(2.5);
    EXPECT_NEAR(s.GetChord(), 2.5, kTol);
}

TEST(SectionUnitTest, SetTwistAngle) {
    Section s;
    s.SetTwistAngle(3.0);
    EXPECT_NEAR(s.GetTwistAngle(), 3.0, kTol);
}

TEST(SectionUnitTest, TwistAngleRadians) {
    Section s;
    s.SetTwistAngle(90.0);
    EXPECT_NEAR(s.GetTwistAngleRadians(), kPi / 2.0, 1e-6);
}

TEST(SectionUnitTest, ClAlphaFactor) {
    Section s;
    EXPECT_NEAR(s.GetClAlphaFactor(), 1.0, kTol);
    s.SetClAlphaFactor(0.8);
    EXPECT_NEAR(s.GetClAlphaFactor(), 0.8, kTol);
}

TEST(SectionUnitTest, TwistOffsetDefaultsToZero) {
    Section s;
    EXPECT_NEAR(s.GetTwistOffset(), 0.0, kTol);
    EXPECT_NEAR(s.GetEffectiveTwistAngle(), s.GetTwistAngle(), kTol);
}

TEST(SectionUnitTest, TwistOffsetAddsToGeometricTwist) {
    Section s;
    s.SetTwistAngle(2.0);
    s.SetTwistOffset(0.5);
    EXPECT_NEAR(s.GetTwistAngle(), 2.0, kTol);  // geometric preserved
    EXPECT_NEAR(s.GetTwistOffset(), 0.5, kTol);
    EXPECT_NEAR(s.GetEffectiveTwistAngle(), 2.5, kTol);
    EXPECT_NEAR(s.GetEffectiveTwistAngleRadians(), 2.5 * kDegreesToRadians,
                1e-12);
}

// ============================================================================
// Control Surfaces
// ============================================================================

TEST(SectionUnitTest, AddControlSurface) {
    Section s;
    EXPECT_EQ(s.GetNumControlSurfaces(), 0u);

    ControlSurfaceDefinition ctrl;
    ctrl.name = "aileron";
    ctrl.hinge_x_fraction = 0.75;
    ctrl.hinge_vector = Vec3(0, 1, 0);
    ctrl.sign = 1.0;
    ctrl.gain = 1.0;
    s.AddControlSurface(ctrl);

    EXPECT_EQ(s.GetNumControlSurfaces(), 1u);
    EXPECT_EQ(s.GetControlSurfaces()[0].name, "aileron");
    EXPECT_NEAR(s.GetControlSurfaces()[0].hinge_x_fraction, 0.75, kTol);
}

// ============================================================================
// Design Variables
// ============================================================================

TEST(SectionUnitTest, AddDesignVariable) {
    Section s;
    EXPECT_TRUE(s.GetDesignVariables().empty());

    DesignVariableDefinition dv;
    dv.name = "twist";
    dv.weight = 1.5;
    s.AddDesignVariable(dv);

    EXPECT_EQ(s.GetDesignVariables().size(), 1u);
    EXPECT_EQ(s.GetDesignVariables()[0].name, "twist");
    EXPECT_NEAR(s.GetDesignVariables()[0].weight, 1.5, kTol);
}

// ============================================================================
// Paneling Overrides
// ============================================================================

TEST(SectionUnitTest, SpanwisePanelCount) {
    Section s;
    EXPECT_EQ(s.GetSpanwisePanelCount(), 0);
    s.SetSpanwisePanelCount(8);
    EXPECT_EQ(s.GetSpanwisePanelCount(), 8);
}

TEST(SectionUnitTest, SpanwiseSpacingEnum) {
    Section s;
    EXPECT_EQ(s.GetSpanwiseSpacing(), SpacingType::kNone);
    s.SetSpanwiseSpacing(SpacingType::kCosine);
    EXPECT_EQ(s.GetSpanwiseSpacing(), SpacingType::kCosine);
}

TEST(SectionUnitTest, SpanwiseSpacingNumericUniform) {
    Section s;
    s.SetSpanwiseSpacing(0.0);  // near-zero -> uniform
    EXPECT_EQ(s.GetSpanwiseSpacing(), SpacingType::kUniform);
}

TEST(SectionUnitTest, SpanwiseSpacingNumericCosine) {
    Section s;
    s.SetSpanwiseSpacing(1.0);  // positive -> cosine
    EXPECT_EQ(s.GetSpanwiseSpacing(), SpacingType::kCosine);
}

TEST(SectionUnitTest, SpanwiseSpacingNumericSine) {
    Section s;
    s.SetSpanwiseSpacing(-1.0);  // negative -> sine
    EXPECT_EQ(s.GetSpanwiseSpacing(), SpacingType::kSine);
}

// ============================================================================
// Viscous Polar
// ============================================================================

TEST(SectionUnitTest, ViscousPolar) {
    Section s;
    auto default_polar = s.GetViscousPolar();
    for (auto v : default_polar) {
        EXPECT_NEAR(v, 0.0, kTol);
    }

    std::array<Scalar, 6> polar = {0.01, 0.0, 0.005, 0.0, 0.0, 0.0};
    s.SetViscousPolar(polar);
    EXPECT_NEAR(s.GetViscousPolar()[0], 0.01, kTol);
    EXPECT_NEAR(s.GetViscousPolar()[2], 0.005, kTol);
}

}  // namespace
}  // namespace pegasus
