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

#include "pegasus/aero/stability_derivatives.h"

#include "pegasus/geometry/aircraft.h"
#include "pegasus/geometry/section.h"
#include "pegasus/geometry/surface.h"
#include "pegasus/solver/vortex_lattice_solver.h"

#include <cmath>

#include <gtest/gtest.h>

namespace pegasus {
namespace {

constexpr double kTestTolerance = 1.0e-6;
constexpr double kDerivTolerance = 0.15;  // 15% tolerance for derivatives

// ============================================================================
// StabilityDerivatives Tests
// ============================================================================

class StabilityDerivativesTest : public ::testing::Test {
protected:
    Aircraft CreateRectangularWing(Scalar span, Scalar chord, int num_chordwise,
                                   int num_spanwise) {
        Aircraft aircraft("Test Wing");

        Surface wing("Wing");
        wing.SetChordwisePanelCount(num_chordwise);
        wing.SetSpanwisePanelCount(num_spanwise);
        wing.SetChordwiseSpacing(SpacingType::kCosine);
        wing.SetSpanwiseSpacing(SpacingType::kCosine);

        // Full wing
        wing.AddSection(Section(Vec3(0, -span / 2, 0), chord, 0.0));
        wing.AddSection(Section(Vec3(0, span / 2, 0), chord, 0.0));

        aircraft.AddSurface(wing);
        aircraft.SetReferenceArea(span * chord);
        aircraft.SetReferenceChord(chord);
        aircraft.SetReferenceSpan(span);

        return aircraft;
    }

    Aircraft CreateTaperedWing(Scalar span, Scalar root_chord, Scalar tip_chord,
                               int num_chordwise, int num_spanwise) {
        Aircraft aircraft("Tapered Wing");

        Surface wing("Wing");
        wing.SetChordwisePanelCount(num_chordwise);
        wing.SetSpanwisePanelCount(num_spanwise);
        wing.SetChordwiseSpacing(SpacingType::kCosine);
        wing.SetSpanwiseSpacing(SpacingType::kCosine);

        wing.AddSection(Section(Vec3(0, -span / 2, 0), tip_chord, 0.0));
        wing.AddSection(Section(Vec3(0, 0, 0), root_chord, 0.0));
        wing.AddSection(Section(Vec3(0, span / 2, 0), tip_chord, 0.0));

        Scalar avg_chord = (root_chord + tip_chord) / 2.0;
        aircraft.AddSurface(wing);
        aircraft.SetReferenceArea(span * avg_chord);
        aircraft.SetReferenceChord(avg_chord);
        aircraft.SetReferenceSpan(span);

        return aircraft;
    }
};

TEST_F(StabilityDerivativesTest, LiftCurveSlopePositive) {
    // Lift curve slope should be positive for any wing
    Aircraft aircraft = CreateRectangularWing(10.0, 1.0, 4, 16);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition condition;
    condition.alpha = 5.0 * kDegreesToRadians;

    solver.Solve(condition);

    StabilityDerivativeCalculator calculator;
    StabilityDerivatives derivs =
        calculator.ComputeDerivatives(aircraft, solver, condition);

    EXPECT_TRUE(derivs.is_valid);
    EXPECT_GT(derivs.CL_alpha, 0.0);
}

TEST_F(StabilityDerivativesTest, LiftCurveSlopeReasonable) {
    // For finite wing: CL_alpha ≈ 2*pi / (1 + 2/AR) (Helmbold approximation)
    // For AR=10: CL_alpha ≈ 5.24 rad^-1
    Scalar span = 10.0;
    Scalar chord = 1.0;
    Aircraft aircraft = CreateRectangularWing(span, chord, 4, 16);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition condition;
    condition.alpha = 5.0 * kDegreesToRadians;

    solver.Solve(condition);

    StabilityDerivativeCalculator calculator;
    StabilityDerivatives derivs =
        calculator.ComputeDerivatives(aircraft, solver, condition);

    Scalar AR = span / chord;
    Scalar CL_alpha_theory = 2.0 * kPi / (1.0 + 2.0 / AR);

    EXPECT_TRUE(derivs.is_valid);
    EXPECT_NEAR(derivs.CL_alpha, CL_alpha_theory, 0.3 * CL_alpha_theory);
}

TEST_F(StabilityDerivativesTest, ZeroAlphaSymmetry) {
    // At zero alpha, symmetric wing should have zero sideslip-related
    // derivatives
    Aircraft aircraft = CreateRectangularWing(10.0, 1.0, 4, 16);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition condition;
    condition.alpha = 0.0;
    condition.beta = 0.0;

    solver.Solve(condition);

    StabilityDerivativeCalculator calculator;
    StabilityDerivatives derivs =
        calculator.ComputeDerivatives(aircraft, solver, condition);

    EXPECT_TRUE(derivs.is_valid);
    // Side force derivative w.r.t. alpha should be small for symmetric wing
    EXPECT_NEAR(derivs.CY_alpha, 0.0, kTestTolerance);
}

TEST_F(StabilityDerivativesTest, RollDampingNegative) {
    // Roll damping Cl_p should be negative (opposes roll rate)
    Aircraft aircraft = CreateRectangularWing(10.0, 1.0, 4, 16);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition condition;
    condition.alpha = 5.0 * kDegreesToRadians;

    solver.Solve(condition);

    StabilityDerivativeCalculator calculator;
    StabilityDerivatives derivs =
        calculator.ComputeDerivatives(aircraft, solver, condition);

    EXPECT_TRUE(derivs.is_valid);
    EXPECT_LT(derivs.Cl_p, 0.0);
}

TEST_F(StabilityDerivativesTest, PitchDampingFinite) {
    // Pitch damping Cm_q for a wing alone (no tail) is typically small
    // For a complete aircraft with tail, Cm_q should be negative (opposes pitch
    // rate)
    Aircraft aircraft = CreateRectangularWing(10.0, 1.0, 4, 16);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition condition;
    condition.alpha = 5.0 * kDegreesToRadians;

    solver.Solve(condition);

    StabilityDerivativeCalculator calculator;
    StabilityDerivatives derivs =
        calculator.ComputeDerivatives(aircraft, solver, condition);

    EXPECT_TRUE(derivs.is_valid);
    // For wing-alone, Cm_q may be positive or small; just check it's finite
    EXPECT_TRUE(std::isfinite(derivs.Cm_q));
}

TEST_F(StabilityDerivativesTest, YawDampingNonPositive) {
    // Yaw damping Cn_r should be non-positive (opposes yaw rate)
    // For a planar wing without vertical surfaces, Cn_r is approximately zero.
    // Actual negative yaw damping comes from vertical tail in real aircraft.
    Aircraft aircraft = CreateRectangularWing(10.0, 1.0, 4, 16);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition condition;
    condition.alpha = 5.0 * kDegreesToRadians;

    solver.Solve(condition);

    StabilityDerivativeCalculator calculator;
    StabilityDerivatives derivs =
        calculator.ComputeDerivatives(aircraft, solver, condition);

    EXPECT_TRUE(derivs.is_valid);
    // For planar wing, Cn_r should be approximately zero or slightly negative
    EXPECT_LE(derivs.Cn_r, kTestTolerance);
}

TEST_F(StabilityDerivativesTest, BodyAxisDerivativesValid) {
    Aircraft aircraft = CreateRectangularWing(10.0, 1.0, 4, 16);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition condition;
    condition.alpha = 5.0 * kDegreesToRadians;

    solver.Solve(condition);

    StabilityDerivativeCalculator calculator;
    BodyAxisDerivatives derivs =
        calculator.ComputeBodyAxisDerivatives(aircraft, solver, condition);

    EXPECT_TRUE(derivs.is_valid);
    // CZ_alpha in Z-up convention: increasing alpha increases lift (positive Z
    // force) Therefore CZ_alpha should be positive
    EXPECT_GT(derivs.CZ_alpha, 0.0);
}

TEST_F(StabilityDerivativesTest, DerivativesAtDifferentAlpha) {
    // Derivatives should be relatively constant over small alpha range
    Aircraft aircraft = CreateRectangularWing(10.0, 1.0, 4, 16);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    StabilityDerivativeCalculator calculator;

    FlightCondition condition1;
    condition1.alpha = 2.0 * kDegreesToRadians;
    solver.Solve(condition1);
    StabilityDerivatives derivs1 =
        calculator.ComputeDerivatives(aircraft, solver, condition1);

    FlightCondition condition2;
    condition2.alpha = 8.0 * kDegreesToRadians;
    solver.Solve(condition2);
    StabilityDerivatives derivs2 =
        calculator.ComputeDerivatives(aircraft, solver, condition2);

    // Lift slope should be similar (within 20%)
    EXPECT_NEAR(derivs1.CL_alpha, derivs2.CL_alpha, 0.2 * derivs1.CL_alpha);

    // Roll damping should be similar
    EXPECT_NEAR(derivs1.Cl_p, derivs2.Cl_p, 0.2 * std::abs(derivs1.Cl_p));
}

TEST_F(StabilityDerivativesTest, StaticMarginComputed) {
    // For a wing with CG at quarter chord, static margin should be computed
    Aircraft aircraft = CreateRectangularWing(10.0, 1.0, 4, 16);
    aircraft.SetMomentReferencePoint(Vec3(0.25, 0.0, 0.0));

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition condition;
    condition.alpha = 5.0 * kDegreesToRadians;

    solver.Solve(condition);

    StabilityDerivativeCalculator calculator;
    StabilityDerivatives derivs =
        calculator.ComputeDerivatives(aircraft, solver, condition);

    EXPECT_TRUE(derivs.is_valid);
    // Static margin should be finite
    EXPECT_TRUE(std::isfinite(derivs.static_margin));
}

TEST_F(StabilityDerivativesTest, ReferenceQuantitiesPreserved) {
    Scalar span = 8.0;
    Scalar chord = 2.0;
    Aircraft aircraft = CreateRectangularWing(span, chord, 4, 8);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition condition;
    condition.alpha = 5.0 * kDegreesToRadians;

    solver.Solve(condition);

    StabilityDerivativeCalculator calculator;
    StabilityDerivatives derivs =
        calculator.ComputeDerivatives(aircraft, solver, condition);

    EXPECT_NEAR(derivs.reference_area, span * chord, kTestTolerance);
    EXPECT_NEAR(derivs.reference_chord, chord, kTestTolerance);
    EXPECT_NEAR(derivs.reference_span, span, kTestTolerance);
}

TEST_F(StabilityDerivativesTest, HighAspectRatioRollDamping) {
    // Higher aspect ratio should give more roll damping (more negative Cl_p)
    Aircraft aircraft_low = CreateRectangularWing(4.0, 1.0, 4, 8);
    Aircraft aircraft_high = CreateRectangularWing(16.0, 1.0, 4, 32);

    VortexLatticeSolver solver_low, solver_high;
    solver_low.SetAircraft(aircraft_low);
    solver_high.SetAircraft(aircraft_high);

    FlightCondition condition;
    condition.alpha = 5.0 * kDegreesToRadians;

    solver_low.Solve(condition);
    solver_high.Solve(condition);

    StabilityDerivativeCalculator calculator;
    StabilityDerivatives derivs_low =
        calculator.ComputeDerivatives(aircraft_low, solver_low, condition);
    StabilityDerivatives derivs_high =
        calculator.ComputeDerivatives(aircraft_high, solver_high, condition);

    // Higher AR wing should have more negative roll damping
    EXPECT_LT(derivs_high.Cl_p, derivs_low.Cl_p);
}

}  // namespace
}  // namespace pegasus
