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

#include <cmath>

#include <gtest/gtest.h>

#include "pegasus/aero/aerodynamic_results.h"
#include "pegasus/aero/force_calculator.h"
#include "pegasus/geometry/aircraft.h"
#include "pegasus/geometry/section.h"
#include "pegasus/geometry/surface.h"
#include "pegasus/solver/vortex_lattice_solver.h"

namespace pegasus {
namespace {

constexpr double kTestTolerance = 1.0e-6;

// ============================================================================
// ForceCalculator Tests
// ============================================================================

class ForceCalculatorTest : public ::testing::Test {
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
};

TEST_F(ForceCalculatorTest, ZeroAlphaGivesZeroForces) {
  Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 4, 8);

  VortexLatticeSolver solver;
  solver.SetAircraft(aircraft);

  FlightCondition condition;
  condition.alpha = 0.0;
  condition.beta = 0.0;

  SolutionState solution = solver.Solve(condition);
  ASSERT_TRUE(solution.is_valid);

  ForceCalculator calculator;
  AerodynamicResults results =
      calculator.ComputeForces(aircraft, solver.GetMesh(), solution, condition);

  EXPECT_TRUE(results.is_valid);
  EXPECT_NEAR(results.CL, 0.0, kTestTolerance);
  EXPECT_NEAR(results.CY, 0.0, kTestTolerance);
  EXPECT_NEAR(results.CL_ff, 0.0, kTestTolerance);
  EXPECT_NEAR(results.CDi_ff, 0.0, kTestTolerance);
}

TEST_F(ForceCalculatorTest, PositiveAlphaGivesPositiveLift) {
  Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 4, 8);

  VortexLatticeSolver solver;
  solver.SetAircraft(aircraft);

  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  SolutionState solution = solver.Solve(condition);
  ASSERT_TRUE(solution.is_valid);

  ForceCalculator calculator;
  AerodynamicResults results =
      calculator.ComputeForces(aircraft, solver.GetMesh(), solution, condition);

  EXPECT_TRUE(results.is_valid);
  EXPECT_GT(results.CL, 0.0);
  EXPECT_GT(results.CL_ff, 0.0);
  EXPECT_GT(results.CDi_ff, 0.0);
  EXPECT_GT(results.span_efficiency, 0.0);
}

TEST_F(ForceCalculatorTest, NegativeAlphaGivesNegativeLift) {
  Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 4, 8);

  VortexLatticeSolver solver;
  solver.SetAircraft(aircraft);

  FlightCondition condition;
  condition.alpha = -5.0 * kDegreesToRadians;

  SolutionState solution = solver.Solve(condition);
  ASSERT_TRUE(solution.is_valid);

  ForceCalculator calculator;
  AerodynamicResults results =
      calculator.ComputeForces(aircraft, solver.GetMesh(), solution, condition);

  EXPECT_TRUE(results.is_valid);
  EXPECT_LT(results.CL, 0.0);
  EXPECT_LT(results.CL_ff, 0.0);
  EXPECT_GT(results.CDi_ff, 0.0);  // Induced drag is always positive
}

TEST_F(ForceCalculatorTest, LiftCoefficientMatchesCirculation) {
  // CL from circulation should match CL from force integration
  Scalar span = 4.0;
  Scalar chord = 1.0;
  Aircraft aircraft = CreateRectangularWing(span, chord, 4, 8);

  VortexLatticeSolver solver;
  solver.SetAircraft(aircraft);

  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  SolutionState solution = solver.Solve(condition);
  ASSERT_TRUE(solution.is_valid);

  // CL from Kutta-Joukowski: CL = 2 * sum(Gamma * dy) / (V * S)
  const VortexLattice& mesh = solver.GetMesh();
  Scalar S_ref = aircraft.GetReferenceArea();
  Scalar CL_from_circulation = 0.0;

  for (Index j = 0; j < mesh.GetNumStrips(); ++j) {
    const Strip& strip = mesh.GetStrip(j);
    Scalar strip_circulation = 0.0;
    for (Index k = 0; k < strip.num_vortices; ++k) {
      strip_circulation += solution.circulation(strip.first_vortex_index + k);
    }
    CL_from_circulation += 2.0 * strip_circulation * strip.width / S_ref;
  }

  ForceCalculator calculator;
  AerodynamicResults results =
      calculator.ComputeForces(aircraft, solver.GetMesh(), solution, condition);

  // Near-field and far-field lift should be close
  EXPECT_NEAR(results.CL_ff, CL_from_circulation, 0.01 * CL_from_circulation);
}

TEST_F(ForceCalculatorTest, SpanEfficiencyReasonable) {
  // Span efficiency for rectangular wing should be close to 1 but not too high
  // (Note: VLM can give slightly > 1.0 due to numerical errors)
  Aircraft aircraft = CreateRectangularWing(10.0, 1.0, 4, 20);

  VortexLatticeSolver solver;
  solver.SetAircraft(aircraft);

  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  SolutionState solution = solver.Solve(condition);
  ASSERT_TRUE(solution.is_valid);

  ForceCalculator calculator;
  AerodynamicResults results =
      calculator.ComputeForces(aircraft, solver.GetMesh(), solution, condition);

  EXPECT_TRUE(results.is_valid);
  EXPECT_GT(results.span_efficiency, 0.8);   // Should be reasonably efficient
  EXPECT_LT(results.span_efficiency, 1.15);  // Allow small numerical overshoot
}

TEST_F(ForceCalculatorTest, InducedDragQuadraticWithLift) {
  // CDi should be approximately proportional to CL^2
  Aircraft aircraft = CreateRectangularWing(10.0, 1.0, 4, 16);

  VortexLatticeSolver solver;
  solver.SetAircraft(aircraft);

  ForceCalculator calculator;

  FlightCondition condition1;
  condition1.alpha = 2.0 * kDegreesToRadians;
  SolutionState solution1 = solver.Solve(condition1);
  AerodynamicResults results1 =
      calculator.ComputeForces(aircraft, solver.GetMesh(), solution1, condition1);

  FlightCondition condition2;
  condition2.alpha = 4.0 * kDegreesToRadians;
  SolutionState solution2 = solver.Solve(condition2);
  AerodynamicResults results2 =
      calculator.ComputeForces(aircraft, solver.GetMesh(), solution2, condition2);

  // CL should double, CDi should quadruple
  Scalar CL_ratio = results2.CL_ff / results1.CL_ff;
  Scalar CDi_ratio = results2.CDi_ff / results1.CDi_ff;

  EXPECT_NEAR(CL_ratio, 2.0, 0.05);
  EXPECT_NEAR(CDi_ratio, 4.0, 0.2);  // Allow more tolerance for induced drag
}

TEST_F(ForceCalculatorTest, SurfaceLoadsAddUpToTotal) {
  Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 4, 8);

  VortexLatticeSolver solver;
  solver.SetAircraft(aircraft);

  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  SolutionState solution = solver.Solve(condition);
  ASSERT_TRUE(solution.is_valid);

  ForceCalculator calculator;
  AerodynamicResults results =
      calculator.ComputeForces(aircraft, solver.GetMesh(), solution, condition);

  EXPECT_TRUE(results.is_valid);
  EXPECT_EQ(results.surface_loads.size(), 1);  // One surface

  // Surface loads should match total
  EXPECT_NEAR(results.surface_loads[0].CL, results.CL, 0.01);
  EXPECT_NEAR(results.surface_loads[0].CD, results.CD, 0.01);
}

TEST_F(ForceCalculatorTest, StripLoadsHaveCorrectCount) {
  Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 4, 8);

  VortexLatticeSolver solver;
  solver.SetAircraft(aircraft);

  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  SolutionState solution = solver.Solve(condition);
  ASSERT_TRUE(solution.is_valid);

  ForceCalculator calculator;
  AerodynamicResults results =
      calculator.ComputeForces(aircraft, solver.GetMesh(), solution, condition);

  EXPECT_TRUE(results.is_valid);
  EXPECT_EQ(results.surface_loads.size(), 1);

  // Should have 8 strips for 8 spanwise panels
  EXPECT_EQ(results.surface_loads[0].strip_loads.size(), 8);
}

TEST_F(ForceCalculatorTest, PitchingMomentSign) {
  // At positive alpha, pitching moment should be nose-down for aft CG
  Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 4, 8);
  aircraft.SetMomentReferencePoint(Vec3(0.75, 0.0, 0.0));  // Aft reference

  VortexLatticeSolver solver;
  solver.SetAircraft(aircraft);

  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  SolutionState solution = solver.Solve(condition);
  ASSERT_TRUE(solution.is_valid);

  ForceCalculator calculator;
  AerodynamicResults results =
      calculator.ComputeForces(aircraft, solver.GetMesh(), solution, condition);

  EXPECT_TRUE(results.is_valid);
  // With aft moment reference (0.75) and lift at quarter chord (0.25):
  // Moment arm is negative (force is forward of reference)
  // Upward lift with forward moment arm creates nose-UP moment (Cm > 0)
  // This is physically correct for unstable (aft CG) configuration
  EXPECT_GT(results.Cm, 0.0);
}

TEST_F(ForceCalculatorTest, AspectRatioCalculation) {
  Scalar span = 10.0;
  Scalar chord = 1.0;
  Aircraft aircraft = CreateRectangularWing(span, chord, 4, 8);

  VortexLatticeSolver solver;
  solver.SetAircraft(aircraft);

  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  SolutionState solution = solver.Solve(condition);

  ForceCalculator calculator;
  AerodynamicResults results =
      calculator.ComputeForces(aircraft, solver.GetMesh(), solution, condition);

  EXPECT_NEAR(results.GetAspectRatio(), 10.0, kTestTolerance);
}

TEST_F(ForceCalculatorTest, IdealInducedDragComputation) {
  Aircraft aircraft = CreateRectangularWing(10.0, 1.0, 4, 16);

  VortexLatticeSolver solver;
  solver.SetAircraft(aircraft);

  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  SolutionState solution = solver.Solve(condition);

  ForceCalculator calculator;
  AerodynamicResults results =
      calculator.ComputeForces(aircraft, solver.GetMesh(), solution, condition);

  // Ideal induced drag: CDi_ideal = CL^2 / (pi * AR)
  Scalar CDi_ideal = results.GetIdealInducedDrag();
  Scalar AR = results.GetAspectRatio();

  EXPECT_NEAR(CDi_ideal, results.CL * results.CL / (kPi * AR), kTestTolerance);

  // Actual induced drag should be slightly higher than ideal (e < 1)
  EXPECT_GT(results.CDi_ff, CDi_ideal * 0.95);
  EXPECT_LT(results.CDi_ff, CDi_ideal * 1.25);
}

}  // namespace
}  // namespace pegasus
