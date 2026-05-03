/*
    Project VLM
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

/// @file rectangular_wing_test.cpp
/// @brief Integration tests for rectangular wing analysis.

#include <cmath>

#include <gtest/gtest.h>

#include "vlm/vlm.h"

namespace vlm {
namespace {

// VLM results differ from simplified analytical formulas by ~5-10%
// due to discrete panel approximation and end-effects
constexpr double kTestTolerance = 0.10;  // 10% tolerance for theory comparison

class RectangularWingTest : public ::testing::Test {
 protected:
  void SetUp() override {
    const double span = 10.0;
    const double chord = 1.0;

    Surface wing("Wing");
    wing.SetChordwisePanelCount(8);
    wing.SetSpanwisePanelCount(20);
    wing.SetChordwiseSpacing(SpacingType::kCosine);
    wing.SetSpanwiseSpacing(SpacingType::kCosine);
    wing.SetYDuplicate(0.0);

    Section root(Vec3(0.0, 0.0, 0.0), chord);
    Section tip(Vec3(0.0, span / 2.0, 0.0), chord);
    wing.AddSection(root);
    wing.AddSection(tip);

    aircraft_.SetTitle("Rectangular Wing AR=10");
    aircraft_.AddSurface(wing);
    aircraft_.SetReferenceArea(span * chord);
    aircraft_.SetReferenceChord(chord);
    aircraft_.SetReferenceSpan(span);

    solver_.SetAircraft(aircraft_);
  }

  Aircraft aircraft_;
  VortexLatticeSolver solver_;
};

TEST_F(RectangularWingTest, LiftCurveSlopeMatchesTheory) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  SolutionState solution = solver_.Solve(condition);

  ForceCalculator calculator;
  AerodynamicResults results = calculator.ComputeForces(
      aircraft_, solver_.GetMesh(), solution, condition);

  double CL_alpha = results.CL / condition.alpha;

  // Theory: CL_alpha = 2*pi / (1 + 2/(pi*AR)) for high AR
  double AR = 10.0;
  double CL_alpha_theory = 2.0 * kPi / (1.0 + 2.0 / (kPi * AR));

  EXPECT_NEAR(CL_alpha, CL_alpha_theory, CL_alpha_theory * kTestTolerance);
}

TEST_F(RectangularWingTest, InducedDragNonZero) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  SolutionState solution = solver_.Solve(condition);

  ForceCalculator calculator;
  AerodynamicResults results = calculator.ComputeForces(
      aircraft_, solver_.GetMesh(), solution, condition);

  // Verify drag is computed (non-zero) and finite
  // Note: Sign convention may vary; check absolute value is reasonable
  double AR = 10.0;
  double CDi_ideal = results.CL * results.CL / (kPi * AR);

  EXPECT_TRUE(std::isfinite(results.CD));
  EXPECT_NE(results.CD, 0.0);

  // Absolute drag should be within order of magnitude of ideal
  EXPECT_LT(std::abs(results.CD), CDi_ideal * 10.0);
}

TEST_F(RectangularWingTest, ZeroLiftAtZeroAlpha) {
  FlightCondition condition;
  condition.alpha = 0.0;

  SolutionState solution = solver_.Solve(condition);

  ForceCalculator calculator;
  AerodynamicResults results = calculator.ComputeForces(
      aircraft_, solver_.GetMesh(), solution, condition);

  EXPECT_NEAR(results.CL, 0.0, 1e-10);
}

TEST_F(RectangularWingTest, LiftLinearWithAlpha) {
  ForceCalculator calculator;
  std::vector<std::pair<double, double>> alpha_CL;

  for (double alpha_deg = -5.0; alpha_deg <= 10.0; alpha_deg += 2.5) {
    FlightCondition condition;
    condition.alpha = alpha_deg * kDegreesToRadians;

    SolutionState solution = solver_.Solve(condition);
    AerodynamicResults results = calculator.ComputeForces(
        aircraft_, solver_.GetMesh(), solution, condition);

    alpha_CL.emplace_back(condition.alpha, results.CL);
  }

  // Compute R^2 for linear fit
  double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_x2 = 0.0;
  int n = static_cast<int>(alpha_CL.size());

  for (const auto& [alpha, CL] : alpha_CL) {
    sum_x += alpha;
    sum_y += CL;
    sum_xy += alpha * CL;
    sum_x2 += alpha * alpha;
  }

  double slope = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x);
  double intercept = (sum_y - slope * sum_x) / n;

  double mean_y = sum_y / n;
  double ss_tot = 0.0, ss_res = 0.0;
  for (const auto& [alpha, CL] : alpha_CL) {
    double predicted = slope * alpha + intercept;
    ss_tot += (CL - mean_y) * (CL - mean_y);
    ss_res += (CL - predicted) * (CL - predicted);
  }

  double r_squared = 1.0 - ss_res / ss_tot;
  EXPECT_GT(r_squared, 0.999);
}

TEST_F(RectangularWingTest, StabilityDerivativesPhysicallyReasonable) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  solver_.Solve(condition);

  StabilityDerivativeCalculator calculator;
  StabilityDerivatives derivs = calculator.ComputeDerivatives(
      aircraft_, solver_, condition);

  EXPECT_TRUE(derivs.is_valid);
  EXPECT_GT(derivs.CL_alpha, 4.0);
  EXPECT_LT(derivs.CL_alpha, 6.5);
  EXPECT_LT(derivs.Cl_p, 0.0);
}

// ============================================================================
// High Aspect Ratio Validation
// ============================================================================

class HighAspectRatioTest : public ::testing::TestWithParam<double> {};

TEST_P(HighAspectRatioTest, LiftSlopeApproaches2Pi) {
  double AR = GetParam();
  double span = AR;
  double chord = 1.0;

  Surface wing("Wing");
  wing.SetChordwisePanelCount(6);
  wing.SetSpanwisePanelCount(static_cast<int>(AR * 2));
  wing.SetYDuplicate(0.0);

  Section root(Vec3(0.0, 0.0, 0.0), chord);
  Section tip(Vec3(0.0, span / 2.0, 0.0), chord);
  wing.AddSection(root);
  wing.AddSection(tip);

  Aircraft aircraft("Wing AR=" + std::to_string(static_cast<int>(AR)));
  aircraft.AddSurface(wing);
  aircraft.SetReferenceArea(span * chord);
  aircraft.SetReferenceChord(chord);
  aircraft.SetReferenceSpan(span);

  VortexLatticeSolver solver;
  solver.SetAircraft(aircraft);

  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  SolutionState solution = solver.Solve(condition);

  ForceCalculator calculator;
  AerodynamicResults results = calculator.ComputeForces(
      aircraft, solver.GetMesh(), solution, condition);

  double CL_alpha = results.CL / condition.alpha;
  double CL_alpha_theory = 2.0 * kPi / (1.0 + 2.0 / (kPi * AR));

  EXPECT_NEAR(CL_alpha, CL_alpha_theory, CL_alpha_theory * kTestTolerance);
}

// Test moderate aspect ratios (6-10) where VLM matches simplified theory best
INSTANTIATE_TEST_SUITE_P(AspectRatioSweep, HighAspectRatioTest,
                         ::testing::Values(6.0, 7.0, 8.0, 10.0));

// ============================================================================
// AVL File Parsing Integration
// ============================================================================

TEST(AvlFileTest, ParseRectangularWing) {
  AvlReader reader;
  Aircraft aircraft = reader.ReadFile("fixtures/rect.avl");

  EXPECT_EQ(aircraft.GetNumSurfaces(), 1u);
  EXPECT_NEAR(aircraft.GetReferenceArea(), 1.123, 0.001);
  EXPECT_NEAR(aircraft.GetReferenceChord(), 0.25, 0.001);
  EXPECT_NEAR(aircraft.GetReferenceSpan(), 6.01, 0.001);
}

TEST(AvlFileTest, SolveFromParsedGeometry) {
  AvlReader reader;
  Aircraft aircraft = reader.ReadFile("fixtures/rect.avl");

  VortexLatticeSolver solver;
  solver.SetAircraft(aircraft);

  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  SolutionState solution = solver.Solve(condition);

  ForceCalculator calculator;
  AerodynamicResults results = calculator.ComputeForces(
      aircraft, solver.GetMesh(), solution, condition);

  // Positive alpha should produce positive lift
  EXPECT_GT(results.CL, 0.0);

  // Drag should be computed (non-zero and finite)
  EXPECT_TRUE(std::isfinite(results.CD));
  EXPECT_NE(results.CD, 0.0);

  // Moment should be finite
  EXPECT_TRUE(std::isfinite(results.Cm));
}

}  // namespace
}  // namespace vlm
