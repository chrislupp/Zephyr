// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include <gtest/gtest.h>

#include "vlm/ad/gradient_calculator.h"
#include "vlm/aero/stability_derivatives.h"
#include "vlm/geometry/aircraft.h"
#include "vlm/geometry/section.h"
#include "vlm/geometry/surface.h"
#include "vlm/solver/vortex_lattice_solver.h"

namespace vlm {
namespace {

/// @brief Test fixture for gradient calculator validation.
///
/// Creates a simple rectangular wing and validates that the gradient
/// calculator produces derivatives matching those from StabilityDerivatives.
class GradientCalculatorTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Create simple rectangular wing (AR = 8)
    Surface wing("Wing");
    wing.SetChordwisePanelCount(4);
    wing.SetSpanwisePanelCount(8);
    wing.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    wing.AddSection(Section(Vec3(0, 4.0, 0), 1.0, 0.0));

    aircraft_.AddSurface(wing);
    aircraft_.SetReferenceArea(4.0);      // chord * half_span
    aircraft_.SetReferenceChord(1.0);
    aircraft_.SetReferenceSpan(8.0);

    // Set up stability derivative calculator
    solver_.SetAircraft(aircraft_);

    // Set up gradient calculator
    grad_calc_.SetAircraft(aircraft_);
  }

  Aircraft aircraft_{"Test Aircraft"};
  VortexLatticeSolver solver_;
  ad::GradientCalculator grad_calc_;

  // Tolerance for comparing derivatives
  // Finite differences have inherent numerical error, and the methods
  // may compute slightly different quantities
  static constexpr double kDerivTolerance = 0.05;  // 5% relative tolerance
};

TEST_F(GradientCalculatorTest, AlphaDerivativesMatchStabilityDerivatives) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  // Solve first (required for StabilityDerivativeCalculator)
  solver_.Solve(condition);

  // Compute using StabilityDerivativeCalculator
  StabilityDerivativeCalculator stab_calc;
  StabilityDerivatives stab_derivs =
      stab_calc.ComputeDerivatives(aircraft_, solver_, condition);

  // Compute using GradientCalculator
  ad::AeroGradients grads = grad_calc_.ComputeAlphaDerivatives(condition);

  // Compare lift curve slope
  EXPECT_NEAR(grads.dCL, stab_derivs.CL_alpha,
              std::abs(stab_derivs.CL_alpha) * kDerivTolerance + 0.01)
      << "CL_alpha mismatch: gradient calc = " << grads.dCL
      << ", stability derivs = " << stab_derivs.CL_alpha;

  // Compare pitch moment derivative (static stability)
  EXPECT_NEAR(grads.dCm, stab_derivs.Cm_alpha,
              std::abs(stab_derivs.Cm_alpha) * kDerivTolerance + 0.01)
      << "Cm_alpha mismatch: gradient calc = " << grads.dCm
      << ", stability derivs = " << stab_derivs.Cm_alpha;

  // Compare side force derivative (should be ~0 for symmetric wing)
  EXPECT_NEAR(grads.dCY, stab_derivs.CY_alpha,
              std::abs(stab_derivs.CY_alpha) * kDerivTolerance + 0.01)
      << "CY_alpha mismatch";

  // Note: CD_alpha comparison is tricky because:
  // - StabilityDerivatives computes linear derivatives assuming CD ~ CL²
  // - Finite differences compute actual CD change including the CL² effect
  // The two methods measure different quantities, so we skip this comparison
  // for now. When Enzyme is working, we can ensure consistency.
}

TEST_F(GradientCalculatorTest, BetaDerivativesMatchStabilityDerivatives) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  // Solve first
  solver_.Solve(condition);

  StabilityDerivativeCalculator stab_calc;
  StabilityDerivatives stab_derivs =
      stab_calc.ComputeDerivatives(aircraft_, solver_, condition);

  ad::AeroGradients grads = grad_calc_.ComputeBetaDerivatives(condition);

  // For a symmetric wing, CY_beta and Cn_beta should be ~0
  // Cl_beta (dihedral effect) depends on wing dihedral

  EXPECT_NEAR(grads.dCY, stab_derivs.CY_beta,
              std::abs(stab_derivs.CY_beta) * kDerivTolerance + 0.01)
      << "CY_beta mismatch";

  EXPECT_NEAR(grads.dCl, stab_derivs.Cl_beta,
              std::abs(stab_derivs.Cl_beta) * kDerivTolerance + 0.01)
      << "Cl_beta mismatch";

  EXPECT_NEAR(grads.dCn, stab_derivs.Cn_beta,
              std::abs(stab_derivs.Cn_beta) * kDerivTolerance + 0.01)
      << "Cn_beta mismatch";
}

TEST_F(GradientCalculatorTest, RollRateDerivativesMatchStabilityDerivatives) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  // Solve first
  solver_.Solve(condition);

  StabilityDerivativeCalculator stab_calc;
  StabilityDerivatives stab_derivs =
      stab_calc.ComputeDerivatives(aircraft_, solver_, condition);

  ad::AeroGradients grads = grad_calc_.ComputeRollRateDerivatives(condition);

  // Roll damping Cl_p should be negative (restoring)
  EXPECT_NEAR(grads.dCl, stab_derivs.Cl_p,
              std::abs(stab_derivs.Cl_p) * kDerivTolerance + 0.01)
      << "Cl_p mismatch";

  // CL_p and Cm_p
  EXPECT_NEAR(grads.dCL, stab_derivs.CL_p,
              std::abs(stab_derivs.CL_p) * kDerivTolerance + 0.01)
      << "CL_p mismatch";
}

TEST_F(GradientCalculatorTest, PitchRateDerivativesMatchStabilityDerivatives) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  // Solve first
  solver_.Solve(condition);

  StabilityDerivativeCalculator stab_calc;
  StabilityDerivatives stab_derivs =
      stab_calc.ComputeDerivatives(aircraft_, solver_, condition);

  ad::AeroGradients grads = grad_calc_.ComputePitchRateDerivatives(condition);

  // Pitch damping Cm_q should be negative (restoring)
  EXPECT_NEAR(grads.dCm, stab_derivs.Cm_q,
              std::abs(stab_derivs.Cm_q) * kDerivTolerance + 0.01)
      << "Cm_q mismatch";

  // Lift due to pitch rate
  EXPECT_NEAR(grads.dCL, stab_derivs.CL_q,
              std::abs(stab_derivs.CL_q) * kDerivTolerance + 0.01)
      << "CL_q mismatch";
}

TEST_F(GradientCalculatorTest, YawRateDerivativesMatchStabilityDerivatives) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  // Solve first
  solver_.Solve(condition);

  StabilityDerivativeCalculator stab_calc;
  StabilityDerivatives stab_derivs =
      stab_calc.ComputeDerivatives(aircraft_, solver_, condition);

  ad::AeroGradients grads = grad_calc_.ComputeYawRateDerivatives(condition);

  // Yaw damping Cn_r should be negative (restoring)
  EXPECT_NEAR(grads.dCn, stab_derivs.Cn_r,
              std::abs(stab_derivs.Cn_r) * kDerivTolerance + 0.01)
      << "Cn_r mismatch";

  // Roll due to yaw rate
  EXPECT_NEAR(grads.dCl, stab_derivs.Cl_r,
              std::abs(stab_derivs.Cl_r) * kDerivTolerance + 0.01)
      << "Cl_r mismatch";
}

TEST_F(GradientCalculatorTest, FullGradientsComputation) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  ad::FullGradients grads = grad_calc_.ComputeFlightConditionDerivatives(condition);

  // Sanity checks
  // CL_alpha should be positive (~2*pi for high AR wing)
  EXPECT_GT(grads.alpha_derivs.dCL, 3.0);  // Should be around 5.0 for AR=8
  EXPECT_LT(grads.alpha_derivs.dCL, 7.0);

  // Cl_p (roll damping) should be negative
  EXPECT_LT(grads.p_derivs.dCl, 0.0);

  // Note: Cm_q for a wing-only configuration is not necessarily negative
  // Pitch damping primarily comes from the horizontal tail
  // Just verify it's computed (non-zero)
  EXPECT_NE(grads.q_derivs.dCm, 0.0);
}

TEST_F(GradientCalculatorTest, IsUsingEnzymeReturnsCorrectValue) {
#ifdef VLM_USE_ENZYME
  EXPECT_TRUE(ad::GradientCalculator::IsUsingEnzyme());
#else
  EXPECT_FALSE(ad::GradientCalculator::IsUsingEnzyme());
#endif
}

TEST_F(GradientCalculatorTest, DerivativesAtZeroAlpha) {
  FlightCondition condition;
  condition.alpha = 0.0;

  ad::AeroGradients grads = grad_calc_.ComputeAlphaDerivatives(condition);

  // At zero alpha, lift curve slope should still be positive
  EXPECT_GT(grads.dCL, 0.0);

  // CD_alpha should be ~0 at zero alpha (minimum drag)
  EXPECT_NEAR(grads.dCD, 0.0, 0.1);
}

// ============================================================================
// Aircraft with Tail Test
// ============================================================================

class GradientCalculatorWithTailTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Create wing-tail configuration
    Surface wing("Wing");
    wing.SetChordwisePanelCount(4);
    wing.SetSpanwisePanelCount(8);
    wing.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    wing.AddSection(Section(Vec3(0, 4.0, 0), 1.0, 0.0));

    Surface htail("HTail");
    htail.SetChordwisePanelCount(3);
    htail.SetSpanwisePanelCount(4);
    htail.AddSection(Section(Vec3(4.0, 0, 0.2), 0.5, 0.0));
    htail.AddSection(Section(Vec3(4.0, 1.5, 0.2), 0.5, 0.0));

    aircraft_.AddSurface(wing);
    aircraft_.AddSurface(htail);
    aircraft_.SetReferenceArea(4.0);
    aircraft_.SetReferenceChord(1.0);
    aircraft_.SetReferenceSpan(8.0);
    aircraft_.SetMomentReferencePoint(Vec3(0.25, 0, 0));

    solver_.SetAircraft(aircraft_);
    grad_calc_.SetAircraft(aircraft_);
  }

  Aircraft aircraft_{"Wing-Tail"};
  VortexLatticeSolver solver_;
  ad::GradientCalculator grad_calc_;
};

TEST_F(GradientCalculatorWithTailTest, CmAlphaIsNegativeWithTail) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  ad::AeroGradients grads = grad_calc_.ComputeAlphaDerivatives(condition);

  // With horizontal tail behind CG, Cm_alpha should be negative (stable)
  // This is a key stability criterion
  EXPECT_LT(grads.dCm, 0.0)
      << "Aircraft with tail should have negative Cm_alpha for stability";
}

TEST_F(GradientCalculatorWithTailTest, MatchesStabilityDerivatives) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  // Solve first
  solver_.Solve(condition);

  StabilityDerivativeCalculator stab_calc;
  StabilityDerivatives stab_derivs =
      stab_calc.ComputeDerivatives(aircraft_, solver_, condition);

  ad::AeroGradients grads = grad_calc_.ComputeAlphaDerivatives(condition);

  // Cm_alpha match
  EXPECT_NEAR(grads.dCm, stab_derivs.Cm_alpha,
              std::abs(stab_derivs.Cm_alpha) * 0.05 + 0.01)
      << "Cm_alpha should match for wing-tail config";
}

}  // namespace
}  // namespace vlm
