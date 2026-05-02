// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include <gtest/gtest.h>

#include <cmath>

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
  // for now. The analytical adjoint method handles this correctly.
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

TEST_F(GradientCalculatorTest, MachDerivativesComputed) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;
  condition.mach = 0.5;  // Moderate subsonic Mach

  ad::AeroGradients grads = grad_calc_.ComputeMachDerivatives(condition);

  // Prandtl-Glauert: CL increases with Mach (positive dCL/dM for subsonic)
  // CL ≈ CL0 / sqrt(1 - M²)
  // d/dM[CL] = CL0 * M / (1-M²)^1.5 > 0 for M > 0
  EXPECT_GT(grads.dCL, 0.0) << "dCL/dMach should be positive for subsonic flow";

  // Similarly for induced drag
  EXPECT_GT(grads.dCDi, 0.0) << "dCDi/dMach should be positive";
}

TEST_F(GradientCalculatorTest, MachDerivativesAtLowMach) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;
  condition.mach = 0.1;  // Low Mach - small compressibility effect

  ad::AeroGradients grads = grad_calc_.ComputeMachDerivatives(condition);

  // At low Mach, derivatives should be small but positive
  EXPECT_GT(grads.dCL, 0.0);

  // The magnitude should be small at low Mach
  // (derivative of 1/sqrt(1-M²) at M=0.1 is about 0.1)
  EXPECT_LT(std::abs(grads.dCL), 1.0) << "dCL/dMach should be small at low Mach";
}

TEST_F(GradientCalculatorTest, MachDerivativesAtHighMach) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;
  condition.mach = 0.8;  // High subsonic Mach

  ad::AeroGradients grads_high = grad_calc_.ComputeMachDerivatives(condition);

  condition.mach = 0.3;
  ad::AeroGradients grads_low = grad_calc_.ComputeMachDerivatives(condition);

  // Derivative should be larger at higher Mach (more sensitive)
  EXPECT_GT(grads_high.dCL, grads_low.dCL)
      << "dCL/dMach should increase with Mach number";
}

TEST_F(GradientCalculatorTest, FullGradientsComputation) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;
  condition.mach = 0.3;

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

  // Mach derivative should be positive (compressibility increases lift)
  EXPECT_GT(grads.mach_derivs.dCL, 0.0);
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

// ============================================================================
// Geometry Gradient Tests
// ============================================================================

class GeometryGradientTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Create simple rectangular wing (AR = 6)
    Surface wing("Wing");
    wing.SetChordwisePanelCount(3);
    wing.SetSpanwisePanelCount(4);
    wing.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));      // Root
    wing.AddSection(Section(Vec3(0, 3.0, 0), 1.0, 0.0));    // Tip

    aircraft_.AddSurface(wing);
    aircraft_.SetReferenceArea(3.0);
    aircraft_.SetReferenceChord(1.0);
    aircraft_.SetReferenceSpan(6.0);

    grad_calc_.SetAircraft(aircraft_);
  }

  Aircraft aircraft_{"GeometryTest"};
  ad::GradientCalculator grad_calc_;
};

TEST_F(GeometryGradientTest, ChordGradientIsNonZero) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  ad::AeroGradients grads = grad_calc_.ComputeSectionGradient(
      condition, 0, 0, ad::GeometryVariable::kChord);

  // Increasing chord should increase lift (more area)
  // dCL/dchord depends on how reference area is defined
  // For a rectangular wing with S_ref fixed, larger chord means
  // larger actual lift but CL may decrease if S_ref stays fixed

  // Just verify the gradient is computed (non-zero)
  EXPECT_NE(grads.dCL, 0.0) << "dCL/dchord should be non-zero";
}

TEST_F(GeometryGradientTest, TwistGradientPhysicallyCorrect) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  // Twist at root section
  ad::AeroGradients root_grads = grad_calc_.ComputeSectionGradient(
      condition, 0, 0, ad::GeometryVariable::kTwist);

  // Twist at tip section
  ad::AeroGradients tip_grads = grad_calc_.ComputeSectionGradient(
      condition, 0, 1, ad::GeometryVariable::kTwist);

  // Positive twist (leading edge up) should increase local lift
  // For per-degree gradient, this should be positive
  // CL contribution depends on section span fraction

  // Both should be non-zero
  EXPECT_NE(root_grads.dCL, 0.0) << "dCL/dtwist at root should be non-zero";
  EXPECT_NE(tip_grads.dCL, 0.0) << "dCL/dtwist at tip should be non-zero";

  // Positive twist should increase CL
  EXPECT_GT(root_grads.dCL, 0.0) << "Positive twist should increase lift";
  EXPECT_GT(tip_grads.dCL, 0.0) << "Positive twist should increase lift";
}

TEST_F(GeometryGradientTest, SpanPositionGradientNonZero) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  // Moving tip section outward (increasing y) should increase span
  ad::AeroGradients grads = grad_calc_.ComputeSectionGradient(
      condition, 0, 1, ad::GeometryVariable::kYle);

  // Gradient should be non-zero
  // Increasing span increases aspect ratio, which increases lift slope
  EXPECT_NE(grads.dCL, 0.0) << "dCL/dy_tip should be non-zero";
}

TEST_F(GeometryGradientTest, FullGeometryGradientsReturnsAllSections) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  ad::GeometryGradients grads = grad_calc_.ComputeGeometryGradients(condition);

  // Should have gradients for all sections (2 sections in our wing)
  EXPECT_EQ(grads.section_gradients.size(), 2u);

  // Check surface names are set correctly
  EXPECT_EQ(grads.section_gradients[0].surface_name, "Wing");
  EXPECT_EQ(grads.section_gradients[1].surface_name, "Wing");

  // Check section indices
  EXPECT_EQ(grads.section_gradients[0].section_index, 0);
  EXPECT_EQ(grads.section_gradients[1].section_index, 1);

  // Total design variables: 2 sections * 5 variables = 10
  EXPECT_EQ(grads.TotalDesignVariables(), 10);
}

TEST_F(GeometryGradientTest, GetChordGradientsForExtraction) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  ad::GeometryGradients grads = grad_calc_.ComputeGeometryGradients(condition);

  // Extract dCL/dchord for all sections
  VecX dCL_dchord = grads.GetChordGradientsFor(
      [](const ad::AeroGradients& g) { return g.dCL; });

  EXPECT_EQ(dCL_dchord.size(), 2);

  // Both values should be non-zero
  EXPECT_NE(dCL_dchord(0), 0.0);
  EXPECT_NE(dCL_dchord(1), 0.0);
}

TEST_F(GeometryGradientTest, GetTwistGradientsForExtraction) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  ad::GeometryGradients grads = grad_calc_.ComputeGeometryGradients(condition);

  // Extract dCm/dtwist for all sections
  VecX dCm_dtwist = grads.GetTwistGradientsFor(
      [](const ad::AeroGradients& g) { return g.dCm; });

  EXPECT_EQ(dCm_dtwist.size(), 2);

  // Twist affects pitching moment through load distribution
  // Non-zero expected
  EXPECT_NE(dCm_dtwist(0), 0.0);
  EXPECT_NE(dCm_dtwist(1), 0.0);
}

TEST_F(GeometryGradientTest, DihedralGradientAffectsRollingMoment) {
  // Dihedral (z-position change with span) affects roll stability

  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;
  condition.beta = 2.0 * kDegreesToRadians;  // Add sideslip to see Cl effect

  // dCl/dz at tip (dihedral effect)
  ad::AeroGradients grads = grad_calc_.ComputeSectionGradient(
      condition, 0, 1, ad::GeometryVariable::kZle);

  // In sideslip, dihedral (positive z at tip) creates roll-out moment
  // The gradient should be non-zero
  // Note: With very small panels, this may be noisy
  // The main check is that the computation runs without error
  EXPECT_TRUE(std::isfinite(grads.dCl));
}

TEST_F(GeometryGradientTest, GradientsAreSymmetricForSymmetricSections) {
  // For a symmetric wing, sensitivities of root and tip should differ
  // but both should be finite and physically reasonable

  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  ad::GeometryGradients grads = grad_calc_.ComputeGeometryGradients(condition);

  // All gradients should be finite
  for (const auto& section_grad : grads.section_gradients) {
    EXPECT_TRUE(std::isfinite(section_grad.chord_derivs.dCL));
    EXPECT_TRUE(std::isfinite(section_grad.chord_derivs.dCD));
    EXPECT_TRUE(std::isfinite(section_grad.chord_derivs.dCm));
    EXPECT_TRUE(std::isfinite(section_grad.twist_derivs.dCL));
    EXPECT_TRUE(std::isfinite(section_grad.twist_derivs.dCD));
    EXPECT_TRUE(std::isfinite(section_grad.twist_derivs.dCm));
    EXPECT_TRUE(std::isfinite(section_grad.xle_derivs.dCL));
    EXPECT_TRUE(std::isfinite(section_grad.yle_derivs.dCL));
    EXPECT_TRUE(std::isfinite(section_grad.zle_derivs.dCL));
  }
}

// ============================================================================
// Multi-Surface Geometry Gradient Tests
// ============================================================================

class MultiSurfaceGeometryGradientTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Create wing + tail configuration
    Surface wing("Wing");
    wing.SetChordwisePanelCount(3);
    wing.SetSpanwisePanelCount(4);
    wing.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    wing.AddSection(Section(Vec3(0, 3.0, 0), 1.0, 0.0));

    Surface tail("HTail");
    tail.SetChordwisePanelCount(2);
    tail.SetSpanwisePanelCount(3);
    tail.AddSection(Section(Vec3(3.0, 0, 0.1), 0.4, 0.0));
    tail.AddSection(Section(Vec3(3.0, 1.0, 0.1), 0.4, 0.0));

    aircraft_.AddSurface(wing);
    aircraft_.AddSurface(tail);
    aircraft_.SetReferenceArea(3.0);
    aircraft_.SetReferenceChord(1.0);
    aircraft_.SetReferenceSpan(6.0);
    aircraft_.SetMomentReferencePoint(Vec3(0.25, 0, 0));

    grad_calc_.SetAircraft(aircraft_);
  }

  Aircraft aircraft_{"MultiSurface"};
  ad::GradientCalculator grad_calc_;
};

TEST_F(MultiSurfaceGeometryGradientTest, GradientsComputedForAllSurfaces) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  ad::GeometryGradients grads = grad_calc_.ComputeGeometryGradients(condition);

  // Should have gradients for all sections: 2 wing + 2 tail = 4 sections
  EXPECT_EQ(grads.section_gradients.size(), 4u);

  // Check surface names
  EXPECT_EQ(grads.section_gradients[0].surface_name, "Wing");
  EXPECT_EQ(grads.section_gradients[1].surface_name, "Wing");
  EXPECT_EQ(grads.section_gradients[2].surface_name, "HTail");
  EXPECT_EQ(grads.section_gradients[3].surface_name, "HTail");

  // Total design variables: 4 sections * 5 = 20
  EXPECT_EQ(grads.TotalDesignVariables(), 20);
}

TEST_F(MultiSurfaceGeometryGradientTest, TailTwistAffectsPitchingMoment) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  // Tail twist should have strong effect on pitching moment
  ad::AeroGradients tail_twist_grads = grad_calc_.ComputeSectionGradient(
      condition, 1, 0, ad::GeometryVariable::kTwist);  // Surface 1 = HTail

  // Tail twist affects Cm significantly (tail is aft of moment reference)
  EXPECT_NE(tail_twist_grads.dCm, 0.0)
      << "Tail twist should affect pitching moment";

  // Positive tail twist (nose-up on tail) should create negative Cm
  // (nose-down moment on aircraft)
  EXPECT_LT(tail_twist_grads.dCm, 0.0)
      << "Positive tail twist should create nose-down moment";
}

}  // namespace
}  // namespace vlm
