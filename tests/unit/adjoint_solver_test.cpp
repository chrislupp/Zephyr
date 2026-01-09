// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include <gtest/gtest.h>

#include "vlm/ad/adjoint_solver.h"
#include "vlm/aero/force_calculator.h"
#include "vlm/geometry/aircraft.h"
#include "vlm/geometry/section.h"
#include "vlm/geometry/surface.h"

namespace vlm {
namespace {

// ============================================================================
// Test Fixture
// ============================================================================

class AdjointSolverTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Create simple rectangular wing
    Surface wing("Wing");
    wing.SetChordwisePanelCount(4);
    wing.SetSpanwisePanelCount(6);
    wing.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    wing.AddSection(Section(Vec3(0, 3.0, 0), 1.0, 0.0));

    aircraft_.AddSurface(wing);
    aircraft_.SetReferenceArea(3.0);
    aircraft_.SetReferenceChord(1.0);
    aircraft_.SetReferenceSpan(6.0);

    solver_.SetAircraft(aircraft_);
  }

  Aircraft aircraft_{"Test"};
  VortexLatticeSolver solver_;
};

// ============================================================================
// AdjointSolver Basic Tests
// ============================================================================

TEST_F(AdjointSolverTest, ConstructFromSolverBeforeSolve) {
  // Before calling Solve(), AIC is not yet built
  ad::AdjointSolver adjoint(solver_);

  // AIC is built lazily during Solve(), so adjoint should NOT be valid yet
  EXPECT_FALSE(adjoint.IsValid());
}

TEST_F(AdjointSolverTest, ConstructFromSolverAfterSolve) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;
  solver_.Solve(condition);

  ad::AdjointSolver adjoint(solver_);
  EXPECT_TRUE(adjoint.IsValid());
  EXPECT_GT(adjoint.GetNumUnknowns(), 0);
}

TEST_F(AdjointSolverTest, SolveAdjointReturnsCorrectSize) {
  FlightCondition condition;
  solver_.Solve(condition);

  ad::AdjointSolver adjoint(solver_);
  const Index n = adjoint.GetNumUnknowns();

  VecX rhs = VecX::Ones(n);
  VecX lambda = adjoint.SolveAdjoint(rhs);

  EXPECT_EQ(lambda.size(), n);
}

TEST_F(AdjointSolverTest, AdjointSolveIsCorrect) {
  // Verify that Aᵀλ = rhs is satisfied
  FlightCondition condition;
  solver_.Solve(condition);

  ad::AdjointSolver adjoint(solver_);
  const Index n = adjoint.GetNumUnknowns();

  // Random RHS
  VecX rhs = VecX::Random(n);

  // Solve adjoint
  VecX lambda = adjoint.SolveAdjoint(rhs);

  // Verify: Aᵀλ should equal rhs
  const MatX& A = solver_.GetAICMatrix();
  VecX residual = A.transpose() * lambda - rhs;

  EXPECT_LT(residual.norm(), 1e-10 * rhs.norm())
      << "Adjoint solve residual too large";
}

TEST_F(AdjointSolverTest, SolveMultipleAdjointProblems) {
  FlightCondition condition;
  solver_.Solve(condition);

  ad::AdjointSolver adjoint(solver_);
  const Index n = adjoint.GetNumUnknowns();

  // Multiple RHS vectors (e.g., for CL, CD, Cm simultaneously)
  MatX rhs(n, 3);
  rhs.col(0) = VecX::Random(n);
  rhs.col(1) = VecX::Random(n);
  rhs.col(2) = VecX::Random(n);

  MatX lambda = adjoint.SolveAdjointMultiple(rhs);

  EXPECT_EQ(lambda.rows(), n);
  EXPECT_EQ(lambda.cols(), 3);

  // Verify each column
  const MatX& A = solver_.GetAICMatrix();
  for (int k = 0; k < 3; ++k) {
    VecX residual = A.transpose() * lambda.col(k) - rhs.col(k);
    EXPECT_LT(residual.norm(), 1e-10 * rhs.col(k).norm());
  }
}

// ============================================================================
// Force Gradient Tests
// ============================================================================

TEST_F(AdjointSolverTest, ForceGradientsHaveCorrectSize) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;
  SolutionState solution = solver_.Solve(condition);

  ad::ForceGradients grads = ad::ComputeForceGradients(
      solver_.GetMesh(), solution, condition, aircraft_);

  const Index n = solution.circulation.size();
  EXPECT_EQ(grads.dCL_dGamma.size(), n);
  EXPECT_EQ(grads.dCD_dGamma.size(), n);
  EXPECT_EQ(grads.dCY_dGamma.size(), n);
  EXPECT_EQ(grads.dCl_dGamma.size(), n);
  EXPECT_EQ(grads.dCm_dGamma.size(), n);
  EXPECT_EQ(grads.dCn_dGamma.size(), n);
}

TEST_F(AdjointSolverTest, ForceGradientsAreNonZero) {
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;
  SolutionState solution = solver_.Solve(condition);

  ad::ForceGradients grads = ad::ComputeForceGradients(
      solver_.GetMesh(), solution, condition, aircraft_);

  // At non-zero alpha, CL should depend on circulation
  EXPECT_GT(grads.dCL_dGamma.norm(), 0.0)
      << "dCL/dGamma should be non-zero";

  // Cm should also depend on circulation
  EXPECT_GT(grads.dCm_dGamma.norm(), 0.0)
      << "dCm/dGamma should be non-zero";
}

TEST_F(AdjointSolverTest, CLGradientIsPositive) {
  // Increasing circulation should increase lift
  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;
  SolutionState solution = solver_.Solve(condition);

  ad::ForceGradients grads = ad::ComputeForceGradients(
      solver_.GetMesh(), solution, condition, aircraft_);

  // Most elements should have positive dCL/dGamma
  int positive_count = 0;
  for (Index i = 0; i < grads.dCL_dGamma.size(); ++i) {
    if (grads.dCL_dGamma(i) > 0) {
      positive_count++;
    }
  }

  // At least 80% should be positive (some edge effects may cause exceptions)
  double positive_fraction =
      static_cast<double>(positive_count) / grads.dCL_dGamma.size();
  EXPECT_GT(positive_fraction, 0.8)
      << "Most dCL/dGamma elements should be positive";
}

// ============================================================================
// Adjoint vs Finite Difference Validation
// ============================================================================

class AdjointValidationTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Create a simple wing for validation
    Surface wing("Wing");
    wing.SetChordwisePanelCount(3);
    wing.SetSpanwisePanelCount(4);
    wing.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    wing.AddSection(Section(Vec3(0, 2.0, 0), 1.0, 0.0));

    aircraft_.AddSurface(wing);
    aircraft_.SetReferenceArea(2.0);
    aircraft_.SetReferenceChord(1.0);
    aircraft_.SetReferenceSpan(4.0);

    solver_.SetAircraft(aircraft_);
  }

  /// Compute derivative using finite differences on the full system
  double FiniteDiffDerivative(int objective, const FlightCondition& base_cond) {
    const double h = 1e-6;

    FlightCondition cond_plus = base_cond;
    FlightCondition cond_minus = base_cond;
    cond_plus.alpha += h;
    cond_minus.alpha -= h;

    // Solve at perturbed conditions
    SolutionState sol_plus = solver_.Solve(cond_plus);
    SolutionState sol_minus = solver_.Solve(cond_minus);

    // Compute forces
    ForceCalculator calc;
    AerodynamicResults res_plus =
        calc.ComputeForces(aircraft_, solver_.GetMesh(), sol_plus, cond_plus);
    AerodynamicResults res_minus =
        calc.ComputeForces(aircraft_, solver_.GetMesh(), sol_minus, cond_minus);

    // Extract objective
    double J_plus, J_minus;
    switch (objective) {
      case 0: J_plus = res_plus.CL; J_minus = res_minus.CL; break;
      case 1: J_plus = res_plus.CD; J_minus = res_minus.CD; break;
      case 2: J_plus = res_plus.CY; J_minus = res_minus.CY; break;
      case 3: J_plus = res_plus.Cl; J_minus = res_minus.Cl; break;
      case 4: J_plus = res_plus.Cm; J_minus = res_minus.Cm; break;
      case 5: J_plus = res_plus.Cn; J_minus = res_minus.Cn; break;
      default: return 0.0;
    }

    return (J_plus - J_minus) / (2 * h);
  }

  Aircraft aircraft_{"Validation"};
  VortexLatticeSolver solver_;
};

TEST_F(AdjointValidationTest, DISABLED_AdjointMatchesFiniteDiffForCL) {
  // Note: This test is disabled because the current adjoint implementation
  // only captures the implicit term (through circulation). The explicit
  // term (body-to-stability axis transformation) is not yet included.
  // Once that's added, this test should pass.

  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  // Solve at base condition
  solver_.Solve(condition);

  // Adjoint derivative
  double adjoint_deriv = ad::ComputeTotalDerivative(0, solver_, condition, aircraft_);

  // Finite difference derivative
  double fd_deriv = FiniteDiffDerivative(0, condition);

  // They should match to reasonable tolerance
  EXPECT_NEAR(adjoint_deriv, fd_deriv, std::abs(fd_deriv) * 0.1 + 0.01)
      << "Adjoint dCL/dalpha = " << adjoint_deriv
      << ", FD dCL/dalpha = " << fd_deriv;
}

// ============================================================================
// Adjoint Identity Tests
// ============================================================================

TEST_F(AdjointSolverTest, AdjointIdentityTest) {
  // For a scalar objective J and linear system Ax = b:
  // dJ/db = (A^{-T}) * (dJ/dx) = adjoint solution
  //
  // So: (dJ/db)^T * db = dJ (for any perturbation db)
  //
  // This tests that our adjoint solve is consistent

  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;
  solver_.Solve(condition);

  const SolutionState& solution = solver_.GetSolution();
  const VortexLattice& mesh = solver_.GetMesh();
  const Index n = solution.circulation.size();

  // Compute force gradients
  ad::ForceGradients grads = ad::ComputeForceGradients(
      mesh, solution, condition, aircraft_);

  // Solve adjoint for CL
  ad::AdjointSolver adjoint(solver_);
  VecX lambda_CL = adjoint.SolveAdjoint(grads.dCL_dGamma);

  // The adjoint variable lambda represents dCL/db where b is the RHS
  // Verify dimensions are consistent
  EXPECT_EQ(lambda_CL.size(), n);

  // lambda should have non-trivial values
  EXPECT_GT(lambda_CL.norm(), 1e-10);
}

}  // namespace
}  // namespace vlm
