// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include "vlm/ad/adjoint_solver.h"

#include <cmath>

#include "vlm/aero/force_calculator.h"

namespace vlm {
namespace ad {

// ============================================================================
// AdjointSolver Implementation
// ============================================================================

AdjointSolver::AdjointSolver(const VortexLatticeSolver& solver)
    : solver_(&solver), is_valid_(false), num_unknowns_(0) {
  if (solver.IsAICValid()) {
    is_valid_ = true;
    num_unknowns_ = solver.GetAICMatrix().rows();
  }
}

VecX AdjointSolver::SolveAdjoint(const VecX& rhs) const {
  if (!is_valid_) {
    return VecX::Zero(rhs.size());
  }

  // Solve Aᵀ·λ = rhs
  // Since we have LU factorization of A, we solve using the transpose.
  // For PA = LU, we have Aᵀ = UᵀLᵀPᵀ
  // Eigen's PartialPivLU supports transpose solve via .transpose().solve()

  // Get the AIC matrix and solve the transpose system
  // Note: We need to re-solve since PartialPivLU doesn't store transpose
  // In production, we'd want to cache this or use a different factorization
  const MatX& A = solver_->GetAICMatrix();

  // Solve Aᵀλ = rhs by solving (Aᵀ)λ = rhs
  // This is equivalent to λ = (Aᵀ)⁻¹ rhs = (A⁻¹)ᵀ rhs
  // Using the identity: (A⁻¹)ᵀ = (Aᵀ)⁻¹

  // For efficiency, we use the LU factorization of Aᵀ
  // In the future, this could be optimized by directly using the LU factors
  Eigen::PartialPivLU<MatX> lu_transpose(A.transpose());
  return lu_transpose.solve(rhs);
}

MatX AdjointSolver::SolveAdjointMultiple(const MatX& rhs) const {
  if (!is_valid_) {
    return MatX::Zero(rhs.rows(), rhs.cols());
  }

  const MatX& A = solver_->GetAICMatrix();
  Eigen::PartialPivLU<MatX> lu_transpose(A.transpose());
  return lu_transpose.solve(rhs);
}

// ============================================================================
// Force Gradient Computation
// ============================================================================

namespace {

/// Finite difference step for gradient computation
constexpr double kFiniteDiffStep = 1e-7;

/// Compute forces for a perturbed circulation
AerodynamicResults ComputeForcesWithPerturbedCirculation(
    const VortexLattice& mesh,
    const SolutionState& base_solution,
    const FlightCondition& condition,
    const Aircraft& aircraft,
    Index perturb_index,
    double perturbation) {

  // Create perturbed solution
  SolutionState perturbed = base_solution;
  perturbed.circulation(perturb_index) += perturbation;

  // Recompute induced velocities for perturbed circulation
  // v_induced = V_influence * Gamma
  // For simplicity, we approximate by scaling the induced velocity contribution
  // A more accurate approach would recompute v_ind = sum_j V_ij * Gamma_j
  // but this is expensive. The induced velocity effect is often small.

  // For now, we'll use the base induced velocities
  // This is an approximation that works well for flight condition derivatives
  // For geometry derivatives, we'd need the full recomputation

  ForceCalculator calc;
  return calc.ComputeForces(aircraft, mesh, perturbed, condition);
}

}  // namespace

ForceGradients ComputeForceGradients(
    const VortexLattice& mesh,
    const SolutionState& solution,
    const FlightCondition& condition,
    const Aircraft& aircraft) {

  const Index n = solution.circulation.size();

  ForceGradients grads;
  grads.dCL_dGamma = VecX::Zero(n);
  grads.dCD_dGamma = VecX::Zero(n);
  grads.dCY_dGamma = VecX::Zero(n);
  grads.dCl_dGamma = VecX::Zero(n);
  grads.dCm_dGamma = VecX::Zero(n);
  grads.dCn_dGamma = VecX::Zero(n);

  if (!solution.is_valid || n == 0) {
    return grads;
  }

  // Compute base forces
  ForceCalculator calc;
  AerodynamicResults base_results =
      calc.ComputeForces(aircraft, mesh, solution, condition);

  // Compute gradient via finite differences for each circulation element
  // dJ/dGamma_i = (J(Gamma + e_i*h) - J(Gamma - e_i*h)) / (2h)

  const double h = kFiniteDiffStep;
  const double inv_2h = 0.5 / h;

  for (Index i = 0; i < n; ++i) {
    // Positive perturbation
    AerodynamicResults res_plus = ComputeForcesWithPerturbedCirculation(
        mesh, solution, condition, aircraft, i, h);

    // Negative perturbation
    AerodynamicResults res_minus = ComputeForcesWithPerturbedCirculation(
        mesh, solution, condition, aircraft, i, -h);

    // Central difference gradients
    grads.dCL_dGamma(i) = (res_plus.CL - res_minus.CL) * inv_2h;
    grads.dCD_dGamma(i) = (res_plus.CD - res_minus.CD) * inv_2h;
    grads.dCY_dGamma(i) = (res_plus.CY - res_minus.CY) * inv_2h;
    grads.dCl_dGamma(i) = (res_plus.Cl - res_minus.Cl) * inv_2h;
    grads.dCm_dGamma(i) = (res_plus.Cm - res_minus.Cm) * inv_2h;
    grads.dCn_dGamma(i) = (res_plus.Cn - res_minus.Cn) * inv_2h;
  }

  return grads;
}

// ============================================================================
// Total Derivative Computation
// ============================================================================

double ComputeTotalDerivative(
    int objective,
    const VortexLatticeSolver& solver,
    const FlightCondition& condition,
    const Aircraft& aircraft) {

  // Get the forward solution
  const SolutionState& solution = solver.GetSolution();
  const VortexLattice& mesh = solver.GetMesh();

  if (!solution.is_valid) {
    return 0.0;
  }

  // Step 1: Compute force gradients ∂J/∂Γ
  ForceGradients force_grads = ComputeForceGradients(
      mesh, solution, condition, aircraft);

  // Select the appropriate gradient based on objective
  VecX dJ_dGamma;
  switch (objective) {
    case 0: dJ_dGamma = force_grads.dCL_dGamma; break;
    case 1: dJ_dGamma = force_grads.dCD_dGamma; break;
    case 2: dJ_dGamma = force_grads.dCY_dGamma; break;
    case 3: dJ_dGamma = force_grads.dCl_dGamma; break;
    case 4: dJ_dGamma = force_grads.dCm_dGamma; break;
    case 5: dJ_dGamma = force_grads.dCn_dGamma; break;
    default: return 0.0;
  }

  // Step 2: Solve adjoint system Aᵀ·λ = ∂J/∂Γ
  AdjointSolver adjoint(solver);
  VecX lambda = adjoint.SolveAdjoint(dJ_dGamma);

  // Step 3: Compute total derivative
  // For flight condition derivatives (α, β, p, q, r):
  //   dJ/dα = ∂J/∂α|_Γ + λᵀ·(∂b/∂α)
  //
  // The term ∂J/∂α|_Γ captures the explicit dependence of forces on α
  // (e.g., body-to-stability axis transformation).
  //
  // The term λᵀ·(∂b/∂α) captures the implicit dependence through circulation.
  // Since A·Γ = b, and b depends on freestream direction:
  //   ∂b/∂α = RHS sensitivity to angle of attack
  //
  // For now, we return the full derivative computed via the adjoint.
  // The ∂b/∂α term requires differentiating the RHS construction.

  // Compute ∂b/∂α using finite differences on the RHS
  // This is a simplified version - full implementation would use AD

  const double h = 1e-6;
  FlightCondition cond_plus = condition;
  FlightCondition cond_minus = condition;
  cond_plus.alpha += h;
  cond_minus.alpha -= h;

  // Get unit responses and compute how they change with alpha
  // For the unit freestream response, the RHS is -n·V_inf
  // V_inf = (cos(α)cos(β), sin(β), sin(α)cos(β))
  // ∂V_inf/∂α = (-sin(α)cos(β), 0, cos(α)cos(β))

  Vec3 v_inf_plus = VortexLatticeSolver::ComputeFreestreamDirection(cond_plus);
  Vec3 v_inf_minus = VortexLatticeSolver::ComputeFreestreamDirection(cond_minus);
  Vec3 dV_dalpha = (v_inf_plus - v_inf_minus) / (2 * h);

  // The RHS for the VLM system is: b_i = -n_i · V_inf
  // So: ∂b_i/∂α = -n_i · (∂V_inf/∂α)

  const Index n = solution.circulation.size();
  VecX db_dalpha(n);
  for (Index i = 0; i < n; ++i) {
    Vec3 normal = mesh.GetVortex(i).normal_at_control_point;
    db_dalpha(i) = -normal.dot(dV_dalpha);
  }

  // Total derivative: dJ/dα = λᵀ · (∂b/∂α)
  // Note: The explicit term ∂J/∂α|_Γ is captured in the force calculation
  // through the body-to-stability axis transformation, which depends on α.
  // For a complete treatment, we'd need to add that term separately.

  double implicit_term = lambda.dot(db_dalpha);

  // For now, return just the implicit term
  // The explicit term would require differentiating BodyToStabilityAxis
  return implicit_term;
}

}  // namespace ad
}  // namespace vlm
