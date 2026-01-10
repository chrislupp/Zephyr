// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include "vlm/ad/adjoint_solver.h"

#include <cmath>
#include <vector>

#include "vlm/aero/force_calculator.h"
#include "vlm/ad/enzyme_forces.h"

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

/// Extract mesh data into raw arrays for Enzyme kernels
void ExtractMeshData(
    const VortexLattice& mesh,
    const SolutionState& solution,
    const Aircraft& aircraft,
    std::vector<double>& bound_left,
    std::vector<double>& bound_right,
    std::vector<double>& control_points,
    std::vector<double>& moment_arm,
    std::vector<double>& induced_velocity) {

  const Index n = mesh.GetNumVortices();
  bound_left.resize(3 * n);
  bound_right.resize(3 * n);
  control_points.resize(3 * n);
  moment_arm.resize(3 * n);
  induced_velocity.resize(3 * n);

  Vec3 moment_ref = aircraft.GetMomentReferencePoint();

  for (Index i = 0; i < n; ++i) {
    const HorseshoeVortex& vor = mesh.GetVortex(i);

    bound_left[3*i + 0] = vor.bound_vortex_left.x();
    bound_left[3*i + 1] = vor.bound_vortex_left.y();
    bound_left[3*i + 2] = vor.bound_vortex_left.z();

    bound_right[3*i + 0] = vor.bound_vortex_right.x();
    bound_right[3*i + 1] = vor.bound_vortex_right.y();
    bound_right[3*i + 2] = vor.bound_vortex_right.z();

    control_points[3*i + 0] = vor.control_point.x();
    control_points[3*i + 1] = vor.control_point.y();
    control_points[3*i + 2] = vor.control_point.z();

    moment_arm[3*i + 0] = vor.control_point.x() - moment_ref.x();
    moment_arm[3*i + 1] = vor.control_point.y() - moment_ref.y();
    moment_arm[3*i + 2] = vor.control_point.z() - moment_ref.z();

    if (solution.induced_velocity.rows() > i) {
      induced_velocity[3*i + 0] = solution.induced_velocity(i, 0);
      induced_velocity[3*i + 1] = solution.induced_velocity(i, 1);
      induced_velocity[3*i + 2] = solution.induced_velocity(i, 2);
    } else {
      induced_velocity[3*i + 0] = 0.0;
      induced_velocity[3*i + 1] = 0.0;
      induced_velocity[3*i + 2] = 0.0;
    }
  }
}

#ifdef VLM_USE_ENZYME

/// Compute force gradients using Enzyme reverse-mode AD
ForceGradients ComputeForceGradientsEnzyme(
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

  // Extract mesh data into raw arrays
  std::vector<double> bound_left, bound_right, control_points, moment_arm, induced_velocity;
  ExtractMeshData(mesh, solution, aircraft, bound_left, bound_right,
                  control_points, moment_arm, induced_velocity);

  // Copy circulation to contiguous array
  std::vector<double> circulation(n);
  for (Index i = 0; i < n; ++i) {
    circulation[i] = solution.circulation(i);
  }

  double reference_area = aircraft.GetReferenceArea();

  // Compute gradients for each output using Enzyme reverse-mode
  // Output indices: 0=CL, 1=CD, 2=CY, 3=Cl, 4=Cm, 5=Cn
  std::vector<double> d_circulation(n);

  // CL gradient
  enzyme_force_circulation_gradient(
      condition.alpha, condition.beta,
      circulation.data(), induced_velocity.data(),
      bound_left.data(), bound_right.data(),
      control_points.data(), moment_arm.data(),
      static_cast<int>(n), reference_area, 0, d_circulation.data());
  for (Index i = 0; i < n; ++i) {
    grads.dCL_dGamma(i) = d_circulation[i];
  }

  // CD gradient
  enzyme_force_circulation_gradient(
      condition.alpha, condition.beta,
      circulation.data(), induced_velocity.data(),
      bound_left.data(), bound_right.data(),
      control_points.data(), moment_arm.data(),
      static_cast<int>(n), reference_area, 1, d_circulation.data());
  for (Index i = 0; i < n; ++i) {
    grads.dCD_dGamma(i) = d_circulation[i];
  }

  // CY gradient
  enzyme_force_circulation_gradient(
      condition.alpha, condition.beta,
      circulation.data(), induced_velocity.data(),
      bound_left.data(), bound_right.data(),
      control_points.data(), moment_arm.data(),
      static_cast<int>(n), reference_area, 2, d_circulation.data());
  for (Index i = 0; i < n; ++i) {
    grads.dCY_dGamma(i) = d_circulation[i];
  }

  // Cl gradient
  enzyme_force_circulation_gradient(
      condition.alpha, condition.beta,
      circulation.data(), induced_velocity.data(),
      bound_left.data(), bound_right.data(),
      control_points.data(), moment_arm.data(),
      static_cast<int>(n), reference_area, 3, d_circulation.data());
  for (Index i = 0; i < n; ++i) {
    grads.dCl_dGamma(i) = d_circulation[i];
  }

  // Cm gradient
  enzyme_force_circulation_gradient(
      condition.alpha, condition.beta,
      circulation.data(), induced_velocity.data(),
      bound_left.data(), bound_right.data(),
      control_points.data(), moment_arm.data(),
      static_cast<int>(n), reference_area, 4, d_circulation.data());
  for (Index i = 0; i < n; ++i) {
    grads.dCm_dGamma(i) = d_circulation[i];
  }

  // Cn gradient
  enzyme_force_circulation_gradient(
      condition.alpha, condition.beta,
      circulation.data(), induced_velocity.data(),
      bound_left.data(), bound_right.data(),
      control_points.data(), moment_arm.data(),
      static_cast<int>(n), reference_area, 5, d_circulation.data());
  for (Index i = 0; i < n; ++i) {
    grads.dCn_dGamma(i) = d_circulation[i];
  }

  return grads;
}

#endif  // VLM_USE_ENZYME

/// Compute force gradients using finite differences (fallback)
ForceGradients ComputeForceGradientsFD(
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

  // Compute gradient via finite differences for each circulation element
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

ForceGradients ComputeForceGradients(
    const VortexLattice& mesh,
    const SolutionState& solution,
    const FlightCondition& condition,
    const Aircraft& aircraft) {
#ifdef VLM_USE_ENZYME
  return ComputeForceGradientsEnzyme(mesh, solution, condition, aircraft);
#else
  return ComputeForceGradientsFD(mesh, solution, condition, aircraft);
#endif
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
