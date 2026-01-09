// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include "vlm/ad/gradient_calculator.h"

#include <cmath>

#include "vlm/ad/adjoint_solver.h"
#include "vlm/aero/force_calculator.h"

namespace vlm {
namespace ad {

// ============================================================================
// Implementation Class
// ============================================================================

class GradientCalculator::Impl {
 public:
  Aircraft aircraft_{"GradientCalc"};
  VortexLatticeSolver solver_;
  ForceCalculator force_calc_;
  bool aircraft_set_{false};
  Scalar fd_step_{1e-6};
  bool use_adjoint_{true};  // Use adjoint method by default

  // Cached solution for adjoint computations
  SolutionState cached_solution_;
  FlightCondition cached_condition_;
  bool solution_cached_{false};

  /// Solve and cache the solution at a flight condition
  void SolveAndCache(const FlightCondition& condition) {
    cached_solution_ = solver_.Solve(condition);
    cached_condition_ = condition;
    solution_cached_ = true;
  }

  /// Compute aerodynamic results at a flight condition (full solve)
  AerodynamicResults ComputeResults(const FlightCondition& condition) {
    SolutionState solution = solver_.Solve(condition);
    return force_calc_.ComputeForces(aircraft_, solver_.GetMesh(),
                                     solution, condition);
  }

  /// Compute aerodynamic results using cached circulation but different condition
  /// This gives us ∂J/∂condition at fixed Γ (the explicit term)
  AerodynamicResults ComputeResultsExplicit(const FlightCondition& condition) {
    // Use cached circulation, but evaluate forces at new condition
    // This captures how forces change when only the axis transformation changes
    return force_calc_.ComputeForces(aircraft_, solver_.GetMesh(),
                                     cached_solution_, condition);
  }

  /// Compute RHS sensitivity ∂b/∂(parameter)
  /// The RHS is b_i = -n_i · V_eff where V_eff depends on flight condition
  VecX ComputeRHSSensitivity(const FlightCondition& condition, int param_index) {
    const VortexLattice& mesh = solver_.GetMesh();
    const Index n = cached_solution_.circulation.size();
    VecX db_dparam(n);

    const Scalar h = 1e-7;

    // Perturbed conditions
    FlightCondition cond_plus = condition;
    FlightCondition cond_minus = condition;

    switch (param_index) {
      case 0:  // alpha
        cond_plus.alpha += h;
        cond_minus.alpha -= h;
        break;
      case 1:  // beta
        cond_plus.beta += h;
        cond_minus.beta -= h;
        break;
      case 2:  // roll rate
        cond_plus.roll_rate += h;
        cond_minus.roll_rate -= h;
        break;
      case 3:  // pitch rate
        cond_plus.pitch_rate += h;
        cond_minus.pitch_rate -= h;
        break;
      case 4:  // yaw rate
        cond_plus.yaw_rate += h;
        cond_minus.yaw_rate -= h;
        break;
    }

    // Compute freestream directions
    Vec3 v_plus = VortexLatticeSolver::ComputeFreestreamDirection(cond_plus);
    Vec3 v_minus = VortexLatticeSolver::ComputeFreestreamDirection(cond_minus);
    Vec3 dV_dparam = (v_plus - v_minus) / (2 * h);

    // RHS sensitivity: ∂b_i/∂param = -n_i · (∂V/∂param)
    for (Index i = 0; i < n; ++i) {
      Vec3 normal = mesh.GetVortex(i).normal_at_control_point;
      db_dparam(i) = -normal.dot(dV_dparam);
    }

    // For rotation rates, there's also a contribution from the rotation velocity
    // The RHS includes -n · (ω × r) where r is the control point position
    if (param_index >= 2 && param_index <= 4) {
      const Vec3& ref = solver_.GetReferencePoint();

      for (Index i = 0; i < n; ++i) {
        Vec3 normal = mesh.GetVortex(i).normal_at_control_point;
        Vec3 cp = mesh.GetVortex(i).control_point;
        Vec3 r = cp - ref;

        Vec3 omega_dir;
        switch (param_index) {
          case 2: omega_dir = Vec3(1, 0, 0); break;  // p (roll)
          case 3: omega_dir = Vec3(0, 1, 0); break;  // q (pitch)
          case 4: omega_dir = Vec3(0, 0, 1); break;  // r (yaw)
        }

        // Rotation velocity: V_rot = ω × r
        // ∂V_rot/∂ω = (direction) × r
        Vec3 dVrot = omega_dir.cross(r);
        db_dparam(i) -= normal.dot(dVrot);
      }
    }

    return db_dparam;
  }

  /// Compute explicit term ∂J/∂param at fixed circulation
  /// This uses finite diff on forces only (no re-solve)
  AeroGradients ComputeExplicitDerivatives(int param_index) {
    const Scalar h = fd_step_;
    FlightCondition cond_plus = cached_condition_;
    FlightCondition cond_minus = cached_condition_;

    switch (param_index) {
      case 0:
        cond_plus.alpha += h;
        cond_minus.alpha -= h;
        break;
      case 1:
        cond_plus.beta += h;
        cond_minus.beta -= h;
        break;
      case 2:
        cond_plus.roll_rate += h;
        cond_minus.roll_rate -= h;
        break;
      case 3:
        cond_plus.pitch_rate += h;
        cond_minus.pitch_rate -= h;
        break;
      case 4:
        cond_plus.yaw_rate += h;
        cond_minus.yaw_rate -= h;
        break;
    }

    // Compute forces at perturbed conditions using SAME circulation
    AerodynamicResults res_plus = ComputeResultsExplicit(cond_plus);
    AerodynamicResults res_minus = ComputeResultsExplicit(cond_minus);

    Scalar inv_2h = 0.5 / h;

    AeroGradients grads;
    grads.dCL = (res_plus.CL - res_minus.CL) * inv_2h;
    grads.dCD = (res_plus.CD - res_minus.CD) * inv_2h;
    grads.dCY = (res_plus.CY - res_minus.CY) * inv_2h;
    grads.dCl = (res_plus.Cl - res_minus.Cl) * inv_2h;
    grads.dCm = (res_plus.Cm - res_minus.Cm) * inv_2h;
    grads.dCn = (res_plus.Cn - res_minus.Cn) * inv_2h;
    grads.dCDi = (res_plus.CDi_ff - res_minus.CDi_ff) * inv_2h;

    return grads;
  }

  /// Compute derivatives using adjoint method
  /// dJ/d(param) = ∂J/∂param|_Γ + λᵀ·(∂b/∂param)
  AeroGradients AdjointDerivatives(const FlightCondition& condition,
                                   int param_index) {
    // Step 1: Solve at base condition and cache
    if (!solution_cached_ ||
        condition.alpha != cached_condition_.alpha ||
        condition.beta != cached_condition_.beta) {
      SolveAndCache(condition);
    }

    // Step 2: Compute force gradients ∂J/∂Γ for all outputs
    ForceGradients force_grads = ComputeForceGradients(
        solver_.GetMesh(), cached_solution_, condition, aircraft_);

    // Step 3: Solve adjoint for each output
    AdjointSolver adjoint(solver_);
    if (!adjoint.IsValid()) {
      // Fall back to finite differences
      return FiniteDiffDerivatives(condition, param_index);
    }

    // Solve all adjoints at once
    const Index n = cached_solution_.circulation.size();
    MatX rhs(n, 6);  // CL, CD, CY, Cl, Cm, Cn
    rhs.col(0) = force_grads.dCL_dGamma;
    rhs.col(1) = force_grads.dCD_dGamma;
    rhs.col(2) = force_grads.dCY_dGamma;
    rhs.col(3) = force_grads.dCl_dGamma;
    rhs.col(4) = force_grads.dCm_dGamma;
    rhs.col(5) = force_grads.dCn_dGamma;

    MatX lambda = adjoint.SolveAdjointMultiple(rhs);

    // Step 4: Compute RHS sensitivity ∂b/∂param
    VecX db_dparam = ComputeRHSSensitivity(condition, param_index);

    // Step 5: Compute implicit terms λᵀ·(∂b/∂param)
    AeroGradients implicit_grads;
    implicit_grads.dCL = lambda.col(0).dot(db_dparam);
    implicit_grads.dCD = lambda.col(1).dot(db_dparam);
    implicit_grads.dCY = lambda.col(2).dot(db_dparam);
    implicit_grads.dCl = lambda.col(3).dot(db_dparam);
    implicit_grads.dCm = lambda.col(4).dot(db_dparam);
    implicit_grads.dCn = lambda.col(5).dot(db_dparam);

    // Step 6: Compute explicit terms ∂J/∂param|_Γ
    AeroGradients explicit_grads = ComputeExplicitDerivatives(param_index);

    // Step 7: Total derivative = explicit + implicit
    AeroGradients total;
    total.dCL = explicit_grads.dCL + implicit_grads.dCL;
    total.dCD = explicit_grads.dCD + implicit_grads.dCD;
    total.dCY = explicit_grads.dCY + implicit_grads.dCY;
    total.dCl = explicit_grads.dCl + implicit_grads.dCl;
    total.dCm = explicit_grads.dCm + implicit_grads.dCm;
    total.dCn = explicit_grads.dCn + implicit_grads.dCn;
    total.dCDi = explicit_grads.dCDi;  // CDi is computed differently

    return total;
  }

  /// Compute derivatives using central finite differences (original method)
  AeroGradients FiniteDiffDerivatives(const FlightCondition& base_condition,
                                      int param_index) {
    FlightCondition cond_plus = base_condition;
    FlightCondition cond_minus = base_condition;

    switch (param_index) {
      case 0:
        cond_plus.alpha += fd_step_;
        cond_minus.alpha -= fd_step_;
        break;
      case 1:
        cond_plus.beta += fd_step_;
        cond_minus.beta -= fd_step_;
        break;
      case 2:
        cond_plus.roll_rate += fd_step_;
        cond_minus.roll_rate -= fd_step_;
        break;
      case 3:
        cond_plus.pitch_rate += fd_step_;
        cond_minus.pitch_rate -= fd_step_;
        break;
      case 4:
        cond_plus.yaw_rate += fd_step_;
        cond_minus.yaw_rate -= fd_step_;
        break;
    }

    AerodynamicResults res_plus = ComputeResults(cond_plus);
    AerodynamicResults res_minus = ComputeResults(cond_minus);

    Scalar inv_2step = 0.5 / fd_step_;

    AeroGradients grads;
    grads.dCL = (res_plus.CL - res_minus.CL) * inv_2step;
    grads.dCD = (res_plus.CD - res_minus.CD) * inv_2step;
    grads.dCY = (res_plus.CY - res_minus.CY) * inv_2step;
    grads.dCl = (res_plus.Cl - res_minus.Cl) * inv_2step;
    grads.dCm = (res_plus.Cm - res_minus.Cm) * inv_2step;
    grads.dCn = (res_plus.Cn - res_minus.Cn) * inv_2step;
    grads.dCDi = (res_plus.CDi_ff - res_minus.CDi_ff) * inv_2step;

    return grads;
  }
};

// ============================================================================
// GradientCalculator Public Interface
// ============================================================================

GradientCalculator::GradientCalculator() : impl_(std::make_unique<Impl>()) {}

GradientCalculator::~GradientCalculator() = default;

void GradientCalculator::SetAircraft(const Aircraft& aircraft) {
  impl_->aircraft_ = aircraft;
  impl_->solver_.SetAircraft(aircraft);
  impl_->aircraft_set_ = true;
  impl_->solution_cached_ = false;  // Invalidate cache
}

void GradientCalculator::SetReferencePoint(const Vec3& ref) {
  impl_->solver_.SetReferencePoint(ref);
}

bool GradientCalculator::IsUsingEnzyme() {
#ifdef VLM_USE_ENZYME
  return true;
#else
  return false;
#endif
}

void GradientCalculator::SetFiniteDiffStep(Scalar step) {
  fd_step_ = step;
  impl_->fd_step_ = step;
}

// ============================================================================
// Derivative Computation Methods
// ============================================================================

AeroGradients GradientCalculator::ComputeAlphaDerivatives(
    const FlightCondition& condition) {
  if (impl_->use_adjoint_) {
    return impl_->AdjointDerivatives(condition, 0);
  }
  return impl_->FiniteDiffDerivatives(condition, 0);
}

AeroGradients GradientCalculator::ComputeBetaDerivatives(
    const FlightCondition& condition) {
  if (impl_->use_adjoint_) {
    return impl_->AdjointDerivatives(condition, 1);
  }
  return impl_->FiniteDiffDerivatives(condition, 1);
}

AeroGradients GradientCalculator::ComputeRollRateDerivatives(
    const FlightCondition& condition) {
  if (impl_->use_adjoint_) {
    return impl_->AdjointDerivatives(condition, 2);
  }
  return impl_->FiniteDiffDerivatives(condition, 2);
}

AeroGradients GradientCalculator::ComputePitchRateDerivatives(
    const FlightCondition& condition) {
  if (impl_->use_adjoint_) {
    return impl_->AdjointDerivatives(condition, 3);
  }
  return impl_->FiniteDiffDerivatives(condition, 3);
}

AeroGradients GradientCalculator::ComputeYawRateDerivatives(
    const FlightCondition& condition) {
  if (impl_->use_adjoint_) {
    return impl_->AdjointDerivatives(condition, 4);
  }
  return impl_->FiniteDiffDerivatives(condition, 4);
}

FullGradients GradientCalculator::ComputeFlightConditionDerivatives(
    const FlightCondition& condition) {
  FullGradients grads;

  grads.alpha_derivs = ComputeAlphaDerivatives(condition);
  grads.beta_derivs = ComputeBetaDerivatives(condition);
  grads.p_derivs = ComputeRollRateDerivatives(condition);
  grads.q_derivs = ComputePitchRateDerivatives(condition);
  grads.r_derivs = ComputeYawRateDerivatives(condition);

  return grads;
}

}  // namespace ad
}  // namespace vlm
