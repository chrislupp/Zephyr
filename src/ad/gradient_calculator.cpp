// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include "vlm/ad/gradient_calculator.h"

#include <cmath>

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

  /// Compute aerodynamic results at a flight condition
  AerodynamicResults ComputeResults(const FlightCondition& condition) {
    SolutionState solution = solver_.Solve(condition);
    return force_calc_.ComputeForces(aircraft_, solver_.GetMesh(),
                                     solution, condition);
  }

  /// Compute derivatives using central finite differences
  AeroGradients FiniteDiffDerivatives(
      const FlightCondition& base_condition,
      int param_index) {  // 0=alpha, 1=beta, 2=p, 3=q, 4=r

    FlightCondition cond_plus = base_condition;
    FlightCondition cond_minus = base_condition;

    // Perturb the appropriate parameter
    switch (param_index) {
      case 0:  // alpha
        cond_plus.alpha += fd_step_;
        cond_minus.alpha -= fd_step_;
        break;
      case 1:  // beta
        cond_plus.beta += fd_step_;
        cond_minus.beta -= fd_step_;
        break;
      case 2:  // roll rate
        cond_plus.roll_rate += fd_step_;
        cond_minus.roll_rate -= fd_step_;
        break;
      case 3:  // pitch rate
        cond_plus.pitch_rate += fd_step_;
        cond_minus.pitch_rate -= fd_step_;
        break;
      case 4:  // yaw rate
        cond_plus.yaw_rate += fd_step_;
        cond_minus.yaw_rate -= fd_step_;
        break;
    }

    // Compute results at perturbed conditions
    AerodynamicResults res_plus = ComputeResults(cond_plus);
    AerodynamicResults res_minus = ComputeResults(cond_minus);

    // Central difference
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
// Enzyme Implementation
// ============================================================================

#ifdef VLM_USE_ENZYME

// The Enzyme-based implementation will be added here when we have verified
// the basic infrastructure works. For now, we use finite differences even
// when Enzyme is available, to allow incremental testing.

// TODO: Implement enzyme_compute_forces_kernel and use __enzyme_autodiff
// to compute gradients analytically.

#endif  // VLM_USE_ENZYME

// ============================================================================
// Derivative Computation Methods
// ============================================================================

AeroGradients GradientCalculator::ComputeAlphaDerivatives(
    const FlightCondition& condition) {
  // TODO: When Enzyme is working, use AD here instead of finite differences
  return impl_->FiniteDiffDerivatives(condition, 0);
}

AeroGradients GradientCalculator::ComputeBetaDerivatives(
    const FlightCondition& condition) {
  return impl_->FiniteDiffDerivatives(condition, 1);
}

AeroGradients GradientCalculator::ComputeRollRateDerivatives(
    const FlightCondition& condition) {
  return impl_->FiniteDiffDerivatives(condition, 2);
}

AeroGradients GradientCalculator::ComputePitchRateDerivatives(
    const FlightCondition& condition) {
  return impl_->FiniteDiffDerivatives(condition, 3);
}

AeroGradients GradientCalculator::ComputeYawRateDerivatives(
    const FlightCondition& condition) {
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
