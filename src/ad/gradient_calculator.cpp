// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include "vlm/ad/gradient_calculator.h"

#include <cmath>
#include <vector>

#include "vlm/ad/adjoint_solver.h"
#include "vlm/ad/enzyme_forces.h"
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

    // Compute freestream direction gradient analytically
    Vec3 dV_dparam = Vec3::Zero();

    switch (param_index) {
      case 0: {  // alpha - analytical gradient
        double dv_dalpha[3];
        compute_freestream_alpha_gradient(condition.alpha, condition.beta, dv_dalpha);
        dV_dparam = Vec3(dv_dalpha[0], dv_dalpha[1], dv_dalpha[2]);
        break;
      }
      case 1: {  // beta - analytical gradient
        double dv_dbeta[3];
        compute_freestream_beta_gradient(condition.alpha, condition.beta, dv_dbeta);
        dV_dparam = Vec3(dv_dbeta[0], dv_dbeta[1], dv_dbeta[2]);
        break;
      }
      case 2:  // roll rate - freestream direction doesn't depend on p
      case 3:  // pitch rate - freestream direction doesn't depend on q
      case 4:  // yaw rate - freestream direction doesn't depend on r
        // dV_dparam stays zero; the rotation contribution is handled below
        break;
    }

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

  /// Extract mesh data into raw arrays for Enzyme kernels
  void ExtractMeshData(
      std::vector<double>& bound_left,
      std::vector<double>& bound_right,
      std::vector<double>& control_points,
      std::vector<double>& moment_arm,
      std::vector<double>& induced_velocity) const {

    const VortexLattice& mesh = solver_.GetMesh();
    const Index n = mesh.GetNumVortices();
    bound_left.resize(3 * n);
    bound_right.resize(3 * n);
    control_points.resize(3 * n);
    moment_arm.resize(3 * n);
    induced_velocity.resize(3 * n);

    Vec3 moment_ref = aircraft_.GetMomentReferencePoint();

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

      if (cached_solution_.induced_velocity.rows() > i) {
        induced_velocity[3*i + 0] = cached_solution_.induced_velocity(i, 0);
        induced_velocity[3*i + 1] = cached_solution_.induced_velocity(i, 1);
        induced_velocity[3*i + 2] = cached_solution_.induced_velocity(i, 2);
      } else {
        induced_velocity[3*i + 0] = 0.0;
        induced_velocity[3*i + 1] = 0.0;
        induced_velocity[3*i + 2] = 0.0;
      }
    }
  }

  /// Compute explicit term ∂J/∂param at fixed circulation using FD
  AeroGradients ComputeExplicitDerivativesFD(int param_index) {
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

#ifdef VLM_USE_ENZYME
  /// Compute explicit term ∂J/∂param at fixed circulation using Enzyme
  AeroGradients ComputeExplicitDerivativesEnzyme(int param_index) {
    // Enzyme only handles alpha and beta directly
    // For rotation rates, fall back to FD
    if (param_index >= 2) {
      return ComputeExplicitDerivativesFD(param_index);
    }

    const Index n = cached_solution_.circulation.size();
    if (n == 0) {
      return AeroGradients{};
    }

    // Extract mesh data
    std::vector<double> bound_left, bound_right, control_points, moment_arm, induced_velocity;
    ExtractMeshData(bound_left, bound_right, control_points, moment_arm, induced_velocity);

    // Copy circulation to contiguous array
    std::vector<double> circulation(n);
    for (Index i = 0; i < n; ++i) {
      circulation[i] = cached_solution_.circulation(i);
    }

    double reference_area = aircraft_.GetReferenceArea();
    double d_result[6];

    AeroGradients grads;

    if (param_index == 0) {
      // Alpha derivative
      enzyme_force_alpha_gradient(
          cached_condition_.alpha, cached_condition_.beta,
          circulation.data(), induced_velocity.data(),
          bound_left.data(), bound_right.data(),
          control_points.data(), moment_arm.data(),
          static_cast<int>(n), reference_area, d_result);
    } else {
      // Beta derivative
      enzyme_force_beta_gradient(
          cached_condition_.alpha, cached_condition_.beta,
          circulation.data(), induced_velocity.data(),
          bound_left.data(), bound_right.data(),
          control_points.data(), moment_arm.data(),
          static_cast<int>(n), reference_area, d_result);
    }

    grads.dCL = d_result[0];
    grads.dCD = d_result[1];
    grads.dCY = d_result[2];
    grads.dCl = d_result[3];
    grads.dCm = d_result[4];
    grads.dCn = d_result[5];

    // CDi_ff (Trefftz plane drag) still needs FD - it's not in the kernel
    // Fall back to FD for CDi_ff
    AeroGradients fd_grads = ComputeExplicitDerivativesFD(param_index);
    grads.dCDi = fd_grads.dCDi;

    return grads;
  }
#endif

  /// Compute explicit term ∂J/∂param at fixed circulation
  AeroGradients ComputeExplicitDerivatives(int param_index) {
#ifdef VLM_USE_ENZYME
    return ComputeExplicitDerivativesEnzyme(param_index);
#else
    return ComputeExplicitDerivativesFD(param_index);
#endif
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

AeroGradients GradientCalculator::ComputeMachDerivatives(
    const FlightCondition& condition) {
  // Mach changes the Prandtl-Glauert factor, which affects the AIC matrix.
  // We must re-solve at perturbed Mach numbers (can't use adjoint shortcut).
  const Scalar h = 0.001;  // Small Mach perturbation

  FlightCondition cond_plus = condition;
  FlightCondition cond_minus = condition;
  cond_plus.mach = condition.mach + h;
  cond_minus.mach = condition.mach - h;

  // Clamp Mach to valid range
  if (cond_plus.mach >= 1.0) {
    cond_plus.mach = 0.999;
  }
  if (cond_minus.mach < 0.0) {
    cond_minus.mach = 0.0;
  }

  AerodynamicResults res_plus = impl_->ComputeResults(cond_plus);
  AerodynamicResults res_minus = impl_->ComputeResults(cond_minus);

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

FullGradients GradientCalculator::ComputeFlightConditionDerivatives(
    const FlightCondition& condition) {
  FullGradients grads;

  grads.alpha_derivs = ComputeAlphaDerivatives(condition);
  grads.beta_derivs = ComputeBetaDerivatives(condition);
  grads.p_derivs = ComputeRollRateDerivatives(condition);
  grads.q_derivs = ComputePitchRateDerivatives(condition);
  grads.r_derivs = ComputeYawRateDerivatives(condition);
  grads.mach_derivs = ComputeMachDerivatives(condition);

  return grads;
}

// ============================================================================
// Geometry Gradient Computation
// ============================================================================

namespace {

/// Create a perturbed copy of an aircraft with one geometry parameter changed
Aircraft PerturbAircraft(const Aircraft& base, Index surface_idx,
                         Index section_idx, GeometryVariable variable,
                         Scalar perturbation) {
  Aircraft perturbed = base;  // Deep copy

  if (surface_idx >= static_cast<Index>(perturbed.GetNumSurfaces())) {
    return perturbed;
  }

  Surface& surface = perturbed.GetSurfacesMutable()[surface_idx];
  if (section_idx >= static_cast<Index>(surface.GetNumSections())) {
    return perturbed;
  }

  Section& section = surface.GetSectionMutable(section_idx);

  switch (variable) {
    case GeometryVariable::kChord:
      section.SetChord(section.GetChord() + perturbation);
      break;
    case GeometryVariable::kTwist:
      section.SetTwistAngle(section.GetTwistAngle() + perturbation);
      break;
    case GeometryVariable::kXle: {
      Vec3 pos = section.GetLeadingEdgePosition();
      pos.x() += perturbation;
      section.SetLeadingEdgePosition(pos);
      break;
    }
    case GeometryVariable::kYle: {
      Vec3 pos = section.GetLeadingEdgePosition();
      pos.y() += perturbation;
      section.SetLeadingEdgePosition(pos);
      break;
    }
    case GeometryVariable::kZle: {
      Vec3 pos = section.GetLeadingEdgePosition();
      pos.z() += perturbation;
      section.SetLeadingEdgePosition(pos);
      break;
    }
  }

  return perturbed;
}

/// Compute forces for a perturbed aircraft configuration
AerodynamicResults ComputeResultsForAircraft(const Aircraft& aircraft,
                                              const FlightCondition& condition) {
  VortexLatticeSolver solver;
  ForceCalculator force_calc;
  solver.SetAircraft(aircraft);
  SolutionState solution = solver.Solve(condition);
  return force_calc.ComputeForces(aircraft, solver.GetMesh(), solution, condition);
}

/// Compute gradients for one geometry variable at one section using FD
AeroGradients ComputeSectionGradientFD(const Aircraft& aircraft,
                                       const FlightCondition& condition,
                                       Index surface_idx,
                                       Index section_idx,
                                       GeometryVariable variable,
                                       Scalar fd_step) {
  // Step size depends on variable type
  Scalar h = fd_step;
  if (variable == GeometryVariable::kTwist) {
    // Twist is stored in degrees. Use 0.1 degree step.
    // Note: For wings in the XY plane, twist currently only affects
    // stored values, not actual panel geometry (known limitation).
    h = 0.1;  // degrees
  }

  // Forward and backward perturbations
  Aircraft aircraft_plus = PerturbAircraft(aircraft, surface_idx, section_idx,
                                           variable, h);
  Aircraft aircraft_minus = PerturbAircraft(aircraft, surface_idx, section_idx,
                                            variable, -h);

  // Compute forces at perturbed configurations
  AerodynamicResults res_plus = ComputeResultsForAircraft(aircraft_plus, condition);
  AerodynamicResults res_minus = ComputeResultsForAircraft(aircraft_minus, condition);

  // Central difference
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

}  // namespace

AeroGradients GradientCalculator::ComputeSectionGradient(
    const FlightCondition& condition,
    Index surface_index,
    Index section_index,
    GeometryVariable variable) {
  if (!impl_->aircraft_set_) {
    return AeroGradients{};
  }

  return ComputeSectionGradientFD(impl_->aircraft_, condition,
                                  surface_index, section_index,
                                  variable, fd_step_);
}

GeometryGradients GradientCalculator::ComputeGeometryGradients(
    const FlightCondition& condition) {
  GeometryGradients grads;

  if (!impl_->aircraft_set_) {
    return grads;
  }

  const Aircraft& aircraft = impl_->aircraft_;

  // Iterate over all surfaces and sections
  for (Index surf_idx = 0;
       surf_idx < static_cast<Index>(aircraft.GetNumSurfaces());
       ++surf_idx) {
    const Surface& surface = aircraft.GetSurface(surf_idx);

    for (Index sect_idx = 0;
         sect_idx < static_cast<Index>(surface.GetNumSections());
         ++sect_idx) {
      SectionGeometryGradients section_grads;
      section_grads.surface_name = surface.GetName();
      section_grads.section_index = sect_idx;

      // Compute gradients for each geometry variable
      section_grads.chord_derivs = ComputeSectionGradient(
          condition, surf_idx, sect_idx, GeometryVariable::kChord);
      section_grads.twist_derivs = ComputeSectionGradient(
          condition, surf_idx, sect_idx, GeometryVariable::kTwist);
      section_grads.xle_derivs = ComputeSectionGradient(
          condition, surf_idx, sect_idx, GeometryVariable::kXle);
      section_grads.yle_derivs = ComputeSectionGradient(
          condition, surf_idx, sect_idx, GeometryVariable::kYle);
      section_grads.zle_derivs = ComputeSectionGradient(
          condition, surf_idx, sect_idx, GeometryVariable::kZle);

      grads.section_gradients.push_back(section_grads);
    }
  }

  return grads;
}

}  // namespace ad
}  // namespace vlm
