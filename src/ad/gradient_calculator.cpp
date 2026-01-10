// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include "vlm/ad/gradient_calculator.h"

#include <cmath>
#include <vector>

#include "vlm/ad/adjoint_solver.h"
#include "vlm/ad/enzyme_forces.h"
#include "vlm/ad/enzyme_biot_savart.h"
#include "vlm/ad/enzyme_mesh.h"
#include "vlm/ad/enzyme_aic.h"
#include "vlm/ad/mesh_sensitivity.h"
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

// Forward declaration for FD function used by adjoint implementation
AeroGradients ComputeSectionGradientFD(const Aircraft& aircraft,
                                       const FlightCondition& condition,
                                       Index surface_idx,
                                       Index section_idx,
                                       GeometryVariable variable,
                                       Scalar fd_step);

#ifdef VLM_USE_ENZYME

/// Adjoint-based geometry gradient computation
///
/// For a geometry parameter θ at section s, the total gradient is:
/// dJ/dθ = ∂J/∂θ|_Γ + λᵀ·(∂b/∂θ) - λᵀ·(∂A/∂θ)·Γ
///
/// where:
/// - λ is the adjoint vector (solves Aᵀλ = ∂J/∂Γ)
/// - ∂b/∂θ is the RHS sensitivity (how normal vectors change)
/// - ∂A/∂θ is the AIC matrix sensitivity (how induced velocities change)
///
/// This implementation uses analytical mesh gradients and Biot-Savart gradients.

/// Compute the gradient of the RHS with respect to a geometry parameter.
/// The RHS is b_i = -n_i · V_inf, so ∂b_i/∂θ = -(∂n_i/∂θ) · V_inf
///
/// For geometry changes, only the normal vector n changes (through twist).
/// For chord/position changes, the normal also changes through the bound vortex direction.
VecX ComputeRHSGeometrySensitivity(
    const VortexLattice& mesh,
    const FlightCondition& condition,
    Index surface_idx,
    [[maybe_unused]] Index section_idx,  // TODO: Use for localized sensitivity
    GeometryVariable variable) {

  const Index n = mesh.GetNumVortices();
  VecX db_dtheta = VecX::Zero(n);

  // Compute freestream direction
  double v_inf[3];
  compute_freestream_direction(condition.alpha, condition.beta, v_inf);

  // For twist changes, the normal vector rotates
  // For chord changes, bound vortex positions shift, affecting the normal
  // For LE position changes, everything shifts, affecting the normal

  // Get surface mesh info
  if (surface_idx >= mesh.GetNumSurfaces()) return db_dtheta;
  const SurfaceMeshInfo& surf_info = mesh.GetSurfaceInfo(surface_idx);

  // Identify affected vortices (those in strips adjacent to this section)
  // This is a simplified model - in practice we'd need the exact section-to-strip mapping
  // For now, assume uniform strip distribution and affect all strips on this surface
  // TODO: Use section_idx to only affect panels adjacent to this section

  // For twist perturbation, the normal change is approximately:
  // ∂n/∂twist ≈ (chord_dir_perp) where chord_dir_perp is perpendicular to current normal
  // The exact derivative depends on the twist axis and current orientation

  if (variable == GeometryVariable::kTwist) {
    // Twist affects normals directly
    // ∂n/∂twist = rotation of normal about spanwise (bound vortex) axis
    // For small angles: ∂n/∂θ ≈ (spanwise × n)

    for (Index i = surf_info.first_vortex_index;
         i < surf_info.first_vortex_index + surf_info.num_vortices; ++i) {
      const HorseshoeVortex& vor = mesh.GetVortex(i);

      // Check if this vortex is affected by this section
      // Use strip's surface_index and position to determine influence
      // This is a simplification - exact mapping would require tracking section intervals

      // For now, apply to all vortices in the surface with a weight based on position
      // A proper implementation would track which strips correspond to which sections

      // Compute bound vortex (spanwise) direction - this is the twist rotation axis
      Vec3 bound_dir = (vor.bound_vortex_right - vor.bound_vortex_left).normalized();

      Vec3 n_vec = vor.normal_at_control_point;

      // Twist rotates the normal about the spanwise axis
      // ∂n/∂twist = spanwise × n
      Vec3 dn_dtwist = bound_dir.cross(n_vec);

      // Convert from radians to degrees (twist is stored in degrees)
      dn_dtwist *= kDegreesToRadians;

      // ∂b/∂θ = -∂n/∂θ · V_inf
      Scalar db = -(dn_dtwist.x() * v_inf[0] +
                    dn_dtwist.y() * v_inf[1] +
                    dn_dtwist.z() * v_inf[2]);

      db_dtheta(i) = db;
    }
  }
  // For chord and position changes, the effect on normals is through
  // bound vortex direction changes - this is more complex and often small
  // For now, return zero for these (falling back to FD)

  return db_dtheta;
}

/// Compute λᵀ·(∂A/∂θ)·Γ for a geometry parameter
/// This is the contribution from how the AIC matrix changes with geometry.
///
/// The AIC sensitivity uses the adjoint approach:
/// - Row contributions (control point changes): enzyme_aic_row_gradient_rc with Γ as weights
/// - Column contributions (bound vortex changes): enzyme_aic_column_gradient_rv with λ as weights
///
/// Chain rule: λᵀ·(∂A/∂θ)·Γ = Σ_mesh_vars (∂(λᵀ·A·Γ)/∂mesh_var) · (∂mesh_var/∂θ)
/// Helper function to compute AIC sensitivity for one interval
Scalar ComputeAICGeometrySensitivityForInterval(
    const VortexLattice& mesh,
    const MeshArrays& arrays,
    const VecX& lambda,
    const VecX& gamma,
    Index first_strip,
    Index num_strips,
    Index num_chordwise,
    double core_chord_frac,
    double core_width_frac,
    double beta,
    GeometryVariable variable,
    bool is_left_section) {  // true if section is left of interval, false if right

  Scalar result = 0.0;
  const Index n = arrays.num_vortices;

  // Loop over strips in this section interval
  for (Index istrip = 0; istrip < num_strips; ++istrip) {
    Index strip_idx = first_strip + istrip;
    const Strip& strip = mesh.GetStrip(strip_idx);

    // Compute interpolation fractions
    double f1, f2;
    ComputeStripEdgeFractions(istrip, num_strips, &f1, &f2);
    double fc = ComputeStripInterpolationFraction(istrip, num_strips);

    // Section contribution factors:
    // - Left section uses (1-f) factors
    // - Right section uses f factors
    double w1 = is_left_section ? (1.0 - f1) : f1;
    double w2 = is_left_section ? (1.0 - f2) : f2;
    double wc = is_left_section ? (1.0 - fc) : fc;

    // Loop over vortex elements in this strip
    for (Index ichord = 0; ichord < num_chordwise; ++ichord) {
      Index vor_idx = strip.first_vortex_index + ichord;

      // Compute chordwise fraction
      double dx = 1.0 / static_cast<double>(num_chordwise);
      double x1 = ichord * dx;
      double xvr = x1 + 0.25 * dx;  // Vortex at 1/4 panel
      double xcp = x1 + 0.75 * dx;  // Control at 3/4 panel

      // =========================================================
      // Compute mesh position sensitivities ∂mesh/∂θ
      // =========================================================
      double d_rv1[3] = {0, 0, 0};
      double d_rv2[3] = {0, 0, 0};
      double d_rc[3] = {0, 0, 0};
      double d_normal[3] = {0, 0, 0};

      switch (variable) {
        case GeometryVariable::kChord:
          // Chord affects x-positions: RV1.x = RLE1.x + XVR*CHORD1
          d_rv1[0] = w1 * xvr;
          d_rv2[0] = w2 * xvr;
          d_rc[0] = wc * xcp;
          break;

        case GeometryVariable::kXle:
          // LE.x affects x-positions directly
          d_rv1[0] = w1;
          d_rv2[0] = w2;
          d_rc[0] = wc;
          break;

        case GeometryVariable::kYle:
          // LE.y affects y-positions
          d_rv1[1] = w1;
          d_rv2[1] = w2;
          d_rc[1] = wc;
          break;

        case GeometryVariable::kZle:
          // LE.z affects z-positions
          d_rv1[2] = w1;
          d_rv2[2] = w2;
          d_rc[2] = wc;
          break;

        case GeometryVariable::kTwist: {
          // Twist primarily affects the normal vector
          // The effect on positions is secondary (through the bound vortex)
          // For now, focus on the normal contribution
          const HorseshoeVortex& vor = mesh.GetVortex(vor_idx);

          // Compute bound direction
          double dx_b = vor.bound_vortex_right.x() - vor.bound_vortex_left.x();
          double dy_b = vor.bound_vortex_right.y() - vor.bound_vortex_left.y();
          double dz_b = vor.bound_vortex_right.z() - vor.bound_vortex_left.z();
          double bound_len = std::sqrt(dx_b*dx_b + dy_b*dy_b + dz_b*dz_b);
          double bound_dir[3] = {0, 1, 0};
          if (bound_len > 1e-12) {
            bound_dir[0] = dx_b / bound_len;
            bound_dir[1] = dy_b / bound_len;
            bound_dir[2] = dz_b / bound_len;
          }

          // Compute d(normal)/d(twist) analytically
          double twist = strip.twist_angle;
          double camber = vor.camber_slope_at_control_point;
          double angle = twist - std::atan(camber);
          double ca = std::cos(angle);
          double sa = std::sin(angle);

          // Chord direction: EC = (cos, -sin*ES_y, -sin*ES_z)
          // d(EC)/d(twist) = (-sin, -cos*ES_y, -cos*ES_z)
          double dec[3];
          dec[0] = -sa;
          dec[1] = -ca * strip.trefftz_normal_y;
          dec[2] = -ca * strip.trefftz_normal_z;

          // Normal = normalize(EC × EB)
          double ec[3] = {ca, -sa * strip.trefftz_normal_y, -sa * strip.trefftz_normal_z};
          double ecxb[3];
          ecxb[0] = ec[1] * bound_dir[2] - ec[2] * bound_dir[1];
          ecxb[1] = ec[2] * bound_dir[0] - ec[0] * bound_dir[2];
          ecxb[2] = ec[0] * bound_dir[1] - ec[1] * bound_dir[0];
          double emag = std::sqrt(ecxb[0]*ecxb[0] + ecxb[1]*ecxb[1] + ecxb[2]*ecxb[2]);

          if (emag > 1e-12) {
            // d(ECXB)/d(twist) = d(EC)/d(twist) × EB
            double decxb[3];
            decxb[0] = dec[1] * bound_dir[2] - dec[2] * bound_dir[1];
            decxb[1] = dec[2] * bound_dir[0] - dec[0] * bound_dir[2];
            decxb[2] = dec[0] * bound_dir[1] - dec[1] * bound_dir[0];

            double normal[3] = {ecxb[0]/emag, ecxb[1]/emag, ecxb[2]/emag};
            double demag = normal[0]*decxb[0] + normal[1]*decxb[1] + normal[2]*decxb[2];

            d_normal[0] = wc * (decxb[0] - normal[0]*demag) / emag;
            d_normal[1] = wc * (decxb[1] - normal[1]*demag) / emag;
            d_normal[2] = wc * (decxb[2] - normal[2]*demag) / emag;
          }
          break;
        }
      }

      // =========================================================
      // Compute AIC gradient contributions
      // =========================================================

      // Row contribution: How does changing the control point/normal at row vor_idx
      // affect the AIC row, weighted by Γ?
      // The result is d(AIC_row · Γ)/d(rc) and d(AIC_row · Γ)/d(normal)
      if (d_rc[0] != 0.0 || d_rc[1] != 0.0 || d_rc[2] != 0.0 ||
          d_normal[0] != 0.0 || d_normal[1] != 0.0 || d_normal[2] != 0.0) {

        double d_rc_out[3] = {0, 0, 0};
        double d_normal_out[3] = {0, 0, 0};

        // Use enzyme to get d(AIC_row · Γ)/d(rc, normal)
        // Note: enzyme_aic_row_gradient_rc computes gradient of weighted sum
#ifdef VLM_USE_ENZYME
        enzyme_aic_row_gradient_rc(
            arrays.rv1.data(),
            arrays.rv2.data(),
            &arrays.rc[3 * vor_idx],
            &arrays.normal[3 * vor_idx],
            arrays.chord.data(),
            arrays.surface_idx.data(),
            static_cast<int>(n),
            arrays.surface_idx[vor_idx],  // Surface index of this control point
            core_chord_frac,
            core_width_frac,
            beta,
            gamma.data(),  // Use Γ as weights
            d_rc_out,
            d_normal_out);
#else
        (void)core_chord_frac;
        (void)core_width_frac;
        (void)beta;
#endif

        // Chain rule: contribution = d(AIC·Γ)/d(rc) · d(rc)/d(θ) + d(AIC·Γ)/d(normal) · d(normal)/d(θ)
        // This gives us the row's contribution to Σ_j A_ij·Γ_j, weighted by λ_i
        Scalar row_contrib = d_rc_out[0] * d_rc[0] + d_rc_out[1] * d_rc[1] + d_rc_out[2] * d_rc[2]
                           + d_normal_out[0] * d_normal[0] + d_normal_out[1] * d_normal[1]
                           + d_normal_out[2] * d_normal[2];

        result += lambda(vor_idx) * row_contrib;
      }

      // Column contribution: How does changing the bound vortex at column vor_idx
      // affect all AIC rows, weighted by λ?
      if (d_rv1[0] != 0.0 || d_rv1[1] != 0.0 || d_rv1[2] != 0.0 ||
          d_rv2[0] != 0.0 || d_rv2[1] != 0.0 || d_rv2[2] != 0.0) {

        double d_rv1_out[3] = {0, 0, 0};
        double d_rv2_out[3] = {0, 0, 0};

#ifdef VLM_USE_ENZYME
        enzyme_aic_column_gradient_rv(
            arrays.rv1.data(),
            arrays.rv2.data(),
            arrays.rc.data(),
            arrays.normal.data(),
            arrays.chord.data(),
            arrays.surface_idx.data(),
            static_cast<int>(n),
            static_cast<int>(vor_idx),
            core_chord_frac,
            core_width_frac,
            beta,
            lambda.data(),  // Use λ as weights
            d_rv1_out,
            d_rv2_out);
#endif

        // Chain rule: contribution = d(λᵀ·A_col)/d(rv1) · d(rv1)/d(θ) + ...
        // This gives us the column's contribution, already weighted by λ
        Scalar col_contrib = d_rv1_out[0] * d_rv1[0] + d_rv1_out[1] * d_rv1[1] + d_rv1_out[2] * d_rv1[2]
                           + d_rv2_out[0] * d_rv2[0] + d_rv2_out[1] * d_rv2[1] + d_rv2_out[2] * d_rv2[2];

        result += gamma(vor_idx) * col_contrib;
      }
    }
  }

  return result;
}

/// Compute AIC sensitivity for a section, considering both adjacent intervals.
/// A section contributes to:
/// - Interval [section_idx, section_idx+1] as the "left" section (if exists)
/// - Interval [section_idx-1, section_idx] as the "right" section (if exists)
/// This matches OptVL's approach where each section gets contributions from both sides.
Scalar ComputeAICGeometrySensitivity(
    const VortexLattice& mesh,
    const Aircraft& aircraft,
    const VecX& lambda,
    const VecX& gamma,
    const FlightCondition& condition,
    Index surface_idx,
    Index section_idx,
    GeometryVariable variable) {

  const Index n = mesh.GetNumVortices();
  if (n == 0) return 0.0;

  // Build section-strip mapping
  SectionStripMapping mapping = BuildSectionStripMapping(aircraft, mesh, surface_idx);
  if (mapping.num_sections == 0) return 0.0;

  const auto& surface_info = mesh.GetSurfaceInfo(surface_idx);
  Index num_chordwise = surface_info.num_chordwise;

  // Extract mesh arrays for enzyme functions
  MeshArrays arrays = ExtractMeshArrays(mesh);

  // Compressibility factor
  double beta = std::sqrt(std::max(0.01, 1.0 - condition.mach * condition.mach));

  // Core radius fractions
  double core_chord_frac = mesh.GetVortexCoreRadiusChordFraction();
  double core_width_frac = mesh.GetVortexCoreRadiusWidthFraction();

  Scalar result = 0.0;

  // Contribution as "left" section of interval [section_idx, section_idx+1]
  if (section_idx < mapping.num_sections - 1) {
    auto [first_strip, num_strips] = mapping.GetStripRange(section_idx);
    if (num_strips > 0) {
      result += ComputeAICGeometrySensitivityForInterval(
          mesh, arrays, lambda, gamma,
          first_strip, num_strips, num_chordwise,
          core_chord_frac, core_width_frac, beta,
          variable, true);  // is_left_section = true
    }
  }

  // Contribution as "right" section of interval [section_idx-1, section_idx]
  if (section_idx > 0) {
    auto [first_strip, num_strips] = mapping.GetStripRange(section_idx - 1);
    if (num_strips > 0) {
      result += ComputeAICGeometrySensitivityForInterval(
          mesh, arrays, lambda, gamma,
          first_strip, num_strips, num_chordwise,
          core_chord_frac, core_width_frac, beta,
          variable, false);  // is_left_section = false
    }
  }

  return result;
}

/// Compute geometry gradient using the adjoint method (partial implementation)
/// Falls back to FD for terms that are not yet analytically implemented.
AeroGradients ComputeSectionGradientAdjoint(
    const Aircraft& aircraft,
    const FlightCondition& condition,
    Index surface_idx,
    Index section_idx,
    GeometryVariable variable,
    Scalar fd_step) {

  // Create solver and compute solution
  VortexLatticeSolver solver;
  solver.SetAircraft(aircraft);
  SolutionState solution = solver.Solve(condition);

  const VortexLattice& mesh = solver.GetMesh();
  const Index n = solution.circulation.size();

  if (n == 0) {
    return AeroGradients{};
  }

  // Step 1: Compute force gradients ∂J/∂Γ
  ForceGradients force_grads = ComputeForceGradients(mesh, solution, condition, aircraft);

  // Step 2: Solve adjoint system Aᵀλ = ∂J/∂Γ
  AdjointSolver adjoint(solver);
  if (!adjoint.IsValid()) {
    // Fall back to full FD
    return ComputeSectionGradientFD(aircraft, condition, surface_idx,
                                    section_idx, variable, fd_step);
  }

  // Solve all adjoints at once
  MatX rhs(n, 6);  // CL, CD, CY, Cl, Cm, Cn
  rhs.col(0) = force_grads.dCL_dGamma;
  rhs.col(1) = force_grads.dCD_dGamma;
  rhs.col(2) = force_grads.dCY_dGamma;
  rhs.col(3) = force_grads.dCl_dGamma;
  rhs.col(4) = force_grads.dCm_dGamma;
  rhs.col(5) = force_grads.dCn_dGamma;

  MatX lambda = adjoint.SolveAdjointMultiple(rhs);

  // Step 3: Compute RHS sensitivity ∂b/∂θ
  VecX db_dtheta = ComputeRHSGeometrySensitivity(
      mesh, condition, surface_idx, section_idx, variable);

  // Step 4: Compute implicit term λᵀ·(∂b/∂θ)
  AeroGradients implicit_rhs;
  implicit_rhs.dCL = lambda.col(0).dot(db_dtheta);
  implicit_rhs.dCD = lambda.col(1).dot(db_dtheta);
  implicit_rhs.dCY = lambda.col(2).dot(db_dtheta);
  implicit_rhs.dCl = lambda.col(3).dot(db_dtheta);
  implicit_rhs.dCm = lambda.col(4).dot(db_dtheta);
  implicit_rhs.dCn = lambda.col(5).dot(db_dtheta);
  implicit_rhs.dCDi = 0.0;  // CDi computed separately via Trefftz plane

  // Step 5: Compute AIC sensitivity term -λᵀ·(∂A/∂θ)·Γ
  // This uses the enzyme_aic gradients combined with mesh sensitivity
  AeroGradients aic_sens;
  aic_sens.dCL = -ComputeAICGeometrySensitivity(mesh, aircraft, lambda.col(0),
      solution.circulation, condition, surface_idx, section_idx, variable);
  aic_sens.dCD = -ComputeAICGeometrySensitivity(mesh, aircraft, lambda.col(1),
      solution.circulation, condition, surface_idx, section_idx, variable);
  aic_sens.dCY = -ComputeAICGeometrySensitivity(mesh, aircraft, lambda.col(2),
      solution.circulation, condition, surface_idx, section_idx, variable);
  aic_sens.dCl = -ComputeAICGeometrySensitivity(mesh, aircraft, lambda.col(3),
      solution.circulation, condition, surface_idx, section_idx, variable);
  aic_sens.dCm = -ComputeAICGeometrySensitivity(mesh, aircraft, lambda.col(4),
      solution.circulation, condition, surface_idx, section_idx, variable);
  aic_sens.dCn = -ComputeAICGeometrySensitivity(mesh, aircraft, lambda.col(5),
      solution.circulation, condition, surface_idx, section_idx, variable);
  aic_sens.dCDi = 0.0;  // CDi handled via Trefftz plane

  // Step 6: Compute explicit term ∂J/∂θ|_Γ (forces change at fixed circulation)
  // For geometry changes, the explicit force term changes through:
  // - Bound vortex direction changes
  // - Reference point changes for moments
  // This is a secondary effect compared to the AIC/RHS terms
  // For now, assume it's small and set to zero (could add later)
  AeroGradients explicit_term;
  explicit_term.dCL = 0.0;
  explicit_term.dCD = 0.0;
  explicit_term.dCY = 0.0;
  explicit_term.dCl = 0.0;
  explicit_term.dCm = 0.0;
  explicit_term.dCn = 0.0;
  explicit_term.dCDi = 0.0;

  // Total gradient = explicit + implicit_rhs + aic_sens
  AeroGradients result;
  result.dCL = explicit_term.dCL + implicit_rhs.dCL + aic_sens.dCL;
  result.dCD = explicit_term.dCD + implicit_rhs.dCD + aic_sens.dCD;
  result.dCY = explicit_term.dCY + implicit_rhs.dCY + aic_sens.dCY;
  result.dCl = explicit_term.dCl + implicit_rhs.dCl + aic_sens.dCl;
  result.dCm = explicit_term.dCm + implicit_rhs.dCm + aic_sens.dCm;
  result.dCn = explicit_term.dCn + implicit_rhs.dCn + aic_sens.dCn;
  result.dCDi = 0.0;  // TODO: Add Trefftz plane gradient

  return result;
}

#endif  // VLM_USE_ENZYME

// Forward declaration for FD function
AeroGradients ComputeSectionGradientFD(const Aircraft& aircraft,
                                       const FlightCondition& condition,
                                       Index surface_idx,
                                       Index section_idx,
                                       GeometryVariable variable,
                                       Scalar fd_step);

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

#ifdef VLM_USE_ENZYME
  // Use adjoint-based approach when Enzyme is available
  // Note: This currently falls back to FD for the AIC sensitivity term,
  // but provides the infrastructure for full analytical gradients
  return ComputeSectionGradientAdjoint(impl_->aircraft_, condition,
                                       surface_index, section_index,
                                       variable, fd_step_);
#else
  return ComputeSectionGradientFD(impl_->aircraft_, condition,
                                  surface_index, section_index,
                                  variable, fd_step_);
#endif
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
