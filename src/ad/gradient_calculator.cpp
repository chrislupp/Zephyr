// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include "vlm/ad/gradient_calculator.h"

#include <cmath>
#include <vector>

#include "vlm/ad/adjoint_solver.h"
#include "vlm/ad/analytical_forces.h"
#include "vlm/ad/analytical_aic.h"
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
  bool use_adjoint_{true};

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

  /// Extract mesh data into raw arrays for analytical kernels
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

  /// Compute explicit term ∂J/∂param at fixed circulation using analytical derivatives
  AeroGradients ComputeExplicitDerivatives(int param_index) {
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
      analytical_force_alpha_gradient(
          cached_condition_.alpha, cached_condition_.beta,
          circulation.data(), induced_velocity.data(),
          bound_left.data(), bound_right.data(),
          control_points.data(), moment_arm.data(),
          static_cast<int>(n), reference_area, d_result);
    } else if (param_index == 1) {
      // Beta derivative
      analytical_force_beta_gradient(
          cached_condition_.alpha, cached_condition_.beta,
          circulation.data(), induced_velocity.data(),
          bound_left.data(), bound_right.data(),
          control_points.data(), moment_arm.data(),
          static_cast<int>(n), reference_area, d_result);
    } else {
      // Rotation rates (p, q, r): The force kernel only uses V_inf + V_induced
      // and does not include V_rotation. The rotation effect on forces is
      // entirely captured through the implicit term (how Γ changes when the
      // RHS changes with rotation rate). The explicit force term is zero.
      for (int out_idx = 0; out_idx < 6; ++out_idx) {
        d_result[out_idx] = 0.0;
      }
    }

    grads.dCL = d_result[0];
    grads.dCD = d_result[1];
    grads.dCY = d_result[2];
    grads.dCl = d_result[3];
    grads.dCm = d_result[4];
    grads.dCn = d_result[5];

    // CDi_ff (Trefftz plane drag) doesn't have an explicit alpha/beta/rotation term
    // because it's computed from circulation and geometry only
    grads.dCDi = 0.0;

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
      return AeroGradients{};
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
    total.dCDi = explicit_grads.dCDi;

    return total;
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
  impl_->solution_cached_ = false;
}

void GradientCalculator::SetReferencePoint(const Vec3& ref) {
  impl_->solver_.SetReferencePoint(ref);
}

void GradientCalculator::SetFiniteDiffStep(Scalar) {
  // No-op: analytical derivatives don't use finite differences
}

// ============================================================================
// Derivative Computation Methods
// ============================================================================

AeroGradients GradientCalculator::ComputeAlphaDerivatives(
    const FlightCondition& condition) {
  return impl_->AdjointDerivatives(condition, 0);
}

AeroGradients GradientCalculator::ComputeBetaDerivatives(
    const FlightCondition& condition) {
  return impl_->AdjointDerivatives(condition, 1);
}

AeroGradients GradientCalculator::ComputeRollRateDerivatives(
    const FlightCondition& condition) {
  return impl_->AdjointDerivatives(condition, 2);
}

AeroGradients GradientCalculator::ComputePitchRateDerivatives(
    const FlightCondition& condition) {
  return impl_->AdjointDerivatives(condition, 3);
}

AeroGradients GradientCalculator::ComputeYawRateDerivatives(
    const FlightCondition& condition) {
  return impl_->AdjointDerivatives(condition, 4);
}

AeroGradients GradientCalculator::ComputeMachDerivatives(
    const FlightCondition& condition) {
  // Mach changes the Prandtl-Glauert factor, which affects the AIC matrix.
  // This requires re-solving the system at perturbed Mach numbers.
  // TODO: Implement analytical ∂A/∂M via adjoint with PG chain rule
  const Scalar h = 0.001;

  FlightCondition cond_plus = condition;
  FlightCondition cond_minus = condition;
  cond_plus.mach = condition.mach + h;
  cond_minus.mach = condition.mach - h;

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

/// Compute the gradient of the RHS with respect to a geometry parameter.
/// The RHS is b_i = -n_i · V_inf, so ∂b_i/∂θ = -(∂n_i/∂θ) · V_inf
VecX ComputeRHSGeometrySensitivity(
    const VortexLattice& mesh,
    const Aircraft& aircraft,
    const FlightCondition& condition,
    Index surface_idx,
    Index section_idx,
    GeometryVariable variable) {

  const Index n = mesh.GetNumVortices();
  VecX db_dtheta = VecX::Zero(n);

  double v_inf[3];
  compute_freestream_direction(condition.alpha, condition.beta, v_inf);

  if (surface_idx >= mesh.GetNumSurfaces()) return db_dtheta;
  const SurfaceMeshInfo& surf_info = mesh.GetSurfaceInfo(surface_idx);

  if (variable == GeometryVariable::kTwist) {
    SectionStripMapping mapping = BuildSectionStripMapping(
        aircraft, mesh, surface_idx);
    if (mapping.num_sections == 0) return db_dtheta;

    Index num_chordwise = surf_info.num_chordwise;

    auto compute_twist_contribution = [&](Index first_strip, Index num_strips, bool is_left_section) {
      for (Index istrip = 0; istrip < num_strips; ++istrip) {
        Index strip_idx = first_strip + istrip;
        const Strip& strip = mesh.GetStrip(strip_idx);

        double f1, f2;
        ComputeStripEdgeFractions(istrip, num_strips, &f1, &f2);

        double fc = (f1 + f2) / 2.0;
        double weight = is_left_section ? (1.0 - fc) : fc;

        for (Index ichord = 0; ichord < num_chordwise; ++ichord) {
          Index vor_idx = strip.first_vortex_index + ichord;
          const HorseshoeVortex& vor = mesh.GetVortex(vor_idx);

          Vec3 bound_dir = (vor.bound_vortex_right - vor.bound_vortex_left).normalized();
          Vec3 n_vec = vor.normal_at_control_point;

          Vec3 dn_dtwist = bound_dir.cross(n_vec);
          dn_dtwist *= kDegreesToRadians;
          dn_dtwist *= weight;

          Scalar db = -(dn_dtwist.x() * v_inf[0] +
                        dn_dtwist.y() * v_inf[1] +
                        dn_dtwist.z() * v_inf[2]);

          db_dtheta(vor_idx) += db;
        }
      }
    };

    if (section_idx < mapping.num_sections - 1) {
      auto [first_strip, num_strips] = mapping.GetStripRange(section_idx);
      if (num_strips > 0) {
        compute_twist_contribution(first_strip, num_strips, true);
      }
    }

    if (section_idx > 0) {
      auto [first_strip, num_strips] = mapping.GetStripRange(section_idx - 1);
      if (num_strips > 0) {
        compute_twist_contribution(first_strip, num_strips, false);
      }
    }
  } else if (variable == GeometryVariable::kChord ||
             variable == GeometryVariable::kXle ||
             variable == GeometryVariable::kYle ||
             variable == GeometryVariable::kZle) {
    SectionStripMapping mapping = BuildSectionStripMapping(
        aircraft, mesh, surface_idx);
    if (mapping.num_sections == 0) return db_dtheta;

    Index num_chordwise = surf_info.num_chordwise;

    auto compute_interval_contribution = [&](Index first_strip, Index num_strips, bool is_left_section) {
      for (Index istrip = 0; istrip < num_strips; ++istrip) {
        Index strip_idx = first_strip + istrip;
        const Strip& strip = mesh.GetStrip(strip_idx);

        double f1, f2;
        ComputeStripEdgeFractions(istrip, num_strips, &f1, &f2);

        double w1 = is_left_section ? (1.0 - f1) : f1;
        double w2 = is_left_section ? (1.0 - f2) : f2;

        for (Index ichord = 0; ichord < num_chordwise; ++ichord) {
          Index vor_idx = strip.first_vortex_index + ichord;
          const HorseshoeVortex& vor = mesh.GetVortex(vor_idx);

          double dx = 1.0 / static_cast<double>(num_chordwise);
          double x1 = ichord * dx;
          double xvr = x1 + 0.25 * dx;

          double d_rv1[3] = {0, 0, 0};
          double d_rv2[3] = {0, 0, 0};

          switch (variable) {
            case GeometryVariable::kChord:
              d_rv1[0] = w1 * xvr;
              d_rv2[0] = w2 * xvr;
              break;
            case GeometryVariable::kXle:
              d_rv1[0] = w1;
              d_rv2[0] = w2;
              break;
            case GeometryVariable::kYle:
              d_rv1[1] = w1;
              d_rv2[1] = w2;
              break;
            case GeometryVariable::kZle:
              d_rv1[2] = w1;
              d_rv2[2] = w2;
              break;
            default:
              break;
          }

          Vec3 bound_vec = vor.bound_vortex_right - vor.bound_vortex_left;
          double bound_len = bound_vec.norm();
          if (bound_len < 1e-12) continue;
          Vec3 bound_dir = bound_vec / bound_len;

          double d_rv[3] = {d_rv2[0] - d_rv1[0], d_rv2[1] - d_rv1[1], d_rv2[2] - d_rv1[2]};

          double bd_dot_drv = bound_dir.x() * d_rv[0] + bound_dir.y() * d_rv[1] + bound_dir.z() * d_rv[2];
          Vec3 d_bound_dir;
          d_bound_dir.x() = d_rv[0] / bound_len - bound_dir.x() * bd_dot_drv / bound_len;
          d_bound_dir.y() = d_rv[1] / bound_len - bound_dir.y() * bd_dot_drv / bound_len;
          d_bound_dir.z() = d_rv[2] / bound_len - bound_dir.z() * bd_dot_drv / bound_len;

          Vec3 normal = vor.normal_at_control_point;
          Vec3 chord_dir = bound_dir.cross(normal);
          double chord_len = chord_dir.norm();
          if (chord_len < 1e-12) {
            chord_dir = Vec3(1, 0, 0);
          } else {
            chord_dir /= chord_len;
          }

          Vec3 d_normal_unnorm = chord_dir.cross(d_bound_dir);

          Vec3 n_unnorm = chord_dir.cross(bound_dir);
          double n_unnorm_len = n_unnorm.norm();
          if (n_unnorm_len < 1e-12) continue;

          double n_dot_dnu = normal.x() * d_normal_unnorm.x() +
                             normal.y() * d_normal_unnorm.y() +
                             normal.z() * d_normal_unnorm.z();
          Vec3 d_normal;
          d_normal.x() = d_normal_unnorm.x() / n_unnorm_len - normal.x() * n_dot_dnu / n_unnorm_len;
          d_normal.y() = d_normal_unnorm.y() / n_unnorm_len - normal.y() * n_dot_dnu / n_unnorm_len;
          d_normal.z() = d_normal_unnorm.z() / n_unnorm_len - normal.z() * n_dot_dnu / n_unnorm_len;

          Scalar db = -(d_normal.x() * v_inf[0] +
                        d_normal.y() * v_inf[1] +
                        d_normal.z() * v_inf[2]);

          db_dtheta(vor_idx) += db;
        }
      }
    };

    if (section_idx < mapping.num_sections - 1) {
      auto [first_strip, num_strips] = mapping.GetStripRange(section_idx);
      if (num_strips > 0) {
        compute_interval_contribution(first_strip, num_strips, true);
      }
    }

    if (section_idx > 0) {
      auto [first_strip, num_strips] = mapping.GetStripRange(section_idx - 1);
      if (num_strips > 0) {
        compute_interval_contribution(first_strip, num_strips, false);
      }
    }
  }

  return db_dtheta;
}

/// Compute λᵀ·(∂A/∂θ)·Γ for a geometry parameter using analytical AIC gradients.
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
    bool is_left_section) {

  Scalar result = 0.0;
  const Index n = arrays.num_vortices;

  for (Index istrip = 0; istrip < num_strips; ++istrip) {
    Index strip_idx = first_strip + istrip;
    const Strip& strip = mesh.GetStrip(strip_idx);

    double f1, f2;
    ComputeStripEdgeFractions(istrip, num_strips, &f1, &f2);
    double fc = ComputeStripInterpolationFraction(istrip, num_strips);

    double w1 = is_left_section ? (1.0 - f1) : f1;
    double w2 = is_left_section ? (1.0 - f2) : f2;
    double wc = is_left_section ? (1.0 - fc) : fc;

    for (Index ichord = 0; ichord < num_chordwise; ++ichord) {
      Index vor_idx = strip.first_vortex_index + ichord;

      double dx = 1.0 / static_cast<double>(num_chordwise);
      double x1 = ichord * dx;
      double xvr = x1 + 0.25 * dx;
      double xcp = x1 + 0.75 * dx;

      double d_rv1[3] = {0, 0, 0};
      double d_rv2[3] = {0, 0, 0};
      double d_rc[3] = {0, 0, 0};
      double d_normal[3] = {0, 0, 0};

      switch (variable) {
        case GeometryVariable::kChord:
          d_rv1[0] = w1 * xvr;
          d_rv2[0] = w2 * xvr;
          d_rc[0] = wc * xcp;
          break;

        case GeometryVariable::kXle:
          d_rv1[0] = w1;
          d_rv2[0] = w2;
          d_rc[0] = wc;
          break;

        case GeometryVariable::kYle:
          d_rv1[1] = w1;
          d_rv2[1] = w2;
          d_rc[1] = wc;
          break;

        case GeometryVariable::kZle:
          d_rv1[2] = w1;
          d_rv2[2] = w2;
          d_rc[2] = wc;
          break;

        case GeometryVariable::kTwist: {
          const HorseshoeVortex& vor = mesh.GetVortex(vor_idx);

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

          double twist = strip.twist_angle;
          double camber = vor.camber_slope_at_control_point;
          double angle = twist - std::atan(camber);
          double ca = std::cos(angle);
          double sa = std::sin(angle);

          double dec[3];
          dec[0] = -sa;
          dec[1] = -ca * strip.trefftz_normal_y;
          dec[2] = -ca * strip.trefftz_normal_z;

          double ec[3] = {ca, -sa * strip.trefftz_normal_y, -sa * strip.trefftz_normal_z};
          double ecxb[3];
          ecxb[0] = ec[1] * bound_dir[2] - ec[2] * bound_dir[1];
          ecxb[1] = ec[2] * bound_dir[0] - ec[0] * bound_dir[2];
          ecxb[2] = ec[0] * bound_dir[1] - ec[1] * bound_dir[0];
          double emag = std::sqrt(ecxb[0]*ecxb[0] + ecxb[1]*ecxb[1] + ecxb[2]*ecxb[2]);

          if (emag > 1e-12) {
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

      // Row contribution: d(AIC_row · Γ)/d(rc, normal)
      if (d_rc[0] != 0.0 || d_rc[1] != 0.0 || d_rc[2] != 0.0 ||
          d_normal[0] != 0.0 || d_normal[1] != 0.0 || d_normal[2] != 0.0) {

        double d_rc_out[3] = {0, 0, 0};
        double d_normal_out[3] = {0, 0, 0};

        analytical_aic_row_gradient_rc(
            arrays.rv1.data(),
            arrays.rv2.data(),
            &arrays.rc[3 * vor_idx],
            &arrays.normal[3 * vor_idx],
            arrays.chord.data(),
            arrays.surface_idx.data(),
            static_cast<int>(n),
            arrays.surface_idx[vor_idx],
            core_chord_frac,
            core_width_frac,
            beta,
            gamma.data(),
            d_rc_out,
            d_normal_out);

        Scalar row_contrib = d_rc_out[0] * d_rc[0] + d_rc_out[1] * d_rc[1] + d_rc_out[2] * d_rc[2]
                           + d_normal_out[0] * d_normal[0] + d_normal_out[1] * d_normal[1]
                           + d_normal_out[2] * d_normal[2];

        result += lambda(vor_idx) * row_contrib;
      }

      // Column contribution: d(λᵀ·A_col)/d(rv1, rv2)
      if (d_rv1[0] != 0.0 || d_rv1[1] != 0.0 || d_rv1[2] != 0.0 ||
          d_rv2[0] != 0.0 || d_rv2[1] != 0.0 || d_rv2[2] != 0.0) {

        double d_rv1_out[3] = {0, 0, 0};
        double d_rv2_out[3] = {0, 0, 0};

        analytical_aic_column_gradient_rv(
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
            lambda.data(),
            d_rv1_out,
            d_rv2_out);

        Scalar col_contrib = d_rv1_out[0] * d_rv1[0] + d_rv1_out[1] * d_rv1[1] + d_rv1_out[2] * d_rv1[2]
                           + d_rv2_out[0] * d_rv2[0] + d_rv2_out[1] * d_rv2[1] + d_rv2_out[2] * d_rv2[2];

        result += gamma(vor_idx) * col_contrib;
      }
    }
  }

  return result;
}

/// Compute AIC sensitivity for a section, considering both adjacent intervals.
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

  SectionStripMapping mapping = BuildSectionStripMapping(aircraft, mesh, surface_idx);
  if (mapping.num_sections == 0) return 0.0;

  const auto& surface_info = mesh.GetSurfaceInfo(surface_idx);
  Index num_chordwise = surface_info.num_chordwise;

  MeshArrays arrays = ExtractMeshArrays(mesh);

  double beta = std::sqrt(std::max(0.01, 1.0 - condition.mach * condition.mach));

  double core_chord_frac = mesh.GetVortexCoreRadiusChordFraction();
  double core_width_frac = mesh.GetVortexCoreRadiusWidthFraction();

  Scalar result = 0.0;

  if (section_idx < mapping.num_sections - 1) {
    auto [first_strip, num_strips] = mapping.GetStripRange(section_idx);
    if (num_strips > 0) {
      result += ComputeAICGeometrySensitivityForInterval(
          mesh, arrays, lambda, gamma,
          first_strip, num_strips, num_chordwise,
          core_chord_frac, core_width_frac, beta,
          variable, true);
    }
  }

  if (section_idx > 0) {
    auto [first_strip, num_strips] = mapping.GetStripRange(section_idx - 1);
    if (num_strips > 0) {
      result += ComputeAICGeometrySensitivityForInterval(
          mesh, arrays, lambda, gamma,
          first_strip, num_strips, num_chordwise,
          core_chord_frac, core_width_frac, beta,
          variable, false);
    }
  }

  return result;
}

/// Compute velocity kernel geometry sensitivity for one interval.
Scalar ComputeVelocityKernelGeometrySensitivityForInterval(
    const VortexLattice& mesh,
    const MeshArrays& arrays,
    const std::vector<double>& dJ_dV,
    const VecX& gamma,
    Index first_strip,
    Index num_strips,
    Index num_chordwise,
    double core_chord_frac,
    double core_width_frac,
    double beta,
    GeometryVariable variable,
    bool is_left_section) {

  Scalar result = 0.0;
  const Index n = arrays.num_vortices;

  for (Index istrip = 0; istrip < num_strips; ++istrip) {
    Index strip_idx = first_strip + istrip;
    const Strip& strip = mesh.GetStrip(strip_idx);

    double f1, f2;
    ComputeStripEdgeFractions(istrip, num_strips, &f1, &f2);
    double fc = ComputeStripInterpolationFraction(istrip, num_strips);

    double w1 = is_left_section ? (1.0 - f1) : f1;
    double w2 = is_left_section ? (1.0 - f2) : f2;
    double wc = is_left_section ? (1.0 - fc) : fc;

    for (Index ichord = 0; ichord < num_chordwise; ++ichord) {
      Index vor_idx = strip.first_vortex_index + ichord;

      double dx = 1.0 / static_cast<double>(num_chordwise);
      double x1 = ichord * dx;
      double xvr = x1 + 0.25 * dx;
      double xcp = x1 + 0.75 * dx;

      double d_rv1[3] = {0, 0, 0};
      double d_rv2[3] = {0, 0, 0};
      double d_rc[3] = {0, 0, 0};

      switch (variable) {
        case GeometryVariable::kChord:
          d_rv1[0] = w1 * xvr;
          d_rv2[0] = w2 * xvr;
          d_rc[0] = wc * xcp;
          break;
        case GeometryVariable::kXle:
          d_rv1[0] = w1;
          d_rv2[0] = w2;
          d_rc[0] = wc;
          break;
        case GeometryVariable::kYle:
          d_rv1[1] = w1;
          d_rv2[1] = w2;
          d_rc[1] = wc;
          break;
        case GeometryVariable::kZle:
          d_rv1[2] = w1;
          d_rv2[2] = w2;
          d_rc[2] = wc;
          break;
        case GeometryVariable::kTwist:
          break;
      }

      // Row contribution: when control point moves
      if (d_rc[0] != 0.0 || d_rc[1] != 0.0 || d_rc[2] != 0.0) {
        double d_rc_out[3] = {0, 0, 0};
        const double* dJ_dV_i = &dJ_dV[3 * vor_idx];

        analytical_velocity_row_gradient_rc(
            arrays.rv1.data(),
            arrays.rv2.data(),
            &arrays.rc[3 * vor_idx],
            arrays.chord.data(),
            arrays.surface_idx.data(),
            static_cast<int>(n),
            arrays.surface_idx[vor_idx],
            core_chord_frac,
            core_width_frac,
            beta,
            gamma.data(),
            dJ_dV_i,
            d_rc_out);

        result += d_rc_out[0] * d_rc[0] + d_rc_out[1] * d_rc[1] + d_rc_out[2] * d_rc[2];
      }

      // Column contribution: when bound vortex moves
      if (d_rv1[0] != 0.0 || d_rv1[1] != 0.0 || d_rv1[2] != 0.0 ||
          d_rv2[0] != 0.0 || d_rv2[1] != 0.0 || d_rv2[2] != 0.0) {
        double d_rv1_out[3] = {0, 0, 0};
        double d_rv2_out[3] = {0, 0, 0};

        analytical_velocity_column_gradient_rv(
            arrays.rv1.data(),
            arrays.rv2.data(),
            arrays.rc.data(),
            arrays.chord.data(),
            arrays.surface_idx.data(),
            static_cast<int>(n),
            static_cast<int>(vor_idx),
            core_chord_frac,
            core_width_frac,
            beta,
            gamma(vor_idx),
            dJ_dV.data(),
            d_rv1_out,
            d_rv2_out);

        result += d_rv1_out[0] * d_rv1[0] + d_rv1_out[1] * d_rv1[1] + d_rv1_out[2] * d_rv1[2];
        result += d_rv2_out[0] * d_rv2[0] + d_rv2_out[1] * d_rv2[1] + d_rv2_out[2] * d_rv2[2];
      }
    }
  }

  return result;
}

/// Compute velocity kernel geometry sensitivity for a section.
Scalar ComputeVelocityKernelGeometrySensitivity(
    const VortexLattice& mesh,
    const Aircraft& aircraft,
    const std::vector<double>& dJ_dV,
    const VecX& gamma,
    const FlightCondition& condition,
    Index surface_idx,
    Index section_idx,
    GeometryVariable variable) {

  SectionStripMapping mapping = BuildSectionStripMapping(aircraft, mesh, surface_idx);
  if (mapping.num_sections == 0) return 0.0;

  const auto& surface_info = mesh.GetSurfaceInfo(surface_idx);
  Index num_chordwise = surface_info.num_chordwise;

  MeshArrays arrays = ExtractMeshArrays(mesh);
  double beta = std::sqrt(std::max(0.01, 1.0 - condition.mach * condition.mach));
  double core_chord_frac = mesh.GetVortexCoreRadiusChordFraction();
  double core_width_frac = mesh.GetVortexCoreRadiusWidthFraction();

  Scalar result = 0.0;

  if (section_idx < mapping.num_sections - 1) {
    auto [first_strip, num_strips] = mapping.GetStripRange(section_idx);
    if (num_strips > 0) {
      result += ComputeVelocityKernelGeometrySensitivityForInterval(
          mesh, arrays, dJ_dV, gamma,
          first_strip, num_strips, num_chordwise,
          core_chord_frac, core_width_frac, beta,
          variable, true);
    }
  }

  if (section_idx > 0) {
    auto [first_strip, num_strips] = mapping.GetStripRange(section_idx - 1);
    if (num_strips > 0) {
      result += ComputeVelocityKernelGeometrySensitivityForInterval(
          mesh, arrays, dJ_dV, gamma,
          first_strip, num_strips, num_chordwise,
          core_chord_frac, core_width_frac, beta,
          variable, false);
    }
  }

  return result;
}

/// Compute explicit force term ∂J/∂θ|_Γ,V for a geometry parameter.
AeroGradients ComputeExplicitForceTerm(
    const VortexLattice& mesh,
    const Aircraft& aircraft,
    const SolutionState& solution,
    const FlightCondition& condition,
    Index surface_idx,
    Index section_idx,
    GeometryVariable variable) {

  const Index n = mesh.GetNumVortices();
  if (n == 0) return AeroGradients{};

  SectionStripMapping mapping = BuildSectionStripMapping(aircraft, mesh, surface_idx);
  if (mapping.num_sections == 0) return AeroGradients{};

  const auto& surface_info = mesh.GetSurfaceInfo(surface_idx);
  Index num_chordwise = surface_info.num_chordwise;

  Vec3 moment_ref = aircraft.GetMomentReferencePoint();
  std::vector<double> bound_left(3 * n);
  std::vector<double> bound_right(3 * n);
  std::vector<double> moment_arm(3 * n);
  std::vector<double> induced_velocity(3 * n);

  for (Index i = 0; i < n; ++i) {
    const HorseshoeVortex& vor = mesh.GetVortex(i);
    bound_left[3*i + 0] = vor.bound_vortex_left.x();
    bound_left[3*i + 1] = vor.bound_vortex_left.y();
    bound_left[3*i + 2] = vor.bound_vortex_left.z();
    bound_right[3*i + 0] = vor.bound_vortex_right.x();
    bound_right[3*i + 1] = vor.bound_vortex_right.y();
    bound_right[3*i + 2] = vor.bound_vortex_right.z();
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

  std::vector<double> circulation(n);
  for (Index i = 0; i < n; ++i) {
    circulation[i] = solution.circulation(i);
  }

  // Initialize perturbation arrays to zero
  std::vector<double> d_bound_left(3 * n, 0.0);
  std::vector<double> d_bound_right(3 * n, 0.0);
  std::vector<double> d_moment_arm(3 * n, 0.0);

  auto compute_interval_perturbations = [&](
      Index first_strip, Index num_strips, bool is_left_section) {
    for (Index istrip = 0; istrip < num_strips; ++istrip) {
      Index strip_idx = first_strip + istrip;
      const Strip& strip = mesh.GetStrip(strip_idx);

      double f1, f2;
      ComputeStripEdgeFractions(istrip, num_strips, &f1, &f2);
      double fc = ComputeStripInterpolationFraction(istrip, num_strips);

      double w1 = is_left_section ? (1.0 - f1) : f1;
      double w2 = is_left_section ? (1.0 - f2) : f2;
      double wc = is_left_section ? (1.0 - fc) : fc;

      for (Index ichord = 0; ichord < num_chordwise; ++ichord) {
        Index vor_idx = strip.first_vortex_index + ichord;

        double dx = 1.0 / static_cast<double>(num_chordwise);
        double x1 = ichord * dx;
        double xvr = x1 + 0.25 * dx;
        double xcp = x1 + 0.75 * dx;

        switch (variable) {
          case GeometryVariable::kChord:
            d_bound_left[3*vor_idx + 0] += w1 * xvr;
            d_bound_right[3*vor_idx + 0] += w2 * xvr;
            d_moment_arm[3*vor_idx + 0] += wc * xcp;
            break;

          case GeometryVariable::kXle:
            d_bound_left[3*vor_idx + 0] += w1;
            d_bound_right[3*vor_idx + 0] += w2;
            d_moment_arm[3*vor_idx + 0] += wc;
            break;

          case GeometryVariable::kYle:
            d_bound_left[3*vor_idx + 1] += w1;
            d_bound_right[3*vor_idx + 1] += w2;
            d_moment_arm[3*vor_idx + 1] += wc;
            break;

          case GeometryVariable::kZle:
            d_bound_left[3*vor_idx + 2] += w1;
            d_bound_right[3*vor_idx + 2] += w2;
            d_moment_arm[3*vor_idx + 2] += wc;
            break;

          case GeometryVariable::kTwist:
            break;
        }
      }
    }
  };

  if (section_idx < mapping.num_sections - 1) {
    auto [first_strip, num_strips] = mapping.GetStripRange(section_idx);
    if (num_strips > 0) {
      compute_interval_perturbations(first_strip, num_strips, true);
    }
  }

  if (section_idx > 0) {
    auto [first_strip, num_strips] = mapping.GetStripRange(section_idx - 1);
    if (num_strips > 0) {
      compute_interval_perturbations(first_strip, num_strips, false);
    }
  }

  double d_result[6] = {0};
  double ref_area = aircraft.GetReferenceArea();

  analytical_force_geometry_gradient(
      condition.alpha, condition.beta,
      circulation.data(), induced_velocity.data(),
      bound_left.data(), bound_right.data(),
      moment_arm.data(),
      static_cast<int>(n), ref_area,
      d_bound_left.data(), d_bound_right.data(),
      d_moment_arm.data(),
      d_result);

  AeroGradients grads;
  grads.dCL = d_result[0];
  grads.dCD = d_result[1];
  grads.dCY = d_result[2];
  grads.dCl = d_result[3];
  grads.dCm = d_result[4];
  grads.dCn = d_result[5];
  grads.dCDi = 0.0;

  return grads;
}

/// Compute geometry gradient using the adjoint method with analytical derivatives.
AeroGradients ComputeSectionGradientAdjoint(
    const Aircraft& aircraft,
    const FlightCondition& condition,
    Index surface_idx,
    Index section_idx,
    GeometryVariable variable) {

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
    return AeroGradients{};
  }

  MatX rhs(n, 6);
  rhs.col(0) = force_grads.dCL_dGamma;
  rhs.col(1) = force_grads.dCD_dGamma;
  rhs.col(2) = force_grads.dCY_dGamma;
  rhs.col(3) = force_grads.dCl_dGamma;
  rhs.col(4) = force_grads.dCm_dGamma;
  rhs.col(5) = force_grads.dCn_dGamma;

  MatX lambda = adjoint.SolveAdjointMultiple(rhs);

  // Step 3: Compute RHS sensitivity ∂b/∂θ
  VecX db_dtheta = ComputeRHSGeometrySensitivity(
      mesh, aircraft, condition, surface_idx, section_idx, variable);

  // Step 4: Compute implicit term λᵀ·(∂b/∂θ)
  AeroGradients implicit_rhs;
  implicit_rhs.dCL = lambda.col(0).dot(db_dtheta);
  implicit_rhs.dCD = lambda.col(1).dot(db_dtheta);
  implicit_rhs.dCY = lambda.col(2).dot(db_dtheta);
  implicit_rhs.dCl = lambda.col(3).dot(db_dtheta);
  implicit_rhs.dCm = lambda.col(4).dot(db_dtheta);
  implicit_rhs.dCn = lambda.col(5).dot(db_dtheta);
  implicit_rhs.dCDi = 0.0;

  // Step 5: Compute AIC sensitivity term -λᵀ·(∂A/∂θ)·Γ
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
  aic_sens.dCDi = 0.0;

  // Step 6: Compute explicit term ∂J/∂θ|_Γ,V
  AeroGradients explicit_term = ComputeExplicitForceTerm(
      mesh, aircraft, solution, condition, surface_idx, section_idx, variable);

  // Step 7: Compute velocity kernel geometry sensitivity (∂J/∂V)·(∂K/∂θ)·Γ
  AeroGradients velocity_kernel_sens;
  velocity_kernel_sens.dCL = 0.0;
  velocity_kernel_sens.dCD = 0.0;
  velocity_kernel_sens.dCY = 0.0;
  velocity_kernel_sens.dCl = 0.0;
  velocity_kernel_sens.dCm = 0.0;
  velocity_kernel_sens.dCn = 0.0;
  velocity_kernel_sens.dCDi = 0.0;

  {
    Vec3 moment_ref = aircraft.GetMomentReferencePoint();
    std::vector<double> bound_left(3 * n);
    std::vector<double> bound_right(3 * n);
    std::vector<double> moment_arm(3 * n);
    std::vector<double> induced_velocity(3 * n);
    std::vector<double> circulation(n);

    for (Index i = 0; i < n; ++i) {
      const HorseshoeVortex& vor = mesh.GetVortex(i);
      bound_left[3*i + 0] = vor.bound_vortex_left.x();
      bound_left[3*i + 1] = vor.bound_vortex_left.y();
      bound_left[3*i + 2] = vor.bound_vortex_left.z();
      bound_right[3*i + 0] = vor.bound_vortex_right.x();
      bound_right[3*i + 1] = vor.bound_vortex_right.y();
      bound_right[3*i + 2] = vor.bound_vortex_right.z();
      moment_arm[3*i + 0] = vor.control_point.x() - moment_ref.x();
      moment_arm[3*i + 1] = vor.control_point.y() - moment_ref.y();
      moment_arm[3*i + 2] = vor.control_point.z() - moment_ref.z();
      if (solution.induced_velocity.rows() > i) {
        induced_velocity[3*i + 0] = solution.induced_velocity(i, 0);
        induced_velocity[3*i + 1] = solution.induced_velocity(i, 1);
        induced_velocity[3*i + 2] = solution.induced_velocity(i, 2);
      }
      circulation[i] = solution.circulation(i);
    }

    double ref_area = aircraft.GetReferenceArea();

    auto compute_velocity_kernel_term = [&](int output_index) -> Scalar {
      std::vector<double> dJ_dV(3 * n, 0.0);
      analytical_force_induced_velocity_gradient(
          condition.alpha, condition.beta,
          circulation.data(), induced_velocity.data(),
          bound_left.data(), bound_right.data(),
          nullptr, moment_arm.data(),
          static_cast<int>(n), ref_area,
          output_index, dJ_dV.data());

      return ComputeVelocityKernelGeometrySensitivity(
          mesh, aircraft, dJ_dV, solution.circulation, condition,
          surface_idx, section_idx, variable);
    };

    velocity_kernel_sens.dCL = compute_velocity_kernel_term(0);
    velocity_kernel_sens.dCD = compute_velocity_kernel_term(1);
    velocity_kernel_sens.dCY = compute_velocity_kernel_term(2);
    velocity_kernel_sens.dCl = compute_velocity_kernel_term(3);
    velocity_kernel_sens.dCm = compute_velocity_kernel_term(4);
    velocity_kernel_sens.dCn = compute_velocity_kernel_term(5);
  }

  // Total gradient = explicit + velocity_kernel + implicit_rhs + aic_sens
  AeroGradients result;
  result.dCL = explicit_term.dCL + velocity_kernel_sens.dCL + implicit_rhs.dCL + aic_sens.dCL;
  result.dCD = explicit_term.dCD + velocity_kernel_sens.dCD + implicit_rhs.dCD + aic_sens.dCD;
  result.dCY = explicit_term.dCY + velocity_kernel_sens.dCY + implicit_rhs.dCY + aic_sens.dCY;
  result.dCl = explicit_term.dCl + velocity_kernel_sens.dCl + implicit_rhs.dCl + aic_sens.dCl;
  result.dCm = explicit_term.dCm + velocity_kernel_sens.dCm + implicit_rhs.dCm + aic_sens.dCm;
  result.dCn = explicit_term.dCn + velocity_kernel_sens.dCn + implicit_rhs.dCn + aic_sens.dCn;
  result.dCDi = 0.0;  // TODO: Add Trefftz plane gradient

#ifndef NDEBUG
  if (std::getenv("VLM_DEBUG_ADJOINT")) {
    fprintf(stderr, "Adjoint terms for var=%d sec=%lld:\n",
            static_cast<int>(variable), static_cast<long long>(section_idx));
    fprintf(stderr, "  explicit:      dCL=%+.6f dCD=%+.6f dCm=%+.6f\n",
            explicit_term.dCL, explicit_term.dCD, explicit_term.dCm);
    fprintf(stderr, "  vel_kernel:    dCL=%+.6f dCD=%+.6f dCm=%+.6f\n",
            velocity_kernel_sens.dCL, velocity_kernel_sens.dCD, velocity_kernel_sens.dCm);
    fprintf(stderr, "  implicit_rhs:  dCL=%+.6f dCD=%+.6f dCm=%+.6f\n",
            implicit_rhs.dCL, implicit_rhs.dCD, implicit_rhs.dCm);
    fprintf(stderr, "  aic_sens:      dCL=%+.6f dCD=%+.6f dCm=%+.6f\n",
            aic_sens.dCL, aic_sens.dCD, aic_sens.dCm);
    fprintf(stderr, "  TOTAL:         dCL=%+.6f dCD=%+.6f dCm=%+.6f\n",
            result.dCL, result.dCD, result.dCm);
  }
#endif

  return result;
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

  return ComputeSectionGradientAdjoint(impl_->aircraft_, condition,
                                       surface_index, section_index,
                                       variable);
}

GeometryGradients GradientCalculator::ComputeGeometryGradients(
    const FlightCondition& condition) {
  GeometryGradients grads;

  if (!impl_->aircraft_set_) {
    return grads;
  }

  const Aircraft& aircraft = impl_->aircraft_;

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
