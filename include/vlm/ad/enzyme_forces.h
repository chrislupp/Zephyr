// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#pragma once

/// @file enzyme_forces.h
/// @brief Pure kernel functions for force computations that Enzyme can differentiate.
///
/// These functions compute aerodynamic forces from circulation and mesh data.
/// They are designed as pure functions with raw pointer interfaces suitable
/// for Enzyme automatic differentiation.
///
/// The key equations are:
/// - Kutta-Joukowski: F = ρ * Γ * (V × dl)
/// - Body-to-stability: CL = CZ*cos(α) - CX*sin(α), CD = CX*cos(α) + CZ*sin(α)

#include <cmath>

#include "vlm/ad/enzyme_interface.h"

namespace vlm {
namespace ad {

// ============================================================================
// Freestream Direction Computation
// ============================================================================

/// @brief Compute freestream velocity direction from alpha and beta.
///
/// @param alpha Angle of attack (radians)
/// @param beta Sideslip angle (radians)
/// @param v_inf Output freestream direction [3]
inline void compute_freestream_direction(
    double alpha,
    double beta,
    double* v_inf) {
  double ca = std::cos(alpha);
  double sa = std::sin(alpha);
  double cb = std::cos(beta);
  double sb = std::sin(beta);

  v_inf[0] = ca * cb;
  v_inf[1] = sb;
  v_inf[2] = sa * cb;
}

/// @brief Compute analytical gradient of freestream direction w.r.t. alpha.
///
/// ∂V/∂α = (-sin(α)cos(β), 0, cos(α)cos(β))
///
/// @param alpha Angle of attack (radians)
/// @param beta Sideslip angle (radians)
/// @param dv_dalpha Output gradient [3]
inline void compute_freestream_alpha_gradient(
    double alpha,
    double beta,
    double* dv_dalpha) {
  double ca = std::cos(alpha);
  double sa = std::sin(alpha);
  double cb = std::cos(beta);

  dv_dalpha[0] = -sa * cb;
  dv_dalpha[1] = 0.0;
  dv_dalpha[2] = ca * cb;
}

/// @brief Compute analytical gradient of freestream direction w.r.t. beta.
///
/// ∂V/∂β = (-cos(α)sin(β), cos(β), -sin(α)sin(β))
///
/// @param alpha Angle of attack (radians)
/// @param beta Sideslip angle (radians)
/// @param dv_dbeta Output gradient [3]
inline void compute_freestream_beta_gradient(
    double alpha,
    double beta,
    double* dv_dbeta) {
  double ca = std::cos(alpha);
  double sa = std::sin(alpha);
  double cb = std::cos(beta);
  double sb = std::sin(beta);

  dv_dbeta[0] = -ca * sb;
  dv_dbeta[1] = cb;
  dv_dbeta[2] = -sa * sb;
}

// ============================================================================
// Body-to-Stability Axis Transformation
// ============================================================================

/// @brief Transform body-axis forces to stability-axis (CL, CD).
///
/// @param cf_x Body-axis X force coefficient (positive forward)
/// @param cf_z Body-axis Z force coefficient (positive down in body, but we use up)
/// @param alpha Angle of attack (radians)
/// @param CL Output lift coefficient
/// @param CD Output drag coefficient
inline void body_to_stability(
    double cf_x,
    double cf_z,
    double alpha,
    double* CL,
    double* CD) {
  double ca = std::cos(alpha);
  double sa = std::sin(alpha);

  // Lift: perpendicular to freestream (positive up)
  *CL = cf_z * ca - cf_x * sa;

  // Drag: along freestream (positive aft)
  *CD = cf_x * ca + cf_z * sa;
}

// ============================================================================
// Single Vortex Force Computation
// ============================================================================

/// @brief Compute force from a single horseshoe vortex using Kutta-Joukowski.
///
/// @param gamma Vortex circulation strength
/// @param v_eff Effective velocity at bound vortex [3]
/// @param bound_left Left bound vortex endpoint [3]
/// @param bound_right Right bound vortex endpoint [3]
/// @param force Output force vector [3]
inline void compute_vortex_force(
    double gamma,
    const double* v_eff,
    const double* bound_left,
    const double* bound_right,
    double* force) {
  // Bound vortex vector (from left to right)
  double dl[3];
  dl[0] = bound_right[0] - bound_left[0];
  dl[1] = bound_right[1] - bound_left[1];
  dl[2] = bound_right[2] - bound_left[2];

  // Kutta-Joukowski: F = 2 * Γ * (V × dl) for unit dynamic pressure
  // Cross product: V × dl
  double cross[3];
  cross[0] = v_eff[1] * dl[2] - v_eff[2] * dl[1];
  cross[1] = v_eff[2] * dl[0] - v_eff[0] * dl[2];
  cross[2] = v_eff[0] * dl[1] - v_eff[1] * dl[0];

  force[0] = 2.0 * gamma * cross[0];
  force[1] = 2.0 * gamma * cross[1];
  force[2] = 2.0 * gamma * cross[2];
}

// ============================================================================
// Complete Force Coefficient Computation
// ============================================================================

/// @brief Compute total force coefficients from all vortices.
///
/// This is the main kernel for force computation that Enzyme will differentiate.
///
/// @param alpha Angle of attack (radians)
/// @param beta Sideslip angle (radians)
/// @param circulation Vortex circulation values [num_vortices]
/// @param induced_velocity Induced velocities at vortices [3 * num_vortices]
/// @param bound_left Left bound vortex endpoints [3 * num_vortices]
/// @param bound_right Right bound vortex endpoints [3 * num_vortices]
/// @param num_vortices Number of vortices
/// @param reference_area Reference area for normalization
/// @param result Output [CL, CD, CY, Cl, Cm, Cn] (6 values)
inline void compute_force_coefficients_kernel(
    double alpha,
    double beta,
    const double* circulation,
    const double* induced_velocity,
    const double* bound_left,
    const double* bound_right,
    const double* control_points,
    const double* moment_arm,  // control_point - moment_ref for each vortex [3*N]
    int num_vortices,
    double reference_area,
    double* result) {
  // Initialize outputs
  double cf_x = 0.0, cf_y = 0.0, cf_z = 0.0;
  double cm_x = 0.0, cm_y = 0.0, cm_z = 0.0;

  // Compute freestream direction
  double v_inf[3];
  compute_freestream_direction(alpha, beta, v_inf);

  // Accumulate forces from all vortices
  for (int i = 0; i < num_vortices; ++i) {
    // Effective velocity = freestream + induced
    double v_eff[3];
    v_eff[0] = v_inf[0] + induced_velocity[3 * i + 0];
    v_eff[1] = v_inf[1] + induced_velocity[3 * i + 1];
    v_eff[2] = v_inf[2] + induced_velocity[3 * i + 2];

    // Compute force from this vortex
    double force[3];
    compute_vortex_force(
        circulation[i],
        v_eff,
        &bound_left[3 * i],
        &bound_right[3 * i],
        force);

    // Accumulate body-axis forces
    cf_x += force[0];
    cf_y += force[1];
    cf_z += force[2];

    // Accumulate moments about reference point
    // M = r × F
    const double* r = &moment_arm[3 * i];
    cm_x += r[1] * force[2] - r[2] * force[1];
    cm_y += r[2] * force[0] - r[0] * force[2];
    cm_z += r[0] * force[1] - r[1] * force[0];
  }

  // Normalize by reference area
  cf_x /= reference_area;
  cf_y /= reference_area;
  cf_z /= reference_area;
  cm_x /= reference_area;
  cm_y /= reference_area;
  cm_z /= reference_area;

  // Transform to stability axes
  double CL, CD;
  body_to_stability(cf_x, cf_z, alpha, &CL, &CD);

  // Output results
  result[0] = CL;
  result[1] = CD;
  result[2] = cf_y;  // CY (side force)
  result[3] = -cm_x; // Cl (roll moment, sign convention)
  result[4] = cm_y;  // Cm (pitch moment)
  result[5] = -cm_z; // Cn (yaw moment, sign convention)
}

// ============================================================================
// Enzyme Gradient Functions
// ============================================================================

#ifdef VLM_USE_ENZYME

/// @brief Compute gradients of force coefficients w.r.t. alpha using Enzyme.
///
/// @param alpha Angle of attack
/// @param beta Sideslip angle
/// @param circulation Circulation values
/// @param induced_velocity Induced velocities
/// @param bound_left Left bound vortex points
/// @param bound_right Right bound vortex points
/// @param control_points Control point positions
/// @param moment_arm Moment arm vectors
/// @param num_vortices Number of vortices
/// @param reference_area Reference area
/// @param d_result Output gradients d(result)/d(alpha) [6]
void enzyme_force_alpha_gradient(
    double alpha,
    double beta,
    const double* circulation,
    const double* induced_velocity,
    const double* bound_left,
    const double* bound_right,
    const double* control_points,
    const double* moment_arm,
    int num_vortices,
    double reference_area,
    double* d_result);

/// @brief Compute gradients of force coefficients w.r.t. beta using Enzyme.
void enzyme_force_beta_gradient(
    double alpha,
    double beta,
    const double* circulation,
    const double* induced_velocity,
    const double* bound_left,
    const double* bound_right,
    const double* control_points,
    const double* moment_arm,
    int num_vortices,
    double reference_area,
    double* d_result);

/// @brief Compute gradients of force coefficients w.r.t. circulation using Enzyme.
///
/// This computes dJ/dΓ for use in the adjoint method.
///
/// @param alpha Angle of attack
/// @param ... other inputs same as above
/// @param output_seed Which output to differentiate (0=CL, 1=CD, etc.)
/// @param d_circulation Output: d(output)/d(Γ) for each vortex [num_vortices]
void enzyme_force_circulation_gradient(
    double alpha,
    double beta,
    const double* circulation,
    const double* induced_velocity,
    const double* bound_left,
    const double* bound_right,
    const double* control_points,
    const double* moment_arm,
    int num_vortices,
    double reference_area,
    int output_index,
    double* d_circulation);

#endif  // VLM_USE_ENZYME

// ============================================================================
// Finite Difference Fallback
// ============================================================================

/// @brief Compute force alpha gradient using finite differences.
void fd_force_alpha_gradient(
    double alpha,
    double beta,
    const double* circulation,
    const double* induced_velocity,
    const double* bound_left,
    const double* bound_right,
    const double* control_points,
    const double* moment_arm,
    int num_vortices,
    double reference_area,
    double* d_result,
    double h = 1e-7);

/// @brief Compute force beta gradient using finite differences.
void fd_force_beta_gradient(
    double alpha,
    double beta,
    const double* circulation,
    const double* induced_velocity,
    const double* bound_left,
    const double* bound_right,
    const double* control_points,
    const double* moment_arm,
    int num_vortices,
    double reference_area,
    double* d_result,
    double h = 1e-7);

}  // namespace ad
}  // namespace vlm
