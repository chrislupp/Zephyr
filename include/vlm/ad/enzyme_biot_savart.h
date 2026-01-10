// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#pragma once

/// @file enzyme_biot_savart.h
/// @brief Pure kernel functions for Biot-Savart velocity computation that Enzyme can differentiate.
///
/// These functions compute induced velocities from horseshoe vortices.
/// They are designed as pure functions with raw pointer interfaces suitable
/// for Enzyme automatic differentiation.
///
/// The key equations are:
/// - Biot-Savart: V = (Γ/4π) * ∫(dl × r)/|r|³
/// - Leishman core model: regularization using r⁴ + rcore⁴

#include <cmath>

#include "vlm/ad/enzyme_interface.h"

namespace vlm {
namespace ad {

// 1/(4*pi)
constexpr double kPi4InvBiotSavart = 0.079577471545947667884;

// ============================================================================
// Core Biot-Savart Kernel
// ============================================================================

/// @brief Compute induced velocity from a horseshoe vortex at a field point.
///
/// This is the core Biot-Savart computation that Enzyme will differentiate.
/// It computes the velocity induced by:
/// - A bound vortex segment from bound_left to bound_right
/// - Two trailing legs extending to +∞ in the x-direction
///
/// @param field_point Field point where velocity is computed [3]
/// @param bound_left Left endpoint of bound vortex [3]
/// @param bound_right Right endpoint of bound vortex [3]
/// @param prandtl_glauert Prandtl-Glauert compressibility factor sqrt(1-M²)
/// @param core_radius Vortex core radius for regularization
/// @param include_bound Whether to include bound vortex contribution
/// @param velocity Output induced velocity [3] (for unit circulation)
inline void compute_horseshoe_velocity_kernel(
    const double* field_point,
    const double* bound_left,
    const double* bound_right,
    double prandtl_glauert,
    double core_radius,
    bool include_bound,
    double* velocity) {

  constexpr double kEps = 1e-12;

  // Prandtl-Glauert transformed vectors from vortex endpoints to field point
  // A = bound_left - field_point (transformed)
  // B = bound_right - field_point (transformed)
  double Ax = (bound_left[0] - field_point[0]) / prandtl_glauert;
  double Ay = bound_left[1] - field_point[1];
  double Az = bound_left[2] - field_point[2];

  double Bx = (bound_right[0] - field_point[0]) / prandtl_glauert;
  double By = bound_right[1] - field_point[1];
  double Bz = bound_right[2] - field_point[2];

  double Asq = Ax*Ax + Ay*Ay + Az*Az;
  double Bsq = Bx*Bx + By*By + Bz*Bz;
  double Amag = std::sqrt(Asq);
  double Bmag = std::sqrt(Bsq);

  double core_sq = core_radius * core_radius;
  double core_4 = core_sq * core_sq;

  double U = 0.0, V = 0.0, W = 0.0;

  // Contribution from the transverse bound leg
  if (include_bound && Amag * Bmag > kEps) {
    // Cross product A × B
    double AxBx = Ay * Bz - Az * By;
    double AxBy = Az * Bx - Ax * Bz;
    double AxBz = Ax * By - Ay * Bx;
    double AxBsq = AxBx*AxBx + AxBy*AxBy + AxBz*AxBz;

    if (AxBsq > kEps) {
      double AdB = Ax*Bx + Ay*By + Az*Bz;
      double ALsq = Asq + Bsq - 2.0 * AdB;

      // Leishman R^4 core model
      double T =
          ((Bsq - AdB) / std::sqrt(std::sqrt(Bsq * Bsq + core_4)) +
           (Asq - AdB) / std::sqrt(std::sqrt(Asq * Asq + core_4))) /
          std::sqrt(AxBsq * AxBsq + ALsq * ALsq * core_4);

      U = AxBx * T;
      V = AxBy * T;
      W = AxBz * T;
    }
  }

  // Contribution from trailing leg attached to A (semi-infinite along +x)
  if (Amag > kEps) {
    double Axisq = Ay * Ay + Az * Az;
    double Adx = Ax;

    // Leishman core model for trailing leg
    double T = -(1.0 - Adx / Amag) / std::sqrt(Axisq * Axisq + core_4);

    V += Az * T;
    W -= Ay * T;
  }

  // Contribution from trailing leg attached to B (semi-infinite along +x)
  if (Bmag > kEps) {
    double Bxisq = By * By + Bz * Bz;
    double Bdx = Bx;

    // Leishman core model for trailing leg
    double T = (1.0 - Bdx / Bmag) / std::sqrt(Bxisq * Bxisq + core_4);

    V += Bz * T;
    W -= By * T;
  }

  // Apply 1/(4*pi) factor and Prandtl-Glauert correction to U
  velocity[0] = U * kPi4InvBiotSavart / prandtl_glauert;
  velocity[1] = V * kPi4InvBiotSavart;
  velocity[2] = W * kPi4InvBiotSavart;
}

/// @brief Wrapper without include_bound flag (always includes bound vortex)
inline void compute_horseshoe_velocity_kernel_with_bound(
    const double* field_point,
    const double* bound_left,
    const double* bound_right,
    double prandtl_glauert,
    double core_radius,
    double* velocity) {
  compute_horseshoe_velocity_kernel(field_point, bound_left, bound_right,
                                    prandtl_glauert, core_radius, true, velocity);
}

/// @brief Wrapper without include_bound flag (excludes bound vortex - for self-influence)
inline void compute_horseshoe_velocity_kernel_no_bound(
    const double* field_point,
    const double* bound_left,
    const double* bound_right,
    double prandtl_glauert,
    double core_radius,
    double* velocity) {
  compute_horseshoe_velocity_kernel(field_point, bound_left, bound_right,
                                    prandtl_glauert, core_radius, false, velocity);
}

// ============================================================================
// AIC Matrix Element Computation
// ============================================================================

/// @brief Compute single AIC matrix element (normalwash influence coefficient).
///
/// AIC(i,j) = n_i · V_ij where V_ij is velocity at i induced by vortex j.
///
/// @param field_point Control point position [3]
/// @param normal Normal vector at control point [3]
/// @param bound_left Left endpoint of bound vortex [3]
/// @param bound_right Right endpoint of bound vortex [3]
/// @param prandtl_glauert Prandtl-Glauert factor
/// @param core_radius Core radius
/// @param include_bound Whether to include bound vortex
/// @return AIC matrix element
inline double compute_aic_element_kernel(
    const double* field_point,
    const double* normal,
    const double* bound_left,
    const double* bound_right,
    double prandtl_glauert,
    double core_radius,
    bool include_bound) {

  double velocity[3];
  compute_horseshoe_velocity_kernel(field_point, bound_left, bound_right,
                                    prandtl_glauert, core_radius, include_bound,
                                    velocity);

  // Dot product with normal
  return velocity[0] * normal[0] + velocity[1] * normal[1] + velocity[2] * normal[2];
}

// ============================================================================
// Enzyme Gradient Functions
// ============================================================================

#ifdef VLM_USE_ENZYME

/// @brief Compute gradient of AIC element w.r.t. field point using Enzyme.
///
/// @param field_point Control point position [3]
/// @param normal Normal vector at control point [3]
/// @param bound_left Left endpoint of bound vortex [3]
/// @param bound_right Right endpoint of bound vortex [3]
/// @param prandtl_glauert Prandtl-Glauert factor
/// @param core_radius Core radius
/// @param include_bound Whether to include bound vortex
/// @param d_field_point Output gradient d(AIC)/d(field_point) [3]
void enzyme_aic_field_point_gradient(
    const double* field_point,
    const double* normal,
    const double* bound_left,
    const double* bound_right,
    double prandtl_glauert,
    double core_radius,
    bool include_bound,
    double* d_field_point);

/// @brief Compute gradient of AIC element w.r.t. bound vortex endpoints using Enzyme.
///
/// @param field_point Control point position [3]
/// @param normal Normal vector at control point [3]
/// @param bound_left Left endpoint of bound vortex [3]
/// @param bound_right Right endpoint of bound vortex [3]
/// @param prandtl_glauert Prandtl-Glauert factor
/// @param core_radius Core radius
/// @param include_bound Whether to include bound vortex
/// @param d_bound_left Output gradient d(AIC)/d(bound_left) [3]
/// @param d_bound_right Output gradient d(AIC)/d(bound_right) [3]
void enzyme_aic_bound_gradient(
    const double* field_point,
    const double* normal,
    const double* bound_left,
    const double* bound_right,
    double prandtl_glauert,
    double core_radius,
    bool include_bound,
    double* d_bound_left,
    double* d_bound_right);

/// @brief Compute gradient of AIC element w.r.t. normal vector using Enzyme.
void enzyme_aic_normal_gradient(
    const double* field_point,
    const double* normal,
    const double* bound_left,
    const double* bound_right,
    double prandtl_glauert,
    double core_radius,
    bool include_bound,
    double* d_normal);

#endif  // VLM_USE_ENZYME

// ============================================================================
// Finite Difference Fallbacks
// ============================================================================

/// @brief Compute gradient of AIC element w.r.t. field point using FD.
void fd_aic_field_point_gradient(
    const double* field_point,
    const double* normal,
    const double* bound_left,
    const double* bound_right,
    double prandtl_glauert,
    double core_radius,
    bool include_bound,
    double* d_field_point,
    double h = 1e-7);

/// @brief Compute gradient of AIC element w.r.t. bound endpoints using FD.
void fd_aic_bound_gradient(
    const double* field_point,
    const double* normal,
    const double* bound_left,
    const double* bound_right,
    double prandtl_glauert,
    double core_radius,
    bool include_bound,
    double* d_bound_left,
    double* d_bound_right,
    double h = 1e-7);

}  // namespace ad
}  // namespace vlm
