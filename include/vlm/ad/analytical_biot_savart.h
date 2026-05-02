// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#pragma once

/// @file analytical_biot_savart.h
/// @brief Analytical Biot-Savart kernel and its vector-Jacobian products (VJPs).
///
/// Provides the primal horseshoe vortex velocity computation and hand-derived
/// analytical derivatives. Uses VJP form (seed^T * Jacobian) because callers
/// always contract with a seed vector (panel normal for AIC, force direction
/// for forces).
///
/// Core equations:
/// - Biot-Savart with Leishman R^4 core: V = (AxB)*T / (4pi)
/// - Trailing legs: semi-infinite along +x with R^4 regularization

#include <cmath>

namespace vlm {
namespace ad {

constexpr double kPi4Inv = 0.079577471545947667884;  // 1/(4*pi)
constexpr double kBSEps = 1e-12;

// ============================================================================
// Primal Kernel
// ============================================================================

/// @brief Compute induced velocity from a horseshoe vortex at a field point.
///
/// Uses Leishman R^4 core model for regularization and Prandtl-Glauert
/// compressibility correction on the x-component.
///
/// @param field_point Control point position [3]
/// @param bound_left Left endpoint of bound vortex [3]
/// @param bound_right Right endpoint of bound vortex [3]
/// @param prandtl_glauert sqrt(1-M^2)
/// @param core_radius Vortex core radius
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

  // Bound vortex contribution
  if (include_bound && Amag * Bmag > kBSEps) {
    double AxBx = Ay * Bz - Az * By;
    double AxBy = Az * Bx - Ax * Bz;
    double AxBz = Ax * By - Ay * Bx;
    double AxBsq = AxBx*AxBx + AxBy*AxBy + AxBz*AxBz;

    if (AxBsq > kBSEps) {
      double AdB = Ax*Bx + Ay*By + Az*Bz;
      double ALsq = Asq + Bsq - 2.0 * AdB;

      double T =
          ((Bsq - AdB) / std::sqrt(std::sqrt(Bsq * Bsq + core_4)) +
           (Asq - AdB) / std::sqrt(std::sqrt(Asq * Asq + core_4))) /
          std::sqrt(AxBsq * AxBsq + ALsq * ALsq * core_4);

      U = AxBx * T;
      V = AxBy * T;
      W = AxBz * T;
    }
  }

  // Trailing leg A (semi-infinite along +x)
  if (Amag > kBSEps) {
    double Axisq = Ay * Ay + Az * Az;
    double T = -(1.0 - Ax / Amag) / std::sqrt(Axisq * Axisq + core_4);
    V += Az * T;
    W -= Ay * T;
  }

  // Trailing leg B (semi-infinite along +x)
  if (Bmag > kBSEps) {
    double Bxisq = By * By + Bz * Bz;
    double T = (1.0 - Bx / Bmag) / std::sqrt(Bxisq * Bxisq + core_4);
    V += Bz * T;
    W -= By * T;
  }

  velocity[0] = U * kPi4Inv / prandtl_glauert;
  velocity[1] = V * kPi4Inv;
  velocity[2] = W * kPi4Inv;
}

/// @brief Compute AIC element: AIC = V(field_point, rv1, rv2) . normal
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
  return velocity[0]*normal[0] + velocity[1]*normal[1] + velocity[2]*normal[2];
}

// ============================================================================
// Vector-Jacobian Products (VJPs)
// ============================================================================

/// @brief Compute VJP of horseshoe velocity w.r.t. field point.
///
/// Returns: result[k] = sum_m seed[m] * dV[m]/d(field_point[k])
///
/// @param field_point Control point position [3]
/// @param bound_left Left bound vortex endpoint [3]
/// @param bound_right Right bound vortex endpoint [3]
/// @param prandtl_glauert sqrt(1-M^2)
/// @param core_radius Vortex core radius
/// @param include_bound Whether to include bound vortex
/// @param seed Seed vector [3] (e.g., panel normal)
/// @param result Output VJP [3]
void horseshoe_velocity_vjp_field_point(
    const double* field_point,
    const double* bound_left,
    const double* bound_right,
    double prandtl_glauert,
    double core_radius,
    bool include_bound,
    const double* seed,
    double* result);

/// @brief Compute VJP of horseshoe velocity w.r.t. bound_left.
///
/// Returns: result[k] = sum_m seed[m] * dV[m]/d(bound_left[k])
void horseshoe_velocity_vjp_bound_left(
    const double* field_point,
    const double* bound_left,
    const double* bound_right,
    double prandtl_glauert,
    double core_radius,
    bool include_bound,
    const double* seed,
    double* result);

/// @brief Compute VJP of horseshoe velocity w.r.t. bound_right.
///
/// Returns: result[k] = sum_m seed[m] * dV[m]/d(bound_right[k])
void horseshoe_velocity_vjp_bound_right(
    const double* field_point,
    const double* bound_left,
    const double* bound_right,
    double prandtl_glauert,
    double core_radius,
    bool include_bound,
    const double* seed,
    double* result);

}  // namespace ad
}  // namespace vlm
