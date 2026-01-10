// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#pragma once

/// @file enzyme_aic.h
/// @brief Pure kernel functions for AIC matrix computation that Enzyme can differentiate.
///
/// These functions compute the Aerodynamic Influence Coefficient (AIC) matrix,
/// which relates vortex circulations to induced velocities at control points.
/// The kernels are designed for Enzyme automatic differentiation.
///
/// The AIC matrix element A_ij represents the normal component of velocity
/// induced at control point i by a unit-circulation horseshoe vortex j.
///
/// Key equations (from OptVL/AVL):
/// - Biot-Savart law with Lamb-Oseen core regularization
/// - Horseshoe vortex = bound vortex + two semi-infinite trailing vortices
/// - Prandtl-Glauert compressibility correction

#include <cmath>

#include "vlm/ad/enzyme_interface.h"

namespace vlm {
namespace ad {

// ============================================================================
// Constants
// ============================================================================

constexpr double kPiAIC = 3.14159265358979323846;
constexpr double kFourPiInv = 0.07957747154594767;  // 1/(4*pi)
constexpr double kEpsAIC = 1e-12;

// ============================================================================
// Core Biot-Savart Kernel (matches OptVL's VVOR)
// ============================================================================

/// @brief Compute induced velocity from a single vortex segment.
///
/// This implements the Biot-Savart law for a finite vortex segment
/// with Lamb-Oseen core regularization, matching OptVL's implementation.
///
/// @param r1 Start point of vortex segment [3]
/// @param r2 End point of vortex segment [3]
/// @param rc Control point [3]
/// @param gamma Circulation strength
/// @param core_radius Core radius for regularization
/// @param velocity Output induced velocity [3] (accumulated)
inline void compute_aic_segment_velocity(
    const double* r1,
    const double* r2,
    const double* rc,
    double gamma,
    double core_radius,
    double* velocity) {

  // Vectors from vortex endpoints to control point
  double a1 = rc[0] - r1[0];
  double a2 = rc[1] - r1[1];
  double a3 = rc[2] - r1[2];

  double b1 = rc[0] - r2[0];
  double b2 = rc[1] - r2[1];
  double b3 = rc[2] - r2[2];

  // Magnitudes
  double a_mag_sq = a1*a1 + a2*a2 + a3*a3;
  double b_mag_sq = b1*b1 + b2*b2 + b3*b3;
  double a_mag = std::sqrt(a_mag_sq);
  double b_mag = std::sqrt(b_mag_sq);

  // Cross product a x b
  double axb1 = a2*b3 - a3*b2;
  double axb2 = a3*b1 - a1*b3;
  double axb3 = a1*b2 - a2*b1;

  double axb_sq = axb1*axb1 + axb2*axb2 + axb3*axb3;

  // Core regularization (Lamb-Oseen model)
  double core_sq = core_radius * core_radius;
  double denom = axb_sq + core_sq * a_mag_sq * b_mag_sq / (a_mag_sq + b_mag_sq + kEpsAIC);

  if (denom < kEpsAIC || a_mag < kEpsAIC || b_mag < kEpsAIC) {
    return;  // Singular, no contribution
  }

  // Vortex segment direction
  double s1 = r2[0] - r1[0];
  double s2 = r2[1] - r1[1];
  double s3 = r2[2] - r1[2];

  // Dot products for the Biot-Savart factor
  double a_dot_s = a1*s1 + a2*s2 + a3*s3;
  double b_dot_s = b1*s1 + b2*s2 + b3*s3;

  // Biot-Savart factor: (1/4π) * (a·s/|a| - b·s/|b|) / |a×b|²
  double factor = kFourPiInv * gamma * (a_dot_s/a_mag - b_dot_s/b_mag) / denom;

  // Accumulate velocity
  velocity[0] += factor * axb1;
  velocity[1] += factor * axb2;
  velocity[2] += factor * axb3;
}

/// @brief Compute induced velocity from a semi-infinite trailing vortex.
///
/// A trailing vortex extends from a point to infinity in the +x direction.
/// This uses the standard Biot-Savart formula for a semi-infinite vortex.
///
/// @param r_start Start point of trailing vortex [3]
/// @param rc Control point [3]
/// @param gamma Circulation strength (positive = +x trailing direction)
/// @param core_radius Core radius for regularization
/// @param velocity Output induced velocity [3] (accumulated)
inline void compute_aic_trailing_velocity(
    const double* r_start,
    const double* rc,
    double gamma,
    double core_radius,
    double* velocity) {

  // Vector from vortex start to control point
  double r1 = rc[0] - r_start[0];
  double r2 = rc[1] - r_start[1];
  double r3 = rc[2] - r_start[2];

  // Distance in the y-z plane (perpendicular to trailing direction)
  double r_perp_sq = r2*r2 + r3*r3;
  double core_sq = core_radius * core_radius;

  // Regularized denominator
  double denom = r_perp_sq + core_sq;
  if (denom < kEpsAIC) {
    return;
  }

  // For semi-infinite vortex in +x direction:
  // V = (Γ/4π) * (1 + x/|r|) * (r_perp × x̂) / |r_perp|²
  double r_mag = std::sqrt(r1*r1 + r_perp_sq);
  if (r_mag < kEpsAIC) {
    return;
  }

  double factor = kFourPiInv * gamma * (1.0 + r1/r_mag) / denom;

  // r_perp × x̂ = (0, -r3, r2)
  velocity[1] -= factor * r3;
  velocity[2] += factor * r2;
}

/// @brief Compute induced velocity from a complete horseshoe vortex.
///
/// A horseshoe vortex consists of:
/// - A bound vortex from rv1 to rv2
/// - A trailing vortex from rv1 to +∞ (negative circulation)
/// - A trailing vortex from rv2 to +∞ (positive circulation)
///
/// @param rv1 Left bound vortex endpoint [3]
/// @param rv2 Right bound vortex endpoint [3]
/// @param rc Control point [3]
/// @param gamma Circulation strength
/// @param core_radius Core radius for regularization
/// @param include_bound Whether to include bound vortex contribution
/// @param velocity Output induced velocity [3]
inline void compute_aic_horseshoe_velocity(
    const double* rv1,
    const double* rv2,
    const double* rc,
    double gamma,
    double core_radius,
    bool include_bound,
    double* velocity) {

  velocity[0] = 0.0;
  velocity[1] = 0.0;
  velocity[2] = 0.0;

  // Bound vortex contribution (rv1 → rv2)
  if (include_bound) {
    compute_aic_segment_velocity(rv1, rv2, rc, gamma, core_radius, velocity);
  }

  // Left trailing vortex (rv1 → +∞, negative circulation)
  compute_aic_trailing_velocity(rv1, rc, -gamma, core_radius, velocity);

  // Right trailing vortex (rv2 → +∞, positive circulation)
  compute_aic_trailing_velocity(rv2, rc, gamma, core_radius, velocity);
}

/// @brief Compute one element of the AIC matrix.
///
/// This computes the normal component of velocity induced at control point i
/// by a unit-circulation horseshoe vortex j.
///
/// @param rv1 Left bound vortex endpoint [3]
/// @param rv2 Right bound vortex endpoint [3]
/// @param rc Control point [3]
/// @param normal Panel normal at control point [3]
/// @param core_radius Core radius for regularization
/// @param beta Prandtl-Glauert factor sqrt(1-M²)
/// @param include_bound Whether to include bound vortex
/// @return AIC matrix element (normalwash per unit circulation)
inline double compute_aic_element(
    const double* rv1,
    const double* rv2,
    const double* rc,
    const double* normal,
    double core_radius,
    double beta,
    bool include_bound) {

  // Compute induced velocity from unit circulation
  double velocity[3];
  compute_aic_horseshoe_velocity(rv1, rv2, rc, 1.0, core_radius,
                                  include_bound, velocity);

  // Apply Prandtl-Glauert compressibility correction
  // (stretch y,z coordinates by 1/beta)
  if (std::abs(beta - 1.0) > kEpsAIC) {
    double beta_inv = 1.0 / beta;
    velocity[1] *= beta_inv;
    velocity[2] *= beta_inv;
  }

  // Return normal component
  return velocity[0]*normal[0] + velocity[1]*normal[1] + velocity[2]*normal[2];
}

// ============================================================================
// AIC Matrix Row Computation (for Enzyme differentiation)
// ============================================================================

/// @brief Compute one row of the AIC matrix.
///
/// This computes all influence coefficients for one control point.
/// Designed for Enzyme differentiation w.r.t. geometry parameters.
///
/// @param rv1_all Left bound vortex endpoints [3*num_vortices]
/// @param rv2_all Right bound vortex endpoints [3*num_vortices]
/// @param rc Control point [3]
/// @param normal Panel normal [3]
/// @param chord_all Chord values for each vortex [num_vortices]
/// @param surface_idx_all Surface index for each vortex [num_vortices]
/// @param num_vortices Number of vortices
/// @param control_surface_idx Surface index of the control point
/// @param core_chord_fraction Core radius as fraction of chord
/// @param core_width_fraction Core radius as fraction of strip width
/// @param beta Prandtl-Glauert factor
/// @param aic_row Output AIC row [num_vortices]
inline void compute_aic_row_kernel(
    const double* rv1_all,
    const double* rv2_all,
    const double* rc,
    const double* normal,
    const double* chord_all,
    const int* surface_idx_all,
    int num_vortices,
    int control_surface_idx,
    double core_chord_fraction,
    double core_width_fraction,
    double beta,
    double* aic_row) {

  for (int j = 0; j < num_vortices; ++j) {
    const double* rv1 = &rv1_all[3*j];
    const double* rv2 = &rv2_all[3*j];

    // Compute strip width for core radius
    double dy = rv2[1] - rv1[1];
    double dz = rv2[2] - rv1[2];
    double strip_width = std::sqrt(dy*dy + dz*dz);

    // Core radius (different for same/different surfaces)
    double core_radius = 0.0;
    if (surface_idx_all[j] != control_surface_idx) {
      double chord_core = core_chord_fraction * chord_all[j];
      double width_core = core_width_fraction * strip_width;
      core_radius = (chord_core > width_core) ? chord_core : width_core;
    }

    // Include bound vortex only for different surfaces
    bool include_bound = (surface_idx_all[j] != control_surface_idx);

    aic_row[j] = compute_aic_element(rv1, rv2, rc, normal, core_radius,
                                      beta, include_bound);
  }
}

/// @brief Compute the full AIC matrix.
///
/// This is the main kernel for AIC computation, matching OptVL's VVOR.
///
/// @param rv1_all Left bound vortex endpoints [3*num_vortices]
/// @param rv2_all Right bound vortex endpoints [3*num_vortices]
/// @param rc_all Control points [3*num_vortices]
/// @param normal_all Panel normals [3*num_vortices]
/// @param chord_all Chord values [num_vortices]
/// @param surface_idx_all Surface index for each vortex [num_vortices]
/// @param num_vortices Number of vortices
/// @param core_chord_fraction Core radius as fraction of chord
/// @param core_width_fraction Core radius as fraction of strip width
/// @param beta Prandtl-Glauert factor
/// @param aic_matrix Output AIC matrix [num_vortices * num_vortices] (row-major)
inline void compute_aic_matrix_kernel(
    const double* rv1_all,
    const double* rv2_all,
    const double* rc_all,
    const double* normal_all,
    const double* chord_all,
    const int* surface_idx_all,
    int num_vortices,
    double core_chord_fraction,
    double core_width_fraction,
    double beta,
    double* aic_matrix) {

  for (int i = 0; i < num_vortices; ++i) {
    const double* rc = &rc_all[3*i];
    const double* normal = &normal_all[3*i];
    double* aic_row = &aic_matrix[i * num_vortices];

    compute_aic_row_kernel(
        rv1_all, rv2_all, rc, normal, chord_all, surface_idx_all,
        num_vortices, surface_idx_all[i],
        core_chord_fraction, core_width_fraction, beta, aic_row);
  }
}

// ============================================================================
// Enzyme Gradient Declarations
// ============================================================================

#ifdef VLM_USE_ENZYME

/// @brief Compute gradient of AIC element w.r.t. vortex positions using Enzyme.
///
/// @param rv1 Left bound vortex endpoint [3]
/// @param rv2 Right bound vortex endpoint [3]
/// @param rc Control point [3]
/// @param normal Panel normal [3]
/// @param core_radius Core radius
/// @param beta Prandtl-Glauert factor
/// @param include_bound Include bound vortex
/// @param d_rv1 Output gradient w.r.t. rv1 [3]
/// @param d_rv2 Output gradient w.r.t. rv2 [3]
/// @param d_rc Output gradient w.r.t. rc [3]
/// @param d_normal Output gradient w.r.t. normal [3]
void enzyme_aic_element_gradient(
    const double* rv1,
    const double* rv2,
    const double* rc,
    const double* normal,
    double core_radius,
    double beta,
    bool include_bound,
    double* d_rv1,
    double* d_rv2,
    double* d_rc,
    double* d_normal);

/// @brief Compute gradient of AIC row w.r.t. control point position using Enzyme.
///
/// @param rv1_all Left bound vortex endpoints [3*num_vortices]
/// @param rv2_all Right bound vortex endpoints [3*num_vortices]
/// @param rc Control point [3]
/// @param normal Panel normal [3]
/// @param chord_all Chord values [num_vortices]
/// @param surface_idx_all Surface indices [num_vortices]
/// @param num_vortices Number of vortices
/// @param control_surface_idx Control point surface index
/// @param core_chord_fraction Core radius chord fraction
/// @param core_width_fraction Core radius width fraction
/// @param beta Prandtl-Glauert factor
/// @param d_aic_row Seed gradient for AIC row [num_vortices]
/// @param d_rc Output gradient w.r.t. control point [3]
/// @param d_normal Output gradient w.r.t. normal [3]
void enzyme_aic_row_gradient_rc(
    const double* rv1_all,
    const double* rv2_all,
    const double* rc,
    const double* normal,
    const double* chord_all,
    const int* surface_idx_all,
    int num_vortices,
    int control_surface_idx,
    double core_chord_fraction,
    double core_width_fraction,
    double beta,
    const double* d_aic_row,
    double* d_rc,
    double* d_normal);

/// @brief Compute gradient of AIC matrix column w.r.t. vortex positions.
///
/// For geometry sensitivity, we need ∂A/∂θ where θ affects vortex positions.
/// This computes the contribution from one vortex j to all control points.
///
/// @param rv1_all Left bound vortex endpoints [3*num_vortices]
/// @param rv2_all Right bound vortex endpoints [3*num_vortices]
/// @param rc_all Control points [3*num_vortices]
/// @param normal_all Panel normals [3*num_vortices]
/// @param chord_all Chord values [num_vortices]
/// @param surface_idx_all Surface indices [num_vortices]
/// @param num_vortices Number of vortices
/// @param vortex_idx Index of vortex being differentiated
/// @param core_chord_fraction Core radius chord fraction
/// @param core_width_fraction Core radius width fraction
/// @param beta Prandtl-Glauert factor
/// @param d_aic_column Seed gradient for AIC column [num_vortices]
/// @param d_rv1 Output gradient w.r.t. rv1 of vortex [3]
/// @param d_rv2 Output gradient w.r.t. rv2 of vortex [3]
void enzyme_aic_column_gradient_rv(
    const double* rv1_all,
    const double* rv2_all,
    const double* rc_all,
    const double* normal_all,
    const double* chord_all,
    const int* surface_idx_all,
    int num_vortices,
    int vortex_idx,
    double core_chord_fraction,
    double core_width_fraction,
    double beta,
    const double* d_aic_column,
    double* d_rv1,
    double* d_rv2);

#endif  // VLM_USE_ENZYME

// ============================================================================
// Finite Difference Fallback
// ============================================================================

/// @brief Compute AIC element gradient using finite differences (fallback).
void fd_aic_element_gradient(
    const double* rv1,
    const double* rv2,
    const double* rc,
    const double* normal,
    double core_radius,
    double beta,
    bool include_bound,
    double* d_rv1,
    double* d_rv2,
    double* d_rc,
    double* d_normal,
    double h = 1e-7);

}  // namespace ad
}  // namespace vlm
