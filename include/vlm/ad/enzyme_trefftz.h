// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#pragma once

/// @file enzyme_trefftz.h
/// @brief Pure kernel functions for Trefftz plane drag computation that Enzyme can differentiate.
///
/// These functions compute far-field induced drag using the Trefftz plane method.
/// The Trefftz plane is located infinitely far downstream and contains the wake.
/// Drag is computed from the kinetic energy of the crossflow in this plane.
///
/// Key equations:
/// - 2D Biot-Savart: V = (Γ/2π) * r⊥/|r|² (with core regularization)
/// - Far-field lift: CL_ff = Σ 2*Γ*dy / S_ref
/// - Far-field drag: CDi_ff = Σ Γ*(dz*vy - dy*vz) / S_ref
/// - Span efficiency: e = CL² / (π * AR * CDi)

#include <cmath>

#include "vlm/ad/enzyme_interface.h"

namespace vlm {
namespace ad {

// 1/(2*pi) for 2D Biot-Savart
constexpr double kHalfPiInvTrefftz = 0.15915494309189533577;
constexpr double kPiTrefftz = 3.14159265358979323846;
constexpr double kEpsTrefftz = 1e-12;

// ============================================================================
// Trefftz Plane Drag Kernel
// ============================================================================

/// @brief Compute 2D Biot-Savart induced velocity in Trefftz plane.
///
/// @param gamma Circulation strength
/// @param yv Y-coordinate of vortex
/// @param zv Z-coordinate of vortex
/// @param yc Y-coordinate of control point
/// @param zc Z-coordinate of control point
/// @param core_sq Square of core radius
/// @param vy Output Y velocity component
/// @param vz Output Z velocity component
inline void compute_trefftz_velocity_point(
    double gamma,
    double yv, double zv,
    double yc, double zc,
    double core_sq,
    double* vy, double* vz) {
  double dy = yc - yv;
  double dz = zc - zv;
  double rsq = dy * dy + dz * dz;

  // Core regularization (R^2 model)
  double denom = rsq + core_sq;

  if (denom > kEpsTrefftz) {
    // V = (Γ/2π) * r_perp / |r|²
    // where r_perp = (-dz, dy) for circulation pointing in +x
    *vy += kHalfPiInvTrefftz * gamma * dz / denom;
    *vz -= kHalfPiInvTrefftz * gamma * dy / denom;
  }
}

/// @brief Compute Trefftz plane drag coefficients.
///
/// This is the main kernel for Trefftz plane drag computation.
///
/// @param strip_gamma Circulation per strip [num_strips]
/// @param strip_te_left Left trailing edge points [3*num_strips]
/// @param strip_te_right Right trailing edge points [3*num_strips]
/// @param strip_chord Strip chord lengths [num_strips]
/// @param strip_surface_index Surface index per strip [num_strips]
/// @param num_strips Number of strips
/// @param reference_area Reference area
/// @param reference_span Reference span
/// @param core_chord_fraction Core radius as fraction of chord
/// @param core_width_fraction Core radius as fraction of strip width
/// @param CL_ff Output far-field lift coefficient
/// @param CY_ff Output far-field side force coefficient
/// @param CDi_ff Output far-field induced drag coefficient
inline void compute_trefftz_drag_kernel(
    const double* strip_gamma,
    const double* strip_te_left,
    const double* strip_te_right,
    const double* strip_chord,
    const int* strip_surface_index,
    int num_strips,
    double reference_area,
    double reference_span,
    double core_chord_fraction,
    double core_width_fraction,
    double* CL_ff,
    double* CY_ff,
    double* CDi_ff) {

  *CL_ff = 0.0;
  *CY_ff = 0.0;
  *CDi_ff = 0.0;

  // Precompute strip centers and widths
  for (int jc = 0; jc < num_strips; ++jc) {
    // Strip geometry
    double y_left = strip_te_left[3 * jc + 1];
    double z_left = strip_te_left[3 * jc + 2];
    double y_right = strip_te_right[3 * jc + 1];
    double z_right = strip_te_right[3 * jc + 2];

    double dy_c = y_right - y_left;
    double dz_c = z_right - z_left;
    double ds_c = std::sqrt(dy_c * dy_c + dz_c * dz_c);

    if (ds_c < kEpsTrefftz) continue;

    // Control point (strip center)
    double yc = 0.5 * (y_left + y_right);
    double zc = 0.5 * (z_left + z_right);

    // Compute induced velocity from all other strips
    double vy = 0.0, vz = 0.0;

    for (int jv = 0; jv < num_strips; ++jv) {
      double y_left_v = strip_te_left[3 * jv + 1];
      double z_left_v = strip_te_left[3 * jv + 2];
      double y_right_v = strip_te_right[3 * jv + 1];
      double z_right_v = strip_te_right[3 * jv + 2];

      double dy_v = y_right_v - y_left_v;
      double dz_v = z_right_v - z_left_v;
      double ds_v = std::sqrt(dy_v * dy_v + dz_v * dz_v);

      if (ds_v < kEpsTrefftz) continue;

      // Core radius (only for different surfaces)
      double core_radius = 0.0;
      if (strip_surface_index[jc] != strip_surface_index[jv]) {
        double chord_core = core_chord_fraction * strip_chord[jv];
        double width_core = core_width_fraction * ds_v;
        core_radius = (chord_core > width_core) ? chord_core : width_core;
      }
      double core_sq = core_radius * core_radius;

      // Velocity from left endpoint
      compute_trefftz_velocity_point(
          strip_gamma[jv], y_left_v, z_left_v, yc, zc, core_sq, &vy, &vz);

      // Velocity from right endpoint (negative circulation)
      compute_trefftz_velocity_point(
          -strip_gamma[jv], y_right_v, z_right_v, yc, zc, core_sq, &vy, &vz);
    }

    // Far-field lift (Kutta-Joukowski)
    *CL_ff += 2.0 * strip_gamma[jc] * dy_c / reference_area;

    // Far-field side force
    *CY_ff -= 2.0 * strip_gamma[jc] * dz_c / reference_area;

    // Far-field drag (kinetic energy of crossflow)
    *CDi_ff += strip_gamma[jc] * (dz_c * vy - dy_c * vz) / reference_area;
  }
}

/// @brief Compute span efficiency from Trefftz plane results.
///
/// @param CL_ff Far-field lift coefficient
/// @param CY_ff Far-field side force coefficient
/// @param CDi_ff Far-field induced drag coefficient
/// @param reference_area Reference area
/// @param reference_span Reference span
/// @return Span efficiency factor
inline double compute_span_efficiency(
    double CL_ff, double CY_ff, double CDi_ff,
    double reference_area, double reference_span) {
  double AR = reference_span * reference_span / reference_area;
  if (std::abs(CDi_ff) > kEpsTrefftz) {
    return (CL_ff * CL_ff + CY_ff * CY_ff) / (kPiTrefftz * AR * CDi_ff);
  }
  return 0.0;
}

// ============================================================================
// Enzyme Gradient Functions
// ============================================================================

#ifdef VLM_USE_ENZYME

/// @brief Compute gradient of Trefftz plane drag w.r.t. circulation using Enzyme.
///
/// @param strip_gamma Circulation per strip [num_strips]
/// @param strip_te_left Left trailing edge points [3*num_strips]
/// @param strip_te_right Right trailing edge points [3*num_strips]
/// @param strip_chord Strip chord lengths [num_strips]
/// @param strip_surface_index Surface index per strip [num_strips]
/// @param num_strips Number of strips
/// @param reference_area Reference area
/// @param reference_span Reference span
/// @param core_chord_fraction Core radius as fraction of chord
/// @param core_width_fraction Core radius as fraction of strip width
/// @param output_index Which output to differentiate (0=CL_ff, 1=CY_ff, 2=CDi_ff)
/// @param d_gamma Output gradient d(output)/d(gamma) [num_strips]
void enzyme_trefftz_circulation_gradient(
    const double* strip_gamma,
    const double* strip_te_left,
    const double* strip_te_right,
    const double* strip_chord,
    const int* strip_surface_index,
    int num_strips,
    double reference_area,
    double reference_span,
    double core_chord_fraction,
    double core_width_fraction,
    int output_index,
    double* d_gamma);

#endif  // VLM_USE_ENZYME

// ============================================================================
// Finite Difference Fallback
// ============================================================================

/// @brief Compute gradient of Trefftz drag w.r.t. circulation using FD.
void fd_trefftz_circulation_gradient(
    const double* strip_gamma,
    const double* strip_te_left,
    const double* strip_te_right,
    const double* strip_chord,
    const int* strip_surface_index,
    int num_strips,
    double reference_area,
    double reference_span,
    double core_chord_fraction,
    double core_width_fraction,
    int output_index,
    double* d_gamma,
    double h = 1e-7);

}  // namespace ad
}  // namespace vlm
