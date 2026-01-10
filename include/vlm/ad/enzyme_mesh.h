// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#pragma once

/// @file enzyme_mesh.h
/// @brief Pure kernel functions for mesh geometry computations that Enzyme can differentiate.
///
/// These functions compute panel positions from section parameters.
/// They are designed as pure functions with raw pointer interfaces suitable
/// for Enzyme automatic differentiation.
///
/// The mesh geometry depends on section parameters through linear interpolation:
/// - Panel chord at spanwise fraction eta: c = (1-eta)*c_left + eta*c_right
/// - Panel position: p = (1-eta)*p_left + eta*p_right + xi*chord*direction
///
/// Key parameters:
/// - eta: Spanwise fraction within section interval [0, 1]
/// - xi: Chordwise fraction [0, 1]
/// - Twist angles: Rotate normal vector about spanwise axis

#include <cmath>

#include "vlm/ad/enzyme_interface.h"

namespace vlm {
namespace ad {

// ============================================================================
// Panel Geometry Kernels
// ============================================================================

/// @brief Interpolate scalar value between two sections.
///
/// @param value_left Value at left section
/// @param value_right Value at right section
/// @param eta Spanwise interpolation fraction [0, 1]
/// @return Interpolated value
inline double interpolate_scalar(double value_left, double value_right, double eta) {
  return (1.0 - eta) * value_left + eta * value_right;
}

/// @brief Compute panel leading edge position.
///
/// @param le_left Left section leading edge [3]
/// @param le_right Right section leading edge [3]
/// @param eta Spanwise interpolation fraction [0, 1]
/// @param le_out Output leading edge position [3]
inline void compute_panel_leading_edge(
    const double* le_left,
    const double* le_right,
    double eta,
    double* le_out) {
  for (int i = 0; i < 3; ++i) {
    le_out[i] = (1.0 - eta) * le_left[i] + eta * le_right[i];
  }
}

/// @brief Compute bound vortex endpoint position.
///
/// The bound vortex is placed at x_vr fraction along the chord from the leading edge.
///
/// @param le Leading edge position [3]
/// @param chord Local chord length
/// @param x_vr Chordwise fraction for vortex position [0, 1]
/// @param bound_out Output bound vortex position [3]
inline void compute_bound_vortex_position(
    const double* le,
    double chord,
    double x_vr,
    double* bound_out) {
  bound_out[0] = le[0] + x_vr * chord;
  bound_out[1] = le[1];
  bound_out[2] = le[2];
}

/// @brief Compute control point position.
///
/// The control point is placed at x_cp fraction along the chord from the leading edge.
///
/// @param le Leading edge position [3]
/// @param chord Local chord length
/// @param x_cp Chordwise fraction for control point [0, 1]
/// @param cp_out Output control point position [3]
inline void compute_control_point_position(
    const double* le,
    double chord,
    double x_cp,
    double* cp_out) {
  cp_out[0] = le[0] + x_cp * chord;
  cp_out[1] = le[1];
  cp_out[2] = le[2];
}

/// @brief Compute panel normal vector from twist and bound vortex direction.
///
/// The normal is computed as: n = (chord_direction) x (bound_direction)
/// where chord_direction is rotated by the twist angle.
///
/// @param twist_angle Twist angle (radians)
/// @param trefftz_normal_y Trefftz plane normal y component
/// @param trefftz_normal_z Trefftz plane normal z component
/// @param bound_left Left bound vortex endpoint [3]
/// @param bound_right Right bound vortex endpoint [3]
/// @param camber_slope Camber slope at point (dz/dx)
/// @param normal_out Output normal vector [3]
inline void compute_panel_normal(
    double twist_angle,
    double trefftz_normal_y,
    double trefftz_normal_z,
    const double* bound_left,
    const double* bound_right,
    double camber_slope,
    double* normal_out) {
  // Bound vortex direction
  double dx = bound_right[0] - bound_left[0];
  double dy = bound_right[1] - bound_left[1];
  double dz = bound_right[2] - bound_left[2];
  double len = std::sqrt(dx * dx + dy * dy + dz * dz);

  double bound_dir[3] = {dx / len, dy / len, dz / len};

  // Effective angle = twist - atan(camber_slope)
  double ang = twist_angle - std::atan(camber_slope);

  // Chord direction (rotated about Trefftz normal)
  double chord_dir[3];
  chord_dir[0] = std::cos(ang);
  chord_dir[1] = -std::sin(ang) * trefftz_normal_y;
  chord_dir[2] = -std::sin(ang) * trefftz_normal_z;

  // Normal = chord_dir x bound_dir
  double nx = chord_dir[1] * bound_dir[2] - chord_dir[2] * bound_dir[1];
  double ny = chord_dir[2] * bound_dir[0] - chord_dir[0] * bound_dir[2];
  double nz = chord_dir[0] * bound_dir[1] - chord_dir[1] * bound_dir[0];

  double nmag = std::sqrt(nx * nx + ny * ny + nz * nz);
  if (nmag > 1e-12) {
    normal_out[0] = nx / nmag;
    normal_out[1] = ny / nmag;
    normal_out[2] = nz / nmag;
  } else {
    // Fallback to Trefftz normal
    normal_out[0] = 0.0;
    normal_out[1] = trefftz_normal_y;
    normal_out[2] = trefftz_normal_z;
  }
}

/// @brief Compute Trefftz plane normal from leading edge positions.
///
/// The Trefftz normal is perpendicular to the span direction in the y-z plane.
///
/// @param le_left Left leading edge [3]
/// @param le_right Right leading edge [3]
/// @param trefftz_ny Output Trefftz normal y component
/// @param trefftz_nz Output Trefftz normal z component
inline void compute_trefftz_normal(
    const double* le_left,
    const double* le_right,
    double* trefftz_ny,
    double* trefftz_nz) {
  double dy = le_right[1] - le_left[1];
  double dz = le_right[2] - le_left[2];
  double yz_mag = std::sqrt(dy * dy + dz * dz);

  if (yz_mag > 1e-12) {
    *trefftz_ny = -dz / yz_mag;
    *trefftz_nz = dy / yz_mag;
  } else {
    *trefftz_ny = 0.0;
    *trefftz_nz = 1.0;
  }
}

/// @brief Compute chord-weighted twist interpolation.
///
/// Twist is interpolated using chord-weighted sin/cos to avoid discontinuities.
///
/// @param chord_left Left section chord
/// @param chord_right Right section chord
/// @param twist_left Left section twist (radians)
/// @param twist_right Right section twist (radians)
/// @param eta Spanwise interpolation fraction
/// @return Interpolated twist angle (radians)
inline double interpolate_twist(
    double chord_left, double chord_right,
    double twist_left, double twist_right,
    double eta) {
  double sin_left = chord_left * std::sin(twist_left);
  double cos_left = chord_left * std::cos(twist_left);
  double sin_right = chord_right * std::sin(twist_right);
  double cos_right = chord_right * std::cos(twist_right);

  double sin_interp = (1.0 - eta) * sin_left + eta * sin_right;
  double cos_interp = (1.0 - eta) * cos_left + eta * cos_right;

  return std::atan2(sin_interp, cos_interp);
}

// ============================================================================
// Combined Panel Geometry Kernel
// ============================================================================

/// @brief Compute complete panel geometry from section parameters.
///
/// This kernel computes all geometry needed for AIC matrix computation:
/// - Control point position
/// - Bound vortex endpoints (left and right)
/// - Normal vector at control point
///
/// @param chord_left Left section chord
/// @param chord_right Right section chord
/// @param twist_left Left section twist (radians)
/// @param twist_right Right section twist (radians)
/// @param le_left Left section leading edge [3]
/// @param le_right Right section leading edge [3]
/// @param eta Spanwise interpolation fraction [0, 1]
/// @param x_vr Chordwise fraction for bound vortex
/// @param x_cp Chordwise fraction for control point
/// @param camber_slope Camber slope at control point
/// @param control_point Output control point [3]
/// @param bound_left Output left bound vortex endpoint [3]
/// @param bound_right Output right bound vortex endpoint [3]
/// @param normal Output normal vector [3]
inline void compute_panel_geometry_kernel(
    double chord_left,
    double chord_right,
    double twist_left,
    double twist_right,
    const double* le_left,
    const double* le_right,
    double eta,
    double x_vr,
    double x_cp,
    double camber_slope,
    double* control_point,
    double* bound_left,
    double* bound_right,
    double* normal) {
  // Interpolate chord
  double chord_center = interpolate_scalar(chord_left, chord_right, eta);

  // Interpolate leading edge positions
  double le_center[3];
  compute_panel_leading_edge(le_left, le_right, eta, le_center);

  // For bound vortex endpoints, we need LE at panel edges (eta_left, eta_right)
  // For simplicity, assume eta_left = eta - delta, eta_right = eta + delta
  // In practice, the caller provides the appropriate eta values

  // Compute bound vortex left endpoint (at left panel edge)
  // Use left section values directly (eta = 0 for left edge of first panel)
  compute_bound_vortex_position(le_left, chord_left, x_vr, bound_left);

  // Compute bound vortex right endpoint (at right panel edge)
  compute_bound_vortex_position(le_right, chord_right, x_vr, bound_right);

  // Compute control point (at panel center)
  compute_control_point_position(le_center, chord_center, x_cp, control_point);

  // Compute Trefftz normal
  double trefftz_ny, trefftz_nz;
  compute_trefftz_normal(le_left, le_right, &trefftz_ny, &trefftz_nz);

  // Interpolate twist angle
  double twist_center = interpolate_twist(chord_left, chord_right,
                                          twist_left, twist_right, eta);

  // Compute normal
  compute_panel_normal(twist_center, trefftz_ny, trefftz_nz,
                       bound_left, bound_right, camber_slope, normal);
}

// ============================================================================
// Gradient Functions (Analytical)
// ============================================================================

/// @brief Compute gradient of panel chord w.r.t. section chords.
///
/// Since chord_center = (1-eta)*chord_left + eta*chord_right:
/// - d(chord_center)/d(chord_left) = 1 - eta
/// - d(chord_center)/d(chord_right) = eta
///
/// @param eta Spanwise interpolation fraction
/// @param d_chord_left Output gradient w.r.t. left chord
/// @param d_chord_right Output gradient w.r.t. right chord
inline void compute_chord_gradient(
    double eta,
    double* d_chord_left,
    double* d_chord_right) {
  *d_chord_left = 1.0 - eta;
  *d_chord_right = eta;
}

/// @brief Compute gradient of panel leading edge w.r.t. section leading edges.
///
/// Since le_center = (1-eta)*le_left + eta*le_right:
/// - d(le_center)/d(le_left) = (1 - eta) * I (identity)
/// - d(le_center)/d(le_right) = eta * I
///
/// @param eta Spanwise interpolation fraction
/// @param d_le_left Output gradient w.r.t. left LE (scalar multiplier for identity)
/// @param d_le_right Output gradient w.r.t. right LE (scalar multiplier for identity)
inline void compute_leading_edge_gradient(
    double eta,
    double* d_le_left,
    double* d_le_right) {
  *d_le_left = 1.0 - eta;
  *d_le_right = eta;
}

/// @brief Compute gradient of control point x-coordinate w.r.t. chord.
///
/// cp.x = le.x + x_cp * chord
/// d(cp.x)/d(chord) = x_cp
///
/// @param x_cp Chordwise fraction for control point
/// @return Gradient of cp.x w.r.t. chord
inline double compute_control_point_chord_gradient(double x_cp) {
  return x_cp;
}

/// @brief Compute gradient of bound vortex x-coordinate w.r.t. chord.
///
/// bound.x = le.x + x_vr * chord
/// d(bound.x)/d(chord) = x_vr
///
/// @param x_vr Chordwise fraction for bound vortex
/// @return Gradient of bound.x w.r.t. chord
inline double compute_bound_vortex_chord_gradient(double x_vr) {
  return x_vr;
}

}  // namespace ad
}  // namespace vlm
