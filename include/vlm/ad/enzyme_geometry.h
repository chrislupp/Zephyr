// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#pragma once

/// @file enzyme_geometry.h
/// @brief Pure kernel functions for geometry computations that Enzyme can differentiate.
///
/// These functions compute mesh quantities (positions, normals) from section
/// geometry parameters. They are designed as pure functions with raw pointer
/// interfaces suitable for Enzyme automatic differentiation.
///
/// The key mesh quantities are:
/// - Control point position: where flow tangency is enforced
/// - Bound vortex endpoints: where the horseshoe vortex is located
/// - Panel normal: direction of flow tangency constraint
///
/// These depend on section parameters:
/// - Leading edge position (x, y, z)
/// - Chord length
/// - Twist angle
/// - Dihedral angle (implicit in LE positions)

#include <cmath>

#include "vlm/ad/enzyme_interface.h"

// ============================================================================
// Custom atan2 for Enzyme compatibility
// ============================================================================

namespace vlm {
namespace ad {
namespace detail {

/// Custom atan2 that Enzyme can differentiate.
/// We use this instead of std::atan2 which may compile to llvm.atan2 intrinsic.
inline double enzyme_safe_atan2(double y, double x) {
  // Use explicit computation instead of intrinsic
  if (x > 0.0) {
    return std::atan(y / x);
  } else if (x < 0.0 && y >= 0.0) {
    return std::atan(y / x) + 3.14159265358979323846;
  } else if (x < 0.0 && y < 0.0) {
    return std::atan(y / x) - 3.14159265358979323846;
  } else if (x == 0.0 && y > 0.0) {
    return 3.14159265358979323846 / 2.0;
  } else if (x == 0.0 && y < 0.0) {
    return -3.14159265358979323846 / 2.0;
  } else {
    return 0.0;  // x == 0 && y == 0
  }
}

}  // namespace detail
}  // namespace ad
}  // namespace vlm

namespace vlm {
namespace ad {

// ============================================================================
// Pure Kernel Functions for Mesh Generation
// ============================================================================

/// @brief Compute interpolated leading edge position for a strip.
///
/// Interpolates between two section LE positions based on spanwise fraction.
///
/// @param le_left Leading edge position of left section [3]
/// @param le_right Leading edge position of right section [3]
/// @param frac Interpolation fraction (0 = left, 1 = right)
/// @param le_out Output interpolated LE position [3]
inline void compute_interpolated_le(
    const double* le_left,
    const double* le_right,
    double frac,
    double* le_out) {
  for (int i = 0; i < 3; ++i) {
    le_out[i] = (1.0 - frac) * le_left[i] + frac * le_right[i];
  }
}

/// @brief Compute interpolated chord for a strip.
///
/// @param chord_left Chord at left section
/// @param chord_right Chord at right section
/// @param frac Interpolation fraction
/// @return Interpolated chord
inline double compute_interpolated_chord(
    double chord_left,
    double chord_right,
    double frac) {
  return (1.0 - frac) * chord_left + frac * chord_right;
}

/// @brief Compute interpolated twist angle using chord-weighted sin/cos.
///
/// This matches the VLM implementation which interpolates sin and cos
/// of twist weighted by chord, then recovers the angle.
///
/// @param twist_left Twist angle at left section (radians)
/// @param twist_right Twist angle at right section (radians)
/// @param chord_left Chord at left section
/// @param chord_right Chord at right section
/// @param frac Interpolation fraction
/// @return Interpolated twist angle (radians)
inline double compute_interpolated_twist(
    double twist_left,
    double twist_right,
    double chord_left,
    double chord_right,
    double frac) {
  // Chord-weighted sin/cos interpolation
  double sin_left = chord_left * std::sin(twist_left);
  double cos_left = chord_left * std::cos(twist_left);
  double sin_right = chord_right * std::sin(twist_right);
  double cos_right = chord_right * std::cos(twist_right);

  double sin_interp = (1.0 - frac) * sin_left + frac * sin_right;
  double cos_interp = (1.0 - frac) * cos_left + frac * cos_right;

  return detail::enzyme_safe_atan2(sin_interp, cos_interp);
}

/// @brief Compute Trefftz plane normal components (Y-Z plane normal).
///
/// The Trefftz normal is perpendicular to the strip span direction in the
/// Y-Z plane. For a wing lying in the XY plane, this is (0, 0, 1).
///
/// @param le_left Left edge leading edge position [3]
/// @param le_right Right edge leading edge position [3]
/// @param normal_y Output: Y-component of Trefftz normal
/// @param normal_z Output: Z-component of Trefftz normal
inline void compute_trefftz_normal(
    const double* le_left,
    const double* le_right,
    double* normal_y,
    double* normal_z) {
  double dy = le_right[1] - le_left[1];
  double dz = le_right[2] - le_left[2];
  double yz_mag = std::sqrt(dy * dy + dz * dz);

  if (yz_mag > 1e-12) {
    *normal_y = -dz / yz_mag;
    *normal_z = dy / yz_mag;
  } else {
    *normal_y = 0.0;
    *normal_z = 1.0;
  }
}

/// @brief Compute control point position from strip geometry.
///
/// Control point is at the 3/4 chord position of the panel.
///
/// @param le Leading edge position [3]
/// @param chord Chord length
/// @param x_frac Chordwise fraction for control point (typically 0.75)
/// @param cp_out Output control point position [3]
inline void compute_control_point(
    const double* le,
    double chord,
    double x_frac,
    double* cp_out) {
  cp_out[0] = le[0] + x_frac * chord;
  cp_out[1] = le[1];
  cp_out[2] = le[2];
}

/// @brief Compute bound vortex endpoint from strip edge geometry.
///
/// @param le Leading edge position [3]
/// @param chord Chord length
/// @param x_frac Chordwise fraction for vortex (typically 0.25)
/// @param bv_out Output bound vortex endpoint [3]
inline void compute_bound_vortex_point(
    const double* le,
    double chord,
    double x_frac,
    double* bv_out) {
  bv_out[0] = le[0] + x_frac * chord;
  bv_out[1] = le[1];
  bv_out[2] = le[2];
}

/// @brief Compute panel normal vector from twist and Trefftz normal.
///
/// The chord direction is computed using the twist angle and Trefftz normal:
///   chord_dir = (cos(twist), -sin(twist)*trefftz_y, -sin(twist)*trefftz_z)
///
/// The normal is then chord_dir × bound_dir, normalized.
///
/// @param twist Twist angle (radians)
/// @param trefftz_y Y-component of Trefftz plane normal
/// @param trefftz_z Z-component of Trefftz plane normal
/// @param bound_dir Bound vortex direction (unit vector) [3]
/// @param normal_out Output normal vector [3]
inline void compute_panel_normal(
    double twist,
    double trefftz_y,
    double trefftz_z,
    const double* bound_dir,
    double* normal_out) {
  // Chord direction using Trefftz normal (matches OptVL)
  double cos_t = std::cos(twist);
  double sin_t = std::sin(twist);
  double chord_x = cos_t;
  double chord_y = -sin_t * trefftz_y;
  double chord_z = -sin_t * trefftz_z;

  // Normal = chord_dir × bound_dir
  double nx = chord_y * bound_dir[2] - chord_z * bound_dir[1];
  double ny = chord_z * bound_dir[0] - chord_x * bound_dir[2];
  double nz = chord_x * bound_dir[1] - chord_y * bound_dir[0];

  // Normalize
  double mag = std::sqrt(nx * nx + ny * ny + nz * nz);
  if (mag > 1e-12) {
    normal_out[0] = nx / mag;
    normal_out[1] = ny / mag;
    normal_out[2] = nz / mag;
  } else {
    // Fallback
    normal_out[0] = 0.0;
    normal_out[1] = trefftz_y;
    normal_out[2] = trefftz_z;
  }
}

/// @brief Compute bound vortex direction (unit vector).
///
/// @param bv_left Left bound vortex endpoint [3]
/// @param bv_right Right bound vortex endpoint [3]
/// @param dir_out Output direction unit vector [3]
inline void compute_bound_direction(
    const double* bv_left,
    const double* bv_right,
    double* dir_out) {
  double dx = bv_right[0] - bv_left[0];
  double dy = bv_right[1] - bv_left[1];
  double dz = bv_right[2] - bv_left[2];
  double mag = std::sqrt(dx * dx + dy * dy + dz * dz);

  if (mag > 1e-12) {
    dir_out[0] = dx / mag;
    dir_out[1] = dy / mag;
    dir_out[2] = dz / mag;
  } else {
    dir_out[0] = 0.0;
    dir_out[1] = 1.0;
    dir_out[2] = 0.0;
  }
}

// ============================================================================
// Complete Panel Geometry Computation
// ============================================================================

/// @brief Compute all panel geometry from section parameters.
///
/// This is the main kernel for computing panel positions and normal from
/// section geometry. It computes all quantities needed for the VLM solve.
///
/// @param le_left_sec Left section leading edge [3]
/// @param le_right_sec Right section leading edge [3]
/// @param chord_left Left section chord
/// @param chord_right Right section chord
/// @param twist_left Left section twist (radians)
/// @param twist_right Right section twist (radians)
/// @param frac_left Interpolation fraction for left panel edge
/// @param frac_center Interpolation fraction for panel center
/// @param frac_right Interpolation fraction for right panel edge
/// @param x_vortex Chordwise fraction for bound vortex (typically 0.25)
/// @param x_control Chordwise fraction for control point (typically 0.75)
/// @param bv_left_out Output left bound vortex point [3]
/// @param bv_right_out Output right bound vortex point [3]
/// @param cp_out Output control point [3]
/// @param normal_out Output panel normal [3]
inline void compute_panel_geometry(
    const double* le_left_sec,
    const double* le_right_sec,
    double chord_left,
    double chord_right,
    double twist_left,
    double twist_right,
    double frac_left,
    double frac_center,
    double frac_right,
    double x_vortex,
    double x_control,
    double* bv_left_out,
    double* bv_right_out,
    double* cp_out,
    double* normal_out) {
  // Interpolate LE positions
  double le_left[3], le_center[3], le_right[3];
  compute_interpolated_le(le_left_sec, le_right_sec, frac_left, le_left);
  compute_interpolated_le(le_left_sec, le_right_sec, frac_center, le_center);
  compute_interpolated_le(le_left_sec, le_right_sec, frac_right, le_right);

  // Interpolate chords
  double c_left = compute_interpolated_chord(chord_left, chord_right, frac_left);
  double c_center = compute_interpolated_chord(chord_left, chord_right, frac_center);
  double c_right = compute_interpolated_chord(chord_left, chord_right, frac_right);

  // Interpolate twist
  double twist = compute_interpolated_twist(
      twist_left, twist_right, chord_left, chord_right, frac_center);

  // Compute Trefftz normal from panel edge positions
  double trefftz_y, trefftz_z;
  compute_trefftz_normal(le_left, le_right, &trefftz_y, &trefftz_z);

  // Compute bound vortex endpoints
  compute_bound_vortex_point(le_left, c_left, x_vortex, bv_left_out);
  compute_bound_vortex_point(le_right, c_right, x_vortex, bv_right_out);

  // Compute control point
  compute_control_point(le_center, c_center, x_control, cp_out);

  // Compute bound direction
  double bound_dir[3];
  compute_bound_direction(bv_left_out, bv_right_out, bound_dir);

  // Compute normal
  compute_panel_normal(twist, trefftz_y, trefftz_z, bound_dir, normal_out);
}

// ============================================================================
// Enzyme-Differentiated Geometry Gradients
// ============================================================================

#ifdef VLM_USE_ENZYME

/// @brief Compute gradients of panel geometry w.r.t. section chord.
///
/// Uses Enzyme forward mode to compute:
///   d(bv_left)/d(chord), d(bv_right)/d(chord), d(cp)/d(chord), d(normal)/d(chord)
///
/// @param le_left_sec Left section leading edge [3]
/// @param le_right_sec Right section leading edge [3]
/// @param chord_left Left section chord
/// @param chord_right Right section chord
/// @param twist_left Left section twist (radians)
/// @param twist_right Right section twist (radians)
/// @param frac_left Interpolation fraction for left edge
/// @param frac_center Interpolation fraction for center
/// @param frac_right Interpolation fraction for right edge
/// @param x_vortex Chordwise fraction for bound vortex
/// @param x_control Chordwise fraction for control point
/// @param which_section 0 = left section, 1 = right section
/// @param d_bv_left Output: d(bv_left)/d(chord) [3]
/// @param d_bv_right Output: d(bv_right)/d(chord) [3]
/// @param d_cp Output: d(cp)/d(chord) [3]
/// @param d_normal Output: d(normal)/d(chord) [3]
void enzyme_panel_geometry_chord_gradient(
    const double* le_left_sec,
    const double* le_right_sec,
    double chord_left,
    double chord_right,
    double twist_left,
    double twist_right,
    double frac_left,
    double frac_center,
    double frac_right,
    double x_vortex,
    double x_control,
    int which_section,
    double* d_bv_left,
    double* d_bv_right,
    double* d_cp,
    double* d_normal);

/// @brief Compute gradients of panel geometry w.r.t. section twist.
void enzyme_panel_geometry_twist_gradient(
    const double* le_left_sec,
    const double* le_right_sec,
    double chord_left,
    double chord_right,
    double twist_left,
    double twist_right,
    double frac_left,
    double frac_center,
    double frac_right,
    double x_vortex,
    double x_control,
    int which_section,
    double* d_bv_left,
    double* d_bv_right,
    double* d_cp,
    double* d_normal);

/// @brief Compute gradients of panel geometry w.r.t. section LE position.
///
/// @param component 0=x, 1=y, 2=z
void enzyme_panel_geometry_le_gradient(
    const double* le_left_sec,
    const double* le_right_sec,
    double chord_left,
    double chord_right,
    double twist_left,
    double twist_right,
    double frac_left,
    double frac_center,
    double frac_right,
    double x_vortex,
    double x_control,
    int which_section,
    int component,
    double* d_bv_left,
    double* d_bv_right,
    double* d_cp,
    double* d_normal);

#endif  // VLM_USE_ENZYME

// ============================================================================
// Finite Difference Fallback
// ============================================================================

/// @brief Compute gradients of panel geometry using finite differences.
///
/// This is the fallback when Enzyme is not available.
void fd_panel_geometry_chord_gradient(
    const double* le_left_sec,
    const double* le_right_sec,
    double chord_left,
    double chord_right,
    double twist_left,
    double twist_right,
    double frac_left,
    double frac_center,
    double frac_right,
    double x_vortex,
    double x_control,
    int which_section,
    double* d_bv_left,
    double* d_bv_right,
    double* d_cp,
    double* d_normal,
    double h = 1e-7);

void fd_panel_geometry_twist_gradient(
    const double* le_left_sec,
    const double* le_right_sec,
    double chord_left,
    double chord_right,
    double twist_left,
    double twist_right,
    double frac_left,
    double frac_center,
    double frac_right,
    double x_vortex,
    double x_control,
    int which_section,
    double* d_bv_left,
    double* d_bv_right,
    double* d_cp,
    double* d_normal,
    double h = 1e-6);

void fd_panel_geometry_le_gradient(
    const double* le_left_sec,
    const double* le_right_sec,
    double chord_left,
    double chord_right,
    double twist_left,
    double twist_right,
    double frac_left,
    double frac_center,
    double frac_right,
    double x_vortex,
    double x_control,
    int which_section,
    int component,
    double* d_bv_left,
    double* d_bv_right,
    double* d_cp,
    double* d_normal,
    double h = 1e-7);

}  // namespace ad
}  // namespace vlm
