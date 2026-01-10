// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include "vlm/ad/enzyme_geometry.h"

#include <cmath>
#include <cstring>

namespace vlm {
namespace ad {

// ============================================================================
// Wrapper function for Enzyme differentiation
// ============================================================================

// This wrapper allows Enzyme to differentiate with respect to specific inputs
// by taking the inputs we want to vary as separate arguments.

namespace {

// Wrapper that takes chord as explicit parameter for differentiation
void panel_geometry_chord_wrapper(
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
  compute_panel_geometry(
      le_left_sec, le_right_sec,
      chord_left, chord_right,
      twist_left, twist_right,
      frac_left, frac_center, frac_right,
      x_vortex, x_control,
      bv_left_out, bv_right_out, cp_out, normal_out);
}

// Wrapper with LE position for differentiation
// For this we need a mutable copy of the LE arrays
void panel_geometry_le_wrapper(
    double* le_left_sec,  // mutable for Enzyme
    double* le_right_sec,  // mutable for Enzyme
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
  compute_panel_geometry(
      le_left_sec, le_right_sec,
      chord_left, chord_right,
      twist_left, twist_right,
      frac_left, frac_center, frac_right,
      x_vortex, x_control,
      bv_left_out, bv_right_out, cp_out, normal_out);
}

}  // namespace

// ============================================================================
// Enzyme Implementations
// ============================================================================

#ifdef VLM_USE_ENZYME

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
    double* d_normal) {
  // Use forward mode: seed the chord we're differentiating, get output tangents
  double bv_left[3], bv_right[3], cp[3], normal[3];

  // Zero shadow outputs
  std::memset(d_bv_left, 0, 3 * sizeof(double));
  std::memset(d_bv_right, 0, 3 * sizeof(double));
  std::memset(d_cp, 0, 3 * sizeof(double));
  std::memset(d_normal, 0, 3 * sizeof(double));

  // Forward mode differentiation
  double d_chord_left = (which_section == 0) ? 1.0 : 0.0;
  double d_chord_right = (which_section == 1) ? 1.0 : 0.0;

  __enzyme_fwddiff<void>(
      reinterpret_cast<void*>(panel_geometry_chord_wrapper),
      enzyme_const, le_left_sec,
      enzyme_const, le_right_sec,
      enzyme_dup, chord_left, d_chord_left,
      enzyme_dup, chord_right, d_chord_right,
      enzyme_const, twist_left,
      enzyme_const, twist_right,
      enzyme_const, frac_left,
      enzyme_const, frac_center,
      enzyme_const, frac_right,
      enzyme_const, x_vortex,
      enzyme_const, x_control,
      enzyme_dup, bv_left, d_bv_left,
      enzyme_dup, bv_right, d_bv_right,
      enzyme_dup, cp, d_cp,
      enzyme_dup, normal, d_normal);
}

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
    double* d_normal) {
  double bv_left[3], bv_right[3], cp[3], normal[3];

  std::memset(d_bv_left, 0, 3 * sizeof(double));
  std::memset(d_bv_right, 0, 3 * sizeof(double));
  std::memset(d_cp, 0, 3 * sizeof(double));
  std::memset(d_normal, 0, 3 * sizeof(double));

  double d_twist_left = (which_section == 0) ? 1.0 : 0.0;
  double d_twist_right = (which_section == 1) ? 1.0 : 0.0;

  __enzyme_fwddiff<void>(
      reinterpret_cast<void*>(panel_geometry_chord_wrapper),
      enzyme_const, le_left_sec,
      enzyme_const, le_right_sec,
      enzyme_const, chord_left,
      enzyme_const, chord_right,
      enzyme_dup, twist_left, d_twist_left,
      enzyme_dup, twist_right, d_twist_right,
      enzyme_const, frac_left,
      enzyme_const, frac_center,
      enzyme_const, frac_right,
      enzyme_const, x_vortex,
      enzyme_const, x_control,
      enzyme_dup, bv_left, d_bv_left,
      enzyme_dup, bv_right, d_bv_right,
      enzyme_dup, cp, d_cp,
      enzyme_dup, normal, d_normal);
}

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
    double* d_normal) {
  // Make mutable copies of LE arrays
  double le_left_mut[3] = {le_left_sec[0], le_left_sec[1], le_left_sec[2]};
  double le_right_mut[3] = {le_right_sec[0], le_right_sec[1], le_right_sec[2]};

  double bv_left[3], bv_right[3], cp[3], normal[3];

  std::memset(d_bv_left, 0, 3 * sizeof(double));
  std::memset(d_bv_right, 0, 3 * sizeof(double));
  std::memset(d_cp, 0, 3 * sizeof(double));
  std::memset(d_normal, 0, 3 * sizeof(double));

  // Shadow inputs for LE
  double d_le_left[3] = {0.0, 0.0, 0.0};
  double d_le_right[3] = {0.0, 0.0, 0.0};

  if (which_section == 0) {
    d_le_left[component] = 1.0;
  } else {
    d_le_right[component] = 1.0;
  }

  __enzyme_fwddiff<void>(
      reinterpret_cast<void*>(panel_geometry_le_wrapper),
      enzyme_dup, le_left_mut, d_le_left,
      enzyme_dup, le_right_mut, d_le_right,
      enzyme_const, chord_left,
      enzyme_const, chord_right,
      enzyme_const, twist_left,
      enzyme_const, twist_right,
      enzyme_const, frac_left,
      enzyme_const, frac_center,
      enzyme_const, frac_right,
      enzyme_const, x_vortex,
      enzyme_const, x_control,
      enzyme_dup, bv_left, d_bv_left,
      enzyme_dup, bv_right, d_bv_right,
      enzyme_dup, cp, d_cp,
      enzyme_dup, normal, d_normal);
}

#endif  // VLM_USE_ENZYME

// ============================================================================
// Finite Difference Implementations
// ============================================================================

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
    double h) {
  double bv_left_plus[3], bv_right_plus[3], cp_plus[3], normal_plus[3];
  double bv_left_minus[3], bv_right_minus[3], cp_minus[3], normal_minus[3];

  double c_left_p = chord_left + (which_section == 0 ? h : 0.0);
  double c_right_p = chord_right + (which_section == 1 ? h : 0.0);
  double c_left_m = chord_left - (which_section == 0 ? h : 0.0);
  double c_right_m = chord_right - (which_section == 1 ? h : 0.0);

  compute_panel_geometry(
      le_left_sec, le_right_sec, c_left_p, c_right_p,
      twist_left, twist_right,
      frac_left, frac_center, frac_right, x_vortex, x_control,
      bv_left_plus, bv_right_plus, cp_plus, normal_plus);

  compute_panel_geometry(
      le_left_sec, le_right_sec, c_left_m, c_right_m,
      twist_left, twist_right,
      frac_left, frac_center, frac_right, x_vortex, x_control,
      bv_left_minus, bv_right_minus, cp_minus, normal_minus);

  for (int i = 0; i < 3; ++i) {
    d_bv_left[i] = (bv_left_plus[i] - bv_left_minus[i]) / (2.0 * h);
    d_bv_right[i] = (bv_right_plus[i] - bv_right_minus[i]) / (2.0 * h);
    d_cp[i] = (cp_plus[i] - cp_minus[i]) / (2.0 * h);
    d_normal[i] = (normal_plus[i] - normal_minus[i]) / (2.0 * h);
  }
}

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
    double h) {
  double bv_left_plus[3], bv_right_plus[3], cp_plus[3], normal_plus[3];
  double bv_left_minus[3], bv_right_minus[3], cp_minus[3], normal_minus[3];

  double t_left_p = twist_left + (which_section == 0 ? h : 0.0);
  double t_right_p = twist_right + (which_section == 1 ? h : 0.0);
  double t_left_m = twist_left - (which_section == 0 ? h : 0.0);
  double t_right_m = twist_right - (which_section == 1 ? h : 0.0);

  compute_panel_geometry(
      le_left_sec, le_right_sec, chord_left, chord_right,
      t_left_p, t_right_p,
      frac_left, frac_center, frac_right, x_vortex, x_control,
      bv_left_plus, bv_right_plus, cp_plus, normal_plus);

  compute_panel_geometry(
      le_left_sec, le_right_sec, chord_left, chord_right,
      t_left_m, t_right_m,
      frac_left, frac_center, frac_right, x_vortex, x_control,
      bv_left_minus, bv_right_minus, cp_minus, normal_minus);

  for (int i = 0; i < 3; ++i) {
    d_bv_left[i] = (bv_left_plus[i] - bv_left_minus[i]) / (2.0 * h);
    d_bv_right[i] = (bv_right_plus[i] - bv_right_minus[i]) / (2.0 * h);
    d_cp[i] = (cp_plus[i] - cp_minus[i]) / (2.0 * h);
    d_normal[i] = (normal_plus[i] - normal_minus[i]) / (2.0 * h);
  }
}

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
    double h) {
  double bv_left_plus[3], bv_right_plus[3], cp_plus[3], normal_plus[3];
  double bv_left_minus[3], bv_right_minus[3], cp_minus[3], normal_minus[3];

  double le_left_p[3] = {le_left_sec[0], le_left_sec[1], le_left_sec[2]};
  double le_right_p[3] = {le_right_sec[0], le_right_sec[1], le_right_sec[2]};
  double le_left_m[3] = {le_left_sec[0], le_left_sec[1], le_left_sec[2]};
  double le_right_m[3] = {le_right_sec[0], le_right_sec[1], le_right_sec[2]};

  if (which_section == 0) {
    le_left_p[component] += h;
    le_left_m[component] -= h;
  } else {
    le_right_p[component] += h;
    le_right_m[component] -= h;
  }

  compute_panel_geometry(
      le_left_p, le_right_p, chord_left, chord_right,
      twist_left, twist_right,
      frac_left, frac_center, frac_right, x_vortex, x_control,
      bv_left_plus, bv_right_plus, cp_plus, normal_plus);

  compute_panel_geometry(
      le_left_m, le_right_m, chord_left, chord_right,
      twist_left, twist_right,
      frac_left, frac_center, frac_right, x_vortex, x_control,
      bv_left_minus, bv_right_minus, cp_minus, normal_minus);

  for (int i = 0; i < 3; ++i) {
    d_bv_left[i] = (bv_left_plus[i] - bv_left_minus[i]) / (2.0 * h);
    d_bv_right[i] = (bv_right_plus[i] - bv_right_minus[i]) / (2.0 * h);
    d_cp[i] = (cp_plus[i] - cp_minus[i]) / (2.0 * h);
    d_normal[i] = (normal_plus[i] - normal_minus[i]) / (2.0 * h);
  }
}

}  // namespace ad
}  // namespace vlm
