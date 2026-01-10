// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include "vlm/ad/enzyme_aic.h"

namespace vlm {
namespace ad {

// ============================================================================
// Enzyme Gradient Implementations
// ============================================================================

#ifdef VLM_USE_ENZYME

// Wrapper for single AIC element (scalar output for Enzyme)
static double aic_element_wrapper(
    const double* rv1,
    const double* rv2,
    const double* rc,
    const double* normal,
    double core_radius,
    double beta,
    int include_bound) {
  return compute_aic_element(rv1, rv2, rc, normal, core_radius, beta,
                              include_bound != 0);
}

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
    double* d_normal) {

  // Initialize gradients
  for (int i = 0; i < 3; ++i) {
    d_rv1[i] = 0.0;
    d_rv2[i] = 0.0;
    d_rc[i] = 0.0;
    d_normal[i] = 0.0;
  }

  // Shadow arrays for Enzyme
  double d_rv1_local[3] = {0.0, 0.0, 0.0};
  double d_rv2_local[3] = {0.0, 0.0, 0.0};
  double d_rc_local[3] = {0.0, 0.0, 0.0};
  double d_normal_local[3] = {0.0, 0.0, 0.0};

  int ib = include_bound ? 1 : 0;

  // Reverse mode autodiff
  __enzyme_autodiff<void>(
      (void*)aic_element_wrapper,
      enzyme_dup, rv1, d_rv1_local,
      enzyme_dup, rv2, d_rv2_local,
      enzyme_dup, rc, d_rc_local,
      enzyme_dup, normal, d_normal_local,
      enzyme_const, core_radius,
      enzyme_const, beta,
      enzyme_const, ib);

  // Copy results
  for (int i = 0; i < 3; ++i) {
    d_rv1[i] = d_rv1_local[i];
    d_rv2[i] = d_rv2_local[i];
    d_rc[i] = d_rc_local[i];
    d_normal[i] = d_normal_local[i];
  }
}

// Wrapper for AIC row with scalar output (sum of weighted elements)
static double aic_row_dot_wrapper(
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
    const double* weights) {

  double result = 0.0;
  for (int j = 0; j < num_vortices; ++j) {
    const double* rv1 = &rv1_all[3*j];
    const double* rv2 = &rv2_all[3*j];

    // Compute strip width
    double dy = rv2[1] - rv1[1];
    double dz = rv2[2] - rv1[2];
    double strip_width = std::sqrt(dy*dy + dz*dz);

    // Core radius
    double core_radius = 0.0;
    if (surface_idx_all[j] != control_surface_idx) {
      double chord_core = core_chord_fraction * chord_all[j];
      double width_core = core_width_fraction * strip_width;
      core_radius = (chord_core > width_core) ? chord_core : width_core;
    }

    bool include_bound = (surface_idx_all[j] != control_surface_idx);
    double aic_elem = compute_aic_element(rv1, rv2, rc, normal, core_radius,
                                           beta, include_bound);
    result += weights[j] * aic_elem;
  }
  return result;
}

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
    double* d_normal) {

  // Initialize gradients
  for (int i = 0; i < 3; ++i) {
    d_rc[i] = 0.0;
    d_normal[i] = 0.0;
  }

  double d_rc_local[3] = {0.0, 0.0, 0.0};
  double d_normal_local[3] = {0.0, 0.0, 0.0};

  // Use the weighted sum wrapper - the weights are the seed gradients
  __enzyme_autodiff<void>(
      (void*)aic_row_dot_wrapper,
      enzyme_const, rv1_all,
      enzyme_const, rv2_all,
      enzyme_dup, rc, d_rc_local,
      enzyme_dup, normal, d_normal_local,
      enzyme_const, chord_all,
      enzyme_const, surface_idx_all,
      enzyme_const, num_vortices,
      enzyme_const, control_surface_idx,
      enzyme_const, core_chord_fraction,
      enzyme_const, core_width_fraction,
      enzyme_const, beta,
      enzyme_const, d_aic_row);

  for (int i = 0; i < 3; ++i) {
    d_rc[i] = d_rc_local[i];
    d_normal[i] = d_normal_local[i];
  }
}

// Wrapper for AIC column contribution (fixed vortex j, all control points)
static double aic_column_dot_wrapper(
    const double* rv1_j,
    const double* rv2_j,
    const double* rc_all,
    const double* normal_all,
    double chord_j,
    int surface_idx_j,
    const int* surface_idx_all,
    int num_vortices,
    double core_chord_fraction,
    double core_width_fraction,
    double beta,
    const double* weights) {

  double result = 0.0;

  // Compute strip width for this vortex
  double dy = rv2_j[1] - rv1_j[1];
  double dz = rv2_j[2] - rv1_j[2];
  double strip_width = std::sqrt(dy*dy + dz*dz);

  for (int i = 0; i < num_vortices; ++i) {
    const double* rc = &rc_all[3*i];
    const double* normal = &normal_all[3*i];

    // Core radius (depends on whether same surface)
    double core_radius = 0.0;
    if (surface_idx_j != surface_idx_all[i]) {
      double chord_core = core_chord_fraction * chord_j;
      double width_core = core_width_fraction * strip_width;
      core_radius = (chord_core > width_core) ? chord_core : width_core;
    }

    bool include_bound = (surface_idx_j != surface_idx_all[i]);
    double aic_elem = compute_aic_element(rv1_j, rv2_j, rc, normal,
                                           core_radius, beta, include_bound);
    result += weights[i] * aic_elem;
  }
  return result;
}

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
    double* d_rv2) {

  // Initialize gradients
  for (int i = 0; i < 3; ++i) {
    d_rv1[i] = 0.0;
    d_rv2[i] = 0.0;
  }

  const double* rv1_j = &rv1_all[3 * vortex_idx];
  const double* rv2_j = &rv2_all[3 * vortex_idx];

  // Make local copies for Enzyme
  double rv1_local[3] = {rv1_j[0], rv1_j[1], rv1_j[2]};
  double rv2_local[3] = {rv2_j[0], rv2_j[1], rv2_j[2]};
  double d_rv1_local[3] = {0.0, 0.0, 0.0};
  double d_rv2_local[3] = {0.0, 0.0, 0.0};

  __enzyme_autodiff<void>(
      (void*)aic_column_dot_wrapper,
      enzyme_dup, rv1_local, d_rv1_local,
      enzyme_dup, rv2_local, d_rv2_local,
      enzyme_const, rc_all,
      enzyme_const, normal_all,
      enzyme_const, chord_all[vortex_idx],
      enzyme_const, surface_idx_all[vortex_idx],
      enzyme_const, surface_idx_all,
      enzyme_const, num_vortices,
      enzyme_const, core_chord_fraction,
      enzyme_const, core_width_fraction,
      enzyme_const, beta,
      enzyme_const, d_aic_column);

  for (int i = 0; i < 3; ++i) {
    d_rv1[i] = d_rv1_local[i];
    d_rv2[i] = d_rv2_local[i];
  }
}

#endif  // VLM_USE_ENZYME

// ============================================================================
// Finite Difference Fallback
// ============================================================================

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
    double h) {

  double rv1_pert[3], rv2_pert[3], rc_pert[3], normal_pert[3];

  // Gradient w.r.t. rv1
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) rv1_pert[j] = rv1[j];
    rv1_pert[i] += h;
    double plus = compute_aic_element(rv1_pert, rv2, rc, normal,
                                       core_radius, beta, include_bound);
    rv1_pert[i] = rv1[i] - h;
    double minus = compute_aic_element(rv1_pert, rv2, rc, normal,
                                        core_radius, beta, include_bound);
    d_rv1[i] = (plus - minus) / (2.0 * h);
  }

  // Gradient w.r.t. rv2
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) rv2_pert[j] = rv2[j];
    rv2_pert[i] += h;
    double plus = compute_aic_element(rv1, rv2_pert, rc, normal,
                                       core_radius, beta, include_bound);
    rv2_pert[i] = rv2[i] - h;
    double minus = compute_aic_element(rv1, rv2_pert, rc, normal,
                                        core_radius, beta, include_bound);
    d_rv2[i] = (plus - minus) / (2.0 * h);
  }

  // Gradient w.r.t. rc
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) rc_pert[j] = rc[j];
    rc_pert[i] += h;
    double plus = compute_aic_element(rv1, rv2, rc_pert, normal,
                                       core_radius, beta, include_bound);
    rc_pert[i] = rc[i] - h;
    double minus = compute_aic_element(rv1, rv2, rc_pert, normal,
                                        core_radius, beta, include_bound);
    d_rc[i] = (plus - minus) / (2.0 * h);
  }

  // Gradient w.r.t. normal
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) normal_pert[j] = normal[j];
    normal_pert[i] += h;
    double plus = compute_aic_element(rv1, rv2, rc, normal_pert,
                                       core_radius, beta, include_bound);
    normal_pert[i] = normal[i] - h;
    double minus = compute_aic_element(rv1, rv2, rc, normal_pert,
                                        core_radius, beta, include_bound);
    d_normal[i] = (plus - minus) / (2.0 * h);
  }
}

}  // namespace ad
}  // namespace vlm
