// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include "vlm/ad/enzyme_biot_savart.h"

namespace vlm {
namespace ad {

// ============================================================================
// Enzyme Gradient Implementations
// ============================================================================

#ifdef VLM_USE_ENZYME

// Wrapper function for Enzyme that computes scalar AIC element
// (needed because compute_aic_element_kernel returns double, not void)
static void aic_element_wrapper(
    const double* field_point,
    const double* normal,
    const double* bound_left,
    const double* bound_right,
    double prandtl_glauert,
    double core_radius,
    int include_bound_int,
    double* result) {
  bool include_bound = (include_bound_int != 0);
  *result = compute_aic_element_kernel(field_point, normal, bound_left, bound_right,
                                       prandtl_glauert, core_radius, include_bound);
}

void enzyme_aic_field_point_gradient(
    const double* field_point,
    const double* normal,
    const double* bound_left,
    const double* bound_right,
    double prandtl_glauert,
    double core_radius,
    bool include_bound,
    double* d_field_point) {

  // Initialize gradient output
  d_field_point[0] = 0.0;
  d_field_point[1] = 0.0;
  d_field_point[2] = 0.0;

  double result = 0.0;
  double d_result = 1.0;  // Seed for reverse mode

  int include_bound_int = include_bound ? 1 : 0;

  // Reverse mode: compute d(result)/d(field_point)
  __enzyme_autodiff<void>(
      (void*)aic_element_wrapper,
      enzyme_dup, field_point, d_field_point,  // differentiate w.r.t. field_point
      enzyme_const, normal,
      enzyme_const, bound_left,
      enzyme_const, bound_right,
      enzyme_const, prandtl_glauert,
      enzyme_const, core_radius,
      enzyme_const, include_bound_int,
      enzyme_dupnoneed, &result, &d_result);
}

void enzyme_aic_bound_gradient(
    const double* field_point,
    const double* normal,
    const double* bound_left,
    const double* bound_right,
    double prandtl_glauert,
    double core_radius,
    bool include_bound,
    double* d_bound_left,
    double* d_bound_right) {

  // Initialize gradient outputs
  for (int i = 0; i < 3; ++i) {
    d_bound_left[i] = 0.0;
    d_bound_right[i] = 0.0;
  }

  double result = 0.0;
  double d_result = 1.0;

  int include_bound_int = include_bound ? 1 : 0;

  // Reverse mode: compute d(result)/d(bound_left) and d(result)/d(bound_right)
  __enzyme_autodiff<void>(
      (void*)aic_element_wrapper,
      enzyme_const, field_point,
      enzyme_const, normal,
      enzyme_dup, bound_left, d_bound_left,
      enzyme_dup, bound_right, d_bound_right,
      enzyme_const, prandtl_glauert,
      enzyme_const, core_radius,
      enzyme_const, include_bound_int,
      enzyme_dupnoneed, &result, &d_result);
}

void enzyme_aic_normal_gradient(
    const double* field_point,
    const double* normal,
    const double* bound_left,
    const double* bound_right,
    double prandtl_glauert,
    double core_radius,
    bool include_bound,
    double* d_normal) {

  // For normal gradient, the AIC is linear in normal:
  // AIC = V · n, so d(AIC)/d(n) = V
  // We can compute this directly without AD
  double velocity[3];
  compute_horseshoe_velocity_kernel(field_point, bound_left, bound_right,
                                    prandtl_glauert, core_radius, include_bound,
                                    velocity);
  d_normal[0] = velocity[0];
  d_normal[1] = velocity[1];
  d_normal[2] = velocity[2];
}

#endif  // VLM_USE_ENZYME

// ============================================================================
// Finite Difference Fallbacks
// ============================================================================

void fd_aic_field_point_gradient(
    const double* field_point,
    const double* normal,
    const double* bound_left,
    const double* bound_right,
    double prandtl_glauert,
    double core_radius,
    bool include_bound,
    double* d_field_point,
    double h) {

  double fp_pert[3];

  for (int dim = 0; dim < 3; ++dim) {
    // Positive perturbation
    for (int i = 0; i < 3; ++i) fp_pert[i] = field_point[i];
    fp_pert[dim] += h;
    double aic_plus = compute_aic_element_kernel(fp_pert, normal, bound_left, bound_right,
                                                  prandtl_glauert, core_radius, include_bound);

    // Negative perturbation
    for (int i = 0; i < 3; ++i) fp_pert[i] = field_point[i];
    fp_pert[dim] -= h;
    double aic_minus = compute_aic_element_kernel(fp_pert, normal, bound_left, bound_right,
                                                   prandtl_glauert, core_radius, include_bound);

    d_field_point[dim] = (aic_plus - aic_minus) / (2.0 * h);
  }
}

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
    double h) {

  double bl_pert[3], br_pert[3];

  // Gradient w.r.t. bound_left
  for (int dim = 0; dim < 3; ++dim) {
    for (int i = 0; i < 3; ++i) bl_pert[i] = bound_left[i];
    bl_pert[dim] += h;
    double aic_plus = compute_aic_element_kernel(field_point, normal, bl_pert, bound_right,
                                                  prandtl_glauert, core_radius, include_bound);

    for (int i = 0; i < 3; ++i) bl_pert[i] = bound_left[i];
    bl_pert[dim] -= h;
    double aic_minus = compute_aic_element_kernel(field_point, normal, bl_pert, bound_right,
                                                   prandtl_glauert, core_radius, include_bound);

    d_bound_left[dim] = (aic_plus - aic_minus) / (2.0 * h);
  }

  // Gradient w.r.t. bound_right
  for (int dim = 0; dim < 3; ++dim) {
    for (int i = 0; i < 3; ++i) br_pert[i] = bound_right[i];
    br_pert[dim] += h;
    double aic_plus = compute_aic_element_kernel(field_point, normal, bound_left, br_pert,
                                                  prandtl_glauert, core_radius, include_bound);

    for (int i = 0; i < 3; ++i) br_pert[i] = bound_right[i];
    br_pert[dim] -= h;
    double aic_minus = compute_aic_element_kernel(field_point, normal, bound_left, br_pert,
                                                   prandtl_glauert, core_radius, include_bound);

    d_bound_right[dim] = (aic_plus - aic_minus) / (2.0 * h);
  }
}

}  // namespace ad
}  // namespace vlm
