// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include "vlm/ad/enzyme_forces.h"

namespace vlm {
namespace ad {

// ============================================================================
// Enzyme Gradient Implementations
// ============================================================================

#ifdef VLM_USE_ENZYME

// Helper to compute forces and accumulate into result
// This is the primal function that Enzyme differentiates
static void force_kernel_wrapper(
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
    double* result) {
  compute_force_coefficients_kernel(
      alpha, beta, circulation, induced_velocity,
      bound_left, bound_right, control_points, moment_arm,
      num_vortices, reference_area, result);
}

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
    double* d_result) {
  // Use forward mode: seed alpha with 1.0, get d(result)/d(alpha)
  double result[6] = {0};
  double d_result_local[6] = {0};

  __enzyme_fwddiff<void>(
      (void*)force_kernel_wrapper,
      enzyme_dup, alpha, 1.0,           // differentiate w.r.t. alpha
      enzyme_const, beta,               // beta is constant
      enzyme_const, circulation,
      enzyme_const, induced_velocity,
      enzyme_const, bound_left,
      enzyme_const, bound_right,
      enzyme_const, control_points,
      enzyme_const, moment_arm,
      enzyme_const, num_vortices,
      enzyme_const, reference_area,
      enzyme_dup, result, d_result_local);  // output and its tangent

  // Copy results
  for (int i = 0; i < 6; ++i) {
    d_result[i] = d_result_local[i];
  }
}

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
    double* d_result) {
  // Use forward mode: seed beta with 1.0, get d(result)/d(beta)
  double result[6] = {0};
  double d_result_local[6] = {0};

  __enzyme_fwddiff<void>(
      (void*)force_kernel_wrapper,
      enzyme_const, alpha,              // alpha is constant
      enzyme_dup, beta, 1.0,            // differentiate w.r.t. beta
      enzyme_const, circulation,
      enzyme_const, induced_velocity,
      enzyme_const, bound_left,
      enzyme_const, bound_right,
      enzyme_const, control_points,
      enzyme_const, moment_arm,
      enzyme_const, num_vortices,
      enzyme_const, reference_area,
      enzyme_dup, result, d_result_local);

  for (int i = 0; i < 6; ++i) {
    d_result[i] = d_result_local[i];
  }
}

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
    double* d_circulation) {
  // Use reverse mode to get gradient of one output w.r.t. all circulation values
  // This gives ∂J/∂Γ needed for the adjoint RHS

  // Zero the gradient output
  for (int i = 0; i < num_vortices; ++i) {
    d_circulation[i] = 0.0;
  }

  // Create seed vector for the output we want
  double result[6] = {0};
  double d_result[6] = {0};
  d_result[output_index] = 1.0;  // Seed the output we're differentiating

  // Reverse mode differentiation
  __enzyme_autodiff<void>(
      (void*)force_kernel_wrapper,
      enzyme_const, alpha,
      enzyme_const, beta,
      enzyme_dup, circulation, d_circulation,  // input and its gradient
      enzyme_const, induced_velocity,
      enzyme_const, bound_left,
      enzyme_const, bound_right,
      enzyme_const, control_points,
      enzyme_const, moment_arm,
      enzyme_const, num_vortices,
      enzyme_const, reference_area,
      enzyme_dupnoneed, result, d_result);  // output with seed
}

#endif  // VLM_USE_ENZYME

// ============================================================================
// Finite Difference Fallbacks
// ============================================================================

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
    double h) {
  double result_plus[6], result_minus[6];

  compute_force_coefficients_kernel(
      alpha + h, beta, circulation, induced_velocity,
      bound_left, bound_right, control_points, moment_arm,
      num_vortices, reference_area, result_plus);

  compute_force_coefficients_kernel(
      alpha - h, beta, circulation, induced_velocity,
      bound_left, bound_right, control_points, moment_arm,
      num_vortices, reference_area, result_minus);

  for (int i = 0; i < 6; ++i) {
    d_result[i] = (result_plus[i] - result_minus[i]) / (2.0 * h);
  }
}

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
    double h) {
  double result_plus[6], result_minus[6];

  compute_force_coefficients_kernel(
      alpha, beta + h, circulation, induced_velocity,
      bound_left, bound_right, control_points, moment_arm,
      num_vortices, reference_area, result_plus);

  compute_force_coefficients_kernel(
      alpha, beta - h, circulation, induced_velocity,
      bound_left, bound_right, control_points, moment_arm,
      num_vortices, reference_area, result_minus);

  for (int i = 0; i < 6; ++i) {
    d_result[i] = (result_plus[i] - result_minus[i]) / (2.0 * h);
  }
}

}  // namespace ad
}  // namespace vlm
