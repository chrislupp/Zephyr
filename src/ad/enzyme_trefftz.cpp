// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include "vlm/ad/enzyme_trefftz.h"

namespace vlm {
namespace ad {

// ============================================================================
// Enzyme Gradient Implementations
// ============================================================================

#ifdef VLM_USE_ENZYME

// Wrapper function for Enzyme that outputs to a specific result index
static void trefftz_kernel_wrapper(
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
    double* result) {  // result[3] = {CL_ff, CY_ff, CDi_ff}

  compute_trefftz_drag_kernel(
      strip_gamma, strip_te_left, strip_te_right, strip_chord,
      strip_surface_index, num_strips, reference_area, reference_span,
      core_chord_fraction, core_width_fraction,
      &result[0], &result[1], &result[2]);
}

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
    double* d_gamma) {

  // Initialize gradient output
  for (int i = 0; i < num_strips; ++i) {
    d_gamma[i] = 0.0;
  }

  // Allocate shadow for gamma
  double* d_gamma_local = new double[num_strips];
  for (int i = 0; i < num_strips; ++i) {
    d_gamma_local[i] = 0.0;
  }

  double result[3] = {0.0, 0.0, 0.0};
  double d_result[3] = {0.0, 0.0, 0.0};
  d_result[output_index] = 1.0;  // Seed for the specific output

  // Reverse mode: compute d(result)/d(gamma)
  __enzyme_autodiff<void>(
      (void*)trefftz_kernel_wrapper,
      enzyme_dup, strip_gamma, d_gamma_local,
      enzyme_const, strip_te_left,
      enzyme_const, strip_te_right,
      enzyme_const, strip_chord,
      enzyme_const, strip_surface_index,
      enzyme_const, num_strips,
      enzyme_const, reference_area,
      enzyme_const, reference_span,
      enzyme_const, core_chord_fraction,
      enzyme_const, core_width_fraction,
      enzyme_dupnoneed, result, d_result);

  // Copy result
  for (int i = 0; i < num_strips; ++i) {
    d_gamma[i] = d_gamma_local[i];
  }

  delete[] d_gamma_local;
}

#endif  // VLM_USE_ENZYME

// ============================================================================
// Finite Difference Fallback
// ============================================================================

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
    double h) {

  // Allocate perturbed gamma array
  double* gamma_pert = new double[num_strips];

  for (int i = 0; i < num_strips; ++i) {
    // Copy original values
    for (int j = 0; j < num_strips; ++j) {
      gamma_pert[j] = strip_gamma[j];
    }

    // Positive perturbation
    gamma_pert[i] += h;
    double CL_plus, CY_plus, CDi_plus;
    compute_trefftz_drag_kernel(
        gamma_pert, strip_te_left, strip_te_right, strip_chord,
        strip_surface_index, num_strips, reference_area, reference_span,
        core_chord_fraction, core_width_fraction,
        &CL_plus, &CY_plus, &CDi_plus);

    // Negative perturbation
    gamma_pert[i] = strip_gamma[i] - h;
    double CL_minus, CY_minus, CDi_minus;
    compute_trefftz_drag_kernel(
        gamma_pert, strip_te_left, strip_te_right, strip_chord,
        strip_surface_index, num_strips, reference_area, reference_span,
        core_chord_fraction, core_width_fraction,
        &CL_minus, &CY_minus, &CDi_minus);

    // Compute finite difference
    switch (output_index) {
      case 0:
        d_gamma[i] = (CL_plus - CL_minus) / (2.0 * h);
        break;
      case 1:
        d_gamma[i] = (CY_plus - CY_minus) / (2.0 * h);
        break;
      case 2:
        d_gamma[i] = (CDi_plus - CDi_minus) / (2.0 * h);
        break;
      default:
        d_gamma[i] = 0.0;
    }
  }

  delete[] gamma_pert;
}

}  // namespace ad
}  // namespace vlm
