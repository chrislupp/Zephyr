/*
    Project Pegasus
    Version 0.2.0

    Copyright 2024-2026, Christopher A. Lupp

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#pragma once

/// @file analytical_aic.h
/// @brief Analytical AIC matrix derivatives using Biot-Savart VJPs.
///
/// Provides row and column gradient functions for the AIC matrix that are
/// used by the adjoint method for geometry sensitivity computation.

#include "pegasus/gradients/analytical_biot_savart.h"

#include <cmath>

namespace pegasus::gradients {

// ============================================================================
// AIC Element Gradient
// ============================================================================

/// @brief Compute gradient of single AIC element w.r.t. all geometric inputs.
///
/// AIC = V(rc, rv1, rv2) . normal
/// d(AIC)/d(rc) = (dV/d(rc))^T * normal  [via VJP with seed=normal]
/// d(AIC)/d(rv1) = (dV/d(rv1))^T * normal
/// d(AIC)/d(rv2) = (dV/d(rv2))^T * normal
/// d(AIC)/d(normal) = V  [linear in normal]
void analytical_aic_element_gradient(const double* rv1, const double* rv2,
                                     const double* rc, const double* normal,
                                     double core_radius, double prandtl_glauert,
                                     bool include_bound, double* d_rv1,
                                     double* d_rv2, double* d_rc,
                                     double* d_normal);

// ============================================================================
// AIC Row/Column Gradients (for adjoint method)
// ============================================================================

/// @brief Row gradient: d(sum_j AIC[i,j]*weights[j]) / d(rc_i, normal_i)
///
/// Used in adjoint: lambda[i] * d(AIC_row . Gamma)/d(rc, normal)
///
/// @param rv1_all Left bound vortex endpoints [3*num_vortices]
/// @param rv2_all Right bound vortex endpoints [3*num_vortices]
/// @param rc Control point for this row [3]
/// @param normal Panel normal for this row [3]
/// @param chord_all Chord values [num_vortices]
/// @param surface_idx_all Surface index for each vortex [num_vortices]
/// @param num_vortices Total number of vortices
/// @param row_idx Index of this row (for self-influence exclusion)
/// @param control_surface_idx Surface index of the control point
/// @param core_chord_fraction Core radius as fraction of chord
/// @param core_width_fraction Core radius as fraction of strip width
/// @param beta Prandtl-Glauert factor
/// @param weights Seed weights (typically Gamma) [num_vortices]
/// @param d_rc Output gradient w.r.t. control point [3]
/// @param d_normal Output gradient w.r.t. normal [3]
void analytical_aic_row_gradient_rc(
    const double* rv1_all, const double* rv2_all, const double* rc,
    const double* normal, const double* chord_all, const int* component_idx_all,
    int num_vortices, int row_idx, int control_component_idx,
    double core_chord_fraction, double core_width_fraction, double beta,
    const double* weights, double* d_rc, double* d_normal);

/// @brief Column gradient: d(sum_i weights[i]*AIC[i,j]) / d(rv1_j, rv2_j)
///
/// Used in adjoint: d(lambda^T . AIC_col)/d(rv1, rv2)
///
/// @param rv1_all Left bound vortex endpoints [3*num_vortices]
/// @param rv2_all Right bound vortex endpoints [3*num_vortices]
/// @param rc_all Control points [3*num_vortices]
/// @param normal_all Panel normals [3*num_vortices]
/// @param chord_all Chord values [num_vortices]
/// @param surface_idx_all Surface indices [num_vortices]
/// @param num_vortices Total number of vortices
/// @param vortex_idx Index of vortex being differentiated
/// @param core_chord_fraction Core radius chord fraction
/// @param core_width_fraction Core radius width fraction
/// @param beta Prandtl-Glauert factor
/// @param weights Seed weights (typically lambda) [num_vortices]
/// @param d_rv1 Output gradient w.r.t. rv1 [3]
/// @param d_rv2 Output gradient w.r.t. rv2 [3]
void analytical_aic_column_gradient_rv(
    const double* rv1_all, const double* rv2_all, const double* rc_all,
    const double* normal_all, const double* chord_all,
    const int* component_idx_all, int num_vortices, int vortex_idx,
    double core_chord_fraction, double core_width_fraction, double beta,
    const double* weights, double* d_rv1, double* d_rv2);

// ============================================================================
// Velocity Kernel Gradients (for explicit force geometry term)
// ============================================================================

/// @brief Velocity row gradient: d(sum_j dJ_dV . V[i,j] * gamma_j) / d(rc_i)
void analytical_velocity_row_gradient_rc(
    const double* rv1_all, const double* rv2_all, const double* rc,
    const double* chord_all, const int* component_idx_all, int num_vortices,
    int row_idx, int control_component_idx, double core_chord_fraction,
    double core_width_fraction, double beta, const double* gamma,
    const double* dJ_dV, double* d_rc);

/// @brief Velocity column gradient: d(sum_i dJ_dV_i . V[i,j] * gamma_j) /
/// d(rv1_j, rv2_j)
void analytical_velocity_column_gradient_rv(
    const double* rv1_all, const double* rv2_all, const double* rc_all,
    const double* chord_all, const int* component_idx_all, int num_vortices,
    int vortex_idx, double core_chord_fraction, double core_width_fraction,
    double beta, double gamma_j, const double* dJ_dV_all, double* d_rv1,
    double* d_rv2);

}  // namespace pegasus::gradients
