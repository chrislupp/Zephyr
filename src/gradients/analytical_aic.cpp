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

#include "pegasus/gradients/analytical_aic.h"

#include <cmath>

namespace pegasus::gradients {

namespace {

/// Compute core radius for a vortex, matching solver logic.
/// same_component: true when both vortices are on the same aerodynamic
/// component. For same-component, uses tiny core (just prevents division by
/// zero). For different-component, uses larger core to regularize near-field
/// interactions.
inline double compute_core_radius(double strip_width, bool same_component,
                                  double chord, double core_chord_fraction,
                                  double core_width_fraction) {
    if (same_component) {
        return 0.0001 * strip_width;
    }
    double chord_core = core_chord_fraction * chord;
    double width_core = core_width_fraction * strip_width;
    return (chord_core > width_core) ? chord_core : width_core;
}

inline double compute_strip_width(const double* rv1, const double* rv2) {
    double dy = rv2[1] - rv1[1];
    double dz = rv2[2] - rv1[2];
    return std::sqrt(dy * dy + dz * dz);
}

}  // namespace

// ============================================================================
// AIC Element Gradient
// ============================================================================

void analytical_aic_element_gradient(const double* rv1, const double* rv2,
                                     const double* rc, const double* normal,
                                     double core_radius, double prandtl_glauert,
                                     bool include_bound, double* d_rv1,
                                     double* d_rv2, double* d_rc,
                                     double* d_normal) {
    // AIC = V . normal, where V = horseshoe_velocity(rc, rv1, rv2)
    // d(AIC)/d(rc) = VJP(seed=normal) w.r.t. field_point=rc
    // d(AIC)/d(rv1) = VJP(seed=normal) w.r.t. bound_left=rv1
    // d(AIC)/d(rv2) = VJP(seed=normal) w.r.t. bound_right=rv2
    // d(AIC)/d(normal) = V (linear)

    horseshoe_velocity_vjp_field_point(rc, rv1, rv2, prandtl_glauert,
                                       core_radius, include_bound, normal,
                                       d_rc);
    horseshoe_velocity_vjp_bound_left(rc, rv1, rv2, prandtl_glauert,
                                      core_radius, include_bound, normal,
                                      d_rv1);
    horseshoe_velocity_vjp_bound_right(rc, rv1, rv2, prandtl_glauert,
                                       core_radius, include_bound, normal,
                                       d_rv2);

    // Normal gradient is just the velocity
    double velocity[3];
    compute_horseshoe_velocity_kernel(rc, rv1, rv2, prandtl_glauert,
                                      core_radius, include_bound, velocity);
    d_normal[0] = velocity[0];
    d_normal[1] = velocity[1];
    d_normal[2] = velocity[2];
}

// ============================================================================
// AIC Row Gradient
// ============================================================================

void analytical_aic_row_gradient_rc(
    const double* rv1_all, const double* rv2_all, const double* rc,
    const double* normal, const double* chord_all, const int* component_idx_all,
    int num_vortices, int row_idx, int control_component_idx,
    double core_chord_fraction, double core_width_fraction, double beta,
    const double* weights, double* d_rc, double* d_normal) {
    d_rc[0] = d_rc[1] = d_rc[2] = 0.0;
    d_normal[0] = d_normal[1] = d_normal[2] = 0.0;

    for (int j = 0; j < num_vortices; ++j) {
        if (weights[j] == 0.0)
            continue;

        const double* rv1 = &rv1_all[3 * j];
        const double* rv2 = &rv2_all[3 * j];

        double strip_width = compute_strip_width(rv1, rv2);
        bool same_component = (component_idx_all[j] == control_component_idx);
        double core_radius =
            compute_core_radius(strip_width, same_component, chord_all[j],
                                core_chord_fraction, core_width_fraction);
        // Match solver: only exclude bound vortex for self-influence (i == j)
        bool include_bound = (j != row_idx);

        // d(AIC[i,j])/d(rc): VJP of horseshoe velocity w.r.t. field point,
        // seeded with normal
        double vjp_rc[3];
        horseshoe_velocity_vjp_field_point(rc, rv1, rv2, beta, core_radius,
                                           include_bound, normal, vjp_rc);

        // Accumulate weighted by gamma_j (or whatever weights are)
        d_rc[0] += weights[j] * vjp_rc[0];
        d_rc[1] += weights[j] * vjp_rc[1];
        d_rc[2] += weights[j] * vjp_rc[2];

        // d(AIC[i,j])/d(normal) = V[i,j] (velocity induced at rc by vortex j)
        double velocity[3];
        compute_horseshoe_velocity_kernel(rc, rv1, rv2, beta, core_radius,
                                          include_bound, velocity);
        d_normal[0] += weights[j] * velocity[0];
        d_normal[1] += weights[j] * velocity[1];
        d_normal[2] += weights[j] * velocity[2];
    }
}

// ============================================================================
// AIC Column Gradient
// ============================================================================

void analytical_aic_column_gradient_rv(
    const double* rv1_all, const double* rv2_all, const double* rc_all,
    const double* normal_all, const double* chord_all,
    const int* component_idx_all, int num_vortices, int vortex_idx,
    double core_chord_fraction, double core_width_fraction, double beta,
    const double* weights, double* d_rv1, double* d_rv2) {
    d_rv1[0] = d_rv1[1] = d_rv1[2] = 0.0;
    d_rv2[0] = d_rv2[1] = d_rv2[2] = 0.0;

    const double* rv1_j = &rv1_all[3 * vortex_idx];
    const double* rv2_j = &rv2_all[3 * vortex_idx];
    int component_j = component_idx_all[vortex_idx];
    double strip_width_j = compute_strip_width(rv1_j, rv2_j);

    for (int i = 0; i < num_vortices; ++i) {
        if (weights[i] == 0.0)
            continue;

        const double* rc = &rc_all[3 * i];
        const double* normal = &normal_all[3 * i];

        bool same_component = (component_j == component_idx_all[i]);
        double core_radius = compute_core_radius(
            strip_width_j, same_component, chord_all[vortex_idx],
            core_chord_fraction, core_width_fraction);
        // Match solver: only exclude bound vortex for self-influence (i == j)
        bool include_bound = (i != vortex_idx);

        // d(AIC[i,j])/d(rv1_j): VJP seeded with normal_i
        double vjp_rv1[3], vjp_rv2[3];
        horseshoe_velocity_vjp_bound_left(rc, rv1_j, rv2_j, beta, core_radius,
                                          include_bound, normal, vjp_rv1);
        horseshoe_velocity_vjp_bound_right(rc, rv1_j, rv2_j, beta, core_radius,
                                           include_bound, normal, vjp_rv2);

        d_rv1[0] += weights[i] * vjp_rv1[0];
        d_rv1[1] += weights[i] * vjp_rv1[1];
        d_rv1[2] += weights[i] * vjp_rv1[2];
        d_rv2[0] += weights[i] * vjp_rv2[0];
        d_rv2[1] += weights[i] * vjp_rv2[1];
        d_rv2[2] += weights[i] * vjp_rv2[2];
    }
}

// ============================================================================
// Velocity Row Gradient
// ============================================================================

void analytical_velocity_row_gradient_rc(
    const double* rv1_all, const double* rv2_all, const double* rc,
    const double* chord_all, const int* component_idx_all, int num_vortices,
    int row_idx, int control_component_idx, double core_chord_fraction,
    double core_width_fraction, double beta, const double* gamma,
    const double* dJ_dV, double* d_rc) {
    d_rc[0] = d_rc[1] = d_rc[2] = 0.0;

    // We want: d(sum_j dJ_dV . V[i,j] * gamma_j) / d(rc_i)
    // = sum_j gamma_j * dJ_dV^T * (dV[i,j]/d(rc_i))
    // = sum_j gamma_j * VJP(seed=dJ_dV) w.r.t. field_point

    for (int j = 0; j < num_vortices; ++j) {
        if (gamma[j] == 0.0)
            continue;

        const double* rv1 = &rv1_all[3 * j];
        const double* rv2 = &rv2_all[3 * j];

        double strip_width = compute_strip_width(rv1, rv2);
        bool same_component = (component_idx_all[j] == control_component_idx);
        double core_radius =
            compute_core_radius(strip_width, same_component, chord_all[j],
                                core_chord_fraction, core_width_fraction);
        // Match solver: only exclude bound vortex for self-influence (i == j)
        bool include_bound = (j != row_idx);

        double vjp[3];
        horseshoe_velocity_vjp_field_point(rc, rv1, rv2, beta, core_radius,
                                           include_bound, dJ_dV, vjp);

        d_rc[0] += gamma[j] * vjp[0];
        d_rc[1] += gamma[j] * vjp[1];
        d_rc[2] += gamma[j] * vjp[2];
    }
}

// ============================================================================
// Velocity Column Gradient
// ============================================================================

void analytical_velocity_column_gradient_rv(
    const double* rv1_all, const double* rv2_all, const double* rc_all,
    const double* chord_all, const int* component_idx_all, int num_vortices,
    int vortex_idx, double core_chord_fraction, double core_width_fraction,
    double beta, double gamma_j, const double* dJ_dV_all, double* d_rv1,
    double* d_rv2) {
    d_rv1[0] = d_rv1[1] = d_rv1[2] = 0.0;
    d_rv2[0] = d_rv2[1] = d_rv2[2] = 0.0;

    if (gamma_j == 0.0)
        return;

    const double* rv1_j = &rv1_all[3 * vortex_idx];
    const double* rv2_j = &rv2_all[3 * vortex_idx];
    int component_j = component_idx_all[vortex_idx];
    double strip_width_j = compute_strip_width(rv1_j, rv2_j);

    for (int i = 0; i < num_vortices; ++i) {
        const double* rc = &rc_all[3 * i];
        const double* dJ_dV = &dJ_dV_all[3 * i];

        bool same_component = (component_j == component_idx_all[i]);
        double core_radius = compute_core_radius(
            strip_width_j, same_component, chord_all[vortex_idx],
            core_chord_fraction, core_width_fraction);
        // Match solver: only exclude bound vortex for self-influence (i == j)
        bool include_bound = (i != vortex_idx);

        double vjp_rv1[3], vjp_rv2[3];
        horseshoe_velocity_vjp_bound_left(rc, rv1_j, rv2_j, beta, core_radius,
                                          include_bound, dJ_dV, vjp_rv1);
        horseshoe_velocity_vjp_bound_right(rc, rv1_j, rv2_j, beta, core_radius,
                                           include_bound, dJ_dV, vjp_rv2);

        d_rv1[0] += gamma_j * vjp_rv1[0];
        d_rv1[1] += gamma_j * vjp_rv1[1];
        d_rv1[2] += gamma_j * vjp_rv1[2];
        d_rv2[0] += gamma_j * vjp_rv2[0];
        d_rv2[1] += gamma_j * vjp_rv2[1];
        d_rv2[2] += gamma_j * vjp_rv2[2];
    }
}

}  // namespace pegasus::gradients
