/*
    Project VLM
    Version 0.0.0

    Copyright 2024, Christopher A. Lupp

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

/// @file analytical_forces.h
/// @brief Analytical derivatives of force coefficient computation.
///
/// The Kutta-Joukowski force F = 2*Gamma*(V_eff x dl) is bilinear, so all
/// derivatives are straightforward cross-product algebra.

#include <cmath>

namespace vlm {
namespace ad {

// ============================================================================
// Freestream Direction and its Gradients (already analytical)
// ============================================================================

inline void compute_freestream_direction(double alpha, double beta,
                                         double* v_inf) {
    double ca = std::cos(alpha);
    double sa = std::sin(alpha);
    double cb = std::cos(beta);
    double sb = std::sin(beta);
    v_inf[0] = ca * cb;
    v_inf[1] = sb;
    v_inf[2] = sa * cb;
}

inline void compute_freestream_alpha_gradient(double alpha, double beta,
                                              double* dv_dalpha) {
    double ca = std::cos(alpha);
    double sa = std::sin(alpha);
    double cb = std::cos(beta);
    dv_dalpha[0] = -sa * cb;
    dv_dalpha[1] = 0.0;
    dv_dalpha[2] = ca * cb;
}

inline void compute_freestream_beta_gradient(double alpha, double beta,
                                             double* dv_dbeta) {
    double ca = std::cos(alpha);
    double sa = std::sin(alpha);
    double cb = std::cos(beta);
    double sb = std::sin(beta);
    dv_dbeta[0] = -ca * sb;
    dv_dbeta[1] = cb;
    dv_dbeta[2] = -sa * sb;
}

// ============================================================================
// Body-to-Stability Transformation
// ============================================================================

inline void body_to_stability(double cf_x, double cf_z, double alpha,
                              double* CL, double* CD) {
    double ca = std::cos(alpha);
    double sa = std::sin(alpha);
    *CL = cf_z * ca - cf_x * sa;
    *CD = cf_x * ca + cf_z * sa;
}

// ============================================================================
// Force Coefficient Kernel (primal)
// ============================================================================

inline void compute_force_coefficients_kernel(
    double alpha, double beta, const double* circulation,
    const double* induced_velocity, const double* bound_left,
    const double* bound_right, [[maybe_unused]] const double* control_points,
    const double* moment_arm, int num_vortices, double reference_area,
    double* result) {
    double cf_x = 0.0, cf_y = 0.0, cf_z = 0.0;
    double cm_x = 0.0, cm_y = 0.0, cm_z = 0.0;

    double v_inf[3];
    compute_freestream_direction(alpha, beta, v_inf);

    for (int i = 0; i < num_vortices; ++i) {
        double v_eff[3] = {v_inf[0] + induced_velocity[3 * i + 0],
                           v_inf[1] + induced_velocity[3 * i + 1],
                           v_inf[2] + induced_velocity[3 * i + 2]};

        double dl[3] = {bound_right[3 * i + 0] - bound_left[3 * i + 0],
                        bound_right[3 * i + 1] - bound_left[3 * i + 1],
                        bound_right[3 * i + 2] - bound_left[3 * i + 2]};

        double cross[3] = {v_eff[1] * dl[2] - v_eff[2] * dl[1],
                           v_eff[2] * dl[0] - v_eff[0] * dl[2],
                           v_eff[0] * dl[1] - v_eff[1] * dl[0]};

        double f[3] = {2.0 * circulation[i] * cross[0],
                       2.0 * circulation[i] * cross[1],
                       2.0 * circulation[i] * cross[2]};

        cf_x += f[0];
        cf_y += f[1];
        cf_z += f[2];

        const double* r = &moment_arm[3 * i];
        cm_x += r[1] * f[2] - r[2] * f[1];
        cm_y += r[2] * f[0] - r[0] * f[2];
        cm_z += r[0] * f[1] - r[1] * f[0];
    }

    cf_x /= reference_area;
    cf_y /= reference_area;
    cf_z /= reference_area;
    cm_x /= reference_area;
    cm_y /= reference_area;
    cm_z /= reference_area;

    double CL, CD;
    body_to_stability(cf_x, cf_z, alpha, &CL, &CD);

    result[0] = CL;
    result[1] = CD;
    result[2] = cf_y;
    result[3] = -cm_x;
    result[4] = cm_y;
    result[5] = -cm_z;
}

// ============================================================================
// Analytical Gradient Functions
// ============================================================================

/// @brief Analytical gradient of force coefficients w.r.t. alpha.
void analytical_force_alpha_gradient(double alpha, double beta,
                                     const double* circulation,
                                     const double* induced_velocity,
                                     const double* bound_left,
                                     const double* bound_right,
                                     const double* control_points,
                                     const double* moment_arm, int num_vortices,
                                     double reference_area, double* d_result);

/// @brief Analytical gradient of force coefficients w.r.t. beta.
void analytical_force_beta_gradient(double alpha, double beta,
                                    const double* circulation,
                                    const double* induced_velocity,
                                    const double* bound_left,
                                    const double* bound_right,
                                    const double* control_points,
                                    const double* moment_arm, int num_vortices,
                                    double reference_area, double* d_result);

/// @brief Analytical gradient of one output w.r.t. circulation (reverse mode).
void analytical_force_circulation_gradient(
    double alpha, double beta, const double* circulation,
    const double* induced_velocity, const double* bound_left,
    const double* bound_right, const double* control_points,
    const double* moment_arm, int num_vortices, double reference_area,
    int output_index, double* d_circulation);

/// @brief Analytical gradient of one output w.r.t. induced velocity (reverse
/// mode).
void analytical_force_induced_velocity_gradient(
    double alpha, double beta, const double* circulation,
    const double* induced_velocity, const double* bound_left,
    const double* bound_right, const double* control_points,
    const double* moment_arm, int num_vortices, double reference_area,
    int output_index, double* d_induced_velocity);

/// @brief Analytical gradient of force coefficients w.r.t. mesh geometry
/// (forward mode).
void analytical_force_geometry_gradient(
    double alpha, double beta, const double* circulation,
    const double* induced_velocity, const double* bound_left,
    const double* bound_right, const double* moment_arm, int num_vortices,
    double reference_area, const double* d_bound_left,
    const double* d_bound_right, const double* d_moment_arm, double* d_result);

}  // namespace ad
}  // namespace vlm
