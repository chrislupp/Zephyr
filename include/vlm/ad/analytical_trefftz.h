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

/// @file analytical_trefftz.h
/// @brief Analytical derivatives of Trefftz plane drag computation.
///
/// The Trefftz plane drag uses 2D Biot-Savart (linear in Gamma) plus a
/// quadratic coupling (Gamma_i * K_ij * Gamma_j). The analytical derivative
/// is straightforward.

#include <cmath>

namespace vlm {
namespace ad {

constexpr double kHalfPiInvTrefftz = 0.15915494309189533577;  // 1/(2*pi)
constexpr double kEpsTrefftz = 1e-12;

// ============================================================================
// Primal Kernel
// ============================================================================

inline void compute_trefftz_velocity_point(double gamma, double yv, double zv,
                                           double yc, double zc, double core_sq,
                                           double* vy, double* vz) {
    double dy = yc - yv;
    double dz = zc - zv;
    double rsq = dy * dy + dz * dz;
    double denom = rsq + core_sq;
    if (denom > kEpsTrefftz) {
        *vy += kHalfPiInvTrefftz * gamma * dz / denom;
        *vz -= kHalfPiInvTrefftz * gamma * dy / denom;
    }
}

inline void compute_trefftz_drag_kernel(
    const double* strip_gamma, const double* strip_te_left,
    const double* strip_te_right, const double* strip_chord,
    const int* strip_surface_index, int num_strips, double reference_area,
    double reference_span, double core_chord_fraction,
    double core_width_fraction, double* CL_ff, double* CY_ff, double* CDi_ff) {
    (void)reference_span;
    *CL_ff = 0.0;
    *CY_ff = 0.0;
    *CDi_ff = 0.0;

    for (int jc = 0; jc < num_strips; ++jc) {
        double y_left = strip_te_left[3 * jc + 1];
        double z_left = strip_te_left[3 * jc + 2];
        double y_right = strip_te_right[3 * jc + 1];
        double z_right = strip_te_right[3 * jc + 2];

        double dy_c = y_right - y_left;
        double dz_c = z_right - z_left;
        double ds_c = std::sqrt(dy_c * dy_c + dz_c * dz_c);
        if (ds_c < kEpsTrefftz)
            continue;

        double yc = 0.5 * (y_left + y_right);
        double zc = 0.5 * (z_left + z_right);

        double vy = 0.0, vz = 0.0;
        for (int jv = 0; jv < num_strips; ++jv) {
            double y_left_v = strip_te_left[3 * jv + 1];
            double z_left_v = strip_te_left[3 * jv + 2];
            double y_right_v = strip_te_right[3 * jv + 1];
            double z_right_v = strip_te_right[3 * jv + 2];

            double dy_v = y_right_v - y_left_v;
            double dz_v = z_right_v - z_left_v;
            double ds_v = std::sqrt(dy_v * dy_v + dz_v * dz_v);
            if (ds_v < kEpsTrefftz)
                continue;

            double core_radius = 0.0;
            if (strip_surface_index[jc] != strip_surface_index[jv]) {
                double chord_core = core_chord_fraction * strip_chord[jv];
                double width_core = core_width_fraction * ds_v;
                core_radius =
                    (chord_core > width_core) ? chord_core : width_core;
            }
            double core_sq = core_radius * core_radius;

            compute_trefftz_velocity_point(strip_gamma[jv], y_left_v, z_left_v,
                                           yc, zc, core_sq, &vy, &vz);
            compute_trefftz_velocity_point(-strip_gamma[jv], y_right_v,
                                           z_right_v, yc, zc, core_sq, &vy,
                                           &vz);
        }

        *CL_ff += 2.0 * strip_gamma[jc] * dy_c / reference_area;
        *CY_ff -= 2.0 * strip_gamma[jc] * dz_c / reference_area;
        *CDi_ff += strip_gamma[jc] * (dz_c * vy - dy_c * vz) / reference_area;
    }
}

inline double compute_span_efficiency(double CL_ff, double CY_ff, double CDi_ff,
                                      double reference_area,
                                      double reference_span) {
    constexpr double kPi = 3.14159265358979323846;
    double AR = reference_span * reference_span / reference_area;
    if (std::abs(CDi_ff) > kEpsTrefftz) {
        return (CL_ff * CL_ff + CY_ff * CY_ff) / (kPi * AR * CDi_ff);
    }
    return 0.0;
}

// ============================================================================
// Analytical Gradient
// ============================================================================

/// @brief Analytical gradient of Trefftz outputs w.r.t. circulation.
void analytical_trefftz_circulation_gradient(
    const double* strip_gamma, const double* strip_te_left,
    const double* strip_te_right, const double* strip_chord,
    const int* strip_surface_index, int num_strips, double reference_area,
    double reference_span, double core_chord_fraction,
    double core_width_fraction, int output_index, double* d_gamma);

}  // namespace ad
}  // namespace vlm
