// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include "vlm/ad/analytical_trefftz.h"

namespace vlm {
namespace ad {

void analytical_trefftz_circulation_gradient(
    const double* strip_gamma, const double* strip_te_left,
    const double* strip_te_right, const double* strip_chord,
    const int* strip_surface_index, int num_strips, double reference_area,
    double reference_span, double core_chord_fraction,
    double core_width_fraction, int output_index, double* d_gamma) {
    (void)reference_span;

    // Initialize output
    for (int i = 0; i < num_strips; ++i) {
        d_gamma[i] = 0.0;
    }

    // The Trefftz plane computation:
    // For each control strip jc:
    //   vy_jc = sum_jv K_y(jc, jv) * gamma_jv   (2D Biot-Savart, linear in
    //   gamma) vz_jc = sum_jv K_z(jc, jv) * gamma_jv CL_ff += 2 * gamma_jc *
    //   dy_jc / S_ref CY_ff -= 2 * gamma_jc * dz_jc / S_ref CDi_ff += gamma_jc
    //   * (dz_jc * vy_jc - dy_jc * vz_jc) / S_ref
    //
    // Derivatives:
    // d(CL_ff)/d(gamma_k) = 2 * dy_k / S_ref  [direct, only when jc=k]
    // d(CY_ff)/d(gamma_k) = -2 * dz_k / S_ref [direct, only when jc=k]
    // d(CDi_ff)/d(gamma_k) = [direct term] + [indirect through induced
    // velocity]
    //   Direct: (dz_k * vy_k - dy_k * vz_k) / S_ref  [when jc=k]
    //   Indirect: sum_jc gamma_jc * (dz_jc * dvy_jc/dgamma_k - dy_jc *
    //   dvz_jc/dgamma_k) / S_ref Since vy_jc = sum_jv K_y(jc,jv)*gamma_jv,
    //   dvy_jc/dgamma_k = K_y(jc,k)

    double inv_S = 1.0 / reference_area;

    // Precompute strip geometry
    struct StripGeom {
        double y_left, z_left, y_right, z_right;
        double dy, dz, ds, yc, zc;
    };

    auto* geom = new StripGeom[num_strips];
    for (int i = 0; i < num_strips; ++i) {
        geom[i].y_left = strip_te_left[3 * i + 1];
        geom[i].z_left = strip_te_left[3 * i + 2];
        geom[i].y_right = strip_te_right[3 * i + 1];
        geom[i].z_right = strip_te_right[3 * i + 2];
        geom[i].dy = geom[i].y_right - geom[i].y_left;
        geom[i].dz = geom[i].z_right - geom[i].z_left;
        geom[i].ds =
            std::sqrt(geom[i].dy * geom[i].dy + geom[i].dz * geom[i].dz);
        geom[i].yc = 0.5 * (geom[i].y_left + geom[i].y_right);
        geom[i].zc = 0.5 * (geom[i].z_left + geom[i].z_right);
    }

    if (output_index == 0) {
        // d(CL_ff)/d(gamma_k) = 2 * dy_k / S_ref
        for (int k = 0; k < num_strips; ++k) {
            if (geom[k].ds < kEpsTrefftz)
                continue;
            d_gamma[k] = 2.0 * geom[k].dy * inv_S;
        }
    } else if (output_index == 1) {
        // d(CY_ff)/d(gamma_k) = -2 * dz_k / S_ref
        for (int k = 0; k < num_strips; ++k) {
            if (geom[k].ds < kEpsTrefftz)
                continue;
            d_gamma[k] = -2.0 * geom[k].dz * inv_S;
        }
    } else if (output_index == 2) {
        // d(CDi_ff)/d(gamma_k)
        // = [direct] + [indirect]

        // First compute induced velocities (needed for direct term at k)
        auto* vy_arr = new double[num_strips];
        auto* vz_arr = new double[num_strips];

        for (int jc = 0; jc < num_strips; ++jc) {
            vy_arr[jc] = 0.0;
            vz_arr[jc] = 0.0;
            if (geom[jc].ds < kEpsTrefftz)
                continue;

            for (int jv = 0; jv < num_strips; ++jv) {
                if (geom[jv].ds < kEpsTrefftz)
                    continue;
                double core_radius = 0.0;
                if (strip_surface_index[jc] != strip_surface_index[jv]) {
                    double chord_core = core_chord_fraction * strip_chord[jv];
                    double width_core = core_width_fraction * geom[jv].ds;
                    core_radius =
                        (chord_core > width_core) ? chord_core : width_core;
                }
                double core_sq = core_radius * core_radius;

                compute_trefftz_velocity_point(strip_gamma[jv], geom[jv].y_left,
                                               geom[jv].z_left, geom[jc].yc,
                                               geom[jc].zc, core_sq,
                                               &vy_arr[jc], &vz_arr[jc]);
                compute_trefftz_velocity_point(
                    -strip_gamma[jv], geom[jv].y_right, geom[jv].z_right,
                    geom[jc].yc, geom[jc].zc, core_sq, &vy_arr[jc],
                    &vz_arr[jc]);
            }
        }

        for (int k = 0; k < num_strips; ++k) {
            if (geom[k].ds < kEpsTrefftz)
                continue;

            // Direct term: (dz_k * vy_k - dy_k * vz_k) / S_ref
            double direct =
                (geom[k].dz * vy_arr[k] - geom[k].dy * vz_arr[k]) * inv_S;

            // Indirect term: sum_jc gamma_jc * (dz_jc * K_y(jc,k) - dy_jc *
            // K_z(jc,k)) / S_ref where K_y(jc,k) = velocity_y at jc from unit
            // gamma at k K is from both the left and right TE points of strip
            // k:
            //   K_y(jc,k) = dvy/dgamma_k at strip jc
            //   For left TE: +kHalfPiInv * dz / denom
            //   For right TE: -kHalfPiInv * dz / denom (negative gamma)
            double indirect = 0.0;
            for (int jc = 0; jc < num_strips; ++jc) {
                if (geom[jc].ds < kEpsTrefftz)
                    continue;

                double core_radius = 0.0;
                if (strip_surface_index[jc] != strip_surface_index[k]) {
                    double chord_core = core_chord_fraction * strip_chord[k];
                    double width_core = core_width_fraction * geom[k].ds;
                    core_radius =
                        (chord_core > width_core) ? chord_core : width_core;
                }
                double core_sq = core_radius * core_radius;

                // Velocity kernel at jc from unit gamma at k (left TE)
                double dy_l = geom[jc].yc - geom[k].y_left;
                double dz_l = geom[jc].zc - geom[k].z_left;
                double rsq_l = dy_l * dy_l + dz_l * dz_l;
                double denom_l = rsq_l + core_sq;

                double ky_l = 0.0, kz_l = 0.0;
                if (denom_l > kEpsTrefftz) {
                    ky_l = kHalfPiInvTrefftz * dz_l / denom_l;
                    kz_l = -kHalfPiInvTrefftz * dy_l / denom_l;
                }

                // Right TE (negative circulation)
                double dy_r = geom[jc].yc - geom[k].y_right;
                double dz_r = geom[jc].zc - geom[k].z_right;
                double rsq_r = dy_r * dy_r + dz_r * dz_r;
                double denom_r = rsq_r + core_sq;

                double ky_r = 0.0, kz_r = 0.0;
                if (denom_r > kEpsTrefftz) {
                    ky_r = -kHalfPiInvTrefftz * dz_r / denom_r;
                    kz_r = kHalfPiInvTrefftz * dy_r / denom_r;
                }

                double ky_total = ky_l + ky_r;
                double kz_total = kz_l + kz_r;

                indirect += strip_gamma[jc] *
                            (geom[jc].dz * ky_total - geom[jc].dy * kz_total);
            }
            indirect *= inv_S;

            d_gamma[k] = direct + indirect;
        }

        delete[] vy_arr;
        delete[] vz_arr;
    }

    delete[] geom;
}

}  // namespace ad
}  // namespace vlm
