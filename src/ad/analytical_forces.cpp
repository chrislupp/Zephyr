// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include "vlm/ad/analytical_forces.h"

namespace vlm {
namespace ad {

// ============================================================================
// Helper: Body-axis force/moment → output coefficients mapping
// ============================================================================

namespace {

/// Given body-axis forces (cf_x, cf_y, cf_z) and moments (cm_x, cm_y, cm_z),
/// compute the 6 output coefficients and the Jacobian of outputs w.r.t.
/// body-axis values.
///
/// Outputs: [CL, CD, CY, Cl, Cm, Cn]
/// CL = cf_z*cos(a) - cf_x*sin(a)
/// CD = cf_x*cos(a) + cf_z*sin(a)
/// CY = cf_y
/// Cl = -cm_x
/// Cm = cm_y
/// Cn = -cm_z
///
/// The "reverse-mode" operation is: given d_output[6], produce d_body[6].
struct BodyToOutputMap {
    double ca, sa;

    BodyToOutputMap(double alpha) : ca(std::cos(alpha)), sa(std::sin(alpha)) {}

    // Forward: body → output
    void forward(double cf_x, double cf_y, double cf_z, double cm_x,
                 double cm_y, double cm_z, double* out) const {
        out[0] = cf_z * ca - cf_x * sa;  // CL
        out[1] = cf_x * ca + cf_z * sa;  // CD
        out[2] = cf_y;                   // CY
        out[3] = -cm_x;                  // Cl
        out[4] = cm_y;                   // Cm
        out[5] = -cm_z;                  // Cn
    }

    // Reverse: d_output → d_body (d_cf_x, d_cf_y, d_cf_z, d_cm_x, d_cm_y,
    // d_cm_z)
    void reverse(const double* d_out, double* d_body) const {
        // d_cf_x = d_CL*(-sa) + d_CD*(ca)
        d_body[0] = d_out[0] * (-sa) + d_out[1] * ca;
        // d_cf_y = d_CY
        d_body[1] = d_out[2];
        // d_cf_z = d_CL*(ca) + d_CD*(sa)
        d_body[2] = d_out[0] * ca + d_out[1] * sa;
        // d_cm_x = d_Cl*(-1)
        d_body[3] = -d_out[3];
        // d_cm_y = d_Cm
        d_body[4] = d_out[4];
        // d_cm_z = d_Cn*(-1)
        d_body[5] = -d_out[5];
    }
};

}  // anonymous namespace

// ============================================================================
// Alpha Gradient
// ============================================================================

void analytical_force_alpha_gradient(double alpha, double beta,
                                     const double* circulation,
                                     const double* induced_velocity,
                                     const double* bound_left,
                                     const double* bound_right,
                                     const double* /*control_points*/,
                                     const double* moment_arm, int num_vortices,
                                     double reference_area, double* d_result) {
    // Alpha affects:
    // 1. Freestream direction: dV_inf/dalpha
    // 2. Body-to-stability rotation: d(CL,CD)/dalpha

    double v_inf[3], dv_dalpha[3];
    compute_freestream_direction(alpha, beta, v_inf);
    compute_freestream_alpha_gradient(alpha, beta, dv_dalpha);

    double ca = std::cos(alpha);
    double sa = std::sin(alpha);

    // Accumulate body-axis force/moment derivatives w.r.t. alpha
    double d_cf_x = 0.0, d_cf_y = 0.0, d_cf_z = 0.0;
    double d_cm_x = 0.0, d_cm_y = 0.0, d_cm_z = 0.0;

    // Also need the base forces for the rotation derivative
    double cf_x = 0.0, cf_z = 0.0;

    for (int i = 0; i < num_vortices; ++i) {
        double v_eff[3] = {v_inf[0] + induced_velocity[3 * i + 0],
                           v_inf[1] + induced_velocity[3 * i + 1],
                           v_inf[2] + induced_velocity[3 * i + 2]};

        double dl[3] = {bound_right[3 * i + 0] - bound_left[3 * i + 0],
                        bound_right[3 * i + 1] - bound_left[3 * i + 1],
                        bound_right[3 * i + 2] - bound_left[3 * i + 2]};

        // Base force for this vortex (for rotation derivative)
        double cross[3] = {v_eff[1] * dl[2] - v_eff[2] * dl[1],
                           v_eff[2] * dl[0] - v_eff[0] * dl[2],
                           v_eff[0] * dl[1] - v_eff[1] * dl[0]};
        cf_x += 2.0 * circulation[i] * cross[0];
        cf_z += 2.0 * circulation[i] * cross[2];

        // Derivative of force w.r.t. alpha (through V_inf change):
        // dF/dalpha = 2*Gamma * (dV_inf/dalpha x dl)
        double d_cross[3] = {dv_dalpha[1] * dl[2] - dv_dalpha[2] * dl[1],
                             dv_dalpha[2] * dl[0] - dv_dalpha[0] * dl[2],
                             dv_dalpha[0] * dl[1] - dv_dalpha[1] * dl[0]};

        double df[3] = {2.0 * circulation[i] * d_cross[0],
                        2.0 * circulation[i] * d_cross[1],
                        2.0 * circulation[i] * d_cross[2]};

        d_cf_x += df[0];
        d_cf_y += df[1];
        d_cf_z += df[2];

        const double* r = &moment_arm[3 * i];
        d_cm_x += r[1] * df[2] - r[2] * df[1];
        d_cm_y += r[2] * df[0] - r[0] * df[2];
        d_cm_z += r[0] * df[1] - r[1] * df[0];
    }

    double inv_S = 1.0 / reference_area;
    d_cf_x *= inv_S;
    d_cf_y *= inv_S;
    d_cf_z *= inv_S;
    d_cm_x *= inv_S;
    d_cm_y *= inv_S;
    d_cm_z *= inv_S;
    cf_x *= inv_S;
    cf_z *= inv_S;

    // Body-to-stability derivative has two parts:
    // CL = cf_z*cos(a) - cf_x*sin(a)
    // dCL/da = d_cf_z*cos(a) - d_cf_x*sin(a)          [force change]
    //        + cf_z*(-sin(a)) - cf_x*cos(a)            [rotation change]
    // CD = cf_x*cos(a) + cf_z*sin(a)
    // dCD/da = d_cf_x*cos(a) + d_cf_z*sin(a)
    //        + cf_x*(-sin(a)) + cf_z*cos(a)

    d_result[0] =
        (d_cf_z * ca - d_cf_x * sa) + (-cf_z * sa - cf_x * ca);  // dCL/da
    d_result[1] =
        (d_cf_x * ca + d_cf_z * sa) + (-cf_x * sa + cf_z * ca);  // dCD/da
    d_result[2] = d_cf_y;                                        // dCY/da
    d_result[3] = -d_cm_x;                                       // dCl/da
    d_result[4] = d_cm_y;                                        // dCm/da
    d_result[5] = -d_cm_z;                                       // dCn/da
}

// ============================================================================
// Beta Gradient
// ============================================================================

void analytical_force_beta_gradient(double alpha, double beta,
                                    const double* circulation,
                                    const double* induced_velocity,
                                    const double* bound_left,
                                    const double* bound_right,
                                    const double* /*control_points*/,
                                    const double* moment_arm, int num_vortices,
                                    double reference_area, double* d_result) {
    double v_inf[3], dv_dbeta[3];
    compute_freestream_direction(alpha, beta, v_inf);
    compute_freestream_beta_gradient(alpha, beta, dv_dbeta);

    double ca = std::cos(alpha);
    double sa = std::sin(alpha);

    double d_cf_x = 0.0, d_cf_y = 0.0, d_cf_z = 0.0;
    double d_cm_x = 0.0, d_cm_y = 0.0, d_cm_z = 0.0;

    for (int i = 0; i < num_vortices; ++i) {
        double dl[3] = {bound_right[3 * i + 0] - bound_left[3 * i + 0],
                        bound_right[3 * i + 1] - bound_left[3 * i + 1],
                        bound_right[3 * i + 2] - bound_left[3 * i + 2]};

        double d_cross[3] = {dv_dbeta[1] * dl[2] - dv_dbeta[2] * dl[1],
                             dv_dbeta[2] * dl[0] - dv_dbeta[0] * dl[2],
                             dv_dbeta[0] * dl[1] - dv_dbeta[1] * dl[0]};

        double df[3] = {2.0 * circulation[i] * d_cross[0],
                        2.0 * circulation[i] * d_cross[1],
                        2.0 * circulation[i] * d_cross[2]};

        d_cf_x += df[0];
        d_cf_y += df[1];
        d_cf_z += df[2];

        const double* r = &moment_arm[3 * i];
        d_cm_x += r[1] * df[2] - r[2] * df[1];
        d_cm_y += r[2] * df[0] - r[0] * df[2];
        d_cm_z += r[0] * df[1] - r[1] * df[0];
    }

    double inv_S = 1.0 / reference_area;
    d_cf_x *= inv_S;
    d_cf_y *= inv_S;
    d_cf_z *= inv_S;
    d_cm_x *= inv_S;
    d_cm_y *= inv_S;
    d_cm_z *= inv_S;

    // Body-to-stability: beta doesn't affect the rotation (only alpha does)
    d_result[0] = d_cf_z * ca - d_cf_x * sa;  // dCL/dbeta
    d_result[1] = d_cf_x * ca + d_cf_z * sa;  // dCD/dbeta
    d_result[2] = d_cf_y;
    d_result[3] = -d_cm_x;
    d_result[4] = d_cm_y;
    d_result[5] = -d_cm_z;
}

// ============================================================================
// Circulation Gradient (reverse mode)
// ============================================================================

void analytical_force_circulation_gradient(
    double alpha, double beta, const double* /*circulation*/,
    const double* induced_velocity, const double* bound_left,
    const double* bound_right, const double* /*control_points*/,
    const double* moment_arm, int num_vortices, double reference_area,
    int output_index, double* d_circulation) {
    // Force is linear in circulation: F_i = 2*Gamma_i*(V_eff_i x dl_i)
    // So dJ/dGamma_i = d(output)/d(Gamma_i) = (output sensitivity to Gamma_i)
    //
    // We need: d(output)/d(Gamma_i) which is the contribution of vortex i
    // to the output, divided by Gamma_i (since it's linear).

    double v_inf[3];
    compute_freestream_direction(alpha, beta, v_inf);
    BodyToOutputMap map(alpha);

    // Create seed: which output we're differentiating
    double d_out[6] = {0, 0, 0, 0, 0, 0};
    d_out[output_index] = 1.0;

    // Reverse through body-to-output mapping
    double d_body[6];
    map.reverse(d_out, d_body);

    double inv_S = 1.0 / reference_area;

    for (int i = 0; i < num_vortices; ++i) {
        double v_eff[3] = {v_inf[0] + induced_velocity[3 * i + 0],
                           v_inf[1] + induced_velocity[3 * i + 1],
                           v_inf[2] + induced_velocity[3 * i + 2]};

        double dl[3] = {bound_right[3 * i + 0] - bound_left[3 * i + 0],
                        bound_right[3 * i + 1] - bound_left[3 * i + 1],
                        bound_right[3 * i + 2] - bound_left[3 * i + 2]};

        // F_i = 2*Gamma_i * (V_eff x dl) / S_ref  (body axis)
        // dF_i/dGamma_i = 2*(V_eff x dl) / S_ref
        double cross[3] = {v_eff[1] * dl[2] - v_eff[2] * dl[1],
                           v_eff[2] * dl[0] - v_eff[0] * dl[2],
                           v_eff[0] * dl[1] - v_eff[1] * dl[0]};

        double dfi_dgi[3] = {2.0 * cross[0] * inv_S, 2.0 * cross[1] * inv_S,
                             2.0 * cross[2] * inv_S};

        // Moment contribution: dM_i/dGamma_i = r x dF_i/dGamma_i
        const double* r = &moment_arm[3 * i];
        double dmi[3] = {r[1] * dfi_dgi[2] - r[2] * dfi_dgi[1],
                         r[2] * dfi_dgi[0] - r[0] * dfi_dgi[2],
                         r[0] * dfi_dgi[1] - r[1] * dfi_dgi[0]};

        // Contract with the adjoint seed (d_body)
        d_circulation[i] = d_body[0] * dfi_dgi[0] + d_body[1] * dfi_dgi[1] +
                           d_body[2] * dfi_dgi[2] + d_body[3] * dmi[0] +
                           d_body[4] * dmi[1] + d_body[5] * dmi[2];
    }
}

// ============================================================================
// Induced Velocity Gradient (reverse mode)
// ============================================================================

void analytical_force_induced_velocity_gradient(
    double alpha, double beta, const double* circulation,
    const double* /*induced_velocity*/, const double* bound_left,
    const double* bound_right, const double* /*control_points*/,
    const double* moment_arm, int num_vortices, double reference_area,
    int output_index, double* d_induced_velocity) {
    // Force depends on V_induced through V_eff = V_inf + V_induced:
    // F_i = 2*Gamma_i * ((V_inf + V_ind_i) x dl_i) / S_ref
    // dF_i/dV_ind_i = 2*Gamma_i * (• x dl_i) / S_ref
    // This is a 3x3 matrix: dF_i[m]/dV_ind_i[k] = 2*Gamma_i * eps_{mkj}*dl_j /
    // S_ref
    //   = 2*Gamma_i * (-skew(dl))_{mk} / S_ref

    BodyToOutputMap map(alpha);
    (void)beta;

    double d_out[6] = {0, 0, 0, 0, 0, 0};
    d_out[output_index] = 1.0;

    double d_body[6];
    map.reverse(d_out, d_body);

    double inv_S = 1.0 / reference_area;

    for (int i = 0; i < num_vortices; ++i) {
        double dl[3] = {bound_right[3 * i + 0] - bound_left[3 * i + 0],
                        bound_right[3 * i + 1] - bound_left[3 * i + 1],
                        bound_right[3 * i + 2] - bound_left[3 * i + 2]};

        double g2 = 2.0 * circulation[i] * inv_S;

        // dF_i/dV_k: The force cross product (V x dl) differentiated w.r.t. V
        // gives: d(V x dl)/dV_k = e_k x dl (the k-th column is the cross
        // product of basis vector k with dl) So dF_i[m]/dV_k = g2 * (e_k x
        // dl)_m = g2 * eps_{mkj} dl_j

        // We need the adjoint: sum_m d_body_force_m * dF_i[m]/dV_k + moment
        // terms = g2 * sum_m d_body_force_m * eps_{mkj} dl_j = g2 *
        // (d_body_force x dl)_k? Let me verify: sum_m eps_{mkj} * s_m * dl_j
        // where s = d_body_force This is -(s x dl)_k? Let me be careful with
        // index: eps_{mkj} s_m dl_j = (we sum over m and j) This is the
        // k-component of (dl x s): (dl x s)_k = eps_{kij} dl_i s_j Hmm no. Let
        // me just compute component by component: k=0: sum_{m,j} eps_{m0j} s_m
        // dl_j = eps_{10j}s_1 dl_j + eps_{20j}s_2 dl_j
        //     = eps_{102}s_1 dl_2 + eps_{201}s_2 dl_1 = s_1*dl_2 - s_2*dl_1 =
        //     (s x dl)_0... wait (s x dl)_0 = s_1*dl_2 - s_2*dl_1. Yes!
        // So the force contribution is: d_V_k = g2 * (d_body_force x dl)_k...
        // let me double check. Actually I need to be more careful. eps_{mkj}
        // with free index k, summed over m and j. Let's use: sum_m eps_{mkj}
        // s_m = (component of cross product): For fixed k,j: sum_m eps_{mkj}
        // s_m is the (k,j) element of the matrix that gives s x • So: sum_{m,j}
        // eps_{mkj} s_m dl_j = (s x dl)_k? No: sum_j (sum_m eps_{mkj} s_m) dl_j
        // Hmm, eps_{mkj} = eps with indices (m,k,j). The cross product (a x
        // b)_i = eps_{ijk} a_j b_k. So sum_{m,j} eps_{mkj} s_m dl_j: let me
        // rename i=m, then sum_{i,j} eps_{ikj} s_i dl_j = -sum_{i,j} eps_{ijk}
        // s_i dl_j = -(s x dl)_k Hmm wait: eps_{ikj} = -eps_{ijk} (swap second
        // two indices). So it's -( s x dl )_k.
        //
        // OK so: force adjoint contribution = g2 * (-(d_body_force x dl))_k =
        // g2 * (dl x d_body_force)_k

        double df[3];  // d_body_force vector (just the force part)
        df[0] = d_body[0];
        df[1] = d_body[1];
        df[2] = d_body[2];

        // dl x df
        double dlxdf[3] = {dl[1] * df[2] - dl[2] * df[1],
                           dl[2] * df[0] - dl[0] * df[2],
                           dl[0] * df[1] - dl[1] * df[0]};

        // Moment contribution: M = r x F, so dM/dV = r x (dF/dV)
        // dM[m]/dV_k = sum_p eps_{mpq} r_p * dF_q/dV_k = r x (column k of
        // dF/dV) The adjoint: sum_m d_body_moment_m * dM[m]/dV_k = sum_m
        // d_body_moment_m * (r x (dF/dV column k))_m We know dF/dV column k
        // contribution to force adjoint is g2*(dl x •) applied to basis k
        // Actually let's just compute directly:
        // The total adjoint for V_ind_i is:
        //   d_V[k] = g2 * (dl x df)_k + moment_contribution_k

        // For moment: dM/dV_k = r x (2*Gamma/S * (e_k x dl))
        // Adjoint: sum_m dm_m * (r x (e_k x dl))_m where dm is d_body[3..5]
        // = dm . (r x (e_k x dl))
        // Using BAC-CAB: r x (e_k x dl) = e_k*(r.dl) - dl*(r.e_k)
        //                                 = (r.dl)*e_k - dl*(r_k)
        // So: dm . (r x (e_k x dl)) = (r.dl)*dm_k - (dm.dl)*r_k
        // Scaled by g2.

        const double* r = &moment_arm[3 * i];
        double dm[3] = {d_body[3], d_body[4], d_body[5]};
        double r_dot_dl = r[0] * dl[0] + r[1] * dl[1] + r[2] * dl[2];
        double dm_dot_dl = dm[0] * dl[0] + dm[1] * dl[1] + dm[2] * dl[2];

        d_induced_velocity[3 * i + 0] =
            g2 * (dlxdf[0] + r_dot_dl * dm[0] - dm_dot_dl * r[0]);
        d_induced_velocity[3 * i + 1] =
            g2 * (dlxdf[1] + r_dot_dl * dm[1] - dm_dot_dl * r[1]);
        d_induced_velocity[3 * i + 2] =
            g2 * (dlxdf[2] + r_dot_dl * dm[2] - dm_dot_dl * r[2]);
    }
}

// ============================================================================
// Geometry Gradient (forward mode)
// ============================================================================

void analytical_force_geometry_gradient(
    double alpha, double beta, const double* circulation,
    const double* induced_velocity, const double* bound_left,
    const double* bound_right, const double* moment_arm, int num_vortices,
    double reference_area, const double* d_bound_left,
    const double* d_bound_right, const double* d_moment_arm, double* d_result) {
    // Forward mode: seed mesh perturbations, get force sensitivity.
    // F_i = 2*Gamma_i * (V_eff x dl) / S_ref
    // dl = bound_right - bound_left
    // d_dl = d_bound_right - d_bound_left
    // dF_i = 2*Gamma_i * (V_eff x d_dl) / S_ref  (V_eff held fixed)
    // dM_i = d_r x F_i + r x dF_i

    double v_inf[3];
    compute_freestream_direction(alpha, beta, v_inf);

    double d_cf_x = 0.0, d_cf_y = 0.0, d_cf_z = 0.0;
    double d_cm_x = 0.0, d_cm_y = 0.0, d_cm_z = 0.0;

    double ca = std::cos(alpha);
    double sa = std::sin(alpha);

    for (int i = 0; i < num_vortices; ++i) {
        double v_eff[3] = {v_inf[0] + induced_velocity[3 * i + 0],
                           v_inf[1] + induced_velocity[3 * i + 1],
                           v_inf[2] + induced_velocity[3 * i + 2]};

        double dl[3] = {bound_right[3 * i + 0] - bound_left[3 * i + 0],
                        bound_right[3 * i + 1] - bound_left[3 * i + 1],
                        bound_right[3 * i + 2] - bound_left[3 * i + 2]};

        double d_dl[3] = {d_bound_right[3 * i + 0] - d_bound_left[3 * i + 0],
                          d_bound_right[3 * i + 1] - d_bound_left[3 * i + 1],
                          d_bound_right[3 * i + 2] - d_bound_left[3 * i + 2]};

        // Base force
        double cross[3] = {v_eff[1] * dl[2] - v_eff[2] * dl[1],
                           v_eff[2] * dl[0] - v_eff[0] * dl[2],
                           v_eff[0] * dl[1] - v_eff[1] * dl[0]};
        double f[3] = {2.0 * circulation[i] * cross[0],
                       2.0 * circulation[i] * cross[1],
                       2.0 * circulation[i] * cross[2]};

        // Force derivative from dl change
        double d_cross[3] = {v_eff[1] * d_dl[2] - v_eff[2] * d_dl[1],
                             v_eff[2] * d_dl[0] - v_eff[0] * d_dl[2],
                             v_eff[0] * d_dl[1] - v_eff[1] * d_dl[0]};
        double df[3] = {2.0 * circulation[i] * d_cross[0],
                        2.0 * circulation[i] * d_cross[1],
                        2.0 * circulation[i] * d_cross[2]};

        d_cf_x += df[0];
        d_cf_y += df[1];
        d_cf_z += df[2];

        // Moment: M = r x F, dM = d_r x F + r x dF
        const double* r = &moment_arm[3 * i];
        const double* d_r = &d_moment_arm[3 * i];

        // d_r x F
        d_cm_x += d_r[1] * f[2] - d_r[2] * f[1];
        d_cm_y += d_r[2] * f[0] - d_r[0] * f[2];
        d_cm_z += d_r[0] * f[1] - d_r[1] * f[0];

        // r x dF
        d_cm_x += r[1] * df[2] - r[2] * df[1];
        d_cm_y += r[2] * df[0] - r[0] * df[2];
        d_cm_z += r[0] * df[1] - r[1] * df[0];
    }

    double inv_S = 1.0 / reference_area;
    d_cf_x *= inv_S;
    d_cf_y *= inv_S;
    d_cf_z *= inv_S;
    d_cm_x *= inv_S;
    d_cm_y *= inv_S;
    d_cm_z *= inv_S;

    d_result[0] = d_cf_z * ca - d_cf_x * sa;  // dCL
    d_result[1] = d_cf_x * ca + d_cf_z * sa;  // dCD
    d_result[2] = d_cf_y;
    d_result[3] = -d_cm_x;
    d_result[4] = d_cm_y;
    d_result[5] = -d_cm_z;
}

}  // namespace ad
}  // namespace vlm
