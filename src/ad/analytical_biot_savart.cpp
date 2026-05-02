// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include "vlm/ad/analytical_biot_savart.h"

namespace vlm {
namespace ad {

// ============================================================================
// Internal helpers for computing derivatives of sub-expressions
// ============================================================================

namespace {

/// Compute the bound vortex contribution and its VJP w.r.t. A and B vectors.
///
/// The bound vortex velocity is:
///   V_bound = AxB * T * kPi4Inv
/// where:
///   T = T_num / T_den
///   T_num = (Bsq - AdB)/sqrt(sqrt(Bsq^2 + c4)) + (Asq - AdB)/sqrt(sqrt(Asq^2 + c4))
///   T_den = sqrt(AxBsq^2 + ALsq^2 * c4)
///   ALsq = Asq + Bsq - 2*AdB
///
/// The VJP computation: given seed s, compute s^T * dV/dA and s^T * dV/dB
/// Using product rule: V = AxB * T
///   dV/dA_k = d(AxB)/dA_k * T + AxB * dT/dA_k
///   s^T * dV/dA_k = s^T * d(AxB)/dA_k * T + (s . AxB) * dT/dA_k
///
/// d(AxB)/dA_k is: d(A x B)/dA = -skew(B), so (d(AxB)/dA_k)_m = -eps_{mk*}B_*
/// Actually: (AxB)_m = eps_{mij} A_i B_j
/// d(AxB)_m / dA_k = eps_{mkj} B_j = (e_k x B)_m = -(B x e_k)_m
/// So seed^T * d(AxB)/dA_k = (seed x B)_k... no wait.
/// Let's be careful:
///   seed . d(AxB)/dA_k = sum_m seed_m * eps_{mkj} B_j = (B x seed)_k
///
/// Similarly: seed . d(AxB)/dB_k = sum_m seed_m * eps_{mik} A_i...
/// d(AxB)_m / dB_k = eps_{mik} A_i... no.
/// (AxB)_m = eps_{mij} A_i B_j, d/dB_k = eps_{mik} A_i
/// seed . d(AxB)/dB_k = sum_m seed_m * eps_{mik} A_i = -(A x seed)_k...
/// Actually: sum_m eps_{mik} seed_m A_i = (seed x A)_k... let me verify:
/// eps_{mik} = -eps_{imk} = -eps_{imk}. Let me just use the formula:
///   d(AxB)/dA = -[B]_x  (skew matrix of B, applied on the right gives cross product)
///   d(AxB)/dB = [A]_x   (skew matrix of A)
/// So:
///   seed^T * d(AxB)/dA = -seed^T * [B]_x = -(seed x B)? No...
///   [B]_x * v = B x v, so v^T * [B]_x = (B x v)^T? No.
///   Actually [B]_x is skew-symmetric, so v^T * [B]_x = -(B x v)^T = (v x B)^T
///   Wait: [B]_x is such that [B]_x * w = B x w.
///   v^T * [B]_x = v^T * [B]_x. Since [B]_x is antisymmetric, v^T [B]_x = -([B]_x v)^T = -(B x v)^T
///   Hmm, let me just verify component by component.
///   (AxB)_0 = A_1 B_2 - A_2 B_1
///   d(AxB)_0/dA_0 = 0, d(AxB)_0/dA_1 = B_2, d(AxB)_0/dA_2 = -B_1
///   seed^T * d(AxB)/dA_k = s_0 * d(AxB)_0/dA_k + s_1 * d(AxB)_1/dA_k + s_2 * d(AxB)_2/dA_k
///   For k=0: s_0*0 + s_1*(-B_2) + s_2*B_1 = -(s_1*B_2 - s_2*B_1) = -(s x B)_0... wait no.
///   (s x B)_0 = s_1*B_2 - s_2*B_1. So for k=0: result = -(s x B)_0? Let me check:
///   s_1*d(AxB)_1/dA_0 + s_2*d(AxB)_2/dA_0.
///   (AxB)_1 = A_2 B_0 - A_0 B_2 => d/dA_0 = -B_2
///   (AxB)_2 = A_0 B_1 - A_1 B_0 => d/dA_0 = B_1
///   So: s_0*0 + s_1*(-B_2) + s_2*B_1 = s_2*B_1 - s_1*B_2 = -(s_1*B_2 - s_2*B_1) = -(sxB)_0
///   Hmm: (BxS)_0 = B_1*s_2 - B_2*s_1 = s_2*B_1 - s_1*B_2. That's exactly what we got.
///   So result_k = (B x s)_k for the A derivative of the cross product part.
///
/// Summary:
///   seed^T * d(AxB)/dA = B x seed
///   seed^T * d(AxB)/dB = seed x A  (= -(A x seed))
///
/// For the T derivative (scalar), dT/dA_k needs chain rule through T_num and T_den.
void bound_vortex_vjp(
    double Ax, double Ay, double Az,
    double Bx, double By, double Bz,
    double Asq, double Bsq,
    double core_4,
    const double* seed,  // [3] already scaled by kPi4Inv and PG factors
    double* vjp_A,       // [3] output: seed^T * dV_bound/dA
    double* vjp_B) {     // [3] output: seed^T * dV_bound/dB

  double AxBx = Ay * Bz - Az * By;
  double AxBy = Az * Bx - Ax * Bz;
  double AxBz = Ax * By - Ay * Bx;
  double AxBsq = AxBx*AxBx + AxBy*AxBy + AxBz*AxBz;

  if (AxBsq <= kBSEps) {
    vjp_A[0] = vjp_A[1] = vjp_A[2] = 0.0;
    vjp_B[0] = vjp_B[1] = vjp_B[2] = 0.0;
    return;
  }

  double AdB = Ax*Bx + Ay*By + Az*Bz;
  double ALsq = Asq + Bsq - 2.0 * AdB;

  // Compute T components
  double Asq_sq = Asq * Asq;
  double Bsq_sq = Bsq * Bsq;
  double inv_r4A = 1.0 / std::sqrt(std::sqrt(Asq_sq + core_4));  // (Asq^2 + c4)^{-1/4}
  double inv_r4B = 1.0 / std::sqrt(std::sqrt(Bsq_sq + core_4));  // (Bsq^2 + c4)^{-1/4}

  double numA = Asq - AdB;  // Contribution from A leg
  double numB = Bsq - AdB;  // Contribution from B leg
  double T_num = numB * inv_r4B + numA * inv_r4A;

  double den_arg = AxBsq * AxBsq + ALsq * ALsq * core_4;
  double inv_T_den = 1.0 / std::sqrt(den_arg);
  double T = T_num * inv_T_den;

  // Compute seed . AxB (scalar)
  double s_dot_AxB = seed[0]*AxBx + seed[1]*AxBy + seed[2]*AxBz;

  // --- Cross product part: seed^T * d(AxB)/dA_k * T ---
  // = (B x seed)_k * T
  double BxS[3], SxA[3];
  BxS[0] = By*seed[2] - Bz*seed[1];
  BxS[1] = Bz*seed[0] - Bx*seed[2];
  BxS[2] = Bx*seed[1] - By*seed[0];

  SxA[0] = seed[1]*Az - seed[2]*Ay;
  SxA[1] = seed[2]*Ax - seed[0]*Az;
  SxA[2] = seed[0]*Ay - seed[1]*Ax;

  // --- Scalar part: (seed . AxB) * dT/dA_k ---
  // T = T_num * inv_T_den
  // dT/dA_k = dT_num/dA_k * inv_T_den + T_num * d(inv_T_den)/dA_k

  // dT_num/dA:
  //   d(numA)/dA_k = d(Asq - AdB)/dA_k = 2*A_k - B_k
  //   d(numB)/dA_k = d(Bsq - AdB)/dA_k = -B_k  (Bsq doesn't depend on A)
  //   d(inv_r4A)/dA_k = d_fourth_root_inv(Asq^2, c4) * d(Asq^2)/dA_k...
  //   Actually inv_r4A = (Asq^2 + c4)^{-1/4}, and Asq = A.A
  //   d(inv_r4A)/dA_k = -1/4 * (Asq^2+c4)^{-5/4} * 2*Asq * 2*A_k
  //                    = -A_k * Asq / ((Asq^2+c4)^{5/4})
  //   Let's define: d_r4A_dAsq = -0.25 * (Asq^2+c4)^{-5/4} * 2*Asq = -0.5*Asq / (Asq^2+c4)^{5/4}
  //   Then d(inv_r4A)/dA_k = d_r4A_dAsq * 2*A_k

  // Actually let me be more careful. Let g = (x+c4)^{-1/4} where x = Asq^2.
  // dg/dx = -1/4 * (x+c4)^{-5/4}
  // dx/d(Asq) = 2*Asq
  // d(Asq)/dA_k = 2*A_k
  // So dg/dA_k = -1/4 * (Asq^2+c4)^{-5/4} * 2*Asq * 2*A_k = -Asq*A_k / (Asq^2+c4)^{5/4}

  double base_A5 = std::sqrt(Asq_sq + core_4);  // (Asq^2+c4)^{1/2}
  base_A5 *= base_A5;                            // (Asq^2+c4)
  double base_A5_full = base_A5 * std::sqrt(std::sqrt(Asq_sq + core_4));  // (Asq^2+c4)^{5/4}
  double d_inv_r4A_dAk_coeff = -Asq / base_A5_full;  // multiply by A_k to get d(inv_r4A)/dA_k

  double base_B5_full = std::sqrt(Bsq_sq + core_4);
  base_B5_full *= base_B5_full * std::sqrt(std::sqrt(Bsq_sq + core_4));  // (Bsq^2+c4)^{5/4}
  double d_inv_r4B_dBk_coeff = -Bsq / base_B5_full;

  // dT_num/dA_k = (2*A_k - B_k) * inv_r4A + numA * d_inv_r4A_dAk_coeff * A_k
  //             + (-B_k) * inv_r4B  [numB doesn't have A dependence on inv_r4B]
  // Wait: numB * inv_r4B + numA * inv_r4A
  // d/dA_k(numB * inv_r4B) = d(numB)/dA_k * inv_r4B = (-B_k) * inv_r4B
  // d/dA_k(numA * inv_r4A) = d(numA)/dA_k * inv_r4A + numA * d(inv_r4A)/dA_k
  //                         = (2*A_k - B_k) * inv_r4A + numA * d_inv_r4A_dAk_coeff * A_k

  // d(inv_T_den)/dA_k:
  // inv_T_den = (den_arg)^{-1/2}
  // d(inv_T_den)/dA_k = -0.5 * den_arg^{-3/2} * d(den_arg)/dA_k
  // den_arg = AxBsq^2 + ALsq^2 * c4
  // d(den_arg)/dA_k = 2*AxBsq * d(AxBsq)/dA_k + 2*ALsq * d(ALsq)/dA_k * c4
  //
  // d(AxBsq)/dA_k = 2*(AxBx * d(AxBx)/dA_k + AxBy * d(AxBy)/dA_k + AxBz * d(AxBz)/dA_k)
  //               = 2 * AxB . d(AxB)/dA_k = 2 * AxB . (e_k x B)...
  // Actually from above: d(AxB)_m/dA_k comes from the cross product Jacobian.
  // d(AxB)/dA is the matrix where column k gives d(AxB)/dA_k.
  // We showed: (seed^T * d(AxB)/dA)_k = (B x seed)_k
  // So: (AxB^T * d(AxB)/dA)_k = (B x AxB)_k
  // Thus d(AxBsq)/dA_k = 2 * (B x AxB)_k
  // Wait: d(AxBsq)/dA_k = 2 * sum_m AxB_m * d(AxB_m)/dA_k = 2*(AxB^T * Jac_A)_k = 2*(B x AxB)_k

  double BxAxB[3];
  BxAxB[0] = By*AxBz - Bz*AxBy;
  BxAxB[1] = Bz*AxBx - Bx*AxBz;
  BxAxB[2] = Bx*AxBy - By*AxBx;

  double AxBxA[3];
  AxBxA[0] = AxBy*Az - AxBz*Ay;
  AxBxA[1] = AxBz*Ax - AxBx*Az;
  AxBxA[2] = AxBx*Ay - AxBy*Ax;

  // d(ALsq)/dA_k = d(Asq + Bsq - 2*AdB)/dA_k = 2*A_k - 2*B_k
  // d(den_arg)/dA_k = 2*AxBsq * 2*(B x AxB)_k + 2*ALsq * (2*A_k - 2*B_k) * c4
  //                 = 4*AxBsq*(BxAxB)_k + 4*ALsq*(A_k - B_k)*c4

  double den_arg_inv32 = -0.5 * inv_T_den * inv_T_den * inv_T_den;  // -0.5 * den_arg^{-3/2}

  // Now assemble the full VJP for A (and similarly for B):
  for (int k = 0; k < 3; ++k) {
    double Ak = (k==0) ? Ax : ((k==1) ? Ay : Az);
    double Bk = (k==0) ? Bx : ((k==1) ? By : Bz);

    // dT_num/dA_k
    double dTnum_dAk = (-Bk) * inv_r4B +
                       (2.0*Ak - Bk) * inv_r4A +
                       numA * d_inv_r4A_dAk_coeff * Ak;

    // d(den_arg)/dA_k
    double dden_dAk = 4.0*AxBsq * BxAxB[k] + 4.0*ALsq*(Ak - Bk)*core_4;

    // d(inv_T_den)/dA_k
    double dinvTden_dAk = den_arg_inv32 * dden_dAk;

    // dT/dA_k
    double dT_dAk = dTnum_dAk * inv_T_den + T_num * dinvTden_dAk;

    // VJP for A_k: (B x seed)_k * T + s_dot_AxB * dT_dAk
    vjp_A[k] = BxS[k] * T + s_dot_AxB * dT_dAk;
  }

  // Same for B
  // d(numB)/dB_k = 2*B_k - A_k (from d(Bsq)/dB_k = 2*B_k, d(AdB)/dB_k = A_k)
  // d(numA)/dB_k = -A_k (from d(AdB)/dB_k = A_k)
  // d(inv_r4B)/dB_k = d_inv_r4B_dBk_coeff * B_k
  // d(inv_r4A)/dB_k = 0 (inv_r4A only depends on Asq)
  // d(AxBsq)/dB_k = 2 * (AxB x A)_k  (seed x A with seed = AxB)
  //   From: seed^T * d(AxB)/dB = seed x A (derived above with seed=AxB)
  //   So: d(AxBsq)/dB_k = 2*(AxB x A)_k = -2*(A x AxB)_k = 2*AxBxA[k]...
  //   Wait: (AxB x A)_k. Let me compute directly.
  //   (AxB x A)_0 = AxBy*Az - AxBz*Ay = AxBxA[0]? No:
  //   AxBxA was defined as AxBy*Az - AxBz*Ay etc. Yes that's (AxB x A).
  //   Actually let me redefine: let's call AxBcrossA = AxB x A.
  //   AxBcrossA[0] = AxBy*Az - AxBz*Ay
  //   Hmm that's exactly AxBxA. OK good.

  // d(ALsq)/dB_k = 2*B_k - 2*A_k

  for (int k = 0; k < 3; ++k) {
    double Ak = (k==0) ? Ax : ((k==1) ? Ay : Az);
    double Bk = (k==0) ? Bx : ((k==1) ? By : Bz);

    // dT_num/dB_k
    double dTnum_dBk = (2.0*Bk - Ak) * inv_r4B +
                       numB * d_inv_r4B_dBk_coeff * Bk +
                       (-Ak) * inv_r4A;

    // d(den_arg)/dB_k = 4*AxBsq*(AxBxA)_k + 4*ALsq*(B_k - A_k)*c4
    double dden_dBk = 4.0*AxBsq * AxBxA[k] + 4.0*ALsq*(Bk - Ak)*core_4;

    double dinvTden_dBk = den_arg_inv32 * dden_dBk;

    double dT_dBk = dTnum_dBk * inv_T_den + T_num * dinvTden_dBk;

    // VJP for B_k: (seed x A)_k * T + s_dot_AxB * dT_dBk
    vjp_B[k] = SxA[k] * T + s_dot_AxB * dT_dBk;
  }
}

/// Trailing leg VJP for leg at point P (A or B).
/// Trailing leg velocity from point P:
///   sign * (Pz * T, -Py * T) in (V, W) components
/// where T = sign_factor / sqrt(Pxi^2 + c4)  and  sign_factor depends on (1 - Px/|P|)
///
/// For leg A: sign = -1, velocity contribution to V += Az*T_A, W -= Ay*T_A
///   T_A = -(1 - Ax/Amag) / sqrt(Axisq^2 + c4)  where Axisq = Ay^2 + Az^2
///
/// VJP w.r.t. P (= A for left leg, = B for right leg):
///   The output is [0, Pz*T, -Py*T] (no x-component from trailing legs before PG scaling)
///   seed^T * dV/dP_k = seed[1] * d(Pz*T)/dP_k + seed[2] * d(-Py*T)/dP_k
void trailing_leg_vjp(
    double Px, double Py, double Pz,
    double, double Pmag,
    double core_4,
    double leg_sign,  // -1 for A leg, +1 for B leg
    const double* seed,  // [3] after kPi4Inv scaling
    double* vjp_P) {     // [3] output (accumulated)

  if (Pmag <= kBSEps) return;

  double Pxisq = Py*Py + Pz*Pz;
  double denom = std::sqrt(Pxisq * Pxisq + core_4);
  double factor = leg_sign * (1.0 - Px / Pmag);
  double T = factor / denom;  // Note: for A leg, leg_sign=-1 matches the -(1-Ax/Amag)

  // The trailing leg produces: dV = (0, Pz*T, -Py*T)
  // VJP seed:  s[1] * d(Pz*T)/dP_k + s[2] * d(-Py*T)/dP_k

  // We need dT/dP_k:
  // T = leg_sign * (1 - Px/Pmag) / sqrt(Pxisq^2 + c4)
  // Let F = leg_sign * (1 - Px/Pmag), G = 1/sqrt(Pxisq^2 + c4)
  // T = F * G
  // dT/dPx = dF/dPx * G + F * dG/dPx
  // dF/dPx = leg_sign * (-(Pmag - Px*(Px/Pmag))/Pmag^2) = leg_sign * (-(Psq - Px^2)/(Pmag^3))
  //        = -leg_sign * Pxisq / (Pmag^3)
  // dG/dPx = -0.5 * (Pxisq^2+c4)^{-3/2} * d(Pxisq^2)/dPx = 0 (Pxisq doesn't depend on Px)
  // So dT/dPx = -leg_sign * Pxisq / (Pmag^3 * denom)

  double dT_dPx = -leg_sign * Pxisq / (Pmag * Pmag * Pmag * denom);

  // dT/dPy:
  // dF/dPy = leg_sign * (Px * Py / Pmag^3) = leg_sign * Px*Py/Pmag^3
  // dG/dPy = -0.5 * (Pxisq^2+c4)^{-3/2} * 2*Pxisq * 2*Py = -2*Pxisq*Py / (denom^3)
  double denom3 = denom * denom * denom;
  double dF_dPy = leg_sign * Px * Py / (Pmag * Pmag * Pmag);
  double dG_dPy = -2.0 * Pxisq * Py / denom3;
  double G = 1.0 / denom;
  double F = factor;
  double dT_dPy = dF_dPy * G + F * dG_dPy;

  // dT/dPz: same structure as dT/dPy with Pz
  double dF_dPz = leg_sign * Px * Pz / (Pmag * Pmag * Pmag);
  double dG_dPz = -2.0 * Pxisq * Pz / denom3;
  double dT_dPz = dF_dPz * G + F * dG_dPz;

  // VJP accumulation:
  // V_trail = Pz * T, W_trail = -Py * T
  // d(V_trail)/dPx = Pz * dT/dPx
  // d(V_trail)/dPy = Pz * dT/dPy
  // d(V_trail)/dPz = T + Pz * dT/dPz
  // d(W_trail)/dPx = -Py * dT/dPx
  // d(W_trail)/dPy = -T - Py * dT/dPy
  // d(W_trail)/dPz = -Py * dT/dPz

  // Result_k = seed[1] * d(V_trail)/dP_k + seed[2] * d(W_trail)/dP_k
  vjp_P[0] += seed[1] * Pz * dT_dPx + seed[2] * (-Py * dT_dPx);
  vjp_P[1] += seed[1] * Pz * dT_dPy + seed[2] * (-T - Py * dT_dPy);
  vjp_P[2] += seed[1] * (T + Pz * dT_dPz) + seed[2] * (-Py * dT_dPz);
}

}  // anonymous namespace

// ============================================================================
// Public VJP Functions
// ============================================================================

void horseshoe_velocity_vjp_field_point(
    const double* field_point,
    const double* bound_left,
    const double* bound_right,
    double prandtl_glauert,
    double core_radius,
    bool include_bound,
    const double* seed,
    double* result) {

  result[0] = 0.0;
  result[1] = 0.0;
  result[2] = 0.0;

  double pg_inv = 1.0 / prandtl_glauert;

  // Transform to PG coordinates
  double Ax = (bound_left[0] - field_point[0]) * pg_inv;
  double Ay = bound_left[1] - field_point[1];
  double Az = bound_left[2] - field_point[2];
  double Bx = (bound_right[0] - field_point[0]) * pg_inv;
  double By = bound_right[1] - field_point[1];
  double Bz = bound_right[2] - field_point[2];

  double Asq = Ax*Ax + Ay*Ay + Az*Az;
  double Bsq = Bx*Bx + By*By + Bz*Bz;
  double Amag = std::sqrt(Asq);
  double Bmag = std::sqrt(Bsq);

  double core_sq = core_radius * core_radius;
  double core_4 = core_sq * core_sq;

  // The output velocity is:
  //   V[0] = U * kPi4Inv / prandtl_glauert = U * kPi4Inv * pg_inv
  //   V[1] = V_comp * kPi4Inv
  //   V[2] = W * kPi4Inv
  // So the effective seed for the (U, V_comp, W) quantities is:
  double eff_seed[3] = {
    seed[0] * kPi4Inv * pg_inv,
    seed[1] * kPi4Inv,
    seed[2] * kPi4Inv
  };

  // A = bound_left - field_point (in PG coords for x)
  // dA/d(field_point) = [-pg_inv, -1, -1] (diagonal)
  // So vjp_field_point[k] = vjp_A[k] * dA_k/d(fp_k) + vjp_B[k] * dB_k/d(fp_k)
  //                       = -(vjp_A[k] + vjp_B[k]) * scale_k
  // where scale = [pg_inv, 1, 1]

  // Bound vortex contribution
  if (include_bound && Amag * Bmag > kBSEps) {
    double vjp_A[3] = {0, 0, 0};
    double vjp_B[3] = {0, 0, 0};
    bound_vortex_vjp(Ax, Ay, Az, Bx, By, Bz, Asq, Bsq, core_4, eff_seed, vjp_A, vjp_B);

    // Chain rule: dV/d(fp) through A and B
    result[0] -= (vjp_A[0] + vjp_B[0]) * pg_inv;
    result[1] -= (vjp_A[1] + vjp_B[1]);
    result[2] -= (vjp_A[2] + vjp_B[2]);
  }

  // Trailing leg A: A depends on field_point: dA/d(fp) = [-pg_inv, -1, -1]
  if (Amag > kBSEps) {
    double vjp_A_trail[3] = {0, 0, 0};
    trailing_leg_vjp(Ax, Ay, Az, Asq, Amag, core_4, -1.0, eff_seed, vjp_A_trail);
    result[0] -= vjp_A_trail[0] * pg_inv;
    result[1] -= vjp_A_trail[1];
    result[2] -= vjp_A_trail[2];
  }

  // Trailing leg B: B depends on field_point: dB/d(fp) = [-pg_inv, -1, -1]
  if (Bmag > kBSEps) {
    double vjp_B_trail[3] = {0, 0, 0};
    trailing_leg_vjp(Bx, By, Bz, Bsq, Bmag, core_4, 1.0, eff_seed, vjp_B_trail);
    result[0] -= vjp_B_trail[0] * pg_inv;
    result[1] -= vjp_B_trail[1];
    result[2] -= vjp_B_trail[2];
  }
}

void horseshoe_velocity_vjp_bound_left(
    const double* field_point,
    const double* bound_left,
    const double* bound_right,
    double prandtl_glauert,
    double core_radius,
    bool include_bound,
    const double* seed,
    double* result) {

  result[0] = 0.0;
  result[1] = 0.0;
  result[2] = 0.0;

  double pg_inv = 1.0 / prandtl_glauert;

  double Ax = (bound_left[0] - field_point[0]) * pg_inv;
  double Ay = bound_left[1] - field_point[1];
  double Az = bound_left[2] - field_point[2];
  double Bx = (bound_right[0] - field_point[0]) * pg_inv;
  double By = bound_right[1] - field_point[1];
  double Bz = bound_right[2] - field_point[2];

  double Asq = Ax*Ax + Ay*Ay + Az*Az;
  double Bsq = Bx*Bx + By*By + Bz*Bz;
  double Amag = std::sqrt(Asq);
  double Bmag = std::sqrt(Bsq);

  double core_sq = core_radius * core_radius;
  double core_4 = core_sq * core_sq;

  double eff_seed[3] = {
    seed[0] * kPi4Inv * pg_inv,
    seed[1] * kPi4Inv,
    seed[2] * kPi4Inv
  };

  // A = bound_left - field_point (PG x)
  // dA/d(bound_left) = [pg_inv, 1, 1] (diagonal)
  // B doesn't depend on bound_left

  // Bound vortex
  if (include_bound && Amag * Bmag > kBSEps) {
    double vjp_A[3] = {0, 0, 0};
    double vjp_B[3] = {0, 0, 0};
    bound_vortex_vjp(Ax, Ay, Az, Bx, By, Bz, Asq, Bsq, core_4, eff_seed, vjp_A, vjp_B);

    result[0] += vjp_A[0] * pg_inv;
    result[1] += vjp_A[1];
    result[2] += vjp_A[2];
  }

  // Trailing leg A depends on bound_left
  if (Amag > kBSEps) {
    double vjp_A_trail[3] = {0, 0, 0};
    trailing_leg_vjp(Ax, Ay, Az, Asq, Amag, core_4, -1.0, eff_seed, vjp_A_trail);
    result[0] += vjp_A_trail[0] * pg_inv;
    result[1] += vjp_A_trail[1];
    result[2] += vjp_A_trail[2];
  }

  // Trailing leg B does NOT depend on bound_left
}

void horseshoe_velocity_vjp_bound_right(
    const double* field_point,
    const double* bound_left,
    const double* bound_right,
    double prandtl_glauert,
    double core_radius,
    bool include_bound,
    const double* seed,
    double* result) {

  result[0] = 0.0;
  result[1] = 0.0;
  result[2] = 0.0;

  double pg_inv = 1.0 / prandtl_glauert;

  double Ax = (bound_left[0] - field_point[0]) * pg_inv;
  double Ay = bound_left[1] - field_point[1];
  double Az = bound_left[2] - field_point[2];
  double Bx = (bound_right[0] - field_point[0]) * pg_inv;
  double By = bound_right[1] - field_point[1];
  double Bz = bound_right[2] - field_point[2];

  double Asq = Ax*Ax + Ay*Ay + Az*Az;
  double Bsq = Bx*Bx + By*By + Bz*Bz;
  double Amag = std::sqrt(Asq);
  double Bmag = std::sqrt(Bsq);

  double core_sq = core_radius * core_radius;
  double core_4 = core_sq * core_sq;

  double eff_seed[3] = {
    seed[0] * kPi4Inv * pg_inv,
    seed[1] * kPi4Inv,
    seed[2] * kPi4Inv
  };

  // B = bound_right - field_point (PG x)
  // dB/d(bound_right) = [pg_inv, 1, 1]
  // A doesn't depend on bound_right

  // Bound vortex
  if (include_bound && Amag * Bmag > kBSEps) {
    double vjp_A[3] = {0, 0, 0};
    double vjp_B[3] = {0, 0, 0};
    bound_vortex_vjp(Ax, Ay, Az, Bx, By, Bz, Asq, Bsq, core_4, eff_seed, vjp_A, vjp_B);

    result[0] += vjp_B[0] * pg_inv;
    result[1] += vjp_B[1];
    result[2] += vjp_B[2];
  }

  // Trailing leg B depends on bound_right
  if (Bmag > kBSEps) {
    double vjp_B_trail[3] = {0, 0, 0};
    trailing_leg_vjp(Bx, By, Bz, Bsq, Bmag, core_4, 1.0, eff_seed, vjp_B_trail);
    result[0] += vjp_B_trail[0] * pg_inv;
    result[1] += vjp_B_trail[1];
    result[2] += vjp_B_trail[2];
  }

  // Trailing leg A does NOT depend on bound_right
}

}  // namespace ad
}  // namespace vlm
