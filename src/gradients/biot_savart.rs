/*
    Project Zephyr

    Copyright 2024-2026, Christopher A. Lupp

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This work has been cleared for public release, distribution unlimited, case
    number: AFRL-2026-2198.

    The views expressed are those of the authors and do not reflect the
    official guidance or position of the United States Government, the
    Department of Defense or of the United States Air Force.

    Statement from DoD: The Appearance of external hyperlinks does not
    constitute endorsement by the United States Department of Defense (DoD) of
    the linked websites, of the information, products, or services contained
    therein. The DoD does not exercise any editorial, security, or other
    control over the information you may find at these locations.
*/

//! Exact VJPs of the regularized horseshoe-vortex kernel.
//!
//! The derivatives are hand-derived closed forms, matching the structure of the
//! v0.2.0 oracle's `analytical_biot_savart.cpp`. VJP rather than JVP form
//! because every caller contracts with a seed vector — the panel normal for the
//! AIC, a force direction for the loads — so the seed can be pushed through the
//! chain instead of carrying a full 3x3 Jacobian.
//!
//! The primal is [`horseshoe_vortex_velocity`]:
//!
//! ```text
//! V = (a x b) * T / (4 pi)   +   trailing legs
//! T = [ (b.b - a.b) / (|b|^4 + c^4)^(1/4)
//!     + (a.a - a.b) / (|a|^4 + c^4)^(1/4) ] / sqrt(|a x b|^4 + L^4 c^4)
//! ```
//!
//! with `a`, `b` the vectors from the field point to the bound endpoints in
//! Prandtl-Glauert coordinates, `L^2 = a.a + b.b - 2 a.b`, and `c` the Leishman
//! R^4 core radius.
//!
//! # Degeneracy guards
//!
//! The guards below are the primal's, not the oracle's. `horseshoe_vortex_velocity`
//! replaced the C++ absolute thresholds (`|a||b| > 1e-12`, `|a x b|^2 > 1e-12`)
//! with `a_sq > 0`, `b_sq > 0` and the dimensionless `cross_sq > EPSILON * a_sq *
//! b_sq`, because comparing a length^4 quantity against an absolute epsilon made
//! the kernel depend on the unit the model is expressed in. A VJP that kept the
//! absolute thresholds would drop the bound term below roughly millimetre scale
//! while the primal still included it, so the derivative would silently stop
//! describing the function it differentiates. These branches must therefore
//! track the primal's exactly.

use crate::core::{EPSILON, Mat3, PI, Scalar, Vec3};
use crate::solver::horseshoe_vortex_velocity;

const INV_4_PI: Scalar = 1.0 / (4.0 * PI);

/// The horseshoe velocity's three geometry Jacobians.
///
/// Each block is `d(velocity)/d(point)`, so entry `(m, k)` is the derivative of
/// velocity component `m` with respect to component `k` of that point. Left
/// multiplying by a seed's transpose recovers the matching VJP, which is why
/// `Mat3::tr_mul` is how callers contract these.
#[derive(Clone, Copy, Debug, Default)]
pub struct HorseshoeVelocityJacobians {
    /// Derivative with respect to the velocity evaluation point.
    pub field_point: Mat3,
    /// Derivative with respect to the bound-vortex left endpoint.
    pub bound_left: Mat3,
    /// Derivative with respect to the bound-vortex right endpoint.
    pub bound_right: Mat3,
}

/// All three geometry Jacobians of the horseshoe velocity in one pass.
///
/// The VJP kernels are the right shape when a caller has one seed per field
/// point, which is the case for a single output. A reverse sweep carrying `k`
/// simultaneous seeds — one per column of a coupled Jacobian — would otherwise
/// evaluate the kernel `k` times per vortex pair for the same geometry. This
/// builds the dense blocks once, from the same bound and trailing-leg
/// derivatives the VJPs use, so the cost per pair stops growing with the number
/// of seeds.
#[must_use]
pub fn horseshoe_velocity_jacobians(
    field_point: &Vec3,
    bound_left: &Vec3,
    bound_right: &Vec3,
    prandtl_glauert_factor: Scalar,
    core_radius: Scalar,
    include_bound: bool,
) -> HorseshoeVelocityJacobians {
    let geometry = Geometry::new(field_point, bound_left, bound_right, prandtl_glauert_factor);
    let core_4 = core_four(core_radius);
    let bound_active = include_bound && geometry.bound_is_active();

    let mut result = HorseshoeVelocityJacobians::default();
    for component in 0..3 {
        let mut basis = Vec3::zeros();
        basis[component] = 1.0;
        let seed = effective_seed(&basis, prandtl_glauert_factor);

        let (mut left, mut right) = (Vec3::zeros(), Vec3::zeros());
        if bound_active {
            let (vjp_a, vjp_b) = geometry.bound_vjp(core_4, &seed);
            left += scale_x(vjp_a, prandtl_glauert_factor);
            right += scale_x(vjp_b, prandtl_glauert_factor);
        }
        if geometry.a_sq > 0.0 {
            let vjp = trailing_leg_vjp(&geometry.a, geometry.a_mag, core_4, -1.0, &seed);
            left += scale_x(vjp, prandtl_glauert_factor);
        }
        if geometry.b_sq > 0.0 {
            let vjp = trailing_leg_vjp(&geometry.b, geometry.b_mag, core_4, 1.0, &seed);
            right += scale_x(vjp, prandtl_glauert_factor);
        }

        // `a` and `b` are both measured from the field point, so the field-point
        // row is minus the sum of the two endpoint rows, exactly as
        // `horseshoe_velocity_vjp_field_point` assembles it.
        result.bound_left.set_row(component, &left.transpose());
        result.bound_right.set_row(component, &right.transpose());
        result
            .field_point
            .set_row(component, &(-(left + right)).transpose());
    }
    result
}

/// VJP with respect to the velocity evaluation point.
#[must_use]
pub fn horseshoe_velocity_vjp_field_point(
    field_point: &Vec3,
    bound_left: &Vec3,
    bound_right: &Vec3,
    prandtl_glauert_factor: Scalar,
    core_radius: Scalar,
    include_bound: bool,
    velocity_cotangent: &Vec3,
) -> Vec3 {
    let geometry = Geometry::new(field_point, bound_left, bound_right, prandtl_glauert_factor);
    let seed = effective_seed(velocity_cotangent, prandtl_glauert_factor);
    let core_4 = core_four(core_radius);

    // `a` and `b` are both measured from the field point, so both legs and the
    // bound contribution respond to it: da/d(fp) = db/d(fp) = -diag(1/pg, 1, 1).
    let mut result = Vec3::zeros();
    if include_bound && geometry.bound_is_active() {
        let (vjp_a, vjp_b) = geometry.bound_vjp(core_4, &seed);
        result -= scale_x(vjp_a + vjp_b, prandtl_glauert_factor);
    }
    if geometry.a_sq > 0.0 {
        let vjp = trailing_leg_vjp(&geometry.a, geometry.a_mag, core_4, -1.0, &seed);
        result -= scale_x(vjp, prandtl_glauert_factor);
    }
    if geometry.b_sq > 0.0 {
        let vjp = trailing_leg_vjp(&geometry.b, geometry.b_mag, core_4, 1.0, &seed);
        result -= scale_x(vjp, prandtl_glauert_factor);
    }
    result
}

/// VJP with respect to the bound-vortex left endpoint.
#[must_use]
pub fn horseshoe_velocity_vjp_bound_left(
    field_point: &Vec3,
    bound_left: &Vec3,
    bound_right: &Vec3,
    prandtl_glauert_factor: Scalar,
    core_radius: Scalar,
    include_bound: bool,
    velocity_cotangent: &Vec3,
) -> Vec3 {
    let geometry = Geometry::new(field_point, bound_left, bound_right, prandtl_glauert_factor);
    let seed = effective_seed(velocity_cotangent, prandtl_glauert_factor);
    let core_4 = core_four(core_radius);

    // Only `a` depends on the left endpoint: da/d(bl) = diag(1/pg, 1, 1), and
    // db/d(bl) = 0, so the right trailing leg contributes nothing here.
    let mut result = Vec3::zeros();
    if include_bound && geometry.bound_is_active() {
        let (vjp_a, _) = geometry.bound_vjp(core_4, &seed);
        result += scale_x(vjp_a, prandtl_glauert_factor);
    }
    if geometry.a_sq > 0.0 {
        let vjp = trailing_leg_vjp(&geometry.a, geometry.a_mag, core_4, -1.0, &seed);
        result += scale_x(vjp, prandtl_glauert_factor);
    }
    result
}

/// VJP with respect to the bound-vortex right endpoint.
#[must_use]
pub fn horseshoe_velocity_vjp_bound_right(
    field_point: &Vec3,
    bound_left: &Vec3,
    bound_right: &Vec3,
    prandtl_glauert_factor: Scalar,
    core_radius: Scalar,
    include_bound: bool,
    velocity_cotangent: &Vec3,
) -> Vec3 {
    let geometry = Geometry::new(field_point, bound_left, bound_right, prandtl_glauert_factor);
    let seed = effective_seed(velocity_cotangent, prandtl_glauert_factor);
    let core_4 = core_four(core_radius);

    // Mirror of the left case: db/d(br) = diag(1/pg, 1, 1), da/d(br) = 0.
    let mut result = Vec3::zeros();
    if include_bound && geometry.bound_is_active() {
        let (_, vjp_b) = geometry.bound_vjp(core_4, &seed);
        result += scale_x(vjp_b, prandtl_glauert_factor);
    }
    if geometry.b_sq > 0.0 {
        let vjp = trailing_leg_vjp(&geometry.b, geometry.b_mag, core_4, 1.0, &seed);
        result += scale_x(vjp, prandtl_glauert_factor);
    }
    result
}

/// Directional derivative of the horseshoe velocity with respect to the
/// Prandtl–Glauert factor.
#[must_use]
pub fn horseshoe_velocity_pg_derivative(
    field_point: &Vec3,
    bound_left: &Vec3,
    bound_right: &Vec3,
    prandtl_glauert_factor: Scalar,
    core_radius: Scalar,
    include_bound: bool,
) -> Vec3 {
    if prandtl_glauert_factor.abs() <= EPSILON {
        return Vec3::zeros();
    }
    let left_direction = Vec3::new(
        -(bound_left.x - field_point.x) / prandtl_glauert_factor,
        0.0,
        0.0,
    );
    let right_direction = Vec3::new(
        -(bound_right.x - field_point.x) / prandtl_glauert_factor,
        0.0,
        0.0,
    );
    let velocity = horseshoe_vortex_velocity(
        field_point,
        bound_left,
        bound_right,
        prandtl_glauert_factor,
        core_radius,
        include_bound,
    );
    Vec3::from_fn(|component, _| {
        let mut seed = Vec3::zeros();
        seed[component] = 1.0;
        let endpoint_derivative = horseshoe_velocity_vjp_bound_left(
            field_point,
            bound_left,
            bound_right,
            prandtl_glauert_factor,
            core_radius,
            include_bound,
            &seed,
        )
        .dot(&left_direction)
            + horseshoe_velocity_vjp_bound_right(
                field_point,
                bound_left,
                bound_right,
                prandtl_glauert_factor,
                core_radius,
                include_bound,
                &seed,
            )
            .dot(&right_direction);
        endpoint_derivative
            - if component == 0 {
                velocity.x / prandtl_glauert_factor
            } else {
                0.0
            }
    })
}

fn core_four(core_radius: Scalar) -> Scalar {
    let core_sq = core_radius * core_radius;
    core_sq * core_sq
}

/// The cotangent of the unscaled `(u, v, w)` accumulator.
///
/// The primal returns `(u / (4 pi pg), v / (4 pi), w / (4 pi))`, so folding the
/// constant factors into the seed once keeps them out of every downstream term.
fn effective_seed(velocity_cotangent: &Vec3, prandtl_glauert_factor: Scalar) -> Vec3 {
    Vec3::new(
        velocity_cotangent.x * INV_4_PI / prandtl_glauert_factor,
        velocity_cotangent.y * INV_4_PI,
        velocity_cotangent.z * INV_4_PI,
    )
}

/// Chain through the Prandtl-Glauert stretch, which scales only the x row.
fn scale_x(vector: Vec3, prandtl_glauert_factor: Scalar) -> Vec3 {
    Vec3::new(vector.x / prandtl_glauert_factor, vector.y, vector.z)
}

/// Field-point-relative endpoint vectors in Prandtl-Glauert coordinates, with
/// the norms every branch needs.
struct Geometry {
    a: Vec3,
    b: Vec3,
    a_sq: Scalar,
    b_sq: Scalar,
    a_mag: Scalar,
    b_mag: Scalar,
    cross_ab: Vec3,
    cross_sq: Scalar,
}

impl Geometry {
    fn new(
        field_point: &Vec3,
        bound_left: &Vec3,
        bound_right: &Vec3,
        prandtl_glauert_factor: Scalar,
    ) -> Self {
        let a = Vec3::new(
            (bound_left.x - field_point.x) / prandtl_glauert_factor,
            bound_left.y - field_point.y,
            bound_left.z - field_point.z,
        );
        let b = Vec3::new(
            (bound_right.x - field_point.x) / prandtl_glauert_factor,
            bound_right.y - field_point.y,
            bound_right.z - field_point.z,
        );
        let a_sq = a.norm_squared();
        let b_sq = b.norm_squared();
        let cross_ab = a.cross(&b);
        Self {
            a,
            b,
            a_sq,
            b_sq,
            a_mag: a_sq.sqrt(),
            b_mag: b_sq.sqrt(),
            cross_ab,
            cross_sq: cross_ab.norm_squared(),
        }
    }

    /// The primal's bound-leg branch, reproduced exactly. `cross_sq > EPSILON *
    /// a_sq * b_sq` is `sin^2(theta) > EPSILON`, invariant under a uniform
    /// rescaling of the geometry.
    fn bound_is_active(&self) -> bool {
        self.a_sq > 0.0 && self.b_sq > 0.0 && self.cross_sq > EPSILON * self.a_sq * self.b_sq
    }

    /// Seed-contracted derivatives of the bound contribution `(a x b) * T` with
    /// respect to `a` and `b`.
    ///
    /// Product rule on `V = (a x b) T`:
    ///
    /// ```text
    /// s . dV/da_k = (b x s)_k T + (s . (a x b)) dT/da_k
    /// s . dV/db_k = (s x a)_k T + (s . (a x b)) dT/db_k
    /// ```
    ///
    /// The cross-product Jacobians follow from `(a x b)_m = eps_mij a_i b_j`:
    /// `s^T d(a x b)/da = b x s` and `s^T d(a x b)/db = s x a`.
    fn bound_vjp(&self, core_4: Scalar, seed: &Vec3) -> (Vec3, Vec3) {
        let dot_ab = self.a.dot(&self.b);
        let length_sq = self.a_sq + self.b_sq - 2.0 * dot_ab;

        let a_sq_sq = self.a_sq * self.a_sq;
        let b_sq_sq = self.b_sq * self.b_sq;
        let root_a = (a_sq_sq + core_4).sqrt().sqrt();
        let root_b = (b_sq_sq + core_4).sqrt().sqrt();
        let inv_r4_a = 1.0 / root_a;
        let inv_r4_b = 1.0 / root_b;

        let num_a = self.a_sq - dot_ab;
        let num_b = self.b_sq - dot_ab;
        let t_num = num_b * inv_r4_b + num_a * inv_r4_a;

        let den_arg = self.cross_sq * self.cross_sq + length_sq * length_sq * core_4;
        let inv_t_den = 1.0 / den_arg.sqrt();
        let t = t_num * inv_t_den;

        // d/da_k of (a.a^2 + c^4)^(-1/4) is -a.a a_k / (a.a^2 + c^4)^(5/4).
        let d_inv_r4_a = -self.a_sq / ((a_sq_sq + core_4) * root_a);
        let d_inv_r4_b = -self.b_sq / ((b_sq_sq + core_4) * root_b);

        let seed_dot_cross = seed.dot(&self.cross_ab);
        let b_cross_seed = self.b.cross(seed);
        let seed_cross_a = seed.cross(&self.a);

        // d(|a x b|^2)/da_k = 2 (b x (a x b))_k, the same Jacobian identity
        // applied with the cross product itself as the seed.
        let b_cross_ab = self.b.cross(&self.cross_ab);
        let ab_cross_a = self.cross_ab.cross(&self.a);

        // d(den_arg)^(-1/2)/dx = -1/2 den_arg^(-3/2) d(den_arg)/dx
        let d_inv_t_den = -0.5 * inv_t_den * inv_t_den * inv_t_den;

        let mut vjp_a = Vec3::zeros();
        let mut vjp_b = Vec3::zeros();
        for k in 0..3 {
            let a_k = self.a[k];
            let b_k = self.b[k];

            let d_t_num_a =
                (-b_k) * inv_r4_b + (2.0 * a_k - b_k) * inv_r4_a + num_a * d_inv_r4_a * a_k;
            let d_den_a =
                4.0 * self.cross_sq * b_cross_ab[k] + 4.0 * length_sq * (a_k - b_k) * core_4;
            let d_t_a = d_t_num_a * inv_t_den + t_num * d_inv_t_den * d_den_a;
            vjp_a[k] = b_cross_seed[k] * t + seed_dot_cross * d_t_a;

            let d_t_num_b =
                (2.0 * b_k - a_k) * inv_r4_b + num_b * d_inv_r4_b * b_k + (-a_k) * inv_r4_a;
            let d_den_b =
                4.0 * self.cross_sq * ab_cross_a[k] + 4.0 * length_sq * (b_k - a_k) * core_4;
            let d_t_b = d_t_num_b * inv_t_den + t_num * d_inv_t_den * d_den_b;
            vjp_b[k] = seed_cross_a[k] * t + seed_dot_cross * d_t_b;
        }
        (vjp_a, vjp_b)
    }
}

/// Seed-contracted derivative of one semi-infinite trailing leg.
///
/// The leg from endpoint `p` contributes `(0, p_z T, -p_y T)` with
/// `T = sign (1 - p_x / |p|) / sqrt((p_y^2 + p_z^2)^2 + c^4)`, so `sign` is
/// `-1` for the left leg and `+1` for the right. Only the y and z rows of the
/// seed are reachable: a trailing leg induces no x velocity before the
/// Prandtl-Glauert stretch.
fn trailing_leg_vjp(p: &Vec3, p_mag: Scalar, core_4: Scalar, sign: Scalar, seed: &Vec3) -> Vec3 {
    let axis_sq = p.y * p.y + p.z * p.z;
    let denom = (axis_sq * axis_sq + core_4).sqrt();
    let factor = sign * (1.0 - p.x / p_mag);
    let t = factor / denom;

    let inv_denom = 1.0 / denom;
    let denom_cubed = denom * denom * denom;
    let p_mag_cubed = p_mag * p_mag * p_mag;

    // d(1 - p_x/|p|)/dp_x = -(|p|^2 - p_x^2)/|p|^3 = -axis_sq/|p|^3, and the
    // denominator has no p_x dependence at all.
    let d_t_x = -sign * axis_sq / (p_mag_cubed * denom);
    let d_t_y =
        sign * p.x * p.y / p_mag_cubed * inv_denom + factor * (-2.0 * axis_sq * p.y / denom_cubed);
    let d_t_z =
        sign * p.x * p.z / p_mag_cubed * inv_denom + factor * (-2.0 * axis_sq * p.z / denom_cubed);

    // v = p_z T and w = -p_y T, so the p_z and p_y rows pick up the extra term
    // from differentiating the prefactor rather than T.
    Vec3::new(
        seed.y * p.z * d_t_x + seed.z * (-p.y * d_t_x),
        seed.y * p.z * d_t_y + seed.z * (-t - p.y * d_t_y),
        seed.y * (t + p.z * d_t_z) + seed.z * (-p.y * d_t_z),
    )
}
