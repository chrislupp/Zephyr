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
//! The primal is [`horseshoe_vortex_velocity`](crate::solver::horseshoe_vortex_velocity):
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

//!
//! # Shared per-pair terms
//!
//! Everything above except the seed is a property of the field point and the
//! horseshoe: `T`, its gradients with respect to `a` and `b`, and the trailing
//! legs' scalar factors and their gradients. `HorseshoeDerivatives` computes
//! those once per pair and then contracts them with as many seeds as a caller
//! needs, so a full Jacobian, a single-seed VJP and the Prandtl-Glauert
//! derivative all reuse one evaluation of the geometry. The public functions in
//! this module are built on the same pieces and return the same bits they
//! always did; the Prandtl-Glauert derivative may differ only in the sign of an
//! entry that is exactly zero.

use crate::core::{EPSILON, Mat3, PI, Scalar, Vec3};

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

impl HorseshoeVelocityJacobians {
    /// The seed-independent terms of one (field point, horseshoe) pair, for
    /// crate code outside this module, which is private to `gradients`.
    /// Arguments are those of
    /// [`horseshoe_vortex_velocity`](crate::solver::horseshoe_vortex_velocity);
    /// the result offers single-seed endpoint VJPs, the Jacobians and the
    /// Prandtl-Glauert derivative, each bit-identical to the public kernels.
    #[inline]
    pub(crate) fn shared_terms(
        field_point: &Vec3,
        bound_left: &Vec3,
        bound_right: &Vec3,
        prandtl_glauert_factor: Scalar,
        core_radius: Scalar,
        include_bound: bool,
    ) -> HorseshoeDerivatives {
        HorseshoeDerivatives::new(
            field_point,
            bound_left,
            bound_right,
            prandtl_glauert_factor,
            core_radius,
            include_bound,
        )
    }
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
    HorseshoeDerivatives::new(
        field_point,
        bound_left,
        bound_right,
        prandtl_glauert_factor,
        core_radius,
        include_bound,
    )
    .jacobians()
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
        let (vjp_a, vjp_b) = BoundTerms::new(&geometry, core_4).vjp(&geometry, &seed);
        result -= scale_x(vjp_a + vjp_b, prandtl_glauert_factor);
    }
    if geometry.a_sq > 0.0 {
        let vjp = LegTerms::new(&geometry.a, geometry.a_mag, core_4, -1.0).vjp(&seed);
        result -= scale_x(vjp, prandtl_glauert_factor);
    }
    if geometry.b_sq > 0.0 {
        let vjp = LegTerms::new(&geometry.b, geometry.b_mag, core_4, 1.0).vjp(&seed);
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
        let (vjp_a, _) = BoundTerms::new(&geometry, core_4).vjp(&geometry, &seed);
        result += scale_x(vjp_a, prandtl_glauert_factor);
    }
    if geometry.a_sq > 0.0 {
        let vjp = LegTerms::new(&geometry.a, geometry.a_mag, core_4, -1.0).vjp(&seed);
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
        let (_, vjp_b) = BoundTerms::new(&geometry, core_4).vjp(&geometry, &seed);
        result += scale_x(vjp_b, prandtl_glauert_factor);
    }
    if geometry.b_sq > 0.0 {
        let vjp = LegTerms::new(&geometry.b, geometry.b_mag, core_4, 1.0).vjp(&seed);
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
    HorseshoeDerivatives::new(
        field_point,
        bound_left,
        bound_right,
        prandtl_glauert_factor,
        core_radius,
        include_bound,
    )
    .pg_derivative()
}

/// Seed-independent derivative terms of one horseshoe at one field point.
///
/// Built once per (field point, horseshoe) pair; every method then contracts
/// the stored terms with seeds without touching the geometry again. The public
/// kernels in this module are thin wrappers over the same pieces, so
/// switching a caller between them changes nothing but speed.
pub(crate) struct HorseshoeDerivatives {
    geometry: Geometry,
    prandtl_glauert_factor: Scalar,
    /// `None` when the primal drops the bound segment.
    bound: Option<BoundTerms>,
    /// `None` when the left endpoint coincides with the field point.
    leg_a: Option<LegTerms>,
    /// `None` when the right endpoint coincides with the field point.
    leg_b: Option<LegTerms>,
}

impl HorseshoeDerivatives {
    /// Evaluate the shared terms. Arguments are those of
    /// [`horseshoe_vortex_velocity`](crate::solver::horseshoe_vortex_velocity).
    #[inline]
    pub(crate) fn new(
        field_point: &Vec3,
        bound_left: &Vec3,
        bound_right: &Vec3,
        prandtl_glauert_factor: Scalar,
        core_radius: Scalar,
        include_bound: bool,
    ) -> Self {
        let geometry = Geometry::new(field_point, bound_left, bound_right, prandtl_glauert_factor);
        let core_4 = core_four(core_radius);
        let bound = (include_bound && geometry.bound_is_active())
            .then(|| BoundTerms::new(&geometry, core_4));
        let leg_a =
            (geometry.a_sq > 0.0).then(|| LegTerms::new(&geometry.a, geometry.a_mag, core_4, -1.0));
        let leg_b =
            (geometry.b_sq > 0.0).then(|| LegTerms::new(&geometry.b, geometry.b_mag, core_4, 1.0));
        Self {
            geometry,
            prandtl_glauert_factor,
            bound,
            leg_a,
            leg_b,
        }
    }

    /// Single-seed VJP with respect to the left and right bound endpoints.
    ///
    /// The field-point VJP is `-(left + right)`: both endpoint vectors are
    /// measured from the field point. Equal to the Jacobians contracted with the
    /// same seed in exact arithmetic, at roughly the cost of one of their rows.
    #[inline]
    pub(crate) fn endpoint_vjp(&self, velocity_cotangent: &Vec3) -> (Vec3, Vec3) {
        let pg = self.prandtl_glauert_factor;
        let seed = effective_seed(velocity_cotangent, pg);
        let (mut left, mut right) = (Vec3::zeros(), Vec3::zeros());
        if let Some(bound) = &self.bound {
            let (vjp_a, vjp_b) = bound.vjp(&self.geometry, &seed);
            left += scale_x(vjp_a, pg);
            right += scale_x(vjp_b, pg);
        }
        if let Some(leg) = &self.leg_a {
            left += scale_x(leg.vjp(&seed), pg);
        }
        if let Some(leg) = &self.leg_b {
            right += scale_x(leg.vjp(&seed), pg);
        }
        (left, right)
    }

    /// All three geometry Jacobians, identical to
    /// [`horseshoe_velocity_jacobians`].
    #[inline]
    pub(crate) fn jacobians(&self) -> HorseshoeVelocityJacobians {
        let mut result = HorseshoeVelocityJacobians::default();
        for component in 0..3 {
            let mut basis = Vec3::zeros();
            basis[component] = 1.0;
            let (left, right) = self.endpoint_vjp(&basis);
            // `a` and `b` are both measured from the field point, so the
            // field-point row is minus the sum of the two endpoint rows.
            result.bound_left.set_row(component, &left.transpose());
            result.bound_right.set_row(component, &right.transpose());
            result
                .field_point
                .set_row(component, &(-(left + right)).transpose());
        }
        result
    }

    /// Derivative of the velocity with respect to the Prandtl-Glauert factor.
    ///
    /// The factor enters through the x-components of `a` and `b`, which scale
    /// as `1/pg`, and through the primal's final `u / pg`. Only the x-columns
    /// of the endpoint Jacobians are reachable, and each basis contraction
    /// reduces to its x-entry.
    #[inline]
    pub(crate) fn pg_derivative(&self) -> Vec3 {
        let pg = self.prandtl_glauert_factor;
        // d a_x / d pg = -(bound_left.x - field_point.x) / pg^2, which chains
        // through the 1/pg x-row already folded into the endpoint VJPs.
        let left_direction = -self.geometry.a.x;
        let right_direction = -self.geometry.b.x;
        let velocity_x = self
            .bound
            .as_ref()
            .map_or(0.0, |bound| bound.primal_u(&self.geometry) * INV_4_PI / pg);
        Vec3::from_fn(|component, _| {
            let mut basis = Vec3::zeros();
            basis[component] = 1.0;
            let (left, right) = self.endpoint_vjp(&basis);
            // Only the x-entries are needed, so the compiler drops the y and z
            // lanes of each VJP. A full dot product with `(dx, 0, 0)` would add
            // `y * 0 + z * 0`, which changes nothing but the sign of a result
            // that is exactly zero.
            let endpoint_derivative = left.x * left_direction + right.x * right_direction;
            endpoint_derivative - if component == 0 { velocity_x / pg } else { 0.0 }
        })
    }
}

fn core_four(core_radius: Scalar) -> Scalar {
    let core_sq = core_radius * core_radius;
    core_sq * core_sq
}

/// The cotangent of the unscaled `(u, v, w)` accumulator.
///
/// The primal returns `(u / (4 pi pg), v / (4 pi), w / (4 pi))`, so folding the
/// constant factors into the seed once keeps them out of every downstream term.
#[inline]
fn effective_seed(velocity_cotangent: &Vec3, prandtl_glauert_factor: Scalar) -> Vec3 {
    Vec3::new(
        velocity_cotangent.x * INV_4_PI / prandtl_glauert_factor,
        velocity_cotangent.y * INV_4_PI,
        velocity_cotangent.z * INV_4_PI,
    )
}

/// Chain through the Prandtl-Glauert stretch, which scales only the x row.
#[inline]
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
    #[inline]
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
    #[inline]
    fn bound_is_active(&self) -> bool {
        self.a_sq > 0.0 && self.b_sq > 0.0 && self.cross_sq > EPSILON * self.a_sq * self.b_sq
    }
}

/// Seed-independent part of the bound contribution `(a x b) * T`: `T` and its
/// gradients with respect to `a` and `b`.
struct BoundTerms {
    t: Scalar,
    d_t_a: Vec3,
    d_t_b: Vec3,
    /// Pieces of `T` kept so [`Self::primal_u`] can reproduce the primal's own
    /// association of the same quantity.
    num_a: Scalar,
    num_b: Scalar,
    root_a: Scalar,
    root_b: Scalar,
    t_den: Scalar,
}

impl BoundTerms {
    #[inline]
    fn new(geometry: &Geometry, core_4: Scalar) -> Self {
        let (a, b) = (&geometry.a, &geometry.b);
        let (a_sq, b_sq, cross_sq) = (geometry.a_sq, geometry.b_sq, geometry.cross_sq);
        let dot_ab = a.dot(b);
        let length_sq = a_sq + b_sq - 2.0 * dot_ab;

        let a_sq_sq = a_sq * a_sq;
        let b_sq_sq = b_sq * b_sq;
        let root_a = (a_sq_sq + core_4).sqrt().sqrt();
        let root_b = (b_sq_sq + core_4).sqrt().sqrt();
        let inv_r4_a = 1.0 / root_a;
        let inv_r4_b = 1.0 / root_b;

        let num_a = a_sq - dot_ab;
        let num_b = b_sq - dot_ab;
        let t_num = num_b * inv_r4_b + num_a * inv_r4_a;

        let den_arg = cross_sq * cross_sq + length_sq * length_sq * core_4;
        let t_den = den_arg.sqrt();
        let inv_t_den = 1.0 / t_den;
        let t = t_num * inv_t_den;

        // d/da_k of (a.a^2 + c^4)^(-1/4) is -a.a a_k / (a.a^2 + c^4)^(5/4).
        let d_inv_r4_a = -a_sq / ((a_sq_sq + core_4) * root_a);
        let d_inv_r4_b = -b_sq / ((b_sq_sq + core_4) * root_b);

        // d(|a x b|^2)/da_k = 2 (b x (a x b))_k, the cross-product Jacobian
        // identity below applied with the cross product itself as the seed.
        let b_cross_ab = b.cross(&geometry.cross_ab);
        let ab_cross_a = geometry.cross_ab.cross(a);

        // d(den_arg)^(-1/2)/dx = -1/2 den_arg^(-3/2) d(den_arg)/dx
        let d_inv_t_den = -0.5 * inv_t_den * inv_t_den * inv_t_den;

        let mut d_t_a = Vec3::zeros();
        let mut d_t_b = Vec3::zeros();
        for k in 0..3 {
            let a_k = a[k];
            let b_k = b[k];

            let d_t_num_a =
                (-b_k) * inv_r4_b + (2.0 * a_k - b_k) * inv_r4_a + num_a * d_inv_r4_a * a_k;
            let d_den_a = 4.0 * cross_sq * b_cross_ab[k] + 4.0 * length_sq * (a_k - b_k) * core_4;
            d_t_a[k] = d_t_num_a * inv_t_den + t_num * d_inv_t_den * d_den_a;

            let d_t_num_b =
                (2.0 * b_k - a_k) * inv_r4_b + num_b * d_inv_r4_b * b_k + (-a_k) * inv_r4_a;
            let d_den_b = 4.0 * cross_sq * ab_cross_a[k] + 4.0 * length_sq * (b_k - a_k) * core_4;
            d_t_b[k] = d_t_num_b * inv_t_den + t_num * d_inv_t_den * d_den_b;
        }
        Self {
            t,
            d_t_a,
            d_t_b,
            num_a,
            num_b,
            root_a,
            root_b,
            t_den,
        }
    }

    /// Seed-contracted derivatives of `(a x b) * T` with respect to `a` and `b`.
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
    #[inline]
    fn vjp(&self, geometry: &Geometry, seed: &Vec3) -> (Vec3, Vec3) {
        let seed_dot_cross = seed.dot(&geometry.cross_ab);
        let b_cross_seed = geometry.b.cross(seed);
        let seed_cross_a = seed.cross(&geometry.a);
        let mut vjp_a = Vec3::zeros();
        let mut vjp_b = Vec3::zeros();
        for k in 0..3 {
            vjp_a[k] = b_cross_seed[k] * self.t + seed_dot_cross * self.d_t_a[k];
            vjp_b[k] = seed_cross_a[k] * self.t + seed_dot_cross * self.d_t_b[k];
        }
        (vjp_a, vjp_b)
    }

    /// The primal's unscaled `u`, associated exactly as
    /// [`horseshoe_vortex_velocity`](crate::solver::horseshoe_vortex_velocity)
    /// associates it (divisions rather than the reciprocals `T` uses), so the
    /// Prandtl-Glauert derivative sees the primal's value bit for bit.
    #[inline]
    fn primal_u(&self, geometry: &Geometry) -> Scalar {
        let factor = (self.num_b / self.root_b + self.num_a / self.root_a) / self.t_den;
        geometry.cross_ab.x * factor
    }
}

/// Seed-independent part of one semi-infinite trailing leg.
///
/// The leg from endpoint `p` contributes `(0, p_z T, -p_y T)` with
/// `T = sign (1 - p_x / |p|) / sqrt((p_y^2 + p_z^2)^2 + c^4)`, so `sign` is
/// `-1` for the left leg and `+1` for the right.
struct LegTerms {
    p_y: Scalar,
    p_z: Scalar,
    t: Scalar,
    d_t: Vec3,
}

impl LegTerms {
    #[inline]
    fn new(p: &Vec3, p_mag: Scalar, core_4: Scalar, sign: Scalar) -> Self {
        let axis_sq = p.y * p.y + p.z * p.z;
        let denom = (axis_sq * axis_sq + core_4).sqrt();
        let factor = sign * (1.0 - p.x / p_mag);
        let t = factor / denom;

        let inv_denom = 1.0 / denom;
        let denom_cubed = denom * denom * denom;
        let p_mag_cubed = p_mag * p_mag * p_mag;

        // d(1 - p_x/|p|)/dp_x = -(|p|^2 - p_x^2)/|p|^3 = -axis_sq/|p|^3, and
        // the denominator has no p_x dependence at all.
        let d_t_x = -sign * axis_sq / (p_mag_cubed * denom);
        let d_t_y = sign * p.x * p.y / p_mag_cubed * inv_denom
            + factor * (-2.0 * axis_sq * p.y / denom_cubed);
        let d_t_z = sign * p.x * p.z / p_mag_cubed * inv_denom
            + factor * (-2.0 * axis_sq * p.z / denom_cubed);
        Self {
            p_y: p.y,
            p_z: p.z,
            t,
            d_t: Vec3::new(d_t_x, d_t_y, d_t_z),
        }
    }

    /// Seed-contracted derivative of the leg with respect to its endpoint
    /// vector. Only the y and z rows of the seed are reachable: a trailing leg
    /// induces no x velocity before the Prandtl-Glauert stretch.
    #[inline]
    fn vjp(&self, seed: &Vec3) -> Vec3 {
        let (p_y, p_z, t) = (self.p_y, self.p_z, self.t);
        let d_t = &self.d_t;
        // v = p_z T and w = -p_y T, so the p_z and p_y rows pick up the extra
        // term from differentiating the prefactor rather than T.
        Vec3::new(
            seed.y * p_z * d_t.x + seed.z * (-p_y * d_t.x),
            seed.y * p_z * d_t.y + seed.z * (-t - p_y * d_t.y),
            seed.y * (t + p_z * d_t.z) + seed.z * (-p_y * d_t.z),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solver::horseshoe_vortex_velocity;

    /// One kernel evaluation: field point, endpoints, pg, core, include_bound.
    #[derive(Clone, Copy)]
    struct Case {
        field: Vec3,
        left: Vec3,
        right: Vec3,
        pg: Scalar,
        core: Scalar,
        include_bound: bool,
    }

    impl Case {
        fn velocity(&self, field: &Vec3, left: &Vec3, right: &Vec3, pg: Scalar) -> Vec3 {
            horseshoe_vortex_velocity(field, left, right, pg, self.core, self.include_bound)
        }

        fn kernel(&self) -> HorseshoeDerivatives {
            HorseshoeDerivatives::new(
                &self.field,
                &self.left,
                &self.right,
                self.pg,
                self.core,
                self.include_bound,
            )
        }
    }

    const SEED: Vec3 = Vec3::new(0.3, -0.7, 0.9);

    fn base() -> Case {
        Case {
            field: Vec3::new(0.3, 0.2, 0.4),
            left: Vec3::new(0.0, -0.5, 0.0),
            right: Vec3::new(0.1, 0.6, 0.05),
            pg: 1.0,
            core: 1.0e-3,
            include_bound: true,
        }
    }

    /// The image `mach.rs` evaluates for a y-symmetry plane at `y = 0`:
    /// reflected endpoints, swapped so the circulation sense is preserved.
    fn y_image(case: Case) -> Case {
        Case {
            left: Vec3::new(case.right.x, -case.right.y, case.right.z),
            right: Vec3::new(case.left.x, -case.left.y, case.left.z),
            ..case
        }
    }

    /// Cases on which the kernel is smooth, so central differences apply.
    fn smooth_cases() -> Vec<Case> {
        let generic = base();
        let midpoint = 0.5 * (generic.left + generic.right);
        let span = generic.right - generic.left;
        let off_axis = span.cross(&Vec3::new(1.0, 0.0, 0.0)).normalize();
        vec![
            generic,
            Case { pg: 0.7, ..generic },
            Case {
                include_bound: false,
                ..generic
            },
            Case {
                include_bound: false,
                pg: 0.6,
                ..generic
            },
            // Inside the core: a third of a (large) core radius off the segment.
            Case {
                field: midpoint + 0.3 * 0.05 * off_axis,
                core: 0.05,
                pg: 0.8,
                ..generic
            },
            y_image(Case { pg: 0.9, ..generic }),
            // Field point well downstream, between the trailing legs.
            Case {
                field: Vec3::new(5.0, 0.05, -0.1),
                pg: 0.75,
                ..generic
            },
        ]
    }

    /// Cases on a branch boundary of the primal, where only agreement with the
    /// public kernels (not a difference quotient) is meaningful.
    fn degenerate_cases() -> Vec<Case> {
        let generic = base();
        let midpoint = 0.5 * (generic.left + generic.right);
        vec![
            // On the bound segment: the collinearity guard drops the bound leg.
            Case {
                field: midpoint,
                pg: 0.8,
                ..generic
            },
            // Self-influence at the bound centre with the bound leg excluded,
            // as the bound-force sweep evaluates it.
            Case {
                field: midpoint,
                include_bound: false,
                ..generic
            },
            // Coincident with an endpoint: that trailing leg and the bound leg
            // both vanish.
            Case {
                field: generic.left,
                ..generic
            },
            Case {
                field: generic.right,
                pg: 0.5,
                ..generic
            },
            y_image(Case {
                field: Vec3::new(0.05, 0.025, 0.025),
                ..generic
            }),
        ]
    }

    fn all_cases() -> Vec<Case> {
        let mut cases = smooth_cases();
        cases.extend(degenerate_cases());
        cases
    }

    fn assert_close(actual: Vec3, expected: Vec3, relative: Scalar, absolute: Scalar) {
        let scale = expected.abs().max().max(actual.abs().max());
        let error = (actual - expected).abs().max();
        assert!(
            error <= relative * scale + absolute,
            "{actual:?} vs {expected:?}: error {error:e}, scale {scale:e}"
        );
    }

    fn assert_same_bits(actual: Vec3, expected: Vec3) {
        for (a, e) in actual.iter().zip(expected.iter()) {
            assert_eq!(a.to_bits(), e.to_bits(), "{actual:?} vs {expected:?}");
        }
    }

    /// Central difference of `seed . V` along each coordinate of one point.
    fn difference(case: &Case, perturb: impl Fn(&mut Case, usize, Scalar)) -> Vec3 {
        let step = 1.0e-7;
        Vec3::from_fn(|k, _| {
            let mut plus = *case;
            let mut minus = *case;
            perturb(&mut plus, k, step);
            perturb(&mut minus, k, -step);
            let value = |c: &Case| SEED.dot(&c.velocity(&c.field, &c.left, &c.right, c.pg));
            (value(&plus) - value(&minus)) / (2.0 * step)
        })
    }

    #[test]
    fn the_endpoint_vjp_is_the_public_endpoint_vjps_bit_for_bit() {
        for case in all_cases() {
            let (left, right) = case.kernel().endpoint_vjp(&SEED);
            let args = (&case.field, &case.left, &case.right, case.pg, case.core);
            assert_same_bits(
                left,
                horseshoe_velocity_vjp_bound_left(
                    args.0,
                    args.1,
                    args.2,
                    args.3,
                    args.4,
                    case.include_bound,
                    &SEED,
                ),
            );
            assert_same_bits(
                right,
                horseshoe_velocity_vjp_bound_right(
                    args.0,
                    args.1,
                    args.2,
                    args.3,
                    args.4,
                    case.include_bound,
                    &SEED,
                ),
            );
            // The field-point VJP sums the same terms in another order.
            let field = horseshoe_velocity_vjp_field_point(
                args.0,
                args.1,
                args.2,
                args.3,
                args.4,
                case.include_bound,
                &SEED,
            );
            assert_close(-(left + right), field, 1.0e-14, 1.0e-300);
            assert!(
                left.iter()
                    .chain(right.iter())
                    .all(|value| value.is_finite())
            );
        }
    }

    #[test]
    fn the_endpoint_vjp_matches_the_jacobians_contracted_with_the_seed() {
        for case in all_cases() {
            let kernel = case.kernel();
            let (left, right) = kernel.endpoint_vjp(&SEED);
            let jacobians = kernel.jacobians();
            assert_close(left, jacobians.bound_left.tr_mul(&SEED), 1.0e-14, 1.0e-300);
            assert_close(
                right,
                jacobians.bound_right.tr_mul(&SEED),
                1.0e-14,
                1.0e-300,
            );
            assert_close(
                -(left + right),
                jacobians.field_point.tr_mul(&SEED),
                1.0e-14,
                1.0e-300,
            );
        }
    }

    #[test]
    fn the_endpoint_vjp_matches_central_differences() {
        for case in smooth_cases() {
            let (left, right) = case.kernel().endpoint_vjp(&SEED);
            let by_left = difference(&case, |c, k, h| c.left[k] += h);
            let by_right = difference(&case, |c, k, h| c.right[k] += h);
            let by_field = difference(&case, |c, k, h| c.field[k] += h);
            assert_close(left, by_left, 1.0e-6, 1.0e-8);
            assert_close(right, by_right, 1.0e-6, 1.0e-8);
            assert_close(-(left + right), by_field, 1.0e-6, 1.0e-8);
        }
    }

    #[test]
    fn the_jacobians_match_central_differences_column_by_column() {
        let step = 1.0e-7;
        for case in smooth_cases() {
            let jacobians = case.kernel().jacobians();
            for k in 0..3 {
                let column = |shift: fn(&mut Case, usize, Scalar)| {
                    let mut plus = case;
                    let mut minus = case;
                    shift(&mut plus, k, step);
                    shift(&mut minus, k, -step);
                    (plus.velocity(&plus.field, &plus.left, &plus.right, plus.pg)
                        - minus.velocity(&minus.field, &minus.left, &minus.right, minus.pg))
                        / (2.0 * step)
                };
                assert_close(
                    jacobians.bound_left.column(k).into(),
                    column(|c, k, h| c.left[k] += h),
                    1.0e-6,
                    1.0e-8,
                );
                assert_close(
                    jacobians.bound_right.column(k).into(),
                    column(|c, k, h| c.right[k] += h),
                    1.0e-6,
                    1.0e-8,
                );
                assert_close(
                    jacobians.field_point.column(k).into(),
                    column(|c, k, h| c.field[k] += h),
                    1.0e-6,
                    1.0e-8,
                );
            }
        }
    }

    #[test]
    fn the_pg_derivative_matches_central_differences() {
        let step = 1.0e-7;
        for case in smooth_cases() {
            let analytic = case.kernel().pg_derivative();
            let finite = (case.velocity(&case.field, &case.left, &case.right, case.pg + step)
                - case.velocity(&case.field, &case.left, &case.right, case.pg - step))
                / (2.0 * step);
            assert_close(analytic, finite, 1.0e-6, 1.0e-8);
        }
    }

    #[test]
    fn the_pg_derivative_agrees_with_the_public_kernel_and_its_guard() {
        for case in all_cases() {
            let public = horseshoe_velocity_pg_derivative(
                &case.field,
                &case.left,
                &case.right,
                case.pg,
                case.core,
                case.include_bound,
            );
            assert_eq!(case.kernel().pg_derivative(), public);
            assert!(public.iter().all(|value| value.is_finite()));
        }
        let case = base();
        let degenerate = horseshoe_velocity_pg_derivative(
            &case.field,
            &case.left,
            &case.right,
            0.0,
            case.core,
            true,
        );
        assert_eq!(degenerate, Vec3::zeros());
    }

    #[test]
    fn a_missing_bound_leg_leaves_only_the_trailing_legs() {
        // Coincident with the left endpoint: only the right trailing leg moves
        // with the geometry, so the left-endpoint VJP is exactly zero.
        let case = Case {
            field: base().left,
            ..base()
        };
        let kernel = case.kernel();
        assert!(kernel.bound.is_none() && kernel.leg_a.is_none() && kernel.leg_b.is_some());
        let (left, right) = kernel.endpoint_vjp(&SEED);
        assert_eq!(left, Vec3::zeros());
        assert!(right.norm() > 0.0);
        // On the segment the collinearity guard drops the bound leg.
        let on_segment = Case {
            field: 0.5 * (base().left + base().right),
            ..base()
        };
        assert!(on_segment.kernel().bound.is_none());
    }
}
