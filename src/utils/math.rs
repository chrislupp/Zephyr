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

//! Small mathematical kernels shared by the solver and gradient code.

use crate::core::{EPSILON, Mat3, PI, Scalar, SpacingType, Vec3};

/// Prandtl–Glauert factor `sqrt(1-M²)`.
#[must_use]
pub fn prandtl_glauert_factor(mach: Scalar) -> Scalar {
    (1.0 - mach * mach).sqrt()
}

/// Velocity induced by a finite unit-circulation vortex segment.
#[must_use]
pub fn vortex_segment_velocity(
    field: &Vec3,
    start: &Vec3,
    end: &Vec3,
    core_radius: Scalar,
) -> Vec3 {
    let r1 = field - start;
    let r2 = field - end;
    let r0 = end - start;
    let cross = r1.cross(&r2);
    let cross_sq = cross.norm_squared().max(core_radius * core_radius);
    let r1_norm = r1.norm();
    let r2_norm = r2.norm();
    if r1_norm < EPSILON || r2_norm < EPSILON {
        return Vec3::zeros();
    }
    let factor = r0.dot(&(r1 / r1_norm - r2 / r2_norm)) / (4.0 * PI * cross_sq);
    factor * cross
}

/// Velocity induced by a semi-infinite unit-circulation vortex along +X.
#[must_use]
pub fn trailing_vortex_velocity(field: &Vec3, start: &Vec3, core_radius: Scalar) -> Vec3 {
    let r = field - start;
    let cross = Vec3::new(0.0, r.z, -r.y);
    let cross_sq = cross.norm_squared().max(core_radius * core_radius);
    let norm = r.norm();
    if norm < EPSILON {
        return Vec3::zeros();
    }
    ((1.0 + r.x / norm) / (4.0 * PI * cross_sq)) * cross
}

/// Uniform points in `[0,1]`.
#[must_use]
pub fn uniform_spacing(count: usize) -> Vec<Scalar> {
    if count < 2 {
        return vec![0.0; count];
    }
    (0..count)
        .map(|i| i as Scalar / (count - 1) as Scalar)
        .collect()
}

/// Cosine-clustered points in `[0,1]`.
#[must_use]
pub fn cosine_spacing(count: usize) -> Vec<Scalar> {
    if count < 2 {
        return vec![0.0; count];
    }
    (0..count)
        .map(|i| {
            let theta = PI * i as Scalar / (count - 1) as Scalar;
            0.5 * (1.0 - theta.cos())
        })
        .collect()
}

/// Sine-clustered points in `[0,1]`.
#[must_use]
pub fn sine_spacing(count: usize) -> Vec<Scalar> {
    if count < 2 {
        return vec![0.0; count];
    }
    (0..count)
        .map(|i| (0.5 * PI * i as Scalar / (count - 1) as Scalar).sin())
        .collect()
}

/// Linear blend of uniform and cosine spacing.
#[must_use]
pub fn blended_spacing(count: usize, blend: Scalar) -> Vec<Scalar> {
    uniform_spacing(count)
        .into_iter()
        .zip(cosine_spacing(count))
        .map(|(uniform, cosine)| (1.0 - blend) * uniform + blend * cosine)
        .collect()
}

/// Generate the requested spacing distribution.
#[must_use]
pub fn spacing(count: usize, distribution: SpacingType) -> Vec<Scalar> {
    match distribution {
        SpacingType::None | SpacingType::Uniform => uniform_spacing(count),
        SpacingType::Cosine => cosine_spacing(count),
        SpacingType::Sine => sine_spacing(count),
        SpacingType::Blend => blended_spacing(count, 0.5),
    }
}

/// Rotate around X.
#[must_use]
pub fn rotate_about_x(v: &Vec3, angle: Scalar) -> Vec3 {
    let (s, c) = angle.sin_cos();
    Vec3::new(v.x, c * v.y - s * v.z, s * v.y + c * v.z)
}

/// Rotate around Y.
#[must_use]
pub fn rotate_about_y(v: &Vec3, angle: Scalar) -> Vec3 {
    let (s, c) = angle.sin_cos();
    Vec3::new(c * v.x + s * v.z, v.y, -s * v.x + c * v.z)
}

/// Rotate around Z.
#[must_use]
pub fn rotate_about_z(v: &Vec3, angle: Scalar) -> Vec3 {
    let (s, c) = angle.sin_cos();
    Vec3::new(c * v.x - s * v.y, s * v.x + c * v.y, v.z)
}

/// Wind-to-body rotation matrix.
#[must_use]
pub fn wind_to_body_rotation(alpha: Scalar, beta: Scalar) -> Mat3 {
    let (sa, ca) = alpha.sin_cos();
    let (sb, cb) = beta.sin_cos();
    Mat3::new(ca * cb, -ca * sb, -sa, sb, cb, 0.0, sa * cb, -sa * sb, ca)
}

/// Body-to-stability rotation matrix.
#[must_use]
pub fn body_to_stability_rotation(alpha: Scalar) -> Mat3 {
    let (sa, ca) = alpha.sin_cos();
    Mat3::new(ca, 0.0, sa, 0.0, 1.0, 0.0, -sa, 0.0, ca)
}
