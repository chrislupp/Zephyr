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

//! Exact rotation tangents and pullbacks shared by control and mesh gradients.

use super::{ControlSurfaceInfluence, FlightCondition};
use crate::core::{EPSILON, Scalar, Vec3};

pub(crate) fn rotate(vector: Vec3, axis: &Vec3, angle: Scalar) -> Vec3 {
    let (sin, cos) = angle.sin_cos();
    vector * cos + axis.cross(&vector) * sin + axis * axis.dot(&vector) * (1.0 - cos)
}

#[cfg(test)]
fn axis_jvp(vector: Vec3, axis: Vec3, direction: Vec3, angle: Scalar) -> Vec3 {
    let (sin, cos) = angle.sin_cos();
    direction.cross(&vector) * sin
        + (direction * axis.dot(&vector) + axis * direction.dot(&vector)) * (1.0 - cos)
}

/// Derivative with respect to the base normal, normalized axes and gains.
///
/// The geometry adjoints run in reverse through [`normal_vjp_with_gains`];
/// this forward form is kept to pin that transpose in the tests below.
#[cfg(test)]
fn normal_jvp(
    mut normal: Vec3,
    mut tangent: Vec3,
    influences: &[ControlSurfaceInfluence],
    axes: &[Vec3],
    gains: &[Scalar],
    condition: &FlightCondition,
) -> Vec3 {
    for (name, deflection) in &condition.control_deflections {
        for (k, influence) in influences
            .iter()
            .enumerate()
            .filter(|(_, v)| v.name == *name)
        {
            if influence.hinge_axis.norm_squared() <= EPSILON * EPSILON {
                continue;
            }
            let axis = influence.hinge_axis.normalize();
            let angle = deflection * influence.gain;
            tangent = rotate(tangent, &axis, angle) + axis_jvp(normal, axis, axes[k], angle);
            normal = rotate(normal, &axis, angle);
            tangent += deflection * gains[k] * axis.cross(&normal);
        }
    }
    tangent
}

/// Pull a normal cotangent back through every rotation, in reverse order.
/// Axis results are with respect to the normalized influence axes.
pub(crate) fn normal_vjp(
    normal: Vec3,
    seed: Vec3,
    influences: &[ControlSurfaceInfluence],
    condition: &FlightCondition,
) -> (Vec3, Vec<Vec3>) {
    let (seed, axes, _) = normal_vjp_with_gains(normal, seed, influences, condition);
    (seed, axes)
}

/// The transpose of [`normal_jvp`]: base-normal, normalized-axis and gain
/// cotangents of one rotated normal.
pub(crate) fn normal_vjp_with_gains(
    mut normal: Vec3,
    mut seed: Vec3,
    influences: &[ControlSurfaceInfluence],
    condition: &FlightCondition,
) -> (Vec3, Vec<Vec3>, Vec<Scalar>) {
    let mut tape = Vec::new();
    for (name, deflection) in &condition.control_deflections {
        for (k, influence) in influences
            .iter()
            .enumerate()
            .filter(|(_, v)| v.name == *name)
        {
            if influence.hinge_axis.norm_squared() <= EPSILON * EPSILON {
                continue;
            }
            let axis = influence.hinge_axis.normalize();
            let angle = deflection * influence.gain;
            tape.push((k, normal, axis, angle, *deflection));
            normal = rotate(normal, &axis, angle);
        }
    }
    let mut axes = vec![Vec3::zeros(); influences.len()];
    let mut gains = vec![0.0; influences.len()];
    for (k, before, axis, angle, deflection) in tape.into_iter().rev() {
        let (sin, cos) = angle.sin_cos();
        let after = rotate(before, &axis, angle);
        gains[k] += deflection * axis.cross(&after).dot(&seed);
        axes[k] += sin * before.cross(&seed)
            + (1.0 - cos) * (axis.dot(&before) * seed + axis.dot(&seed) * before);
        seed = rotate(seed, &axis, -angle);
    }
    (seed, axes, gains)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `<seed, J t> == <J^T seed, t>` for every tangent channel, with repeated
    /// names and two influences on one panel, so the reverse pullback the
    /// geometry adjoints use is exactly the transpose of the forward rotation.
    #[test]
    fn the_normal_pullback_is_the_transpose_of_the_rotation_tangent() {
        let influences = vec![
            ControlSurfaceInfluence {
                name: "flap".into(),
                gain: 0.8,
                hinge_axis: Vec3::new(0.1, 1.0, 0.05).normalize(),
                hinge_point: Vec3::zeros(),
                reflection_sign: 1.0,
            },
            ControlSurfaceInfluence {
                name: "aileron".into(),
                gain: -1.3,
                hinge_axis: Vec3::new(-0.2, 0.9, 0.3).normalize(),
                hinge_point: Vec3::zeros(),
                reflection_sign: 1.0,
            },
        ];
        let condition = FlightCondition {
            control_deflections: vec![
                ("flap".into(), 0.21),
                ("aileron".into(), -0.17),
                ("flap".into(), 0.05),
            ],
            ..FlightCondition::default()
        };
        let normal = Vec3::new(0.05, -0.1, 1.0).normalize();
        let tangent = Vec3::new(0.3, -0.7, 0.2);
        let axes = [Vec3::new(0.4, 0.1, -0.6), Vec3::new(-0.2, 0.5, 0.3)];
        let gains = [0.7, -0.4];
        let seed = Vec3::new(-0.6, 0.25, 0.9);

        let forward = normal_jvp(normal, tangent, &influences, &axes, &gains, &condition);
        let (by_normal, by_axes, by_gains) =
            normal_vjp_with_gains(normal, seed, &influences, &condition);
        let reverse = by_normal.dot(&tangent)
            + by_axes
                .iter()
                .zip(&axes)
                .map(|(c, t)| c.dot(t))
                .sum::<Scalar>()
            + by_gains
                .iter()
                .zip(&gains)
                .map(|(c, t)| c * t)
                .sum::<Scalar>();
        assert!((seed.dot(&forward) - reverse).abs() < 1e-14);

        let (plain_normal, plain_axes) = normal_vjp(normal, seed, &influences, &condition);
        assert_eq!(plain_normal, by_normal);
        assert_eq!(plain_axes, by_axes);
    }
}
