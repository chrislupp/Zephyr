//! Exact rotation tangents and pullbacks shared by control and mesh gradients.

use super::{ControlSurfaceInfluence, FlightCondition};
use crate::core::{EPSILON, Scalar, Vec3};

pub(crate) fn rotate(vector: Vec3, axis: &Vec3, angle: Scalar) -> Vec3 {
    let (sin, cos) = angle.sin_cos();
    vector * cos + axis.cross(&vector) * sin + axis * axis.dot(&vector) * (1.0 - cos)
}

fn axis_jvp(vector: Vec3, axis: Vec3, direction: Vec3, angle: Scalar) -> Vec3 {
    let (sin, cos) = angle.sin_cos();
    direction.cross(&vector) * sin
        + (direction * axis.dot(&vector) + axis * direction.dot(&vector)) * (1.0 - cos)
}

/// Derivative with respect to the base normal, normalized axes and gains.
pub(crate) fn normal_jvp(
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
    mut normal: Vec3,
    mut seed: Vec3,
    influences: &[ControlSurfaceInfluence],
    condition: &FlightCondition,
) -> (Vec3, Vec<Vec3>) {
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
            tape.push((k, normal, axis, angle));
            normal = rotate(normal, &axis, angle);
        }
    }
    let mut axes = vec![Vec3::zeros(); influences.len()];
    for (k, before, axis, angle) in tape.into_iter().rev() {
        let (sin, cos) = angle.sin_cos();
        axes[k] += sin * before.cross(&seed)
            + (1.0 - cos) * (axis.dot(&before) * seed + axis.dot(&seed) * before);
        seed = rotate(seed, &axis, -angle);
    }
    (seed, axes)
}
