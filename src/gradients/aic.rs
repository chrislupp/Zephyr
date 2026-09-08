//! Safe AIC and velocity-influence geometry derivatives.

use crate::core::{Scalar, Vec3, VecX};
use crate::solver::{VortexLattice, horseshoe_vortex_velocity};

use super::{
    horseshoe_velocity_vjp_bound_left, horseshoe_velocity_vjp_bound_right,
    horseshoe_velocity_vjp_field_point,
};

/// Gradient of one scalar AIC entry.
#[derive(Clone, Copy, Debug, Default)]
pub struct AicElementGradient {
    /// Derivative with respect to the source left endpoint.
    pub bound_left: Vec3,
    /// Derivative with respect to the source right endpoint.
    pub bound_right: Vec3,
    /// Derivative with respect to the target control point.
    pub control_point: Vec3,
    /// Derivative with respect to the target normal.
    pub normal: Vec3,
}

/// Differentiate `normal · horseshoe_velocity(control_point, left, right)`.
#[must_use]
pub fn aic_element_gradient(
    bound_left: &Vec3,
    bound_right: &Vec3,
    control_point: &Vec3,
    normal: &Vec3,
    core_radius: Scalar,
    prandtl_glauert_factor: Scalar,
    include_bound: bool,
) -> AicElementGradient {
    AicElementGradient {
        bound_left: horseshoe_velocity_vjp_bound_left(
            control_point,
            bound_left,
            bound_right,
            prandtl_glauert_factor,
            core_radius,
            include_bound,
            normal,
        ),
        bound_right: horseshoe_velocity_vjp_bound_right(
            control_point,
            bound_left,
            bound_right,
            prandtl_glauert_factor,
            core_radius,
            include_bound,
            normal,
        ),
        control_point: horseshoe_velocity_vjp_field_point(
            control_point,
            bound_left,
            bound_right,
            prandtl_glauert_factor,
            core_radius,
            include_bound,
            normal,
        ),
        normal: horseshoe_vortex_velocity(
            control_point,
            bound_left,
            bound_right,
            prandtl_glauert_factor,
            core_radius,
            include_bound,
        ),
    }
}

fn core_radius(
    mesh: &VortexLattice,
    source: usize,
    target_component: usize,
    chord_fraction: Scalar,
    width_fraction: Scalar,
) -> Scalar {
    let vortex = &mesh.vortices()[source];
    let width = (vortex.bound_vortex_right.y - vortex.bound_vortex_left.y)
        .hypot(vortex.bound_vortex_right.z - vortex.bound_vortex_left.z);
    if vortex.component_index == target_component {
        0.0001 * width
    } else {
        (chord_fraction * vortex.local_chord).max(width_fraction * width)
    }
}

/// Differentiate a weighted AIC row with respect to its target control point
/// and normal.
#[must_use]
pub fn aic_row_gradient(
    mesh: &VortexLattice,
    row: usize,
    weights: &VecX,
    prandtl_glauert_factor: Scalar,
    core_chord_fraction: Scalar,
    core_width_fraction: Scalar,
) -> (Vec3, Vec3) {
    let target = &mesh.vortices()[row];
    let mut control_point = Vec3::zeros();
    let mut normal = Vec3::zeros();
    for (source_index, (source, weight)) in mesh.vortices().iter().zip(weights).enumerate() {
        if *weight == 0.0 {
            continue;
        }
        let core = core_radius(
            mesh,
            source_index,
            target.component_index,
            core_chord_fraction,
            core_width_fraction,
        );
        let gradient = aic_element_gradient(
            &source.bound_vortex_left,
            &source.bound_vortex_right,
            &target.control_point,
            &target.normal_at_control_point,
            core,
            prandtl_glauert_factor,
            true,
        );
        control_point += *weight * gradient.control_point;
        normal += *weight * gradient.normal;
    }
    (control_point, normal)
}

/// Differentiate a weighted AIC column with respect to its source endpoints.
#[must_use]
pub fn aic_column_gradient(
    mesh: &VortexLattice,
    column: usize,
    weights: &VecX,
    prandtl_glauert_factor: Scalar,
    core_chord_fraction: Scalar,
    core_width_fraction: Scalar,
) -> (Vec3, Vec3) {
    let source = &mesh.vortices()[column];
    let mut left = Vec3::zeros();
    let mut right = Vec3::zeros();
    for (target, weight) in mesh.vortices().iter().zip(weights) {
        if *weight == 0.0 {
            continue;
        }
        let core = core_radius(
            mesh,
            column,
            target.component_index,
            core_chord_fraction,
            core_width_fraction,
        );
        let gradient = aic_element_gradient(
            &source.bound_vortex_left,
            &source.bound_vortex_right,
            &target.control_point,
            &target.normal_at_control_point,
            core,
            prandtl_glauert_factor,
            true,
        );
        left += *weight * gradient.bound_left;
        right += *weight * gradient.bound_right;
    }
    (left, right)
}

/// Differentiate one weighted force-velocity row with respect to its bound-force point.
#[must_use]
pub fn velocity_row_gradient(
    mesh: &VortexLattice,
    row: usize,
    circulation: &VecX,
    velocity_cotangent: &Vec3,
    prandtl_glauert_factor: Scalar,
    core_chord_fraction: Scalar,
    core_width_fraction: Scalar,
) -> Vec3 {
    let target = &mesh.vortices()[row];
    let mut result = Vec3::zeros();
    for (column, (source, gamma)) in mesh.vortices().iter().zip(circulation).enumerate() {
        if *gamma == 0.0 {
            continue;
        }
        let core = core_radius(
            mesh,
            column,
            target.component_index,
            core_chord_fraction,
            core_width_fraction,
        );
        result += *gamma
            * horseshoe_velocity_vjp_field_point(
                &target.bound_vortex_center,
                &source.bound_vortex_left,
                &source.bound_vortex_right,
                prandtl_glauert_factor,
                core,
                row != column,
                velocity_cotangent,
            );
    }
    result
}

/// Differentiate one velocity column with respect to its source endpoints.
#[must_use]
pub fn velocity_column_gradient(
    mesh: &VortexLattice,
    column: usize,
    gamma: Scalar,
    velocity_cotangents: &[Vec3],
    prandtl_glauert_factor: Scalar,
    core_chord_fraction: Scalar,
    core_width_fraction: Scalar,
) -> (Vec3, Vec3) {
    if gamma == 0.0 {
        return (Vec3::zeros(), Vec3::zeros());
    }
    let source = &mesh.vortices()[column];
    let mut left = Vec3::zeros();
    let mut right = Vec3::zeros();
    for (row, (target, cotangent)) in mesh.vortices().iter().zip(velocity_cotangents).enumerate() {
        let core = core_radius(
            mesh,
            column,
            target.component_index,
            core_chord_fraction,
            core_width_fraction,
        );
        left += gamma
            * horseshoe_velocity_vjp_bound_left(
                &target.bound_vortex_center,
                &source.bound_vortex_left,
                &source.bound_vortex_right,
                prandtl_glauert_factor,
                core,
                row != column,
                cotangent,
            );
        right += gamma
            * horseshoe_velocity_vjp_bound_right(
                &target.bound_vortex_center,
                &source.bound_vortex_left,
                &source.bound_vortex_right,
                prandtl_glauert_factor,
                core,
                row != column,
                cotangent,
            );
    }
    (left, right)
}
