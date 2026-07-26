//! Exact VJPs of the regularized horseshoe-vortex kernel.

use crate::core::{EPSILON, PI, Scalar, Vec3};
use crate::solver::horseshoe_vortex_velocity;

use super::dual::Dual;

type D3 = Dual<3>;
type DVec3 = [D3; 3];

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
    let field = variables(field_point);
    let left = constants(bound_left);
    let right = constants(bound_right);
    pullback(
        kernel(
            field,
            left,
            right,
            prandtl_glauert_factor,
            core_radius,
            include_bound,
        ),
        velocity_cotangent,
    )
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
    pullback(
        kernel(
            constants(field_point),
            variables(bound_left),
            constants(bound_right),
            prandtl_glauert_factor,
            core_radius,
            include_bound,
        ),
        velocity_cotangent,
    )
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
    pullback(
        kernel(
            constants(field_point),
            constants(bound_left),
            variables(bound_right),
            prandtl_glauert_factor,
            core_radius,
            include_bound,
        ),
        velocity_cotangent,
    )
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

fn constants(vector: &Vec3) -> DVec3 {
    [
        D3::constant(vector.x),
        D3::constant(vector.y),
        D3::constant(vector.z),
    ]
}

fn variables(vector: &Vec3) -> DVec3 {
    [
        D3::variable(vector.x, 0),
        D3::variable(vector.y, 1),
        D3::variable(vector.z, 2),
    ]
}

fn pullback(output: DVec3, cotangent: &Vec3) -> Vec3 {
    Vec3::from_fn(|index, _| {
        output[0].derivative[index] * cotangent.x
            + output[1].derivative[index] * cotangent.y
            + output[2].derivative[index] * cotangent.z
    })
}

fn dot(left: DVec3, right: DVec3) -> D3 {
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
}

fn cross(left: DVec3, right: DVec3) -> DVec3 {
    [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]
}

fn kernel(
    field_point: DVec3,
    bound_left: DVec3,
    bound_right: DVec3,
    pg: Scalar,
    core_radius: Scalar,
    include_bound: bool,
) -> DVec3 {
    let pg = D3::constant(pg);
    let a = [
        (bound_left[0] - field_point[0]) / pg,
        bound_left[1] - field_point[1],
        bound_left[2] - field_point[2],
    ];
    let b = [
        (bound_right[0] - field_point[0]) / pg,
        bound_right[1] - field_point[1],
        bound_right[2] - field_point[2],
    ];
    let a_sq = dot(a, a);
    let b_sq = dot(b, b);
    let a_mag = a_sq.sqrt();
    let b_mag = b_sq.sqrt();
    let core_4 = D3::constant(core_radius.powi(4));
    let mut output = [D3::constant(0.0); 3];
    if include_bound && a_mag.value * b_mag.value > EPSILON {
        let cross_ab = cross(a, b);
        let cross_sq = dot(cross_ab, cross_ab);
        if cross_sq.value > EPSILON {
            let dot_ab = dot(a, b);
            let length_sq = a_sq + b_sq - D3::constant(2.0) * dot_ab;
            let factor = ((b_sq - dot_ab) / (b_sq * b_sq + core_4).sqrt().sqrt()
                + (a_sq - dot_ab) / (a_sq * a_sq + core_4).sqrt().sqrt())
                / (cross_sq * cross_sq + length_sq * length_sq * core_4).sqrt();
            for component in 0..3 {
                output[component] = cross_ab[component] * factor;
            }
        }
    }
    if a_mag.value > EPSILON {
        let axis_sq = a[1] * a[1] + a[2] * a[2];
        let factor = -(D3::constant(1.0) - a[0] / a_mag) / (axis_sq * axis_sq + core_4).sqrt();
        output[1] = output[1] + a[2] * factor;
        output[2] = output[2] - a[1] * factor;
    }
    if b_mag.value > EPSILON {
        let axis_sq = b[1] * b[1] + b[2] * b[2];
        let factor = (D3::constant(1.0) - b[0] / b_mag) / (axis_sq * axis_sq + core_4).sqrt();
        output[1] = output[1] + b[2] * factor;
        output[2] = output[2] - b[1] * factor;
    }
    let factor = D3::constant(1.0 / (4.0 * PI));
    [
        output[0] * factor / pg,
        output[1] * factor,
        output[2] * factor,
    ]
}
