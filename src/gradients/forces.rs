//! Analytical derivatives of the Kutta–Joukowski force kernel.

use crate::core::{Scalar, Vec3, VecX};
use crate::error::{PegasusError, Result};
use crate::geometry::ReferenceQuantities;
use crate::solver::{FlightCondition, SolutionState, VortexLattice, VortexLatticeSolver};

use super::{AeroOutput, OutputSelection};

/// Circulation gradients for the six primary coefficients.
#[derive(Clone, Debug, Default)]
pub struct ForceCoefficientGradients {
    pub lift: VecX,
    pub drag: VecX,
    pub side: VecX,
    pub roll: VecX,
    pub pitch: VecX,
    pub yaw: VecX,
}

impl ForceCoefficientGradients {
    #[must_use]
    pub fn get(&self, output: AeroOutput) -> Option<&VecX> {
        match output {
            AeroOutput::Lift => Some(&self.lift),
            AeroOutput::Drag => Some(&self.drag),
            AeroOutput::Side => Some(&self.side),
            AeroOutput::Roll => Some(&self.roll),
            AeroOutput::Pitch => Some(&self.pitch),
            AeroOutput::Yaw => Some(&self.yaw),
            AeroOutput::InducedDrag => None,
        }
    }
}

/// Fixed-circulation explicit derivative with respect to angle of attack.
#[must_use]
pub fn force_alpha_gradient(
    mesh: &VortexLattice,
    solution: &SolutionState,
    condition: &FlightCondition,
    reference: &ReferenceQuantities,
) -> [Scalar; 6] {
    let (sin_alpha, cos_alpha) = condition.alpha.sin_cos();
    let (_, cos_beta) = condition.beta.sin_cos();
    let freestream = VortexLatticeSolver::freestream_direction(condition);
    let derivative = Vec3::new(-sin_alpha * cos_beta, 0.0, cos_alpha * cos_beta);
    let (base_force, _) = base_force_and_moment(mesh, solution, reference, freestream);
    let (force_derivative, moment_derivative) =
        force_and_moment(mesh, solution, reference, derivative);
    [
        force_derivative.z * cos_alpha
            - force_derivative.x * sin_alpha
            - base_force.z * sin_alpha
            - base_force.x * cos_alpha,
        force_derivative.x * cos_alpha + force_derivative.z * sin_alpha - base_force.x * sin_alpha
            + base_force.z * cos_alpha,
        force_derivative.y,
        -moment_derivative.x,
        moment_derivative.y,
        -moment_derivative.z,
    ]
}

/// Fixed-circulation explicit derivative with respect to sideslip.
#[must_use]
pub fn force_beta_gradient(
    mesh: &VortexLattice,
    solution: &SolutionState,
    condition: &FlightCondition,
    reference: &ReferenceQuantities,
) -> [Scalar; 6] {
    let (sin_alpha, cos_alpha) = condition.alpha.sin_cos();
    let (sin_beta, cos_beta) = condition.beta.sin_cos();
    let derivative = Vec3::new(-cos_alpha * sin_beta, cos_beta, -sin_alpha * sin_beta);
    let (force, moment) = force_and_moment(mesh, solution, reference, derivative);
    [
        force.z * cos_alpha - force.x * sin_alpha,
        force.x * cos_alpha + force.z * sin_alpha,
        force.y,
        -moment.x,
        moment.y,
        -moment.z,
    ]
}

/// Total partial `∂J/∂Γ`, including induced-velocity dependence.
#[must_use]
pub fn force_circulation_gradient(
    solver: &VortexLatticeSolver,
    condition: &FlightCondition,
    output: AeroOutput,
) -> VecX {
    let mesh = solver.mesh();
    let solution = solver.solution();
    let reference = solver.reference_quantities();
    let n = mesh.vortex_count();
    if !solution.is_valid || matches!(output, AeroOutput::InducedDrag) {
        return VecX::zeros(n);
    }
    let (force_seed, moment_seed) = body_output_cotangent(output, condition.alpha);
    let freestream = VortexLatticeSolver::freestream_direction(condition);
    let mut direct = VecX::zeros(n);
    let mut velocity_cotangent = VecX::zeros(3 * n);
    for i in 0..n {
        let vortex = &mesh.vortices()[i];
        let velocity = freestream
            + if solution.induced_velocity.nrows() > i {
                Vec3::new(
                    solution.induced_velocity[(i, 0)],
                    solution.induced_velocity[(i, 1)],
                    solution.induced_velocity[(i, 2)],
                )
            } else {
                Vec3::zeros()
            };
        let line = vortex.bound_vortex_right - vortex.bound_vortex_left;
        let arm = vortex.control_point - reference.moment_reference_point;
        let combined_seed = force_seed + moment_seed.cross(&arm) / reference.chord;
        direct[i] = 2.0 / reference.area * combined_seed.dot(&velocity.cross(&line));
        let cotangent = 2.0 * solution.circulation[i] / reference.area * line.cross(&combined_seed);
        velocity_cotangent[3 * i] = cotangent.x;
        velocity_cotangent[3 * i + 1] = cotangent.y;
        velocity_cotangent[3 * i + 2] = cotangent.z;
    }
    direct + solver.velocity_matrix().transpose() * velocity_cotangent
}

/// Selected circulation gradients, with unselected vectors left at zero.
#[must_use]
pub fn selected_force_circulation_gradients(
    solver: &VortexLatticeSolver,
    condition: &FlightCondition,
    selection: OutputSelection,
) -> ForceCoefficientGradients {
    let n = solver.mesh().vortex_count();
    let mut result = ForceCoefficientGradients {
        lift: VecX::zeros(n),
        drag: VecX::zeros(n),
        side: VecX::zeros(n),
        roll: VecX::zeros(n),
        pitch: VecX::zeros(n),
        yaw: VecX::zeros(n),
    };
    for output in [
        AeroOutput::Lift,
        AeroOutput::Drag,
        AeroOutput::Side,
        AeroOutput::Roll,
        AeroOutput::Pitch,
        AeroOutput::Yaw,
    ] {
        if !selection.contains(output) {
            continue;
        }
        let gradient = force_circulation_gradient(solver, condition, output);
        match output {
            AeroOutput::Lift => result.lift = gradient,
            AeroOutput::Drag => result.drag = gradient,
            AeroOutput::Side => result.side = gradient,
            AeroOutput::Roll => result.roll = gradient,
            AeroOutput::Pitch => result.pitch = gradient,
            AeroOutput::Yaw => result.yaw = gradient,
            AeroOutput::InducedDrag => {}
        }
    }
    result
}

/// Fixed-circulation reverse derivative with respect to the induced velocity
/// at each control point.
#[must_use]
pub fn force_induced_velocity_gradient(
    mesh: &VortexLattice,
    solution: &SolutionState,
    condition: &FlightCondition,
    reference: &ReferenceQuantities,
    output: AeroOutput,
) -> Vec<Vec3> {
    if !solution.is_valid || output == AeroOutput::InducedDrag {
        return vec![Vec3::zeros(); mesh.vortex_count()];
    }
    let (force_seed, moment_seed) = body_output_cotangent(output, condition.alpha);
    mesh.vortices()
        .iter()
        .enumerate()
        .map(|(index, vortex)| {
            let line = vortex.bound_vortex_right - vortex.bound_vortex_left;
            let arm = vortex.control_point - reference.moment_reference_point;
            let combined_seed = force_seed + moment_seed.cross(&arm) / reference.chord;
            2.0 * solution.circulation[index] / reference.area * line.cross(&combined_seed)
        })
        .collect()
}

/// Forward geometry derivative of all six primary coefficients while holding
/// circulation and induced velocity fixed.
///
/// Endpoint and control-point direction slices must contain exactly one vector
/// per vortex; a mismatch is an error rather than a silent truncation, which
/// would return a partial derivative indistinguishable from a correct one. The
/// control-point direction is also the moment-arm direction; callers moving the
/// reference point should subtract that direction first.
pub fn force_geometry_jvp(
    mesh: &VortexLattice,
    solution: &SolutionState,
    condition: &FlightCondition,
    reference: &ReferenceQuantities,
    bound_left_direction: &[Vec3],
    bound_right_direction: &[Vec3],
    control_point_direction: &[Vec3],
) -> Result<[Scalar; 6]> {
    let count = mesh.vortex_count();
    if bound_left_direction.len() != count
        || bound_right_direction.len() != count
        || control_point_direction.len() != count
        || solution.circulation.len() != count
    {
        return Err(PegasusError::InvalidGeometry(
            "geometry direction slices must have one entry per vortex".into(),
        ));
    }
    let freestream = VortexLatticeSolver::freestream_direction(condition);
    let mut force_derivative = Vec3::zeros();
    let mut moment_derivative = Vec3::zeros();
    for index in 0..count {
        let vortex = &mesh.vortices()[index];
        let induced = if solution.induced_velocity.nrows() > index {
            Vec3::new(
                solution.induced_velocity[(index, 0)],
                solution.induced_velocity[(index, 1)],
                solution.induced_velocity[(index, 2)],
            )
        } else {
            Vec3::zeros()
        };
        let velocity = freestream + induced;
        let line = vortex.bound_vortex_right - vortex.bound_vortex_left;
        let line_direction = bound_right_direction[index] - bound_left_direction[index];
        let panel_force = 2.0 * solution.circulation[index] * velocity.cross(&line);
        let panel_force_direction =
            2.0 * solution.circulation[index] * velocity.cross(&line_direction);
        let arm = vortex.control_point - reference.moment_reference_point;
        force_derivative += panel_force_direction;
        moment_derivative +=
            control_point_direction[index].cross(&panel_force) + arm.cross(&panel_force_direction);
    }
    force_derivative /= reference.area;
    moment_derivative /= reference.area * reference.chord;
    let (sin, cos) = condition.alpha.sin_cos();
    Ok([
        force_derivative.z * cos - force_derivative.x * sin,
        force_derivative.x * cos + force_derivative.z * sin,
        force_derivative.y,
        -moment_derivative.x,
        moment_derivative.y,
        -moment_derivative.z,
    ])
}

fn force_and_moment(
    mesh: &VortexLattice,
    solution: &SolutionState,
    reference: &ReferenceQuantities,
    velocity: Vec3,
) -> (Vec3, Vec3) {
    let mut force = Vec3::zeros();
    let mut moment = Vec3::zeros();
    for (index, vortex) in mesh.vortices().iter().enumerate() {
        let panel_force = 2.0
            * solution.circulation[index]
            * velocity.cross(&(vortex.bound_vortex_right - vortex.bound_vortex_left));
        force += panel_force;
        moment += (vortex.control_point - reference.moment_reference_point).cross(&panel_force);
    }
    (
        force / reference.area,
        moment / (reference.area * reference.chord),
    )
}

fn base_force_and_moment(
    mesh: &VortexLattice,
    solution: &SolutionState,
    reference: &ReferenceQuantities,
    freestream: Vec3,
) -> (Vec3, Vec3) {
    let mut force = Vec3::zeros();
    let mut moment = Vec3::zeros();
    for (index, vortex) in mesh.vortices().iter().enumerate() {
        let induced = if solution.induced_velocity.nrows() > index {
            Vec3::new(
                solution.induced_velocity[(index, 0)],
                solution.induced_velocity[(index, 1)],
                solution.induced_velocity[(index, 2)],
            )
        } else {
            Vec3::zeros()
        };
        let panel_force = 2.0
            * solution.circulation[index]
            * (freestream + induced).cross(&(vortex.bound_vortex_right - vortex.bound_vortex_left));
        force += panel_force;
        moment += (vortex.control_point - reference.moment_reference_point).cross(&panel_force);
    }
    (
        force / reference.area,
        moment / (reference.area * reference.chord),
    )
}

fn body_output_cotangent(output: AeroOutput, alpha: Scalar) -> (Vec3, Vec3) {
    let (sin, cos) = alpha.sin_cos();
    match output {
        AeroOutput::Lift => (Vec3::new(-sin, 0.0, cos), Vec3::zeros()),
        AeroOutput::Drag => (Vec3::new(cos, 0.0, sin), Vec3::zeros()),
        AeroOutput::Side => (Vec3::y(), Vec3::zeros()),
        AeroOutput::Roll => (Vec3::zeros(), Vec3::new(-1.0, 0.0, 0.0)),
        AeroOutput::Pitch => (Vec3::zeros(), Vec3::y()),
        AeroOutput::Yaw => (Vec3::zeros(), Vec3::new(0.0, 0.0, -1.0)),
        AeroOutput::InducedDrag => (Vec3::zeros(), Vec3::zeros()),
    }
}
