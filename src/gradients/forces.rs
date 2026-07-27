//! Analytical derivatives of the Kutta–Joukowski force kernel.

use crate::aero::ForceCalculator;
use crate::aero::forces::drag_polar_value_and_slope;
use crate::core::{EPSILON, Scalar, Vec3, VecX};
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
    let (mut base_force, _) = base_force_and_moment(mesh, solution, reference, freestream);
    let (mut force_derivative, mut moment_derivative) =
        force_and_moment(mesh, solution, reference, derivative);
    // Section drag polars act along the stability drag axis, which both rotates
    // with alpha and carries a polar whose argument moves with alpha.
    let viscous = viscous_strips(mesh, solution, condition, reference);
    if !viscous.is_empty() {
        let (viscous_force, _) = viscous_force_and_moment(&viscous, condition.alpha, reference);
        base_force += viscous_force;
        let (viscous_force_derivative, viscous_moment_derivative) =
            viscous_force_and_moment_derivative(
                mesh,
                solution,
                &viscous,
                condition.alpha,
                reference,
                derivative,
                true,
            );
        force_derivative += viscous_force_derivative;
        moment_derivative += viscous_moment_derivative;
    }
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
    let (mut force, mut moment) = force_and_moment(mesh, solution, reference, derivative);
    // Sideslip moves the strip lift coefficients, and so the polars, but the
    // stability axes themselves do not turn with beta.
    let viscous = viscous_strips(mesh, solution, condition, reference);
    if !viscous.is_empty() {
        let (viscous_force, viscous_moment) = viscous_force_and_moment_derivative(
            mesh,
            solution,
            &viscous,
            condition.alpha,
            reference,
            derivative,
            false,
        );
        force += viscous_force;
        moment += viscous_moment;
    }
    [
        force.z * cos_alpha - force.x * sin_alpha,
        force.x * cos_alpha + force.z * sin_alpha,
        force.y,
        -moment.x,
        moment.y,
        -moment.z,
    ]
}

/// Primal viscous state of one strip, cached for the derivative kernels.
struct ViscousStrip {
    first_vortex_index: usize,
    num_vortices: usize,
    area: Scalar,
    /// Moment arm of the strip's viscous force, `quarter_chord - reference`.
    arm: Vec3,
    /// `cd_v` after the primal's `max(0.0)` clamp.
    viscous_drag: Scalar,
    /// `dCd/dCl` at the strip's lift coefficient, zero where the clamp is
    /// active, since the primal is flat there.
    slope: Scalar,
    /// Strip drag coefficient *excluding* the viscous term.
    inviscid_drag: Scalar,
}

/// Per-strip viscous state, or an empty vector when no strip carries a polar.
///
/// Runs the primal force calculation rather than re-deriving the strip lift
/// coefficient, so the polar cannot be evaluated at a different point here than
/// it is in `strip_loads`.
fn viscous_strips(
    mesh: &VortexLattice,
    solution: &SolutionState,
    condition: &FlightCondition,
    reference: &ReferenceQuantities,
) -> Vec<ViscousStrip> {
    if !solution.is_valid || !mesh.strips().iter().any(|strip| strip.has_viscous_data) {
        return Vec::new();
    }
    let results =
        ForceCalculator::new().compute_with_references(reference, mesh, solution, condition);
    let mut viscous = Vec::new();
    for surface in &results.surface_loads {
        for loads in &surface.strip_loads {
            let strip = &mesh.strips()[loads.strip_index];
            if !strip.has_viscous_data || loads.area <= EPSILON {
                continue;
            }
            let slope = if loads.viscous_drag > 0.0 {
                drag_polar_value_and_slope(&strip.viscous_polar, loads.coefficients.lift).1
            } else {
                0.0
            };
            viscous.push(ViscousStrip {
                first_vortex_index: strip.first_vortex_index,
                num_vortices: strip.num_vortices,
                area: loads.area,
                arm: loads.quarter_chord_position - reference.moment_reference_point,
                viscous_drag: loads.viscous_drag,
                slope,
                inviscid_drag: loads.coefficients.drag - loads.viscous_drag,
            });
        }
    }
    viscous
}

/// Viscous contribution to the aircraft body force and moment, normalized the
/// same way as `force_and_moment`.
fn viscous_force_and_moment(
    strips: &[ViscousStrip],
    alpha: Scalar,
    reference: &ReferenceQuantities,
) -> (Vec3, Vec3) {
    let (sin_alpha, cos_alpha) = alpha.sin_cos();
    let drag_axis = Vec3::new(cos_alpha, 0.0, sin_alpha);
    let mut force = Vec3::zeros();
    let mut moment = Vec3::zeros();
    for strip in strips {
        let contribution = strip.area * strip.viscous_drag * drag_axis;
        force += contribution;
        moment += strip.arm.cross(&contribution);
    }
    (
        force / reference.area,
        moment / (reference.area * reference.chord),
    )
}

/// Fixed-circulation derivative of the viscous force and moment along a
/// freestream direction change.
///
/// `alpha_rotation` additionally accounts for the stability axes turning with
/// the parameter, which applies to alpha and not to sideslip: the drag axis the
/// viscous force acts along rotates, and so does the lift axis the polar is
/// evaluated against.
fn viscous_force_and_moment_derivative(
    mesh: &VortexLattice,
    solution: &SolutionState,
    strips: &[ViscousStrip],
    alpha: Scalar,
    reference: &ReferenceQuantities,
    velocity_derivative: Vec3,
    alpha_rotation: bool,
) -> (Vec3, Vec3) {
    let (sin_alpha, cos_alpha) = alpha.sin_cos();
    let drag_axis = Vec3::new(cos_alpha, 0.0, sin_alpha);
    let lift_axis = Vec3::new(-sin_alpha, 0.0, cos_alpha);
    let mut force = Vec3::zeros();
    let mut moment = Vec3::zeros();
    for strip in strips {
        // The strip force changes because the freestream direction does.
        let mut strip_force_derivative = Vec3::zeros();
        for offset in 0..strip.num_vortices {
            let index = strip.first_vortex_index + offset;
            let vortex = &mesh.vortices()[index];
            let line = vortex.bound_vortex_right - vortex.bound_vortex_left;
            strip_force_derivative +=
                2.0 * solution.circulation[index] * velocity_derivative.cross(&line);
        }
        // `d(lift_axis)/d(alpha) = -drag_axis`, so its contraction with the
        // strip force is exactly minus the strip's inviscid drag coefficient.
        let axis_term = if alpha_rotation {
            -strip.inviscid_drag
        } else {
            0.0
        };
        let lift_derivative = axis_term + lift_axis.dot(&strip_force_derivative) / strip.area;
        let drag_derivative = strip.slope * lift_derivative;

        let mut contribution = strip.area * drag_derivative * drag_axis;
        if alpha_rotation {
            // `d(drag_axis)/d(alpha) = lift_axis`.
            contribution += strip.area * strip.viscous_drag * lift_axis;
        }
        force += contribution;
        moment += strip.arm.cross(&contribution);
    }
    (
        force / reference.area,
        moment / (reference.area * reference.chord),
    )
}

/// Per-strip sensitivity of one output to that strip's viscous drag
/// coefficient, divided by the strip lift coefficient's own scaling.
///
/// `strip_loads` adds `cd_v * (cos a, 0, sin a)` to the strip's *body force*,
/// so the term reaches the totals as `area_s * cd_v,s * u` and the moments as
/// `arm_s x (area_s * cd_v,s * u)`. Contracting with the output cotangent,
///
/// ```text
/// dJ/d(cd_v,s) = (area_s / S_ref) * (force_seed + moment_seed x arm_s / c_ref) . u
/// ```
///
/// and since `cl_s` is itself `s_hat . F_s / area_s`, the areas cancel in the
/// chain and never appear in the returned weight.
///
/// The entry is zero for strips with no polar and for strips where the
/// quadratic went negative, since the primal clamps that branch to zero and
/// its subderivative there is zero. Returns `None` when no strip is viscous,
/// so the common inviscid path costs nothing.
fn viscous_strip_weights(
    mesh: &VortexLattice,
    solution: &SolutionState,
    condition: &FlightCondition,
    reference: &ReferenceQuantities,
    output: AeroOutput,
) -> Option<Vec<Scalar>> {
    if !mesh.strips().iter().any(|strip| strip.has_viscous_data) {
        return None;
    }
    let (force_seed, moment_seed) = body_output_cotangent(output, condition.alpha);
    let (sin_alpha, cos_alpha) = condition.alpha.sin_cos();
    let drag_axis = Vec3::new(cos_alpha, 0.0, sin_alpha);

    // Reuse the primal rather than re-deriving the strip lift coefficient, so
    // the two cannot disagree about what the polar is evaluated at.
    let results =
        ForceCalculator::new().compute_with_references(reference, mesh, solution, condition);
    let mut weights = vec![0.0; mesh.strip_count()];
    for surface in &results.surface_loads {
        for loads in &surface.strip_loads {
            let strip = &mesh.strips()[loads.strip_index];
            if !strip.has_viscous_data || loads.viscous_drag <= 0.0 {
                continue;
            }
            let (_, slope) =
                drag_polar_value_and_slope(&strip.viscous_polar, loads.coefficients.lift);
            let arm = loads.quarter_chord_position - reference.moment_reference_point;
            let combined = force_seed + moment_seed.cross(&arm) / reference.chord;
            weights[loads.strip_index] = combined.dot(&drag_axis) * slope / reference.area;
        }
    }
    Some(weights)
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

    // Section drag polars make the strip's viscous drag a function of that
    // strip's lift coefficient, and hence of every circulation in the mesh
    // through the induced velocity. Omitting this left drag, roll, pitch and
    // yaw wrong by the polar slope whenever a `CDCL` block was present.
    if let Some(weights) = viscous_strip_weights(mesh, solution, condition, reference, output) {
        let lift_axis = Vec3::new(-condition.alpha.sin(), 0.0, condition.alpha.cos());
        for (strip_index, strip) in mesh.strips().iter().enumerate() {
            let weight = weights[strip_index];
            if weight == 0.0 {
                continue;
            }
            if strip.chord_center * strip.width <= EPSILON {
                continue;
            }
            // `cl_s = lift_axis . F_s / area_s`, so each panel in the strip
            // contributes both a direct Kutta-Joukowski term and, through its
            // own induced velocity, a non-local one. The strip area already
            // cancelled when the weight was formed and must not be applied
            // again here.
            let scale = weight;
            for offset in 0..strip.num_vortices {
                let i = strip.first_vortex_index + offset;
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
                direct[i] += 2.0 * scale * lift_axis.dot(&velocity.cross(&line));
                let cotangent = 2.0 * solution.circulation[i] * scale * line.cross(&lift_axis);
                velocity_cotangent[3 * i] += cotangent.x;
                velocity_cotangent[3 * i + 1] += cotangent.y;
                velocity_cotangent[3 * i + 2] += cotangent.z;
            }
        }
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
