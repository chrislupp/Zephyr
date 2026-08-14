//! Ported from the frozen v0.2.0 `analytical_forces_test.cpp`.
//!
//! Near-field force and moment derivatives with respect to alpha, beta,
//! circulation, induced velocity and geometry.
//!
//! The C++ differenced its own `compute_force_coefficients_kernel`. The Rust
//! gradient functions take a lattice and do not expose that kernel, so the
//! primal is written out here instead — deliberately, since a difference taken
//! against the same code the gradient calls would not be an independent check.
//! It is the Kutta-Joukowski sum the library documents:
//!
//! ```text
//! F = sum 2 gamma_i (v_inf + v_ind,i) x (br_i - bl_i) / S
//! M = sum (rc_i - rmrp) x F_i / (S c)
//! ```
//!
//! rotated into stability axes as `[CL, CD, CY, Cl, Cm, Cn]`.
//!
//! The fixture reproduces the frozen arrays: two panels with bound vortices at
//! `x = 0.25` spanning `y` in `[0, 1]`, control points at `x = 0.75`, moment
//! reference at `(0.25, 0, 0)` so the arms come out as the frozen values.

use approx::assert_abs_diff_eq;
use pegasus_vlm::core::{PI, SpacingType, Vec3, VecX};
use pegasus_vlm::geometry::{Aircraft, ReferenceQuantities, Section, Surface};
use pegasus_vlm::gradients::{
    AeroOutput, force_alpha_gradient, force_beta_gradient, force_circulation_gradient,
    force_geometry_jvp, force_induced_velocity_gradient,
};
use pegasus_vlm::solver::{FlightCondition, SolutionState, VortexLattice, VortexLatticeSolver};
use pegasus_vlm::utils::body_to_stability_rotation;

const FD_STEP: f64 = 1.0e-7;
const FD_TOL: f64 = 1.0e-4;
const REFERENCE_AREA: f64 = 10.0;
const ALPHA: f64 = 5.0 * PI / 180.0;
const CIRCULATION: [f64; 2] = [0.1, 0.08];
const INDUCED: [[f64; 3]; 2] = [[0.0, 0.0, -0.005], [0.0, 0.0, -0.004]];

fn reference() -> ReferenceQuantities {
    ReferenceQuantities {
        area: REFERENCE_AREA,
        chord: 1.0,
        span: 1.0,
        moment_reference_point: Vec3::new(0.25, 0.0, 0.0),
    }
}

/// Two panels over `y` in `[0, 1]`, one chordwise.
fn lattice() -> VortexLattice {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(1);
    wing.set_spanwise_panel_count(2);
    wing.set_chordwise_spacing(SpacingType::Uniform);
    wing.set_spanwise_spacing(SpacingType::Uniform);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 1.0, 0.0), 1.0, 0.0));
    let mut aircraft = Aircraft::new("Force Fixture");
    aircraft.add_surface(wing);
    let mut mesh = VortexLattice::new();
    mesh.generate_from_aircraft(&aircraft).unwrap();
    mesh
}

fn solution(circulation: &[f64], induced: &[[f64; 3]]) -> SolutionState {
    let mut matrix = pegasus_vlm::core::MatX::zeros(induced.len(), 3);
    for (row, values) in induced.iter().enumerate() {
        for (column, value) in values.iter().enumerate() {
            matrix[(row, column)] = *value;
        }
    }
    SolutionState {
        circulation: VecX::from_column_slice(circulation),
        induced_velocity: matrix,
        is_valid: true,
        ..SolutionState::default()
    }
}

fn condition_at(alpha: f64, beta: f64) -> FlightCondition {
    FlightCondition {
        alpha,
        beta,
        ..FlightCondition::default()
    }
}

/// The six stability-axis coefficients, written out independently of the
/// library's gradient code.
fn primal(
    mesh: &VortexLattice,
    alpha: f64,
    beta: f64,
    circulation: &[f64],
    induced: &[[f64; 3]],
    bound_left: &[Vec3],
    bound_right: &[Vec3],
) -> [f64; 6] {
    let reference = reference();
    let freestream = VortexLatticeSolver::freestream_direction(&condition_at(alpha, beta));
    let mut force = Vec3::zeros();
    let mut moment = Vec3::zeros();
    for (index, vortex) in mesh.vortices().iter().enumerate() {
        let velocity =
            freestream + Vec3::new(induced[index][0], induced[index][1], induced[index][2]);
        let line = bound_right[index] - bound_left[index];
        let panel = 2.0 * circulation[index] * velocity.cross(&line);
        force += panel;
        moment += (vortex.control_point - reference.moment_reference_point).cross(&panel);
    }
    force /= reference.area;
    moment /= reference.area * reference.chord;

    let (sin_alpha, cos_alpha) = alpha.sin_cos();
    [
        force.z * cos_alpha - force.x * sin_alpha,
        force.x * cos_alpha + force.z * sin_alpha,
        force.y,
        -moment.x,
        moment.y,
        -moment.z,
    ]
}

fn endpoints(mesh: &VortexLattice) -> (Vec<Vec3>, Vec<Vec3>) {
    (
        mesh.vortices()
            .iter()
            .map(|v| v.bound_vortex_left)
            .collect(),
        mesh.vortices()
            .iter()
            .map(|v| v.bound_vortex_right)
            .collect(),
    )
}

fn assert_close(exact: f64, difference: f64) {
    assert_abs_diff_eq!(
        exact,
        difference,
        epsilon = FD_TOL * (difference.abs() + 1.0e-8)
    );
}

/// Ports `PrimalKernelComputesAllOutputs`.
#[test]
fn the_primal_produces_finite_positive_lift() {
    let mesh = lattice();
    let (left, right) = endpoints(&mesh);
    let result = primal(&mesh, ALPHA, 0.0, &CIRCULATION, &INDUCED, &left, &right);
    assert!(result[0] > 0.0);
    assert!(result.iter().all(|value| value.is_finite()));
}

/// Ports `FreestreamDirection`.
#[test]
fn the_freestream_direction_resolves_alpha() {
    let direction = VortexLatticeSolver::freestream_direction(&condition_at(0.0, 0.0));
    assert_abs_diff_eq!(direction, Vec3::new(1.0, 0.0, 0.0), epsilon = 1.0e-12);

    let angle = PI / 4.0;
    let direction = VortexLatticeSolver::freestream_direction(&condition_at(angle, 0.0));
    assert_abs_diff_eq!(direction.x, angle.cos(), epsilon = 1.0e-12);
    assert_abs_diff_eq!(direction.z, angle.sin(), epsilon = 1.0e-12);
}

/// Ports `FreestreamAlphaGradient`.
///
/// The closed form the force gradients use internally is
/// `(-sin a cos b, 0, cos a cos b)`; this differences the direction itself to
/// confirm it.
#[test]
fn the_freestream_alpha_derivative_matches_a_difference() {
    let beta = 0.0f64;
    let (sin_alpha, cos_alpha) = ALPHA.sin_cos();
    let (_, cos_beta) = beta.sin_cos();
    let exact = Vec3::new(-sin_alpha * cos_beta, 0.0, cos_alpha * cos_beta);

    let plus = VortexLatticeSolver::freestream_direction(&condition_at(ALPHA + FD_STEP, beta));
    let minus = VortexLatticeSolver::freestream_direction(&condition_at(ALPHA - FD_STEP, beta));
    let difference = (plus - minus) / (2.0 * FD_STEP);
    assert_abs_diff_eq!(exact, difference, epsilon = 1.0e-8);
}

/// Ports `FreestreamBetaGradient`.
#[test]
fn the_freestream_beta_derivative_matches_a_difference() {
    let beta = 3.0 * PI / 180.0;
    let (sin_alpha, cos_alpha) = ALPHA.sin_cos();
    let (sin_beta, cos_beta) = beta.sin_cos();
    let exact = Vec3::new(-cos_alpha * sin_beta, cos_beta, -sin_alpha * sin_beta);

    let plus = VortexLatticeSolver::freestream_direction(&condition_at(ALPHA, beta + FD_STEP));
    let minus = VortexLatticeSolver::freestream_direction(&condition_at(ALPHA, beta - FD_STEP));
    let difference = (plus - minus) / (2.0 * FD_STEP);
    assert_abs_diff_eq!(exact, difference, epsilon = 1.0e-8);
}

/// Ports `BodyToStability`.
///
/// At zero incidence the body and stability axes coincide, so body x is pure
/// drag and body z is pure lift.
#[test]
fn the_body_to_stability_rotation_is_identity_at_zero_alpha() {
    let rotation = body_to_stability_rotation(0.0);
    let drag_axis = rotation * Vec3::new(1.0, 0.0, 0.0);
    assert_abs_diff_eq!(drag_axis[0], 1.0, epsilon = 1.0e-12);
    assert_abs_diff_eq!(drag_axis[2], 0.0, epsilon = 1.0e-12);

    let lift_axis = rotation * Vec3::new(0.0, 0.0, 1.0);
    assert_abs_diff_eq!(lift_axis[0], 0.0, epsilon = 1.0e-12);
    assert_abs_diff_eq!(lift_axis[2], 1.0, epsilon = 1.0e-12);
}

/// Ports `AlphaGradientMatchesFD`.
///
/// All six outputs at once. Alpha enters twice — through the freestream vector
/// and through the body-to-stability rotation — so this is the case a
/// single-path derivative gets wrong.
#[test]
fn the_alpha_gradient_matches_finite_differences() {
    let mesh = lattice();
    let (left, right) = endpoints(&mesh);
    let state = solution(&CIRCULATION, &INDUCED);
    let exact = force_alpha_gradient(&mesh, &state, &condition_at(ALPHA, 0.0), &reference());

    let plus = primal(
        &mesh,
        ALPHA + FD_STEP,
        0.0,
        &CIRCULATION,
        &INDUCED,
        &left,
        &right,
    );
    let minus = primal(
        &mesh,
        ALPHA - FD_STEP,
        0.0,
        &CIRCULATION,
        &INDUCED,
        &left,
        &right,
    );
    for ((exact, plus), minus) in exact.iter().zip(&plus).zip(&minus) {
        assert_close(*exact, (plus - minus) / (2.0 * FD_STEP));
    }
}

/// Ports `BetaGradientMatchesFD`.
#[test]
fn the_beta_gradient_matches_finite_differences() {
    let beta = 3.0 * PI / 180.0;
    let mesh = lattice();
    let (left, right) = endpoints(&mesh);
    let state = solution(&CIRCULATION, &INDUCED);
    let exact = force_beta_gradient(&mesh, &state, &condition_at(ALPHA, beta), &reference());

    let plus = primal(
        &mesh,
        ALPHA,
        beta + FD_STEP,
        &CIRCULATION,
        &INDUCED,
        &left,
        &right,
    );
    let minus = primal(
        &mesh,
        ALPHA,
        beta - FD_STEP,
        &CIRCULATION,
        &INDUCED,
        &left,
        &right,
    );
    for ((exact, plus), minus) in exact.iter().zip(&plus).zip(&minus) {
        assert_close(*exact, (plus - minus) / (2.0 * FD_STEP));
    }
}

/// Ports `CirculationGradientMatchesFD`, over all six outputs and both panels.
///
/// The frozen kernel took induced velocity as an independent input, so its
/// circulation gradient was a partial derivative at fixed induced velocity.
/// `force_circulation_gradient` is documented as the *total* partial including
/// the induced-velocity dependence, and ends with `direct + V^T c`, because in
/// a real solve the induced velocity is `V Gamma` and moves with the
/// circulation. The difference here therefore recomputes induced velocity from
/// the velocity influence matrix at each perturbed circulation. Holding it
/// fixed instead disagrees by about 0.8% on this fixture — small, systematic,
/// and exactly the induced contribution.
#[test]
fn the_circulation_gradient_matches_finite_differences() {
    let mut aircraft = Aircraft::new("Force Fixture");
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(1);
    wing.set_spanwise_panel_count(2);
    wing.set_chordwise_spacing(SpacingType::Uniform);
    wing.set_spanwise_spacing(SpacingType::Uniform);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 1.0, 0.0), 1.0, 0.0));
    aircraft.add_surface(wing);
    aircraft.set_reference_area(REFERENCE_AREA);
    aircraft.set_reference_chord(1.0);
    aircraft.set_reference_span(1.0);
    aircraft.set_moment_reference_point(Vec3::new(0.25, 0.0, 0.0));

    let condition = condition_at(ALPHA, 0.0);
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    solver.solve(&condition).unwrap();

    let mesh = solver.mesh().clone();
    let (left, right) = endpoints(&mesh);
    let solved: Vec<f64> = solver.solution().circulation.iter().copied().collect();
    let velocity_matrix = solver.velocity_matrix().clone();

    // Induced velocity implied by a circulation, v = V Gamma.
    let induced_for = |circulation: &[f64]| -> Vec<[f64; 3]> {
        let gamma = VecX::from_column_slice(circulation);
        let flat = &velocity_matrix * gamma;
        (0..circulation.len())
            .map(|i| [flat[3 * i], flat[3 * i + 1], flat[3 * i + 2]])
            .collect()
    };

    for (index, output) in [
        (0, AeroOutput::Lift),
        (1, AeroOutput::Drag),
        (2, AeroOutput::Side),
        (3, AeroOutput::Roll),
        (4, AeroOutput::Pitch),
        (5, AeroOutput::Yaw),
    ] {
        let exact = force_circulation_gradient(&solver, &condition, output);
        for panel in 0..solved.len() {
            let mut plus = solved.clone();
            let mut minus = solved.clone();
            plus[panel] += FD_STEP;
            minus[panel] -= FD_STEP;
            let difference = (primal(&mesh, ALPHA, 0.0, &plus, &induced_for(&plus), &left, &right)
                [index]
                - primal(
                    &mesh,
                    ALPHA,
                    0.0,
                    &minus,
                    &induced_for(&minus),
                    &left,
                    &right,
                )[index])
                / (2.0 * FD_STEP);
            assert_close(exact[panel], difference);
        }
    }
}

/// Ports `InducedVelocityGradientMatchesFD`, over all six outputs and all three
/// components of both panels.
#[test]
fn the_induced_velocity_gradient_matches_finite_differences() {
    let mesh = lattice();
    let (left, right) = endpoints(&mesh);
    let state = solution(&CIRCULATION, &INDUCED);
    let condition = condition_at(ALPHA, 0.0);

    for (index, output) in [
        (0, AeroOutput::Lift),
        (1, AeroOutput::Drag),
        (2, AeroOutput::Side),
        (3, AeroOutput::Roll),
        (4, AeroOutput::Pitch),
        (5, AeroOutput::Yaw),
    ] {
        let exact =
            force_induced_velocity_gradient(&mesh, &state, &condition, &reference(), output);
        for panel in 0..2 {
            for component in 0..3 {
                let mut plus = INDUCED;
                let mut minus = INDUCED;
                plus[panel][component] += FD_STEP;
                minus[panel][component] -= FD_STEP;
                let difference = (primal(&mesh, ALPHA, 0.0, &CIRCULATION, &plus, &left, &right)
                    [index]
                    - primal(&mesh, ALPHA, 0.0, &CIRCULATION, &minus, &left, &right)[index])
                    / (2.0 * FD_STEP);
                assert_close(exact[panel][component], difference);
            }
        }
    }
}

/// Ports `GeometryGradientMatchesFD`.
///
/// The frozen case seeds a unit perturbation of the first panel's left endpoint
/// in `y` and differences all six outputs along it.
#[test]
fn the_geometry_jvp_matches_a_directional_difference() {
    let mesh = lattice();
    let (left, right) = endpoints(&mesh);
    let state = solution(&CIRCULATION, &INDUCED);
    let condition = condition_at(ALPHA, 0.0);

    let mut left_direction = vec![Vec3::zeros(); 2];
    left_direction[0].y = 1.0;
    let right_direction = vec![Vec3::zeros(); 2];
    let control_direction = vec![Vec3::zeros(); 2];

    let exact = force_geometry_jvp(
        &mesh,
        &state,
        &condition,
        &reference(),
        &left_direction,
        &right_direction,
        &control_direction,
    )
    .unwrap();

    let mut plus = left.clone();
    let mut minus = left.clone();
    plus[0].y += FD_STEP;
    minus[0].y -= FD_STEP;
    let forward = primal(&mesh, ALPHA, 0.0, &CIRCULATION, &INDUCED, &plus, &right);
    let backward = primal(&mesh, ALPHA, 0.0, &CIRCULATION, &INDUCED, &minus, &right);
    for ((exact, forward), backward) in exact.iter().zip(&forward).zip(&backward) {
        assert_close(*exact, (forward - backward) / (2.0 * FD_STEP));
    }
}

/// Ports `ZeroAlpha`.
#[test]
fn the_alpha_gradient_stays_finite_at_zero_incidence() {
    let mesh = lattice();
    let state = solution(&CIRCULATION, &INDUCED);
    let exact = force_alpha_gradient(&mesh, &state, &condition_at(0.0, 0.0), &reference());
    assert!(exact.iter().all(|value| value.is_finite()));
}
