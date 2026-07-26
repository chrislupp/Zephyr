mod common;

use approx::assert_abs_diff_eq;
use pegasus_vlm::aero::ForceCalculator;
use pegasus_vlm::core::{DEGREES_TO_RADIANS, Vec3, VecX};
use pegasus_vlm::gradients::{
    AeroOutput, GeometryVariable, GradientCalculator, OutputSelection, TotalDerivative,
    aic_element_gradient, horseshoe_velocity_pg_derivative, horseshoe_velocity_vjp_bound_left,
    horseshoe_velocity_vjp_bound_right, horseshoe_velocity_vjp_field_point,
    trefftz_circulation_gradients,
};
use pegasus_vlm::solver::{
    FlightCondition, SolutionState, VortexLatticeSolver, horseshoe_vortex_velocity,
};

use common::{rectangular_mesh, rectangular_wing};

fn coefficient(
    aircraft: &pegasus_vlm::geometry::Aircraft,
    condition: &FlightCondition,
    output: AeroOutput,
) -> f64 {
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(aircraft).unwrap();
    let solution = solver.solve(condition).unwrap();
    let result =
        ForceCalculator::new().compute_for_aircraft(aircraft, solver.mesh(), &solution, condition);
    match output {
        AeroOutput::Lift => result.coefficients.lift,
        AeroOutput::Drag => result.coefficients.drag,
        AeroOutput::Side => result.coefficients.side,
        AeroOutput::Roll => result.coefficients.roll,
        AeroOutput::Pitch => result.coefficients.pitch,
        AeroOutput::Yaw => result.coefficients.yaw,
        AeroOutput::InducedDrag => result.far_field_induced_drag,
    }
}

fn mesh_coefficient(
    mesh: &pegasus_vlm::geometry::MeshInput,
    condition: &FlightCondition,
    output: AeroOutput,
) -> f64 {
    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(mesh).unwrap();
    let solution = solver.solve(condition).unwrap();
    let result = ForceCalculator::new().compute_with_references(
        &mesh.reference,
        solver.mesh(),
        &solution,
        condition,
    );
    match output {
        AeroOutput::Lift => result.coefficients.lift,
        AeroOutput::Drag => result.coefficients.drag,
        AeroOutput::Side => result.coefficients.side,
        AeroOutput::Roll => result.coefficients.roll,
        AeroOutput::Pitch => result.coefficients.pitch,
        AeroOutput::Yaw => result.coefficients.yaw,
        AeroOutput::InducedDrag => result.far_field_induced_drag,
    }
}

#[test]
fn biot_savart_vjps_match_centered_differences() {
    let field = Vec3::new(0.72, 0.31, 0.18);
    let left = Vec3::new(0.25, -0.4, 0.02);
    let right = Vec3::new(0.28, 0.55, -0.01);
    let seed = Vec3::new(0.4, -0.7, 1.2);
    let step = 1.0e-7;
    let analytical = [
        horseshoe_velocity_vjp_field_point(&field, &left, &right, 0.92, 0.01, true, &seed),
        horseshoe_velocity_vjp_bound_left(&field, &left, &right, 0.92, 0.01, true, &seed),
        horseshoe_velocity_vjp_bound_right(&field, &left, &right, 0.92, 0.01, true, &seed),
    ];
    for (which, base) in [field, left, right].into_iter().enumerate() {
        for component in 0..3 {
            let mut plus = [field, left, right];
            let mut minus = plus;
            plus[which] = base;
            minus[which] = base;
            plus[which][component] += step;
            minus[which][component] -= step;
            let value = |points: [Vec3; 3]| {
                horseshoe_vortex_velocity(&points[0], &points[1], &points[2], 0.92, 0.01, true)
                    .dot(&seed)
            };
            let finite = (value(plus) - value(minus)) / (2.0 * step);
            assert_abs_diff_eq!(analytical[which][component], finite, epsilon = 2.0e-7);
        }
    }
}

#[test]
fn biot_savart_pg_derivative_matches_centered_difference() {
    let field = Vec3::new(0.72, 0.31, 0.18);
    let left = Vec3::new(0.25, -0.4, 0.02);
    let right = Vec3::new(0.28, 0.55, -0.01);
    let pg = 0.91;
    let step = 1.0e-7;
    let analytical = horseshoe_velocity_pg_derivative(&field, &left, &right, pg, 0.01, true);
    let finite = (horseshoe_vortex_velocity(&field, &left, &right, pg + step, 0.01, true)
        - horseshoe_vortex_velocity(&field, &left, &right, pg - step, 0.01, true))
        / (2.0 * step);
    assert_abs_diff_eq!(analytical, finite, epsilon = 3.0e-8);
}

#[test]
fn aic_element_gradient_matches_centered_differences() {
    let field = Vec3::new(0.72, 0.31, 0.18);
    let left = Vec3::new(0.25, -0.4, 0.02);
    let right = Vec3::new(0.28, 0.55, -0.01);
    let normal = Vec3::new(0.1, -0.2, 0.97).normalize();
    let gradient = aic_element_gradient(&left, &right, &field, &normal, 0.02, 0.9, true);
    let step = 1.0e-7;
    for (which, analytical) in [
        gradient.bound_left,
        gradient.bound_right,
        gradient.control_point,
        gradient.normal,
    ]
    .into_iter()
    .enumerate()
    {
        for component in 0..3 {
            let mut plus = [left, right, field, normal];
            let mut minus = plus;
            plus[which][component] += step;
            minus[which][component] -= step;
            let value = |points: [Vec3; 4]| {
                horseshoe_vortex_velocity(&points[2], &points[0], &points[1], 0.9, 0.02, true)
                    .dot(&points[3])
            };
            let finite = (value(plus) - value(minus)) / (2.0 * step);
            assert_abs_diff_eq!(analytical[component], finite, epsilon = 3.0e-7);
        }
    }
}

#[test]
fn trefftz_circulation_gradient_matches_direct_perturbations() {
    let aircraft = rectangular_wing(3, 6);
    let condition = FlightCondition {
        alpha: 5.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    let solution = solver.solve(&condition).unwrap();
    let gradient =
        trefftz_circulation_gradients(solver.mesh(), &solution, 4.0, 0.01, 0.1).induced_drag;
    let step = 1.0e-7;
    for index in [
        0,
        solution.circulation.len() / 2,
        solution.circulation.len() - 1,
    ] {
        let mut plus = solution.clone();
        let mut minus = solution.clone();
        plus.circulation[index] += step;
        minus.circulation[index] -= step;
        let calculator = ForceCalculator::new();
        let plus_value = calculator.trefftz_plane(solver.mesh(), &plus, 4.0, 4.0).2;
        let minus_value = calculator.trefftz_plane(solver.mesh(), &minus, 4.0, 4.0).2;
        assert_abs_diff_eq!(
            gradient[index],
            (plus_value - minus_value) / (2.0 * step),
            epsilon = 1.0e-8
        );
    }
}

#[test]
fn adjoint_induced_drag_matches_full_resolves() {
    let aircraft = rectangular_wing(3, 6);
    let condition = FlightCondition {
        alpha: 5.0 * DEGREES_TO_RADIANS,
        beta: 1.5 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    solver.solve(&condition).unwrap();
    let derivative =
        TotalDerivative::compute(&solver, &condition, AeroOutput::InducedDrag).unwrap();
    let step = 1.0e-6;
    let mut plus = condition.clone();
    let mut minus = condition.clone();
    plus.alpha += step;
    minus.alpha -= step;
    let finite_alpha = (coefficient(&aircraft, &plus, AeroOutput::InducedDrag)
        - coefficient(&aircraft, &minus, AeroOutput::InducedDrag))
        / (2.0 * step);
    assert_abs_diff_eq!(derivative.by_alpha, finite_alpha, epsilon = 2.0e-6);
}

#[test]
fn gradient_calculator_flight_derivatives_match_resolves() {
    let aircraft = rectangular_wing(3, 6);
    let condition = FlightCondition {
        alpha: 4.0 * DEGREES_TO_RADIANS,
        beta: 2.0 * DEGREES_TO_RADIANS,
        mach: 0.2,
        ..FlightCondition::default()
    };
    let mut calculator = GradientCalculator::new();
    calculator.set_aircraft(&aircraft);
    let gradients = calculator
        .compute_flight_condition_gradients(&condition)
        .unwrap();
    let step = 1.0e-6;
    for output in [
        AeroOutput::Lift,
        AeroOutput::Drag,
        AeroOutput::Side,
        AeroOutput::Roll,
        AeroOutput::Pitch,
        AeroOutput::Yaw,
        AeroOutput::InducedDrag,
    ] {
        let mut plus = condition.clone();
        let mut minus = condition.clone();
        plus.alpha += step;
        minus.alpha -= step;
        let finite = (coefficient(&aircraft, &plus, output)
            - coefficient(&aircraft, &minus, output))
            / (2.0 * step);
        assert_abs_diff_eq!(gradients.by_alpha.get(output), finite, epsilon = 3.0e-5);

        plus = condition.clone();
        minus = condition.clone();
        plus.mach += step;
        minus.mach -= step;
        let finite_mach = (coefficient(&aircraft, &plus, output)
            - coefficient(&aircraft, &minus, output))
            / (2.0 * step);
        assert_abs_diff_eq!(gradients.by_mach.get(output), finite_mach, epsilon = 3.0e-5);
    }
}

#[test]
fn output_selection_leaves_disabled_outputs_exactly_zero() {
    let aircraft = rectangular_wing(2, 4);
    let mut calculator = GradientCalculator::new();
    calculator.set_aircraft(&aircraft);
    calculator.set_output_selection(OutputSelection::only(AeroOutput::Pitch));
    let gradients = calculator
        .compute_flight_condition_gradients(&FlightCondition::default())
        .unwrap();
    assert_eq!(gradients.by_alpha.coefficients.lift, 0.0);
    assert_eq!(gradients.by_alpha.coefficients.drag, 0.0);
    assert_eq!(gradients.by_alpha.induced_drag, 0.0);
    assert!(gradients.by_alpha.coefficients.pitch.is_finite());
}

#[test]
fn section_geometry_gradient_has_expected_shape_and_nonzero_chord_effect() {
    let aircraft = rectangular_wing(2, 4);
    let mut calculator = GradientCalculator::new();
    calculator.set_aircraft(&aircraft);
    calculator.set_output_selection(OutputSelection::only(AeroOutput::Lift));
    let condition = FlightCondition {
        alpha: 5.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let gradients = calculator.compute_geometry_gradients(&condition).unwrap();
    assert_eq!(gradients.sections.len(), 2);
    assert_eq!(gradients.design_variable_count(), 10);
    assert!(
        gradients.sections[0]
            .get(GeometryVariable::Chord)
            .coefficients
            .lift
            .abs()
            > 1.0e-6
    );
    assert_eq!(gradients.sections[0].by_chord.coefficients.drag, 0.0);
}

#[test]
fn direct_mesh_node_gradient_matches_external_resolve() {
    let mesh = rectangular_mesh(1, 2);
    let mut calculator = GradientCalculator::new();
    calculator.set_mesh(&mesh).unwrap();
    calculator.set_output_selection(OutputSelection::only(AeroOutput::Lift));
    let condition = FlightCondition {
        alpha: 5.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let gradients = calculator.compute_mesh_node_gradients(&condition).unwrap();
    assert_eq!(gradients.len(), 1);
    assert_eq!(gradients[0].lift.nodes.len(), 6);
    assert_eq!(gradients[0].drag.nodes, vec![Vec3::zeros(); 6]);

    let step = 1.0e-6;
    let mut plus = mesh.clone();
    let mut minus = mesh.clone();
    plus.surfaces[0].nodes[0].z += step;
    minus.surfaces[0].nodes[0].z -= step;
    let solve = |input: &pegasus_vlm::geometry::MeshInput| {
        let mut solver = VortexLatticeSolver::new();
        solver.set_mesh(input).unwrap();
        let solution = solver.solve(&condition).unwrap();
        ForceCalculator::new()
            .compute_with_references(&input.reference, solver.mesh(), &solution, &condition)
            .coefficients
            .lift
    };
    assert_abs_diff_eq!(
        gradients[0].lift.nodes[0].z,
        (solve(&plus) - solve(&minus)) / (2.0 * step),
        epsilon = 1.0e-7
    );
}

#[test]
fn direct_mesh_node_adjoint_matches_all_outputs() {
    let mesh = rectangular_mesh(2, 3);
    let condition = FlightCondition {
        alpha: 4.0 * DEGREES_TO_RADIANS,
        beta: 1.5 * DEGREES_TO_RADIANS,
        roll_rate: 0.03,
        pitch_rate: -0.02,
        yaw_rate: 0.01,
        ..FlightCondition::default()
    };
    let mut calculator = GradientCalculator::new();
    calculator.set_mesh(&mesh).unwrap();
    let gradients = calculator.compute_mesh_node_gradients(&condition).unwrap();
    let node_i = 1;
    let node_j = 1;
    let node_index = node_i * (mesh.surfaces[0].ny + 1) + node_j;
    let component = 2;
    let step = 1.0e-6;
    let mut plus = mesh.clone();
    let mut minus = mesh.clone();
    plus.surfaces[0].nodes[node_index][component] += step;
    minus.surfaces[0].nodes[node_index][component] -= step;
    for (output, analytical) in [
        (
            AeroOutput::Lift,
            gradients[0].lift.nodes[node_index][component],
        ),
        (
            AeroOutput::Drag,
            gradients[0].drag.nodes[node_index][component],
        ),
        (
            AeroOutput::Side,
            gradients[0].side.nodes[node_index][component],
        ),
        (
            AeroOutput::Roll,
            gradients[0].roll.nodes[node_index][component],
        ),
        (
            AeroOutput::Pitch,
            gradients[0].pitch.nodes[node_index][component],
        ),
        (
            AeroOutput::Yaw,
            gradients[0].yaw.nodes[node_index][component],
        ),
        (
            AeroOutput::InducedDrag,
            gradients[0].induced_drag.nodes[node_index][component],
        ),
    ] {
        let finite = (mesh_coefficient(&plus, &condition, output)
            - mesh_coefficient(&minus, &condition, output))
            / (2.0 * step);
        assert_abs_diff_eq!(analytical, finite, epsilon = 2.0e-6);
    }
}

#[test]
fn adjoint_batched_transpose_identity() {
    let aircraft = rectangular_wing(2, 4);
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    solver.solve(&FlightCondition::default()).unwrap();
    let rhs = VecX::from_element(solver.mesh().vortex_count(), 0.25);
    let lambda = pegasus_vlm::gradients::AdjointSolver::new(&solver)
        .unwrap()
        .solve(&rhs)
        .unwrap();
    assert_abs_diff_eq!(
        solver.linearized_aic().transpose() * lambda,
        rhs,
        epsilon = 1.0e-10
    );
}

#[test]
fn invalid_solution_is_rejected_by_adjoint() {
    let solver = VortexLatticeSolver::new();
    assert!(
        TotalDerivative::compute(&solver, &FlightCondition::default(), AeroOutput::Lift).is_err()
    );
    let state = SolutionState::default();
    assert!(!state.is_valid);
}
