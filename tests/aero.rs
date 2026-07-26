mod common;

use approx::assert_abs_diff_eq;
use pegasus_vlm::aero::{ForceCalculator, StabilityDerivativeCalculator};
use pegasus_vlm::core::{DEGREES_TO_RADIANS, Vec3};
use pegasus_vlm::solver::{FlightCondition, VortexLatticeSolver};

use common::{controlled_wing, rectangular_wing};

#[test]
fn named_control_deflection_changes_the_primal_solution() {
    let aircraft = controlled_wing();
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    let base = solver.solve(&FlightCondition::default()).unwrap();
    let controlled = solver
        .solve(&FlightCondition {
            control_deflections: vec![("elevator".into(), 5.0 * DEGREES_TO_RADIANS)],
            ..FlightCondition::default()
        })
        .unwrap();
    assert!((&controlled.circulation - &base.circulation).norm() > 1.0e-5);
    assert_eq!(
        solver.solution_aic().shape(),
        (solver.mesh().vortex_count(), solver.mesh().vortex_count())
    );
}

#[test]
fn analytical_control_circulation_derivative_matches_resolves() {
    let aircraft = controlled_wing();
    let condition = FlightCondition {
        alpha: 3.0 * DEGREES_TO_RADIANS,
        control_deflections: vec![("elevator".into(), 2.0 * DEGREES_TO_RADIANS)],
        ..FlightCondition::default()
    };
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    solver.solve(&condition).unwrap();
    let analytical = solver
        .control_circulation_derivative(&condition, "elevator")
        .unwrap();
    let step = 1.0e-6;
    let mut plus = condition.clone();
    let mut minus = condition.clone();
    plus.control_deflections[0].1 += step;
    minus.control_deflections[0].1 -= step;
    let mut plus_solver = VortexLatticeSolver::new();
    let mut minus_solver = VortexLatticeSolver::new();
    plus_solver.set_aircraft(&aircraft).unwrap();
    minus_solver.set_aircraft(&aircraft).unwrap();
    let plus_solution = plus_solver.solve(&plus).unwrap();
    let minus_solution = minus_solver.solve(&minus).unwrap();
    let finite = (plus_solution.circulation - minus_solution.circulation) / (2.0 * step);
    assert_abs_diff_eq!(analytical, finite, epsilon = 2.0e-8);
}

#[test]
fn absent_control_entry_has_a_zero_deflection_derivative() {
    let aircraft = controlled_wing();
    let condition = FlightCondition {
        alpha: 3.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    solver.solve(&condition).unwrap();
    let derivative = solver
        .control_circulation_derivative(&condition, "elevator")
        .unwrap();
    assert!(derivative.norm() > 1.0e-5);
}

#[test]
fn stability_results_include_named_control_derivatives() {
    let aircraft = controlled_wing();
    let condition = FlightCondition {
        alpha: 3.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    let derivatives = StabilityDerivativeCalculator
        .compute(&aircraft, &mut solver, &condition)
        .unwrap();
    assert_eq!(derivatives.control_derivatives.len(), 1);
    assert_eq!(derivatives.control_derivatives[0].name, "elevator");
    assert!(derivatives.control_derivatives[0].coefficients.lift.abs() > 1.0e-4);
}

#[test]
fn repeated_control_entries_compose_in_order() {
    let aircraft = controlled_wing();
    let one = FlightCondition {
        control_deflections: vec![("elevator".into(), 0.04)],
        ..FlightCondition::default()
    };
    let repeated = FlightCondition {
        control_deflections: vec![("elevator".into(), 0.01), ("elevator".into(), 0.03)],
        ..FlightCondition::default()
    };
    let solve = |condition: &FlightCondition| {
        let mut solver = VortexLatticeSolver::new();
        solver.set_aircraft(&aircraft).unwrap();
        solver.solve(condition).unwrap().circulation
    };
    assert_abs_diff_eq!(solve(&one), solve(&repeated), epsilon = 1.0e-12);
}

#[test]
fn aircraft_profile_drag_is_added_without_changing_induced_drag() {
    let mut aircraft = rectangular_wing(3, 6);
    aircraft.set_profile_drag_coefficient(0.017);
    let condition = FlightCondition {
        alpha: 4.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    let solution = solver.solve(&condition).unwrap();
    let result = ForceCalculator::new().compute_for_aircraft(
        &aircraft,
        solver.mesh(),
        &solution,
        &condition,
    );
    assert_abs_diff_eq!(result.viscous_drag, 0.017, epsilon = 1.0e-14);
    assert_abs_diff_eq!(
        result.coefficients.drag,
        result.near_field_induced_drag + result.viscous_drag,
        epsilon = 1.0e-12
    );
}

#[test]
fn section_drag_polar_contributes_strip_and_total_viscous_drag() {
    let mut aircraft = rectangular_wing(3, 4);
    let polar = [-1.0, 0.03, 0.0, 0.01, 1.0, 0.03];
    for section in aircraft.surfaces_mut()[0].sections_mut() {
        section.set_viscous_polar(polar);
    }
    let condition = FlightCondition {
        alpha: 5.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    let solution = solver.solve(&condition).unwrap();
    let result = ForceCalculator::new().compute_for_aircraft(
        &aircraft,
        solver.mesh(),
        &solution,
        &condition,
    );
    assert!(result.viscous_drag >= 0.01);
    assert!(
        result.surface_loads[0]
            .strip_loads
            .iter()
            .all(|strip| strip.viscous_drag > 0.0)
    );
}

#[test]
fn force_results_reject_invalid_solution_without_nan_values() {
    let aircraft = rectangular_wing(2, 3);
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    let result = ForceCalculator::new().compute_for_aircraft(
        &aircraft,
        solver.mesh(),
        &pegasus_vlm::solver::SolutionState::default(),
        &FlightCondition::default(),
    );
    assert!(!result.is_valid);
    assert_eq!(result.body_force, Vec3::zeros());
}
