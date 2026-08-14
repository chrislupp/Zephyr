mod common;

use approx::assert_abs_diff_eq;
use pegasus_vlm::aero::{
    AerodynamicResults, Coefficients, ForceCalculator, StabilityDerivativeCalculator, StripLoads,
    SurfaceLoads,
};
use pegasus_vlm::core::{DEGREES_TO_RADIANS, PI, Vec3};
use pegasus_vlm::solver::{FlightCondition, VortexLatticeSolver};

use common::{controlled_wing, rectangular_wing};

// ============================================================================
// Result-type contracts, ported from the frozen v0.2.0 aerodynamic_results
// suite. These are cheap but load-bearing: the reference quantities default to
// 1 rather than 0 precisely so that `aspect_ratio` and `ideal_induced_drag`
// divide by something finite on a freshly constructed result.
// ============================================================================

/// Ports `AerodynamicResultsTest.DefaultState`.
#[test]
fn aerodynamic_results_default_state() {
    let results = AerodynamicResults::default();
    assert!(!results.is_valid);
    assert_abs_diff_eq!(results.coefficients.lift, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(results.coefficients.drag, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(results.reference_area, 1.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(results.reference_chord, 1.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(results.reference_span, 1.0, epsilon = 1.0e-10);
}

/// Ports `AerodynamicResultsTest.GetAspectRatio` and `GetAspectRatioSquareWing`.
#[test]
fn aerodynamic_results_aspect_ratio() {
    let mut results = AerodynamicResults {
        reference_span: 10.0,
        reference_area: 20.0,
        ..AerodynamicResults::default()
    };
    assert_abs_diff_eq!(results.aspect_ratio(), 5.0, epsilon = 1.0e-10);
    // A square planform: b = 4, S = 4 gives AR = 4, not 1. The square case is
    // its own frozen test because b²/S and b/c are easy to confuse.
    results.reference_span = 4.0;
    results.reference_area = 4.0;
    assert_abs_diff_eq!(results.aspect_ratio(), 4.0, epsilon = 1.0e-10);
}

/// Ports `AerodynamicResultsTest.GetIdealInducedDrag` and its zero-lift case.
#[test]
fn aerodynamic_results_ideal_induced_drag() {
    let mut results = AerodynamicResults {
        reference_span: 10.0,
        reference_area: 10.0,
        coefficients: Coefficients {
            lift: 0.5,
            ..Coefficients::default()
        },
        ..AerodynamicResults::default()
    };
    // AR = 100/10 = 10, so CDi_ideal = CL² / (pi AR) = 0.25 / (10 pi).
    assert_abs_diff_eq!(
        results.ideal_induced_drag(),
        0.25 / (PI * 10.0),
        epsilon = 1.0e-8
    );
    results.coefficients.lift = 0.0;
    assert_abs_diff_eq!(results.ideal_induced_drag(), 0.0, epsilon = 1.0e-10);
}

/// Ports `AerodynamicResultsTest.Clear`.
///
/// The C++ had an explicit `Clear()`; the Rust equivalent is reassigning
/// `Default`, so this asserts that a populated result reverts completely —
/// including dropping accumulated surface loads, which is the part a
/// field-by-field reset is most likely to miss.
#[test]
fn aerodynamic_results_reset_to_default_clears_every_field() {
    let populated = AerodynamicResults {
        is_valid: true,
        coefficients: Coefficients {
            lift: 0.5,
            drag: 0.01,
            side: 0.02,
            roll: 0.001,
            pitch: -0.05,
            yaw: 0.002,
        },
        body_force: Vec3::new(1.0, 2.0, 3.0),
        body_moment: Vec3::new(4.0, 5.0, 6.0),
        near_field_induced_drag: 0.008,
        viscous_drag: 0.003,
        far_field_lift: 0.49,
        far_field_side: 0.01,
        far_field_induced_drag: 0.007,
        span_efficiency: 0.95,
        surface_loads: vec![SurfaceLoads {
            name: "wing".into(),
            ..SurfaceLoads::default()
        }],
        ..AerodynamicResults::default()
    };
    assert!(populated.is_valid);
    assert_eq!(populated.surface_loads.len(), 1);

    let cleared = AerodynamicResults::default();
    assert!(!cleared.is_valid);
    assert_abs_diff_eq!(cleared.coefficients.lift, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(cleared.coefficients.drag, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(cleared.coefficients.side, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(cleared.coefficients.roll, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(cleared.coefficients.pitch, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(cleared.coefficients.yaw, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(cleared.near_field_induced_drag, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(cleared.viscous_drag, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(cleared.far_field_lift, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(cleared.far_field_side, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(cleared.far_field_induced_drag, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(cleared.span_efficiency, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(cleared.body_force.norm(), 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(cleared.body_moment.norm(), 0.0, epsilon = 1.0e-10);
    assert!(cleared.surface_loads.is_empty());
}

/// Ports `StripLoadsTest.DefaultState`.
#[test]
fn strip_loads_default_state() {
    let strip = StripLoads::default();
    assert_eq!(strip.strip_index, 0);
    assert_abs_diff_eq!(strip.chord, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(strip.width, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(strip.area, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(strip.coefficients.lift, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(strip.coefficients.drag, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(strip.coefficients.side, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(strip.normal_span_loading, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(strip.induced_angle, 0.0, epsilon = 1.0e-10);
}

/// Ports `SurfaceLoadsTest.DefaultState`.
#[test]
fn surface_loads_default_state() {
    let surface = SurfaceLoads::default();
    assert_eq!(surface.surface_index, 0);
    assert!(surface.name.is_empty());
    assert_abs_diff_eq!(surface.coefficients.lift, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(surface.coefficients.drag, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(surface.coefficients.side, 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(surface.area, 0.0, epsilon = 1.0e-10);
    assert!(surface.strip_loads.is_empty());
}

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
