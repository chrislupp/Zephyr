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

mod common;

use approx::assert_abs_diff_eq;
use serde_json::Value;
use zephyr_vlm::aero::{
    Coefficients, ForceCalculator, StabilityDerivativeCalculator, StabilityDerivatives,
};
use zephyr_vlm::core::{DEGREES_TO_RADIANS, MatX, Vec3};
use zephyr_vlm::dynamics::{EigenmodeAnalysis, FlightState, MassProperties};
use zephyr_vlm::gradients::GradientCalculator;
use zephyr_vlm::solver::{FlightCondition, VortexLatticeSolver};
use zephyr_vlm::trim::{TrimPresets, TrimSolver};

fn array(value: &Value) -> Vec<f64> {
    value
        .as_array()
        .expect("fixture array")
        .iter()
        .map(|entry| entry.as_f64().expect("fixture scalar"))
        .collect()
}

fn assert_slice_close(actual: &[f64], expected: &[f64], epsilon: f64) {
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.iter().zip(expected) {
        assert_abs_diff_eq!(actual, expected, epsilon = epsilon);
    }
}

#[test]
fn rectangular_wing_matches_v020_oracle() {
    let oracle: Value =
        serde_json::from_str(include_str!("fixtures/parity/rectangular_wing_v0.2.0.json"))
            .expect("valid committed oracle fixture");
    assert_eq!(oracle["provenance"]["tag"], "v0.2.0");
    assert_eq!(oracle["provenance"]["commit"], "556b0ea");

    let tolerances = &oracle["tolerances"];
    let geometry_tolerance = tolerances["geometry"].as_f64().unwrap();
    let matrix_tolerance = tolerances["matrix"].as_f64().unwrap();
    let solution_tolerance = tolerances["solution"].as_f64().unwrap();
    let aero_tolerance = tolerances["aero"].as_f64().unwrap();
    let stability_tolerance = tolerances["stability"].as_f64().unwrap();

    let aircraft = common::rectangular_wing(2, 4);
    let condition = FlightCondition {
        alpha: 4.0 * DEGREES_TO_RADIANS,
        beta: 2.0 * DEGREES_TO_RADIANS,
        mach: 0.2,
        ..FlightCondition::default()
    };
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    let solution = solver.solve(&condition).unwrap();
    let mesh = solver.mesh();
    assert_eq!(
        mesh.vortex_count(),
        oracle["dimensions"]["vortices"].as_u64().unwrap() as usize
    );
    assert_eq!(
        mesh.strip_count(),
        oracle["dimensions"]["strips"].as_u64().unwrap() as usize
    );

    let first = &mesh.vortices()[0];
    assert_slice_close(
        first.bound_vortex_left.as_slice(),
        &array(&oracle["geometry"]["first_bound_left"]),
        geometry_tolerance,
    );
    assert_slice_close(
        first.bound_vortex_right.as_slice(),
        &array(&oracle["geometry"]["first_bound_right"]),
        geometry_tolerance,
    );
    assert_slice_close(
        first.control_point.as_slice(),
        &array(&oracle["geometry"]["first_control_point"]),
        geometry_tolerance,
    );
    assert_slice_close(
        first.normal_at_control_point.as_slice(),
        &array(&oracle["geometry"]["first_normal"]),
        geometry_tolerance,
    );

    for (row, expected_row) in oracle["aic"].as_array().unwrap().iter().enumerate() {
        let actual_row: Vec<_> = solver.aic().row(row).iter().copied().collect();
        assert_slice_close(&actual_row, &array(expected_row), matrix_tolerance);
    }
    for (row, expected_row) in oracle["velocity"].as_array().unwrap().iter().enumerate() {
        let actual_row: Vec<_> = solver.velocity_matrix().row(row).iter().copied().collect();
        assert_slice_close(&actual_row, &array(expected_row), matrix_tolerance);
    }
    assert_slice_close(
        solution.circulation.as_slice(),
        &array(&oracle["circulation"]),
        solution_tolerance,
    );
    for (name, actual) in [
        ("u", &solution.circulation_per_u),
        ("v", &solution.circulation_per_v),
        ("w", &solution.circulation_per_w),
        ("p", &solution.circulation_per_p),
        ("q", &solution.circulation_per_q),
        ("r", &solution.circulation_per_r),
    ] {
        assert_slice_close(
            actual.as_slice(),
            &array(&oracle["unit_responses"][name]),
            solution_tolerance,
        );
    }

    let aero = ForceCalculator::new().compute_for_aircraft(&aircraft, mesh, &solution, &condition);
    for (actual, name) in [
        (aero.coefficients.lift, "CL"),
        (aero.coefficients.drag, "CD"),
        (aero.coefficients.side, "CY"),
        (aero.coefficients.roll, "Cl"),
        (aero.coefficients.pitch, "Cm"),
        (aero.coefficients.yaw, "Cn"),
        (aero.far_field_lift, "CL_ff"),
        (aero.far_field_side, "CY_ff"),
        (aero.far_field_induced_drag, "CDi_ff"),
        (aero.span_efficiency, "span_efficiency"),
    ] {
        assert_abs_diff_eq!(
            actual,
            oracle["aero"][name].as_f64().unwrap(),
            epsilon = aero_tolerance
        );
    }

    let stability = StabilityDerivativeCalculator
        .compute(&aircraft, &mut solver, &condition)
        .unwrap();
    for (actual, name) in [
        (stability.by_alpha.lift, "CL_alpha"),
        (stability.by_alpha.drag, "CD_alpha"),
        (stability.by_alpha.pitch, "Cm_alpha"),
        (stability.by_beta.side, "CY_beta"),
        (stability.by_roll_rate.roll, "Cl_p"),
        (stability.by_pitch_rate.pitch, "Cm_q"),
        (stability.by_yaw_rate.yaw, "Cn_r"),
    ] {
        assert_abs_diff_eq!(
            actual,
            oracle["stability"][name].as_f64().unwrap(),
            epsilon = stability_tolerance
        );
    }
}

fn assert_matrix_close(actual: &MatX, expected: &Value, epsilon: f64) {
    let rows = expected.as_array().unwrap();
    assert_eq!(actual.nrows(), rows.len());
    assert_eq!(actual.ncols(), rows[0].as_array().unwrap().len());
    for (row, expected_row) in rows.iter().enumerate() {
        for (column, expected) in expected_row.as_array().unwrap().iter().enumerate() {
            assert_abs_diff_eq!(
                actual[(row, column)],
                expected.as_f64().unwrap(),
                epsilon = epsilon
            );
        }
    }
}

fn mode_fixture() -> (MassProperties, StabilityDerivatives, FlightState) {
    (
        MassProperties {
            mass: 1_000.0,
            center_of_gravity: Vec3::new(0.25, 0.0, 0.0),
            i_xx: 900.0,
            i_yy: 1_500.0,
            i_zz: 2_000.0,
            i_xy: 10.0,
            i_yz: 20.0,
            i_xz: 30.0,
        },
        StabilityDerivatives {
            is_valid: true,
            reference_area: 4.0,
            reference_chord: 1.0,
            reference_span: 4.0,
            velocity: 50.0,
            by_alpha: Coefficients {
                lift: 4.8,
                drag: 0.15,
                pitch: -0.8,
                ..Coefficients::default()
            },
            by_beta: Coefficients {
                side: -0.6,
                roll: -0.12,
                yaw: 0.18,
                ..Coefficients::default()
            },
            by_roll_rate: Coefficients {
                side: -0.1,
                roll: -0.55,
                yaw: -0.04,
                ..Coefficients::default()
            },
            by_pitch_rate: Coefficients {
                lift: 3.0,
                pitch: -8.0,
                ..Coefficients::default()
            },
            by_yaw_rate: Coefficients {
                side: 0.15,
                roll: 0.08,
                yaw: -0.25,
                ..Coefficients::default()
            },
            ..StabilityDerivatives::default()
        },
        FlightState {
            velocity: 50.0,
            density: 1.225,
            gravity: 9.81,
            ..FlightState::default()
        },
    )
}

fn assert_eigenvalues_match(actual: &[(f64, f64)], expected: &Value, epsilon: f64) {
    let mut actual = actual.to_vec();
    actual.sort_by(|left, right| {
        left.0
            .total_cmp(&right.0)
            .then_with(|| left.1.total_cmp(&right.1))
    });
    let mut expected: Vec<_> = expected
        .as_array()
        .unwrap()
        .iter()
        .map(|pair| (pair[0].as_f64().unwrap(), pair[1].as_f64().unwrap().abs()))
        .collect();
    expected.sort_by(|left, right| {
        left.0
            .total_cmp(&right.0)
            .then_with(|| left.1.total_cmp(&right.1))
    });
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.iter().zip(expected) {
        assert_abs_diff_eq!(actual.0, expected.0, epsilon = epsilon);
        assert_abs_diff_eq!(actual.1, expected.1, epsilon = epsilon);
    }
}

#[test]
fn trim_gradients_and_modes_match_v020_oracle_with_documented_transitions() {
    let oracle: Value = serde_json::from_str(include_str!(
        "fixtures/parity/trim_eigen_gradients_v0.2.0.json"
    ))
    .unwrap();
    assert_eq!(oracle["provenance"]["commit"], "556b0ea");
    assert_eq!(oracle["known_transitions"]["rust_full_state_count"], 12);
    let aircraft = common::rectangular_wing(2, 4);

    let mut trim_solver = VortexLatticeSolver::new();
    let trim = TrimSolver::new().solve(
        &aircraft,
        &mut trim_solver,
        &TrimPresets::level_flight(0.35),
        None,
    );
    let trim_tolerance = oracle["tolerances"]["trim"].as_f64().unwrap();
    assert!(trim.converged);
    assert_eq!(
        trim.iterations,
        oracle["trim"]["iterations"].as_u64().unwrap() as usize
    );
    assert_abs_diff_eq!(
        trim.state.alpha,
        oracle["trim"]["alpha"].as_f64().unwrap(),
        epsilon = trim_tolerance
    );
    assert_abs_diff_eq!(
        trim.aerodynamic_results.coefficients.lift,
        oracle["trim"]["CL"].as_f64().unwrap(),
        epsilon = trim_tolerance
    );

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
    let gradient_tolerance = oracle["tolerances"]["gradient"].as_f64().unwrap();
    for (actual, name) in [
        (gradients.by_alpha.coefficients.lift, "CL_alpha"),
        (gradients.by_alpha.coefficients.drag, "CD_alpha"),
        (gradients.by_alpha.coefficients.pitch, "Cm_alpha"),
        (gradients.by_beta.coefficients.side, "CY_beta"),
        (gradients.by_roll_rate.coefficients.roll, "Cl_p"),
        (gradients.by_pitch_rate.coefficients.pitch, "Cm_q"),
        (gradients.by_yaw_rate.coefficients.yaw, "Cn_r"),
        (gradients.by_mach.coefficients.lift, "CL_mach"),
    ] {
        assert_abs_diff_eq!(
            actual,
            oracle["gradients"][name].as_f64().unwrap(),
            epsilon = gradient_tolerance
        );
    }

    let (mass, derivatives, flight) = mode_fixture();
    let analysis = EigenmodeAnalysis::new();
    let longitudinal = analysis
        .analyze_longitudinal(&aircraft, &mass, &derivatives, &flight)
        .unwrap();
    let lateral = analysis
        .analyze_lateral(&aircraft, &mass, &derivatives, &flight)
        .unwrap();
    let full = analysis
        .analyze_full(&aircraft, &mass, &derivatives, &flight)
        .unwrap();
    let matrix_tolerance = oracle["tolerances"]["system_matrix"].as_f64().unwrap();
    assert_matrix_close(
        &longitudinal.system_matrix,
        &oracle["dynamics"]["longitudinal_matrix"],
        matrix_tolerance,
    );
    assert_matrix_close(
        &lateral.system_matrix,
        &oracle["dynamics"]["lateral_matrix"],
        matrix_tolerance,
    );
    assert_eq!(full.system_matrix.shape(), (12, 12));
    let eigenvalue_tolerance = oracle["tolerances"]["eigenvalue"].as_f64().unwrap();
    assert_eigenvalues_match(
        &longitudinal
            .modes
            .iter()
            .map(|mode| (mode.eigenvalue.re, mode.eigenvalue.im.abs()))
            .collect::<Vec<_>>(),
        &oracle["dynamics"]["longitudinal_eigenvalues"],
        eigenvalue_tolerance,
    );
    assert_eigenvalues_match(
        &lateral
            .modes
            .iter()
            .map(|mode| (mode.eigenvalue.re, mode.eigenvalue.im.abs()))
            .collect::<Vec<_>>(),
        &oracle["dynamics"]["lateral_eigenvalues"],
        eigenvalue_tolerance,
    );
}
