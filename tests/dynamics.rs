mod common;

use approx::assert_abs_diff_eq;
use num_complex::Complex64;
use pegasus_vlm::aero::{Coefficients, StabilityDerivatives};
use pegasus_vlm::core::{Mat3, Vec3};
use pegasus_vlm::dynamics::{
    ApparentMass, Eigenmode, EigenmodeAnalysis, FlightState, MassProperties, ModeType,
};

use common::rectangular_wing;

fn mass_properties() -> MassProperties {
    MassProperties {
        mass: 1_000.0,
        center_of_gravity: Vec3::new(0.25, 0.0, 0.0),
        i_xx: 900.0,
        i_yy: 1_500.0,
        i_zz: 2_000.0,
        i_xy: 10.0,
        i_yz: 20.0,
        i_xz: 30.0,
    }
}

fn derivatives() -> StabilityDerivatives {
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
    }
}

fn flight() -> FlightState {
    FlightState {
        velocity: 50.0,
        density: 1.225,
        gravity: 9.81,
        ..FlightState::default()
    }
}

#[test]
fn mass_property_validity_and_inertia_tensor_signs() {
    assert!(!MassProperties::default().is_valid());
    let mass = mass_properties();
    assert!(mass.is_valid());
    let tensor = mass.inertia_tensor();
    assert_eq!(tensor[(0, 0)], 900.0);
    assert_eq!(tensor[(0, 1)], -10.0);
    assert_eq!(tensor[(0, 2)], -30.0);
    assert_eq!(tensor[(1, 2)], -20.0);
    assert_abs_diff_eq!(tensor, tensor.transpose());
}

#[test]
fn inertia_tensor_setter_and_inverse_are_consistent() {
    let tensor = Mat3::new(2.0, -0.1, -0.2, -0.1, 3.0, -0.3, -0.2, -0.3, 4.0);
    let mut mass = MassProperties {
        mass: 1.0,
        ..MassProperties::default()
    };
    mass.set_inertia_tensor(&tensor);
    assert_abs_diff_eq!(mass.inertia_tensor(), tensor);
    let inverse = mass.inverse_inertia_tensor().unwrap();
    assert_abs_diff_eq!(tensor * inverse, Mat3::identity(), epsilon = 1.0e-12);
}

#[test]
fn flight_state_body_velocity_and_dynamic_pressure() {
    let state = FlightState {
        velocity: 100.0,
        alpha: 0.1,
        beta: -0.05,
        density: 1.0,
        ..FlightState::default()
    };
    assert_abs_diff_eq!(state.body_velocity().norm(), 100.0, epsilon = 1.0e-12);
    assert!(state.body_velocity().y < 0.0);
    assert!(state.body_velocity().z > 0.0);
    assert_eq!(state.dynamic_pressure(), 5_000.0);
}

#[test]
fn eigenmode_helpers_handle_oscillatory_and_real_modes() {
    let oscillatory = Eigenmode {
        eigenvalue: Complex64::new(-0.5, 2.0),
        mode_type: ModeType::DutchRoll,
        ..Eigenmode::default()
    };
    assert!(oscillatory.is_stable());
    assert!(oscillatory.is_oscillatory());
    assert_abs_diff_eq!(oscillatory.frequency(), 2.0);
    assert_abs_diff_eq!(oscillatory.period(), std::f64::consts::PI);
    assert!(oscillatory.damping_ratio() > 0.0);
    assert!(oscillatory.time_constant() > 0.0);
    let real = Eigenmode {
        eigenvalue: Complex64::new(0.2, 0.0),
        ..Eigenmode::default()
    };
    assert!(!real.is_stable());
    assert_eq!(real.period(), 0.0);
}

#[test]
fn mode_type_names_are_deterministic() {
    assert_eq!(ModeType::ShortPeriod.name(), "Short Period");
    assert_eq!(ModeType::Phugoid.name(), "Phugoid");
    assert_eq!(ModeType::DutchRoll.name(), "Dutch Roll");
    assert_eq!(ModeType::RollSubsidence.name(), "Roll Subsidence");
    assert_eq!(ModeType::Spiral.name(), "Spiral");
}

#[test]
fn longitudinal_matrix_and_eigenpairs_have_expected_dimensions() {
    let aircraft = rectangular_wing(2, 4);
    let analysis = EigenmodeAnalysis::new();
    let result = analysis
        .analyze_longitudinal(&aircraft, &mass_properties(), &derivatives(), &flight())
        .unwrap();
    assert_eq!(result.system_matrix.shape(), (4, 4));
    assert_eq!(result.state_count, 4);
    assert!(!result.modes.is_empty());
    assert!(result.modes.iter().all(|mode| {
        (mode
            .eigenvector
            .iter()
            .map(|value| value.norm_sqr())
            .sum::<f64>()
            - 1.0)
            .abs()
            < 1.0e-10
    }));
}

#[test]
fn lateral_matrix_and_eigenpairs_have_expected_dimensions() {
    let aircraft = rectangular_wing(2, 4);
    let analysis = EigenmodeAnalysis::new();
    let result = analysis
        .analyze_lateral(&aircraft, &mass_properties(), &derivatives(), &flight())
        .unwrap();
    assert_eq!(result.system_matrix.shape(), (4, 4));
    assert_eq!(result.state_count, 4);
    assert!(!result.modes.is_empty());
    assert!(result.modes.iter().any(|mode| {
        matches!(
            mode.mode_type,
            ModeType::DutchRoll | ModeType::RollSubsidence | ModeType::Spiral
        )
    }));
}

#[test]
fn full_analysis_builds_a_true_twelve_state_matrix() {
    let aircraft = rectangular_wing(2, 4);
    let analysis = EigenmodeAnalysis::new();
    let result = analysis
        .analyze_full(&aircraft, &mass_properties(), &derivatives(), &flight())
        .unwrap();
    assert_eq!(result.system_matrix.shape(), (12, 12));
    assert_eq!(result.state_count, 12);
    assert_eq!(result.reference_velocity, 50.0);
}

#[test]
fn eigensolver_rejects_invalid_mass_derivatives_and_flight_state() {
    let aircraft = rectangular_wing(2, 4);
    let analysis = EigenmodeAnalysis::new();
    assert!(
        analysis
            .analyze_longitudinal(
                &aircraft,
                &MassProperties::default(),
                &derivatives(),
                &flight()
            )
            .is_err()
    );
    assert!(
        analysis
            .analyze_longitudinal(
                &aircraft,
                &mass_properties(),
                &StabilityDerivatives::default(),
                &flight()
            )
            .is_err()
    );
    assert!(
        analysis
            .analyze_longitudinal(
                &aircraft,
                &mass_properties(),
                &derivatives(),
                &FlightState::default()
            )
            .is_err()
    );
}

#[test]
fn apparent_mass_changes_system_matrices() {
    let aircraft = rectangular_wing(2, 4);
    let base = EigenmodeAnalysis::new().longitudinal_matrix(
        &aircraft,
        &mass_properties(),
        &derivatives(),
        &flight(),
    );
    let mut analysis = EigenmodeAnalysis::new();
    analysis.set_apparent_mass(ApparentMass {
        mass_matrix: Mat3::from_diagonal_element(100.0),
        inertia_matrix: Mat3::from_diagonal_element(50.0),
        include_effects: true,
    });
    let adjusted =
        analysis.longitudinal_matrix(&aircraft, &mass_properties(), &derivatives(), &flight());
    assert!((&base - adjusted).norm() > 1.0e-8);
}

#[test]
fn mode_filters_and_stability_helpers_are_safe_when_empty() {
    let result = pegasus_vlm::dynamics::EigenmodeResults::default();
    assert!(result.short_period().is_none());
    assert!(result.phugoid().is_none());
    assert!(result.dutch_roll().is_none());
    assert!(result.is_stable());
}
