//! Ported from the frozen v0.2.0 `eigenmode_test.cpp`.
//!
//! Mass properties, flight state, and the rigid-body eigenmode analysis. The
//! frozen derivative set describes a conventionally stable aircraft — negative
//! `Cm_alpha`, `Cm_q`, `Cl_p` and `Cn_r`, positive `Cn_beta` — so every mode it
//! produces must have non-positive damping.
//!
//! Two cases from this file are not ported: `ApparentMassTest.DefaultState` and
//! `EigenmodeAnalysisTest.SetApparentMass`. `ApparentMass` was carried by
//! v0.2.0 but never read, so it has been removed rather than kept as dead API;
//! `tests/parity_cases.rs` records both as unportable.

use approx::assert_abs_diff_eq;
use num_complex::Complex64;
use pegasus_vlm::aero::{Coefficients, StabilityDerivatives};
use pegasus_vlm::core::{DEGREES_TO_RADIANS, Mat3, PI, Vec3};
use pegasus_vlm::dynamics::{
    Eigenmode, EigenmodeAnalysis, EigenmodeResults, FlightState, MassProperties, ModeType,
};
use pegasus_vlm::geometry::Aircraft;

const TOL: f64 = 1.0e-6;

fn typical_mass() -> MassProperties {
    MassProperties {
        mass: 1_000.0,
        center_of_gravity: Vec3::zeros(),
        i_xx: 1_000.0,
        i_yy: 2_000.0,
        i_zz: 3_000.0,
        i_xy: 0.0,
        i_yz: 0.0,
        i_xz: 100.0,
    }
}

fn typical_flight() -> FlightState {
    FlightState {
        velocity: 50.0,
        alpha: 5.0 * DEGREES_TO_RADIANS,
        beta: 0.0,
        density: 1.225,
        gravity: 9.81,
        theta: 5.0 * DEGREES_TO_RADIANS,
        ..FlightState::default()
    }
}

/// The frozen derivative set: a conventionally stable configuration.
fn typical_derivatives() -> StabilityDerivatives {
    StabilityDerivatives {
        is_valid: true,
        reference_area: 16.0,
        reference_chord: 1.5,
        reference_span: 10.0,
        by_alpha: Coefficients {
            lift: 5.0,
            drag: 0.1,
            pitch: -1.0,
            ..Coefficients::default()
        },
        by_beta: Coefficients {
            side: -0.5,
            roll: -0.1,
            yaw: 0.1,
            ..Coefficients::default()
        },
        by_roll_rate: Coefficients {
            side: 0.0,
            roll: -0.5,
            yaw: -0.05,
            ..Coefficients::default()
        },
        by_pitch_rate: Coefficients {
            lift: 5.0,
            pitch: -15.0,
            ..Coefficients::default()
        },
        by_yaw_rate: Coefficients {
            side: 0.2,
            roll: 0.1,
            yaw: -0.1,
            ..Coefficients::default()
        },
        ..StabilityDerivatives::default()
    }
}

fn aircraft_for(derivs: &StabilityDerivatives) -> Aircraft {
    let mut aircraft = Aircraft::new("Test");
    aircraft.set_reference_area(derivs.reference_area);
    aircraft.set_reference_chord(derivs.reference_chord);
    aircraft.set_reference_span(derivs.reference_span);
    aircraft
}

fn mode_with(real: f64, imaginary: f64) -> Eigenmode {
    Eigenmode {
        eigenvalue: Complex64::new(real, imaginary),
        ..Eigenmode::default()
    }
}

// ============================================================================
// MassProperties
// ============================================================================

/// Ports `DefaultIsInvalid` and `ValidWithPositiveValues`.
#[test]
fn mass_properties_validity() {
    assert!(!MassProperties::default().is_valid());
    assert!(typical_mass().is_valid());
}

/// Ports `InertiaTensorSymmetric`.
///
/// Products of inertia enter the tensor negated, and the tensor must come out
/// symmetric — the sign convention is the easy thing to get wrong.
#[test]
fn inertia_tensor_is_symmetric_with_negated_products() {
    let mass = MassProperties {
        i_xx: 100.0,
        i_yy: 200.0,
        i_zz: 300.0,
        i_xy: 10.0,
        i_yz: 20.0,
        i_xz: 30.0,
        ..MassProperties::default()
    };
    let tensor = mass.inertia_tensor();
    assert_abs_diff_eq!(tensor[(0, 1)], tensor[(1, 0)], epsilon = TOL);
    assert_abs_diff_eq!(tensor[(1, 2)], tensor[(2, 1)], epsilon = TOL);
    assert_abs_diff_eq!(tensor[(0, 2)], tensor[(2, 0)], epsilon = TOL);
    assert_abs_diff_eq!(tensor[(0, 0)], 100.0, epsilon = TOL);
    assert_abs_diff_eq!(tensor[(1, 1)], 200.0, epsilon = TOL);
    assert_abs_diff_eq!(tensor[(2, 2)], 300.0, epsilon = TOL);
}

/// Ports `InverseInertiaCorrect`.
#[test]
fn inverse_inertia_tensor_is_a_true_inverse() {
    let mass = MassProperties {
        mass: 1_000.0,
        i_xx: 100.0,
        i_yy: 200.0,
        i_zz: 300.0,
        ..MassProperties::default()
    };
    let product = mass.inertia_tensor() * mass.inverse_inertia_tensor().unwrap();
    assert_abs_diff_eq!(product, Mat3::identity(), epsilon = TOL);
}

/// Ports `SetInertiaTensor`.
///
/// Reading a tensor back into scalars flips the sign of the off-diagonals
/// again, so a tensor with `-5` at `(0,1)` yields `i_xy = +5`.
#[test]
fn setting_the_inertia_tensor_recovers_positive_products() {
    let mut mass = MassProperties::default();
    mass.set_inertia_tensor(&Mat3::new(
        100.0, -5.0, -2.0, -5.0, 200.0, -3.0, -2.0, -3.0, 300.0,
    ));
    assert_abs_diff_eq!(mass.i_xx, 100.0, epsilon = TOL);
    assert_abs_diff_eq!(mass.i_yy, 200.0, epsilon = TOL);
    assert_abs_diff_eq!(mass.i_zz, 300.0, epsilon = TOL);
    assert_abs_diff_eq!(mass.i_xy, 5.0, epsilon = TOL);
    assert_abs_diff_eq!(mass.i_yz, 3.0, epsilon = TOL);
    assert_abs_diff_eq!(mass.i_xz, 2.0, epsilon = TOL);
}

// ============================================================================
// FlightState
// ============================================================================

/// Ports `FlightStateTest.BodyVelocityAtZeroAngles` and `BodyVelocityWithAlpha`.
#[test]
fn body_velocity_resolves_alpha_into_x_and_z() {
    let state = FlightState {
        velocity: 100.0,
        ..FlightState::default()
    };
    assert_abs_diff_eq!(
        state.body_velocity(),
        Vec3::new(100.0, 0.0, 0.0),
        epsilon = TOL
    );

    let state = FlightState {
        velocity: 100.0,
        alpha: 5.0 * DEGREES_TO_RADIANS,
        ..FlightState::default()
    };
    let velocity = state.body_velocity();
    assert_abs_diff_eq!(velocity.x, 100.0 * state.alpha.cos(), epsilon = TOL);
    assert_abs_diff_eq!(velocity.y, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(velocity.z, 100.0 * state.alpha.sin(), epsilon = TOL);
}

/// Ports `FlightStateTest.BodyVelocityWithBeta`.
#[test]
fn positive_sideslip_gives_positive_lateral_velocity() {
    let state = FlightState {
        velocity: 50.0,
        beta: 0.1,
        ..FlightState::default()
    };
    let velocity = state.body_velocity();
    assert!(velocity.y > 0.0);
    assert!(velocity.iter().all(|value| value.is_finite()));
}

/// Ports `FlightStateTest.DynamicPressure`.
#[test]
fn dynamic_pressure_is_half_rho_v_squared() {
    let state = FlightState {
        velocity: 50.0,
        density: 1.225,
        ..FlightState::default()
    };
    assert_abs_diff_eq!(
        state.dynamic_pressure(),
        0.5 * 1.225 * 50.0 * 50.0,
        epsilon = TOL
    );
}

/// Ports `FlightStateTest.DefaultValues`.
///
/// The atmosphere defaults to sea level rather than to zero, so an analysis
/// with only a velocity set is still physically meaningful.
#[test]
fn flight_state_defaults_to_sea_level() {
    let state = FlightState::default();
    assert_abs_diff_eq!(state.velocity, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(state.alpha, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(state.beta, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(state.density, 1.225, epsilon = TOL);
    assert_abs_diff_eq!(state.gravity, 9.81, epsilon = TOL);
}

// ============================================================================
// Eigenmode scalar helpers
// ============================================================================

/// Ports `EigenmodeGettersWork`.
#[test]
fn eigenmode_getters() {
    let mode = mode_with(-0.5, 2.0);
    assert_abs_diff_eq!(mode.damping(), -0.5, epsilon = TOL);
    assert_abs_diff_eq!(mode.frequency(), 2.0, epsilon = TOL);
    assert!(mode.is_stable());
    assert!(mode.is_oscillatory());
    // Period = 2 pi / omega = pi for omega = 2.
    assert_abs_diff_eq!(mode.period(), PI, epsilon = 0.01);
}

/// Ports `NonOscillatoryModeHasZeroPeriod`.
#[test]
fn a_real_eigenvalue_has_no_period() {
    let mode = mode_with(-1.0, 0.0);
    assert!(!mode.is_oscillatory());
    assert_abs_diff_eq!(mode.period(), 0.0, epsilon = TOL);
}

/// Ports `DampingRatioComputation` and `DampingRatioZeroFrequency`.
#[test]
fn damping_ratio_recovers_zeta() {
    // Construct an eigenvalue from a known zeta and natural frequency, then
    // recover zeta from it.
    let zeta = 0.5;
    let omega_n = 2.0;
    let mode = mode_with(-zeta * omega_n, omega_n * (1.0 - zeta * zeta).sqrt());
    assert_abs_diff_eq!(mode.damping_ratio(), zeta, epsilon = 0.01);
    // A zero eigenvalue has no natural frequency to divide by.
    assert_abs_diff_eq!(mode_with(0.0, 0.0).damping_ratio(), 0.0, epsilon = TOL);
}

/// Ports `GetTimeConstant` and `GetTimeConstantZeroDamping`.
///
/// The time constant is the time to half amplitude, `ln 2 / |sigma|`, not
/// `1/|sigma|`.
#[test]
fn time_constant_is_the_time_to_half_amplitude() {
    let mode = mode_with(-0.5, 2.0);
    assert_abs_diff_eq!(mode.time_constant(), 2.0f64.ln() * 2.0, epsilon = 0.01);
    assert_abs_diff_eq!(mode_with(0.0, 2.0).time_constant(), 0.0, epsilon = TOL);
}

/// Ports `UnstableModeDetected`.
#[test]
fn a_positive_real_part_is_unstable() {
    let mode = mode_with(0.5, 1.0);
    assert!(!mode.is_stable());
    assert!(mode.damping() > 0.0);
}

/// Ports `ModeTypeNames`.
///
/// Both `Unknown` and `Other` display as "Unknown": a mode that failed
/// classification and one that classified into no named category read the same
/// to a user.
#[test]
fn mode_type_names() {
    assert_eq!(ModeType::ShortPeriod.name(), "Short Period");
    assert_eq!(ModeType::Phugoid.name(), "Phugoid");
    assert_eq!(ModeType::DutchRoll.name(), "Dutch Roll");
    assert_eq!(ModeType::RollSubsidence.name(), "Roll Subsidence");
    assert_eq!(ModeType::Spiral.name(), "Spiral");
    assert_eq!(ModeType::Unknown.name(), "Unknown");
    assert_eq!(ModeType::Other.name(), "Unknown");
}

// ============================================================================
// EigenmodeResults queries
// ============================================================================

/// Ports `FindModesByTypeEmpty` and `FindModesByTypeWithMatch`.
#[test]
fn modes_can_be_filtered_by_type() {
    let empty = EigenmodeResults::default();
    assert_eq!(empty.modes_of_type(ModeType::ShortPeriod).count(), 0);

    let results = EigenmodeResults {
        modes: vec![
            Eigenmode {
                mode_type: ModeType::ShortPeriod,
                eigenvalue: Complex64::new(-1.0, 3.0),
                ..Eigenmode::default()
            },
            Eigenmode {
                mode_type: ModeType::Phugoid,
                eigenvalue: Complex64::new(-0.1, 0.5),
                ..Eigenmode::default()
            },
        ],
        ..EigenmodeResults::default()
    };
    let short_period: Vec<_> = results.modes_of_type(ModeType::ShortPeriod).collect();
    assert_eq!(short_period.len(), 1);
    assert_abs_diff_eq!(short_period[0].frequency(), 3.0, epsilon = TOL);
}

/// Ports `GetDutchRollNotFound`, `GetPhugoidNotFound` and
/// `GetShortPeriodNotFound`.
#[test]
fn named_mode_lookups_return_none_when_absent() {
    let results = EigenmodeResults::default();
    assert!(results.dutch_roll().is_none());
    assert!(results.phugoid().is_none());
    assert!(results.short_period().is_none());
}

/// Ports `IsStableAllStable`, `IsStableOneUnstable` and `IsStableEmpty`.
///
/// Stability is a conjunction, so an empty result is vacuously stable — the
/// convention a caller has to know.
#[test]
fn overall_stability_is_a_conjunction_over_modes() {
    let stable = EigenmodeResults {
        modes: vec![mode_with(-1.0, 2.0), mode_with(-0.5, 0.0)],
        ..EigenmodeResults::default()
    };
    assert!(stable.is_stable());

    let mixed = EigenmodeResults {
        modes: vec![mode_with(-1.0, 2.0), mode_with(0.1, 0.5)],
        ..EigenmodeResults::default()
    };
    assert!(!mixed.is_stable());

    assert!(EigenmodeResults::default().is_stable());
}

// ============================================================================
// Full analyses
// ============================================================================

/// Ports `LongitudinalAnalysisProducesModes` and `SystemMatrixHasCorrectSize`.
#[test]
fn longitudinal_analysis_produces_a_four_state_system() {
    let derivs = typical_derivatives();
    let results = EigenmodeAnalysis::new()
        .analyze_longitudinal(
            &aircraft_for(&derivs),
            &typical_mass(),
            &derivs,
            &typical_flight(),
        )
        .unwrap();
    assert_eq!(results.state_count, 4);
    assert!(!results.modes.is_empty());
    assert_eq!(results.system_matrix.shape(), (4, 4));
}

/// Ports `LateralDirectionalAnalysisProducesModes`.
#[test]
fn lateral_analysis_produces_a_four_state_system() {
    let derivs = typical_derivatives();
    let results = EigenmodeAnalysis::new()
        .analyze_lateral(
            &aircraft_for(&derivs),
            &typical_mass(),
            &derivs,
            &typical_flight(),
        )
        .unwrap();
    assert_eq!(results.state_count, 4);
    assert!(!results.modes.is_empty());
}

/// Ports `StableAircraftHasNegativeDamping`.
///
/// The frozen derivative set is conventionally stable, so every longitudinal
/// eigenvalue must sit in the closed left half plane.
#[test]
fn a_stable_configuration_produces_no_diverging_mode() {
    let derivs = typical_derivatives();
    let results = EigenmodeAnalysis::new()
        .analyze_longitudinal(
            &aircraft_for(&derivs),
            &typical_mass(),
            &derivs,
            &typical_flight(),
        )
        .unwrap();
    for mode in &results.modes {
        assert!(
            mode.damping() <= 0.0,
            "mode with eigenvalue {} is unstable",
            mode.eigenvalue
        );
    }
}

/// Ports `ShortPeriodHigherFrequencyThanPhugoid`.
///
/// The defining separation of the two longitudinal modes. Guarded because
/// classification may not name both on every configuration, which is how the
/// frozen case was written.
#[test]
fn the_short_period_is_faster_than_the_phugoid() {
    let derivs = typical_derivatives();
    let results = EigenmodeAnalysis::new()
        .analyze_longitudinal(
            &aircraft_for(&derivs),
            &typical_mass(),
            &derivs,
            &typical_flight(),
        )
        .unwrap();
    if let (Some(short_period), Some(phugoid)) = (results.short_period(), results.phugoid()) {
        assert!(short_period.frequency() > phugoid.frequency());
    }
}

/// Ports `InvalidMassPropertiesFails` and `InvalidDerivativesFails`.
///
/// The C++ returned an invalid result; the Rust returns `Err`, which cannot be
/// read past by accident.
#[test]
fn invalid_inputs_are_rejected() {
    let derivs = typical_derivatives();
    let aircraft = aircraft_for(&derivs);
    assert!(
        EigenmodeAnalysis::new()
            .analyze_longitudinal(
                &aircraft,
                &MassProperties::default(),
                &derivs,
                &typical_flight()
            )
            .is_err()
    );
    assert!(
        EigenmodeAnalysis::new()
            .analyze_longitudinal(
                &aircraft,
                &typical_mass(),
                &StabilityDerivatives::default(),
                &typical_flight()
            )
            .is_err()
    );
}

/// Ports `FullAnalysisProducesMoreModes`.
#[test]
fn the_full_analysis_carries_more_states_than_the_longitudinal_one() {
    let derivs = typical_derivatives();
    let aircraft = aircraft_for(&derivs);
    let full = EigenmodeAnalysis::new()
        .analyze_full(&aircraft, &typical_mass(), &derivs, &typical_flight())
        .unwrap();
    let longitudinal = EigenmodeAnalysis::new()
        .analyze_longitudinal(&aircraft, &typical_mass(), &derivs, &typical_flight())
        .unwrap();
    assert!(full.state_count > longitudinal.state_count);
}
