//! Ported from the frozen v0.2.0 `rectangular_wing_test.cpp`, the C++ suite's
//! only integration file.
//!
//! Everything here is compared against closed-form finite-wing theory rather
//! than against a recorded number, so these are the cases that would catch a
//! result that is internally consistent but physically wrong.

use approx::assert_abs_diff_eq;
use pegasus_vlm::aero::{AerodynamicResults, ForceCalculator, StabilityDerivativeCalculator};
use pegasus_vlm::core::{DEGREES_TO_RADIANS, PI, SpacingType, Vec3};
use pegasus_vlm::geometry::{Aircraft, Section, Surface};
use pegasus_vlm::io::AvlReader;
use pegasus_vlm::solver::{FlightCondition, VortexLatticeSolver};

/// The frozen suite compared against theory at 10%.
const THEORY_TOLERANCE: f64 = 0.10;

/// Helmbold's finite-wing lift slope. The v0.2.0 fixture had an erroneous
/// extra pi in its denominator, overpredicting the finite-span lift slope.
fn lift_slope_theory(aspect_ratio: f64) -> f64 {
    2.0 * PI * aspect_ratio / (2.0 + (4.0 + aspect_ratio * aspect_ratio).sqrt())
}

/// A y-duplicated rectangular wing of the given aspect ratio.
fn wing(aspect_ratio: f64, chordwise: usize, spanwise: usize, cosine: bool) -> Aircraft {
    let span = aspect_ratio;
    let chord = 1.0;
    let mut surface = Surface::new("Wing");
    surface.set_chordwise_panel_count(chordwise);
    surface.set_spanwise_panel_count(spanwise);
    if cosine {
        surface.set_chordwise_spacing(SpacingType::Cosine);
        surface.set_spanwise_spacing(SpacingType::Cosine);
    }
    surface.set_y_duplicate(0.0);
    surface.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), chord, 0.0));
    surface.add_section(Section::new(Vec3::new(0.0, span / 2.0, 0.0), chord, 0.0));

    let mut aircraft = Aircraft::new("Rectangular Wing");
    aircraft.add_surface(surface);
    aircraft.set_reference_area(span * chord);
    aircraft.set_reference_chord(chord);
    aircraft.set_reference_span(span);
    aircraft
}

fn standard_wing() -> Aircraft {
    wing(10.0, 8, 20, true)
}

fn solve_at(aircraft: &Aircraft, alpha_radians: f64) -> AerodynamicResults {
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(aircraft).unwrap();
    let condition = FlightCondition {
        alpha: alpha_radians,
        ..FlightCondition::default()
    };
    let solution = solver.solve(&condition).unwrap();
    ForceCalculator::new().compute_for_aircraft(aircraft, solver.mesh(), &solution, &condition)
}

/// Ports `LiftCurveSlopeMatchesTheory`.
#[test]
fn the_lift_curve_slope_matches_theory_at_aspect_ratio_ten() {
    let aircraft = standard_wing();
    let alpha = 5.0 * DEGREES_TO_RADIANS;
    let results = solve_at(&aircraft, alpha);
    let slope = results.coefficients.lift / alpha;
    let theory = lift_slope_theory(10.0);
    assert_abs_diff_eq!(slope, theory, epsilon = theory * THEORY_TOLERANCE);
}

/// Ports `AspectRatioSweep/HighAspectRatioTest.LiftSlopeApproaches2Pi` at
/// aspect ratios 6, 7, 8 and 10.
///
/// The sweep is the real content: agreement at a single aspect ratio could come
/// from a coincidence, agreement across four cannot.
#[test]
fn the_lift_slope_tracks_theory_across_an_aspect_ratio_sweep() {
    for aspect_ratio in [6.0, 7.0, 8.0, 10.0] {
        let spanwise = (aspect_ratio * 2.0) as usize;
        let aircraft = wing(aspect_ratio, 6, spanwise, false);
        let alpha = 5.0 * DEGREES_TO_RADIANS;
        let results = solve_at(&aircraft, alpha);
        let slope = results.coefficients.lift / alpha;
        let theory = lift_slope_theory(aspect_ratio);
        assert_abs_diff_eq!(slope, theory, epsilon = theory * THEORY_TOLERANCE);
    }
}

/// Ports `InducedDragNonZero`.
#[test]
fn induced_drag_is_nonzero_and_of_the_right_order() {
    let aircraft = standard_wing();
    let results = solve_at(&aircraft, 5.0 * DEGREES_TO_RADIANS);
    let ideal = results.coefficients.lift.powi(2) / (PI * 10.0);
    assert!(results.coefficients.drag.is_finite());
    assert!(results.coefficients.drag != 0.0);
    assert!(results.coefficients.drag.abs() < ideal * 10.0);
}

/// Ports `ZeroLiftAtZeroAlpha`.
#[test]
fn a_flat_wing_produces_no_lift_at_zero_incidence() {
    let aircraft = standard_wing();
    let results = solve_at(&aircraft, 0.0);
    assert_abs_diff_eq!(results.coefficients.lift, 0.0, epsilon = 1.0e-10);
}

/// Ports `LiftLinearWithAlpha`.
///
/// A least-squares fit over alpha from -5 to +10 degrees must explain better
/// than 99.9% of the variance. This is a stronger statement than any two-point
/// ratio: it would catch a curvature that a doubling check cannot see.
#[test]
fn lift_is_linear_in_alpha_across_the_sweep() {
    let aircraft = standard_wing();
    let mut samples = Vec::new();
    let mut degrees = -5.0;
    while degrees <= 10.0 + 1.0e-9 {
        let alpha = degrees * DEGREES_TO_RADIANS;
        samples.push((alpha, solve_at(&aircraft, alpha).coefficients.lift));
        degrees += 2.5;
    }

    let n = samples.len() as f64;
    let sum_x: f64 = samples.iter().map(|(x, _)| x).sum();
    let sum_y: f64 = samples.iter().map(|(_, y)| y).sum();
    let sum_xy: f64 = samples.iter().map(|(x, y)| x * y).sum();
    let sum_x2: f64 = samples.iter().map(|(x, _)| x * x).sum();
    let slope = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x);
    let intercept = (sum_y - slope * sum_x) / n;
    let mean_y = sum_y / n;

    let mut total = 0.0;
    let mut residual = 0.0;
    for (x, y) in &samples {
        let predicted = slope * x + intercept;
        total += (y - mean_y).powi(2);
        residual += (y - predicted).powi(2);
    }
    let r_squared = 1.0 - residual / total;
    assert!(r_squared > 0.999, "R^2 was {r_squared}");
}

/// Ports `StabilityDerivativesPhysicallyReasonable`.
#[test]
fn the_stability_derivatives_are_physically_reasonable() {
    let aircraft = standard_wing();
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    let condition = FlightCondition {
        alpha: 5.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    solver.solve(&condition).unwrap();
    let derivs = StabilityDerivativeCalculator
        .compute(&aircraft, &mut solver, &condition)
        .unwrap();

    assert!(derivs.is_valid);
    assert!(derivs.by_alpha.lift > 4.0);
    assert!(derivs.by_alpha.lift < 6.5);
    // Negative roll damping is what makes the wing laterally stable.
    assert!(derivs.by_roll_rate.roll < 0.0);
}

/// Ports `AvlFileTest.ParseRectangularWing`.
#[test]
fn the_repository_avl_fixture_parses_with_its_reference_quantities() {
    let aircraft = AvlReader::default()
        .read_path("tests/fixtures/rect.avl")
        .unwrap();
    assert_eq!(aircraft.surfaces().len(), 1);
    assert_abs_diff_eq!(aircraft.reference_area(), 1.123, epsilon = 0.001);
    assert_abs_diff_eq!(aircraft.reference_chord(), 0.25, epsilon = 0.001);
    assert_abs_diff_eq!(aircraft.reference_span(), 6.01, epsilon = 0.001);
}

/// Ports `AvlFileTest.SolveFromParsedGeometry`.
///
/// Closes the loop: geometry read from a file, not built in code, solves to
/// finite and sensible loads.
#[test]
fn the_repository_avl_fixture_solves_to_sensible_loads() {
    let aircraft = AvlReader::default()
        .read_path("tests/fixtures/rect.avl")
        .unwrap();
    let results = solve_at(&aircraft, 5.0 * DEGREES_TO_RADIANS);
    assert!(results.coefficients.lift > 0.0);
    assert!(results.coefficients.drag.is_finite());
    assert!(results.coefficients.drag != 0.0);
    assert!(results.coefficients.pitch.is_finite());
}
