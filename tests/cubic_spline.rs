//! Ported from the frozen v0.2.0 `cubic_spline_test.cpp`.
//!
//! The spline backs airfoil camber lines, so its interpolation and both
//! derivatives are load-bearing rather than incidental. Where the C++ threw,
//! the Rust returns a typed `Result`; those cases assert the error instead of
//! the panic.

use approx::assert_abs_diff_eq;
use pegasus_vlm::core::PI;
use pegasus_vlm::utils::{CubicSpline, SplineBoundaryCondition};

const TOL: f64 = 1.0e-10;

/// Ports `LinearDataInterpolation`.
///
/// A cubic through collinear data must reproduce the line exactly: constant
/// first derivative and an identically zero second derivative. Any spurious
/// curvature from the tridiagonal solve shows up here first.
#[test]
fn linear_data_interpolation() {
    let s = [0.0, 1.0, 2.0, 3.0, 4.0];
    let x = [0.0, 2.0, 4.0, 6.0, 8.0];
    let spline = CubicSpline::build(&s, &x).unwrap();
    for (at, expected) in [(0.5, 1.0), (1.5, 3.0), (2.5, 5.0), (3.5, 7.0)] {
        assert_abs_diff_eq!(spline.evaluate(at).unwrap(), expected, epsilon = TOL);
    }
    for at in [0.5, 2.0, 3.5] {
        assert_abs_diff_eq!(spline.evaluate_derivative(at).unwrap(), 2.0, epsilon = TOL);
    }
    for at in [1.0, 2.5] {
        assert_abs_diff_eq!(
            spline.evaluate_second_derivative(at).unwrap(),
            0.0,
            epsilon = TOL
        );
    }
}

/// Ports `QuadraticDataWithMatchingBoundaryConditions`.
///
/// A cubic spline reproduces a quadratic exactly only when the end slopes are
/// supplied; with natural ends the curvature is forced to zero there and the
/// fit degrades. Here `x = s²` on `[0, 4]` is pinned with `x'(0) = 0` and
/// `x'(4) = 8`.
#[test]
fn quadratic_data_with_matching_boundary_conditions() {
    let s = [0.0, 1.0, 2.0, 3.0, 4.0];
    let x: Vec<f64> = s.iter().map(|value| value * value).collect();
    let spline = CubicSpline::build_with_boundaries(
        &s,
        &x,
        SplineBoundaryCondition::SpecifiedSlope,
        SplineBoundaryCondition::SpecifiedSlope,
        0.0,
        8.0,
    )
    .unwrap();
    for (at, expected) in [(0.5, 0.25), (1.5, 2.25), (2.5, 6.25)] {
        assert_abs_diff_eq!(spline.evaluate(at).unwrap(), expected, epsilon = TOL);
    }
    for (at, expected) in [(0.5, 1.0), (1.5, 3.0), (2.5, 5.0)] {
        assert_abs_diff_eq!(
            spline.evaluate_derivative(at).unwrap(),
            expected,
            epsilon = TOL
        );
    }
    for at in [1.0, 2.5] {
        assert_abs_diff_eq!(
            spline.evaluate_second_derivative(at).unwrap(),
            2.0,
            epsilon = TOL
        );
    }
}

/// Ports `SinusoidalData`.
#[test]
fn sinusoidal_data_interpolates_to_fourth_order() {
    let count = 21;
    let s: Vec<f64> = (0..count)
        .map(|i| 2.0 * PI * i as f64 / (count - 1) as f64)
        .collect();
    let x: Vec<f64> = s.iter().map(|value| value.sin()).collect();
    let spline = CubicSpline::build(&s, &x).unwrap();
    for window in s.windows(2) {
        let midpoint = 0.5 * (window[0] + window[1]);
        assert_abs_diff_eq!(
            spline.evaluate(midpoint).unwrap(),
            midpoint.sin(),
            epsilon = 1.0e-4
        );
    }
    let at = PI / 4.0;
    assert_abs_diff_eq!(
        spline.evaluate_derivative(at).unwrap(),
        at.cos(),
        epsilon = 1.0e-3
    );
}

/// Ports `EvaluateAllConsistency`.
#[test]
fn evaluate_all_agrees_with_the_individual_accessors() {
    let s = [0.0, 1.0, 2.0, 3.0];
    let x = [1.0, 2.0, 0.5, 3.0];
    let spline = CubicSpline::build(&s, &x).unwrap();
    let at = 1.7;
    let (value, first, second) = spline.evaluate_all(at).unwrap();
    assert_abs_diff_eq!(value, spline.evaluate(at).unwrap(), epsilon = TOL);
    assert_abs_diff_eq!(
        first,
        spline.evaluate_derivative(at).unwrap(),
        epsilon = TOL
    );
    assert_abs_diff_eq!(
        second,
        spline.evaluate_second_derivative(at).unwrap(),
        epsilon = TOL
    );
}

/// Ports `EndpointInterpolation`.
#[test]
fn the_spline_passes_through_every_knot() {
    let s = [0.0, 1.0, 2.0, 3.0, 4.0];
    let x = [1.5, 2.3, 0.7, 4.1, 2.9];
    let spline = CubicSpline::build(&s, &x).unwrap();
    for (parameter, value) in s.iter().zip(&x) {
        assert_abs_diff_eq!(spline.evaluate(*parameter).unwrap(), value, epsilon = TOL);
    }
}

/// Ports `SpecifiedSlopeBoundaryCondition`.
#[test]
fn specified_slopes_are_honored_at_both_ends() {
    let s = [0.0, 1.0, 2.0];
    let x = [0.0, 1.0, 0.0];
    let spline = CubicSpline::build_with_boundaries(
        &s,
        &x,
        SplineBoundaryCondition::SpecifiedSlope,
        SplineBoundaryCondition::SpecifiedSlope,
        2.0,
        -2.0,
    )
    .unwrap();
    assert_abs_diff_eq!(spline.evaluate_derivative(0.0).unwrap(), 2.0, epsilon = TOL);
    assert_abs_diff_eq!(
        spline.evaluate_derivative(2.0).unwrap(),
        -2.0,
        epsilon = TOL
    );
    for (parameter, value) in s.iter().zip(&x) {
        assert_abs_diff_eq!(spline.evaluate(*parameter).unwrap(), value, epsilon = TOL);
    }
}

/// Ports `TwoPointSpline`.
///
/// Two knots is the minimum a spline accepts, and the interior tridiagonal
/// system is empty there — the result has to degrade to the straight line.
#[test]
fn a_two_point_spline_is_the_straight_line() {
    let spline = CubicSpline::build(&[0.0, 1.0], &[0.0, 2.0]).unwrap();
    assert_abs_diff_eq!(spline.evaluate(0.5).unwrap(), 1.0, epsilon = TOL);
    assert_abs_diff_eq!(spline.evaluate_derivative(0.5).unwrap(), 2.0, epsilon = TOL);
}

/// Ports `NonUniformSpacing`.
///
/// Unequal intervals are the case a naive uniform-step tridiagonal gets wrong,
/// so this checks knot interpolation and that curvature is continuous across
/// an interior knot where the interval length jumps by 4x.
#[test]
fn non_uniform_spacing_interpolates_and_stays_curvature_continuous() {
    let s = [0.0, 0.1, 0.5, 2.0, 4.0];
    let x = [1.0, 1.2, 2.0, 3.0, 2.5];
    let spline = CubicSpline::build(&s, &x).unwrap();
    for (parameter, value) in s.iter().zip(&x) {
        assert_abs_diff_eq!(spline.evaluate(*parameter).unwrap(), value, epsilon = TOL);
    }
    let left = spline.evaluate_second_derivative(0.499).unwrap();
    let right = spline.evaluate_second_derivative(0.501).unwrap();
    assert_abs_diff_eq!(left, right, epsilon = 0.01);
}

/// Ports `ExtrapolationBehavior`.
#[test]
fn evaluation_outside_the_knot_range_stays_finite() {
    let spline = CubicSpline::build(&[0.0, 1.0, 2.0, 3.0], &[0.0, 1.0, 4.0, 9.0]).unwrap();
    assert!(spline.evaluate(-0.5).unwrap().is_finite());
    assert!(spline.evaluate(3.5).unwrap().is_finite());
}

/// Ports `InvalidInputThrows`, as a typed error rather than an exception.
#[test]
fn mismatched_or_too_short_input_is_rejected() {
    assert!(CubicSpline::build(&[0.0, 1.0, 2.0], &[0.0, 1.0]).is_err());
    assert!(CubicSpline::build(&[0.0], &[1.0]).is_err());
}

/// Ports `EvaluateBeforeBuildThrows`.
///
/// The Rust type cannot be observed unbuilt through `build`, which returns
/// `Result`, so the unbuilt state is reached through `Default` — and every
/// evaluator must refuse it rather than index an empty knot vector.
#[test]
fn evaluating_an_unbuilt_spline_is_rejected() {
    let spline = CubicSpline::default();
    assert!(!spline.is_valid());
    assert!(spline.evaluate(0.5).is_err());
    assert!(spline.evaluate_derivative(0.5).is_err());
    assert!(spline.evaluate_second_derivative(0.5).is_err());
    assert!(spline.evaluate_all(0.5).is_err());
}

/// Ports `ValidityAndSizeQueries`.
#[test]
fn validity_size_and_parameter_range_queries() {
    let spline = CubicSpline::default();
    assert!(!spline.is_valid());
    assert_eq!(spline.len(), 0);
    assert!(spline.is_empty());

    let spline = CubicSpline::build(&[0.0, 1.0, 2.0, 3.0], &[1.0, 2.0, 3.0, 4.0]).unwrap();
    assert!(spline.is_valid());
    assert_eq!(spline.len(), 4);
    assert!(!spline.is_empty());
    assert_abs_diff_eq!(spline.min_parameter(), 0.0, epsilon = TOL);
    assert_abs_diff_eq!(spline.max_parameter(), 3.0, epsilon = TOL);
}
