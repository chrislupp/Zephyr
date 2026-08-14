//! Ported from the frozen v0.2.0 `analytical_biot_savart_test.cpp`.
//!
//! Every case here validates a hand-derived VJP against a centered difference
//! of the primal. That is the strongest form of gradient test available: the
//! two sides share no code, so agreement cannot come from a common error.
//!
//! These now guard the closed-form kernel that replaced the forward-mode dual
//! this port previously carried, and they cover the branches the dual version
//! got wrong — with and without the bound segment, compressible and
//! incompressible, small and large core.

use approx::assert_abs_diff_eq;
use pegasus_vlm::core::Vec3;
use pegasus_vlm::gradients::{
    horseshoe_velocity_vjp_bound_left, horseshoe_velocity_vjp_bound_right,
    horseshoe_velocity_vjp_field_point,
};
use pegasus_vlm::solver::horseshoe_vortex_velocity;

const FD_STEP: f64 = 1.0e-6;
/// The frozen relative tolerance; the C++ compared against
/// `kFDTol * (|fd| + floor)`.
const FD_TOL: f64 = 1.0e-4;

const FIELD: Vec3 = Vec3::new(0.5, 0.5, 0.0);
const BOUND_LEFT: Vec3 = Vec3::new(0.25, 0.0, 0.0);
const BOUND_RIGHT: Vec3 = Vec3::new(0.25, 1.0, 0.0);

/// Which of the three inputs a difference perturbs.
#[derive(Clone, Copy)]
enum Wrt {
    FieldPoint,
    BoundLeft,
    BoundRight,
}

fn primal(field: Vec3, left: Vec3, right: Vec3, pg: f64, core: f64, include_bound: bool) -> Vec3 {
    horseshoe_vortex_velocity(&field, &left, &right, pg, core, include_bound)
}

/// Centered difference of `seed . velocity` with respect to one input.
fn finite_difference(wrt: Wrt, pg: f64, core: f64, include_bound: bool, seed: Vec3) -> Vec3 {
    Vec3::from_fn(|component, _| {
        let mut plus = [FIELD, BOUND_LEFT, BOUND_RIGHT];
        let mut minus = plus;
        let slot = match wrt {
            Wrt::FieldPoint => 0,
            Wrt::BoundLeft => 1,
            Wrt::BoundRight => 2,
        };
        plus[slot][component] += FD_STEP;
        minus[slot][component] -= FD_STEP;
        let evaluate =
            |points: [Vec3; 3]| primal(points[0], points[1], points[2], pg, core, include_bound);
        (evaluate(plus) - evaluate(minus)).dot(&seed) / (2.0 * FD_STEP)
    })
}

fn analytical(wrt: Wrt, pg: f64, core: f64, include_bound: bool, seed: Vec3) -> Vec3 {
    let function = match wrt {
        Wrt::FieldPoint => horseshoe_velocity_vjp_field_point,
        Wrt::BoundLeft => horseshoe_velocity_vjp_bound_left,
        Wrt::BoundRight => horseshoe_velocity_vjp_bound_right,
    };
    function(
        &FIELD,
        &BOUND_LEFT,
        &BOUND_RIGHT,
        pg,
        core,
        include_bound,
        &seed,
    )
}

/// Compares one VJP against its difference at the frozen tolerance.
///
/// `floor` is the absolute term the C++ added so a component that differences
/// to nearly zero is not held to a relative bound it cannot meet.
fn assert_vjp_matches(wrt: Wrt, pg: f64, core: f64, include_bound: bool, seed: Vec3, floor: f64) {
    let exact = analytical(wrt, pg, core, include_bound, seed);
    let difference = finite_difference(wrt, pg, core, include_bound, seed);
    for component in 0..3 {
        assert_abs_diff_eq!(
            exact[component],
            difference[component],
            epsilon = FD_TOL * (difference[component].abs() + floor)
        );
    }
}

fn unit_z() -> Vec3 {
    Vec3::new(0.0, 0.0, 1.0)
}

/// Ports `PrimalNonZero`.
///
/// The differences below are meaningless if the primal is zero at this
/// configuration, so this establishes there is something to differentiate.
#[test]
fn the_primal_is_finite_and_nonzero_at_the_test_point() {
    let velocity = primal(FIELD, BOUND_LEFT, BOUND_RIGHT, 1.0, 0.01, true);
    assert!(velocity.iter().all(|value| value.is_finite()));
    assert!(velocity.norm() > 1.0e-10);
}

/// Ports `FieldPointVJPMatchesFD`.
#[test]
fn the_field_point_vjp_matches_a_centered_difference() {
    assert_vjp_matches(Wrt::FieldPoint, 1.0, 0.01, true, unit_z(), 1.0e-8);
}

/// Ports `BoundLeftVJPMatchesFD`.
#[test]
fn the_bound_left_vjp_matches_a_centered_difference() {
    assert_vjp_matches(Wrt::BoundLeft, 1.0, 0.01, true, unit_z(), 1.0e-8);
}

/// Ports `BoundRightVJPMatchesFD`.
#[test]
fn the_bound_right_vjp_matches_a_centered_difference() {
    assert_vjp_matches(Wrt::BoundRight, 1.0, 0.01, true, unit_z(), 1.0e-8);
}

/// Ports `WithoutBoundVortex`.
///
/// The self-influence path: with the bound segment excluded only the two
/// trailing legs remain, so this isolates their derivative from the bound
/// term's.
#[test]
fn all_three_vjps_match_with_the_bound_segment_excluded() {
    for wrt in [Wrt::FieldPoint, Wrt::BoundLeft, Wrt::BoundRight] {
        assert_vjp_matches(wrt, 1.0, 0.01, false, unit_z(), 1.0e-8);
    }
}

/// Ports `PrandtlGlauertScaling`.
///
/// A compressible factor stretches only the x coordinate, so it exercises the
/// asymmetric chain through the Prandtl-Glauert transform.
#[test]
fn all_three_vjps_match_under_compressibility() {
    for wrt in [Wrt::FieldPoint, Wrt::BoundLeft, Wrt::BoundRight] {
        assert_vjp_matches(wrt, 0.8, 0.01, true, unit_z(), 1.0e-4);
    }
}

/// Ports `LargerCoreRadius`.
///
/// A ten-fold core radius moves the regularization terms from negligible to
/// dominant, which is where the `core^4` chain has to be right.
#[test]
fn the_field_point_vjp_matches_with_a_large_core() {
    assert_vjp_matches(Wrt::FieldPoint, 1.0, 0.1, true, unit_z(), 1.0e-8);
}

/// Ports `ArbitrarySeedDirection`.
///
/// The unit-z seed used above would hide an error confined to the x or y rows
/// of the Jacobian; a general seed contracts all three.
#[test]
fn the_field_point_vjp_matches_for_an_arbitrary_seed() {
    assert_vjp_matches(
        Wrt::FieldPoint,
        1.0,
        0.01,
        true,
        Vec3::new(1.0, -0.5, 0.3),
        1.0e-8,
    );
}

/// Ports `BiotSavartAICTest.AICElementKernel`.
///
/// An AIC entry is exactly the induced velocity contracted with the target
/// normal. The C++ had a dedicated kernel for this; the Rust composes it from
/// the primal, so the identity is asserted directly and exactly.
#[test]
fn an_aic_element_is_the_velocity_contracted_with_the_normal() {
    let field = Vec3::new(0.75, 0.5, 0.0);
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let velocity = primal(field, BOUND_LEFT, BOUND_RIGHT, 1.0, 0.01, true);
    let element = velocity.dot(&normal);
    assert!(element.is_finite());
    assert!(element != 0.0);
    assert_eq!(element, velocity.z);
}

/// Ports `BiotSavartAICTest.AICElementKernelWithoutBound`.
#[test]
fn excluding_the_bound_segment_changes_the_aic_element() {
    let field = Vec3::new(0.75, 0.5, 0.0);
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let with_bound = primal(field, BOUND_LEFT, BOUND_RIGHT, 1.0, 0.01, true).dot(&normal);
    let without_bound = primal(field, BOUND_LEFT, BOUND_RIGHT, 1.0, 0.01, false).dot(&normal);
    assert!(without_bound.is_finite());
    assert!(with_bound != without_bound);
}
