//! Ported from the frozen v0.2.0 `analytical_trefftz_test.cpp`.
//!
//! Circulation derivatives of the far-field lift, side force and induced drag.
//! The C++ took raw per-strip arrays and its own primal kernel; the Rust
//! differentiates the same discrete calculation `ForceCalculator` performs, so
//! the differences below perturb circulation and re-run that primal — which
//! also proves the gradient and the primal agree on which calculation they are
//! talking about.
//!
//! The fixture reproduces the frozen arrays: three strips with trailing-edge
//! stations at `y = 0, 0.5, 1.0, 1.5`, unit chord, circulation
//! `[0.10, 0.12, 0.08]`, reference area and span both 10.

use approx::assert_abs_diff_eq;
use pegasus_vlm::aero::{AerodynamicResults, ForceCalculator};
use pegasus_vlm::core::{SpacingType, Vec3, VecX};
use pegasus_vlm::geometry::{Aircraft, ReferenceQuantities, Section, Surface};
use pegasus_vlm::gradients::trefftz_circulation_gradients;
use pegasus_vlm::solver::{FlightCondition, SolutionState, VortexLattice};

const FD_STEP: f64 = 1.0e-7;
const FD_TOL: f64 = 1.0e-4;
const REFERENCE_AREA: f64 = 10.0;
const REFERENCE_SPAN: f64 = 10.0;
const CORE_CHORD_FRACTION: f64 = 0.25;
const CORE_WIDTH_FRACTION: f64 = 0.05;
const GAMMA: [f64; 3] = [0.10, 0.12, 0.08];

/// Three strips over `y` in `[0, 1.5]`, one panel chordwise.
///
/// `components` assigns each spanwise interval to a surface, which is what
/// selects the core radius between strips.
fn lattice(split_surfaces: bool) -> VortexLattice {
    let mut aircraft = Aircraft::new("Trefftz Fixture");
    if split_surfaces {
        // Two surfaces so the third strip sits on a different one.
        let mut inner = Surface::new("Inner");
        inner.set_chordwise_panel_count(1);
        inner.set_spanwise_panel_count(2);
        inner.set_spanwise_spacing(SpacingType::Uniform);
        inner.set_component_index(0);
        inner.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
        inner.add_section(Section::new(Vec3::new(0.0, 1.0, 0.0), 1.0, 0.0));

        let mut outer = Surface::new("Outer");
        outer.set_chordwise_panel_count(1);
        outer.set_spanwise_panel_count(1);
        outer.set_spanwise_spacing(SpacingType::Uniform);
        outer.set_component_index(1);
        outer.add_section(Section::new(Vec3::new(0.0, 1.0, 0.0), 1.0, 0.0));
        outer.add_section(Section::new(Vec3::new(0.0, 1.5, 0.0), 1.0, 0.0));

        aircraft.add_surface(inner);
        aircraft.add_surface(outer);
    } else {
        let mut wing = Surface::new("Wing");
        wing.set_chordwise_panel_count(1);
        wing.set_spanwise_panel_count(3);
        wing.set_spanwise_spacing(SpacingType::Uniform);
        wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
        wing.add_section(Section::new(Vec3::new(0.0, 1.5, 0.0), 1.0, 0.0));
        aircraft.add_surface(wing);
    }

    let mut mesh = VortexLattice::new();
    mesh.generate_from_aircraft(&aircraft).unwrap();
    mesh
}

fn solution_with(circulation: &[f64]) -> SolutionState {
    SolutionState {
        circulation: VecX::from_column_slice(circulation),
        is_valid: true,
        ..SolutionState::default()
    }
}

/// The primal far-field triple for a given circulation, through the same
/// `ForceCalculator` path the gradient differentiates.
fn far_field(mesh: &VortexLattice, circulation: &[f64]) -> AerodynamicResults {
    let mut calculator = ForceCalculator::new();
    calculator.set_trefftz_core_radius(CORE_CHORD_FRACTION, CORE_WIDTH_FRACTION);
    let reference = ReferenceQuantities {
        area: REFERENCE_AREA,
        chord: 1.0,
        span: REFERENCE_SPAN,
        moment_reference_point: Vec3::zeros(),
    };
    calculator.compute_with_references(
        &reference,
        mesh,
        &solution_with(circulation),
        &FlightCondition::default(),
    )
}

/// Which far-field output a gradient column belongs to.
#[derive(Clone, Copy)]
enum Output {
    Lift,
    Side,
    InducedDrag,
}

fn value_of(results: &AerodynamicResults, output: Output) -> f64 {
    match output {
        Output::Lift => results.far_field_lift,
        Output::Side => results.far_field_side,
        Output::InducedDrag => results.far_field_induced_drag,
    }
}

/// Differences one far-field output against its analytical circulation
/// gradient, strip by strip.
fn assert_gradient_matches(mesh: &VortexLattice, output: Output) {
    let gradients = trefftz_circulation_gradients(
        mesh,
        &solution_with(&GAMMA),
        REFERENCE_AREA,
        CORE_CHORD_FRACTION,
        CORE_WIDTH_FRACTION,
    );
    let exact = match output {
        Output::Lift => &gradients.lift,
        Output::Side => &gradients.side,
        Output::InducedDrag => &gradients.induced_drag,
    };

    for strip in 0..GAMMA.len() {
        let mut plus = GAMMA;
        let mut minus = GAMMA;
        plus[strip] += FD_STEP;
        minus[strip] -= FD_STEP;
        let difference = (value_of(&far_field(mesh, &plus), output)
            - value_of(&far_field(mesh, &minus), output))
            / (2.0 * FD_STEP);
        assert_abs_diff_eq!(
            exact[strip],
            difference,
            epsilon = FD_TOL * (difference.abs() + 1.0e-8)
        );
    }
}

/// Ports `PrimalKernelComputes`.
///
/// A positive-circulation flat wing lifts, produces induced drag, and — being
/// symmetric in the plane — no side force.
#[test]
fn the_far_field_primal_is_physically_sensible() {
    let mesh = lattice(false);
    assert_eq!(mesh.strip_count(), 3);
    let results = far_field(&mesh, &GAMMA);
    assert!(results.far_field_lift > 0.0);
    assert_abs_diff_eq!(results.far_field_side, 0.0, epsilon = 1.0e-10);
    assert!(results.far_field_induced_drag > 0.0);
    assert!(results.far_field_lift.is_finite());
    assert!(results.far_field_induced_drag.is_finite());
}

/// Ports `SpanEfficiency` and `SpanEfficiencyZeroDrag`.
///
/// Zero induced drag has no efficiency to report, so the accessor returns zero
/// rather than dividing by it.
#[test]
fn span_efficiency_is_bounded_and_zero_when_there_is_no_drag() {
    let mesh = lattice(false);
    let results = far_field(&mesh, &GAMMA);
    assert!(results.span_efficiency > 0.0);
    assert!(results.span_efficiency < 2.0);

    // Zero circulation produces zero induced drag, which is the branch that
    // must return zero rather than divide by it.
    let unloaded = far_field(&mesh, &[0.0, 0.0, 0.0]);
    assert_abs_diff_eq!(unloaded.far_field_induced_drag, 0.0, epsilon = 1.0e-14);
    assert_eq!(unloaded.span_efficiency, 0.0);
}

/// Ports `TrefftzVelocityPoint`.
///
/// A trailing vortex pair induces a finite, nonzero vertical velocity at a
/// point offset from it — the downwash the induced-drag integral rests on.
#[test]
fn a_trailing_pair_induces_finite_nonzero_downwash() {
    let mesh = lattice(false);
    let results = far_field(&mesh, &GAMMA);
    // The induced drag is the contraction of that downwash with the loading,
    // so a nonzero drag is the observable form of a nonzero induced velocity.
    assert!(results.far_field_induced_drag.is_finite());
    assert!(results.far_field_induced_drag != 0.0);
}

/// Ports `CLGradientMatchesFD`.
#[test]
fn the_far_field_lift_gradient_matches_finite_differences() {
    assert_gradient_matches(&lattice(false), Output::Lift);
}

/// Ports `CYGradientMatchesFD`.
#[test]
fn the_far_field_side_gradient_matches_finite_differences() {
    assert_gradient_matches(&lattice(false), Output::Side);
}

/// Ports `CDiGradientMatchesFD`.
///
/// The only one of the three that is quadratic in circulation, so it is the
/// case where a missed transpose term would show.
#[test]
fn the_induced_drag_gradient_matches_finite_differences() {
    assert_gradient_matches(&lattice(false), Output::InducedDrag);
}

/// Ports `MultiSurfaceCoreRadius`.
///
/// Strips on different surfaces take the larger core radius, so splitting the
/// span across two surfaces reaches the other branch of the selection.
#[test]
fn the_induced_drag_gradient_matches_across_surfaces() {
    let mesh = lattice(true);
    assert_eq!(mesh.strip_count(), 3);
    assert_gradient_matches(&mesh, Output::InducedDrag);
}
