//! Ported from the frozen v0.2.0 `solver_test.cpp`.
//!
//! The primal solve: AIC assembly, unit responses, and the circulation field
//! itself. Several of these assert physics a linear method must satisfy —
//! circulation proportional to alpha, a lift slope near finite-wing theory,
//! and a spanwise loading close to elliptic on a rectangular planform.

use approx::assert_abs_diff_eq;
use pegasus_vlm::core::{DEGREES_TO_RADIANS, PI, SpacingType, SymmetryType, Vec3};
use pegasus_vlm::geometry::{Aircraft, Section, Surface};
use pegasus_vlm::solver::{FlightCondition, SymmetryConfig, VortexLatticeSolver};

const TOL: f64 = 1.0e-6;

/// The frozen fixture. `duplicated` builds a half wing mirrored about `y = 0`;
/// otherwise the full span is meshed directly.
fn rectangular_wing(
    span: f64,
    chord: f64,
    chordwise: usize,
    spanwise: usize,
    duplicated: bool,
) -> Aircraft {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(chordwise);
    wing.set_spanwise_panel_count(spanwise);
    wing.set_chordwise_spacing(SpacingType::Cosine);
    wing.set_spanwise_spacing(SpacingType::Cosine);
    if duplicated {
        wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), chord, 0.0));
        wing.add_section(Section::new(Vec3::new(0.0, span / 2.0, 0.0), chord, 0.0));
        wing.set_y_duplicate(0.0);
    } else {
        wing.add_section(Section::new(Vec3::new(0.0, -span / 2.0, 0.0), chord, 0.0));
        wing.add_section(Section::new(Vec3::new(0.0, span / 2.0, 0.0), chord, 0.0));
    }

    let mut aircraft = Aircraft::new("Test Wing");
    aircraft.add_surface(wing);
    aircraft.set_reference_area(span * chord);
    aircraft.set_reference_chord(chord);
    aircraft.set_reference_span(span);
    aircraft
}

fn solver_for(aircraft: &Aircraft) -> VortexLatticeSolver {
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(aircraft).unwrap();
    solver
}

fn at_alpha(degrees: f64) -> FlightCondition {
    FlightCondition {
        alpha: degrees * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    }
}

/// Sum of circulation over each strip, paired with the strip's spanwise station.
fn spanwise_circulation(
    solver: &VortexLatticeSolver,
    solution: &pegasus_vlm::solver::SolutionState,
) -> Vec<(f64, f64)> {
    solver
        .mesh()
        .strips()
        .iter()
        .map(|strip| {
            let circulation: f64 = (0..strip.num_vortices)
                .map(|k| solution.circulation[strip.first_vortex_index + k])
                .sum();
            (strip.leading_edge_center.y, circulation)
        })
        .collect()
}

/// Ports `SetAircraftGeneratesMesh`.
#[test]
fn setting_an_aircraft_generates_the_mesh() {
    let solver = solver_for(&rectangular_wing(4.0, 1.0, 4, 4, false));
    assert!(solver.mesh().vortex_count() > 0);
    assert!(solver.mesh().strip_count() > 0);
}

/// Ports `RebuildAICCreatesValidMatrix` and `GetMatrices`.
///
/// The AIC is square in the panel count; the velocity influence matrix is
/// `3N x N` because it carries a vector per panel pair.
#[test]
fn the_assembled_matrices_have_the_expected_shapes() {
    let mut solver = solver_for(&rectangular_wing(4.0, 1.0, 4, 4, false));
    solver.rebuild_aic().unwrap();
    assert!(solver.is_aic_valid());
    let n = solver.mesh().vortex_count();
    assert_eq!(solver.aic().shape(), (n, n));
    assert_eq!(solver.velocity_matrix().shape(), (3 * n, n));
}

/// Ports `SolveUnitResponsesComputes6Columns`.
#[test]
fn unit_responses_become_valid_once_solved() {
    let mut solver = solver_for(&rectangular_wing(4.0, 1.0, 4, 4, false));
    solver.solve_unit_responses().unwrap();
    assert!(solver.are_unit_responses_valid());
}

/// Ports `InitialStateInvalid`.
#[test]
fn a_fresh_solver_has_neither_an_aic_nor_unit_responses() {
    let solver = VortexLatticeSolver::new();
    assert!(!solver.is_aic_valid());
    assert!(!solver.are_unit_responses_valid());
}

/// Ports `SolveAtZeroAlphaGivesZeroCirculation`.
#[test]
fn a_flat_wing_at_zero_alpha_has_zero_circulation() {
    let mut solver = solver_for(&rectangular_wing(4.0, 1.0, 4, 4, false));
    let solution = solver.solve(&at_alpha(0.0)).unwrap();
    assert!(solution.is_valid);
    for value in solution.circulation.iter() {
        assert_abs_diff_eq!(*value, 0.0, epsilon = TOL);
    }
}

/// Ports `SolveAtPositiveAlphaGivesPositiveCirculation`.
#[test]
fn positive_alpha_gives_positive_total_circulation() {
    let mut solver = solver_for(&rectangular_wing(4.0, 1.0, 4, 4, false));
    let solution = solver.solve(&at_alpha(5.0)).unwrap();
    assert!(solution.is_valid);
    assert!(solution.circulation.sum() > 0.0);
}

/// Ports `CirculationLinearWithAlpha`.
///
/// The method is linear in the boundary condition, so doubling alpha doubles
/// the circulation to within a percent — the sharpest statement that the solve
/// is a linear system and not something else.
#[test]
fn circulation_is_linear_in_alpha() {
    let mut solver = solver_for(&rectangular_wing(4.0, 1.0, 6, 8, false));
    let low = solver.solve(&at_alpha(2.0)).unwrap().circulation.sum();
    let high = solver.solve(&at_alpha(4.0)).unwrap().circulation.sum();
    assert_abs_diff_eq!(high / low, 2.0, epsilon = 0.01);
}

/// Ports `FreestreamDirectionComputation` and `FreestreamDirectionWithBeta`.
#[test]
fn the_freestream_direction_follows_alpha_and_beta() {
    let direction = VortexLatticeSolver::freestream_direction(&FlightCondition::default());
    assert_abs_diff_eq!(direction, Vec3::new(1.0, 0.0, 0.0), epsilon = TOL);

    let direction = VortexLatticeSolver::freestream_direction(&FlightCondition {
        alpha: PI / 2.0,
        ..FlightCondition::default()
    });
    assert_abs_diff_eq!(direction, Vec3::new(0.0, 0.0, 1.0), epsilon = TOL);

    let direction = VortexLatticeSolver::freestream_direction(&FlightCondition {
        beta: PI / 2.0,
        ..FlightCondition::default()
    });
    assert_abs_diff_eq!(direction, Vec3::new(0.0, 1.0, 0.0), epsilon = TOL);

    // Positive sideslip must give a positive y component, not merely a nonzero
    // one.
    let direction = VortexLatticeSolver::freestream_direction(&FlightCondition {
        beta: 10.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    });
    assert!(direction.y > 0.0);
    assert!(direction.iter().all(|value| value.is_finite()));
}

/// Ports `MachNumberChangeInvalidatesAIC`.
#[test]
fn changing_mach_rebuilds_and_still_solves() {
    let mut solver = solver_for(&rectangular_wing(4.0, 1.0, 4, 4, false));
    solver
        .solve(&FlightCondition {
            mach: 0.0,
            ..at_alpha(5.0)
        })
        .unwrap();
    assert!(solver.is_aic_valid());

    let solution = solver
        .solve(&FlightCondition {
            mach: 0.5,
            ..at_alpha(5.0)
        })
        .unwrap();
    assert!(solution.is_valid);
}

/// Ports `HighAspectRatioWingLiftSlope`.
///
/// The lift slope is recomputed here from the circulation field by
/// Kutta-Joukowski and compared against `2 pi AR / (AR + 2)`. Both the relative
/// bound and the absolute bracket from the frozen case are kept: the absolute
/// one is what catches a result that is self-consistent but physically wrong.
#[test]
fn a_high_aspect_ratio_wing_matches_finite_wing_lift_slope() {
    let (span, chord) = (10.0, 1.0);
    let aspect_ratio = span / chord;
    let aircraft = rectangular_wing(span, chord, 8, 20, false);
    let mut solver = solver_for(&aircraft);
    let alpha = 1.0 * DEGREES_TO_RADIANS;
    let solution = solver.solve(&at_alpha(1.0)).unwrap();

    let reference_area = aircraft.reference_area();
    let lift: f64 = solver
        .mesh()
        .strips()
        .iter()
        .map(|strip| {
            let circulation: f64 = (0..strip.num_vortices)
                .map(|k| solution.circulation[strip.first_vortex_index + k])
                .sum();
            2.0 * circulation * strip.width / reference_area
        })
        .sum();
    let lift_slope = lift / alpha;

    let theory = 2.0 * PI * aspect_ratio / (aspect_ratio + 2.0);
    assert_abs_diff_eq!(lift_slope, theory, epsilon = 0.25 * theory);
    assert!(lift_slope > 4.5);
    assert!(lift_slope < 7.0);
}

/// Ports `EllipticLiftDistributionShape`.
///
/// A rectangular wing is not elliptic, but its loading stays within a factor
/// of two of the elliptic shape over the mid-span band the frozen case checked.
/// The tips are excluded because that is exactly where a rectangular planform
/// departs from elliptic.
#[test]
fn spanwise_loading_is_roughly_elliptic() {
    let span = 10.0;
    let aircraft = rectangular_wing(span, 1.0, 4, 20, false);
    let mut solver = solver_for(&aircraft);
    let solution = solver.solve(&at_alpha(5.0)).unwrap();

    let stations = spanwise_circulation(&solver, &solution);
    let peak = stations
        .iter()
        .map(|(_, circulation)| *circulation)
        .fold(0.0, f64::max);

    let mut checked = 0;
    for (y, circulation) in &stations {
        if y.abs() > span * 0.1 && y.abs() < span * 0.4 {
            let eta = 2.0 * y / span;
            let elliptic = peak * (1.0 - eta * eta).sqrt();
            assert!(*circulation > 0.5 * elliptic);
            assert!(*circulation < 1.5 * elliptic);
            checked += 1;
        }
    }
    assert!(checked > 0, "the mid-span band contained no strips");
}

/// Ports `YDuplicateProducesSymmetricCirculation`.
#[test]
fn a_duplicated_wing_carries_equal_circulation_on_both_halves() {
    let aircraft = rectangular_wing(4.0, 1.0, 4, 8, true);
    let mut solver = solver_for(&aircraft);
    let solution = solver.solve(&at_alpha(5.0)).unwrap();
    assert_eq!(solver.mesh().surface_count(), 2);

    let total = |index: usize| -> f64 {
        let info = &solver.mesh().surface_info()[index];
        (0..info.num_vortices)
            .map(|offset| solution.circulation[info.first_vortex_index + offset])
            .sum()
    };
    let (left, right) = (total(0), total(1));
    assert_abs_diff_eq!(left, right, epsilon = 0.01 * left.abs());
}

/// Ports `InducedVelocityAtControlPoints`.
///
/// Downwash: a lifting wing induces a mean negative `w` at its own control
/// points, which is the mechanism that produces induced drag.
#[test]
fn a_lifting_wing_induces_mean_downwash() {
    let mut solver = solver_for(&rectangular_wing(4.0, 1.0, 4, 4, false));
    let solution = solver.solve(&at_alpha(5.0)).unwrap();
    let panels = solver.mesh().vortex_count();
    assert_eq!(solution.induced_velocity.shape(), (panels, 3));

    let mean_w: f64 = (0..panels)
        .map(|row| solution.induced_velocity[(row, 2)])
        .sum::<f64>()
        / panels as f64;
    assert!(mean_w < 0.0);
}

/// Ports `GetMeshMutable`.
///
/// The C++ handed out a mutable reference to the lattice. The Rust scopes the
/// mutation through a closure so the caches derived from the mesh are
/// invalidated with it, rather than leaving a stale AIC behind a live
/// reference.
#[test]
fn the_mesh_can_be_mutated_through_a_scoped_closure() {
    let mut solver = solver_for(&rectangular_wing(4.0, 1.0, 4, 4, false));
    solver.rebuild_aic().unwrap();
    assert!(solver.is_aic_valid());

    let mut observed = 0;
    solver.mutate_mesh(|mesh| observed = mesh.vortex_count());
    assert!(observed > 0);
    // Mutating through the closure drops the derived state rather than letting
    // it go stale.
    assert!(!solver.is_aic_valid());
}

/// Ports `GetReferenceQuantities`.
#[test]
fn reference_quantities_default_to_unity() {
    let solver = VortexLatticeSolver::new();
    let reference = solver.reference_quantities();
    assert_abs_diff_eq!(reference.area, 1.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(reference.chord, 1.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(reference.span, 1.0, epsilon = 1.0e-10);
}

/// Ports `SetReferencePoint`.
#[test]
fn the_reference_point_round_trips() {
    let mut solver = VortexLatticeSolver::new();
    solver.set_reference_point(Vec3::new(0.25, 0.0, 0.0));
    assert_abs_diff_eq!(
        *solver.reference_point(),
        Vec3::new(0.25, 0.0, 0.0),
        epsilon = 1.0e-10
    );
}

/// Ports `SetSymmetry`.
#[test]
fn the_symmetry_configuration_round_trips() {
    let mut solver = VortexLatticeSolver::new();
    solver.set_symmetry(SymmetryConfig {
        y_symmetry: SymmetryType::Symmetric,
        ..SymmetryConfig::default()
    });
    assert_eq!(solver.symmetry().y_symmetry, SymmetryType::Symmetric);
}

/// Ports `GetSolution`.
#[test]
fn the_solver_retains_its_last_solution() {
    let mut solver = solver_for(&rectangular_wing(4.0, 1.0, 4, 4, false));
    solver.solve(&at_alpha(5.0)).unwrap();
    let solution = solver.solution();
    assert!(solution.is_valid);
    assert!(!solution.circulation.is_empty());
}
