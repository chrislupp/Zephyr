//! Ported from the frozen v0.2.0 `influence_matrix_test.cpp`.
//!
//! The Biot-Savart kernel itself and the matrices assembled from it. The kernel
//! cases pin the symmetry and sign structure that everything downstream
//! inherits; the matrix cases pin shapes, conditioning, and that every optional
//! path — compressibility, symmetry images, multi-component core radii — stays
//! finite.

use approx::assert_abs_diff_eq;
use pegasus_vlm::core::{MatX, SpacingType, SymmetryType, Vec3, VecX};
use pegasus_vlm::geometry::{Aircraft, Section, Surface};
use pegasus_vlm::solver::{
    InfluenceMatrixBuilder, SymmetryConfig, VortexLattice, horseshoe_vortex_velocity,
};

const TOL: f64 = 1.0e-6;

const BOUND_LEFT: Vec3 = Vec3::new(-0.5, 0.0, 0.0);
const BOUND_RIGHT: Vec3 = Vec3::new(0.5, 0.0, 0.0);

fn velocity(field: Vec3, pg: f64, core: f64, include_bound: bool) -> Vec3 {
    horseshoe_vortex_velocity(&field, &BOUND_LEFT, &BOUND_RIGHT, pg, core, include_bound)
}

// ============================================================================
// Biot-Savart kernel
// ============================================================================

/// Ports `VelocityAboveBoundVortex`.
#[test]
fn velocity_directly_above_the_bound_vortex_has_no_spanwise_component() {
    let result = velocity(Vec3::new(0.0, 0.0, 1.0), 1.0, 0.001, true);
    assert!(result.iter().all(|value| value.is_finite()));
    assert_abs_diff_eq!(result.y, 0.0, epsilon = 0.01);
}

/// Ports `VelocitySymmetry`.
///
/// Reflecting the field point in `y` flips `x` and `z` but preserves `y`. That
/// parity is what makes a mirrored panel's influence derivable by reflection
/// rather than requiring its own evaluation.
#[test]
fn the_kernel_has_the_expected_parity_under_a_y_reflection() {
    let positive = velocity(Vec3::new(0.0, 0.5, 0.5), 1.0, 0.001, true);
    let negative = velocity(Vec3::new(0.0, -0.5, 0.5), 1.0, 0.001, true);
    assert_abs_diff_eq!(positive.x, -negative.x, epsilon = TOL);
    assert_abs_diff_eq!(positive.y, negative.y, epsilon = TOL);
    assert_abs_diff_eq!(positive.z, -negative.z, epsilon = TOL);
}

/// Ports `CoreRegularization`.
///
/// A larger core spreads the vorticity further, so it must weaken the induced
/// velocity near the filament rather than merely keep it finite.
#[test]
fn a_larger_core_radius_weakens_the_induced_velocity() {
    let field = Vec3::new(1.0, 0.3, 0.02);
    let small = velocity(field, 1.0, 0.01, false);
    let large = velocity(field, 1.0, 0.5, false);
    assert!(small.norm().is_finite());
    assert!(large.norm().is_finite());
    assert!(small.norm() > 0.01);
    assert!(large.norm() < small.norm());
}

/// Ports `PrandtlGlauertCorrection`.
#[test]
fn the_prandtl_glauert_factor_changes_the_kernel() {
    let field = Vec3::new(1.0, 0.5, 0.5);
    let incompressible = velocity(field, 1.0, 0.001, true);
    let compressible = velocity(field, 0.75f64.sqrt(), 0.001, true);
    assert!((compressible.norm() - incompressible.norm()).abs() > 0.0);
}

/// Ports `ExcludeBoundVortex`.
#[test]
fn excluding_the_bound_segment_changes_the_result() {
    let field = Vec3::new(0.0, 0.3, 0.5);
    let with_bound = velocity(field, 1.0, 0.001, true);
    let without_bound = velocity(field, 1.0, 0.001, false);
    assert!(with_bound.norm().is_finite());
    assert!(without_bound.norm().is_finite());
    assert!((with_bound.norm() - without_bound.norm()).abs() > 0.0);
}

/// Ports `TrailingVortexInDownstreamField`.
///
/// Far downstream, between the two trailing legs, the induced velocity is
/// downward — the two-dimensional limit that gives rise to downwash.
#[test]
fn the_field_far_downstream_between_the_legs_is_downwash() {
    let result = velocity(Vec3::new(10.0, 0.5, 0.0), 1.0, 0.001, true);
    assert!(result.z < 0.0);
}

// ============================================================================
// Assembled matrices
// ============================================================================

fn rectangular_wing(span: f64, chord: f64, chordwise: usize, spanwise: usize) -> Aircraft {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(chordwise);
    wing.set_spanwise_panel_count(spanwise);
    wing.set_chordwise_spacing(SpacingType::Uniform);
    wing.set_spanwise_spacing(SpacingType::Uniform);
    wing.add_section(Section::new(Vec3::new(0.0, -span / 2.0, 0.0), chord, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, span / 2.0, 0.0), chord, 0.0));

    let mut aircraft = Aircraft::new("Test Wing");
    aircraft.add_surface(wing);
    aircraft.set_reference_area(span * chord);
    aircraft.set_reference_chord(chord);
    aircraft.set_reference_span(span);
    aircraft
}

fn lattice(aircraft: &Aircraft) -> VortexLattice {
    let mut mesh = VortexLattice::new();
    mesh.generate_from_aircraft(aircraft).unwrap();
    mesh
}

fn all_finite(matrix: &MatX) -> bool {
    matrix.iter().all(|value| value.is_finite())
}

/// Ports `VelocityMatrixDimensions` and `AICMatrixDimensions`.
#[test]
fn assembled_matrices_have_the_expected_shapes() {
    let aircraft = rectangular_wing(4.0, 1.0, 4, 4);
    let mesh = lattice(&aircraft);
    let builder = InfluenceMatrixBuilder::new(&mesh);
    let n = mesh.vortex_count();
    assert_eq!(builder.velocity_influence_matrix().shape(), (3 * n, n));
    assert_eq!(builder.normalwash_matrix().shape(), (n, n));
}

/// Ports `AICDiagonalDominance`.
///
/// The frozen case only asserted a nonzero diagonal, not true dominance; a
/// panel always influences itself, which is what keeps the system solvable.
#[test]
fn every_aic_diagonal_entry_is_nonzero() {
    let aircraft = rectangular_wing(4.0, 1.0, 4, 4);
    let mesh = lattice(&aircraft);
    let aic = InfluenceMatrixBuilder::new(&mesh).normalwash_matrix();
    for index in 0..mesh.vortex_count() {
        assert!(aic[(index, index)].abs() > 0.0);
    }
}

/// Ports `AICMatrixNonsingular`.
#[test]
fn the_aic_is_nonsingular() {
    let aircraft = rectangular_wing(4.0, 1.0, 4, 4);
    let mesh = lattice(&aircraft);
    let aic = InfluenceMatrixBuilder::new(&mesh).normalwash_matrix();
    let determinant = aic.determinant();
    assert!(determinant.is_finite());
    assert!(determinant != 0.0);
}

/// Ports `MachNumberEffect`.
#[test]
fn compressibility_moves_the_assembled_aic() {
    let aircraft = rectangular_wing(4.0, 1.0, 2, 2);
    let mesh = lattice(&aircraft);
    let mut builder = InfluenceMatrixBuilder::new(&mesh);
    builder.set_mach_number(0.0);
    let incompressible = builder.normalwash_matrix();
    builder.set_mach_number(0.5);
    let compressible = builder.normalwash_matrix();
    assert!((compressible - incompressible).norm() > 0.0);
}

/// Ports `SupersonicMachClamped`.
///
/// The builder clamps the Prandtl-Glauert factor to one at and above Mach 1
/// rather than taking the square root of a negative number. `solve` rejects
/// supersonic input outright; this is the lower-level guard that keeps the
/// assembled matrix finite regardless.
#[test]
fn a_supersonic_mach_still_assembles_a_finite_matrix() {
    let aircraft = rectangular_wing(4.0, 1.0, 2, 2);
    let mesh = lattice(&aircraft);
    let mut builder = InfluenceMatrixBuilder::new(&mesh);
    builder.set_mach_number(1.5);
    assert_abs_diff_eq!(builder.prandtl_glauert_factor(), 1.0, epsilon = TOL);
    assert!(all_finite(&builder.normalwash_matrix()));
}

/// Ports `InducedVelocityFromMatrix`.
///
/// Every panel's self-influence block must be nonzero, or that panel would
/// contribute nothing to its own normalwash.
#[test]
fn every_panel_has_a_nonzero_self_influence_block() {
    let aircraft = rectangular_wing(4.0, 1.0, 2, 2);
    let mesh = lattice(&aircraft);
    let (_, velocity_matrix) = InfluenceMatrixBuilder::new(&mesh).normalwash_and_velocity();
    assert!(all_finite(&velocity_matrix));
    for index in 0..mesh.vortex_count() {
        let block = Vec3::new(
            velocity_matrix[(3 * index, index)],
            velocity_matrix[(3 * index + 1, index)],
            velocity_matrix[(3 * index + 2, index)],
        );
        assert!(block.norm() > 0.0);
    }
}

/// Ports `ComputeInducedVelocity`.
#[test]
fn a_uniformly_loaded_lattice_induces_downwash_above_itself() {
    let aircraft = rectangular_wing(4.0, 1.0, 2, 4);
    let mesh = lattice(&aircraft);
    let builder = InfluenceMatrixBuilder::new(&mesh);
    let circulation = VecX::from_element(mesh.vortex_count(), 0.1);
    let result = builder.induced_velocity(&Vec3::new(0.5, 0.0, 0.5), &circulation);
    assert!(result.iter().all(|value| value.is_finite()));
    assert!(result.z < 0.0);
}

/// Ports `ComputeSingleVortexVelocity`.
#[test]
fn a_single_panel_contribution_is_finite() {
    let aircraft = rectangular_wing(4.0, 1.0, 2, 4);
    let mesh = lattice(&aircraft);
    let builder = InfluenceMatrixBuilder::new(&mesh);
    let result = builder.single_vortex_velocity(&Vec3::new(0.5, 0.0, 1.0), 0, true);
    assert!(result.iter().all(|value| value.is_finite()));
}

/// Builds the symmetric-wing fixture shared by the two symmetry cases.
fn symmetric_wing() -> Aircraft {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(2);
    wing.set_spanwise_panel_count(4);
    wing.set_y_duplicate(0.0);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 2.0, 0.0), 1.0, 0.0));

    let mut aircraft = Aircraft::new("Symmetric Wing");
    aircraft.add_surface(wing);
    aircraft.set_reference_area(4.0);
    aircraft.set_reference_chord(1.0);
    aircraft.set_reference_span(4.0);
    aircraft.set_y_symmetry(SymmetryType::Symmetric);
    aircraft
}

fn symmetry_of(aircraft: &Aircraft) -> SymmetryConfig {
    SymmetryConfig {
        y_symmetry: aircraft.y_symmetry(),
        z_symmetry: aircraft.z_symmetry(),
        z_symmetry_plane: aircraft.z_symmetry_plane(),
        ..SymmetryConfig::default()
    }
}

/// Ports `YSymmetryInfluence`.
#[test]
fn y_symmetry_images_assemble_a_finite_matrix() {
    let aircraft = symmetric_wing();
    let mesh = lattice(&aircraft);
    let mut builder = InfluenceMatrixBuilder::new(&mesh);
    builder.set_symmetry(symmetry_of(&aircraft));
    assert!(all_finite(&builder.normalwash_matrix()));
}

/// Ports `InducedVelocityWithYSymmetry`.
#[test]
fn y_symmetry_images_give_a_finite_induced_velocity() {
    let aircraft = symmetric_wing();
    let mesh = lattice(&aircraft);
    let mut builder = InfluenceMatrixBuilder::new(&mesh);
    builder.set_symmetry(symmetry_of(&aircraft));
    let circulation = VecX::from_element(mesh.vortex_count(), 0.1);
    let result = builder.induced_velocity(&Vec3::new(0.5, 0.0, 0.5), &circulation);
    assert!(result.iter().all(|value| value.is_finite()));
}

/// Ports `DifferentComponentCoreRadius`.
///
/// Panels on different components use a larger core radius than panels on the
/// same one, so a two-component aircraft exercises the other branch of the
/// core-radius selection.
#[test]
fn a_two_component_aircraft_assembles_a_finite_matrix() {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(2);
    wing.set_spanwise_panel_count(2);
    wing.set_component_index(0);
    wing.add_section(Section::new(Vec3::new(0.0, -2.0, 0.0), 1.0, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 2.0, 0.0), 1.0, 0.0));

    let mut tail = Surface::new("Tail");
    tail.set_chordwise_panel_count(2);
    tail.set_spanwise_panel_count(2);
    tail.set_component_index(1);
    tail.add_section(Section::new(Vec3::new(4.0, -1.0, 0.0), 0.5, 0.0));
    tail.add_section(Section::new(Vec3::new(4.0, 1.0, 0.0), 0.5, 0.0));

    let mut aircraft = Aircraft::new("Two Component");
    aircraft.add_surface(wing);
    aircraft.add_surface(tail);
    aircraft.set_reference_area(4.0);
    aircraft.set_reference_chord(1.0);
    aircraft.set_reference_span(4.0);

    let mesh = lattice(&aircraft);
    assert!(all_finite(
        &InfluenceMatrixBuilder::new(&mesh).normalwash_matrix()
    ));
}
