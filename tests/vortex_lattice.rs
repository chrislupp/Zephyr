//! Ported from the frozen v0.2.0 `vortex_lattice_test.cpp`.
//!
//! Mesh generation: panel counts, quarter-chord and three-quarter-chord
//! placement, normals, mirroring, and the index bookkeeping that ties strips to
//! vortices to surfaces.
//!
//! Two cases from this file are not ported. `CoreRadiusFractions` exercised
//! `VortexLattice::{Get,Set}VortexCoreRadius*`, a second pair of core-radius
//! knobs that was never passed to assembly — setting them left the AIC
//! bit-identical while appearing to tune the regularization — and has been
//! removed. `MutableVortexAccess` asserted that a panel can be mutated in place
//! through a public reference; the Rust keeps `vortices_mut` crate-private so
//! the lattice cannot be desynchronised from the caches derived from it.

use approx::assert_abs_diff_eq;
use pegasus_vlm::core::{SpacingType, Vec3};
use pegasus_vlm::geometry::{Aircraft, Section, Surface};
use pegasus_vlm::solver::VortexLattice;

const TOL: f64 = 1.0e-6;

/// The frozen fixture: a half wing mirrored about `y = 0`, cosine spaced.
fn rectangular_wing(span: f64, chord: f64, chordwise: usize, spanwise: usize) -> Aircraft {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(chordwise);
    wing.set_spanwise_panel_count(spanwise);
    wing.set_chordwise_spacing(SpacingType::Cosine);
    wing.set_spanwise_spacing(SpacingType::Cosine);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), chord, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, span / 2.0, 0.0), chord, 0.0));
    wing.set_y_duplicate(0.0);

    let mut aircraft = Aircraft::new("Rectangular Wing");
    aircraft.add_surface(wing);
    aircraft.set_reference_area(span * chord);
    aircraft.set_reference_chord(chord);
    aircraft.set_reference_span(span);
    aircraft
}

fn tapered_wing(
    span: f64,
    root_chord: f64,
    taper_ratio: f64,
    chordwise: usize,
    spanwise: usize,
) -> Aircraft {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(chordwise);
    wing.set_spanwise_panel_count(spanwise);
    wing.set_chordwise_spacing(SpacingType::Cosine);
    wing.set_spanwise_spacing(SpacingType::Cosine);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), root_chord, 0.0));
    wing.add_section(Section::new(
        Vec3::new(0.0, span / 2.0, 0.0),
        root_chord * taper_ratio,
        0.0,
    ));
    wing.set_y_duplicate(0.0);

    let mut aircraft = Aircraft::new("Tapered Wing");
    aircraft.add_surface(wing);
    aircraft
}

/// A single undivided panel: one section pair, one chordwise, one spanwise.
fn single_panel_wing(span: f64, chord: f64, chordwise: usize, spanwise: usize) -> Aircraft {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(chordwise);
    wing.set_spanwise_panel_count(spanwise);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), chord, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, span, 0.0), chord, 0.0));
    let mut aircraft = Aircraft::new("Test");
    aircraft.add_surface(wing);
    aircraft
}

fn lattice(aircraft: &Aircraft) -> VortexLattice {
    let mut mesh = VortexLattice::new();
    mesh.generate_from_aircraft(aircraft).unwrap();
    mesh
}

/// Ports `EmptyAircraft`.
#[test]
fn an_empty_aircraft_meshes_to_nothing() {
    let mesh = lattice(&Aircraft::default());
    assert_eq!(mesh.vortex_count(), 0);
    assert_eq!(mesh.strip_count(), 0);
    assert_eq!(mesh.surface_count(), 0);
}

/// Ports `RectangularWingBasicCounts`.
///
/// A duplicated surface produces two entries, and every strip carries exactly
/// the chordwise count regardless of which half it belongs to.
#[test]
fn a_duplicated_wing_produces_two_surfaces_with_uniform_strips() {
    let chordwise = 6;
    let mesh = lattice(&rectangular_wing(4.0, 1.0, chordwise, 8));
    assert_eq!(mesh.surface_count(), 2);
    assert!(mesh.strip_count() > 0);
    assert!(mesh.vortex_count() > 0);
    for strip in mesh.strips() {
        assert_eq!(strip.num_vortices, chordwise);
    }
}

/// Ports `RectangularWingGeometry`.
#[test]
fn every_panel_lies_inside_the_planform() {
    let (span, chord) = (4.0, 1.0);
    let mesh = lattice(&rectangular_wing(span, chord, 4, 4));
    for strip in mesh.strips() {
        assert_abs_diff_eq!(strip.chord_center, chord, epsilon = TOL);
    }
    for vortex in mesh.vortices() {
        assert!(vortex.control_point.x >= 0.0);
        assert!(vortex.control_point.x <= chord);
        assert!(vortex.control_point.y >= -span / 2.0);
        assert!(vortex.control_point.y <= span / 2.0);
    }
}

/// Ports `VortexQuarterChord`.
///
/// The defining placement of the method: on a single undivided panel the bound
/// vortex sits at the quarter chord and the control point at three quarters.
#[test]
fn a_single_panel_places_the_vortex_and_control_point_correctly() {
    let chord = 1.0;
    let mesh = lattice(&single_panel_wing(2.0, chord, 1, 1));
    assert_eq!(mesh.vortex_count(), 1);
    let vortex = &mesh.vortices()[0];
    assert_abs_diff_eq!(vortex.bound_vortex_center.x, 0.25 * chord, epsilon = 0.1);
    assert_abs_diff_eq!(vortex.control_point.x, 0.75 * chord, epsilon = 0.1);
}

/// Ports `ControlPointThreeQuarterChord`.
///
/// Panel by panel, the control point must lie aft of its own bound vortex.
#[test]
fn every_control_point_lies_aft_of_its_bound_vortex() {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(4);
    wing.set_spanwise_panel_count(2);
    wing.set_chordwise_spacing(SpacingType::Uniform);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 2.0, 0.0), 1.0, 0.0));
    let mut aircraft = Aircraft::new("Test");
    aircraft.add_surface(wing);

    let mesh = lattice(&aircraft);
    for vortex in mesh.vortices() {
        assert!(vortex.control_point.x > vortex.bound_vortex_center.x);
    }
}

/// Ports `NormalVectorsFlatWing`.
#[test]
fn a_flat_wing_has_vertical_normals() {
    let mesh = lattice(&single_panel_wing(2.0, 1.0, 4, 4));
    for vortex in mesh.vortices() {
        assert_abs_diff_eq!(vortex.normal_at_control_point.z, 1.0, epsilon = 0.1);
    }
}

/// Ports `TwistedWingNormals`.
///
/// Twist is interpolated along the span, so a root-twisted wing has more twist
/// at the root strip than at the tip.
#[test]
fn twist_decays_from_root_to_tip() {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(4);
    wing.set_spanwise_panel_count(4);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 5.0));
    wing.add_section(Section::new(Vec3::new(0.0, 2.0, 0.0), 1.0, 0.0));
    let mut aircraft = Aircraft::new("Test");
    aircraft.add_surface(wing);

    let mesh = lattice(&aircraft);
    let root = &mesh.strips()[0];
    let tip = &mesh.strips()[mesh.strip_count() - 1];
    assert!(root.twist_angle.abs() > tip.twist_angle.abs());
}

/// Ports `YDuplicateSymmetry`.
///
/// The mirrored half must be a true reflection: same `x` and `z`, negated `y`,
/// strip for strip.
#[test]
fn the_mirrored_half_reflects_strip_for_strip() {
    let mesh = lattice(&rectangular_wing(4.0, 1.0, 4, 4));
    assert_eq!(mesh.surface_count(), 2);
    let original = &mesh.surface_info()[0];
    let mirror = &mesh.surface_info()[1];
    assert_eq!(original.num_vortices, mirror.num_vortices);
    assert_eq!(original.num_strips, mirror.num_strips);

    for offset in 0..original.num_strips {
        let left = &mesh.strips()[original.first_strip_index + offset];
        let right = &mesh.strips()[mirror.first_strip_index + offset];
        assert_abs_diff_eq!(
            left.leading_edge_center.y,
            -right.leading_edge_center.y,
            epsilon = TOL
        );
        assert_abs_diff_eq!(
            left.leading_edge_center.x,
            right.leading_edge_center.x,
            epsilon = TOL
        );
        assert_abs_diff_eq!(
            left.leading_edge_center.z,
            right.leading_edge_center.z,
            epsilon = TOL
        );
    }
}

/// Ports `TaperedWingChordDistribution`.
#[test]
fn a_tapered_wing_has_monotonically_shrinking_strips() {
    let (root_chord, taper_ratio) = (1.0, 0.5);
    let mesh = lattice(&tapered_wing(4.0, root_chord, taper_ratio, 4, 8));
    let info = &mesh.surface_info()[0];

    let mut previous = root_chord;
    for offset in 0..info.num_strips {
        let strip = &mesh.strips()[info.first_strip_index + offset];
        assert!(strip.chord_center <= previous + 0.01);
        previous = strip.chord_center;
    }
    let tip = &mesh.strips()[info.first_strip_index + info.num_strips - 1];
    assert_abs_diff_eq!(tip.chord_center, root_chord * taper_ratio, epsilon = 0.1);
}

/// Ports `SurfaceInfoConsistency`.
///
/// Every surface's vortex count must factor exactly as chordwise x strips, and
/// its index ranges must stay inside the arrays they point into.
#[test]
fn surface_info_ranges_are_internally_consistent() {
    let (chordwise, spanwise) = (6, 10);
    let mesh = lattice(&rectangular_wing(4.0, 1.0, chordwise, spanwise));
    for info in mesh.surface_info() {
        assert_eq!(info.num_chordwise, chordwise);
        assert_eq!(info.num_spanwise, spanwise);
        assert_eq!(info.num_vortices, info.num_chordwise * info.num_strips);
        assert!(info.first_vortex_index + info.num_vortices <= mesh.vortex_count());
        assert!(info.first_strip_index + info.num_strips <= mesh.strip_count());
    }
}

/// Ports `StripVortexIndicesConsistent`.
///
/// Each panel must name the strip that owns it, which is what lets strip loads
/// be accumulated by a single pass over the panels.
#[test]
fn every_panel_names_the_strip_that_owns_it() {
    let chordwise = 4;
    let mesh = lattice(&rectangular_wing(4.0, 1.0, chordwise, 5));
    for (index, strip) in mesh.strips().iter().enumerate() {
        assert_eq!(strip.num_vortices, chordwise);
        assert!(strip.first_vortex_index + strip.num_vortices <= mesh.vortex_count());
        for offset in 0..strip.num_vortices {
            assert_eq!(
                mesh.vortices()[strip.first_vortex_index + offset].strip_index,
                index
            );
        }
    }
}

/// Ports `ClearAndRegenerate`.
#[test]
fn clearing_and_regenerating_reproduces_the_same_lattice() {
    let aircraft = rectangular_wing(4.0, 1.0, 4, 4);
    let mut mesh = VortexLattice::new();
    mesh.generate_from_aircraft(&aircraft).unwrap();
    let (vortices, strips) = (mesh.vortex_count(), mesh.strip_count());

    mesh.clear();
    assert_eq!(mesh.vortex_count(), 0);
    assert_eq!(mesh.strip_count(), 0);

    mesh.generate_from_aircraft(&aircraft).unwrap();
    assert_eq!(mesh.vortex_count(), vortices);
    assert_eq!(mesh.strip_count(), strips);
}

/// Ports `GetMeshInputNullByDefault`.
///
/// The stored `MeshInput` is `Some` only for the direct-mesh path, so a lattice
/// built from an `Aircraft` — or not built at all — reports `None`.
#[test]
fn a_lattice_has_no_mesh_input_unless_built_from_one() {
    assert!(VortexLattice::new().mesh_input().is_none());
    assert!(
        lattice(&rectangular_wing(4.0, 1.0, 2, 2))
            .mesh_input()
            .is_none()
    );
}

/// Ports `MutableStripAccess`, `GetVorticesMutable` and `GetStripsMutable`.
///
/// The C++ reached these through mutable accessors; the assertions were about
/// content, not mutability, so they port onto the immutable slices.
#[test]
fn generated_vortex_and_strip_collections_are_populated() {
    let mesh = lattice(&rectangular_wing(4.0, 1.0, 2, 2));
    assert!(!mesh.vortices().is_empty());
    assert!(!mesh.strips().is_empty());
    assert!(mesh.strips()[0].num_vortices >= 1);
}

/// Ports `ControlSurfaceNamesEmpty` and `ControlSurfaceInfluencesSize`.
///
/// A wing with no controls still gets one influence row per panel — an empty
/// row, not a missing one.
#[test]
fn an_uncontrolled_wing_has_empty_but_present_influence_rows() {
    let mesh = lattice(&rectangular_wing(4.0, 1.0, 2, 2));
    assert!(mesh.control_surface_names().is_empty());
    assert_eq!(mesh.control_influences().len(), mesh.vortex_count());
}

/// Ports `RecalculateNormals`.
///
/// Recalculating from unchanged geometry has to be a no-op to the last bit, or
/// repeated calls would walk the normals.
#[test]
fn recalculating_normals_on_unchanged_geometry_is_idempotent() {
    let aircraft = rectangular_wing(4.0, 1.0, 2, 2);
    let mut mesh = VortexLattice::new();
    mesh.generate_from_aircraft(&aircraft).unwrap();
    let before = mesh.vortices()[0].normal_at_control_point;
    mesh.recalculate_normals();
    assert_abs_diff_eq!(
        (mesh.vortices()[0].normal_at_control_point - before).norm(),
        0.0,
        epsilon = 1.0e-10
    );
}
