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

//! Ported from the frozen v0.2.0 `mesh_sensitivity_test.cpp`.
//!
//! The C++ built section-geometry gradients out of free functions operating on
//! flattened coordinate arrays: `ExtractMeshArrays`, `BuildSectionStripMapping`,
//! `MeshGradients`, `PropagateLeadingEdgeGradient` and so on. The Rust reaches
//! the same result through `section_geometry` and typed structs, so those
//! helpers have no one-to-one counterpart.
//!
//! Each case below therefore reproduces the *assertion* against the Rust
//! construct that carries the same information, and the doc comment names the
//! substitution. Where the C++ helper's behaviour was replaced rather than
//! relocated — `PropagateTwistGradient`, whose small-angle linearization the
//! Rust deliberately does not use — the case says so and asserts the
//! replacement's own property instead.

use approx::assert_abs_diff_eq;
use zephyr_vlm::core::{DEGREES_TO_RADIANS, SpacingType, Vec3};
use zephyr_vlm::geometry::{Aircraft, Section, Surface};
use zephyr_vlm::gradients::{GradientCalculator, interpolated_twist_jvp, panel_normal_jvp};
use zephyr_vlm::solver::{FlightCondition, VortexLattice};
use zephyr_vlm::utils::{cosine_spacing, uniform_spacing};

/// The frozen fixture: a 4x8 duplicated wing of span 10.
fn aircraft() -> Aircraft {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(4);
    wing.set_spanwise_panel_count(8);
    wing.set_y_duplicate(0.0);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 5.0, 0.0), 1.0, 0.0));

    let mut aircraft = Aircraft::new("Test Wing");
    aircraft.add_surface(wing);
    aircraft.set_reference_area(10.0);
    aircraft.set_reference_chord(1.0);
    aircraft.set_reference_span(10.0);
    aircraft
}

fn lattice(aircraft: &Aircraft) -> VortexLattice {
    let mut mesh = VortexLattice::new();
    mesh.generate_from_aircraft(aircraft).unwrap();
    mesh
}

/// Ports `ExtractMeshArrays`.
///
/// The C++ flattened the lattice into parallel `3N` coordinate arrays for its
/// kernels; the Rust passes the lattice itself. The assertion is the same: one
/// finite entry per panel for every quantity the gradients read.
#[test]
fn every_panel_carries_finite_geometry_and_an_owning_surface() {
    let aircraft = aircraft();
    let mesh = lattice(&aircraft);
    assert!(mesh.vortex_count() > 0);
    for vortex in mesh.vortices() {
        for value in [
            vortex.bound_vortex_left,
            vortex.bound_vortex_right,
            vortex.control_point,
            vortex.normal_at_control_point,
        ] {
            assert!(value.iter().all(|component| component.is_finite()));
        }
        assert!(vortex.local_chord.is_finite());
        assert!(vortex.surface_index < mesh.surface_count());
    }
}

/// Ports `MeshSensitivityMultiSurfaceTest.TwoSurfaces`.
///
/// Panels must be attributable to both surfaces, and each surface must own a
/// nonempty strip range — the mapping the C++ built explicitly.
#[test]
fn a_two_surface_aircraft_attributes_panels_and_strips_to_both() {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(4);
    wing.set_spanwise_panel_count(6);
    wing.set_y_duplicate(0.0);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 5.0, 0.0), 1.0, 0.0));

    let mut tail = Surface::new("Tail");
    tail.set_chordwise_panel_count(4);
    tail.set_spanwise_panel_count(4);
    tail.set_y_duplicate(0.0);
    tail.add_section(Section::new(Vec3::new(4.0, 0.0, 0.0), 0.5, 0.0));
    tail.add_section(Section::new(Vec3::new(4.0, 1.5, 0.0), 0.3, 0.0));

    let mut aircraft = Aircraft::new("Wing-Tail");
    aircraft.add_surface(wing);
    aircraft.add_surface(tail);
    aircraft.set_reference_area(10.0);
    aircraft.set_reference_chord(1.0);
    aircraft.set_reference_span(10.0);

    let mesh = lattice(&aircraft);
    let surfaces: std::collections::HashSet<usize> =
        mesh.vortices().iter().map(|v| v.surface_index).collect();
    assert!(surfaces.contains(&0));
    assert!(surfaces.contains(&1));
    for info in mesh.surface_info() {
        assert!(info.num_strips > 0);
    }
}

/// Ports `BuildSectionStripMappingSingleSurface` and `...YDuplicate`.
///
/// Both the original surface and its mirror own strips, and the first surface's
/// strip range starts at zero.
#[test]
fn each_generated_surface_owns_a_contiguous_strip_range() {
    let aircraft = aircraft();
    let mesh = lattice(&aircraft);
    // The wing is duplicated, so the lattice carries the original and its
    // mirror.
    assert_eq!(mesh.surface_count(), 2);
    assert_eq!(mesh.surface_info()[0].first_strip_index, 0);

    let mut expected_first = 0;
    for info in mesh.surface_info() {
        assert!(info.num_strips > 0);
        assert_eq!(info.first_strip_index, expected_first);
        expected_first += info.num_strips;
    }
    assert_eq!(expected_first, mesh.strip_count());
    // Two sections define the surface, whichever half a strip belongs to.
    assert_eq!(aircraft.surfaces()[0].sections().len(), 2);
}

/// Ports `SectionStripMappingGetStripRange`.
///
/// The C++ derived a per-section strip range from a `strip_at_section` table.
/// The Rust equivalent is the `SurfaceMeshInfo` ranges, which must tile the
/// strip array exactly once with no gap and no overlap — the property the
/// index arithmetic existed to provide.
#[test]
fn surface_strip_ranges_tile_the_strip_array_exactly() {
    let aircraft = aircraft();
    let mesh = lattice(&aircraft);
    let mut covered = vec![0usize; mesh.strip_count()];
    for info in mesh.surface_info() {
        for offset in 0..info.num_strips {
            covered[info.first_strip_index + offset] += 1;
        }
    }
    assert!(covered.iter().all(|count| *count == 1));
}

/// Ports `BuildDefaultChordwiseFractions`.
///
/// The C++ helper produced the per-panel quarter-chord and three-quarter-chord
/// fractions. The Rust bakes them into the mesher, so this asserts the property
/// on the generated panels instead: within every strip the bound vortex is
/// ahead of its control point, both lie inside the chord, and the first panel's
/// vortex sits in the front 20%.
#[test]
fn chordwise_placement_puts_each_vortex_ahead_of_its_control_point() {
    let aircraft = aircraft();
    let mesh = lattice(&aircraft);
    for strip in mesh.strips() {
        let leading_edge = strip.leading_edge_center.x;
        let chord = strip.chord_center;
        for offset in 0..strip.num_vortices {
            let vortex = &mesh.vortices()[strip.first_vortex_index + offset];
            let vortex_fraction = (vortex.bound_vortex_center.x - leading_edge) / chord;
            let control_fraction = (vortex.control_point.x - leading_edge) / chord;
            assert!(vortex_fraction < control_fraction);
            assert!(vortex_fraction >= 0.0);
            assert!(control_fraction <= 1.0 + 1.0e-12);
            if offset == 0 {
                assert!(vortex_fraction < 0.2);
            }
        }
    }
}

/// Ports `StripInterpolationFraction`.
///
/// A single strip interpolates at its midpoint, and across eight strips the
/// first sits near the root end and the last near the tip.
#[test]
fn strip_interpolation_fractions_span_the_interval() {
    // One strip: the centre of a uniform division of [0, 1].
    let single = uniform_spacing(2);
    assert_abs_diff_eq!(0.5 * (single[0] + single[1]), 0.5, epsilon = 1.0e-12);

    let stations = uniform_spacing(9);
    let centre = |index: usize| 0.5 * (stations[index] + stations[index + 1]);
    assert!(centre(0) < 0.15);
    assert!(centre(7) > 0.85);
}

/// Ports `StripEdgeFractions`.
#[test]
fn strip_edge_fractions_are_ordered_and_reach_the_tip() {
    let stations = uniform_spacing(5);
    assert!(stations[0] >= 0.0);
    assert!(stations[1] <= 1.0);
    assert!(stations[0] < stations[1]);
    // The last strip's outboard edge is the tip.
    assert_abs_diff_eq!(stations[4], 1.0, epsilon = 1.0e-10);
}

/// Ports `BuildSpanwiseFractions`.
///
/// Cosine spacing still spans `[0, 1]`, and every strip centre lies between its
/// own two edges.
#[test]
fn cosine_spanwise_fractions_bracket_their_centres() {
    let stations = cosine_spacing(5);
    assert_eq!(stations.len(), 5);
    assert_abs_diff_eq!(stations[0], 0.0, epsilon = 1.0e-10);
    assert_abs_diff_eq!(stations[4], 1.0, epsilon = 1.0e-10);
    for index in 0..4 {
        let centre = 0.5 * (stations[index] + stations[index + 1]);
        assert!(centre >= stations[index]);
        assert!(centre <= stations[index + 1]);
    }
}

/// Ports `MeshGradientsResizeAndZero`.
///
/// The C++ carried a mutable scratch buffer that had to be sized and cleared by
/// hand. The Rust returns typed gradient structs, so the equivalent property is
/// that a default-constructed set is fully zero — there is no separate `Zero()`
/// step that could be forgotten.
#[test]
fn default_constructed_section_gradients_are_zero() {
    use zephyr_vlm::gradients::SectionGeometryGradients;
    let gradients = SectionGeometryGradients::default();
    assert!(gradients.surface_name.is_empty());
    assert_eq!(gradients.section_index, 0);
    for group in [
        &gradients.by_chord,
        &gradients.by_twist,
        &gradients.by_leading_edge_x,
        &gradients.by_leading_edge_y,
        &gradients.by_leading_edge_z,
    ] {
        assert_abs_diff_eq!(group.coefficients.lift, 0.0);
        assert_abs_diff_eq!(group.coefficients.drag, 0.0);
        assert_abs_diff_eq!(group.coefficients.side, 0.0);
        assert_abs_diff_eq!(group.coefficients.roll, 0.0);
        assert_abs_diff_eq!(group.coefficients.pitch, 0.0);
        assert_abs_diff_eq!(group.coefficients.yaw, 0.0);
    }
}

/// Ports `PropagateLeadingEdgeGradient`.
///
/// A quantity interpolated at fraction `f` between two section stations
/// distributes its sensitivity as `(1 - f)` inboard and `f` outboard. The Rust
/// applies that inside `section_geometry`; `interpolated_twist_jvp` exposes the
/// same partition, and for equal chords it reduces to exactly those weights.
#[test]
fn interval_interpolation_splits_sensitivity_as_one_minus_f_and_f() {
    let fraction = 0.3;
    // Equal chords and equal twists: the chord weighting drops out and the
    // split is purely linear in the interpolation fraction.
    // (fraction, chord_l, chord_r, twist_l, twist_r, dchord_l, dchord_r,
    //  dtwist_l, dtwist_r) — perturb one end's twist at a time.
    let inboard = interpolated_twist_jvp(fraction, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    let outboard = interpolated_twist_jvp(fraction, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    assert_abs_diff_eq!(inboard, 1.0 - fraction, epsilon = 1.0e-12);
    assert_abs_diff_eq!(outboard, fraction, epsilon = 1.0e-12);
    // The two halves account for the whole perturbation.
    assert_abs_diff_eq!(inboard + outboard, 1.0, epsilon = 1.0e-12);
}

/// Ports `PropagateChordGradient`.
///
/// Both bounding sections' chords must influence a strip that lies between
/// them, so neither derivative may come out zero.
#[test]
fn both_bounding_section_chords_influence_an_interior_strip() {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(2);
    wing.set_spanwise_panel_count(4);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 4.0, 0.0), 0.5, 0.0));
    let mut aircraft = Aircraft::new("Tapered");
    aircraft.add_surface(wing);
    aircraft.set_reference_area(3.0);
    aircraft.set_reference_chord(0.75);
    aircraft.set_reference_span(4.0);

    let mut calculator = GradientCalculator::new();
    calculator.set_aircraft(&aircraft);
    let gradients = calculator
        .compute_geometry_gradients(&FlightCondition {
            alpha: 5.0 * DEGREES_TO_RADIANS,
            ..FlightCondition::default()
        })
        .unwrap();

    assert_eq!(gradients.sections.len(), 2);
    for section in &gradients.sections {
        assert!(section.by_chord.coefficients.lift != 0.0);
        assert!(section.by_chord.coefficients.lift.is_finite());
    }
}

/// Ports `ComputeSectionGradientFromMesh`.
///
/// The full per-section gradient: every design variable produces a finite
/// derivative for every output.
#[test]
fn section_gradients_are_finite_for_every_design_variable() {
    let aircraft = aircraft();
    let mut calculator = GradientCalculator::new();
    calculator.set_aircraft(&aircraft);
    let gradients = calculator
        .compute_geometry_gradients(&FlightCondition {
            alpha: 5.0 * DEGREES_TO_RADIANS,
            ..FlightCondition::default()
        })
        .unwrap();

    assert!(!gradients.sections.is_empty());
    for section in &gradients.sections {
        for group in [
            &section.by_chord,
            &section.by_twist,
            &section.by_leading_edge_x,
            &section.by_leading_edge_y,
            &section.by_leading_edge_z,
        ] {
            for value in [
                group.coefficients.lift,
                group.coefficients.drag,
                group.coefficients.side,
                group.coefficients.roll,
                group.coefficients.pitch,
                group.coefficients.yaw,
            ] {
                assert!(value.is_finite());
            }
        }
    }
}

/// Ports `PropagateTwistGradient`, whose formula the Rust deliberately does not
/// reproduce.
///
/// v0.2.0 propagated twist through a small-angle linearization that is exact
/// only when the chords and twists match at both ends of an interval. The Rust
/// uses the exact chord-weighted form and exports `panel_normal_jvp`, so this
/// asserts the replacement's own property: a twist perturbation rotates the
/// panel normal by a finite amount, orthogonally to the normal itself, since a
/// rotation cannot change its length.
#[test]
fn a_twist_perturbation_rotates_the_panel_normal_orthogonally() {
    let left = Vec3::new(0.25, 0.0, 0.0);
    let right = Vec3::new(0.25, 1.0, 0.0);
    let derivative = panel_normal_jvp(
        &left,
        &right,
        &Vec3::zeros(),
        &Vec3::zeros(),
        0.05,
        1.0,
        0.0,
        (0.0, 1.0),
        (0.0, 0.0),
    );
    assert!(derivative.iter().all(|value| value.is_finite()));
    assert!(derivative.norm() > 0.0);

    // The unperturbed normal, built exactly as the primal does: the chord
    // vector is `(cos a, -sin a * n_y, -sin a * n_z)` for incidence `a` and
    // Trefftz normal `(n_y, n_z)`, and the panel normal is `chord x bound`.
    let angle = 0.05f64;
    let unit_bound = (right - left).normalize();
    let chord = Vec3::new(angle.cos(), -angle.sin() * 0.0, -angle.sin() * 1.0);
    let normal = chord.cross(&unit_bound).normalize();
    // The normal is a unit vector, so its derivative is orthogonal to it: a
    // rotation cannot change its length.
    assert_abs_diff_eq!(normal.norm(), 1.0, epsilon = 1.0e-12);
    assert_abs_diff_eq!(derivative.dot(&normal), 0.0, epsilon = 1.0e-9);
}

/// Ports `MeshSensitivityTest` spacing coverage for the non-uniform case.
#[test]
fn cosine_spacing_clusters_more_tightly_than_uniform_at_the_root() {
    let uniform = uniform_spacing(9);
    let cosine = cosine_spacing(9);
    assert!(cosine[1] - cosine[0] < uniform[1] - uniform[0]);
    assert_eq!(
        zephyr_vlm::utils::spacing(9, SpacingType::Cosine).len(),
        cosine.len()
    );
}
