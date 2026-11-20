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

//! Ported from the frozen v0.2.0 `analytical_aic_test.cpp`.
//!
//! Gradients of a single AIC entry and of weighted rows and columns. The C++
//! took raw coordinate arrays; the Rust row and column functions take a real
//! lattice, so the fixture here is the geometry that produces exactly the
//! frozen arrays: a chord-1 flat wing spanning `y` in `[0, 1]`, one panel
//! chordwise and two spanwise, which puts bound vortices at `x = 0.25` and
//! control points at `x = 0.75`.

use approx::assert_abs_diff_eq;
use zephyr_vlm::core::{SpacingType, Vec3, VecX};
use zephyr_vlm::geometry::{Aircraft, Section, Surface};
use zephyr_vlm::gradients::{
    aic_column_gradient, aic_element_gradient, aic_row_gradient, velocity_column_gradient,
    velocity_row_gradient,
};
use zephyr_vlm::solver::{VortexLattice, horseshoe_vortex_velocity};

const FD_STEP: f64 = 1.0e-6;
const FD_TOL: f64 = 1.0e-4;

const BOUND_LEFT: Vec3 = Vec3::new(0.25, 0.0, 0.0);
const BOUND_RIGHT: Vec3 = Vec3::new(0.25, 1.0, 0.0);
const CONTROL: Vec3 = Vec3::new(0.75, 0.5, 0.0);
const NORMAL: Vec3 = Vec3::new(0.0, 0.0, 1.0);
const CORE: f64 = 0.01;

/// One AIC entry: the induced velocity contracted with the target normal.
fn aic_element(
    control: Vec3,
    normal: Vec3,
    left: Vec3,
    right: Vec3,
    pg: f64,
    core: f64,
    include_bound: bool,
) -> f64 {
    horseshoe_vortex_velocity(&control, &left, &right, pg, core, include_bound).dot(&normal)
}

/// Compares an exact derivative against its centered difference at the frozen
/// tolerance, which is relative with a small absolute floor so a component that
/// differences to nearly zero is not held to a bound it cannot meet.
fn assert_close(exact: f64, difference: f64) {
    assert_abs_diff_eq!(
        exact,
        difference,
        epsilon = FD_TOL * (difference.abs() + 1.0e-8)
    );
}

/// Ports `GradientsMatchFD`.
///
/// All four gradient blocks at once. The normal block is exact rather than
/// differenced: the entry is linear in the normal, so its derivative is the
/// induced velocity itself.
#[test]
fn the_aic_element_gradient_matches_finite_differences() {
    let pg = 1.0;
    let gradient =
        aic_element_gradient(&BOUND_LEFT, &BOUND_RIGHT, &CONTROL, &NORMAL, CORE, pg, true);

    for component in 0..3 {
        let mut plus = CONTROL;
        let mut minus = CONTROL;
        plus[component] += FD_STEP;
        minus[component] -= FD_STEP;
        let difference = (aic_element(plus, NORMAL, BOUND_LEFT, BOUND_RIGHT, pg, CORE, true)
            - aic_element(minus, NORMAL, BOUND_LEFT, BOUND_RIGHT, pg, CORE, true))
            / (2.0 * FD_STEP);
        assert_close(gradient.control_point[component], difference);
    }

    for component in 0..3 {
        let mut plus = BOUND_LEFT;
        let mut minus = BOUND_LEFT;
        plus[component] += FD_STEP;
        minus[component] -= FD_STEP;
        let difference = (aic_element(CONTROL, NORMAL, plus, BOUND_RIGHT, pg, CORE, true)
            - aic_element(CONTROL, NORMAL, minus, BOUND_RIGHT, pg, CORE, true))
            / (2.0 * FD_STEP);
        assert_close(gradient.bound_left[component], difference);
    }

    for component in 0..3 {
        let mut plus = BOUND_RIGHT;
        let mut minus = BOUND_RIGHT;
        plus[component] += FD_STEP;
        minus[component] -= FD_STEP;
        let difference = (aic_element(CONTROL, NORMAL, BOUND_LEFT, plus, pg, CORE, true)
            - aic_element(CONTROL, NORMAL, BOUND_LEFT, minus, pg, CORE, true))
            / (2.0 * FD_STEP);
        assert_close(gradient.bound_right[component], difference);
    }

    // d(n . v)/dn = v exactly.
    let velocity = horseshoe_vortex_velocity(&CONTROL, &BOUND_LEFT, &BOUND_RIGHT, pg, CORE, true);
    assert_abs_diff_eq!(gradient.normal, velocity, epsilon = 1.0e-12);
}

/// Ports `WithoutBound`.
#[test]
fn the_aic_element_gradient_matches_with_the_bound_segment_excluded() {
    let gradient = aic_element_gradient(
        &BOUND_LEFT,
        &BOUND_RIGHT,
        &CONTROL,
        &NORMAL,
        CORE,
        1.0,
        false,
    );
    for component in 0..3 {
        let mut plus = CONTROL;
        let mut minus = CONTROL;
        plus[component] += FD_STEP;
        minus[component] -= FD_STEP;
        let difference = (aic_element(plus, NORMAL, BOUND_LEFT, BOUND_RIGHT, 1.0, CORE, false)
            - aic_element(minus, NORMAL, BOUND_LEFT, BOUND_RIGHT, 1.0, CORE, false))
            / (2.0 * FD_STEP);
        assert_close(gradient.control_point[component], difference);
    }
}

/// Ports `PrandtlGlauert`.
#[test]
fn the_aic_element_gradient_matches_under_compressibility() {
    let pg = 0.8;
    let gradient =
        aic_element_gradient(&BOUND_LEFT, &BOUND_RIGHT, &CONTROL, &NORMAL, CORE, pg, true);
    for component in 0..3 {
        let mut plus = CONTROL;
        let mut minus = CONTROL;
        plus[component] += FD_STEP;
        minus[component] -= FD_STEP;
        let difference = (aic_element(plus, NORMAL, BOUND_LEFT, BOUND_RIGHT, pg, CORE, true)
            - aic_element(minus, NORMAL, BOUND_LEFT, BOUND_RIGHT, pg, CORE, true))
            / (2.0 * FD_STEP);
        assert_close(gradient.control_point[component], difference);
    }
}

// ============================================================================
// Weighted rows and columns
// ============================================================================

/// The geometry that reproduces the frozen two-panel arrays.
fn two_panel_lattice() -> VortexLattice {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(1);
    wing.set_spanwise_panel_count(2);
    wing.set_chordwise_spacing(SpacingType::Uniform);
    wing.set_spanwise_spacing(SpacingType::Uniform);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 1.0, 0.0), 1.0, 0.0));
    let mut aircraft = Aircraft::new("Two Panel");
    aircraft.add_surface(wing);

    let mut mesh = VortexLattice::new();
    mesh.generate_from_aircraft(&aircraft).unwrap();
    mesh
}

/// Same-component core radius, matching what the row and column gradients use.
fn same_component_core(mesh: &VortexLattice, source: usize) -> f64 {
    let vortex = &mesh.vortices()[source];
    0.0001 * (vortex.bound_vortex_right - vortex.bound_vortex_left).norm()
}

/// Ports `RowGradientMatchesFD`.
///
/// The weighted row is differenced with respect to the target control point,
/// summing over both sources with the frozen weights `[1.0, 0.5]`. The
/// self-influence term drops its bound segment, so the sum exercises both
/// branches in one difference.
#[test]
fn the_weighted_row_gradient_matches_finite_differences() {
    let mesh = two_panel_lattice();
    assert_eq!(mesh.vortex_count(), 2);
    let weights = VecX::from_vec(vec![1.0, 0.5]);
    let row = 0;
    let (control_gradient, _) = aic_row_gradient(&mesh, row, &weights, 1.0, 0.25, 0.05);

    let target = &mesh.vortices()[row];
    let weighted_row = |control: Vec3| -> f64 {
        mesh.vortices()
            .iter()
            .enumerate()
            .map(|(source_index, source)| {
                weights[source_index]
                    * aic_element(
                        control,
                        target.normal_at_control_point,
                        source.bound_vortex_left,
                        source.bound_vortex_right,
                        1.0,
                        same_component_core(&mesh, source_index),
                        row != source_index,
                    )
            })
            .sum()
    };

    for component in 0..3 {
        let mut plus = target.control_point;
        let mut minus = target.control_point;
        plus[component] += FD_STEP;
        minus[component] -= FD_STEP;
        let difference = (weighted_row(plus) - weighted_row(minus)) / (2.0 * FD_STEP);
        assert_close(control_gradient[component], difference);
    }
}

/// Ports `ColumnGradientFinite`.
#[test]
fn the_weighted_column_gradient_is_finite() {
    let mesh = two_panel_lattice();
    let weights = VecX::from_vec(vec![1.0, 0.5]);
    let (left, right) = aic_column_gradient(&mesh, 0, &weights, 1.0, 0.25, 0.05);
    assert!(left.iter().all(|value| value.is_finite()));
    assert!(right.iter().all(|value| value.is_finite()));
}

/// Ports `DifferentComponents`.
///
/// Panels on different components take the larger core radius, so this reaches
/// the other branch of the selection.
#[test]
fn the_row_gradient_is_finite_across_components() {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(1);
    wing.set_spanwise_panel_count(1);
    wing.set_component_index(0);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 0.5, 0.0), 1.0, 0.0));

    let mut tail = Surface::new("Tail");
    tail.set_chordwise_panel_count(1);
    tail.set_spanwise_panel_count(1);
    tail.set_component_index(1);
    tail.add_section(Section::new(Vec3::new(0.0, 0.5, 0.0), 1.0, 0.0));
    tail.add_section(Section::new(Vec3::new(0.0, 1.0, 0.0), 1.0, 0.0));

    let mut aircraft = Aircraft::new("Two Component");
    aircraft.add_surface(wing);
    aircraft.add_surface(tail);
    let mut mesh = VortexLattice::new();
    mesh.generate_from_aircraft(&aircraft).unwrap();

    let weights = VecX::from_vec(vec![1.0, 0.5]);
    let (control, normal) = aic_row_gradient(&mesh, 0, &weights, 1.0, 0.25, 0.05);
    assert!(control.iter().all(|value| value.is_finite()));
    assert!(normal.iter().all(|value| value.is_finite()));
}

/// Ports `VelocityRowGradientFinite`.
#[test]
fn the_velocity_row_gradient_is_finite() {
    let mesh = two_panel_lattice();
    let circulation = VecX::from_vec(vec![1.0, 0.5]);
    let cotangent = Vec3::new(0.1, -0.2, 0.5);
    let gradient = velocity_row_gradient(&mesh, 0, &circulation, &cotangent, 1.0, 0.25, 0.05);
    assert!(gradient.iter().all(|value| value.is_finite()));
}

/// Ports `VelocityColumnGradientFinite`.
#[test]
fn the_velocity_column_gradient_is_finite() {
    let mesh = two_panel_lattice();
    let cotangents = vec![Vec3::new(0.0, 0.0, 1.0), Vec3::new(0.0, 0.0, 0.5)];
    let (left, right) = velocity_column_gradient(&mesh, 0, 1.0, &cotangents, 1.0, 0.25, 0.05);
    assert!(left.iter().all(|value| value.is_finite()));
    assert!(right.iter().all(|value| value.is_finite()));
}

/// Ports `VelocityColumnGradientZeroGamma`.
///
/// A column's contribution scales with its own circulation, so a zero
/// circulation must give an exactly zero gradient — not merely a small one.
#[test]
fn a_zero_circulation_column_has_an_exactly_zero_gradient() {
    let mesh = two_panel_lattice();
    let cotangents = vec![Vec3::new(1.0, 1.0, 1.0), Vec3::new(1.0, 1.0, 1.0)];
    let (left, right) = velocity_column_gradient(&mesh, 0, 0.0, &cotangents, 1.0, 0.25, 0.05);
    assert_eq!(left, Vec3::zeros());
    assert_eq!(right, Vec3::zeros());
}
