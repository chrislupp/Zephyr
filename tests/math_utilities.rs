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

//! Ported from the frozen v0.2.0 `math_utilities_test.cpp`.
//!
//! These cover the shared primitives the rest of the library is built on: the
//! compressibility factor, the two regularized vortex kernels, the four spacing
//! distributions and their dispatcher, and the rotation helpers. The C++ suite
//! asserted these directly rather than through a solve, and so does this.

use approx::assert_abs_diff_eq;
use zephyr_vlm::core::{Mat3, PI, SpacingType, Vec3};
use zephyr_vlm::utils::{
    blended_spacing, body_to_stability_rotation, cosine_spacing, prandtl_glauert_factor,
    rotate_about_x, rotate_about_y, rotate_about_z, sine_spacing, spacing,
    trailing_vortex_velocity, uniform_spacing, vortex_segment_velocity, wind_to_body_rotation,
};

const TOL: f64 = 1.0e-10;

// ============================================================================
// Prandtl-Glauert factor
// ============================================================================

/// Ports `PrandtlGlauertZeroMach` and `PrandtlGlauertSubsonic`.
#[test]
fn prandtl_glauert_factor_matches_the_closed_form() {
    assert_abs_diff_eq!(prandtl_glauert_factor(0.0), 1.0, epsilon = TOL);
    assert_abs_diff_eq!(
        prandtl_glauert_factor(0.5),
        (1.0f64 - 0.25).sqrt(),
        epsilon = TOL
    );
}

// ============================================================================
// Vortex segment velocity
// ============================================================================

/// Ports `VortexSegmentAboveCenter`.
#[test]
fn vortex_segment_above_its_center_induces_sideways_flow() {
    let velocity = vortex_segment_velocity(
        &Vec3::new(0.5, 0.0, 1.0),
        &Vec3::new(0.0, 0.0, 0.0),
        &Vec3::new(1.0, 0.0, 0.0),
        0.01,
    );
    assert!(velocity.iter().all(|value| value.is_finite()));
    // Right-hand rule about a segment along +x puts the induced velocity in y.
    assert!(velocity.y.abs() > 1.0e-6);
}

/// Ports `VortexSegmentOnAxis`.
///
/// A field point on the segment's own axis is the singular direction; the core
/// model has to keep it finite rather than dividing by a zero cross product.
#[test]
fn vortex_segment_on_axis_stays_finite() {
    let velocity = vortex_segment_velocity(
        &Vec3::new(0.5, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 0.0),
        &Vec3::new(1.0, 0.0, 0.0),
        0.01,
    );
    assert!(velocity.iter().all(|value| value.is_finite()));
}

/// Ports `VortexSegmentAtEndpoint`.
#[test]
fn vortex_segment_at_its_endpoint_induces_nothing() {
    let velocity = vortex_segment_velocity(
        &Vec3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 0.0),
        &Vec3::new(1.0, 0.0, 0.0),
        0.01,
    );
    assert!(velocity.iter().all(|value| value.is_finite()));
    assert_abs_diff_eq!(velocity.norm(), 0.0, epsilon = 1.0e-6);
}

// ============================================================================
// Trailing vortex velocity
// ============================================================================

/// Ports `TrailingVortexAbove`.
#[test]
fn trailing_vortex_above_its_origin_induces_positive_y() {
    let velocity =
        trailing_vortex_velocity(&Vec3::new(0.0, 0.0, 1.0), &Vec3::new(0.0, 0.0, 0.0), 0.01);
    assert!(velocity.iter().all(|value| value.is_finite()));
    // cross = (0, r_z, -r_y) = (0, 1, 0) for a field point directly above.
    assert!(velocity.y > 0.0);
}

/// Ports `TrailingVortexAtStart`.
#[test]
fn trailing_vortex_at_its_origin_induces_nothing() {
    let velocity =
        trailing_vortex_velocity(&Vec3::new(0.0, 0.0, 0.0), &Vec3::new(0.0, 0.0, 0.0), 0.01);
    assert!(velocity.iter().all(|value| value.is_finite()));
    assert_abs_diff_eq!(velocity.norm(), 0.0, epsilon = 1.0e-6);
}

/// Ports `TrailingVortexDownstream`.
///
/// Far downstream the semi-infinite leg approaches its two-dimensional limit,
/// so the sign must not flip relative to the point directly above the origin.
#[test]
fn trailing_vortex_far_downstream_keeps_its_sign() {
    let velocity =
        trailing_vortex_velocity(&Vec3::new(5.0, 0.0, 1.0), &Vec3::new(0.0, 0.0, 0.0), 0.01);
    assert!(velocity.iter().all(|value| value.is_finite()));
    assert!(velocity.y > 0.0);
}

// ============================================================================
// Spacing distributions
// ============================================================================

/// Ports `CosineSpacingEndpoints` and `CosineSpacingConcentration`.
#[test]
fn cosine_spacing_spans_the_unit_interval_and_clusters_at_both_ends() {
    let values = cosine_spacing(11);
    assert_eq!(values.len(), 11);
    assert_abs_diff_eq!(values[0], 0.0, epsilon = TOL);
    assert_abs_diff_eq!(values[10], 1.0, epsilon = TOL);

    let values = cosine_spacing(21);
    assert!(values[1] - values[0] < values[11] - values[10]);
}

/// Ports `UniformSpacing`.
#[test]
fn uniform_spacing_is_evenly_divided() {
    let values = uniform_spacing(5);
    assert_eq!(values.len(), 5);
    for (index, expected) in [0.0, 0.25, 0.5, 0.75, 1.0].into_iter().enumerate() {
        assert_abs_diff_eq!(values[index], expected, epsilon = TOL);
    }
}

/// Ports `SineSpacing`.
///
/// Sine clustering is one-sided: it packs points toward the *far* end, so the
/// first gap is the largest. That asymmetry is what distinguishes it from
/// cosine, whose gaps shrink at both ends.
#[test]
fn sine_spacing_spans_the_unit_interval_and_clusters_at_one_end() {
    let values = sine_spacing(11);
    assert_eq!(values.len(), 11);
    assert_abs_diff_eq!(values[0], 0.0, epsilon = TOL);
    assert_abs_diff_eq!(values[10], 1.0, epsilon = TOL);
    assert!(values[1] - values[0] > values[10] - values[9]);
}

/// Ports `BlendedSpacingUniform` and `BlendedSpacingCosine`.
#[test]
fn blended_spacing_interpolates_between_uniform_and_cosine() {
    let uniform = uniform_spacing(5);
    for (blended, expected) in blended_spacing(5, 0.0).into_iter().zip(&uniform) {
        assert_abs_diff_eq!(blended, expected, epsilon = TOL);
    }
    let cosine = cosine_spacing(5);
    for (blended, expected) in blended_spacing(5, 1.0).into_iter().zip(&cosine) {
        assert_abs_diff_eq!(blended, expected, epsilon = TOL);
    }
}

/// Ports the five `SpacingDispatcher*` cases.
///
/// `None` is the one worth stating explicitly: it means "caller has no
/// preference" and resolves to uniform, not to an empty or degenerate list.
#[test]
fn spacing_dispatcher_routes_each_distribution() {
    let uniform = spacing(5, SpacingType::Uniform);
    assert_eq!(uniform.len(), 5);
    assert_abs_diff_eq!(uniform[2], 0.5, epsilon = TOL);

    let cosine = spacing(5, SpacingType::Cosine);
    assert_eq!(cosine.len(), 5);
    assert_abs_diff_eq!(cosine[0], 0.0, epsilon = TOL);
    assert_abs_diff_eq!(cosine[4], 1.0, epsilon = TOL);

    assert_eq!(spacing(5, SpacingType::Sine).len(), 5);
    assert_eq!(spacing(5, SpacingType::Blend).len(), 5);

    let none = spacing(5, SpacingType::None);
    assert_eq!(none.len(), 5);
    assert_abs_diff_eq!(none[2], 0.5, epsilon = TOL);
}

// ============================================================================
// Rotations
// ============================================================================

/// Ports `RotateAboutXIdentity` and `RotateAboutX90`.
#[test]
fn rotation_about_x_is_identity_at_zero_and_maps_y_to_z_at_ninety() {
    let v = Vec3::new(1.0, 2.0, 3.0);
    assert_abs_diff_eq!(rotate_about_x(&v, 0.0), v, epsilon = TOL);
    let rotated = rotate_about_x(&Vec3::new(0.0, 1.0, 0.0), PI / 2.0);
    assert_abs_diff_eq!(rotated, Vec3::new(0.0, 0.0, 1.0), epsilon = TOL);
}

/// Ports `RotateAboutYIdentity` and `RotateAboutY90`.
#[test]
fn rotation_about_y_is_identity_at_zero_and_maps_x_to_negative_z_at_ninety() {
    let v = Vec3::new(1.0, 2.0, 3.0);
    assert_abs_diff_eq!(rotate_about_y(&v, 0.0), v, epsilon = TOL);
    let rotated = rotate_about_y(&Vec3::new(1.0, 0.0, 0.0), PI / 2.0);
    assert_abs_diff_eq!(rotated, Vec3::new(0.0, 0.0, -1.0), epsilon = TOL);
}

/// Ports `RotateAboutZIdentity` and `RotateAboutZ90`.
#[test]
fn rotation_about_z_is_identity_at_zero_and_maps_x_to_y_at_ninety() {
    let v = Vec3::new(1.0, 2.0, 3.0);
    assert_abs_diff_eq!(rotate_about_z(&v, 0.0), v, epsilon = TOL);
    let rotated = rotate_about_z(&Vec3::new(1.0, 0.0, 0.0), PI / 2.0);
    assert_abs_diff_eq!(rotated, Vec3::new(0.0, 1.0, 0.0), epsilon = TOL);
}

// ============================================================================
// Coordinate transformations
// ============================================================================

/// Ports `WindToBodyIdentity`, `WindToBodyOrthogonality` and `WindToBodyAlphaOnly`.
#[test]
fn wind_to_body_rotation_is_orthogonal_with_the_expected_alpha_row() {
    assert_abs_diff_eq!(
        wind_to_body_rotation(0.0, 0.0),
        Mat3::identity(),
        epsilon = TOL
    );

    let rotation = wind_to_body_rotation(0.1, 0.05);
    assert!((rotation.transpose() * rotation - Mat3::identity()).norm() < 1.0e-10);

    let alpha = PI / 6.0;
    let rotation = wind_to_body_rotation(alpha, 0.0);
    assert_abs_diff_eq!(rotation[(0, 0)], alpha.cos(), epsilon = TOL);
    assert_abs_diff_eq!(rotation[(0, 2)], -alpha.sin(), epsilon = TOL);
}

/// Ports `BodyToStabilityIdentity`, `BodyToStabilityOrthogonality` and
/// `BodyToStabilityElements`.
#[test]
fn body_to_stability_rotation_is_orthogonal_with_the_expected_elements() {
    assert_abs_diff_eq!(
        body_to_stability_rotation(0.0),
        Mat3::identity(),
        epsilon = TOL
    );

    let rotation = body_to_stability_rotation(0.2);
    assert!((rotation.transpose() * rotation - Mat3::identity()).norm() < 1.0e-10);

    let alpha = PI / 4.0;
    let rotation = body_to_stability_rotation(alpha);
    assert_abs_diff_eq!(rotation[(0, 0)], alpha.cos(), epsilon = TOL);
    assert_abs_diff_eq!(rotation[(0, 2)], alpha.sin(), epsilon = TOL);
    assert_abs_diff_eq!(rotation[(1, 1)], 1.0, epsilon = TOL);
    assert_abs_diff_eq!(rotation[(2, 0)], -alpha.sin(), epsilon = TOL);
    assert_abs_diff_eq!(rotation[(2, 2)], alpha.cos(), epsilon = TOL);
}
