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

//! Ported from the frozen v0.2.0 `body_test.cpp`.
//!
//! `Body` carries no aerodynamics of its own — it is geometry plus a transform
//! that the mesher and the VTK exporter both read. The defaults are the part
//! worth pinning: a unit scale and a zero translation mean an untouched body
//! passes through both unchanged.

use approx::assert_abs_diff_eq;
use zephyr_vlm::core::Vec3;
use zephyr_vlm::geometry::{Body, BodySection};

/// Ports `DefaultConstruction`, `DefaultScale` and `DefaultTranslation`.
#[test]
fn body_default_construction() {
    let body = Body::default();
    assert!(body.name().is_empty());
    assert_eq!(body.sections().len(), 0);
    assert_eq!(body.node_count(), 20);
    assert_abs_diff_eq!(body.node_spacing(), 1.0);
    assert!(!body.is_y_duplicated());
    // Identity transform: scale one, translation zero.
    assert_abs_diff_eq!(*body.scale(), Vec3::new(1.0, 1.0, 1.0));
    assert_abs_diff_eq!(*body.translation(), Vec3::zeros());
}

/// Ports `NamedConstruction` and `SetName`.
#[test]
fn body_name_construction_and_mutation() {
    let body = Body::new("fuselage");
    assert_eq!(body.name(), "fuselage");
    let mut body = Body::new("old_name");
    body.set_name("new_name");
    assert_eq!(body.name(), "new_name");
}

/// Ports `AddAndGetSections`.
#[test]
fn body_sections_are_appended_in_order() {
    let mut body = Body::new("fuselage");
    body.add_section(BodySection {
        center: Vec3::new(0.0, 0.0, 0.0),
        width: 0.5,
        height: 0.3,
    });
    body.add_section(BodySection {
        center: Vec3::new(1.0, 0.0, 0.0),
        width: 0.8,
        height: 0.6,
    });
    assert_eq!(body.sections().len(), 2);
    assert_abs_diff_eq!(body.sections()[0].width, 0.5);
    assert_abs_diff_eq!(body.sections()[1].width, 0.8);
}

/// Ports `NodeCount` and `NodeSpacing`.
#[test]
fn body_node_count_and_spacing_round_trip() {
    let mut body = Body::default();
    body.set_node_count(30);
    assert_eq!(body.node_count(), 30);
    body.set_node_spacing(2.5);
    assert_abs_diff_eq!(body.node_spacing(), 2.5);
}

/// Ports `Translation` and `Scale`.
#[test]
fn body_transform_round_trips() {
    let mut body = Body::default();
    body.set_translation(Vec3::new(1.0, 2.0, 3.0));
    assert_abs_diff_eq!(*body.translation(), Vec3::new(1.0, 2.0, 3.0));
    body.set_scale(Vec3::new(2.0, 3.0, 4.0));
    assert_abs_diff_eq!(*body.scale(), Vec3::new(2.0, 3.0, 4.0));
}

/// Ports `YDuplicate`.
#[test]
fn body_y_duplicate_toggles_both_ways() {
    let mut body = Body::default();
    assert!(!body.is_y_duplicated());
    body.set_y_duplicated(true);
    assert!(body.is_y_duplicated());
    body.set_y_duplicated(false);
    assert!(!body.is_y_duplicated());
}
