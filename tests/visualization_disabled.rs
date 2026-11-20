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

#![cfg(not(feature = "vtk"))]

use zephyr_vlm::error::ZephyrError;
use zephyr_vlm::solver::{SolutionState, VortexLattice};
use zephyr_vlm::visualization::{VtkExportOptions, VtkExporter, is_available};

#[test]
fn vtk_availability_is_false_without_the_feature() {
    assert!(!is_available());
}

#[test]
fn visualization_api_returns_typed_feature_disabled_error() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("disabled.vtu");
    let error = VtkExporter::export_surface_mesh(
        &path,
        &VortexLattice::default(),
        &SolutionState::default(),
        VtkExportOptions::default(),
    )
    .unwrap_err();
    assert!(matches!(error, ZephyrError::FeatureDisabled("vtk")));
    assert!(!path.exists());
}
