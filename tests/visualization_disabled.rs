#![cfg(not(feature = "vtk"))]

use pegasus_vlm::error::PegasusError;
use pegasus_vlm::solver::{SolutionState, VortexLattice};
use pegasus_vlm::visualization::{VtkExportOptions, VtkExporter, is_available};

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
    assert!(matches!(error, PegasusError::FeatureDisabled("vtk")));
    assert!(!path.exists());
}
