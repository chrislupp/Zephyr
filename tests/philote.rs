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

//! Tests for the Philote-MDO disciplines.
//!
//! The load-bearing cases difference each discipline's own `compute` and check
//! the declared partials against it. A discipline that reports a wrong Jacobian
//! is worse than one that reports none: the optimizer will follow it.
//!
//! Comparisons are relative with an absolute floor of [`NOISE_FLOOR`]. A
//! centered difference of step `h` on an O(1) quantity carries roundoff of
//! about `eps / h`, which is near `1e-9` at the step used here, so a tighter
//! floor would compare noise against an exact zero rather than checking
//! anything.

#![cfg(feature = "philote")]

use std::collections::HashMap;

use approx::assert_abs_diff_eq;
use ndarray::{ArrayD, IxDyn};
use philote_mdo::traits::{Discipline, ExplicitDiscipline};
use philote_mdo::{ArrayMap, PartialMap};
use zephyr_vlm::philote::{MeshDiscipline, SurfaceOption, VlmDiscipline};

/// Absolute tolerance floor: the roundoff of a centered difference at the step
/// these tests use, `eps / h` with `h = 1e-7`.
const NOISE_FLOOR: f64 = 1.0e-8;

fn array(values: &[f64]) -> ArrayD<f64> {
    ArrayD::from_shape_vec(IxDyn(&[values.len()]), values.to_vec()).unwrap()
}

fn options(pairs: &[(&str, serde_json::Value)]) -> HashMap<String, serde_json::Value> {
    pairs
        .iter()
        .map(|(key, value)| ((*key).to_string(), value.clone()))
        .collect()
}

/// Centered difference of every output with respect to one element of one
/// input, run through the discipline's own `compute`.
async fn difference<D: ExplicitDiscipline>(
    discipline: &D,
    inputs: &ArrayMap,
    name: &str,
    index: usize,
    step: f64,
) -> HashMap<String, Vec<f64>> {
    let mut plus = inputs.clone();
    let mut minus = inputs.clone();
    plus.get_mut(name).unwrap().as_slice_mut().unwrap()[index] += step;
    minus.get_mut(name).unwrap().as_slice_mut().unwrap()[index] -= step;

    let forward = discipline.compute(&plus).await.unwrap();
    let backward = discipline.compute(&minus).await.unwrap();
    forward
        .iter()
        .map(|(output, values)| {
            let other = &backward[output];
            let derivative = values
                .iter()
                .zip(other.iter())
                .map(|(a, b)| (a - b) / (2.0 * step))
                .collect();
            (output.clone(), derivative)
        })
        .collect()
}

/// Compare every declared partial against a difference, input element by input
/// element.
async fn assert_partials_match<D: ExplicitDiscipline>(
    discipline: &D,
    inputs: &ArrayMap,
    partials: &PartialMap,
    step: f64,
    tolerance: f64,
) {
    let mut checked = 0usize;
    for (name, values) in inputs {
        for index in 0..values.len() {
            let differences = difference(discipline, inputs, name, index, step).await;
            for (output, derivative) in &differences {
                let Some(block) = partials.get(&(output.clone(), name.clone())) else {
                    continue;
                };
                let stride = values.len();
                for (row, expected) in derivative.iter().enumerate() {
                    let exact = block.as_slice().unwrap()[row * stride + index];
                    assert!(
                        (exact - expected).abs() <= tolerance * expected.abs() + NOISE_FLOOR,
                        "d{output}/d{name}[{index}], row {row}: {exact} != {expected}"
                    );
                    checked += 1;
                }
            }
        }
    }
    assert!(checked > 0, "no partials were compared");
}

// ============================================================================
// MeshDiscipline
// ============================================================================

fn mesh_discipline(nx: usize, ny: usize, duplicated: bool) -> MeshDiscipline {
    let mut discipline = MeshDiscipline::default();
    discipline.initialize().unwrap();
    discipline
        .set_options(&options(&[
            ("nx", nx.into()),
            ("ny", ny.into()),
            ("y_duplicated", duplicated.into()),
        ]))
        .unwrap();
    discipline.setup().unwrap();
    discipline.setup_partials().unwrap();
    discipline
}

fn mesh_inputs(stations: usize) -> ArrayMap {
    let eta: Vec<f64> = (0..stations)
        .map(|i| i as f64 / (stations - 1) as f64)
        .collect();
    let chord: Vec<f64> = eta.iter().map(|e| 1.0 - 0.4 * e).collect();
    let twist: Vec<f64> = eta.iter().map(|e| 0.05 - 0.08 * e).collect();
    let xle: Vec<f64> = eta.iter().map(|e| 0.3 * e).collect();
    let zle: Vec<f64> = eta.iter().map(|e| 0.15 * e).collect();

    let mut inputs = HashMap::new();
    inputs.insert("chord".to_string(), array(&chord));
    inputs.insert("twist".to_string(), array(&twist));
    inputs.insert("xle".to_string(), array(&xle));
    inputs.insert("zle".to_string(), array(&zle));
    inputs.insert("eta".to_string(), array(&eta));
    inputs.insert("halfspan".to_string(), array(&[4.0]));
    inputs
}

#[tokio::test]
async fn mesh_discipline_declares_its_variables_with_the_right_shapes() {
    let discipline = mesh_discipline(4, 6, false);
    let names: Vec<&str> = discipline
        .registry()
        .var_meta()
        .iter()
        .map(|meta| meta.name.as_str())
        .collect();
    for expected in [
        "chord",
        "twist",
        "xle",
        "zle",
        "halfspan",
        "eta",
        "mesh",
        "reference_area",
        "reference_span",
        "reference_chord",
    ] {
        assert!(names.contains(&expected), "missing `{expected}`");
    }
    let mesh = discipline
        .registry()
        .var_meta()
        .iter()
        .find(|meta| meta.name == "mesh")
        .unwrap();
    // (nx + 1) * (ny + 1) * 3
    assert_eq!(mesh.shape, vec![5 * 7 * 3]);
}

#[tokio::test]
async fn mesh_discipline_places_nodes_and_reference_quantities() {
    let discipline = mesh_discipline(2, 4, false);
    let inputs = mesh_inputs(5);
    let outputs = discipline.compute(&inputs).await.unwrap();

    let mesh = outputs["mesh"].as_slice().unwrap();
    assert_eq!(mesh.len(), 3 * 5 * 3);

    // The root leading edge: chordwise fraction zero, so it sits exactly at
    // (xle, 0, zle) regardless of chord or twist.
    assert_abs_diff_eq!(mesh[0], 0.0, epsilon = 1.0e-12);
    assert_abs_diff_eq!(mesh[1], 0.0, epsilon = 1.0e-12);
    assert_abs_diff_eq!(mesh[2], 0.0, epsilon = 1.0e-12);

    // A rectangular untwisted planform has a known area and MAC, so check
    // against that rather than against the code's own output.
    let mut rectangular = mesh_inputs(5);
    rectangular.insert("chord".to_string(), array(&[1.0; 5]));
    rectangular.insert("twist".to_string(), array(&[0.0; 5]));
    rectangular.insert("zle".to_string(), array(&[0.0; 5]));
    let outputs = discipline.compute(&rectangular).await.unwrap();
    assert_abs_diff_eq!(outputs["reference_area"][0], 4.0, epsilon = 1.0e-12);
    assert_abs_diff_eq!(outputs["reference_span"][0], 4.0, epsilon = 1.0e-12);
    assert_abs_diff_eq!(outputs["reference_chord"][0], 1.0, epsilon = 1.0e-12);
}

#[tokio::test]
async fn mesh_discipline_doubles_area_and_span_when_duplicated() {
    let inputs = {
        let mut inputs = mesh_inputs(5);
        inputs.insert("chord".to_string(), array(&[1.0; 5]));
        inputs.insert("twist".to_string(), array(&[0.0; 5]));
        inputs.insert("zle".to_string(), array(&[0.0; 5]));
        inputs
    };
    let half = mesh_discipline(2, 4, false).compute(&inputs).await.unwrap();
    let full = mesh_discipline(2, 4, true).compute(&inputs).await.unwrap();

    assert_abs_diff_eq!(
        full["reference_area"][0],
        2.0 * half["reference_area"][0],
        epsilon = 1.0e-12
    );
    assert_abs_diff_eq!(
        full["reference_span"][0],
        2.0 * half["reference_span"][0],
        epsilon = 1.0e-12
    );
    // The mean aerodynamic chord is an average, so mirroring leaves it alone.
    assert_abs_diff_eq!(
        full["reference_chord"][0],
        half["reference_chord"][0],
        epsilon = 1.0e-12
    );
}

#[tokio::test]
async fn mesh_discipline_partials_match_finite_differences() {
    let discipline = mesh_discipline(2, 4, true);
    let inputs = mesh_inputs(5);
    let partials = discipline.compute_partials(&inputs).await.unwrap();
    assert_partials_match(&discipline, &inputs, &partials, 1.0e-7, 1.0e-4).await;
}

// ============================================================================
// VlmDiscipline
// ============================================================================

fn vlm_discipline(nx: usize, ny: usize) -> VlmDiscipline {
    let mut discipline = VlmDiscipline::default();
    discipline.initialize().unwrap();
    discipline
        .set_options(&options(&[
            ("num_surfaces", 1.into()),
            ("surface_names", serde_json::json!(["Wing"])),
            ("surface_nx", serde_json::json!([nx])),
            ("surface_ny", serde_json::json!([ny])),
        ]))
        .unwrap();
    discipline.setup().unwrap();
    discipline.setup_partials().unwrap();
    discipline
}

/// A flat rectangular wing, plus a five-degree flight condition.
fn vlm_inputs(nx: usize, ny: usize, chord: f64, span: f64) -> ArrayMap {
    let mut nodes = Vec::with_capacity((nx + 1) * (ny + 1) * 3);
    for i in 0..=nx {
        for j in 0..=ny {
            nodes.push(chord * i as f64 / nx as f64);
            nodes.push(span * j as f64 / ny as f64);
            nodes.push(0.0);
        }
    }
    let mut inputs = HashMap::new();
    inputs.insert("Wing".to_string(), array(&nodes));
    inputs.insert("reference_area".to_string(), array(&[chord * span]));
    inputs.insert("reference_chord".to_string(), array(&[chord]));
    inputs.insert("reference_span".to_string(), array(&[span]));
    inputs.insert("moment_ref_x".to_string(), array(&[0.25 * chord]));
    inputs.insert("moment_ref_y".to_string(), array(&[0.0]));
    inputs.insert("moment_ref_z".to_string(), array(&[0.0]));
    inputs.insert("alpha".to_string(), array(&[5.0_f64.to_radians()]));
    inputs.insert("beta".to_string(), array(&[0.0]));
    inputs.insert("mach".to_string(), array(&[0.0]));
    inputs.insert("dynamic_pressure".to_string(), array(&[1.0]));
    inputs.insert("roll_rate".to_string(), array(&[0.0]));
    inputs.insert("pitch_rate".to_string(), array(&[0.0]));
    inputs.insert("yaw_rate".to_string(), array(&[0.0]));
    inputs
}

#[tokio::test]
async fn vlm_discipline_declares_the_six_coefficient_outputs() {
    let discipline = vlm_discipline(2, 4);
    let outputs: Vec<&str> = discipline
        .registry()
        .var_meta()
        .iter()
        .filter(|meta| ["CL", "CD", "CY", "Cl", "Cm", "Cn"].contains(&meta.name.as_str()))
        .map(|meta| meta.name.as_str())
        .collect();
    assert_eq!(outputs.len(), 6);
    assert_eq!(discipline.surfaces().len(), 1);
    assert_eq!(discipline.surfaces()[0].name, "Wing");
}

#[tokio::test]
async fn vlm_discipline_produces_positive_lift_at_positive_alpha() {
    let discipline = vlm_discipline(2, 4);
    let outputs = discipline
        .compute(&vlm_inputs(2, 4, 1.0, 4.0))
        .await
        .unwrap();
    assert!(outputs["CL"][0] > 0.0);
    assert!(outputs["CD"][0] > 0.0);
    // A symmetric planform at zero sideslip has no lateral loads.
    assert_abs_diff_eq!(outputs["CY"][0], 0.0, epsilon = 1.0e-10);
}

#[tokio::test]
async fn vlm_discipline_node_and_flight_partials_match_finite_differences() {
    let discipline = vlm_discipline(2, 3);
    let inputs = vlm_inputs(2, 3, 1.0, 3.0);
    let partials = discipline.compute_partials(&inputs).await.unwrap();

    // Node coordinates and the two flight variables the fixture actually
    // responds to. Rates and Mach are zero here, so differencing them would
    // compare noise against noise.
    for name in ["Wing", "alpha", "beta"] {
        let length = inputs[name].len();
        for index in 0..length {
            let differences = difference(&discipline, &inputs, name, index, 1.0e-7).await;
            for (output, derivative) in &differences {
                let block = &partials[&(output.clone(), name.to_string())];
                let exact = block.as_slice().unwrap()[index];
                assert_abs_diff_eq!(
                    exact,
                    derivative[0],
                    epsilon = 1.0e-4 * derivative[0].abs() + NOISE_FLOOR
                );
            }
        }
    }
}

#[tokio::test]
async fn vlm_discipline_reference_area_partials_are_the_normalisation_derivative() {
    let discipline = vlm_discipline(2, 3);
    let inputs = vlm_inputs(2, 3, 1.0, 3.0);
    let outputs = discipline.compute(&inputs).await.unwrap();
    let partials = discipline.compute_partials(&inputs).await.unwrap();
    let area = inputs["reference_area"][0];

    // The coefficients are loads divided by the reference area, so the
    // derivative is exactly minus the coefficient over that area.
    for name in ["CL", "CD", "CY", "Cl", "Cm", "Cn"] {
        let exact = partials[&(name.to_string(), "reference_area".to_string())][0];
        assert_abs_diff_eq!(exact, -outputs[name][0] / area, epsilon = 1.0e-12);
    }
}

#[tokio::test]
async fn vlm_discipline_rejects_a_wrongly_sized_node_array() {
    let discipline = vlm_discipline(2, 3);
    let mut inputs = vlm_inputs(2, 3, 1.0, 3.0);
    inputs.insert("Wing".to_string(), array(&[0.0; 9]));
    assert!(discipline.compute(&inputs).await.is_err());
}

// ============================================================================
// Options
// ============================================================================

#[tokio::test]
async fn surface_options_accept_arrays_and_comma_separated_strings() {
    let configure = |values: Vec<(&str, serde_json::Value)>| {
        let mut discipline = VlmDiscipline::default();
        discipline.initialize().unwrap();
        discipline.set_options(&options(&values)).unwrap();
        discipline
    };

    let from_arrays = configure(vec![
        ("num_surfaces", 2.into()),
        ("surface_names", serde_json::json!(["Wing", "Tail"])),
        ("surface_nx", serde_json::json!([4, 2])),
        ("surface_ny", serde_json::json!([8, 3])),
        ("surface_y_duplicated", serde_json::json!([true, false])),
    ]);
    // The spelling the v0.2.0 C++ discipline used.
    let from_strings = configure(vec![
        ("num_surfaces", 2.into()),
        ("surface_names", "Wing,Tail".into()),
        ("surface_nx", "4,2".into()),
        ("surface_ny", "8,3".into()),
        ("surface_y_duplicated", "true,false".into()),
    ]);

    for discipline in [&from_arrays, &from_strings] {
        let surfaces = discipline.surfaces();
        assert_eq!(surfaces.len(), 2);
        assert_eq!(surfaces[0].name, "Wing");
        assert_eq!(surfaces[0].nx, 4);
        assert_eq!(surfaces[0].ny, 8);
        assert!(surfaces[0].y_duplicated);
        assert_eq!(surfaces[1].name, "Tail");
        assert_eq!(surfaces[1].nx, 2);
        assert_eq!(surfaces[1].ny, 3);
        assert!(!surfaces[1].y_duplicated);
    }
}

/// The constructor path, for a caller building a discipline in process rather
/// than configuring one over the wire.
#[tokio::test]
async fn a_discipline_can_be_built_directly_from_surface_options() {
    let mut discipline = VlmDiscipline::new(vec![SurfaceOption {
        name: "Panel".into(),
        nx: 2,
        ny: 3,
        ..SurfaceOption::default()
    }]);
    discipline.setup().unwrap();
    discipline.setup_partials().unwrap();

    let names: Vec<&str> = discipline
        .registry()
        .var_meta()
        .iter()
        .map(|meta| meta.name.as_str())
        .collect();
    assert!(names.contains(&"Panel"));
    // (nx + 1) * (ny + 1) * 3
    let panel = discipline
        .registry()
        .var_meta()
        .iter()
        .find(|meta| meta.name == "Panel")
        .unwrap();
    assert_eq!(panel.shape, vec![3 * 4 * 3]);
}

#[tokio::test]
async fn a_surface_name_count_mismatch_is_rejected() {
    let mut discipline = VlmDiscipline::default();
    discipline.initialize().unwrap();
    let error = discipline
        .set_options(&options(&[
            ("num_surfaces", 2.into()),
            ("surface_names", serde_json::json!(["Wing"])),
        ]))
        .unwrap_err();
    // Names are the input variable names, so a mismatch would produce a
    // discipline whose inputs a client cannot address.
    assert!(error.to_string().contains("surface_names"));
}

// ============================================================================
// The two disciplines chained
// ============================================================================

#[tokio::test]
async fn the_mesh_discipline_output_drives_the_vlm_discipline() {
    let (nx, ny) = (2, 4);
    let mesh = mesh_discipline(nx, ny, false);
    let mut geometry = mesh_inputs(ny + 1);
    geometry.insert("chord".to_string(), array(&[1.0; 5]));
    geometry.insert("twist".to_string(), array(&[0.0; 5]));
    geometry.insert("zle".to_string(), array(&[0.0; 5]));
    geometry.insert("xle".to_string(), array(&[0.0; 5]));
    let generated = mesh.compute(&geometry).await.unwrap();

    let vlm = vlm_discipline(nx, ny);
    let mut inputs = vlm_inputs(nx, ny, 1.0, 4.0);
    inputs.insert("Wing".to_string(), generated["mesh"].clone());
    inputs.insert(
        "reference_area".to_string(),
        generated["reference_area"].clone(),
    );
    inputs.insert(
        "reference_chord".to_string(),
        generated["reference_chord"].clone(),
    );
    inputs.insert(
        "reference_span".to_string(),
        generated["reference_span"].clone(),
    );

    let outputs = vlm.compute(&inputs).await.unwrap();
    assert!(outputs["CL"][0] > 0.0);
    assert!(outputs["CL"][0].is_finite());
}

// ============================================================================
// The three Philote cases the frozen manifest deferred
//
// `ctest-manifest.json` records 491 C++ cases, of which 3 were excluded as
// `deferred_philote_cases` because the disciplines did not exist in the Rust
// tree. They do now, so those three are ported here and the C++ suite is
// accounted for end to end.
// ============================================================================

/// Ports `MeshDisciplineRefTest.RectangularWingArea`.
///
/// A duplicated unit-chord wing of half span 5: area 10, span 10, MAC 1. The
/// numbers are closed-form, so this checks the reference integrals rather than
/// the code's own output.
#[tokio::test]
async fn deferred_rectangular_wing_reference_quantities() {
    let (nx, ny) = (2, 4);
    let stations = ny + 1;
    let discipline = mesh_discipline(nx, ny, true);

    let eta: Vec<f64> = (0..stations).map(|j| j as f64 / ny as f64).collect();
    let mut inputs = HashMap::new();
    inputs.insert("chord".to_string(), array(&vec![1.0; stations]));
    inputs.insert("twist".to_string(), array(&vec![0.0; stations]));
    inputs.insert("xle".to_string(), array(&vec![0.0; stations]));
    inputs.insert("zle".to_string(), array(&vec![0.0; stations]));
    inputs.insert("eta".to_string(), array(&eta));
    inputs.insert("halfspan".to_string(), array(&[5.0]));

    let outputs = discipline.compute(&inputs).await.unwrap();
    assert_abs_diff_eq!(outputs["reference_area"][0], 10.0, epsilon = 1.0e-12);
    assert_abs_diff_eq!(outputs["reference_span"][0], 10.0, epsilon = 1.0e-12);
    assert_abs_diff_eq!(outputs["reference_chord"][0], 1.0, epsilon = 1.0e-12);
}

/// Ports `MeshDisciplineRefTest.TaperedWingArea`.
///
/// Chord runs 2 to 1 over a half span of 4, duplicated: area 12, span 8. The
/// MAC is recomputed here from the chord-weighted definition rather than
/// hard-coded, so it checks the weighting and not just a number.
#[tokio::test]
async fn deferred_tapered_wing_reference_quantities() {
    let (nx, ny) = (1, 4);
    let stations = ny + 1;
    let half_span = 4.0;
    let discipline = mesh_discipline(nx, ny, true);

    let eta: Vec<f64> = (0..stations).map(|j| j as f64 / ny as f64).collect();
    let chord: Vec<f64> = eta.iter().map(|e| 2.0 - e).collect();
    let mut inputs = HashMap::new();
    inputs.insert("chord".to_string(), array(&chord));
    inputs.insert("twist".to_string(), array(&vec![0.0; stations]));
    inputs.insert("xle".to_string(), array(&vec![0.0; stations]));
    inputs.insert("zle".to_string(), array(&vec![0.0; stations]));
    inputs.insert("eta".to_string(), array(&eta));
    inputs.insert("halfspan".to_string(), array(&[half_span]));

    let outputs = discipline.compute(&inputs).await.unwrap();
    assert_abs_diff_eq!(outputs["reference_area"][0], 12.0, epsilon = 1.0e-12);
    assert_abs_diff_eq!(outputs["reference_span"][0], 8.0, epsilon = 1.0e-12);

    let segment = half_span / ny as f64;
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    for interval in 0..ny {
        let (left, right) = (chord[interval], chord[interval + 1]);
        numerator += 0.5 * (left * left + right * right) * segment;
        denominator += 0.5 * (left + right) * segment;
    }
    assert_abs_diff_eq!(
        outputs["reference_chord"][0],
        numerator / denominator,
        epsilon = 1.0e-12
    );
}

/// Ports `MeshDisciplineRefTest.PartialsFiniteDifference`.
///
/// The reference-quantity partials against differences of `compute`, on the
/// tapered fixture the frozen case used. `mesh_discipline_partials_match_finite_differences`
/// covers every partial on a different geometry; this pins the same check to
/// the frozen one.
#[tokio::test]
async fn deferred_reference_quantity_partials_match_finite_differences() {
    let (nx, ny) = (1, 4);
    let stations = ny + 1;
    let discipline = mesh_discipline(nx, ny, true);

    let eta: Vec<f64> = (0..stations).map(|j| j as f64 / ny as f64).collect();
    let chord: Vec<f64> = eta.iter().map(|e| 2.0 - e).collect();
    // A nonzero dihedral so the `zle` partials are not trivially zero.
    let zle: Vec<f64> = eta.iter().map(|e| 0.2 * e).collect();
    let mut inputs = HashMap::new();
    inputs.insert("chord".to_string(), array(&chord));
    inputs.insert("twist".to_string(), array(&vec![0.0; stations]));
    inputs.insert("xle".to_string(), array(&vec![0.0; stations]));
    inputs.insert("zle".to_string(), array(&zle));
    inputs.insert("eta".to_string(), array(&eta));
    inputs.insert("halfspan".to_string(), array(&[4.0]));

    let partials = discipline.compute_partials(&inputs).await.unwrap();
    let mut checked = 0usize;
    for name in ["chord", "zle", "halfspan", "eta"] {
        for index in 0..inputs[name].len() {
            let differences = difference(&discipline, &inputs, name, index, 1.0e-7).await;
            for output in ["reference_area", "reference_chord"] {
                let block = &partials[&(output.to_string(), name.to_string())];
                let expected = differences[output][0];
                assert_abs_diff_eq!(
                    block.as_slice().unwrap()[index],
                    expected,
                    epsilon = 1.0e-5 * expected.abs() + NOISE_FLOOR
                );
                checked += 1;
            }
        }
    }
    assert!(checked > 0);
}

/// Every entry of the VLM Jacobian, against differences of `compute`.
///
/// `vlm_discipline_node_and_flight_partials_match_finite_differences` covers the
/// node, alpha and beta columns. This covers all thirteen scalar inputs at a
/// condition where nothing is zero by symmetry — sideslip, compressibility, all
/// three body rates, and an offset reference point — so no column can be
/// silently zero and still pass.
#[tokio::test]
async fn vlm_discipline_every_scalar_partial_matches_finite_differences() {
    let discipline = vlm_discipline(3, 4);
    let mut inputs = vlm_inputs(3, 4, 1.0, 4.0);
    inputs.insert("beta".to_string(), array(&[3.0_f64.to_radians()]));
    inputs.insert("mach".to_string(), array(&[0.3]));
    inputs.insert("roll_rate".to_string(), array(&[0.02]));
    inputs.insert("pitch_rate".to_string(), array(&[0.03]));
    inputs.insert("yaw_rate".to_string(), array(&[0.04]));
    inputs.insert("moment_ref_z".to_string(), array(&[0.1]));

    let partials = discipline.compute_partials(&inputs).await.unwrap();
    for name in [
        "reference_area",
        "reference_chord",
        "reference_span",
        "moment_ref_x",
        "moment_ref_y",
        "moment_ref_z",
        "alpha",
        "beta",
        "mach",
        "dynamic_pressure",
        "roll_rate",
        "pitch_rate",
        "yaw_rate",
    ] {
        let differences = difference(&discipline, &inputs, name, 0, 1.0e-6).await;
        for output in ["CL", "CD", "CY", "Cl", "Cm", "Cn"] {
            let exact = partials[&(output.to_string(), name.to_string())]
                .as_slice()
                .unwrap()[0];
            let expected = differences[output][0];
            assert!(
                (exact - expected).abs() <= 1.0e-5 * expected.abs() + NOISE_FLOOR,
                "{output}/{name}: analytic {exact}, finite {expected}"
            );
        }
    }
}

/// A partial that is absent is a partial an optimizer reads as zero, so an
/// undeclared pair whose difference is nonzero is as wrong as a declared pair
/// carrying the wrong number. This sweeps both disciplines for that case.
#[tokio::test]
async fn disciplines_declare_every_partial_that_is_not_identically_zero() {
    let vlm = vlm_discipline(3, 4);
    let mut inputs = vlm_inputs(3, 4, 1.0, 4.0);
    inputs.insert("beta".to_string(), array(&[3.0_f64.to_radians()]));
    inputs.insert("mach".to_string(), array(&[0.3]));
    inputs.insert("roll_rate".to_string(), array(&[0.02]));
    inputs.insert("pitch_rate".to_string(), array(&[0.03]));
    inputs.insert("yaw_rate".to_string(), array(&[0.04]));
    let declared = vlm.compute_partials(&inputs).await.unwrap();
    let names: Vec<String> = inputs.keys().cloned().collect();
    for name in &names {
        // The node array is long; a few entries establish the column exists.
        for index in 0..inputs[name].len().min(4) {
            let differences = difference(&vlm, &inputs, name, index, 1.0e-6).await;
            for output in ["CL", "CD", "CY", "Cl", "Cm", "Cn"] {
                let key = (output.to_string(), name.clone());
                assert!(
                    declared.contains_key(&key) || differences[output][0].abs() <= 1.0e-9,
                    "d{output}/d{name}[{index}] is {:.3e} but is not declared",
                    differences[output][0]
                );
            }
        }
    }

    let mesh = mesh_discipline(2, 4, true);
    let inputs = mesh_inputs(5);
    let declared = mesh.compute_partials(&inputs).await.unwrap();
    let names: Vec<String> = inputs.keys().cloned().collect();
    for name in &names {
        for index in 0..inputs[name].len() {
            let differences = difference(&mesh, &inputs, name, index, 1.0e-7).await;
            for (output, values) in &differences {
                let peak = values.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
                let key = (output.clone(), name.clone());
                assert!(
                    declared.contains_key(&key) || peak <= 1.0e-8,
                    "d{output}/d{name}[{index}] peaks at {peak:.3e} but is not declared"
                );
            }
        }
    }
}

#[tokio::test]
async fn selected_analysis_declares_and_computes_only_requested_outputs_and_partials() {
    let mut discipline = VlmDiscipline::default();
    discipline.initialize().unwrap();
    discipline
        .set_options(&options(&[
            ("surface_names", serde_json::json!(["Wing"])),
            ("surface_nx", serde_json::json!([2])),
            ("surface_ny", serde_json::json!([3])),
            ("surface_y_duplicated", serde_json::json!([true])),
            ("analysis_mode", "half_wing".into()),
            ("analysis_outputs", serde_json::json!(["CL", "CDi"])),
            (
                "analysis_gradients",
                serde_json::json!({"CL":["alpha"],"CDi":["Wing"]}),
            ),
        ]))
        .unwrap();
    discipline.setup().unwrap();
    discipline.setup_partials().unwrap();
    let inputs = vlm_inputs(2, 3, 1.0, 3.0);
    let outputs = discipline.compute(&inputs).await.unwrap();
    assert_eq!(outputs.len(), 2);
    assert!(outputs.contains_key("CL"));
    assert!(outputs.contains_key("CDi"));
    let partials = discipline.compute_partials(&inputs).await.unwrap();
    assert_eq!(partials.len(), 2);
    assert!(partials.contains_key(&("CL".into(), "alpha".into())));
    assert!(partials.contains_key(&("CDi".into(), "Wing".into())));
    let alpha_difference = difference(&discipline, &inputs, "alpha", 0, 1e-6).await;
    assert_abs_diff_eq!(
        partials[&("CL".into(), "alpha".into())][0],
        alpha_difference["CL"][0],
        epsilon = 1e-7
    );
    // An interior node preserves the half-space restriction in both directions.
    let node_index = 4;
    let derivative = difference(&discipline, &inputs, "Wing", node_index, 1e-6).await;
    assert_abs_diff_eq!(
        partials[&("CDi".into(), "Wing".into())][node_index],
        derivative["CDi"][0],
        epsilon = 1e-7
    );
}

#[tokio::test]
async fn explicit_empty_analysis_and_outputs_only_do_not_return_partials() {
    for selected in [serde_json::json!([]), serde_json::json!(["CL"])] {
        let mut discipline = VlmDiscipline::default();
        discipline.initialize().unwrap();
        discipline
            .set_options(&options(&[
                ("surface_names", serde_json::json!(["Wing"])),
                ("surface_nx", serde_json::json!([2])),
                ("surface_ny", serde_json::json!([3])),
                ("analysis_outputs", selected.clone()),
            ]))
            .unwrap();
        discipline.setup().unwrap();
        discipline.setup_partials().unwrap();
        let inputs = vlm_inputs(2, 3, 1.0, 3.0);
        assert_eq!(
            discipline.compute(&inputs).await.unwrap().len(),
            selected.as_array().unwrap().len()
        );
        assert!(
            discipline
                .compute_partials(&inputs)
                .await
                .unwrap()
                .is_empty()
        );
    }
}

fn selected_vlm() -> VlmDiscipline {
    let mut discipline = VlmDiscipline::default();
    discipline.initialize().unwrap();
    discipline.set_options(&options(&[
        ("surface_names", "Wing".into()),
        ("surface_nx", "2".into()), ("surface_ny", "3".into()),
        ("consistent_forces", true.into()),
        ("outputs", "CL,CD,Cm".into()),
        ("fixed_inputs", "beta,mach,roll_rate,pitch_rate,yaw_rate,moment_ref_x,moment_ref_y,moment_ref_z,reference_area,reference_chord,reference_span,dynamic_pressure".into()),
    ])).unwrap();
    discipline.setup().unwrap();
    discipline.setup_partials().unwrap();
    discipline
}

#[tokio::test]
async fn selected_vlm_partials_match_full_derivatives_and_independent_oracle() {
    let selected = selected_vlm();
    let mut full = vlm_discipline(2, 3);
    full.set_options(&options(&[
        ("consistent_forces", true.into()),
        ("surface_names", "Wing".into()),
        ("surface_nx", "2".into()),
        ("surface_ny", "3".into()),
    ]))
    .unwrap();
    let inputs = vlm_inputs(2, 3, 1.0, 3.0);
    let forward = selected.compute(&inputs).await.unwrap();
    assert_eq!(forward.len(), 3);
    let partials = selected.compute_partials(&inputs).await.unwrap();
    let reference = full.compute_partials(&inputs).await.unwrap();
    assert_eq!(partials.len(), 6); // Three outputs, mesh and alpha only.
    for (key, values) in &partials {
        for (a, b) in values.iter().zip(reference[key].iter()) {
            assert_abs_diff_eq!(a, b, epsilon = 1e-12);
        }
    }
    for input in ["Wing", "alpha"] {
        for (index, _) in inputs[input].iter().enumerate() {
            let differences = difference(&selected, &inputs, input, index, 1e-6).await;
            for (output, values) in differences {
                let analytic = partials[&(output, input.into())][index];
                assert_abs_diff_eq!(analytic, values[0], epsilon = 1e-7);
            }
        }
    }
}

#[tokio::test]
async fn vlm_cache_tracks_all_inputs_and_option_changes() {
    let mut cached = selected_vlm();
    let inputs = vlm_inputs(2, 3, 1.0, 3.0);
    let original = cached.compute(&inputs).await.unwrap();
    assert_eq!(original, cached.compute(&inputs).await.unwrap());
    for name in ["Wing", "alpha", "mach", "reference_area", "moment_ref_x"] {
        let mut changed = inputs.clone();
        changed.get_mut(name).unwrap()[0] += 0.01;
        let expected = selected_vlm().compute(&changed).await.unwrap();
        assert_eq!(expected, cached.compute(&changed).await.unwrap());
        assert_eq!(original, cached.compute(&inputs).await.unwrap());
    }
    // The legacy v0.2 force convention was removed: asking for it is an error.
    let error = cached
        .set_options(&options(&[
            ("consistent_forces", false.into()),
            ("surface_names", "Wing".into()),
            ("surface_nx", "2".into()),
            ("surface_ny", "3".into()),
        ]))
        .unwrap_err();
    assert!(
        error
            .to_string()
            .contains("legacy v0.2 force convention was removed")
    );
    // Explicit `true` is the default: after the option change clears the cache,
    // the result matches a discipline that never set it.
    cached
        .set_options(&options(&[
            ("consistent_forces", true.into()),
            ("surface_names", "Wing".into()),
            ("surface_nx", "2".into()),
            ("surface_ny", "3".into()),
        ]))
        .unwrap();
    let reference = vlm_discipline(2, 3).compute(&inputs).await.unwrap();
    let changed = cached.compute(&inputs).await.unwrap();
    for (name, value) in changed {
        assert_eq!(value, reference[&name]);
    }
}

#[test]
fn disciplines_report_zephyr_names() {
    assert_eq!(MeshDiscipline::default().name(), "ZephyrMesh");
    assert_eq!(VlmDiscipline::default().name(), "ZephyrVLM");
}

#[test]
fn vlm_rejects_unknown_output_and_fixed_input_names() {
    let mut discipline = VlmDiscipline::default();
    discipline.initialize().unwrap();
    assert!(
        discipline
            .set_options(&options(&[("outputs", "CD,typo".into())]))
            .is_err()
    );
    discipline
        .set_options(&options(&[("fixed_inputs", "typo".into())]))
        .unwrap();
    assert!(discipline.setup().is_err());
}

#[test]
fn selective_analysis_rejects_options_it_would_otherwise_ignore() {
    let control = serde_json::json!([
        {"surface": "Wing", "name": "flap", "span_start": 0, "span_end": 3, "hinge": 1.2}
    ]);
    let trim = serde_json::json!({
        "enabled": true, "variables": ["alpha"],
        "balances": [{"quantity": "Fz", "scale": 1000.0}]
    });
    let conflicting: [(&str, serde_json::Value); 4] = [
        ("controls", control),
        ("fixed_inputs", "mach".into()),
        ("outputs", "CL,CD".into()),
        ("trim", trim),
    ];
    for (name, value) in conflicting {
        let mut discipline = VlmDiscipline::default();
        discipline.initialize().unwrap();
        discipline
            .set_options(&options(&[
                ("surface_names", serde_json::json!(["Wing"])),
                ("surface_nx", serde_json::json!([2])),
                ("surface_ny", serde_json::json!([3])),
                ("surface_y_duplicated", serde_json::json!([true])),
                ("analysis_outputs", serde_json::json!(["CL"])),
                (name, value),
            ]))
            .unwrap_or_else(|e| panic!("{name} should be accepted as an option: {e}"));
        let error = discipline
            .setup()
            .expect_err(&format!("{name} must be rejected with analysis options"));
        assert!(
            error.to_string().contains("cannot be combined"),
            "{name}: {error}"
        );
    }
}

#[tokio::test]
async fn mesh_controls_have_exact_philote_partials_and_invalidate_cache() {
    let mut discipline = VlmDiscipline::default();
    discipline.initialize().unwrap();
    let controls = serde_json::json!([
        {"surface": "Wing", "name": "elevator", "span_start": 0, "span_end": 3,
         "hinge": [0.6, 0.8, 0.9, 1.1]},
        {"surface": "Wing", "name": "aileron", "span_start": 1, "span_end": 3,
         "hinge": 1.2, "gain": -0.7, "duplicate_sign": -1.0}
    ]);
    let configuration = options(&[
        ("surface_names", "Wing".into()),
        ("surface_nx", "2".into()),
        ("surface_ny", "3".into()),
        ("surface_y_duplicated", "true".into()),
        ("consistent_forces", true.into()),
        ("controls", controls.to_string().into()),
    ]);
    discipline.set_options(&configuration).unwrap();
    discipline.setup().unwrap();
    discipline.setup_partials().unwrap();
    let mut inputs = vlm_inputs(2, 3, 1.0, 3.0);
    inputs.insert("mach".into(), array(&[0.3]));
    inputs.insert("pitch_rate".into(), array(&[0.005]));
    inputs.insert("control_elevator".into(), array(&[0.13]));
    inputs.insert("control_aileron".into(), array(&[-0.05]));
    let original = discipline.compute(&inputs).await.unwrap();
    let partials = discipline.compute_partials(&inputs).await.unwrap();
    assert!(partials.contains_key(&("CL".into(), "control_elevator".into())));
    assert_partials_match(&discipline, &inputs, &partials, 1e-7, 2e-5).await;
    inputs.get_mut("control_elevator").unwrap()[0] += 0.04;
    assert_ne!(original, discipline.compute(&inputs).await.unwrap());
    inputs.get_mut("control_elevator").unwrap()[0] = 0.13;
    assert_eq!(original, discipline.compute(&inputs).await.unwrap());

    let mut fixed = VlmDiscipline::default();
    fixed.initialize().unwrap();
    let mut configuration = configuration.clone();
    configuration.insert(
        "fixed_inputs".into(),
        "control_elevator,control_aileron".into(),
    );
    fixed.set_options(&configuration).unwrap();
    fixed.setup().unwrap();
    fixed.setup_partials().unwrap();
    let partials = fixed.compute_partials(&inputs).await.unwrap();
    assert!(
        partials
            .keys()
            .all(|(_, input)| !input.starts_with("control_"))
    );
    assert_eq!(original, fixed.compute(&inputs).await.unwrap());
}

#[test]
fn philote_rejects_invalid_control_metadata() {
    for definition in [
        serde_json::json!({"surface":"Wing", "name":"e", "span_start":0, "span_end":1, "hinge":0.5, "typo":1}),
        serde_json::json!({"surface":"missing", "name":"e", "span_start":0, "span_end":1, "hinge":0.5}),
        serde_json::json!({"surface":"Wing", "name":"e", "span_start":0, "span_end":2, "hinge":0.5}),
        serde_json::json!({"surface":"Wing", "name":"e", "span_start":0, "span_end":1, "hinge":2.0}),
    ] {
        let mut discipline = VlmDiscipline::default();
        discipline.initialize().unwrap();
        assert!(
            discipline
                .set_options(&options(&[
                    ("surface_names", "Wing".into()),
                    ("controls", serde_json::json!([definition]))
                ]))
                .is_err()
        );
    }
}

// ============================================================================
// Opt-in solver accelerations: `mirror_reuse` and `symmetric_factorization`
// ============================================================================

/// Relative bound on the symmetric factorization's departure from the full LU.
const SYMMETRIC_BOUND: f64 = 1.0e-12;

/// Every entry of two output or partial maps has the same bits.
fn assert_bits<K: std::fmt::Debug + Eq + std::hash::Hash>(
    label: &str,
    actual: &HashMap<K, ArrayD<f64>>,
    expected: &HashMap<K, ArrayD<f64>>,
) {
    assert_eq!(actual.len(), expected.len(), "{label}: key count");
    for (key, values) in expected {
        let other = actual
            .get(key)
            .unwrap_or_else(|| panic!("{label}: missing {key:?}"));
        assert_eq!(values.len(), other.len(), "{label}: {key:?} length");
        for (a, b) in other.iter().zip(values.iter()) {
            assert_eq!(a.to_bits(), b.to_bits(), "{label}: {key:?}: {a} != {b}");
        }
    }
}

/// Every entry agrees to [`SYMMETRIC_BOUND`] relative to the largest magnitude
/// in its group, and returns the largest such relative deviation.
///
/// `group` maps a key to its group. Outputs are grouped by themselves; partials
/// by output, i.e. by Jacobian row, because a derivative that vanishes by
/// symmetry (`dCm/droll_rate` on a symmetric wing) is a cancellation of terms
/// as large as the rest of its row, and its rounding is of their size, not of
/// its own.
fn assert_near<K: std::fmt::Debug + Eq + std::hash::Hash>(
    label: &str,
    actual: &HashMap<K, ArrayD<f64>>,
    expected: &HashMap<K, ArrayD<f64>>,
    group: impl Fn(&K) -> String,
) -> f64 {
    assert_eq!(actual.len(), expected.len(), "{label}: key count");
    let mut scales: HashMap<String, f64> = HashMap::new();
    for (key, values) in expected {
        let scale = scales.entry(group(key)).or_default();
        *scale = values.iter().fold(*scale, |m, v| m.max(v.abs()));
    }
    let mut largest = 0.0_f64;
    for (key, values) in expected {
        let other = actual
            .get(key)
            .unwrap_or_else(|| panic!("{label}: missing {key:?}"));
        assert_eq!(values.len(), other.len(), "{label}: {key:?} length");
        let scale = scales[&group(key)];
        for (a, b) in other.iter().zip(values.iter()) {
            let deviation = (a - b).abs();
            assert!(
                deviation <= SYMMETRIC_BOUND * scale,
                "{label}: {key:?}: {a} vs {b} (scale {scale})"
            );
            if scale > 0.0 {
                largest = largest.max(deviation / scale);
            }
        }
    }
    largest
}

/// A y-duplicated wing with a symmetric elevator and an antisymmetric aileron.
fn accelerated_configuration(
    extra: &[(&str, serde_json::Value)],
) -> HashMap<String, serde_json::Value> {
    let controls = serde_json::json!([
        {"surface": "Wing", "name": "elevator", "span_start": 0, "span_end": 3,
         "hinge": [0.6, 0.8, 0.9, 1.1]},
        {"surface": "Wing", "name": "aileron", "span_start": 1, "span_end": 3,
         "hinge": 1.2, "duplicate_sign": -1.0}
    ]);
    let mut configuration = options(&[
        ("surface_names", "Wing".into()),
        ("surface_nx", "2".into()),
        ("surface_ny", "3".into()),
        ("surface_y_duplicated", "true".into()),
        ("controls", controls),
    ]);
    configuration.extend(options(extra));
    configuration
}

/// Groups output values by themselves.
fn by_output(key: &str) -> String {
    key.to_owned()
}

/// Groups partials by Jacobian row.
fn by_row(key: &(String, String)) -> String {
    key.0.clone()
}

fn configured(configuration: &HashMap<String, serde_json::Value>) -> VlmDiscipline {
    let mut discipline = VlmDiscipline::default();
    discipline.initialize().unwrap();
    discipline.set_options(configuration).unwrap();
    discipline.setup().unwrap();
    discipline.setup_partials().unwrap();
    discipline
}

/// Nonzero sideslip, rates and Mach so every output and partial is populated;
/// none of them enters the system matrix.
fn accelerated_inputs(aileron: f64) -> ArrayMap {
    let mut inputs = vlm_inputs(2, 3, 1.0, 3.0);
    inputs.insert("beta".into(), array(&[0.03]));
    inputs.insert("mach".into(), array(&[0.3]));
    inputs.insert("roll_rate".into(), array(&[0.01]));
    inputs.insert("pitch_rate".into(), array(&[0.005]));
    inputs.insert("yaw_rate".into(), array(&[-0.004]));
    inputs.insert("control_elevator".into(), array(&[0.13]));
    inputs.insert("control_aileron".into(), array(&[aileron]));
    inputs
}

/// Outputs, partials and the applies flag of one discipline at one point.
async fn evaluate(
    configuration: &HashMap<String, serde_json::Value>,
    inputs: &ArrayMap,
) -> (ArrayMap, PartialMap, Option<bool>) {
    let discipline = configured(configuration);
    let partials = discipline.compute_partials(inputs).await.unwrap();
    let outputs = discipline.compute(inputs).await.unwrap();
    let applies = discipline.symmetric_factorization_applies().unwrap();
    (outputs, partials, applies)
}

#[test]
fn acceleration_options_are_declared_as_booleans() {
    let mut discipline = VlmDiscipline::default();
    discipline.initialize().unwrap();
    let declared = discipline.get_available_options().unwrap();
    for name in ["mirror_reuse", "symmetric_factorization"] {
        assert_eq!(
            declared.get(name).map(String::as_str),
            Some("bool"),
            "{name}"
        );
    }
}

#[tokio::test]
async fn acceleration_options_default_off_and_parse_boolean_spellings() {
    let inputs = accelerated_inputs(0.0);
    // Off by default: the mesh qualifies, yet the full LU is used.
    let default = configured(&accelerated_configuration(&[]));
    assert_eq!(default.symmetric_factorization_applies().unwrap(), None);
    default.compute(&inputs).await.unwrap();
    assert_eq!(
        default.symmetric_factorization_applies().unwrap(),
        Some(false)
    );

    let spellings: [(serde_json::Value, bool); 10] = [
        (true.into(), true),
        (false.into(), false),
        ("true".into(), true),
        ("False".into(), false),
        (" yes ".into(), true),
        ("no".into(), false),
        ("1".into(), true),
        ("0".into(), false),
        (1.into(), true),
        (0.into(), false),
    ];
    for (value, expected) in spellings {
        let discipline = configured(&accelerated_configuration(&[(
            "symmetric_factorization",
            value.clone(),
        )]));
        discipline.compute(&inputs).await.unwrap();
        assert_eq!(
            discipline.symmetric_factorization_applies().unwrap(),
            Some(expected),
            "{value}"
        );
        // Mirror reuse never changes results, so a parse is all that shows.
        configured(&accelerated_configuration(&[("mirror_reuse", value)]));
    }
}

#[test]
fn acceleration_options_reject_values_that_are_not_booleans() {
    for name in ["mirror_reuse", "symmetric_factorization"] {
        for value in [
            serde_json::json!("maybe"),
            serde_json::json!(""),
            serde_json::json!(2),
            serde_json::json!(0.5),
            serde_json::json!(null),
            serde_json::json!([true]),
            serde_json::json!({"enabled": true}),
        ] {
            let mut discipline = VlmDiscipline::default();
            discipline.initialize().unwrap();
            let error = discipline
                .set_options(&accelerated_configuration(&[(name, value.clone())]))
                .expect_err(&format!("{name} = {value} must be rejected"));
            let message = error.to_string();
            assert!(
                message.contains(name) && message.contains("true or false"),
                "{name} = {value}: {message}"
            );
        }
    }
}

#[tokio::test]
async fn mirror_reuse_is_bit_identical_on_the_ordinary_path() {
    let all_outputs = serde_json::json!([
        "CL", "CD", "CY", "Cl", "Cm", "Cn", "Fx", "Fy", "Fz", "Mx", "My", "Mz"
    ]);
    // Coefficients take the adjoint branch, dimensional outputs the
    // linearization branch; aileron deflection takes the deflected system.
    for (extra, aileron) in [
        (vec![], 0.0),
        (vec![], -0.05),
        (vec![("outputs", all_outputs.clone())], 0.0),
        (vec![("outputs", all_outputs)], -0.05),
    ] {
        let inputs = accelerated_inputs(aileron);
        let (outputs, partials, _) = evaluate(&accelerated_configuration(&extra), &inputs).await;
        let mut on = extra.clone();
        on.push(("mirror_reuse", true.into()));
        let (fast, fast_partials, _) = evaluate(&accelerated_configuration(&on), &inputs).await;
        assert_bits("outputs", &fast, &outputs);
        assert_bits("partials", &fast_partials, &partials);
        // Also with the symmetric factorization on: mirror reuse changes no
        // bit of what that produces either.
        let mut both = on.clone();
        both.push(("symmetric_factorization", true.into()));
        let mut symmetric = extra.clone();
        symmetric.push(("symmetric_factorization", true.into()));
        let (a, a_partials, _) = evaluate(&accelerated_configuration(&both), &inputs).await;
        let (b, b_partials, _) = evaluate(&accelerated_configuration(&symmetric), &inputs).await;
        assert_bits("outputs with symmetric", &a, &b);
        assert_bits("partials with symmetric", &a_partials, &b_partials);
    }
}

#[tokio::test]
async fn symmetric_factorization_matches_the_full_lu_on_the_ordinary_path() {
    let all_outputs = serde_json::json!([
        "CL", "CD", "CY", "Cl", "Cm", "Cn", "Fx", "Fy", "Fz", "Mx", "My", "Mz"
    ]);
    let mut largest = 0.0_f64;
    for extra in [vec![], vec![("outputs", all_outputs)]] {
        let inputs = accelerated_inputs(0.0);
        let (outputs, partials, applies) =
            evaluate(&accelerated_configuration(&extra), &inputs).await;
        assert_eq!(applies, Some(false));
        let mut on = extra.clone();
        on.push(("symmetric_factorization", true.into()));
        let (fast, fast_partials, applies) =
            evaluate(&accelerated_configuration(&on), &inputs).await;
        assert_eq!(applies, Some(true));
        largest = largest.max(assert_near("outputs", &fast, &outputs, |k| by_output(k)));
        largest = largest.max(assert_near("partials", &fast_partials, &partials, by_row));
        // Control partials are present, so the deflected system was used.
        assert!(fast_partials.contains_key(&("CL".into(), "control_elevator".into())));
    }
    println!("ordinary path: largest symmetric-factorization deviation {largest:e}");
}

#[tokio::test]
async fn symmetric_factorization_falls_back_bit_for_bit_where_it_cannot_apply() {
    // An antisymmetric aileron deflection breaks the block structure.
    let inputs = accelerated_inputs(-0.05);
    let (outputs, partials, _) = evaluate(&accelerated_configuration(&[]), &inputs).await;
    let (fast, fast_partials, applies) = evaluate(
        &accelerated_configuration(&[("symmetric_factorization", true.into())]),
        &inputs,
    )
    .await;
    assert_eq!(applies, Some(false));
    assert_bits("antisymmetric outputs", &fast, &outputs);
    assert_bits("antisymmetric partials", &fast_partials, &partials);

    // A surface that is not y-duplicated has no mirror pairs.
    let plain = |extra: &[(&str, serde_json::Value)]| {
        let mut configuration = options(&[
            ("surface_names", "Wing".into()),
            ("surface_nx", "2".into()),
            ("surface_ny", "3".into()),
        ]);
        configuration.extend(options(extra));
        configuration
    };
    let inputs = vlm_inputs(2, 3, 1.0, 3.0);
    let (outputs, partials, _) = evaluate(&plain(&[]), &inputs).await;
    let (fast, fast_partials, applies) = evaluate(
        &plain(&[
            ("symmetric_factorization", true.into()),
            ("mirror_reuse", true.into()),
        ]),
        &inputs,
    )
    .await;
    assert_eq!(applies, Some(false));
    assert_bits("unduplicated outputs", &fast, &outputs);
    assert_bits("unduplicated partials", &fast_partials, &partials);
}

/// The selective interface on a y-duplicated wing in the given mode.
fn selective_configuration(
    mode: &str,
    extra: &[(&str, serde_json::Value)],
) -> HashMap<String, serde_json::Value> {
    let mut configuration = options(&[
        ("surface_names", serde_json::json!(["Wing"])),
        ("surface_nx", serde_json::json!([2])),
        ("surface_ny", serde_json::json!([3])),
        ("surface_y_duplicated", serde_json::json!([true])),
        ("analysis_mode", mode.into()),
        (
            "analysis_outputs",
            serde_json::json!(["CL", "CD", "Cm", "CDi", "e"]),
        ),
        (
            "analysis_gradients",
            serde_json::json!({
                "CL": ["alpha", "Wing", "mach", "reference_area"],
                "CD": ["Wing", "pitch_rate"],
                "Cm": ["alpha", "Wing", "moment_ref_x"],
                "CDi": ["Wing", "alpha"]
            }),
        ),
    ]);
    configuration.extend(options(extra));
    configuration
}

fn selective_inputs() -> ArrayMap {
    let mut inputs = vlm_inputs(2, 3, 1.0, 3.0);
    inputs.insert("mach".into(), array(&[0.3]));
    inputs.insert("pitch_rate".into(), array(&[0.005]));
    inputs
}

#[tokio::test]
async fn acceleration_options_reach_the_selective_path() {
    let inputs = selective_inputs();
    let mut largest = 0.0_f64;
    for mode in ["full", "half_wing"] {
        let (outputs, partials, applies) =
            evaluate(&selective_configuration(mode, &[]), &inputs).await;
        assert_eq!(applies, Some(false), "{mode}");
        assert!(!partials.is_empty());

        let (fast, fast_partials, _) = evaluate(
            &selective_configuration(mode, &[("mirror_reuse", true.into())]),
            &inputs,
        )
        .await;
        assert_bits(&format!("{mode} mirror outputs"), &fast, &outputs);
        assert_bits(
            &format!("{mode} mirror partials"),
            &fast_partials,
            &partials,
        );

        let symmetric = selective_configuration(mode, &[("symmetric_factorization", true.into())]);
        let (fast, fast_partials, applies) = evaluate(&symmetric, &inputs).await;
        let both = selective_configuration(
            mode,
            &[
                ("symmetric_factorization", true.into()),
                ("mirror_reuse", true.into()),
            ],
        );
        let (both_outputs, both_partials, _) = evaluate(&both, &inputs).await;
        assert_bits(&format!("{mode} both outputs"), &both_outputs, &fast);
        assert_bits(
            &format!("{mode} both partials"),
            &both_partials,
            &fast_partials,
        );
        if mode == "full" {
            assert_eq!(applies, Some(true));
            largest = largest.max(assert_near("selective outputs", &fast, &outputs, |k| {
                by_output(k)
            }));
            largest = largest.max(assert_near(
                "selective partials",
                &fast_partials,
                &partials,
                by_row,
            ));
        } else {
            // Half-wing mode solves one semispan, which has no mirror pairs.
            assert_eq!(applies, Some(false));
            assert_bits("half-wing outputs", &fast, &outputs);
            assert_bits("half-wing partials", &fast_partials, &partials);
        }
    }
    println!("selective full path: largest symmetric-factorization deviation {largest:e}");
}

#[tokio::test]
async fn acceleration_option_changes_invalidate_the_solved_state() {
    let inputs = accelerated_inputs(0.0);
    let off = accelerated_configuration(&[("symmetric_factorization", false.into())]);
    let on = accelerated_configuration(&[("symmetric_factorization", true.into())]);
    let mut discipline = configured(&off);
    let original = discipline.compute(&inputs).await.unwrap();
    let original_partials = discipline.compute_partials(&inputs).await.unwrap();
    assert_eq!(
        discipline.symmetric_factorization_applies().unwrap(),
        Some(false)
    );

    // After setup without trim, options are re-sent in full, as for every
    // other option; the cached state goes with them.
    discipline.set_options(&on).unwrap();
    assert_eq!(discipline.symmetric_factorization_applies().unwrap(), None);
    let partials = discipline.compute_partials(&inputs).await.unwrap();
    let changed = discipline.compute(&inputs).await.unwrap();
    assert_eq!(
        discipline.symmetric_factorization_applies().unwrap(),
        Some(true)
    );
    let (fresh, fresh_partials, _) = evaluate(&on, &inputs).await;
    assert_bits("flipped on", &changed, &fresh);
    assert_bits("flipped on partials", &partials, &fresh_partials);

    discipline.set_options(&off).unwrap();
    assert_eq!(discipline.symmetric_factorization_applies().unwrap(), None);
    assert_bits(
        "flipped back",
        &discipline.compute(&inputs).await.unwrap(),
        &original,
    );
    assert_bits(
        "flipped back partials",
        &discipline.compute_partials(&inputs).await.unwrap(),
        &original_partials,
    );
    assert_eq!(
        discipline.symmetric_factorization_applies().unwrap(),
        Some(false)
    );

    // Mirror reuse changes no result bit, so a stale state could only be seen
    // through the cache itself: the flag reads `None` once it is cleared.
    discipline
        .set_options(&accelerated_configuration(&[("mirror_reuse", true.into())]))
        .unwrap();
    assert_eq!(discipline.symmetric_factorization_applies().unwrap(), None);
    assert_bits(
        "mirror reuse",
        &discipline.compute(&inputs).await.unwrap(),
        &original,
    );
    // Options absent from a later request keep their value, like `outputs`.
    discipline
        .set_options(&accelerated_configuration(&[(
            "symmetric_factorization",
            true.into(),
        )]))
        .unwrap();
    discipline.compute(&inputs).await.unwrap();
    discipline
        .set_options(&accelerated_configuration(&[]))
        .unwrap();
    discipline.compute(&inputs).await.unwrap();
    assert_eq!(
        discipline.symmetric_factorization_applies().unwrap(),
        Some(true)
    );
}

#[tokio::test]
async fn acceleration_option_changes_invalidate_the_selective_analysis() {
    let inputs = selective_inputs();
    let off = selective_configuration("full", &[]);
    let on = selective_configuration("full", &[("symmetric_factorization", true.into())]);
    let mut discipline = configured(&off);
    let original = discipline.compute(&inputs).await.unwrap();
    assert_eq!(
        discipline.symmetric_factorization_applies().unwrap(),
        Some(false)
    );
    discipline.set_options(&on).unwrap();
    assert_eq!(discipline.symmetric_factorization_applies().unwrap(), None);
    let changed = discipline.compute(&inputs).await.unwrap();
    let partials = discipline.compute_partials(&inputs).await.unwrap();
    assert_eq!(
        discipline.symmetric_factorization_applies().unwrap(),
        Some(true)
    );
    let (fresh, fresh_partials, _) = evaluate(&on, &inputs).await;
    assert_bits("selective flipped on", &changed, &fresh);
    assert_bits("selective flipped on partials", &partials, &fresh_partials);

    discipline
        .set_options(&selective_configuration(
            "full",
            &[
                ("symmetric_factorization", false.into()),
                ("mirror_reuse", true.into()),
            ],
        ))
        .unwrap();
    assert_eq!(discipline.symmetric_factorization_applies().unwrap(), None);
    assert_bits(
        "selective flipped back",
        &discipline.compute(&inputs).await.unwrap(),
        &original,
    );
    assert_eq!(
        discipline.symmetric_factorization_applies().unwrap(),
        Some(false)
    );
}
