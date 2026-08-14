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
use pegasus_vlm::philote::{MeshDiscipline, SurfaceOption, VlmDiscipline};
use philote_mdo::traits::{Discipline, ExplicitDiscipline};
use philote_mdo::{ArrayMap, PartialMap};

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
                    assert_abs_diff_eq!(
                        exact,
                        expected,
                        epsilon = tolerance * expected.abs() + NOISE_FLOOR
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
