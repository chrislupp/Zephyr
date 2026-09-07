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
