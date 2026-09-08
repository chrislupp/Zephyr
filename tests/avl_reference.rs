//! Independent OptVL 2.1.1 results from the committed AVL geometry files.
//! These cases catch the v0.2.0 control/self-bound and force-location errors.
use approx::assert_abs_diff_eq;
use pegasus_vlm::aero::ForceCalculator;
use pegasus_vlm::analysis::{
    Analysis, AnalysisRequest, Details, GradientRequest, GradientValue, Input, Output, SolverMode,
};
use pegasus_vlm::core::{DEGREES_TO_RADIANS, SymmetryType};
use pegasus_vlm::io::AvlReader;
use pegasus_vlm::solver::{FlightCondition, VortexLatticeSolver};
use serde_json::Value;

#[test]
fn sections_match_optvl_circulation_loads_and_alpha_derivatives() {
    let reference: Value =
        serde_json::from_str(include_str!("fixtures/avl_reference/optvl.json")).unwrap();
    for case in reference["cases"].as_array().unwrap() {
        let name = case["name"].as_str().unwrap();
        let mut aircraft = AvlReader::default()
            .read_path(format!("tests/fixtures/avl_reference/{name}.avl"))
            .unwrap();
        let condition = FlightCondition {
            alpha: case["alpha_deg"].as_f64().unwrap() * DEGREES_TO_RADIANS,
            beta: case["beta_deg"].as_f64().unwrap() * DEGREES_TO_RADIANS,
            mach: case["mach"].as_f64().unwrap(),
            ..FlightCondition::default()
        };
        let symmetric = case["iysym"] == 1;
        if symmetric {
            // Represent AVL's image wing using Pegasus's explicit duplication
            // metadata; HalfWing then solves only the original panels.
            aircraft.set_y_symmetry(SymmetryType::None);
            for surface in aircraft.surfaces_mut() {
                surface.set_y_duplicate(0.0);
            }
        }
        let n = case["vortices"].as_u64().unwrap() as usize;
        let outputs = [
            (Output::Lift, "CL"),
            (Output::Drag, "CD"),
            (Output::Pitch, "Cm"),
        ];
        // The finest half-wing case exercises parallel assembly when enabled.
        let modes = if symmetric && n < 100 {
            vec![SolverMode::Full, SolverMode::HalfWing]
        } else if symmetric {
            vec![SolverMode::HalfWing]
        } else {
            vec![SolverMode::Full]
        };
        for mode in modes {
            let result = Analysis::from_aircraft(&aircraft)
                .evaluate(
                    &condition,
                    &AnalysisRequest {
                        mode,
                        outputs: outputs.iter().map(|(o, _)| *o).collect(),
                        gradients: outputs
                            .iter()
                            .map(|(o, _)| GradientRequest {
                                output: *o,
                                input: Input::Alpha,
                            })
                            .collect(),
                        details: Details {
                            circulation: true,
                            ..Details::default()
                        },
                    },
                )
                .unwrap();
            for (output, key) in outputs {
                let actual = result.outputs[&output];
                let expected = case["coefficients"][key].as_f64().unwrap();
                assert!(
                    (actual - expected).abs() < 2e-8,
                    "{name} {mode:?} {key}: {actual} vs OptVL {expected}"
                );
                let GradientValue::Scalar(actual) = result
                    .gradients
                    .iter()
                    .find(|g| g.request.output == output)
                    .unwrap()
                    .value
                else {
                    panic!()
                };
                let expected = case["alpha_derivatives"][key].as_f64().unwrap();
                assert!(
                    (actual - expected).abs() < 2e-7,
                    "{name} {mode:?} d{key}/dalpha: {actual} vs OptVL {expected}"
                );
            }
            for (actual, expected) in result
                .circulation
                .as_ref()
                .unwrap()
                .iter()
                .take(n)
                .zip(case["circulation"].as_array().unwrap())
            {
                assert_abs_diff_eq!(*actual, expected.as_f64().unwrap(), epsilon = 2e-8);
            }
        }
        if n < 100 {
            // Exercise the legacy solver/force path independently of Analysis.
            let mut solver = VortexLatticeSolver::new();
            solver.set_aircraft(&aircraft).unwrap();
            let solution = solver.solve(&condition).unwrap();
            let loads = ForceCalculator::new().compute_for_aircraft(
                &aircraft,
                solver.mesh(),
                &solution,
                &condition,
            );
            for (actual, key) in [
                (loads.coefficients.lift, "CL"),
                (loads.coefficients.drag, "CD"),
                (loads.coefficients.pitch, "Cm"),
            ] {
                assert_abs_diff_eq!(
                    actual,
                    case["coefficients"][key].as_f64().unwrap(),
                    epsilon = 2e-8
                );
            }
            if let Some(geometry) = case["geometry"].as_object() {
                for (i, vortex) in solver.mesh().vortices().iter().take(n).enumerate() {
                    for (key, actual) in [
                        ("RC", vortex.control_point),
                        ("RV", vortex.bound_vortex_center),
                        ("RV1", vortex.bound_vortex_left),
                        ("RV2", vortex.bound_vortex_right),
                        ("ENC", vortex.normal_at_control_point),
                    ] {
                        for axis in 0..3 {
                            assert_abs_diff_eq!(
                                actual[axis],
                                geometry[key][i][axis].as_f64().unwrap(),
                                epsilon = 2e-8
                            );
                        }
                    }
                }
            }
        }
    }
}
