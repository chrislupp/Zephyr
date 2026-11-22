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

//! Opt-in performance settings of the selective analysis path, and the
//! bit-identity of its threaded kernels across thread counts.

use zephyr_vlm::analysis::{
    Analysis, AnalysisRequest, AnalysisResult, Details, GradientRequest, GradientValue, Input,
    Output, SolverMode,
};
use zephyr_vlm::core::Vec3;
use zephyr_vlm::geometry::{ControlSurfaceDefinition, MeshInput, SurfaceMesh};
use zephyr_vlm::solver::FlightCondition;

mod common;

fn mesh(nx: usize, ny: usize) -> MeshInput {
    let mut s = SurfaceMesh::new("wing", nx, ny);
    s.y_duplicated = true;
    for i in 0..=nx {
        for j in 0..=ny {
            let eta = j as f64 / ny as f64;
            let chord = i as f64 / nx as f64;
            s.set_node(
                i,
                j,
                Vec3::new(
                    0.7 * eta + chord * (1.0 - 0.3 * eta),
                    4.0 * eta,
                    0.4 * eta - 0.08 * chord * eta,
                ),
            )
            .unwrap();
        }
    }
    let mut m = MeshInput::default();
    m.surfaces.push(s);
    m.reference.area = 6.8;
    m.reference.chord = 0.9;
    m.reference.span = 8.0;
    m.reference.moment_reference_point = Vec3::new(0.25, 0.3, 0.12);
    m
}

fn condition() -> FlightCondition {
    FlightCondition {
        alpha: 0.09,
        mach: 0.3,
        pitch_rate: 0.002,
        ..FlightCondition::default()
    }
}

/// Every output, every detail array, and alpha, Mach and mesh-node gradients
/// of four near-field outputs (the Jacobian branch of the load term) and one
/// far-field output.
fn request(mode: SolverMode) -> AnalysisRequest {
    let gradient_outputs = [
        Output::Lift,
        Output::Drag,
        Output::Pitch,
        Output::NearFieldInducedDrag,
        Output::FarFieldInducedDrag,
    ];
    AnalysisRequest {
        mode,
        outputs: Output::ALL.to_vec(),
        gradients: gradient_outputs
            .into_iter()
            .flat_map(|output| {
                [Input::Alpha, Input::Mach, Input::MeshNodes(0)]
                    .into_iter()
                    .map(move |input| GradientRequest { output, input })
            })
            .collect(),
        details: Details {
            circulation: true,
            induced_velocity: true,
            surface_loads: true,
        },
    }
}

/// Primal results: outputs, circulation, induced velocity and surface loads.
fn primal(r: &AnalysisResult) -> Vec<u64> {
    let mut bits: Vec<_> = r.outputs.values().map(|v| v.to_bits()).collect();
    bits.extend(r.circulation.iter().flatten().map(|v| v.to_bits()));
    bits.extend(r.induced_velocity.iter().flatten().map(|v| v.to_bits()));
    for surface in r.surface_loads.iter().flatten() {
        // `Debug` prints the shortest round-tripping decimal of every float,
        // so equal strings mean equal bits.
        bits.extend(format!("{surface:?}").bytes().map(u64::from));
    }
    bits
}

fn gradients(r: &AnalysisResult) -> Vec<f64> {
    let mut values = Vec::new();
    for g in &r.gradients {
        match &g.value {
            GradientValue::Scalar(v) => values.push(*v),
            GradientValue::MeshNodes(nodes) => {
                values.extend(nodes.iter().flat_map(|v| v.iter().copied()));
            }
        }
    }
    values
}

fn all_bits(r: &AnalysisResult) -> Vec<u64> {
    let mut bits = primal(r);
    bits.extend(gradients(r).iter().map(|v| v.to_bits()));
    bits
}

fn evaluate(m: &MeshInput, mode: SolverMode, mirror: bool, symmetric: bool) -> AnalysisResult {
    let mut analysis = Analysis::from_mesh(m);
    analysis.set_mirror_reuse(mirror);
    analysis.set_symmetric_factorization(symmetric);
    analysis.evaluate(&condition(), &request(mode)).unwrap()
}

/// Largest deviation of each gradient relative to its own largest entry.
fn max_gradient_deviation(actual: &AnalysisResult, expected: &AnalysisResult) -> f64 {
    let mut worst = 0.0_f64;
    for (a, e) in actual.gradients.iter().zip(&expected.gradients) {
        assert_eq!(a.request, e.request);
        let (a, e) = match (&a.value, &e.value) {
            (GradientValue::Scalar(a), GradientValue::Scalar(e)) => (vec![*a], vec![*e]),
            (GradientValue::MeshNodes(a), GradientValue::MeshNodes(e)) => (
                a.iter().flat_map(|v| v.iter().copied()).collect(),
                e.iter().flat_map(|v| v.iter().copied()).collect(),
            ),
            _ => panic!("gradient kinds differ"),
        };
        let scale = e.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        let error = a
            .iter()
            .zip(&e)
            .fold(0.0_f64, |m, (x, y)| m.max((x - y).abs()));
        if scale > 0.0 {
            worst = worst.max(error / scale);
        } else {
            assert_eq!(error, 0.0);
        }
    }
    worst
}

#[test]
fn settings_default_off_and_survive_a_source_change() {
    let m = mesh(2, 4);
    let mut analysis = Analysis::from_mesh(&m);
    assert!(!analysis.mirror_reuse());
    assert!(!analysis.symmetric_factorization());
    assert!(!analysis.symmetric_factorization_applies());
    analysis.set_mirror_reuse(true);
    analysis.set_symmetric_factorization(true);
    analysis.set_mesh(&mesh(3, 4));
    assert!(analysis.mirror_reuse());
    assert!(analysis.symmetric_factorization());
    analysis.set_aircraft(&common::tapered_twisted_wing_mirrored());
    assert!(analysis.mirror_reuse());
    assert!(analysis.symmetric_factorization());
    analysis.set_mirror_reuse(false);
    analysis.set_symmetric_factorization(false);
    assert!(!analysis.mirror_reuse());
    assert!(!analysis.symmetric_factorization());
}

/// Mirror reuse changes nothing but the work: every primal value and every
/// gradient is the same bits, in both modes. Full mode evaluates half the
/// kernels; a half-wing lattice has no mirrored panels, so there it is a no-op.
#[test]
fn mirror_reuse_is_bit_identical_and_halves_full_mode_assembly() {
    for (nx, ny) in [(2, 5), (6, 24)] {
        let m = mesh(nx, ny);
        for mode in [SolverMode::Full, SolverMode::HalfWing] {
            let plain = evaluate(&m, mode, false, false);
            let reused = evaluate(&m, mode, true, false);
            assert_eq!(all_bits(&reused), all_bits(&plain), "{mode:?}");
            let expected = if mode == SolverMode::Full {
                plain.work.influence_evaluations / 2
            } else {
                plain.work.influence_evaluations
            };
            assert_eq!(reused.work.influence_evaluations, expected);
        }
    }
}

/// Full mode factorizes the mirror-symmetric system by blocks and agrees
/// with the full LU to rounding. Half-wing mode has no mirror pairs, so the
/// setting leaves every bit unchanged there.
#[test]
fn symmetric_factorization_applies_in_full_mode_only() {
    for (nx, ny) in [(2, 5), (6, 24)] {
        let m = mesh(nx, ny);
        let full = evaluate(&m, SolverMode::Full, false, false);
        for mirror in [false, true] {
            let mut analysis = Analysis::from_mesh(&m);
            analysis.set_mirror_reuse(mirror);
            analysis.set_symmetric_factorization(true);
            let block = analysis
                .evaluate(&condition(), &request(SolverMode::Full))
                .unwrap();
            assert!(analysis.symmetric_factorization_applies());
            for (output, value) in &full.outputs {
                let actual = block.outputs[output];
                assert!(
                    (actual - value).abs() <= 1e-13 * value.abs().max(1e-3),
                    "{output:?}: {actual:e} != {value:e}"
                );
            }
            let gamma = full.circulation.as_ref().unwrap();
            let scale = gamma.amax();
            assert!((block.circulation.as_ref().unwrap() - gamma).amax() <= 1e-13 * scale);
            assert!(max_gradient_deviation(&block, &full) <= 1e-12);

            let half = evaluate(&m, SolverMode::HalfWing, false, false);
            let half_block = analysis
                .evaluate(&condition(), &request(SolverMode::HalfWing))
                .unwrap();
            assert!(!analysis.symmetric_factorization_applies());
            assert_eq!(all_bits(&half_block), all_bits(&half));
        }
    }
}

/// A setting change discards the cached operator, factorization and results;
/// setting the current value again keeps them.
#[test]
fn changing_a_setting_invalidates_the_cached_solve() {
    let m = mesh(3, 6);
    let c = condition();
    let req = request(SolverMode::Full);
    let mut analysis = Analysis::from_mesh(&m);
    let first = analysis.evaluate(&c, &req).unwrap();
    assert!(!analysis.symmetric_factorization_applies());

    analysis.set_symmetric_factorization(false);
    analysis.set_mirror_reuse(false);
    let cached = analysis.evaluate(&c, &req).unwrap();
    assert_eq!(cached.work.influence_evaluations, 0);
    assert_eq!(cached.work.factorizations, 0);
    assert_eq!(cached.work.primal_solves, 0);
    assert_eq!(cached.work.adjoint_solves, 0);
    assert_eq!(all_bits(&cached), all_bits(&first));

    analysis.set_symmetric_factorization(true);
    assert!(!analysis.symmetric_factorization_applies());
    let block = analysis.evaluate(&c, &req).unwrap();
    assert!(analysis.symmetric_factorization_applies());
    assert!(block.work.influence_evaluations > 0);
    assert_eq!(block.work.factorizations, 1);
    assert_eq!(block.work.primal_solves, 1);
    assert!(block.work.adjoint_solves > 0);

    analysis.set_mirror_reuse(true);
    let reused = analysis.evaluate(&c, &req).unwrap();
    assert_eq!(
        reused.work.influence_evaluations,
        block.work.influence_evaluations / 2
    );
    assert_eq!(reused.work.factorizations, 1);
    assert_eq!(all_bits(&reused), all_bits(&block));

    analysis.set_symmetric_factorization(false);
    analysis.set_mirror_reuse(false);
    let restored = analysis.evaluate(&c, &req).unwrap();
    assert!(!analysis.symmetric_factorization_applies());
    assert_eq!(restored.work.factorizations, 1);
    assert_eq!(all_bits(&restored), all_bits(&first));

    // A setting change before the first solve has nothing to invalidate.
    let mut fresh = Analysis::from_mesh(&m);
    fresh.set_symmetric_factorization(true);
    fresh.evaluate(&c, &req).unwrap();
    assert!(fresh.symmetric_factorization_applies());
}

/// Control deflections rotate the normals the full-mode operator contracts
/// with. Symmetric deflections keep the block structure; mirror reuse stays
/// exact either way.
#[test]
fn deflected_controls_keep_both_settings_consistent() {
    let mut aircraft = common::tapered_twisted_wing_mirrored();
    for section in aircraft.surfaces_mut()[0].sections_mut() {
        section.add_control_surface(ControlSurfaceDefinition {
            name: "flap".into(),
            hinge_x_fraction: 0.7,
            hinge_vector: Vec3::y(),
            sign: 1.0,
            gain: 1.0,
        });
    }
    let mut c = condition();
    c.control_deflections = vec![("flap".into(), 0.07)];
    let req = AnalysisRequest {
        outputs: vec![Output::Lift, Output::Drag, Output::Pitch],
        gradients: vec![
            GradientRequest {
                output: Output::Lift,
                input: Input::Control("flap".into()),
            },
            GradientRequest {
                output: Output::Pitch,
                input: Input::Mach,
            },
        ],
        details: Details {
            circulation: true,
            ..Details::default()
        },
        ..AnalysisRequest::default()
    };
    let run = |mirror: bool, symmetric: bool| {
        let mut analysis = Analysis::from_aircraft(&aircraft);
        analysis.set_mirror_reuse(mirror);
        analysis.set_symmetric_factorization(symmetric);
        let result = analysis.evaluate(&c, &req).unwrap();
        (result, analysis.symmetric_factorization_applies())
    };
    let (plain, plain_applies) = run(false, false);
    let (reused, _) = run(true, false);
    let (block, block_applies) = run(true, true);
    assert!(!plain_applies);
    assert!(block_applies);
    assert_eq!(all_bits(&reused), all_bits(&plain));
    for (output, value) in &plain.outputs {
        assert!((block.outputs[output] - value).abs() <= 1e-13 * value.abs());
    }
    assert!(max_gradient_deviation(&block, &plain) <= 1e-12);
}

/// At Mach zero the Mach derivative needs no lattice work in either mode.
#[test]
fn incompressible_mach_gradients_skip_the_lattice() {
    let m = mesh(2, 5);
    let mut c = condition();
    c.mach = 0.0;
    for mode in [SolverMode::Full, SolverMode::HalfWing] {
        let result = Analysis::from_mesh(&m)
            .evaluate(
                &c,
                &AnalysisRequest {
                    mode,
                    gradients: vec![GradientRequest {
                        output: Output::Lift,
                        input: Input::Mach,
                    }],
                    ..AnalysisRequest::default()
                },
            )
            .unwrap();
        assert_eq!(result.work.mach_interactions, 0);
        let GradientValue::Scalar(value) = result.gradients[0].value else {
            panic!("scalar expected")
        };
        assert_eq!(value, 0.0);
    }
}

/// Every threaded kernel of the selective path -- full-mode assembly with and
/// without mirror reuse, half-wing assembly, the Mach products, the geometry
/// reverse sweep -- gives the same bits at 1, 2, 3 and 8 threads, and the
/// same bits as the serial build gives. Both meshes put each mode above the
/// threading threshold, the first only in full mode.
#[cfg(feature = "parallel")]
#[test]
fn threaded_kernels_are_independent_of_thread_count() {
    for (nx, ny) in [(4, 20), (6, 24)] {
        let m = mesh(nx, ny);
        for mode in [SolverMode::Full, SolverMode::HalfWing] {
            for (mirror, symmetric) in [(false, false), (true, false), (true, true)] {
                let mut reference = None;
                for threads in [1, 2, 3, 8] {
                    let pool = rayon::ThreadPoolBuilder::new()
                        .num_threads(threads)
                        .build()
                        .unwrap();
                    let actual = pool.install(|| all_bits(&evaluate(&m, mode, mirror, symmetric)));
                    if let Some(expected) = &reference {
                        assert_eq!(&actual, expected, "{mode:?} at {threads} threads");
                    } else {
                        reference = Some(actual);
                    }
                }
            }
        }
    }
}
