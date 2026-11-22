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

//! Gradients that run through the fused horseshoe kernels.
//!
//! The geometry sweep contracts the bound-centre kernel either as one VJP per
//! output seed or through full Jacobians, depending on how many outputs need
//! it; the Mach derivative is evaluated as a matrix-free product. These tests
//! pin both sweeps to each other, to re-solves, and to the requirement that
//! threading never changes a bit.

use zephyr_vlm::aero::ForceCalculator;
use zephyr_vlm::core::{DEGREES_TO_RADIANS, Vec3};
use zephyr_vlm::geometry::{MeshInput, ReferenceQuantities, SurfaceMesh};
use zephyr_vlm::gradients::{
    FullGradients, GradientCalculator, OutputSelection, SurfaceMeshNodeGradients,
};
use zephyr_vlm::solver::{FlightCondition, VortexLatticeSolver};

/// A swept, dihedralled, cosine-spaced wing mirrored about `y = 0`, with
/// `2 * nx * ny` panels.
fn wing(nx: usize, ny: usize) -> MeshInput {
    let mut surface = SurfaceMesh::new("Wing", nx, ny);
    surface.y_duplicated = true;
    for i in 0..=nx {
        let x = 0.5 * (1.0 - (std::f64::consts::PI * i as f64 / nx as f64).cos());
        for j in 0..=ny {
            let y = 5.0 * j as f64 / ny as f64;
            let chord = 1.0 - 0.06 * y;
            surface
                .set_node(i, j, Vec3::new(chord * x + 0.2 * y, y, 0.05 * y))
                .unwrap();
        }
    }
    MeshInput {
        surfaces: vec![surface],
        reference: ReferenceQuantities {
            area: 10.0,
            chord: 1.0,
            span: 10.0,
            moment_reference_point: Vec3::new(0.25, 0.0, 0.0),
        },
    }
}

fn condition(mach: f64) -> FlightCondition {
    FlightCondition {
        alpha: 4.0 * DEGREES_TO_RADIANS,
        beta: 1.5 * DEGREES_TO_RADIANS,
        mach,
        roll_rate: 0.01,
        pitch_rate: 0.02,
        yaw_rate: -0.01,
        ..FlightCondition::default()
    }
}

/// Lift and induced drag only: the near-field sweep then carries one seed.
fn lift_and_induced_drag() -> OutputSelection {
    let mut selection = OutputSelection::none();
    selection.lift = true;
    selection.induced_drag = true;
    selection
}

fn node_gradients(
    input: &MeshInput,
    condition: &FlightCondition,
    selection: OutputSelection,
) -> Vec<SurfaceMeshNodeGradients> {
    let mut calculator = GradientCalculator::new();
    calculator.set_mesh(input).unwrap();
    calculator.set_output_selection(selection);
    calculator.compute_mesh_node_gradients(condition).unwrap()
}

fn lift(input: &MeshInput, condition: &FlightCondition) -> f64 {
    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(input).unwrap();
    let solution = solver.solve(condition).unwrap();
    ForceCalculator::new()
        .compute_with_references(&input.reference, solver.mesh(), &solution, condition)
        .coefficients
        .lift
}

fn max_abs(nodes: &[Vec3]) -> f64 {
    nodes.iter().fold(0.0, |m, v| m.max(v.abs().max()))
}

fn max_difference(a: &[Vec3], b: &[Vec3]) -> f64 {
    a.iter()
        .zip(b)
        .fold(0.0, |m, (x, y)| m.max((x - y).abs().max()))
}

#[test]
fn single_seed_and_jacobian_sweeps_agree_to_rounding() {
    let input = wing(4, 8);
    for mach in [0.0, 0.4] {
        let condition = condition(mach);
        let single = node_gradients(&input, &condition, lift_and_induced_drag());
        let all = node_gradients(&input, &condition, OutputSelection::all());
        for (one, every) in [
            (&single[0].lift.nodes, &all[0].lift.nodes),
            (&single[0].induced_drag.nodes, &all[0].induced_drag.nodes),
        ] {
            let scale = max_abs(every);
            assert!(scale > 0.0);
            assert!(
                max_difference(one, every) <= 1.0e-13 * scale,
                "Mach {mach}: {} vs scale {scale}",
                max_difference(one, every)
            );
        }
    }
}

#[test]
fn single_seed_node_gradients_match_a_directional_difference_at_mach() {
    let input = wing(3, 5);
    let condition = condition(0.5);
    let gradients = node_gradients(&input, &condition, lift_and_induced_drag());
    let direction: Vec<Vec3> = (0..input.surfaces[0].node_count())
        .map(|i| {
            let t = i as f64;
            Vec3::new(0.3 * t.sin(), 0.2 * (1.3 * t).cos(), 0.1 * (0.7 * t).sin())
        })
        .collect();
    let analytic: f64 = gradients[0]
        .lift
        .nodes
        .iter()
        .zip(&direction)
        .map(|(g, d)| g.dot(d))
        .sum();
    let step = 1.0e-6;
    let shifted = |sign: f64| {
        let mut moved = input.clone();
        for (node, d) in moved.surfaces[0].nodes.iter_mut().zip(&direction) {
            *node += sign * step * d;
        }
        lift(&moved, &condition)
    };
    let finite = (shifted(1.0) - shifted(-1.0)) / (2.0 * step);
    assert!(
        (analytic - finite).abs() <= 1.0e-6 * finite.abs() + 1.0e-9,
        "{analytic} vs {finite}"
    );
}

fn flight_bits(gradients: &FullGradients) -> Vec<u64> {
    let mut groups = vec![
        gradients.by_alpha,
        gradients.by_beta,
        gradients.by_roll_rate,
        gradients.by_pitch_rate,
        gradients.by_yaw_rate,
        gradients.by_mach,
    ];
    groups.extend(gradients.by_moment_reference);
    groups
        .iter()
        .flat_map(|g| {
            let c = g.coefficients;
            [
                c.lift,
                c.drag,
                c.side,
                c.roll,
                c.pitch,
                c.yaw,
                g.induced_drag,
            ]
        })
        .map(f64::to_bits)
        .collect()
}

fn node_bits(gradients: &[SurfaceMeshNodeGradients]) -> Vec<u64> {
    gradients
        .iter()
        .flat_map(|s| {
            [
                &s.lift,
                &s.drag,
                &s.side,
                &s.roll,
                &s.pitch,
                &s.yaw,
                &s.induced_drag,
            ]
        })
        .flat_map(|g| g.nodes.iter().flat_map(|v| [v.x, v.y, v.z]))
        .map(f64::to_bits)
        .collect()
}

/// Mach flight-condition gradients and both geometry-sweep branches, above the
/// parallel threshold, as raw bits.
fn fingerprint() -> Vec<u64> {
    let input = wing(4, 20);
    let condition = condition(0.3);
    let mut calculator = GradientCalculator::new();
    calculator.set_mesh(&input).unwrap();
    let mut bits = flight_bits(
        &calculator
            .compute_flight_condition_gradients(&condition)
            .unwrap(),
    );
    bits.extend(node_bits(&node_gradients(
        &input,
        &condition,
        lift_and_induced_drag(),
    )));
    bits.extend(node_bits(&node_gradients(
        &input,
        &condition,
        OutputSelection::all(),
    )));
    bits
}

#[test]
fn threaded_mach_and_geometry_gradients_are_bit_identical() {
    let evaluate = |threads: usize| {
        #[cfg(feature = "parallel")]
        {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(fingerprint)
        }
        #[cfg(not(feature = "parallel"))]
        {
            let _ = threads;
            fingerprint()
        }
    };
    let serial = evaluate(1);
    assert!(serial.iter().any(|bits| *bits != 0));
    // Without the feature there is only one code path to compare.
    #[cfg(feature = "parallel")]
    for threads in [2, 3, 8] {
        assert_eq!(evaluate(threads), serial, "{threads} threads");
    }
}
