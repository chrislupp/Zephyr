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

//! The strip-load adjoint's pairwise sweep, by functional count and thread
//! count.
//!
//! The sweep contracts the per-pair kernel terms with all functionals in one
//! dense product per group of source panels. A product of a handful of columns
//! and one of many take different code paths through the matrix kernel, so
//! these check both against re-solve differences and against each other, on a
//! mirrored, deflected, compressible case whose block layout is ragged.

use zephyr_vlm::aero::{AerodynamicResults, ForceCalculator};
use zephyr_vlm::core::Vec3;
use zephyr_vlm::geometry::{MeshControlSurface, MeshInput, SurfaceMesh};
use zephyr_vlm::gradients::{
    GradientCalculator, MeshNodeCotangent, OutputSelection, StripLoadAdjoint, StripLoadCotangent,
};
use zephyr_vlm::solver::{FlightCondition, VortexLatticeSolver};

/// A swept, tapered, cambered surface with two deflected controls, mirrored.
fn mesh(nx: usize, ny: usize) -> MeshInput {
    let mut surface = SurfaceMesh::new("Wing", nx, ny);
    surface.y_duplicated = true;
    for i in 0..=nx {
        for j in 0..=ny {
            let eta = j as f64 / ny as f64;
            let x = i as f64 / nx as f64;
            surface
                .set_node(
                    i,
                    j,
                    Vec3::new(
                        0.25 * eta + (1.0 - 0.35 * eta) * x,
                        0.1 + 3.0 * eta,
                        0.12 * eta * eta + 0.05 * x * (1.0 - x),
                    ),
                )
                .unwrap();
        }
    }
    let hinge = nx as f64 * 0.7;
    surface.controls.push(MeshControlSurface::trailing_edge(
        "flap",
        0..ny / 2,
        vec![hinge; ny + 1],
    ));
    let mut aileron = MeshControlSurface::trailing_edge("aileron", ny / 2..ny, vec![hinge; ny + 1]);
    aileron.duplicate_sign = -1.0;
    aileron.gain = -0.8;
    surface.controls.push(aileron);
    let mut input = MeshInput::default();
    input.surfaces.push(surface);
    input.reference.area = 5.5;
    input.reference.chord = 0.8;
    input.reference.span = 6.2;
    input.reference.moment_reference_point = Vec3::new(0.3, 0.05, -0.02);
    input
}

fn flight() -> FlightCondition {
    FlightCondition {
        alpha: 0.07,
        beta: -0.03,
        mach: 0.35,
        roll_rate: 0.004,
        pitch_rate: -0.006,
        yaw_rate: 0.003,
        control_deflections: vec![("flap".into(), 0.09), ("aileron".into(), -0.06)],
        ..FlightCondition::default()
    }
}

fn points(strip_count: usize) -> Vec<Vec3> {
    (0..strip_count)
        .map(|strip| {
            let fraction = (strip as f64 + 0.5) / strip_count as f64;
            Vec3::new(0.3 + 0.1 * fraction, 6.0 * fraction - 3.0, 0.04)
        })
        .collect()
}

/// `count` functionals with no shared structure between columns.
fn columns(strip_count: usize, count: usize) -> Vec<StripLoadCotangent> {
    (0..count)
        .map(|column| {
            let mut cotangent = StripLoadCotangent::zeros(strip_count);
            for strip in 0..strip_count {
                let t = 0.37 * (strip + 1) as f64 + 1.13 * column as f64;
                cotangent.force[strip] = Vec3::new(t.sin(), (1.3 * t).cos(), 0.5 + (0.7 * t).sin());
                cotangent.moment[strip] =
                    Vec3::new((0.9 * t).cos(), 0.3 * t.sin(), (1.7 * t).sin());
            }
            cotangent
        })
        .collect()
}

fn solve(
    input: &MeshInput,
    condition: &FlightCondition,
) -> (VortexLatticeSolver, AerodynamicResults) {
    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(input).unwrap();
    let solution = solver.solve(condition).unwrap();
    let result = ForceCalculator::new().compute_with_references(
        &input.reference,
        solver.mesh(),
        &solution,
        condition,
    );
    (solver, result)
}

/// The functional as a coupled solver evaluates it from the reported strip
/// loads, redimensionalized and transported to the functional's moment point.
fn functional(result: &AerodynamicResults, points: &[Vec3], column: &StripLoadCotangent) -> f64 {
    result
        .surface_loads
        .iter()
        .flat_map(|surface| &surface.strip_loads)
        .map(|loads| {
            let force = loads.area * loads.body_force;
            let moment = loads.area
                * loads.chord
                * Vec3::new(
                    -loads.body_moment.x,
                    loads.body_moment.y,
                    -loads.body_moment.z,
                )
                + (loads.quarter_chord_position - points[loads.strip_index]).cross(&force);
            column.force[loads.strip_index].dot(&force)
                + column.moment[loads.strip_index].dot(&moment)
        })
        .sum()
}

fn flatten(nodes: &[MeshNodeCotangent]) -> Vec<f64> {
    nodes
        .iter()
        .flat_map(|node| node.surfaces.iter().flatten())
        .flat_map(|value| value.iter().copied())
        .collect()
}

/// Every public result the sweep feeds, flattened to bits.
fn all_results(input: &MeshInput, condition: &FlightCondition, count: usize) -> Vec<u64> {
    let (solver, _) = solve(input, condition);
    let strip_count = solver.mesh().strip_count();
    let points = points(strip_count);
    let columns = columns(strip_count, count);
    let adjoint = StripLoadAdjoint::new(&solver, condition, &points).unwrap();
    let mut values = flatten(&adjoint.mesh_node_cotangents(input, &columns).unwrap());
    for total in adjoint.flight_condition_cotangents(&columns).unwrap() {
        values.extend([total.by_alpha, total.by_beta]);
    }
    for per_column in adjoint.control_cotangents(&columns).unwrap() {
        values.extend(per_column.iter().map(|(_, value)| *value));
    }
    let mut calculator = GradientCalculator::new();
    calculator.set_mesh(input).unwrap();
    let mut selection = OutputSelection::all();
    selection.induced_drag = false;
    calculator.set_output_selection(selection);
    for surface in calculator
        .mesh_node_gradients_from(&solver, condition)
        .unwrap()
    {
        for output in [
            &surface.lift,
            &surface.drag,
            &surface.side,
            &surface.roll,
            &surface.pitch,
            &surface.yaw,
        ] {
            values.extend(output.nodes.iter().flat_map(|value| value.iter().copied()));
        }
    }
    values.into_iter().map(f64::to_bits).collect()
}

/// Many functionals at once against re-solve differences: the dense-product
/// path, with deflected controls, compressibility and a mirrored surface.
#[test]
fn many_functionals_match_resolve_differences() {
    let input = mesh(2, 3);
    let condition = flight();
    let (solver, _) = solve(&input, &condition);
    let strip_count = solver.mesh().strip_count();
    let points = points(strip_count);
    let columns = columns(strip_count, 7);
    let adjoint = StripLoadAdjoint::new(&solver, &condition, &points).unwrap();
    let nodes = adjoint.mesh_node_cotangents(&input, &columns).unwrap();
    assert_eq!(nodes.len(), columns.len());

    let step = 1.0e-7;
    let mut worst = 0.0f64;
    for node in 0..input.surfaces[0].node_count() {
        for direction in 0..3 {
            let mut plus = input.clone();
            let mut minus = input.clone();
            plus.surfaces[0].nodes[node][direction] += step;
            minus.surfaces[0].nodes[node][direction] -= step;
            let (_, high) = solve(&plus, &condition);
            let (_, low) = solve(&minus, &condition);
            for (column, cotangent) in columns.iter().enumerate() {
                let difference = (functional(&high, &points, cotangent)
                    - functional(&low, &points, cotangent))
                    / (2.0 * step);
                let analytical = nodes[column].surfaces[0][node][direction];
                let error = (analytical - difference).abs()
                    / (1.0 + difference.abs().max(analytical.abs()));
                worst = worst.max(error);
                assert!(
                    error < 1.0e-6,
                    "column {column}, node {node}, direction {direction}: \
                     adjoint {analytical:e} vs difference {difference:e}"
                );
            }
        }
    }
    assert!(worst > 0.0, "every component was identically zero");
}

/// A batch of functionals and the same functionals one at a time take
/// different paths through the dense product; they must agree to rounding, and
/// no column may leak into another.
#[test]
fn a_batch_agrees_with_one_functional_at_a_time() {
    // 7 x 10 x 2 = 140 panels: above the parallel threshold, and neither a
    // multiple of the source grouping nor of the block count.
    let input = mesh(7, 10);
    let condition = flight();
    let (solver, _) = solve(&input, &condition);
    assert_eq!(solver.mesh().vortex_count(), 140);
    let strip_count = solver.mesh().strip_count();
    let points = points(strip_count);
    let columns = columns(strip_count, 8);
    let adjoint = StripLoadAdjoint::new(&solver, &condition, &points).unwrap();
    let batch = adjoint.mesh_node_cotangents(&input, &columns).unwrap();
    let batch_controls = adjoint.control_cotangents(&columns).unwrap();
    let batch_angles = adjoint.flight_condition_cotangents(&columns).unwrap();
    for (column, cotangent) in columns.iter().enumerate() {
        let single = std::slice::from_ref(cotangent);
        let alone = flatten(&adjoint.mesh_node_cotangents(&input, single).unwrap());
        let together = flatten(std::slice::from_ref(&batch[column]));
        let scale = alone.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        assert!(scale > 0.0);
        for (a, b) in alone.iter().zip(&together) {
            assert!(
                (a - b).abs() <= 1.0e-12 * scale,
                "column {column}: {a:e} alone vs {b:e} in the batch"
            );
        }
        let controls = &adjoint.control_cotangents(single).unwrap()[0];
        for ((name, a), (other, b)) in controls.iter().zip(&batch_controls[column]) {
            assert_eq!(name, other);
            assert!(
                (a - b).abs() <= 1.0e-12 * (1.0 + a.abs()),
                "{name}: {a:e} vs {b:e}"
            );
        }
        let angles = adjoint.flight_condition_cotangents(single).unwrap()[0];
        for (a, b) in [
            (angles.by_alpha, batch_angles[column].by_alpha),
            (angles.by_beta, batch_angles[column].by_beta),
        ] {
            assert!((a - b).abs() <= 1.0e-12 * (1.0 + a.abs()), "{a:e} vs {b:e}");
        }
    }
}

/// Every public result is the same bits at every pool size, for one
/// functional, a few, and many.
#[test]
fn every_result_is_independent_of_the_thread_count() {
    let input = mesh(7, 10);
    let condition = flight();
    for count in [1, 7] {
        let evaluate = |threads: usize| {
            #[cfg(feature = "parallel")]
            {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap()
                    .install(|| all_results(&input, &condition, count))
            }
            #[cfg(not(feature = "parallel"))]
            {
                let _ = threads;
                all_results(&input, &condition, count)
            }
        };
        let serial = evaluate(1);
        assert!(!serial.is_empty());
        // Without the parallel feature there is no pool to vary; one repeat
        // still checks the result is deterministic.
        let pools: &[usize] = if cfg!(feature = "parallel") {
            &[2, 3, 8]
        } else {
            &[1]
        };
        for &threads in pools {
            assert_eq!(
                serial,
                evaluate(threads),
                "{count} functionals differ at {threads} threads"
            );
        }
    }
}
