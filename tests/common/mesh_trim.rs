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

#![allow(dead_code)]
use zephyr_vlm::core::Vec3;
use zephyr_vlm::geometry::{MeshControlSurface, MeshInput, SurfaceMesh};
use zephyr_vlm::solver::{FlightCondition, VortexLatticeSolver};
use zephyr_vlm::trim::{
    AeroQuantity, MeshTrimProblem, TrimBalance, TrimVariable, aerodynamic_values,
};

pub fn mesh() -> MeshInput {
    let mut mesh = MeshInput::default();
    for (name, xle, span, chord) in [("Wing", 0.0, 3.0, 1.0), ("Tail", 3.0, 1.3, 0.65)] {
        let mut s = SurfaceMesh::new(name, 2, 3);
        s.y_duplicated = true;
        for i in 0..=2 {
            for j in 0..=3 {
                let eta = j as f64 / 3.0;
                let x = i as f64 / 2.0;
                s.set_node(
                    i,
                    j,
                    Vec3::new(
                        xle + 0.12 * eta + chord * x * (1.0 - 0.2 * eta),
                        span * eta,
                        0.08 * eta * eta + 0.015 * x * eta + if name == "Tail" { 0.2 } else { 0.0 },
                    ),
                )
                .unwrap();
            }
        }
        s.controls.push(MeshControlSurface::trailing_edge(
            if name == "Wing" { "flap" } else { "elevator" },
            0..3,
            vec![1.2; 4],
        ));
        if name == "Tail" {
            s.controls
                .push(MeshControlSurface::trailing_edge("tab", 1..3, vec![1.6; 4]));
        }
        mesh.surfaces.push(s);
    }
    mesh.reference.area = 7.0;
    mesh.reference.chord = 0.9;
    mesh.reference.span = 6.0;
    mesh.reference.moment_reference_point = Vec3::new(0.3, 0.0, -0.05);
    mesh
}

pub fn condition() -> FlightCondition {
    FlightCondition {
        alpha: 0.075,
        beta: 0.012,
        mach: 0.25,
        dynamic_pressure: 650.0,
        roll_rate: 0.002,
        pitch_rate: -0.003,
        yaw_rate: 0.001,
        control_deflections: vec![
            ("elevator".into(), 0.06),
            ("flap".into(), 0.025),
            ("tab".into(), -0.015),
        ],
    }
}

pub fn problem() -> MeshTrimProblem {
    let mut alpha = TrimVariable::new("alpha", 0.03);
    alpha.bounds = [-0.3, 0.3];
    alpha.scale = 0.1;
    let mut elevator = TrimVariable::new("control_elevator", 0.0);
    elevator.bounds = [-0.4, 0.4];
    elevator.scale = 0.2;
    MeshTrimProblem {
        variables: vec![alpha, elevator],
        balances: vec![
            TrimBalance::new(AeroQuantity::Fz, 1000.0),
            TrimBalance::new(AeroQuantity::My, 1000.0),
        ],
        options: zephyr_vlm::trim::MeshTrimOptions {
            absolute_tolerance: 1e-12,
            relative_tolerance: 0.0,
            ..zephyr_vlm::trim::MeshTrimOptions::default()
        },
        ..MeshTrimProblem::default()
    }
}

pub fn targets(
    mesh: &MeshInput,
    condition: &FlightCondition,
    problem: &MeshTrimProblem,
    external: &[f64; 6],
) -> Vec<f64> {
    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(mesh).unwrap();
    solver.solve(condition).unwrap();
    let v = aerodynamic_values(mesh, &solver, condition);
    problem
        .balances
        .iter()
        .map(|b| {
            v[b.quantity as usize]
                + if b.quantity as usize >= 6 {
                    external[b.quantity as usize - 6]
                } else {
                    0.0
                }
        })
        .collect()
}

pub fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 5e-7 + 3e-5 * expected.abs(),
        "{actual} != {expected}"
    );
}
