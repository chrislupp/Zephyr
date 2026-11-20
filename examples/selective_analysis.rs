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

//! Run with `cargo run --release --example selective_analysis`.
use zephyr_vlm::analysis::{Analysis, AnalysisRequest, GradientRequest, Input, Output, SolverMode};
use zephyr_vlm::core::Vec3;
use zephyr_vlm::geometry::{MeshInput, SurfaceMesh};
use zephyr_vlm::solver::FlightCondition;

fn main() -> zephyr_vlm::Result<()> {
    let mut wing = SurfaceMesh::new("wing", 4, 16);
    wing.y_duplicated = true;
    for i in 0..=wing.nx {
        for j in 0..=wing.ny {
            wing.set_node(i, j, Vec3::new(i as f64 / 4.0, 5.0 * j as f64 / 16.0, 0.0))?;
        }
    }
    let mut input = MeshInput::default();
    input.surfaces.push(wing);
    input.reference.area = 10.0;
    input.reference.span = 10.0;
    let condition = FlightCondition {
        alpha: 5.0_f64.to_radians(),
        ..FlightCondition::default()
    };
    let request = AnalysisRequest {
        mode: SolverMode::HalfWing,
        outputs: vec![Output::Lift, Output::FarFieldInducedDrag],
        gradients: vec![GradientRequest {
            output: Output::Lift,
            input: Input::Alpha,
        }],
        ..AnalysisRequest::default()
    };
    let mut analysis = Analysis::from_mesh(&input);
    let result = analysis.evaluate(&condition, &request)?;
    println!(
        "CL = {:.10}, CDi = {:.10}",
        result.outputs[&Output::Lift],
        result.outputs[&Output::FarFieldInducedDrag]
    );
    println!("dCL/dalpha: {:?}", result.gradients[0].value);
    println!("Work: {:?}", result.work);
    println!(
        "Repeated request: {:?}",
        analysis.evaluate(&condition, &request)?.work
    );
    Ok(())
}
