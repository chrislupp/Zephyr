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

use std::path::PathBuf;

use zephyr_vlm::aero::ForceCalculator;
use zephyr_vlm::core::{DEGREES_TO_RADIANS, Vec3};
use zephyr_vlm::geometry::{Aircraft, Section, Surface};
use zephyr_vlm::solver::{FlightCondition, VortexLatticeSolver};
use zephyr_vlm::visualization::{VtkExportOptions, VtkExporter, is_available};

fn main() -> zephyr_vlm::Result<()> {
    if !is_available() {
        println!("VTK support is disabled in this build");
        return Ok(());
    }
    let output = std::env::args_os()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("vtk_output"));

    let mut wing = Surface::new("Tapered wing");
    wing.set_chordwise_panel_count(6);
    wing.set_spanwise_panel_count(16);
    wing.set_y_duplicate(0.0);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 2.0));
    wing.add_section(Section::new(Vec3::new(0.25, 5.0, 0.0), 0.5, -1.0));
    let mut aircraft = Aircraft::new("VTK example");
    aircraft.set_reference_area(7.5);
    aircraft.set_reference_chord(0.75);
    aircraft.set_reference_span(10.0);
    aircraft.add_surface(wing);

    let condition = FlightCondition {
        alpha: 5.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft)?;
    let solution = solver.solve(&condition)?;
    let results = ForceCalculator::new().compute_for_aircraft(
        &aircraft,
        solver.mesh(),
        &solution,
        &condition,
    );
    VtkExporter::export_all(
        &output,
        solver.mesh(),
        &solution,
        &results,
        &condition,
        aircraft.reference_chord(),
        VtkExportOptions::default(),
    )?;
    println!("wrote VTK products to {}", output.display());
    Ok(())
}
