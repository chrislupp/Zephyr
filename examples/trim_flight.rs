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

use zephyr_vlm::core::Vec3;
use zephyr_vlm::geometry::{Aircraft, Section, Surface};
use zephyr_vlm::solver::VortexLatticeSolver;
use zephyr_vlm::trim::{TrimPresets, TrimSolver};

fn main() -> zephyr_vlm::Result<()> {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(5);
    wing.set_spanwise_panel_count(16);
    wing.add_section(Section::new(Vec3::new(0.0, -5.0, 0.0), 1.5, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 5.0, 0.0), 1.5, 0.0));

    let mut tail = Surface::new("Horizontal tail");
    tail.set_chordwise_panel_count(3);
    tail.set_spanwise_panel_count(8);
    tail.add_section(Section::new(Vec3::new(4.0, -1.5, 0.0), 0.6, 0.0));
    tail.add_section(Section::new(Vec3::new(4.0, 1.5, 0.0), 0.4, 0.0));

    let mut aircraft = Aircraft::new("Wing-tail trim example");
    aircraft.set_reference_area(15.0);
    aircraft.set_reference_chord(1.5);
    aircraft.set_reference_span(10.0);
    aircraft.set_moment_reference_point(Vec3::new(0.375, 0.0, 0.0));
    aircraft.add_surface(wing);
    aircraft.add_surface(tail);

    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft)?;
    let result = TrimSolver::new().solve(
        &aircraft,
        &mut solver,
        &TrimPresets::level_flight(0.5),
        None,
    );
    println!("status: {:?}", result.status);
    println!("alpha: {:.4} deg", result.state.alpha.to_degrees());
    println!("CL: {:.6}", result.aerodynamic_results.coefficients.lift);
    println!("residual: {:.3e}", result.residual_norm);
    Ok(())
}
