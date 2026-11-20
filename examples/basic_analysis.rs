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

use zephyr_vlm::aero::ForceCalculator;
use zephyr_vlm::core::{DEGREES_TO_RADIANS, Vec3};
use zephyr_vlm::geometry::{Aircraft, Section, Surface};
use zephyr_vlm::solver::{FlightCondition, VortexLatticeSolver};

fn main() -> zephyr_vlm::Result<()> {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(6);
    wing.set_spanwise_panel_count(20);
    wing.add_section(Section::new(Vec3::new(0.0, -5.0, 0.0), 1.0, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 5.0, 0.0), 1.0, 0.0));

    let mut aircraft = Aircraft::new("Rectangular wing");
    aircraft.set_reference_area(10.0);
    aircraft.set_reference_chord(1.0);
    aircraft.set_reference_span(10.0);
    aircraft.set_moment_reference_point(Vec3::new(0.25, 0.0, 0.0));
    aircraft.add_surface(wing);

    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft)?;
    println!("alpha_deg,CL,CDi,e");
    for alpha_degrees in [-2.0_f64, 0.0, 2.0, 4.0, 6.0, 8.0] {
        let condition = FlightCondition {
            alpha: alpha_degrees * DEGREES_TO_RADIANS,
            ..FlightCondition::default()
        };
        let solution = solver.solve(&condition)?;
        let result = ForceCalculator::new().compute_for_aircraft(
            &aircraft,
            solver.mesh(),
            &solution,
            &condition,
        );
        println!(
            "{alpha_degrees:.1},{:.6},{:.6},{:.4}",
            result.coefficients.lift, result.far_field_induced_drag, result.span_efficiency
        );
    }
    Ok(())
}
