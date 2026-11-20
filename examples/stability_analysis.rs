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

use zephyr_vlm::aero::StabilityDerivativeCalculator;
use zephyr_vlm::core::{DEGREES_TO_RADIANS, Vec3};
use zephyr_vlm::dynamics::{EigenmodeAnalysis, FlightState, MassProperties};
use zephyr_vlm::geometry::{Aircraft, Section, Surface};
use zephyr_vlm::solver::{FlightCondition, VortexLatticeSolver};

fn main() -> zephyr_vlm::Result<()> {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(4);
    wing.set_spanwise_panel_count(12);
    wing.add_section(Section::new(Vec3::new(0.0, -4.0, 0.0), 1.2, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 4.0, 0.0), 1.2, 0.0));
    let mut aircraft = Aircraft::new("Stability example");
    aircraft.set_reference_area(9.6);
    aircraft.set_reference_chord(1.2);
    aircraft.set_reference_span(8.0);
    aircraft.set_moment_reference_point(Vec3::new(0.3, 0.0, 0.0));
    aircraft.add_surface(wing);

    let condition = FlightCondition {
        alpha: 4.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft)?;
    let derivatives = StabilityDerivativeCalculator.compute(&aircraft, &mut solver, &condition)?;
    println!("CL_alpha = {:.5} /rad", derivatives.by_alpha.lift);
    println!("Cm_alpha = {:.5} /rad", derivatives.by_alpha.pitch);
    println!("Cl_p     = {:.5}", derivatives.by_roll_rate.roll);

    let mass = MassProperties {
        mass: 900.0,
        i_xx: 800.0,
        i_yy: 1_400.0,
        i_zz: 1_900.0,
        ..MassProperties::default()
    };
    let flight = FlightState {
        velocity: 50.0,
        ..FlightState::default()
    };
    let modes =
        EigenmodeAnalysis::new().analyze_longitudinal(&aircraft, &mass, &derivatives, &flight)?;
    for mode in modes.modes {
        println!(
            "{}: lambda={:.5}{:+.5}i, zeta={:.3}",
            mode.mode_type.name(),
            mode.eigenvalue.re,
            mode.eigenvalue.im,
            mode.damping_ratio()
        );
    }
    Ok(())
}
