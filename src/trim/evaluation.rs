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

//! Selected aerodynamic quantities and their fixed-operating-state partials.

use std::collections::BTreeMap;

use crate::aero::ForceCalculator;
use crate::core::MatX;
use crate::error::{Result, ZephyrError};
use crate::geometry::MeshInput;
use crate::gradients::{GradientCalculator, OutputSelection};
use crate::solver::{FlightCondition, VortexLatticeSolver};

/// An exposed scalar function of the aerodynamic solution.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(usize)]
pub enum AeroQuantity {
    /// Stability-axis lift coefficient.
    CL,
    /// Near-field drag coefficient.
    CD,
    /// Body side-force coefficient.
    CY,
    /// Conventional rolling-moment coefficient.
    Cl,
    /// Pitching-moment coefficient.
    Cm,
    /// Conventional yawing-moment coefficient.
    Cn,
    /// Dimensional force along the mesh x axis, positive with drag at zero alpha.
    Fx,
    /// Dimensional force along the mesh y axis.
    Fy,
    /// Dimensional force along the mesh z axis, positive with lift at zero alpha.
    Fz,
    /// Mathematical right-hand moment about the mesh x axis.
    Mx,
    /// Mathematical right-hand moment about the mesh y axis.
    My,
    /// Mathematical right-hand moment about the mesh z axis.
    Mz,
}

impl AeroQuantity {
    /// All supported quantities, in coefficient then dimensional-load order.
    pub const ALL: [Self; 12] = [
        Self::CL,
        Self::CD,
        Self::CY,
        Self::Cl,
        Self::Cm,
        Self::Cn,
        Self::Fx,
        Self::Fy,
        Self::Fz,
        Self::Mx,
        Self::My,
        Self::Mz,
    ];

    /// Parse an exact, case-sensitive public quantity name.
    #[must_use]
    pub fn from_name(name: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|q| q.name() == name)
    }

    /// Public variable name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::CL => "CL",
            Self::CD => "CD",
            Self::CY => "CY",
            Self::Cl => "Cl",
            Self::Cm => "Cm",
            Self::Cn => "Cn",
            Self::Fx => "Fx",
            Self::Fy => "Fy",
            Self::Fz => "Fz",
            Self::Mx => "Mx",
            Self::My => "My",
            Self::Mz => "Mz",
        }
    }

    /// Philote-compatible units.
    #[must_use]
    pub const fn units(self) -> &'static str {
        match self {
            Self::Fx | Self::Fy | Self::Fz => "N",
            Self::Mx | Self::My | Self::Mz => "N*m",
            _ => "",
        }
    }

    pub(crate) fn load_index(self) -> Option<usize> {
        (self as usize).checked_sub(6)
    }
}

/// Values and input Jacobian blocks, with rows indexed by [`AeroQuantity`].
///
/// Circulation has already been eliminated. Operating variables are held
/// fixed while differentiating other inputs; these are not trimmed totals.
#[derive(Clone, Debug)]
pub struct AeroLinearization {
    /// All twelve aerodynamic values.
    pub values: [f64; 12],
    /// Each block has twelve rows and one column per scalar input coordinate.
    pub partials: BTreeMap<String, MatX>,
}

/// Evaluate coefficients and dimensional, right-hand body resultants.
///
/// Moments are about `mesh.reference.moment_reference_point`. Zephyr's Cl/Cn
/// sign flips are removed for Mx/Mz; all moment coefficients use reference chord.
#[must_use]
pub fn aerodynamic_values(
    mesh: &MeshInput,
    solver: &VortexLatticeSolver,
    condition: &FlightCondition,
) -> [f64; 12] {
    let result = ForceCalculator::new().compute_with_references(
        &mesh.reference,
        solver.mesh(),
        solver.solution(),
        condition,
    );
    let c = result.coefficients;
    let qs = condition.dynamic_pressure * mesh.reference.area;
    let qsc = qs * mesh.reference.chord;
    [
        c.lift,
        c.drag,
        c.side,
        c.roll,
        c.pitch,
        c.yaw,
        qs * result.body_force.x,
        qs * result.body_force.y,
        qs * result.body_force.z,
        -qsc * c.roll,
        qsc * c.pitch,
        -qsc * c.yaw,
    ]
}

fn selection(quantities: &[AeroQuantity]) -> OutputSelection {
    let mut s = OutputSelection::none();
    for q in quantities {
        match q {
            AeroQuantity::CL => s.lift = true,
            AeroQuantity::CD => s.drag = true,
            AeroQuantity::CY | AeroQuantity::Fy => s.side = true,
            AeroQuantity::Cl | AeroQuantity::Mx => s.roll = true,
            AeroQuantity::Cm | AeroQuantity::My => s.pitch = true,
            AeroQuantity::Cn | AeroQuantity::Mz => s.yaw = true,
            AeroQuantity::Fx | AeroQuantity::Fz => {
                s.lift = true;
                s.drag = true;
            }
        }
    }
    s
}

/// Linearize selected quantities at a matching solved aerodynamic state.
///
/// Set `nodes` to false inside Newton iterations to avoid mesh adjoint work.
/// Unselected coefficient derivative rows are zero. Dimensional quantities
/// requested through `quantities` have their necessary coefficient rows selected.
pub fn linearize_aerodynamics(
    mesh: &MeshInput,
    solver: &VortexLatticeSolver,
    condition: &FlightCondition,
    quantities: &[AeroQuantity],
    nodes: bool,
) -> Result<AeroLinearization> {
    linearize(mesh, solver, condition, quantities, nodes, false)
}

// Trim varies only alpha, beta, and controls. The common longitudinal case
// bypasses prescribed Mach/rate/reference derivative assembly inside Newton.
pub(crate) fn linearize_trim_states(
    mesh: &MeshInput,
    solver: &VortexLatticeSolver,
    condition: &FlightCondition,
    quantities: &[AeroQuantity],
    solve_beta: bool,
) -> Result<AeroLinearization> {
    linearize(mesh, solver, condition, quantities, false, !solve_beta)
}

fn linearize(
    mesh: &MeshInput,
    solver: &VortexLatticeSolver,
    condition: &FlightCondition,
    quantities: &[AeroQuantity],
    nodes: bool,
    alpha_only: bool,
) -> Result<AeroLinearization> {
    solver.check_solved_condition(condition)?;
    if solver.mesh().mesh_input() != Some(mesh) {
        return Err(ZephyrError::InvalidTrim(
            "linearization mesh must match the solved aerodynamic mesh".into(),
        ));
    }
    let values = aerodynamic_values(mesh, solver, condition);
    let mut calculator = GradientCalculator::new();
    calculator.set_output_selection(selection(quantities));
    let flight = calculator.flight_gradients_at_solution(solver, condition, alpha_only)?;
    let mut partials = BTreeMap::new();
    for (name, gradient) in [
        ("alpha", flight.by_alpha),
        ("beta", flight.by_beta),
        ("mach", flight.by_mach),
        ("roll_rate", flight.by_roll_rate),
        ("pitch_rate", flight.by_pitch_rate),
        ("yaw_rate", flight.by_yaw_rate),
        ("moment_ref_x", flight.by_moment_reference[0]),
        ("moment_ref_y", flight.by_moment_reference[1]),
        ("moment_ref_z", flight.by_moment_reference[2]),
    ] {
        let c = gradient.coefficients;
        let v = [c.lift, c.drag, c.side, c.roll, c.pitch, c.yaw];
        partials.insert(
            name.into(),
            MatX::from_fn(12, 1, |r, _| v.get(r).copied().unwrap_or(0.0)),
        );
    }
    for control in calculator.control_gradients_at_solution(solver, condition)? {
        let c = control.by_deflection.coefficients;
        let v = [c.lift, c.drag, c.side, c.roll, c.pitch, c.yaw];
        partials.insert(
            format!("control_{}", control.name),
            MatX::from_fn(12, 1, |r, _| v.get(r).copied().unwrap_or(0.0)),
        );
    }
    if nodes {
        for surface in calculator.node_gradients_at_solution(mesh, solver, condition)? {
            let blocks = [
                &surface.lift,
                &surface.drag,
                &surface.side,
                &surface.roll,
                &surface.pitch,
                &surface.yaw,
            ];
            let n = surface.lift.nodes.len() * 3;
            partials.insert(
                surface.surface_name,
                MatX::from_fn(12, n, |r, c| {
                    if r < 6 {
                        blocks[r].nodes[c / 3][c % 3]
                    } else {
                        0.0
                    }
                }),
            );
        }
    }
    for name in [
        "reference_area",
        "reference_chord",
        "reference_span",
        "dynamic_pressure",
    ] {
        let mut block = MatX::zeros(12, 1);
        for r in 0..6 {
            block[(r, 0)] = match name {
                "reference_area" => -values[r] / mesh.reference.area,
                "reference_chord" if r >= 3 => -values[r] / mesh.reference.chord,
                _ => 0.0,
            };
        }
        partials.insert(name.into(), block);
    }
    let (sin, cos) = condition.alpha.sin_cos();
    let s = mesh.reference.area;
    let chord = mesh.reference.chord;
    let q = condition.dynamic_pressure;
    for (name, block) in &mut partials {
        for col in 0..block.ncols() {
            let d_alpha = f64::from(name == "alpha");
            let d_qs = match name.as_str() {
                "dynamic_pressure" => s,
                "reference_area" => q,
                _ => 0.0,
            };
            let d_qsc = d_qs * chord
                + if name == "reference_chord" {
                    q * s
                } else {
                    0.0
                };
            block[(6, col)] = q
                * s
                * (cos * block[(1, col)]
                    - sin * block[(0, col)]
                    - (sin * values[1] + cos * values[0]) * d_alpha)
                + d_qs * (cos * values[1] - sin * values[0]);
            block[(7, col)] = q * s * block[(2, col)] + d_qs * values[2];
            block[(8, col)] = q
                * s
                * (sin * block[(1, col)]
                    + cos * block[(0, col)]
                    + (cos * values[1] - sin * values[0]) * d_alpha)
                + d_qs * (sin * values[1] + cos * values[0]);
            block[(9, col)] = -q * s * chord * block[(3, col)] - d_qsc * values[3];
            block[(10, col)] = q * s * chord * block[(4, col)] + d_qsc * values[4];
            block[(11, col)] = -q * s * chord * block[(5, col)] - d_qsc * values[5];
        }
    }
    Ok(AeroLinearization { values, partials })
}
