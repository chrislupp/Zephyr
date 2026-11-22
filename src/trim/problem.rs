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

//! Transport-independent definition of a direct-mesh trim problem.

use std::collections::{BTreeMap, BTreeSet};

use crate::core::{MatX, VecX};
use crate::error::{Result, ZephyrError};
use crate::geometry::MeshInput;
use crate::solver::FlightCondition;

use super::{AeroLinearization, AeroQuantity};

/// One independently solved operating variable, in radians.
#[derive(Clone, Debug)]
pub struct TrimVariable {
    /// `alpha`, `beta`, or `control_<name>` (a physical or mixer command).
    pub name: String,
    /// Reproducible initial guess.
    pub initial_guess: f64,
    /// Inclusive lower and upper admissible values.
    pub bounds: [f64; 2],
    /// Positive variable scale for the Newton solve and rank check.
    pub scale: f64,
}

impl TrimVariable {
    /// An unbounded radian variable with unit scaling.
    #[must_use]
    pub fn new(name: impl Into<String>, initial_guess: f64) -> Self {
        Self {
            name: name.into(),
            initial_guess,
            bounds: [f64::NEG_INFINITY, f64::INFINITY],
            scale: 1.0,
        }
    }
}

/// One coefficient target or dimensional resultant balance.
#[derive(Clone, Debug)]
pub struct TrimBalance {
    /// Selected aerodynamic quantity.
    pub quantity: AeroQuantity,
    /// Public physical input carrying this row's target.
    pub target_input: String,
    /// Fixed positive residual scale, in the quantity's units.
    pub scale: f64,
}

impl TrimBalance {
    /// A row with the standard `trim_target_<quantity>` input name.
    #[must_use]
    pub fn new(quantity: AeroQuantity, scale: f64) -> Self {
        Self {
            quantity,
            target_input: format!("trim_target_{}", quantity.name()),
            scale,
        }
    }
}

/// A fixed linkage `delta = offsets + matrix * commands`.
///
/// Command names use `control_<name>` and must all be solved trim variables.
/// Surface names are physical mesh control names without the `control_` prefix.
#[derive(Clone, Debug, Default)]
pub struct ControlMixer {
    /// Independent solved command names, in matrix-column order.
    pub commands: Vec<String>,
    /// Physical control names, in matrix-row order.
    pub surfaces: Vec<String>,
    /// Per-surface weights of the independent commands.
    pub matrix: Vec<Vec<f64>>,
    /// Per-surface offsets in radians.
    pub offsets: Vec<f64>,
    /// Optional limits on realized physical surface commands, in radians.
    pub bounds: BTreeMap<String, [f64; 2]>,
}

/// Numerical controls for internal trim; never physical design inputs.
#[derive(Clone, Debug)]
pub struct MeshTrimOptions {
    /// Maximum accepted Newton steps.
    pub max_iterations: usize,
    /// Per-row absolute tolerance on scaled residuals.
    pub absolute_tolerance: f64,
    /// Relative tolerance referenced to each row's initial scaled residual.
    pub relative_tolerance: f64,
    /// Minimum scaled singular-value ratio admitted for a trim Jacobian.
    pub rank_tolerance: f64,
    /// Maximum backtracking trials per Newton step.
    pub max_backtracks: usize,
}

impl Default for MeshTrimOptions {
    fn default() -> Self {
        Self {
            max_iterations: 30,
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-10,
            rank_tolerance: 1e-12,
            max_backtracks: 25,
        }
    }
}

/// Square residual problem whose locally converged solution is an explicit map.
#[derive(Clone, Debug, Default)]
pub struct MeshTrimProblem {
    /// Independently solved states in Jacobian-column order.
    pub variables: Vec<TrimVariable>,
    /// Balance equations in Jacobian-row order.
    pub balances: Vec<TrimBalance>,
    /// Optional linkage for redundant controls.
    pub mixer: Option<ControlMixer>,
    /// Numerical solver settings.
    pub options: MeshTrimOptions,
}

fn invalid(message: impl Into<String>) -> ZephyrError {
    ZephyrError::InvalidTrim(message.into())
}

impl MeshTrimProblem {
    /// Validate configuration, physical control ownership, and the square system.
    pub fn validate(&self, mesh: &MeshInput) -> Result<()> {
        mesh.validate()?;
        if [
            mesh.reference.area,
            mesh.reference.chord,
            mesh.reference.span,
        ]
        .iter()
        .chain(mesh.reference.moment_reference_point.iter())
        .any(|v| !v.is_finite())
        {
            return Err(invalid("trim reference quantities must be finite"));
        }
        if self.variables.is_empty() || self.variables.len() != self.balances.len() {
            return Err(invalid(
                "one independent variable is required per balance row",
            ));
        }
        let opts = &self.options;
        if opts.max_iterations == 0
            || opts.max_backtracks == 0
            || !opts.absolute_tolerance.is_finite()
            || opts.absolute_tolerance <= 0.0
            || !opts.relative_tolerance.is_finite()
            || opts.relative_tolerance < 0.0
            || !opts.rank_tolerance.is_finite()
            || opts.rank_tolerance <= 0.0
            || opts.rank_tolerance >= 1.0
        {
            return Err(invalid(
                "invalid iteration limits, convergence tolerances, or rank tolerance",
            ));
        }
        let controls: BTreeSet<_> = mesh
            .surfaces
            .iter()
            .flat_map(|s| s.controls.iter().map(|c| c.name.as_str()))
            .collect();
        let mut surface_names = BTreeSet::new();
        for surface in &mesh.surfaces {
            if !surface_names.insert(&surface.name)
                || AeroQuantity::from_name(&surface.name).is_some()
                || [
                    "alpha",
                    "beta",
                    "mach",
                    "dynamic_pressure",
                    "roll_rate",
                    "pitch_rate",
                    "yaw_rate",
                    "reference_area",
                    "reference_chord",
                    "reference_span",
                    "moment_ref_x",
                    "moment_ref_y",
                    "moment_ref_z",
                ]
                .contains(&surface.name.as_str())
                || surface.name.starts_with("control_")
                || surface.name.starts_with("external_")
            {
                return Err(invalid(format!(
                    "mesh input name {} conflicts with another mesh or operating/output variable",
                    surface.name
                )));
            }
        }
        let commands: BTreeSet<_> = self
            .mixer
            .as_ref()
            .map(|m| m.commands.iter().map(String::as_str).collect())
            .unwrap_or_default();
        let mut names = BTreeSet::new();
        for v in &self.variables {
            if !names.insert(v.name.as_str()) {
                return Err(invalid(format!("duplicate variable {}", v.name)));
            }
            let physical = v
                .name
                .strip_prefix("control_")
                .is_some_and(|n| controls.contains(n));
            if !matches!(v.name.as_str(), "alpha" | "beta")
                && !physical
                && !commands.contains(v.name.as_str())
            {
                return Err(invalid(format!("unknown trim variable {}", v.name)));
            }
            if !v.initial_guess.is_finite()
                || !v.scale.is_finite()
                || v.scale <= 0.0
                || v.bounds[0].is_nan()
                || v.bounds[1].is_nan()
                || v.bounds[0] >= v.bounds[1]
                || v.initial_guess < v.bounds[0]
                || v.initial_guess > v.bounds[1]
            {
                return Err(invalid(format!(
                    "invalid guess, bounds, or scale for {}",
                    v.name
                )));
            }
        }
        let mut targets = BTreeSet::new();
        let mut quantities = BTreeSet::new();
        for b in &self.balances {
            if AeroQuantity::from_name(&b.target_input).is_some()
                || [
                    "alpha",
                    "beta",
                    "mach",
                    "dynamic_pressure",
                    "roll_rate",
                    "pitch_rate",
                    "yaw_rate",
                    "reference_area",
                    "reference_chord",
                    "reference_span",
                    "moment_ref_x",
                    "moment_ref_y",
                    "moment_ref_z",
                ]
                .contains(&b.target_input.as_str())
                || b.target_input.starts_with("control_")
                || b.target_input.starts_with("external_")
                || mesh.surfaces.iter().any(|s| s.name == b.target_input)
            {
                return Err(invalid(format!(
                    "target input {} conflicts with an operating, mesh, or load variable",
                    b.target_input
                )));
            }
            if b.target_input.is_empty()
                || !targets.insert(b.target_input.as_str())
                || !quantities.insert(b.quantity as usize)
                || !b.scale.is_finite()
                || b.scale <= 0.0
            {
                return Err(invalid(
                    "balance names and quantities must be unique, with finite positive scales",
                ));
            }
        }
        if let Some(m) = &self.mixer {
            if m.commands.is_empty()
                || m.surfaces.is_empty()
                || commands.len() != m.commands.len()
                || m.matrix.len() != m.surfaces.len()
                || m.offsets.len() != m.surfaces.len()
                || m.matrix
                    .iter()
                    .any(|r| r.len() != m.commands.len() || r.iter().any(|v| !v.is_finite()))
                || m.offsets.iter().any(|v| !v.is_finite())
            {
                return Err(invalid(
                    "invalid mixer dimensions, duplicate commands, or nonfinite weights",
                ));
            }
            for c in &m.commands {
                if c.strip_prefix("control_")
                    .is_some_and(|n| controls.contains(n))
                {
                    return Err(invalid(format!(
                        "independent mixer command {c} conflicts with a physical control"
                    )));
                }
                if !c.starts_with("control_") || c.len() == 8 || !names.contains(c.as_str()) {
                    return Err(invalid(format!(
                        "mixer command {c} must be a solved control variable"
                    )));
                }
            }
            let mut surfaces = BTreeSet::new();
            for s in &m.surfaces {
                if !controls.contains(s.as_str()) || !surfaces.insert(s.as_str()) {
                    return Err(invalid(format!("unknown or duplicate mixer surface {s}")));
                }
                let command = format!("control_{s}");
                if names.contains(command.as_str()) && !commands.contains(command.as_str()) {
                    return Err(invalid(format!(
                        "mixed surface {s} is also independently solved"
                    )));
                }
            }
            for (s, b) in &m.bounds {
                if !surfaces.contains(s.as_str()) || b[0].is_nan() || b[1].is_nan() || b[0] >= b[1]
                {
                    return Err(invalid(format!("invalid mixer bounds for {s}")));
                }
            }
        }
        Ok(())
    }

    /// Whether this variable is supplied by the internal solve or mixer.
    #[must_use]
    pub fn owns(&self, name: &str) -> bool {
        self.variables.iter().any(|v| v.name == name)
            || self
                .mixer
                .as_ref()
                .is_some_and(|m| m.surfaces.iter().any(|s| name == format!("control_{s}")))
    }

    /// Apply trial states and the linkage, preserving prescribed flight data.
    pub fn condition(&self, base: &FlightCondition, z: &VecX) -> Result<FlightCondition> {
        if z.len() != self.variables.len() || z.iter().any(|v| !v.is_finite()) {
            return Err(invalid("invalid trim state vector"));
        }
        let mut condition = base.clone();
        for (v, &value) in self.variables.iter().zip(z.iter()) {
            if value < v.bounds[0] || value > v.bounds[1] {
                return Err(invalid(format!("{} is outside its bounds", v.name)));
            }
            match v.name.as_str() {
                "alpha" => condition.alpha = value,
                "beta" => condition.beta = value,
                _ if self
                    .mixer
                    .as_ref()
                    .is_some_and(|m| m.commands.contains(&v.name)) => {}
                _ => {
                    let name = v
                        .name
                        .strip_prefix("control_")
                        .ok_or_else(|| invalid(format!("unknown trim variable {}", v.name)))?;
                    set_control(&mut condition, name, value);
                }
            }
        }
        if let Some(m) = &self.mixer {
            for (r, surface) in m.surfaces.iter().enumerate() {
                let value = m.offsets[r]
                    + m.commands
                        .iter()
                        .enumerate()
                        .map(|(c, name)| {
                            let index = self
                                .variables
                                .iter()
                                .position(|v| v.name == *name)
                                .expect("validated command");
                            m.matrix[r][c] * z[index]
                        })
                        .sum::<f64>();
                if !value.is_finite()
                    || m.bounds
                        .get(surface)
                        .is_some_and(|b| value < b[0] || value > b[1])
                {
                    return Err(invalid(format!(
                        "mixed control {surface} is outside its bounds"
                    )));
                }
                set_control(&mut condition, surface, value);
            }
        }
        // Stable ordering is part of the normal-rotation model.
        condition.control_deflections.sort_by(|a, b| a.0.cmp(&b.0));
        Ok(condition)
    }

    /// Evaluate scaled residual rows from aerodynamic values and prescribed loads.
    pub fn residual(
        &self,
        values: &[f64; 12],
        targets: &[f64],
        external: &[f64; 6],
    ) -> Result<VecX> {
        if targets.len() != self.balances.len()
            || targets.iter().chain(external).any(|v| !v.is_finite())
        {
            return Err(invalid(
                "nonfinite external loads or missing/nonfinite balance targets",
            ));
        }
        Ok(VecX::from_iterator(
            self.balances.len(),
            self.balances.iter().zip(targets).map(|(b, target)| {
                (values[b.quantity as usize] + b.quantity.load_index().map_or(0.0, |i| external[i])
                    - target)
                    / b.scale
            }),
        ))
    }

    /// Complete aerodynamic output Jacobian with respect to independent states.
    pub fn state_partials(&self, aero: &AeroLinearization) -> MatX {
        MatX::from_fn(12, self.variables.len(), |r, c| {
            let name = &self.variables[c].name;
            if let Some((m, col)) = self
                .mixer
                .as_ref()
                .and_then(|m| m.commands.iter().position(|n| n == name).map(|c| (m, c)))
            {
                m.surfaces
                    .iter()
                    .enumerate()
                    .map(|(i, s)| m.matrix[i][col] * aero.partials[&format!("control_{s}")][(r, 0)])
                    .sum()
            } else {
                aero.partials[name][(r, 0)]
            }
        })
    }

    pub(crate) fn residual_jacobian(&self, state: &MatX) -> MatX {
        MatX::from_fn(self.balances.len(), self.variables.len(), |r, c| {
            state[(self.balances[r].quantity as usize, c)] / self.balances[r].scale
        })
    }
}

fn set_control(condition: &mut FlightCondition, name: &str, value: f64) {
    condition.control_deflections.retain(|(n, _)| n != name);
    condition.control_deflections.push((name.into(), value));
}
