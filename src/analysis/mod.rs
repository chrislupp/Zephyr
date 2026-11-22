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

//! Selective aerodynamic analysis with optional half-wing symmetry.
//!
//! [`AnalysisRequest`] selects both ordinary outputs and exact derivative pairs.
//! The legacy solver and force APIs retain their complete-result behavior.

mod derivatives;
mod loads;
pub(crate) mod operator;
mod request;
mod trefftz;
pub use request::*;

use crate::core::{Scalar, SymmetryType};
use crate::error::{Result, ZephyrError};
use crate::geometry::{Aircraft, MeshInput};
use crate::solver::{FlightCondition, VortexLatticeSolver};
use std::collections::BTreeMap;

#[derive(Clone, Debug)]
pub(crate) enum Source {
    Aircraft(Aircraft),
    Mesh(MeshInput),
}

impl Source {
    fn half(&self) -> Result<Self> {
        let mut result = self.clone();
        match &mut result {
            Self::Mesh(input) => {
                for surface in &mut input.surfaces {
                    if !surface.y_duplicated || surface.y_duplicate_plane != 0.0 {
                        return Err(unsupported(
                            "half-wing mode requires every surface to be duplicated about y = 0",
                        ));
                    }
                    surface.y_duplicated = false;
                }
            }
            Self::Aircraft(aircraft) => {
                if aircraft.y_symmetry() != SymmetryType::None
                    || aircraft.z_symmetry() != SymmetryType::None
                {
                    return Err(unsupported(
                        "half-wing mode does not support additional image planes",
                    ));
                }
                for surface in aircraft.surfaces_mut() {
                    if surface.y_duplicate() != Some(0.0) {
                        return Err(unsupported(
                            "half-wing mode requires every surface to be duplicated about y = 0",
                        ));
                    }
                    surface.clear_y_duplicate();
                }
            }
        }
        Ok(result)
    }
    pub(crate) fn profile_drag(&self) -> Scalar {
        match self {
            Self::Aircraft(a) => a.profile_drag_coefficient(),
            Self::Mesh(_) => 0.0,
        }
    }
    fn solver(&self, settings: Settings) -> Result<VortexLatticeSolver> {
        let mut solver = VortexLatticeSolver::new();
        solver.set_mirror_reuse(settings.mirror_reuse);
        solver.set_symmetric_factorization(settings.symmetric_factorization);
        match self {
            Self::Aircraft(a) => solver.set_aircraft(a)?,
            Self::Mesh(m) => solver.set_mesh(m)?,
        }
        let reference = solver.reference_quantities();
        if ![reference.area, reference.chord, reference.span]
            .iter()
            .all(|v| v.is_finite() && *v > 0.0)
            || !reference
                .moment_reference_point
                .iter()
                .all(|v| v.is_finite())
        {
            return Err(ZephyrError::InvalidGeometry(
                "reference quantities must be finite with positive area, chord and span".into(),
            ));
        }
        if solver.mesh().vortex_count() == 0 {
            return Err(ZephyrError::InvalidMesh("empty lattice".into()));
        }
        Ok(solver)
    }
}

pub(crate) fn unsupported(message: &str) -> ZephyrError {
    ZephyrError::InvalidGeometry(format!(
        "{message}; use full mode for unrestricted existing functionality"
    ))
}

/// Opt-in performance settings, carried across source changes.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Settings {
    mirror_reuse: bool,
    symmetric_factorization: bool,
}

/// Reusable analysis context. Source changes invalidate all caches. Requests can
/// grow or shrink without rebuilding an unchanged circulation system.
#[derive(Clone, Debug)]
pub struct Analysis {
    source: Source,
    settings: Settings,
    active: Option<Source>,
    solver: Option<VortexLatticeSolver>,
    mode: SolverMode,
    condition: Option<FlightCondition>,
    cached_outputs: BTreeMap<Output, Scalar>,
    cached_gradients: Vec<GradientResult>,
    cached_adjoints: BTreeMap<Output, crate::core::VecX>,
    cached_loads: Option<Vec<crate::aero::SurfaceLoads>>,
}

impl Analysis {
    #[must_use]
    pub fn from_mesh(input: &MeshInput) -> Self {
        Self::new(Source::Mesh(input.clone()))
    }
    #[must_use]
    pub fn from_aircraft(aircraft: &Aircraft) -> Self {
        Self::new(Source::Aircraft(aircraft.clone()))
    }
    fn new(source: Source) -> Self {
        Self {
            source,
            settings: Settings::default(),
            active: None,
            solver: None,
            mode: SolverMode::Full,
            condition: None,
            cached_outputs: BTreeMap::new(),
            cached_gradients: Vec::new(),
            cached_adjoints: BTreeMap::new(),
            cached_loads: None,
        }
    }
    /// Replaces the geometry and drops every cache. The opt-in settings below
    /// are kept, as [`VortexLatticeSolver::set_mesh`] keeps its own.
    pub fn set_mesh(&mut self, input: &MeshInput) {
        let settings = self.settings;
        *self = Self::from_mesh(input);
        self.settings = settings;
    }
    /// Aircraft counterpart of [`Self::set_mesh`].
    pub fn set_aircraft(&mut self, aircraft: &Aircraft) {
        let settings = self.settings;
        *self = Self::from_aircraft(aircraft);
        self.settings = settings;
    }

    /// Lets full-mode operator assembly derive the mirrored half of a
    /// `y_duplicated` lattice by reflection instead of evaluating it, as
    /// [`VortexLatticeSolver::set_mirror_reuse`] does for the solver.
    ///
    /// Off by default. Results are bit-identical either way: reflection is
    /// used only where it is proven exact, and a lattice that does not qualify
    /// falls back to evaluating every entry. Both full-mode passes use it, the
    /// control-point pass behind the circulation system and the bound-centre
    /// pass behind near-field loads, roughly halving their kernel evaluations.
    ///
    /// [`SolverMode::HalfWing`] already solves only one semispan, whose
    /// lattice has no mirrored panels, so there the setting has no effect.
    ///
    /// Changing the setting discards the assembled operators, factorization
    /// and cached results; setting the current value again keeps them.
    pub fn set_mirror_reuse(&mut self, enabled: bool) {
        if self.settings.mirror_reuse != enabled {
            self.settings.mirror_reuse = enabled;
            if let Some(solver) = self.solver.as_mut() {
                solver.set_mirror_reuse(enabled);
            }
            self.condition = None;
        }
    }

    #[must_use]
    pub fn mirror_reuse(&self) -> bool {
        self.settings.mirror_reuse
    }

    /// Lets full mode factorize the circulation system of a mirror-symmetric
    /// lattice as two half-size systems, as
    /// [`VortexLatticeSolver::set_symmetric_factorization`] does for the
    /// solver. The primal solve and every adjoint solve then use the block
    /// form.
    ///
    /// Off by default, and with it off results are unchanged bit for bit.
    /// When on, the assembled system is checked for the exact block structure
    /// `[[P, Q], [Q, P]]`; if it has it, `P + Q` and `P - Q` are factorized in
    /// place of the whole matrix at about a quarter of the cost, and results
    /// agree with the full LU to rounding of the elimination (a few `1e-15`
    /// relative on the meshes in the test suite). Antisymmetric deflections
    /// and lattices that are not fully duplicated fail the check and silently
    /// use the full LU.
    ///
    /// In [`SolverMode::HalfWing`] the lattice is one semispan, which has no
    /// mirror pairs, so the full LU is always used and results are those with
    /// the setting off, bit for bit. [`Self::symmetric_factorization_applies`]
    /// reports which factorization the last solve used.
    ///
    /// Changing the setting discards the assembled operators, factorization
    /// and cached results; setting the current value again keeps them.
    pub fn set_symmetric_factorization(&mut self, enabled: bool) {
        if self.settings.symmetric_factorization != enabled {
            self.settings.symmetric_factorization = enabled;
            if let Some(solver) = self.solver.as_mut() {
                solver.set_symmetric_factorization(enabled);
            }
            self.condition = None;
        }
    }

    #[must_use]
    pub fn symmetric_factorization(&self) -> bool {
        self.settings.symmetric_factorization
    }

    /// Whether the system behind the last evaluation was factorized in the
    /// half-size block form of [`Self::set_symmetric_factorization`]. `false`
    /// before the first solve, after a setting change, with the setting off,
    /// in half-wing mode, and whenever the block check fell back to the full
    /// LU.
    #[must_use]
    pub fn symmetric_factorization_applies(&self) -> bool {
        self.condition.is_some()
            && self
                .solver
                .as_ref()
                .is_some_and(VortexLatticeSolver::symmetric_factorization_applies)
    }

    pub fn evaluate(
        &mut self,
        condition: &FlightCondition,
        request: &AnalysisRequest,
    ) -> Result<AnalysisResult> {
        validate_condition(condition)?;
        let plan = request.plan();
        let mut result = AnalysisResult::default();
        if !plan.solve {
            for gradient in &request.gradients {
                if !result.gradients.iter().any(|g| g.request == *gradient) {
                    result.gradients.push(GradientResult {
                        request: gradient.clone(),
                        value: GradientValue::Scalar(0.0),
                    });
                }
            }
            return Ok(result);
        }
        if self.solver.is_none() || self.mode != request.mode {
            let source = if request.mode == SolverMode::HalfWing {
                self.source.half()?
            } else {
                self.source.clone()
            };
            let solver = source.solver(self.settings)?;
            self.active = Some(source);
            self.solver = Some(solver);
            self.mode = request.mode;
            self.condition = None;
        }
        let half = request.mode == SolverMode::HalfWing;
        let solver = self.solver.as_mut().expect("prepared solver");
        // `VortexLatticeSolver::solve` rejects these; the selective solve must
        // not silently ignore a misspelled control instead.
        for (name, _) in &condition.control_deflections {
            if !solver.mesh().control_surface_names().contains(name) {
                return Err(ZephyrError::InvalidFlightCondition(format!(
                    "unknown control '{name}'"
                )));
            }
        }
        if half {
            validate_half(solver, condition)?;
        }
        derivatives::validate(self.active.as_ref().unwrap(), solver, condition, request)?;
        let reuse = self.condition.as_ref() == Some(condition);
        if !reuse {
            self.cached_outputs.clear();
            self.cached_gradients.clear();
            self.cached_adjoints.clear();
            self.cached_loads = None;
        }
        result.work.unknowns = solver.mesh().vortex_count();
        // If evaluation fails, the next call must never reuse its partial state.
        self.condition = None;
        solver.solve_selective(condition, half, plan.velocity, reuse, &mut result.work)?;
        result.work.operator_bytes = std::mem::size_of::<Scalar>()
            * (solver.aic().len()
                + solver.control_velocity_matrix().len()
                + solver.velocity_matrix().len());
        let missing: Vec<_> = request
            .outputs
            .iter()
            .copied()
            .filter(|o| !self.cached_outputs.contains_key(o))
            .collect();
        if !missing.is_empty() || (request.details.surface_loads && self.cached_loads.is_none()) {
            let output_request = AnalysisRequest {
                mode: request.mode,
                outputs: missing,
                details: Details {
                    surface_loads: request.details.surface_loads && self.cached_loads.is_none(),
                    ..Details::default()
                },
                ..AnalysisRequest::default()
            };
            let (values, details) = loads::evaluate(
                solver,
                condition,
                half,
                self.active.as_ref().unwrap().profile_drag(),
                &output_request,
                &mut result.work,
            )?;
            self.cached_outputs.extend(values);
            if details.is_some() {
                self.cached_loads = details;
            }
        }
        let mut missing_gradients = Vec::new();
        for gradient in &request.gradients {
            if !self.cached_gradients.iter().any(|g| g.request == *gradient)
                && !missing_gradients.contains(gradient)
            {
                missing_gradients.push(gradient.clone());
            }
        }
        if !missing_gradients.is_empty() {
            self.cached_gradients.extend(derivatives::evaluate(
                self.active.as_ref().unwrap(),
                solver,
                condition,
                half,
                &missing_gradients,
                &mut result.work,
                &mut self.cached_adjoints,
            )?);
        }
        for output in &request.outputs {
            result.outputs.insert(*output, self.cached_outputs[output]);
        }
        for gradient in &request.gradients {
            if !result.gradients.iter().any(|g| g.request == *gradient) {
                result.gradients.push(
                    self.cached_gradients
                        .iter()
                        .find(|g| g.request == *gradient)
                        .unwrap()
                        .clone(),
                );
            }
        }
        if request.details.circulation {
            result.circulation = Some(loads::expand_circulation(solver, half));
        }
        if request.details.induced_velocity {
            result.induced_velocity = Some(loads::expand_velocity(solver, half));
        }
        if request.details.surface_loads {
            result.surface_loads = self.cached_loads.clone();
        }
        self.condition = Some(condition.clone());
        Ok(result)
    }
}

fn validate_condition(c: &FlightCondition) -> Result<()> {
    if ![
        c.alpha,
        c.beta,
        c.roll_rate,
        c.pitch_rate,
        c.yaw_rate,
        c.mach,
        c.dynamic_pressure,
    ]
    .iter()
    .all(|v| v.is_finite())
        || c.control_deflections.iter().any(|(_, v)| !v.is_finite())
        || c.mach < 0.0
        || c.mach >= 1.0
        || c.dynamic_pressure < 0.0
    {
        return Err(ZephyrError::InvalidFlightCondition(
            "finite subsonic flight conditions and nonnegative pressure are required".into(),
        ));
    }
    Ok(())
}

fn validate_half(solver: &VortexLatticeSolver, c: &FlightCondition) -> Result<()> {
    if c.beta != 0.0 || c.roll_rate != 0.0 || c.yaw_rate != 0.0 {
        return Err(unsupported(
            "half-wing mode requires zero sideslip, roll rate, and yaw rate",
        ));
    }
    for info in solver.mesh().surface_info() {
        let mut side = 0.0_f64;
        for v in &solver.mesh().vortices()
            [info.first_vortex_index..info.first_vortex_index + info.num_vortices]
        {
            for p in [v.bound_vortex_left, v.bound_vortex_right, v.control_point] {
                if p.y != 0.0 {
                    if side != 0.0 && side != p.y.signum() {
                        return Err(unsupported("a half-wing surface cannot cross y = 0"));
                    }
                    side = p.y.signum();
                }
            }
        }
        if side == 0.0 {
            return Err(unsupported(
                "surfaces lying in the symmetry plane are not supported in half-wing mode",
            ));
        }
    }
    for influences in solver.mesh().control_influences() {
        for influence in influences {
            if influence.reflection_sign != 1.0
                && c.control_deflections
                    .iter()
                    .any(|(name, value)| name == &influence.name && *value != 0.0)
            {
                return Err(unsupported(
                    "half-wing mode requires symmetric control deflections",
                ));
            }
        }
    }
    Ok(())
}
