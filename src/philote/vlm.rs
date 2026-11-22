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

//! Node mesh and flight condition to force and moment coefficients, with
//! adjoint partials.
//!
//! Ports the v0.2.0 C++ `VlmDiscipline`. Node partials come from
//! [`GradientCalculator::compute_mesh_node_gradients`], so the cost is one
//! transposed solve per output rather than one primal solve per node
//! coordinate — the reason to run this through an adjoint at all.

use crate::analysis::{
    Analysis, AnalysisRequest, GradientRequest, GradientValue, Input, Output, SolverMode,
};
use std::collections::HashMap;
use std::sync::{Mutex, MutexGuard};

use async_trait::async_trait;
use philote_mdo::{
    ArrayMap, PartialMap, Result, impl_registry,
    registry::VariableRegistry,
    traits::{Discipline, ExplicitDiscipline},
};

use super::aero_evaluation::SolvedState;
use super::{
    flat, invalid, option_flag, option_tokens, option_usize, parse_bool, rejected, scalar,
    scalar_input, sized_input,
};
use crate::aero::ForceCalculator;
use crate::core::Vec3;
use crate::geometry::{MeshControlSurface, MeshInput, SurfaceMesh};
use crate::gradients::GradientCalculator;
use crate::solver::{FlightCondition, VortexLatticeSolver};
use crate::trim::{
    AeroQuantity, MeshTrimProblem, MeshTrimSolver, aerodynamic_values, linearize_aerodynamics,
};

/// Per-surface configuration supplied through discipline options.
#[derive(Clone, Debug)]
pub struct SurfaceOption {
    /// Controls attached to material locations in this surface mesh.
    pub controls: Vec<MeshControlSurface>,
    /// Variable name carrying this surface's nodes.
    pub name: String,
    /// Chordwise panel count.
    pub nx: usize,
    /// Spanwise panel count.
    pub ny: usize,
    /// Whether the surface trails a wake.
    pub has_wake: bool,
    /// Whether the surface is mirrored about a spanwise plane.
    pub y_duplicated: bool,
    /// The mirror plane, when duplicated.
    pub y_duplicate_plane: f64,
    /// Component index, which selects the inter-panel core radius.
    pub component_index: usize,
    /// Whether the freestream acts on this surface.
    pub affected_by_freestream: bool,
}

impl Default for SurfaceOption {
    fn default() -> Self {
        Self {
            controls: Vec::new(),
            name: "Wing".into(),
            nx: 1,
            ny: 1,
            has_wake: true,
            y_duplicated: false,
            y_duplicate_plane: 0.0,
            component_index: 0,
            affected_by_freestream: true,
        }
    }
}

impl SurfaceOption {
    /// Number of node coordinates this surface contributes.
    fn node_values(&self) -> usize {
        (self.nx + 1) * (self.ny + 1) * 3
    }
}

/// The six coefficient outputs, in the order the C++ discipline used.
const OUTPUTS: [&str; 6] = ["CL", "CD", "CY", "Cl", "Cm", "Cn"];

/// Options that change neither the registry nor the trim configuration, so
/// they may still be set after setup when trim is configured.
const RUNTIME_OPTIONS: [&str; 3] = [
    "consistent_forces",
    "mirror_reuse",
    "symmetric_factorization",
];

/// The opt-in solver accelerations, applied to every solver and [`Analysis`]
/// the discipline builds.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct SolverSettings {
    mirror_reuse: bool,
    symmetric_factorization: bool,
}

impl SolverSettings {
    /// These settings with any present in `options` replacing them. Absent
    /// options keep their current value, as `outputs` and `fixed_inputs` do.
    fn updated(self, options: &HashMap<String, serde_json::Value>) -> Result<Self> {
        let mut settings = self;
        if let Some(value) = options.get("mirror_reuse") {
            settings.mirror_reuse = option_flag(value, "mirror_reuse")?;
        }
        if let Some(value) = options.get("symmetric_factorization") {
            settings.symmetric_factorization = option_flag(value, "symmetric_factorization")?;
        }
        Ok(settings)
    }

    /// An empty solver carrying these settings.
    fn solver(self) -> VortexLatticeSolver {
        let mut solver = VortexLatticeSolver::new();
        solver.set_mirror_reuse(self.mirror_reuse);
        solver.set_symmetric_factorization(self.symmetric_factorization);
        solver
    }

    /// A selective analysis of `mesh` carrying these settings.
    fn analysis(self, mesh: &MeshInput) -> Analysis {
        let mut analysis = Analysis::from_mesh(mesh);
        analysis.set_mirror_reuse(self.mirror_reuse);
        analysis.set_symmetric_factorization(self.symmetric_factorization);
        analysis
    }
}

/// Explicit discipline mapping a node mesh and flight condition to loads.
///
/// `outputs` selects coefficient names (comma-separated, default all six).
/// `fixed_inputs` names inputs held constant by the caller's differentiation
/// problem (default none). Their partials are not declared or returned; remove
/// a name before using that input as a design variable. Forward evaluations
/// still honor every supplied value, including fixed inputs.
///
/// One solved state is cached using exact equality of all input arrays. Option
/// changes clear the cache. Function and derivative requests share that state;
/// no approximate matching or numerical differentiation is used.
///
/// `trim` optionally configures an internal square root solve. Selected operating
/// variables become outputs; coefficient and dimensional load partials include
/// their converged trim response. Configure trim before setup. `Fx,Fy,Fz` are
/// dimensional mesh-axis forces; `Mx,My,Mz` are right-hand mesh-axis moments.
///
/// `mirror_reuse` and `symmetric_factorization` (booleans, both default
/// `false`) switch on [`VortexLatticeSolver::set_mirror_reuse`] and
/// [`VortexLatticeSolver::set_symmetric_factorization`] for every solve the
/// discipline makes: the ordinary solved state, the internal trim iteration,
/// and the selective `analysis_*` [`Analysis`]. Mirror reuse leaves results
/// bit-identical; the symmetric factorization agrees with the full LU to
/// rounding where it applies and is bit-identical where it falls back. They
/// combine with every other option and, unlike the trim and registry options,
/// may still be changed after setup when trim is configured. Changing either
/// clears the cached state; [`Self::symmetric_factorization_applies`] reports
/// whether the last solve used the half-size factorization.
#[derive(Debug)]
pub struct VlmDiscipline {
    registry: VariableRegistry,
    surfaces: Vec<SurfaceOption>,
    settings: SolverSettings,
    analysis_request: Option<AnalysisRequest>,
    analysis_cache: Mutex<Option<(MeshInput, Analysis)>>,
    outputs: Vec<String>,
    fixed_inputs: Vec<String>,
    trim: Option<MeshTrimProblem>,
    setup_complete: bool,
    cache: Mutex<Option<SolvedState>>,
}

impl Default for VlmDiscipline {
    fn default() -> Self {
        Self {
            registry: VariableRegistry::new(false),
            surfaces: vec![SurfaceOption::default()],
            settings: SolverSettings::default(),
            analysis_request: None,
            analysis_cache: Mutex::new(None),
            outputs: OUTPUTS.iter().map(|s| (*s).into()).collect(),
            fixed_inputs: Vec::new(),
            trim: None,
            setup_complete: false,
            cache: Mutex::new(None),
        }
    }
}

impl VlmDiscipline {
    /// A discipline configured for the given surfaces.
    #[must_use]
    pub fn new(surfaces: Vec<SurfaceOption>) -> Self {
        Self {
            registry: VariableRegistry::new(false),
            surfaces,
            settings: SolverSettings::default(),
            analysis_request: None,
            analysis_cache: Mutex::new(None),
            outputs: OUTPUTS.iter().map(|s| (*s).into()).collect(),
            fixed_inputs: Vec::new(),
            trim: None,
            setup_complete: false,
            cache: Mutex::new(None),
        }
    }

    fn solved(&self, inputs: &ArrayMap) -> Result<MutexGuard<'_, Option<SolvedState>>> {
        let mut cache = self
            .cache
            .lock()
            .map_err(|_| rejected("cache", "poisoned solution cache"))?;
        if cache.as_ref().is_none_or(|state| state.inputs != *inputs) {
            *cache = None;
            if self.trim.is_some() {
                for name in self.physical_inputs() {
                    if !self.surfaces.iter().any(|s| s.name == name) {
                        let value = sized_input(inputs, &name, 1)?[0];
                        if !value.is_finite() {
                            return Err(rejected("trim_input", format!("{name} must be finite")));
                        }
                    }
                }
            }
            let mesh = self.build_mesh(inputs)?;
            let condition = self.build_condition(inputs)?;
            let mut solver = self.settings.solver();
            let trim = if let Some(problem) = &self.trim {
                Some(
                    MeshTrimSolver
                        .solve_with_solver(
                            self.settings.solver(),
                            &mesh,
                            &condition,
                            problem,
                            &super::trim::targets(inputs, problem)?,
                            &super::trim::external(inputs, problem)?,
                        )
                        .map_err(|e| rejected("trim", e.to_string()))?,
                )
            } else {
                solver
                    .set_mesh(&mesh)
                    .map_err(|e| rejected("set_mesh", e.to_string()))?;
                solver
                    .solve(&condition)
                    .map_err(|e| rejected("solve", e.to_string()))?;
                None
            };
            let condition = trim.as_ref().map_or(condition, |t| t.condition.clone());
            *cache = Some(SolvedState {
                inputs: inputs.clone(),
                mesh,
                condition,
                solver,
                trim,
                partials: None,
            });
        }
        Ok(cache)
    }

    fn active(&self, input: &str) -> bool {
        !self.fixed_inputs.iter().any(|name| name == input)
    }

    fn owned(&self, name: &str) -> bool {
        self.trim.as_ref().is_some_and(|p| p.owns(name))
    }

    fn output_names(&self) -> Vec<String> {
        let mut names = self.outputs.clone();
        if let Some(problem) = &self.trim {
            names.extend(problem.variables.iter().map(|v| v.name.clone()));
            if let Some(mixer) = &problem.mixer {
                names.extend(mixer.surfaces.iter().map(|s| format!("control_{s}")));
            }
        }
        names.sort();
        names.dedup();
        names
    }

    fn physical_inputs(&self) -> Vec<String> {
        let mut names: Vec<_> = self.surfaces.iter().map(|s| s.name.clone()).collect();
        names.extend(
            [
                "reference_area",
                "reference_chord",
                "reference_span",
                "moment_ref_x",
                "moment_ref_y",
                "moment_ref_z",
                "dynamic_pressure",
                "alpha",
                "beta",
                "mach",
                "roll_rate",
                "pitch_rate",
                "yaw_rate",
            ]
            .map(String::from),
        );
        names.extend(
            self.control_names()
                .into_iter()
                .map(|n| format!("control_{n}")),
        );
        names.retain(|n| !self.owned(n));
        if let Some(problem) = &self.trim {
            names.extend(problem.balances.iter().map(|b| b.target_input.clone()));
            if problem.balances.iter().any(|b| b.quantity as usize >= 6) {
                names.extend(
                    AeroQuantity::ALL[6..]
                        .iter()
                        .map(|q| format!("external_{}", q.name())),
                );
            }
        }
        names
    }

    /// Diagnostics from the cached successful trim, if trim is enabled and solved.
    pub fn trim_diagnostics(&self) -> Result<Option<crate::trim::TrimDiagnostics>> {
        let cache = self
            .cache
            .lock()
            .map_err(|_| rejected("cache", "poisoned solution cache"))?;
        Ok(cache
            .as_ref()
            .and_then(|s| s.trim.as_ref())
            .map(|t| t.diagnostics.clone()))
    }

    /// Whether the cached solve factorized its system in the half-size block
    /// form of the `symmetric_factorization` option.
    ///
    /// `None` when nothing is cached: before the first evaluation, after an
    /// option change, or after a failed solve. `Some(false)` when the option is
    /// off or the system fell back to the full LU (an antisymmetric control
    /// deflection, a surface that is not y-duplicated, half-wing analysis
    /// mode). With trim, this describes the converged trim solve; with the
    /// `analysis_*` options, the cached selective analysis's last solve.
    pub fn symmetric_factorization_applies(&self) -> Result<Option<bool>> {
        if self.analysis_request.is_some() {
            let cache = self
                .analysis_cache
                .lock()
                .map_err(|_| invalid("analysis cache is poisoned"))?;
            return Ok(cache
                .as_ref()
                .map(|(_, analysis)| analysis.symmetric_factorization_applies()));
        }
        let cache = self
            .cache
            .lock()
            .map_err(|_| rejected("cache", "poisoned solution cache"))?;
        Ok(cache
            .as_ref()
            .map(|state| state.solver().symmetric_factorization_applies()))
    }

    /// The configured surfaces.
    #[must_use]
    pub fn surfaces(&self) -> &[SurfaceOption] {
        &self.surfaces
    }

    fn control_names(&self) -> Vec<String> {
        self.surfaces
            .iter()
            .flat_map(|s| s.controls.iter().map(|c| c.name.clone()))
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .collect()
    }

    /// Assemble a [`MeshInput`] from the incoming node arrays.
    fn build_mesh(&self, inputs: &ArrayMap) -> Result<MeshInput> {
        let mut mesh = MeshInput::default();
        for surface in &self.surfaces {
            let nodes = sized_input(inputs, &surface.name, surface.node_values())?;
            let mut generated = SurfaceMesh::new(&surface.name, surface.nx, surface.ny);
            generated.controls = surface.controls.clone();
            generated.component_index = surface.component_index;
            generated.affected_by_freestream = surface.affected_by_freestream;
            generated.has_wake = surface.has_wake;
            generated.y_duplicated = surface.y_duplicated;
            generated.y_duplicate_plane = surface.y_duplicate_plane;
            for (index, chunk) in nodes.chunks_exact(3).enumerate() {
                generated.nodes[index] = Vec3::new(chunk[0], chunk[1], chunk[2]);
            }
            mesh.surfaces.push(generated);
        }
        mesh.reference.area = scalar_input(inputs, "reference_area")?;
        mesh.reference.chord = scalar_input(inputs, "reference_chord")?;
        mesh.reference.span = scalar_input(inputs, "reference_span")?;
        mesh.reference.moment_reference_point = Vec3::new(
            scalar_input(inputs, "moment_ref_x")?,
            scalar_input(inputs, "moment_ref_y")?,
            scalar_input(inputs, "moment_ref_z")?,
        );
        Ok(mesh)
    }

    /// Assemble a [`FlightCondition`] from the incoming scalars.
    fn build_condition(&self, inputs: &ArrayMap) -> Result<FlightCondition> {
        let value = |name: &str| -> Result<f64> {
            if self.owned(name) {
                Ok(0.0)
            } else {
                scalar_input(inputs, name)
            }
        };
        Ok(FlightCondition {
            alpha: value("alpha")?,
            beta: value("beta")?,
            mach: scalar_input(inputs, "mach")?,
            dynamic_pressure: scalar_input(inputs, "dynamic_pressure")?,
            roll_rate: scalar_input(inputs, "roll_rate")?,
            pitch_rate: scalar_input(inputs, "pitch_rate")?,
            yaw_rate: scalar_input(inputs, "yaw_rate")?,
            control_deflections: self
                .control_names()
                .into_iter()
                .map(|name| {
                    let command = value(&format!("control_{name}"))?;
                    Ok((name, command))
                })
                .collect::<Result<_>>()?,
        })
    }
}

/// Read one per-surface option list, padding shorter lists with the default.
fn apply_surface_option<T>(
    options: &HashMap<String, serde_json::Value>,
    key: &str,
    surfaces: &mut [SurfaceOption],
    mut parse: impl FnMut(&str) -> Option<T>,
    mut assign: impl FnMut(&mut SurfaceOption, T),
) -> Result<()> {
    let Some(value) = options.get(key) else {
        return Ok(());
    };
    let tokens = option_tokens(value, key)?;
    for (surface, token) in surfaces.iter_mut().zip(&tokens) {
        let parsed = parse(token)
            .ok_or_else(|| invalid(format!("option `{key}` has a bad value `{token}`")))?;
        assign(surface, parsed);
    }
    Ok(())
}

impl Discipline for VlmDiscipline {
    impl_registry!(registry);

    fn name(&self) -> &str {
        "ZephyrVLM"
    }

    fn is_differentiable(&self) -> bool {
        true
    }

    fn provides_gradients(&self) -> bool {
        true
    }

    fn initialize(&mut self) -> Result<()> {
        self.add_option("num_surfaces", "int")?;
        self.add_option("analysis_outputs", "str")?;
        self.add_option("analysis_gradients", "dict")?;
        self.add_option("analysis_mode", "str")?;
        self.add_option("consistent_forces", "bool")?;
        self.add_option("mirror_reuse", "bool")?;
        self.add_option("symmetric_factorization", "bool")?;
        self.add_option("outputs", "str")?;
        self.add_option("fixed_inputs", "str")?;
        self.add_option("controls", "str")?;
        self.add_option("trim", "str")?;
        for key in [
            "surface_names",
            "surface_nx",
            "surface_ny",
            "surface_has_wake",
            "surface_y_duplicated",
            "surface_y_duplicate_plane",
            "surface_component_index",
            "surface_affected_by_freestream",
        ] {
            // philote-mdo's option types are bool/int/float/str/dict; the
            // per-surface lists arrive as `str` even when written as JSON
            // arrays, since the array form is decoded from the same field.
            self.add_option(key, "str")?;
        }
        Ok(())
    }

    fn set_options(&mut self, options: &HashMap<String, serde_json::Value>) -> Result<()> {
        // Forces always act at bound-vortex centers. The option stays declared
        // so clients that pass `true` keep working; `false` is refused rather
        // than silently ignored.
        if let Some(value) = options.get("consistent_forces") {
            require_consistent_forces(value)?;
        }
        let settings = self.settings.updated(options)?;
        if self.setup_complete && self.trim.is_some() {
            if options
                .keys()
                .any(|name| !RUNTIME_OPTIONS.contains(&name.as_str()))
            {
                return Err(invalid(
                    "trim configuration and registry options must be set before setup; create a new discipline to reconfigure",
                ));
            }
            // Only the solver settings can change here. Keep the cached trim
            // unless they did, as a repeated `consistent_forces` keeps it.
            if settings != self.settings {
                self.settings = settings;
                *self
                    .cache
                    .get_mut()
                    .map_err(|_| rejected("cache", "poisoned solution cache"))? = None;
            }
            return Ok(());
        }
        self.settings = settings;
        *self
            .cache
            .get_mut()
            .map_err(|_| rejected("cache", "poisoned solution cache"))? = None;
        if options.contains_key("trim") {
            if self.setup_complete {
                return Err(invalid(
                    "configure trim before setup; create a new discipline to rebuild variable ownership",
                ));
            }
            self.trim = super::trim::selected(options)?;
        }
        if let Some(value) = options.get("outputs") {
            let outputs = option_tokens(value, "outputs")?;
            if outputs.is_empty() || outputs.iter().any(|s| AeroQuantity::from_name(s).is_none()) {
                return Err(invalid(
                    "outputs must select coefficients or Fx,Fy,Fz,Mx,My,Mz",
                ));
            }
            self.outputs = outputs;
        }
        if let Some(value) = options.get("fixed_inputs") {
            self.fixed_inputs = option_tokens(value, "fixed_inputs")?;
        }
        let count = match options.get("num_surfaces") {
            Some(value) => option_usize(value, "num_surfaces")?,
            None => self.surfaces.len(),
        };
        if count == 0 {
            return Err(invalid("`num_surfaces` must be at least 1"));
        }

        // Names are required whenever more than the default single surface is
        // requested: they are the input variable names, so a wrong count here
        // produces a discipline whose inputs the client cannot address.
        let names = match options.get("surface_names") {
            Some(value) => option_tokens(value, "surface_names")?,
            None => (0..count).map(|index| format!("surface_{index}")).collect(),
        };
        if names.len() != count {
            return Err(invalid(format!(
                "`surface_names` has {} entries, expected {count}",
                names.len()
            )));
        }

        self.surfaces = names
            .into_iter()
            .map(|name| SurfaceOption {
                name,
                ..SurfaceOption::default()
            })
            .collect();

        for (key, assign) in [
            ("surface_nx", 0usize),
            ("surface_ny", 1),
            ("surface_component_index", 2),
        ] {
            apply_surface_option(
                options,
                key,
                &mut self.surfaces,
                |token| token.parse::<usize>().ok(),
                |surface, value| match assign {
                    0 => surface.nx = value,
                    1 => surface.ny = value,
                    _ => surface.component_index = value,
                },
            )?;
        }
        apply_surface_option(
            options,
            "surface_y_duplicate_plane",
            &mut self.surfaces,
            |token| token.parse::<f64>().ok(),
            |surface, value| surface.y_duplicate_plane = value,
        )?;
        for (key, slot) in [
            ("surface_has_wake", 0usize),
            ("surface_y_duplicated", 1),
            ("surface_affected_by_freestream", 2),
        ] {
            apply_surface_option(
                options,
                key,
                &mut self.surfaces,
                |token| Some(parse_bool(token)),
                |surface, value| match slot {
                    0 => surface.has_wake = value,
                    1 => surface.y_duplicated = value,
                    _ => surface.affected_by_freestream = value,
                },
            )?;
        }
        if let Some(value) = options.get("controls") {
            let decoded = if let Some(text) = value.as_str() {
                serde_json::from_str(text)
                    .map_err(|e| invalid(format!("invalid controls JSON: {e}")))?
            } else {
                value.clone()
            };
            let definitions = decoded
                .as_array()
                .ok_or_else(|| invalid("controls must be a JSON array"))?;
            for definition in definitions {
                let object = definition
                    .as_object()
                    .ok_or_else(|| invalid("control must be an object"))?;
                for key in object.keys() {
                    if ![
                        "surface",
                        "name",
                        "span_start",
                        "span_end",
                        "hinge",
                        "gain",
                        "duplicate_sign",
                        "leading_edge",
                    ]
                    .contains(&key.as_str())
                    {
                        return Err(invalid(format!("unknown control field '{key}'")));
                    }
                }
                let text = |key| {
                    object
                        .get(key)
                        .and_then(|v| v.as_str())
                        .ok_or_else(|| invalid(format!("control requires string '{key}'")))
                };
                let surface_name = text("surface")?;
                let surface = self
                    .surfaces
                    .iter_mut()
                    .find(|s| s.name == surface_name)
                    .ok_or_else(|| invalid(format!("unknown control surface '{surface_name}'")))?;
                let index = |key: &str| -> Result<usize> {
                    let value = object
                        .get(key)
                        .and_then(|v| v.as_u64())
                        .ok_or_else(|| invalid(format!("control requires integer '{key}'")))?;
                    usize::try_from(value).map_err(|_| invalid("control index is too large"))
                };
                let hinge = object
                    .get("hinge")
                    .ok_or_else(|| invalid("control requires hinge"))?;
                let hinge = if let Some(value) = hinge.as_f64() {
                    vec![value; surface.ny + 1]
                } else {
                    hinge
                        .as_array()
                        .ok_or_else(|| invalid("hinge must be a number or array"))?
                        .iter()
                        .map(|v| {
                            v.as_f64()
                                .ok_or_else(|| invalid("hinge must contain numbers"))
                        })
                        .collect::<Result<Vec<_>>>()?
                };
                let mut control = MeshControlSurface::trailing_edge(
                    text("name")?,
                    index("span_start")?..index("span_end")?,
                    hinge,
                );
                for (key, target) in [
                    ("gain", &mut control.gain),
                    ("duplicate_sign", &mut control.duplicate_sign),
                ] {
                    if let Some(value) = object.get(key) {
                        *target = value
                            .as_f64()
                            .ok_or_else(|| invalid(format!("{key} must be numeric")))?;
                    }
                }
                if let Some(value) = object.get("leading_edge") {
                    control.leading_edge = value
                        .as_bool()
                        .ok_or_else(|| invalid("leading_edge must be boolean"))?;
                }
                control
                    .validate_layout(surface.nx, surface.ny)
                    .map_err(|e| invalid(e.to_string()))?;
                surface.controls.push(control);
            }
        }
        if ["analysis_outputs", "analysis_gradients", "analysis_mode"]
            .iter()
            .any(|key| options.contains_key(*key))
        {
            let mut request = self
                .analysis_request
                .clone()
                .unwrap_or_else(|| AnalysisRequest {
                    outputs: vec![
                        Output::Lift,
                        Output::Drag,
                        Output::Side,
                        Output::Roll,
                        Output::Pitch,
                        Output::Yaw,
                    ],
                    ..AnalysisRequest::default()
                });
            if let Some(value) = options.get("analysis_mode") {
                request.mode = match value.as_str() {
                    Some("full") => SolverMode::Full,
                    Some("half_wing") => SolverMode::HalfWing,
                    _ => return Err(invalid("analysis_mode must be 'full' or 'half_wing'")),
                };
            }
            if let Some(value) = options.get("analysis_outputs") {
                request.outputs = option_tokens(value, "analysis_outputs")?
                    .iter()
                    .map(|name| analysis_output(name))
                    .collect::<Result<_>>()?;
            }
            if let Some(value) = options.get("analysis_gradients") {
                let entries = value.as_object().ok_or_else(|| {
                    invalid("analysis_gradients must map output names to input-name lists")
                })?;
                request.gradients.clear();
                for (name, inputs) in entries {
                    let output = analysis_output(name)?;
                    for input in option_tokens(inputs, "analysis_gradients")? {
                        request.gradients.push(GradientRequest {
                            output,
                            input: self.analysis_input(&input)?,
                        });
                    }
                }
            }
            self.analysis_request = Some(request);
        }
        *self
            .analysis_cache
            .get_mut()
            .map_err(|_| invalid("analysis cache is poisoned"))? = None;
        Ok(())
    }

    fn setup(&mut self) -> Result<()> {
        if self.setup_complete && !self.registry.var_meta().is_empty() {
            return Err(invalid("discipline is already set up"));
        }
        let mut layout = MeshInput::default();
        for surface in &self.surfaces {
            let mut mesh = SurfaceMesh::new(&surface.name, surface.nx, surface.ny);
            for i in 0..=surface.nx {
                for j in 0..=surface.ny {
                    mesh.set_node(i, j, Vec3::new(i as f64, j as f64, 0.0))
                        .map_err(|e| invalid(e.to_string()))?;
                }
            }
            mesh.controls = surface.controls.clone();
            for control in &mesh.controls {
                control
                    .validate_layout(surface.nx, surface.ny)
                    .map_err(|e| invalid(e.to_string()))?;
            }
            layout.surfaces.push(mesh);
        }
        if let Some(problem) = &self.trim {
            problem
                .validate(&layout)
                .map_err(|e| invalid(e.to_string()))?;
        }
        if self.analysis_request.is_some()
            && (self.trim.is_some()
                || !self.fixed_inputs.is_empty()
                || self.outputs.iter().map(String::as_str).ne(OUTPUTS)
                || self.surfaces.iter().any(|s| !s.controls.is_empty()))
        {
            // The selective path has no control inputs and its geometry
            // adjoint requires an empty deflection list, so declared controls
            // would make every mesh-node partial fail at evaluation time.
            return Err(invalid(
                "analysis_outputs, analysis_gradients and analysis_mode select their own outputs \
                 and partials; they cannot be combined with trim, outputs, fixed_inputs or \
                 controls",
            ));
        }
        let inputs = self.physical_inputs();
        let outputs = match &self.analysis_request {
            Some(request) => {
                let mut outputs = request.outputs.clone();
                // A gradient output must exist in the discipline registry, even if
                // its value is only an internal dependency of the core request.
                outputs.extend(request.gradients.iter().map(|g| g.output));
                outputs.sort();
                outputs.dedup();
                outputs
                    .into_iter()
                    .map(|output| output_name(output).to_owned())
                    .collect()
            }
            None => self.output_names(),
        };
        let mut unique = std::collections::BTreeSet::new();
        for name in inputs.iter().chain(&outputs) {
            if !unique.insert(name) {
                return Err(invalid(format!("input/output name collision '{name}'")));
            }
        }
        for input in &self.fixed_inputs {
            if !inputs.contains(input) {
                return Err(invalid(format!(
                    "unknown or internally solved fixed input {input}"
                )));
            }
        }
        for name in inputs {
            if let Some(surface) = self.surfaces.iter().find(|s| s.name == name) {
                self.add_input(&name, &[surface.node_values()], "m")?;
                continue;
            }
            let units = match name.as_str() {
                "reference_area" => "m**2",
                "reference_chord" | "reference_span" | "moment_ref_x" | "moment_ref_y"
                | "moment_ref_z" => "m",
                "alpha" | "beta" => "rad",
                "dynamic_pressure" => "Pa",
                n if n.starts_with("control_") => "rad",
                _ => self
                    .trim
                    .as_ref()
                    .and_then(|p| p.balances.iter().find(|b| b.target_input == name))
                    .map(|b| b.quantity.units())
                    .or_else(|| {
                        name.strip_prefix("external_")
                            .and_then(AeroQuantity::from_name)
                            .map(AeroQuantity::units)
                    })
                    .unwrap_or(""),
            };
            self.add_input(&name, &[1], units)?;
        }
        for name in outputs {
            let units = if self.analysis_request.is_some() {
                ""
            } else {
                AeroQuantity::from_name(&name).map_or("rad", AeroQuantity::units)
            };
            self.add_output(&name, &[1], units)?;
        }
        self.setup_complete = true;
        Ok(())
    }

    fn setup_partials(&mut self) -> Result<()> {
        if let Some(request) = &self.analysis_request {
            let pairs: Vec<_> = request
                .gradients
                .iter()
                .map(|g| (output_name(g.output), self.input_name(&g.input)))
                .collect();
            for (output, input) in pairs {
                self.declare_partials(output, &input)?;
            }
            return Ok(());
        }
        for output in self.output_names() {
            for input in self.physical_inputs() {
                if self.active(&input) {
                    self.declare_partials(&output, &input)?;
                }
            }
        }
        Ok(())
    }
}

#[async_trait]
impl ExplicitDiscipline for VlmDiscipline {
    async fn compute(&self, inputs: &ArrayMap) -> Result<ArrayMap> {
        if let Some(request) = &self.analysis_request {
            let mut values_request = request.clone();
            values_request
                .outputs
                .extend(request.gradients.iter().map(|g| g.output));
            values_request.gradients.clear();
            let result = self.evaluate_selected(inputs, &values_request)?;
            return Ok(result
                .outputs
                .into_iter()
                .map(|(o, v)| (output_name(o).to_owned(), scalar(v)))
                .collect());
        }
        let cache = self.solved(inputs)?;
        let state = cache.as_ref().expect("solved state");
        let values = state.trim.as_ref().map_or_else(
            || aerodynamic_values(&state.mesh, state.solver(), &state.condition),
            |trim| trim.values,
        );
        let mut outputs = HashMap::new();
        for name in &self.outputs {
            let q = AeroQuantity::from_name(name).expect("validated output");
            outputs.insert(name.clone(), scalar(values[q as usize]));
        }
        if let Some(trim) = &state.trim {
            for (variable, value) in trim.problem.variables.iter().zip(trim.states.iter()) {
                outputs.insert(variable.name.clone(), scalar(*value));
            }
            if let Some(mixer) = &trim.problem.mixer {
                for surface in &mixer.surfaces {
                    let value = trim
                        .condition
                        .control_deflections
                        .iter()
                        .find(|(n, _)| n == surface)
                        .expect("realized control")
                        .1;
                    outputs.insert(format!("control_{surface}"), scalar(value));
                }
            }
        }
        Ok(outputs)
    }

    async fn compute_partials(&self, inputs: &ArrayMap) -> Result<PartialMap> {
        if let Some(request) = &self.analysis_request {
            let mut derivative_request = request.clone();
            derivative_request.outputs.clear();
            let result = self.evaluate_selected(inputs, &derivative_request)?;
            return Ok(result
                .gradients
                .into_iter()
                .map(|gradient| {
                    let key = (
                        output_name(gradient.request.output).to_owned(),
                        self.input_name(&gradient.request.input),
                    );
                    let values = match gradient.value {
                        GradientValue::Scalar(v) => vec![v],
                        GradientValue::MeshNodes(nodes) => {
                            nodes.iter().flat_map(|node| node.iter().copied()).collect()
                        }
                    };
                    (key, flat(values))
                })
                .collect());
        }
        if self.trim.is_some()
            || self
                .outputs
                .iter()
                .any(|n| AeroQuantity::from_name(n).is_some_and(|q| q as usize >= 6))
        {
            let mut cache = self.solved(inputs)?;
            let state = cache.as_mut().expect("solved state");
            if let Some(partials) = &state.partials {
                return Ok(partials.clone());
            }
            let quantities: Vec<_> = self
                .outputs
                .iter()
                .map(|n| AeroQuantity::from_name(n).expect("validated output"))
                .collect();
            let nodes = self.surfaces.iter().any(|s| self.active(&s.name));
            let partials: PartialMap = if let Some(trim) = &state.trim {
                trim.compute_partials(&state.mesh, &quantities, nodes)
                    .map_err(|e| rejected("trim_partials", e.to_string()))?
                    .into_iter()
                    .filter(|((_, input), _)| self.active(input))
                    .map(|(key, values)| (key, flat(values)))
                    .collect()
            } else {
                let aero = linearize_aerodynamics(
                    &state.mesh,
                    state.solver(),
                    &state.condition,
                    &quantities,
                    nodes,
                )
                .map_err(|e| rejected("load_partials", e.to_string()))?;
                let mut partials = HashMap::new();
                for (input, block) in aero.partials {
                    if self.active(&input) {
                        for output in &quantities {
                            let values = (0..block.ncols())
                                .map(|c| block[(*output as usize, c)])
                                .collect();
                            partials.insert((output.name().into(), input.clone()), flat(values));
                        }
                    }
                }
                partials
            };
            state.partials = Some(partials.clone());
            return Ok(partials);
        }
        let cache = self.solved(inputs)?;
        let state = cache.as_ref().expect("solved state");
        let (mesh, condition, solver) = (&state.mesh, &state.condition, state.solver());
        let solution = solver.solution();
        let results = ForceCalculator::new().compute_with_references(
            &mesh.reference,
            solver.mesh(),
            solution,
            condition,
        );

        // Both gradient sets run against the solution above rather than
        // solving again: the factorization is the expensive part of an
        // operating point, and one linearize used to pay for it three times.
        let mut calculator = GradientCalculator::new();
        calculator.set_output_selection(crate::gradients::OutputSelection {
            induced_drag: false,
            lift: self.outputs.iter().any(|s| s == "CL"),
            drag: self.outputs.iter().any(|s| s == "CD"),
            side: self.outputs.iter().any(|s| s == "CY"),
            roll: self.outputs.iter().any(|s| s == "Cl"),
            pitch: self.outputs.iter().any(|s| s == "Cm"),
            yaw: self.outputs.iter().any(|s| s == "Cn"),
        });
        let node_gradients = calculator
            .node_gradients_at_solution(mesh, solver, condition)
            .map_err(|error| rejected("compute_mesh_node_gradients", error.to_string()))?;
        let alpha_only = [
            "beta",
            "mach",
            "roll_rate",
            "pitch_rate",
            "yaw_rate",
            "moment_ref_x",
            "moment_ref_y",
            "moment_ref_z",
        ]
        .iter()
        .all(|name| !self.active(name));
        let flight = calculator
            .flight_gradients_at_solution(solver, condition, alpha_only)
            .map_err(|error| rejected("compute_flight_condition_gradients", error.to_string()))?;

        let mut partials: PartialMap = HashMap::new();
        let mut insert = |output: &str, input: &str, values: Vec<f64>| {
            if self.outputs.iter().any(|s| s == output) && self.active(input) {
                partials.insert((output.to_string(), input.to_string()), flat(values));
            }
        };

        // Node partials, surface by surface.
        for (surface, gradients) in self.surfaces.iter().zip(&node_gradients) {
            for (output, block) in [
                ("CL", &gradients.lift),
                ("CD", &gradients.drag),
                ("CY", &gradients.side),
                ("Cl", &gradients.roll),
                ("Cm", &gradients.pitch),
                ("Cn", &gradients.yaw),
            ] {
                let flattened: Vec<f64> = block
                    .nodes
                    .iter()
                    .flat_map(|node| [node.x, node.y, node.z])
                    .collect();
                insert(output, &surface.name, flattened);
            }
        }

        // Flight-condition partials.
        let columns = [
            ("alpha", &flight.by_alpha),
            ("beta", &flight.by_beta),
            ("mach", &flight.by_mach),
            ("roll_rate", &flight.by_roll_rate),
            ("pitch_rate", &flight.by_pitch_rate),
            ("yaw_rate", &flight.by_yaw_rate),
            ("moment_ref_x", &flight.by_moment_reference[0]),
            ("moment_ref_y", &flight.by_moment_reference[1]),
            ("moment_ref_z", &flight.by_moment_reference[2]),
        ];
        for (input, gradients) in columns {
            for (output, value) in [
                ("CL", gradients.coefficients.lift),
                ("CD", gradients.coefficients.drag),
                ("CY", gradients.coefficients.side),
                ("Cl", gradients.coefficients.roll),
                ("Cm", gradients.coefficients.pitch),
                ("Cn", gradients.coefficients.yaw),
            ] {
                insert(output, input, vec![value]);
            }
        }

        if self
            .control_names()
            .iter()
            .any(|name| self.active(&format!("control_{name}")))
        {
            let controls = calculator
                .control_gradients_at_solution(solver, condition)
                .map_err(|error| rejected("compute_control_gradients", error.to_string()))?;
            for control in controls {
                let input = format!("control_{}", control.name);
                let values = control.by_deflection.coefficients;
                for (output, value) in [
                    ("CL", values.lift),
                    ("CD", values.drag),
                    ("CY", values.side),
                    ("Cl", values.roll),
                    ("Cm", values.pitch),
                    ("Cn", values.yaw),
                ] {
                    insert(output, &input, vec![value]);
                }
            }
        }

        // The coefficients are already normalized by dynamic pressure, so it
        // scales out entirely.
        for output in OUTPUTS {
            insert(output, "dynamic_pressure", vec![0.0]);
        }

        // Reference quantities enter only as normalizing divisors, so each
        // derivative is minus the coefficient over that reference. Span does
        // not appear in any of the six, and chord only in the moments.
        let area = mesh.reference.area;
        let chord = mesh.reference.chord;
        for (output, value) in [
            ("CL", results.coefficients.lift),
            ("CD", results.coefficients.drag),
            ("CY", results.coefficients.side),
            ("Cl", results.coefficients.roll),
            ("Cm", results.coefficients.pitch),
            ("Cn", results.coefficients.yaw),
        ] {
            insert(output, "reference_area", vec![-value / area]);
            insert(output, "reference_span", vec![0.0]);
            let by_chord = if matches!(output, "Cl" | "Cm" | "Cn") {
                -value / chord
            } else {
                0.0
            };
            insert(output, "reference_chord", vec![by_chord]);
        }

        Ok(partials)
    }
}

/// Accept the retained `consistent_forces` option only in its supported form.
///
/// Forces act at bound-vortex centers, with the full collocation self
/// influence. The historical v0.2 convention it once switched to is gone.
fn require_consistent_forces(value: &serde_json::Value) -> Result<()> {
    if super::option_bool(value) {
        return Ok(());
    }
    Err(invalid(
        "consistent_forces must be true: the legacy v0.2 force convention was removed, \
         and forces always act at bound-vortex centers",
    ))
}

fn analysis_output(name: &str) -> Result<Output> {
    match name {
        "CL" => Ok(Output::Lift),
        "CD" => Ok(Output::Drag),
        "CY" => Ok(Output::Side),
        "Cl" => Ok(Output::Roll),
        "Cm" => Ok(Output::Pitch),
        "Cn" => Ok(Output::Yaw),
        "CDi_near" => Ok(Output::NearFieldInducedDrag),
        "CDv" => Ok(Output::ViscousDrag),
        "CL_far" => Ok(Output::FarFieldLift),
        "CY_far" => Ok(Output::FarFieldSide),
        "CDi" => Ok(Output::FarFieldInducedDrag),
        "e" => Ok(Output::SpanEfficiency),
        _ => Err(invalid(format!("unknown analysis output '{name}'"))),
    }
}
fn output_name(output: Output) -> &'static str {
    match output {
        Output::Lift => "CL",
        Output::Drag => "CD",
        Output::Side => "CY",
        Output::Roll => "Cl",
        Output::Pitch => "Cm",
        Output::Yaw => "Cn",
        Output::NearFieldInducedDrag => "CDi_near",
        Output::ViscousDrag => "CDv",
        Output::FarFieldLift => "CL_far",
        Output::FarFieldSide => "CY_far",
        Output::FarFieldInducedDrag => "CDi",
        Output::SpanEfficiency => "e",
    }
}
impl VlmDiscipline {
    fn analysis_input(&self, name: &str) -> Result<Input> {
        if let Some(index) = self.surfaces.iter().position(|s| s.name == name) {
            return Ok(Input::MeshNodes(index));
        }
        match name {
            "alpha" => Ok(Input::Alpha),
            "beta" => Ok(Input::Beta),
            "mach" => Ok(Input::Mach),
            "roll_rate" => Ok(Input::RollRate),
            "pitch_rate" => Ok(Input::PitchRate),
            "yaw_rate" => Ok(Input::YawRate),
            "moment_ref_x" => Ok(Input::MomentReference(0)),
            "moment_ref_y" => Ok(Input::MomentReference(1)),
            "moment_ref_z" => Ok(Input::MomentReference(2)),
            "reference_area" => Ok(Input::ReferenceArea),
            "reference_chord" => Ok(Input::ReferenceChord),
            "reference_span" => Ok(Input::ReferenceSpan),
            "dynamic_pressure" => Ok(Input::DynamicPressure),
            _ => Err(invalid(format!("unknown analysis input '{name}'"))),
        }
    }
    fn input_name(&self, input: &Input) -> String {
        match input {
            Input::Alpha => "alpha",
            Input::Beta => "beta",
            Input::Mach => "mach",
            Input::RollRate => "roll_rate",
            Input::PitchRate => "pitch_rate",
            Input::YawRate => "yaw_rate",
            Input::MomentReference(0) => "moment_ref_x",
            Input::MomentReference(1) => "moment_ref_y",
            Input::MomentReference(_) => "moment_ref_z",
            Input::ReferenceArea => "reference_area",
            Input::ReferenceChord => "reference_chord",
            Input::ReferenceSpan => "reference_span",
            Input::DynamicPressure => "dynamic_pressure",
            Input::MeshNodes(i) => return self.surfaces[*i].name.clone(),
            _ => unreachable!("only structured-mesh and flight inputs are parsed"),
        }
        .to_owned()
    }
    fn evaluate_selected(
        &self,
        inputs: &ArrayMap,
        request: &AnalysisRequest,
    ) -> Result<crate::analysis::AnalysisResult> {
        let mesh = self.build_mesh(inputs)?;
        let condition = self.build_condition(inputs)?;
        let mut cache = self
            .analysis_cache
            .lock()
            .map_err(|_| invalid("analysis cache is poisoned"))?;
        if cache.as_ref().is_none_or(|(previous, _)| previous != &mesh) {
            let analysis = self.settings.analysis(&mesh);
            *cache = Some((mesh, analysis));
        }
        cache
            .as_mut()
            .unwrap()
            .1
            .evaluate(&condition, request)
            .map_err(|e| rejected("selective analysis", e.to_string()))
    }
}
