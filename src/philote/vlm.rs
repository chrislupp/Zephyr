//! Node mesh and flight condition to force and moment coefficients, with
//! adjoint partials.
//!
//! Ports the v0.2.0 C++ `VlmDiscipline`. Node partials come from
//! [`GradientCalculator::compute_mesh_node_gradients`], so the cost is one
//! transposed solve per output rather than one primal solve per node
//! coordinate — the reason to run this through an adjoint at all.

use std::{
    collections::HashMap,
    sync::{Mutex, MutexGuard},
};

use async_trait::async_trait;
use philote_mdo::{
    ArrayMap, PartialMap, Result, impl_registry,
    registry::VariableRegistry,
    traits::{Discipline, ExplicitDiscipline},
};

use super::{
    flat, invalid, option_tokens, option_usize, parse_bool, rejected, scalar, scalar_input,
    sized_input,
};
use crate::aero::ForceCalculator;
use crate::core::Vec3;
use crate::geometry::{MeshControlSurface, MeshInput, SurfaceMesh};
use crate::gradients::GradientCalculator;
use crate::solver::{FlightCondition, VortexLatticeSolver};

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

/// Flight-condition inputs that carry a nonzero derivative.
const FLIGHT_INPUTS: [&str; 6] = [
    "alpha",
    "beta",
    "mach",
    "roll_rate",
    "pitch_rate",
    "yaw_rate",
];

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
#[derive(Debug)]
pub struct VlmDiscipline {
    registry: VariableRegistry,
    surfaces: Vec<SurfaceOption>,
    consistent_forces: bool,
    outputs: Vec<String>,
    fixed_inputs: Vec<String>,
    cache: Mutex<Option<SolvedState>>,
}

#[derive(Debug)]
struct SolvedState {
    inputs: ArrayMap,
    mesh: MeshInput,
    condition: FlightCondition,
    solver: VortexLatticeSolver,
}

impl Default for VlmDiscipline {
    fn default() -> Self {
        Self {
            registry: VariableRegistry::new(false),
            surfaces: vec![SurfaceOption::default()],
            consistent_forces: false,
            outputs: OUTPUTS.iter().map(|s| (*s).into()).collect(),
            fixed_inputs: Vec::new(),
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
            consistent_forces: false,
            outputs: OUTPUTS.iter().map(|s| (*s).into()).collect(),
            fixed_inputs: Vec::new(),
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
            let mesh = self.build_mesh(inputs)?;
            let condition = self.build_condition(inputs)?;
            let mut solver = VortexLatticeSolver::new();
            solver
                .set_mesh(&mesh)
                .map_err(|e| rejected("set_mesh", e.to_string()))?;
            solver
                .solve(&condition)
                .map_err(|e| rejected("solve", e.to_string()))?;
            *cache = Some(SolvedState {
                inputs: inputs.clone(),
                mesh,
                condition,
                solver,
            });
        }
        Ok(cache)
    }

    fn active(&self, input: &str) -> bool {
        !self.fixed_inputs.iter().any(|name| name == input)
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
        let mut mesh = MeshInput {
            consistent_forces: self.consistent_forces,
            ..MeshInput::default()
        };
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
        Ok(FlightCondition {
            alpha: scalar_input(inputs, "alpha")?,
            beta: scalar_input(inputs, "beta")?,
            mach: scalar_input(inputs, "mach")?,
            dynamic_pressure: scalar_input(inputs, "dynamic_pressure")?,
            roll_rate: scalar_input(inputs, "roll_rate")?,
            pitch_rate: scalar_input(inputs, "pitch_rate")?,
            yaw_rate: scalar_input(inputs, "yaw_rate")?,
            control_deflections: self
                .control_names()
                .into_iter()
                .map(|name| {
                    let value = scalar_input(inputs, &format!("control_{name}"))?;
                    Ok((name, value))
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
        "PegasusVLM"
    }

    fn is_differentiable(&self) -> bool {
        true
    }

    fn provides_gradients(&self) -> bool {
        true
    }

    fn initialize(&mut self) -> Result<()> {
        self.add_option("num_surfaces", "int")?;
        self.add_option("consistent_forces", "bool")?;
        self.add_option("outputs", "str")?;
        self.add_option("fixed_inputs", "str")?;
        self.add_option("controls", "str")?;
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
        *self
            .cache
            .get_mut()
            .map_err(|_| rejected("cache", "poisoned solution cache"))? = None;
        if let Some(value) = options.get("outputs") {
            let outputs = option_tokens(value, "outputs")?;
            if outputs.is_empty() || outputs.iter().any(|s| !OUTPUTS.contains(&s.as_str())) {
                return Err(invalid("outputs must select from CL,CD,CY,Cl,Cm,Cn"));
            }
            self.outputs = outputs;
        }
        if let Some(value) = options.get("fixed_inputs") {
            self.fixed_inputs = option_tokens(value, "fixed_inputs")?;
        }
        if let Some(value) = options.get("consistent_forces") {
            self.consistent_forces = super::option_bool(value);
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
        Ok(())
    }

    fn setup(&mut self) -> Result<()> {
        for surface in self.surfaces.clone() {
            for control in &surface.controls {
                control
                    .validate_layout(surface.nx, surface.ny)
                    .map_err(|e| invalid(e.to_string()))?;
            }
            self.add_input(&surface.name, &[surface.node_values()], "m")?;
        }
        self.add_input("reference_area", &[1], "m**2")?;
        self.add_input("reference_chord", &[1], "m")?;
        self.add_input("reference_span", &[1], "m")?;
        self.add_input("moment_ref_x", &[1], "m")?;
        self.add_input("moment_ref_y", &[1], "m")?;
        self.add_input("moment_ref_z", &[1], "m")?;

        self.add_input("alpha", &[1], "rad")?;
        self.add_input("beta", &[1], "rad")?;
        self.add_input("mach", &[1], "")?;
        self.add_input("dynamic_pressure", &[1], "Pa")?;
        self.add_input("roll_rate", &[1], "")?;
        self.add_input("pitch_rate", &[1], "")?;
        self.add_input("yaw_rate", &[1], "")?;

        for name in self.control_names() {
            let input = format!("control_{name}");
            if self.surfaces.iter().any(|surface| surface.name == input) {
                return Err(invalid(format!(
                    "control input '{input}' conflicts with a mesh input"
                )));
            }
            self.add_input(&input, &[1], "rad")?;
        }
        for input in &self.fixed_inputs {
            let known = self
                .control_names()
                .iter()
                .any(|name| format!("control_{name}") == *input)
                || self.surfaces.iter().any(|s| s.name == *input)
                || FLIGHT_INPUTS.contains(&input.as_str())
                || [
                    "reference_area",
                    "reference_chord",
                    "reference_span",
                    "moment_ref_x",
                    "moment_ref_y",
                    "moment_ref_z",
                    "dynamic_pressure",
                ]
                .contains(&input.as_str());
            if !known {
                return Err(invalid(format!("unknown fixed input {input}")));
            }
        }
        for output in self.outputs.clone() {
            self.add_output(&output, &[1], "")?;
        }
        Ok(())
    }

    fn setup_partials(&mut self) -> Result<()> {
        let names: Vec<String> = self
            .surfaces
            .iter()
            .map(|s| s.name.clone())
            .chain(
                self.control_names()
                    .into_iter()
                    .map(|name| format!("control_{name}")),
            )
            .collect();
        for output in self.outputs.clone() {
            for name in &names {
                if self.active(name) {
                    self.declare_partials(&output, name)?;
                }
            }
            for input in [
                "reference_area",
                "reference_chord",
                "reference_span",
                "moment_ref_x",
                "moment_ref_y",
                "moment_ref_z",
                "dynamic_pressure",
            ] {
                if self.active(input) {
                    self.declare_partials(&output, input)?;
                }
            }
            for input in FLIGHT_INPUTS {
                if self.active(input) {
                    self.declare_partials(&output, input)?;
                }
            }
        }
        Ok(())
    }
}

#[async_trait]
impl ExplicitDiscipline for VlmDiscipline {
    async fn compute(&self, inputs: &ArrayMap) -> Result<ArrayMap> {
        let cache = self.solved(inputs)?;
        let state = cache.as_ref().expect("solved state");
        let (mesh, condition, solver) = (&state.mesh, &state.condition, &state.solver);
        let solution = solver.solution();
        let results = ForceCalculator::new().compute_with_references(
            &mesh.reference,
            solver.mesh(),
            solution,
            condition,
        );

        let mut outputs = HashMap::new();
        for (name, value) in [
            ("CL", results.coefficients.lift),
            ("CD", results.coefficients.drag),
            ("CY", results.coefficients.side),
            ("Cl", results.coefficients.roll),
            ("Cm", results.coefficients.pitch),
            ("Cn", results.coefficients.yaw),
        ] {
            if self.outputs.iter().any(|s| s == name) {
                outputs.insert(name.to_string(), scalar(value));
            }
        }
        Ok(outputs)
    }

    async fn compute_partials(&self, inputs: &ArrayMap) -> Result<PartialMap> {
        let cache = self.solved(inputs)?;
        let state = cache.as_ref().expect("solved state");
        let (mesh, condition, solver) = (&state.mesh, &state.condition, &state.solver);
        let solution = solver.solution();
        let results = ForceCalculator::new().compute_with_references(
            &mesh.reference,
            solver.mesh(),
            solution,
            condition,
        );

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
