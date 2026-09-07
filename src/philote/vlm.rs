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
use std::sync::Mutex;

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
use crate::geometry::{MeshInput, SurfaceMesh};
use crate::gradients::GradientCalculator;
use crate::solver::{FlightCondition, VortexLatticeSolver};

/// Per-surface configuration supplied through discipline options.
#[derive(Clone, Debug)]
pub struct SurfaceOption {
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
#[derive(Debug)]
pub struct VlmDiscipline {
    registry: VariableRegistry,
    surfaces: Vec<SurfaceOption>,
    analysis_request: Option<AnalysisRequest>,
    analysis_cache: Mutex<Option<(MeshInput, Analysis)>>,
}

impl Default for VlmDiscipline {
    fn default() -> Self {
        Self {
            registry: VariableRegistry::new(false),
            surfaces: vec![SurfaceOption::default()],
            analysis_request: None,
            analysis_cache: Mutex::new(None),
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
            analysis_request: None,
            analysis_cache: Mutex::new(None),
        }
    }

    /// The configured surfaces.
    #[must_use]
    pub fn surfaces(&self) -> &[SurfaceOption] {
        &self.surfaces
    }

    /// Assemble a [`MeshInput`] from the incoming node arrays.
    fn build_mesh(&self, inputs: &ArrayMap) -> Result<MeshInput> {
        let mut mesh = MeshInput::default();
        for surface in &self.surfaces {
            let nodes = sized_input(inputs, &surface.name, surface.node_values())?;
            let mut generated = SurfaceMesh::new(&surface.name, surface.nx, surface.ny);
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
            ..FlightCondition::default()
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
        self.add_option("analysis_outputs", "str")?;
        self.add_option("analysis_gradients", "dict")?;
        self.add_option("analysis_mode", "str")?;
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
        for surface in self.surfaces.clone() {
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

        if let Some(request) = &self.analysis_request {
            let mut outputs = request.outputs.clone();
            // A gradient output must exist in the discipline registry, even if
            // its value is only an internal dependency of the core request.
            outputs.extend(request.gradients.iter().map(|g| g.output));
            outputs.sort();
            outputs.dedup();
            for output in outputs {
                self.add_output(output_name(output), &[1], "")?;
            }
        } else {
            for output in OUTPUTS {
                self.add_output(output, &[1], "")?;
            }
        }
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
        let names: Vec<String> = self.surfaces.iter().map(|s| s.name.clone()).collect();
        for output in OUTPUTS {
            for name in &names {
                self.declare_partials(output, name)?;
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
                self.declare_partials(output, input)?;
            }
            for input in FLIGHT_INPUTS {
                self.declare_partials(output, input)?;
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
        let mesh = self.build_mesh(inputs)?;
        let condition = self.build_condition(inputs)?;

        let mut solver = VortexLatticeSolver::new();
        solver
            .set_mesh(&mesh)
            .map_err(|error| rejected("set_mesh", error.to_string()))?;
        let solution = solver
            .solve(&condition)
            .map_err(|error| rejected("solve", error.to_string()))?;
        let results = ForceCalculator::new().compute_with_references(
            &mesh.reference,
            solver.mesh(),
            &solution,
            &condition,
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
            outputs.insert(name.to_string(), scalar(value));
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
        let mesh = self.build_mesh(inputs)?;
        let condition = self.build_condition(inputs)?;

        let mut solver = VortexLatticeSolver::new();
        solver
            .set_mesh(&mesh)
            .map_err(|error| rejected("set_mesh", error.to_string()))?;
        let solution = solver
            .solve(&condition)
            .map_err(|error| rejected("solve", error.to_string()))?;
        let results = ForceCalculator::new().compute_with_references(
            &mesh.reference,
            solver.mesh(),
            &solution,
            &condition,
        );

        // Both gradient sets run against the solution above rather than
        // solving again: the factorization is the expensive part of an
        // operating point, and one linearize used to pay for it three times.
        let mut calculator = GradientCalculator::new();
        calculator
            .set_mesh(&mesh)
            .map_err(|error| rejected("set_mesh", error.to_string()))?;
        let node_gradients = calculator
            .mesh_node_gradients_from(&solver, &condition)
            .map_err(|error| rejected("compute_mesh_node_gradients", error.to_string()))?;
        let flight = calculator
            .flight_condition_gradients_from(&solver, &condition)
            .map_err(|error| rejected("compute_flight_condition_gradients", error.to_string()))?;

        let mut partials: PartialMap = HashMap::new();
        let mut insert = |output: &str, input: &str, values: Vec<f64>| {
            partials.insert((output.to_string(), input.to_string()), flat(values));
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
            let analysis = Analysis::from_mesh(&mesh);
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
