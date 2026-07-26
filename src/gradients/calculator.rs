//! High-level aerodynamic, geometry, and structured-node gradients.

use crate::aero::Coefficients;
use crate::core::{EPSILON, Scalar, Vec3, VecX};
use crate::error::{PegasusError, Result};
use crate::geometry::{Aircraft, MeshInput};
use crate::solver::{FlightCondition, VortexLatticeSolver};

use super::geometry_adjoint::{
    GeometryAdjointContext, LatticeDirections, direct_mesh_node_direction,
};
use super::mach::velocity_influence_mach_derivative;
use super::{
    AdjointSolver, AeroOutput, OutputSelection, TotalDerivative, force_circulation_gradient,
    force_induced_velocity_gradient, trefftz_circulation_gradients,
};

/// Derivatives of all aerodynamic outputs with respect to one input.
#[derive(Clone, Copy, Debug, Default)]
pub struct AeroGradients {
    /// Primary force and moment coefficient derivatives.
    pub coefficients: Coefficients,
    /// Far-field induced-drag derivative.
    pub induced_drag: Scalar,
}

impl AeroGradients {
    /// Read a selected output.
    #[must_use]
    pub fn get(self, output: AeroOutput) -> Scalar {
        match output {
            AeroOutput::Lift => self.coefficients.lift,
            AeroOutput::Drag => self.coefficients.drag,
            AeroOutput::Side => self.coefficients.side,
            AeroOutput::Roll => self.coefficients.roll,
            AeroOutput::Pitch => self.coefficients.pitch,
            AeroOutput::Yaw => self.coefficients.yaw,
            AeroOutput::InducedDrag => self.induced_drag,
        }
    }

    pub(crate) fn set(&mut self, output: AeroOutput, value: Scalar) {
        match output {
            AeroOutput::Lift => self.coefficients.lift = value,
            AeroOutput::Drag => self.coefficients.drag = value,
            AeroOutput::Side => self.coefficients.side = value,
            AeroOutput::Roll => self.coefficients.roll = value,
            AeroOutput::Pitch => self.coefficients.pitch = value,
            AeroOutput::Yaw => self.coefficients.yaw = value,
            AeroOutput::InducedDrag => self.induced_drag = value,
        }
    }
}

/// Complete flight-condition derivative set.
#[derive(Clone, Copy, Debug, Default)]
pub struct FullGradients {
    /// Derivatives per radian of angle of attack.
    pub by_alpha: AeroGradients,
    /// Derivatives per radian of sideslip.
    pub by_beta: AeroGradients,
    /// Derivatives per unit normalized roll rate.
    pub by_roll_rate: AeroGradients,
    /// Derivatives per unit normalized pitch rate.
    pub by_pitch_rate: AeroGradients,
    /// Derivatives per unit normalized yaw rate.
    pub by_yaw_rate: AeroGradients,
    /// Derivatives per unit Mach number.
    pub by_mach: AeroGradients,
}

/// Section design variable.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum GeometryVariable {
    /// Section chord.
    Chord,
    /// Section twist in degrees.
    Twist,
    /// Leading-edge x position.
    LeadingEdgeX,
    /// Leading-edge y position.
    LeadingEdgeY,
    /// Leading-edge z position.
    LeadingEdgeZ,
}

/// Gradients for all five geometric variables on one section.
#[derive(Clone, Debug, Default)]
pub struct SectionGeometryGradients {
    /// Owning surface name.
    pub surface_name: String,
    /// Zero-based section index.
    pub section_index: usize,
    /// Chord derivatives.
    pub by_chord: AeroGradients,
    /// Twist derivatives, per degree.
    pub by_twist: AeroGradients,
    /// Leading-edge x derivatives.
    pub by_leading_edge_x: AeroGradients,
    /// Leading-edge y derivatives.
    pub by_leading_edge_y: AeroGradients,
    /// Leading-edge z derivatives.
    pub by_leading_edge_z: AeroGradients,
}

impl SectionGeometryGradients {
    /// Read one geometric derivative group.
    #[must_use]
    pub fn get(&self, variable: GeometryVariable) -> AeroGradients {
        match variable {
            GeometryVariable::Chord => self.by_chord,
            GeometryVariable::Twist => self.by_twist,
            GeometryVariable::LeadingEdgeX => self.by_leading_edge_x,
            GeometryVariable::LeadingEdgeY => self.by_leading_edge_y,
            GeometryVariable::LeadingEdgeZ => self.by_leading_edge_z,
        }
    }
}

/// Geometry gradients for every aircraft section.
#[derive(Clone, Debug, Default)]
pub struct GeometryGradients {
    /// Surface-major section gradients.
    pub sections: Vec<SectionGeometryGradients>,
}

impl GeometryGradients {
    /// Number of scalar section design variables represented.
    #[must_use]
    pub fn design_variable_count(&self) -> usize {
        self.sections.len() * 5
    }

    /// Extract all chord derivatives for one aerodynamic output.
    #[must_use]
    pub fn chord_gradients(&self, output: AeroOutput) -> Vec<Scalar> {
        self.sections
            .iter()
            .map(|section| section.by_chord.get(output))
            .collect()
    }

    /// Extract all twist derivatives for one aerodynamic output.
    #[must_use]
    pub fn twist_gradients(&self, output: AeroOutput) -> Vec<Scalar> {
        self.sections
            .iter()
            .map(|section| section.by_twist.get(output))
            .collect()
    }
}

/// Gradient of one scalar output with respect to every structured mesh node.
#[derive(Clone, Debug)]
pub struct MeshNodeGradients {
    /// Chordwise panel count.
    pub nx: usize,
    /// Spanwise panel count.
    pub ny: usize,
    /// Row-major node-coordinate gradients.
    pub nodes: Vec<Vec3>,
}

impl MeshNodeGradients {
    fn zeros(nx: usize, ny: usize) -> Self {
        Self {
            nx,
            ny,
            nodes: vec![Vec3::zeros(); (nx + 1) * (ny + 1)],
        }
    }

    /// Checked gradient access.
    #[must_use]
    pub fn node(&self, i: usize, j: usize) -> Option<&Vec3> {
        (i <= self.nx && j <= self.ny).then(|| &self.nodes[i * (self.ny + 1) + j])
    }
}

/// Node-coordinate gradients of every primary aerodynamic coefficient on one
/// input surface.
#[derive(Clone, Debug)]
pub struct SurfaceMeshNodeGradients {
    /// Surface name.
    pub surface_name: String,
    /// Lift gradients.
    pub lift: MeshNodeGradients,
    /// Drag gradients.
    pub drag: MeshNodeGradients,
    /// Side-force gradients.
    pub side: MeshNodeGradients,
    /// Rolling-moment gradients.
    pub roll: MeshNodeGradients,
    /// Pitching-moment gradients.
    pub pitch: MeshNodeGradients,
    /// Yawing-moment gradients.
    pub yaw: MeshNodeGradients,
    /// Far-field induced-drag gradients.
    pub induced_drag: MeshNodeGradients,
}

impl SurfaceMeshNodeGradients {
    fn zeros(name: String, nx: usize, ny: usize) -> Self {
        Self {
            surface_name: name,
            lift: MeshNodeGradients::zeros(nx, ny),
            drag: MeshNodeGradients::zeros(nx, ny),
            side: MeshNodeGradients::zeros(nx, ny),
            roll: MeshNodeGradients::zeros(nx, ny),
            pitch: MeshNodeGradients::zeros(nx, ny),
            yaw: MeshNodeGradients::zeros(nx, ny),
            induced_drag: MeshNodeGradients::zeros(nx, ny),
        }
    }

    fn output_mut(&mut self, output: AeroOutput) -> &mut MeshNodeGradients {
        match output {
            AeroOutput::Lift => &mut self.lift,
            AeroOutput::Drag => &mut self.drag,
            AeroOutput::Side => &mut self.side,
            AeroOutput::Roll => &mut self.roll,
            AeroOutput::Pitch => &mut self.pitch,
            AeroOutput::Yaw => &mut self.yaw,
            AeroOutput::InducedDrag => &mut self.induced_drag,
        }
    }
}

#[derive(Clone, Debug)]
enum GeometrySource {
    Aircraft(Aircraft),
    Mesh(MeshInput),
}

/// High-level gradient interface for aircraft and direct structured meshes.
#[derive(Clone, Debug, Default)]
pub struct GradientCalculator {
    source: Option<GeometrySource>,
    output_selection: OutputSelection,
}

const OUTPUTS: [AeroOutput; 7] = [
    AeroOutput::Lift,
    AeroOutput::Drag,
    AeroOutput::Side,
    AeroOutput::Roll,
    AeroOutput::Pitch,
    AeroOutput::Yaw,
    AeroOutput::InducedDrag,
];

impl GradientCalculator {
    /// Construct an empty calculator.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Use aircraft geometry as the design source.
    pub fn set_aircraft(&mut self, aircraft: &Aircraft) {
        self.source = Some(GeometrySource::Aircraft(aircraft.clone()));
    }

    /// Use a direct structured mesh as the design source.
    pub fn set_mesh(&mut self, mesh: &MeshInput) -> Result<()> {
        mesh.validate()?;
        self.source = Some(GeometrySource::Mesh(mesh.clone()));
        Ok(())
    }

    /// Change the moment reference point without changing other geometry.
    pub fn set_reference_point(&mut self, point: Vec3) -> Result<()> {
        match self.source.as_mut() {
            Some(GeometrySource::Aircraft(aircraft)) => {
                aircraft.set_moment_reference_point(point);
                Ok(())
            }
            Some(GeometrySource::Mesh(mesh)) => {
                mesh.reference.moment_reference_point = point;
                Ok(())
            }
            None => Err(PegasusError::InvalidGeometry(
                "set aircraft or mesh before setting the reference point".into(),
            )),
        }
    }

    /// Select outputs. Unselected output arrays remain exactly zero.
    pub fn set_output_selection(&mut self, selection: OutputSelection) {
        self.output_selection = selection;
    }

    fn configured_solver(&self) -> Result<VortexLatticeSolver> {
        let mut solver = VortexLatticeSolver::new();
        match self.source.as_ref() {
            Some(GeometrySource::Aircraft(aircraft)) => solver.set_aircraft(aircraft)?,
            Some(GeometrySource::Mesh(mesh)) => solver.set_mesh(mesh)?,
            None => {
                return Err(PegasusError::InvalidGeometry(
                    "no aircraft or mesh has been set".into(),
                ));
            }
        }
        Ok(solver)
    }

    /// Compute analytical angle, rate, and Mach derivatives.
    pub fn compute_flight_condition_gradients(
        &self,
        condition: &FlightCondition,
    ) -> Result<FullGradients> {
        let mut solver = self.configured_solver()?;
        solver.solve(condition)?;
        let mut result = FullGradients::default();
        let velocity_mach = velocity_influence_mach_derivative(&solver);
        let induced_velocity_mach = &velocity_mach * &solver.solution().circulation;
        let residual_mach = VecX::from_fn(solver.mesh().vortex_count(), |row, _| {
            -solver.linearized_normals()[row].dot(&Vec3::new(
                induced_velocity_mach[3 * row],
                induced_velocity_mach[3 * row + 1],
                induced_velocity_mach[3 * row + 2],
            ))
        });
        for output in OUTPUTS {
            if !self.output_selection.contains(output) {
                continue;
            }
            let total = TotalDerivative::compute(&solver, condition, output)?;
            result.by_alpha.set(output, total.by_alpha);
            result.by_beta.set(output, total.by_beta);
            let circulation_gradient = if output == AeroOutput::InducedDrag {
                trefftz_circulation_gradients(
                    solver.mesh(),
                    solver.solution(),
                    solver.reference_quantities().area,
                    0.01,
                    0.1,
                )
                .induced_drag
            } else {
                force_circulation_gradient(&solver, condition, output)
            };
            result.by_roll_rate.set(
                output,
                circulation_gradient.dot(&solver.solution().circulation_per_p),
            );
            result.by_pitch_rate.set(
                output,
                circulation_gradient.dot(&solver.solution().circulation_per_q),
            );
            result.by_yaw_rate.set(
                output,
                circulation_gradient.dot(&solver.solution().circulation_per_r),
            );
            let explicit_mach = if output == AeroOutput::InducedDrag {
                0.0
            } else {
                force_induced_velocity_gradient(
                    solver.mesh(),
                    solver.solution(),
                    condition,
                    solver.reference_quantities(),
                    output,
                )
                .iter()
                .enumerate()
                .map(|(row, seed)| {
                    seed.dot(&Vec3::new(
                        induced_velocity_mach[3 * row],
                        induced_velocity_mach[3 * row + 1],
                        induced_velocity_mach[3 * row + 2],
                    ))
                })
                .sum()
            };
            let lambda = AdjointSolver::new(&solver)?.solve(&circulation_gradient)?;
            result
                .by_mach
                .set(output, explicit_mach + lambda.dot(&residual_mach));
        }
        Ok(result)
    }

    /// Compute derivatives with respect to angle of attack.
    pub fn compute_alpha_gradients(&self, condition: &FlightCondition) -> Result<AeroGradients> {
        Ok(self.compute_flight_condition_gradients(condition)?.by_alpha)
    }

    /// Compute derivatives with respect to sideslip.
    pub fn compute_beta_gradients(&self, condition: &FlightCondition) -> Result<AeroGradients> {
        Ok(self.compute_flight_condition_gradients(condition)?.by_beta)
    }

    /// Compute all section geometry derivatives.
    ///
    /// Aerodynamic derivatives use one primal solve plus batched adjoint
    /// contractions. Section-to-lattice directions are isolated from the
    /// aerodynamic solve and evaluated symmetrically to preserve the legacy
    /// mesher's interpolation conventions.
    pub fn compute_geometry_gradients(
        &self,
        condition: &FlightCondition,
    ) -> Result<GeometryGradients> {
        let Some(GeometrySource::Aircraft(aircraft)) = self.source.as_ref() else {
            return Err(PegasusError::InvalidGeometry(
                "section gradients require an aircraft source".into(),
            ));
        };
        let mut solver = self.configured_solver()?;
        solver.solve(condition)?;
        let context = GeometryAdjointContext::new(&solver, condition, self.output_selection)?;
        let mut result = GeometryGradients::default();
        for (surface_index, surface) in aircraft.surfaces().iter().enumerate() {
            for section_index in 0..surface.sections().len() {
                let mut section_result = SectionGeometryGradients {
                    surface_name: surface.name().into(),
                    section_index,
                    ..SectionGeometryGradients::default()
                };
                for variable in [
                    GeometryVariable::Chord,
                    GeometryVariable::Twist,
                    GeometryVariable::LeadingEdgeX,
                    GeometryVariable::LeadingEdgeY,
                    GeometryVariable::LeadingEdgeZ,
                ] {
                    let gradient = self.compute_section_gradient_with_context(
                        aircraft,
                        &context,
                        surface_index,
                        section_index,
                        variable,
                    )?;
                    match variable {
                        GeometryVariable::Chord => section_result.by_chord = gradient,
                        GeometryVariable::Twist => section_result.by_twist = gradient,
                        GeometryVariable::LeadingEdgeX => {
                            section_result.by_leading_edge_x = gradient;
                        }
                        GeometryVariable::LeadingEdgeY => {
                            section_result.by_leading_edge_y = gradient;
                        }
                        GeometryVariable::LeadingEdgeZ => {
                            section_result.by_leading_edge_z = gradient;
                        }
                    }
                }
                result.sections.push(section_result);
            }
        }
        Ok(result)
    }

    /// Compute one section geometry derivative group.
    pub fn compute_section_gradient(
        &self,
        condition: &FlightCondition,
        surface_index: usize,
        section_index: usize,
        variable: GeometryVariable,
    ) -> Result<AeroGradients> {
        let Some(GeometrySource::Aircraft(aircraft)) = self.source.as_ref() else {
            return Err(PegasusError::InvalidGeometry(
                "section gradients require an aircraft source".into(),
            ));
        };
        let mut solver = self.configured_solver()?;
        solver.solve(condition)?;
        let context = GeometryAdjointContext::new(&solver, condition, self.output_selection)?;
        self.compute_section_gradient_with_context(
            aircraft,
            &context,
            surface_index,
            section_index,
            variable,
        )
    }

    fn compute_section_gradient_with_context(
        &self,
        aircraft: &Aircraft,
        context: &GeometryAdjointContext<'_>,
        surface_index: usize,
        section_index: usize,
        variable: GeometryVariable,
    ) -> Result<AeroGradients> {
        let section = aircraft
            .surfaces()
            .get(surface_index)
            .and_then(|surface| surface.section(section_index))
            .ok_or_else(|| PegasusError::InvalidGeometry("section index is out of range".into()))?;
        let base_value = match variable {
            GeometryVariable::Chord => section.chord(),
            GeometryVariable::Twist => section.twist_degrees(),
            GeometryVariable::LeadingEdgeX => section.leading_edge().x,
            GeometryVariable::LeadingEdgeY => section.leading_edge().y,
            GeometryVariable::LeadingEdgeZ => section.leading_edge().z,
        };
        let step = match variable {
            GeometryVariable::Twist => 1.0e-4,
            _ => 1.0e-6 * base_value.abs().max(1.0),
        };
        let mut plus = aircraft.clone();
        let mut minus = aircraft.clone();
        perturb_section(&mut plus, surface_index, section_index, variable, step)?;
        perturb_section(&mut minus, surface_index, section_index, variable, -step)?;
        let mut plus_solver = VortexLatticeSolver::new();
        plus_solver.set_aircraft(&plus)?;
        let mut minus_solver = VortexLatticeSolver::new();
        minus_solver.set_aircraft(&minus)?;
        let direction =
            LatticeDirections::difference(plus_solver.mesh(), minus_solver.mesh(), 2.0 * step)?;
        let mut gradient = AeroGradients::default();
        for (output, value) in context.evaluate(&direction)? {
            gradient.set(output, value);
        }
        Ok(gradient)
    }

    /// Compute gradients with respect to every direct-mesh node coordinate.
    pub fn compute_mesh_node_gradients(
        &self,
        condition: &FlightCondition,
    ) -> Result<Vec<SurfaceMeshNodeGradients>> {
        let Some(GeometrySource::Mesh(mesh)) = self.source.as_ref() else {
            return Err(PegasusError::InvalidMesh(
                "node gradients require a direct mesh source".into(),
            ));
        };
        let mut solver = self.configured_solver()?;
        solver.solve(condition)?;
        let context = GeometryAdjointContext::new(&solver, condition, self.output_selection)?;
        let mut result: Vec<_> = mesh
            .surfaces
            .iter()
            .map(|surface| {
                SurfaceMeshNodeGradients::zeros(surface.name.clone(), surface.nx, surface.ny)
            })
            .collect();
        for (surface_index, surface_result) in result.iter_mut().enumerate() {
            let surface = &mesh.surfaces[surface_index];
            for i in 0..=surface.nx {
                for j in 0..=surface.ny {
                    let node_index = i * (surface.ny + 1) + j;
                    for component in 0..3 {
                        let direction = direct_mesh_node_direction(
                            mesh,
                            solver.mesh(),
                            surface_index,
                            i,
                            j,
                            component,
                        )?;
                        for (output, value) in context.evaluate(&direction)? {
                            surface_result.output_mut(output).nodes[node_index][component] = value;
                        }
                    }
                }
            }
        }
        Ok(result)
    }

    /// Compatibility no-op. Analytical kernels do not consume a finite
    /// difference step.
    #[deprecated(note = "analytical derivatives do not use a configurable finite-difference step")]
    pub fn set_finite_difference_step(&mut self, _step: Scalar) {}
}

fn perturb_section(
    aircraft: &mut Aircraft,
    surface_index: usize,
    section_index: usize,
    variable: GeometryVariable,
    delta: Scalar,
) -> Result<()> {
    let section = aircraft
        .surfaces_mut()
        .get_mut(surface_index)
        .and_then(|surface| surface.section_mut(section_index))
        .ok_or_else(|| PegasusError::InvalidGeometry("section index is out of range".into()))?;
    match variable {
        GeometryVariable::Chord => {
            let chord = section.chord() + delta;
            if chord <= EPSILON {
                return Err(PegasusError::InvalidGeometry(
                    "chord perturbation is nonpositive".into(),
                ));
            }
            section.set_chord(chord);
        }
        GeometryVariable::Twist => section.set_twist_degrees(section.twist_degrees() + delta),
        GeometryVariable::LeadingEdgeX
        | GeometryVariable::LeadingEdgeY
        | GeometryVariable::LeadingEdgeZ => {
            let mut point = *section.leading_edge();
            let component = match variable {
                GeometryVariable::LeadingEdgeX => 0,
                GeometryVariable::LeadingEdgeY => 1,
                GeometryVariable::LeadingEdgeZ => 2,
                GeometryVariable::Chord | GeometryVariable::Twist => unreachable!(),
            };
            point[component] += delta;
            section.set_leading_edge(point);
        }
    }
    Ok(())
}
