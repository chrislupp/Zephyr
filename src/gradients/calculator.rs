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

//! High-level aerodynamic, geometry, and structured-node gradients.

use crate::aero::Coefficients;
use crate::core::{MatX, Scalar, Vec3, VecX};
use crate::error::{Result, ZephyrError};
use crate::geometry::{Aircraft, MeshInput};
use crate::solver::{FlightCondition, VortexLatticeSolver};

use super::geometry_adjoint::{GeometryAdjointContext, scatter_mesh_node_cotangents};
use super::mach::velocity_influence_mach_derivative;
use super::section_geometry::section_direction;
use super::{
    AdjointSolver, AeroOutput, OutputSelection, force_alpha_gradient, force_beta_gradient,
    force_circulation_gradient, force_geometry_jvp, force_induced_velocity_gradient,
    trefftz_circulation_gradients,
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
    /// Derivatives per unit movement of the moment reference point, indexed
    /// x, y, z.
    ///
    /// Moving the reference point does two independent things: it changes every
    /// moment arm, and — because the point is also the center the body rates
    /// rotate about — it shifts the rotational onset velocity, which changes the
    /// circulation. The second term vanishes when all three rates are zero, so
    /// at a non-rotating condition only the moment outputs respond.
    pub by_moment_reference: [AeroGradients; 3],
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
            None => Err(ZephyrError::InvalidGeometry(
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
                return Err(ZephyrError::InvalidGeometry(
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
        self.flight_condition_gradients_from(&solver, condition)
    }

    /// As [`Self::compute_flight_condition_gradients`], reusing a solution the
    /// caller already has.
    ///
    /// `solver` must hold a valid solution of `condition` on this calculator's
    /// geometry. A caller that needs several gradient sets from one operating
    /// point pays for the factorization once instead of once per set.
    pub fn flight_condition_gradients_from(
        &self,
        solver: &VortexLatticeSolver,
        condition: &FlightCondition,
    ) -> Result<FullGradients> {
        let mut result = FullGradients::default();
        let velocity_mach = velocity_influence_mach_derivative(solver, true);
        let induced_velocity_mach = &velocity_mach * &solver.solution().circulation;
        let control_velocity_mach =
            velocity_influence_mach_derivative(solver, false) * &solver.solution().circulation;
        let residual_mach = VecX::from_fn(solver.mesh().vortex_count(), |row, _| {
            -solver.linearized_normals()[row].dot(&Vec3::new(
                control_velocity_mach[3 * row],
                control_velocity_mach[3 * row + 1],
                control_velocity_mach[3 * row + 2],
            ))
        });
        // Every selected output needs the same adjoint vector for alpha, beta
        // and Mach, so the circulation gradients are assembled once and solved
        // in a single batched transpose solve. Calling `TotalDerivative` here
        // and then rebuilding the gradient for the rate terms did both twice —
        // 2K gradient assemblies and 2K solves for K outputs.
        let selected: Vec<AeroOutput> = OUTPUTS
            .into_iter()
            .filter(|output| self.output_selection.contains(*output))
            .collect();
        if selected.is_empty() {
            return Ok(result);
        }
        let unknowns = solver.mesh().vortex_count();
        let mut circulation_gradients = MatX::zeros(unknowns, selected.len());
        for (column, &output) in selected.iter().enumerate() {
            let gradient = if output == AeroOutput::InducedDrag {
                trefftz_circulation_gradients(
                    solver.mesh(),
                    solver.solution(),
                    solver.reference_quantities().area,
                    0.01,
                    0.1,
                )
                .induced_drag
            } else {
                force_circulation_gradient(solver, condition, output)
            };
            circulation_gradients.set_column(column, &gradient);
        }
        let lambdas = AdjointSolver::new(solver)?.solve_batched(&circulation_gradients)?;

        // Freestream sensitivities of the residual, shared by every output.
        let (sin_alpha, cos_alpha) = condition.alpha.sin_cos();
        let (sin_beta, cos_beta) = condition.beta.sin_cos();
        let velocity_alpha = Vec3::new(-sin_alpha * cos_beta, 0.0, cos_alpha * cos_beta);
        let velocity_beta = Vec3::new(-cos_alpha * sin_beta, cos_beta, -sin_alpha * sin_beta);
        let mut rhs_alpha = VecX::zeros(unknowns);
        let mut rhs_beta = VecX::zeros(unknowns);
        for (index, (vortex, normal)) in solver
            .mesh()
            .vortices()
            .iter()
            .zip(solver.linearized_normals())
            .enumerate()
        {
            if vortex.affected_by_freestream {
                rhs_alpha[index] = -normal.dot(&velocity_alpha);
                rhs_beta[index] = -normal.dot(&velocity_beta);
            }
        }
        // Moment-reference sensitivities, also shared by every output.
        //
        // The reference point enters the residual only as the center the body
        // rates rotate about: the onset term is `omega x (r_cp - r_ref)`, whose
        // derivative along `e_k` is the *uniform* field `-(omega x e_k)`. That
        // makes the residual sensitivity `(omega x e_k) . n`, the same shape as
        // the alpha and beta rows above, and it is identically zero when the
        // aircraft is not rotating.
        let omega = Vec3::new(
            condition.roll_rate,
            condition.pitch_rate,
            condition.yaw_rate,
        );
        let mut rhs_moment_reference = [
            VecX::zeros(unknowns),
            VecX::zeros(unknowns),
            VecX::zeros(unknowns),
        ];
        if omega != Vec3::zeros() {
            for axis in 0..3 {
                let mut unit = Vec3::zeros();
                unit[axis] = 1.0;
                let sensitivity = omega.cross(&unit);
                for (index, (vortex, normal)) in solver
                    .mesh()
                    .vortices()
                    .iter()
                    .zip(solver.linearized_normals())
                    .enumerate()
                {
                    if vortex.affected_by_freestream {
                        rhs_moment_reference[axis][index] = normal.dot(&sensitivity);
                    }
                }
            }
        }
        // The explicit half is the moment-arm change at frozen circulation.
        // `force_geometry_jvp` already expresses exactly that: moving the
        // reference point by `+e_k` is moving every control point by `-e_k`
        // relative to it, with the vortex lines held still. Reusing it rather
        // than rewriting the contraction keeps this tied to the same force
        // convention and normalization the primal uses.
        let still = vec![Vec3::zeros(); unknowns];
        let mut explicit_moment_reference = [[0.0; 6]; 3];
        for axis in 0..3 {
            let mut shift = Vec3::zeros();
            shift[axis] = -1.0;
            explicit_moment_reference[axis] = force_geometry_jvp(
                solver.mesh(),
                solver.solution(),
                condition,
                solver.reference_quantities(),
                &still,
                &still,
                &vec![shift; unknowns],
            )?;
        }

        // Both return all six components at once, so they are built once
        // rather than once per output.
        let explicit_alpha = force_alpha_gradient(
            solver.mesh(),
            solver.solution(),
            condition,
            solver.reference_quantities(),
        );
        let explicit_beta = force_beta_gradient(
            solver.mesh(),
            solver.solution(),
            condition,
            solver.reference_quantities(),
        );

        for (column, &output) in selected.iter().enumerate() {
            let circulation_gradient = circulation_gradients.column(column).into_owned();
            let lambda = lambdas.column(column).into_owned();
            let output_index = match output {
                AeroOutput::Lift => 0,
                AeroOutput::Drag => 1,
                AeroOutput::Side => 2,
                AeroOutput::Roll => 3,
                AeroOutput::Pitch => 4,
                AeroOutput::Yaw => 5,
                AeroOutput::InducedDrag => 0,
            };
            let (alpha_term, beta_term) = if output == AeroOutput::InducedDrag {
                (0.0, 0.0)
            } else {
                (explicit_alpha[output_index], explicit_beta[output_index])
            };
            result
                .by_alpha
                .set(output, alpha_term + lambda.dot(&rhs_alpha));
            result
                .by_beta
                .set(output, beta_term + lambda.dot(&rhs_beta));
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
            result
                .by_mach
                .set(output, explicit_mach + lambda.dot(&residual_mach));
            for axis in 0..3 {
                // Induced drag comes from the Trefftz plane, which knows nothing
                // of the moment reference, so it has no explicit term.
                let explicit = if output == AeroOutput::InducedDrag {
                    0.0
                } else {
                    explicit_moment_reference[axis][output_index]
                };
                result.by_moment_reference[axis]
                    .set(output, explicit + lambda.dot(&rhs_moment_reference[axis]));
            }
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
            return Err(ZephyrError::InvalidGeometry(
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
            return Err(ZephyrError::InvalidGeometry(
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
        let direction = section_direction(
            aircraft,
            context.mesh(),
            surface_index,
            section_index,
            variable,
        )?;
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
        let mut solver = self.configured_solver()?;
        solver.solve(condition)?;
        self.mesh_node_gradients_from(&solver, condition)
    }

    /// As [`Self::compute_mesh_node_gradients`], reusing a solution the caller
    /// already has.
    ///
    /// `solver` must hold a valid solution of `condition` on this calculator's
    /// geometry.
    pub fn mesh_node_gradients_from(
        &self,
        solver: &VortexLatticeSolver,
        condition: &FlightCondition,
    ) -> Result<Vec<SurfaceMeshNodeGradients>> {
        let Some(GeometrySource::Mesh(mesh)) = self.source.as_ref() else {
            return Err(ZephyrError::InvalidMesh(
                "node gradients require a direct mesh source".into(),
            ));
        };
        let context = GeometryAdjointContext::new(solver, condition, self.output_selection)?;
        let mut result: Vec<_> = mesh
            .surfaces
            .iter()
            .map(|surface| {
                SurfaceMeshNodeGradients::zeros(surface.name.clone(), surface.nx, surface.ny)
            })
            .collect();
        for (surface_index, surface_result) in result.iter_mut().enumerate() {
            for (output, cotangent) in context.cotangents() {
                let nodes =
                    scatter_mesh_node_cotangents(mesh, solver.mesh(), surface_index, cotangent)?;
                surface_result.output_mut(output).nodes = nodes;
            }
        }
        Ok(result)
    }
}
