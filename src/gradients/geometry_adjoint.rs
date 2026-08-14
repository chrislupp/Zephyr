//! Adjoint contraction of lattice-geometry directions.

use crate::core::{EPSILON, MatX, Scalar, SymmetryType, Vec3, VecX};
use crate::error::{PegasusError, Result};
use crate::geometry::MeshInput;
use crate::solver::{FlightCondition, VortexLattice, VortexLatticeSolver};

use super::{
    AdjointSolver, AeroOutput, OutputSelection, aic_element_gradient, force_circulation_gradient,
    force_geometry_jvp, force_induced_velocity_gradient, horseshoe_velocity_vjp_bound_left,
    horseshoe_velocity_vjp_bound_right, horseshoe_velocity_vjp_field_point,
    trefftz_circulation_gradients, trefftz_geometry_jvp,
};

#[derive(Clone, Debug)]
pub(crate) struct LatticeDirections {
    pub bound_left: Vec<Vec3>,
    pub bound_right: Vec<Vec3>,
    pub control_point: Vec<Vec3>,
    pub normal: Vec<Vec3>,
    pub local_chord: Vec<Scalar>,
    pub strip_chord: Vec<Scalar>,
}

impl LatticeDirections {
    pub(crate) fn zeros(mesh: &VortexLattice) -> Self {
        Self {
            bound_left: vec![Vec3::zeros(); mesh.vortex_count()],
            bound_right: vec![Vec3::zeros(); mesh.vortex_count()],
            control_point: vec![Vec3::zeros(); mesh.vortex_count()],
            normal: vec![Vec3::zeros(); mesh.vortex_count()],
            local_chord: vec![0.0; mesh.vortex_count()],
            strip_chord: vec![0.0; mesh.strip_count()],
        }
    }
}

pub(crate) struct GeometryAdjointContext<'a> {
    solver: &'a VortexLatticeSolver,
    condition: &'a FlightCondition,
    outputs: Vec<AeroOutput>,
    adjoints: MatX,
}

impl<'a> GeometryAdjointContext<'a> {
    pub(crate) fn new(
        solver: &'a VortexLatticeSolver,
        condition: &'a FlightCondition,
        selection: OutputSelection,
    ) -> Result<Self> {
        if !solver.solution().is_valid {
            return Err(PegasusError::InvalidSolution);
        }
        if !condition.control_deflections.is_empty() {
            return Err(PegasusError::InvalidFlightCondition(
                "geometry adjoints with deflected control normals are not supported".into(),
            ));
        }
        if solver.symmetry().y_symmetry != SymmetryType::None
            || solver.symmetry().z_symmetry != SymmetryType::None
        {
            return Err(PegasusError::InvalidGeometry(
                "geometry adjoints require explicit mirrored lattice geometry".into(),
            ));
        }
        // Flight-condition gradients carry the viscous chain, but the geometry
        // path cannot yet: a section polar makes the strip's viscous drag a
        // function of its lift coefficient, and hence of the strip *area* and
        // the moment arm of its quarter chord. `LatticeDirections` has a
        // tangent for the strip chord but none for the strip width or the
        // quarter-chord position, so the area and arm sensitivities are
        // unavailable. Returning an inviscid-only geometry gradient for a
        // viscous mesh would be a silently wrong derivative, so this is an
        // error until those channels exist.
        if solver
            .mesh()
            .strips()
            .iter()
            .any(|strip| strip.has_viscous_data)
        {
            return Err(PegasusError::InvalidGeometry(
                "geometry adjoints do not yet support section drag polars; the strip \
                 area and quarter-chord sensitivities they need are not represented"
                    .into(),
            ));
        }
        let outputs: Vec<_> = [
            AeroOutput::Lift,
            AeroOutput::Drag,
            AeroOutput::Side,
            AeroOutput::Roll,
            AeroOutput::Pitch,
            AeroOutput::Yaw,
            AeroOutput::InducedDrag,
        ]
        .into_iter()
        .filter(|output| selection.contains(*output))
        .collect();
        let n = solver.mesh().vortex_count();
        let mut rhs = MatX::zeros(n, outputs.len());
        for (column, output) in outputs.iter().enumerate() {
            let gradient = if *output == AeroOutput::InducedDrag {
                trefftz_circulation_gradients(
                    solver.mesh(),
                    solver.solution(),
                    solver.reference_quantities().area,
                    0.01,
                    0.1,
                )
                .induced_drag
            } else {
                force_circulation_gradient(solver, condition, *output)
            };
            rhs.set_column(column, &gradient);
        }
        let adjoints = if outputs.is_empty() {
            MatX::zeros(n, 0)
        } else {
            AdjointSolver::new(solver)?.solve_batched(&rhs)?
        };
        Ok(Self {
            solver,
            condition,
            outputs,
            adjoints,
        })
    }

    /// Lattice the context was built against.
    pub(crate) fn mesh(&self) -> &VortexLattice {
        self.solver.mesh()
    }

    pub(crate) fn evaluate(
        &self,
        direction: &LatticeDirections,
    ) -> Result<Vec<(AeroOutput, Scalar)>> {
        let mesh = self.solver.mesh();
        if direction.bound_left.len() != mesh.vortex_count()
            || direction.bound_right.len() != mesh.vortex_count()
            || direction.control_point.len() != mesh.vortex_count()
            || direction.normal.len() != mesh.vortex_count()
        {
            return Err(PegasusError::InvalidGeometry(
                "lattice direction dimensions do not match the solved mesh".into(),
            ));
        }
        if self.outputs.is_empty() {
            return Ok(Vec::new());
        }
        let residual_direction = self.residual_direction(direction);
        let force_geometry = force_geometry_jvp(
            mesh,
            self.solver.solution(),
            self.condition,
            self.solver.reference_quantities(),
            &direction.bound_left,
            &direction.bound_right,
            &direction.control_point,
        )?;
        let trefftz_geometry = trefftz_geometry_jvp(
            mesh,
            self.solver.solution(),
            self.solver.reference_quantities().area,
            0.01,
            0.1,
            &direction.bound_left,
            &direction.bound_right,
            &direction.strip_chord,
        )?;
        let mut result = Vec::with_capacity(self.outputs.len());
        for (column, output) in self.outputs.iter().copied().enumerate() {
            let explicit = if output == AeroOutput::InducedDrag {
                trefftz_geometry[2]
            } else {
                let output_index = match output {
                    AeroOutput::Lift => 0,
                    AeroOutput::Drag => 1,
                    AeroOutput::Side => 2,
                    AeroOutput::Roll => 3,
                    AeroOutput::Pitch => 4,
                    AeroOutput::Yaw => 5,
                    AeroOutput::InducedDrag => unreachable!(),
                };
                let velocity_seed = force_induced_velocity_gradient(
                    mesh,
                    self.solver.solution(),
                    self.condition,
                    self.solver.reference_quantities(),
                    output,
                );
                force_geometry[output_index]
                    + induced_velocity_geometry_contraction(
                        mesh,
                        self.solver.solution().circulation.as_slice(),
                        &velocity_seed,
                        direction,
                        prandtl_glauert_factor(self.solver.current_mach()),
                    )
            };
            result.push((
                output,
                explicit + self.adjoints.column(column).dot(&residual_direction),
            ));
        }
        Ok(result)
    }

    fn residual_direction(&self, direction: &LatticeDirections) -> VecX {
        let mesh = self.solver.mesh();
        let solution = self.solver.solution();
        let pg = prandtl_glauert_factor(self.solver.current_mach());
        let freestream = VortexLatticeSolver::freestream_direction(self.condition);
        let omega = Vec3::new(
            self.condition.roll_rate,
            self.condition.pitch_rate,
            self.condition.yaw_rate,
        );
        VecX::from_fn(mesh.vortex_count(), |row, _| {
            let target = &mesh.vortices()[row];
            let local_velocity =
                freestream + omega.cross(&(target.control_point - self.solver.reference_point()));
            let rhs_direction = if target.affected_by_freestream {
                -omega
                    .cross(&direction.control_point[row])
                    .dot(&target.normal_at_control_point)
                    - local_velocity.dot(&direction.normal[row])
            } else {
                0.0
            };
            let mut matrix_product = 0.0;
            for (column, source) in mesh.vortices().iter().enumerate() {
                let core = core_radius(mesh, column, target.component_index);
                let gradient = aic_element_gradient(
                    &source.bound_vortex_left,
                    &source.bound_vortex_right,
                    &target.control_point,
                    &target.normal_at_control_point,
                    core,
                    pg,
                    row != column,
                );
                let entry_direction = gradient.bound_left.dot(&direction.bound_left[column])
                    + gradient.bound_right.dot(&direction.bound_right[column])
                    + gradient.control_point.dot(&direction.control_point[row])
                    + gradient.normal.dot(&direction.normal[row]);
                matrix_product += entry_direction * solution.circulation[column];
            }
            rhs_direction - matrix_product
        })
    }
}

fn prandtl_glauert_factor(mach: Scalar) -> Scalar {
    if mach < 1.0 {
        (1.0 - mach * mach).sqrt()
    } else {
        1.0
    }
}

fn core_radius(mesh: &VortexLattice, source_index: usize, target_component: usize) -> Scalar {
    let source = &mesh.vortices()[source_index];
    let width = (source.bound_vortex_right - source.bound_vortex_left).norm();
    if source.component_index == target_component {
        0.0001 * width
    } else {
        (0.0001 * source.local_chord).max(0.0001 * width)
    }
}

fn induced_velocity_geometry_contraction(
    mesh: &VortexLattice,
    circulation: &[Scalar],
    velocity_seeds: &[Vec3],
    direction: &LatticeDirections,
    pg: Scalar,
) -> Scalar {
    let mut result = 0.0;
    for (row, (target, seed)) in mesh.vortices().iter().zip(velocity_seeds).enumerate() {
        for (column, (source, gamma)) in mesh.vortices().iter().zip(circulation).enumerate() {
            if gamma.abs() <= EPSILON {
                continue;
            }
            let core = core_radius(mesh, column, target.component_index);
            let field_gradient = horseshoe_velocity_vjp_field_point(
                &target.control_point,
                &source.bound_vortex_left,
                &source.bound_vortex_right,
                pg,
                core,
                row != column,
                seed,
            );
            let left_gradient = horseshoe_velocity_vjp_bound_left(
                &target.control_point,
                &source.bound_vortex_left,
                &source.bound_vortex_right,
                pg,
                core,
                row != column,
                seed,
            );
            let right_gradient = horseshoe_velocity_vjp_bound_right(
                &target.control_point,
                &source.bound_vortex_left,
                &source.bound_vortex_right,
                pg,
                core,
                row != column,
                seed,
            );
            result += gamma
                * (field_gradient.dot(&direction.control_point[row])
                    + left_gradient.dot(&direction.bound_left[column])
                    + right_gradient.dot(&direction.bound_right[column]));
        }
    }
    result
}

pub(crate) fn direct_mesh_node_direction(
    input: &MeshInput,
    lattice: &VortexLattice,
    surface_index: usize,
    node_i: usize,
    node_j: usize,
    component: usize,
) -> Result<LatticeDirections> {
    let surface = input
        .surfaces
        .get(surface_index)
        .ok_or_else(|| PegasusError::InvalidMesh("surface index is out of range".into()))?;
    if node_i > surface.nx || node_j > surface.ny || component >= 3 {
        return Err(PegasusError::InvalidMesh(
            "mesh-node direction index is out of range".into(),
        ));
    }
    let surface_slot: usize = input.surfaces[..surface_index]
        .iter()
        .map(|surface| 1 + usize::from(surface.y_duplicated))
        .sum();
    let info = lattice
        .surface_info()
        .get(surface_slot)
        .ok_or_else(|| PegasusError::InvalidMesh("surface lattice is missing".into()))?;
    let mut result = LatticeDirections::zeros(lattice);
    let basis = Vec3::from_fn(|row, _| if row == component { 1.0 } else { 0.0 });
    let node = |i, j| *surface.node(i, j).expect("validated mesh index");

    for strip_offset in 0..surface.ny {
        let chord_left = node(surface.nx, strip_offset) - node(0, strip_offset);
        let chord_right = node(surface.nx, strip_offset + 1) - node(0, strip_offset + 1);
        let left_direction = if node_j == strip_offset {
            if node_i == surface.nx {
                basis
            } else if node_i == 0 {
                -basis
            } else {
                Vec3::zeros()
            }
        } else {
            Vec3::zeros()
        };
        let right_direction = if node_j == strip_offset + 1 {
            if node_i == surface.nx {
                basis
            } else if node_i == 0 {
                -basis
            } else {
                Vec3::zeros()
            }
        } else {
            Vec3::zeros()
        };
        let chord_direction = 0.5
            * (normalized_dot(&chord_left, &left_direction)
                + normalized_dot(&chord_right, &right_direction));
        let strip_index = info.first_strip_index + strip_offset;
        result.strip_chord[strip_index] = chord_direction;

        for chord_offset in 0..surface.nx {
            let n00 = node(chord_offset, strip_offset);
            let n10 = node(chord_offset + 1, strip_offset);
            let n01 = node(chord_offset, strip_offset + 1);
            let n11 = node(chord_offset + 1, strip_offset + 1);
            let d00 = node_direction(node_i, node_j, chord_offset, strip_offset, component);
            let d10 = node_direction(node_i, node_j, chord_offset + 1, strip_offset, component);
            let d01 = node_direction(node_i, node_j, chord_offset, strip_offset + 1, component);
            let d11 = node_direction(
                node_i,
                node_j,
                chord_offset + 1,
                strip_offset + 1,
                component,
            );
            let vortex_index = info.first_vortex_index + strip_offset * surface.nx + chord_offset;
            result.bound_left[vortex_index] = 0.75 * d00 + 0.25 * d10;
            result.bound_right[vortex_index] = 0.75 * d01 + 0.25 * d11;
            result.control_point[vortex_index] =
                0.125 * d00 + 0.375 * d10 + 0.125 * d01 + 0.375 * d11;
            let diagonal_a = n11 - n00;
            let diagonal_b = n10 - n01;
            let cross = diagonal_a.cross(&diagonal_b);
            let cross_direction = (d11 - d00).cross(&diagonal_b) + diagonal_a.cross(&(d10 - d01));
            result.normal[vortex_index] = normalized_direction(&cross, &cross_direction);
            result.local_chord[vortex_index] = chord_direction;
        }
    }

    if surface.y_duplicated {
        let duplicate = &lattice.surface_info()[surface_slot + 1];
        for strip_offset in 0..info.num_strips {
            let original_strip = info.first_strip_index + strip_offset;
            let duplicate_strip = duplicate.first_strip_index + strip_offset;
            result.strip_chord[duplicate_strip] = result.strip_chord[original_strip];
            for chord_offset in 0..info.num_chordwise {
                let original_index =
                    info.first_vortex_index + strip_offset * info.num_chordwise + chord_offset;
                let duplicate_index =
                    duplicate.first_vortex_index + strip_offset * info.num_chordwise + chord_offset;
                result.bound_left[duplicate_index] =
                    reflected_direction(result.bound_right[original_index]);
                result.bound_right[duplicate_index] =
                    reflected_direction(result.bound_left[original_index]);
                result.control_point[duplicate_index] =
                    reflected_direction(result.control_point[original_index]);
                result.normal[duplicate_index] = reflected_direction(result.normal[original_index]);
                result.local_chord[duplicate_index] = result.local_chord[original_index];
            }
        }
    }
    Ok(result)
}

fn node_direction(target_i: usize, target_j: usize, i: usize, j: usize, component: usize) -> Vec3 {
    if target_i == i && target_j == j {
        Vec3::from_fn(|row, _| if row == component { 1.0 } else { 0.0 })
    } else {
        Vec3::zeros()
    }
}

pub(crate) fn normalized_direction(vector: &Vec3, direction: &Vec3) -> Vec3 {
    let norm = vector.norm();
    if norm <= EPSILON {
        Vec3::zeros()
    } else {
        let unit = *vector / norm;
        (direction - unit * unit.dot(direction)) / norm
    }
}

fn normalized_dot(vector: &Vec3, direction: &Vec3) -> Scalar {
    let norm = vector.norm();
    if norm <= EPSILON {
        0.0
    } else {
        vector.dot(direction) / norm
    }
}

pub(crate) fn reflected_direction(direction: Vec3) -> Vec3 {
    Vec3::new(direction.x, -direction.y, direction.z)
}
