//! Hand-derived mesh JVPs and VJPs about one immutable solved state.
use super::coupled::LatticeCotangents;
use super::geometry_adjoint::{LatticeDirections, normalized_direction, reflected_direction};
use super::{
    AdjointSolver, AeroGradients, AeroOutput, MeshNodeCotangent, OutputSelection, StripLoadAdjoint,
    StripLoadCotangent, force_circulation_gradient, force_geometry_jvp,
    force_induced_velocity_gradient, horseshoe_velocity_jacobians, trefftz_circulation_gradients,
    trefftz_geometry_jvp, trefftz_geometry_vjp,
};
use crate::core::{MatX, Vec3};
use crate::error::{PegasusError, Result};
use crate::geometry::MeshInput;
use crate::solver::{FlightCondition, VortexLatticeSolver, horseshoe_vortex_velocity};

/// Canonical output/seed order. Moments use Pegasus's reference chord for all
/// three axes. Reference quantities and the flight condition are held fixed.
pub const COEFFICIENT_ORDER: [AeroOutput; 7] = [
    AeroOutput::Lift,
    AeroOutput::Drag,
    AeroOutput::Side,
    AeroOutput::Roll,
    AeroOutput::Pitch,
    AeroOutput::Yaw,
    AeroOutput::InducedDrag,
];

/// A tangent has the same surface/node/XYZ storage as a cotangent.
pub type MeshNodeDirection = MeshNodeCotangent;

/// Linearization of direct-mesh aerodynamic coefficients. Borrows a solved
/// state and its original mesh and flight condition, preventing stale-state
/// mixing. Both products reuse the primal LU. No mesh-coordinate Jacobian is
/// formed. Mirrors through `y_duplicated` are supported; symmetry planes,
/// deflected control normals, and section polars retain the adjoint restrictions.
pub struct MeshLinearization<'a> {
    solver: &'a VortexLatticeSolver,
    condition: &'a FlightCondition,
    input: &'a MeshInput,
    loads: StripLoadAdjoint<'a>,
}
impl<'a> MeshLinearization<'a> {
    pub fn new(solver: &'a VortexLatticeSolver) -> Result<Self> {
        let condition = solver
            .solution_condition()
            .ok_or(PegasusError::InvalidSolution)?;
        let input = solver.mesh().mesh_input().ok_or_else(|| {
            PegasusError::InvalidMesh("mesh products require a direct mesh source".into())
        })?;
        for surface in &input.surfaces {
            for j in 0..=surface.ny {
                if (surface.nodes[surface.nx * (surface.ny + 1) + j] - surface.nodes[j]).norm()
                    <= crate::core::EPSILON
                {
                    return Err(PegasusError::InvalidGeometry(
                        "mesh products require noncollapsed chord stations".into(),
                    ));
                }
            }
        }
        let points = vec![*solver.reference_point(); solver.mesh().strip_count()];
        let loads = StripLoadAdjoint::new(solver, condition, &points)?;
        Ok(Self {
            solver,
            condition,
            input,
            loads,
        })
    }

    /// Weighted total coefficient gradients with respect to input mesh nodes.
    /// Each seed is an arbitrary linear combination of the seven outputs.
    /// One batched transpose solve and one vortex-pair sweep serve all seeds.
    pub fn vjp(&self, seeds: &[AeroGradients]) -> Result<Vec<MeshNodeCotangent>> {
        if seeds
            .iter()
            .any(|s| COEFFICIENT_ORDER.iter().any(|&o| !s.get(o).is_finite()))
        {
            return Err(PegasusError::InvalidGeometry(
                "nonfinite coefficient seed".into(),
            ));
        }
        if seeds.is_empty() {
            return Ok(Vec::new());
        }
        let mesh = self.solver.mesh();
        let reference = self.solver.reference_quantities();
        let (sin, cos) = self.condition.alpha.sin_cos();
        let columns: Vec<_> = seeds
            .iter()
            .map(|s| {
                let c = s.coefficients;
                let force = Vec3::new(
                    -sin * c.lift + cos * c.drag,
                    c.side,
                    cos * c.lift + sin * c.drag,
                ) / reference.area;
                let moment =
                    Vec3::new(-c.roll, c.pitch, -c.yaw) / (reference.area * reference.chord);
                StripLoadCotangent {
                    force: vec![force; mesh.strip_count()],
                    moment: vec![moment; mesh.strip_count()],
                }
            })
            .collect();
        let (panel_seed, velocity_seed) = self.loads.panel_cotangents(&columns);
        let mut rhs = self
            .loads
            .circulation_gradients(&panel_seed, &velocity_seed);
        let mut bars = LatticeCotangents::zeros(mesh, seeds.len());
        self.loads
            .accumulate_explicit(&columns, &panel_seed, &mut bars);
        if seeds.iter().any(|s| s.induced_drag != 0.0) {
            let gamma = trefftz_circulation_gradients(
                mesh,
                self.solver.solution(),
                reference.area,
                0.01,
                0.1,
            )
            .induced_drag;
            let geometry = trefftz_geometry_vjp(
                mesh,
                self.solver.solution(),
                reference.area,
                0.01,
                0.1,
                [0.0, 0.0, 1.0],
            )?;
            for (col, seed) in seeds.iter().enumerate() {
                for i in 0..mesh.vortex_count() {
                    rhs[(i, col)] += seed.induced_drag * gamma[i];
                    for k in 0..3 {
                        bars.bound_left[(3 * i + k, col)] +=
                            seed.induced_drag * geometry.bound_left[i][k];
                        bars.bound_right[(3 * i + k, col)] +=
                            seed.induced_drag * geometry.bound_right[i][k];
                    }
                }
                for i in 0..mesh.strip_count() {
                    bars.strip_chord[(i, col)] += seed.induced_drag * geometry.strip_chord[i];
                }
            }
        }
        let adjoints = AdjointSolver::new(self.solver)?.solve_batched(&rhs)?;
        self.loads
            .accumulate_pairwise(&velocity_seed, &adjoints, &mut bars);
        self.loads.to_mesh_nodes(self.input, &bars)
    }

    /// True tangent solve: A dGamma = db - dA Gamma. Cost follows the number
    /// of supplied directions, independent of the number of mesh coordinates.
    /// Unselected outputs are zero. All directions share one kernel sweep and LU.
    pub fn jvp(
        &self,
        directions: &[MeshNodeDirection],
        selection: OutputSelection,
    ) -> Result<Vec<AeroGradients>> {
        let directions = directions
            .iter()
            .map(|d| self.lattice_direction(d))
            .collect::<Result<Vec<_>>>()?;
        if directions.is_empty() {
            return Ok(Vec::new());
        }
        if selection == OutputSelection::none() {
            return Ok(vec![AeroGradients::default(); directions.len()]);
        }
        let mesh = self.solver.mesh();
        let condition = self.condition;
        let reference = self.solver.reference_quantities();
        let solution = self.solver.solution();
        let (gamma, velocity) = self.tangent_state(&directions)?;
        let mut result = vec![AeroGradients::default(); directions.len()];
        // Output gradients and velocity seeds are independent of direction.
        for (output_index, output) in COEFFICIENT_ORDER
            .into_iter()
            .enumerate()
            .filter(|(_, o)| selection.contains(*o))
        {
            let circulation = if output == AeroOutput::InducedDrag {
                trefftz_circulation_gradients(mesh, solution, reference.area, 0.01, 0.1)
                    .induced_drag
            } else {
                force_circulation_gradient(self.solver, condition, output)
            };
            let velocity_seed = if output != AeroOutput::InducedDrag {
                force_induced_velocity_gradient(mesh, solution, condition, reference, output)
            } else {
                Vec::new()
            };
            for (k, d) in directions.iter().enumerate() {
                let explicit = if output == AeroOutput::InducedDrag {
                    trefftz_geometry_jvp(
                        mesh,
                        solution,
                        reference.area,
                        0.01,
                        0.1,
                        &d.bound_left,
                        &d.bound_right,
                        &d.strip_chord,
                    )?[2]
                } else {
                    force_geometry_jvp(
                        mesh,
                        solution,
                        condition,
                        reference,
                        &d.bound_left,
                        &d.bound_right,
                        &d.control_point,
                    )?[output_index]
                        + velocity_seed
                            .iter()
                            .zip(&velocity[k])
                            .map(|(a, b)| a.dot(b))
                            .sum::<f64>()
                };
                result[k].set(output, explicit + circulation.dot(&gamma.column(k)));
            }
        }
        Ok(result)
    }

    fn tangent_state(&self, directions: &[LatticeDirections]) -> Result<(MatX, Vec<Vec<Vec3>>)> {
        let mesh = self.solver.mesh();
        let condition = self.condition;
        let solution = self.solver.solution();
        let n = mesh.vortex_count();
        let pg = if self.solver.current_mach() < 1.0 {
            (1.0 - self.solver.current_mach().powi(2)).sqrt()
        } else {
            1.0
        };
        let freestream = VortexLatticeSolver::freestream_direction(condition);
        let omega = Vec3::new(
            condition.roll_rate,
            condition.pitch_rate,
            condition.yaw_rate,
        );
        let mut rhs = MatX::zeros(n, directions.len());
        let mut velocity = vec![vec![Vec3::zeros(); n]; directions.len()];
        for (i, target) in mesh.vortices().iter().enumerate() {
            let normal = target.normal_at_control_point;
            let local =
                freestream + omega.cross(&(target.control_point - self.solver.reference_point()));
            for (k, d) in directions.iter().enumerate() {
                if target.affected_by_freestream {
                    rhs[(i, k)] =
                        -omega.cross(&d.control_point[i]).dot(&normal) - local.dot(&d.normal[i]);
                }
            }
            for (j, source) in mesh.vortices().iter().enumerate() {
                let width = (source.bound_vortex_right - source.bound_vortex_left).norm();
                let core = if source.component_index == target.component_index {
                    0.0001 * width
                } else {
                    0.0001 * source.local_chord.max(width)
                };
                let jac = horseshoe_velocity_jacobians(
                    &target.control_point,
                    &source.bound_vortex_left,
                    &source.bound_vortex_right,
                    pg,
                    core,
                    i != j,
                );
                let v = horseshoe_vortex_velocity(
                    &target.control_point,
                    &source.bound_vortex_left,
                    &source.bound_vortex_right,
                    pg,
                    core,
                    i != j,
                );
                for (k, d) in directions.iter().enumerate() {
                    let core_direction = 0.0001
                        * if source.component_index != target.component_index
                            && source.local_chord >= width
                        {
                            d.local_chord[j]
                        } else {
                            (source.bound_vortex_right - source.bound_vortex_left)
                                .dot(&(d.bound_right[j] - d.bound_left[j]))
                                / width
                        };
                    let dv = jac.field_point * d.control_point[i]
                        + jac.bound_left * d.bound_left[j]
                        + jac.bound_right * d.bound_right[j]
                        + jac.core_radius * core_direction;
                    rhs[(i, k)] -=
                        solution.circulation[j] * (normal.dot(&dv) + v.dot(&d.normal[i]));
                    velocity[k][i] += solution.circulation[j] * dv;
                }
            }
        }
        let gamma = self.solver.solve_linearized_batched(&rhs)?;
        Ok((gamma, velocity))
    }

    /// Tangent strip force and moment resultants for load-rich coupled systems.
    /// Values have the same raw Kutta-Joukowski normalization as
    /// `StripLoadCotangent`: physical dynamic pressure is supplied by the caller.
    /// Moments are about `moment_points`; matching point directions include the
    /// derivative of moving beam stations. An empty point-direction batch means
    /// all moment points are fixed. All other batches must match node directions.
    pub fn strip_load_jvp(
        &self,
        directions: &[MeshNodeDirection],
        moment_points: &[Vec3],
        point_directions: &[Vec<Vec3>],
    ) -> Result<Vec<StripLoadCotangent>> {
        let mesh = self.solver.mesh();
        let ns = mesh.strip_count();
        if moment_points.len() != ns
            || moment_points
                .iter()
                .any(|v| v.iter().any(|x| !x.is_finite()))
            || (!point_directions.is_empty()
                && (point_directions.len() != directions.len()
                    || point_directions.iter().any(|d| {
                        d.len() != ns || d.iter().any(|v| v.iter().any(|x| !x.is_finite()))
                    })))
        {
            return Err(PegasusError::InvalidGeometry(
                "invalid strip moment-point directions".into(),
            ));
        }
        let directions = directions
            .iter()
            .map(|d| self.lattice_direction(d))
            .collect::<Result<Vec<_>>>()?;
        if directions.is_empty() {
            return Ok(Vec::new());
        }
        let (gamma, velocity) = self.tangent_state(&directions)?;
        let induced = self.solver.velocity_matrix() * &gamma;
        let solution = self.solver.solution();
        let freestream = VortexLatticeSolver::freestream_direction(self.condition);
        let mut result = vec![StripLoadCotangent::zeros(ns); directions.len()];
        for (i, vortex) in mesh.vortices().iter().enumerate() {
            let s = vortex.strip_index;
            let line = vortex.bound_vortex_right - vortex.bound_vortex_left;
            let v = freestream
                + Vec3::new(
                    solution.induced_velocity[(i, 0)],
                    solution.induced_velocity[(i, 1)],
                    solution.induced_velocity[(i, 2)],
                );
            let force = 2.0 * solution.circulation[i] * v.cross(&line);
            let arm = vortex.control_point - moment_points[s];
            for (k, d) in directions.iter().enumerate() {
                let dv = velocity[k][i]
                    + Vec3::new(
                        induced[(3 * i, k)],
                        induced[(3 * i + 1, k)],
                        induced[(3 * i + 2, k)],
                    );
                let df = 2.0
                    * (gamma[(i, k)] * v.cross(&line)
                        + solution.circulation[i]
                            * (dv.cross(&line) + v.cross(&(d.bound_right[i] - d.bound_left[i]))));
                let dp = if point_directions.is_empty() {
                    Vec3::zeros()
                } else {
                    point_directions[k][s]
                };
                result[k].force[s] += df;
                result[k].moment[s] += (d.control_point[i] - dp).cross(&force) + arm.cross(&df);
            }
        }
        Ok(result)
    }

    fn lattice_direction(&self, d: &MeshNodeDirection) -> Result<LatticeDirections> {
        if d.surfaces.len() != self.input.surfaces.len()
            || d.surfaces.iter().zip(&self.input.surfaces).any(|(d, s)| {
                d.len() != s.nodes.len() || d.iter().any(|v| v.iter().any(|x| !x.is_finite()))
            })
        {
            return Err(PegasusError::InvalidGeometry(
                "invalid mesh tangent layout or values".into(),
            ));
        }
        let mesh = self.solver.mesh();
        let mut out = LatticeDirections::zeros(mesh);
        let mut slot = 0;
        for (s, nodes) in self.input.surfaces.iter().zip(&d.surfaces) {
            let info = &mesh.surface_info()[slot];
            let index = |i: usize, j: usize| i * (s.ny + 1) + j;
            for j in 0..s.ny {
                let chord_direction: f64 = [j, j + 1]
                    .into_iter()
                    .map(|j| {
                        let chord = s.nodes[index(s.nx, j)] - s.nodes[index(0, j)];
                        0.5 * chord.dot(&(nodes[index(s.nx, j)] - nodes[index(0, j)]))
                            / chord.norm()
                    })
                    .sum();
                out.strip_chord[info.first_strip_index + j] = chord_direction;
                for i in 0..s.nx {
                    let ids = [
                        index(i, j),
                        index(i + 1, j),
                        index(i, j + 1),
                        index(i + 1, j + 1),
                    ];
                    let [p00, p10, p01, p11] = ids.map(|k| s.nodes[k]);
                    let [d00, d10, d01, d11] = ids.map(|k| nodes[k]);
                    let k = info.first_vortex_index + j * s.nx + i;
                    out.bound_left[k] = 0.75 * d00 + 0.25 * d10;
                    out.bound_right[k] = 0.75 * d01 + 0.25 * d11;
                    out.control_point[k] = 0.125 * d00 + 0.375 * d10 + 0.125 * d01 + 0.375 * d11;
                    let a = p11 - p00;
                    let b = p10 - p01;
                    out.normal[k] = normalized_direction(
                        &a.cross(&b),
                        &((d11 - d00).cross(&b) + a.cross(&(d10 - d01))),
                    );
                    out.local_chord[k] = chord_direction;
                }
            }
            if s.y_duplicated {
                let mirror = &mesh.surface_info()[slot + 1];
                for j in 0..s.ny {
                    out.strip_chord[mirror.first_strip_index + j] =
                        out.strip_chord[info.first_strip_index + j];
                    for i in 0..s.nx {
                        let a = info.first_vortex_index + j * s.nx + i;
                        let b = mirror.first_vortex_index + j * s.nx + i;
                        out.bound_left[b] = reflected_direction(out.bound_right[a]);
                        out.bound_right[b] = reflected_direction(out.bound_left[a]);
                        out.control_point[b] = reflected_direction(out.control_point[a]);
                        out.normal[b] = reflected_direction(out.normal[a]);
                        out.local_chord[b] = out.local_chord[a];
                    }
                }
            }
            slot += 1 + usize::from(s.y_duplicated);
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::super::geometry_adjoint::{GeometryAdjointContext, direct_mesh_node_direction};
    use super::*;
    use crate::geometry::SurfaceMesh;
    #[test]
    fn reverse_rows_match_independent_coordinate_adjoint_oracle() {
        let mut surface = SurfaceMesh::new("test", 1, 2);
        for i in 0..=1 {
            for j in 0..=2 {
                surface
                    .set_node(
                        i,
                        j,
                        Vec3::new(
                            i as f64 + 0.1 * j as f64,
                            j as f64,
                            0.07 * i as f64 * j as f64,
                        ),
                    )
                    .unwrap();
            }
        }
        let mut input = MeshInput::default();
        input.surfaces.push(surface);
        let mut solver = VortexLatticeSolver::new();
        solver.set_mesh(&input).unwrap();
        let condition = FlightCondition {
            alpha: 0.1,
            ..Default::default()
        };
        solver.solve(&condition).unwrap();
        let oracle =
            GeometryAdjointContext::new(&solver, &condition, OutputSelection::all()).unwrap();
        let seeds: Vec<_> = COEFFICIENT_ORDER
            .into_iter()
            .map(|o| {
                let mut s = AeroGradients::default();
                s.set(o, 1.0);
                s
            })
            .collect();
        let fast = MeshLinearization::new(&solver)
            .unwrap()
            .vjp(&seeds)
            .unwrap();
        for i in 0..=1 {
            for j in 0..=2 {
                for k in 0..3 {
                    let d = direct_mesh_node_direction(&input, solver.mesh(), 0, i, j, k).unwrap();
                    for (o, v) in oracle.evaluate(&d).unwrap() {
                        let row = COEFFICIENT_ORDER.iter().position(|&x| x == o).unwrap();
                        assert!((fast[row].surfaces[0][i * 3 + j][k] - v).abs() < 2e-12);
                    }
                }
            }
        }
    }
}
