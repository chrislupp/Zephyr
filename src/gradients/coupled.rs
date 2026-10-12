//! Reverse-mode sensitivities of per-strip loads to the lattice geometry.
//!
//! The aerodynamic coefficient adjoints behind [`GradientCalculator`] answer
//! "how does one scalar coefficient move when the mesh moves". An aeroelastic
//! coupling
//! asks something else: it needs the derivative of every *strip's* resultant
//! force and moment with respect to every mesh node, because those resultants
//! are what the structure feels.
//!
//! Rather than hand out that full tensor, this module takes the contraction the
//! caller wants. A [`StripLoadCotangent`] is one linear functional of the strip
//! loads; the adjoint returns the matching cotangent on the surface mesh nodes.
//! A coupled Jacobian is `ndof` such functionals, one per structural residual
//! row, and they are all carried through the same batched adjoint solve and the
//! same single sweep over vortex pairs.
//!
//! Everything here is reverse mode. There is no tangent solve and no forward
//! geometry JVP on this path — the batched transpose solve costs the same as a
//! forward one would, since the seed is `nvortex x ndof` either way.
//!
//! [`GradientCalculator`]: super::GradientCalculator

use crate::core::{MatX, Scalar, SymmetryType, Vec3};
use crate::error::{PegasusError, Result};
use crate::geometry::MeshInput;
use crate::solver::{
    FlightCondition, VortexLattice, VortexLatticeSolver, horseshoe_vortex_velocity,
};

use super::geometry_adjoint::{normalized_direction, reflected_direction};
use super::{AdjointSolver, TotalDerivative, horseshoe_velocity_jacobians};

/// One linear functional of a lattice's strip loads.
///
/// The two vectors are the cotangents of each strip's resultant force and of
/// its resultant moment, so the functional is
///
/// ```text
/// phi = sum_s ( force[s] . F_s + moment[s] . M_s )
/// ```
///
/// with `F_s` the sum of the strip's panel forces and `M_s` their moment about
/// that strip's moment point. Both are the raw Kutta–Joukowski sums, not the
/// area- and chord-normalized coefficients [`StripLoads`] reports, so a caller
/// working in physical units folds its dynamic pressure into the cotangents.
///
/// [`StripLoads`]: crate::aero::StripLoads
#[derive(Clone, Debug, Default)]
pub struct StripLoadCotangent {
    /// Cotangent of each strip's resultant force.
    pub force: Vec<Vec3>,
    /// Cotangent of each strip's resultant moment.
    pub moment: Vec<Vec3>,
}

impl StripLoadCotangent {
    /// A zeroed cotangent for a lattice with `strip_count` strips.
    #[must_use]
    pub fn zeros(strip_count: usize) -> Self {
        Self {
            force: vec![Vec3::zeros(); strip_count],
            moment: vec![Vec3::zeros(); strip_count],
        }
    }
}

/// Cotangents of a mesh's node coordinates.
#[derive(Clone, Debug, Default)]
pub struct MeshNodeCotangent {
    /// One entry per input surface, laid out like [`SurfaceMesh::nodes`] —
    /// row major, `i * (ny + 1) + j`.
    ///
    /// [`SurfaceMesh::nodes`]: crate::geometry::SurfaceMesh::nodes
    pub surfaces: Vec<Vec<Vec3>>,
}

/// Cotangents of the lattice's own geometry, three components per vortex.
///
/// Columns are the caller's functionals; rows are laid out `3 * vortex + d`,
/// matching the velocity-influence matrix.
struct LatticeCotangents {
    bound_left: MatX,
    bound_right: MatX,
    control_point: MatX,
    normal: MatX,
}

impl LatticeCotangents {
    fn zeros(vortex_count: usize, columns: usize) -> Self {
        Self {
            bound_left: MatX::zeros(3 * vortex_count, columns),
            bound_right: MatX::zeros(3 * vortex_count, columns),
            control_point: MatX::zeros(3 * vortex_count, columns),
            normal: MatX::zeros(3 * vortex_count, columns),
        }
    }
}

/// Reverse-mode adjoint of a solved lattice's per-strip loads.
///
/// Built against a converged solution, it turns cotangents of the strip
/// resultants into cotangents of the surface mesh nodes, carrying both the
/// explicit Kutta–Joukowski geometry dependence and the implicit one through the
/// circulation.
///
/// # Restrictions
///
/// The same ones the mesh-node coefficient gradients carry, for the same
/// reasons: symmetry planes and section drag polars each add a geometry path
/// this module does not represent, and returning a
/// derivative that silently omits one is worse than refusing. Mirror a wing
/// explicitly rather than through a symmetry plane.
pub struct StripLoadAdjoint<'a> {
    solver: &'a VortexLatticeSolver,
    condition: &'a FlightCondition,
    /// The point each strip's moment cotangent is taken about.
    moment_points: Vec<Vec3>,
    /// Total velocity at each control point, freestream plus induced.
    velocity: Vec<Vec3>,
    /// Kutta–Joukowski force on each panel.
    panel_force: Vec<Vec3>,
    prandtl_glauert: Scalar,
}

impl<'a> StripLoadAdjoint<'a> {
    /// Builds the adjoint against a solved lattice.
    ///
    /// # Arguments
    ///
    /// * `solver` — Solver holding a valid solution and its AIC factorization.
    /// * `condition` — The flight condition the solution was produced at.
    /// * `moment_points` — One point per strip, the moment cotangents are taken
    ///   about these.
    ///
    /// # Errors
    ///
    /// Returns [`PegasusError::InvalidSolution`] when the solution is stale, and
    /// [`PegasusError::InvalidGeometry`] or
    /// [`PegasusError::InvalidFlightCondition`] for the unsupported cases listed
    /// above or a moment-point count that does not match the lattice.
    pub fn new(
        solver: &'a VortexLatticeSolver,
        condition: &'a FlightCondition,
        moment_points: &[Vec3],
    ) -> Result<Self> {
        if !solver.solution().is_valid {
            return Err(PegasusError::InvalidSolution);
        }
        solver.check_solved_condition(condition)?;
        if solver.symmetry().y_symmetry != SymmetryType::None
            || solver.symmetry().z_symmetry != SymmetryType::None
        {
            return Err(PegasusError::InvalidGeometry(
                "strip load adjoints require explicit mirrored lattice geometry".into(),
            ));
        }
        if solver
            .mesh()
            .strips()
            .iter()
            .any(|strip| strip.has_viscous_data)
        {
            return Err(PegasusError::InvalidGeometry(
                "strip load adjoints do not yet support section drag polars; the strip \
                 area and lift-coefficient sensitivities they need are not represented"
                    .into(),
            ));
        }
        let mesh = solver.mesh();
        if moment_points.len() != mesh.strip_count() {
            return Err(PegasusError::InvalidGeometry(format!(
                "{} moment points for {} strips",
                moment_points.len(),
                mesh.strip_count()
            )));
        }

        let solution = solver.solution();
        let freestream = VortexLatticeSolver::freestream_direction(condition);
        let count = mesh.vortex_count();
        let mut velocity = Vec::with_capacity(count);
        let mut panel_force = Vec::with_capacity(count);
        for (index, vortex) in mesh.vortices().iter().enumerate() {
            let induced = if solution.induced_velocity.nrows() > index {
                Vec3::new(
                    solution.induced_velocity[(index, 0)],
                    solution.induced_velocity[(index, 1)],
                    solution.induced_velocity[(index, 2)],
                )
            } else {
                Vec3::zeros()
            };
            let total = freestream + induced;
            let line = vortex.bound_vortex_right - vortex.bound_vortex_left;
            velocity.push(total);
            panel_force.push(2.0 * solution.circulation[index] * total.cross(&line));
        }

        Ok(Self {
            solver,
            condition,
            moment_points: moment_points.to_vec(),
            velocity,
            panel_force,
            prandtl_glauert: prandtl_glauert_factor(solver.current_mach()),
        })
    }

    /// The lattice the adjoint was built against.
    #[must_use]
    pub fn mesh(&self) -> &VortexLattice {
        self.solver.mesh()
    }

    /// Cotangents of the surface mesh nodes, one per input functional.
    ///
    /// The result mirrors `input.surfaces`: entry `column` holds, for each
    /// surface, one [`Vec3`] per mesh node. A duplicated surface contributes to
    /// its original's nodes, since that is where its lattice came from.
    ///
    /// # Arguments
    ///
    /// * `input` — The mesh the lattice was built from.
    /// * `columns` — One strip-load cotangent per functional.
    ///
    /// # Errors
    ///
    /// Returns [`PegasusError::InvalidGeometry`] if a cotangent's length does not
    /// match the strip count or if the mesh does not describe this lattice, and
    /// propagates adjoint solve failures.
    pub fn mesh_node_cotangents(
        &self,
        input: &MeshInput,
        columns: &[StripLoadCotangent],
    ) -> Result<Vec<MeshNodeCotangent>> {
        if self.mesh().mesh_input() != Some(input) {
            return Err(PegasusError::InvalidMesh(
                "mesh input does not describe the solved lattice".into(),
            ));
        }
        let mesh = self.mesh();
        for column in columns {
            if column.force.len() != mesh.strip_count() || column.moment.len() != mesh.strip_count()
            {
                return Err(PegasusError::InvalidGeometry(
                    "strip load cotangent length does not match the solved lattice".into(),
                ));
            }
        }
        if columns.is_empty() {
            return Ok(Vec::new());
        }

        let (panel_seed, velocity_seed) = self.panel_cotangents(columns);
        let circulation_gradients = self.circulation_gradients(&panel_seed, &velocity_seed);
        let adjoints = AdjointSolver::new(self.solver)?.solve_batched(&circulation_gradients)?;

        let mut cotangents = LatticeCotangents::zeros(mesh.vortex_count(), columns.len());
        self.accumulate_explicit(columns, &panel_seed, &mut cotangents);
        self.accumulate_pairwise(&velocity_seed, &adjoints, &mut cotangents);
        self.to_mesh_nodes(input, &cotangents)
    }

    /// Derivatives of the same functionals with respect to the flight
    /// condition, one [`TotalDerivative`] per input functional.
    ///
    /// [`mesh_node_cotangents`](Self::mesh_node_cotangents) answers for the
    /// geometry; this answers for the angles. An aeroelastic coupling needs
    /// both, because the trim state is a design variable alongside the shape,
    /// and it needs them per *strip* rather than per coefficient — which is
    /// what [`GradientCalculator::compute_alpha_gradients`] gives and what a
    /// structure cannot use.
    ///
    /// The angles reach a strip's load twice, and both are here. The freestream
    /// is part of the total velocity the Kutta–Joukowski force is built from,
    /// which is explicit at frozen circulation. And it is the flow-tangency
    /// right-hand side, which moves the circulation itself; that term is the
    /// adjoint's, and it costs the same batched transpose solve
    /// `mesh_node_cotangents` already pays, with no sweep over vortex pairs
    /// after it.
    ///
    /// # Arguments
    ///
    /// * `columns` — One strip-load cotangent per functional.
    ///
    /// # Errors
    ///
    /// Returns [`PegasusError::InvalidGeometry`] if a cotangent's length does
    /// not match the strip count, and propagates adjoint solve failures.
    ///
    /// [`GradientCalculator::compute_alpha_gradients`]:
    ///     super::GradientCalculator::compute_alpha_gradients
    pub fn flight_condition_cotangents(
        &self,
        columns: &[StripLoadCotangent],
    ) -> Result<Vec<TotalDerivative>> {
        let mesh = self.mesh();
        for column in columns {
            if column.force.len() != mesh.strip_count() || column.moment.len() != mesh.strip_count()
            {
                return Err(PegasusError::InvalidGeometry(
                    "strip load cotangent length does not match the solved lattice".into(),
                ));
            }
        }
        if columns.is_empty() {
            return Ok(Vec::new());
        }

        let (panel_seed, velocity_seed) = self.panel_cotangents(columns);
        let circulation_gradients = self.circulation_gradients(&panel_seed, &velocity_seed);
        let adjoints = AdjointSolver::new(self.solver)?.solve_batched(&circulation_gradients)?;

        // The freestream is a unit direction; the dynamic pressure that
        // redimensionalizes it is the caller's, exactly as it is for the loads.
        let (sin_alpha, cos_alpha) = self.condition.alpha.sin_cos();
        let (sin_beta, cos_beta) = self.condition.beta.sin_cos();
        let by_alpha = Vec3::new(-sin_alpha * cos_beta, 0.0, cos_alpha * cos_beta);
        let by_beta = Vec3::new(-cos_alpha * sin_beta, cos_beta, -sin_alpha * sin_beta);

        let mut out = vec![TotalDerivative::default(); columns.len()];
        for (index, (vortex, normal)) in mesh
            .vortices()
            .iter()
            .zip(self.solver.linearized_normals())
            .enumerate()
        {
            for (column, total) in out.iter_mut().enumerate() {
                let seed = Vec3::new(
                    velocity_seed[(3 * index, column)],
                    velocity_seed[(3 * index + 1, column)],
                    velocity_seed[(3 * index + 2, column)],
                );
                // Explicit: the force reads the total velocity directly.
                total.by_alpha += seed.dot(&by_alpha);
                total.by_beta += seed.dot(&by_beta);

                // Implicit: `rhs = -n . V`, so the residual's own freestream
                // derivative carries the circulation's change with it. A panel
                // the freestream does not reach contributes nothing.
                if vortex.affected_by_freestream {
                    let lambda = adjoints[(index, column)];
                    total.by_alpha -= lambda * normal.dot(&by_alpha);
                    total.by_beta -= lambda * normal.dot(&by_beta);
                }
            }
        }
        Ok(out)
    }

    /// Derivatives of each strip-load functional with respect to every named
    /// command, in radians. Names are sorted like the lattice's control list.
    /// Includes induced velocity and uses the solved, deflected AIC.
    pub fn control_cotangents(
        &self,
        columns: &[StripLoadCotangent],
    ) -> Result<Vec<Vec<(String, Scalar)>>> {
        for column in columns {
            if column.force.len() != self.mesh().strip_count()
                || column.moment.len() != self.mesh().strip_count()
            {
                return Err(PegasusError::InvalidGeometry(
                    "strip load cotangent length does not match the solved lattice".into(),
                ));
            }
        }
        if columns.is_empty() {
            return Ok(Vec::new());
        }
        let (panel_seed, velocity_seed) = self.panel_cotangents(columns);
        let gradients = self.circulation_gradients(&panel_seed, &velocity_seed);
        let adjoints = AdjointSolver::new(self.solver)?.solve_batched(&gradients)?;
        let mut result = vec![Vec::new(); columns.len()];
        for name in self.mesh().control_surface_names() {
            let rhs = self
                .solver
                .control_residual_derivative(self.condition, name)?;
            for (column, values) in result.iter_mut().enumerate() {
                values.push((name.clone(), adjoints.column(column).dot(&rhs)));
            }
        }
        Ok(result)
    }

    // ========================================================================
    // Load model
    // ========================================================================

    /// Per-panel cotangents of the force and of the total velocity.
    ///
    /// A strip's resultant moment about `p_s` is `sum_i (cp_i - p_s) x F_i`, so
    /// contracting the functional strip by strip leaves one force cotangent per
    /// *panel*,
    ///
    /// ```text
    /// c_i = force[s] + moment[s] x (cp_i - p_s)
    /// ```
    ///
    /// and `phi = sum_i c_i . F_i`. The velocity cotangent follows from
    /// `F_i = 2 gamma_i V_i x l_i`, which is what carries the induced-velocity
    /// dependence into both the circulation gradient and the geometry sweep.
    ///
    /// Both are returned laid out `3 * vortex + d` by column.
    fn panel_cotangents(&self, columns: &[StripLoadCotangent]) -> (MatX, MatX) {
        let mesh = self.mesh();
        let solution = self.solver.solution();
        let count = mesh.vortex_count();
        let mut panel_seed = MatX::zeros(3 * count, columns.len());
        let mut velocity_seed = MatX::zeros(3 * count, columns.len());

        for (strip_index, strip) in mesh.strips().iter().enumerate() {
            let point = self.moment_points[strip_index];
            for offset in 0..strip.num_vortices {
                let index = strip.first_vortex_index + offset;
                let vortex = &mesh.vortices()[index];
                let arm = mesh.force_point(vortex) - point;
                let line = vortex.bound_vortex_right - vortex.bound_vortex_left;
                for (column, cotangent) in columns.iter().enumerate() {
                    let seed =
                        cotangent.force[strip_index] + cotangent.moment[strip_index].cross(&arm);
                    let velocity = 2.0 * solution.circulation[index] * line.cross(&seed);
                    for d in 0..3 {
                        panel_seed[(3 * index + d, column)] = seed[d];
                        velocity_seed[(3 * index + d, column)] = velocity[d];
                    }
                }
            }
        }
        (panel_seed, velocity_seed)
    }

    /// Total `d(phi)/d(circulation)`, including the induced-velocity path.
    ///
    /// The direct term is the panel force's own linearity in its circulation;
    /// the second is every *other* panel's circulation reaching this one through
    /// the velocity influence matrix. Mirrors
    /// [`force_circulation_gradient`](super::force_circulation_gradient), which
    /// does the same thing for a single aircraft-level coefficient.
    fn circulation_gradients(&self, panel_seed: &MatX, velocity_seed: &MatX) -> MatX {
        let mesh = self.mesh();
        let mut direct = MatX::zeros(mesh.vortex_count(), panel_seed.ncols());
        for (index, vortex) in mesh.vortices().iter().enumerate() {
            let line = vortex.bound_vortex_right - vortex.bound_vortex_left;
            let rate = 2.0 * self.velocity[index].cross(&line);
            for column in 0..panel_seed.ncols() {
                let seed = Vec3::new(
                    panel_seed[(3 * index, column)],
                    panel_seed[(3 * index + 1, column)],
                    panel_seed[(3 * index + 2, column)],
                );
                direct[(index, column)] = rate.dot(&seed);
            }
        }
        direct + self.solver.velocity_matrix().transpose() * velocity_seed
    }

    /// The geometry terms that need no sweep over vortex pairs.
    ///
    /// Two of them: the bound-vortex line the panel force is proportional to,
    /// and the moment arm from the strip's moment point to the panel's control
    /// point.
    fn accumulate_explicit(
        &self,
        columns: &[StripLoadCotangent],
        panel_seed: &MatX,
        out: &mut LatticeCotangents,
    ) {
        let mesh = self.mesh();
        let solution = self.solver.solution();
        for (strip_index, strip) in mesh.strips().iter().enumerate() {
            for offset in 0..strip.num_vortices {
                let index = strip.first_vortex_index + offset;
                for (column, cotangent) in columns.iter().enumerate() {
                    let seed = Vec3::new(
                        panel_seed[(3 * index, column)],
                        panel_seed[(3 * index + 1, column)],
                        panel_seed[(3 * index + 2, column)],
                    );
                    // `phi` picks up `c_i . 2 gamma_i (V_i x l_i)`, and
                    // `c . (V x l) = l . (c x V)`.
                    let line =
                        2.0 * solution.circulation[index] * seed.cross(&self.velocity[index]);
                    let arm = self.panel_force[index].cross(&cotangent.moment[strip_index]);
                    for d in 0..3 {
                        out.bound_right[(3 * index + d, column)] += line[d];
                        out.bound_left[(3 * index + d, column)] -= line[d];
                        if mesh.consistent_forces() {
                            out.bound_left[(3 * index + d, column)] += 0.5 * arm[d];
                            out.bound_right[(3 * index + d, column)] += 0.5 * arm[d];
                        } else {
                            out.control_point[(3 * index + d, column)] += arm[d];
                        }
                    }
                }
            }
        }
    }

    /// The geometry terms that couple every vortex to every other.
    ///
    /// One sweep serves both: the residual's dependence on geometry, contracted
    /// with the adjoints, and the induced velocity's, contracted with the
    /// velocity seeds. They share the horseshoe kernel — an AIC entry is just
    /// the induced velocity projected onto the target normal — so evaluating the
    /// dense Jacobian blocks once per pair covers both, for every column at
    /// once. Evaluating a seeded VJP per column instead would re-derive the same
    /// geometry `ndof` times per pair.
    fn accumulate_pairwise(
        &self,
        velocity_seed: &MatX,
        adjoints: &MatX,
        out: &mut LatticeCotangents,
    ) {
        let n = self.mesh().vortex_count();
        if !self.mesh().consistent_forces() || n < 128 {
            self.accumulate_pairwise_rows(velocity_seed, adjoints, out, 0..n);
            return;
        }
        // Fixed blocks and ordered reduction make the arithmetic independent
        // of the pool size (and of whether the parallel feature is enabled).
        // Bound scratch storage to 16 O(N * outputs) blocks, not O(N^2).
        let block_size = n.div_ceil(16);
        let blocks: Vec<_> = (0..n)
            .step_by(block_size)
            .map(|start| start..(start + block_size).min(n))
            .collect();
        let evaluate = |rows: &std::ops::Range<usize>| {
            let mut local = LatticeCotangents::zeros(n, adjoints.ncols());
            self.accumulate_pairwise_rows(velocity_seed, adjoints, &mut local, rows.clone());
            local
        };
        #[cfg(feature = "parallel")]
        let partials: Vec<_> = {
            use rayon::prelude::*;
            blocks.par_iter().map(evaluate).collect()
        };
        #[cfg(not(feature = "parallel"))]
        let partials: Vec<_> = blocks.iter().map(evaluate).collect();
        for local in partials {
            out.bound_left += local.bound_left;
            out.bound_right += local.bound_right;
            out.control_point += local.control_point;
            out.normal += local.normal;
        }
    }

    fn accumulate_pairwise_rows(
        &self,
        velocity_seed: &MatX,
        adjoints: &MatX,
        out: &mut LatticeCotangents,
        rows: std::ops::Range<usize>,
    ) {
        let mesh = self.mesh();
        let solution = self.solver.solution();
        let freestream = VortexLatticeSolver::freestream_direction(self.condition);
        let omega = Vec3::new(
            self.condition.roll_rate,
            self.condition.pitch_rate,
            self.condition.yaw_rate,
        );

        for row in rows {
            let target = &mesh.vortices()[row];
            let normal = self.solver.linearized_normals()[row];
            let seed_rows = velocity_seed.rows(3 * row, 3);
            let adjoint_row = adjoints.row(row);

            for (column, source) in mesh.vortices().iter().enumerate() {
                let gamma = solution.circulation[column];
                if gamma == 0.0 {
                    continue;
                }
                let core = core_radius(mesh, column, target.component_index);
                let include_bound = mesh.consistent_forces() || row != column;
                let jacobians = horseshoe_velocity_jacobians(
                    &target.control_point,
                    &source.bound_vortex_left,
                    &source.bound_vortex_right,
                    self.prandtl_glauert,
                    core,
                    include_bound,
                );

                // The AIC entry is `normal . velocity`, so its geometry gradient
                // is these same blocks contracted with the normal, plus the
                // velocity itself for the normal's own derivative. This is
                // exactly `aic_element_gradient`, reusing work the sweep has
                // already done.
                let gradient_left = jacobians.bound_left.tr_mul(&normal);
                let gradient_right = jacobians.bound_right.tr_mul(&normal);
                let gradient_control = jacobians.field_point.tr_mul(&normal);
                let gradient_normal = horseshoe_vortex_velocity(
                    &target.control_point,
                    &source.bound_vortex_left,
                    &source.bound_vortex_right,
                    self.prandtl_glauert,
                    core,
                    include_bound,
                );

                // Residual path. `R = A gamma - rhs`, and the coupled Jacobian
                // consumes `-psi^T dR/dg`, hence the negative weight.
                out.bound_left.rows_mut(3 * column, 3).gemm(
                    -gamma,
                    &gradient_left,
                    &adjoint_row,
                    1.0,
                );
                out.bound_right.rows_mut(3 * column, 3).gemm(
                    -gamma,
                    &gradient_right,
                    &adjoint_row,
                    1.0,
                );
                out.control_point.rows_mut(3 * row, 3).gemm(
                    -gamma,
                    &gradient_control,
                    &adjoint_row,
                    1.0,
                );
                out.normal
                    .rows_mut(3 * row, 3)
                    .gemm(-gamma, &gradient_normal, &adjoint_row, 1.0);

                // Load path: the induced velocity at this control point moves
                // when either the target or the source moves.
                let jacobians = if mesh.consistent_forces() {
                    horseshoe_velocity_jacobians(
                        &mesh.force_point(target),
                        &source.bound_vortex_left,
                        &source.bound_vortex_right,
                        self.prandtl_glauert,
                        core,
                        row != column,
                    )
                } else {
                    jacobians
                };
                out.bound_left.rows_mut(3 * column, 3).gemm_tr(
                    gamma,
                    &jacobians.bound_left,
                    &seed_rows,
                    1.0,
                );
                out.bound_right.rows_mut(3 * column, 3).gemm_tr(
                    gamma,
                    &jacobians.bound_right,
                    &seed_rows,
                    1.0,
                );
                if mesh.consistent_forces() {
                    out.bound_left.rows_mut(3 * row, 3).gemm_tr(
                        0.5 * gamma,
                        &jacobians.field_point,
                        &seed_rows,
                        1.0,
                    );
                    out.bound_right.rows_mut(3 * row, 3).gemm_tr(
                        0.5 * gamma,
                        &jacobians.field_point,
                        &seed_rows,
                        1.0,
                    );
                } else {
                    out.control_point.rows_mut(3 * row, 3).gemm_tr(
                        gamma,
                        &jacobians.field_point,
                        &seed_rows,
                        1.0,
                    );
                }
            }

            // The right-hand side is the freestream and rotation normalwash,
            // which depends on the target's control point and normal alone.
            if target.affected_by_freestream {
                let local_velocity = freestream
                    + omega.cross(&(target.control_point - self.solver.reference_point()));
                out.control_point.rows_mut(3 * row, 3).gemm(
                    -1.0,
                    &normal.cross(&omega),
                    &adjoint_row,
                    1.0,
                );
                out.normal
                    .rows_mut(3 * row, 3)
                    .gemm(-1.0, &local_velocity, &adjoint_row, 1.0);
            }
        }
    }

    // ========================================================================
    // Lattice to mesh
    // ========================================================================

    /// Pulls lattice cotangents back onto the mesh nodes they were built from.
    ///
    /// The transpose of the lattice construction: the bound vortex sits at the
    /// quarter chord of its panel's two chordwise edges, the control point at
    /// three quarters, and the normal comes from the panel's diagonals.
    fn to_mesh_nodes(
        &self,
        input: &MeshInput,
        cotangents: &LatticeCotangents,
    ) -> Result<Vec<MeshNodeCotangent>> {
        let lattice = self.mesh();
        let columns = cotangents.bound_left.ncols();
        let mut result = vec![
            MeshNodeCotangent {
                surfaces: input
                    .surfaces
                    .iter()
                    .map(|surface| vec![Vec3::zeros(); surface.node_count()])
                    .collect(),
            };
            columns
        ];

        let mut slot = 0usize;
        for (surface_index, surface) in input.surfaces.iter().enumerate() {
            let info = lattice
                .surface_info()
                .get(slot)
                .ok_or_else(|| PegasusError::InvalidMesh("surface lattice is missing".into()))?;
            if info.num_strips != surface.ny || info.num_chordwise != surface.nx {
                return Err(PegasusError::InvalidMesh(
                    "mesh input does not describe the solved lattice".into(),
                ));
            }

            for (column, per_column) in result.iter_mut().enumerate() {
                let nodes = &mut per_column.surfaces[surface_index];
                for strip_offset in 0..surface.ny {
                    for chord_offset in 0..surface.nx {
                        let index =
                            info.first_vortex_index + strip_offset * surface.nx + chord_offset;
                        let mut bound_left = read(&cotangents.bound_left, index, column);
                        let mut bound_right = read(&cotangents.bound_right, index, column);
                        let mut control_point = read(&cotangents.control_point, index, column);
                        let (mut normal, mut axes) = crate::solver::control::normal_vjp(
                            lattice.vortices()[index].normal_at_control_point,
                            read(&cotangents.normal, index, column),
                            &lattice.control_influences()[index],
                            self.condition,
                        );

                        // A duplicated surface has no nodes of its own; its
                        // lattice is a reflection of this one's, so its
                        // cotangents come back here with the reflection undone
                        // and the two bound endpoints swapped.
                        if surface.y_duplicated {
                            let duplicate =
                                lattice.surface_info().get(slot + 1).ok_or_else(|| {
                                    PegasusError::InvalidMesh("duplicate lattice is missing".into())
                                })?;
                            let mirrored = duplicate.first_vortex_index
                                + strip_offset * surface.nx
                                + chord_offset;
                            bound_left += reflected_direction(read(
                                &cotangents.bound_right,
                                mirrored,
                                column,
                            ));
                            bound_right +=
                                reflected_direction(read(&cotangents.bound_left, mirrored, column));
                            control_point += reflected_direction(read(
                                &cotangents.control_point,
                                mirrored,
                                column,
                            ));
                            let (mirror_normal, mirror_axes) = crate::solver::control::normal_vjp(
                                lattice.vortices()[mirrored].normal_at_control_point,
                                read(&cotangents.normal, mirrored, column),
                                &lattice.control_influences()[mirrored],
                                self.condition,
                            );
                            normal += reflected_direction(mirror_normal);
                            for (axis, mirrored_axis) in axes.iter_mut().zip(mirror_axes) {
                                *axis += reflected_direction(mirrored_axis);
                            }
                        }

                        // Hinges follow the incoming deformed mesh. Pull the
                        // normalized axis seed back onto its material nodes.
                        let mut influence = 0;
                        for control in &surface.controls {
                            if control.panel_gain(chord_offset, strip_offset).abs()
                                <= crate::core::EPSILON
                            {
                                continue;
                            }
                            let seed = normalized_direction(
                                &control.hinge_vector(surface, strip_offset),
                                &axes[influence],
                            );
                            for (j, sign) in [(strip_offset, -1.0), (strip_offset + 1, 1.0)] {
                                for (k, weight) in control.weights(surface, j) {
                                    nodes[k] += sign * weight * seed;
                                }
                            }
                            influence += 1;
                        }

                        let node = |i: usize, j: usize| -> usize { i * (surface.ny + 1) + j };
                        let n00 = node(chord_offset, strip_offset);
                        let n10 = node(chord_offset + 1, strip_offset);
                        let n01 = node(chord_offset, strip_offset + 1);
                        let n11 = node(chord_offset + 1, strip_offset + 1);

                        nodes[n00] += 0.75 * bound_left;
                        nodes[n10] += 0.25 * bound_left;
                        nodes[n01] += 0.75 * bound_right;
                        nodes[n11] += 0.25 * bound_right;

                        nodes[n00] += 0.125 * control_point;
                        nodes[n10] += 0.375 * control_point;
                        nodes[n01] += 0.125 * control_point;
                        nodes[n11] += 0.375 * control_point;

                        let position = |index: usize| surface.nodes[index];
                        let diagonal_a = position(n11) - position(n00);
                        let diagonal_b = position(n10) - position(n01);
                        let cross = diagonal_a.cross(&diagonal_b);
                        // Normalization is self-adjoint — it projects out the
                        // component along the unit vector and scales — so the
                        // forward helper is its own transpose here.
                        let bar_cross = normalized_direction(&cross, &normal);
                        let bar_a = diagonal_b.cross(&bar_cross);
                        let bar_b = bar_cross.cross(&diagonal_a);
                        nodes[n11] += bar_a;
                        nodes[n00] -= bar_a;
                        nodes[n10] += bar_b;
                        nodes[n01] -= bar_b;
                    }
                }
            }
            slot += 1 + usize::from(surface.y_duplicated);
        }
        Ok(result)
    }
}

/// Reads one vortex's three cotangent components out of a stacked matrix.
fn read(matrix: &MatX, vortex: usize, column: usize) -> Vec3 {
    Vec3::new(
        matrix[(3 * vortex, column)],
        matrix[(3 * vortex + 1, column)],
        matrix[(3 * vortex + 2, column)],
    )
}

fn prandtl_glauert_factor(mach: Scalar) -> Scalar {
    if mach < 1.0 {
        (1.0 - mach * mach).sqrt()
    } else {
        1.0
    }
}

/// The Leishman core radius the influence matrix builder uses, at its defaults.
fn core_radius(mesh: &VortexLattice, source_index: usize, target_component: usize) -> Scalar {
    let source = &mesh.vortices()[source_index];
    let width = (source.bound_vortex_right - source.bound_vortex_left).norm();
    if source.component_index == target_component {
        0.0001 * width
    } else {
        (0.0001 * source.local_chord).max(0.0001 * width)
    }
}

#[cfg(test)]
mod tests {
    use super::super::aic_element_gradient;
    use super::*;

    /// The pairwise sweep derives the AIC entry's geometry gradient from the
    /// dense velocity blocks rather than calling `aic_element_gradient`. This
    /// pins the two to each other, so a change to either kernel cannot silently
    /// split the coupled Jacobian off from the aerodynamic adjoints.
    #[test]
    fn the_swept_blocks_reproduce_the_aic_element_gradient() {
        let bound_left = Vec3::new(0.1, -0.4, 0.05);
        let bound_right = Vec3::new(0.15, 0.6, -0.02);
        let control_point = Vec3::new(0.9, 0.2, 0.3);
        let normal = Vec3::new(0.1, -0.2, 0.97).normalize();
        let core = 1.0e-4;
        let pg = 0.9;

        let reference = aic_element_gradient(
            &bound_left,
            &bound_right,
            &control_point,
            &normal,
            core,
            pg,
            true,
        );
        let jacobians =
            horseshoe_velocity_jacobians(&control_point, &bound_left, &bound_right, pg, core, true);

        for (swept, expected) in [
            (jacobians.bound_left.tr_mul(&normal), reference.bound_left),
            (jacobians.bound_right.tr_mul(&normal), reference.bound_right),
            (
                jacobians.field_point.tr_mul(&normal),
                reference.control_point,
            ),
        ] {
            assert!(
                (swept - expected).norm() <= 1.0e-12 * (1.0 + expected.norm()),
                "swept block {swept:?} does not match {expected:?}"
            );
        }
    }
}
