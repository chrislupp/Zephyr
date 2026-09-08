//! Cached primal vortex-lattice solver.

use faer::{
    MatMut, MatRef,
    linalg::solvers::{PartialPivLu, Solve},
};

use crate::core::{EPSILON, MatX, Scalar, Vec3, VecX};
use crate::error::{PegasusError, Result};
use crate::geometry::{Aircraft, MeshInput, ReferenceQuantities};

use super::{
    FlightCondition, InfluenceMatrixBuilder, SolutionState, SymmetryConfig, VortexLattice,
};

type DynamicLu = PartialPivLu<Scalar>;

fn factorize(matrix: &MatX) -> DynamicLu {
    MatRef::from_column_major_slice(matrix.as_slice(), matrix.nrows(), matrix.ncols())
        .partial_piv_lu()
}

fn solve_matrix(factorization: &DynamicLu, rhs: &MatX) -> Result<MatX> {
    let mut solution = rhs.clone();
    let (nrows, ncols) = solution.shape();
    factorization.solve_in_place(MatMut::from_column_major_slice_mut(
        solution.as_mut_slice(),
        nrows,
        ncols,
    ));
    if solution.iter().all(|value| value.is_finite()) {
        Ok(solution)
    } else {
        Err(PegasusError::SingularSystem)
    }
}

fn solve_vector(factorization: &DynamicLu, rhs: &VecX) -> Result<VecX> {
    let mut solution = rhs.clone();
    let nrows = solution.nrows();
    factorization.solve_in_place(MatMut::from_column_major_slice_mut(
        solution.as_mut_slice(),
        nrows,
        1,
    ));
    if solution.iter().all(|value| value.is_finite()) {
        Ok(solution)
    } else {
        Err(PegasusError::SingularSystem)
    }
}

fn solve_transpose_matrix(factorization: &DynamicLu, rhs: &MatX) -> Result<MatX> {
    let mut solution = rhs.clone();
    let (nrows, ncols) = solution.shape();
    factorization.solve_transpose_in_place(MatMut::from_column_major_slice_mut(
        solution.as_mut_slice(),
        nrows,
        ncols,
    ));
    if solution.iter().all(|value| value.is_finite()) {
        Ok(solution)
    } else {
        Err(PegasusError::SingularSystem)
    }
}

fn solve_transpose_vector(factorization: &DynamicLu, rhs: &VecX) -> Result<VecX> {
    let mut solution = rhs.clone();
    let nrows = solution.nrows();
    factorization.solve_transpose_in_place(MatMut::from_column_major_slice_mut(
        solution.as_mut_slice(),
        nrows,
        1,
    ));
    if solution.iter().all(|value| value.is_finite()) {
        Ok(solution)
    } else {
        Err(PegasusError::SingularSystem)
    }
}

/// Main VLM solver with cached AIC and six unit responses.
#[derive(Clone, Debug)]
pub struct VortexLatticeSolver {
    mesh: VortexLattice,
    symmetry: SymmetryConfig,
    mirror_reuse: bool,
    reference_point: Vec3,
    reference_quantities: ReferenceQuantities,
    aic_matrix: MatX,
    velocity_influence_matrix: MatX,
    control_influence_matrix: MatX,
    solution_aic_matrix: Option<MatX>,
    aic_lu: Option<DynamicLu>,
    solution_lu: Option<DynamicLu>,
    solution_normals: Vec<Vec3>,
    aic_valid: bool,
    unit_circulation: MatX,
    unit_responses_valid: bool,
    solution: SolutionState,
    current_mach: Scalar,
}

impl Default for VortexLatticeSolver {
    fn default() -> Self {
        Self {
            mesh: VortexLattice::default(),
            symmetry: SymmetryConfig::default(),
            mirror_reuse: false,
            reference_point: Vec3::zeros(),
            reference_quantities: ReferenceQuantities::default(),
            aic_matrix: MatX::zeros(0, 0),
            velocity_influence_matrix: MatX::zeros(0, 0),
            control_influence_matrix: MatX::zeros(0, 0),
            solution_aic_matrix: None,
            aic_lu: None,
            solution_lu: None,
            solution_normals: Vec::new(),
            aic_valid: false,
            unit_circulation: MatX::zeros(0, 6),
            unit_responses_valid: false,
            solution: SolutionState::default(),
            current_mach: 0.0,
        }
    }
}

impl VortexLatticeSolver {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_aircraft(&mut self, aircraft: &Aircraft) -> Result<()> {
        self.mesh.generate_from_aircraft(aircraft)?;
        self.reference_point = *aircraft.moment_reference_point();
        self.reference_quantities = ReferenceQuantities::from(aircraft);
        self.current_mach = aircraft.mach_number();
        self.invalidate();
        Ok(())
    }

    pub fn set_mesh(&mut self, input: &MeshInput) -> Result<()> {
        self.mesh.generate_from_mesh(input)?;
        self.reference_quantities = input.reference.clone();
        self.reference_point = input.reference.moment_reference_point;
        self.invalidate();
        Ok(())
    }

    pub fn set_symmetry(&mut self, symmetry: SymmetryConfig) {
        self.symmetry = symmetry;
        self.invalidate();
    }

    /// Lets influence assembly derive mirrored entries by reflection rather
    /// than evaluating them, for lattices built from `y_duplicated` surfaces.
    ///
    /// Off by default. Results are unaffected: assembly takes the fast path
    /// only when it can prove the mesh satisfies the reflection identity
    /// exactly, and the derived entries are bit-identical to evaluated ones.
    /// Meshes that do not qualify -- symmetry images in use, partially
    /// duplicated geometry, or a mirror plane that is not a pure sign flip --
    /// fall back to evaluating every entry.
    pub fn set_mirror_reuse(&mut self, enabled: bool) {
        self.mirror_reuse = enabled;
        self.invalidate();
    }

    #[must_use]
    pub fn mirror_reuse(&self) -> bool {
        self.mirror_reuse
    }

    pub fn set_reference_point(&mut self, point: Vec3) {
        self.reference_point = point;
        self.unit_responses_valid = false;
    }

    /// Apply a scoped mesh mutation and invalidate all dependent caches.
    pub fn mutate_mesh(&mut self, mutation: impl FnOnce(&mut VortexLattice)) {
        mutation(&mut self.mesh);
        self.invalidate();
    }

    fn invalidate(&mut self) {
        self.aic_valid = false;
        self.unit_responses_valid = false;
        self.solution_aic_matrix = None;
        self.aic_lu = None;
        self.solution_lu = None;
        self.solution_normals.clear();
        self.solution.is_valid = false;
    }

    pub fn rebuild_aic(&mut self) -> Result<()> {
        if self.mesh.vortex_count() == 0 {
            self.aic_valid = false;
            return Err(PegasusError::InvalidMesh(
                "cannot build an AIC for an empty lattice".into(),
            ));
        }
        let mut builder = InfluenceMatrixBuilder::new(&self.mesh);
        builder.set_mach_number(self.current_mach);
        builder.set_symmetry(self.symmetry);
        builder.set_mirror_reuse(self.mirror_reuse);
        // The finiteness flag is accumulated during assembly, so rejecting a
        // singular system no longer costs a second full pass over the matrix
        // (10.6 MB at 1152 panels).
        // Assemble into locals and publish only once the result is known to be
        // finite. Writing the fields first left `aic()`, `velocity_matrix()`,
        // `linearized_aic()` and `solution_aic()` handing out NaN-laden
        // matrices after a failed rebuild, paired with a stale `aic_lu` from
        // whatever preceded them.
        let (aic_matrix, control_influence_matrix, velocity_influence_matrix, all_finite) =
            builder.normalwash_control_and_velocity_checked();
        if !all_finite {
            self.aic_valid = false;
            self.unit_responses_valid = false;
            self.aic_lu = None;
            return Err(PegasusError::SingularSystem);
        }
        self.aic_matrix = aic_matrix;
        self.control_influence_matrix = control_influence_matrix;
        self.velocity_influence_matrix = velocity_influence_matrix;
        self.aic_lu = Some(factorize(&self.aic_matrix));
        self.aic_valid = true;
        self.unit_responses_valid = false;
        Ok(())
    }

    fn unit_freestream_rhs(&self, component: usize) -> VecX {
        VecX::from_iterator(
            self.mesh.vortex_count(),
            self.mesh.vortices().iter().map(|vortex| {
                if vortex.affected_by_freestream {
                    -vortex.normal_at_control_point[component]
                } else {
                    0.0
                }
            }),
        )
    }

    fn unit_rotation_rhs(&self, component: usize) -> VecX {
        let mut omega = Vec3::zeros();
        omega[component] = 1.0;
        VecX::from_iterator(
            self.mesh.vortex_count(),
            self.mesh.vortices().iter().map(|vortex| {
                if vortex.affected_by_freestream {
                    -omega
                        .cross(&(vortex.control_point - self.reference_point))
                        .dot(&vortex.normal_at_control_point)
                } else {
                    0.0
                }
            }),
        )
    }

    pub fn solve_unit_responses(&mut self) -> Result<()> {
        if !self.aic_valid {
            self.rebuild_aic()?;
        }
        let n = self.mesh.vortex_count();
        let mut rhs = MatX::zeros(n, 6);
        for component in 0..3 {
            rhs.set_column(component, &self.unit_freestream_rhs(component));
            rhs.set_column(component + 3, &self.unit_rotation_rhs(component));
        }
        // No residual check here: `rebuild_aic` already rejects a non-finite
        // AIC and `solve_matrix` already rejects a non-finite solution, which
        // is what a singular system yields. Forming `A x - b` only to test it
        // for finiteness cost a full N²×6 product and an N×6 allocation
        // without testing the residual's magnitude, so it never established
        // any accuracy property the other two checks did not already cover.
        self.unit_circulation = solve_matrix(
            self.aic_lu.as_ref().ok_or(PegasusError::SingularSystem)?,
            &rhs,
        )?;
        self.unit_responses_valid = true;
        Ok(())
    }

    #[must_use]
    pub fn freestream_direction(condition: &FlightCondition) -> Vec3 {
        let (sin_alpha, cos_alpha) = condition.alpha.sin_cos();
        let (sin_beta, cos_beta) = condition.beta.sin_cos();
        Vec3::new(cos_alpha * cos_beta, sin_beta, sin_alpha * cos_beta)
    }

    fn combined_circulation(&self, condition: &FlightCondition) -> VecX {
        let velocity = Self::freestream_direction(condition);
        self.unit_circulation.column(0) * velocity.x
            + self.unit_circulation.column(1) * velocity.y
            + self.unit_circulation.column(2) * velocity.z
            + self.unit_circulation.column(3) * condition.roll_rate
            + self.unit_circulation.column(4) * condition.pitch_rate
            + self.unit_circulation.column(5) * condition.yaw_rate
    }

    fn rotate_about_axis(vector: Vec3, axis: &Vec3, angle: Scalar) -> Vec3 {
        let (sin, cos) = angle.sin_cos();
        vector * cos + axis.cross(&vector) * sin + axis * axis.dot(&vector) * (1.0 - cos)
    }

    pub(crate) fn effective_normals(&self, condition: &FlightCondition) -> Vec<Vec3> {
        self.mesh
            .vortices()
            .iter()
            .zip(self.mesh.control_influences())
            .map(|(vortex, influences)| {
                let mut normal = vortex.normal_at_control_point;
                // Entries are intentionally processed in caller-provided order.
                // Repeated names therefore compose repeated rotations.
                for (name, deflection) in &condition.control_deflections {
                    for influence in influences.iter().filter(|item| item.name == *name) {
                        if influence.hinge_axis.norm_squared() > EPSILON * EPSILON {
                            normal = Self::rotate_about_axis(
                                normal,
                                &influence.hinge_axis.normalize(),
                                deflection * influence.gain,
                            );
                        }
                    }
                }
                normal
            })
            .collect()
    }

    pub(crate) fn effective_normals_with_derivative(
        &self,
        condition: &FlightCondition,
        control_name: &str,
    ) -> (Vec<Vec3>, Vec<Vec3>) {
        self.mesh
            .vortices()
            .iter()
            .zip(self.mesh.control_influences())
            .map(|(vortex, influences)| {
                let mut normal = vortex.normal_at_control_point;
                let mut derivative = Vec3::zeros();
                let mut target_present = false;
                for (name, deflection) in &condition.control_deflections {
                    target_present |= name == control_name;
                    for influence in influences.iter().filter(|item| item.name == *name) {
                        if influence.hinge_axis.norm_squared() <= EPSILON * EPSILON {
                            continue;
                        }
                        let axis = influence.hinge_axis.normalize();
                        let angle = deflection * influence.gain;
                        normal = Self::rotate_about_axis(normal, &axis, angle);
                        derivative = Self::rotate_about_axis(derivative, &axis, angle);
                        if name == control_name {
                            derivative += influence.gain * axis.cross(&normal);
                        }
                    }
                }
                if !target_present {
                    for influence in influences.iter().filter(|item| item.name == control_name) {
                        if influence.hinge_axis.norm_squared() > EPSILON * EPSILON {
                            let axis = influence.hinge_axis.normalize();
                            derivative += influence.gain * axis.cross(&normal);
                        }
                    }
                }
                (normal, derivative)
            })
            .unzip()
    }

    fn matrix_for_normals(&self, normals: &[Vec3]) -> MatX {
        let n = self.mesh.vortex_count();
        MatX::from_fn(n, n, |row, column| {
            normals[row].dot(&Vec3::new(
                self.control_influence_matrix[(3 * row, column)],
                self.control_influence_matrix[(3 * row + 1, column)],
                self.control_influence_matrix[(3 * row + 2, column)],
            ))
        })
    }

    pub(crate) fn rhs_for_normals(&self, condition: &FlightCondition, normals: &[Vec3]) -> VecX {
        let freestream = Self::freestream_direction(condition);
        let omega = Vec3::new(
            condition.roll_rate,
            condition.pitch_rate,
            condition.yaw_rate,
        );
        VecX::from_iterator(
            self.mesh.vortex_count(),
            self.mesh
                .vortices()
                .iter()
                .zip(normals)
                .map(|(vortex, normal)| {
                    if vortex.affected_by_freestream {
                        -(freestream + omega.cross(&(vortex.control_point - self.reference_point)))
                            .dot(normal)
                    } else {
                        0.0
                    }
                }),
        )
    }

    fn unit_responses_for_normals(
        &self,
        factorization: &DynamicLu,
        normals: &[Vec3],
    ) -> Result<MatX> {
        let n = self.mesh.vortex_count();
        let mut rhs = MatX::zeros(n, 6);
        for (row, (vortex, normal)) in self.mesh.vortices().iter().zip(normals).enumerate() {
            if !vortex.affected_by_freestream {
                continue;
            }
            for component in 0..3 {
                rhs[(row, component)] = -normal[component];
                let mut omega = Vec3::zeros();
                omega[component] = 1.0;
                rhs[(row, component + 3)] = -omega
                    .cross(&(vortex.control_point - self.reference_point))
                    .dot(normal);
            }
        }
        solve_matrix(factorization, &rhs)
    }

    /// Analytical circulation derivative for one named control at the current
    /// flight condition.
    pub fn control_circulation_derivative(
        &self,
        condition: &FlightCondition,
        control_name: &str,
    ) -> Result<VecX> {
        if !self.solution.is_valid {
            return Err(PegasusError::InvalidSolution);
        }
        let (normals, normal_derivatives) =
            self.effective_normals_with_derivative(condition, control_name);
        let n = self.mesh.vortex_count();

        // `matrix_derivative` was materialized as a dense N×N only to be used
        // once, in `matrix_derivative * circulation`. Contracting as we go
        // keeps the same O(N²) work without the allocation.
        let circulation = &self.solution.circulation;
        let mut matrix_derivative_times_circulation = VecX::zeros(n);
        for row in 0..n {
            let mut sum = 0.0;
            for column in 0..n {
                let entry = normal_derivatives[row].dot(&Vec3::new(
                    self.control_influence_matrix[(3 * row, column)],
                    self.control_influence_matrix[(3 * row + 1, column)],
                    self.control_influence_matrix[(3 * row + 2, column)],
                ));
                sum += entry * circulation[column];
            }
            matrix_derivative_times_circulation[row] = sum;
        }
        let freestream = Self::freestream_direction(condition);
        let omega = Vec3::new(
            condition.roll_rate,
            condition.pitch_rate,
            condition.yaw_rate,
        );
        let rhs_derivative = VecX::from_iterator(
            n,
            self.mesh
                .vortices()
                .iter()
                .zip(normal_derivatives)
                .map(|(vortex, derivative)| {
                    if vortex.affected_by_freestream {
                        -(freestream + omega.cross(&(vortex.control_point - self.reference_point)))
                            .dot(&derivative)
                    } else {
                        0.0
                    }
                }),
        );
        let rhs = rhs_derivative - matrix_derivative_times_circulation;

        // The system here is the one the primal already solved: this function
        // requires a valid solution and reads its circulation, so it is only
        // meaningful at the solved condition, where `effective_normals_with_derivative`
        // reproduces exactly the normals `solve` used. When that holds the
        // cached factorization applies and a fresh O(N³) decomposition is pure
        // waste. The normals are compared rather than assumed, so a caller who
        // reaches here with a stale cache still gets a correct answer.
        if self.solution_normals == normals
            && let Some(lu) = self.solution_lu.as_ref().or(self.aic_lu.as_ref())
        {
            return solve_vector(lu, &rhs);
        }
        let matrix = self.matrix_for_normals(&normals);
        solve_vector(&factorize(&matrix), &rhs)
    }

    /// Internal selective path. Its solver is owned exclusively by `Analysis`,
    /// so partial response availability never escapes through the legacy API.
    pub(crate) fn solve_selective(
        &mut self,
        condition: &FlightCondition,
        half: bool,
        velocity: bool,
        reuse_solution: bool,
        work: &mut crate::analysis::Work,
    ) -> Result<()> {
        let normals = self.effective_normals(condition);
        let n = self.mesh.vortex_count();
        let operator_changed = !self.aic_valid || self.current_mach != condition.mach;
        let matrix_changed = operator_changed || self.solution_normals != normals;
        if matrix_changed {
            self.solution.is_valid = false;
            self.current_mach = condition.mach;
            let (matrix, control, velocities) =
                if !operator_changed && self.control_influence_matrix.nrows() == 3 * n {
                    // A control changes the normal projection, not the velocity
                    // operator. Reuse its coefficients when they are available.
                    (
                        self.matrix_for_normals(&normals),
                        std::mem::take(&mut self.control_influence_matrix),
                        std::mem::take(&mut self.velocity_influence_matrix),
                    )
                } else {
                    crate::analysis::operator::assemble(
                        &self.mesh,
                        condition.mach,
                        &normals,
                        half,
                        true,
                        velocity,
                        work,
                    )
                };
            if !matrix.iter().all(|v| v.is_finite()) {
                self.invalidate();
                return Err(PegasusError::SingularSystem);
            }
            self.aic_matrix = matrix;
            self.control_influence_matrix = control;
            self.velocity_influence_matrix = velocities;
            self.aic_lu = Some(factorize(&self.aic_matrix));
            self.solution_lu = None;
            self.solution_aic_matrix = None;
            self.aic_valid = true;
            self.solution_normals = normals;
            work.factorizations += 1;
        } else if velocity && self.velocity_influence_matrix.nrows() == 0 {
            self.velocity_influence_matrix = crate::analysis::operator::assemble(
                &self.mesh,
                condition.mach,
                &normals,
                half,
                false,
                true,
                work,
            )
            .2;
        }
        if matrix_changed || !reuse_solution || !self.solution.is_valid {
            let rhs = self.rhs_for_normals(condition, &self.solution_normals);
            let gamma = solve_vector(
                self.aic_lu.as_ref().ok_or(PegasusError::InvalidSolution)?,
                &rhs,
            )?;
            self.solution = SolutionState {
                circulation: gamma,
                is_valid: true,
                ..SolutionState::default()
            };
            work.primal_solves += 1;
        }
        if velocity && self.solution.induced_velocity.nrows() == 0 {
            let flat = &self.velocity_influence_matrix * &self.solution.circulation;
            self.solution.induced_velocity = MatX::from_fn(n, 3, |i, j| flat[3 * i + j]);
            work.velocity_products += 1;
        }
        Ok(())
    }

    pub fn solve(&mut self, condition: &FlightCondition) -> Result<SolutionState> {
        if self.mesh.vortex_count() == 0 {
            return Err(PegasusError::InvalidMesh(
                "cannot solve an empty lattice".into(),
            ));
        }
        // Every scalar that reaches the right-hand side is checked, including
        // the three body rates. They were previously unchecked, so a NaN rate
        // produced an all-NaN circulation reported as a valid solution.
        //
        // `dynamic_pressure` is checked for finiteness before its sign: every
        // ordered comparison against NaN is false, so `x < 0.0` alone let NaN
        // straight through.
        if !condition.alpha.is_finite()
            || !condition.beta.is_finite()
            || !condition.mach.is_finite()
            || !condition.roll_rate.is_finite()
            || !condition.pitch_rate.is_finite()
            || !condition.yaw_rate.is_finite()
            || condition
                .control_deflections
                .iter()
                .any(|(_, value)| !value.is_finite())
            || condition.mach < 0.0
            || !condition.dynamic_pressure.is_finite()
            || condition.dynamic_pressure < 0.0
        {
            return Err(PegasusError::InvalidFlightCondition(
                "angles, body rates, Mach, and dynamic pressure must be finite; \
                 Mach and dynamic pressure must be nonnegative"
                    .into(),
            ));
        }
        // Prandtl-Glauert has no meaning at or above Mach 1. Returning the
        // incompressible answer silently, as clamping beta to 1 did, hands back
        // plausible numbers for a condition the method cannot model.
        if condition.mach >= 1.0 {
            return Err(PegasusError::InvalidFlightCondition(format!(
                "Mach {} is not subsonic; the vortex-lattice method requires Mach < 1",
                condition.mach
            )));
        }
        if (condition.mach - self.current_mach).abs() > EPSILON {
            self.current_mach = condition.mach;
            self.invalidate();
        }
        if !self.unit_responses_valid {
            self.solve_unit_responses()?;
        }
        let (circulation, solution_units) = if condition.control_deflections.is_empty() {
            self.solution_aic_matrix = None;
            self.solution_lu = None;
            self.solution_normals = self
                .mesh
                .vortices()
                .iter()
                .map(|vortex| vortex.normal_at_control_point)
                .collect();
            (
                self.combined_circulation(condition),
                self.unit_circulation.clone(),
            )
        } else {
            let normals = self.effective_normals(condition);
            let matrix = self.matrix_for_normals(&normals);
            let rhs = self.rhs_for_normals(condition, &normals);
            let factorization = factorize(&matrix);
            let circulation = solve_vector(&factorization, &rhs)?;
            let units = self.unit_responses_for_normals(&factorization, &normals)?;
            self.solution_lu = Some(factorization);
            self.solution_aic_matrix = Some(matrix);
            self.solution_normals = normals;
            (circulation, units)
        };
        let n = self.mesh.vortex_count();
        let flat_velocity = &self.velocity_influence_matrix * &circulation;
        let mut induced = MatX::zeros(n, 3);
        for i in 0..n {
            induced[(i, 0)] = flat_velocity[3 * i];
            induced[(i, 1)] = flat_velocity[3 * i + 1];
            induced[(i, 2)] = flat_velocity[3 * i + 2];
        }
        for (vortex, gamma) in self.mesh.vortices_mut().iter_mut().zip(circulation.iter()) {
            vortex.circulation = *gamma;
        }
        let state = SolutionState {
            circulation,
            circulation_per_u: solution_units.column(0).into_owned(),
            circulation_per_v: solution_units.column(1).into_owned(),
            circulation_per_w: solution_units.column(2).into_owned(),
            circulation_per_p: solution_units.column(3).into_owned(),
            circulation_per_q: solution_units.column(4).into_owned(),
            circulation_per_r: solution_units.column(5).into_owned(),
            induced_velocity: induced,
            is_valid: true,
        };
        self.solution = state.clone();
        Ok(state)
    }

    #[must_use]
    pub fn mesh(&self) -> &VortexLattice {
        &self.mesh
    }

    #[must_use]
    pub fn symmetry(&self) -> SymmetryConfig {
        self.symmetry
    }

    #[must_use]
    pub fn reference_point(&self) -> &Vec3 {
        &self.reference_point
    }

    #[must_use]
    pub fn reference_quantities(&self) -> &ReferenceQuantities {
        &self.reference_quantities
    }

    #[must_use]
    pub fn aic(&self) -> &MatX {
        &self.aic_matrix
    }

    /// Linear system used by the current solution, including control-normal
    /// rotations. Before the first solve this is empty.
    #[must_use]
    pub fn solution_aic(&self) -> &MatX {
        self.solution_aic_matrix
            .as_ref()
            .unwrap_or(&self.aic_matrix)
    }

    /// AIC to use for linearization: the current control-adjusted matrix when
    /// available, otherwise the cached undeflected matrix.
    #[must_use]
    pub fn linearized_aic(&self) -> &MatX {
        self.solution_aic_matrix
            .as_ref()
            .unwrap_or(&self.aic_matrix)
    }

    pub(crate) fn linearized_normals(&self) -> &[Vec3] {
        &self.solution_normals
    }

    pub(crate) fn solve_linearized_transpose(&self, rhs: &VecX) -> Result<VecX> {
        let factorization = self
            .solution_lu
            .as_ref()
            .or(self.aic_lu.as_ref())
            .ok_or(PegasusError::InvalidSolution)?;
        solve_transpose_vector(factorization, rhs)
    }

    pub(crate) fn solve_linearized_transpose_batched(&self, rhs: &MatX) -> Result<MatX> {
        let factorization = self
            .solution_lu
            .as_ref()
            .or(self.aic_lu.as_ref())
            .ok_or(PegasusError::InvalidSolution)?;
        solve_transpose_matrix(factorization, rhs)
    }

    /// Velocity influence at boundary-condition control points (all bound legs included).
    #[must_use]
    pub fn control_velocity_matrix(&self) -> &MatX {
        &self.control_influence_matrix
    }

    /// Velocity influence at bound-force points, excluding each self-bound leg.
    #[must_use]
    pub fn velocity_matrix(&self) -> &MatX {
        &self.velocity_influence_matrix
    }

    #[must_use]
    pub fn solution(&self) -> &SolutionState {
        &self.solution
    }

    #[must_use]
    pub fn is_aic_valid(&self) -> bool {
        self.aic_valid
    }

    #[must_use]
    pub fn are_unit_responses_valid(&self) -> bool {
        self.unit_responses_valid
    }

    #[must_use]
    pub fn current_mach(&self) -> Scalar {
        self.current_mach
    }
}
