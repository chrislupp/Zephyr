//! Cached primal vortex-lattice solver.

use nalgebra::{Dyn, linalg::LU};

use crate::core::{EPSILON, MatX, Scalar, Vec3, VecX};
use crate::error::{PegasusError, Result};
use crate::geometry::{Aircraft, MeshInput, ReferenceQuantities};

use super::{
    FlightCondition, InfluenceMatrixBuilder, SolutionState, SymmetryConfig, VortexLattice,
};

type DynamicLu = LU<Scalar, Dyn, Dyn>;

/// Main VLM solver with cached AIC and six unit responses.
#[derive(Clone, Debug)]
pub struct VortexLatticeSolver {
    mesh: VortexLattice,
    symmetry: SymmetryConfig,
    reference_point: Vec3,
    reference_quantities: ReferenceQuantities,
    aic_matrix: MatX,
    velocity_influence_matrix: MatX,
    solution_aic_matrix: MatX,
    aic_lu: Option<DynamicLu>,
    aic_transpose_lu: Option<DynamicLu>,
    solution_transpose_lu: Option<DynamicLu>,
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
            reference_point: Vec3::zeros(),
            reference_quantities: ReferenceQuantities::default(),
            aic_matrix: MatX::zeros(0, 0),
            velocity_influence_matrix: MatX::zeros(0, 0),
            solution_aic_matrix: MatX::zeros(0, 0),
            aic_lu: None,
            aic_transpose_lu: None,
            solution_transpose_lu: None,
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
        self.solution_aic_matrix = MatX::zeros(0, 0);
        self.aic_lu = None;
        self.aic_transpose_lu = None;
        self.solution_transpose_lu = None;
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
        (self.aic_matrix, self.velocity_influence_matrix) = builder.normalwash_and_velocity();
        if !self.aic_matrix.iter().all(|value| value.is_finite()) {
            return Err(PegasusError::SingularSystem);
        }
        self.aic_lu = Some(self.aic_matrix.clone().lu());
        self.aic_transpose_lu = Some(self.aic_matrix.transpose().lu());
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
        self.unit_circulation = self
            .aic_lu
            .as_ref()
            .ok_or(PegasusError::SingularSystem)?
            .solve(&rhs)
            .ok_or(PegasusError::SingularSystem)?;
        let residual = &self.aic_matrix * &self.unit_circulation - rhs;
        if !residual.iter().all(|value| value.is_finite()) {
            return Err(PegasusError::SingularSystem);
        }
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

    fn effective_normals(&self, condition: &FlightCondition) -> Vec<Vec3> {
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

    fn effective_normals_with_derivative(
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
                self.velocity_influence_matrix[(3 * row, column)],
                self.velocity_influence_matrix[(3 * row + 1, column)],
                self.velocity_influence_matrix[(3 * row + 2, column)],
            ))
        })
    }

    fn rhs_for_normals(&self, condition: &FlightCondition, normals: &[Vec3]) -> VecX {
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
        factorization
            .solve(&rhs)
            .ok_or(PegasusError::SingularSystem)
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
        let matrix = self.matrix_for_normals(&normals);
        let n = self.mesh.vortex_count();
        let mut matrix_derivative = MatX::zeros(n, n);
        for row in 0..n {
            for column in 0..n {
                matrix_derivative[(row, column)] = normal_derivatives[row].dot(&Vec3::new(
                    self.velocity_influence_matrix[(3 * row, column)],
                    self.velocity_influence_matrix[(3 * row + 1, column)],
                    self.velocity_influence_matrix[(3 * row + 2, column)],
                ));
            }
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
        let rhs = rhs_derivative - matrix_derivative * &self.solution.circulation;
        matrix.lu().solve(&rhs).ok_or(PegasusError::SingularSystem)
    }

    pub fn solve(&mut self, condition: &FlightCondition) -> Result<SolutionState> {
        if self.mesh.vortex_count() == 0 {
            return Err(PegasusError::InvalidMesh(
                "cannot solve an empty lattice".into(),
            ));
        }
        if !condition.alpha.is_finite()
            || !condition.beta.is_finite()
            || !condition.mach.is_finite()
            || condition
                .control_deflections
                .iter()
                .any(|(_, value)| !value.is_finite())
            || condition.mach < 0.0
            || condition.dynamic_pressure < 0.0
        {
            return Err(PegasusError::InvalidFlightCondition(
                "angles, Mach, and dynamic pressure must be finite; Mach and \
                 dynamic pressure must be nonnegative"
                    .into(),
            ));
        }
        if (condition.mach - self.current_mach).abs() > EPSILON {
            self.current_mach = condition.mach;
            self.invalidate();
        }
        if !self.unit_responses_valid {
            self.solve_unit_responses()?;
        }
        let (circulation, solution_units) = if condition.control_deflections.is_empty() {
            self.solution_aic_matrix = self.aic_matrix.clone();
            self.solution_transpose_lu = self.aic_transpose_lu.clone();
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
            let factorization = matrix.clone().lu();
            let circulation = factorization
                .solve(&rhs)
                .ok_or(PegasusError::SingularSystem)?;
            let units = self.unit_responses_for_normals(&factorization, &normals)?;
            self.solution_transpose_lu = Some(matrix.transpose().lu());
            self.solution_aic_matrix = matrix;
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
        &self.solution_aic_matrix
    }

    /// AIC to use for linearization: the current control-adjusted matrix when
    /// available, otherwise the cached undeflected matrix.
    #[must_use]
    pub fn linearized_aic(&self) -> &MatX {
        if self.solution_aic_matrix.shape() == self.aic_matrix.shape()
            && self.solution_aic_matrix.nrows() > 0
        {
            &self.solution_aic_matrix
        } else {
            &self.aic_matrix
        }
    }

    pub(crate) fn linearized_normals(&self) -> &[Vec3] {
        &self.solution_normals
    }

    pub(crate) fn solve_linearized_transpose(&self, rhs: &VecX) -> Result<VecX> {
        self.solution_transpose_lu
            .as_ref()
            .or(self.aic_transpose_lu.as_ref())
            .ok_or(PegasusError::InvalidSolution)?
            .solve(rhs)
            .ok_or(PegasusError::SingularSystem)
    }

    pub(crate) fn solve_linearized_transpose_batched(&self, rhs: &MatX) -> Result<MatX> {
        self.solution_transpose_lu
            .as_ref()
            .or(self.aic_transpose_lu.as_ref())
            .ok_or(PegasusError::InvalidSolution)?
            .solve(rhs)
            .ok_or(PegasusError::SingularSystem)
    }

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
