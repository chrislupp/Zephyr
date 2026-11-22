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

//! Cached primal vortex-lattice solver.

use faer::{
    MatMut, MatRef,
    linalg::solvers::{PartialPivLu, Solve},
};

use crate::core::{EPSILON, MatX, Scalar, Vec3, VecX};
use crate::error::{Result, ZephyrError};
use crate::geometry::{Aircraft, MeshInput, ReferenceQuantities};

use super::{
    FlightCondition, InfluenceMatrixBuilder, SolutionState, SymmetryConfig, VortexLattice,
    influence::mirror_pairs,
};

type DynamicLu = PartialPivLu<Scalar>;

fn factorize(matrix: &MatX) -> DynamicLu {
    MatRef::from_column_major_slice(matrix.as_slice(), matrix.nrows(), matrix.ncols())
        .partial_piv_lu()
}

/// Block-diagonalized LU of a mirror-symmetric system.
///
/// With `o_k` the original and `m_k` the mirrored panel of pair `k`, grouping
/// each panel with its image puts the system matrix in the form
///
/// ```text
/// A = [[P, Q],      P[k, l] = A[o_k, o_l] = A[m_k, m_l]
///      [Q, P]]      Q[k, l] = A[o_k, m_l] = A[m_k, o_l]
/// ```
///
/// which the orthogonal change of variables `s = x1 + x2`, `d = x1 - x2`
/// splits into the two half-size systems `(P + Q) s = b1 + b2` and
/// `(P - Q) d = b1 - b2`. Then `x1 = (s + d) / 2` and `x2 = (s - d) / 2`.
/// `Aᵀ` has the same structure with `Pᵀ` and `Qᵀ`, so a transpose solve uses
/// the transposed half systems. Two half-size factorizations cost a quarter
/// of the full one.
///
/// The eigenvalues of `A` are exactly those of `P + Q` together with those of
/// `P - Q`, so neither half is worse conditioned than `A` itself.
#[derive(Clone, Debug)]
struct BlockLu {
    originals: Vec<usize>,
    mirrors: Vec<usize>,
    sum: DynamicLu,
    difference: DynamicLu,
}

impl BlockLu {
    /// Factorizes `matrix` by blocks if, and only if, it is exactly
    /// block-symmetric under `pairs`.
    ///
    /// The criterion is exact IEEE equality, `A[m_k, m_l] == A[o_k, o_l]` and
    /// `A[m_k, o_l] == A[o_k, m_l]` for every `k, l`; no tolerance. Under it
    /// the block matrix the half systems represent *is* `A`, entry for entry,
    /// so the only difference from the full LU is the order of rounding in
    /// the elimination, never the system being solved. Any unequal pair --
    /// including a NaN, which compares unequal to itself -- rejects the block
    /// form and the caller falls back to the full LU. The scan is `O(N²)`
    /// next to an `O(N³)` factorization, and it is fused with forming
    /// `P + Q` and `P - Q`, which needs the same reads.
    fn try_new(matrix: &MatX, pairs: &[(usize, usize)]) -> Option<Self> {
        let h = pairs.len();
        if h == 0 || matrix.nrows() != 2 * h || matrix.ncols() != 2 * h {
            return None;
        }
        let mut sum = MatX::zeros(h, h);
        let mut difference = MatX::zeros(h, h);
        for (l, &(original_l, mirror_l)) in pairs.iter().enumerate() {
            for (k, &(original_k, mirror_k)) in pairs.iter().enumerate() {
                let p = matrix[(original_k, original_l)];
                let q = matrix[(original_k, mirror_l)];
                if matrix[(mirror_k, mirror_l)] != p || matrix[(mirror_k, original_l)] != q {
                    return None;
                }
                sum[(k, l)] = p + q;
                difference[(k, l)] = p - q;
            }
        }
        let (sum, difference) = factorize_pair(&sum, &difference);
        Some(Self {
            originals: pairs.iter().map(|pair| pair.0).collect(),
            mirrors: pairs.iter().map(|pair| pair.1).collect(),
            sum,
            difference,
        })
    }

    /// Solves `A X = B` (or `Aᵀ X = B`) in place for a column-major `B` with
    /// `ncols` columns.
    fn solve_in_place(&self, data: &mut [Scalar], ncols: usize, transpose: bool) {
        let h = self.originals.len();
        let n = 2 * h;
        let mut sum = vec![0.0; h * ncols];
        let mut difference = vec![0.0; h * ncols];
        for column in 0..ncols {
            let rhs = &data[column * n..(column + 1) * n];
            for (k, (&original, &mirror)) in self.originals.iter().zip(&self.mirrors).enumerate() {
                sum[column * h + k] = rhs[original] + rhs[mirror];
                difference[column * h + k] = rhs[original] - rhs[mirror];
            }
        }
        for (lu, values) in [(&self.sum, &mut sum), (&self.difference, &mut difference)] {
            let view = MatMut::from_column_major_slice_mut(values.as_mut_slice(), h, ncols);
            if transpose {
                lu.solve_transpose_in_place(view);
            } else {
                lu.solve_in_place(view);
            }
        }
        for column in 0..ncols {
            let solution = &mut data[column * n..(column + 1) * n];
            for (k, (&original, &mirror)) in self.originals.iter().zip(&self.mirrors).enumerate() {
                let s = sum[column * h + k];
                let d = difference[column * h + k];
                solution[original] = 0.5 * (s + d);
                solution[mirror] = 0.5 * (s - d);
            }
        }
    }
}

/// Factorizes the two half systems, concurrently under `parallel`.
///
/// Each LU is computed by the same serial faer call either way and the two
/// share nothing, so the results are bit-identical between builds and
/// independent of the thread count.
#[cfg(feature = "parallel")]
fn factorize_pair(sum: &MatX, difference: &MatX) -> (DynamicLu, DynamicLu) {
    rayon::join(|| factorize(sum), || factorize(difference))
}

#[cfg(not(feature = "parallel"))]
fn factorize_pair(sum: &MatX, difference: &MatX) -> (DynamicLu, DynamicLu) {
    (factorize(sum), factorize(difference))
}

/// The one factorization type every solve in this solver goes through.
///
/// `Full` is the plain partial-pivot LU and is what every path used before
/// the block form existed; with symmetric factorization off it is the only
/// variant ever built, so default results are unchanged bit for bit.
#[derive(Clone, Debug)]
enum Factorization {
    Full(DynamicLu),
    Block(Box<BlockLu>),
}

impl Factorization {
    /// Full LU, or the block form when `pairs` is given and `matrix` is
    /// exactly block-symmetric under it.
    fn new(matrix: &MatX, pairs: Option<&[(usize, usize)]>) -> Self {
        pairs
            .and_then(|pairs| BlockLu::try_new(matrix, pairs))
            .map_or_else(
                || Self::Full(factorize(matrix)),
                |block| Self::Block(Box::new(block)),
            )
    }

    fn is_block(&self) -> bool {
        matches!(self, Self::Block(_))
    }

    fn solve_in_place(&self, data: &mut [Scalar], nrows: usize, ncols: usize, transpose: bool) {
        match self {
            Self::Full(lu) => {
                let view = MatMut::from_column_major_slice_mut(data, nrows, ncols);
                if transpose {
                    lu.solve_transpose_in_place(view);
                } else {
                    lu.solve_in_place(view);
                }
            }
            Self::Block(block) => block.solve_in_place(data, ncols, transpose),
        }
    }
}

/// A partial-pivot LU of a singular matrix still "succeeds"; the breakdown
/// only shows up as non-finite entries in what it solves for.
fn ensure_finite(solution: &[Scalar]) -> Result<()> {
    if solution.iter().all(|value| value.is_finite()) {
        Ok(())
    } else {
        Err(ZephyrError::SingularSystem)
    }
}

fn solve_matrix(factorization: &Factorization, rhs: &MatX) -> Result<MatX> {
    let mut solution = rhs.clone();
    let (nrows, ncols) = solution.shape();
    factorization.solve_in_place(solution.as_mut_slice(), nrows, ncols, false);
    ensure_finite(solution.as_slice())?;
    Ok(solution)
}

fn solve_vector(factorization: &Factorization, rhs: &VecX) -> Result<VecX> {
    let mut solution = rhs.clone();
    let nrows = solution.nrows();
    factorization.solve_in_place(solution.as_mut_slice(), nrows, 1, false);
    ensure_finite(solution.as_slice())?;
    Ok(solution)
}

fn solve_transpose_matrix(factorization: &Factorization, rhs: &MatX) -> Result<MatX> {
    let mut solution = rhs.clone();
    let (nrows, ncols) = solution.shape();
    factorization.solve_in_place(solution.as_mut_slice(), nrows, ncols, true);
    ensure_finite(solution.as_slice())?;
    Ok(solution)
}

fn solve_transpose_vector(factorization: &Factorization, rhs: &VecX) -> Result<VecX> {
    let mut solution = rhs.clone();
    let nrows = solution.nrows();
    factorization.solve_in_place(solution.as_mut_slice(), nrows, 1, true);
    ensure_finite(solution.as_slice())?;
    Ok(solution)
}

/// Main VLM solver with cached AIC and six unit responses.
#[derive(Clone, Debug)]
pub struct VortexLatticeSolver {
    mesh: VortexLattice,
    symmetry: SymmetryConfig,
    mirror_reuse: bool,
    symmetric_factorization: bool,
    reference_point: Vec3,
    reference_quantities: ReferenceQuantities,
    aic_matrix: MatX,
    velocity_influence_matrix: MatX,
    control_influence_matrix: MatX,
    solution_aic_matrix: Option<MatX>,
    aic_lu: Option<Factorization>,
    solution_lu: Option<Factorization>,
    /// Unit responses of the deflected system `solution_lu` factorizes. They
    /// depend only on that factorization, the normals and the reference
    /// point, so a repeated solve at unchanged deflections reuses them.
    solution_units: Option<MatX>,
    solution_normals: Vec<Vec3>,
    aic_valid: bool,
    unit_circulation: MatX,
    unit_responses_valid: bool,
    solution: SolutionState,
    current_mach: Scalar,
    solved_condition: Option<FlightCondition>,
}

impl Default for VortexLatticeSolver {
    fn default() -> Self {
        Self {
            mesh: VortexLattice::default(),
            symmetry: SymmetryConfig::default(),
            mirror_reuse: false,
            symmetric_factorization: false,
            reference_point: Vec3::zeros(),
            reference_quantities: ReferenceQuantities::default(),
            aic_matrix: MatX::zeros(0, 0),
            velocity_influence_matrix: MatX::zeros(0, 0),
            control_influence_matrix: MatX::zeros(0, 0),
            solution_aic_matrix: None,
            aic_lu: None,
            solution_lu: None,
            solution_units: None,
            solution_normals: Vec::new(),
            aic_valid: false,
            unit_circulation: MatX::zeros(0, 6),
            unit_responses_valid: false,
            solution: SolutionState::default(),
            current_mach: 0.0,
            solved_condition: None,
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
    /// Both assembly passes use it: the control-point pass that builds the
    /// AIC and the bound-centre pass that builds the force-velocity operator.
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

    /// Lets the solver factorize a mirror-symmetric system as two half-size
    /// systems instead of one full one, for lattices built from
    /// `y_duplicated` surfaces.
    ///
    /// Off by default, and with it off every result is exactly what it was
    /// before this option existed. When on, each system the solver factorizes
    /// -- the undeflected AIC, a control-deflected system, the operator of the
    /// selective [`crate::analysis::Analysis`] path -- is checked for the
    /// block form `[[P, Q], [Q, P]]` under the lattice's original/mirror
    /// pairing, and if it has it, `P + Q` and `P - Q` are factorized in place
    /// of the whole matrix, costing about a quarter of the full LU. Every
    /// solve that uses the factorization, including the transpose solves of
    /// the adjoints, goes through the same block form.
    ///
    /// The check is exact: every entry must equal its mirror-image entry
    /// bit for bit, so the block system is the assembled matrix itself, not
    /// an approximation of it. Results then agree with the full LU to within
    /// rounding of the elimination, not bit for bit: relative differences of a
    /// few `1e-15` in circulations, forces and gradients on the meshes in the
    /// test suite, with residuals as small as the full LU's. Serial and
    /// `parallel` builds give identical bits.
    ///
    /// The structure holds for an undeflected mirrored lattice at any flight
    /// condition, since angles and rates only enter the right-hand side, and
    /// for symmetrically deflected controls. It does not hold for
    /// antisymmetric deflections such as ailerons, for meshes that are not
    /// fully duplicated, or for half-wing models; those silently use the full
    /// LU. [`Self::symmetric_factorization_applies`] reports which happened.
    pub fn set_symmetric_factorization(&mut self, enabled: bool) {
        self.symmetric_factorization = enabled;
        self.invalidate();
    }

    #[must_use]
    pub fn symmetric_factorization(&self) -> bool {
        self.symmetric_factorization
    }

    /// Whether the factorization behind the current solution is the
    /// half-size block form of [`Self::set_symmetric_factorization`].
    ///
    /// This is the factorization the linearized solves -- control
    /// derivatives and adjoints -- use: the control-deflected system's when
    /// the last solve had deflections, otherwise the undeflected AIC's.
    /// `false` before anything has been factorized, when the setting is off,
    /// and whenever the system failed the block-structure check and fell back
    /// to the full LU.
    #[must_use]
    pub fn symmetric_factorization_applies(&self) -> bool {
        self.active_factorization()
            .is_some_and(Factorization::is_block)
    }

    /// The factorization linearized solves use: the deflected system's when
    /// one is cached, else the undeflected AIC's.
    fn active_factorization(&self) -> Option<&Factorization> {
        self.solution_lu.as_ref().or(self.aic_lu.as_ref())
    }

    /// The single place every factorization is built.
    fn factorize_system(&self, matrix: &MatX) -> Factorization {
        let pairs = if self.symmetric_factorization {
            mirror_pairs(&self.mesh)
        } else {
            None
        };
        Factorization::new(matrix, pairs.as_deref())
    }

    pub fn set_reference_point(&mut self, point: Vec3) {
        self.reference_point = point;
        self.unit_responses_valid = false;
        self.solution_units = None;
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
        self.solution_units = None;
        self.solution_normals.clear();
        self.solution.is_valid = false;
        self.solved_condition = None;
    }

    pub fn rebuild_aic(&mut self) -> Result<()> {
        if self.mesh.vortex_count() == 0 {
            self.aic_valid = false;
            return Err(ZephyrError::InvalidMesh(
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
            return Err(ZephyrError::SingularSystem);
        }
        self.aic_matrix = aic_matrix;
        self.control_influence_matrix = control_influence_matrix;
        self.velocity_influence_matrix = velocity_influence_matrix;
        self.aic_lu = Some(self.factorize_system(&self.aic_matrix));
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
            self.aic_lu.as_ref().ok_or(ZephyrError::SingularSystem)?,
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
                            normal = super::control::rotate(
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
                        normal = super::control::rotate(normal, &axis, angle);
                        derivative = super::control::rotate(derivative, &axis, angle);
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
        factorization: &Factorization,
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
        let rhs = self.control_residual_derivative(condition, control_name)?;
        let lu = self
            .active_factorization()
            .ok_or(ZephyrError::InvalidSolution)?;
        solve_vector(lu, &rhs)
    }

    /// Check that derivative inputs describe the solved state.
    pub(crate) fn check_solved_condition(&self, condition: &FlightCondition) -> Result<()> {
        if !self.solution.is_valid || self.solved_condition.as_ref() != Some(condition) {
            return Err(ZephyrError::InvalidSolution);
        }
        Ok(())
    }

    /// `db/d(delta) - dA/d(delta) Gamma`, for adjoint contractions.
    pub(crate) fn control_residual_derivative(
        &self,
        condition: &FlightCondition,
        control_name: &str,
    ) -> Result<VecX> {
        self.check_solved_condition(condition)?;
        if !self
            .mesh
            .control_surface_names()
            .iter()
            .any(|name| name == control_name)
        {
            return Err(ZephyrError::InvalidFlightCondition(format!(
                "unknown control '{control_name}'"
            )));
        }
        let (_, derivatives) = self.effective_normals_with_derivative(condition, control_name);
        let induced = &self.control_influence_matrix * &self.solution.circulation;
        let freestream = Self::freestream_direction(condition);
        let omega = Vec3::new(
            condition.roll_rate,
            condition.pitch_rate,
            condition.yaw_rate,
        );
        Ok(VecX::from_fn(self.mesh.vortex_count(), |row, _| {
            let vortex = &self.mesh.vortices()[row];
            let mut velocity =
                Vec3::new(induced[3 * row], induced[3 * row + 1], induced[3 * row + 2]);
            if vortex.affected_by_freestream {
                velocity +=
                    freestream + omega.cross(&(vortex.control_point - self.reference_point));
            }
            -derivatives[row].dot(&velocity)
        }))
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
                        self.mirror_reuse,
                        work,
                    )
                };
            if !matrix.iter().all(|v| v.is_finite()) {
                self.invalidate();
                return Err(ZephyrError::SingularSystem);
            }
            self.aic_matrix = matrix;
            self.control_influence_matrix = control;
            self.velocity_influence_matrix = velocities;
            self.aic_lu = Some(self.factorize_system(&self.aic_matrix));
            self.solution_lu = None;
            self.solution_units = None;
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
                self.mirror_reuse,
                work,
            )
            .2;
        }
        if matrix_changed || !reuse_solution || !self.solution.is_valid {
            let rhs = self.rhs_for_normals(condition, &self.solution_normals);
            let gamma = solve_vector(
                self.aic_lu.as_ref().ok_or(ZephyrError::InvalidSolution)?,
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
        // Derivative entry points check the condition a solution belongs to;
        // a reused solution is reused only for an identical condition.
        self.solved_condition = Some(condition.clone());
        Ok(())
    }

    pub fn solve(&mut self, condition: &FlightCondition) -> Result<SolutionState> {
        if self.mesh.vortex_count() == 0 {
            return Err(ZephyrError::InvalidMesh(
                "cannot solve an empty lattice".into(),
            ));
        }
        for (name, _) in &condition.control_deflections {
            if !self.mesh.control_surface_names().contains(name) {
                return Err(ZephyrError::InvalidFlightCondition(format!(
                    "unknown control '{name}'"
                )));
            }
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
            return Err(ZephyrError::InvalidFlightCondition(
                "angles, body rates, Mach, and dynamic pressure must be finite; \
                 Mach and dynamic pressure must be nonnegative"
                    .into(),
            ));
        }
        // Prandtl-Glauert has no meaning at or above Mach 1. Returning the
        // incompressible answer silently, as clamping beta to 1 did, hands back
        // plausible numbers for a condition the method cannot model.
        if condition.mach >= 1.0 {
            return Err(ZephyrError::InvalidFlightCondition(format!(
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
            self.solution_units = None;
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
            // Geometry/Mach changes invalidate this cache. Angle and rate
            // updates at fixed controls can reuse the deflected factorization.
            if self.solution_lu.is_none() || self.solution_normals != normals {
                let matrix = self.matrix_for_normals(&normals);
                self.solution_lu = Some(self.factorize_system(&matrix));
                self.solution_units = None;
                self.solution_aic_matrix = Some(matrix);
            }
            let rhs = self.rhs_for_normals(condition, &normals);
            let factorization = self.solution_lu.as_ref().expect("deflected factorization");
            let circulation = solve_vector(factorization, &rhs)?;
            // The unit responses depend on the factorization, the normals and
            // the reference point only, none of which changed if they are
            // cached, so reusing them returns the same bits as recomputing.
            let units = match self.solution_units.take() {
                Some(units) => units,
                None => self.unit_responses_for_normals(factorization, &normals)?,
            };
            self.solution_units = Some(units.clone());
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
        self.solved_condition = Some(condition.clone());
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
            .active_factorization()
            .ok_or(ZephyrError::InvalidSolution)?;
        solve_transpose_vector(factorization, rhs)
    }

    pub(crate) fn solve_linearized_transpose_batched(&self, rhs: &MatX) -> Result<MatX> {
        let factorization = self
            .active_factorization()
            .ok_or(ZephyrError::InvalidSolution)?;
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
