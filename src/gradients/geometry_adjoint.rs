//! Reverse-mode geometry sensitivities of the lattice.
//!
//! Every geometry derivative here is adjoint-based. The context solves the flow
//! adjoint once per output, then sweeps the lattice **once** in reverse to build
//! a cotangent -- the gradient of each output with respect to every lattice
//! quantity (bound-vortex endpoints, control points, normals, strip chords).
//! A directional derivative is then an inner product, and a mesh-node gradient
//! is a scatter of the cotangent onto the nodes that formed each panel.
//!
//! The distinction matters: the reverse sweep costs one pass over the influence
//! matrix regardless of how many design variables there are, whereas seeding
//! forward once per variable multiplies that same pass by their count.

use crate::core::{EPSILON, MatX, PI, Scalar, SymmetryType, Vec3};
use crate::error::{PegasusError, Result};
use crate::geometry::{MeshInput, ReferenceQuantities};
use crate::solver::{FlightCondition, SolutionState, VortexLattice, VortexLatticeSolver};

const INV_2_PI: Scalar = 1.0 / (2.0 * PI);

use super::forces::body_output_cotangent;
use super::{
    AdjointSolver, AeroOutput, OutputSelection, force_circulation_gradient,
    force_induced_velocity_gradient, horseshoe_velocity_jacobians, trefftz_circulation_gradients,
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

/// Gradient of one scalar output with respect to every lattice quantity.
///
/// The transpose of [`LatticeDirections`]: where a direction carries a
/// perturbation into the lattice, a cotangent carries a sensitivity back out.
#[derive(Clone, Debug)]
pub(crate) struct LatticeCotangents {
    pub bound_left: Vec<Vec3>,
    pub bound_right: Vec<Vec3>,
    pub control_point: Vec<Vec3>,
    pub normal: Vec<Vec3>,
    pub local_chord: Vec<Scalar>,
    pub strip_chord: Vec<Scalar>,
}

impl LatticeCotangents {
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

    /// Directional derivative along `direction`.
    pub(crate) fn dot(&self, direction: &LatticeDirections) -> Scalar {
        let vectors = self
            .bound_left
            .iter()
            .zip(&direction.bound_left)
            .chain(self.bound_right.iter().zip(&direction.bound_right))
            .chain(self.control_point.iter().zip(&direction.control_point))
            .chain(self.normal.iter().zip(&direction.normal))
            .map(|(cotangent, tangent)| cotangent.dot(tangent))
            .sum::<Scalar>();
        let scalars = self
            .local_chord
            .iter()
            .zip(&direction.local_chord)
            .chain(self.strip_chord.iter().zip(&direction.strip_chord))
            .map(|(cotangent, tangent)| cotangent * tangent)
            .sum::<Scalar>();
        vectors + scalars
    }
}

pub(crate) struct GeometryAdjointContext<'a> {
    solver: &'a VortexLatticeSolver,
    outputs: Vec<AeroOutput>,
    cotangents: Vec<LatticeCotangents>,
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
        let cotangents = geometry_cotangents(solver, condition, &outputs, &adjoints);
        Ok(Self {
            solver,
            outputs,
            cotangents,
        })
    }

    /// Lattice the context was built against.
    pub(crate) fn mesh(&self) -> &VortexLattice {
        self.solver.mesh()
    }

    /// Cotangent of each selected output, paired with the output it belongs to.
    pub(crate) fn cotangents(&self) -> impl Iterator<Item = (AeroOutput, &LatticeCotangents)> {
        self.outputs.iter().copied().zip(&self.cotangents)
    }

    /// Directional derivative of every selected output along `direction`.
    ///
    /// An inner product against the cotangents built when the context was
    /// created, so evaluating many directions costs one pass over the lattice
    /// each rather than one pass over the influence matrix each.
    pub(crate) fn evaluate(
        &self,
        direction: &LatticeDirections,
    ) -> Result<Vec<(AeroOutput, Scalar)>> {
        let mesh = self.solver.mesh();
        if direction.bound_left.len() != mesh.vortex_count()
            || direction.bound_right.len() != mesh.vortex_count()
            || direction.control_point.len() != mesh.vortex_count()
            || direction.normal.len() != mesh.vortex_count()
            || direction.local_chord.len() != mesh.vortex_count()
            || direction.strip_chord.len() != mesh.strip_count()
        {
            return Err(PegasusError::InvalidGeometry(
                "lattice direction dimensions do not match the solved mesh".into(),
            ));
        }
        Ok(self
            .cotangents()
            .map(|(output, cotangent)| (output, cotangent.dot(direction)))
            .collect())
    }
}

/// Reverse sweep: gradient of every selected output with respect to the lattice.
///
/// One pass over the influence matrix serves every output and every downstream
/// design variable. The implicit part rides the adjoints solved by the caller;
/// the explicit near-field and Trefftz terms are accumulated alongside it.
fn geometry_cotangents(
    solver: &VortexLatticeSolver,
    condition: &FlightCondition,
    outputs: &[AeroOutput],
    adjoints: &MatX,
) -> Vec<LatticeCotangents> {
    let mesh = solver.mesh();
    let solution = solver.solution();
    let reference = solver.reference_quantities();
    let mut result = vec![LatticeCotangents::zeros(mesh); outputs.len()];
    if outputs.is_empty() {
        return result;
    }

    // Seeds for the near-field force path: how each output responds to the
    // induced velocity at a control point. Identically zero for far-field drag,
    // which takes its geometry sensitivity from the Trefftz sweep instead.
    let velocity_seeds: Vec<Vec<Vec3>> = outputs
        .iter()
        .map(|output| {
            force_induced_velocity_gradient(mesh, solution, condition, reference, *output)
        })
        .collect();

    let velocity_influence = solver.velocity_matrix();
    let pg = prandtl_glauert_factor(solver.current_mach());
    let freestream = VortexLatticeSolver::freestream_direction(condition);
    let omega = Vec3::new(
        condition.roll_rate,
        condition.pitch_rate,
        condition.yaw_rate,
    );
    let sweep = SweepInputs {
        mesh,
        solution,
        adjoints,
        velocity_influence,
        velocity_seeds: &velocity_seeds,
        reference_point: *solver.reference_point(),
        freestream,
        omega,
        pg,
        outputs: outputs.len(),
    };

    // Rows are swept in fixed-size chunks whose partial column sums are added
    // back in chunk order. Both builds run the same chunking, so the parallel
    // result is bit-identical to the serial one and independent of the thread
    // count -- the property the assembly driver in `solver::influence` holds to.
    let count = mesh.vortex_count();
    let rows_per_chunk = count.div_ceil(SWEEP_CHUNKS).max(1);
    let ranges: Vec<(usize, usize)> = (0..count)
        .step_by(rows_per_chunk)
        .map(|first| (first, (first + rows_per_chunk).min(count)))
        .collect();
    for chunk in drive_row_chunks(&ranges, &sweep) {
        chunk.accumulate_into(&mut result);
    }

    for (index, output) in outputs.iter().enumerate() {
        if *output == AeroOutput::InducedDrag {
            trefftz_geometry_cotangents(
                mesh,
                solution,
                reference.area,
                0.01,
                0.1,
                &mut result[index],
            );
        } else {
            force_geometry_cotangents(
                mesh,
                solution,
                condition,
                reference,
                *output,
                &mut result[index],
            );
        }
    }
    result
}

/// Row chunks the geometry sweep is split into.
///
/// Fixed rather than derived from the thread count: the chunk partials are
/// summed in chunk order, so the count has to be a property of the problem for
/// the result to be reproducible.
const SWEEP_CHUNKS: usize = 32;

/// Panel count below which the sweep stays on one thread.
///
/// Matches the assembly driver's threshold: the sweep is `O(N^2)` while handing
/// out chunks is `O(N)`, so small lattices lose more to distribution than they
/// gain. Affects speed only -- both branches run the same per-chunk code.
#[cfg(feature = "parallel")]
const PARALLEL_MIN_PANELS: usize = 128;

/// Everything the row sweep reads. Shared immutably across chunks.
struct SweepInputs<'a> {
    mesh: &'a VortexLattice,
    solution: &'a SolutionState,
    adjoints: &'a MatX,
    velocity_influence: &'a MatX,
    velocity_seeds: &'a [Vec<Vec3>],
    reference_point: Vec3,
    freestream: Vec3,
    omega: Vec3,
    pg: Scalar,
    outputs: usize,
}

/// One chunk's contribution: exclusive rows, and partial sums over columns.
struct SweepChunk {
    first_row: usize,
    control_point: Vec<Vec<Vec3>>,
    normal: Vec<Vec<Vec3>>,
    bound_left: Vec<Vec<Vec3>>,
    bound_right: Vec<Vec<Vec3>>,
}

impl SweepChunk {
    fn accumulate_into(self, result: &mut [LatticeCotangents]) {
        for (index, cotangent) in result.iter_mut().enumerate() {
            for (offset, value) in self.control_point[index].iter().enumerate() {
                cotangent.control_point[self.first_row + offset] += value;
            }
            for (offset, value) in self.normal[index].iter().enumerate() {
                cotangent.normal[self.first_row + offset] += value;
            }
            for (column, value) in self.bound_left[index].iter().enumerate() {
                cotangent.bound_left[column] += value;
            }
            for (column, value) in self.bound_right[index].iter().enumerate() {
                cotangent.bound_right[column] += value;
            }
        }
    }
}

/// Sweep one range of rows. Pure in its inputs, which is what lets the chunks
/// run in any order.
fn sweep_row_range(range: (usize, usize), sweep: &SweepInputs<'_>) -> SweepChunk {
    let (first_row, last_row) = range;
    let count = sweep.mesh.vortex_count();
    let rows = last_row - first_row;
    let mut chunk = SweepChunk {
        first_row,
        control_point: vec![vec![Vec3::zeros(); rows]; sweep.outputs],
        normal: vec![vec![Vec3::zeros(); rows]; sweep.outputs],
        bound_left: vec![vec![Vec3::zeros(); count]; sweep.outputs],
        bound_right: vec![vec![Vec3::zeros(); count]; sweep.outputs],
    };

    for row in first_row..last_row {
        let local = row - first_row;
        let target = &sweep.mesh.vortices()[row];

        // Right-hand side: the boundary condition moves with the control point
        // and the normal.
        if target.affected_by_freestream {
            let local_velocity = sweep.freestream
                + sweep
                    .omega
                    .cross(&(target.control_point - sweep.reference_point));
            let control_point_seed = -target.normal_at_control_point.cross(&sweep.omega);
            for index in 0..sweep.outputs {
                let psi = sweep.adjoints[(row, index)];
                chunk.control_point[index][local] += psi * control_point_seed;
                chunk.normal[index][local] -= psi * local_velocity;
            }
        }

        let normal = target.normal_at_control_point;
        for (column, source) in sweep.mesh.vortices().iter().enumerate() {
            let core = core_radius(sweep.mesh, column, target.component_index);
            let include_bound = row != column;
            let gamma = sweep.solution.circulation[column];

            // One kernel evaluation serves both geometry terms. The influence
            // coefficient is `normal . velocity`, so its gradient is these same
            // Jacobians contracted with the target normal -- evaluating
            // `aic_element_gradient` as well would rebuild the horseshoe
            // geometry three more times for an answer already in hand.
            let jacobians = horseshoe_velocity_jacobians(
                &target.control_point,
                &source.bound_vortex_left,
                &source.bound_vortex_right,
                sweep.pg,
                core,
                include_bound,
            );

            // Influence-matrix term: psi^T (dA/dX) Gamma. The normal derivative
            // of `normal . velocity` is the velocity itself, which the solver
            // already assembled for the primal.
            if gamma != 0.0 {
                let by_bound_left = jacobians.bound_left.tr_mul(&normal);
                let by_bound_right = jacobians.bound_right.tr_mul(&normal);
                let by_control_point = jacobians.field_point.tr_mul(&normal);
                let by_normal = Vec3::new(
                    sweep.velocity_influence[(3 * row, column)],
                    sweep.velocity_influence[(3 * row + 1, column)],
                    sweep.velocity_influence[(3 * row + 2, column)],
                );
                for index in 0..sweep.outputs {
                    let weight = -sweep.adjoints[(row, index)] * gamma;
                    if weight == 0.0 {
                        continue;
                    }
                    chunk.bound_left[index][column] += weight * by_bound_left;
                    chunk.bound_right[index][column] += weight * by_bound_right;
                    chunk.control_point[index][local] += weight * by_control_point;
                    chunk.normal[index][local] += weight * by_normal;
                }
            }

            // Near-field term: the induced velocity that sets each panel force
            // moves with the same geometry. Contracting the shared Jacobians
            // with each output's seed keeps the kernel count independent of the
            // number of outputs.
            if gamma.abs() <= EPSILON {
                continue;
            }
            for index in 0..sweep.outputs {
                let seed = sweep.velocity_seeds[index][row];
                if seed.norm_squared() <= 0.0 {
                    continue;
                }
                let scaled = gamma * seed;
                chunk.control_point[index][local] += jacobians.field_point.tr_mul(&scaled);
                chunk.bound_left[index][column] += jacobians.bound_left.tr_mul(&scaled);
                chunk.bound_right[index][column] += jacobians.bound_right.tr_mul(&scaled);
            }
        }
    }
    chunk
}

/// Sweep every row chunk, serially.
#[cfg(not(feature = "parallel"))]
fn drive_row_chunks(ranges: &[(usize, usize)], sweep: &SweepInputs<'_>) -> Vec<SweepChunk> {
    ranges
        .iter()
        .map(|range| sweep_row_range(*range, sweep))
        .collect()
}

/// Sweep every row chunk, one rayon task per chunk.
///
/// Chunks share nothing mutable and are collected in index order, so the
/// partial sums are added back in the same order the serial driver uses.
#[cfg(feature = "parallel")]
fn drive_row_chunks(ranges: &[(usize, usize)], sweep: &SweepInputs<'_>) -> Vec<SweepChunk> {
    use rayon::prelude::*;

    if sweep.mesh.vortex_count() < PARALLEL_MIN_PANELS {
        return ranges
            .iter()
            .map(|range| sweep_row_range(*range, sweep))
            .collect();
    }
    ranges
        .par_iter()
        .map(|range| sweep_row_range(*range, sweep))
        .collect()
}

/// Reverse of the fixed-circulation near-field force sweep.
///
/// Panel force depends on geometry through the bound-vortex line, and the
/// moment additionally through the control point that carries the moment arm.
fn force_geometry_cotangents(
    mesh: &VortexLattice,
    solution: &SolutionState,
    condition: &FlightCondition,
    reference: &ReferenceQuantities,
    output: AeroOutput,
    cotangent: &mut LatticeCotangents,
) {
    if !solution.is_valid || reference.area <= 0.0 || reference.chord <= 0.0 {
        return;
    }
    let (force_seed, moment_seed) = body_output_cotangent(output, condition.alpha);
    let freestream = VortexLatticeSolver::freestream_direction(condition);
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
        let velocity = freestream + induced;
        let line = vortex.bound_vortex_right - vortex.bound_vortex_left;
        let arm = vortex.control_point - reference.moment_reference_point;
        let panel_force = 2.0 * solution.circulation[index] * velocity.cross(&line);

        // Moving the control point swings the moment arm.
        cotangent.control_point[index] +=
            panel_force.cross(&moment_seed) / (reference.area * reference.chord);

        // Moving either endpoint changes the bound-vortex line. `combined` is
        // the same force/moment seed the induced-velocity path uses.
        let combined = force_seed / reference.area
            + moment_seed.cross(&arm) / (reference.area * reference.chord);
        let line_cotangent = 2.0 * solution.circulation[index] * combined.cross(&velocity);
        cotangent.bound_right[index] += line_cotangent;
        cotangent.bound_left[index] -= line_cotangent;
    }
}

/// Reverse of the fixed-circulation Trefftz-plane induced-drag sweep.
///
/// Mirrors `trefftz_geometry_jvp` term for term; the strip trailing-edge
/// coordinates and the strip chord (through the core radius) are the only
/// geometry it sees.
fn trefftz_geometry_cotangents(
    mesh: &VortexLattice,
    solution: &SolutionState,
    reference_area: Scalar,
    core_chord_fraction: Scalar,
    core_width_fraction: Scalar,
    cotangent: &mut LatticeCotangents,
) {
    let strip_count = mesh.strip_count();
    if !solution.is_valid || strip_count == 0 || reference_area <= 0.0 {
        return;
    }
    let mut gamma = vec![0.0; strip_count];
    let mut y_left = vec![0.0; strip_count];
    let mut z_left = vec![0.0; strip_count];
    let mut y_right = vec![0.0; strip_count];
    let mut z_right = vec![0.0; strip_count];
    for (index, strip) in mesh.strips().iter().enumerate() {
        gamma[index] = solution
            .circulation
            .rows(strip.first_vortex_index, strip.num_vortices)
            .sum();
        let vortex = &mesh.vortices()[strip.first_vortex_index + strip.num_vortices - 1];
        y_left[index] = vortex.bound_vortex_left.y;
        z_left[index] = vortex.bound_vortex_left.z;
        y_right[index] = vortex.bound_vortex_right.y;
        z_right[index] = vortex.bound_vortex_right.z;
    }
    let y_center = |i: usize| 0.5 * (y_left[i] + y_right[i]);
    let z_center = |i: usize| 0.5 * (z_left[i] + z_right[i]);

    let mut c_y_left = vec![0.0; strip_count];
    let mut c_z_left = vec![0.0; strip_count];
    let mut c_y_right = vec![0.0; strip_count];
    let mut c_z_right = vec![0.0; strip_count];
    let mut c_chord = vec![0.0; strip_count];

    for target in 0..strip_count {
        let dy = y_right[target] - y_left[target];
        let dz = z_right[target] - z_left[target];
        if dy.hypot(dz) < EPSILON {
            continue;
        }

        // Primal induced velocity at this strip, needed by the width terms.
        let (mut velocity_y, mut velocity_z) = (0.0, 0.0);
        for (source, source_gamma) in gamma.iter().enumerate() {
            let Some(kernel) = trefftz_kernel(
                mesh,
                target,
                source,
                &y_left,
                &z_left,
                &y_right,
                &z_right,
                y_center(target),
                z_center(target),
                core_chord_fraction,
                core_width_fraction,
            ) else {
                continue;
            };
            velocity_y += INV_2_PI
                * source_gamma
                * (kernel.dz1 / kernel.radius1 - kernel.dz2 / kernel.radius2);
            velocity_z += INV_2_PI
                * source_gamma
                * (-kernel.dy1 / kernel.radius1 + kernel.dy2 / kernel.radius2);
        }

        // Explicit width terms: 2 Gamma dy / S and -2 Gamma dz / S carry the
        // strip width straight into the far-field coefficients.
        let target_gamma = gamma[target];
        let by_dz = target_gamma * velocity_y / reference_area;
        let by_dy = -target_gamma * velocity_z / reference_area;
        c_z_right[target] += by_dz;
        c_z_left[target] -= by_dz;
        c_y_right[target] += by_dy;
        c_y_left[target] -= by_dy;

        // Reverse of the induced-velocity accumulation.
        let weight_y = target_gamma * dz / reference_area;
        let weight_z = -target_gamma * dy / reference_area;
        for source in 0..strip_count {
            let Some(kernel) = trefftz_kernel(
                mesh,
                target,
                source,
                &y_left,
                &z_left,
                &y_right,
                &z_right,
                y_center(target),
                z_center(target),
                core_chord_fraction,
                core_width_fraction,
            ) else {
                continue;
            };
            let scale = INV_2_PI * gamma[source];
            let (radius1, radius2) = (kernel.radius1, kernel.radius2);
            let mut c_dz1 = scale * weight_y / radius1;
            let mut c_dy1 = -scale * weight_z / radius1;
            let mut c_dz2 = -scale * weight_y / radius2;
            let mut c_dy2 = scale * weight_z / radius2;
            let c_radius1 =
                scale * (weight_z * kernel.dy1 - weight_y * kernel.dz1) / (radius1 * radius1);
            let c_radius2 =
                scale * (weight_y * kernel.dz2 - weight_z * kernel.dy2) / (radius2 * radius2);

            // radius = sqrt(q^2 + core^4), so dr = (q dq + 0.5 d(core^4)) / r
            // and dq = 2 (dy d_dy + dz d_dz).
            let c_q1 = c_radius1 * kernel.q1 / radius1;
            let c_q2 = c_radius2 * kernel.q2 / radius2;
            c_dy1 += 2.0 * kernel.dy1 * c_q1;
            c_dz1 += 2.0 * kernel.dz1 * c_q1;
            c_dy2 += 2.0 * kernel.dy2 * c_q2;
            c_dz2 += 2.0 * kernel.dz2 * c_q2;

            let c_core =
                4.0 * kernel.core.powi(3) * 0.5 * (c_radius1 / radius1 + c_radius2 / radius2);
            match kernel.core_source {
                CoreSource::None => {}
                CoreSource::Chord => c_chord[source] += core_chord_fraction * c_core,
                CoreSource::Width => {
                    let source_dy = y_right[source] - y_left[source];
                    let source_dz = z_right[source] - z_left[source];
                    let width = source_dy.hypot(source_dz);
                    if width > 0.0 {
                        let by_width = core_width_fraction * c_core / width;
                        c_y_right[source] += by_width * source_dy;
                        c_y_left[source] -= by_width * source_dy;
                        c_z_right[source] += by_width * source_dz;
                        c_z_left[source] -= by_width * source_dz;
                    }
                }
            }

            // The offsets are taken from the target centre to each source
            // endpoint, so the centre collects half of both.
            let center_y = 0.5 * (c_dy1 + c_dy2);
            let center_z = 0.5 * (c_dz1 + c_dz2);
            c_y_left[target] += center_y;
            c_y_right[target] += center_y;
            c_z_left[target] += center_z;
            c_z_right[target] += center_z;
            c_y_left[source] -= c_dy1;
            c_y_right[source] -= c_dy2;
            c_z_left[source] -= c_dz1;
            c_z_right[source] -= c_dz2;
        }
    }

    for (index, strip) in mesh.strips().iter().enumerate() {
        let vortex_index = strip.first_vortex_index + strip.num_vortices - 1;
        cotangent.bound_left[vortex_index].y += c_y_left[index];
        cotangent.bound_left[vortex_index].z += c_z_left[index];
        cotangent.bound_right[vortex_index].y += c_y_right[index];
        cotangent.bound_right[vortex_index].z += c_z_right[index];
        cotangent.strip_chord[index] += c_chord[index];
    }
}

/// Which quantity set the source strip's core radius, so the reverse sweep
/// sends the sensitivity to the same place the forward one took it from.
#[derive(Clone, Copy, PartialEq, Eq)]
enum CoreSource {
    None,
    Chord,
    Width,
}

/// Primal Trefftz-plane geometry shared by the velocity and reverse loops.
struct TrefftzKernel {
    dy1: Scalar,
    dy2: Scalar,
    dz1: Scalar,
    dz2: Scalar,
    q1: Scalar,
    q2: Scalar,
    radius1: Scalar,
    radius2: Scalar,
    core: Scalar,
    core_source: CoreSource,
}

#[allow(clippy::too_many_arguments)]
fn trefftz_kernel(
    mesh: &VortexLattice,
    target: usize,
    source: usize,
    y_left: &[Scalar],
    z_left: &[Scalar],
    y_right: &[Scalar],
    z_right: &[Scalar],
    y_center: Scalar,
    z_center: Scalar,
    core_chord_fraction: Scalar,
    core_width_fraction: Scalar,
) -> Option<TrefftzKernel> {
    let target_strip = &mesh.strips()[target];
    let source_strip = &mesh.strips()[source];
    let source_width = (y_right[source] - y_left[source]).hypot(z_right[source] - z_left[source]);
    let (core, core_source) = if target_strip.surface_index == source_strip.surface_index {
        (0.0, CoreSource::None)
    } else {
        let chord = core_chord_fraction * source_strip.chord_center;
        let width = core_width_fraction * source_width;
        if chord >= width {
            (chord, CoreSource::Chord)
        } else {
            (width, CoreSource::Width)
        }
    };
    let core_4 = core.powi(4);
    let dy1 = y_center - y_left[source];
    let dy2 = y_center - y_right[source];
    let dz1 = z_center - z_left[source];
    let dz2 = z_center - z_right[source];
    let q1 = dy1 * dy1 + dz1 * dz1;
    let q2 = dy2 * dy2 + dz2 * dz2;
    let radius1 = (q1 * q1 + core_4).sqrt();
    let radius2 = (q2 * q2 + core_4).sqrt();
    (radius1 > EPSILON && radius2 > EPSILON).then_some(TrefftzKernel {
        dy1,
        dy2,
        dz1,
        dz2,
        q1,
        q2,
        radius1,
        radius2,
        core,
        core_source,
    })
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

pub(crate) fn normalized_direction(vector: &Vec3, direction: &Vec3) -> Vec3 {
    let norm = vector.norm();
    if norm <= EPSILON {
        Vec3::zeros()
    } else {
        let unit = *vector / norm;
        (direction - unit * unit.dot(direction)) / norm
    }
}

pub(crate) fn reflected_direction(direction: Vec3) -> Vec3 {
    Vec3::new(direction.x, -direction.y, direction.z)
}

/// Scatter a lattice cotangent onto the mesh nodes that formed the panels.
///
/// The transpose of [`direct_mesh_node_direction`]: instead of seeding one node
/// coordinate and pushing it through the lattice, every panel hands its
/// sensitivity back to its four corners in a single pass. That is what makes a
/// node gradient cost the same whether the mesh has ten nodes or ten thousand.
pub(crate) fn scatter_mesh_node_cotangents(
    input: &MeshInput,
    lattice: &VortexLattice,
    surface_index: usize,
    cotangent: &LatticeCotangents,
) -> Result<Vec<Vec3>> {
    let surface = input
        .surfaces
        .get(surface_index)
        .ok_or_else(|| PegasusError::InvalidMesh("surface index is out of range".into()))?;
    let surface_slot: usize = input.surfaces[..surface_index]
        .iter()
        .map(|surface| 1 + usize::from(surface.y_duplicated))
        .sum();
    let info = lattice
        .surface_info()
        .get(surface_slot)
        .ok_or_else(|| PegasusError::InvalidMesh("surface lattice is missing".into()))?;

    // Fold the mirrored half back onto the original. The reflection is its own
    // transpose, and the duplicate's left endpoint came from the original's
    // right, so the two swap on the way back.
    let mut bound_left = cotangent.bound_left.clone();
    let mut bound_right = cotangent.bound_right.clone();
    let mut control_point = cotangent.control_point.clone();
    let mut normal = cotangent.normal.clone();
    let mut local_chord = cotangent.local_chord.clone();
    let mut strip_chord = cotangent.strip_chord.clone();
    if surface.y_duplicated {
        let duplicate = lattice
            .surface_info()
            .get(surface_slot + 1)
            .ok_or_else(|| {
                PegasusError::InvalidMesh("mirrored surface lattice is missing".into())
            })?;
        for strip_offset in 0..info.num_strips {
            let original_strip = info.first_strip_index + strip_offset;
            let duplicate_strip = duplicate.first_strip_index + strip_offset;
            strip_chord[original_strip] += strip_chord[duplicate_strip];
            for chord_offset in 0..info.num_chordwise {
                let original =
                    info.first_vortex_index + strip_offset * info.num_chordwise + chord_offset;
                let mirrored =
                    duplicate.first_vortex_index + strip_offset * info.num_chordwise + chord_offset;
                bound_right[original] += reflected_direction(bound_left[mirrored]);
                bound_left[original] += reflected_direction(bound_right[mirrored]);
                let mirrored_control = control_point[mirrored];
                let mirrored_normal = normal[mirrored];
                control_point[original] += reflected_direction(mirrored_control);
                normal[original] += reflected_direction(mirrored_normal);
                local_chord[original] += local_chord[mirrored];
            }
        }
    }

    let node = |i: usize, j: usize| *surface.node(i, j).expect("validated mesh index");
    let index_of = |i: usize, j: usize| i * (surface.ny + 1) + j;
    let mut nodes = vec![Vec3::zeros(); (surface.nx + 1) * (surface.ny + 1)];

    for strip_offset in 0..surface.ny {
        let strip_index = info.first_strip_index + strip_offset;
        let mut chord_cotangent = strip_chord[strip_index];
        for chord_offset in 0..surface.nx {
            chord_cotangent +=
                local_chord[info.first_vortex_index + strip_offset * surface.nx + chord_offset];
        }
        // The strip chord is the mean of the two edge chord lengths, so its
        // sensitivity runs along each edge's unit vector.
        for j in [strip_offset, strip_offset + 1] {
            let chord = node(surface.nx, j) - node(0, j);
            let norm = chord.norm();
            if norm <= EPSILON {
                continue;
            }
            let unit = chord / norm;
            nodes[index_of(surface.nx, j)] += 0.5 * chord_cotangent * unit;
            nodes[index_of(0, j)] -= 0.5 * chord_cotangent * unit;
        }

        for chord_offset in 0..surface.nx {
            let vortex_index = info.first_vortex_index + strip_offset * surface.nx + chord_offset;
            let (i00, i10) = (
                index_of(chord_offset, strip_offset),
                index_of(chord_offset + 1, strip_offset),
            );
            let (i01, i11) = (
                index_of(chord_offset, strip_offset + 1),
                index_of(chord_offset + 1, strip_offset + 1),
            );

            // Quarter-chord bound vortex and three-quarter-chord control point.
            let left = bound_left[vortex_index];
            let right = bound_right[vortex_index];
            let control = control_point[vortex_index];
            nodes[i00] += 0.75 * left + 0.125 * control;
            nodes[i10] += 0.25 * left + 0.375 * control;
            nodes[i01] += 0.75 * right + 0.125 * control;
            nodes[i11] += 0.25 * right + 0.375 * control;

            // Normal: the panel diagonals cross to an unnormalized normal, so
            // the cotangent passes through the same projection before it is
            // split between the corners.
            let diagonal_a =
                node(chord_offset + 1, strip_offset + 1) - node(chord_offset, strip_offset);
            let diagonal_b =
                node(chord_offset + 1, strip_offset) - node(chord_offset, strip_offset + 1);
            let cross = diagonal_a.cross(&diagonal_b);
            let projected = normalized_direction(&cross, &normal[vortex_index]);
            let by_diagonal_a = diagonal_b.cross(&projected);
            let by_diagonal_b = projected.cross(&diagonal_a);
            nodes[i11] += by_diagonal_a;
            nodes[i00] -= by_diagonal_a;
            nodes[i10] += by_diagonal_b;
            nodes[i01] -= by_diagonal_b;
        }
    }
    Ok(nodes)
}
