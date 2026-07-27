//! Aerodynamic influence-coefficient kernels.

use crate::core::{EPSILON, MatX, PI, Scalar, SymmetryType, Vec3, VecX};

use super::VortexLattice;

const INV_4_PI: Scalar = 1.0 / (4.0 * PI);

/// Symmetry planes used while assembling influence matrices.
#[derive(Clone, Copy, Debug)]
pub struct SymmetryConfig {
    pub y_symmetry: SymmetryType,
    pub y_symmetry_plane: Scalar,
    pub z_symmetry: SymmetryType,
    pub z_symmetry_plane: Scalar,
}

impl Default for SymmetryConfig {
    fn default() -> Self {
        Self {
            y_symmetry: SymmetryType::None,
            y_symmetry_plane: 0.0,
            z_symmetry: SymmetryType::None,
            z_symmetry_plane: 0.0,
        }
    }
}

/// Velocity induced by a unit horseshoe using the legacy Leishman `R⁴` core.
///
/// `#[inline]` rather than `#[inline(always)]`: this is the O(N²) inner kernel
/// of influence assembly and must inline there, but it is also called from the
/// gradient kernels where forcing it would only add code bloat.
#[must_use]
#[inline]
pub fn horseshoe_vortex_velocity(
    field_point: &Vec3,
    bound_left: &Vec3,
    bound_right: &Vec3,
    prandtl_glauert_factor: Scalar,
    core_radius: Scalar,
    include_bound: bool,
) -> Vec3 {
    let a = Vec3::new(
        (bound_left.x - field_point.x) / prandtl_glauert_factor,
        bound_left.y - field_point.y,
        bound_left.z - field_point.z,
    );
    let b = Vec3::new(
        (bound_right.x - field_point.x) / prandtl_glauert_factor,
        bound_right.y - field_point.y,
        bound_right.z - field_point.z,
    );
    let a_sq = a.norm_squared();
    let b_sq = b.norm_squared();
    let a_mag = a_sq.sqrt();
    let b_mag = b_sq.sqrt();
    let core_sq = core_radius * core_radius;
    let core_4 = core_sq * core_sq;

    let (mut u, mut v, mut w) = (0.0, 0.0, 0.0);
    if include_bound && a_mag * b_mag > EPSILON {
        let cross = a.cross(&b);
        let cross_sq = cross.norm_squared();
        if cross_sq > EPSILON {
            let dot = a.dot(&b);
            let length_sq = a_sq + b_sq - 2.0 * dot;
            let factor = ((b_sq - dot) / (b_sq * b_sq + core_4).sqrt().sqrt()
                + (a_sq - dot) / (a_sq * a_sq + core_4).sqrt().sqrt())
                / (cross_sq * cross_sq + length_sq * length_sq * core_4).sqrt();
            u = cross.x * factor;
            v = cross.y * factor;
            w = cross.z * factor;
        }
    }
    // The two trailing legs are summed with each other first, then added to the
    // bound contribution, rather than accumulated into it one leg at a time.
    //
    // This makes the kernel equivariant under a y-reflection in floating point,
    // not merely in exact arithmetic. A mirrored panel is evaluated with its
    // bound endpoints swapped, so `(a, b)` become `(M·b, M·a)`; accumulating
    // leg-at-a-time then associates the same three terms as
    // `((bound + b_leg) + a_leg)` instead of `((bound + a_leg) + b_leg)`, and
    // addition is commutative but not associative. Summing the legs pairwise
    // removes the asymmetry, because `a_leg + b_leg` and `b_leg + a_leg` agree
    // exactly.
    //
    // `0.0 - x` rather than `-x` for the `w` legs: the two are identical for
    // every finite non-zero value, but unary negation maps `+0.0` to `-0.0`,
    // which would leave the reflected and directly evaluated forms differing in
    // the sign of zero.
    let (mut leg_a_v, mut leg_a_w) = (0.0, 0.0);
    if a_mag > EPSILON {
        let axis_sq = a.y * a.y + a.z * a.z;
        let factor = -(1.0 - a.x / a_mag) / (axis_sq * axis_sq + core_4).sqrt();
        leg_a_v = a.z * factor;
        leg_a_w = 0.0 - a.y * factor;
    }
    let (mut leg_b_v, mut leg_b_w) = (0.0, 0.0);
    if b_mag > EPSILON {
        let axis_sq = b.y * b.y + b.z * b.z;
        let factor = (1.0 - b.x / b_mag) / (axis_sq * axis_sq + core_4).sqrt();
        leg_b_v = b.z * factor;
        leg_b_w = 0.0 - b.y * factor;
    }
    v += leg_a_v + leg_b_v;
    w += leg_a_w + leg_b_w;
    Vec3::new(
        u * INV_4_PI / prandtl_glauert_factor,
        v * INV_4_PI,
        w * INV_4_PI,
    )
}

/// Target-side data packed to 56 bytes from the 224-byte `HorseshoeVortex`.
/// Assembly sweeps all N targets once per source column, so narrowing the
/// stride keeps that swept working set within L1.
struct Target {
    control_point: Vec3,
    normal: Vec3,
    component_index: usize,
}

/// Output of one assembly pass.
struct Assembly {
    velocity: MatX,
    /// Present only when assembled with `WITH_AIC`.
    aic: Option<MatX>,
    /// Whether every normalwash entry is finite. Trivially `true` when the
    /// normalwash matrix was not requested.
    aic_all_finite: bool,
}

/// Panel count below which the `parallel` feature keeps to the serial driver.
///
/// Assembly is `O(N²)` while the cost of distributing the columns is `O(N)`, so
/// there is a size under which spawning does not pay for itself. Measured on an
/// M1 Pro, 96 panels was ~5% slower parallel while 144 was already 1.3x faster;
/// the threshold sits between them. Both sides of the branch run the same
/// per-column closure, so this affects speed only, never results.
#[cfg(feature = "parallel")]
const PARALLEL_MIN_PANELS: usize = 128;

/// Runs `fill` over every source column, serially.
///
/// The accumulation is deliberately written as a non-short-circuiting `&=`
/// rather than `Iterator::all`. `fill` is called for its side effect of writing
/// the column, so short-circuiting on the first non-finite entry would leave
/// the remaining columns unwritten.
fn drive_columns_serial<F>(
    velocity: &mut [Scalar],
    aic: &mut [Scalar],
    n: usize,
    with_aic: bool,
    fill: F,
) -> bool
where
    F: Fn(usize, &mut [Scalar], &mut [Scalar]) -> bool + Sync,
{
    let mut all_finite = true;
    if with_aic {
        for (j, (velocity_column, aic_column)) in velocity
            .chunks_mut(3 * n)
            .zip(aic.chunks_mut(n))
            .enumerate()
        {
            all_finite &= fill(j, velocity_column, aic_column);
        }
    } else {
        for (j, velocity_column) in velocity.chunks_mut(3 * n).enumerate() {
            all_finite &= fill(j, velocity_column, &mut []);
        }
    }
    all_finite
}

/// Runs `fill` over every source column, serially. The `parallel` feature
/// swaps in the rayon variant below; both hand `fill` the same arguments, so
/// the two builds evaluate identical code per column and differ only in the
/// order independent columns are visited.
#[cfg(not(feature = "parallel"))]
fn drive_columns<F>(
    velocity: &mut [Scalar],
    aic: &mut [Scalar],
    n: usize,
    with_aic: bool,
    fill: F,
) -> bool
where
    F: Fn(usize, &mut [Scalar], &mut [Scalar]) -> bool + Sync,
{
    drive_columns_serial(velocity, aic, n, with_aic, fill)
}

/// Runs `fill` over every source column, one rayon task per column.
///
/// Each column writes a disjoint slice of both outputs and reads nothing that
/// another column writes, so the result is bit-identical to the serial driver
/// and independent of the thread count. The reduction is a boolean `and`, which
/// is order-independent.
#[cfg(feature = "parallel")]
fn drive_columns<F>(
    velocity: &mut [Scalar],
    aic: &mut [Scalar],
    n: usize,
    with_aic: bool,
    fill: F,
) -> bool
where
    F: Fn(usize, &mut [Scalar], &mut [Scalar]) -> bool + Sync,
{
    use rayon::prelude::*;

    if n < PARALLEL_MIN_PANELS {
        return drive_columns_serial(velocity, aic, n, with_aic, fill);
    }

    if with_aic {
        velocity
            .par_chunks_mut(3 * n)
            .zip(aic.par_chunks_mut(n))
            .enumerate()
            .map(|(j, (velocity_column, aic_column))| fill(j, velocity_column, aic_column))
            .reduce(|| true, |all, finite| all && finite)
    } else {
        velocity
            .par_chunks_mut(3 * n)
            .enumerate()
            .map(|(j, velocity_column)| fill(j, velocity_column, &mut []))
            .reduce(|| true, |all, finite| all && finite)
    }
}

/// Builds velocity and normalwash influence matrices for a fixed mesh.
pub struct InfluenceMatrixBuilder<'a> {
    mesh: &'a VortexLattice,
    mach_number: Scalar,
    prandtl_glauert_factor: Scalar,
    symmetry: SymmetryConfig,
    core_radius_chord_fraction: Scalar,
    core_radius_width_fraction: Scalar,
}

impl<'a> InfluenceMatrixBuilder<'a> {
    #[must_use]
    pub fn new(mesh: &'a VortexLattice) -> Self {
        Self {
            mesh,
            mach_number: 0.0,
            prandtl_glauert_factor: 1.0,
            symmetry: SymmetryConfig::default(),
            core_radius_chord_fraction: 0.0001,
            core_radius_width_fraction: 0.0001,
        }
    }

    pub fn set_mach_number(&mut self, mach: Scalar) {
        self.mach_number = mach;
        self.prandtl_glauert_factor = if mach < 1.0 {
            (1.0 - mach * mach).sqrt()
        } else {
            1.0
        };
    }

    #[must_use]
    pub fn mach_number(&self) -> Scalar {
        self.mach_number
    }

    #[must_use]
    pub fn prandtl_glauert_factor(&self) -> Scalar {
        self.prandtl_glauert_factor
    }

    pub fn set_symmetry(&mut self, symmetry: SymmetryConfig) {
        self.symmetry = symmetry;
    }

    #[must_use]
    pub fn symmetry(&self) -> SymmetryConfig {
        self.symmetry
    }

    pub fn set_core_radius_chord_fraction(&mut self, fraction: Scalar) {
        self.core_radius_chord_fraction = fraction;
    }

    pub fn set_core_radius_width_fraction(&mut self, fraction: Scalar) {
        self.core_radius_width_fraction = fraction;
    }

    fn core_radius(&self, vortex_index: usize, same_component: bool) -> Scalar {
        let vortex = &self.mesh.vortices()[vortex_index];
        let panel_width = (vortex.bound_vortex_right - vortex.bound_vortex_left).norm();
        let mut core = 0.0001 * panel_width;
        if !same_component {
            core = (self.core_radius_chord_fraction * vortex.local_chord)
                .max(self.core_radius_width_fraction * panel_width);
        }
        core
    }

    #[must_use]
    pub fn single_vortex_velocity(
        &self,
        field_point: &Vec3,
        vortex_index: usize,
        include_bound: bool,
    ) -> Vec3 {
        let vortex = &self.mesh.vortices()[vortex_index];
        horseshoe_vortex_velocity(
            field_point,
            &vortex.bound_vortex_left,
            &vortex.bound_vortex_right,
            self.prandtl_glauert_factor,
            self.core_radius(vortex_index, true),
            include_bound,
        )
    }

    fn y_image_velocity(&self, field_point: &Vec3, vortex_index: usize, core: Scalar) -> Vec3 {
        let vortex = &self.mesh.vortices()[vortex_index];
        let offset = 2.0 * self.symmetry.y_symmetry_plane;
        let left = Vec3::new(
            vortex.bound_vortex_right.x,
            offset - vortex.bound_vortex_right.y,
            vortex.bound_vortex_right.z,
        );
        let right = Vec3::new(
            vortex.bound_vortex_left.x,
            offset - vortex.bound_vortex_left.y,
            vortex.bound_vortex_left.z,
        );
        let include_bound = !(self.symmetry.y_symmetry == SymmetryType::Symmetric
            && (field_point - 0.5 * (left + right)).norm_squared() < EPSILON);
        self.symmetry.y_symmetry.factor()
            * horseshoe_vortex_velocity(
                field_point,
                &left,
                &right,
                self.prandtl_glauert_factor,
                core,
                include_bound,
            )
    }

    fn z_image_velocity(&self, field_point: &Vec3, vortex_index: usize, core: Scalar) -> Vec3 {
        let vortex = &self.mesh.vortices()[vortex_index];
        let offset = 2.0 * self.symmetry.z_symmetry_plane;
        let left = Vec3::new(
            vortex.bound_vortex_right.x,
            vortex.bound_vortex_right.y,
            offset - vortex.bound_vortex_right.z,
        );
        let right = Vec3::new(
            vortex.bound_vortex_left.x,
            vortex.bound_vortex_left.y,
            offset - vortex.bound_vortex_left.z,
        );
        self.symmetry.z_symmetry.factor()
            * horseshoe_vortex_velocity(
                field_point,
                &left,
                &right,
                self.prandtl_glauert_factor,
                core,
                true,
            )
    }

    /// Column-oriented influence assembly, optionally fused with the
    /// normalwash contraction.
    ///
    /// `MatX` is column-major, so iterating source panels (`j`) on the outside
    /// and target panels (`i`) on the inside makes both outputs stream
    /// sequentially: velocity at `j*3n + 3i + k`, normalwash at `j*n + i`. The
    /// previous row-major traversal strode `24n` bytes per element and then
    /// needed a second full pass over the velocity matrix (31.8 MB at N=1152)
    /// to contract it against the panel normals.
    ///
    /// `WITH_AIC` is a const generic so the normalwash half compiles out
    /// entirely for callers that only want the velocity matrix.
    ///
    /// Bit-exactness: every output element is produced by the same operations
    /// on the same operands in the same order as the two-pass version. Only
    /// the traversal order over mutually independent elements changes, plus
    /// hoisting of values that depend solely on the source panel.
    fn assemble<const WITH_AIC: bool>(&self) -> Assembly {
        let vortices = self.mesh.vortices();
        let n = vortices.len();
        if n == 0 {
            return Assembly {
                velocity: MatX::zeros(0, 0),
                aic: WITH_AIC.then(|| MatX::zeros(0, 0)),
                aic_all_finite: true,
            };
        }

        let targets: Vec<Target> = vortices
            .iter()
            .map(|vortex| Target {
                control_point: vortex.control_point,
                normal: vortex.normal_at_control_point,
                component_index: vortex.component_index,
            })
            .collect();

        let mut velocity_data = vec![0.0; 3 * n * n];
        let mut aic_data = if WITH_AIC {
            vec![0.0; n * n]
        } else {
            Vec::new()
        };

        // One source column. Returns whether every normalwash entry it wrote is
        // finite, which `rebuild_aic` previously established with a separate
        // full pass over the assembled matrix.
        //
        // Taking the column slices as arguments rather than indexing captured
        // buffers is what lets the parallel driver hand each invocation a
        // disjoint borrow.
        let fill_column =
            |j: usize, velocity_column: &mut [Scalar], aic_column: &mut [Scalar]| -> bool {
                let source = &vortices[j];
                let left = source.bound_vortex_left;
                let right = source.bound_vortex_right;

                // `core_radius` recomputed this width (and its `sqrt`) on every one
                // of the N² elements; it depends only on the source panel. The
                // 0.0001 literal is the original same-component constant, which is
                // deliberately not `core_radius_chord_fraction`.
                let panel_width = (right - left).norm();
                let core_same = 0.0001 * panel_width;
                let core_other = (self.core_radius_chord_fraction * source.local_chord)
                    .max(self.core_radius_width_fraction * panel_width);

                // Image geometry likewise depends only on the source panel. Note
                // the left/right swap in the y and z images and its deliberate
                // absence in the yz image, both preserved from the originals.
                let y_image = (self.symmetry.y_symmetry != SymmetryType::None).then(|| {
                    let offset = 2.0 * self.symmetry.y_symmetry_plane;
                    let image_left = Vec3::new(right.x, offset - right.y, right.z);
                    let image_right = Vec3::new(left.x, offset - left.y, left.z);
                    let midpoint = 0.5 * (image_left + image_right);
                    (
                        image_left,
                        image_right,
                        midpoint,
                        self.symmetry.y_symmetry.factor(),
                    )
                });
                let z_image = (self.symmetry.z_symmetry != SymmetryType::None).then(|| {
                    let offset = 2.0 * self.symmetry.z_symmetry_plane;
                    (
                        Vec3::new(right.x, right.y, offset - right.z),
                        Vec3::new(left.x, left.y, offset - left.z),
                        self.symmetry.z_symmetry.factor(),
                    )
                });
                let yz_image = (self.symmetry.z_symmetry != SymmetryType::None
                    && self.symmetry.y_symmetry != SymmetryType::None)
                    .then(|| {
                        let y_offset = 2.0 * self.symmetry.y_symmetry_plane;
                        let z_offset = 2.0 * self.symmetry.z_symmetry_plane;
                        (
                            Vec3::new(left.x, y_offset - left.y, z_offset - left.z),
                            Vec3::new(right.x, y_offset - right.y, z_offset - right.z),
                            // `a * b * v` parses as `(a * b) * v`, so folding the
                            // two scalar factors is exact.
                            self.symmetry.y_symmetry.factor() * self.symmetry.z_symmetry.factor(),
                        )
                    });

                let mut all_finite = true;
                for (i, target) in targets.iter().enumerate() {
                    let field_point = target.control_point;
                    let core = if target.component_index == source.component_index {
                        core_same
                    } else {
                        core_other
                    };

                    let mut velocity = horseshoe_vortex_velocity(
                        &field_point,
                        &left,
                        &right,
                        self.prandtl_glauert_factor,
                        core,
                        i != j,
                    );
                    if let Some((image_left, image_right, midpoint, factor)) = y_image {
                        let include_bound = !(self.symmetry.y_symmetry == SymmetryType::Symmetric
                            && (field_point - midpoint).norm_squared() < EPSILON);
                        velocity += factor
                            * horseshoe_vortex_velocity(
                                &field_point,
                                &image_left,
                                &image_right,
                                self.prandtl_glauert_factor,
                                core,
                                include_bound,
                            );
                    }
                    if let Some((image_left, image_right, factor)) = z_image {
                        velocity += factor
                            * horseshoe_vortex_velocity(
                                &field_point,
                                &image_left,
                                &image_right,
                                self.prandtl_glauert_factor,
                                core,
                                true,
                            );
                        if let Some((image_left, image_right, factor)) = yz_image {
                            velocity += factor
                                * horseshoe_vortex_velocity(
                                    &field_point,
                                    &image_left,
                                    &image_right,
                                    self.prandtl_glauert_factor,
                                    core,
                                    true,
                                );
                        }
                    }

                    velocity_column[3 * i] = velocity.x;
                    velocity_column[3 * i + 1] = velocity.y;
                    velocity_column[3 * i + 2] = velocity.z;
                    if WITH_AIC {
                        let normal = target.normal;
                        let entry =
                            velocity.x * normal.x + velocity.y * normal.y + velocity.z * normal.z;
                        aic_column[i] = entry;
                        all_finite &= entry.is_finite();
                    }
                }
                all_finite
            };

        let aic_all_finite =
            drive_columns(&mut velocity_data, &mut aic_data, n, WITH_AIC, fill_column);

        Assembly {
            velocity: MatX::from_vec(3 * n, n, velocity_data),
            aic: WITH_AIC.then(|| MatX::from_vec(n, n, aic_data)),
            aic_all_finite,
        }
    }

    #[must_use]
    pub fn velocity_influence_matrix(&self) -> MatX {
        self.assemble::<false>().velocity
    }

    #[must_use]
    pub fn normalwash_matrix(&self) -> MatX {
        self.normalwash_and_velocity().0
    }

    #[must_use]
    pub fn normalwash_and_velocity(&self) -> (MatX, MatX) {
        let (aic, velocity, _) = self.normalwash_and_velocity_checked();
        (aic, velocity)
    }

    /// As [`Self::normalwash_and_velocity`], plus whether every normalwash
    /// entry is finite.
    ///
    /// The flag is accumulated while the matrix is written, so callers that
    /// need it — [`super::VortexLatticeSolver::rebuild_aic`] does, to reject a
    /// singular system — avoid a second full pass over the assembled matrix.
    #[must_use]
    pub fn normalwash_and_velocity_checked(&self) -> (MatX, MatX, bool) {
        let assembly = self.assemble::<true>();
        (
            assembly
                .aic
                .expect("assemble::<true> always produces the normalwash matrix"),
            assembly.velocity,
            assembly.aic_all_finite,
        )
    }

    #[must_use]
    pub fn induced_velocity(&self, field_point: &Vec3, circulation: &VecX) -> Vec3 {
        let mut total = Vec3::zeros();
        for (j, source) in self.mesh.vortices().iter().enumerate() {
            let core = self.core_radius(j, false);
            let mut velocity = horseshoe_vortex_velocity(
                field_point,
                &source.bound_vortex_left,
                &source.bound_vortex_right,
                self.prandtl_glauert_factor,
                core,
                true,
            );
            if self.symmetry.y_symmetry != SymmetryType::None {
                velocity += self.y_image_velocity(field_point, j, core);
            }
            if self.symmetry.z_symmetry != SymmetryType::None {
                velocity += self.z_image_velocity(field_point, j, core);
            }
            total += circulation[j] * velocity;
        }
        total
    }
}
