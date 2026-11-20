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

    // Degeneracy tests are dimensionless ratios, not comparisons against an
    // absolute length.
    //
    // `cross_sq` is a length⁴ quantity and `a_mag * b_mag` a length², so
    // testing either against a fixed `EPSILON` makes the kernel depend on the
    // unit the model is expressed in: at millimetre scale `|a×b|²` falls below
    // 1e-12 for ordinary panels and the bound leg is silently dropped, which
    // moved every AIC entry by several percent with no error reported.
    //
    // `cross_sq > EPSILON * a_sq * b_sq` is `sin²θ > EPSILON`, the collinearity
    // test that was intended. It is invariant under a uniform rescaling of the
    // geometry and agrees with the previous threshold at unit scale, where
    // `a_sq * b_sq ≈ 1`. `core_radius` is itself derived from panel size, so
    // the `core_4` terms already scale correctly and need no change.
    let (mut u, mut v, mut w) = (0.0, 0.0, 0.0);
    if include_bound && a_sq > 0.0 && b_sq > 0.0 {
        let cross = a.cross(&b);
        let cross_sq = cross.norm_squared();
        if cross_sq > EPSILON * a_sq * b_sq {
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
    // A trailing leg only has to avoid `0 / 0`. `a.x / a_mag` is a direction
    // cosine bounded by 1 for any non-zero `a`, and the denominator is bounded
    // below by `core_sq`, so a merely small `a_mag` is well behaved and only an
    // exactly coincident endpoint is degenerate. Testing `a_sq > 0.0` rather
    // than a fixed length keeps the legs scale-invariant too.
    let (mut leg_a_v, mut leg_a_w) = (0.0, 0.0);
    if a_sq > 0.0 {
        let axis_sq = a.y * a.y + a.z * a.z;
        let factor = -(1.0 - a.x / a_mag) / (axis_sq * axis_sq + core_4).sqrt();
        leg_a_v = a.z * factor;
        leg_a_w = 0.0 - a.y * factor;
    }
    let (mut leg_b_v, mut leg_b_w) = (0.0, 0.0);
    if b_sq > 0.0 {
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

/// A contiguous run of source panels immediately followed by its y-mirrored
/// copy. Columns `[first, first + count)` are the originals and
/// `[first + count, first + 2 * count)` their mirrors, paired by offset.
struct MirrorBlock {
    first: usize,
    count: usize,
}

/// Pairing that lets assembly derive half its entries by reflection.
struct MirrorPlan {
    blocks: Vec<MirrorBlock>,
    /// Partner index of every panel; an involution over `0..n`.
    partner: Vec<usize>,
    /// Whether each panel is an original (rather than a mirrored copy).
    is_original: Vec<bool>,
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

/// Fills the mirrored rows of a column pair by reflection.
///
/// `fill_column` has already written the original rows of both columns. With
/// `mi` a mirrored panel and `i = partner[mi]` its original:
///
/// ```text
/// V(mi, mj) = M · V(i, j)        and        V(mi, j) = M · V(i, mj)
/// ```
///
/// so each column's mirrored half is the reflection of its sibling's original
/// half. `M = diag(1, -1, 1)` negates only `y`, which is exact.
///
/// The normalwash entries need no reflection at all: contracting a reflected
/// velocity with a reflected normal gives `vx*nx + (-vy)*(-ny) + vz*nz`, the
/// same three products summed in the same order as the unreflected entry, so
/// `aic(mi, mj) == aic(i, j)` bit-for-bit.
fn derive_mirrored_rows(
    velocity_j: &mut [Scalar],
    velocity_mj: &mut [Scalar],
    aic_j: &mut [Scalar],
    aic_mj: &mut [Scalar],
    plan: &MirrorPlan,
) {
    for (mirror, &original) in plan.partner.iter().enumerate() {
        if plan.is_original[mirror] {
            continue;
        }
        // `0.0 - y` rather than `-y`: both negate every finite non-zero value
        // identically, but unary negation maps `+0.0` to `-0.0` while direct
        // evaluation produces `+0.0`. Subtracting from zero normalizes the
        // sign, keeping derived and evaluated entries bit-identical.
        velocity_mj[3 * mirror] = velocity_j[3 * original];
        velocity_mj[3 * mirror + 1] = 0.0 - velocity_j[3 * original + 1];
        velocity_mj[3 * mirror + 2] = velocity_j[3 * original + 2];

        velocity_j[3 * mirror] = velocity_mj[3 * original];
        velocity_j[3 * mirror + 1] = 0.0 - velocity_mj[3 * original + 1];
        velocity_j[3 * mirror + 2] = velocity_mj[3 * original + 2];

        aic_mj[mirror] = aic_j[original];
        aic_j[mirror] = aic_mj[original];
    }
}

/// Assembles a mesh that admits reflection reuse, halving kernel evaluations.
///
/// Columns are visited in `(original, mirror)` pairs. Each pair evaluates the
/// kernel only for original target rows in both of its columns and derives the
/// mirrored rows, so the pair costs `n` evaluations instead of `2n`.
///
/// Pairing across columns is why this cannot reuse [`drive_columns`]: a
/// mirrored entry depends on the sibling column, not on its own. Pairs remain
/// mutually independent, so the decomposition stays safe to parallelize; the
/// blocks tile the column space, so `split_at_mut` can hand each pair disjoint
/// borrows of both outputs.
fn drive_mirrored_columns<F>(
    velocity: &mut [Scalar],
    aic: &mut [Scalar],
    n: usize,
    plan: &MirrorPlan,
    fill: F,
) -> bool
where
    F: Fn(usize, &mut [Scalar], &mut [Scalar]) -> bool + Sync,
{
    #[cfg(feature = "parallel")]
    use rayon::prelude::*;

    let mut all_finite = true;
    let mut velocity_rest = velocity;
    let mut aic_rest = aic;
    let mut consumed = 0usize;

    for block in &plan.blocks {
        debug_assert_eq!(block.first, consumed, "mirror blocks must tile the columns");
        let (velocity_block, velocity_tail) = velocity_rest.split_at_mut(2 * block.count * 3 * n);
        let (aic_block, aic_tail) = aic_rest.split_at_mut(2 * block.count * n);
        velocity_rest = velocity_tail;
        aic_rest = aic_tail;

        let (velocity_a, velocity_b) = velocity_block.split_at_mut(block.count * 3 * n);
        let (aic_a, aic_b) = aic_block.split_at_mut(block.count * n);

        let first = block.first;
        let count = block.count;
        let pair = |offset: usize,
                    velocity_j: &mut [Scalar],
                    velocity_mj: &mut [Scalar],
                    aic_j: &mut [Scalar],
                    aic_mj: &mut [Scalar]|
         -> bool {
            let j = first + offset;
            let mj = first + count + offset;
            let finite = fill(j, velocity_j, aic_j) & fill(mj, velocity_mj, aic_mj);
            derive_mirrored_rows(velocity_j, velocity_mj, aic_j, aic_mj, plan);
            finite
        };

        #[cfg(feature = "parallel")]
        let block_finite = if count >= PARALLEL_MIN_PANELS {
            velocity_a
                .par_chunks_mut(3 * n)
                .zip(velocity_b.par_chunks_mut(3 * n))
                .zip(aic_a.par_chunks_mut(n).zip(aic_b.par_chunks_mut(n)))
                .enumerate()
                .map(|(offset, ((velocity_j, velocity_mj), (aic_j, aic_mj)))| {
                    pair(offset, velocity_j, velocity_mj, aic_j, aic_mj)
                })
                .reduce(|| true, |all, finite| all && finite)
        } else {
            serial_pairs(velocity_a, velocity_b, aic_a, aic_b, n, pair)
        };

        #[cfg(not(feature = "parallel"))]
        let block_finite = serial_pairs(velocity_a, velocity_b, aic_a, aic_b, n, pair);

        all_finite &= block_finite;
        consumed += 2 * block.count;
    }

    all_finite
}

/// Serial column-pair sweep shared by both drivers.
fn serial_pairs<P>(
    velocity_a: &mut [Scalar],
    velocity_b: &mut [Scalar],
    aic_a: &mut [Scalar],
    aic_b: &mut [Scalar],
    n: usize,
    pair: P,
) -> bool
where
    P: Fn(usize, &mut [Scalar], &mut [Scalar], &mut [Scalar], &mut [Scalar]) -> bool + Sync,
{
    let mut all_finite = true;
    for (offset, ((velocity_j, velocity_mj), (aic_j, aic_mj))) in velocity_a
        .chunks_mut(3 * n)
        .zip(velocity_b.chunks_mut(3 * n))
        .zip(aic_a.chunks_mut(n).zip(aic_b.chunks_mut(n)))
        .enumerate()
    {
        all_finite &= pair(offset, velocity_j, velocity_mj, aic_j, aic_mj);
    }
    all_finite
}

/// Builds velocity and normalwash influence matrices for a fixed mesh.
pub struct InfluenceMatrixBuilder<'a> {
    mesh: &'a VortexLattice,
    mach_number: Scalar,
    prandtl_glauert_factor: Scalar,
    symmetry: SymmetryConfig,
    core_radius_chord_fraction: Scalar,
    core_radius_width_fraction: Scalar,
    mirror_reuse: bool,
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
            mirror_reuse: false,
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

    /// Enables deriving mirrored influence entries by reflection instead of
    /// evaluating them, for meshes built with `y_duplicated` surfaces.
    ///
    /// Off by default. Enabling it never changes results: the fast path is
    /// taken only when the mesh is verified to satisfy the reflection identity
    /// exactly, and the derived entries are then bit-identical to evaluated
    /// ones. A mesh that does not qualify — symmetry images in use, partially
    /// duplicated geometry, or a mirror plane that is not a pure sign flip —
    /// falls back to evaluating every entry.
    ///
    /// Use [`Self::mirror_reuse_applies`] to tell which of the two will happen.
    pub fn set_mirror_reuse(&mut self, enabled: bool) {
        self.mirror_reuse = enabled;
    }

    #[must_use]
    pub fn mirror_reuse(&self) -> bool {
        self.mirror_reuse
    }

    /// Whether reflection reuse is both enabled and applicable to this mesh.
    ///
    /// [`Self::set_mirror_reuse`] is a request, not a guarantee — a mesh that
    /// cannot satisfy the reflection identity exactly falls back to evaluating
    /// every entry. This reports which of the two will happen, so callers can
    /// tell whether enabling it actually bought anything.
    #[must_use]
    pub fn mirror_reuse_applies(&self) -> bool {
        self.mirror_plan().is_some()
    }

    /// Detects whether this mesh admits reflection reuse, and how.
    ///
    /// Assembly can halve its kernel evaluations when the lattice is built from
    /// `y_duplicated` surfaces, because the induced velocity is equivariant
    /// under the reflection `M = diag(1, -1, 1)`:
    ///
    /// ```text
    /// V(mirror(i), mirror(j)) = M · V(i, j)
    /// ```
    ///
    /// That identity holds *exactly* in floating point only under conditions
    /// this function verifies rather than assumes:
    ///
    /// - **No symmetry images.** The image kernels add terms whose own
    ///   equivariance depends on the image plane, so they are excluded outright.
    /// - **Every panel is paired**, in contiguous original/mirror blocks. The
    ///   lattice appends a duplicate immediately after its original, so blocks
    ///   are contiguous; a partially duplicated mesh is rejected.
    /// - **The reflection is a pure sign flip.** `reflect_swapped` maps `y` to
    ///   `2*plane - y`, which rounds for a non-zero plane and would make the
    ///   derived entries disagree with evaluated ones. Rather than trusting the
    ///   plane metadata, this checks the resulting geometry directly: each
    ///   mirrored panel must satisfy the reflection relation exactly, including
    ///   the deliberate left/right endpoint swap.
    ///
    /// The whole check is `O(N)` against an `O(N²)` assembly, so validating is
    /// far cheaper than the work it guards.
    fn mirror_plan(&self) -> Option<MirrorPlan> {
        if !self.mirror_reuse
            || self.symmetry.y_symmetry != SymmetryType::None
            || self.symmetry.z_symmetry != SymmetryType::None
        {
            return None;
        }

        let vortices = self.mesh.vortices();
        let n = vortices.len();
        if n == 0 || !n.is_multiple_of(2) {
            return None;
        }

        let surfaces = self.mesh.surface_info();
        let mut blocks = Vec::new();
        let mut partner = vec![usize::MAX; n];
        let mut is_original = vec![false; n];

        for surface in surfaces {
            if !surface.is_y_duplicated {
                continue;
            }
            let original = surfaces.get(surface.original_surface_index)?;
            // The duplicate must directly follow its original and match it in
            // size, which is how `duplicate_surface` lays them out.
            if original.num_vortices != surface.num_vortices
                || original.first_vortex_index + original.num_vortices != surface.first_vortex_index
            {
                return None;
            }
            blocks.push(MirrorBlock {
                first: original.first_vortex_index,
                count: original.num_vortices,
            });
            for offset in 0..original.num_vortices {
                let a = original.first_vortex_index + offset;
                let b = surface.first_vortex_index + offset;
                if a >= n || b >= n || partner[a] != usize::MAX || partner[b] != usize::MAX {
                    return None;
                }
                partner[a] = b;
                partner[b] = a;
                is_original[a] = true;
            }
        }

        // Every panel must belong to a pair; a mesh mixing duplicated and
        // plain surfaces has no usable pairing and falls back.
        if blocks.is_empty() || partner.contains(&usize::MAX) {
            return None;
        }

        // The blocks must tile `0..n` contiguously and in ascending order.
        // `drive_mirrored_columns` walks them with `split_at_mut`, so a gap or
        // an out-of-order block would either panic there or, worse, write a
        // column into the wrong block. That held only because both generators
        // append each duplicate immediately after its original; checking it
        // here makes the guarantee total instead of incidental, and costs
        // `O(blocks)` against an `O(N²)` assembly.
        let mut consumed = 0usize;
        for block in &blocks {
            if block.first != consumed {
                return None;
            }
            consumed += 2 * block.count;
        }
        if consumed != vortices.len() {
            return None;
        }

        // Verify the reflection relation holds exactly for every pair.
        for (a, &b) in partner.iter().enumerate() {
            if !is_original[a] {
                continue;
            }
            let original = &vortices[a];
            let mirrored = &vortices[b];
            let reflected = |point: Vec3| Vec3::new(point.x, -point.y, point.z);
            if mirrored.control_point != reflected(original.control_point)
                || mirrored.normal_at_control_point != reflected(original.normal_at_control_point)
                // `duplicate_surface` swaps the bound endpoints so the mirrored
                // panel keeps a consistent circulation orientation.
                || mirrored.bound_vortex_left != reflected(original.bound_vortex_right)
                || mirrored.bound_vortex_right != reflected(original.bound_vortex_left)
                || mirrored.local_chord != original.local_chord
                || mirrored.component_index != original.component_index
            {
                return None;
            }
        }

        Some(MirrorPlan {
            blocks,
            partner,
            is_original,
        })
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

    /// Velocity from the image reflected in *both* symmetry planes.
    ///
    /// Present only when y and z symmetry are both active. Unlike the single
    /// reflections, the bound endpoints are *not* swapped: reflecting twice
    /// restores the original orientation, so swapping would reverse it. The
    /// two scalar factors are folded before scaling the vector, matching the
    /// assembly path exactly.
    fn yz_image_velocity(&self, field_point: &Vec3, vortex_index: usize, core: Scalar) -> Vec3 {
        let vortex = &self.mesh.vortices()[vortex_index];
        let y_offset = 2.0 * self.symmetry.y_symmetry_plane;
        let z_offset = 2.0 * self.symmetry.z_symmetry_plane;
        let left = Vec3::new(
            vortex.bound_vortex_left.x,
            y_offset - vortex.bound_vortex_left.y,
            z_offset - vortex.bound_vortex_left.z,
        );
        let right = Vec3::new(
            vortex.bound_vortex_right.x,
            y_offset - vortex.bound_vortex_right.y,
            z_offset - vortex.bound_vortex_right.z,
        );
        self.symmetry.y_symmetry.factor()
            * self.symmetry.z_symmetry.factor()
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

        // Reflection reuse applies only to the fused path. The velocity-only
        // path is not hot -- it serves the gradient kernels and tests, not
        // `rebuild_aic` -- so it is not worth a second set of drivers.
        let plan = if WITH_AIC { self.mirror_plan() } else { None };
        let originals: Option<&[bool]> = plan.as_ref().map(|plan| plan.is_original.as_slice());

        // One source column. Returns whether every normalwash entry it wrote is
        // finite, which `rebuild_aic` previously established with a separate
        // full pass over the assembled matrix.
        //
        // Taking the column slices as arguments rather than indexing captured
        // buffers is what lets the parallel driver hand each invocation a
        // disjoint borrow.
        //
        // Under a mirror plan this fills only the rows of original panels; the
        // mirrored rows are then derived by reflection.
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
                    if let Some(originals) = originals
                        && !originals[i]
                    {
                        continue;
                    }
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

        let aic_all_finite = match &plan {
            Some(plan) => {
                drive_mirrored_columns(&mut velocity_data, &mut aic_data, n, plan, fill_column)
            }
            None => drive_columns(&mut velocity_data, &mut aic_data, n, WITH_AIC, fill_column),
        };

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
                // Assembly adds the doubly-reflected image inside its z branch;
                // omitting it here made field-point velocities disagree with
                // the influence matrix by several percent, including a sign
                // flip in y, whenever both planes were active.
                if self.symmetry.y_symmetry != SymmetryType::None {
                    velocity += self.yz_image_velocity(field_point, j, core);
                }
            }
            total += circulation[j] * velocity;
        }
        total
    }
}
