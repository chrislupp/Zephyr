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

//! Operators evaluated at stored targets, with reflected source images on demand.
//!
//! Full mode hands assembly to the solver's tuned [`InfluenceMatrixBuilder`],
//! which evaluates exactly the kernel calls this module used to make, in the
//! same order per entry, so its matrices are the same bits. Half-wing mode
//! adds each source's reflection about `y = 0`, which the builder has no notion
//! of, so it keeps its own column driver here.
use super::Work;
use crate::core::{MatX, Scalar, Vec3};
use crate::solver::{InfluenceMatrixBuilder, VortexLattice, horseshoe_vortex_velocity};

pub(crate) fn reflect(v: Vec3) -> Vec3 {
    Vec3::new(v.x, -v.y, v.z)
}

/// Per-source data every `N²` sweep reads, hoisted out of its inner loop and
/// packed away from the 224-byte `HorseshoeVortex`.
///
/// The core radii are the legacy per-pair rule evaluated once per source: a
/// target on the same component sees `0.0001 * width`, any other target
/// `max(0.0001 * chord, 0.0001 * width)`. Hoisting changes no bits; the
/// width's `hypot` only ever depended on the source.
#[derive(Clone, Copy, Debug)]
pub(crate) struct SourcePanel {
    pub left: Vec3,
    pub right: Vec3,
    /// Endpoints of the reflected image about `y = 0`. Reflection reverses the
    /// circulation sense, so the image's left endpoint is the reflected right.
    pub image_left: Vec3,
    pub image_right: Vec3,
    core_same: Scalar,
    core_other: Scalar,
    component: usize,
}

impl SourcePanel {
    /// Core radius of this source seen from a target on `component`.
    #[inline]
    pub fn core(&self, component: usize) -> Scalar {
        if component == self.component {
            self.core_same
        } else {
            self.core_other
        }
    }

    /// Bound endpoints of the direct source (`image == 0`) or its reflection.
    #[inline]
    pub fn endpoints(&self, image: usize) -> (Vec3, Vec3) {
        if image == 0 {
            (self.left, self.right)
        } else {
            (self.image_left, self.image_right)
        }
    }
}

pub(crate) fn source_panels(mesh: &VortexLattice) -> Vec<SourcePanel> {
    mesh.vortices()
        .iter()
        .map(|source| {
            let left = source.bound_vortex_left;
            let right = source.bound_vortex_right;
            let width = (right.y - left.y).hypot(right.z - left.z);
            SourcePanel {
                left,
                right,
                image_left: reflect(right),
                image_right: reflect(left),
                core_same: 0.0001 * width,
                core_other: (0.0001 * source.local_chord).max(0.0001 * width),
                component: source.component_index,
            }
        })
        .collect()
}

/// Target-side data for one assembly pass.
struct Target {
    point: Vec3,
    normal: Vec3,
    component: usize,
}

/// Assembles the selective operators: the normalwash matrix and control-point
/// velocity operator (`aic`), and the bound-force velocity operator
/// (`velocity`), each only when requested.
///
/// `mirror_reuse` lets full-mode assembly derive mirrored entries by
/// reflection where the builder proves that exact; it never changes a bit.
/// Half-wing lattices have no mirrored panels, so there it does nothing.
#[allow(clippy::too_many_arguments)]
pub(crate) fn assemble(
    mesh: &VortexLattice,
    mach: Scalar,
    normals: &[Vec3],
    half: bool,
    aic: bool,
    velocity: bool,
    mirror_reuse: bool,
    work: &mut Work,
) -> (MatX, MatX, MatX) {
    let n = mesh.vortex_count();
    let empty = || MatX::zeros(0, 0);
    let mut evaluations_per_pass = n * n;
    let (matrix, control, velocities) = if half {
        let pg = (1.0 - mach * mach).sqrt();
        let sources = source_panels(mesh);
        let (matrix, control) = if aic {
            assemble_half(mesh, &sources, pg, normals, false)
        } else {
            (empty(), empty())
        };
        let velocities = if velocity {
            assemble_half(mesh, &sources, pg, normals, true).1
        } else {
            empty()
        };
        evaluations_per_pass *= 2;
        (matrix, control, velocities)
    } else {
        // Default symmetry: the selective path never assembled image planes.
        let mut builder = InfluenceMatrixBuilder::new(mesh);
        builder.set_mach_number(mach);
        builder.set_mirror_reuse(mirror_reuse);
        if builder.mirror_reuse_applies() {
            evaluations_per_pass /= 2;
        }
        let (matrix, control) = if aic {
            let (matrix, control, _) = builder.normalwash_and_control_checked();
            (
                contract_if_rotated(mesh, matrix, &control, normals),
                control,
            )
        } else {
            (empty(), empty())
        };
        let velocities = if velocity {
            builder.velocity_influence_matrix()
        } else {
            empty()
        };
        (matrix, control, velocities)
    };
    work.influence_evaluations += (usize::from(aic) + usize::from(velocity)) * evaluations_per_pass;
    (matrix, control, velocities)
}

/// The builder contracts with the lattice normals. When the requested normals
/// differ (rotated by a control deflection), recontract the stored
/// control-point velocities with them, associated exactly as the fused
/// contraction is, `vx*nx + vy*ny + vz*nz`, so each entry is the bits a fused
/// pass with those normals would have written.
fn contract_if_rotated(
    mesh: &VortexLattice,
    matrix: MatX,
    control: &MatX,
    normals: &[Vec3],
) -> MatX {
    let same = |a: &Vec3, b: &Vec3| {
        a.iter()
            .zip(b.iter())
            .all(|(x, y)| x.to_bits() == y.to_bits())
    };
    if mesh
        .vortices()
        .iter()
        .zip(normals)
        .all(|(vortex, normal)| same(&vortex.normal_at_control_point, normal))
    {
        return matrix;
    }
    let n = mesh.vortex_count();
    let mut data = vec![0.0; n * n];
    let velocity = control.as_slice();
    for (j, column) in data.chunks_mut(n.max(1)).enumerate() {
        let velocity_column = &velocity[3 * n * j..3 * n * (j + 1)];
        for (i, entry) in column.iter_mut().enumerate() {
            let normal = normals[i];
            *entry = velocity_column[3 * i] * normal.x
                + velocity_column[3 * i + 1] * normal.y
                + velocity_column[3 * i + 2] * normal.z;
        }
    }
    MatX::from_vec(n, n, data)
}

/// One half-wing pass: control points with the normalwash contraction
/// (`at_bound == false`) or bound centres without it, each source paired with
/// its reflected image. Per entry this makes the same two kernel calls in the
/// same order as before -- direct source first, image added to it -- with the
/// same core radius and bound flags, so the entries are unchanged bit for bit.
fn assemble_half(
    mesh: &VortexLattice,
    sources: &[SourcePanel],
    pg: Scalar,
    normals: &[Vec3],
    at_bound: bool,
) -> (MatX, MatX) {
    let n = mesh.vortex_count();
    let targets: Vec<Target> = mesh
        .vortices()
        .iter()
        .zip(normals)
        .map(|(vortex, normal)| Target {
            point: if at_bound {
                vortex.bound_vortex_center
            } else {
                vortex.control_point
            },
            normal: *normal,
            component: vortex.component_index,
        })
        .collect();
    let fill = |j: usize, column: &mut [Scalar], velocity_column: &mut [Scalar]| {
        let source = &sources[j];
        for (i, target) in targets.iter().enumerate() {
            let core = source.core(target.component);
            let mut v = horseshoe_vortex_velocity(
                &target.point,
                &source.left,
                &source.right,
                pg,
                core,
                !at_bound || i != j,
            );
            v += horseshoe_vortex_velocity(
                &target.point,
                &source.image_left,
                &source.image_right,
                pg,
                core,
                true,
            );
            if !column.is_empty() {
                let normal = target.normal;
                column[i] = v.x * normal.x + v.y * normal.y + v.z * normal.z;
            }
            velocity_column[3 * i] = v.x;
            velocity_column[3 * i + 1] = v.y;
            velocity_column[3 * i + 2] = v.z;
        }
    };
    let mut matrix = if at_bound {
        MatX::zeros(0, 0)
    } else {
        MatX::zeros(n, n)
    };
    let mut velocity = MatX::zeros(3 * n, n);
    if n > 0 {
        drive_columns(matrix.as_mut_slice(), velocity.as_mut_slice(), n, &fill);
    }
    (matrix, velocity)
}

// Match the legacy assembly threshold; small systems stay on one thread.
#[cfg(feature = "parallel")]
pub(super) const PARALLEL_MIN_PANELS: usize = 128;

/// Runs `fill` over every source column. `velocity` always holds `3n` rows
/// per column; `matrix` holds `n` per column, or is empty for a pass without
/// the normalwash contraction. Each column writes only its own slices and reads
/// nothing another column writes, so the parallel driver is bit-identical to
/// the serial one at every thread count.
fn drive_columns<F>(matrix: &mut [Scalar], velocity: &mut [Scalar], n: usize, fill: &F)
where
    F: Fn(usize, &mut [Scalar], &mut [Scalar]) + Sync,
{
    #[cfg(feature = "parallel")]
    if n >= PARALLEL_MIN_PANELS {
        use rayon::prelude::*;
        if matrix.is_empty() {
            velocity
                .par_chunks_mut(3 * n)
                .enumerate()
                .for_each(|(j, v)| fill(j, &mut [], v));
        } else {
            matrix
                .par_chunks_mut(n)
                .zip(velocity.par_chunks_mut(3 * n))
                .enumerate()
                .for_each(|(j, (a, v))| fill(j, a, v));
        }
        return;
    }
    if matrix.is_empty() {
        for (j, v) in velocity.chunks_mut(3 * n).enumerate() {
            fill(j, &mut [], v);
        }
    } else {
        for (j, (a, v)) in matrix
            .chunks_mut(n)
            .zip(velocity.chunks_mut(3 * n))
            .enumerate()
        {
            fill(j, a, v);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{MeshInput, SurfaceMesh};

    /// A wing and a tail on separate components, so both core-radius branches
    /// are taken, each duplicated about `y = 0` unless `half`.
    fn lattice(half: bool, ny: usize) -> VortexLattice {
        let surface = |name: &str, nx: usize, ny: usize, component: usize, offset: Vec3| {
            let mut surface = SurfaceMesh::new(name, nx, ny);
            surface.y_duplicated = !half;
            surface.component_index = component;
            for i in 0..=nx {
                for j in 0..=ny {
                    let eta = j as Scalar / ny as Scalar;
                    let chord = i as Scalar / nx as Scalar;
                    surface
                        .set_node(
                            i,
                            j,
                            offset
                                + Vec3::new(
                                    0.6 * eta + chord * (1.0 - 0.4 * eta),
                                    0.05 + 3.0 * eta,
                                    0.3 * eta - 0.05 * chord * eta,
                                ),
                        )
                        .unwrap();
                }
            }
            surface
        };
        let mut lattice = VortexLattice::default();
        lattice
            .generate_from_mesh(&MeshInput {
                surfaces: vec![
                    surface("wing", 3, ny, 0, Vec3::zeros()),
                    surface("tail", 2, 3, 1, Vec3::new(4.0, 0.0, 0.5)),
                ],
                ..MeshInput::default()
            })
            .unwrap();
        lattice
    }

    /// The assembly as it stood at the 0.5 performance baseline: one kernel
    /// call per entry with the core radius recomputed per pair.
    fn reference(
        mesh: &VortexLattice,
        mach: Scalar,
        normals: &[Vec3],
        half: bool,
        at_bound: bool,
    ) -> (Vec<Scalar>, Vec<Scalar>) {
        let n = mesh.vortex_count();
        let pg = (1.0 - mach * mach).sqrt();
        let core = |row: usize, column: usize| {
            let source = &mesh.vortices()[column];
            let width = (source.bound_vortex_right.y - source.bound_vortex_left.y)
                .hypot(source.bound_vortex_right.z - source.bound_vortex_left.z);
            if mesh.vortices()[row].component_index == source.component_index {
                0.0001 * width
            } else {
                (0.0001 * source.local_chord).max(0.0001 * width)
            }
        };
        let mut matrix = vec![0.0; n * n];
        let mut velocity = vec![0.0; 3 * n * n];
        for (j, source) in mesh.vortices().iter().enumerate() {
            for (i, target) in mesh.vortices().iter().enumerate() {
                let point = if at_bound {
                    target.bound_vortex_center
                } else {
                    target.control_point
                };
                let radius = core(i, j);
                let mut v = horseshoe_vortex_velocity(
                    &point,
                    &source.bound_vortex_left,
                    &source.bound_vortex_right,
                    pg,
                    radius,
                    !at_bound || i != j,
                );
                if half {
                    v += horseshoe_vortex_velocity(
                        &point,
                        &reflect(source.bound_vortex_right),
                        &reflect(source.bound_vortex_left),
                        pg,
                        radius,
                        true,
                    );
                }
                matrix[j * n + i] = v.x * normals[i].x + v.y * normals[i].y + v.z * normals[i].z;
                for k in 0..3 {
                    velocity[j * 3 * n + 3 * i + k] = v[k];
                }
            }
        }
        (matrix, velocity)
    }

    fn bits(values: &[Scalar]) -> Vec<u64> {
        values.iter().map(|value| value.to_bits()).collect()
    }

    /// Every operator configuration reproduces the baseline assembly bit for
    /// bit: full and half-wing, lattice and rotated normals, mirror reuse on
    /// and off, and each pass on its own or both together.
    #[test]
    fn assembly_matches_the_baseline_bit_for_bit() {
        for half in [false, true] {
            // Above the threading threshold in both modes.
            let mesh = lattice(half, if half { 44 } else { 22 });
            let n = mesh.vortex_count();
            assert!(n >= 128);
            let lattice_normals: Vec<_> = mesh
                .vortices()
                .iter()
                .map(|v| v.normal_at_control_point)
                .collect();
            let rotated: Vec<_> = lattice_normals
                .iter()
                .enumerate()
                .map(|(i, normal)| {
                    (normal + Vec3::new(0.02 * (i as Scalar).sin(), 0.0, 0.01)).normalize()
                })
                .collect();
            for normals in [&lattice_normals, &rotated] {
                let (matrix, control) = reference(&mesh, 0.35, normals, half, false);
                let (_, bound) = reference(&mesh, 0.35, normals, half, true);
                for mirror_reuse in [false, true] {
                    for (aic, velocity) in [(true, true), (true, false), (false, true)] {
                        let mut work = Work::default();
                        let (a, c, v) = assemble(
                            &mesh,
                            0.35,
                            normals,
                            half,
                            aic,
                            velocity,
                            mirror_reuse,
                            &mut work,
                        );
                        if aic {
                            assert_eq!(bits(a.as_slice()), bits(&matrix));
                            assert_eq!(bits(c.as_slice()), bits(&control));
                        } else {
                            assert_eq!((a.len(), c.len()), (0, 0));
                        }
                        if velocity {
                            assert_eq!(bits(v.as_slice()), bits(&bound));
                        } else {
                            assert_eq!(v.len(), 0);
                        }
                        // Mirror reuse halves the full-mode kernel evaluations
                        // and cannot apply to a half-wing lattice.
                        let per_pass = if half {
                            2 * n * n
                        } else if mirror_reuse {
                            n * n / 2
                        } else {
                            n * n
                        };
                        assert_eq!(
                            work.influence_evaluations,
                            (usize::from(aic) + usize::from(velocity)) * per_pass
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn source_panels_hoist_the_per_pair_core_rule() {
        let mesh = lattice(false, 4);
        let sources = source_panels(&mesh);
        for (j, (source, vortex)) in sources.iter().zip(mesh.vortices()).enumerate() {
            assert_eq!(
                source.endpoints(0),
                (vortex.bound_vortex_left, vortex.bound_vortex_right)
            );
            assert_eq!(
                source.endpoints(1),
                (
                    reflect(vortex.bound_vortex_right),
                    reflect(vortex.bound_vortex_left)
                )
            );
            for target in mesh.vortices() {
                let width = (vortex.bound_vortex_right.y - vortex.bound_vortex_left.y)
                    .hypot(vortex.bound_vortex_right.z - vortex.bound_vortex_left.z);
                let expected = if target.component_index == vortex.component_index {
                    0.0001 * width
                } else {
                    (0.0001 * vortex.local_chord).max(0.0001 * width)
                };
                assert_eq!(
                    source.core(target.component_index).to_bits(),
                    expected.to_bits(),
                    "{j}"
                );
            }
        }
    }
}
