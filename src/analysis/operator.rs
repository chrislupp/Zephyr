//! Operators evaluated at stored targets, with reflected source images on demand.
use super::Work;
use crate::core::{MatX, Scalar, Vec3};
use crate::solver::{VortexLattice, horseshoe_vortex_velocity};

pub(crate) fn reflect(v: Vec3) -> Vec3 {
    Vec3::new(v.x, -v.y, v.z)
}

pub(crate) fn core(mesh: &VortexLattice, row: usize, column: usize) -> Scalar {
    let source = &mesh.vortices()[column];
    let width = (source.bound_vortex_right - source.bound_vortex_left).norm();
    if mesh.vortices()[row].component_index == source.component_index {
        0.0001 * width
    } else {
        (0.0001 * source.local_chord).max(0.0001 * width)
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn assemble(
    mesh: &VortexLattice,
    mach: Scalar,
    normals: &[Vec3],
    half: bool,
    aic: bool,
    velocity: bool,
    work: &mut Work,
) -> (MatX, MatX) {
    let n = mesh.vortex_count();
    let pg = (1.0 - mach * mach).sqrt();
    let mut matrix = if aic {
        MatX::zeros(n, n)
    } else {
        MatX::zeros(0, 0)
    };
    let mut velocities = if velocity {
        MatX::zeros(3 * n, n)
    } else {
        MatX::zeros(0, 0)
    };
    let fill = |j: usize, column: &mut [Scalar], velocity_column: &mut [Scalar]| {
        let source = &mesh.vortices()[j];
        for (i, target) in mesh.vortices().iter().enumerate() {
            let radius = core(mesh, i, j);
            let mut v = horseshoe_vortex_velocity(
                &target.control_point,
                &source.bound_vortex_left,
                &source.bound_vortex_right,
                pg,
                radius,
                i != j,
            );
            if half {
                // This is an explicit duplicate's influence: its bound leg is
                // included even for a coincident image midpoint, as in full mode.
                v += horseshoe_vortex_velocity(
                    &target.control_point,
                    &reflect(source.bound_vortex_right),
                    &reflect(source.bound_vortex_left),
                    pg,
                    radius,
                    true,
                );
            }
            if aic {
                column[i] = v.x * normals[i].x + v.y * normals[i].y + v.z * normals[i].z;
            }
            if velocity {
                for k in 0..3 {
                    velocity_column[3 * i + k] = v[k];
                }
            }
        }
    };
    if n > 0 {
        drive_columns(matrix.as_mut_slice(), velocities.as_mut_slice(), n, &fill);
    }
    work.influence_evaluations += n * n * if half { 2 } else { 1 };
    (matrix, velocities)
}

// Match the legacy assembly threshold; small systems stay on one thread.
#[cfg(feature = "parallel")]
pub(super) const PARALLEL_MIN_PANELS: usize = 128;

fn drive_columns_serial<F>(matrix: &mut [Scalar], velocity: &mut [Scalar], n: usize, fill: &F)
where
    F: Fn(usize, &mut [Scalar], &mut [Scalar]) + Sync,
{
    for j in 0..n {
        let column = if matrix.is_empty() {
            &mut []
        } else {
            &mut matrix[j * n..(j + 1) * n]
        };
        let velocities = if velocity.is_empty() {
            &mut []
        } else {
            &mut velocity[j * 3 * n..(j + 1) * 3 * n]
        };
        fill(j, column, velocities);
    }
}

fn drive_columns<F>(matrix: &mut [Scalar], velocity: &mut [Scalar], n: usize, fill: &F)
where
    F: Fn(usize, &mut [Scalar], &mut [Scalar]) + Sync,
{
    #[cfg(feature = "parallel")]
    if n >= PARALLEL_MIN_PANELS {
        use rayon::prelude::*;
        match (matrix.is_empty(), velocity.is_empty()) {
            (false, false) => matrix
                .par_chunks_mut(n)
                .zip(velocity.par_chunks_mut(3 * n))
                .enumerate()
                .for_each(|(j, (a, v))| fill(j, a, v)),
            (false, true) => matrix
                .par_chunks_mut(n)
                .enumerate()
                .for_each(|(j, a)| fill(j, a, &mut [])),
            (true, false) => velocity
                .par_chunks_mut(3 * n)
                .enumerate()
                .for_each(|(j, v)| fill(j, &mut [], v)),
            (true, true) => (0..n)
                .into_par_iter()
                .for_each(|j| fill(j, &mut [], &mut [])),
        }
        return;
    }
    drive_columns_serial(matrix, velocity, n, fill);
}
