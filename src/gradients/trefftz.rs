//! Analytical circulation derivatives of Trefftz-plane outputs.

use crate::core::{EPSILON, PI, Scalar, Vec3, VecX};
use crate::error::{PegasusError, Result};
use crate::solver::{SolutionState, VortexLattice};

const INV_2_PI: Scalar = 1.0 / (2.0 * PI);

/// Circulation derivatives of the far-field lift, side force, and induced
/// drag coefficients.
#[derive(Clone, Debug, Default)]
pub struct TrefftzCirculationGradients {
    /// `dCL_far/dGamma`.
    pub lift: VecX,
    /// `dCY_far/dGamma`.
    pub side: VecX,
    /// `dCDi/dGamma`.
    pub induced_drag: VecX,
}

/// Differentiate the same discrete Trefftz calculation used by
/// `ForceCalculator`.
#[must_use]
pub fn trefftz_circulation_gradients(
    mesh: &VortexLattice,
    solution: &SolutionState,
    reference_area: Scalar,
    core_chord_fraction: Scalar,
    core_width_fraction: Scalar,
) -> TrefftzCirculationGradients {
    let vortex_count = mesh.vortex_count();
    let strip_count = mesh.strip_count();
    let mut result = TrefftzCirculationGradients {
        lift: VecX::zeros(vortex_count),
        side: VecX::zeros(vortex_count),
        induced_drag: VecX::zeros(vortex_count),
    };
    if !solution.is_valid || strip_count == 0 || reference_area <= 0.0 {
        return result;
    }

    let mut strip_gamma = vec![0.0; strip_count];
    let mut left_yz = vec![(0.0, 0.0); strip_count];
    let mut right_yz = vec![(0.0, 0.0); strip_count];
    let mut center_yz = vec![(0.0, 0.0); strip_count];
    for (index, strip) in mesh.strips().iter().enumerate() {
        strip_gamma[index] = solution
            .circulation
            .rows(strip.first_vortex_index, strip.num_vortices)
            .sum();
        let vortex = &mesh.vortices()[strip.first_vortex_index + strip.num_vortices - 1];
        left_yz[index] = (vortex.bound_vortex_left.y, vortex.bound_vortex_left.z);
        right_yz[index] = (vortex.bound_vortex_right.y, vortex.bound_vortex_right.z);
        center_yz[index] = (
            0.5 * (left_yz[index].0 + right_yz[index].0),
            0.5 * (left_yz[index].1 + right_yz[index].1),
        );
    }

    // `kernel[target, source]` maps source strip circulation to the target
    // strip's unnormalized induced-drag velocity contraction.
    let mut kernel = vec![0.0; strip_count * strip_count];
    for target in 0..strip_count {
        let target_strip = &mesh.strips()[target];
        let dy = right_yz[target].0 - left_yz[target].0;
        let dz = right_yz[target].1 - left_yz[target].1;
        // `<` rather than `<=`, matching the primal far-field loop. The two
        // disagreed only at exactly `EPSILON`, where the primal summed a strip
        // whose circulation derivative this loop left at zero.
        if dy.hypot(dz) < EPSILON {
            continue;
        }
        for source in 0..strip_count {
            let source_strip = &mesh.strips()[source];
            let source_width = (right_yz[source].0 - left_yz[source].0)
                .hypot(right_yz[source].1 - left_yz[source].1);
            let core = if target_strip.surface_index == source_strip.surface_index {
                0.0
            } else {
                (core_chord_fraction * source_strip.chord_center)
                    .max(core_width_fraction * source_width)
            };
            let core_4 = core.powi(4);
            let dy1 = center_yz[target].0 - left_yz[source].0;
            let dy2 = center_yz[target].0 - right_yz[source].0;
            let dz1 = center_yz[target].1 - left_yz[source].1;
            let dz2 = center_yz[target].1 - right_yz[source].1;
            let radius1 = ((dy1 * dy1 + dz1 * dz1).powi(2) + core_4).sqrt();
            let radius2 = ((dy2 * dy2 + dz2 * dz2).powi(2) + core_4).sqrt();
            if radius1 <= EPSILON || radius2 <= EPSILON {
                continue;
            }
            let velocity_y = INV_2_PI * (dz1 / radius1 - dz2 / radius2);
            let velocity_z = INV_2_PI * (-dy1 / radius1 + dy2 / radius2);
            kernel[target * strip_count + source] = dz * velocity_y - dy * velocity_z;
        }
    }

    let mut strip_lift = vec![0.0; strip_count];
    let mut strip_side = vec![0.0; strip_count];
    let mut strip_drag = vec![0.0; strip_count];
    for target in 0..strip_count {
        let dy = right_yz[target].0 - left_yz[target].0;
        let dz = right_yz[target].1 - left_yz[target].1;
        strip_lift[target] = 2.0 * dy / reference_area;
        strip_side[target] = -2.0 * dz / reference_area;
        let direct: Scalar = (0..strip_count)
            .map(|source| kernel[target * strip_count + source] * strip_gamma[source])
            .sum();
        let transpose: Scalar = (0..strip_count)
            .map(|source| kernel[source * strip_count + target] * strip_gamma[source])
            .sum();
        strip_drag[target] = (direct + transpose) / reference_area;
    }

    for (strip_index, strip) in mesh.strips().iter().enumerate() {
        for vortex_index in strip.first_vortex_index..strip.first_vortex_index + strip.num_vortices
        {
            result.lift[vortex_index] = strip_lift[strip_index];
            result.side[vortex_index] = strip_side[strip_index];
            result.induced_drag[vortex_index] = strip_drag[strip_index];
        }
    }
    result
}

/// Fixed-circulation directional derivative of far-field lift, side force,
/// and induced drag with respect to lattice geometry.
#[allow(clippy::too_many_arguments)]
pub fn trefftz_geometry_jvp(
    mesh: &VortexLattice,
    solution: &SolutionState,
    reference_area: Scalar,
    core_chord_fraction: Scalar,
    core_width_fraction: Scalar,
    bound_left_direction: &[Vec3],
    bound_right_direction: &[Vec3],
    strip_chord_direction: &[Scalar],
) -> Result<[Scalar; 3]> {
    // Zero-filling a short slice returns a partial derivative that looks
    // exactly like a correct one. `strip_chord_direction` in particular was
    // never checked at all.
    if bound_left_direction.len() != mesh.vortex_count()
        || bound_right_direction.len() != mesh.vortex_count()
        || strip_chord_direction.len() != mesh.strip_count()
    {
        return Err(PegasusError::InvalidGeometry(
            "geometry direction slices must match the lattice dimensions".into(),
        ));
    }
    if !solution.is_valid || mesh.strip_count() == 0 || reference_area <= 0.0 {
        return Ok([0.0; 3]);
    }
    let strip_count = mesh.strip_count();
    let mut gamma = vec![0.0; strip_count];
    // Trailing-edge station coordinates and their directional rates. `d_*` is
    // the derivative of the preceding quantity along the perturbation
    // direction; the circulation is held fixed, so its rate is identically zero
    // and never appears below.
    let mut y_left = vec![0.0; strip_count];
    let mut z_left = vec![0.0; strip_count];
    let mut y_right = vec![0.0; strip_count];
    let mut z_right = vec![0.0; strip_count];
    let mut d_y_left = vec![0.0; strip_count];
    let mut d_z_left = vec![0.0; strip_count];
    let mut d_y_right = vec![0.0; strip_count];
    let mut d_z_right = vec![0.0; strip_count];
    for (index, strip) in mesh.strips().iter().enumerate() {
        gamma[index] = solution
            .circulation
            .rows(strip.first_vortex_index, strip.num_vortices)
            .sum();
        let vortex_index = strip.first_vortex_index + strip.num_vortices - 1;
        let vortex = &mesh.vortices()[vortex_index];
        y_left[index] = vortex.bound_vortex_left.y;
        z_left[index] = vortex.bound_vortex_left.z;
        y_right[index] = vortex.bound_vortex_right.y;
        z_right[index] = vortex.bound_vortex_right.z;
        d_y_left[index] = bound_left_direction[vortex_index].y;
        d_z_left[index] = bound_left_direction[vortex_index].z;
        d_y_right[index] = bound_right_direction[vortex_index].y;
        d_z_right[index] = bound_right_direction[vortex_index].z;
    }
    let y_center = |i: usize| 0.5 * (y_left[i] + y_right[i]);
    let z_center = |i: usize| 0.5 * (z_left[i] + z_right[i]);
    let d_y_center = |i: usize| 0.5 * (d_y_left[i] + d_y_right[i]);
    let d_z_center = |i: usize| 0.5 * (d_z_left[i] + d_z_right[i]);

    let mut d_lift = 0.0;
    let mut d_side = 0.0;
    let mut d_drag = 0.0;
    for target in 0..strip_count {
        let target_strip = &mesh.strips()[target];
        let dy = y_right[target] - y_left[target];
        let dz = z_right[target] - z_left[target];
        // Same test as the primal far-field loop in `ForceCalculator`, so the
        // derivative includes exactly the strips the primal sums over.
        if dy.hypot(dz) < EPSILON {
            continue;
        }
        let d_dy = d_y_right[target] - d_y_left[target];
        let d_dz = d_z_right[target] - d_z_left[target];

        let (mut velocity_y, mut velocity_z) = (0.0, 0.0);
        let (mut d_velocity_y, mut d_velocity_z) = (0.0, 0.0);
        for source in 0..strip_count {
            let source_strip = &mesh.strips()[source];
            let source_dy = y_right[source] - y_left[source];
            let source_dz = z_right[source] - z_left[source];
            let d_source_dy = d_y_right[source] - d_y_left[source];
            let d_source_dz = d_z_right[source] - d_z_left[source];
            let source_width = source_dy.hypot(source_dz);
            // d|w|/dt = (w . dw/dt) / |w|, undefined only for a zero-width
            // strip, which contributes no core anyway.
            let d_source_width = if source_width > 0.0 {
                (source_dy * d_source_dy + source_dz * d_source_dz) / source_width
            } else {
                0.0
            };

            // `max` is differentiable away from the tie; the branch is taken on
            // the values so the rate follows whichever term is active, matching
            // the primal's selection exactly.
            let (core, d_core) = if target_strip.surface_index == source_strip.surface_index {
                (0.0, 0.0)
            } else {
                let chord = core_chord_fraction * source_strip.chord_center;
                let d_chord = core_chord_fraction * strip_chord_direction[source];
                let width = core_width_fraction * source_width;
                let d_width = core_width_fraction * d_source_width;
                if chord >= width {
                    (chord, d_chord)
                } else {
                    (width, d_width)
                }
            };
            let core_4 = core.powi(4);
            let d_core_4 = 4.0 * core.powi(3) * d_core;

            let dy1 = y_center(target) - y_left[source];
            let dy2 = y_center(target) - y_right[source];
            let dz1 = z_center(target) - z_left[source];
            let dz2 = z_center(target) - z_right[source];
            let d_dy1 = d_y_center(target) - d_y_left[source];
            let d_dy2 = d_y_center(target) - d_y_right[source];
            let d_dz1 = d_z_center(target) - d_z_left[source];
            let d_dz2 = d_z_center(target) - d_z_right[source];

            // r = sqrt(q^2 + c^4) with q = dy^2 + dz^2, so
            // dr/dt = (q dq/dt + c^3 dc/dt * 2) / r, written out through the
            // intermediates to keep each step readable.
            let q1 = dy1 * dy1 + dz1 * dz1;
            let q2 = dy2 * dy2 + dz2 * dz2;
            let d_q1 = 2.0 * (dy1 * d_dy1 + dz1 * d_dz1);
            let d_q2 = 2.0 * (dy2 * d_dy2 + dz2 * d_dz2);
            let radius1 = (q1 * q1 + core_4).sqrt();
            let radius2 = (q2 * q2 + core_4).sqrt();
            if radius1 > EPSILON && radius2 > EPSILON {
                let d_radius1 = 0.5 * (2.0 * q1 * d_q1 + d_core_4) / radius1;
                let d_radius2 = 0.5 * (2.0 * q2 * d_q2 + d_core_4) / radius2;
                let source_gamma = gamma[source];
                velocity_y += INV_2_PI * source_gamma * (dz1 / radius1 - dz2 / radius2);
                velocity_z += INV_2_PI * source_gamma * (-dy1 / radius1 + dy2 / radius2);
                // Quotient rule on each of the four dz/r and dy/r terms.
                d_velocity_y += INV_2_PI
                    * source_gamma
                    * ((d_dz1 * radius1 - dz1 * d_radius1) / (radius1 * radius1)
                        - (d_dz2 * radius2 - dz2 * d_radius2) / (radius2 * radius2));
                d_velocity_z += INV_2_PI
                    * source_gamma
                    * (-(d_dy1 * radius1 - dy1 * d_radius1) / (radius1 * radius1)
                        + (d_dy2 * radius2 - dy2 * d_radius2) / (radius2 * radius2));
            }
        }
        // The reference area is not a design variable here, so it scales the
        // rate the same way it scales the value.
        let target_gamma = gamma[target];
        d_lift += 2.0 * target_gamma * d_dy / reference_area;
        d_side -= 2.0 * target_gamma * d_dz / reference_area;
        d_drag += target_gamma
            * (d_dz * velocity_y + dz * d_velocity_y - d_dy * velocity_z - dy * d_velocity_z)
            / reference_area;
    }
    Ok([d_lift, d_side, d_drag])
}

/// Cotangents of the fixed-circulation Trefftz geometry. Endpoint arrays follow
/// vortex order; only each strip's last chordwise vortex receives a seed.
#[derive(Clone, Debug)]
pub struct TrefftzGeometryCotangent {
    pub bound_left: Vec<Vec3>,
    pub bound_right: Vec<Vec3>,
    pub strip_chord: Vec<Scalar>,
}

/// Hand-derived transpose of [`trefftz_geometry_jvp`]. Seeds are ordered
/// `[CL_far, CY_far, CDi]`; reference area and core fractions are constants.
/// The active core-radius branch matches the primal (`chord >= width` at a
/// tie). At that nonsmooth tie this is the chosen branch derivative only.
/// Degenerate-strip and coincident-point guards match the primal exactly.
#[allow(clippy::too_many_arguments)]
pub fn trefftz_geometry_vjp(
    mesh: &VortexLattice,
    solution: &SolutionState,
    reference_area: Scalar,
    core_chord_fraction: Scalar,
    core_width_fraction: Scalar,
    seed: [Scalar; 3],
) -> Result<TrefftzGeometryCotangent> {
    if seed.iter().any(|x| !x.is_finite()) {
        return Err(PegasusError::InvalidGeometry(
            "nonfinite Trefftz seed".into(),
        ));
    }
    let mut out = TrefftzGeometryCotangent {
        bound_left: vec![Vec3::zeros(); mesh.vortex_count()],
        bound_right: vec![Vec3::zeros(); mesh.vortex_count()],
        strip_chord: vec![0.0; mesh.strip_count()],
    };
    if !solution.is_valid || mesh.strip_count() == 0 || reference_area <= 0.0 {
        return Ok(out);
    }
    let index: Vec<_> = mesh
        .strips()
        .iter()
        .map(|s| s.first_vortex_index + s.num_vortices - 1)
        .collect();
    let gamma: Vec<_> = mesh
        .strips()
        .iter()
        .map(|s| {
            solution
                .circulation
                .rows(s.first_vortex_index, s.num_vortices)
                .sum()
        })
        .collect();
    let yz = |p: Vec3| Vec3::new(0.0, p.y, p.z);
    for (i, target) in mesh.strips().iter().enumerate() {
        let ti = index[i];
        let left = yz(mesh.vortices()[ti].bound_vortex_left);
        let right = yz(mesh.vortices()[ti].bound_vortex_right);
        let width = right - left;
        if width.norm() < EPSILON {
            continue;
        }
        let center = 0.5 * (left + right);
        let direct = Vec3::new(0.0, 2.0 * seed[0], -2.0 * seed[1]) * (gamma[i] / reference_area);
        out.bound_right[ti] += direct;
        out.bound_left[ti] -= direct;
        if seed[2] == 0.0 {
            continue;
        }
        for (j, source) in mesh.strips().iter().enumerate() {
            let sj = index[j];
            let sl = yz(mesh.vortices()[sj].bound_vortex_left);
            let sr = yz(mesh.vortices()[sj].bound_vortex_right);
            let source_width = sr - sl;
            let chord_core = core_chord_fraction * source.chord_center;
            let width_core = core_width_fraction * source_width.norm();
            let core = if target.surface_index == source.surface_index {
                0.0
            } else {
                chord_core.max(width_core)
            };
            let d1 = center - sl;
            let d2 = center - sr;
            let radius = |d: Vec3| (d.norm_squared().powi(2) + core.powi(4)).sqrt();
            let r1 = radius(d1);
            let r2 = radius(d2);
            if r1 <= EPSILON || r2 <= EPSILON {
                continue;
            }
            let factor = seed[2] * gamma[i] * gamma[j] * INV_2_PI / reference_area;
            let mut core_bar = 0.0;
            for (leg, d, r, sign) in [(0, d1, r1, 1.0), (1, d2, r2, -1.0)] {
                // phi = factor * sign * dot(width, d) / sqrt(dot(d,d)^2 + core^4).
                let weight = factor * sign;
                let numerator = width.dot(&d);
                let width_bar = weight * d / r;
                let d_bar =
                    weight * (width / r - 2.0 * numerator * d.norm_squared() * d / r.powi(3));
                out.bound_right[ti] += width_bar + 0.5 * d_bar;
                out.bound_left[ti] += -width_bar + 0.5 * d_bar;
                if leg == 0 {
                    out.bound_left[sj] -= d_bar;
                } else {
                    out.bound_right[sj] -= d_bar;
                }
                core_bar -= 2.0 * weight * numerator * core.powi(3) / r.powi(3);
            }
            if target.surface_index != source.surface_index {
                if chord_core >= width_core {
                    out.strip_chord[j] += core_chord_fraction * core_bar;
                } else if source_width.norm() > 0.0 {
                    let bar = core_width_fraction * core_bar * source_width / source_width.norm();
                    out.bound_right[sj] += bar;
                    out.bound_left[sj] -= bar;
                }
            }
        }
    }
    Ok(out)
}
