//! Analytical circulation derivatives of Trefftz-plane outputs.

use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::core::{EPSILON, PI, Scalar, Vec3, VecX};
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
        if dy.hypot(dz) <= EPSILON {
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

#[derive(Clone, Copy)]
struct Directional {
    value: Scalar,
    derivative: Scalar,
}

impl Directional {
    const fn new(value: Scalar, derivative: Scalar) -> Self {
        Self { value, derivative }
    }

    fn sqrt(self) -> Self {
        let value = self.value.sqrt();
        Self::new(value, 0.5 * self.derivative / value)
    }

    fn powi4(self) -> Self {
        Self::new(
            self.value.powi(4),
            4.0 * self.value.powi(3) * self.derivative,
        )
    }
}

impl Add for Directional {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.value + rhs.value, self.derivative + rhs.derivative)
    }
}

impl Sub for Directional {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.value - rhs.value, self.derivative - rhs.derivative)
    }
}

impl Mul for Directional {
    type Output = Self;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(
            self.value * rhs.value,
            self.derivative * rhs.value + self.value * rhs.derivative,
        )
    }
}

impl Div for Directional {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        Self::new(
            self.value / rhs.value,
            (self.derivative * rhs.value - self.value * rhs.derivative) / rhs.value.powi(2),
        )
    }
}

impl Neg for Directional {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(-self.value, -self.derivative)
    }
}

/// Fixed-circulation directional derivative of far-field lift, side force,
/// and induced drag with respect to lattice geometry.
#[must_use]
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
) -> [Scalar; 3] {
    if !solution.is_valid || mesh.strip_count() == 0 || reference_area <= 0.0 {
        return [0.0; 3];
    }
    let strip_count = mesh.strip_count();
    let mut gamma = vec![0.0; strip_count];
    let mut left = vec![[Directional::new(0.0, 0.0); 2]; strip_count];
    let mut right = left.clone();
    let mut center = left.clone();
    for (index, strip) in mesh.strips().iter().enumerate() {
        gamma[index] = solution
            .circulation
            .rows(strip.first_vortex_index, strip.num_vortices)
            .sum();
        let vortex_index = strip.first_vortex_index + strip.num_vortices - 1;
        let vortex = &mesh.vortices()[vortex_index];
        let left_direction = bound_left_direction
            .get(vortex_index)
            .copied()
            .unwrap_or_else(Vec3::zeros);
        let right_direction = bound_right_direction
            .get(vortex_index)
            .copied()
            .unwrap_or_else(Vec3::zeros);
        left[index] = [
            Directional::new(vortex.bound_vortex_left.y, left_direction.y),
            Directional::new(vortex.bound_vortex_left.z, left_direction.z),
        ];
        right[index] = [
            Directional::new(vortex.bound_vortex_right.y, right_direction.y),
            Directional::new(vortex.bound_vortex_right.z, right_direction.z),
        ];
        center[index] = [
            Directional::new(
                0.5 * (left[index][0].value + right[index][0].value),
                0.5 * (left[index][0].derivative + right[index][0].derivative),
            ),
            Directional::new(
                0.5 * (left[index][1].value + right[index][1].value),
                0.5 * (left[index][1].derivative + right[index][1].derivative),
            ),
        ];
    }

    let mut lift = Directional::new(0.0, 0.0);
    let mut side = Directional::new(0.0, 0.0);
    let mut drag = Directional::new(0.0, 0.0);
    let scale = Directional::new(reference_area, 0.0);
    let two = Directional::new(2.0, 0.0);
    let inv_2_pi = Directional::new(INV_2_PI, 0.0);
    for target in 0..strip_count {
        let target_strip = &mesh.strips()[target];
        let dy = right[target][0] - left[target][0];
        let dz = right[target][1] - left[target][1];
        if dy.value.hypot(dz.value) < EPSILON {
            continue;
        }
        let mut velocity_y = Directional::new(0.0, 0.0);
        let mut velocity_z = Directional::new(0.0, 0.0);
        for source in 0..strip_count {
            let source_strip = &mesh.strips()[source];
            let source_dy = right[source][0] - left[source][0];
            let source_dz = right[source][1] - left[source][1];
            let source_width = (source_dy * source_dy + source_dz * source_dz).sqrt();
            let core = if target_strip.surface_index == source_strip.surface_index {
                Directional::new(0.0, 0.0)
            } else {
                let chord = Directional::new(
                    source_strip.chord_center,
                    strip_chord_direction.get(source).copied().unwrap_or(0.0),
                ) * Directional::new(core_chord_fraction, 0.0);
                let width = source_width * Directional::new(core_width_fraction, 0.0);
                if chord.value >= width.value {
                    chord
                } else {
                    width
                }
            };
            let dy1 = center[target][0] - left[source][0];
            let dy2 = center[target][0] - right[source][0];
            let dz1 = center[target][1] - left[source][1];
            let dz2 = center[target][1] - right[source][1];
            let radius1 = ((dy1 * dy1 + dz1 * dz1) * (dy1 * dy1 + dz1 * dz1) + core.powi4()).sqrt();
            let radius2 = ((dy2 * dy2 + dz2 * dz2) * (dy2 * dy2 + dz2 * dz2) + core.powi4()).sqrt();
            if radius1.value > EPSILON && radius2.value > EPSILON {
                let source_gamma = Directional::new(gamma[source], 0.0);
                velocity_y = velocity_y + inv_2_pi * source_gamma * (dz1 / radius1 - dz2 / radius2);
                velocity_z =
                    velocity_z + inv_2_pi * source_gamma * (-dy1 / radius1 + dy2 / radius2);
            }
        }
        let target_gamma = Directional::new(gamma[target], 0.0);
        lift = lift + two * target_gamma * dy / scale;
        side = side - two * target_gamma * dz / scale;
        drag = drag + target_gamma * (dz * velocity_y - dy * velocity_z) / scale;
    }
    [lift.derivative, side.derivative, drag.derivative]
}
