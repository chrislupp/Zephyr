//! Near-field Kutta–Joukowski and Trefftz-plane loads.

use crate::core::{EPSILON, PI, Scalar, Vec3};
use crate::geometry::{Aircraft, ReferenceQuantities};
use crate::solver::{FlightCondition, SolutionState, VortexLattice, VortexLatticeSolver};

use super::{AerodynamicResults, Coefficients, StripLoads, SurfaceLoads};

const INV_2_PI: Scalar = 1.0 / (2.0 * PI);

/// Aerodynamic force and moment post-processor.
#[derive(Clone, Debug)]
pub struct ForceCalculator {
    include_induced_velocity: bool,
    trefftz_core_chord_fraction: Scalar,
    trefftz_core_width_fraction: Scalar,
}

impl Default for ForceCalculator {
    fn default() -> Self {
        Self {
            include_induced_velocity: true,
            trefftz_core_chord_fraction: 0.01,
            trefftz_core_width_fraction: 0.1,
        }
    }
}

impl ForceCalculator {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_include_induced_velocity(&mut self, include: bool) {
        self.include_induced_velocity = include;
    }

    pub fn set_trefftz_core_radius(&mut self, chord_fraction: Scalar, width_fraction: Scalar) {
        self.trefftz_core_chord_fraction = chord_fraction;
        self.trefftz_core_width_fraction = width_fraction;
    }

    #[must_use]
    pub fn compute_for_aircraft(
        &self,
        aircraft: &Aircraft,
        mesh: &VortexLattice,
        solution: &SolutionState,
        condition: &FlightCondition,
    ) -> AerodynamicResults {
        let mut results = self.compute_with_references(
            &ReferenceQuantities::from(aircraft),
            mesh,
            solution,
            condition,
        );
        if results.is_valid {
            let profile_drag = aircraft.profile_drag_coefficient();
            let (sin, cos) = condition.alpha.sin_cos();
            results.viscous_drag += profile_drag;
            results.coefficients.drag += profile_drag;
            results.body_force += Vec3::new(profile_drag * cos, 0.0, profile_drag * sin);
        }
        results
    }

    #[must_use]
    pub fn compute_with_references(
        &self,
        reference: &ReferenceQuantities,
        mesh: &VortexLattice,
        solution: &SolutionState,
        condition: &FlightCondition,
    ) -> AerodynamicResults {
        if !solution.is_valid || mesh.vortex_count() == 0 {
            return AerodynamicResults::default();
        }
        let mut results = AerodynamicResults {
            reference_area: reference.area,
            reference_chord: reference.chord,
            reference_span: reference.span,
            moment_reference_point: reference.moment_reference_point,
            ..AerodynamicResults::default()
        };
        let mut total_force = Vec3::zeros();
        let mut total_moment = Vec3::zeros();
        let mut viscous_drag_area = 0.0;
        for (surface_index, info) in mesh.surface_info().iter().enumerate() {
            let mut surface = SurfaceLoads {
                surface_index,
                name: info.name.clone(),
                ..SurfaceLoads::default()
            };
            let mut surface_force = Vec3::zeros();
            let mut surface_moment = Vec3::zeros();
            for (strip_index, strip) in mesh.strips().iter().enumerate() {
                if strip.surface_index != surface_index {
                    continue;
                }
                let loads = self.strip_loads(mesh, solution, condition, strip_index);
                surface.area += loads.area;
                viscous_drag_area += loads.viscous_drag * loads.area;
                let force = loads.body_force * loads.area;
                surface_force += force;
                let arm = loads.quarter_chord_position - results.moment_reference_point;
                surface_moment += Vec3::new(
                    -loads.body_moment.x,
                    loads.body_moment.y,
                    -loads.body_moment.z,
                ) * loads.chord
                    * loads.area;
                surface_moment += arm.cross(&force);
                surface.strip_loads.push(loads);
            }
            surface.body_force = surface_force / results.reference_area;
            surface.body_moment = Vec3::new(-surface_moment.x, surface_moment.y, -surface_moment.z)
                / (results.reference_area * results.reference_chord);
            let (lift, drag) = body_to_stability(&surface.body_force, condition.alpha);
            surface.coefficients = Coefficients {
                lift,
                drag,
                side: surface.body_force.y,
                roll: surface.body_moment.x,
                pitch: surface.body_moment.y,
                yaw: surface.body_moment.z,
            };
            total_force += surface_force;
            total_moment += surface_moment;
            results.surface_loads.push(surface);
        }
        results.body_force = total_force / results.reference_area;
        let mathematical_moment = total_moment / (results.reference_area * results.reference_chord);
        results.body_moment = Vec3::new(
            -mathematical_moment.x,
            mathematical_moment.y,
            -mathematical_moment.z,
        );
        let (lift, drag) = body_to_stability(&results.body_force, condition.alpha);
        results.viscous_drag = viscous_drag_area / results.reference_area;
        results.coefficients = Coefficients {
            lift,
            drag,
            side: results.body_force.y,
            roll: results.body_moment.x,
            pitch: results.body_moment.y,
            yaw: results.body_moment.z,
        };
        results.near_field_induced_drag = drag - results.viscous_drag;
        (
            results.far_field_lift,
            results.far_field_side,
            results.far_field_induced_drag,
            results.span_efficiency,
        ) = self.trefftz_plane(
            mesh,
            solution,
            results.reference_area,
            results.reference_span,
        );
        results.is_valid = true;
        results
    }

    fn strip_loads(
        &self,
        mesh: &VortexLattice,
        solution: &SolutionState,
        condition: &FlightCondition,
        strip_index: usize,
    ) -> StripLoads {
        let strip = &mesh.strips()[strip_index];
        let area = strip.chord_center * strip.width;
        let mut loads = StripLoads {
            strip_index,
            chord: strip.chord_center,
            width: strip.width,
            area,
            leading_edge_position: strip.leading_edge_center,
            quarter_chord_position: strip.leading_edge_center
                + Vec3::new(0.25 * strip.chord_center, 0.0, 0.0),
            ..StripLoads::default()
        };
        let freestream = VortexLatticeSolver::freestream_direction(condition);
        let mut force = Vec3::zeros();
        let mut moment = Vec3::zeros();
        let mut normal_loading = 0.0;
        for offset in 0..strip.num_vortices {
            let index = strip.first_vortex_index + offset;
            let vortex = &mesh.vortices()[index];
            let mut velocity = freestream;
            if self.include_induced_velocity
                && solution.induced_velocity.nrows() > index
                && solution.induced_velocity.ncols() == 3
            {
                velocity += Vec3::new(
                    solution.induced_velocity[(index, 0)],
                    solution.induced_velocity[(index, 1)],
                    solution.induced_velocity[(index, 2)],
                );
            }
            let line = vortex.bound_vortex_right - vortex.bound_vortex_left;
            let panel_force = 2.0 * solution.circulation[index] * velocity.cross(&line);
            force += panel_force;
            moment += (vortex.control_point - loads.quarter_chord_position).cross(&panel_force);
            normal_loading += strip.chord_center
                * panel_force.dot(&Vec3::new(
                    0.0,
                    strip.trefftz_normal_y,
                    strip.trefftz_normal_z,
                ));
        }
        if area > EPSILON {
            loads.body_force = force / area;
            loads.body_moment = Vec3::new(-moment.x, moment.y, -moment.z) / (loads.chord * area);
            loads.normal_span_loading = normal_loading / area;
        }
        let (lift, drag) = body_to_stability(&loads.body_force, condition.alpha);
        loads.coefficients = Coefficients {
            lift,
            drag,
            side: loads.body_force.y,
            roll: loads.body_moment.x,
            pitch: loads.body_moment.y,
            yaw: loads.body_moment.z,
        };
        if strip.has_viscous_data {
            loads.viscous_drag =
                interpolate_drag_polar(&strip.viscous_polar, loads.coefficients.lift).max(0.0);
            let (sin, cos) = condition.alpha.sin_cos();
            loads.body_force += Vec3::new(loads.viscous_drag * cos, 0.0, loads.viscous_drag * sin);
            loads.coefficients.drag += loads.viscous_drag;
        }
        loads
    }

    /// Far-field coefficients `(lift, side, induced_drag, efficiency)`.
    #[must_use]
    pub fn trefftz_plane(
        &self,
        mesh: &VortexLattice,
        solution: &SolutionState,
        reference_area: Scalar,
        reference_span: Scalar,
    ) -> (Scalar, Scalar, Scalar, Scalar) {
        if !solution.is_valid || mesh.strip_count() == 0 || reference_area <= 0.0 {
            return (0.0, 0.0, 0.0, 0.0);
        }
        let count = mesh.strip_count();
        let mut gamma = vec![0.0; count];
        let mut left = vec![Vec3::zeros(); count];
        let mut right = vec![Vec3::zeros(); count];
        let mut center = vec![Vec3::zeros(); count];
        for (index, strip) in mesh.strips().iter().enumerate() {
            gamma[index] = solution
                .circulation
                .rows(strip.first_vortex_index, strip.num_vortices)
                .sum();
            let vortex = &mesh.vortices()[strip.first_vortex_index + strip.num_vortices - 1];
            left[index] = vortex.bound_vortex_left;
            right[index] = vortex.bound_vortex_right;
            center[index] = 0.5 * (left[index] + right[index]);
        }
        let (mut lift, mut side, mut drag) = (0.0, 0.0, 0.0);
        for target in 0..count {
            let target_strip = &mesh.strips()[target];
            let dy = right[target].y - left[target].y;
            let dz = right[target].z - left[target].z;
            if dy.hypot(dz) < EPSILON {
                continue;
            }
            let (mut velocity_y, mut velocity_z) = (0.0, 0.0);
            for source in 0..count {
                let source_strip = &mesh.strips()[source];
                let source_width =
                    (right[source].y - left[source].y).hypot(right[source].z - left[source].z);
                let core = if target_strip.surface_index == source_strip.surface_index {
                    0.0
                } else {
                    (self.trefftz_core_chord_fraction * source_strip.chord_center)
                        .max(self.trefftz_core_width_fraction * source_width)
                };
                let core_4 = core.powi(4);
                let dy1 = center[target].y - left[source].y;
                let dy2 = center[target].y - right[source].y;
                let dz1 = center[target].z - left[source].z;
                let dz2 = center[target].z - right[source].z;
                let radius1 = ((dy1 * dy1 + dz1 * dz1).powi(2) + core_4).sqrt();
                let radius2 = ((dy2 * dy2 + dz2 * dz2).powi(2) + core_4).sqrt();
                if radius1 > EPSILON && radius2 > EPSILON {
                    velocity_y += INV_2_PI * gamma[source] * (dz1 / radius1 - dz2 / radius2);
                    velocity_z += INV_2_PI * gamma[source] * (-dy1 / radius1 + dy2 / radius2);
                }
            }
            lift += 2.0 * gamma[target] * dy / reference_area;
            side -= 2.0 * gamma[target] * dz / reference_area;
            drag += gamma[target] * (dz * velocity_y - dy * velocity_z) / reference_area;
        }
        let aspect_ratio = reference_span * reference_span / reference_area;
        let efficiency = if drag.abs() > EPSILON {
            (lift * lift + side * side) / (PI * aspect_ratio * drag)
        } else {
            0.0
        };
        (lift, side, drag, efficiency)
    }
}

fn body_to_stability(force: &Vec3, alpha: Scalar) -> (Scalar, Scalar) {
    let (sin, cos) = alpha.sin_cos();
    (force.z * cos - force.x * sin, force.x * cos + force.z * sin)
}

fn interpolate_drag_polar(polar: &[Scalar; 6], lift: Scalar) -> Scalar {
    let points = [
        (polar[0], polar[1]),
        (polar[2], polar[3]),
        (polar[4], polar[5]),
    ];
    let mut drag = 0.0;
    for (index, (x, y)) in points.iter().enumerate() {
        let mut basis = 1.0;
        for (other_index, (other_x, _)) in points.iter().enumerate() {
            if index == other_index {
                continue;
            }
            let denominator = x - other_x;
            if denominator.abs() <= EPSILON {
                return points
                    .iter()
                    .min_by(|left, right| (left.0 - lift).abs().total_cmp(&(right.0 - lift).abs()))
                    .map_or(0.0, |point| point.1);
            }
            basis *= (lift - other_x) / denominator;
        }
        drag += y * basis;
    }
    drag
}
