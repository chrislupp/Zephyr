//! Selective near-field integration and optional full-span presentation.
use super::operator::reflect;
use super::{AnalysisRequest, Output, Work};
use crate::aero::forces::drag_polar_value_and_slope;
use crate::aero::{Coefficients, StripLoads, SurfaceLoads};
use crate::core::{EPSILON, MatX, Scalar, Vec3, VecX};
use crate::error::Result;
use crate::solver::{FlightCondition, VortexLatticeSolver};
use std::collections::BTreeMap;

pub(crate) fn induced(solver: &VortexLatticeSolver, i: usize) -> Vec3 {
    let v = &solver.solution().induced_velocity;
    if v.nrows() > i {
        Vec3::new(v[(i, 0)], v[(i, 1)], v[(i, 2)])
    } else {
        Vec3::zeros()
    }
}

fn coefficients(force: Vec3, moment: Vec3, alpha: Scalar) -> Coefficients {
    let (s, c) = alpha.sin_cos();
    Coefficients {
        lift: force.z * c - force.x * s,
        drag: force.x * c + force.z * s,
        side: force.y,
        roll: -moment.x,
        pitch: moment.y,
        yaw: -moment.z,
    }
}

type LoadEvaluation = (BTreeMap<Output, Scalar>, Option<Vec<SurfaceLoads>>);

pub(crate) fn evaluate(
    solver: &VortexLatticeSolver,
    condition: &FlightCondition,
    half: bool,
    profile: Scalar,
    request: &AnalysisRequest,
    work: &mut Work,
) -> Result<LoadEvaluation> {
    let plan = request.plan();
    let mesh = solver.mesh();
    let reference = solver.reference_quantities();
    let freestream = VortexLatticeSolver::freestream_direction(condition);
    let mut values = BTreeMap::new();
    let mut details = request.details.surface_loads.then(Vec::new);
    if plan.near {
        let mut total_force = Vec3::zeros();
        let mut total_moment = Vec3::zeros();
        let mut viscous_area = 0.0;
        let (sin, cos) = condition.alpha.sin_cos();
        let drag_axis = Vec3::new(cos, 0.0, sin);
        let need_viscous = request.details.surface_loads
            || plan.outputs.iter().any(|o| {
                matches!(
                    o,
                    Output::Drag
                        | Output::ViscousDrag
                        | Output::NearFieldInducedDrag
                        | Output::Roll
                        | Output::Pitch
                        | Output::Yaw
                )
            });
        for info in mesh.surface_info() {
            let mut forces = [Vec3::zeros(); 2];
            let mut moments = [Vec3::zeros(); 2];
            let mut records = [Vec::new(), Vec::new()];
            let mut area_sum = 0.0;
            for strip in
                &mesh.strips()[info.first_strip_index..info.first_strip_index + info.num_strips]
            {
                let area = strip.chord_center * strip.width;
                area_sum += area;
                let qc = strip.leading_edge_center + Vec3::new(0.25 * strip.chord_center, 0.0, 0.0);
                let mut force = Vec3::zeros();
                let mut moment_qc = Vec3::zeros();
                let mut normal_load = 0.0;
                for i in strip.first_vortex_index..strip.first_vortex_index + strip.num_vortices {
                    let v = &mesh.vortices()[i];
                    let f = 2.0
                        * solver.solution().circulation[i]
                        * (freestream + induced(solver, i))
                            .cross(&(v.bound_vortex_right - v.bound_vortex_left));
                    force += f;
                    if plan.moments {
                        moment_qc += (v.bound_vortex_center - qc).cross(&f);
                    }
                    if details.is_some() {
                        normal_load += strip.chord_center
                            * f.dot(&Vec3::new(
                                0.0,
                                strip.trefftz_normal_y,
                                strip.trefftz_normal_z,
                            ));
                    }
                }
                let inv_area = if area > EPSILON { 1.0 / area } else { 0.0 };
                // Match the force calculator's degenerate-strip normalization.
                if area <= EPSILON {
                    force = Vec3::zeros();
                    moment_qc = Vec3::zeros();
                }
                let viscous = if need_viscous && strip.has_viscous_data {
                    let cl = (force.z * cos - force.x * sin) * inv_area;
                    drag_polar_value_and_slope(&strip.viscous_polar, cl)
                        .0
                        .max(0.0)
                } else {
                    0.0
                };
                force += area * viscous * drag_axis;
                viscous_area += area * viscous * if half { 2.0 } else { 1.0 };
                for image in 0..if half { 2 } else { 1 } {
                    let f = if image == 0 { force } else { reflect(force) };
                    let q = if image == 0 { qc } else { reflect(qc) };
                    let m = if image == 0 {
                        moment_qc
                    } else {
                        -reflect(moment_qc)
                    };
                    forces[image] += f;
                    if plan.moments {
                        moments[image] += m + (q - reference.moment_reference_point).cross(&f);
                    }
                    if details.is_some() {
                        let local_force = f * inv_area;
                        let local_math_moment = m * inv_area / strip.chord_center;
                        let mut coeff =
                            coefficients(local_force, local_math_moment, condition.alpha);
                        // Full records include the same local coefficients as the legacy API.
                        coeff.drag = local_force.x * cos + local_force.z * sin;
                        records[image].push(StripLoads {
                            strip_index: 0,
                            chord: strip.chord_center,
                            width: strip.width,
                            area,
                            leading_edge_position: if image == 0 {
                                strip.leading_edge_center
                            } else {
                                reflect(strip.leading_edge_center)
                            },
                            quarter_chord_position: q,
                            coefficients: coeff,
                            body_force: local_force,
                            body_moment: Vec3::new(
                                -local_math_moment.x,
                                local_math_moment.y,
                                -local_math_moment.z,
                            ),
                            normal_span_loading: normal_load * inv_area,
                            viscous_drag: viscous,
                            ..StripLoads::default()
                        });
                    }
                }
            }
            for image in 0..if half { 2 } else { 1 } {
                total_force += forces[image];
                total_moment += moments[image];
                if let Some(details) = &mut details {
                    let surface_index = details.len();
                    let first_strip: usize = details
                        .iter()
                        .map(|s: &SurfaceLoads| s.strip_loads.len())
                        .sum();
                    let mut strips = std::mem::take(&mut records[image]);
                    for (offset, strip) in strips.iter_mut().enumerate() {
                        strip.strip_index = first_strip + offset;
                    }
                    let f = forces[image] / reference.area;
                    let m = moments[image] / (reference.area * reference.chord);
                    details.push(SurfaceLoads {
                        surface_index,
                        name: if image == 0 {
                            info.name.clone()
                        } else {
                            format!("{} (yduplicate)", info.name)
                        },
                        coefficients: coefficients(f, m, condition.alpha),
                        body_force: f,
                        body_moment: Vec3::new(-m.x, m.y, -m.z),
                        area: area_sum,
                        strip_loads: strips,
                    });
                }
            }
        }
        let mut coeff = coefficients(
            total_force / reference.area,
            total_moment / (reference.area * reference.chord),
            condition.alpha,
        );
        let viscous = viscous_area / reference.area + profile;
        coeff.drag += profile;
        for output in &request.outputs {
            let value = match output {
                Output::Lift => coeff.lift,
                Output::Drag => coeff.drag,
                Output::Side => coeff.side,
                Output::Roll => coeff.roll,
                Output::Pitch => coeff.pitch,
                Output::Yaw => coeff.yaw,
                Output::NearFieldInducedDrag => coeff.drag - viscous,
                Output::ViscousDrag => viscous,
                _ => continue,
            };
            values.insert(*output, value);
        }
        work.near_field_panels += mesh.vortex_count();
        if plan.moments {
            work.moment_panels += mesh.vortex_count();
        }
    }
    if plan.far {
        let far = super::trefftz::evaluate(solver, half, plan.far_drag, work);
        for output in &request.outputs {
            if output.far() {
                values.insert(*output, far.value(*output, reference));
            }
        }
    }
    Ok((values, details))
}

pub(crate) fn expand_circulation(solver: &VortexLatticeSolver, half: bool) -> VecX {
    if !half {
        return solver.solution().circulation.clone();
    }
    let mut data = Vec::with_capacity(2 * solver.mesh().vortex_count());
    for s in solver.mesh().surface_info() {
        for _ in 0..2 {
            data.extend(
                solver
                    .solution()
                    .circulation
                    .rows(s.first_vortex_index, s.num_vortices)
                    .iter()
                    .copied(),
            );
        }
    }
    VecX::from_vec(data)
}
pub(crate) fn expand_velocity(solver: &VortexLatticeSolver, half: bool) -> MatX {
    if !half {
        return solver.solution().induced_velocity.clone();
    }
    let mut rows = Vec::with_capacity(2 * solver.mesh().vortex_count());
    for s in solver.mesh().surface_info() {
        for image in 0..2 {
            for i in s.first_vortex_index..s.first_vortex_index + s.num_vortices {
                let v = induced(solver, i);
                rows.push(if image == 0 { v } else { reflect(v) });
            }
        }
    }
    MatX::from_fn(rows.len(), 3, |i, j| rows[i][j])
}
