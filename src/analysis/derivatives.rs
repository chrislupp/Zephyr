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

//! Selected analytical derivatives. One adjoint per requested scalar output.
use super::operator::{core, reflect};
use super::{
    AnalysisRequest, GradientRequest, GradientResult, GradientValue, Input, Output, SolverMode,
    Source, Work, unsupported,
};
use crate::core::{EPSILON, MatX, Scalar, Vec3, VecX};
use crate::error::{Result, ZephyrError};
use crate::gradients::geometry_adjoint::{LatticeCotangents, scatter_mesh_node_cotangents};
use crate::gradients::section_geometry::section_direction;
use crate::gradients::{horseshoe_velocity_jacobians, horseshoe_velocity_pg_derivative};
use crate::solver::{FlightCondition, VortexLatticeSolver, horseshoe_vortex_velocity};

pub(crate) fn validate(
    source: &Source,
    s: &VortexLatticeSolver,
    c: &FlightCondition,
    r: &AnalysisRequest,
) -> Result<()> {
    for g in &r.gradients {
        if r.mode == SolverMode::HalfWing
            && matches!(g.input, Input::Beta | Input::RollRate | Input::YawRate)
        {
            return Err(unsupported(
                "this derivative requires an antisymmetric response, which half-wing mode does not yet support",
            ));
        }
        if g.input.geometry()
            && (!c.control_deflections.is_empty()
                || s.mesh().strips().iter().any(|v| v.has_viscous_data))
        {
            return Err(ZephyrError::InvalidGeometry(
                "geometry adjoints require undeflected controls and no section drag polars".into(),
            ));
        }
        match (&g.input, source) {
            (Input::MomentReference(axis), _) if *axis > 2 => {
                return Err(ZephyrError::InvalidGeometry(
                    "moment-reference axis must be 0, 1, or 2".into(),
                ));
            }
            (Input::MeshNodes(i), Source::Mesh(m)) if *i < m.surfaces.len() => {}
            (Input::MeshNodes(_), _) => {
                return Err(ZephyrError::InvalidGeometry(
                    "mesh-node derivatives require a valid direct-mesh surface index".into(),
                ));
            }
            (
                Input::Section {
                    surface, section, ..
                },
                Source::Aircraft(a),
            ) if a
                .surfaces()
                .get(*surface)
                .is_some_and(|s| *section < s.sections().len()) => {}
            (Input::Section { .. }, _) => {
                return Err(ZephyrError::InvalidGeometry(
                    "section derivatives require a valid aircraft section".into(),
                ));
            }
            (Input::Control(name), _) => {
                let influences: Vec<_> = s
                    .mesh()
                    .control_influences()
                    .iter()
                    .flatten()
                    .filter(|i| &i.name == name)
                    .collect();
                if influences.is_empty() {
                    return Err(ZephyrError::InvalidGeometry(format!(
                        "unknown control '{name}'"
                    )));
                }
                if r.mode == SolverMode::HalfWing
                    && influences.iter().any(|i| i.reflection_sign != 1.0)
                {
                    return Err(unsupported(
                        "an antisymmetric control derivative is unavailable in half-wing mode",
                    ));
                }
            }
            _ => {}
        }
    }
    Ok(())
}

struct Seeds {
    gamma: VecX,
    velocity: Vec<Vec3>,
    geometry: LatticeCotangents,
    alpha: Scalar,
    beta: Scalar,
    reference: [Scalar; 3],
}

fn axes(o: Output, a: Scalar) -> (Vec3, Vec3, Vec3) {
    let (s, c) = a.sin_cos();
    match o {
        Output::Lift => (Vec3::new(-s, 0.0, c), Vec3::zeros(), Vec3::new(-c, 0.0, -s)),
        Output::Drag | Output::NearFieldInducedDrag => {
            (Vec3::new(c, 0.0, s), Vec3::zeros(), Vec3::new(-s, 0.0, c))
        }
        Output::Side => (Vec3::y(), Vec3::zeros(), Vec3::zeros()),
        Output::Roll => (Vec3::zeros(), -Vec3::x(), Vec3::zeros()),
        Output::Pitch => (Vec3::zeros(), Vec3::y(), Vec3::zeros()),
        Output::Yaw => (Vec3::zeros(), -Vec3::z(), Vec3::zeros()),
        _ => (Vec3::zeros(), Vec3::zeros(), Vec3::zeros()),
    }
}
fn freestream_derivatives(c: &FlightCondition) -> (Vec3, Vec3) {
    let (sa, ca) = c.alpha.sin_cos();
    let (sb, cb) = c.beta.sin_cos();
    (
        Vec3::new(-sa * cb, 0.0, ca * cb),
        Vec3::new(-ca * sb, cb, -sa * sb),
    )
}

#[derive(Clone, Copy)]
struct Needs {
    geometry: bool,
    alpha: bool,
    beta: bool,
    reference: bool,
}

fn near_seeds(
    s: &VortexLatticeSolver,
    c: &FlightCondition,
    half: bool,
    o: Output,
    needs: Needs,
    work: &mut Work,
) -> Seeds {
    let mesh = s.mesh();
    let n = mesh.vortex_count();
    let r = s.reference_quantities();
    let (fs, ms, dfs) = axes(o, c.alpha);
    let u = VortexLatticeSolver::freestream_direction(c);
    let (ua, ub) = freestream_derivatives(c);
    let mut out = Seeds {
        gamma: VecX::zeros(n),
        velocity: vec![Vec3::zeros(); n],
        geometry: LatticeCotangents::zeros(mesh),
        alpha: 0.0,
        beta: 0.0,
        reference: [0.0; 3],
    };
    for (i, v) in mesh.vortices().iter().enumerate() {
        let line = v.bound_vortex_right - v.bound_vortex_left;
        let velocity = u + super::loads::induced(s, i);
        let gamma = s.solution().circulation[i];
        let force = 2.0 * gamma * velocity.cross(&line);
        let mut combined = fs / r.area
            + ms.cross(&(v.control_point - r.moment_reference_point)) / (r.area * r.chord);
        if half {
            combined += reflect(
                fs / r.area
                    + ms.cross(&(reflect(v.control_point) - r.moment_reference_point))
                        / (r.area * r.chord),
            );
        }
        out.gamma[i] = 2.0 * combined.dot(&velocity.cross(&line));
        out.velocity[i] = 2.0 * gamma * line.cross(&combined);
        if needs.alpha {
            out.alpha += out.velocity[i].dot(&ua);
            out.alpha +=
                dfs.dot(&(force + if half { reflect(force) } else { Vec3::zeros() })) / r.area;
        }
        if needs.beta {
            out.beta += out.velocity[i].dot(&ub);
        }
        if needs.reference {
            let reference_seed = -(force + if half { reflect(force) } else { Vec3::zeros() })
                .cross(&ms)
                / (r.area * r.chord);
            for k in 0..3 {
                out.reference[k] += reference_seed[k];
            }
        }
        if needs.geometry {
            let line_seed = 2.0 * gamma * combined.cross(&velocity);
            out.geometry.bound_right[i] += line_seed;
            out.geometry.bound_left[i] -= line_seed;
            out.geometry.control_point[i] += force.cross(&ms) / (r.area * r.chord);
            if half {
                out.geometry.control_point[i] +=
                    reflect(reflect(force).cross(&ms)) / (r.area * r.chord);
            }
        }
    }
    // Profile polar derivatives follow the strip lift and include both the
    // changing onset and the rotating lift/drag axes. Geometry polars are
    // explicitly rejected, matching the existing geometry-adjoint limitation.
    if !matches!(
        o,
        Output::Lift | Output::Side | Output::NearFieldInducedDrag
    ) {
        let (sa, ca) = c.alpha.sin_cos();
        let lift = Vec3::new(-sa, 0.0, ca);
        let dlift = Vec3::new(-ca, 0.0, -sa);
        let drag = Vec3::new(ca, 0.0, sa);
        let ddrag = Vec3::new(-sa, 0.0, ca);
        for strip in mesh.strips() {
            if !strip.has_viscous_data {
                continue;
            }
            let area = strip.chord_center * strip.width;
            if area <= EPSILON {
                continue;
            }
            let mut force = Vec3::zeros();
            let mut da = Vec3::zeros();
            let mut db = Vec3::zeros();
            for i in strip.first_vortex_index..strip.first_vortex_index + strip.num_vortices {
                let v = &mesh.vortices()[i];
                let line = v.bound_vortex_right - v.bound_vortex_left;
                let scale = 2.0 * s.solution().circulation[i];
                force += scale * (u + super::loads::induced(s, i)).cross(&line);
                if needs.alpha {
                    da += scale * ua.cross(&line);
                }
                if needs.beta {
                    db += scale * ub.cross(&line);
                }
            }
            let (cd, slope) = crate::aero::forces::drag_polar_value_and_slope(
                &strip.viscous_polar,
                lift.dot(&force) / area,
            );
            if cd <= 0.0 {
                continue;
            }
            let qc = strip.leading_edge_center + Vec3::new(0.25 * strip.chord_center, 0.0, 0.0);
            let mut weight = 0.0;
            let mut alpha_weight = 0.0;
            for image in 0..if half { 2 } else { 1 } {
                let q = if image == 0 { qc } else { reflect(qc) };
                let seed = if o == Output::ViscousDrag {
                    1.0 / r.area
                } else {
                    (fs + ms.cross(&(q - r.moment_reference_point)) / r.chord).dot(&drag) / r.area
                };
                weight += seed;
                if o != Output::ViscousDrag {
                    if needs.alpha {
                        alpha_weight += (dfs.dot(&drag)
                            + (fs + ms.cross(&(q - r.moment_reference_point)) / r.chord)
                                .dot(&ddrag))
                            / r.area;
                    }
                    if needs.reference {
                        let rs = -area * cd * drag.cross(&ms) / (r.area * r.chord);
                        for k in 0..3 {
                            out.reference[k] += rs[k];
                        }
                    }
                }
            }
            if needs.alpha {
                out.alpha +=
                    weight * slope * (dlift.dot(&force) + lift.dot(&da)) + area * cd * alpha_weight;
            }
            if needs.beta {
                out.beta += weight * slope * lift.dot(&db);
            }
            for i in strip.first_vortex_index..strip.first_vortex_index + strip.num_vortices {
                let v = &mesh.vortices()[i];
                let line = v.bound_vortex_right - v.bound_vortex_left;
                out.gamma[i] += 2.0
                    * weight
                    * slope
                    * lift.dot(&(u + super::loads::induced(s, i)).cross(&line));
                out.velocity[i] +=
                    2.0 * s.solution().circulation[i] * weight * slope * line.cross(&lift);
            }
        }
    }
    let flat = VecX::from_iterator(3 * n, out.velocity.iter().flat_map(|v| v.iter().copied()));
    out.gamma += s.velocity_matrix().transpose() * flat;
    work.force_derivative_panels += n;
    if needs.alpha {
        work.alpha_derivative_panels += n;
    }
    if needs.beta {
        work.beta_derivative_panels += n;
    }
    out
}

pub(crate) fn evaluate(
    source: &Source,
    s: &VortexLatticeSolver,
    c: &FlightCondition,
    half: bool,
    requests: &[GradientRequest],
    work: &mut Work,
    cached_adjoints: &mut std::collections::BTreeMap<Output, VecX>,
) -> Result<Vec<GradientResult>> {
    let implicit: Vec<_> = requests
        .iter()
        .filter(|g| {
            !matches!(
                g.input,
                Input::ReferenceArea
                    | Input::ReferenceChord
                    | Input::ReferenceSpan
                    | Input::DynamicPressure
            )
        })
        .cloned()
        .collect();
    let mut result = if implicit.is_empty() {
        Vec::new()
    } else {
        evaluate_implicit(source, s, c, half, &implicit, work, cached_adjoints)?
    };
    let normalizing: Vec<_> = requests
        .iter()
        .filter(|g| {
            matches!(
                g.input,
                Input::ReferenceArea
                    | Input::ReferenceChord
                    | Input::ReferenceSpan
                    | Input::DynamicPressure
            )
        })
        .collect();
    let req = AnalysisRequest {
        outputs: normalizing
            .iter()
            .filter(|g| !g.input.known_zero(g.output))
            .map(|g| g.output)
            .collect(),
        ..AnalysisRequest::default()
    };
    let values = if req.outputs.is_empty() {
        std::collections::BTreeMap::new()
    } else {
        super::loads::evaluate(s, c, half, source.profile_drag(), &req, work)?.0
    };
    for g in normalizing {
        let value = if g.input.known_zero(g.output) {
            0.0
        } else {
            let v = values[&g.output];
            let reference = s.reference_quantities();
            match g.input {
                Input::ReferenceArea => {
                    -(v - if matches!(g.output, Output::Drag | Output::ViscousDrag) {
                        source.profile_drag()
                    } else {
                        0.0
                    }) / reference.area
                }
                Input::ReferenceChord => -v / reference.chord,
                Input::ReferenceSpan => -2.0 * v / reference.span,
                _ => unreachable!(),
            }
        };
        result.push(GradientResult {
            request: g.clone(),
            value: GradientValue::Scalar(value),
        });
    }
    Ok(result)
}

fn evaluate_implicit(
    source: &Source,
    s: &VortexLatticeSolver,
    c: &FlightCondition,
    half: bool,
    requests: &[GradientRequest],
    work: &mut Work,
    cached_adjoints: &mut std::collections::BTreeMap<Output, VecX>,
) -> Result<Vec<GradientResult>> {
    let mut outputs: Vec<_> = requests.iter().map(|r| r.output).collect();
    outputs.sort();
    outputs.dedup();
    let n = s.mesh().vortex_count();
    let mut seeds = Vec::new();
    let missing: Vec<_> = outputs
        .iter()
        .copied()
        .filter(|o| !cached_adjoints.contains_key(o))
        .collect();
    let mut rhs = MatX::zeros(n, missing.len());
    for &o in &outputs {
        let geometry = requests.iter().any(|r| r.output == o && r.input.geometry());
        let seed = if o.far() {
            let (gamma, geometry) = super::trefftz::reverse(s, half, o, geometry, work);
            Seeds {
                gamma,
                geometry,
                velocity: vec![Vec3::zeros(); n],
                alpha: 0.0,
                beta: 0.0,
                reference: [0.0; 3],
            }
        } else {
            let needs = Needs {
                geometry,
                alpha: requests
                    .iter()
                    .any(|r| r.output == o && r.input == Input::Alpha),
                beta: requests
                    .iter()
                    .any(|r| r.output == o && r.input == Input::Beta),
                reference: requests
                    .iter()
                    .any(|r| r.output == o && matches!(r.input, Input::MomentReference(_))),
            };
            near_seeds(s, c, half, o, needs, work)
        };
        if let Some(column) = missing.iter().position(|output| *output == o) {
            rhs.set_column(column, &seed.gamma);
        }
        seeds.push(seed);
    }
    if !missing.is_empty() {
        let solved = s.solve_linearized_transpose_batched(&rhs)?;
        for (j, output) in missing.iter().enumerate() {
            cached_adjoints.insert(*output, solved.column(j).into_owned());
        }
        work.adjoint_solves += missing.len();
    }
    let adjoints = MatX::from_fn(n, outputs.len(), |i, j| cached_adjoints[&outputs[j]][i]);
    let needs_mach = requests.iter().any(|r| r.input == Input::Mach);
    let mach_velocity = if needs_mach {
        Some(mach_velocity(s, c, half, work))
    } else {
        None
    };
    let geometry_columns: Vec<_> = outputs
        .iter()
        .enumerate()
        .filter_map(|(i, o)| {
            requests
                .iter()
                .any(|r| r.output == *o && r.input.geometry())
                .then_some(i)
        })
        .collect();
    if !geometry_columns.is_empty() {
        geometry_reverse(s, c, half, &adjoints, &geometry_columns, &mut seeds, work);
    }
    let mut result = Vec::new();
    for request in requests {
        let column = outputs.binary_search(&request.output).unwrap();
        let seed = &seeds[column];
        let lambda = adjoints.column(column);
        let implicit = |velocity: Vec3| -> Scalar {
            s.mesh()
                .vortices()
                .iter()
                .enumerate()
                .filter(|(_, v)| v.affected_by_freestream)
                .map(|(i, _)| -lambda[i] * s.linearized_normals()[i].dot(&velocity))
                .sum()
        };
        let value = match &request.input {
            Input::DynamicPressure => GradientValue::Scalar(0.0),
            Input::Alpha => {
                GradientValue::Scalar(seed.alpha + implicit(freestream_derivatives(c).0))
            }
            Input::Beta => GradientValue::Scalar(seed.beta + implicit(freestream_derivatives(c).1)),
            Input::Mach => {
                let v = mach_velocity.as_ref().unwrap();
                let mut value = 0.0;
                for i in 0..n {
                    value += (seed.velocity[i] - lambda[i] * s.linearized_normals()[i]).dot(&v[i]);
                }
                GradientValue::Scalar(value)
            }
            Input::RollRate | Input::PitchRate | Input::YawRate => {
                let axis = match request.input {
                    Input::RollRate => Vec3::x(),
                    Input::PitchRate => Vec3::y(),
                    _ => Vec3::z(),
                };
                let value = s
                    .mesh()
                    .vortices()
                    .iter()
                    .enumerate()
                    .filter(|(_, v)| v.affected_by_freestream)
                    .map(|(i, v)| {
                        -lambda[i]
                            * s.linearized_normals()[i]
                                .dot(&axis.cross(&(v.control_point - s.reference_point())))
                    })
                    .sum();
                GradientValue::Scalar(value)
            }
            Input::MomentReference(axis) => {
                let mut unit = Vec3::zeros();
                unit[*axis] = 1.0;
                let omega = Vec3::new(c.roll_rate, c.pitch_rate, c.yaw_rate);
                GradientValue::Scalar(seed.reference[*axis] + implicit(-omega.cross(&unit)))
            }
            Input::Control(name) => {
                let (normals, dn) = s.effective_normals_with_derivative(c, name);
                let u = VortexLatticeSolver::freestream_direction(c);
                let omega = Vec3::new(c.roll_rate, c.pitch_rate, c.yaw_rate);
                let mut value = 0.0;
                for (i, v) in s.mesh().vortices().iter().enumerate() {
                    let onset = if v.affected_by_freestream {
                        u + omega.cross(&(v.control_point - s.reference_point()))
                    } else {
                        Vec3::zeros()
                    };
                    value -= lambda[i] * dn[i].dot(&(onset + super::loads::induced(s, i)));
                }
                let _ = normals;
                GradientValue::Scalar(value)
            }
            Input::ReferenceArea | Input::ReferenceChord | Input::ReferenceSpan => unreachable!(),
            Input::MeshNodes(index) => {
                let Source::Mesh(input) = source else {
                    unreachable!()
                };
                GradientValue::MeshNodes(scatter_mesh_node_cotangents(
                    input,
                    s.mesh(),
                    *index,
                    &seed.geometry,
                )?)
            }
            Input::Section {
                surface,
                section,
                variable,
            } => {
                let Source::Aircraft(aircraft) = source else {
                    unreachable!()
                };
                let direction =
                    section_direction(aircraft, s.mesh(), *surface, *section, *variable)?;
                GradientValue::Scalar(seed.geometry.dot(&direction))
            }
        };
        result.push(GradientResult {
            request: request.clone(),
            value,
        });
    }
    Ok(result)
}

fn mach_velocity(
    s: &VortexLatticeSolver,
    c: &FlightCondition,
    half: bool,
    work: &mut Work,
) -> Vec<Vec3> {
    let n = s.mesh().vortex_count();
    let pg = (1.0 - c.mach * c.mach).sqrt();
    let mut result = vec![Vec3::zeros(); n];
    if c.mach == 0.0 {
        return result;
    }
    for (i, target) in s.mesh().vortices().iter().enumerate() {
        for (j, source) in s.mesh().vortices().iter().enumerate() {
            for image in 0..if half { 2 } else { 1 } {
                let left = if image == 0 {
                    source.bound_vortex_left
                } else {
                    reflect(source.bound_vortex_right)
                };
                let right = if image == 0 {
                    source.bound_vortex_right
                } else {
                    reflect(source.bound_vortex_left)
                };
                let derivative = horseshoe_velocity_pg_derivative(
                    &target.control_point,
                    &left,
                    &right,
                    pg,
                    core(s.mesh(), i, j),
                    image == 1 || i != j,
                );
                result[i] += (-c.mach / pg) * s.solution().circulation[j] * derivative;
                work.mach_interactions += 1;
            }
        }
    }
    result
}

fn geometry_reverse(
    s: &VortexLatticeSolver,
    c: &FlightCondition,
    half: bool,
    adjoints: &MatX,
    columns: &[usize],
    seeds: &mut [Seeds],
    work: &mut Work,
) {
    let pg = (1.0 - c.mach * c.mach).sqrt();
    let u = VortexLatticeSolver::freestream_direction(c);
    let omega = Vec3::new(c.roll_rate, c.pitch_rate, c.yaw_rate);
    for (i, target) in s.mesh().vortices().iter().enumerate() {
        let normal = s.linearized_normals()[i];
        for &k in columns {
            let lambda = adjoints[(i, k)];
            if target.affected_by_freestream {
                seeds[k].geometry.normal[i] -=
                    lambda * (u + omega.cross(&(target.control_point - s.reference_point())));
                seeds[k].geometry.control_point[i] -= lambda * normal.cross(&omega);
            }
        }
        for (j, source) in s.mesh().vortices().iter().enumerate() {
            for image in 0..if half { 2 } else { 1 } {
                let left = if image == 0 {
                    source.bound_vortex_left
                } else {
                    reflect(source.bound_vortex_right)
                };
                let right = if image == 0 {
                    source.bound_vortex_right
                } else {
                    reflect(source.bound_vortex_left)
                };
                let radius = core(s.mesh(), i, j);
                let bound = image == 1 || i != j;
                let jac = horseshoe_velocity_jacobians(
                    &target.control_point,
                    &left,
                    &right,
                    pg,
                    radius,
                    bound,
                );
                let velocity = horseshoe_vortex_velocity(
                    &target.control_point,
                    &left,
                    &right,
                    pg,
                    radius,
                    bound,
                );
                let gamma = s.solution().circulation[j];
                for &k in columns {
                    let seed = gamma * (seeds[k].velocity[i] - adjoints[(i, k)] * normal);
                    let l = jac.bound_left.tr_mul(&seed);
                    let r = jac.bound_right.tr_mul(&seed);
                    if image == 0 {
                        seeds[k].geometry.bound_left[j] += l;
                        seeds[k].geometry.bound_right[j] += r;
                    } else {
                        seeds[k].geometry.bound_right[j] += reflect(l);
                        seeds[k].geometry.bound_left[j] += reflect(r);
                    }
                    seeds[k].geometry.control_point[i] += jac.field_point.tr_mul(&seed);
                    seeds[k].geometry.normal[i] -= adjoints[(i, k)] * gamma * velocity;
                }
                work.geometry_interactions += 1;
            }
        }
    }
}
