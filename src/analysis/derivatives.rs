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
use super::operator::{SourcePanel, reflect, source_panels};
use super::{
    AnalysisRequest, GradientRequest, GradientResult, GradientValue, Input, Output, SolverMode,
    Source, Work, unsupported,
};
use crate::core::{EPSILON, MatX, Scalar, Vec3, VecX};
use crate::error::{Result, ZephyrError};
use crate::gradients::geometry_adjoint::{LatticeCotangents, scatter_mesh_node_cotangents};
use crate::gradients::section_geometry::section_direction;
use crate::gradients::{HorseshoeVelocityJacobians, horseshoe_velocity_pg_derivative};
use crate::solver::{FlightCondition, VortexLatticeSolver};

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
            + ms.cross(&(v.bound_vortex_center - r.moment_reference_point)) / (r.area * r.chord);
        if half {
            combined += reflect(
                fs / r.area
                    + ms.cross(&(reflect(v.bound_vortex_center) - r.moment_reference_point))
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
            let mut arm = force.cross(&ms) / (r.area * r.chord);
            if half {
                arm += reflect(reflect(force).cross(&ms)) / (r.area * r.chord);
            }
            out.geometry.bound_center[i] += arm;
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
    // `tr_mul` contracts in place; `transpose() *` would first copy the whole
    // `3N x N` operator.
    out.gamma += s.velocity_matrix().tr_mul(&flat);
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
                let (control, bound) = mach_velocity.as_ref().unwrap();
                let mut value = 0.0;
                for i in 0..n {
                    value += seed.velocity[i].dot(&bound[i])
                        - lambda[i] * s.linearized_normals()[i].dot(&control[i]);
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
                    let velocity =
                        s.control_velocity_matrix().rows(3 * i, 3) * &s.solution().circulation;
                    value -= lambda[i]
                        * dn[i].dot(&(onset + Vec3::new(velocity[0], velocity[1], velocity[2])));
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

/// Mach derivative of the induced velocity at fixed circulation, at every
/// control point and bound centre: `(dV_control/dM Γ, dV_bound/dM Γ)`.
///
/// Zero, with no lattice work at all, at Mach `<= EPSILON`, where the
/// derivative of a `sqrt(1 - M²)` factor vanishes. Full mode is exactly the
/// solver's operator (no image planes, the same core radii and bound flags),
/// so it reuses the solver's fused row kernel,
/// `VortexLatticeSolver::velocity_mach_derivative_times`. Half-wing mode
/// adds each source's reflection, which that kernel does not know, so it runs
/// the same row-wise structure here: one self-contained, column-ordered sum
/// per target row, rows threaded under `parallel` and collected in order, so
/// results are bit-identical between builds and across thread counts.
fn mach_velocity(
    s: &VortexLatticeSolver,
    c: &FlightCondition,
    half: bool,
    work: &mut Work,
) -> (Vec<Vec3>, Vec<Vec3>) {
    let n = s.mesh().vortex_count();
    if c.mach <= EPSILON {
        return (vec![Vec3::zeros(); n], vec![Vec3::zeros(); n]);
    }
    let gamma = &s.solution().circulation;
    work.mach_interactions += 2 * n * n * if half { 2 } else { 1 };
    if !half {
        let triples = |flat: VecX| -> Vec<Vec3> {
            (0..n)
                .map(|i| Vec3::new(flat[3 * i], flat[3 * i + 1], flat[3 * i + 2]))
                .collect()
        };
        return (
            triples(s.velocity_mach_derivative_times(false, gamma)),
            triples(s.velocity_mach_derivative_times(true, gamma)),
        );
    }
    let pg = (1.0 - c.mach * c.mach).sqrt();
    let direction = -c.mach / pg;
    let sources = source_panels(s.mesh());
    let vortices = s.mesh().vortices();
    let row = |i: usize| -> (Vec3, Vec3) {
        let target = &vortices[i];
        let (mut control, mut bound) = (Vec3::zeros(), Vec3::zeros());
        for (j, (source, g)) in sources.iter().zip(gamma.iter()).enumerate() {
            let core = source.core(target.component_index);
            let derivative = |point: &Vec3, image: usize, include_bound: bool| {
                let (left, right) = source.endpoints(image);
                horseshoe_velocity_pg_derivative(point, &left, &right, pg, core, include_bound)
            };
            let at_control = derivative(&target.control_point, 0, true)
                + derivative(&target.control_point, 1, true);
            let at_bound = derivative(&target.bound_vortex_center, 0, i != j)
                + derivative(&target.bound_vortex_center, 1, true);
            control += *g * (at_control * direction);
            bound += *g * (at_bound * direction);
        }
        (control, bound)
    };
    drive_rows(n, row).into_iter().unzip()
}

/// Evaluates `row` for every target row and returns the values in row order:
/// one rayon task per row under `parallel` above the threading threshold.
/// Rows share nothing, so the result never depends on the schedule.
fn drive_rows<T: Send>(n: usize, row: impl Fn(usize) -> T + Sync) -> Vec<T> {
    #[cfg(feature = "parallel")]
    if n >= super::operator::PARALLEL_MIN_PANELS {
        use rayon::prelude::*;
        return (0..n).into_par_iter().map(&row).collect();
    }
    (0..n).map(row).collect()
}

/// Seeds at or below this count are contracted as separate single-seed VJPs;
/// above it one set of Jacobians is cheaper. Matches the geometry adjoint.
const MAX_SEPARATE_SEEDS: usize = 3;

/// Adds the endpoint cotangents of a direct source (`image == 0`) or of its
/// reflection, which maps back to the swapped, reflected endpoints.
#[inline]
fn add_endpoints(
    geometry: &mut LatticeCotangents,
    j: usize,
    image: usize,
    left: Vec3,
    right: Vec3,
) {
    if image == 0 {
        geometry.bound_left[j] += left;
        geometry.bound_right[j] += right;
    } else {
        geometry.bound_right[j] += reflect(left);
        geometry.bound_left[j] += reflect(right);
    }
}

/// Geometry cotangents of the residual and near-field loads for each column in
/// `columns`.
///
/// Per (target, source, image) the residual term contracts the control-point
/// Jacobian with the single seed `normal`, shared by every output, so one
/// endpoint VJP replaces the three Jacobians and each output only scales it by
/// `-λ Γ`. The normal's own cotangent is `-λ` times the induced velocity at the
/// control point, `Σ_j Γ_j V_c(i, j)`, which is the stored control-point
/// operator times the circulation; it is formed once rather than re-evaluating
/// every kernel. The load term has one seed per output, contracted as separate
/// VJPs for a few outputs and through one set of Jacobians for many, and is
/// skipped for outputs with no near-field seed at that row (far-field ones).
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
    let n = s.mesh().vortex_count();
    let images = if half { 2 } else { 1 };
    let gamma = &s.solution().circulation;
    let sources: Vec<SourcePanel> = source_panels(s.mesh());
    let control_velocity = s.control_velocity_matrix() * gamma;
    // Fixed ranges preserve reduction order across serial/parallel builds and
    // pool sizes. Workers own their partial cotangents; no shared additions.
    let rows_per_chunk = n.div_ceil(32).max(1);
    let ranges: Vec<_> = (0..n)
        .step_by(rows_per_chunk)
        .map(|first| (first, (first + rows_per_chunk).min(n)))
        .collect();
    let seeds_ref = &*seeds;
    let sweep = |&(first, last): &(usize, usize)| {
        let mut geometry: Vec<_> = columns
            .iter()
            .map(|_| LatticeCotangents::zeros(s.mesh()))
            .collect();
        let mut residual_weights = Vec::with_capacity(columns.len());
        let mut load_seeds = Vec::with_capacity(columns.len());
        for i in first..last {
            let target = &s.mesh().vortices()[i];
            let normal = s.linearized_normals()[i];
            let induced = Vec3::new(
                control_velocity[3 * i],
                control_velocity[3 * i + 1],
                control_velocity[3 * i + 2],
            );
            for (local, &k) in columns.iter().enumerate() {
                let lambda = adjoints[(i, k)];
                if target.affected_by_freestream {
                    geometry[local].normal[i] -=
                        lambda * (u + omega.cross(&(target.control_point - s.reference_point())));
                    geometry[local].control_point[i] -= lambda * normal.cross(&omega);
                }
                geometry[local].normal[i] -= lambda * induced;
            }
            residual_weights.clear();
            residual_weights.extend(
                columns
                    .iter()
                    .enumerate()
                    .map(|(local, &k)| (local, adjoints[(i, k)]))
                    .filter(|(_, lambda)| *lambda != 0.0),
            );
            load_seeds.clear();
            load_seeds.extend(
                columns
                    .iter()
                    .enumerate()
                    .map(|(local, &k)| (local, seeds_ref[k].velocity[i]))
                    .filter(|(_, seed)| seed.norm_squared() > 0.0),
            );
            for (j, source) in sources.iter().enumerate() {
                let g = gamma[j];
                if g == 0.0 {
                    continue;
                }
                let core = source.core(target.component_index);
                for image in 0..images {
                    let (left, right) = source.endpoints(image);
                    if !residual_weights.is_empty() {
                        let (by_left, by_right) = HorseshoeVelocityJacobians::shared_terms(
                            &target.control_point,
                            &left,
                            &right,
                            pg,
                            core,
                            true,
                        )
                        .endpoint_vjp(&normal);
                        let by_point = -(by_left + by_right);
                        for &(local, lambda) in &residual_weights {
                            let weight = -lambda * g;
                            add_endpoints(
                                &mut geometry[local],
                                j,
                                image,
                                weight * by_left,
                                weight * by_right,
                            );
                            geometry[local].control_point[i] += weight * by_point;
                        }
                    }
                    if load_seeds.is_empty() {
                        continue;
                    }
                    let kernel = HorseshoeVelocityJacobians::shared_terms(
                        &target.bound_vortex_center,
                        &left,
                        &right,
                        pg,
                        core,
                        image == 1 || i != j,
                    );
                    if load_seeds.len() <= MAX_SEPARATE_SEEDS {
                        for &(local, seed) in &load_seeds {
                            let (by_left, by_right) = kernel.endpoint_vjp(&(g * seed));
                            geometry[local].bound_center[i] -= by_left + by_right;
                            add_endpoints(&mut geometry[local], j, image, by_left, by_right);
                        }
                    } else {
                        let jacobians = kernel.jacobians();
                        for &(local, seed) in &load_seeds {
                            let scaled = g * seed;
                            geometry[local].bound_center[i] +=
                                jacobians.field_point.tr_mul(&scaled);
                            add_endpoints(
                                &mut geometry[local],
                                j,
                                image,
                                jacobians.bound_left.tr_mul(&scaled),
                                jacobians.bound_right.tr_mul(&scaled),
                            );
                        }
                    }
                }
            }
        }
        geometry
    };
    #[cfg(feature = "parallel")]
    let chunks: Vec<_> = {
        use rayon::prelude::*;
        if n >= super::operator::PARALLEL_MIN_PANELS {
            ranges.par_iter().map(sweep).collect()
        } else {
            ranges.iter().map(sweep).collect()
        }
    };
    #[cfg(not(feature = "parallel"))]
    let chunks: Vec<_> = ranges.iter().map(sweep).collect();
    for chunk in chunks {
        for (&k, part) in columns.iter().zip(chunk) {
            let target = &mut seeds[k].geometry;
            for i in 0..n {
                target.bound_left[i] += part.bound_left[i];
                target.bound_right[i] += part.bound_right[i];
                target.control_point[i] += part.control_point[i];
                target.bound_center[i] += part.bound_center[i];
                target.normal[i] += part.normal[i];
            }
        }
    }
    work.geometry_interactions += 2 * n * n * images;
}

#[cfg(test)]
mod tests {
    use super::super::Settings;
    use super::*;
    use crate::geometry::{MeshInput, SurfaceMesh};
    use crate::gradients::horseshoe_velocity_jacobians;
    use crate::solver::horseshoe_vortex_velocity;

    /// A wing and a tail on separate components, duplicated about `y = 0`.
    fn input(ny: usize) -> MeshInput {
        let surface = |name: &str, nx: usize, ny: usize, component: usize, offset: Vec3| {
            let mut surface = SurfaceMesh::new(name, nx, ny);
            surface.y_duplicated = true;
            surface.component_index = component;
            for i in 0..=nx {
                for j in 0..=ny {
                    let eta = j as Scalar / ny as Scalar;
                    let chord = i as Scalar / nx as Scalar;
                    let point = Vec3::new(
                        0.6 * eta + chord * (1.0 - 0.4 * eta),
                        0.05 + 3.0 * eta,
                        0.3 * eta - 0.05 * chord * eta,
                    );
                    surface.set_node(i, j, offset + point).unwrap();
                }
            }
            surface
        };
        let mut input = MeshInput {
            surfaces: vec![
                surface("wing", 3, ny, 0, Vec3::zeros()),
                surface("tail", 2, 3, 1, Vec3::new(4.0, 0.0, 0.5)),
            ],
            ..MeshInput::default()
        };
        input.reference.area = 6.0;
        input.reference.chord = 0.8;
        input.reference.span = 6.0;
        input.reference.moment_reference_point = Vec3::new(0.3, 0.0, 0.1);
        input
    }

    fn solved(half: bool, mach: Scalar, ny: usize) -> (VortexLatticeSolver, FlightCondition) {
        let source = Source::Mesh(input(ny));
        let source = if half { source.half().unwrap() } else { source };
        let mut solver = source.solver(Settings::default()).unwrap();
        let condition = FlightCondition {
            alpha: 0.08,
            mach,
            pitch_rate: 0.01,
            ..FlightCondition::default()
        };
        solver
            .solve_selective(&condition, half, true, false, &mut Work::default())
            .unwrap();
        (solver, condition)
    }

    fn core(mesh: &crate::solver::VortexLattice, row: usize, column: usize) -> Scalar {
        let source = &mesh.vortices()[column];
        let width = (source.bound_vortex_right.y - source.bound_vortex_left.y)
            .hypot(source.bound_vortex_right.z - source.bound_vortex_left.z);
        if mesh.vortices()[row].component_index == source.component_index {
            0.0001 * width
        } else {
            (0.0001 * source.local_chord).max(0.0001 * width)
        }
    }

    fn image_endpoints(source: &crate::solver::HorseshoeVortex, image: usize) -> (Vec3, Vec3) {
        if image == 0 {
            (source.bound_vortex_left, source.bound_vortex_right)
        } else {
            (
                reflect(source.bound_vortex_right),
                reflect(source.bound_vortex_left),
            )
        }
    }

    /// The baseline Mach derivative: serial, one pair at a time.
    fn reference_mach(
        s: &VortexLatticeSolver,
        c: &FlightCondition,
        half: bool,
    ) -> (Vec<Vec3>, Vec<Vec3>) {
        let n = s.mesh().vortex_count();
        let pg = (1.0 - c.mach * c.mach).sqrt();
        let mut result = (vec![Vec3::zeros(); n], vec![Vec3::zeros(); n]);
        for (i, target) in s.mesh().vortices().iter().enumerate() {
            for (j, source) in s.mesh().vortices().iter().enumerate() {
                for image in 0..if half { 2 } else { 1 } {
                    let (left, right) = image_endpoints(source, image);
                    let radius = core(s.mesh(), i, j);
                    let control = horseshoe_velocity_pg_derivative(
                        &target.control_point,
                        &left,
                        &right,
                        pg,
                        radius,
                        true,
                    );
                    let bound = horseshoe_velocity_pg_derivative(
                        &target.bound_vortex_center,
                        &left,
                        &right,
                        pg,
                        radius,
                        image == 1 || i != j,
                    );
                    result.0[i] += (-c.mach / pg) * s.solution().circulation[j] * control;
                    result.1[i] += (-c.mach / pg) * s.solution().circulation[j] * bound;
                }
            }
        }
        result
    }

    /// The baseline geometry reverse: two full Jacobians and a primal velocity
    /// per pair, contracted per output.
    fn reference_geometry(
        s: &VortexLatticeSolver,
        c: &FlightCondition,
        half: bool,
        adjoints: &MatX,
        columns: &[usize],
        seeds: &mut [Seeds],
    ) {
        let pg = (1.0 - c.mach * c.mach).sqrt();
        let u = VortexLatticeSolver::freestream_direction(c);
        let omega = Vec3::new(c.roll_rate, c.pitch_rate, c.yaw_rate);
        for (i, target) in s.mesh().vortices().iter().enumerate() {
            let normal = s.linearized_normals()[i];
            for &k in columns {
                let lambda = adjoints[(i, k)];
                let geometry = &mut seeds[k].geometry;
                if target.affected_by_freestream {
                    geometry.normal[i] -=
                        lambda * (u + omega.cross(&(target.control_point - s.reference_point())));
                    geometry.control_point[i] -= lambda * normal.cross(&omega);
                }
            }
            for (j, source) in s.mesh().vortices().iter().enumerate() {
                for image in 0..if half { 2 } else { 1 } {
                    let (left, right) = image_endpoints(source, image);
                    let radius = core(s.mesh(), i, j);
                    let jac = horseshoe_velocity_jacobians(
                        &target.control_point,
                        &left,
                        &right,
                        pg,
                        radius,
                        true,
                    );
                    let velocity = horseshoe_vortex_velocity(
                        &target.control_point,
                        &left,
                        &right,
                        pg,
                        radius,
                        true,
                    );
                    let load_jac = horseshoe_velocity_jacobians(
                        &target.bound_vortex_center,
                        &left,
                        &right,
                        pg,
                        radius,
                        image == 1 || i != j,
                    );
                    let gamma = s.solution().circulation[j];
                    for &k in columns {
                        let seed = -gamma * adjoints[(i, k)] * normal;
                        let load_seed = gamma * seeds[k].velocity[i];
                        let l =
                            jac.bound_left.tr_mul(&seed) + load_jac.bound_left.tr_mul(&load_seed);
                        let r =
                            jac.bound_right.tr_mul(&seed) + load_jac.bound_right.tr_mul(&load_seed);
                        let geometry = &mut seeds[k].geometry;
                        geometry.bound_center[i] += load_jac.field_point.tr_mul(&load_seed);
                        add_endpoints(geometry, j, image, l, r);
                        geometry.control_point[i] += jac.field_point.tr_mul(&seed);
                        geometry.normal[i] -= adjoints[(i, k)] * gamma * velocity;
                    }
                }
            }
        }
    }

    fn assert_close(actual: &[Vec3], expected: &[Vec3], tolerance: Scalar) {
        let scale = expected.iter().map(|v| v.amax()).fold(0.0, Scalar::max);
        let error = actual
            .iter()
            .zip(expected)
            .map(|(a, b)| (a - b).amax())
            .fold(0.0, Scalar::max);
        assert!(scale > 0.0);
        assert!(
            error <= tolerance * scale,
            "error {error:e}, scale {scale:e}"
        );
    }

    #[test]
    fn mach_velocity_matches_the_baseline_in_both_modes() {
        for half in [false, true] {
            let (s, c) = solved(half, 0.4, 6);
            let mut work = Work::default();
            let (control, bound) = mach_velocity(&s, &c, half, &mut work);
            let (expected_control, expected_bound) = reference_mach(&s, &c, half);
            assert_close(&control, &expected_control, 1e-13);
            assert_close(&bound, &expected_bound, 1e-13);
            let n = s.mesh().vortex_count();
            assert_eq!(work.mach_interactions, 2 * n * n * if half { 2 } else { 1 });
        }
    }

    #[test]
    fn mach_velocity_does_no_work_when_incompressible() {
        for half in [false, true] {
            for mach in [0.0, 0.5 * EPSILON] {
                let (s, c) = solved(half, mach, 4);
                let mut work = Work::default();
                let (control, bound) = mach_velocity(&s, &c, half, &mut work);
                assert!(control.iter().chain(&bound).all(|v| *v == Vec3::zeros()));
                assert_eq!(work.mach_interactions, 0);
            }
        }
    }

    /// Both load-seed branches (separate VJPs for up to three near-field
    /// outputs, Jacobians above), plus a far-field output with no load seed.
    #[test]
    fn geometry_reverse_matches_the_baseline_in_both_modes() {
        let output_sets: [&[Output]; 2] = [
            &[Output::Lift, Output::FarFieldInducedDrag],
            &[
                Output::Lift,
                Output::Drag,
                Output::Pitch,
                Output::NearFieldInducedDrag,
            ],
        ];
        for half in [false, true] {
            let (s, c) = solved(half, 0.3, 6);
            let n = s.mesh().vortex_count();
            for outputs in output_sets {
                let needs = Needs {
                    geometry: true,
                    alpha: false,
                    beta: false,
                    reference: false,
                };
                let make_seeds = || -> Vec<Seeds> {
                    outputs
                        .iter()
                        .map(|&o| {
                            if o.far() {
                                let (gamma, geometry) = crate::analysis::trefftz::reverse(
                                    &s,
                                    half,
                                    o,
                                    true,
                                    &mut Work::default(),
                                );
                                Seeds {
                                    gamma,
                                    geometry,
                                    velocity: vec![Vec3::zeros(); n],
                                    alpha: 0.0,
                                    beta: 0.0,
                                    reference: [0.0; 3],
                                }
                            } else {
                                near_seeds(&s, &c, half, o, needs, &mut Work::default())
                            }
                        })
                        .collect()
                };
                let adjoints = MatX::from_fn(n, outputs.len(), |i, k| {
                    0.3 + 0.1 * ((i * 7 + k * 3) as Scalar).sin()
                });
                let columns: Vec<_> = (0..outputs.len()).collect();
                let mut actual = make_seeds();
                let mut expected = make_seeds();
                let mut work = Work::default();
                geometry_reverse(&s, &c, half, &adjoints, &columns, &mut actual, &mut work);
                reference_geometry(&s, &c, half, &adjoints, &columns, &mut expected);
                assert_eq!(
                    work.geometry_interactions,
                    2 * n * n * if half { 2 } else { 1 }
                );
                for (a, e) in actual.iter().zip(&expected) {
                    let (a, e) = (&a.geometry, &e.geometry);
                    assert_close(&a.bound_left, &e.bound_left, 1e-13);
                    assert_close(&a.bound_right, &e.bound_right, 1e-13);
                    assert_close(&a.control_point, &e.control_point, 1e-13);
                    assert_close(&a.normal, &e.normal, 1e-13);
                    if e.bound_center.iter().any(|v| *v != Vec3::zeros()) {
                        assert_close(&a.bound_center, &e.bound_center, 1e-13);
                    } else {
                        assert!(a.bound_center.iter().all(|v| *v == Vec3::zeros()));
                    }
                }
            }
        }
    }
}
