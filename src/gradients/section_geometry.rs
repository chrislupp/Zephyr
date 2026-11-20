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

//! Forward-mode lattice tangents for aircraft section design variables.
//!
//! This is the analytical counterpart of `direct_mesh_node_direction` for the
//! section/AVL geometry path. It differentiates `VortexLattice::generate_surface`,
//! `recalculate_normals`, and `duplicate_surface` in closed form, producing the
//! same [`LatticeDirections`] seed the adjoint context already consumes. The
//! entire chain downstream of that seed is unchanged.
//!
//! Several properties of the mesher make this far smaller than it first looks,
//! and each is load-bearing:
//!
//! - **The spanwise interpolation weights are constants.** A strip's `f1`, `f2`
//!   and `fc` are `(p - arc[k]) / (arc[k+1] - arc[k])` where `p` was built as
//!   `arc[k] + length * fractions[..]`, so the arc lengths cancel identically
//!   and the weights reduce to the `spacing()` fractions. Leading-edge y and z
//!   therefore do *not* propagate globally through arc length. (At an interval
//!   seam the cancellation leaves a residue of order machine epsilon, bounded
//!   and orders below any verification tolerance — not an omission.)
//! - **Twist does not move panel vertices.** Panel geometry is built purely as
//!   `leading_edge + fraction * chord` along x, so `d(bound_vortex)/d(twist)`
//!   and `d(control_point)/d(twist)` are *exactly* zero. Twist reaches the
//!   solution only through the normal.
//! - **Camber contributes nothing.** `camber_slope_at_control_point` is
//!   evaluated at `x_control[i]`, which no design variable moves, and the
//!   camber spline is stored in normalized chord fractions and never rebuilt.
//! - **`x_vortex[0]` is the only movable chordwise fraction**, through `claf`,
//!   and only when the bracketing sections disagree on `cl_alpha_factor`.
//!
//! A section perturbation touches the two intervals it bounds — as the right
//! end of the preceding one, then as the left end of the following one. Those
//! intervals own **disjoint** strip ranges, so each writes its own entries and
//! plain assignment is correct. Do not "fix" this to `+=`: it is a no-op today
//! and becomes a silent double-count the moment anything else writes here.

use crate::core::{DEGREES_TO_RADIANS, EPSILON, Scalar, Vec3};
use crate::error::{Result, ZephyrError};
use crate::geometry::{Aircraft, Surface};
use crate::solver::{VortexLattice, chordwise_spacing};

use super::GeometryVariable;
use super::geometry_adjoint::{LatticeDirections, normalized_direction, reflected_direction};

/// Tangent of one section's mesher-space quantities.
///
/// Mesher space is the section after the surface's `scale`, `translation` and
/// `angle_offset` have been applied, which is what the strip loop interpolates.
#[derive(Clone, Copy, Default)]
struct SectionTangent {
    /// `d(leading_edge * scale + translation)`
    leading_edge: Vec3,
    /// `d(scale.x * chord)`
    chord: Scalar,
    /// `d(twist_degrees * pi/180 + angle_offset)`, radians
    twist: Scalar,
}

/// Exact chord-weighted interpolated-twist JVP.
///
/// The mesher sets a strip's twist to `atan2(S, K)` with
///
/// ```text
/// S = (1-f)*c_l*sin(t_l) + f*c_r*sin(t_r)
/// K = (1-f)*c_l*cos(t_l) + f*c_r*cos(t_r)
/// ```
///
/// so twist is interpolated as a *chord-weighted vector*, not linearly in the
/// angle. The derivative is written as `(K dS - S dK) / (S² + K²)` rather than
/// as four expanded partials: algebraically identical, half the terms, and no
/// opportunity to transpose a subscript.
///
/// This is deliberately **not** the small-angle form
/// `dTheta = (1-f)*dt_l + f*dt_r` used by the v0.2.0 C++ `PropagateTwistGradient`.
/// That form is exact only when `t_l == t_r` *and* `c_l == c_r`, so it is wrong
/// for any tapered or twisted wing — on a 1.6/0.5 taper at 6°/-5° it is off by
/// more than 50%.
///
/// Angles are radians here; the per-degree convention of `Section::twist_degrees`
/// is applied by the caller when it builds the seed.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn interpolated_twist_jvp(
    fraction: Scalar,
    chord_left: Scalar,
    chord_right: Scalar,
    twist_left: Scalar,
    twist_right: Scalar,
    chord_left_direction: Scalar,
    chord_right_direction: Scalar,
    twist_left_direction: Scalar,
    twist_right_direction: Scalar,
) -> Scalar {
    let (sin_left, cos_left) = twist_left.sin_cos();
    let (sin_right, cos_right) = twist_right.sin_cos();
    let weight_left = 1.0 - fraction;
    let weight_right = fraction;

    let s = weight_left * chord_left * sin_left + weight_right * chord_right * sin_right;
    let k = weight_left * chord_left * cos_left + weight_right * chord_right * cos_right;
    let radius = s * s + k * k;
    if radius <= EPSILON {
        return 0.0;
    }

    let d_s = weight_left
        * (chord_left_direction * sin_left + chord_left * cos_left * twist_left_direction)
        + weight_right
            * (chord_right_direction * sin_right + chord_right * cos_right * twist_right_direction);
    let d_k = weight_left
        * (chord_left_direction * cos_left - chord_left * sin_left * twist_left_direction)
        + weight_right
            * (chord_right_direction * cos_right - chord_right * sin_right * twist_right_direction);

    (k * d_s - s * d_k) / radius
}

/// JVP of `VortexLattice::recalculate_normals` for one panel.
///
/// Mirrors the primal's guards exactly: the unnormalized-`bound` fallback when
/// the panel has zero width, and the Trefftz-normal fallback when the cross
/// product degenerates. `normalized_direction` is applied to the *unnormalized*
/// cross product, matching the primal's `normal.normalize()`.
///
/// `camber_slope` enters only through `atan`, and its derivative is zero (see
/// the module docs), so the incidence tangent is the twist tangent alone.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn panel_normal_jvp(
    bound_left: &Vec3,
    bound_right: &Vec3,
    bound_left_direction: &Vec3,
    bound_right_direction: &Vec3,
    twist_angle: Scalar,
    twist_angle_direction: Scalar,
    camber_slope: Scalar,
    trefftz_normal: (Scalar, Scalar),
    trefftz_normal_direction: (Scalar, Scalar),
) -> Vec3 {
    let bound = bound_right - bound_left;
    let bound_direction = bound_right_direction - bound_left_direction;
    let (unit_bound, unit_bound_direction) = if bound.norm() > EPSILON {
        (
            bound.normalize(),
            normalized_direction(&bound, &bound_direction),
        )
    } else {
        // The primal returns the unnormalized vector here, so the tangent is
        // the unnormalized tangent. Only reachable for a zero-width panel,
        // which makes the system singular anyway.
        (bound, bound_direction)
    };

    let angle = twist_angle - camber_slope.atan();
    let (sin_angle, cos_angle) = angle.sin_cos();
    let (normal_y, normal_z) = trefftz_normal;
    let (d_normal_y, d_normal_z) = trefftz_normal_direction;

    let chord = Vec3::new(cos_angle, -sin_angle * normal_y, -sin_angle * normal_z);
    let chord_direction = Vec3::new(
        -sin_angle * twist_angle_direction,
        -cos_angle * twist_angle_direction * normal_y - sin_angle * d_normal_y,
        -cos_angle * twist_angle_direction * normal_z - sin_angle * d_normal_z,
    );

    let normal = chord.cross(&unit_bound);
    let normal_direction = chord_direction.cross(&unit_bound) + chord.cross(&unit_bound_direction);
    if normal.norm() > EPSILON {
        normalized_direction(&normal, &normal_direction)
    } else {
        Vec3::new(0.0, d_normal_y, d_normal_z)
    }
}

/// Forward-mode lattice tangent for a unit perturbation of one section
/// design variable.
///
/// Chord and leading-edge seeds are per unit length; twist is **per degree**,
/// matching `Section::twist_degrees` and the crate's documented unit split.
pub(crate) fn section_direction(
    aircraft: &Aircraft,
    lattice: &VortexLattice,
    surface_index: usize,
    section_index: usize,
    variable: GeometryVariable,
) -> Result<LatticeDirections> {
    let surface = aircraft
        .surfaces()
        .get(surface_index)
        .ok_or_else(|| ZephyrError::InvalidGeometry("section index is out of range".to_owned()))?;
    let sections = surface.sections();
    if section_index >= sections.len() {
        return Err(ZephyrError::InvalidGeometry(
            "section index is out of range".to_owned(),
        ));
    }

    let layout = surface.spanwise_layout().ok_or_else(|| {
        ZephyrError::InvalidGeometry(format!("surface '{}' has zero span", surface.name()))
    })?;

    // Each preceding surface occupies one lattice slot, or two when mirrored.
    let surface_slot: usize = aircraft.surfaces()[..surface_index]
        .iter()
        .map(|surface| 1 + usize::from(surface.y_duplicate().is_some()))
        .sum();
    let info = lattice
        .surface_info()
        .get(surface_slot)
        .cloned()
        .ok_or_else(|| ZephyrError::InvalidGeometry("surface lattice is missing".to_owned()))?;

    let scale = *surface.scale();
    let seed = match variable {
        GeometryVariable::Chord => SectionTangent {
            chord: scale.x,
            ..SectionTangent::default()
        },
        GeometryVariable::Twist => SectionTangent {
            twist: DEGREES_TO_RADIANS,
            ..SectionTangent::default()
        },
        GeometryVariable::LeadingEdgeX => SectionTangent {
            leading_edge: Vec3::new(scale.x, 0.0, 0.0),
            ..SectionTangent::default()
        },
        GeometryVariable::LeadingEdgeY => SectionTangent {
            leading_edge: Vec3::new(0.0, scale.y, 0.0),
            ..SectionTangent::default()
        },
        GeometryVariable::LeadingEdgeZ => SectionTangent {
            leading_edge: Vec3::new(0.0, 0.0, scale.z),
            ..SectionTangent::default()
        },
    };

    let mut result = LatticeDirections::zeros(lattice);
    let mut strip_offset = 0usize;
    for interval in 0..sections.len() - 1 {
        let count = layout.per_interval[interval];
        // The perturbed section bounds at most these two intervals, and the
        // two write disjoint strip ranges.
        let left_seed = if interval == section_index {
            seed
        } else {
            SectionTangent::default()
        };
        let right_seed = if interval + 1 == section_index {
            seed
        } else {
            SectionTangent::default()
        };
        if interval == section_index || interval + 1 == section_index {
            write_interval(
                &mut result,
                lattice,
                surface,
                &info,
                &layout,
                interval,
                strip_offset,
                &left_seed,
                &right_seed,
            );
        }
        strip_offset += count;
    }

    if surface.y_duplicate().is_some() {
        mirror_directions(&mut result, lattice, &info, surface_slot)?;
    }
    Ok(result)
}

/// Writes the tangents for every strip of one section interval.
#[allow(clippy::too_many_arguments)]
fn write_interval(
    result: &mut LatticeDirections,
    lattice: &VortexLattice,
    surface: &Surface,
    info: &crate::solver::SurfaceMeshInfo,
    layout: &crate::geometry::SpanwiseLayout,
    interval: usize,
    first_strip_offset: usize,
    left_seed: &SectionTangent,
    right_seed: &SectionTangent,
) {
    let sections = surface.sections();
    let left = &sections[interval];
    let right = &sections[interval + 1];
    let scale = *surface.scale();
    let angle_offset = surface.angle_offset_degrees() * DEGREES_TO_RADIANS;

    let chord_left = scale.x * left.chord();
    let chord_right = scale.x * right.chord();
    let twist_left = left.twist_degrees() * DEGREES_TO_RADIANS + angle_offset;
    let twist_right = right.twist_degrees() * DEGREES_TO_RADIANS + angle_offset;
    let cl_alpha_left = left.cl_alpha_factor();
    let cl_alpha_right = right.cl_alpha_factor();

    let nx = info.num_chordwise;
    let start = layout.arc[interval];
    let span = layout.arc[interval + 1] - start;

    for panel in 0..layout.per_interval[interval] {
        let strip_offset = first_strip_offset + panel;
        let strip_index = info.first_strip_index + strip_offset;
        let strip = &lattice.strips()[strip_index];

        // Constants: the arc lengths cancel out of these ratios exactly.
        let f1 = (layout.panel_edges[strip_offset] - start) / span;
        let f2 = (layout.panel_edges[strip_offset + 1] - start) / span;
        let fc = (layout.controls[strip_offset] - start) / span;

        let interpolate =
            |weight: Scalar, left: Vec3, right: Vec3| (1.0 - weight) * left + weight * right;
        let le_left_direction = interpolate(f1, left_seed.leading_edge, right_seed.leading_edge);
        let le_right_direction = interpolate(f2, left_seed.leading_edge, right_seed.leading_edge);
        let le_center_direction = interpolate(fc, left_seed.leading_edge, right_seed.leading_edge);
        let chord_left_direction = (1.0 - f1) * left_seed.chord + f1 * right_seed.chord;
        let chord_right_direction = (1.0 - f2) * left_seed.chord + f2 * right_seed.chord;
        let chord_center_direction = (1.0 - fc) * left_seed.chord + fc * right_seed.chord;

        result.strip_chord[strip_index] = chord_center_direction;

        // Trefftz normal: a 2-D normalize of the strip's y-z leading-edge span.
        let delta_y = strip.leading_edge_right.y - strip.leading_edge_left.y;
        let delta_z = strip.leading_edge_right.z - strip.leading_edge_left.z;
        let d_delta_y = le_right_direction.y - le_left_direction.y;
        let d_delta_z = le_right_direction.z - le_left_direction.z;
        let magnitude = delta_y.hypot(delta_z);
        let trefftz_direction = if magnitude > EPSILON {
            let d_magnitude = (delta_y * d_delta_y + delta_z * d_delta_z) / magnitude;
            (
                -d_delta_z / magnitude + delta_z * d_magnitude / (magnitude * magnitude),
                d_delta_y / magnitude - delta_y * d_magnitude / (magnitude * magnitude),
            )
        } else {
            (0.0, 0.0)
        };

        // The section-level tangents go in unweighted: `interpolated_twist_jvp`
        // applies the `fc` weighting itself, exactly as the primal does.
        let twist_direction = interpolated_twist_jvp(
            fc,
            chord_left,
            chord_right,
            twist_left,
            twist_right,
            left_seed.chord,
            right_seed.chord,
            left_seed.twist,
            right_seed.twist,
        );

        // `claf` moves `x_vortex[0]` only when the bracketing sections disagree
        // on `cl_alpha_factor`; otherwise it is identically one and this term
        // vanishes.
        let claf = (1.0 - fc) * chord_left / strip.chord_center * cl_alpha_left
            + fc * chord_right / strip.chord_center * cl_alpha_right;
        let (x_vortex, x_control) = chordwise_spacing(nx, surface.chordwise_spacing(), claf);
        let edges = crate::utils::spacing(nx + 1, surface.chordwise_spacing());
        let raw = 0.25 / claf;
        let leading_vortex_direction = if claf != 1.0 && raw > 0.1 && raw < 0.4 {
            // `claf` weights the *section* chords against the strip's centre
            // chord, not the strip's own left/right chords. Using the strip
            // tangents here makes `d_claf` spuriously non-zero even when both
            // sections carry the default unit factor, where `claf` is
            // identically one and this whole term must vanish.
            let d_claf = ((1.0 - fc) * left_seed.chord * cl_alpha_left
                + fc * right_seed.chord * cl_alpha_right
                - claf * chord_center_direction)
                / strip.chord_center;
            -0.25 / (claf * claf) * d_claf * (edges[1] - edges[0])
        } else {
            0.0
        };

        for chord_index in 0..nx {
            let vortex_index = info.first_vortex_index + strip_offset * nx + chord_index;
            let d_x_vortex = if chord_index == 0 {
                leading_vortex_direction
            } else {
                0.0
            };

            let bound_left_direction = Vec3::new(
                le_left_direction.x
                    + x_vortex[chord_index] * chord_left_direction
                    + d_x_vortex * strip.chord_left,
                le_left_direction.y,
                le_left_direction.z,
            );
            let bound_right_direction = Vec3::new(
                le_right_direction.x
                    + x_vortex[chord_index] * chord_right_direction
                    + d_x_vortex * strip.chord_right,
                le_right_direction.y,
                le_right_direction.z,
            );
            result.bound_left[vortex_index] = bound_left_direction;
            result.bound_right[vortex_index] = bound_right_direction;
            result.control_point[vortex_index] = Vec3::new(
                le_center_direction.x + x_control[chord_index] * chord_center_direction,
                le_center_direction.y,
                le_center_direction.z,
            );
            result.local_chord[vortex_index] = chord_center_direction;

            let vortex = &lattice.vortices()[vortex_index];
            result.normal[vortex_index] = panel_normal_jvp(
                &vortex.bound_vortex_left,
                &vortex.bound_vortex_right,
                &bound_left_direction,
                &bound_right_direction,
                strip.twist_angle,
                twist_direction,
                vortex.camber_slope_at_control_point,
                (strip.trefftz_normal_y, strip.trefftz_normal_z),
                trefftz_direction,
            );
        }
    }
}

/// Reflects an original surface's tangents onto its mirrored copy.
///
/// `recalculate_normals` runs *after* `duplicate_surface`, so the mirrored
/// normal is recomputed rather than reflected — but the reflection identity
/// still holds exactly. With `R = diag(1, -1, 1)`:
///
/// - `bound_m = reflect(bl_o) - reflect(br_o) = -R * bound_o`, so
///   `|bound_m| = |bound_o|` and the zero-width guard agrees;
/// - the Trefftz normal picks up `tny_m = -tny_o`, `tnz_m = tnz_o`, and the
///   twist and camber are cloned unchanged, so `chord_m = R * chord_o`;
/// - for an orthogonal reflection `(Ru) x (Rv) = det(R) * R(u x v)`, and
///   `det(R) = -1`, so `normal_m = (R c) x (-R b) = R (c x b) = R * normal_o`,
///   with `|normal_m| = |normal_o|` so the degenerate branch agrees too.
///
/// Hence the tangent mirrors with a plain reflection. The bound endpoints
/// **swap** — `duplicate_surface` swaps them to keep circulation orientation
/// consistent — while the normal, being a single vector, does not.
fn mirror_directions(
    result: &mut LatticeDirections,
    lattice: &VortexLattice,
    info: &crate::solver::SurfaceMeshInfo,
    surface_slot: usize,
) -> Result<()> {
    let duplicate = lattice
        .surface_info()
        .get(surface_slot + 1)
        .cloned()
        .ok_or_else(|| {
            ZephyrError::InvalidGeometry("mirrored surface lattice is missing".to_owned())
        })?;
    for strip_offset in 0..info.num_strips {
        let original_strip = info.first_strip_index + strip_offset;
        let duplicate_strip = duplicate.first_strip_index + strip_offset;
        result.strip_chord[duplicate_strip] = result.strip_chord[original_strip];
        for chord_offset in 0..info.num_chordwise {
            let original =
                info.first_vortex_index + strip_offset * info.num_chordwise + chord_offset;
            let mirrored =
                duplicate.first_vortex_index + strip_offset * info.num_chordwise + chord_offset;
            result.bound_left[mirrored] = reflected_direction(result.bound_right[original]);
            result.bound_right[mirrored] = reflected_direction(result.bound_left[original]);
            result.control_point[mirrored] = reflected_direction(result.control_point[original]);
            result.normal[mirrored] = reflected_direction(result.normal[original]);
            result.local_chord[mirrored] = result.local_chord[original];
        }
    }
    Ok(())
}
