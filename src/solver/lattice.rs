//! Aircraft and structured-grid discretization.

use crate::core::{DEGREES_TO_RADIANS, EPSILON, Scalar, SpacingType, Vec3};
use crate::error::{PegasusError, Result};
use crate::geometry::{Aircraft, MeshInput, Surface, SurfaceMesh};
use crate::utils::spacing;

use super::{ControlSurfaceInfluence, HorseshoeVortex, Strip, SurfaceMeshInfo};

/// Complete panel lattice used by the solver.
#[derive(Clone, Debug)]
pub struct VortexLattice {
    vortices: Vec<HorseshoeVortex>,
    strips: Vec<Strip>,
    surface_info: Vec<SurfaceMeshInfo>,
    control_influences: Vec<Vec<ControlSurfaceInfluence>>,
    control_surface_names: Vec<String>,
    spanwise_axis_fraction: Scalar,
    mesh_input: Option<MeshInput>,
}

impl Default for VortexLattice {
    fn default() -> Self {
        Self {
            vortices: Vec::new(),
            strips: Vec::new(),
            surface_info: Vec::new(),
            control_influences: Vec::new(),
            control_surface_names: Vec::new(),
            spanwise_axis_fraction: 0.25,
            mesh_input: None,
        }
    }
}

impl VortexLattice {
    /// Whether the direct mesh requests physically consistent force locations.
    #[must_use]
    pub fn consistent_forces(&self) -> bool {
        self.mesh_input
            .as_ref()
            .is_some_and(|input| input.consistent_forces)
    }

    /// Location where Kutta-Joukowski forces act for the selected convention.
    #[must_use]
    pub fn force_point(&self, vortex: &HorseshoeVortex) -> Vec3 {
        if self.consistent_forces() {
            vortex.bound_vortex_center
        } else {
            vortex.control_point
        }
    }

    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    pub fn clear(&mut self) {
        self.vortices.clear();
        self.strips.clear();
        self.surface_info.clear();
        self.control_influences.clear();
        self.control_surface_names.clear();
        self.mesh_input = None;
    }

    pub fn generate_from_aircraft(&mut self, aircraft: &Aircraft) -> Result<()> {
        self.clear();
        self.control_surface_names = aircraft.control_surface_names();
        for surface in aircraft.surfaces() {
            self.generate_surface(surface)?;
            if let Some(plane) = surface.y_duplicate() {
                self.duplicate_surface(self.surface_info.len() - 1, plane);
            }
        }
        self.recalculate_normals();
        Ok(())
    }

    fn generate_surface(&mut self, surface: &Surface) -> Result<()> {
        let sections = surface.sections();
        if sections.len() < 2 {
            return Err(PegasusError::InvalidGeometry(format!(
                "surface '{}' needs at least two sections",
                surface.name()
            )));
        }
        validate_sections(surface)?;
        let first_vortex = self.vortices.len();
        let first_strip = self.strips.len();
        let surface_slot = self.surface_info.len();
        let nx = surface.chordwise_panel_count();
        if nx == 0 {
            return Err(PegasusError::InvalidGeometry(format!(
                "surface '{}' needs at least one chordwise panel",
                surface.name()
            )));
        }

        // The layout is computed once, by `Surface`, and shared with everything
        // that must agree with the lattice this produces. Recomputing the arc
        // accumulation and the per-interval rounding in a second place is how
        // the aircraft counting helpers drifted out of step with the mesher.
        let layout = surface.spanwise_layout().ok_or_else(|| {
            PegasusError::InvalidGeometry(format!("surface '{}' has zero span", surface.name()))
        })?;
        let arc = &layout.arc;
        let per_interval = &layout.per_interval;
        let panel_edges = &layout.panel_edges;
        let controls = &layout.controls;
        let ny = layout.strip_count();

        let scale = surface.scale();
        let translation = surface.translation();
        let angle_offset = surface.angle_offset_degrees() * DEGREES_TO_RADIANS;
        let x_edges = spacing(nx + 1, surface.chordwise_spacing());
        let mut strip_index = 0;
        let mut wetted_area = 0.0;
        let mut total_width = 0.0;

        for section_index in 0..sections.len() - 1 {
            let left = &sections[section_index];
            let right = &sections[section_index + 1];
            let le_left = left.leading_edge().component_mul(scale) + translation;
            let le_right = right.leading_edge().component_mul(scale) + translation;
            let chord_left = scale.x * left.chord();
            let chord_right = scale.x * right.chord();
            let twist_left = left.twist_degrees() * DEGREES_TO_RADIANS + angle_offset;
            let twist_right = right.twist_degrees() * DEGREES_TO_RADIANS + angle_offset;
            let section_width = (le_right - le_left).yz().norm();
            let sweep_le = if section_width > 0.0 {
                (le_right.x - le_left.x) / section_width
            } else {
                0.0
            };
            let sweep_te = if section_width > 0.0 {
                (le_right.x + chord_right - le_left.x - chord_left) / section_width
            } else {
                0.0
            };
            let sin_twist_left = chord_left * twist_left.sin();
            let cos_twist_left = chord_left * twist_left.cos();
            let sin_twist_right = chord_right * twist_right.sin();
            let cos_twist_right = chord_right * twist_right.cos();
            let start = arc[section_index];
            let end = arc[section_index + 1];

            for _ in 0..per_interval[section_index] {
                let p1 = panel_edges[strip_index];
                let p2 = panel_edges[strip_index + 1];
                let pc = controls[strip_index];
                let f1 = (p1 - start) / (end - start);
                let f2 = (p2 - start) / (end - start);
                let fc = (pc - start) / (end - start);
                let mut strip = Strip {
                    leading_edge_left: (1.0 - f1) * le_left + f1 * le_right,
                    leading_edge_right: (1.0 - f2) * le_left + f2 * le_right,
                    leading_edge_center: (1.0 - fc) * le_left + fc * le_right,
                    chord_left: (1.0 - f1) * chord_left + f1 * chord_right,
                    chord_right: (1.0 - f2) * chord_left + f2 * chord_right,
                    chord_center: (1.0 - fc) * chord_left + fc * chord_right,
                    width: (f2 - f1).abs() * section_width,
                    leading_edge_sweep_slope: sweep_le,
                    trailing_edge_sweep_slope: sweep_te,
                    twist_angle: ((1.0 - fc) * sin_twist_left + fc * sin_twist_right)
                        .atan2((1.0 - fc) * cos_twist_left + fc * cos_twist_right),
                    twist_angle_left: ((1.0 - f1) * sin_twist_left + f1 * sin_twist_right)
                        .atan2((1.0 - f1) * cos_twist_left + f1 * cos_twist_right),
                    twist_angle_right: ((1.0 - f2) * sin_twist_left + f2 * sin_twist_right)
                        .atan2((1.0 - f2) * cos_twist_left + f2 * cos_twist_right),
                    first_vortex_index: self.vortices.len(),
                    num_vortices: nx,
                    surface_index: surface_slot,
                    ..Strip::default()
                };
                let dy = strip.leading_edge_right.y - strip.leading_edge_left.y;
                let dz = strip.leading_edge_right.z - strip.leading_edge_left.z;
                let yz_mag = dy.hypot(dz);
                if yz_mag > EPSILON {
                    strip.trefftz_normal_y = -dz / yz_mag;
                    strip.trefftz_normal_z = dy / yz_mag;
                    strip.spanwise_unit_vector = Vec3::new(0.0, dy / yz_mag, dz / yz_mag);
                }
                strip.reference_point = strip.leading_edge_center;
                strip.reference_point.x += self.spanwise_axis_fraction * strip.chord_center;
                for polar_index in 0..6 {
                    strip.viscous_polar[polar_index] = (1.0 - fc)
                        * left.viscous_polar()[polar_index]
                        + fc * right.viscous_polar()[polar_index];
                }
                strip.has_viscous_data = strip
                    .viscous_polar
                    .iter()
                    .any(|value| value.abs() > EPSILON);
                if !strip.has_viscous_data && surface.has_viscous_drag() {
                    strip.viscous_polar = *surface.cl_cd_polar();
                    strip.has_viscous_data = true;
                }
                wetted_area += strip.width * strip.chord_center;
                total_width += strip.width;

                let claf = (1.0 - fc) * chord_left / strip.chord_center * left.cl_alpha_factor()
                    + fc * chord_right / strip.chord_center * right.cl_alpha_factor();
                let (x_vortex, x_control) =
                    chordwise_spacing(nx, surface.chordwise_spacing(), claf);
                self.strips.push(strip.clone());

                for chord_index in 0..nx {
                    let mut vortex = HorseshoeVortex {
                        bound_vortex_left: Vec3::new(
                            strip.leading_edge_left.x + x_vortex[chord_index] * strip.chord_left,
                            strip.leading_edge_left.y,
                            strip.leading_edge_left.z,
                        ),
                        bound_vortex_right: Vec3::new(
                            strip.leading_edge_right.x + x_vortex[chord_index] * strip.chord_right,
                            strip.leading_edge_right.y,
                            strip.leading_edge_right.z,
                        ),
                        bound_vortex_center: Vec3::new(
                            strip.leading_edge_center.x
                                + x_vortex[chord_index] * strip.chord_center,
                            strip.leading_edge_center.y,
                            strip.leading_edge_center.z,
                        ),
                        control_point: Vec3::new(
                            strip.leading_edge_center.x
                                + x_control[chord_index] * strip.chord_center,
                            strip.leading_edge_center.y,
                            strip.leading_edge_center.z,
                        ),
                        chordwise_width: (x_edges[chord_index + 1] - x_edges[chord_index])
                            * strip.chord_center,
                        local_chord: strip.chord_center,
                        strip_index: self.strips.len() - 1,
                        surface_index: surface_slot,
                        component_index: surface.component_index(),
                        affected_by_freestream: surface.affected_by_alpha_beta(),
                        enforce_flow_tangency: surface.has_wake(),
                        ..HorseshoeVortex::default()
                    };
                    if left.has_airfoil() {
                        vortex.camber_slope_at_control_point +=
                            (1.0 - fc) * left.airfoil().camber_slope(x_control[chord_index]);
                        vortex.camber_slope_at_vortex +=
                            (1.0 - fc) * left.airfoil().camber_slope(x_vortex[chord_index]);
                    }
                    if right.has_airfoil() {
                        vortex.camber_slope_at_control_point +=
                            fc * right.airfoil().camber_slope(x_control[chord_index]);
                        vortex.camber_slope_at_vortex +=
                            fc * right.airfoil().camber_slope(x_vortex[chord_index]);
                    }
                    self.vortices.push(vortex);

                    let mut influences = Vec::new();
                    for control_left in left.control_surfaces() {
                        let Some(control_right) = right
                            .control_surfaces()
                            .iter()
                            .find(|candidate| candidate.name == control_left.name)
                        else {
                            continue;
                        };
                        let hinge_left = control_left.hinge_x_fraction * chord_left;
                        let hinge_right = control_right.hinge_x_fraction * chord_right;
                        let hinge = (1.0 - fc) * hinge_left + fc * hinge_right;
                        let (control_le, control_te) = if hinge >= 0.0 {
                            (hinge / strip.chord_center, 1.0)
                        } else {
                            (0.0, -hinge / strip.chord_center)
                        };
                        let panel_dx = x_edges[chord_index + 1] - x_edges[chord_index];
                        let frac_le =
                            ((control_le - x_edges[chord_index]) / panel_dx).clamp(0.0, 1.0);
                        let frac_te =
                            ((control_te - x_edges[chord_index]) / panel_dx).clamp(0.0, 1.0);
                        let gain = ((1.0 - fc) * control_left.gain + fc * control_right.gain)
                            * (frac_te - frac_le);
                        if gain.abs() <= EPSILON {
                            continue;
                        }
                        let hinge_vector = (1.0 - fc) * control_left.hinge_vector
                            + fc * control_right.hinge_vector;
                        let hinge_axis = if hinge_vector.norm() > EPSILON {
                            hinge_vector.normalize()
                        } else {
                            strip.spanwise_unit_vector
                        };
                        let mut hinge_point = strip.leading_edge_center;
                        hinge_point.x += hinge.abs();
                        influences.push(ControlSurfaceInfluence {
                            name: control_left.name.clone(),
                            gain,
                            hinge_axis,
                            hinge_point,
                            reflection_sign: control_left.sign,
                        });
                    }
                    self.control_influences.push(influences);
                }
                strip_index += 1;
            }
        }
        self.surface_info.push(SurfaceMeshInfo {
            name: surface.name().into(),
            first_vortex_index: first_vortex,
            num_vortices: self.vortices.len() - first_vortex,
            first_strip_index: first_strip,
            num_strips: self.strips.len() - first_strip,
            num_chordwise: nx,
            num_spanwise: ny,
            wetted_area,
            average_chord: if total_width > 0.0 {
                wetted_area / total_width
            } else {
                0.0
            },
            is_y_duplicated: false,
            // Not a mirror, so it names its own slot.
            original_surface_index: self.surface_info.len(),
        });
        Ok(())
    }

    pub fn generate_from_mesh(&mut self, input: &MeshInput) -> Result<()> {
        input.validate()?;
        self.clear();
        self.mesh_input = Some(input.clone());
        self.control_surface_names = input
            .surfaces
            .iter()
            .flat_map(|surface| surface.controls.iter().map(|control| control.name.clone()))
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .collect();
        for surface in &input.surfaces {
            self.generate_mesh_surface(surface)?;
            if surface.y_duplicated {
                self.duplicate_surface(self.surface_info.len() - 1, surface.y_duplicate_plane);
            }
        }
        Ok(())
    }

    fn generate_mesh_surface(&mut self, mesh: &SurfaceMesh) -> Result<()> {
        mesh.validate()?;
        let nx = mesh.nx;
        let ny = mesh.ny;
        let surface_slot = self.surface_info.len();
        let first_vortex = self.vortices.len();
        let first_strip = self.strips.len();
        let node = |i, j| *mesh.node(i, j).expect("validated mesh index");

        for j in 0..ny {
            let le_left = node(0, j);
            let le_right = node(0, j + 1);
            let le_center = 0.5 * (le_left + le_right);
            let te_left = node(nx, j);
            let te_right = node(nx, j + 1);
            let chord_left = (te_left - le_left).norm();
            let chord_right = (te_right - le_right).norm();
            let chord_center = 0.5 * (chord_left + chord_right);
            let qc_left = 0.75 * node(0, j) + 0.25 * node(1, j);
            let qc_right = 0.75 * node(0, j + 1) + 0.25 * node(1, j + 1);
            let dy = le_right.y - le_left.y;
            let dz = le_right.z - le_left.z;
            let yz_mag = dy.hypot(dz);
            let mut strip = Strip {
                leading_edge_left: le_left,
                leading_edge_right: le_right,
                leading_edge_center: le_center,
                chord_left,
                chord_right,
                chord_center,
                width: (qc_right - qc_left).norm(),
                twist_angle: (0.5 * (te_left + te_right) - le_center)
                    .z
                    .atan2((0.5 * (te_left + te_right) - le_center).x),
                twist_angle_left: (te_left - le_left).z.atan2((te_left - le_left).x),
                twist_angle_right: (te_right - le_right).z.atan2((te_right - le_right).x),
                first_vortex_index: self.vortices.len(),
                num_vortices: nx,
                surface_index: surface_slot,
                ..Strip::default()
            };
            if yz_mag > EPSILON {
                strip.leading_edge_sweep_slope = (le_right.x - le_left.x) / yz_mag;
                strip.trailing_edge_sweep_slope = (te_right.x - te_left.x) / yz_mag;
                strip.trefftz_normal_y = dy / yz_mag;
                strip.trefftz_normal_z = dz / yz_mag;
            }
            strip.reference_point = strip.leading_edge_center;
            strip.reference_point.x += 0.25 * strip.chord_center;

            for i in 0..nx {
                let n00 = node(i, j);
                let n10 = node(i + 1, j);
                let n01 = node(i, j + 1);
                let n11 = node(i + 1, j + 1);
                let left = 0.75 * n00 + 0.25 * n10;
                let right = 0.75 * n01 + 0.25 * n11;
                let cross = (n11 - n00).cross(&(n10 - n01));
                let normal = if cross.norm() > EPSILON {
                    cross.normalize()
                } else {
                    Vec3::z()
                };
                self.vortices.push(HorseshoeVortex {
                    bound_vortex_left: left,
                    bound_vortex_right: right,
                    bound_vortex_center: 0.5 * (left + right),
                    control_point: 0.5 * (0.25 * n00 + 0.75 * n10 + 0.25 * n01 + 0.75 * n11),
                    normal_at_control_point: normal,
                    normal_at_vortex: normal,
                    chordwise_width: 0.5 * ((n10 - n00).norm() + (n11 - n01).norm()),
                    local_chord: chord_center,
                    strip_index: self.strips.len(),
                    surface_index: surface_slot,
                    component_index: mesh.component_index,
                    affected_by_freestream: mesh.affected_by_freestream,
                    enforce_flow_tangency: true,
                    ..HorseshoeVortex::default()
                });
                self.control_influences.push(
                    mesh.controls
                        .iter()
                        .filter_map(|control| {
                            let gain = control.panel_gain(i, j);
                            (gain.abs() > EPSILON).then(|| ControlSurfaceInfluence {
                                name: control.name.clone(),
                                gain,
                                hinge_axis: control.hinge_vector(mesh, j).normalize(),
                                hinge_point: 0.5
                                    * (control.hinge_point(mesh, j)
                                        + control.hinge_point(mesh, j + 1)),
                                reflection_sign: control.duplicate_sign,
                            })
                        })
                        .collect(),
                );
            }
            let direction = self.vortices[strip.first_vortex_index].bound_vortex_right
                - self.vortices[strip.first_vortex_index].bound_vortex_left;
            if direction.norm() > EPSILON {
                strip.spanwise_unit_vector = direction.normalize();
            }
            self.strips.push(strip);
        }

        let mut area = 0.0;
        let mut chord_area = 0.0;
        for strip in &self.strips[first_strip..] {
            let strip_area = strip.chord_center * strip.width;
            area += strip_area;
            chord_area += strip.chord_center * strip_area;
        }
        self.surface_info.push(SurfaceMeshInfo {
            name: mesh.name.clone(),
            first_vortex_index: first_vortex,
            num_vortices: nx * ny,
            first_strip_index: first_strip,
            num_strips: ny,
            num_chordwise: nx,
            num_spanwise: ny,
            wetted_area: area,
            average_chord: if area > EPSILON {
                chord_area / area
            } else {
                0.0
            },
            // Not a mirror, so it names its own slot. Falling back to the
            // `Default` of 0 would have pointed every direct-mesh surface at
            // the first one.
            original_surface_index: surface_slot,
            ..SurfaceMeshInfo::default()
        });
        Ok(())
    }

    fn duplicate_surface(&mut self, original_index: usize, plane: Scalar) {
        let original = self.surface_info[original_index].clone();
        let surface_slot = self.surface_info.len();
        let first_vortex = self.vortices.len();
        let first_strip = self.strips.len();
        let y_offset = 2.0 * plane;
        for strip_offset in 0..original.num_strips {
            let old_strip = self.strips[original.first_strip_index + strip_offset].clone();
            let new_first = self.vortices.len();
            let mut strip = old_strip.clone();
            strip.leading_edge_left = reflect_swapped(&old_strip.leading_edge_right, y_offset);
            strip.leading_edge_right = reflect_swapped(&old_strip.leading_edge_left, y_offset);
            strip.leading_edge_center = reflect_swapped(&old_strip.leading_edge_center, y_offset);
            std::mem::swap(&mut strip.chord_left, &mut strip.chord_right);
            strip.leading_edge_sweep_slope *= -1.0;
            strip.trailing_edge_sweep_slope *= -1.0;
            std::mem::swap(&mut strip.twist_angle_left, &mut strip.twist_angle_right);
            strip.trefftz_normal_y *= -1.0;
            strip.spanwise_unit_vector.y *= -1.0;
            strip.reference_point = reflect_swapped(&old_strip.reference_point, y_offset);
            strip.first_vortex_index = new_first;
            strip.surface_index = surface_slot;
            self.strips.push(strip);

            for vortex_offset in 0..old_strip.num_vortices {
                let old_index = old_strip.first_vortex_index + vortex_offset;
                let old = self.vortices[old_index].clone();
                let mut vortex = old.clone();
                vortex.bound_vortex_left = reflect_swapped(&old.bound_vortex_right, y_offset);
                vortex.bound_vortex_right = reflect_swapped(&old.bound_vortex_left, y_offset);
                vortex.bound_vortex_center = reflect_swapped(&old.bound_vortex_center, y_offset);
                vortex.control_point = reflect_swapped(&old.control_point, y_offset);
                vortex.normal_at_control_point.y *= -1.0;
                vortex.normal_at_vortex.y *= -1.0;
                vortex.strip_index = self.strips.len() - 1;
                vortex.surface_index = surface_slot;
                self.vortices.push(vortex);

                let mirrored = self.control_influences[old_index]
                    .iter()
                    .map(|old| ControlSurfaceInfluence {
                        name: old.name.clone(),
                        gain: -old.gain * old.reflection_sign,
                        hinge_axis: Vec3::new(
                            old.hinge_axis.x,
                            -old.hinge_axis.y,
                            old.hinge_axis.z,
                        ),
                        hinge_point: reflect_swapped(&old.hinge_point, y_offset),
                        reflection_sign: old.reflection_sign,
                    })
                    .collect();
                self.control_influences.push(mirrored);
            }
        }
        self.surface_info.push(SurfaceMeshInfo {
            name: format!("{} (yduplicate)", original.name),
            first_vortex_index: first_vortex,
            num_vortices: original.num_vortices,
            first_strip_index: first_strip,
            num_strips: original.num_strips,
            num_chordwise: original.num_chordwise,
            num_spanwise: original.num_spanwise,
            wetted_area: original.wetted_area,
            average_chord: original.average_chord,
            is_y_duplicated: true,
            original_surface_index: original_index,
        });
    }

    pub fn recalculate_normals(&mut self) {
        for strip_index in 0..self.strips.len() {
            let strip = &self.strips[strip_index];
            if strip.num_vortices == 0 {
                continue;
            }
            let leading = &self.vortices[strip.first_vortex_index];
            let trailing = &self.vortices[strip.first_vortex_index + strip.num_vortices - 1];
            let average = (1.0 - self.spanwise_axis_fraction)
                * (leading.bound_vortex_right - leading.bound_vortex_left)
                + self.spanwise_axis_fraction
                    * (trailing.bound_vortex_right - trailing.bound_vortex_left);
            if average.norm() > EPSILON {
                self.strips[strip_index].spanwise_unit_vector = average.normalize();
            }
            let strip = self.strips[strip_index].clone();
            for offset in 0..strip.num_vortices {
                let vortex = &mut self.vortices[strip.first_vortex_index + offset];
                let bound = vortex.bound_vortex_right - vortex.bound_vortex_left;
                let bound = if bound.norm() > EPSILON {
                    bound.normalize()
                } else {
                    bound
                };
                let angle = strip.twist_angle - vortex.camber_slope_at_control_point.atan();
                let chord = Vec3::new(
                    angle.cos(),
                    -angle.sin() * strip.trefftz_normal_y,
                    -angle.sin() * strip.trefftz_normal_z,
                );
                let normal = chord.cross(&bound);
                vortex.normal_at_control_point = if normal.norm() > EPSILON {
                    normal.normalize()
                } else {
                    Vec3::new(0.0, strip.trefftz_normal_y, strip.trefftz_normal_z)
                };
                let angle = strip.twist_angle - vortex.camber_slope_at_vortex.atan();
                let chord = Vec3::new(
                    angle.cos(),
                    -angle.sin() * strip.trefftz_normal_y,
                    -angle.sin() * strip.trefftz_normal_z,
                );
                let normal = chord.cross(&bound);
                vortex.normal_at_vortex = if normal.norm() > EPSILON {
                    normal.normalize()
                } else {
                    Vec3::new(0.0, strip.trefftz_normal_y, strip.trefftz_normal_z)
                };
            }
        }
    }

    #[must_use]
    pub fn vortices(&self) -> &[HorseshoeVortex] {
        &self.vortices
    }

    #[must_use]
    pub fn strips(&self) -> &[Strip] {
        &self.strips
    }

    #[must_use]
    pub fn surface_info(&self) -> &[SurfaceMeshInfo] {
        &self.surface_info
    }

    #[must_use]
    pub fn control_influences(&self) -> &[Vec<ControlSurfaceInfluence>] {
        &self.control_influences
    }

    #[must_use]
    pub fn control_surface_names(&self) -> &[String] {
        &self.control_surface_names
    }

    #[must_use]
    pub fn mesh_input(&self) -> Option<&MeshInput> {
        self.mesh_input.as_ref()
    }

    #[must_use]
    pub fn vortex_count(&self) -> usize {
        self.vortices.len()
    }

    #[must_use]
    pub fn strip_count(&self) -> usize {
        self.strips.len()
    }

    #[must_use]
    pub fn surface_count(&self) -> usize {
        self.surface_info.len()
    }

    pub(crate) fn vortices_mut(&mut self) -> &mut [HorseshoeVortex] {
        &mut self.vortices
    }
}

/// Relative tolerance for spanwise station separation.
///
/// Applied against the longest interval on the surface rather than as an
/// absolute length, so the check means the same thing whatever unit the model
/// is expressed in.
const MIN_STATION_SEPARATION: Scalar = 1.0e-9;

/// Rejects section geometry the mesher cannot turn into a usable lattice.
///
/// Each of these cases previously generated "successfully" and produced NaN
/// panel data that silently poisoned the whole AIC:
///
/// - a zero chord makes `claf` a `0.0 / 0.0`;
/// - two sections at the same spanwise station leave a zero-length interval,
///   whose interpolation parameter is `0.0 / 0.0`. Repeating a station is a
///   standard AVL idiom for a spacing or control-surface break, so this is
///   reachable from ordinary input;
/// - sections that backtrack pass the arc-length test, because the arc
///   accumulates `|Δ|` and is monotonic regardless of ordering, but produce
///   strips that physically overlap and a near-singular AIC.
fn validate_sections(surface: &Surface) -> Result<()> {
    let sections = surface.sections();
    let name = surface.name();

    for (index, section) in sections.iter().enumerate() {
        let leading_edge = section.leading_edge();
        if !leading_edge.iter().all(|value| value.is_finite())
            || !section.chord().is_finite()
            || !section.twist_degrees().is_finite()
        {
            return Err(PegasusError::InvalidGeometry(format!(
                "surface '{name}' section {index} has non-finite geometry"
            )));
        }
        if section.chord() <= 0.0 {
            return Err(PegasusError::InvalidGeometry(format!(
                "surface '{name}' section {index} has chord {}; chord must be positive",
                section.chord()
            )));
        }
    }

    let deltas: Vec<Vec3> = sections
        .windows(2)
        .map(|pair| pair[1].leading_edge() - pair[0].leading_edge())
        .collect();
    let longest = deltas
        .iter()
        .map(|delta| delta.yz().norm())
        .fold(0.0, Scalar::max);
    if longest <= 0.0 {
        return Err(PegasusError::InvalidGeometry(format!(
            "surface '{name}' has zero span"
        )));
    }

    for (index, delta) in deltas.iter().enumerate() {
        if delta.yz().norm() <= MIN_STATION_SEPARATION * longest {
            return Err(PegasusError::InvalidGeometry(format!(
                "surface '{name}' sections {index} and {} share a spanwise station; \
                 use a single section, or separate them to define a break",
                index + 1
            )));
        }
    }

    // A reversal in the y-z plane means the surface doubles back on itself.
    for (index, pair) in deltas.windows(2).enumerate() {
        if pair[0].yz().dot(&pair[1].yz()) < 0.0 {
            return Err(PegasusError::InvalidGeometry(format!(
                "surface '{name}' reverses spanwise direction at section {}; \
                 sections must be ordered monotonically along the span",
                index + 1
            )));
        }
    }

    Ok(())
}

fn reflect_swapped(point: &Vec3, y_offset: Scalar) -> Vec3 {
    Vec3::new(point.x, y_offset - point.y, point.z)
}

pub(crate) fn chordwise_spacing(
    count: usize,
    distribution: SpacingType,
    cl_alpha_factor: Scalar,
) -> (Vec<Scalar>, Vec<Scalar>) {
    let edges = spacing(count + 1, distribution);
    let mut vortex = vec![0.0; count];
    let mut control = vec![0.0; count];
    for i in 0..count {
        let dx = edges[i + 1] - edges[i];
        let fraction = if cl_alpha_factor != 1.0 && i == 0 {
            (0.25 / cl_alpha_factor).clamp(0.1, 0.4)
        } else {
            0.25
        };
        vortex[i] = edges[i] + fraction * dx;
        control[i] = edges[i] + 0.75 * dx;
    }
    (vortex, control)
}
