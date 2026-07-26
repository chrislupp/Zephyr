//! Complete aircraft definition.

use std::collections::BTreeSet;

use crate::core::{Scalar, SymmetryType, Vec3};

use super::{Body, Surface};

/// Top-level geometry and reference quantities.
#[derive(Clone, Debug)]
pub struct Aircraft {
    title: String,
    reference_area: Scalar,
    reference_chord: Scalar,
    reference_span: Scalar,
    moment_reference_point: Vec3,
    mach_number: Scalar,
    profile_drag_coefficient: Scalar,
    y_symmetry: SymmetryType,
    z_symmetry: SymmetryType,
    z_symmetry_plane: Scalar,
    surfaces: Vec<Surface>,
    bodies: Vec<Body>,
}

impl Default for Aircraft {
    fn default() -> Self {
        Self::new("")
    }
}

impl Aircraft {
    /// Construct a named aircraft.
    #[must_use]
    pub fn new(title: impl Into<String>) -> Self {
        Self {
            title: title.into(),
            reference_area: 1.0,
            reference_chord: 1.0,
            reference_span: 1.0,
            moment_reference_point: Vec3::zeros(),
            mach_number: 0.0,
            profile_drag_coefficient: 0.0,
            y_symmetry: SymmetryType::None,
            z_symmetry: SymmetryType::None,
            z_symmetry_plane: 0.0,
            surfaces: Vec::new(),
            bodies: Vec::new(),
        }
    }

    /// Title.
    #[must_use]
    pub fn title(&self) -> &str {
        &self.title
    }

    /// Change title.
    pub fn set_title(&mut self, value: impl Into<String>) {
        self.title = value.into();
    }

    /// Reference area.
    #[must_use]
    pub fn reference_area(&self) -> Scalar {
        self.reference_area
    }

    /// Change reference area.
    pub fn set_reference_area(&mut self, value: Scalar) {
        self.reference_area = value;
    }

    /// Reference chord.
    #[must_use]
    pub fn reference_chord(&self) -> Scalar {
        self.reference_chord
    }

    /// Change reference chord.
    pub fn set_reference_chord(&mut self, value: Scalar) {
        self.reference_chord = value;
    }

    /// Reference span.
    #[must_use]
    pub fn reference_span(&self) -> Scalar {
        self.reference_span
    }

    /// Change reference span.
    pub fn set_reference_span(&mut self, value: Scalar) {
        self.reference_span = value;
    }

    /// Moment reference point.
    #[must_use]
    pub fn moment_reference_point(&self) -> &Vec3 {
        &self.moment_reference_point
    }

    /// Change moment reference point.
    pub fn set_moment_reference_point(&mut self, value: Vec3) {
        self.moment_reference_point = value;
    }

    /// Reference Mach number.
    #[must_use]
    pub fn mach_number(&self) -> Scalar {
        self.mach_number
    }

    /// Change reference Mach.
    pub fn set_mach_number(&mut self, value: Scalar) {
        self.mach_number = value;
    }

    /// Constant profile-drag coefficient.
    #[must_use]
    pub fn profile_drag_coefficient(&self) -> Scalar {
        self.profile_drag_coefficient
    }

    /// Change profile drag.
    pub fn set_profile_drag_coefficient(&mut self, value: Scalar) {
        self.profile_drag_coefficient = value;
    }

    /// Y symmetry.
    #[must_use]
    pub fn y_symmetry(&self) -> SymmetryType {
        self.y_symmetry
    }

    /// Change Y symmetry.
    pub fn set_y_symmetry(&mut self, value: SymmetryType) {
        self.y_symmetry = value;
    }

    /// Z symmetry.
    #[must_use]
    pub fn z_symmetry(&self) -> SymmetryType {
        self.z_symmetry
    }

    /// Change Z symmetry.
    pub fn set_z_symmetry(&mut self, value: SymmetryType) {
        self.z_symmetry = value;
    }

    /// Z symmetry plane.
    #[must_use]
    pub fn z_symmetry_plane(&self) -> Scalar {
        self.z_symmetry_plane
    }

    /// Change Z symmetry plane.
    pub fn set_z_symmetry_plane(&mut self, value: Scalar) {
        self.z_symmetry_plane = value;
    }

    /// Append a surface.
    pub fn add_surface(&mut self, value: Surface) {
        self.surfaces.push(value);
    }

    /// Surfaces.
    #[must_use]
    pub fn surfaces(&self) -> &[Surface] {
        &self.surfaces
    }

    /// Mutable surfaces.
    pub fn surfaces_mut(&mut self) -> &mut [Surface] {
        &mut self.surfaces
    }

    /// Find a surface by exact name.
    #[must_use]
    pub fn surface_named(&self, name: &str) -> Option<&Surface> {
        self.surfaces.iter().find(|surface| surface.name() == name)
    }

    /// Append a body.
    pub fn add_body(&mut self, value: Body) {
        self.bodies.push(value);
    }

    /// Bodies.
    #[must_use]
    pub fn bodies(&self) -> &[Body] {
        &self.bodies
    }

    /// Sorted unique control names.
    #[must_use]
    pub fn control_surface_names(&self) -> Vec<String> {
        self.surfaces
            .iter()
            .flat_map(Surface::sections)
            .flat_map(super::Section::control_surfaces)
            .map(|control| control.name.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect()
    }

    /// Sorted unique design-variable names.
    #[must_use]
    pub fn design_variable_names(&self) -> Vec<String> {
        self.surfaces
            .iter()
            .flat_map(Surface::sections)
            .flat_map(super::Section::design_variables)
            .map(|variable| variable.name.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect()
    }

    /// Total generated horseshoe count.
    #[must_use]
    pub fn total_vortex_count(&self) -> usize {
        self.surfaces
            .iter()
            .map(|surface| {
                surface.chordwise_panel_count()
                    * surface.spanwise_panel_count()
                    * if surface.y_duplicate().is_some() {
                        2
                    } else {
                        1
                    }
            })
            .sum()
    }

    /// Total generated strip count.
    #[must_use]
    pub fn total_strip_count(&self) -> usize {
        self.surfaces
            .iter()
            .map(|surface| {
                surface.spanwise_panel_count()
                    * if surface.y_duplicate().is_some() {
                        2
                    } else {
                        1
                    }
            })
            .sum()
    }

    /// Sum of surface planform areas.
    #[must_use]
    pub fn total_planform_area(&self) -> Scalar {
        self.surfaces.iter().map(Surface::planform_area).sum()
    }

    /// Maximum surface span.
    #[must_use]
    pub fn total_span(&self) -> Scalar {
        self.surfaces
            .iter()
            .map(Surface::span)
            .fold(0.0, Scalar::max)
    }
}
