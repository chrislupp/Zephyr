//! Lifting-surface planform definition.

use crate::core::{EPSILON, Scalar, SpacingType, Vec3};

use super::Section;

/// Wing, tail, or other lifting component.
#[derive(Clone, Debug)]
pub struct Surface {
    pub(crate) name: String,
    pub(crate) component_index: usize,
    pub(crate) sections: Vec<Section>,
    pub(crate) chordwise_panel_count: usize,
    pub(crate) chordwise_spacing: SpacingType,
    pub(crate) spanwise_panel_count: usize,
    pub(crate) spanwise_spacing: SpacingType,
    pub(crate) y_duplicate: Option<Scalar>,
    pub(crate) angle_offset_degrees: Scalar,
    pub(crate) translation: Vec3,
    pub(crate) scale: Vec3,
    pub(crate) has_wake: bool,
    pub(crate) contributes_to_loads: bool,
    pub(crate) affected_by_alpha_beta: bool,
    pub(crate) cl_cd_polar: [Scalar; 6],
}

impl Default for Surface {
    fn default() -> Self {
        Self::new("")
    }
}

impl Surface {
    /// Construct a named surface with Pegasus defaults.
    #[must_use]
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            component_index: 0,
            sections: Vec::new(),
            chordwise_panel_count: 12,
            chordwise_spacing: SpacingType::Cosine,
            spanwise_panel_count: 20,
            spanwise_spacing: SpacingType::Cosine,
            y_duplicate: None,
            angle_offset_degrees: 0.0,
            translation: Vec3::zeros(),
            scale: Vec3::repeat(1.0),
            has_wake: true,
            contributes_to_loads: true,
            affected_by_alpha_beta: true,
            cl_cd_polar: [0.0; 6],
        }
    }

    /// Name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Change name.
    pub fn set_name(&mut self, value: impl Into<String>) {
        self.name = value.into();
    }

    /// Component index used by core regularization.
    #[must_use]
    pub fn component_index(&self) -> usize {
        self.component_index
    }

    /// Change component index.
    pub fn set_component_index(&mut self, value: usize) {
        self.component_index = value;
    }

    /// Append a spanwise station.
    pub fn add_section(&mut self, section: Section) {
        self.sections.push(section);
    }

    /// Stations.
    #[must_use]
    pub fn sections(&self) -> &[Section] {
        &self.sections
    }

    /// Mutable stations. Callers must rebuild any mesh created from them.
    pub fn sections_mut(&mut self) -> &mut [Section] {
        &mut self.sections
    }

    /// Checked station access.
    #[must_use]
    pub fn section(&self, index: usize) -> Option<&Section> {
        self.sections.get(index)
    }

    /// Checked mutable station access.
    pub fn section_mut(&mut self, index: usize) -> Option<&mut Section> {
        self.sections.get_mut(index)
    }

    /// Chordwise panel count.
    #[must_use]
    pub fn chordwise_panel_count(&self) -> usize {
        self.chordwise_panel_count
    }

    /// Set chordwise panel count.
    pub fn set_chordwise_panel_count(&mut self, value: usize) {
        self.chordwise_panel_count = value;
    }

    /// Chordwise spacing.
    #[must_use]
    pub fn chordwise_spacing(&self) -> SpacingType {
        self.chordwise_spacing
    }

    /// Set chordwise spacing.
    pub fn set_chordwise_spacing(&mut self, value: SpacingType) {
        self.chordwise_spacing = value;
    }

    /// Set chordwise spacing from AVL's signed numeric convention.
    pub fn set_numeric_chordwise_spacing(&mut self, value: Scalar) {
        self.chordwise_spacing = numeric_spacing(value);
    }

    /// Spanwise panel count.
    #[must_use]
    pub fn spanwise_panel_count(&self) -> usize {
        self.spanwise_panel_count
    }

    /// Set spanwise panel count.
    pub fn set_spanwise_panel_count(&mut self, value: usize) {
        self.spanwise_panel_count = value;
    }

    /// Spanwise spacing.
    #[must_use]
    pub fn spanwise_spacing(&self) -> SpacingType {
        self.spanwise_spacing
    }

    /// Set spanwise spacing.
    pub fn set_spanwise_spacing(&mut self, value: SpacingType) {
        self.spanwise_spacing = value;
    }

    /// Set spanwise spacing from AVL's signed numeric convention.
    pub fn set_numeric_spanwise_spacing(&mut self, value: Scalar) {
        self.spanwise_spacing = numeric_spacing(value);
    }

    /// Y-plane used to generate a mirror copy.
    #[must_use]
    pub fn y_duplicate(&self) -> Option<Scalar> {
        self.y_duplicate
    }

    /// Enable Y duplication.
    pub fn set_y_duplicate(&mut self, plane: Scalar) {
        self.y_duplicate = Some(plane);
    }

    /// Disable Y duplication.
    pub fn clear_y_duplicate(&mut self) {
        self.y_duplicate = None;
    }

    /// Surface incidence offset in degrees.
    #[must_use]
    pub fn angle_offset_degrees(&self) -> Scalar {
        self.angle_offset_degrees
    }

    /// Change incidence offset in degrees.
    pub fn set_angle_offset_degrees(&mut self, value: Scalar) {
        self.angle_offset_degrees = value;
    }

    /// Translation.
    #[must_use]
    pub fn translation(&self) -> &Vec3 {
        &self.translation
    }

    /// Change translation.
    pub fn set_translation(&mut self, value: Vec3) {
        self.translation = value;
    }

    /// Axis scale factors.
    #[must_use]
    pub fn scale(&self) -> &Vec3 {
        &self.scale
    }

    /// Change axis scale factors.
    pub fn set_scale(&mut self, value: Vec3) {
        self.scale = value;
    }

    /// Whether this surface sheds a wake.
    #[must_use]
    pub fn has_wake(&self) -> bool {
        self.has_wake
    }

    /// Set wake behavior.
    pub fn set_has_wake(&mut self, value: bool) {
        self.has_wake = value;
    }

    /// Whether surface forces are accumulated.
    #[must_use]
    pub fn contributes_to_loads(&self) -> bool {
        self.contributes_to_loads
    }

    /// Set load behavior.
    pub fn set_contributes_to_loads(&mut self, value: bool) {
        self.contributes_to_loads = value;
    }

    /// Whether alpha and beta affect the surface.
    #[must_use]
    pub fn affected_by_alpha_beta(&self) -> bool {
        self.affected_by_alpha_beta
    }

    /// Set alpha/beta behavior.
    pub fn set_affected_by_alpha_beta(&mut self, value: bool) {
        self.affected_by_alpha_beta = value;
    }

    /// Surface drag-polar coefficients.
    #[must_use]
    pub fn cl_cd_polar(&self) -> &[Scalar; 6] {
        &self.cl_cd_polar
    }

    /// Change surface drag-polar coefficients.
    pub fn set_cl_cd_polar(&mut self, value: [Scalar; 6]) {
        self.cl_cd_polar = value;
    }

    /// Whether any surface drag-polar coefficient is active.
    #[must_use]
    pub fn has_viscous_drag(&self) -> bool {
        self.cl_cd_polar.iter().any(|value| value.abs() > EPSILON)
    }

    /// Trapezoidal planform area.
    #[must_use]
    pub fn planform_area(&self) -> Scalar {
        let area: Scalar = self
            .sections
            .windows(2)
            .map(|pair| {
                let delta = pair[1].leading_edge - pair[0].leading_edge;
                0.5 * (pair[0].chord + pair[1].chord) * delta.yz().norm()
            })
            .sum();
        if self.y_duplicate.is_some() {
            2.0 * area
        } else {
            area
        }
    }

    /// Mean aerodynamic chord using the original trapezoidal integration.
    #[must_use]
    pub fn mean_aerodynamic_chord(&self) -> Scalar {
        if self.sections.len() < 2 {
            return self.sections.first().map_or(0.0, Section::chord);
        }
        let mut numerator = 0.0;
        let mut denominator = 0.0;
        for pair in self.sections.windows(2) {
            let distance = (pair[1].leading_edge - pair[0].leading_edge).yz().norm();
            numerator += 0.5 * (pair[0].chord.powi(2) + pair[1].chord.powi(2)) * distance;
            denominator += 0.5 * (pair[0].chord + pair[1].chord) * distance;
        }
        if denominator > EPSILON {
            numerator / denominator
        } else {
            0.0
        }
    }

    /// Span measured along Y, doubled for a duplicated semisurface.
    #[must_use]
    pub fn span(&self) -> Scalar {
        let Some(first) = self.sections.first() else {
            return 0.0;
        };
        let (minimum, maximum) = self.sections.iter().fold(
            (first.leading_edge.y, first.leading_edge.y),
            |(minimum, maximum), section| {
                (
                    minimum.min(section.leading_edge.y),
                    maximum.max(section.leading_edge.y),
                )
            },
        );
        let span = maximum - minimum;
        if self.y_duplicate.is_some() {
            2.0 * span
        } else {
            span
        }
    }
}

fn numeric_spacing(value: Scalar) -> SpacingType {
    if value.abs() < 0.01 {
        SpacingType::Uniform
    } else if value > 0.0 {
        SpacingType::Cosine
    } else {
        SpacingType::Sine
    }
}
