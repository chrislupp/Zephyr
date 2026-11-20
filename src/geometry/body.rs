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

//! Fuselage-body metadata parsed from AVL files.

use crate::core::{Scalar, Vec3};

/// Elliptic body cross-section.
#[derive(Clone, Debug, Default)]
pub struct BodySection {
    /// Section center.
    pub center: Vec3,
    /// Half-width.
    pub width: Scalar,
    /// Half-height.
    pub height: Scalar,
}

/// Body definition.
#[derive(Clone, Debug)]
pub struct Body {
    name: String,
    sections: Vec<BodySection>,
    node_count: usize,
    node_spacing: Scalar,
    translation: Vec3,
    scale: Vec3,
    y_duplicated: bool,
}

impl Default for Body {
    fn default() -> Self {
        Self::new("")
    }
}

impl Body {
    /// Construct a named body.
    #[must_use]
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            sections: Vec::new(),
            node_count: 20,
            node_spacing: 1.0,
            translation: Vec3::zeros(),
            scale: Vec3::repeat(1.0),
            y_duplicated: false,
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

    /// Append a section.
    pub fn add_section(&mut self, value: BodySection) {
        self.sections.push(value);
    }

    /// Sections.
    #[must_use]
    pub fn sections(&self) -> &[BodySection] {
        &self.sections
    }

    /// Number of source/doublet nodes.
    #[must_use]
    pub fn node_count(&self) -> usize {
        self.node_count
    }

    /// Change node count.
    pub fn set_node_count(&mut self, value: usize) {
        self.node_count = value;
    }

    /// Node-spacing parameter.
    #[must_use]
    pub fn node_spacing(&self) -> Scalar {
        self.node_spacing
    }

    /// Change node spacing.
    pub fn set_node_spacing(&mut self, value: Scalar) {
        self.node_spacing = value;
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

    /// Axis scales.
    #[must_use]
    pub fn scale(&self) -> &Vec3 {
        &self.scale
    }

    /// Change scales.
    pub fn set_scale(&mut self, value: Vec3) {
        self.scale = value;
    }

    /// Whether AVL requests a Y duplicate.
    #[must_use]
    pub fn is_y_duplicated(&self) -> bool {
        self.y_duplicated
    }

    /// Change Y-duplication behavior.
    pub fn set_y_duplicated(&mut self, value: bool) {
        self.y_duplicated = value;
    }
}
