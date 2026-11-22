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

//! Direct structured-surface mesh input.

use crate::core::{Scalar, Vec3};
use crate::error::{Result, ZephyrError};

/// Coefficient-normalization quantities.
#[derive(Clone, Debug, PartialEq)]
pub struct ReferenceQuantities {
    /// Reference area.
    pub area: Scalar,
    /// Reference mean chord.
    pub chord: Scalar,
    /// Reference span.
    pub span: Scalar,
    /// Moment reference point.
    pub moment_reference_point: Vec3,
}

impl Default for ReferenceQuantities {
    fn default() -> Self {
        Self {
            area: 1.0,
            chord: 1.0,
            span: 1.0,
            moment_reference_point: Vec3::zeros(),
        }
    }
}

impl From<&super::Aircraft> for ReferenceQuantities {
    fn from(aircraft: &super::Aircraft) -> Self {
        Self {
            area: aircraft.reference_area(),
            chord: aircraft.reference_chord(),
            span: aircraft.reference_span(),
            moment_reference_point: *aircraft.moment_reference_point(),
        }
    }
}

/// Structured `(nx+1) × (ny+1)` node grid.
#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceMesh {
    /// Controls attached to this mesh before control actuation.
    pub controls: Vec<super::MeshControlSurface>,
    /// Surface name.
    pub name: String,
    /// Chordwise panel count.
    pub nx: usize,
    /// Spanwise panel count.
    pub ny: usize,
    /// Row-major nodes, indexed by `i*(ny+1)+j`.
    pub nodes: Vec<Vec3>,
    /// Component index used by core regularization.
    pub component_index: usize,
    /// Whether freestream alpha and beta affect the surface.
    pub affected_by_freestream: bool,
    /// Whether the surface sheds a wake.
    pub has_wake: bool,
    /// Whether to generate a mirrored copy.
    pub y_duplicated: bool,
    /// Mirror plane.
    pub y_duplicate_plane: Scalar,
}

impl SurfaceMesh {
    /// Allocate a zero-filled node grid.
    #[must_use]
    pub fn new(name: impl Into<String>, nx: usize, ny: usize) -> Self {
        Self {
            controls: Vec::new(),
            name: name.into(),
            nx,
            ny,
            nodes: vec![Vec3::zeros(); (nx + 1) * (ny + 1)],
            component_index: 0,
            affected_by_freestream: true,
            has_wake: true,
            y_duplicated: false,
            y_duplicate_plane: 0.0,
        }
    }

    /// Node count.
    #[must_use]
    pub const fn node_count(&self) -> usize {
        (self.nx + 1) * (self.ny + 1)
    }

    /// Panel count.
    #[must_use]
    pub const fn panel_count(&self) -> usize {
        self.nx * self.ny
    }

    /// Validate dimensions and finite coordinates.
    pub fn validate(&self) -> Result<()> {
        if self.nx == 0 || self.ny == 0 {
            return Err(ZephyrError::InvalidMesh(format!(
                "surface '{}' needs positive nx and ny",
                self.name
            )));
        }
        if self.nodes.len() != self.node_count() {
            return Err(ZephyrError::InvalidMesh(format!(
                "surface '{}' has {} nodes, expected {}",
                self.name,
                self.nodes.len(),
                self.node_count()
            )));
        }
        if self
            .nodes
            .iter()
            .any(|node| !node.iter().all(|x| x.is_finite()))
        {
            return Err(ZephyrError::InvalidMesh(format!(
                "surface '{}' contains non-finite coordinates",
                self.name
            )));
        }
        for control in &self.controls {
            control.validate(self)?;
        }
        Ok(())
    }

    /// Checked node access.
    #[must_use]
    pub fn node(&self, i: usize, j: usize) -> Option<&Vec3> {
        (i <= self.nx && j <= self.ny).then(|| &self.nodes[i * (self.ny + 1) + j])
    }

    /// Set one node.
    pub fn set_node(&mut self, i: usize, j: usize, point: Vec3) -> Result<()> {
        if i > self.nx || j > self.ny {
            return Err(ZephyrError::InvalidMesh(format!(
                "node ({i},{j}) is outside {} × {} surface '{}'",
                self.nx + 1,
                self.ny + 1,
                self.name
            )));
        }
        self.nodes[i * (self.ny + 1) + j] = point;
        Ok(())
    }
}

/// Complete direct-mesh analysis input.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MeshInput {
    /// Surface grids.
    pub surfaces: Vec<SurfaceMesh>,
    /// Normalization quantities.
    pub reference: ReferenceQuantities,
}

impl MeshInput {
    /// Validate all surfaces and reference values.
    pub fn validate(&self) -> Result<()> {
        if self.surfaces.is_empty() {
            return Err(ZephyrError::InvalidMesh(
                "at least one surface is required".into(),
            ));
        }
        for surface in &self.surfaces {
            surface.validate()?;
        }
        if self.reference.area <= 0.0 || self.reference.chord <= 0.0 || self.reference.span <= 0.0 {
            return Err(ZephyrError::InvalidMesh(
                "reference area, chord, and span must be positive".into(),
            ));
        }
        Ok(())
    }
}
