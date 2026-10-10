//! Direct structured-surface mesh input.

use crate::core::{Scalar, Vec3};
use crate::error::{PegasusError, Result};

/// Coefficient-normalization quantities.
#[derive(Clone, Debug)]
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
#[derive(Clone, Debug)]
pub struct SurfaceMesh {
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
            return Err(PegasusError::InvalidMesh(format!(
                "surface '{}' needs positive nx and ny",
                self.name
            )));
        }
        if self.nodes.len() != self.node_count() {
            return Err(PegasusError::InvalidMesh(format!(
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
            return Err(PegasusError::InvalidMesh(format!(
                "surface '{}' contains non-finite coordinates",
                self.name
            )));
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
            return Err(PegasusError::InvalidMesh(format!(
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
#[derive(Clone, Debug, Default)]
pub struct MeshInput {
    /// Use full collocation self influence and bound-vortex force evaluation.
    /// False retains the historical v0.2 coefficient convention.
    pub consistent_forces: bool,
    /// Surface grids.
    pub surfaces: Vec<SurfaceMesh>,
    /// Normalization quantities.
    pub reference: ReferenceQuantities,
}

impl MeshInput {
    /// Validate all surfaces and reference values.
    pub fn validate(&self) -> Result<()> {
        if self.surfaces.is_empty() {
            return Err(PegasusError::InvalidMesh(
                "at least one surface is required".into(),
            ));
        }
        for surface in &self.surfaces {
            surface.validate()?;
        }
        if self.reference.area <= 0.0 || self.reference.chord <= 0.0 || self.reference.span <= 0.0 {
            return Err(PegasusError::InvalidMesh(
                "reference area, chord, and span must be positive".into(),
            ));
        }
        Ok(())
    }
}
