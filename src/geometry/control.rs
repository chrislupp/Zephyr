//! Controls attached to material coordinates of a structured mesh.

use std::ops::Range;

use crate::core::{EPSILON, Scalar, Vec3};
use crate::error::{PegasusError, Result};

use super::SurfaceMesh;

/// A normal-rotation control attached to the incoming, structurally deformed
/// mesh **before** control actuation. No mesh nodes or wake vertices are moved.
///
/// Hinges use chordwise mesh coordinates: `0` is the leading edge, `nx` the
/// trailing edge, and `i + f` interpolates between node rows `i` and `i + 1`.
/// These are material coordinates, not physical chord fractions. Keep them and
/// the span interval fixed while differentiating with respect to mesh nodes.
#[derive(Clone, Debug, PartialEq)]
pub struct MeshControlSurface {
    /// Shared command name; matching definitions on other surfaces are mixed.
    pub name: String,
    /// Active spanwise strips, start inclusive and end exclusive.
    pub span: Range<usize>,
    /// Hinge coordinate at every spanwise node station (`ny + 1` entries).
    pub hinge: Vec<Scalar>,
    /// Rotate panels ahead of the hinge instead of behind it.
    pub leading_edge: bool,
    /// Rotation per radian of command, using the right-hand rule about the
    /// hinge from station `j` to `j + 1`.
    pub gain: Scalar,
    /// Relative deflection on the mirrored copy: +1 symmetric, -1 antisymmetric.
    pub duplicate_sign: Scalar,
}

impl MeshControlSurface {
    /// Define a trailing-edge control with unit gain and symmetric duplication.
    #[must_use]
    pub fn trailing_edge(name: impl Into<String>, span: Range<usize>, hinge: Vec<Scalar>) -> Self {
        Self {
            name: name.into(),
            span,
            hinge,
            leading_edge: false,
            gain: 1.0,
            duplicate_sign: 1.0,
        }
    }

    pub(crate) fn validate_layout(&self, nx: usize, ny: usize) -> Result<()> {
        if self.name.trim().is_empty()
            || self.span.start >= self.span.end
            || self.span.end > ny
            || self.hinge.len() != ny + 1
            || self
                .hinge
                .iter()
                .any(|h| !h.is_finite() || *h < 0.0 || *h > nx as Scalar)
            || !self.gain.is_finite()
            || !self.duplicate_sign.is_finite()
        {
            return Err(PegasusError::InvalidMesh(format!(
                "invalid layout for control '{}'",
                self.name
            )));
        }
        Ok(())
    }

    pub(crate) fn validate(&self, surface: &SurfaceMesh) -> Result<()> {
        self.validate_layout(surface.nx, surface.ny)?;
        for j in self.span.clone() {
            if self.hinge_vector(surface, j).norm() <= EPSILON {
                return Err(PegasusError::InvalidMesh(format!(
                    "control '{}' has a degenerate hinge in strip {j}",
                    self.name
                )));
            }
        }
        Ok(())
    }

    /// The fixed material weights of one hinge station.
    pub(crate) fn weights(&self, surface: &SurfaceMesh, j: usize) -> [(usize, Scalar); 2] {
        let h = self.hinge[j];
        let i = (h.floor() as usize).min(surface.nx - 1);
        let f = h - i as Scalar;
        [
            (i * (surface.ny + 1) + j, 1.0 - f),
            ((i + 1) * (surface.ny + 1) + j, f),
        ]
    }

    pub(crate) fn hinge_point(&self, surface: &SurfaceMesh, j: usize) -> Vec3 {
        self.weights(surface, j)
            .into_iter()
            .map(|(k, w)| w * surface.nodes[k])
            .sum()
    }

    pub(crate) fn hinge_vector(&self, surface: &SurfaceMesh, j: usize) -> Vec3 {
        self.hinge_point(surface, j + 1) - self.hinge_point(surface, j)
    }

    pub(crate) fn panel_gain(&self, i: usize, j: usize) -> Scalar {
        if !self.span.contains(&j) {
            return 0.0;
        }
        let aft = 0.5
            * ((i as Scalar + 1.0 - self.hinge[j]).clamp(0.0, 1.0)
                + (i as Scalar + 1.0 - self.hinge[j + 1]).clamp(0.0, 1.0));
        self.gain * if self.leading_edge { 1.0 - aft } else { aft }
    }
}
