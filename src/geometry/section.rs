//! Spanwise lifting-surface station.

use std::sync::OnceLock;

use crate::core::{DEGREES_TO_RADIANS, Scalar, SpacingType, Vec3};

use super::Airfoil;

/// Control-surface definition at a spanwise station.
#[derive(Clone, Debug)]
pub struct ControlSurfaceDefinition {
    /// Shared control name.
    pub name: String,
    /// Hinge chord fraction.
    pub hinge_x_fraction: Scalar,
    /// Hinge-axis direction.
    pub hinge_vector: Vec3,
    /// Deflection sign.
    pub sign: Scalar,
    /// Deflection gain.
    pub gain: Scalar,
}

impl Default for ControlSurfaceDefinition {
    fn default() -> Self {
        Self {
            name: String::new(),
            hinge_x_fraction: 0.0,
            hinge_vector: Vec3::zeros(),
            sign: 1.0,
            gain: 1.0,
        }
    }
}

/// Optimization design-variable metadata.
#[derive(Clone, Debug, Default)]
pub struct DesignVariableDefinition {
    /// Variable name.
    pub name: String,
    /// Station interpolation weight.
    pub weight: Scalar,
}

/// A spanwise station defining leading edge, chord, twist, and airfoil.
#[derive(Clone, Debug)]
pub struct Section {
    pub(crate) leading_edge: Vec3,
    pub(crate) chord: Scalar,
    pub(crate) twist_degrees: Scalar,
    pub(crate) airfoil: Option<Airfoil>,
    pub(crate) cl_alpha_factor: Scalar,
    pub(crate) control_surfaces: Vec<ControlSurfaceDefinition>,
    pub(crate) design_variables: Vec<DesignVariableDefinition>,
    pub(crate) spanwise_panel_count: usize,
    pub(crate) spanwise_spacing: SpacingType,
    pub(crate) viscous_polar: [Scalar; 6],
}

impl Default for Section {
    fn default() -> Self {
        Self::new(Vec3::zeros(), 1.0, 0.0)
    }
}

impl Section {
    /// Construct a station. Twist is in degrees.
    #[must_use]
    pub fn new(leading_edge: Vec3, chord: Scalar, twist_degrees: Scalar) -> Self {
        Self {
            leading_edge,
            chord,
            twist_degrees,
            airfoil: None,
            cl_alpha_factor: 1.0,
            control_surfaces: Vec::new(),
            design_variables: Vec::new(),
            spanwise_panel_count: 0,
            spanwise_spacing: SpacingType::None,
            viscous_polar: [0.0; 6],
        }
    }

    /// Construct with an owned airfoil.
    #[must_use]
    pub fn with_airfoil(
        leading_edge: Vec3,
        chord: Scalar,
        twist_degrees: Scalar,
        airfoil: Airfoil,
    ) -> Self {
        let mut section = Self::new(leading_edge, chord, twist_degrees);
        section.airfoil = Some(airfoil);
        section
    }

    /// Leading-edge position.
    #[must_use]
    pub fn leading_edge(&self) -> &Vec3 {
        &self.leading_edge
    }

    /// Change leading-edge position.
    pub fn set_leading_edge(&mut self, value: Vec3) {
        self.leading_edge = value;
    }

    /// Chord.
    #[must_use]
    pub fn chord(&self) -> Scalar {
        self.chord
    }

    /// Change chord.
    pub fn set_chord(&mut self, value: Scalar) {
        self.chord = value;
    }

    /// Twist in degrees.
    #[must_use]
    pub fn twist_degrees(&self) -> Scalar {
        self.twist_degrees
    }

    /// Twist in radians.
    #[must_use]
    pub fn twist_radians(&self) -> Scalar {
        self.twist_degrees * DEGREES_TO_RADIANS
    }

    /// Change twist in degrees.
    pub fn set_twist_degrees(&mut self, value: Scalar) {
        self.twist_degrees = value;
    }

    /// Get the custom airfoil or a shared flat plate.
    #[must_use]
    pub fn airfoil(&self) -> &Airfoil {
        self.airfoil.as_ref().unwrap_or_else(|| {
            static FLAT: OnceLock<Airfoil> = OnceLock::new();
            FLAT.get_or_init(Airfoil::flat_plate)
        })
    }

    /// Get a mutable airfoil, creating a flat plate if necessary.
    pub fn airfoil_mut(&mut self) -> &mut Airfoil {
        self.airfoil.get_or_insert_with(Airfoil::flat_plate)
    }

    /// Replace the airfoil.
    pub fn set_airfoil(&mut self, value: Airfoil) {
        self.airfoil = Some(value);
    }

    /// Whether a custom airfoil is present.
    #[must_use]
    pub fn has_airfoil(&self) -> bool {
        self.airfoil.is_some()
    }

    /// Section lift-curve multiplier.
    #[must_use]
    pub fn cl_alpha_factor(&self) -> Scalar {
        self.cl_alpha_factor
    }

    /// Change the section lift-curve multiplier.
    pub fn set_cl_alpha_factor(&mut self, value: Scalar) {
        self.cl_alpha_factor = value;
    }

    /// Add control metadata.
    pub fn add_control_surface(&mut self, control: ControlSurfaceDefinition) {
        self.control_surfaces.push(control);
    }

    /// Controls defined at this station.
    #[must_use]
    pub fn control_surfaces(&self) -> &[ControlSurfaceDefinition] {
        &self.control_surfaces
    }

    /// Add design-variable metadata.
    pub fn add_design_variable(&mut self, variable: DesignVariableDefinition) {
        self.design_variables.push(variable);
    }

    /// Design variables defined at this station.
    #[must_use]
    pub fn design_variables(&self) -> &[DesignVariableDefinition] {
        &self.design_variables
    }

    /// Spanwise panel override (`0` means surface default).
    #[must_use]
    pub fn spanwise_panel_count(&self) -> usize {
        self.spanwise_panel_count
    }

    /// Change the spanwise panel override.
    pub fn set_spanwise_panel_count(&mut self, value: usize) {
        self.spanwise_panel_count = value;
    }

    /// Station spacing override.
    #[must_use]
    pub fn spanwise_spacing(&self) -> SpacingType {
        self.spanwise_spacing
    }

    /// Change the station spacing override.
    pub fn set_spanwise_spacing(&mut self, value: SpacingType) {
        self.spanwise_spacing = value;
    }

    /// Convert the AVL numeric spacing convention.
    pub fn set_numeric_spanwise_spacing(&mut self, value: Scalar) {
        self.spanwise_spacing = if value.abs() < 0.01 {
            SpacingType::Uniform
        } else if value > 0.0 {
            SpacingType::Cosine
        } else {
            SpacingType::Sine
        };
    }

    /// Section viscous-polar coefficients.
    #[must_use]
    pub fn viscous_polar(&self) -> &[Scalar; 6] {
        &self.viscous_polar
    }

    /// Change section viscous-polar coefficients.
    pub fn set_viscous_polar(&mut self, value: [Scalar; 6]) {
        self.viscous_polar = value;
    }
}
