//! Rigid-body mass and flight-state definitions.

use crate::core::{Mat3, Scalar, Vec3};
use crate::error::{PegasusError, Result};

/// Mass, center of gravity, and inertia about the center of gravity.
#[derive(Clone, Debug)]
pub struct MassProperties {
    /// Total mass.
    pub mass: Scalar,
    /// Center of gravity.
    pub center_of_gravity: Vec3,
    /// Roll inertia.
    pub i_xx: Scalar,
    /// Pitch inertia.
    pub i_yy: Scalar,
    /// Yaw inertia.
    pub i_zz: Scalar,
    /// Stored-positive XY product of inertia.
    pub i_xy: Scalar,
    /// Stored-positive YZ product of inertia.
    pub i_yz: Scalar,
    /// Stored-positive XZ product of inertia.
    pub i_xz: Scalar,
}

impl Default for MassProperties {
    fn default() -> Self {
        Self {
            mass: 0.0,
            center_of_gravity: Vec3::zeros(),
            i_xx: 0.0,
            i_yy: 0.0,
            i_zz: 0.0,
            i_xy: 0.0,
            i_yz: 0.0,
            i_xz: 0.0,
        }
    }
}

impl MassProperties {
    /// Full symmetric inertia tensor.
    #[must_use]
    pub fn inertia_tensor(&self) -> Mat3 {
        Mat3::new(
            self.i_xx, -self.i_xy, -self.i_xz, -self.i_xy, self.i_yy, -self.i_yz, -self.i_xz,
            -self.i_yz, self.i_zz,
        )
    }

    /// Set scalar fields from a symmetric inertia tensor.
    pub fn set_inertia_tensor(&mut self, tensor: &Mat3) {
        self.i_xx = tensor[(0, 0)];
        self.i_yy = tensor[(1, 1)];
        self.i_zz = tensor[(2, 2)];
        self.i_xy = -tensor[(0, 1)];
        self.i_yz = -tensor[(1, 2)];
        self.i_xz = -tensor[(0, 2)];
    }

    /// Whether mass and principal inertias are positive.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        self.mass > 0.0 && self.i_xx > 0.0 && self.i_yy > 0.0 && self.i_zz > 0.0
    }

    /// Inverse inertia tensor.
    pub fn inverse_inertia_tensor(&self) -> Result<Mat3> {
        self.inertia_tensor()
            .try_inverse()
            .ok_or(PegasusError::SingularSystem)
    }
}

/// Flight state used by dynamic-mode analysis.
#[derive(Clone, Debug)]
pub struct FlightState {
    /// Airspeed magnitude.
    pub velocity: Scalar,
    /// Angle of attack in radians.
    pub alpha: Scalar,
    /// Sideslip in radians.
    pub beta: Scalar,
    /// Normalized roll rate.
    pub roll_rate: Scalar,
    /// Normalized pitch rate.
    pub pitch_rate: Scalar,
    /// Normalized yaw rate.
    pub yaw_rate: Scalar,
    /// Bank angle.
    pub phi: Scalar,
    /// Pitch attitude.
    pub theta: Scalar,
    /// Heading.
    pub psi: Scalar,
    /// Density.
    pub density: Scalar,
    /// Gravity.
    pub gravity: Scalar,
}

impl Default for FlightState {
    fn default() -> Self {
        Self {
            velocity: 0.0,
            alpha: 0.0,
            beta: 0.0,
            roll_rate: 0.0,
            pitch_rate: 0.0,
            yaw_rate: 0.0,
            phi: 0.0,
            theta: 0.0,
            psi: 0.0,
            density: 1.225,
            gravity: 9.81,
        }
    }
}

impl FlightState {
    /// Body-axis freestream velocity.
    #[must_use]
    pub fn body_velocity(&self) -> Vec3 {
        let (sa, ca) = self.alpha.sin_cos();
        let (sb, cb) = self.beta.sin_cos();
        self.velocity * Vec3::new(ca * cb, sb, sa * cb)
    }

    /// Dynamic pressure.
    #[must_use]
    pub fn dynamic_pressure(&self) -> Scalar {
        0.5 * self.density * self.velocity * self.velocity
    }
}

/// Optional added-mass terms.
#[derive(Clone, Debug)]
pub struct ApparentMass {
    /// Added translational inertia.
    pub mass_matrix: Mat3,
    /// Added rotational inertia.
    pub inertia_matrix: Mat3,
    /// Whether to include these terms.
    pub include_effects: bool,
}

impl Default for ApparentMass {
    fn default() -> Self {
        Self {
            mass_matrix: Mat3::zeros(),
            inertia_matrix: Mat3::zeros(),
            include_effects: false,
        }
    }
}
