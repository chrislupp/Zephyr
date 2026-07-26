//! Fundamental numeric types and constants.

use nalgebra::{
    DMatrix, DVector, Matrix2, Matrix3, Matrix4, RowDVector, RowVector3, Vector2, Vector3, Vector4,
};

/// Scalar used by all aerodynamic calculations.
pub type Scalar = f64;
/// Two-component vector.
pub type Vec2 = Vector2<Scalar>;
/// Three-component vector.
pub type Vec3 = Vector3<Scalar>;
/// Four-component vector.
pub type Vec4 = Vector4<Scalar>;
/// Two-by-two matrix.
pub type Mat2 = Matrix2<Scalar>;
/// Three-by-three matrix.
pub type Mat3 = Matrix3<Scalar>;
/// Four-by-four matrix.
pub type Mat4 = Matrix4<Scalar>;
/// Dynamically sized vector.
pub type VecX = DVector<Scalar>;
/// Dynamically sized matrix.
pub type MatX = DMatrix<Scalar>;
/// Three-component row vector.
pub type RowVec3 = RowVector3<Scalar>;
/// Dynamically sized row vector.
pub type RowVecX = RowDVector<Scalar>;
/// Index type used throughout the crate.
pub type Index = usize;

/// Mathematical constant π.
pub const PI: Scalar = std::f64::consts::PI;
/// Degrees-to-radians conversion.
pub const DEGREES_TO_RADIANS: Scalar = PI / 180.0;
/// Radians-to-degrees conversion.
pub const RADIANS_TO_DEGREES: Scalar = 180.0 / PI;
/// Small comparison epsilon.
pub const EPSILON: Scalar = 1.0e-12;
/// Default iterative-solver tolerance.
pub const DEFAULT_TOLERANCE: Scalar = 1.0e-10;

/// Image-vortex symmetry type.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
#[repr(i8)]
pub enum SymmetryType {
    /// No symmetry image.
    #[default]
    None = 0,
    /// Same-sign image.
    Symmetric = 1,
    /// Opposite-sign image.
    Antisymmetric = -1,
}

impl SymmetryType {
    /// Numeric image strength multiplier.
    #[must_use]
    pub const fn factor(self) -> Scalar {
        self as i8 as Scalar
    }
}

/// Panel-spacing distribution.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum SpacingType {
    /// No override; callers use their default.
    #[default]
    None,
    /// Uniform distribution.
    Uniform,
    /// Cosine clustering at both ends.
    Cosine,
    /// Sine clustering near the trailing end.
    Sine,
    /// Equal blend of uniform and cosine.
    Blend,
}
