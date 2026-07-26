//! Mathematical helpers and interpolation.

mod cubic_spline;
mod math;

pub use cubic_spline::{CubicSpline, SplineBoundaryCondition};
pub use math::{
    blended_spacing, body_to_stability_rotation, cosine_spacing, prandtl_glauert_factor,
    rotate_about_x, rotate_about_y, rotate_about_z, sine_spacing, spacing,
    trailing_vortex_velocity, uniform_spacing, vortex_segment_velocity, wind_to_body_rotation,
};
