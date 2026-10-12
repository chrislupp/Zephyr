//! Analytical kernels, adjoints, and geometry sensitivities.

mod adjoint;
mod aic;
mod biot_savart;
mod calculator;
mod coupled;
mod forces;
mod geometry_adjoint;
mod mach;
mod section_geometry;
mod selection;
mod trefftz;

pub use adjoint::{AdjointSolver, TotalDerivative};
pub use aic::{
    AicElementGradient, aic_column_gradient, aic_element_gradient, aic_row_gradient,
    velocity_column_gradient, velocity_row_gradient,
};
pub use biot_savart::{
    HorseshoeVelocityJacobians, horseshoe_velocity_jacobians, horseshoe_velocity_pg_derivative,
    horseshoe_velocity_vjp_bound_left, horseshoe_velocity_vjp_bound_right,
    horseshoe_velocity_vjp_field_point,
};
pub use calculator::{
    AeroGradients, ControlGradient, FullGradients, GeometryGradients, GeometryVariable,
    GradientCalculator, MeshNodeGradients, SectionGeometryGradients, SurfaceMeshNodeGradients,
};
pub use coupled::{MeshNodeCotangent, StripLoadAdjoint, StripLoadCotangent};
pub use forces::{
    ForceCoefficientGradients, force_alpha_gradient, force_beta_gradient,
    force_circulation_gradient, force_geometry_jvp, force_induced_velocity_gradient,
    selected_force_circulation_gradients,
};
pub use section_geometry::{interpolated_twist_jvp, panel_normal_jvp};
pub use selection::{AeroOutput, OutputSelection};
pub use trefftz::{
    TrefftzCirculationGradients, trefftz_circulation_gradients, trefftz_geometry_jvp,
};
