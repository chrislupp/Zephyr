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

//! Analytical kernels, adjoints, and geometry sensitivities.

mod adjoint;
mod aic;
mod biot_savart;
mod calculator;
mod coupled;
mod forces;
pub(crate) mod geometry_adjoint;
mod mach;
pub(crate) mod section_geometry;
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
    AeroGradients, FullGradients, GeometryGradients, GeometryVariable, GradientCalculator,
    MeshNodeGradients, SectionGeometryGradients, SurfaceMeshNodeGradients,
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
