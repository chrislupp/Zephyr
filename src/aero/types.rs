//! Public aerodynamic result records.

use crate::core::{PI, Scalar, Vec3};

/// Six conventional aerodynamic coefficients.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Coefficients {
    /// Stability-axis lift.
    pub lift: Scalar,
    /// Stability-axis drag.
    pub drag: Scalar,
    /// Side force.
    pub side: Scalar,
    /// Rolling moment.
    pub roll: Scalar,
    /// Pitching moment.
    pub pitch: Scalar,
    /// Yawing moment.
    pub yaw: Scalar,
}

impl Coefficients {
    pub(crate) fn linear_combination(
        u: Self,
        u_scale: Scalar,
        v: Self,
        v_scale: Scalar,
        w: Self,
        w_scale: Scalar,
    ) -> Self {
        Self {
            lift: u_scale * u.lift + v_scale * v.lift + w_scale * w.lift,
            drag: u_scale * u.drag + v_scale * v.drag + w_scale * w.drag,
            side: u_scale * u.side + v_scale * v.side + w_scale * w.side,
            roll: u_scale * u.roll + v_scale * v.roll + w_scale * w.roll,
            pitch: u_scale * u.pitch + v_scale * v.pitch + w_scale * w.pitch,
            yaw: u_scale * u.yaw + v_scale * v.yaw + w_scale * w.yaw,
        }
    }
}

/// Loads on one spanwise strip.
#[derive(Clone, Debug, Default)]
pub struct StripLoads {
    pub strip_index: usize,
    pub chord: Scalar,
    pub width: Scalar,
    pub area: Scalar,
    pub leading_edge_position: Vec3,
    pub quarter_chord_position: Vec3,
    pub coefficients: Coefficients,
    pub body_force: Vec3,
    pub body_moment: Vec3,
    pub normal_span_loading: Scalar,
    pub induced_angle: Scalar,
    pub induced_velocity: Vec3,
    /// Section-profile drag coefficient for this strip.
    pub viscous_drag: Scalar,
}

/// Loads aggregated over one generated surface.
#[derive(Clone, Debug, Default)]
pub struct SurfaceLoads {
    pub surface_index: usize,
    pub name: String,
    pub coefficients: Coefficients,
    pub body_force: Vec3,
    pub body_moment: Vec3,
    pub area: Scalar,
    pub strip_loads: Vec<StripLoads>,
}

/// Complete near- and far-field aerodynamic result.
#[derive(Clone, Debug)]
pub struct AerodynamicResults {
    pub is_valid: bool,
    pub coefficients: Coefficients,
    pub body_force: Vec3,
    pub body_moment: Vec3,
    pub near_field_induced_drag: Scalar,
    pub viscous_drag: Scalar,
    pub far_field_lift: Scalar,
    pub far_field_side: Scalar,
    pub far_field_induced_drag: Scalar,
    pub span_efficiency: Scalar,
    pub reference_area: Scalar,
    pub reference_chord: Scalar,
    pub reference_span: Scalar,
    pub moment_reference_point: Vec3,
    pub surface_loads: Vec<SurfaceLoads>,
}

impl Default for AerodynamicResults {
    fn default() -> Self {
        Self {
            is_valid: false,
            coefficients: Coefficients::default(),
            body_force: Vec3::zeros(),
            body_moment: Vec3::zeros(),
            near_field_induced_drag: 0.0,
            viscous_drag: 0.0,
            far_field_lift: 0.0,
            far_field_side: 0.0,
            far_field_induced_drag: 0.0,
            span_efficiency: 0.0,
            reference_area: 1.0,
            reference_chord: 1.0,
            reference_span: 1.0,
            moment_reference_point: Vec3::zeros(),
            surface_loads: Vec::new(),
        }
    }
}

impl AerodynamicResults {
    /// Geometric aspect ratio.
    #[must_use]
    pub fn aspect_ratio(&self) -> Scalar {
        self.reference_span * self.reference_span / self.reference_area
    }

    /// Ideal elliptic induced drag at the current lift coefficient.
    #[must_use]
    pub fn ideal_induced_drag(&self) -> Scalar {
        self.coefficients.lift.powi(2) / (PI * self.aspect_ratio())
    }
}

/// Derivatives with respect to one named control.
#[derive(Clone, Debug, Default)]
pub struct ControlDerivatives {
    pub name: String,
    pub deflection: Scalar,
    pub gain: Scalar,
    pub coefficients: Coefficients,
    pub body_force: Vec3,
    pub hinge_moment: Scalar,
}

/// Stability-axis derivatives grouped by perturbation input.
#[derive(Clone, Debug)]
pub struct StabilityDerivatives {
    pub is_valid: bool,
    pub reference_area: Scalar,
    pub reference_chord: Scalar,
    pub reference_span: Scalar,
    pub velocity: Scalar,
    pub alpha: Scalar,
    pub beta: Scalar,
    pub by_alpha: Coefficients,
    pub by_beta: Coefficients,
    pub by_roll_rate: Coefficients,
    pub by_pitch_rate: Coefficients,
    pub by_yaw_rate: Coefficients,
    pub control_derivatives: Vec<ControlDerivatives>,
    pub neutral_point_x: Scalar,
    pub static_margin: Scalar,
}

impl Default for StabilityDerivatives {
    fn default() -> Self {
        Self {
            is_valid: false,
            reference_area: 1.0,
            reference_chord: 1.0,
            reference_span: 1.0,
            velocity: 1.0,
            alpha: 0.0,
            beta: 0.0,
            by_alpha: Coefficients::default(),
            by_beta: Coefficients::default(),
            by_roll_rate: Coefficients::default(),
            by_pitch_rate: Coefficients::default(),
            by_yaw_rate: Coefficients::default(),
            control_derivatives: Vec::new(),
            neutral_point_x: 0.0,
            static_margin: 0.0,
        }
    }
}

/// Body-axis force/moment derivative groups.
#[derive(Clone, Debug)]
pub struct BodyAxisDerivatives {
    pub is_valid: bool,
    pub reference_area: Scalar,
    pub reference_chord: Scalar,
    pub reference_span: Scalar,
    pub velocity: Scalar,
    pub alpha: Scalar,
    pub beta: Scalar,
    pub by_alpha_force: Vec3,
    pub by_alpha_moment: Vec3,
    pub by_beta_force: Vec3,
    pub by_beta_moment: Vec3,
    pub by_roll_rate_force: Vec3,
    pub by_roll_rate_moment: Vec3,
    pub by_pitch_rate_force: Vec3,
    pub by_pitch_rate_moment: Vec3,
    pub by_yaw_rate_force: Vec3,
    pub by_yaw_rate_moment: Vec3,
}

impl Default for BodyAxisDerivatives {
    fn default() -> Self {
        Self {
            is_valid: false,
            reference_area: 1.0,
            reference_chord: 1.0,
            reference_span: 1.0,
            velocity: 1.0,
            alpha: 0.0,
            beta: 0.0,
            by_alpha_force: Vec3::zeros(),
            by_alpha_moment: Vec3::zeros(),
            by_beta_force: Vec3::zeros(),
            by_beta_moment: Vec3::zeros(),
            by_roll_rate_force: Vec3::zeros(),
            by_roll_rate_moment: Vec3::zeros(),
            by_pitch_rate_force: Vec3::zeros(),
            by_pitch_rate_moment: Vec3::zeros(),
            by_yaw_rate_force: Vec3::zeros(),
            by_yaw_rate_moment: Vec3::zeros(),
        }
    }
}
