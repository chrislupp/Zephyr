//! Data shared by lattice generation and aerodynamic post-processing.

use crate::core::{MatX, Scalar, Vec3, VecX};

/// One horseshoe-vortex panel.
#[derive(Clone, Debug)]
pub struct HorseshoeVortex {
    pub bound_vortex_left: Vec3,
    pub bound_vortex_right: Vec3,
    pub bound_vortex_center: Vec3,
    pub control_point: Vec3,
    pub normal_at_control_point: Vec3,
    pub normal_at_vortex: Vec3,
    pub camber_slope_at_control_point: Scalar,
    pub camber_slope_at_vortex: Scalar,
    pub chordwise_width: Scalar,
    pub local_chord: Scalar,
    pub strip_index: usize,
    pub surface_index: usize,
    pub component_index: usize,
    pub affected_by_freestream: bool,
    pub enforce_flow_tangency: bool,
    pub circulation: Scalar,
    pub delta_cp: Scalar,
}

impl Default for HorseshoeVortex {
    fn default() -> Self {
        Self {
            bound_vortex_left: Vec3::zeros(),
            bound_vortex_right: Vec3::zeros(),
            bound_vortex_center: Vec3::zeros(),
            control_point: Vec3::zeros(),
            normal_at_control_point: Vec3::z(),
            normal_at_vortex: Vec3::z(),
            camber_slope_at_control_point: 0.0,
            camber_slope_at_vortex: 0.0,
            chordwise_width: 0.0,
            local_chord: 0.0,
            strip_index: 0,
            surface_index: 0,
            component_index: 0,
            affected_by_freestream: true,
            enforce_flow_tangency: true,
            circulation: 0.0,
            delta_cp: 0.0,
        }
    }
}

/// A spanwise strip containing chordwise vortex elements.
#[derive(Clone, Debug)]
pub struct Strip {
    pub leading_edge_left: Vec3,
    pub leading_edge_right: Vec3,
    pub leading_edge_center: Vec3,
    pub chord_left: Scalar,
    pub chord_right: Scalar,
    pub chord_center: Scalar,
    pub width: Scalar,
    pub leading_edge_sweep_slope: Scalar,
    pub trailing_edge_sweep_slope: Scalar,
    pub twist_angle: Scalar,
    pub twist_angle_left: Scalar,
    pub twist_angle_right: Scalar,
    pub trefftz_normal_y: Scalar,
    pub trefftz_normal_z: Scalar,
    pub spanwise_unit_vector: Vec3,
    pub reference_point: Vec3,
    pub first_vortex_index: usize,
    pub num_vortices: usize,
    pub surface_index: usize,
    pub viscous_polar: [Scalar; 6],
    pub has_viscous_data: bool,
    pub is_active: bool,
}

impl Default for Strip {
    fn default() -> Self {
        Self {
            leading_edge_left: Vec3::zeros(),
            leading_edge_right: Vec3::zeros(),
            leading_edge_center: Vec3::zeros(),
            chord_left: 0.0,
            chord_right: 0.0,
            chord_center: 0.0,
            width: 0.0,
            leading_edge_sweep_slope: 0.0,
            trailing_edge_sweep_slope: 0.0,
            twist_angle: 0.0,
            twist_angle_left: 0.0,
            twist_angle_right: 0.0,
            trefftz_normal_y: 0.0,
            trefftz_normal_z: 1.0,
            spanwise_unit_vector: Vec3::y(),
            reference_point: Vec3::zeros(),
            first_vortex_index: 0,
            num_vortices: 0,
            surface_index: 0,
            viscous_polar: [0.0; 6],
            has_viscous_data: false,
            is_active: true,
        }
    }
}

/// Effect of a named control on one panel.
#[derive(Clone, Debug)]
pub struct ControlSurfaceInfluence {
    pub name: String,
    pub gain: Scalar,
    pub hinge_axis: Vec3,
    pub hinge_point: Vec3,
    pub reflection_sign: Scalar,
}

impl Default for ControlSurfaceInfluence {
    fn default() -> Self {
        Self {
            name: String::new(),
            gain: 0.0,
            hinge_axis: Vec3::y(),
            hinge_point: Vec3::zeros(),
            reflection_sign: 1.0,
        }
    }
}

/// Contiguous ranges and aggregate data for one generated surface.
#[derive(Clone, Debug, Default)]
pub struct SurfaceMeshInfo {
    pub name: String,
    pub first_vortex_index: usize,
    pub num_vortices: usize,
    pub first_strip_index: usize,
    pub num_strips: usize,
    pub num_chordwise: usize,
    pub num_spanwise: usize,
    pub wetted_area: Scalar,
    pub average_chord: Scalar,
    pub is_y_duplicated: bool,
    /// Slot in `VortexLattice::surface_info` of the surface this one mirrors.
    ///
    /// Always a `surface_info` index, never an index into
    /// `Aircraft::surfaces()`. A surface that is not a mirrored copy names its
    /// own slot. The two generators previously disagreed about which index
    /// space this was, so for `[wing(duplicated), tail]` the tail's value
    /// pointed at the wing's mirror; any consumer using it as a `surface_info`
    /// index silently got the wrong surface.
    pub original_surface_index: usize,
}

/// Nondimensional flight state. Angles, rates, and controls are radians.
#[derive(Clone, Debug)]
pub struct FlightCondition {
    pub alpha: Scalar,
    pub beta: Scalar,
    pub roll_rate: Scalar,
    pub pitch_rate: Scalar,
    pub yaw_rate: Scalar,
    pub mach: Scalar,
    pub dynamic_pressure: Scalar,
    /// Ordered control entries, intentionally not a map.
    pub control_deflections: Vec<(String, Scalar)>,
}

impl Default for FlightCondition {
    fn default() -> Self {
        Self {
            alpha: 0.0,
            beta: 0.0,
            roll_rate: 0.0,
            pitch_rate: 0.0,
            yaw_rate: 0.0,
            mach: 0.0,
            dynamic_pressure: 1.0,
            control_deflections: Vec::new(),
        }
    }
}

/// Circulation and induced velocity produced by a solve.
#[derive(Clone, Debug, Default)]
pub struct SolutionState {
    pub circulation: VecX,
    pub circulation_per_u: VecX,
    pub circulation_per_v: VecX,
    pub circulation_per_w: VecX,
    pub circulation_per_p: VecX,
    pub circulation_per_q: VecX,
    pub circulation_per_r: VecX,
    /// `N × 3`, in body axes.
    pub induced_velocity: MatX,
    pub is_valid: bool,
}
