//! Aerodynamic influence-coefficient kernels.

use crate::core::{EPSILON, MatX, PI, Scalar, SymmetryType, Vec3, VecX};

use super::VortexLattice;

const INV_4_PI: Scalar = 1.0 / (4.0 * PI);

/// Symmetry planes used while assembling influence matrices.
#[derive(Clone, Copy, Debug)]
pub struct SymmetryConfig {
    pub y_symmetry: SymmetryType,
    pub y_symmetry_plane: Scalar,
    pub z_symmetry: SymmetryType,
    pub z_symmetry_plane: Scalar,
}

impl Default for SymmetryConfig {
    fn default() -> Self {
        Self {
            y_symmetry: SymmetryType::None,
            y_symmetry_plane: 0.0,
            z_symmetry: SymmetryType::None,
            z_symmetry_plane: 0.0,
        }
    }
}

/// Velocity induced by a unit horseshoe using the legacy Leishman `R⁴` core.
#[must_use]
pub fn horseshoe_vortex_velocity(
    field_point: &Vec3,
    bound_left: &Vec3,
    bound_right: &Vec3,
    prandtl_glauert_factor: Scalar,
    core_radius: Scalar,
    include_bound: bool,
) -> Vec3 {
    let a = Vec3::new(
        (bound_left.x - field_point.x) / prandtl_glauert_factor,
        bound_left.y - field_point.y,
        bound_left.z - field_point.z,
    );
    let b = Vec3::new(
        (bound_right.x - field_point.x) / prandtl_glauert_factor,
        bound_right.y - field_point.y,
        bound_right.z - field_point.z,
    );
    let a_sq = a.norm_squared();
    let b_sq = b.norm_squared();
    let a_mag = a_sq.sqrt();
    let b_mag = b_sq.sqrt();
    let core_sq = core_radius * core_radius;
    let core_4 = core_sq * core_sq;

    let (mut u, mut v, mut w) = (0.0, 0.0, 0.0);
    if include_bound && a_mag * b_mag > EPSILON {
        let cross = a.cross(&b);
        let cross_sq = cross.norm_squared();
        if cross_sq > EPSILON {
            let dot = a.dot(&b);
            let length_sq = a_sq + b_sq - 2.0 * dot;
            let factor = ((b_sq - dot) / (b_sq * b_sq + core_4).sqrt().sqrt()
                + (a_sq - dot) / (a_sq * a_sq + core_4).sqrt().sqrt())
                / (cross_sq * cross_sq + length_sq * length_sq * core_4).sqrt();
            u = cross.x * factor;
            v = cross.y * factor;
            w = cross.z * factor;
        }
    }
    if a_mag > EPSILON {
        let axis_sq = a.y * a.y + a.z * a.z;
        let factor = -(1.0 - a.x / a_mag) / (axis_sq * axis_sq + core_4).sqrt();
        v += a.z * factor;
        w -= a.y * factor;
    }
    if b_mag > EPSILON {
        let axis_sq = b.y * b.y + b.z * b.z;
        let factor = (1.0 - b.x / b_mag) / (axis_sq * axis_sq + core_4).sqrt();
        v += b.z * factor;
        w -= b.y * factor;
    }
    Vec3::new(
        u * INV_4_PI / prandtl_glauert_factor,
        v * INV_4_PI,
        w * INV_4_PI,
    )
}

/// Builds velocity and normalwash influence matrices for a fixed mesh.
pub struct InfluenceMatrixBuilder<'a> {
    mesh: &'a VortexLattice,
    mach_number: Scalar,
    prandtl_glauert_factor: Scalar,
    symmetry: SymmetryConfig,
    core_radius_chord_fraction: Scalar,
    core_radius_width_fraction: Scalar,
}

impl<'a> InfluenceMatrixBuilder<'a> {
    #[must_use]
    pub fn new(mesh: &'a VortexLattice) -> Self {
        Self {
            mesh,
            mach_number: 0.0,
            prandtl_glauert_factor: 1.0,
            symmetry: SymmetryConfig::default(),
            core_radius_chord_fraction: 0.0001,
            core_radius_width_fraction: 0.0001,
        }
    }

    pub fn set_mach_number(&mut self, mach: Scalar) {
        self.mach_number = mach;
        self.prandtl_glauert_factor = if mach < 1.0 {
            (1.0 - mach * mach).sqrt()
        } else {
            1.0
        };
    }

    #[must_use]
    pub fn mach_number(&self) -> Scalar {
        self.mach_number
    }

    #[must_use]
    pub fn prandtl_glauert_factor(&self) -> Scalar {
        self.prandtl_glauert_factor
    }

    pub fn set_symmetry(&mut self, symmetry: SymmetryConfig) {
        self.symmetry = symmetry;
    }

    #[must_use]
    pub fn symmetry(&self) -> SymmetryConfig {
        self.symmetry
    }

    pub fn set_core_radius_chord_fraction(&mut self, fraction: Scalar) {
        self.core_radius_chord_fraction = fraction;
    }

    pub fn set_core_radius_width_fraction(&mut self, fraction: Scalar) {
        self.core_radius_width_fraction = fraction;
    }

    fn core_radius(&self, vortex_index: usize, same_component: bool) -> Scalar {
        let vortex = &self.mesh.vortices()[vortex_index];
        let panel_width = (vortex.bound_vortex_right - vortex.bound_vortex_left).norm();
        let mut core = 0.0001 * panel_width;
        if !same_component {
            core = (self.core_radius_chord_fraction * vortex.local_chord)
                .max(self.core_radius_width_fraction * panel_width);
        }
        core
    }

    #[must_use]
    pub fn single_vortex_velocity(
        &self,
        field_point: &Vec3,
        vortex_index: usize,
        include_bound: bool,
    ) -> Vec3 {
        let vortex = &self.mesh.vortices()[vortex_index];
        horseshoe_vortex_velocity(
            field_point,
            &vortex.bound_vortex_left,
            &vortex.bound_vortex_right,
            self.prandtl_glauert_factor,
            self.core_radius(vortex_index, true),
            include_bound,
        )
    }

    fn y_image_velocity(&self, field_point: &Vec3, vortex_index: usize, core: Scalar) -> Vec3 {
        let vortex = &self.mesh.vortices()[vortex_index];
        let offset = 2.0 * self.symmetry.y_symmetry_plane;
        let left = Vec3::new(
            vortex.bound_vortex_right.x,
            offset - vortex.bound_vortex_right.y,
            vortex.bound_vortex_right.z,
        );
        let right = Vec3::new(
            vortex.bound_vortex_left.x,
            offset - vortex.bound_vortex_left.y,
            vortex.bound_vortex_left.z,
        );
        let include_bound = !(self.symmetry.y_symmetry == SymmetryType::Symmetric
            && (field_point - 0.5 * (left + right)).norm_squared() < EPSILON);
        self.symmetry.y_symmetry.factor()
            * horseshoe_vortex_velocity(
                field_point,
                &left,
                &right,
                self.prandtl_glauert_factor,
                core,
                include_bound,
            )
    }

    fn z_image_velocity(&self, field_point: &Vec3, vortex_index: usize, core: Scalar) -> Vec3 {
        let vortex = &self.mesh.vortices()[vortex_index];
        let offset = 2.0 * self.symmetry.z_symmetry_plane;
        let left = Vec3::new(
            vortex.bound_vortex_right.x,
            vortex.bound_vortex_right.y,
            offset - vortex.bound_vortex_right.z,
        );
        let right = Vec3::new(
            vortex.bound_vortex_left.x,
            vortex.bound_vortex_left.y,
            offset - vortex.bound_vortex_left.z,
        );
        self.symmetry.z_symmetry.factor()
            * horseshoe_vortex_velocity(
                field_point,
                &left,
                &right,
                self.prandtl_glauert_factor,
                core,
                true,
            )
    }

    #[must_use]
    pub fn velocity_influence_matrix(&self) -> MatX {
        let n = self.mesh.vortex_count();
        let mut matrix = MatX::zeros(3 * n, n);
        for i in 0..n {
            let target = &self.mesh.vortices()[i];
            for j in 0..n {
                let source = &self.mesh.vortices()[j];
                let same_component = target.component_index == source.component_index;
                let core = self.core_radius(j, same_component);
                let mut velocity = horseshoe_vortex_velocity(
                    &target.control_point,
                    &source.bound_vortex_left,
                    &source.bound_vortex_right,
                    self.prandtl_glauert_factor,
                    core,
                    i != j,
                );
                if self.symmetry.y_symmetry != SymmetryType::None {
                    velocity += self.y_image_velocity(&target.control_point, j, core);
                }
                if self.symmetry.z_symmetry != SymmetryType::None {
                    velocity += self.z_image_velocity(&target.control_point, j, core);
                    if self.symmetry.y_symmetry != SymmetryType::None {
                        let y_offset = 2.0 * self.symmetry.y_symmetry_plane;
                        let z_offset = 2.0 * self.symmetry.z_symmetry_plane;
                        let left = Vec3::new(
                            source.bound_vortex_left.x,
                            y_offset - source.bound_vortex_left.y,
                            z_offset - source.bound_vortex_left.z,
                        );
                        let right = Vec3::new(
                            source.bound_vortex_right.x,
                            y_offset - source.bound_vortex_right.y,
                            z_offset - source.bound_vortex_right.z,
                        );
                        velocity += self.symmetry.y_symmetry.factor()
                            * self.symmetry.z_symmetry.factor()
                            * horseshoe_vortex_velocity(
                                &target.control_point,
                                &left,
                                &right,
                                self.prandtl_glauert_factor,
                                core,
                                true,
                            );
                    }
                }
                matrix[(3 * i, j)] = velocity.x;
                matrix[(3 * i + 1, j)] = velocity.y;
                matrix[(3 * i + 2, j)] = velocity.z;
            }
        }
        matrix
    }

    #[must_use]
    pub fn normalwash_matrix(&self) -> MatX {
        self.normalwash_and_velocity().0
    }

    #[must_use]
    pub fn normalwash_and_velocity(&self) -> (MatX, MatX) {
        let velocity = self.velocity_influence_matrix();
        let n = self.mesh.vortex_count();
        let mut aic = MatX::zeros(n, n);
        for i in 0..n {
            let normal = self.mesh.vortices()[i].normal_at_control_point;
            for j in 0..n {
                aic[(i, j)] = velocity[(3 * i, j)] * normal.x
                    + velocity[(3 * i + 1, j)] * normal.y
                    + velocity[(3 * i + 2, j)] * normal.z;
            }
        }
        (aic, velocity)
    }

    #[must_use]
    pub fn induced_velocity(&self, field_point: &Vec3, circulation: &VecX) -> Vec3 {
        let mut total = Vec3::zeros();
        for (j, source) in self.mesh.vortices().iter().enumerate() {
            let core = self.core_radius(j, false);
            let mut velocity = horseshoe_vortex_velocity(
                field_point,
                &source.bound_vortex_left,
                &source.bound_vortex_right,
                self.prandtl_glauert_factor,
                core,
                true,
            );
            if self.symmetry.y_symmetry != SymmetryType::None {
                velocity += self.y_image_velocity(field_point, j, core);
            }
            if self.symmetry.z_symmetry != SymmetryType::None {
                velocity += self.z_image_velocity(field_point, j, core);
            }
            total += circulation[j] * velocity;
        }
        total
    }
}
