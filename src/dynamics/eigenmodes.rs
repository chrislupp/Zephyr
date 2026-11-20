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

//! Linearized rigid-body matrices and dynamic-mode decomposition.

use faer::Mat;
use num_complex::Complex64;

use crate::aero::StabilityDerivatives;
use crate::core::{EPSILON, MatX, PI, Scalar};
use crate::error::{Result, ZephyrError};
use crate::geometry::Aircraft;

use super::{FlightState, MassProperties};

/// Conventional dynamic-mode classification.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum ModeType {
    Unknown,
    ShortPeriod,
    Phugoid,
    DutchRoll,
    RollSubsidence,
    Spiral,
    Other,
}

impl ModeType {
    /// Display name used by reports.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::ShortPeriod => "Short Period",
            Self::Phugoid => "Phugoid",
            Self::DutchRoll => "Dutch Roll",
            Self::RollSubsidence => "Roll Subsidence",
            Self::Spiral => "Spiral",
            Self::Unknown | Self::Other => "Unknown",
        }
    }
}

/// One normalized complex eigenpair.
#[derive(Clone, Debug)]
pub struct Eigenmode {
    pub eigenvalue: Complex64,
    pub eigenvector: Vec<Complex64>,
    pub mode_type: ModeType,
}

impl Default for Eigenmode {
    fn default() -> Self {
        Self {
            eigenvalue: Complex64::default(),
            eigenvector: Vec::new(),
            mode_type: ModeType::Unknown,
        }
    }
}

impl Eigenmode {
    #[must_use]
    pub fn damping(&self) -> Scalar {
        self.eigenvalue.re
    }

    #[must_use]
    pub fn frequency(&self) -> Scalar {
        self.eigenvalue.im.abs()
    }

    #[must_use]
    pub fn period(&self) -> Scalar {
        let frequency = self.frequency();
        if frequency > 0.0 {
            2.0 * PI / frequency
        } else {
            0.0
        }
    }

    #[must_use]
    pub fn time_constant(&self) -> Scalar {
        let damping = self.damping();
        if damping != 0.0 {
            -Scalar::ln(2.0) / damping
        } else {
            0.0
        }
    }

    #[must_use]
    pub fn damping_ratio(&self) -> Scalar {
        let natural_frequency = self.damping().hypot(self.frequency());
        if natural_frequency > 0.0 {
            -self.damping() / natural_frequency
        } else {
            0.0
        }
    }

    #[must_use]
    pub fn is_stable(&self) -> bool {
        self.damping() < 0.0
    }

    #[must_use]
    pub fn is_oscillatory(&self) -> bool {
        self.frequency() > EPSILON
    }
}

/// Matrix and classified modes returned by an analysis.
#[derive(Clone, Debug, Default)]
pub struct EigenmodeResults {
    pub modes: Vec<Eigenmode>,
    pub system_matrix: MatX,
    pub state_count: usize,
    pub reference_velocity: Scalar,
    pub reference_span: Scalar,
    pub reference_chord: Scalar,
}

impl EigenmodeResults {
    pub fn modes_of_type(&self, mode_type: ModeType) -> impl Iterator<Item = &Eigenmode> {
        self.modes
            .iter()
            .filter(move |mode| mode.mode_type == mode_type)
    }

    #[must_use]
    pub fn short_period(&self) -> Option<&Eigenmode> {
        self.modes_of_type(ModeType::ShortPeriod).next()
    }

    #[must_use]
    pub fn phugoid(&self) -> Option<&Eigenmode> {
        self.modes_of_type(ModeType::Phugoid).next()
    }

    #[must_use]
    pub fn dutch_roll(&self) -> Option<&Eigenmode> {
        self.modes_of_type(ModeType::DutchRoll).next()
    }

    #[must_use]
    pub fn is_stable(&self) -> bool {
        self.modes.iter().all(Eigenmode::is_stable)
    }
}

/// Builds longitudinal, lateral, and full rigid-body systems.
#[derive(Clone, Debug, Default)]
pub struct EigenmodeAnalysis;

impl EigenmodeAnalysis {
    #[must_use]
    pub fn new() -> Self {
        Self
    }

    pub fn analyze_longitudinal(
        &self,
        aircraft: &Aircraft,
        mass: &MassProperties,
        derivatives: &StabilityDerivatives,
        flight: &FlightState,
    ) -> Result<EigenmodeResults> {
        validate_inputs(mass, derivatives, flight)?;
        let matrix = self.longitudinal_matrix(aircraft, mass, derivatives, flight);
        self.analyze_matrix(matrix, derivatives, flight, true)
    }

    pub fn analyze_lateral(
        &self,
        aircraft: &Aircraft,
        mass: &MassProperties,
        derivatives: &StabilityDerivatives,
        flight: &FlightState,
    ) -> Result<EigenmodeResults> {
        validate_inputs(mass, derivatives, flight)?;
        let matrix = self.lateral_matrix(aircraft, mass, derivatives, flight);
        self.analyze_matrix(matrix, derivatives, flight, false)
    }

    pub fn analyze_full(
        &self,
        aircraft: &Aircraft,
        mass: &MassProperties,
        derivatives: &StabilityDerivatives,
        flight: &FlightState,
    ) -> Result<EigenmodeResults> {
        validate_inputs(mass, derivatives, flight)?;
        let matrix = self.full_matrix(aircraft, mass, derivatives, flight);
        self.analyze_matrix(matrix, derivatives, flight, false)
    }

    #[must_use]
    pub fn longitudinal_matrix(
        &self,
        aircraft: &Aircraft,
        mass: &MassProperties,
        derivatives: &StabilityDerivatives,
        flight: &FlightState,
    ) -> MatX {
        let velocity = flight.velocity;
        let effective_mass = mass.mass;
        let pitch_inertia = mass.i_yy;
        let qs = flight.dynamic_pressure() * aircraft.reference_area();
        let force_factor = qs / (effective_mass * velocity);
        let moment_factor = qs * aircraft.reference_chord() / pitch_inertia;
        let mut matrix = MatX::zeros(4, 4);
        matrix[(0, 0)] = -force_factor * 2.0 * derivatives.by_alpha.drag;
        matrix[(0, 1)] = force_factor * derivatives.by_alpha.lift;
        matrix[(0, 3)] = -flight.gravity;
        matrix[(1, 0)] = -force_factor * derivatives.by_alpha.lift;
        matrix[(1, 1)] = -force_factor * derivatives.by_alpha.lift;
        matrix[(1, 2)] = velocity;
        matrix[(2, 1)] = moment_factor * derivatives.by_alpha.pitch / velocity;
        matrix[(2, 2)] =
            moment_factor * derivatives.by_pitch_rate.pitch * aircraft.reference_chord()
                / (2.0 * velocity);
        matrix[(3, 2)] = 1.0;
        matrix
    }

    #[must_use]
    pub fn lateral_matrix(
        &self,
        aircraft: &Aircraft,
        mass: &MassProperties,
        derivatives: &StabilityDerivatives,
        flight: &FlightState,
    ) -> MatX {
        let velocity = flight.velocity;
        let effective_mass = mass.mass;
        let roll_inertia = mass.i_xx;
        let yaw_inertia = mass.i_zz;
        let span = aircraft.reference_span();
        let qs = flight.dynamic_pressure() * aircraft.reference_area();
        let force_factor = qs / (effective_mass * velocity);
        let roll_factor = qs * span / roll_inertia;
        let yaw_factor = qs * span / yaw_inertia;
        let mut matrix = MatX::zeros(4, 4);
        matrix[(0, 0)] = force_factor * derivatives.by_beta.side;
        matrix[(0, 1)] = force_factor * derivatives.by_roll_rate.side * span / (2.0 * velocity);
        matrix[(0, 2)] =
            force_factor * derivatives.by_yaw_rate.side * span / (2.0 * velocity) - velocity;
        matrix[(0, 3)] = flight.gravity * flight.theta.cos();
        matrix[(1, 0)] = roll_factor * derivatives.by_beta.roll / velocity;
        matrix[(1, 1)] = roll_factor * derivatives.by_roll_rate.roll * span / (2.0 * velocity);
        matrix[(1, 2)] = roll_factor * derivatives.by_yaw_rate.roll * span / (2.0 * velocity);
        matrix[(2, 0)] = yaw_factor * derivatives.by_beta.yaw / velocity;
        matrix[(2, 1)] = yaw_factor * derivatives.by_roll_rate.yaw * span / (2.0 * velocity);
        matrix[(2, 2)] = yaw_factor * derivatives.by_yaw_rate.yaw * span / (2.0 * velocity);
        matrix[(3, 1)] = 1.0;
        matrix
    }

    /// Full `[u,v,w,p,q,r,φ,θ,ψ,x,y,z]` linear system.
    #[must_use]
    pub fn full_matrix(
        &self,
        aircraft: &Aircraft,
        mass: &MassProperties,
        derivatives: &StabilityDerivatives,
        flight: &FlightState,
    ) -> MatX {
        let mut matrix = MatX::zeros(12, 12);
        let longitudinal = self.longitudinal_matrix(aircraft, mass, derivatives, flight);
        for (reduced_row, full_row) in [0, 2, 4, 7].iter().enumerate() {
            for (reduced_column, full_column) in [0, 2, 4, 7].iter().enumerate() {
                matrix[(*full_row, *full_column)] = longitudinal[(reduced_row, reduced_column)];
            }
        }
        let lateral = self.lateral_matrix(aircraft, mass, derivatives, flight);
        for (reduced_row, full_row) in [1, 3, 5, 6].iter().enumerate() {
            for (reduced_column, full_column) in [1, 3, 5, 6].iter().enumerate() {
                matrix[(*full_row, *full_column)] = lateral[(reduced_row, reduced_column)];
            }
        }
        matrix[(8, 5)] = 1.0 / flight.theta.cos();
        let (sin_phi, cos_phi) = flight.phi.sin_cos();
        let (sin_theta, cos_theta) = flight.theta.sin_cos();
        let rotation = [
            [cos_theta, sin_phi * sin_theta, cos_phi * sin_theta],
            [0.0, cos_phi, -sin_phi],
            [-sin_theta, sin_phi * cos_theta, cos_phi * cos_theta],
        ];
        for row in 0..3 {
            for column in 0..3 {
                matrix[(9 + row, column)] = rotation[row][column];
            }
        }
        matrix
    }

    fn analyze_matrix(
        &self,
        system_matrix: MatX,
        derivatives: &StabilityDerivatives,
        flight: &FlightState,
        longitudinal: bool,
    ) -> Result<EigenmodeResults> {
        let n = system_matrix.nrows();
        let faer_matrix = Mat::<Scalar>::from_fn(n, n, |row, column| system_matrix[(row, column)]);
        let decomposition = faer_matrix
            .eigen()
            .map_err(|error| ZephyrError::Eigensolver(format!("{error:?}")))?;
        let eigenvalues = decomposition.S().column_vector();
        let eigenvectors = decomposition.U();
        let mut modes = Vec::new();
        for column in 0..n {
            let value = eigenvalues[column];
            let eigenvalue = Complex64::new(value.re, value.im);
            if eigenvalue.im < -EPSILON {
                continue;
            }
            let mut vector: Vec<_> = (0..n)
                .map(|row| {
                    let value = eigenvectors[(row, column)];
                    Complex64::new(value.re, value.im)
                })
                .collect();
            let norm = vector
                .iter()
                .map(|value| value.norm_sqr())
                .sum::<Scalar>()
                .sqrt();
            if norm > EPSILON {
                for value in &mut vector {
                    *value /= norm;
                }
            }
            let mut mode = Eigenmode {
                eigenvalue,
                eigenvector: vector,
                mode_type: ModeType::Unknown,
            };
            mode.mode_type = classify(&mode, longitudinal);
            modes.push(mode);
        }
        modes.sort_by(|left, right| {
            left.damping()
                .total_cmp(&right.damping())
                .then_with(|| right.frequency().total_cmp(&left.frequency()))
                .then_with(|| left.mode_type.cmp(&right.mode_type))
        });
        Ok(EigenmodeResults {
            state_count: n,
            system_matrix,
            reference_velocity: flight.velocity,
            reference_span: derivatives.reference_span,
            reference_chord: derivatives.reference_chord,
            modes,
        })
    }
}

fn validate_inputs(
    mass: &MassProperties,
    derivatives: &StabilityDerivatives,
    flight: &FlightState,
) -> Result<()> {
    if !mass.is_valid() {
        return Err(ZephyrError::Eigensolver(
            "mass and principal inertias must be positive".into(),
        ));
    }
    if !derivatives.is_valid {
        return Err(ZephyrError::Eigensolver(
            "valid stability derivatives are required".into(),
        ));
    }
    if flight.velocity <= 0.0 || !flight.velocity.is_finite() || flight.density <= 0.0 {
        return Err(ZephyrError::Eigensolver(
            "positive finite velocity and density are required".into(),
        ));
    }
    Ok(())
}

fn classify(mode: &Eigenmode, longitudinal: bool) -> ModeType {
    if longitudinal {
        if mode.is_oscillatory() {
            if mode.frequency() > 1.0 && mode.damping_ratio() > 0.1 {
                ModeType::ShortPeriod
            } else if mode.frequency() < 1.0 {
                ModeType::Phugoid
            } else {
                ModeType::Other
            }
        } else {
            ModeType::Other
        }
    } else if mode.is_oscillatory() {
        ModeType::DutchRoll
    } else if mode.damping() < -0.5 {
        ModeType::RollSubsidence
    } else {
        ModeType::Spiral
    }
}
