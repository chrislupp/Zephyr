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

//! Cubic Hermite spline compatible with AVL's `spline.f`.

use crate::core::Scalar;
use crate::error::{Result, ZephyrError};

/// Endpoint boundary condition.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum SplineBoundaryCondition {
    /// Zero second derivative.
    #[default]
    Natural,
    /// Zero third derivative.
    ZeroThirdDerivative,
    /// Prescribed first derivative.
    SpecifiedSlope,
}

/// One-dimensional cubic spline.
#[derive(Clone, Debug, Default)]
pub struct CubicSpline {
    parameters: Vec<Scalar>,
    values: Vec<Scalar>,
    slopes: Vec<Scalar>,
}

impl CubicSpline {
    /// Build a natural spline.
    pub fn build(parameters: &[Scalar], values: &[Scalar]) -> Result<Self> {
        Self::build_with_boundaries(
            parameters,
            values,
            SplineBoundaryCondition::Natural,
            SplineBoundaryCondition::Natural,
            0.0,
            0.0,
        )
    }

    /// Build a spline with explicit endpoint conditions.
    pub fn build_with_boundaries(
        parameters: &[Scalar],
        values: &[Scalar],
        left: SplineBoundaryCondition,
        right: SplineBoundaryCondition,
        left_slope: Scalar,
        right_slope: Scalar,
    ) -> Result<Self> {
        if parameters.len() != values.len() {
            return Err(ZephyrError::InvalidGeometry(
                "spline parameters and values must have equal length".into(),
            ));
        }
        if parameters.len() < 2 {
            return Err(ZephyrError::InvalidGeometry(
                "a spline needs at least two points".into(),
            ));
        }
        if parameters.windows(2).any(|window| window[1] <= window[0]) {
            return Err(ZephyrError::InvalidGeometry(
                "spline parameters must be strictly increasing".into(),
            ));
        }

        let n = parameters.len();
        let mut diagonal = vec![0.0; n];
        let mut sub = vec![0.0; n];
        let mut sup = vec![0.0; n];
        let mut slopes = vec![0.0; n];

        for i in 1..n - 1 {
            let ds_minus = parameters[i] - parameters[i - 1];
            let ds_plus = parameters[i + 1] - parameters[i];
            sub[i] = ds_plus;
            diagonal[i] = 2.0 * (ds_minus + ds_plus);
            sup[i] = ds_minus;
            slopes[i] = 3.0
                * ((values[i + 1] - values[i]) * ds_minus / ds_plus
                    + (values[i] - values[i - 1]) * ds_plus / ds_minus);
        }

        let ds_left = parameters[1] - parameters[0];
        match left {
            SplineBoundaryCondition::Natural => {
                diagonal[0] = 2.0;
                sup[0] = 1.0;
                slopes[0] = 3.0 * (values[1] - values[0]) / ds_left;
            }
            SplineBoundaryCondition::ZeroThirdDerivative => {
                diagonal[0] = 1.0;
                sup[0] = 1.0;
                slopes[0] = 2.0 * (values[1] - values[0]) / ds_left;
            }
            SplineBoundaryCondition::SpecifiedSlope => {
                diagonal[0] = 1.0;
                slopes[0] = left_slope;
            }
        }

        let ds_right = parameters[n - 1] - parameters[n - 2];
        match right {
            SplineBoundaryCondition::Natural => {
                sub[n - 1] = 1.0;
                diagonal[n - 1] = 2.0;
                slopes[n - 1] = 3.0 * (values[n - 1] - values[n - 2]) / ds_right;
            }
            SplineBoundaryCondition::ZeroThirdDerivative => {
                sub[n - 1] = 1.0;
                diagonal[n - 1] = 1.0;
                slopes[n - 1] = 2.0 * (values[n - 1] - values[n - 2]) / ds_right;
            }
            SplineBoundaryCondition::SpecifiedSlope => {
                diagonal[n - 1] = 1.0;
                slopes[n - 1] = right_slope;
            }
        }
        if n == 2
            && left == SplineBoundaryCondition::ZeroThirdDerivative
            && right == SplineBoundaryCondition::ZeroThirdDerivative
        {
            sub[1] = 1.0;
            diagonal[1] = 2.0;
            slopes[1] = 3.0 * (values[1] - values[0]) / ds_right;
        }
        solve_tridiagonal(&mut diagonal, &sub, &mut sup, &mut slopes);

        Ok(Self {
            parameters: parameters.to_vec(),
            values: values.to_vec(),
            slopes,
        })
    }

    /// Whether the spline contains usable data.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        !self.parameters.is_empty()
    }

    /// Number of control points.
    #[must_use]
    pub fn len(&self) -> usize {
        self.parameters.len()
    }

    /// Whether the spline has no control points.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.parameters.is_empty()
    }

    /// Smallest parameter, or zero for an empty spline.
    #[must_use]
    pub fn min_parameter(&self) -> Scalar {
        self.parameters.first().copied().unwrap_or(0.0)
    }

    /// Largest parameter, or zero for an empty spline.
    #[must_use]
    pub fn max_parameter(&self) -> Scalar {
        self.parameters.last().copied().unwrap_or(0.0)
    }

    /// Evaluate value, first derivative, and second derivative.
    pub fn evaluate_all(&self, parameter: Scalar) -> Result<(Scalar, Scalar, Scalar)> {
        if !self.is_valid() {
            return Err(ZephyrError::InvalidGeometry(
                "spline has not been built".into(),
            ));
        }
        let i = self.find_interval(parameter);
        let ds = self.parameters[i] - self.parameters[i - 1];
        let t = (parameter - self.parameters[i - 1]) / ds;
        let f0 = self.values[i - 1];
        let f1 = ds * self.slopes[i - 1];
        let delta = self.values[i] - self.values[i - 1];
        let f2 = -ds * (2.0 * self.slopes[i - 1] + self.slopes[i]) + 3.0 * delta;
        let f3 = ds * (self.slopes[i - 1] + self.slopes[i]) - 2.0 * delta;
        Ok((
            f0 + t * (f1 + t * (f2 + t * f3)),
            (f1 + t * (2.0 * f2 + t * 3.0 * f3)) / ds,
            (2.0 * f2 + t * 6.0 * f3) / (ds * ds),
        ))
    }

    /// Evaluate the interpolated value.
    pub fn evaluate(&self, parameter: Scalar) -> Result<Scalar> {
        self.evaluate_all(parameter).map(|all| all.0)
    }

    /// Evaluate the first derivative.
    pub fn evaluate_derivative(&self, parameter: Scalar) -> Result<Scalar> {
        self.evaluate_all(parameter).map(|all| all.1)
    }

    /// Evaluate the second derivative.
    pub fn evaluate_second_derivative(&self, parameter: Scalar) -> Result<Scalar> {
        self.evaluate_all(parameter).map(|all| all.2)
    }

    /// Index `i` such that the parameter lies in `[parameters[i-1], parameters[i]]`.
    ///
    /// The result is always in `1..=last`, so callers may index `i - 1`
    /// unconditionally. A NaN parameter satisfies neither ordered comparison
    /// and makes `partition_point` return 0, which would underflow that
    /// subtraction; clamping to 1 instead lets the NaN flow through the
    /// polynomial and out as a NaN result, which is what the legacy binary
    /// search did.
    fn find_interval(&self, parameter: Scalar) -> usize {
        if parameter <= self.parameters[0] {
            return 1;
        }
        let last = self.parameters.len() - 1;
        if parameter >= self.parameters[last] {
            return last;
        }
        self.parameters
            .partition_point(|value| *value <= parameter)
            .clamp(1, last)
    }
}

fn solve_tridiagonal(
    diagonal: &mut [Scalar],
    sub: &[Scalar],
    sup: &mut [Scalar],
    rhs: &mut [Scalar],
) {
    for k in 1..diagonal.len() {
        sup[k - 1] /= diagonal[k - 1];
        rhs[k - 1] /= diagonal[k - 1];
        diagonal[k] -= sub[k] * sup[k - 1];
        rhs[k] -= sub[k] * rhs[k - 1];
    }
    let last = diagonal.len() - 1;
    rhs[last] /= diagonal[last];
    for k in (0..last).rev() {
        rhs[k] -= sup[k] * rhs[k + 1];
    }
}
