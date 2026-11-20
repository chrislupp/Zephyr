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

//! Airfoil camber-line representation.

use std::fs;
use std::path::Path;

use crate::core::Scalar;
use crate::error::{Result, ZephyrError};
use crate::utils::CubicSpline;

/// Airfoil camber and thickness distributions used by VLM meshing.
#[derive(Clone, Debug, Default)]
pub struct Airfoil {
    name: String,
    camber_x: Vec<Scalar>,
    camber_y: Vec<Scalar>,
    thickness: Vec<Scalar>,
    camber_spline: CubicSpline,
    thickness_spline: CubicSpline,
    max_camber: Scalar,
    max_camber_position: Scalar,
    max_thickness: Scalar,
}

impl Airfoil {
    /// Construct a zero-camber, zero-thickness flat plate.
    #[must_use]
    pub fn flat_plate() -> Self {
        Self::from_parts(
            "Flat Plate".into(),
            vec![0.0, 1.0],
            vec![0.0, 0.0],
            vec![0.0, 0.0],
        )
        .expect("two-point flat plate is valid")
    }

    /// Construct a four-digit NACA airfoil from a numeric designation.
    pub fn naca4(code: u16, point_count: usize) -> Result<Self> {
        if code > 9999 {
            return Err(ZephyrError::InvalidGeometry(
                "NACA four-digit code must be in 0000..=9999".into(),
            ));
        }
        if point_count < 2 {
            return Err(ZephyrError::InvalidGeometry(
                "a NACA airfoil needs at least two points".into(),
            ));
        }
        let m_digit = code / 1000;
        let p_digit = (code % 1000) / 100;
        let thickness_digits = code % 100;
        let max_camber = Scalar::from(m_digit) / 100.0;
        let camber_position = Scalar::from(p_digit) / 10.0;
        let max_thickness = Scalar::from(thickness_digits) / 100.0;

        let mut x_values = Vec::with_capacity(point_count);
        let mut camber = Vec::with_capacity(point_count);
        let mut thickness = Vec::with_capacity(point_count);
        for i in 0..point_count {
            let x = i as Scalar / (point_count - 1) as Scalar;
            let y_c = if max_camber > 0.0 && camber_position > 0.0 {
                if x < camber_position {
                    max_camber / camber_position.powi(2) * (2.0 * camber_position * x - x * x)
                } else {
                    let q = 1.0 - camber_position;
                    max_camber / q.powi(2)
                        * (1.0 - 2.0 * camber_position + 2.0 * camber_position * x - x * x)
                }
            } else {
                0.0
            };
            let t = 5.0
                * max_thickness
                * (0.2969 * x.sqrt() - 0.1260 * x - 0.3516 * x.powi(2) + 0.2843 * x.powi(3)
                    - 0.1015 * x.powi(4));
            x_values.push(x);
            camber.push(y_c);
            thickness.push(t);
        }
        let mut airfoil = Self::from_parts(format!("NACA {code}"), x_values, camber, thickness)?;
        airfoil.max_camber = max_camber;
        airfoil.max_camber_position = camber_position;
        airfoil.max_thickness = max_thickness;
        Ok(airfoil)
    }

    /// Construct a four-digit NACA airfoil from text such as `"2412"`.
    pub fn naca4_str(code: &str, point_count: usize) -> Result<Self> {
        if code.len() > 4 || code.is_empty() || !code.bytes().all(|byte| byte.is_ascii_digit()) {
            return Err(ZephyrError::InvalidGeometry(format!(
                "invalid NACA four-digit code '{code}'"
            )));
        }
        let parsed = code.parse::<u16>().map_err(|_| {
            ZephyrError::InvalidGeometry(format!("invalid NACA four-digit code '{code}'"))
        })?;
        Self::naca4(parsed, point_count)
    }

    /// Read Selig/Lednicer-style coordinates from a file.
    pub fn from_coordinate_file(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let content = fs::read_to_string(path).map_err(|error| ZephyrError::io(path, error))?;
        let mut name = None;
        let mut x = Vec::new();
        let mut y = Vec::new();
        for (line_index, line) in content.lines().enumerate() {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            let values: Vec<_> = trimmed.split_whitespace().collect();
            if values.len() >= 2
                && let (Ok(x_value), Ok(y_value)) =
                    (values[0].parse::<Scalar>(), values[1].parse::<Scalar>())
            {
                x.push(x_value);
                y.push(y_value);
                continue;
            }
            if name.is_none() && trimmed.chars().any(char::is_alphabetic) && line_index < 2 {
                name = Some(trimmed.to_owned());
            }
        }
        if x.len() < 3 {
            return Err(ZephyrError::InvalidGeometry(format!(
                "insufficient airfoil coordinates in {}",
                path.display()
            )));
        }
        let mut airfoil = Self::from_coordinates(&x, &y, false)?;
        airfoil.name = name.unwrap_or_else(|| path.display().to_string());
        Ok(airfoil)
    }

    /// Construct from either an explicit camber line or closed airfoil
    /// coordinates.
    pub fn from_coordinates(x: &[Scalar], y: &[Scalar], is_camber_line: bool) -> Result<Self> {
        if x.len() != y.len() {
            return Err(ZephyrError::InvalidGeometry(
                "airfoil X and Y arrays must have equal length".into(),
            ));
        }
        if x.len() < 2 {
            return Err(ZephyrError::InvalidGeometry(
                "an airfoil needs at least two coordinates".into(),
            ));
        }
        let (camber_x, camber_y, thickness) = if is_camber_line {
            (x.to_vec(), y.to_vec(), vec![0.0; x.len()])
        } else {
            extract_camber(x, y)?
        };
        Self::from_parts(String::new(), camber_x, camber_y, thickness)
    }

    fn from_parts(
        name: String,
        camber_x: Vec<Scalar>,
        camber_y: Vec<Scalar>,
        thickness: Vec<Scalar>,
    ) -> Result<Self> {
        let camber_spline = CubicSpline::build(&camber_x, &camber_y)?;
        let thickness_spline = CubicSpline::build(&camber_x, &thickness)?;
        let mut max_camber: Scalar = 0.0;
        let mut max_camber_position = 0.0;
        for (&x, &y) in camber_x.iter().zip(&camber_y) {
            if y.abs() > max_camber.abs() {
                max_camber = y;
                max_camber_position = x;
            }
        }
        let max_thickness = thickness.iter().copied().fold(0.0_f64, Scalar::max);
        Ok(Self {
            name,
            camber_x,
            camber_y,
            thickness,
            camber_spline,
            thickness_spline,
            max_camber,
            max_camber_position,
            max_thickness,
        })
    }

    /// Evaluate camber height, clamping the chord fraction to `[0,1]`.
    #[must_use]
    pub fn camber(&self, x: Scalar) -> Scalar {
        self.camber_spline
            .evaluate(x.clamp(0.0, 1.0))
            .unwrap_or(0.0)
    }

    /// Evaluate camber slope.
    #[must_use]
    pub fn camber_slope(&self, x: Scalar) -> Scalar {
        self.camber_spline
            .evaluate_derivative(x.clamp(0.0, 1.0))
            .unwrap_or(0.0)
    }

    /// Derivative of camber slope, including the chord-fraction clamp.
    pub(crate) fn camber_slope_derivative(&self, x: Scalar) -> Scalar {
        if !(0.0..=1.0).contains(&x) {
            return 0.0;
        }
        self.camber_spline
            .evaluate_second_derivative(x)
            .unwrap_or(0.0)
    }

    /// Evaluate thickness.
    #[must_use]
    pub fn thickness_at(&self, x: Scalar) -> Scalar {
        self.thickness_spline
            .evaluate(x.clamp(0.0, 1.0))
            .unwrap_or(0.0)
    }

    /// Whether this airfoil was initialized.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        !self.camber_x.is_empty()
    }

    /// Name or designation.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Set the display name.
    pub fn set_name(&mut self, name: impl Into<String>) {
        self.name = name.into();
    }

    /// Camber-line chord coordinates.
    #[must_use]
    pub fn camber_x(&self) -> &[Scalar] {
        &self.camber_x
    }

    /// Camber-line heights.
    #[must_use]
    pub fn camber_y(&self) -> &[Scalar] {
        &self.camber_y
    }

    /// Thickness values at [`Self::camber_x`] coordinates.
    #[must_use]
    pub fn thickness_distribution(&self) -> &[Scalar] {
        &self.thickness
    }

    /// Maximum signed camber.
    #[must_use]
    pub fn max_camber(&self) -> Scalar {
        self.max_camber
    }

    /// Chord position of maximum absolute camber.
    #[must_use]
    pub fn max_camber_position(&self) -> Scalar {
        self.max_camber_position
    }

    /// Maximum thickness.
    #[must_use]
    pub fn max_thickness(&self) -> Scalar {
        self.max_thickness
    }
}

fn extract_camber(x: &[Scalar], y: &[Scalar]) -> Result<(Vec<Scalar>, Vec<Scalar>, Vec<Scalar>)> {
    let leading_edge = x
        .iter()
        .enumerate()
        .min_by(|(_, left), (_, right)| left.total_cmp(right))
        .map(|(index, _)| index)
        .ok_or_else(|| ZephyrError::InvalidGeometry("empty airfoil coordinates".into()))?;
    let mut upper_x = x[..=leading_edge].to_vec();
    let mut upper_y = y[..=leading_edge].to_vec();
    upper_x.reverse();
    upper_y.reverse();
    let lower_x = x[leading_edge..].to_vec();
    let lower_y = y[leading_edge..].to_vec();
    let upper = if upper_x.len() >= 2 {
        Some(CubicSpline::build(&upper_x, &upper_y)?)
    } else {
        None
    };
    let lower = if lower_x.len() >= 2 {
        Some(CubicSpline::build(&lower_x, &lower_y)?)
    } else {
        None
    };
    let count = 50;
    let mut camber_x = Vec::with_capacity(count);
    let mut camber_y = Vec::with_capacity(count);
    let mut thickness = Vec::with_capacity(count);
    for index in 0..count {
        let chord = index as Scalar / (count - 1) as Scalar;
        let upper_y = upper
            .as_ref()
            .and_then(|spline| spline.evaluate(chord).ok())
            .unwrap_or(0.0);
        let lower_y = lower
            .as_ref()
            .and_then(|spline| spline.evaluate(chord).ok())
            .unwrap_or(0.0);
        camber_x.push(chord);
        camber_y.push(0.5 * (upper_y + lower_y));
        thickness.push(upper_y - lower_y);
    }
    Ok((camber_x, camber_y, thickness))
}
