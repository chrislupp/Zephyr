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

//! AVL `.mass` reader.

use std::fs;
use std::path::Path;

use crate::core::{Mat3, Scalar, Vec3};
use crate::dynamics::MassProperties;
use crate::error::{Result, ZephyrError};

/// One component from an AVL mass file.
#[derive(Clone, Debug, Default)]
pub struct MassComponent {
    /// Component description from its trailing comment.
    pub name: String,
    /// Component mass in file units.
    pub mass: Scalar,
    /// Component position in file units after scale/offset.
    pub position: Vec3,
    /// Roll inertia.
    pub i_xx: Scalar,
    /// Pitch inertia.
    pub i_yy: Scalar,
    /// Yaw inertia.
    pub i_zz: Scalar,
    /// XY product of inertia.
    pub i_xy: Scalar,
    /// XZ product of inertia.
    pub i_xz: Scalar,
    /// YZ product of inertia.
    pub i_yz: Scalar,
}

/// Parsed mass-file data and units.
#[derive(Clone, Debug)]
pub struct MassFileData {
    /// Length-to-SI multiplier.
    pub length_unit: Scalar,
    /// Mass-to-SI multiplier.
    pub mass_unit: Scalar,
    /// Time-to-SI multiplier.
    pub time_unit: Scalar,
    /// Default gravity.
    pub gravity: Scalar,
    /// Default density.
    pub density: Scalar,
    /// Parsed components.
    pub components: Vec<MassComponent>,
}

impl Default for MassFileData {
    fn default() -> Self {
        Self {
            length_unit: 1.0,
            mass_unit: 1.0,
            time_unit: 1.0,
            gravity: 9.81,
            density: 1.225,
            components: Vec::new(),
        }
    }
}

impl MassFileData {
    /// Total mass in SI units.
    #[must_use]
    pub fn total_mass(&self) -> Scalar {
        self.components
            .iter()
            .map(|component| component.mass * self.mass_unit)
            .sum()
    }

    /// Combined center of gravity in SI units.
    #[must_use]
    pub fn center_of_gravity(&self) -> Vec3 {
        let total = self.total_mass();
        if total <= 0.0 {
            return Vec3::zeros();
        }
        self.components
            .iter()
            .fold(Vec3::zeros(), |sum, component| {
                sum + component.mass * self.mass_unit * component.position * self.length_unit
            })
            / total
    }

    /// Combine components with the parallel-axis theorem.
    #[must_use]
    pub fn total_properties(&self) -> MassProperties {
        let total_mass = self.total_mass();
        if total_mass <= 0.0 {
            return MassProperties::default();
        }
        let center_of_gravity = self.center_of_gravity();
        let inertia_unit = self.mass_unit * self.length_unit * self.length_unit;
        let mut inertia = Mat3::zeros();
        for component in &self.components {
            let mass = component.mass * self.mass_unit;
            let position = component.position * self.length_unit;
            let displacement = position - center_of_gravity;
            let own = inertia_unit
                * Mat3::new(
                    component.i_xx,
                    -component.i_xy,
                    -component.i_xz,
                    -component.i_xy,
                    component.i_yy,
                    -component.i_yz,
                    -component.i_xz,
                    -component.i_yz,
                    component.i_zz,
                );
            let parallel = mass
                * (displacement.norm_squared() * Mat3::identity()
                    - displacement * displacement.transpose());
            inertia += own + parallel;
        }
        let mut properties = MassProperties {
            mass: total_mass,
            center_of_gravity,
            ..MassProperties::default()
        };
        properties.set_inertia_tensor(&inertia);
        properties
    }
}

/// Parser for AVL mass-property text.
#[derive(Clone, Copy, Debug, Default)]
pub struct MassReader;

impl MassReader {
    /// Read a mass file.
    pub fn read_path(path: impl AsRef<Path>) -> Result<MassFileData> {
        let path = path.as_ref();
        let content = fs::read_to_string(path).map_err(|error| ZephyrError::io(path, error))?;
        Self::read_str(&content)
    }

    /// Parse mass-file text.
    pub fn read_str(content: &str) -> Result<MassFileData> {
        let mut data = MassFileData::default();
        let mut scale = Vec3::repeat(1.0);
        let mut offset = Vec3::zeros();
        for line in content.lines() {
            parse_mass_line(line, &mut data, &mut scale, &mut offset);
        }
        Ok(data)
    }
}

/// Parses one mass-file line.
///
/// Comments are stripped before any dispatch decision, matching the legacy
/// reader's `RemoveComment`-then-classify order. Getting that order wrong is
/// not cosmetic: an `=` inside a trailing comment (`! main gear (x = 1.2 m)`)
/// would otherwise capture the line as a unit assignment and discard a real
/// mass component, and a number inside a comment would be read as a scale
/// factor.
///
/// Each branch also mirrors the legacy extraction semantics: a unit line with
/// an unrecognized key falls through to component parsing rather than
/// consuming the line, `*`/`+` require three consecutive numeric fields and
/// leave the previous value untouched otherwise, and numeric runs stop at the
/// first token that does not parse instead of skipping it.
fn parse_mass_line(line: &str, data: &mut MassFileData, scale: &mut Vec3, offset: &mut Vec3) {
    let trimmed = remove_comment(line);
    if trimmed.is_empty() {
        return;
    }
    if let Some((key, value)) = trimmed.split_once('=') {
        let key = key.trim().to_ascii_lowercase();
        let parsed = value
            .split_whitespace()
            .next()
            .and_then(|value| value.parse().ok());
        if let Some(value) = parsed {
            // Only a recognized key consumes the line. An unknown key falls
            // through, as `ParseUnitLine` returning false did.
            match key.as_str() {
                "lunit" => data.length_unit = value,
                "munit" => data.mass_unit = value,
                "tunit" => data.time_unit = value,
                "g" => data.gravity = value,
                "rho" => data.density = value,
                _ => {}
            }
            if matches!(key.as_str(), "lunit" | "munit" | "tunit" | "g" | "rho") {
                return;
            }
        }
    }
    if let Some(values) = trimmed.strip_prefix('*') {
        if let Some(parsed) = parse_vec3(values) {
            *scale = parsed;
        }
        return;
    }
    if let Some(values) = trimmed.strip_prefix('+') {
        if let Some(parsed) = parse_vec3(values) {
            *offset = parsed;
        }
        return;
    }
    // Keyword lines begin with a letter and are not components.
    if trimmed.starts_with(|c: char| c.is_alphabetic()) {
        return;
    }
    let values: Vec<Scalar> = trimmed
        .split_whitespace()
        .map_while(|value| value.parse().ok())
        .collect();
    if values.len() < 7 {
        return;
    }
    data.components.push(MassComponent {
        name: extract_comment(line).to_owned(),
        mass: values[0],
        position: Vec3::new(
            values[1] * scale.x + offset.x,
            values[2] * scale.y + offset.y,
            values[3] * scale.z + offset.z,
        ),
        i_xx: values[4],
        i_yy: values[5],
        i_zz: values[6],
        i_xy: values.get(7).copied().unwrap_or(0.0),
        i_xz: values.get(8).copied().unwrap_or(0.0),
        i_yz: values.get(9).copied().unwrap_or(0.0),
    });
}

fn remove_comment(line: &str) -> &str {
    line.split(['!', '#']).next().unwrap_or(line).trim()
}

fn extract_comment(line: &str) -> &str {
    line.find(['!', '#'])
        .map_or("", |index| line[index + 1..].trim())
}

/// Reads three consecutive numeric fields.
///
/// Stops at the first token that does not parse rather than skipping it, so a
/// trailing word cannot promote a later number into one of the three slots.
fn parse_vec3(text: &str) -> Option<Vec3> {
    let mut values = text
        .split_whitespace()
        .map_while(|value| value.parse::<Scalar>().ok());
    Some(Vec3::new(values.next()?, values.next()?, values.next()?))
}
