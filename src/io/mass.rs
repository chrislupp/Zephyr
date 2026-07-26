//! AVL `.mass` reader.

use std::fs;
use std::path::Path;

use crate::core::{Mat3, Scalar, Vec3};
use crate::dynamics::MassProperties;
use crate::error::{PegasusError, Result};

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
        let content = fs::read_to_string(path).map_err(|error| PegasusError::io(path, error))?;
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

fn parse_mass_line(line: &str, data: &mut MassFileData, scale: &mut Vec3, offset: &mut Vec3) {
    let trimmed = line.trim();
    if trimmed.is_empty() || trimmed.starts_with('!') || trimmed.starts_with('#') {
        return;
    }
    if let Some((key, value)) = trimmed.split_once('=') {
        let key = key.trim().to_ascii_lowercase();
        let value = remove_comment(value).split_whitespace().next();
        if let Some(value) = value.and_then(|value| value.parse().ok()) {
            match key.as_str() {
                "lunit" => data.length_unit = value,
                "munit" => data.mass_unit = value,
                "tunit" => data.time_unit = value,
                "g" => data.gravity = value,
                "rho" => data.density = value,
                _ => {}
            }
        }
        return;
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
    let data_part = remove_comment(line);
    let values: Vec<Scalar> = data_part
        .split_whitespace()
        .filter_map(|value| value.parse().ok())
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

fn parse_vec3(text: &str) -> Option<Vec3> {
    let mut values = text
        .split_whitespace()
        .filter_map(|value| value.parse().ok());
    Some(Vec3::new(values.next()?, values.next()?, values.next()?))
}
