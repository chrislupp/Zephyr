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

//! Athena Vortex Lattice geometry reader.

use std::fs;
use std::path::{Path, PathBuf};

use crate::core::{Scalar, SymmetryType, Vec3};
use crate::error::{Result, ZephyrError};
use crate::geometry::{
    Aircraft, Airfoil, Body, BodySection, ControlSurfaceDefinition, DesignVariableDefinition,
    Section, Surface,
};

/// Parser for AVL geometry files.
#[derive(Clone, Copy, Debug, Default)]
pub struct AvlReader {
    verbose: bool,
}

impl AvlReader {
    /// Construct a reader.
    #[must_use]
    pub const fn new() -> Self {
        Self { verbose: false }
    }

    /// Enable or disable progress messages.
    pub fn set_verbose(&mut self, verbose: bool) {
        self.verbose = verbose;
    }

    /// Read an AVL file, resolving `AFILE` references relative to it.
    pub fn read_path(&self, path: impl AsRef<Path>) -> Result<Aircraft> {
        let path = path.as_ref();
        let content = fs::read_to_string(path).map_err(|error| ZephyrError::io(path, error))?;
        self.read_with_base(&content, path.parent().unwrap_or(Path::new(".")))
            .map_err(|error| error.with_source_path(path))
    }

    /// Parse AVL text without a base directory.
    pub fn read_str(&self, content: &str) -> Result<Aircraft> {
        self.read_with_base(content, Path::new(""))
    }

    /// Parse AVL text with a base directory for relative airfoil files.
    pub fn read_with_base(&self, content: &str, base: &Path) -> Result<Aircraft> {
        let lines: Vec<(usize, &str)> = content
            .lines()
            .enumerate()
            .filter_map(|(index, line)| {
                let trimmed = line.trim();
                (!trimmed.is_empty() && !trimmed.starts_with('#') && !trimmed.starts_with('!'))
                    .then_some((index + 1, trimmed))
            })
            .collect();
        let mut cursor = Cursor::new(lines);
        let mut aircraft = Aircraft::new(cursor.next_required("title")?.1);
        if self.verbose {
            eprintln!("Configuration: {}", aircraft.title());
        }
        if let Some(value) = parse_numbers(cursor.next_required("Mach")?.1).first() {
            aircraft.set_mach_number(*value);
        }
        let symmetry = parse_numbers(cursor.next_required("symmetry")?.1);
        if symmetry.len() >= 3 {
            aircraft.set_y_symmetry(symmetry_type(symmetry[0]));
            aircraft.set_z_symmetry(symmetry_type(symmetry[1]));
            aircraft.set_z_symmetry_plane(symmetry[2]);
        }
        let reference = parse_numbers(cursor.next_required("reference")?.1);
        if reference.len() >= 3 {
            aircraft.set_reference_area(reference[0]);
            aircraft.set_reference_chord(reference[1]);
            aircraft.set_reference_span(reference[2]);
        }
        let moment = parse_numbers(cursor.next_required("moment reference")?.1);
        if moment.len() >= 3 {
            aircraft.set_moment_reference_point(Vec3::new(moment[0], moment[1], moment[2]));
        }
        if let Some((_, line)) = cursor.peek()
            && !matches!(keyword(line).as_str(), "SURF" | "BODY")
        {
            if let Some(value) = parse_numbers(line).first() {
                aircraft.set_profile_drag_coefficient(*value);
            }
            cursor.advance();
        }
        while let Some((line_number, line)) = cursor.next() {
            match keyword(line).as_str() {
                "SURF" => {
                    let surface = self.parse_surface(&mut cursor, base, line_number)?;
                    aircraft.add_surface(surface);
                }
                "BODY" => {
                    let body = self.parse_body(&mut cursor, line_number)?;
                    aircraft.add_body(body);
                }
                _ => {}
            }
        }
        Ok(aircraft)
    }

    fn parse_surface(
        &self,
        cursor: &mut Cursor<'_>,
        base: &Path,
        start_line: usize,
    ) -> Result<Surface> {
        let (_, name) = cursor.next_required_after("SURFACE", start_line)?;
        if self.verbose {
            eprintln!("  Surface: {name}");
        }
        let (_, panel_line) = cursor.next_required_after("surface paneling", start_line)?;
        let paneling = parse_numbers(panel_line);
        let mut surface = Surface::new(name);
        if paneling.len() >= 2 {
            surface.set_chordwise_panel_count(count_from_file(
                paneling[0],
                start_line,
                "chordwise panel",
            )?);
            surface.set_numeric_chordwise_spacing(paneling[1]);
        }
        if paneling.len() >= 4 {
            surface.set_spanwise_panel_count(count_from_file(
                paneling[2],
                start_line,
                "spanwise panel",
            )?);
            surface.set_numeric_spanwise_spacing(paneling[3]);
        }
        let mut section: Option<Section> = None;
        while let Some((line_number, line)) = cursor.peek() {
            let key = keyword(line);
            if matches!(key.as_str(), "SURF" | "BODY") {
                break;
            }
            cursor.advance();
            match key.as_str() {
                "SCAL" => {
                    if let Some(value) = next_vec3(cursor, line_number, "SCALE")? {
                        surface.set_scale(value);
                    }
                }
                "TRAN" => {
                    if let Some(value) = next_vec3(cursor, line_number, "TRANSLATE")? {
                        surface.set_translation(value);
                    }
                }
                "ANGL" => {
                    let values = next_numbers(cursor, line_number, "ANGLE")?;
                    if let Some(value) = values.first() {
                        surface.set_angle_offset_degrees(*value);
                    }
                }
                "YDUP" => {
                    let values = next_numbers(cursor, line_number, "YDUPLICATE")?;
                    if let Some(value) = values.first() {
                        surface.set_y_duplicate(*value);
                    }
                }
                "COMP" => {
                    let values = next_numbers(cursor, line_number, "COMPONENT")?;
                    if let Some(value) = values.first() {
                        surface.set_component_index(count_from_file(
                            *value,
                            line_number,
                            "component index",
                        )?);
                    }
                }
                "NOWA" => surface.set_has_wake(false),
                "NOLO" => surface.set_contributes_to_loads(false),
                "NOAL" => surface.set_affected_by_alpha_beta(false),
                "CDCL" => {
                    let values = next_numbers(cursor, line_number, "CDCL")?;
                    let mut polar = [0.0; 6];
                    for (target, source) in polar.iter_mut().zip(values) {
                        *target = source;
                    }
                    surface.set_cl_cd_polar(polar);
                }
                "SECT" => {
                    if let Some(previous) = section.take() {
                        surface.add_section(previous);
                    }
                    let values = next_numbers(cursor, line_number, "SECTION")?;
                    if values.len() < 5 {
                        return Err(ZephyrError::parse(
                            line_number,
                            "SECTION needs Xle Yle Zle chord incidence",
                        ));
                    }
                    let mut current = Section::new(
                        Vec3::new(values[0], values[1], values[2]),
                        values[3],
                        values[4],
                    );
                    if values.len() >= 7 {
                        current.set_spanwise_panel_count(count_from_file(
                            values[5],
                            line_number,
                            "section spanwise panel",
                        )?);
                        current.set_numeric_spanwise_spacing(values[6]);
                    }
                    section = Some(current);
                }
                "NACA" => {
                    let values = next_numbers(cursor, line_number, "NACA")?;
                    if let (Some(current), Some(code)) = (section.as_mut(), values.first()) {
                        current.set_airfoil(Airfoil::naca4(*code as u16, 50)?);
                    }
                }
                "AFIL" | "AIRF" => {
                    // The filename runs to the end of the line, so a trailing
                    // comment would otherwise become part of the path and the
                    // section would silently keep its default flat plate.
                    let (_, relative) = cursor.next_required_after("AFILE", line_number)?;
                    let path = resolve_path(base, strip_comment(relative));
                    match Airfoil::from_coordinate_file(&path) {
                        Ok(airfoil) => {
                            if let Some(current) = section.as_mut() {
                                current.set_airfoil(airfoil);
                            }
                        }
                        Err(cause) => {
                            // Silently degrading to a flat plate changes the
                            // camber the solve runs on while still reporting
                            // success, so this has to surface.
                            return Err(ZephyrError::parse(
                                line_number,
                                format!("could not load AFILE {}: {cause}", path.display()),
                            ));
                        }
                    }
                }
                "CLAF" => {
                    let values = next_numbers(cursor, line_number, "CLAF")?;
                    if let (Some(current), Some(value)) = (section.as_mut(), values.first()) {
                        current.set_cl_alpha_factor(*value);
                    }
                }
                "CONT" => {
                    let (_, definition) = cursor.next_required_after("CONTROL", line_number)?;
                    if let Some(current) = section.as_mut() {
                        current.add_control_surface(parse_control(definition, line_number)?);
                    }
                }
                "DESI" => {
                    let (_, definition) = cursor.next_required_after("DESIGN", line_number)?;
                    let mut fields = definition.split_whitespace();
                    let name = fields.next().unwrap_or_default().to_owned();
                    let weight = fields
                        .next()
                        .and_then(|value| value.parse().ok())
                        .unwrap_or(0.0);
                    if let Some(current) = section.as_mut() {
                        current.add_design_variable(DesignVariableDefinition { name, weight });
                    }
                }
                _ => {}
            }
        }
        if let Some(section) = section {
            surface.add_section(section);
        }
        Ok(surface)
    }

    fn parse_body(&self, cursor: &mut Cursor<'_>, start_line: usize) -> Result<Body> {
        let (_, name) = cursor.next_required_after("BODY", start_line)?;
        let (_, panel_line) = cursor.next_required_after("body paneling", start_line)?;
        let paneling = parse_numbers(panel_line);
        let mut body = Body::new(name);
        if paneling.len() >= 2 {
            body.set_node_count(count_from_file(paneling[0], start_line, "body node")?);
            body.set_node_spacing(paneling[1]);
        }
        while let Some((line_number, line)) = cursor.peek() {
            let key = keyword(line);
            if matches!(key.as_str(), "SURF" | "BODY") {
                break;
            }
            cursor.advance();
            match key.as_str() {
                "SCAL" => {
                    if let Some(value) = next_vec3(cursor, line_number, "SCALE")? {
                        body.set_scale(value);
                    }
                }
                "TRAN" => {
                    if let Some(value) = next_vec3(cursor, line_number, "TRANSLATE")? {
                        body.set_translation(value);
                    }
                }
                "YDUP" => body.set_y_duplicated(true),
                "BSEC" => {
                    let values = next_numbers(cursor, line_number, "BSEC")?;
                    if values.len() >= 5 {
                        body.add_section(BodySection {
                            center: Vec3::new(values[0], values[1], values[2]),
                            width: values[3],
                            height: values[4],
                        });
                    }
                }
                _ => {}
            }
        }
        Ok(body)
    }
}

struct Cursor<'a> {
    lines: Vec<(usize, &'a str)>,
    index: usize,
}

impl<'a> Cursor<'a> {
    fn new(lines: Vec<(usize, &'a str)>) -> Self {
        Self { lines, index: 0 }
    }

    fn peek(&self) -> Option<(usize, &'a str)> {
        self.lines.get(self.index).copied()
    }

    fn advance(&mut self) {
        self.index += 1;
    }

    fn next(&mut self) -> Option<(usize, &'a str)> {
        let line = self.peek();
        if line.is_some() {
            self.advance();
        }
        line
    }

    fn next_required(&mut self, context: &str) -> Result<(usize, &'a str)> {
        self.next().ok_or_else(|| ZephyrError::Parse {
            location: crate::error::ParseLocation::default(),
            message: format!("unexpected end of AVL data while reading {context}"),
        })
    }

    fn next_required_after(&mut self, context: &str, line: usize) -> Result<(usize, &'a str)> {
        self.next().ok_or_else(|| {
            ZephyrError::parse(line, format!("unexpected end of AVL data after {context}"))
        })
    }
}

fn keyword(line: &str) -> String {
    line.chars()
        .take(4)
        .collect::<String>()
        .to_ascii_uppercase()
}

fn parse_numbers(line: &str) -> Vec<Scalar> {
    line.replace(',', " ")
        .split_whitespace()
        .map_while(|value| value.parse().ok())
        .collect()
}

/// Largest panel or node count accepted from a file.
///
/// Well above any practical lattice, but low enough that the product of two
/// counts cannot overflow the index arithmetic in the mesher or exhaust memory
/// through a single allocation.
const MAX_COUNT_FROM_FILE: Scalar = 100_000.0;

/// Converts a file-supplied count to `usize`, rejecting values the mesher
/// cannot represent.
///
/// `value as usize` saturates in release and is undefined for NaN, so an
/// out-of-range count such as `1e30` would silently become `usize::MAX` and
/// later overflow the panel arithmetic. A library must not abort on file
/// input, so this rejects rather than clamps.
fn count_from_file(value: Scalar, line: usize, context: &str) -> Result<usize> {
    if !value.is_finite() || !(0.0..=MAX_COUNT_FROM_FILE).contains(&value) {
        return Err(ZephyrError::parse(
            line,
            format!("{context} count must be between 0 and {MAX_COUNT_FROM_FILE:.0}"),
        ));
    }
    Ok(value as usize)
}

fn symmetry_type(value: Scalar) -> SymmetryType {
    if value > 0.0 {
        SymmetryType::Symmetric
    } else if value < 0.0 {
        SymmetryType::Antisymmetric
    } else {
        SymmetryType::None
    }
}

fn next_numbers(cursor: &mut Cursor<'_>, line: usize, context: &str) -> Result<Vec<Scalar>> {
    Ok(parse_numbers(cursor.next_required_after(context, line)?.1))
}

fn next_vec3(cursor: &mut Cursor<'_>, line: usize, context: &str) -> Result<Option<Vec3>> {
    let values = next_numbers(cursor, line, context)?;
    Ok((values.len() >= 3).then(|| Vec3::new(values[0], values[1], values[2])))
}

/// Drops an AVL trailing comment and surrounding whitespace.
fn strip_comment(line: &str) -> &str {
    line.split(['!', '#']).next().unwrap_or(line).trim()
}

fn resolve_path(base: &Path, relative: &str) -> PathBuf {
    let relative = Path::new(relative.trim());
    if relative.is_absolute() {
        relative.to_owned()
    } else {
        base.join(relative)
    }
}

/// Parses a `CONTROL` declaration.
///
/// The trailing `SgnDup` field is optional: AVL streams the seven fields and
/// leaves the last at its default when a file supplies only six, which
/// hand-written `.avl` files commonly do. Requiring it would reject the whole
/// aircraft over an omission the legacy reader accepted.
fn parse_control(definition: &str, line: usize) -> Result<ControlSurfaceDefinition> {
    let fields: Vec<_> = definition.split_whitespace().collect();
    if fields.len() < 6 {
        return Err(ZephyrError::parse(
            line,
            "CONTROL needs name gain hinge_x hx hy hz [sign]",
        ));
    }
    let number = |index: usize| {
        fields[index]
            .parse::<Scalar>()
            .map_err(|_| ZephyrError::parse(line, "invalid numeric CONTROL field"))
    };
    let sign = match fields.get(6) {
        Some(_) => number(6)?,
        None => 1.0,
    };
    Ok(ControlSurfaceDefinition {
        name: fields[0].to_owned(),
        gain: number(1)?,
        hinge_x_fraction: number(2)?,
        hinge_vector: Vec3::new(number(3)?, number(4)?, number(5)?),
        sign,
    })
}
