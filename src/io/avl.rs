//! Athena Vortex Lattice geometry reader.

use std::fs;
use std::path::{Path, PathBuf};

use crate::core::{Scalar, SymmetryType, Vec3};
use crate::error::{PegasusError, Result};
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
        let content = fs::read_to_string(path).map_err(|error| PegasusError::io(path, error))?;
        self.read_with_base(&content, path.parent().unwrap_or(Path::new(".")))
            .map_err(|error| match error {
                PegasusError::Parse { location, message } => PegasusError::Parse {
                    location: format!(" in {}{location}", path.display()),
                    message,
                },
                other => other,
            })
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
            surface.set_chordwise_panel_count(paneling[0].max(0.0) as usize);
            surface.set_numeric_chordwise_spacing(paneling[1]);
        }
        if paneling.len() >= 4 {
            surface.set_spanwise_panel_count(paneling[2].max(0.0) as usize);
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
                        surface.set_component_index(value.max(0.0) as usize);
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
                        return Err(PegasusError::parse(
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
                        current.set_spanwise_panel_count(values[5].max(0.0) as usize);
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
                    let (_, relative) = cursor.next_required_after("AFILE", line_number)?;
                    let path = resolve_path(base, relative);
                    if let Ok(airfoil) = Airfoil::from_coordinate_file(&path) {
                        if let Some(current) = section.as_mut() {
                            current.set_airfoil(airfoil);
                        }
                    } else if self.verbose {
                        eprintln!("Warning: could not load airfoil {}", path.display());
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
            body.set_node_count(paneling[0].max(0.0) as usize);
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
        self.next().ok_or_else(|| PegasusError::Parse {
            location: String::new(),
            message: format!("unexpected end of AVL data while reading {context}"),
        })
    }

    fn next_required_after(&mut self, context: &str, line: usize) -> Result<(usize, &'a str)> {
        self.next().ok_or_else(|| {
            PegasusError::parse(line, format!("unexpected end of AVL data after {context}"))
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

fn resolve_path(base: &Path, relative: &str) -> PathBuf {
    let relative = Path::new(relative.trim());
    if relative.is_absolute() {
        relative.to_owned()
    } else {
        base.join(relative)
    }
}

fn parse_control(definition: &str, line: usize) -> Result<ControlSurfaceDefinition> {
    let fields: Vec<_> = definition.split_whitespace().collect();
    if fields.len() < 7 {
        return Err(PegasusError::parse(
            line,
            "CONTROL needs name gain hinge_x hx hy hz sign",
        ));
    }
    let number = |index: usize| {
        fields[index]
            .parse::<Scalar>()
            .map_err(|_| PegasusError::parse(line, "invalid numeric CONTROL field"))
    };
    Ok(ControlSurfaceDefinition {
        name: fields[0].to_owned(),
        gain: number(1)?,
        hinge_x_fraction: number(2)?,
        hinge_vector: Vec3::new(number(3)?, number(4)?, number(5)?),
        sign: number(6)?,
    })
}
