//! Atomic `.vtu` and `.vtp` export with optional `vtkio` support.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use crate::aero::AerodynamicResults;
use crate::core::{Scalar, Vec3};
use crate::error::{PegasusError, Result};
use crate::solver::{FlightCondition, SolutionState, VortexLattice, VortexLatticeSolver};

/// Fields and products included in visualization output.
#[derive(Clone, Copy, Debug)]
pub struct VtkExportOptions {
    pub export_cp: bool,
    pub export_circulation: bool,
    pub export_normals: bool,
    pub export_indices: bool,
    pub export_wake: bool,
    pub wake_length_chords: Scalar,
    pub export_induced_velocity: bool,
    pub export_spanwise_distributions: bool,
    pub binary_format: bool,
}

impl Default for VtkExportOptions {
    fn default() -> Self {
        Self {
            export_cp: true,
            export_circulation: true,
            export_normals: false,
            export_indices: false,
            export_wake: true,
            wake_length_chords: 10.0,
            export_induced_velocity: true,
            export_spanwise_distributions: true,
            binary_format: true,
        }
    }
}

/// Whether this build can serialize VTK products.
#[must_use]
pub const fn is_available() -> bool {
    cfg!(feature = "vtk")
}

/// Namespace for the four visualization products.
pub struct VtkExporter;

impl VtkExporter {
    pub fn export_surface_mesh(
        path: impl AsRef<Path>,
        mesh: &VortexLattice,
        solution: &SolutionState,
        options: VtkExportOptions,
    ) -> Result<()> {
        ensure_available()?;
        if !solution.is_valid || mesh.vortex_count() == 0 {
            return Err(PegasusError::InvalidSolution);
        }
        let mut product = Product::unstructured();
        for (surface_index, info) in mesh.surface_info().iter().enumerate() {
            for strip_offset in 0..info.num_strips {
                let strip_index = info.first_strip_index + strip_offset;
                let strip = &mesh.strips()[strip_index];
                for offset in 0..strip.num_vortices {
                    let index = strip.first_vortex_index + offset;
                    let vortex = &mesh.vortices()[index];
                    let (front_left, front_right) = if offset == 0 {
                        (strip.leading_edge_left, strip.leading_edge_right)
                    } else {
                        let prior = &mesh.vortices()[index - 1];
                        (prior.bound_vortex_left, prior.bound_vortex_right)
                    };
                    product.add_cell(
                        &[
                            front_left,
                            front_right,
                            vortex.bound_vortex_right,
                            vortex.bound_vortex_left,
                        ],
                        CellKind::Quad,
                    );
                    product.cell_value("Cp", 1, vortex.delta_cp);
                    product.cell_value("Circulation", 1, solution.circulation[index]);
                    product.cell_value("Surface", 1, surface_index as Scalar);
                    product.cell_value("Strip", 1, strip_index as Scalar);
                    product.cell_vector("Normal", vortex.normal_at_control_point);
                }
            }
        }
        if !options.export_cp {
            product.remove_cell_array("Cp");
        }
        if !options.export_circulation {
            product.remove_cell_array("Circulation");
        }
        if !options.export_indices {
            product.remove_cell_array("Surface");
            product.remove_cell_array("Strip");
        }
        if !options.export_normals {
            product.remove_cell_array("Normal");
        }
        write_product(path.as_ref(), product, options.binary_format)
    }

    pub fn export_spanwise_loads(
        path: impl AsRef<Path>,
        _mesh: &VortexLattice,
        results: &AerodynamicResults,
        options: VtkExportOptions,
    ) -> Result<()> {
        ensure_available()?;
        if !results.is_valid {
            return Err(PegasusError::InvalidSolution);
        }
        let mut product = Product::poly();
        for surface in &results.surface_loads {
            let mut previous = None;
            for strip in &surface.strip_loads {
                let point = product.add_point(strip.quarter_chord_position);
                if let Some(prior) = previous {
                    product.cells.push(vec![prior, point]);
                }
                previous = Some(point);
                product.point_value("CL", 1, strip.coefficients.lift);
                product.point_value("CD", 1, strip.coefficients.drag);
                product.point_value("CY", 1, strip.coefficients.side);
                product.point_value("SpanLoading", 1, strip.normal_span_loading);
                product.point_value("InducedAlpha", 1, strip.induced_angle);
                product.point_vector("CM", strip.body_moment);
                product.point_vector("InducedVelocity", strip.induced_velocity);
                product.point_value("Surface", 1, surface.surface_index as Scalar);
            }
        }
        if !options.export_induced_velocity {
            product.remove_point_array("InducedVelocity");
        }
        if !options.export_indices {
            product.remove_point_array("Surface");
        }
        write_product(path.as_ref(), product, options.binary_format)
    }

    pub fn export_wake(
        path: impl AsRef<Path>,
        mesh: &VortexLattice,
        solution: &SolutionState,
        condition: &FlightCondition,
        reference_chord: Scalar,
        options: VtkExportOptions,
    ) -> Result<()> {
        ensure_available()?;
        if !solution.is_valid || mesh.strip_count() == 0 {
            return Err(PegasusError::InvalidSolution);
        }
        let mut product = Product::poly();
        let direction = VortexLatticeSolver::freestream_direction(condition);
        let length = options.wake_length_chords * reference_chord;
        for (surface_index, info) in mesh.surface_info().iter().enumerate() {
            for strip_offset in 0..info.num_strips {
                let strip = &mesh.strips()[info.first_strip_index + strip_offset];
                if strip.num_vortices == 0 {
                    continue;
                }
                let vortex = &mesh.vortices()[strip.first_vortex_index + strip.num_vortices - 1];
                let start = product.add_point(vortex.bound_vortex_center);
                let end = product.add_point(vortex.bound_vortex_center + length * direction);
                product.cells.push(vec![start, end]);
                product.cell_value(
                    "Circulation",
                    1,
                    solution
                        .circulation
                        .rows(strip.first_vortex_index, strip.num_vortices)
                        .sum(),
                );
                product.cell_value("Surface", 1, surface_index as Scalar);
            }
        }
        if !options.export_circulation {
            product.remove_cell_array("Circulation");
        }
        if !options.export_indices {
            product.remove_cell_array("Surface");
        }
        write_product(path.as_ref(), product, options.binary_format)
    }

    pub fn export_velocity_field(
        path: impl AsRef<Path>,
        mesh: &VortexLattice,
        solution: &SolutionState,
        condition: &FlightCondition,
        options: VtkExportOptions,
    ) -> Result<()> {
        ensure_available()?;
        if !solution.is_valid || mesh.vortex_count() == 0 {
            return Err(PegasusError::InvalidSolution);
        }
        let mut product = Product::poly();
        let freestream = VortexLatticeSolver::freestream_direction(condition);
        for (index, vortex) in mesh.vortices().iter().enumerate() {
            product.add_point(vortex.control_point);
            let induced = Vec3::new(
                solution.induced_velocity[(index, 0)],
                solution.induced_velocity[(index, 1)],
                solution.induced_velocity[(index, 2)],
            );
            let total = freestream + induced;
            product.point_vector("InducedVelocity", induced);
            product.point_vector("TotalVelocity", total);
            product.point_value("VelocityMagnitude", 1, total.norm());
            product.point_value("Surface", 1, vortex.surface_index as Scalar);
        }
        if !options.export_induced_velocity {
            product.remove_point_array("InducedVelocity");
        }
        if !options.export_indices {
            product.remove_point_array("Surface");
        }
        write_product(path.as_ref(), product, options.binary_format)
    }

    pub fn export_all(
        directory: impl AsRef<Path>,
        mesh: &VortexLattice,
        solution: &SolutionState,
        results: &AerodynamicResults,
        condition: &FlightCondition,
        reference_chord: Scalar,
        options: VtkExportOptions,
    ) -> Result<()> {
        ensure_available()?;
        let directory = directory.as_ref();
        std::fs::create_dir_all(directory).map_err(|error| PegasusError::io(directory, error))?;
        Self::export_surface_mesh(directory.join("surface_mesh.vtu"), mesh, solution, options)?;
        if options.export_spanwise_distributions {
            Self::export_spanwise_loads(
                directory.join("spanwise_loads.vtp"),
                mesh,
                results,
                options,
            )?;
        }
        if options.export_wake {
            Self::export_wake(
                directory.join("wake.vtp"),
                mesh,
                solution,
                condition,
                reference_chord,
                options,
            )?;
        }
        if options.export_induced_velocity {
            Self::export_velocity_field(
                directory.join("velocity_field.vtp"),
                mesh,
                solution,
                condition,
                options,
            )?;
        }
        Ok(())
    }
}

fn ensure_available() -> Result<()> {
    if is_available() {
        Ok(())
    } else {
        Err(PegasusError::FeatureDisabled("vtk"))
    }
}

#[derive(Clone, Copy)]
enum CellKind {
    Quad,
}

struct Array {
    name: &'static str,
    components: usize,
    data: Vec<Scalar>,
}

enum ProductKind {
    Unstructured,
    Poly,
}

struct Product {
    kind: ProductKind,
    points: Vec<Scalar>,
    cells: Vec<Vec<u64>>,
    cell_kinds: Vec<CellKind>,
    point_data: Vec<Array>,
    cell_data: Vec<Array>,
}

impl Product {
    fn unstructured() -> Self {
        Self::new(ProductKind::Unstructured)
    }

    fn poly() -> Self {
        Self::new(ProductKind::Poly)
    }

    fn new(kind: ProductKind) -> Self {
        Self {
            kind,
            points: Vec::new(),
            cells: Vec::new(),
            cell_kinds: Vec::new(),
            point_data: Vec::new(),
            cell_data: Vec::new(),
        }
    }

    fn add_point(&mut self, point: Vec3) -> u64 {
        let index = (self.points.len() / 3) as u64;
        self.points.extend([point.x, point.y, point.z]);
        index
    }

    fn add_cell(&mut self, points: &[Vec3], kind: CellKind) {
        let cell = points.iter().map(|point| self.add_point(*point)).collect();
        self.cells.push(cell);
        self.cell_kinds.push(kind);
    }

    fn point_value(&mut self, name: &'static str, components: usize, value: Scalar) {
        array_mut(&mut self.point_data, name, components)
            .data
            .push(value);
    }

    fn point_vector(&mut self, name: &'static str, value: Vec3) {
        array_mut(&mut self.point_data, name, 3)
            .data
            .extend([value.x, value.y, value.z]);
    }

    fn cell_value(&mut self, name: &'static str, components: usize, value: Scalar) {
        array_mut(&mut self.cell_data, name, components)
            .data
            .push(value);
    }

    fn cell_vector(&mut self, name: &'static str, value: Vec3) {
        array_mut(&mut self.cell_data, name, 3)
            .data
            .extend([value.x, value.y, value.z]);
    }

    fn remove_cell_array(&mut self, name: &str) {
        self.cell_data.retain(|array| array.name != name);
    }

    fn remove_point_array(&mut self, name: &str) {
        self.point_data.retain(|array| array.name != name);
    }
}

fn array_mut<'a>(
    arrays: &'a mut Vec<Array>,
    name: &'static str,
    components: usize,
) -> &'a mut Array {
    if let Some(index) = arrays.iter().position(|array| array.name == name) {
        return &mut arrays[index];
    }
    arrays.push(Array {
        name,
        components,
        data: Vec::new(),
    });
    arrays.last_mut().expect("array was just inserted")
}

fn temporary_sibling(path: &Path) -> Result<PathBuf> {
    let parent = path.parent().unwrap_or_else(|| Path::new("."));
    let stem = path
        .file_stem()
        .and_then(|value| value.to_str())
        .ok_or_else(|| {
            PegasusError::InvalidGeometry("VTK output path needs a UTF-8 file name".into())
        })?;
    let extension = path
        .extension()
        .and_then(|value| value.to_str())
        .ok_or_else(|| {
            PegasusError::InvalidGeometry("VTK output path needs .vtu or .vtp extension".into())
        })?;
    Ok(parent.join(format!(".{stem}.{}.tmp.{extension}", std::process::id())))
}

fn write_product(path: &Path, product: Product, binary: bool) -> Result<()> {
    let temporary = temporary_sibling(path)?;
    let result = if binary {
        write_binary(&temporary, product)
    } else {
        let xml = ascii_xml(&product);
        std::fs::write(&temporary, xml).map_err(|error| PegasusError::io(&temporary, error))
    };
    if let Err(error) = result {
        let _ = std::fs::remove_file(&temporary);
        return Err(error);
    }
    std::fs::rename(&temporary, path).map_err(|error| PegasusError::io(path, error))
}

#[cfg(feature = "vtk")]
fn write_binary(path: &Path, product: Product) -> Result<()> {
    use vtkio::model::{
        Attribute, Attributes, ByteOrder, CellType, Cells, DataSet, IOBuffer, PolyDataPiece,
        UnstructuredGridPiece, Version, VertexNumbers, Vtk,
    };

    let attributes = |arrays: Vec<Array>| {
        arrays
            .into_iter()
            .map(|array| {
                Attribute::scalars(array.name, array.components as u32).with_data(array.data)
            })
            .collect()
    };
    let data = Attributes {
        point: attributes(product.point_data),
        cell: attributes(product.cell_data),
    };
    let (connectivity, offsets) = connectivity(&product.cells);
    let data = match product.kind {
        ProductKind::Unstructured => DataSet::inline(UnstructuredGridPiece {
            points: IOBuffer::F64(product.points),
            cells: Cells {
                cell_verts: VertexNumbers::XML {
                    connectivity,
                    offsets,
                },
                types: product
                    .cell_kinds
                    .iter()
                    .map(|kind| match kind {
                        CellKind::Quad => CellType::Quad,
                    })
                    .collect(),
            },
            data,
        }),
        ProductKind::Poly => DataSet::inline(PolyDataPiece {
            points: IOBuffer::F64(product.points),
            lines: (!product.cells.is_empty()).then_some(VertexNumbers::XML {
                connectivity,
                offsets,
            }),
            data,
            ..PolyDataPiece::default()
        }),
    };
    Vtk {
        version: Version::new((1, 0)),
        title: "Pegasus VLM".into(),
        byte_order: ByteOrder::LittleEndian,
        file_path: None,
        data,
    }
    .export(path)
    .map_err(|error| PegasusError::io(path, std::io::Error::other(error)))
}

#[cfg(not(feature = "vtk"))]
fn write_binary(_path: &Path, _product: Product) -> Result<()> {
    Err(PegasusError::FeatureDisabled("vtk"))
}

fn connectivity(cells: &[Vec<u64>]) -> (Vec<u64>, Vec<u64>) {
    let mut connectivity = Vec::new();
    let mut offsets = Vec::with_capacity(cells.len());
    for cell in cells {
        connectivity.extend(cell);
        offsets.push(connectivity.len() as u64);
    }
    (connectivity, offsets)
}

fn ascii_xml(product: &Product) -> String {
    let mut xml = String::new();
    let vtk_type = match product.kind {
        ProductKind::Unstructured => "UnstructuredGrid",
        ProductKind::Poly => "PolyData",
    };
    let _ = write!(
        xml,
        "<?xml version=\"1.0\"?>\n<VTKFile type=\"{vtk_type}\" version=\"1.0\" \
         byte_order=\"LittleEndian\" header_type=\"UInt64\">\n  <{vtk_type}>\n"
    );
    match product.kind {
        ProductKind::Unstructured => {
            let _ = writeln!(
                xml,
                "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">",
                product.points.len() / 3,
                product.cells.len()
            );
        }
        ProductKind::Poly => {
            let _ = writeln!(
                xml,
                "    <Piece NumberOfPoints=\"{}\" NumberOfVerts=\"0\" \
                 NumberOfLines=\"{}\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">",
                product.points.len() / 3,
                product.cells.len()
            );
        }
    }
    write_arrays(&mut xml, "PointData", &product.point_data);
    write_arrays(&mut xml, "CellData", &product.cell_data);
    let _ = writeln!(xml, "      <Points>");
    write_data_array(&mut xml, "", 3, &product.points, 8);
    let _ = writeln!(xml, "      </Points>");
    let (connectivity, offsets) = connectivity(&product.cells);
    match product.kind {
        ProductKind::Unstructured => {
            let _ = writeln!(xml, "      <Cells>");
            write_integer_array(&mut xml, "connectivity", &connectivity, "Int64");
            write_integer_array(&mut xml, "offsets", &offsets, "Int64");
            let types: Vec<_> = product.cell_kinds.iter().map(|_| 9_u64).collect();
            write_integer_array(&mut xml, "types", &types, "UInt8");
            let _ = writeln!(xml, "      </Cells>");
        }
        ProductKind::Poly if !product.cells.is_empty() => {
            let _ = writeln!(xml, "      <Lines>");
            write_integer_array(&mut xml, "connectivity", &connectivity, "Int64");
            write_integer_array(&mut xml, "offsets", &offsets, "Int64");
            let _ = writeln!(xml, "      </Lines>");
        }
        ProductKind::Poly => {}
    }
    let _ = writeln!(xml, "    </Piece>\n  </{vtk_type}>\n</VTKFile>");
    xml
}

fn write_arrays(xml: &mut String, tag: &str, arrays: &[Array]) {
    if arrays.is_empty() {
        let _ = writeln!(xml, "      <{tag}/>");
        return;
    }
    let _ = writeln!(xml, "      <{tag}>");
    for array in arrays {
        write_data_array(xml, array.name, array.components, &array.data, 8);
    }
    let _ = writeln!(xml, "      </{tag}>");
}

fn write_data_array(
    xml: &mut String,
    name: &str,
    components: usize,
    data: &[Scalar],
    indent: usize,
) {
    let spaces = " ".repeat(indent);
    let name = if name.is_empty() {
        String::new()
    } else {
        format!(" Name=\"{name}\"")
    };
    let _ = write!(
        xml,
        "{spaces}<DataArray type=\"Float64\"{name} \
         NumberOfComponents=\"{components}\" format=\"ascii\">"
    );
    for value in data {
        let _ = write!(xml, " {value:.17e}");
    }
    let _ = writeln!(xml, " </DataArray>");
}

fn write_integer_array(xml: &mut String, name: &str, data: &[u64], scalar_type: &str) {
    let _ = write!(
        xml,
        "        <DataArray type=\"{scalar_type}\" Name=\"{name}\" \
         NumberOfComponents=\"1\" format=\"ascii\">"
    );
    for value in data {
        let _ = write!(xml, " {value}");
    }
    let _ = writeln!(xml, " </DataArray>");
}
