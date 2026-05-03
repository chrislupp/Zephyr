/*
    Project Pegasus
    Version 0.0.0

    Copyright 2024, Christopher A. Lupp

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "pegasus/visualization/vtk_export.h"

#include "pegasus/aero/aerodynamic_results.h"
#include "pegasus/solver/horseshoe_vortex.h"
#include "pegasus/solver/vortex_lattice.h"
#include "pegasus/solver/vortex_lattice_solver.h"

#include <filesystem>
#include <map>
#include <tuple>

#ifdef PEGASUS_USE_VTK
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkLine.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkQuad.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#endif

namespace pegasus {

bool VtkExporter::IsAvailable() {
#ifdef PEGASUS_USE_VTK
    return true;
#else
    return false;
#endif
}

bool VtkExporter::ExportSurfaceMesh(const std::string& filename,
                                    const VortexLattice& mesh,
                                    const SolutionState& solution,
                                    const VtkExportOptions& options) {
#ifdef PEGASUS_USE_VTK
    if (!solution.is_valid || mesh.GetNumVortices() == 0) {
        return false;
    }

    vtkNew<vtkUnstructuredGrid> grid;
    vtkNew<vtkPoints> points;
    vtkNew<vtkCellArray> cells;

    // Data arrays
    vtkNew<vtkDoubleArray> cp_array;
    cp_array->SetName("Cp");
    cp_array->SetNumberOfComponents(1);

    vtkNew<vtkDoubleArray> gamma_array;
    gamma_array->SetName("Circulation");
    gamma_array->SetNumberOfComponents(1);

    vtkNew<vtkDoubleArray> surface_array;
    surface_array->SetName("Surface");
    surface_array->SetNumberOfComponents(1);

    vtkNew<vtkDoubleArray> strip_array;
    strip_array->SetName("Strip");
    strip_array->SetNumberOfComponents(1);

    vtkNew<vtkDoubleArray> normal_array;
    normal_array->SetName("Normal");
    normal_array->SetNumberOfComponents(3);

    // Build point list with unique vertices
    // Map from (vortex_index, corner) -> point_id
    // corner: 0=bound_left, 1=bound_right
    std::map<std::pair<Index, int>, vtkIdType> point_map;

    auto get_or_create_point = [&](Index vor_idx, int corner) -> vtkIdType {
        auto key = std::make_pair(vor_idx, corner);
        auto it = point_map.find(key);
        if (it != point_map.end()) {
            return it->second;
        }

        const HorseshoeVortex& vor = mesh.GetVortex(vor_idx);
        Vec3 pos =
            (corner == 0) ? vor.bound_vortex_left : vor.bound_vortex_right;
        vtkIdType pt_id = points->InsertNextPoint(pos.x(), pos.y(), pos.z());
        point_map[key] = pt_id;
        return pt_id;
    };

    // For each surface
    for (Index s = 0; s < mesh.GetNumSurfaces(); ++s) {
        const SurfaceMeshInfo& info = mesh.GetSurfaceInfo(s);

        // For each strip in surface
        for (Index j = 0; j < info.num_strips; ++j) {
            Index strip_idx = info.first_strip_index + j;
            const Strip& strip = mesh.GetStrip(strip_idx);

            // For each vortex (panel) in strip
            for (Index k = 0; k < strip.num_vortices; ++k) {
                Index vor_idx = strip.first_vortex_index + k;
                const HorseshoeVortex& vor = mesh.GetVortex(vor_idx);

                // Get the four corners of the panel
                // Front edge: from previous vortex (or leading edge for first)
                // Back edge: current vortex bound
                vtkIdType pt0, pt1, pt2, pt3;

                if (k == 0) {
                    // First panel - front edge is leading edge
                    // Use strip LE points
                    pt0 = points->InsertNextPoint(strip.leading_edge_left.x(),
                                                  strip.leading_edge_left.y(),
                                                  strip.leading_edge_left.z());
                    pt1 = points->InsertNextPoint(strip.leading_edge_right.x(),
                                                  strip.leading_edge_right.y(),
                                                  strip.leading_edge_right.z());
                } else {
                    // Use previous vortex's bound vortex as front edge
                    Index prev_vor_idx = strip.first_vortex_index + k - 1;
                    pt0 = get_or_create_point(prev_vor_idx, 0);
                    pt1 = get_or_create_point(prev_vor_idx, 1);
                }

                // Back edge is current vortex's bound vortex
                pt2 = get_or_create_point(vor_idx, 1);
                pt3 = get_or_create_point(vor_idx, 0);

                // Create quad cell (order: counter-clockwise from pt0)
                vtkNew<vtkQuad> quad;
                quad->GetPointIds()->SetId(0, pt0);
                quad->GetPointIds()->SetId(1, pt1);
                quad->GetPointIds()->SetId(2, pt2);
                quad->GetPointIds()->SetId(3, pt3);
                cells->InsertNextCell(quad);

                // Add cell data
                if (options.export_cp) {
                    cp_array->InsertNextValue(vor.delta_cp);
                }
                if (options.export_circulation) {
                    gamma_array->InsertNextValue(solution.circulation(vor_idx));
                }
                if (options.export_indices) {
                    surface_array->InsertNextValue(static_cast<double>(s));
                    strip_array->InsertNextValue(
                        static_cast<double>(strip_idx));
                }
                if (options.export_normals) {
                    normal_array->InsertNextTuple3(
                        vor.normal_at_control_point.x(),
                        vor.normal_at_control_point.y(),
                        vor.normal_at_control_point.z());
                }
            }
        }
    }

    grid->SetPoints(points);
    grid->SetCells(VTK_QUAD, cells);

    if (options.export_cp) {
        grid->GetCellData()->AddArray(cp_array);
    }
    if (options.export_circulation) {
        grid->GetCellData()->AddArray(gamma_array);
    }
    if (options.export_indices) {
        grid->GetCellData()->AddArray(surface_array);
        grid->GetCellData()->AddArray(strip_array);
    }
    if (options.export_normals) {
        grid->GetCellData()->AddArray(normal_array);
    }

    vtkNew<vtkXMLUnstructuredGridWriter> writer;
    writer->SetFileName(filename.c_str());
    writer->SetInputData(grid);
    writer->SetDataMode(options.binary_format ? vtkXMLWriter::Binary
                                              : vtkXMLWriter::Ascii);
    return writer->Write() == 1;
#else
    (void)filename;
    (void)mesh;
    (void)solution;
    (void)options;
    return false;
#endif
}

bool VtkExporter::ExportSpanwiseDistributions(const std::string& filename,
                                              const VortexLattice& mesh,
                                              const AerodynamicResults& results,
                                              const VtkExportOptions& options) {
#ifdef PEGASUS_USE_VTK
    if (!results.is_valid) {
        return false;
    }

    vtkNew<vtkPolyData> polydata;
    vtkNew<vtkPoints> points;
    vtkNew<vtkCellArray> lines;

    // Scalar arrays
    vtkNew<vtkDoubleArray> cl_array;
    cl_array->SetName("CL");
    cl_array->SetNumberOfComponents(1);

    vtkNew<vtkDoubleArray> cd_array;
    cd_array->SetName("CD");
    cd_array->SetNumberOfComponents(1);

    vtkNew<vtkDoubleArray> cy_array;
    cy_array->SetName("CY");
    cy_array->SetNumberOfComponents(1);

    vtkNew<vtkDoubleArray> cnc_array;
    cnc_array->SetName("SpanLoading");
    cnc_array->SetNumberOfComponents(1);

    vtkNew<vtkDoubleArray> alpha_i_array;
    alpha_i_array->SetName("InducedAlpha");
    alpha_i_array->SetNumberOfComponents(1);

    // Vector arrays
    vtkNew<vtkDoubleArray> cm_array;
    cm_array->SetName("CM");
    cm_array->SetNumberOfComponents(3);

    vtkNew<vtkDoubleArray> velocity_array;
    velocity_array->SetName("InducedVelocity");
    velocity_array->SetNumberOfComponents(3);

    // Surface index array
    vtkNew<vtkDoubleArray> surface_array;
    surface_array->SetName("Surface");
    surface_array->SetNumberOfComponents(1);

    // Process each surface
    for (const auto& surf_loads : results.surface_loads) {
        vtkIdType prev_point = -1;
        vtkIdType first_point_in_surface = -1;

        for (const auto& strip_load : surf_loads.strip_loads) {
            // Add point at quarter-chord
            vtkIdType pt =
                points->InsertNextPoint(strip_load.quarter_chord_position.x(),
                                        strip_load.quarter_chord_position.y(),
                                        strip_load.quarter_chord_position.z());

            if (first_point_in_surface < 0) {
                first_point_in_surface = pt;
            }

            // Add scalar data
            cl_array->InsertNextValue(strip_load.cl);
            cd_array->InsertNextValue(strip_load.cd);
            cy_array->InsertNextValue(strip_load.cy);
            cnc_array->InsertNextValue(strip_load.cnc);
            alpha_i_array->InsertNextValue(strip_load.alpha_induced);

            // Add vector data
            cm_array->InsertNextTuple3(strip_load.cm.x(), strip_load.cm.y(),
                                       strip_load.cm.z());
            velocity_array->InsertNextTuple3(strip_load.induced_velocity.x(),
                                             strip_load.induced_velocity.y(),
                                             strip_load.induced_velocity.z());

            surface_array->InsertNextValue(
                static_cast<double>(surf_loads.surface_index));

            // Connect to previous point within same surface
            if (prev_point >= 0 && prev_point >= first_point_in_surface) {
                vtkNew<vtkLine> line;
                line->GetPointIds()->SetId(0, prev_point);
                line->GetPointIds()->SetId(1, pt);
                lines->InsertNextCell(line);
            }
            prev_point = pt;
        }
    }

    polydata->SetPoints(points);
    polydata->SetLines(lines);

    polydata->GetPointData()->AddArray(cl_array);
    polydata->GetPointData()->AddArray(cd_array);
    polydata->GetPointData()->AddArray(cy_array);
    polydata->GetPointData()->AddArray(cnc_array);
    polydata->GetPointData()->AddArray(alpha_i_array);
    polydata->GetPointData()->AddArray(cm_array);
    polydata->GetPointData()->AddArray(velocity_array);
    polydata->GetPointData()->AddArray(surface_array);

    vtkNew<vtkXMLPolyDataWriter> writer;
    writer->SetFileName(filename.c_str());
    writer->SetInputData(polydata);
    writer->SetDataMode(options.binary_format ? vtkXMLWriter::Binary
                                              : vtkXMLWriter::Ascii);
    return writer->Write() == 1;
#else
    (void)filename;
    (void)mesh;
    (void)results;
    (void)options;
    return false;
#endif
}

bool VtkExporter::ExportWakeGeometry(const std::string& filename,
                                     const VortexLattice& mesh,
                                     const SolutionState& solution,
                                     const FlightCondition& condition,
                                     Scalar reference_chord,
                                     const VtkExportOptions& options) {
#ifdef PEGASUS_USE_VTK
    if (!solution.is_valid || mesh.GetNumStrips() == 0) {
        return false;
    }

    vtkNew<vtkPolyData> polydata;
    vtkNew<vtkPoints> points;
    vtkNew<vtkCellArray> lines;

    // Cell data for circulation
    vtkNew<vtkDoubleArray> gamma_array;
    gamma_array->SetName("Circulation");
    gamma_array->SetNumberOfComponents(1);

    vtkNew<vtkDoubleArray> surface_array;
    surface_array->SetName("Surface");
    surface_array->SetNumberOfComponents(1);

    // Compute wake direction (downstream from freestream)
    Vec3 freestream =
        VortexLatticeSolver::ComputeFreestreamDirection(condition);
    Vec3 wake_dir = freestream;  // Wake extends in freestream direction
    Scalar wake_length = options.wake_length_chords * reference_chord;

    // For each strip, create wake filament from trailing edge
    for (Index s = 0; s < mesh.GetNumSurfaces(); ++s) {
        const SurfaceMeshInfo& info = mesh.GetSurfaceInfo(s);

        for (Index j = 0; j < info.num_strips; ++j) {
            Index strip_idx = info.first_strip_index + j;
            const Strip& strip = mesh.GetStrip(strip_idx);

            // Get last vortex in strip (trailing edge)
            if (strip.num_vortices == 0)
                continue;

            Index last_vor_idx =
                strip.first_vortex_index + strip.num_vortices - 1;
            const HorseshoeVortex& last_vor = mesh.GetVortex(last_vor_idx);

            // Sum circulation for the strip
            Scalar strip_circulation = 0.0;
            for (Index k = 0; k < strip.num_vortices; ++k) {
                strip_circulation +=
                    solution.circulation(strip.first_vortex_index + k);
            }

            // Create wake filament from center of trailing edge
            Vec3 te_center = last_vor.bound_vortex_center;
            Vec3 wake_end = te_center + wake_length * wake_dir;

            vtkIdType pt0 = points->InsertNextPoint(
                te_center.x(), te_center.y(), te_center.z());
            vtkIdType pt1 = points->InsertNextPoint(wake_end.x(), wake_end.y(),
                                                    wake_end.z());

            vtkNew<vtkLine> line;
            line->GetPointIds()->SetId(0, pt0);
            line->GetPointIds()->SetId(1, pt1);
            lines->InsertNextCell(line);

            gamma_array->InsertNextValue(strip_circulation);
            surface_array->InsertNextValue(static_cast<double>(s));
        }
    }

    polydata->SetPoints(points);
    polydata->SetLines(lines);
    polydata->GetCellData()->AddArray(gamma_array);
    polydata->GetCellData()->AddArray(surface_array);

    vtkNew<vtkXMLPolyDataWriter> writer;
    writer->SetFileName(filename.c_str());
    writer->SetInputData(polydata);
    writer->SetDataMode(options.binary_format ? vtkXMLWriter::Binary
                                              : vtkXMLWriter::Ascii);
    return writer->Write() == 1;
#else
    (void)filename;
    (void)mesh;
    (void)solution;
    (void)condition;
    (void)reference_chord;
    (void)options;
    return false;
#endif
}

bool VtkExporter::ExportVelocityField(const std::string& filename,
                                      const VortexLattice& mesh,
                                      const SolutionState& solution,
                                      const FlightCondition& condition,
                                      const VtkExportOptions& options) {
#ifdef PEGASUS_USE_VTK
    if (!solution.is_valid || mesh.GetNumVortices() == 0) {
        return false;
    }

    vtkNew<vtkPolyData> polydata;
    vtkNew<vtkPoints> points;

    // Vector arrays
    vtkNew<vtkDoubleArray> induced_velocity_array;
    induced_velocity_array->SetName("InducedVelocity");
    induced_velocity_array->SetNumberOfComponents(3);

    vtkNew<vtkDoubleArray> total_velocity_array;
    total_velocity_array->SetName("TotalVelocity");
    total_velocity_array->SetNumberOfComponents(3);

    // Scalar arrays
    vtkNew<vtkDoubleArray> velocity_magnitude_array;
    velocity_magnitude_array->SetName("VelocityMagnitude");
    velocity_magnitude_array->SetNumberOfComponents(1);

    vtkNew<vtkDoubleArray> surface_array;
    surface_array->SetName("Surface");
    surface_array->SetNumberOfComponents(1);

    // Compute freestream velocity
    Vec3 freestream =
        VortexLatticeSolver::ComputeFreestreamDirection(condition);

    // Add control points with velocity data
    for (Index i = 0; i < mesh.GetNumVortices(); ++i) {
        const HorseshoeVortex& vor = mesh.GetVortex(i);

        points->InsertNextPoint(vor.control_point.x(), vor.control_point.y(),
                                vor.control_point.z());

        // Get induced velocity from solution
        Vec3 v_induced = Vec3::Zero();
        if (solution.induced_velocity.rows() > i) {
            v_induced = solution.induced_velocity.row(i);
        }

        Vec3 v_total = freestream + v_induced;

        induced_velocity_array->InsertNextTuple3(v_induced.x(), v_induced.y(),
                                                 v_induced.z());
        total_velocity_array->InsertNextTuple3(v_total.x(), v_total.y(),
                                               v_total.z());
        velocity_magnitude_array->InsertNextValue(v_total.norm());
        surface_array->InsertNextValue(static_cast<double>(vor.surface_index));
    }

    polydata->SetPoints(points);
    polydata->GetPointData()->AddArray(induced_velocity_array);
    polydata->GetPointData()->AddArray(total_velocity_array);
    polydata->GetPointData()->AddArray(velocity_magnitude_array);
    polydata->GetPointData()->AddArray(surface_array);

    vtkNew<vtkXMLPolyDataWriter> writer;
    writer->SetFileName(filename.c_str());
    writer->SetInputData(polydata);
    writer->SetDataMode(options.binary_format ? vtkXMLWriter::Binary
                                              : vtkXMLWriter::Ascii);
    return writer->Write() == 1;
#else
    (void)filename;
    (void)mesh;
    (void)solution;
    (void)condition;
    (void)options;
    return false;
#endif
}

bool VtkExporter::ExportAll(const std::string& directory,
                            const VortexLattice& mesh,
                            const SolutionState& solution,
                            const AerodynamicResults& results,
                            const FlightCondition& condition,
                            Scalar reference_chord,
                            const VtkExportOptions& options) {
#ifdef PEGASUS_USE_VTK
    // Create directory if needed
    std::filesystem::create_directories(directory);

    std::string sep = directory.back() == '/' ? "" : "/";

    bool success = true;

    // Export surface mesh
    success &= ExportSurfaceMesh(directory + sep + "surface_mesh.vtu", mesh,
                                 solution, options);

    // Export spanwise distributions
    if (options.export_spanwise_distributions) {
        success &= ExportSpanwiseDistributions(
            directory + sep + "spanwise_loads.vtp", mesh, results, options);
    }

    // Export wake
    if (options.export_wake) {
        success &=
            ExportWakeGeometry(directory + sep + "wake.vtp", mesh, solution,
                               condition, reference_chord, options);
    }

    // Export velocity field
    if (options.export_induced_velocity) {
        success &= ExportVelocityField(directory + sep + "velocity_field.vtp",
                                       mesh, solution, condition, options);
    }

    return success;
#else
    (void)directory;
    (void)mesh;
    (void)solution;
    (void)results;
    (void)condition;
    (void)reference_chord;
    (void)options;
    return false;
#endif
}

}  // namespace pegasus
