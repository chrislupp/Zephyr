/*
    Project Pegasus
    Version 0.2.0

    Copyright 2024-2026, Christopher A. Lupp

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

#ifndef PEGASUS_VISUALIZATION_VTK_EXPORT_H_
#define PEGASUS_VISUALIZATION_VTK_EXPORT_H_

#include "pegasus/core/types.h"

#include <string>

namespace pegasus {

// Forward declarations
class VortexLattice;
struct SolutionState;
struct AerodynamicResults;
struct FlightCondition;

/// @brief Configuration options for VTK export.
struct VtkExportOptions {
    /// Include pressure coefficient data on panels
    bool export_cp{true};

    /// Include circulation data on panels
    bool export_circulation{true};

    /// Include surface normal vectors
    bool export_normals{false};

    /// Include strip/surface indices as cell data
    bool export_indices{false};

    /// Export wake geometry (trailing vortex filaments)
    bool export_wake{true};

    /// Wake length as multiple of reference chord
    Scalar wake_length_chords{10.0};

    /// Export induced velocity vectors at control points
    bool export_induced_velocity{true};

    /// Export spanwise distributions as separate file
    bool export_spanwise_distributions{true};

    /// Use binary format (more compact) vs ASCII (human-readable)
    bool binary_format{true};
};

/// @brief Exports VLM results to VTK file format for visualization.
///
/// Provides static methods to export various VLM data structures
/// to VTK XML formats:
/// - `.vtu` (UnstructuredGrid) for panel mesh with Cp
/// - `.vtp` (PolyData) for spanwise loads, wake, velocity field
///
/// All methods gracefully no-op when VTK support is not compiled in.
/// Use IsAvailable() to check if VTK support is present.
///
/// Example usage:
/// @code
/// if (VtkExporter::IsAvailable()) {
///     VtkExportOptions options;
///     options.wake_length_chords = 5.0;
///     VtkExporter::ExportAll("./results/", mesh, solution,
///                            results, condition, options);
/// }
/// @endcode
class VtkExporter {
public:
    /// @brief Check if VTK support is available.
    /// @return true if library was compiled with VTK support
    static bool IsAvailable();

    /// @brief Export surface mesh with Cp distribution.
    ///
    /// Exports the panel mesh as VTK unstructured grid (.vtu) with:
    /// - Panel geometry as quadrilateral cells
    /// - Delta Cp as cell data
    /// - Circulation as cell data
    /// - Surface/strip indices as cell data (optional)
    ///
    /// @param filename Output file path (should end in .vtu)
    /// @param mesh The vortex lattice mesh
    /// @param solution Solution state with circulation values
    /// @param options Export options
    /// @return true if export succeeded, false otherwise
    static bool ExportSurfaceMesh(const std::string& filename,
                                  const VortexLattice& mesh,
                                  const SolutionState& solution,
                                  const VtkExportOptions& options = {});

    /// @brief Export spanwise load distributions.
    ///
    /// Exports spanwise distributions as VTK polydata (.vtp) with:
    /// - Line segments along span (quarter-chord positions)
    /// - CL, CD, CY as point data
    /// - Moment coefficients as point data
    /// - Induced velocity as vector point data
    ///
    /// @param filename Output file path (should end in .vtp)
    /// @param mesh The vortex lattice mesh
    /// @param results Aerodynamic results with strip loads
    /// @param options Export options
    /// @return true if export succeeded, false otherwise
    static bool ExportSpanwiseDistributions(
        const std::string& filename, const VortexLattice& mesh,
        const AerodynamicResults& results,
        const VtkExportOptions& options = {});

    /// @brief Export wake geometry.
    ///
    /// Exports trailing vortex filaments as VTK polydata (.vtp) with:
    /// - Line segments representing trailing vortices
    /// - Circulation strength as cell data
    ///
    /// @param filename Output file path (should end in .vtp)
    /// @param mesh The vortex lattice mesh
    /// @param solution Solution state with circulation values
    /// @param condition Flight condition (for wake direction)
    /// @param reference_chord Reference chord for wake length scaling
    /// @param options Export options (wake_length_chords used)
    /// @return true if export succeeded, false otherwise
    static bool ExportWakeGeometry(const std::string& filename,
                                   const VortexLattice& mesh,
                                   const SolutionState& solution,
                                   const FlightCondition& condition,
                                   Scalar reference_chord,
                                   const VtkExportOptions& options = {});

    /// @brief Export velocity field at control points.
    ///
    /// Exports induced velocity vectors as VTK polydata (.vtp) with:
    /// - Points at control point locations
    /// - Induced velocity as vector point data
    /// - Total velocity (freestream + induced) as vector point data
    /// - Velocity magnitude as scalar point data
    ///
    /// @param filename Output file path (should end in .vtp)
    /// @param mesh The vortex lattice mesh
    /// @param solution Solution state with induced velocities
    /// @param condition Flight condition (for freestream velocity)
    /// @param options Export options
    /// @return true if export succeeded, false otherwise
    static bool ExportVelocityField(const std::string& filename,
                                    const VortexLattice& mesh,
                                    const SolutionState& solution,
                                    const FlightCondition& condition,
                                    const VtkExportOptions& options = {});

    /// @brief Export all visualization data to a directory.
    ///
    /// Convenience method that calls all export methods and creates:
    /// - surface_mesh.vtu - Panel mesh with Cp
    /// - spanwise_loads.vtp - Spanwise distributions
    /// - wake.vtp - Wake geometry
    /// - velocity_field.vtp - Induced velocities
    ///
    /// @param directory Output directory path (created if needed)
    /// @param mesh The vortex lattice mesh
    /// @param solution Solution state
    /// @param results Aerodynamic results
    /// @param condition Flight condition
    /// @param reference_chord Reference chord for wake length
    /// @param options Export options
    /// @return true if all exports succeeded, false otherwise
    static bool ExportAll(const std::string& directory,
                          const VortexLattice& mesh,
                          const SolutionState& solution,
                          const AerodynamicResults& results,
                          const FlightCondition& condition,
                          Scalar reference_chord,
                          const VtkExportOptions& options = {});
};

}  // namespace pegasus

#endif  // PEGASUS_VISUALIZATION_VTK_EXPORT_H_
