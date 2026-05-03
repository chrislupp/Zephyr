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

#include "pegasus/aero/force_calculator.h"
#include "pegasus/geometry/aircraft.h"
#include "pegasus/geometry/section.h"
#include "pegasus/geometry/surface.h"
#include "pegasus/solver/vortex_lattice_solver.h"

#include <filesystem>
#include <fstream>

#include <gtest/gtest.h>

namespace pegasus {
namespace {

class VtkExportTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create simple rectangular wing
        Surface wing("Wing");
        wing.SetChordwisePanelCount(4);
        wing.SetSpanwisePanelCount(6);
        wing.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
        wing.AddSection(Section(Vec3(0, 2.0, 0), 1.0, 0.0));

        aircraft_.AddSurface(wing);
        aircraft_.SetReferenceArea(2.0);
        aircraft_.SetReferenceChord(1.0);
        aircraft_.SetReferenceSpan(4.0);

        solver_.SetAircraft(aircraft_);

        FlightCondition cond;
        cond.alpha = 5.0 * kDegreesToRadians;
        solution_ = solver_.Solve(cond);

        ForceCalculator calc;
        results_ =
            calc.ComputeForces(aircraft_, solver_.GetMesh(), solution_, cond);
        condition_ = cond;
    }

    void TearDown() override {
        // Clean up test files
        std::filesystem::remove("test_surface.vtu");
        std::filesystem::remove("test_loads.vtp");
        std::filesystem::remove("test_wake.vtp");
        std::filesystem::remove("test_velocity.vtp");
        std::filesystem::remove_all("test_vtk_output");
    }

    Aircraft aircraft_{"Test Aircraft"};
    VortexLatticeSolver solver_;
    SolutionState solution_;
    AerodynamicResults results_;
    FlightCondition condition_;
};

TEST_F(VtkExportTest, IsAvailableReturnsConsistentValue) {
    // This test just ensures IsAvailable() doesn't crash and returns a bool
    bool available = VtkExporter::IsAvailable();
    // The value depends on whether VTK was found during build
    (void)available;
    SUCCEED();
}

TEST_F(VtkExportTest, ExportSurfaceMeshCreatesFile) {
    if (!VtkExporter::IsAvailable()) {
        GTEST_SKIP() << "VTK not available";
    }

    VtkExportOptions options;
    bool result = VtkExporter::ExportSurfaceMesh(
        "test_surface.vtu", solver_.GetMesh(), solution_, options);

    EXPECT_TRUE(result);
    EXPECT_TRUE(std::filesystem::exists("test_surface.vtu"));

    // Verify file is non-empty
    std::ifstream file("test_surface.vtu");
    EXPECT_TRUE(file.good());
    file.seekg(0, std::ios::end);
    EXPECT_GT(file.tellg(), 0);
}

TEST_F(VtkExportTest, ExportSpanwiseDistributionsCreatesFile) {
    if (!VtkExporter::IsAvailable()) {
        GTEST_SKIP() << "VTK not available";
    }

    VtkExportOptions options;
    bool result = VtkExporter::ExportSpanwiseDistributions(
        "test_loads.vtp", solver_.GetMesh(), results_, options);

    EXPECT_TRUE(result);
    EXPECT_TRUE(std::filesystem::exists("test_loads.vtp"));
}

TEST_F(VtkExportTest, ExportWakeGeometryCreatesFile) {
    if (!VtkExporter::IsAvailable()) {
        GTEST_SKIP() << "VTK not available";
    }

    VtkExportOptions options;
    options.wake_length_chords = 5.0;
    bool result = VtkExporter::ExportWakeGeometry(
        "test_wake.vtp", solver_.GetMesh(), solution_, condition_,
        aircraft_.GetReferenceChord(), options);

    EXPECT_TRUE(result);
    EXPECT_TRUE(std::filesystem::exists("test_wake.vtp"));
}

TEST_F(VtkExportTest, ExportVelocityFieldCreatesFile) {
    if (!VtkExporter::IsAvailable()) {
        GTEST_SKIP() << "VTK not available";
    }

    VtkExportOptions options;
    bool result = VtkExporter::ExportVelocityField(
        "test_velocity.vtp", solver_.GetMesh(), solution_, condition_, options);

    EXPECT_TRUE(result);
    EXPECT_TRUE(std::filesystem::exists("test_velocity.vtp"));
}

TEST_F(VtkExportTest, ExportAllCreatesAllFiles) {
    if (!VtkExporter::IsAvailable()) {
        GTEST_SKIP() << "VTK not available";
    }

    VtkExportOptions options;
    bool result = VtkExporter::ExportAll(
        "test_vtk_output", solver_.GetMesh(), solution_, results_, condition_,
        aircraft_.GetReferenceChord(), options);

    EXPECT_TRUE(result);
    EXPECT_TRUE(std::filesystem::exists("test_vtk_output/surface_mesh.vtu"));
    EXPECT_TRUE(std::filesystem::exists("test_vtk_output/spanwise_loads.vtp"));
    EXPECT_TRUE(std::filesystem::exists("test_vtk_output/wake.vtp"));
    EXPECT_TRUE(std::filesystem::exists("test_vtk_output/velocity_field.vtp"));
}

TEST_F(VtkExportTest, BinaryVsAsciiFormat) {
    if (!VtkExporter::IsAvailable()) {
        GTEST_SKIP() << "VTK not available";
    }

    VtkExportOptions binary_options;
    binary_options.binary_format = true;
    VtkExporter::ExportSurfaceMesh("test_binary.vtu", solver_.GetMesh(),
                                   solution_, binary_options);

    VtkExportOptions ascii_options;
    ascii_options.binary_format = false;
    VtkExporter::ExportSurfaceMesh("test_ascii.vtu", solver_.GetMesh(),
                                   solution_, ascii_options);

    // Binary should generally be smaller than ASCII for same data
    auto binary_size = std::filesystem::file_size("test_binary.vtu");
    auto ascii_size = std::filesystem::file_size("test_ascii.vtu");
    EXPECT_LT(binary_size, ascii_size);

    std::filesystem::remove("test_binary.vtu");
    std::filesystem::remove("test_ascii.vtu");
}

TEST_F(VtkExportTest, InvalidSolutionReturnsFalse) {
    if (!VtkExporter::IsAvailable()) {
        GTEST_SKIP() << "VTK not available";
    }

    SolutionState invalid_solution;
    invalid_solution.is_valid = false;

    VtkExportOptions options;
    bool result = VtkExporter::ExportSurfaceMesh(
        "test_invalid.vtu", solver_.GetMesh(), invalid_solution, options);
    EXPECT_FALSE(result);
    EXPECT_FALSE(std::filesystem::exists("test_invalid.vtu"));
}

TEST_F(VtkExportTest, ExportOptionsDefaultsAreReasonable) {
    VtkExportOptions options;

    EXPECT_TRUE(options.export_cp);
    EXPECT_TRUE(options.export_circulation);
    EXPECT_FALSE(options.export_normals);
    EXPECT_FALSE(options.export_indices);
    EXPECT_TRUE(options.export_wake);
    EXPECT_GT(options.wake_length_chords, 0.0);
    EXPECT_TRUE(options.export_induced_velocity);
    EXPECT_TRUE(options.export_spanwise_distributions);
    EXPECT_TRUE(options.binary_format);
}

}  // namespace
}  // namespace pegasus
