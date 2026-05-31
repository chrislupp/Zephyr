# Project Pegasus
# Version 0.2.0
#
# Copyright 2024-2026, Christopher A. Lupp
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

include(GoogleTest)

# Unit tests
add_executable(pegasus_unit_tests
    unit/cubic_spline_test.cpp
    unit/geometry_test.cpp
    unit/vortex_lattice_test.cpp
    unit/influence_matrix_test.cpp
    unit/solver_test.cpp
    unit/force_calculator_test.cpp
    unit/stability_derivatives_test.cpp
    unit/trim_solver_test.cpp
    unit/eigenmode_test.cpp
    unit/mass_reader_test.cpp
    unit/vtk_export_test.cpp
    unit/control_surface_test.cpp
    unit/gradient_calculator_test.cpp
    unit/adjoint_solver_test.cpp
    unit/optvl_gradient_validation_test.cpp
    unit/mesh_input_test.cpp
    unit/mesh_node_gradient_test.cpp
    unit/body_test.cpp
    unit/analytical_aic_test.cpp
    unit/analytical_biot_savart_test.cpp
    unit/analytical_forces_test.cpp
    unit/analytical_trefftz_test.cpp
    unit/mesh_sensitivity_test.cpp
    unit/math_utilities_test.cpp
    unit/flight_constraints_test.cpp
    unit/section_test.cpp
    unit/aerodynamic_results_test.cpp
)

if(PEGASUS_ENABLE_VISCOUS_COUPLING)
    target_sources(pegasus_unit_tests PRIVATE
        unit/viscous_strip_polar_test.cpp
    )
endif()

target_link_libraries(pegasus_unit_tests
    PRIVATE
        pegasus
        GTest::gtest_main
)

gtest_discover_tests(pegasus_unit_tests)

# Integration tests
add_executable(pegasus_integration_tests
    integration/rectangular_wing_test.cpp
)

target_link_libraries(pegasus_integration_tests
    PRIVATE
        pegasus
        GTest::gtest_main
)

gtest_discover_tests(pegasus_integration_tests
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# Copy test fixtures to build directory
file(COPY fixtures DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
