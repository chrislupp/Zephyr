include(GoogleTest)

# Unit tests
add_executable(vlm_unit_tests
    unit/cubic_spline_test.cpp
    unit/geometry_test.cpp
    unit/vortex_lattice_test.cpp
    unit/influence_matrix_test.cpp
    unit/solver_test.cpp
    unit/force_calculator_test.cpp
    unit/stability_derivatives_test.cpp
    unit/trim_solver_test.cpp
    unit/eigenmode_test.cpp
    unit/autodiff_test.cpp
    unit/mass_reader_test.cpp
    unit/vtk_export_test.cpp
    unit/control_surface_test.cpp
)

target_link_libraries(vlm_unit_tests
    PRIVATE
        vlm
        GTest::gtest_main
)

gtest_discover_tests(vlm_unit_tests)

# Integration tests
add_executable(vlm_integration_tests
    integration/rectangular_wing_test.cpp
)

target_link_libraries(vlm_integration_tests
    PRIVATE
        vlm
        GTest::gtest_main
)

gtest_discover_tests(vlm_integration_tests
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# Copy test fixtures to build directory
file(COPY fixtures DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
