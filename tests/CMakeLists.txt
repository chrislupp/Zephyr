include(GoogleTest)

# Unit tests
add_executable(vlm_unit_tests
    unit/cubic_spline_test.cpp
    unit/geometry_test.cpp
)

target_link_libraries(vlm_unit_tests
    PRIVATE
        vlm
        GTest::gtest_main
)

gtest_discover_tests(vlm_unit_tests)

# Copy test fixtures to build directory
file(COPY fixtures DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
