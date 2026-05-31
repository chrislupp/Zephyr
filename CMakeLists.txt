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

cmake_minimum_required(VERSION 3.16)
project(Pegasus
    VERSION 0.2.0
    DESCRIPTION "Vortex Lattice Method Aerodynamics Library"
    LANGUAGES CXX
)

# C++17 required
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build options
option(PEGASUS_BUILD_TESTS "Build unit tests" ON)
option(PEGASUS_BUILD_EXAMPLES "Build example programs" ON)
option(PEGASUS_ENABLE_OPENMP "Enable OpenMP parallelization" OFF)
option(PEGASUS_ENABLE_VISCOUS_COUPLING
    "Build the viscous-inviscid coupling module (requires Scirocco/Trysil)"
    OFF)
set(PEGASUS_SCIROCCO_SOURCE_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/../ProjectTrysil/scirocco"
    CACHE PATH "Path to the Scirocco/Trysil source tree")

# Export compile commands for IDE support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Find dependencies
find_package(Eigen3 5.0 CONFIG REQUIRED)

# Collect source files (exclude philote, built separately; exclude viscous,
# added conditionally below)
file(GLOB_RECURSE PEGASUS_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
)
list(FILTER PEGASUS_SOURCES EXCLUDE REGEX ".*/src/philote/.*")
list(FILTER PEGASUS_SOURCES EXCLUDE REGEX ".*/src/viscous/.*")

# Main library target
add_library(pegasus ${PEGASUS_SOURCES})

target_include_directories(pegasus
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(pegasus
    PUBLIC
        Eigen3::Eigen
)

# Compiler warnings
target_compile_options(pegasus PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

# OpenMP support
if(PEGASUS_ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(pegasus PUBLIC OpenMP::OpenMP_CXX)
        target_compile_definitions(pegasus PUBLIC PEGASUS_USE_OPENMP)
    endif()
endif()

# VTK support (optional)
option(PEGASUS_ENABLE_VTK "Enable VTK visualization export" ON)
if(PEGASUS_ENABLE_VTK)
    find_package(VTK QUIET COMPONENTS
        CommonCore
        CommonDataModel
        IOXML
    )
    if(VTK_FOUND)
        message(STATUS "VTK found: ${VTK_VERSION}")
        target_link_libraries(pegasus PUBLIC ${VTK_LIBRARIES})
        target_compile_definitions(pegasus PUBLIC PEGASUS_USE_VTK)
        # VTK 9+ requires module autoinit
        if(VTK_VERSION VERSION_GREATER_EQUAL "9.0")
            vtk_module_autoinit(
                TARGETS pegasus
                MODULES ${VTK_LIBRARIES}
            )
        endif()
    else()
        message(STATUS "VTK not found - visualization export disabled")
    endif()
endif()

# Viscous-inviscid coupling with Scirocco/Trysil (optional)
if(PEGASUS_ENABLE_VISCOUS_COUPLING)
    if(NOT EXISTS "${PEGASUS_SCIROCCO_SOURCE_DIR}/CMakeLists.txt")
        message(FATAL_ERROR
            "PEGASUS_ENABLE_VISCOUS_COUPLING=ON but Scirocco was not found "
            "at ${PEGASUS_SCIROCCO_SOURCE_DIR}. Set PEGASUS_SCIROCCO_SOURCE_DIR "
            "to point at a Scirocco/Trysil source tree."
        )
    endif()

    # Suppress Scirocco's tests / examples / Philote / HDF5 / native-march when
    # consumed as a subdirectory; we only need its library target.
    set(SCIROCCO_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(SCIROCCO_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(SCIROCCO_BUILD_PHILOTE OFF CACHE BOOL "" FORCE)
    set(SCIROCCO_HDF5_OUTPUT OFF CACHE BOOL "" FORCE)
    set(SCIROCCO_NATIVE_ARCH OFF CACHE BOOL "" FORCE)
    set(SCIROCCO_BUILD_PROBES OFF CACHE BOOL "" FORCE)
    set(SCIROCCO_BUILD_SENSITIVITY_DEBUG_TESTS OFF CACHE BOOL "" FORCE)

    add_subdirectory(${PEGASUS_SCIROCCO_SOURCE_DIR}
                     ${CMAKE_BINARY_DIR}/_deps/scirocco-build
                     EXCLUDE_FROM_ALL)

    # Collect and add the viscous-coupling sources to the main library
    file(GLOB_RECURSE PEGASUS_VISCOUS_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/viscous/*.cpp
    )
    target_sources(pegasus PRIVATE ${PEGASUS_VISCOUS_SOURCES})
    target_link_libraries(pegasus PUBLIC scirocco)
    target_compile_definitions(pegasus PUBLIC PEGASUS_HAS_VISCOUS_COUPLING)
endif()


# Tests with Google Test
if(PEGASUS_BUILD_TESTS)
    enable_testing()
    include(FetchContent)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.15.2
    )
    # Prevent overriding parent project's compiler/linker settings on Windows
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
    add_subdirectory(tests)
endif()

# Examples
if(PEGASUS_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# Philote MDO discipline server (optional)
option(PEGASUS_BUILD_PHILOTE "Build Philote MDO discipline server" OFF)
if(PEGASUS_BUILD_PHILOTE)
    find_package(Protobuf CONFIG QUIET)
    if(NOT Protobuf_FOUND)
        find_package(Protobuf REQUIRED)
    endif()
    find_package(gRPC CONFIG REQUIRED)
    find_package(PhiloteCpp REQUIRED)
    # Ensure correct protobuf headers are found before any Homebrew system includes
    get_target_property(_proto_inc protobuf::libprotobuf INTERFACE_INCLUDE_DIRECTORIES)

    add_executable(pegasus_philote_server
        src/philote/vlm_discipline.cpp
        src/philote/vlm_server.cpp
    )
    target_link_libraries(pegasus_philote_server PRIVATE pegasus PhiloteCpp::PhiloteCpp)
    if(_proto_inc)
        target_include_directories(pegasus_philote_server BEFORE PRIVATE ${_proto_inc})
    endif()

    add_executable(mesh_philote_server
        src/philote/mesh_discipline.cpp
        src/philote/mesh_server.cpp
    )
    target_link_libraries(mesh_philote_server PRIVATE pegasus PhiloteCpp::PhiloteCpp)
    if(_proto_inc)
        target_include_directories(mesh_philote_server BEFORE PRIVATE ${_proto_inc})
    endif()

    # Philote discipline tests
    if(PEGASUS_BUILD_TESTS)
        add_executable(pegasus_philote_tests
            tests/unit/mesh_discipline_test.cpp
            src/philote/mesh_discipline.cpp
        )
        target_link_libraries(pegasus_philote_tests
            PRIVATE pegasus PhiloteCpp::PhiloteCpp GTest::gtest_main)
        if(_proto_inc)
            target_include_directories(pegasus_philote_tests BEFORE PRIVATE ${_proto_inc})
        endif()
        include(GoogleTest)
        gtest_discover_tests(pegasus_philote_tests)
    endif()
endif()

# Installation
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS pegasus
    EXPORT pegasus-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT pegasus-targets
    FILE PegasusTargets.cmake
    NAMESPACE Pegasus::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Pegasus
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/PegasusConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/PegasusConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Pegasus
)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/PegasusConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/PegasusConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/PegasusConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Pegasus
)
