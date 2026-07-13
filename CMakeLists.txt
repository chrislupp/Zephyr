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
    "Build the viscous-inviscid coupling module (requires at least one 2D backend)"
    OFF)
option(PEGASUS_ENABLE_SCIROCCO
    "Enable the Scirocco/Trysil 2D backend for viscous coupling"
    ON)
option(PEGASUS_ENABLE_BORA
    "Enable the Bora (NeuralFoil) 2D backend for viscous coupling"
    OFF)
option(PEGASUS_ENABLE_SCIROCCO_TESTS
    "Build and register the Scirocco viscous tests (slow: several minutes)"
    OFF)
set(PEGASUS_SCIROCCO_SOURCE_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/../ProjectTrysil/scirocco"
    CACHE PATH "Path to the Scirocco/Trysil source tree")
set(PEGASUS_BORA_SOURCE_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/../airfoils/Bora"
    CACHE PATH
    "Path to the Bora source tree (used if no installed Bora package is found)")

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

# BLAS/LAPACK-backed Eigen dense kernels (optional). Speeds up the LU
# factorizations and large matrix products; results agree with Eigen's
# built-in kernels to machine precision but are not bitwise identical.
option(PEGASUS_ENABLE_BLAS_LAPACK
    "Back Eigen dense factorizations/products with BLAS + LAPACKE"
    OFF)
if(PEGASUS_ENABLE_BLAS_LAPACK)
    find_library(PEGASUS_OPENBLAS_LIBRARY NAMES openblas
        HINTS /opt/homebrew/opt/openblas/lib /usr/local/opt/openblas/lib)
    find_path(PEGASUS_LAPACKE_INCLUDE_DIR lapacke.h
        HINTS /opt/homebrew/opt/openblas/include
              /usr/local/opt/openblas/include)
    if(NOT PEGASUS_OPENBLAS_LIBRARY OR NOT PEGASUS_LAPACKE_INCLUDE_DIR)
        message(FATAL_ERROR
            "PEGASUS_ENABLE_BLAS_LAPACK=ON requires OpenBLAS with LAPACKE "
            "(e.g. `brew install openblas`)."
        )
    endif()
    message(STATUS "Eigen BLAS/LAPACKE backend: ${PEGASUS_OPENBLAS_LIBRARY}")
    # PUBLIC: Eigen's inline kernels must compile identically in every TU
    # that links pegasus (mixing would violate the ODR).
    target_compile_definitions(pegasus PUBLIC EIGEN_USE_BLAS EIGEN_USE_LAPACKE)
    target_include_directories(pegasus PUBLIC ${PEGASUS_LAPACKE_INCLUDE_DIR})
    target_link_libraries(pegasus PUBLIC ${PEGASUS_OPENBLAS_LIBRARY})
endif()

# Apple Accelerate LAPACK for the AIC LU factorization (optional, macOS
# only). Same partial-pivoted LU algorithm as Eigen's, vendor-tuned (AMX):
# ~2.7x faster at 600 panels, ~6.8x at 2000. Backward error is in the same
# machine-precision class; results are not bitwise identical to Eigen's LU.
option(PEGASUS_ENABLE_ACCELERATE
    "Use Apple Accelerate LAPACK for the AIC LU factorization"
    OFF)
if(PEGASUS_ENABLE_ACCELERATE)
    if(NOT APPLE)
        message(FATAL_ERROR "PEGASUS_ENABLE_ACCELERATE requires macOS.")
    endif()
    if(PEGASUS_ENABLE_BLAS_LAPACK)
        message(FATAL_ERROR
            "PEGASUS_ENABLE_ACCELERATE and PEGASUS_ENABLE_BLAS_LAPACK are "
            "mutually exclusive (conflicting BLAS/LAPACK symbols)."
        )
    endif()
    target_compile_definitions(pegasus PUBLIC PEGASUS_USE_ACCELERATE)
    target_link_libraries(pegasus PUBLIC "-framework Accelerate")
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

# Viscous-inviscid coupling (optional). The coupler core is backend-
# agnostic; at least one 2D airfoil backend (Scirocco and/or Bora) must be
# enabled alongside it.
if(PEGASUS_ENABLE_VISCOUS_COUPLING)
    if(NOT PEGASUS_ENABLE_SCIROCCO AND NOT PEGASUS_ENABLE_BORA)
        message(FATAL_ERROR
            "PEGASUS_ENABLE_VISCOUS_COUPLING=ON requires at least one 2D "
            "backend: set PEGASUS_ENABLE_SCIROCCO=ON and/or "
            "PEGASUS_ENABLE_BORA=ON."
        )
    endif()

    # Backend-agnostic coupler core.
    target_sources(pegasus PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/viscous/coupler.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/viscous/coupling_station.cpp
    )
    target_compile_definitions(pegasus PUBLIC PEGASUS_HAS_VISCOUS_COUPLING)

    # Scirocco/Trysil viscous panel + boundary-layer backend.
    if(PEGASUS_ENABLE_SCIROCCO)
        if(NOT EXISTS "${PEGASUS_SCIROCCO_SOURCE_DIR}/CMakeLists.txt")
            message(FATAL_ERROR
                "PEGASUS_ENABLE_SCIROCCO=ON but Scirocco was not found "
                "at ${PEGASUS_SCIROCCO_SOURCE_DIR}. Set "
                "PEGASUS_SCIROCCO_SOURCE_DIR to point at a Scirocco/Trysil "
                "source tree."
            )
        endif()

        # Suppress Scirocco's tests / examples / Philote / HDF5 / native-march
        # when consumed as a subdirectory; we only need its library target.
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

        target_sources(pegasus PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src/viscous/strip_polar.cpp
        )
        target_link_libraries(pegasus PUBLIC scirocco)
        target_compile_definitions(pegasus PUBLIC PEGASUS_HAS_SCIROCCO)
    endif()

    # Bora (NeuralFoil) neural-surrogate backend. Prefer an installed
    # package; fall back to a source-tree add_subdirectory. Bora is a C++20
    # library, so its wrapper translation unit is compiled separately as an
    # object library — the rest of Pegasus (and its public headers) stay
    # C++17.
    if(PEGASUS_ENABLE_BORA)
        find_package(Bora CONFIG QUIET)
        if(NOT Bora_FOUND)
            if(NOT EXISTS "${PEGASUS_BORA_SOURCE_DIR}/CMakeLists.txt")
                message(FATAL_ERROR
                    "PEGASUS_ENABLE_BORA=ON but no installed Bora package "
                    "was found and no source tree exists at "
                    "${PEGASUS_BORA_SOURCE_DIR}. Install Bora or set "
                    "PEGASUS_BORA_SOURCE_DIR."
                )
            endif()
            # Suppress Bora's tests / examples / Philote / native-march when
            # consumed as a subdirectory; we only need its library target.
            set(BUILD_TESTS OFF CACHE BOOL "" FORCE)
            set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
            set(BORA_NATIVE_ARCH OFF CACHE BOOL "" FORCE)
            set(BORA_BUILD_PHILOTE OFF CACHE BOOL "" FORCE)
            add_subdirectory(${PEGASUS_BORA_SOURCE_DIR}
                             ${CMAKE_BINARY_DIR}/_deps/bora-build
                             EXCLUDE_FROM_ALL)
        else()
            message(STATUS "Bora found: ${Bora_VERSION}")
        endif()

        add_library(pegasus_bora_backend OBJECT
            ${CMAKE_CURRENT_SOURCE_DIR}/src/viscous/bora_polar.cpp
        )
        target_compile_features(pegasus_bora_backend PRIVATE cxx_std_20)
        target_include_directories(pegasus_bora_backend PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/include
        )
        target_link_libraries(pegasus_bora_backend PRIVATE
            Bora::Bora
            Eigen3::Eigen
        )
        target_compile_options(pegasus_bora_backend PRIVATE
            $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
            $<$<CXX_COMPILER_ID:MSVC>:/W4>
        )

        target_sources(pegasus PRIVATE $<TARGET_OBJECTS:pegasus_bora_backend>)
        target_link_libraries(pegasus PUBLIC Bora::Bora)
        target_compile_definitions(pegasus PUBLIC PEGASUS_HAS_BORA)
    endif()
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

    # Viscous-VLM discipline server (requires the viscous coupling module).
    if(PEGASUS_ENABLE_VISCOUS_COUPLING)
        add_executable(viscous_vlm_philote_server
            src/philote/viscous_vlm_discipline.cpp
            src/philote/viscous_vlm_server.cpp
        )
        target_link_libraries(viscous_vlm_philote_server
            PRIVATE pegasus PhiloteCpp::PhiloteCpp)
        if(_proto_inc)
            target_include_directories(viscous_vlm_philote_server BEFORE PRIVATE
                ${_proto_inc})
        endif()
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

# Record which optional dependencies are in the installed link interface so
# PegasusConfig.cmake can re-find them for downstream projects
if(PEGASUS_ENABLE_OPENMP AND OpenMP_CXX_FOUND)
    set(PEGASUS_CONFIG_NEEDS_OPENMP TRUE)
else()
    set(PEGASUS_CONFIG_NEEDS_OPENMP FALSE)
endif()
if(PEGASUS_ENABLE_VTK AND VTK_FOUND)
    set(PEGASUS_CONFIG_NEEDS_VTK TRUE)
else()
    set(PEGASUS_CONFIG_NEEDS_VTK FALSE)
endif()
# Bora can only be re-found downstream when it was consumed as an installed
# package (a source-tree add_subdirectory has no package config to find).
if(PEGASUS_ENABLE_VISCOUS_COUPLING AND PEGASUS_ENABLE_BORA AND Bora_FOUND)
    set(PEGASUS_CONFIG_NEEDS_BORA TRUE)
else()
    set(PEGASUS_CONFIG_NEEDS_BORA FALSE)
endif()

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
