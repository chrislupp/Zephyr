cmake_minimum_required(VERSION 3.16)
project(VLM
    VERSION 1.0.0
    DESCRIPTION "Vortex Lattice Method Aerodynamics Library"
    LANGUAGES CXX
)

# C++17 required
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build options
option(VLM_BUILD_TESTS "Build unit tests" ON)
option(VLM_BUILD_EXAMPLES "Build example programs" ON)
option(VLM_ENABLE_OPENMP "Enable OpenMP parallelization" ON)

# Export compile commands for IDE support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Find dependencies
find_package(Eigen3 5.0 CONFIG REQUIRED)

# Collect source files
file(GLOB_RECURSE VLM_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
)

# Main library target
add_library(vlm ${VLM_SOURCES})

target_include_directories(vlm
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(vlm
    PUBLIC
        Eigen3::Eigen
)

# Compiler warnings
target_compile_options(vlm PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

# OpenMP support
if(VLM_ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(vlm PUBLIC OpenMP::OpenMP_CXX)
        target_compile_definitions(vlm PUBLIC VLM_USE_OPENMP)
    endif()
endif()

# VTK support (optional)
option(VLM_ENABLE_VTK "Enable VTK visualization export" ON)
if(VLM_ENABLE_VTK)
    find_package(VTK QUIET COMPONENTS
        CommonCore
        CommonDataModel
        IOXML
    )
    if(VTK_FOUND)
        message(STATUS "VTK found: ${VTK_VERSION}")
        target_link_libraries(vlm PUBLIC ${VTK_LIBRARIES})
        target_compile_definitions(vlm PUBLIC VLM_USE_VTK)
        # VTK 9+ requires module autoinit
        if(VTK_VERSION VERSION_GREATER_EQUAL "9.0")
            vtk_module_autoinit(
                TARGETS vlm
                MODULES ${VTK_LIBRARIES}
            )
        endif()
    else()
        message(STATUS "VTK not found - visualization export disabled")
    endif()
endif()

# Enzyme automatic differentiation support (optional)
option(VLM_ENABLE_ENZYME "Enable Enzyme automatic differentiation" OFF)
if(VLM_ENABLE_ENZYME)
    # Enzyme requires Clang compiler
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        message(FATAL_ERROR "Enzyme requires Clang compiler. "
            "Set CMAKE_CXX_COMPILER to clang++ or use -DVLM_ENABLE_ENZYME=OFF")
    endif()

    # Get LLVM version from Clang
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} --version
        OUTPUT_VARIABLE CLANG_VERSION_OUTPUT
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    string(REGEX MATCH "version ([0-9]+)" _ ${CLANG_VERSION_OUTPUT})
    set(LLVM_VERSION_MAJOR ${CMAKE_MATCH_1})
    message(STATUS "Detected LLVM/Clang major version: ${LLVM_VERSION_MAJOR}")

    # Find Enzyme plugin
    # Search common locations: user-specified, Homebrew, system
    find_library(ENZYME_PLUGIN
        NAMES ClangEnzyme-${LLVM_VERSION_MAJOR} ClangEnzyme
        HINTS
            ${ENZYME_DIR}
            $ENV{ENZYME_DIR}
            /opt/homebrew/lib
            /usr/local/lib
            /usr/lib
        DOC "Path to Enzyme Clang plugin"
    )

    if(ENZYME_PLUGIN)
        message(STATUS "Enzyme found: ${ENZYME_PLUGIN}")

        # Add Enzyme compiler flags
        # -fplugin loads the Enzyme LLVM pass
        # Optimization flags help Enzyme work correctly
        target_compile_options(vlm PUBLIC
            -fplugin=${ENZYME_PLUGIN}
            -O2
            -fno-vectorize
            -fno-slp-vectorize
            -fno-unroll-loops
        )

        # Define VLM_USE_ENZYME so code can conditionally use AD features
        target_compile_definitions(vlm PUBLIC VLM_USE_ENZYME)

        # Note: -fno-exceptions removed to maintain ABI compatibility with gtest
    else()
        message(WARNING "Enzyme not found - automatic differentiation disabled. "
            "Install with: brew install enzyme (macOS) or set ENZYME_DIR")
        set(VLM_ENABLE_ENZYME OFF)
    endif()
endif()

# Tests with Google Test
if(VLM_BUILD_TESTS)
    enable_testing()
    include(FetchContent)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0
    )
    # Prevent overriding parent project's compiler/linker settings on Windows
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
    add_subdirectory(tests)
endif()

# Examples
if(VLM_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# Installation
include(GNUInstallDirs)
install(TARGETS vlm
    EXPORT vlm-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT vlm-targets
    FILE VLMTargets.cmake
    NAMESPACE VLM::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/VLM
)
