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

# Collect source files (exclude philote subdirectory, built separately)
file(GLOB_RECURSE VLM_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
)
list(FILTER VLM_SOURCES EXCLUDE REGEX ".*/src/philote/.*")

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


# Tests with Google Test
if(VLM_BUILD_TESTS)
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
if(VLM_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# Philote MDO discipline server (optional)
option(VLM_BUILD_PHILOTE "Build Philote MDO discipline server" OFF)
if(VLM_BUILD_PHILOTE)
    find_package(Protobuf CONFIG QUIET)
    if(NOT Protobuf_FOUND)
        find_package(Protobuf REQUIRED)
    endif()
    find_package(gRPC CONFIG REQUIRED)
    find_package(PhiloteCpp REQUIRED)
    add_executable(vlm_philote_server
        src/philote/vlm_discipline.cpp
        src/philote/vlm_server.cpp
    )
    target_link_libraries(vlm_philote_server PRIVATE vlm PhiloteCpp::PhiloteCpp)
    # Ensure correct protobuf headers are found before any Homebrew system includes
    get_target_property(_proto_inc protobuf::libprotobuf INTERFACE_INCLUDE_DIRECTORIES)
    if(_proto_inc)
        target_include_directories(vlm_philote_server BEFORE PRIVATE ${_proto_inc})
    endif()
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
