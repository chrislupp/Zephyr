---
sidebar_position: 1
---

# Installation

## Requirements

- **CMake** 3.16 or higher
- **C++17** compiler (GCC 8+, Clang 7+, MSVC 2019+)
- **Eigen** 5.0+ linear algebra library

### Optional Dependencies

- **OpenMP** — parallel computation support (auto-detected)
- **VTK** 8.0+ — visualization export (auto-detected)
- **Protobuf** and **gRPC** — required only for Philote MDO integration
- **PhiloteCpp** — required only for Philote MDO integration

## Installing Eigen 5.0

Eigen 5.0 must be built from source:

```bash
git clone https://gitlab.com/libeigen/eigen.git
cd eigen
git checkout master
cmake -B build -DBUILD_TESTING=OFF -DEIGEN_BUILD_BLAS=OFF -DEIGEN_BUILD_LAPACK=OFF
cmake --install build
```

To install to a custom prefix, add `-DCMAKE_INSTALL_PREFIX=/path/to/prefix` to the configure step.

## Building VLM

```bash
git clone https://github.com/chrislupp/Pegasus.git
cd Pegasus
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

If Eigen is installed to a non-standard location:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/eigen
```

## Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `VLM_BUILD_TESTS` | `ON` | Build unit and integration tests |
| `VLM_BUILD_EXAMPLES` | `ON` | Build example programs |
| `VLM_ENABLE_OPENMP` | `ON` | Enable OpenMP parallelization |
| `VLM_ENABLE_VTK` | `ON` | Enable VTK visualization export |
| `VLM_BUILD_PHILOTE` | `OFF` | Build Philote MDO discipline server |

Example with custom options:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release \
    -DVLM_BUILD_TESTS=OFF \
    -DVLM_ENABLE_VTK=OFF
```

## Running Tests

```bash
cd build
ctest --output-on-failure
```

Or run test executables directly:

```bash
./build/tests/vlm_unit_tests
./build/tests/vlm_integration_tests
```

## Installing

```bash
cmake --install build
```

This installs headers, the library, and CMake package configuration files. Downstream projects can then use:

```cmake
find_package(VLM REQUIRED)
target_link_libraries(my_target PRIVATE VLM::vlm)
```
