# AI Agent Guidelines for Project Pegasus

## Project Overview

Pegasus is a Vortex Lattice Method aerodynamic analysis library written in C++17. It provides VLM solving, force computation, stability derivatives, trim analysis, eigenmode analysis, and analytical gradient computation. Ported from OptVL.

## Build System

### Requirements
- CMake 3.16+
- C++17 compatible compiler
- Eigen 5.0+

### Build Commands
```bash
# Configure (from repository root)
cmake -B build -DCMAKE_BUILD_TYPE=Release

# Build all
cmake --build build

# Build specific target
cmake --build build --target <target_name>

# Run all tests
ctest --test-dir build

# Run specific test executable
./build/tests/pegasus_unit_tests
./build/tests/pegasus_integration_tests
```

### Linting and Formatting
Before pushing any changes to the remote, always run formatting and linting checks:

```bash
# Check formatting (must pass before push)
find include src tests examples -name '*.h' -o -name '*.cpp' | xargs clang-format --dry-run --Werror

# Auto-fix formatting
find include src tests examples -name '*.h' -o -name '*.cpp' | xargs clang-format -i

# Run clang-tidy (requires cmake configure with compile commands)
cmake -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
run-clang-tidy -p build -header-filter='include/pegasus/' 'src/.*\.cpp$'
```

### Key CMake Options
- `PEGASUS_BUILD_TESTS=ON|OFF` - Enable/disable test building (default: ON)
- `PEGASUS_BUILD_EXAMPLES=ON|OFF` - Enable/disable example programs (default: ON)
- `PEGASUS_ENABLE_VTK=ON|OFF` - Enable VTK visualization export (default: ON)
- `PEGASUS_BUILD_PHILOTE=ON|OFF` - Build Philote-Cpp MDO discipline interface (default: OFF)

## Code Organization

```
include/pegasus/          # Public headers
  core/                   # Types, constants
  geometry/               # Airfoil, Section, Surface, Aircraft
  io/                     # AVL/mass file readers
  solver/                 # VLM solver, influence matrix
  aero/                   # Forces, stability derivatives
  trim/                   # Trim solver
  dynamics/               # Eigenmode analysis
  gradients/              # Analytical gradient computation
  utils/                  # Cubic spline, math utilities
src/                      # Implementation files
tests/                    # Unit and integration tests
examples/                 # Example programs
scripts/                  # Build and release scripts
docs/                     # Docusaurus documentation site
```

## CI/CD

### Build Workflow (`.github/workflows/build.yml`)
Runs on every push and PR:
- **Format**: clang-format-18 check on `include/` and `src/`
- **Clang-Tidy**: Static analysis on PRs to `main`/`develop`
- **Build**: GCC 14, Release, CTest
- **Coverage**: Clang 18, llvm-cov LCOV, Codecov upload

### Docs Workflow (`.github/workflows/docs.yml`)
Deploys Docusaurus docs to GitHub Pages on push to `develop` when `docs/**` or `include/**` change.

### Release Workflow (`.github/workflows/release.yml`)
Label-driven release pipeline. Triggered when a labeled PR is merged to `main`.

**Label system** (applied to PRs targeting `main`):

| Purpose | Labels (pick exactly one per row) |
|---------|-----------------------------------|
| Release type | `release` OR `prerelease` |
| Version bump | `major`, `minor`, or `patch` |
| Pre-release qualifier | `alpha`, `beta`, or `rc` (only with `prerelease`) |

### Changelog
`CHANGELOG.md` follows [Keep a Changelog](https://keepachangelog.com/) format. New changes go under `[Unreleased]`. The release workflow converts this to a versioned section automatically.

**Important**: EVERY code change MUST include a corresponding `CHANGELOG.md` update in the same commit. This applies to all additions, modifications, and removals — no exceptions. If something previously added is later removed, update the changelog to reflect that removal. Never commit code changes without updating the changelog.

## Coding Conventions

### Style
- Use 4-space indentation
- Opening braces on same line for control structures, new line for functions/classes
- Member variables use trailing underscore: `member_variable_`
- Prefer `auto` for complex types, explicit types for primitives
- Google-based clang-format with 80-column limit

### Namespace
- Primary namespace: `pegasus`
- Gradient computation sub-namespace: `pegasus::gradients`

### Copyright Header
All source files must include the copyright header:
```
/*
    Project Pegasus
    Version X.Y.Z

    Copyright YYYY, Christopher A. Lupp
    ...GPL-3.0 license block...
*/
```

## Testing

### Test Framework
Uses Google Test (gtest).

### Running Tests
```bash
# All tests
ctest --test-dir build

# Specific test with verbose output
./build/tests/pegasus_unit_tests --gtest_filter="*PatternHere*"
```

## Pull Requests

When creating PRs:
- Target `develop` for feature work, `main` only for releases
- Label with appropriate type: `enhancement`, `bug`, `documentation`, etc.

## Git Workflow

- **`main`**: Stable releases only. PRs from `develop` with release labels.
- **`develop`**: Integration branch. Feature branches merge here.
- **`feature/*`**: Feature branches off `develop`.

## License

GPL-3.0-or-later. All source files must include the license header.
