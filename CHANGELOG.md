# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- A safe Rust 2024 library implementing geometry, structured meshing, the
  primal VLM solver, near- and far-field loads, stability derivatives, trim,
  mass properties, 4-state and 12-state eigenmodes, and analytical adjoints.
- AVL geometry and AVL mass parsers with typed, source-aware errors.
- Pure-Rust `.vtu` and `.vtp` export through the default `vtk` feature,
  including ASCII/binary output and atomic file replacement.
- Native Rust examples, a Criterion scaling benchmark, Rust API guides, and
  Linux/macOS/Windows CI.
- A committed v0.2.0 numerical oracle and a manifest covering all 488
  non-Philote CTest cases.

### Changed

- The package is now `pegasus-vlm`, imported as `pegasus_vlm`, and requires
  Rust 1.93.1.
- The public API uses owned Rust data, checked slices and indices,
  `PegasusError`, grouped `Coefficients`, and `AeroOutput`.
- A single partial-pivot AIC factorization is cached and reused for primal,
  transposed-adjoint, unit-response, and repeated-flight-condition solves,
  avoiding duplicate cubic work while preserving single-threaded execution.
- Flight-condition, aircraft-section, and direct mesh-node gradients use
  analytical influence/force kernels and adjoint contractions; the Mach
  derivative no longer requires perturbed full solves.
- The missing fixed-circulation alpha term in the adjoint force derivative now
  includes freestream and body-to-stability-axis effects. The former disabled
  CL validation is enabled and extended across all six primary coefficients.

### Removed

- The C++ API, C ABI expectations, CMake build, system Eigen/VTK/OpenMP
  discovery, GoogleTest, and C++ examples.
- Philote servers, `MeshDiscipline`, `VlmDiscipline`, and Philote-dependent
  OpenMDAO workflows. These integrations are deferred to a later release.

## [0.2.0] - 2026-05-04

### Added
- `OutputSelection` struct and `GradientCalculator::SetOutputSelection()` to select which aerodynamic outputs (CL, CD, CY, Cl, Cm, Cn, CDi) have gradients computed, reducing adjoint solves and O(n²) velocity-kernel operations for unused outputs
- `MeshDiscipline` now computes `reference_area`, `reference_span`, and `reference_chord` (MAC) as outputs with analytical partials
- Comprehensive unit tests for geometry (Surface, Aircraft, Airfoil, Section, AVL reader), math utilities, influence matrix, trim solver, flight constraints, aerodynamic results, eigenmode analysis, mass properties, vortex lattice, and VLM solver — line coverage increased from ~80% to ~92%
- Philote and VTK/visualization source exclusion from Codecov coverage reporting
- Optimization convergence history plot (CD vs major iteration) in OpenMDAO optimization example

### Changed
- `MeshDiscipline` replaces absolute `yle` input with scalar `halfspan` and normalized `eta` array parameterization, adding wingspan as a design variable for gradient-based optimization with analytical partials
- License changed from AGPL-3.0-or-later to GPL-3.0-or-later

### Removed
- God header `pegasus/pegasus.h` — consumers now include specific headers directly
- Autodiff tape and dual number headers (`gradients/autodiff_tape.h`, `gradients/dual_number.h`)

### Fixed
- CI job names now match branch protection required checks (`format`, `build`)

## [0.1.0] - 2026-05-03

### Added
- Analytical adjoint derivatives replacing Enzyme AD — exact VJPs for Biot-Savart, AIC matrix, force kernel, and Trefftz plane
- Direct mesh input API (`MeshInput`) as alternative to section-based Aircraft geometry
- Analytical mesh node gradients for all aerodynamic coefficients (CL, CD, CY, Cl, Cm, Cn)
- Philote MDO disciplines (`VlmDiscipline`, `MeshDiscipline`) with gRPC-based remote execution and full `ComputePartials`
- OpenMDAO examples for VLM analysis, chord optimization, and performance benchmarks vs OpenAeroStruct
- GitHub Actions CI pipeline with format checking, build, static analysis, and code coverage
- Docusaurus documentation site with getting started guides, user guide, examples, and API reference
- GitHub Actions workflow for deploying docs to GitHub Pages on push to `develop`
- Automated release workflow with label-driven semantic versioning
- `.clang-format` configuration (Google-based, C++17, 4-space indent, 80-column limit)
- `.clang-tidy` static analysis configuration (bugprone, modernize, performance checks)
- Codecov integration for code coverage tracking

### Changed
- Reorganized examples into subdirectories
- Refactored `ForceCalculator` and `AdjointSolver` to accept `ReferenceQuantities` directly

### Removed
- Enzyme automatic differentiation compiler plugin dependency and all AD kernels (~3,200 lines removed)

### Fixed
- Docs homepage redirect to use relative path instead of absolute (avoiding baseUrl bypass)

[0.1.0]: https://github.com/chrislupp/Pegasus/releases/tag/v0.1.0
[Unreleased]: https://github.com/chrislupp/Pegasus/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/chrislupp/Pegasus/releases/tag/v0.2.0
