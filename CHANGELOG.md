# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `PerVortexForces` struct and `ForceCalculator::ComputePerVortexForces()` for raw per-vortex body-frame forces without strip aggregation (Phase 0A of Perseus aerostructural coupling)
- `per_vortex_force_circulation_vjp()` analytical gradient function for the VLM adjoint: computes (dF/dGamma)^T * seed including both direct and induced-velocity contributions (Phase 0B of Perseus coupling)
- `forward_jvp.h`/`forward_jvp.cpp` module with forward-mode Jacobian-vector products through the VLM chain (mesh nodes → vortex geometry → AIC/b → Γ → per-vortex forces) for condensed-Jacobian aerostructural assembly (Phase 0C of Perseus coupling)
- `VortexLatticeSolver::SolveForward()` exposes the cached AIC LU factorization for forward-mode circulation sensitivities used by the JVP chain
- CMake package config (`cmake/PegasusConfig.cmake.in`, `PegasusConfigVersion.cmake`) installed under `lib/cmake/Pegasus` so downstream projects can consume Pegasus via `find_package(Pegasus)`

- `Section::SetTwistOffset()` / `GetTwistOffset()` / `GetEffectiveTwistAngle()` for injecting a Δα correction on top of the geometric twist without overwriting it; consumed by mesh generation as the effective twist. Used by the viscous-inviscid coupler to inject per-strip Newton corrections.
- `PEGASUS_ENABLE_VISCOUS_COUPLING` CMake option (default OFF) that pulls in Scirocco/Trysil as a subdirectory dependency for the viscous-inviscid coupling module. Path is controlled by `PEGASUS_SCIROCCO_SOURCE_DIR` (defaults to a sibling `ProjectTrysil/scirocco` directory). When ON, defines `PEGASUS_HAS_VISCOUS_COUPLING` and links `pegasus` against `scirocco`.
- `pegasus::viscous::StripPolarSolver`: per-strip wrapper around a Scirocco 2D viscous airfoil analysis with warm-start support (caches converged BL state across Newton iterations via Scirocco's `oper.initbl = false`). Exposes Cl/Cd/Cm at radians and per-radian dCl/dα and dCd/dα for the coupling Jacobian.
- `pegasus::viscous::ViscousVlmCoupler`: Newton-based VLM ↔ 2D viscous coupler. At each outer iteration applies a per-strip Δα to `Strip::twist_angle`, rebuilds the AIC and re-solves the VLM, runs Scirocco at each strip's local effective α, and forms the per-strip Cl mismatch residual. v1 builds the Jacobian by finite differences and returns viscous-corrected `AerodynamicResults` including a non-zero `CDv`. End-to-end convergence is exercised on a small rectangular wing test.
- `ViscousVlmCoupler::ComputeAlphaDerivatives()` returns total dCL/dα∞ and dCD/dα∞ for the converged coupled problem. Applies the implicit-function theorem `J du/dα = -∂R/∂α` followed by `dF/dα = ∂F/∂α + ∂F/∂u · du/dα`, using shared FD sweeps over Δα for both `J = ∂R/∂u` and `∂F/∂u` to amortize the cost. Validated against a 2-point FD of fully re-converged solves at α_base ± h within 5% on CL and 10% on CD.

### Fixed
- Qualify bare `size_t` as `std::size_t` in `per_vortex_force_circulation_vjp()` so it builds with stricter compiler frontends (e.g. AppleClang)

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
