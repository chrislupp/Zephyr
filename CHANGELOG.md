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
- `ViscousVlmCoupler::ComputeAlphaDerivatives()` returns total dCL/dα∞ and dCD/dα∞ for the converged coupled problem. Applies the implicit-function theorem `J du/dα = -∂R/∂α` followed by `dF/dα = ∂F/∂α + ∂F/∂u · du/dα`. All four partials (J = ∂R/∂u, ∂F/∂u, ∂R/∂α, ∂F/∂α) are computed **analytically**: a single forward-mode propagation per perturbation direction through the VLM chain (closed-form normal derivative w.r.t. twist → δR_i = δn_i · v_total_i + n_i · δv_inf → AIC LU back-substitution for δΓ → δCl_strip and δα_eff). Scirocco's adjoint contributes the diagonal `dCl_2D/dα` and `dCd_2D/dα` per strip. Newton's Jacobian build cost drops from N+1 VLM + N²+N Scirocco solves per iter (FD) to 1 VLM + N Scirocco solves (analytical). Validated against fully-re-converged outer FD with central step 0.005°: analytical agrees to ~3e-4 relative on both dCL/dα and dCD/dα (test enforces 1e-3).
- `ViscousVlmDiscipline` Philote MDO discipline (`include/pegasus/philote/viscous_vlm_discipline.h`) exposing the coupled VLM ↔ 2D viscous solve as a single monolithic explicit component. Inputs match `VlmDiscipline` (mesh nodes, reference quantities, flight condition); viscous options (`viscous_airfoil_naca`, `viscous_reynolds`, etc.) are set via `SetOptions`. v1 outputs viscous-corrected CL and CD with declared partials `dCL/dα` and `dCD/dα` driven by `ComputeAlphaDerivatives`. Built only when both `PEGASUS_BUILD_PHILOTE` and `PEGASUS_ENABLE_VISCOUS_COUPLING` are ON. A matching `viscous_vlm_philote_server` executable hosts it over gRPC on port 50053 by default.

### Changed
- **viscous coupling v2**: redesigned the coupler around user-defined 2D viscous stations that are independent of the VLM strip grid. Each station carries its own airfoil coordinate matrix (2×N), panel count, and Reynolds number; stations are placed at user-specified eta values. Cubic-spline operators map VLM strip values ↔ station values along eta. Newton unknowns are now per-station Δα (small N). The per-strip α_eff projection uses the strip's intrinsic local frame (chord direction × spanwise tangent × normal), so swept and dihedral wings are absorbed naturally. The previous per-strip-only API is replaced; tests and the Philote discipline are updated to the station form. Strip Cl residual remains in stability-axis convention to preserve the monotonic Newton path with a single physical fixed point. Analytical chain through the spline operator + closed-form δn/δc_hat/δn_hat + cached AIC LU + Scirocco's dCl/dα and dCd/dα keeps the gradient computation analytical end-to-end; FD validation agrees to ~3e-6 relative on both dCL/dα and dCD/dα (limited by the FD reference's central-difference noise floor).

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
