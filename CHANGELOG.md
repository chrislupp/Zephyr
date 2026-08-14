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
- A committed v0.2.0 numerical oracle, and a frozen inventory of the 488
  non-Philote CTest cases in `tests/fixtures/parity/ctest-manifest.json`. The
  manifest records which cases existed and at what tolerance; it carries no
  expected values and is not itself a behavioural oracle. Numerical parity is
  asserted by `tests/oracle_parity.rs` against committed v0.2.0 fixtures.
- Opt-in reflection reuse for lattices built from `y_duplicated` surfaces, via
  `set_mirror_reuse` on the solver and the influence builder. The induced
  velocity is equivariant under the mirror `M = diag(1, -1, 1)`, so half the
  influence entries can be derived rather than evaluated; measured 1.15-1.34x.
  Off by default, and enabling it never changes results. Assembly verifies the
  reflection relation against the actual panel geometry before using it and
  falls back to evaluating every entry otherwise, so symmetry images, partially
  duplicated meshes, and mirror planes away from `y = 0` are all handled
  correctly rather than approximated. `mirror_reuse_applies` reports which path
  a given mesh will take.
- A non-default `parallel` feature that multithreads influence assembly over
  source panels via rayon. Each column writes a disjoint slice and reads
  nothing another column writes, so results are bit-identical to the serial
  build and independent of the thread count; below 128 panels it keeps to the
  serial driver, where distribution costs more than it saves. Measured 1.5-1.8x
  on an 8-performance-core M1 Pro. The default build remains single-threaded.
  `faer/rayon` is deliberately not enabled, since gemm's cache blocking differs
  between its sequential and parallel paths and would perturb the LU's last ULP.

### Changed

- Section-geometry gradients are now fully analytical. `compute_geometry_gradients`
  previously cloned the aircraft, perturbed each section by a hardcoded step,
  regenerated two lattices and central-differenced them; it now builds the
  lattice tangent in closed form through the new `section_geometry` module and
  feeds the same adjoint context. The interpolated-twist derivative is the exact
  chord-weighted form `(K dS - S dK) / (S² + K²)`, **not** the small-angle
  linearization used by the v0.2.0 C++ `PropagateTwistGradient`, which is exact
  only when the chords and twists match at both ends of an interval and is off
  by tens of percent on any tapered, twisted wing. `interpolated_twist_jvp` and
  `panel_normal_jvp` are exported so the kernels can be verified directly.
- `GeometryAdjointContext` rejects meshes carrying section drag polars rather
  than returning an inviscid-only geometry gradient for them. The viscous term
  needs strip-area and quarter-chord sensitivities that `LatticeDirections`
  does not yet represent, and a silently incomplete derivative is worse than a
  refusal. Flight-condition gradients are unaffected and do handle polars.
- `compute_flight_condition_gradients` assembles each output's circulation
  gradient once and performs a single batched transpose solve. It previously
  built every gradient twice and ran two adjoint solves per output — 2K
  assemblies and 2K solves for K outputs.
- `force_geometry_jvp` and `trefftz_geometry_jvp` return `Result` and reject
  direction slices that do not match the lattice. They previously truncated to
  the shortest slice or zero-filled the missing entries, returning a partial
  derivative indistinguishable from a correct one; `strip_chord_direction` was
  never length-checked at all.
- `control_circulation_derivative` reuses the cached factorization when the
  effective normals match the solved ones, instead of taking a fresh O(N³)
  decomposition of a matrix that was already factorized, and contracts the
  normal-derivative matrix into the circulation directly rather than
  materializing a dense N×N to use once.
- Spanwise discretization is computed once by `Surface` and shared with the
  mesher, replacing a second copy of the arc accumulation and per-interval
  rounding. `Surface::spanwise_strip_count` exposes the result.
- `PegasusError::Parse` carries a structured `ParseLocation { path, line }`
  instead of a pre-formatted string, so callers can recover the line number and
  file programmatically. `PegasusError` is now `#[non_exhaustive]`.
- 400 of the 488 frozen C++ cases now have their real assertions reproduced in
  Rust, up from none, with 4 recorded as unportable because they cover
  behaviour that was deliberately removed. Twenty-two C++ files are complete:
  every non-gradient file, plus `analytical_biot_savart_test` (10),
  `analytical_aic_test` (9) and `analytical_trefftz_test` (8). The remaining 84
  are the rest of the gradient suites. `parity_cases.rs` tracks ported,
  unportable and backlog counts separately and asserts all three, so they
  cannot be misreported.
- `tests/parity_cases.rs` no longer generates one test per frozen CTest case.
  The 488 generated tests all dispatched to one of eight category-level smoke
  checks, so every case within a category ran byte-identical assertions
  regardless of the C++ behaviour it was named after — 488 reported passes
  backed by eight assertions, which also inflated the line-coverage gate. The
  file now asserts the manifest's integrity and records the porting gap
  explicitly; real behavioural coverage lives in the dedicated suites.
- Meshes large enough to reach the multithreaded assembly are now covered by
  `tests/mirror_reuse.rs`. Both rayon drivers gate at 128 panels and no test
  mesh had previously exceeded that, so the parallel code was compiled but
  never executed and CI's `RAYON_NUM_THREADS` sweep exercised nothing.
- The package is now `pegasus-vlm`, imported as `pegasus_vlm`, and requires
  Rust 1.93.1.
- The public API uses owned Rust data, checked slices and indices,
  `PegasusError`, grouped `Coefficients`, and `AeroOutput`.
- A single partial-pivot AIC factorization is cached and reused for primal,
  transposed-adjoint, unit-response, and repeated-flight-condition solves,
  avoiding duplicate cubic work while preserving single-threaded execution.
- The `analysis_scaling` example now mirrors the v0.2.0 C++ `benchmark_analysis`
  harness iteration for iteration: one untimed warm-up, 50 timed repetitions,
  an explicit `rebuild_aic`, and the force calculation inside the timed region.
  It previously ran 10 repetitions with no warm-up and omitted the forces,
  which made it measure less work than the C++ harness it is compared against.
- Influence assembly now iterates source panels on the outside and target
  panels on the inside, matching the column-major storage order, and fuses the
  normalwash contraction into the same pass. This replaces a strided traversal
  plus a second full pass over the velocity matrix, and hoists the source panel
  width, core radii, and symmetry-image geometry out of the inner loop. Output
  is bit-for-bit unchanged.
- The influence kernel sums its two trailing-leg contributions with each other
  before adding them to the bound contribution, instead of accumulating them
  one leg at a time. This is a reassociation, not a correction: neither order is
  more accurate. It makes the kernel equivariant under a y-reflection in
  floating point rather than only in exact arithmetic, which is what lets
  mirrored influence entries be derived instead of evaluated.
  **This shifts results.** Influence coefficients move by about one unit in the
  last place; circulations move by that amplified by the AIC condition number,
  measured as 1.5e-16 relative on an 8-panel wing rising to 4.2e-12 on 384
  panels. All parity fixtures pass unchanged, the strictest with roughly three
  orders of magnitude of tolerance to spare.
- Rejecting a singular system no longer costs a separate full pass over the
  assembled normalwash matrix. Finiteness is accumulated while the matrix is
  written and reported through `normalwash_and_velocity_checked`.
- Unit-response solves no longer form the `A x - b` residual. The product was
  used only to test for non-finite values, which `rebuild_aic` and the solver
  already reject, so it cost a full `N²×6` product and an `N×6` allocation
  without establishing any additional property.
- The release profile enables fat LTO and a single codegen unit. Cargo ignores
  library profiles for downstream consumers, so this affects only this
  workspace's examples, benchmarks, and tests.
- Flight-condition, aircraft-section, and direct mesh-node gradients use
  analytical influence/force kernels and adjoint contractions; the Mach
  derivative no longer requires perturbed full solves.
- The missing fixed-circulation alpha term in the adjoint force derivative now
  includes freestream and body-to-stability-axis effects. The former disabled
  CL validation is enabled and extended across all six primary coefficients.
- Section drag polars are now applied to the loads. v0.2.0 parsed `CDCL`,
  interpolated the six coefficients onto each strip and set `has_viscous_data`,
  but no force path ever read them back — every reference to `viscous_polar` in
  the C++ tree is a write. `strip_loads` folds the polar-interpolated term into
  the strip body force, so drag, roll, pitch and yaw now respond to a `CDCL`
  block that the oracle silently ignored. Lift and side are unchanged, since
  the added vector is parallel to the stability drag axis.

### Removed

- The C++ API, C ABI expectations, CMake build, system Eigen/VTK/OpenMP
  discovery, GoogleTest, and C++ examples.
- Philote servers, `MeshDiscipline`, `VlmDiscipline`, and Philote-dependent
  OpenMDAO workflows. These integrations are deferred to a later release.
- `GradientCalculator::set_finite_difference_step`, a deprecated no-op, along
  with the internal finite-difference machinery it documented. No finite
  difference remains anywhere in `src/`.
- `VortexLattice::{vortex_core_radius_chord_fraction,
  set_vortex_core_radius_chord_fraction, vortex_core_radius_width_fraction,
  set_vortex_core_radius_width_fraction}`. These were a second, independent
  pair of core-radius fractions with different defaults (0.75/0.25 against the
  builder's 0.0001/0.0001) that were never passed to assembly, so setting them
  left the AIC bit-identical while appearing to tune the regularization. The
  live knobs are `InfluenceMatrixBuilder::{set_core_radius_chord_fraction,
  set_core_radius_width_fraction}`, which apply between panels on different
  components.
- The `Directional` forward-mode scalar in `gradients/trefftz.rs`, a second
  dual number under another name, along with its five operator impls.
  `trefftz_geometry_jvp` now carries the value and its directional rate as
  plain paired scalars through hand-written closed forms. With this and the
  Biot-Savart change, no operator-overloaded differentiation scalar of any kind
  remains in `src/`.
- `ApparentMass`, `EigenmodeAnalysis::set_apparent_mass`, and the added-mass
  terms in the state matrices. v0.2.0 exposed `SetApparentMass` and stored the
  value, but `eigenmode_analysis.cpp` never read it back, so the setter was
  inert and the struct was carried purely to be assigned to. The Rust had
  revived it with an invented convention — effective mass `m + tr(M) / 3` and a
  matching inertia diagonal — that no oracle behaviour corresponds to. The
  state matrices now use the real mass and inertia directly. The two frozen
  C++ cases covering this, `ApparentMassTest.DefaultState` and
  `EigenmodeAnalysisTest.SetApparentMass`, are correspondingly unportable: the
  first asserted the struct's zero defaults and the second only that the setter
  did not crash.
- The internal forward-mode dual-number scalar that the Biot-Savart VJPs were
  built on. `v0.2.0` had already deleted `gradients/dual_number.h` and
  `gradients/autodiff_tape.h` in favour of hand-derived closed forms, so
  carrying a dual through the kernel put back the machinery the behavioural
  oracle deliberately removed. The three `horseshoe_velocity_vjp_*` functions
  are now closed-form expressions matching the oracle's
  `analytical_biot_savart.cpp`, and no automatic differentiation of any kind
  remains in `src/`.

### Fixed

- `ControlDerivatives::default()` reports a unit gain again. v0.2.0 declared
  `Scalar gain{1.0}` and asserted it; the Rust used a derived `Default`, which
  gave zero. Gain is a multiplier, so the neutral element is one — a
  default-constructed value scaled every derivative it touched to nothing.
  `StabilityDerivativeCalculator` always set the field explicitly, so no
  computed result was affected; the defect was reachable only through the
  public struct literal.
- The Trefftz circulation gradients now skip exactly the strips the primal
  far-field loop skips. `ForceCalculator` drops a strip when its trailing-edge
  width is `< EPSILON` and `trefftz_circulation_gradients` used `<=`, so a strip
  of width exactly `EPSILON` was summed into the primal while its circulation
  derivative was left at zero.
- The Biot-Savart VJPs are now invariant under a uniform rescaling of the
  geometry, matching the primal. `horseshoe_vortex_velocity` tests collinearity
  with the dimensionless `cross_sq > EPSILON * a_sq * b_sq`, but the gradient
  kernel compared the length⁴ quantity `cross_sq` against an absolute
  `EPSILON` — the same class of bug that was fixed in influence assembly, but
  left in place on the derivative side. Below roughly millimetre scale the
  VJP therefore dropped the bound-leg term that the primal still evaluated: at
  1e-4 scale the derivative disagreed with a centered difference by the entire
  bound contribution, including a sign flip on the bound-right x-component,
  with no error reported. A sweep of 12960 derivative components across four
  decades of geometry scale now agrees with centered differences to 8.6e-7,
  the truncation floor of the difference itself.
- Flight-condition gradients now include the section drag polar. `strip_loads`
  folds a polar-interpolated viscous term into the strip's body force, which
  flows into the totals and, through `arm x force`, into the moments — but the
  gradient module contained no reference to it at all, so loading any AVL file
  with a `CDCL` block produced derivatives that were wrong by the polar slope on
  drag, roll, pitch and yaw. Lift and side were unaffected, because the added
  vector is parallel to the stability drag axis and cancels exactly; that
  asymmetry is now asserted. The chain is non-local: a strip's viscous drag
  depends on its lift coefficient and hence, through the induced velocity, on
  every circulation in the mesh.
- `interpolate_drag_polar` and its slope are produced by one function,
  `drag_polar_value_and_slope`, so the gradient kernels cannot drift from the
  primal. The `max(0.0)` clamp and the coincident-breakpoint fallback both
  report a zero slope, matching the flat primal there.
- Degenerate section geometry is rejected at solver entry instead of meshing
  into NaN. Zero or negative chord, coincident spanwise stations, sections that
  reverse spanwise direction, and non-finite coordinates all previously
  generated a lattice "successfully": repeating a station — a standard AVL
  spacing-break idiom — produced 2 non-finite control points out of 18, which
  then poisoned every entry of the AIC. The station-separation test is relative
  to the longest interval, so it means the same thing at any model scale.
- `Aircraft::total_vortex_count` and `total_strip_count` now agree with the
  lattice the mesher builds. They multiplied the surface-level panel count by
  the duplicate factor, ignoring per-section overrides and the per-interval
  rounding, so they reported 20 strips where the mesher produced 19, and 20
  where a section-level override produced 7.
- Influence assembly is now invariant under a uniform rescaling of the
  geometry. The bound-leg guard compared a length⁴ quantity, and the trailing
  legs a length, against an absolute `EPSILON`, so below roughly millimetre
  scale the bound contribution was silently dropped from every panel: `‖AIC‖`
  moved 3.3% at 1e-3 scale and 4.2% at 1e-4, with no error reported. The tests
  are the dimensionless ratios that were intended, and agree with the previous
  thresholds at unit scale.
- `solve` now rejects non-finite body rates and a NaN dynamic pressure. Only
  `alpha`, `beta`, `mach` and the control deflections were validated, so a NaN
  roll rate returned an all-NaN circulation flagged as a valid solution; NaN
  also defeated the `dynamic_pressure < 0.0` test, since every ordered
  comparison against NaN is false.
- `solve` now rejects Mach numbers at or above 1 instead of clamping the
  Prandtl-Glauert factor to 1 and returning the incompressible answer, which
  produced finite, plausible-looking results for a regime the method cannot
  represent while `mach_number()` still reported the supersonic value.
- `InfluenceMatrixBuilder::induced_velocity` now includes the doubly-reflected
  yz image. Assembly adds it whenever both symmetry planes are active, so the
  sampled velocity field disagreed with the influence matrix by roughly 3%,
  including a sign flip in the y component.
- A failed `rebuild_aic` no longer publishes non-finite matrices. The fields
  were written before the finiteness check, so `aic()`, `velocity_matrix()`,
  `linearized_aic()` and `solution_aic()` handed out NaN-laden matrices paired
  with a stale factorization.
- `SurfaceMeshInfo::original_surface_index` is now always a `surface_info`
  slot. The aircraft generator stored an `Aircraft::surfaces()` index and the
  direct-mesh generator left it at the `Default` of 0, so for
  `[wing(duplicated), tail]` the tail's value pointed at the wing's mirror.
- Reflection reuse now verifies that its mirror blocks tile the column space
  contiguously and in order. That held only because both generators happen to
  append each duplicate immediately after its original, and the check was a
  `debug_assert` that vanished in release, where a violation would have written
  columns into the wrong block.
- The AVL reader no longer panics on out-of-range panel counts. `value as usize`
  saturates, so a paneling field of `1e30` became `usize::MAX` and overflowed
  the mesher's panel arithmetic; counts are now range-checked at parse time and
  rejected with a typed error.
- `CubicSpline::evaluate` no longer panics on a NaN parameter. NaN satisfies
  neither interval comparison, so the interval search returned 0 and the
  caller's `i - 1` underflowed. NaN now propagates through as a NaN result, as
  the legacy binary search did. Reachable from `Airfoil::camber`.
- The mass reader no longer discards a line because of its comment. Comments
  are stripped before the line is classified, so an `=` inside a comment such
  as `! main gear (x = 1.2 m)` no longer captures the line as a unit
  assignment and silently drops the component's mass, centre of gravity and
  inertia.
- Mass `*` scale and `+` offset lines no longer read numbers out of their own
  comments, and leave the previous value untouched when fewer than three
  numeric fields are present.
- Mass component fields stop at the first token that does not parse, matching
  the legacy stream extraction, instead of skipping it and promoting later
  numbers into earlier slots.
- A `CONTROL` line without the optional trailing `SgnDup` field is accepted
  with a default sign of 1.0. Requiring it rejected the entire aircraft over an
  omission the legacy reader tolerated.
- A failed `AFILE` load is now a parse error carrying the line number instead
  of being silently swallowed, which left the section on its default flat plate
  and let the solve report success against the wrong camber. The filename is
  comment-stripped.

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
