# Zephyr VLM

Zephyr is a safe, Rust-only vortex-lattice-method library for aircraft
geometry, aerodynamic analysis, stability and control derivatives, trim,
eigenmodes, analytical/adjoint gradients, AVL input, and VTK output.

The package is named `zephyr-vlm` and imported as `zephyr_vlm`. The
unreleased Rust API is a breaking replacement for the former C++ API. The
repository is publication-ready, but this port does not publish it to
crates.io.

## Status

The Rust port preserves the numerical conventions and panel ordering of
Pegasus 0.2.0. The legacy source at tag `v0.2.0` (commit `556b0ea`) is the
behavioral oracle for parity fixtures. Philote servers and Philote-dependent
OpenMDAO workflows are intentionally deferred; the crate has no C, C++, or
Python runtime dependency.

## Units and axes

Zephyr uses double precision throughout. Geometry is expressed in a
right-handed body frame with x aft, y starboard, and z up. The following unit
distinction is deliberate:

- Section twist and surface incidence are in **degrees**.
- Flight angles, angular rates, control deflections, trim variables, and
  derivatives are in **radians**.
- Length, area, mass, density, and dynamic pressure use a consistent unit
  system selected by the caller.
- Aerodynamic coefficients and normalized body rates are dimensionless.

## Add the crate

```toml
[dependencies]
zephyr-vlm = { git = "https://github.com/chrislupp/Zephyr.git" }
```

The default `vtk` feature uses the pure-Rust `vtkio` crate. A minimal build is
available with `default-features = false`; visualization functions then return
`ZephyrError::FeatureDisabled`.

## Basic analysis

```rust
use zephyr_vlm::aero::ForceCalculator;
use zephyr_vlm::core::{DEGREES_TO_RADIANS, Vec3};
use zephyr_vlm::geometry::{Aircraft, Section, Surface};
use zephyr_vlm::solver::{FlightCondition, VortexLatticeSolver};

let mut wing = Surface::new("Wing");
wing.set_chordwise_panel_count(4);
wing.set_spanwise_panel_count(16);
wing.add_section(Section::new(Vec3::new(0.0, -5.0, 0.0), 1.0, 0.0));
wing.add_section(Section::new(Vec3::new(0.0,  5.0, 0.0), 1.0, 0.0));

let mut aircraft = Aircraft::new("Rectangular wing");
aircraft.set_reference_area(10.0);
aircraft.set_reference_chord(1.0);
aircraft.set_reference_span(10.0);
aircraft.set_moment_reference_point(Vec3::new(0.25, 0.0, 0.0));
aircraft.add_surface(wing);

let condition = FlightCondition {
    alpha: 5.0 * DEGREES_TO_RADIANS,
    ..FlightCondition::default()
};
let mut solver = VortexLatticeSolver::new();
solver.set_aircraft(&aircraft)?;
let solution = solver.solve(&condition)?;
let results = ForceCalculator::new().compute_for_aircraft(
    &aircraft,
    solver.mesh(),
    &solution,
    &condition,
);

assert!(results.is_valid);
assert!(results.coefficients.lift > 0.0);
# Ok::<(), zephyr_vlm::ZephyrError>(())
```

## AVL and mass files

`AvlReader::read_path` and `AvlReader::read_str` load AVL geometry. Relative
`AFILE` paths are resolved from the geometry file. `MassReader` provides the
corresponding path and string APIs for AVL mass files. Parse errors retain the
source name, line, and active keyword.

```rust,no_run
use zephyr_vlm::io::{AvlReader, MassReader};

let aircraft = AvlReader::new().read_path("aircraft.avl")?;
let mass = MassReader::read_path("aircraft.mass")?;
# let _ = (aircraft, mass);
# Ok::<(), zephyr_vlm::ZephyrError>(())
```

## Analytical gradients

The adjoint calculation includes both the implicit circulation response and
the explicit fixed-circulation force response. Flight angles, rates, and Mach
use analytical derivatives; section and mesh-node gradients reuse one primal
solve and cached batched transpose solves. Objectives use `AeroOutput` instead
of integer selectors.

```rust,no_run
use zephyr_vlm::gradients::{AeroOutput, TotalDerivative};
use zephyr_vlm::solver::{FlightCondition, VortexLatticeSolver};
# use zephyr_vlm::geometry::Aircraft;
# fn evaluate(aircraft: &Aircraft) -> zephyr_vlm::Result<()> {

let condition = FlightCondition::default();
let mut solver = VortexLatticeSolver::new();
solver.set_aircraft(aircraft)?;
solver.solve(&condition)?;
let lift = TotalDerivative::compute(&solver, &condition, AeroOutput::Lift)?;
println!("dCL/dalpha = {}", lift.by_alpha);
# Ok(())
# }
```

## Selective outputs and half-wing analysis

`analysis::Analysis` accepts explicit output and derivative requests. Unrequested
force calculations, Trefftz interactions, detailed arrays, and derivative stages
are skipped. An optional `SolverMode::HalfWing` solves only the original half's
unknowns while accounting for reflected influences. Full geometry is the default;
existing APIs keep their current behavior.

See [the selective analysis guide](docs/selective-analysis.md) for the request API,
symmetry restrictions, reference conventions, caching, and Philote options. Run
`cargo run --release --example selective_analysis` for a working example.

## VTK visualization

The visualization module writes surface meshes (`.vtu`) and spanwise loads,
wakes, and velocity fields (`.vtp`) in XML ASCII or binary form. Each export is
written to a temporary sibling and renamed only after serialization succeeds,
so failed exports do not leave a partial result.

## Development

Rust 1.93.1 is pinned in `rust-toolchain.toml`.

```text
cargo test --all-features
cargo test --no-default-features
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
RUSTDOCFLAGS="-D warnings" cargo doc --all-features --no-deps
```

The repository also contains native Rust examples, Criterion scaling
benchmarks, committed parity fixtures, and a manifest mapping the non-Philote
0.2.0 CTest cases to Rust replacements.

## License

Zephyr is licensed under GPL-3.0-or-later. See the repository's `LICENSE`
file.

## Public release

This work has been cleared for public release, distribution unlimited, case
number: AFRL-2026-2198.

The views expressed are those of the authors and do not reflect the official
guidance or position of the United States Government, the Department of
Defense or of the United States Air Force.

Statement from DoD: The Appearance of external hyperlinks does not constitute
endorsement by the United States Department of Defense (DoD) of the linked
websites, of the information, products, or services contained therein. The
DoD does not exercise any editorial, security, or other control over the
information you may find at these locations.
