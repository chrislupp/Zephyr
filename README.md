# Pegasus VLM

Pegasus is a safe, Rust-only vortex-lattice-method library for aircraft
geometry, aerodynamic analysis, stability and control derivatives, trim,
eigenmodes, analytical/adjoint gradients, AVL input, and VTK output.

The package is named `pegasus-vlm` and imported as `pegasus_vlm`. The
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

Pegasus uses double precision throughout. Geometry is expressed in a
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
pegasus-vlm = { git = "https://github.com/chrislupp/Pegasus.git" }
```

The default `vtk` feature uses the pure-Rust `vtkio` crate. A minimal build is
available with `default-features = false`; visualization functions then return
`PegasusError::FeatureDisabled`.

## Basic analysis

```rust
use pegasus_vlm::aero::ForceCalculator;
use pegasus_vlm::core::{DEGREES_TO_RADIANS, Vec3};
use pegasus_vlm::geometry::{Aircraft, Section, Surface};
use pegasus_vlm::solver::{FlightCondition, VortexLatticeSolver};

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
# Ok::<(), pegasus_vlm::PegasusError>(())
```

## AVL and mass files

`AvlReader::read_path` and `AvlReader::read_str` load AVL geometry. Relative
`AFILE` paths are resolved from the geometry file. `MassReader` provides the
corresponding path and string APIs for AVL mass files. Parse errors retain the
source name, line, and active keyword.

```rust,no_run
use pegasus_vlm::io::{AvlReader, MassReader};

let aircraft = AvlReader::new().read_path("aircraft.avl")?;
let mass = MassReader::read_path("aircraft.mass")?;
# let _ = (aircraft, mass);
# Ok::<(), pegasus_vlm::PegasusError>(())
```

## Analytical gradients

The adjoint calculation includes both the implicit circulation response and
the explicit fixed-circulation force response. Flight angles, rates, and Mach
use analytical derivatives; section and mesh-node gradients reuse one primal
solve and cached batched transpose solves. Objectives use `AeroOutput` instead
of integer selectors.

```rust,no_run
use pegasus_vlm::gradients::{AeroOutput, TotalDerivative};
use pegasus_vlm::solver::{FlightCondition, VortexLatticeSolver};
# use pegasus_vlm::geometry::Aircraft;
# fn evaluate(aircraft: &Aircraft) -> pegasus_vlm::Result<()> {

let condition = FlightCondition::default();
let mut solver = VortexLatticeSolver::new();
solver.set_aircraft(aircraft)?;
solver.solve(&condition)?;
let lift = TotalDerivative::compute(&solver, &condition, AeroOutput::Lift)?;
println!("dCL/dalpha = {}", lift.by_alpha);
# Ok(())
# }
```

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

Pegasus is licensed under GPL-3.0-or-later. See the repository's `LICENSE`
file.

### Direct-mesh derivative products

After solving, `gradients::MeshLinearization::new(&solver)` provides analytical
`jvp(&directions, selection)` and `vjp(&coefficient_seeds)` methods. Directions
and returned cotangents use the original `MeshInput` surface/node order. A node
index is `i * (ny + 1) + j`, with XYZ in the mesh's physical coordinate units.
VJP seeds are `AeroGradients` values interpreted as weights on lift, near-field
drag, side force, roll, pitch, yaw and Trefftz induced drag. Directions are
`MeshNodeDirection` values; each supplied direction uses one tangent right-hand
side. Each coefficient seed uses one adjoint right-hand side. Both paths batch
right-hand sides and reuse the primal factorization.

For `R = A(x) Gamma - b(x)`, the tangent solves
`A dGamma = db - dA Gamma` and returns `f_x dx + f_Gamma dGamma`.
The reverse product solves `A^T lambda = f_Gamma^T w` and accumulates
`f_x^T w - R_x^T lambda`. Shared nodes and mirrored surfaces accumulate through
the transpose of lattice construction. Coefficient reference quantities and
flight condition are fixed; all moment coefficients retain Pegasus's reference
chord convention. The linearization borrows the solver, original input mesh
and the condition stored with its valid solution, so they cannot be mutated
while applying products.

The reverse Trefftz geometry kernel differentiates the same endpoint geometry
and cross-surface core law as the primal. Its pair term is
`weight * dot(strip_width, displacement) / radius`, with
`radius = sqrt(dot(displacement, displacement)^2 + core^4)`. It accumulates
endpoint, strip-width and core cotangents in one strip-pair sweep, then maps the
active core's chord/width derivative back to nodes. The core maximum uses the
chord branch at a tie (a nonsmooth, branch-specific derivative). The influence
kernel also includes its core-radius dependence. Point-coincidence and strip
collapse branches follow the primal; collapsed chord stations, section polars,
deflected control normals and symmetry planes are rejected by this interface.

Use direct products for few input directions, adjoint products for few weighted
responses, and choose based on the actual system. The coefficient node-gradient
convenience API now uses this reverse accumulation instead of reconstructing
one coordinate direction at a time. No AD, finite differences, differentiation
code generation or optimizer is used in these production paths.

For many strip-load outputs and few input directions,
`MeshLinearization::strip_load_jvp` returns force/moment directions from tangent
solves, including optional motion of strip moment reference points.
`StripLoadAdjoint::mesh_node_cotangents` and `moment_point_cotangents` supply the
corresponding reverse products. Raw resultants use the existing Kutta-Joukowski
normalization; the caller supplies dynamic-pressure scaling. With
`F = 2 Gamma (V cross line)` and `M = (cp-p) cross F`, the explicit point term is
`dM = -dp cross F` and its transpose is `bar_p = bar_M cross F`. These local
operators do not replace a coupled aeroelastic tangent/adjoint solve.
