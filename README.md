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
behavioral oracle for parity fixtures. The optional `philote` feature exposes
mesh and VLM disciplines for external MDO workflows. The crate has no C, C++,
or Python runtime dependency.

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

## Controls in external optimizations and coupled simulations

A `SurfaceMesh` can carry `MeshControlSurface` definitions. Pass the structurally
deformed mesh **before control actuation** and put each command's deflection in
`FlightCondition::control_deflections`, in radians. Pegasus rotates the flow
boundary-condition normals; it does not rotate mesh nodes or the wake. If the
caller already physically deflected the mesh, omit the corresponding control
actuation here.

```rust
use pegasus_vlm::geometry::MeshControlSurface;

// For a 4 x 6 panel mesh: strips 1 through 5, hinge at chordwise node row 3.
let flap = MeshControlSurface::trailing_edge("flap", 1..6, vec![3.0; 7]);
// wing.controls.push(flap);
```

Hinge coordinates are **material mesh coordinates**, not chord fractions:
`0` is the leading edge, `nx` the trailing edge, and `i + f` interpolates between
node rows `i` and `i + 1`. Provide `ny + 1` coordinates, including inactive
stations. The hinge follows mesh bending and twist. Control coverage uses the
mean panel overlap at the two span stations; mesh derivatives hold those
material coordinates and the active span interval fixed. Align panel edges with
the hinge when practical. Hinge placement and span-extent design derivatives are
not exposed.

`gain` specifies rotation per command using the right-hand rule about the hinge
from station `j` to `j + 1`. `duplicate_sign` is +1 for symmetric deflection and
-1 for antisymmetric deflection on a mirrored mesh. Matching names on several
surfaces share a command. Separate meshes and names give independent left/right
commands; multiple definitions can also mix flap/aileron commands. Set
`leading_edge = true` for panels ahead of a hinge. Commands compose in the order
of the flight-condition list, then definition order on each panel. Keep that
order fixed during optimization. Omitted commands mean zero; unknown names are
errors. For legacy repeated names, a named derivative varies all matching
entries together; the derivative of an omitted command appends it at zero.

- `GradientCalculator::compute_control_gradients` returns derivatives per radian
  for all named commands, honoring `OutputSelection`, including near-field and
  Trefftz induced drag.
- `compute_mesh_node_gradients` and section-geometry gradients support nonzero
  control deflections. Mesh gradients include the changing hinge direction.
- `StripLoadAdjoint::control_cotangents` differentiates arbitrary functionals of
  distributed strip forces and moments. Its `mesh_node_cotangents` and
  `flight_condition_cotangents` also work with active controls. Load seeds use
  the existing raw Kutta–Joukowski resultant convention: multiply seeds by
  dynamic pressure when differentiating dimensional loads. Moment points are
  held fixed; callers must differentiate their motion separately.

The external application handles structural equilibrium, load transfer, and
total coupled derivatives. It can supply trim commands or use the optional
internal trim solve below. Geometry/strip
adjoints still require explicit mesh duplication rather than symmetry planes,
and do not support section drag polars. This remains a quasi-steady VLM normal
rotation model, not a separated-flow or unsteady control model.

Run the complete native example with `cargo run --example control_surfaces`.

With the `philote` feature, the VLM discipline accepts a `controls` option as a
JSON array (or a JSON-encoded string over the wire):

```json
[
  {"surface": "Wing", "name": "flap", "span_start": 1, "span_end": 6,
   "hinge": 3.0, "gain": 1.0, "duplicate_sign": 1.0}
]
```

`hinge` accepts a constant or a `ny + 1` array. Optional `leading_edge` defaults
to false. Set surface names and dimensions in the same options request. Each
unique name declares a scalar input `control_<name>` in radians, with analytical
partials for the selected coefficient outputs. Philote orders commands
alphabetically. Supply these inputs on every
evaluation, including zero values. They participate in solution caching and can
be included in `fixed_inputs`; fixing a command suppresses its partials while
retaining its deflection and its effect on mesh derivatives. Native callers can
also populate `SurfaceOption::controls` directly.

## Optional trim in the explicit Philote discipline

`VlmDiscipline` remains an explicit component. Configure its `trim` option
before setup to solve selected `alpha`, `beta`, and control commands internally.
Each solved variable becomes an output, while unselected operating variables
remain inputs. All declared function gradients include the converged trim
response, including derivatives with respect to deformed mesh nodes, targets,
external loads, dynamic pressure, and prescribed commands. Trim is off by default.

```json
{"trim": {
  "enabled": true,
  "variables": ["alpha", "control_elevator"],
  "balances": [{"quantity": "Fz", "scale": 1000.0},
               {"quantity": "My", "scale": 1000.0}],
  "initial_guess": {"alpha": 0.03, "control_elevator": 0.0},
  "bounds": {"alpha": [-0.3, 0.3], "control_elevator": [-0.4, 0.4]}
}}
```

For dimensional balances, supply `external_Fx`, `external_Fy`, `external_Fz`,
`external_Mx`, `external_My`, and `external_Mz`, plus the declared balance target
inputs (by default `trim_target_Fz` and `trim_target_My` in this example).
The residual is aerodynamic load plus prescribed external load minus target.
This example enforces vertical-force and pitching-moment equilibrium only.
Coefficient targets such as `CL` and `Cm` also work and require no external loads.

`outputs` can select `CL,CD,CY,Cl,Cm,Cn` and dimensional `Fx,Fy,Fz,Mx,My,Mz`.
Forces use the mesh/body coordinate axes: x has positive drag and z has positive
lift at zero angle of attack. Dimensional moments follow the right-hand rule
about the configured reference point; Mx/Mz remove the conventional Cl/Cn sign
flips. External loads must use the same axes and reference point.

The core `MeshTrimSolver` is available without Philote. Run
`cargo run --example mesh_trim` for a deformed wing/tail trim and analytical
sensitivities. For the explicit wire workflow, start
`cargo run --features philote --example philote_vlm_server`, then run
`cargo run --features philote --example philote_trim_client` in another terminal.
The [MDO integration guide](docs/docs/guide/mdo-integration.md) describes control
mixers, numerical options, gradient semantics, and external coupling.

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
