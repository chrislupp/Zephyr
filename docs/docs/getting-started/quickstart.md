---
sidebar_position: 2
---

# Quick Start

An analysis has four parts: geometry, a flight condition, a solve, and a force
integration. This page builds the smallest useful case — a rectangular wing at a
fixed angle of attack — and then shows what to reuse when you sweep it.

## A wing and its loads

```rust
use zephyr_vlm::aero::ForceCalculator;
use zephyr_vlm::core::{DEGREES_TO_RADIANS, Vec3};
use zephyr_vlm::geometry::{Aircraft, Section, Surface};
use zephyr_vlm::solver::{FlightCondition, VortexLatticeSolver};

let mut wing = Surface::new("Wing");
wing.set_chordwise_panel_count(4);
wing.set_spanwise_panel_count(16);
wing.add_section(Section::new(Vec3::new(0.0, -5.0, 0.0), 1.0, 0.0));
wing.add_section(Section::new(Vec3::new(0.0, 5.0, 0.0), 1.0, 0.0));

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

println!("CL  = {:.6}", results.coefficients.lift);
println!("CDi = {:.6}", results.far_field_induced_drag);
println!("Cm  = {:.6}", results.coefficients.pitch);
println!("e   = {:.4}", results.span_efficiency);
```

Figure 1 on the contents page is this wing's lattice, drawn from the mesh the
solver builds for it.

## What each step did

`Surface::new` makes an empty lifting surface. Two `Section`s give it a
planform: each carries a leading-edge point, a chord, and a twist in degrees.
The panel counts set the lattice: 4 chordwise by 16 spanwise, cosine spaced in
both directions by default, so 64 horseshoe vortices and 16 spanwise strips.

`Aircraft` owns the surfaces and the quantities that normalize the answer:
reference area, chord, and span, and the point moments are taken about. Nothing
infers them from the planform — a wrong reference area silently rescales every
coefficient, and it is the most common first mistake.

`FlightCondition::default()` is level flight at zero alpha, zero sideslip, zero
rates, zero Mach, and unit dynamic pressure. Spreading one field over the
default, as above, is the usual way to build one.

`solver.set_aircraft` generates the lattice and invalidates any cached state.
`solve` assembles the influence matrices if they are stale, factors the
normalwash matrix, and returns a `SolutionState` holding the circulation.

`ForceCalculator` turns that circulation into loads. `compute_for_aircraft`
takes the reference quantities and the constant profile-drag coefficient from
the `Aircraft`; `compute_with_references` is the entry point for direct meshes.
The result carries near-field forces and moments, Trefftz-plane lift, side
force and induced drag, span efficiency, and per-surface and per-strip records.

## Units

Two conventions coexist, deliberately:

- **Degrees** for section twist and surface incidence, following AVL.
- **Radians** for flight angles, rates, control deflections, trim variables, and
  every derivative.

Lengths, areas, masses, densities, and dynamic pressures are in whatever
consistent system you choose; the library never converts. See
[Conventions](./conventions.md) for axes, signs, and the coefficient
definitions.

## Sweeping a condition

A solver holds its lattice, influence matrices, and factorization. Changing only
alpha, beta, or the body rates reuses all three, so a sweep should keep one
solver alive rather than building a new one per point.

```rust
for alpha_degrees in [-2.0_f64, 0.0, 2.0, 4.0, 6.0, 8.0] {
    let condition = FlightCondition {
        alpha: alpha_degrees * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let solution = solver.solve(&condition)?;
    let results = ForceCalculator::new().compute_for_aircraft(
        aircraft,
        solver.mesh(),
        &solution,
        &condition,
    );
    println!("{alpha_degrees:4.1}  {:.6}", results.coefficients.lift);
}
```

That is [Basic Analysis](../examples/basic-analysis.md), runnable with
`cargo run --example basic_analysis`.

## Errors

Fallible calls return `Result<T, ZephyrError>`. Geometry and mesh problems are
reported before any matrix is assembled; a stale or invalid `SolutionState` is
rejected by anything that would otherwise consume it silently. Parse failures
from the AVL readers keep the file path, the line number, and the keyword being
read.

## Where to go next

- [Conventions](./conventions.md) — axes, signs, units, and panel ordering.
- [Geometry](../guide/geometry.md) — sections, spacing, duplication, airfoils.
- [Loads and Coefficients](../guide/loads.md) — what the near and far field each
  report, and when they disagree.
- [Analytical and Adjoint Gradients](../guide/automatic-differentiation.md) —
  derivatives without finite differences.
- [AVL Input Files](../guide/avl-files.md) — reading an existing model instead of
  writing one.
