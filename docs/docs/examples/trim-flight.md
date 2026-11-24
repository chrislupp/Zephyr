---
sidebar_position: 2
---

# Trim Flight

A wing–tail aircraft trimmed to a lift target by solving for angle of attack.

```bash
cargo run --example trim_flight
```

Source: `examples/trim_flight.rs`.

## What it does

```rust
use zephyr_vlm::core::Vec3;
use zephyr_vlm::geometry::{Aircraft, Section, Surface};
use zephyr_vlm::solver::VortexLatticeSolver;
use zephyr_vlm::trim::{TrimPresets, TrimSolver};

fn main() -> zephyr_vlm::Result<()> {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(5);
    wing.set_spanwise_panel_count(16);
    wing.add_section(Section::new(Vec3::new(0.0, -5.0, 0.0), 1.5, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 5.0, 0.0), 1.5, 0.0));

    let mut tail = Surface::new("Horizontal tail");
    tail.set_chordwise_panel_count(3);
    tail.set_spanwise_panel_count(8);
    tail.add_section(Section::new(Vec3::new(4.0, -1.5, 0.0), 0.6, 0.0));
    tail.add_section(Section::new(Vec3::new(4.0, 1.5, 0.0), 0.4, 0.0));

    let mut aircraft = Aircraft::new("Wing-tail trim example");
    aircraft.set_reference_area(15.0);
    aircraft.set_reference_chord(1.5);
    aircraft.set_reference_span(10.0);
    aircraft.set_moment_reference_point(Vec3::new(0.375, 0.0, 0.0));
    aircraft.add_surface(wing);
    aircraft.add_surface(tail);

    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft)?;

    let result = TrimSolver::new().solve(
        &aircraft,
        &mut solver,
        &TrimPresets::level_flight(0.5),
        None,
    );

    println!("status: {:?}", result.status);
    println!("alpha: {:.4} deg", result.state.alpha.to_degrees());
    println!("CL: {:.6}", result.aerodynamic_results.coefficients.lift);
    println!("residual: {:.3e}", result.residual_norm);
    Ok(())
}
```

`TrimPresets::level_flight(0.5)` is a one-equation problem: solve alpha for
$C_L = 0.5$, with sideslip and all three rates held at zero. The tail is present
aerodynamically — it carries load and downwashes the wing — but nothing is being
solved for it.

The last argument is an optional initial guess; `None` starts from zero.

## Checking the result

Always read `status` before the numbers. Nonconvergence here is a
`TrimStatus`, not an error, and `TrimResult` carries its last state either way —
so a failed trim that is read as a converged one will hand you a plausible alpha
with a large residual behind it.

```rust
if result.status != TrimStatus::Converged {
    eprintln!("trim failed: {} ({:?})", result.diagnostic, result.status);
}
```

`TrimResult` also carries the full `aerodynamic_results` and the
`StabilityDerivatives` at the trimmed point, so a trim is usually the first step
of an analysis rather than its own output.

## Adding the elevator

To actually trim the aircraft — lift target *and* zero pitching moment — add a
control to the tail sections and use the two-equation preset:

```rust
let mut inboard = Section::new(Vec3::new(4.0, -1.5, 0.0), 0.6, 0.0);
let mut outboard = Section::new(Vec3::new(4.0, 1.5, 0.0), 0.4, 0.0);

let elevator = ControlSurfaceDefinition {
    name: "elevator".into(),
    hinge_x_fraction: 0.7,
    hinge_vector: Vec3::zeros(),
    sign: 1.0,
    gain: 1.0,
};
// Both bounding sections need the definition.
inboard.add_control_surface(elevator.clone());
outboard.add_control_surface(elevator);

let constraints = TrimPresets::level_trimmed_flight(0.5, "elevator");
```

The "both sections" requirement is not a formality — a definition on one section
alone is silently inert for that interval. See
[Control Surfaces](../guide/control-surfaces.md).

## When to use the other trim solver

`TrimSolver` solves coefficient targets on a section model and reports a
trimmed state. It gives no sensitivities.

If trim sits inside an optimization or a coupled solve, you want
`MeshTrimSolver` instead: dimensional force and moment balances, prescribed
external loads, and exact derivatives of the converged mapping. See
[Trim](../guide/trim.md) and [Mesh Trim](./mesh-trim.md).
