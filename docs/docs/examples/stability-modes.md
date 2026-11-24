---
sidebar_position: 3
---

# Stability Modes

Stability derivatives from a solved wing, fed with mass properties and a flight
state into the longitudinal eigenproblem.

```bash
cargo run --example stability_analysis
```

Source: `examples/stability_analysis.rs`.

## What it does

```rust
use zephyr_vlm::aero::StabilityDerivativeCalculator;
use zephyr_vlm::core::{DEGREES_TO_RADIANS, Vec3};
use zephyr_vlm::dynamics::{EigenmodeAnalysis, FlightState, MassProperties};
use zephyr_vlm::geometry::{Aircraft, Section, Surface};
use zephyr_vlm::solver::{FlightCondition, VortexLatticeSolver};

fn main() -> zephyr_vlm::Result<()> {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(4);
    wing.set_spanwise_panel_count(12);
    wing.add_section(Section::new(Vec3::new(0.0, -4.0, 0.0), 1.2, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 4.0, 0.0), 1.2, 0.0));

    let mut aircraft = Aircraft::new("Stability example");
    aircraft.set_reference_area(9.6);
    aircraft.set_reference_chord(1.2);
    aircraft.set_reference_span(8.0);
    aircraft.set_moment_reference_point(Vec3::new(0.3, 0.0, 0.0));
    aircraft.add_surface(wing);

    let condition = FlightCondition {
        alpha: 4.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft)?;

    let derivatives =
        StabilityDerivativeCalculator.compute(&aircraft, &mut solver, &condition)?;
    println!("CL_alpha = {:.5} /rad", derivatives.by_alpha.lift);
    println!("Cm_alpha = {:.5} /rad", derivatives.by_alpha.pitch);
    println!("Cl_p     = {:.5}", derivatives.by_roll_rate.roll);

    let mass = MassProperties {
        mass: 900.0,
        i_xx: 800.0,
        i_yy: 1_400.0,
        i_zz: 1_900.0,
        ..MassProperties::default()
    };
    let flight = FlightState { velocity: 50.0, ..FlightState::default() };

    let modes = EigenmodeAnalysis::new()
        .analyze_longitudinal(&aircraft, &mass, &derivatives, &flight)?;
    for mode in modes.modes {
        println!(
            "{}: lambda={:.5}{:+.5}i, zeta={:.3}",
            mode.mode_type.name(),
            mode.eigenvalue.re,
            mode.eigenvalue.im,
            mode.damping_ratio()
        );
    }
    Ok(())
}
```

Note that `StabilityDerivativeCalculator::compute` takes the solver mutably and
solves for itself if needed — but the solution it uses must match the condition
you pass, or it returns `ZephyrError::InvalidSolution`.

## Reading the derivatives

Derivatives are grouped by input, so `by_alpha.lift` is $C_{L\alpha}$ and
`by_alpha.pitch` is $C_{m\alpha}$. A wing-alone model at a moment reference
ahead of the aerodynamic centre gives a negative $C_{m\alpha}$ and therefore a
positive static margin; `derivatives.static_margin` reports it in reference
chords.

$C_{l_p}$ — roll damping — is negative for any sane wing. Its magnitude is the
quantity most sensitive to spanwise resolution here, because it is dominated by
the outboard strips.

## The caveat on mode frequencies

This example is a **wing alone**. It has no tail, so it has no meaningful
pitch-rate damping, and the longitudinal modes it produces are illustrative
rather than physical. The point of the example is the plumbing.

Separately, read the rate-derivative scaling note in
[Dynamic Modes](../guide/dynamic-modes.md) before taking an absolute frequency
from any model. `EigenmodeAnalysis` applies the conventional $b/2V$ and $c/2V$
factors to the rate derivatives, which is a different normalization from the one
the solver's raw rates use. Mode *signs*, and therefore stability conclusions,
are unaffected.

## Variations

- **Lateral modes**: `analyze_lateral` builds the $\beta$, $p$, $r$, $\phi$
  system. Add a vertical fin, or the Dutch roll will be meaningless.
- **The full system**: `analyze_full` gives the 12-state matrix. It is
  block-decoupled by construction; there is no inertial or aerodynamic
  cross-coupling between the longitudinal and lateral blocks.
- **The matrix itself**: `longitudinal_matrix` and friends return the system
  without solving it, and `EigenmodeResults::system_matrix` is the exact matrix
  the reported eigenvalues came from — so you can modify it, add coupling, and
  solve it yourself.
- **Mass from a file**: `MassReader::read_path("aircraft.mass")?.total_properties()`
  gives a `MassProperties` with the parallel-axis transfer already done. See
  [AVL Input Files](../guide/avl-files.md).
