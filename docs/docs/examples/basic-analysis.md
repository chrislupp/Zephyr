---
sidebar_position: 1
---

# Basic Analysis

An angle-of-attack sweep on a rectangular wing, reusing one solver across the
whole sweep. This is the smallest complete program that produces a result you
would put in a report.

```bash
cargo run --example basic_analysis
```

Source: `examples/basic_analysis.rs`.

## What it does

```rust
use zephyr_vlm::aero::ForceCalculator;
use zephyr_vlm::core::{DEGREES_TO_RADIANS, Vec3};
use zephyr_vlm::geometry::{Aircraft, Section, Surface};
use zephyr_vlm::solver::{FlightCondition, VortexLatticeSolver};

fn main() -> zephyr_vlm::Result<()> {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(6);
    wing.set_spanwise_panel_count(20);
    wing.add_section(Section::new(Vec3::new(0.0, -5.0, 0.0), 1.0, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 5.0, 0.0), 1.0, 0.0));

    let mut aircraft = Aircraft::new("Rectangular wing");
    aircraft.set_reference_area(10.0);
    aircraft.set_reference_chord(1.0);
    aircraft.set_reference_span(10.0);
    aircraft.set_moment_reference_point(Vec3::new(0.25, 0.0, 0.0));
    aircraft.add_surface(wing);

    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft)?;

    println!("alpha_deg,CL,CDi,e");
    for alpha_degrees in [-2.0_f64, 0.0, 2.0, 4.0, 6.0, 8.0] {
        let condition = FlightCondition {
            alpha: alpha_degrees * DEGREES_TO_RADIANS,
            ..FlightCondition::default()
        };
        let solution = solver.solve(&condition)?;
        let result = ForceCalculator::new().compute_for_aircraft(
            &aircraft,
            solver.mesh(),
            &solution,
            &condition,
        );
        println!(
            "{alpha_degrees:.1},{:.6},{:.6},{:.4}",
            result.coefficients.lift,
            result.far_field_induced_drag,
            result.span_efficiency
        );
    }
    Ok(())
}
```

The wing spans $y = -5$ to $y = 5$ with unit chord, so the reference area of 10
and span of 10 are the actual planform values. Six chordwise by twenty spanwise
panels, cosine spaced in both directions by default, gives 120 horseshoe
vortices. The moment reference sits at the quarter chord.

Output is CSV on standard output, so it pipes straight into a plot.

## The part worth copying

One solver, built once, used six times. Angle of attack enters only through the
right-hand side, so the lattice, both influence operators, the factorization,
and the six unit responses are all reused: the sweep costs one assembly and one
factorization, not six.

Building a fresh `VortexLatticeSolver` inside the loop would be a correct
program and roughly six times the work. See
[Solver Workflow](../guide/solver-workflow.md) for what else survives a change.

## Reading the numbers

`far_field_induced_drag` is the Trefftz-plane value, which is the one to report
on a lattice this coarse — it needs no chordwise resolution.
`span_efficiency` should come out near unity for a rectangular wing of this
aspect ratio and fall as the lift distribution departs from elliptic.

Compare against `AerodynamicResults::ideal_induced_drag()`, which is
$C_L^2/\pi\mathrm{AR}$, to see how far from elliptic the loading is.

## Variations

- **Taper and twist**: change the tip section's chord and twist. Twist is in
  degrees.
- **A half model**: give the wing `set_y_duplicate(0.0)` and move the root
  section to `y = 0`. The lattice is the same size; what changes is that the
  symmetry accelerations in [Performance](../guide/performance.md) become
  available.
- **Camber**: `root.set_airfoil(Airfoil::naca4(2412, 61)?)`. Expect a nonzero
  $C_L$ at zero alpha.
- **Spanwise loading**: ask for `results.surface_loads[0].strip_loads` and plot
  `normal_span_loading` against the strip's `quarter_chord_position.y`.
