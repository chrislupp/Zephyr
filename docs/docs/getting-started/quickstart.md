---
sidebar_position: 2
---

# Quick Start

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
aircraft.add_surface(wing);

let condition = FlightCondition {
    alpha: 5.0 * DEGREES_TO_RADIANS,
    ..FlightCondition::default()
};
let mut solver = VortexLatticeSolver::new();
solver.set_aircraft(&aircraft)?;
let solution = solver.solve(&condition)?;
let result = ForceCalculator::new().compute_for_aircraft(
    &aircraft, solver.mesh(), &solution, &condition,
);
println!("CL = {}", result.coefficients.lift);
# Ok::<(), zephyr_vlm::ZephyrError>(())
```

Geometry incidence and section twist use degrees. Flight angles, rates,
control deflections, trim variables, and derivatives use radians.
