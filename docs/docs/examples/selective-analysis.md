---
sidebar_position: 7
---

# Selective Analysis

Two outputs and one derivative from a half-wing model, with the work counters
printed so you can see what was skipped — and then the same request again, to
see the caching.

```bash
cargo run --release --example selective_analysis
```

Source: `examples/selective_analysis.rs`. Run it in release mode; the point of
the example is work avoided, and a debug build buries that in overhead.

## What it does

```rust
use zephyr_vlm::analysis::{Analysis, AnalysisRequest, GradientRequest, Input, Output, SolverMode};
use zephyr_vlm::core::Vec3;
use zephyr_vlm::geometry::{MeshInput, SurfaceMesh};
use zephyr_vlm::solver::FlightCondition;

fn main() -> zephyr_vlm::Result<()> {
    let mut wing = SurfaceMesh::new("wing", 4, 16);
    wing.y_duplicated = true;
    for i in 0..=wing.nx {
        for j in 0..=wing.ny {
            wing.set_node(i, j, Vec3::new(i as f64 / 4.0, 5.0 * j as f64 / 16.0, 0.0))?;
        }
    }
    let mut input = MeshInput::default();
    input.surfaces.push(wing);
    input.reference.area = 10.0;
    input.reference.span = 10.0;

    let condition = FlightCondition {
        alpha: 5.0_f64.to_radians(),
        ..FlightCondition::default()
    };
    let request = AnalysisRequest {
        mode: SolverMode::HalfWing,
        outputs: vec![Output::Lift, Output::FarFieldInducedDrag],
        gradients: vec![GradientRequest {
            output: Output::Lift,
            input: Input::Alpha,
        }],
        ..AnalysisRequest::default()
    };

    let mut analysis = Analysis::from_mesh(&input);
    let result = analysis.evaluate(&condition, &request)?;

    println!(
        "CL = {:.10}, CDi = {:.10}",
        result.outputs[&Output::Lift],
        result.outputs[&Output::FarFieldInducedDrag]
    );
    println!("dCL/dalpha: {:?}", result.gradients[0].value);
    println!("Work: {:?}", result.work);
    println!(
        "Repeated request: {:?}",
        analysis.evaluate(&condition, &request)?.work
    );
    Ok(())
}
```

## What the request says

`AnalysisRequest::default()` requests **nothing**. Every field here is an
explicit opt-in:

- `outputs` asks for near-field lift and Trefftz induced drag. No moments, no
  viscous drag, no span efficiency.
- `gradients` asks for exactly one `(output, input)` pair. Requesting a gradient
  does **not** implicitly request its scalar output — they are listed
  separately because you may want only one of them.
- `details` is left at its default, so no circulation vector, no induced-velocity
  matrix, and no surface or strip load records are built.
- `mode: SolverMode::HalfWing` solves only the stored semispan.

`result.outputs` is a map containing only what was asked for. An absent entry
means not computed; zero is reserved for a computed or analytically known zero.

## What gets skipped

The analysis plans dependencies before assembling anything. For this request:

- near-field lift needs local velocities and forces, so the force-velocity
  operator is built;
- Trefftz induced drag needs the strip-circulation interaction loop but no
  near-field forces of its own;
- no moment accumulation, because nothing requested a moment;
- one adjoint, shared by the single gradient request.

`result.work` reports the actual counts — unknowns, kernel evaluations,
factorizations, primal and adjoint solves, Trefftz interactions, integration
panels, and operator bytes. It is for proving that a request avoids work, not
for estimating elapsed time. See the table in
[Performance](../guide/performance.md).

The second `evaluate` of an identical request prints a `Work` of zeros: the
scalar results and the adjoint are cached. Reuse one `Analysis` across a sweep
and an unchanged flight condition reuses the circulation, while adding a
velocity-dependent request builds only the missing operator rather than
refactoring.

## Half-wing mode

The mesh describes one half with `y_duplicated = true` about `y = 0`, and the
reference area and span are the **full aircraft** values — half-wing mode does
not rescale them for you. The reported `CL` is the whole wing's.

The mode builds only half the panels and factors
$A_{\text{same}} + A_{\text{reflected}}$: reflected sources still influence the
stored half, as images, without a second set of unknowns.

It comes with restrictions. Sideslip, roll rate, and yaw rate must be zero;
deflected controls must be symmetric; and beta, roll-rate, yaw-rate, and
antisymmetric-control derivatives are rejected outright rather than silently
falling back to a full solve. The complete list is in
[Selective Analysis](../guide/selective-analysis.md).

## Variations

- **Full geometry**: drop `mode` to get `SolverMode::Full`, where the
  restrictions do not apply and the two symmetry accelerations become available
  instead.
- **Geometry derivatives**: add
  `GradientRequest { output: Output::FarFieldInducedDrag, input: Input::MeshNodes(0) }`.
  The result is a vector of node cotangents, and the geometry variables share
  one reverse sweep.
- **Everything**: `AnalysisRequest::all_outputs()` requests all twelve scalars,
  with no gradients and no detail arrays.
- **Detail arrays**: set `details.surface_loads = true` to get the complete load
  records, expanded to the full span even in half-wing mode.
