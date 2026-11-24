---
sidebar_position: 8
---

# Analysis

`analysis::Analysis` is the request-driven interface: you state which scalar
outputs, which detail arrays, and which individual `(output, input)` derivative
pairs you want, and stages the request does not reach are not executed. It also
offers a half-wing mode that solves one semispan of a symmetric model.

The guide sheet, [Selective Analysis](../guide/selective-analysis.md), covers the
request API, the symmetry restrictions, caching, and the Philote options in
full. This sheet is the type reference.

```rust
use zephyr_vlm::analysis::{Analysis, AnalysisRequest, GradientRequest, Input, Output, SolverMode};

let mut analysis = Analysis::from_mesh(input);
let request = AnalysisRequest {
    mode: SolverMode::HalfWing,
    outputs: vec![Output::Lift, Output::FarFieldInducedDrag],
    gradients: vec![GradientRequest { output: Output::Lift, input: Input::Alpha }],
    ..AnalysisRequest::default()
};
let result = analysis.evaluate(condition, &request)?;
```

## Construction

| Method | Effect |
| --- | --- |
| `from_mesh(&MeshInput)`, `from_aircraft(&Aircraft)` | build a context |
| `set_mesh`, `set_aircraft` | replace the geometry, dropping every cache but keeping the opt-in settings |
| `set_mirror_reuse(bool)`, `set_symmetric_factorization(bool)` | the solver accelerations; changing either discards operators, factorization, and cached results |
| `symmetric_factorization_applies()` | which factorization the last solve used |
| `evaluate(&FlightCondition, &AnalysisRequest)` | run a request |

Reuse one `Analysis` across requests. That is where the caching lives.

## `AnalysisRequest`

| Field | Meaning |
| --- | --- |
| `mode` | `SolverMode::Full` (default) or `HalfWing` |
| `outputs` | the scalars wanted |
| `gradients` | `(output, input)` pairs |
| `details` | optional full-span arrays |

`AnalysisRequest::default()` requests **nothing** and performs no aerodynamic
work. `all_outputs()` requests all twelve scalars with no gradients and no
details. Requesting a gradient does not implicitly request its scalar output.
Duplicate requests are evaluated once.

## `Output`

`Lift`, `Drag`, `Side`, `Roll`, `Pitch`, `Yaw`, `NearFieldInducedDrag`,
`ViscousDrag`, `FarFieldLift`, `FarFieldSide`, `FarFieldInducedDrag`,
`SpanEfficiency`. `Output::ALL` is the array of all twelve.

Every one uses the full-aircraft normalization, half-wing mode included — the
references are not rescaled for you.

## `Input`

`Alpha`, `Beta`, `RollRate`, `PitchRate`, `YawRate`, `Mach`,
`DynamicPressure`, `MomentReference(index)`, `ReferenceArea`,
`ReferenceChord`, `ReferenceSpan`, `Control(name)`, `MeshNodes(surface_index)`,
and `Section { surface, section, variable }`.

Angles and controls are radians; section twist is degrees. Analytically known
zeros — the dynamic-pressure derivative of a coefficient, for instance — are
returned as exact zero without a flow solve when they are the only requests.

## `Details`

`circulation`, `induced_velocity`, and `surface_loads`, each independently
selectable. Scalar-only requests do not build the load records. Returned
ordering matches the lattice, and half-wing mode expands only the arrays you
asked for.

## `AnalysisResult`

| Field | Contents |
| --- | --- |
| `outputs` | a `BTreeMap<Output, Scalar>` of **only** the requested scalars |
| `gradients` | one `GradientResult` per request, each carrying the request and a `GradientValue` |
| `circulation`, `induced_velocity`, `surface_loads` | `Option`s, present only when requested |
| `work` | what this call actually did |

Omitted quantities are **absent**, never zero-filled. Zero is reserved for a
computed or analytically known zero.

`GradientValue` is either a scalar or a vector of node cotangents, depending on
the input.

## `Work`

Counters for the current call, zero on a cache hit: `unknowns`,
`operator_bytes`, `influence_evaluations`, `factorizations`, `primal_solves`,
`adjoint_solves`, `velocity_products`, `near_field_panels`, `moment_panels`,
`force_derivative_panels`, `trefftz_interactions`, `mach_interactions`,
`geometry_interactions`, `alpha_derivative_panels`, `beta_derivative_panels`.

It is for proving a request avoids work, not for estimating elapsed time. See
[Performance](../guide/performance.md).

[`analysis` rustdoc](https://docs.rs/zephyr-vlm/latest/zephyr_vlm/analysis/index.html)
