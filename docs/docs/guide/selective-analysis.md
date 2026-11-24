---
sidebar_position: 11
---

# Selective Analysis

`analysis::Analysis` evaluates an explicit request: which scalar outputs you
want, which detailed arrays, and which individual `(output, input)` derivative
pairs. Stages a request does not reach are not executed, and an optional
half-wing mode solves only one semispan of a symmetric model.

This is the interface to reach for when [the solver
workflow](./solver-workflow.md) computes more than you asked for — a complete
solution, all six unit responses, both influence operators — and you only want
two coefficients and one derivative.

It accepts aircraft sections and direct structured meshes. The solver, force,
gradient, stability, trim, dynamics, and visualization APIs described elsewhere
on this site keep their own behavior; nothing here changes them.

## Request only the quantities you need

```rust
use zephyr_vlm::analysis::{
    Analysis, AnalysisRequest, GradientRequest, Input, Output, SolverMode,
};

// `input` is a MeshInput whose surfaces describe one half and have
// y_duplicated = true, y_duplicate_plane = 0.0; `condition` is a
// FlightCondition with zero sideslip and zero roll and yaw rate.
let mut analysis = Analysis::from_mesh(&input);
let request = AnalysisRequest {
    mode: SolverMode::HalfWing,
    outputs: vec![Output::Lift, Output::FarFieldInducedDrag],
    gradients: vec![
        GradientRequest { output: Output::Lift, input: Input::Alpha },
        GradientRequest { output: Output::FarFieldInducedDrag, input: Input::MeshNodes(0) },
    ],
    ..AnalysisRequest::default()
};
let result = analysis.evaluate(&condition, &request)?;
```

A runnable version is in `examples/selective_analysis.rs`:

```bash
cargo run --release --example selective_analysis
```

See [Selective Analysis (example)](../examples/selective-analysis.md).

`AnalysisRequest::default()` requests nothing. `all_outputs()` requests every
scalar output, with no gradients or detailed arrays. `SolverMode::Full` is the
default. There is no implicit request for six unit responses or a full Jacobian.

`AnalysisResult.outputs` contains only requested scalars. Gradients are returned
as requested `(output, input)` pairs with either a scalar or a vector of node
cotangents. A gradient may be requested without returning its corresponding
scalar output. Duplicate requests are evaluated once. Unrequested results are
absent; zero is reserved for a computed or analytically known zero.

Available scalars:

| Output | Meaning |
| --- | --- |
| `Lift`, `Drag`, `Side` | Near-field force coefficients in the existing conventions |
| `Roll`, `Pitch`, `Yaw` | Existing body-axis moment coefficients |
| `NearFieldInducedDrag` | Near-field drag less profile/viscous drag |
| `ViscousDrag` | Section polar contribution plus aircraft constant profile drag |
| `FarFieldLift`, `FarFieldSide` | Trefftz-plane force coefficients |
| `FarFieldInducedDrag` | Trefftz-plane induced drag |
| `SpanEfficiency` | Existing far-field span-efficiency definition |

Derivative inputs include alpha, beta, Mach, normalized rates, individual moment
reference coordinates, reference area/chord/span, dynamic pressure, named controls,
mesh nodes on a selected input surface, and one section design variable. Flight
angles and controls use radians; section twist retains the existing degrees
convention. Mesh-node derivatives hold reference quantities fixed. If reference
quantities also depend on design variables, the caller supplies that chain rule.

`Details` selects full-span circulation, induced velocity, and complete
surface/strip load records independently. Returned ordering matches the existing
lattice: each original surface followed by its mirror. Scalar-only requests do
not build those output records. Detailed surface loads request the complete load
record, including its moment fields.

## What is skipped

The analysis plans dependencies before assembling matrices. Near-field quantities
need local velocities and forces. Far-field lift alone needs strip circulation,
with no Trefftz interaction loop. Far-field induced drag needs Trefftz interactions,
but no near-field forces or force-velocity matrix. Moments and detailed records
are accumulated only when requested quantities require them.

Each requested scalar gradient uses one adjoint shared across its requested
inputs. Mach interactions and geometry reverse sweeps run only when needed;
geometry variables share a reverse sweep instead of requiring one solve per node.
Normalization derivatives need no adjoint, and known-zero derivatives need no
flow solve when they are the only requests.

`result.work` reports actual counts for the current call: stored unknowns, kernel
interactions, factorizations, primal/adjoint solves, velocity products, force and
moment integration, and selected derivative stages. It is intended for checking
that a request avoids work, rather than estimating elapsed time.

Reuse an `Analysis` across requests. The circulation factorization is reused when
Mach and effective normals are unchanged. An unchanged flight condition reuses
circulation; adding a velocity-dependent request builds the missing velocity
operator without another factorization. Scalar and derivative results are cached, as are adjoints reused when adding new
input derivatives for an already evaluated output. Control changes reuse an
available velocity operator and update its normal projection.
Changing mode or replacing the source via `set_mesh`/`set_aircraft` resets the
relevant context. Failed evaluations cannot publish a reusable partial result.

## Half-wing restrictions and conventions

Every input surface must describe one half, with explicit duplication about
`y = 0`. Both left-half and right-half inputs are accepted, as are multiple
mirrored surfaces and distinct components. Stored surfaces cannot cross the
plane or lie entirely in it. Mixed unmatched surfaces and additional image
planes are outside this mode's initial scope.

Sideslip, roll rate, and yaw rate must be zero. Angle of attack, subsonic Mach,
and pitch rate are supported. Deflected controls must have symmetric reflection
signs. Beta, roll-rate, yaw-rate, and antisymmetric-control derivatives are
rejected because they require an antisymmetric response system.

The half-wing mode builds only half the panels and factorizes
`A_same_half + A_reflected_half`. Reflected sources still influence the stored
half; they are evaluated as images without a second set of unknowns. The source
component and surface distinctions used by the full solver's core regularization
are preserved, including the cross-surface Trefftz core for the reflected half.
Full-span arrays are reconstructed only when requested.

Coefficients use the supplied **full-aircraft** reference area, chord, and span;
these are not automatically rescaled. Moment reference points may lie away from
the symmetry plane. Mirrored force and moment contributions are transformed about
the actual reference point rather than simply doubling every coefficient.

Geometry derivatives tie the original and reflected geometry. Their reverse
chain includes the reflected contribution and the swapped vortex endpoints.
They are not sensitivities of two independently movable wings. Keep centerline
nodes on `y = 0` when parameterizing a joined wing; reported coordinate cotangents
are in the ambient input coordinates, and the caller applies its geometric
constraints. Symmetry-breaking shape perturbations require full mode.

Geometry adjoints retain the existing restrictions on deflected control normals
and section drag polars. Flight-condition and symmetric-control derivatives
support viscous section polars. Unsupported requests return an error; half-wing
mode never silently falls back to a full solve. The unrestricted full solver
remains available through `VortexLatticeSolver` and `GradientCalculator`.

This mode is distinct from `VortexLatticeSolver::set_symmetry`, which adds image
planes to an existing lattice, and from `set_mirror_reuse`, which preserves the
full system while reusing reflected influence entries.

Full mode has two opt-in accelerations for lattices built from y-duplicated
surfaces, both off by default. `Analysis::set_mirror_reuse` derives reflected
influence entries instead of evaluating them and is bit-identical to plain
assembly. `Analysis::set_symmetric_factorization` factorizes two half-size
systems when the system matrix is exactly mirror-symmetric, and otherwise falls
back to the full factorization; `symmetric_factorization_applies` reports which
was used, and results agree with the full factorization to rounding. Neither
setting has an effect in half-wing mode. [Performance](./performance.md) covers
both in context. The Philote discipline exposes them as its `mirror_reuse` and
`symmetric_factorization` options, described below.

## Philote options

With the `philote` feature, `VlmDiscipline` accepts the options below. See
[MDO Integration](./mdo-integration.md) for the rest of the discipline.

```json
{
  "surface_y_duplicated": [true],
  "analysis_mode": "half_wing",
  "analysis_outputs": ["CL", "CDi"],
  "analysis_gradients": {"CL": ["alpha"], "CDi": ["Wing"]}
}
```

Use the configured surface variable name (here `Wing`) for its mesh-node
partials. `analysis_mode` is `full` or `half_wing`. Ordinary output names are
`CL`, `CD`, `CY`, `Cl`, `Cm`, `Cn`, `CDi_near`, `CDv`, `CL_far`, `CY_far`, `CDi`,
and `e`. Scalar partial names are the existing flight/reference input names.

Supplying any analysis option enables the selective interface. Its default
output list is the existing six coefficients; its default partial list is empty.
An explicit empty output list requests none. The discipline declares and evaluates
only selected partial pairs. Because Philote partials refer to registered output
variables, outputs named by `analysis_gradients` are also registered and returned
by the discipline's `compute`; the core Rust API allows gradient-only results.

The analysis options cannot be combined with the discipline's `trim`, `outputs`,
`fixed_inputs` or `controls` options; setup rejects that combination. Without analysis
options, the discipline keeps its ordinary six-output/full-partials behavior,
including those options. The selective discipline caches its `Analysis` across `compute` and
`compute_partials` calls on identical inputs.

The boolean `mirror_reuse` and `symmetric_factorization` options (both default
`false`) call `Analysis::set_mirror_reuse` and
`Analysis::set_symmetric_factorization` on that `Analysis`, and apply the same
settings to the discipline's ordinary and trim solves. They may be combined with
the analysis options and with every other option; changing either discards the
cached `Analysis`. In `half_wing` mode they have no effect, bit for bit.

## Validation

The numerical tests compare the selective path with the existing solver, and
reduced results with explicit full wings. Scalar comparisons use
`abs(error) <= 1e-10 + 1e-8 * abs(reference)`. Gradient tests additionally compare
analytical results against directional finite differences with separate absolute
and relative allowances for finite-difference truncation and roundoff. These
bounds apply to the tested cases, not to arbitrarily ill-conditioned systems.
The selective algorithms may differ in rounding from the legacy solver.

`tests/analysis_requests.rs` covers outputs, detailed arrays, flight and geometry
derivatives, controls, viscous loads, multiple components, left/right half meshes,
cache transitions, restrictions, and counts proving that stages are skipped.
Existing API tests continue to run independently.
