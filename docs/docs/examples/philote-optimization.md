---
sidebar_position: 9
---

# Philote Optimization

A SciPy SLSQP optimizer sizing a wing and tail through the explicit Philote
gRPC discipline, with trim solved inside each aerodynamic evaluation.

This is the end-to-end demonstration of the MDO path: a Python optimizer, a
Rust server, analytical design gradients, and a trim solve the optimizer never
sees inside.

## Running it

Build the bridge and start a fresh server:

```bash
cargo build --features philote --example philote_trim_optimization
cargo run --features philote --example philote_vlm_server
```

Then, in another terminal, the Python driver:

```bash
pip install numpy scipy
python examples/trim_optimization.py
```

Sources: `examples/philote_trim_optimization.rs` (a JSON-line bridge from the
optimizer to the gRPC server), `examples/trim_optimization.py` (the driver), and
`examples/philote_vlm_server.rs` (the server). Start a **fresh** server: the
discipline's registry and trim settings are fixed after setup.

## The problem

Design variables: wing halfspan and tail chord. The objective is near-field
aerodynamic drag plus a simple analytical area cost standing in for profile drag.
Geometry bounds define the sizing range.

Every aerodynamic evaluation solves $F_z + F_{z,\text{external}} = 0$ and
$M_y = 0$ for angle of attack and elevator. Trim is solved inside the
discipline throughout the optimization, so the optimizer sees a trimmed
explicit function of the design variables.

Design-dependent weight, deformed mesh coordinates, and reference dimensions all
enter the analytical design gradient. There is no finite differencing anywhere in
the chain.

## What the script checks

The driver does not merely run — it verifies. It compares the analytical
objective and trim-state gradients against independent perturbed re-solves at
both the initial and the final design, and it requires optimization success,
objective improvement, and converged trim at every evaluation.

Results and the iteration history are written to JSON. A recorded run is
committed at `examples/trim_optimization_results.json`, including the gradient
checks and the solved operating states. That run converged in 7 iterations,
reducing the objective from 32.481917 to 27.833541 N, with maximum $F_z$ and
$M_y$ residuals of 4.55e-13 N and 1.80e-13 N·m.

## What it is not

A sizing demonstration with prescribed deformation and a prescribed centre of
gravity, balancing the selected $F_z$ and $M_y$ components. It is not an
aircraft design study: there is no structural model, no weight build-up beyond
the analytical stand-in, no atmosphere, and no attempt at the remaining four
equilibrium equations.

## A smaller wire example

For the request and response shapes without an optimizer in the loop, start a
fresh server and run the client:

```bash
cargo run --features philote --example philote_vlm_server
# in another terminal
cargo run --features philote --example philote_trim_client
```

The client configures a deformed wing–tail trim, discovers which variables the
discipline solved, requests functions and analytical gradients, and changes an
external load.

## Without transport

The same trimmed mapping and the same derivatives are available natively,
without a gRPC stack or an async runtime, through `MeshTrimSolver` and
`MeshTrimSolution::compute_partials`:

```bash
cargo run --example mesh_trim
```

See [Mesh Trim](./mesh-trim.md), and
[MDO Integration](../guide/mdo-integration.md) for the discipline options,
control mixers, load conventions, and gradient semantics.
