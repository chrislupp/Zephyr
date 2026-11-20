---
sidebar_position: 5
---

# VortexLatticeSolver

Use `set_aircraft` or `set_mesh`, optionally set symmetry, then call `solve`.
The solver exposes read-only mesh, AIC, velocity matrix, reference quantities,
cache validity, and current solution accessors. Scoped mesh mutation
invalidates all dependent state.

See the [`VortexLatticeSolver` rustdoc](https://docs.rs/zephyr-vlm/latest/zephyr_vlm/solver/struct.VortexLatticeSolver.html).
