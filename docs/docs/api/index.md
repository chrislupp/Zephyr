---
sidebar_position: 1
---

# API Reference

The crate is imported as `zephyr_vlm`. Its public modules are:

- `geometry`: owned aircraft, surfaces, sections, airfoils, bodies, and meshes
- `solver`: lattice generation, influence matrices, and cached primal solves
- `aero`: loads, grouped coefficients, and stability derivatives
- `gradients`: analytical kernels, adjoints, and high-level gradients
- `trim` and `dynamics`: flight constraints, trim, mass, and eigenmodes
- `io` and `visualization`: AVL/mass input and VTK XML output

Generated API documentation is available through rustdoc and
[docs.rs](https://docs.rs/zephyr-vlm).
