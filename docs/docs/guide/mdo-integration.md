---
sidebar_position: 5
---

# MDO Integration

The unreleased Zephyr implementation is a Rust library. Geometry and mesh
gradients can be consumed directly by a Rust optimizer through
`GradientCalculator`.

Philote servers, `MeshDiscipline`, `VlmDiscipline`, and Philote-dependent
OpenMDAO workflows are intentionally deferred. The next release does not
provide a C ABI, C++ facade, or Python binding. Standalone comparison scripts
that do not depend on Philote may still be used outside the library.
