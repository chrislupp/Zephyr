# AI Agent Guidelines for Project Pegasus

Pegasus is a safe Rust 2024 vortex-lattice-method library. The crate package is
`pegasus-vlm`, its import name is `pegasus_vlm`, and the pinned toolchain is
Rust 1.93.1.

## Required checks

Run these before handing off code changes:

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-features
cargo test --release --all-features
cargo test --no-default-features
RUSTDOCFLAGS="-D warnings" cargo doc --all-features --no-deps
cargo package
```

Coverage is enforced in CI with `cargo llvm-cov --fail-under-lines 92`.
Docusaurus changes should additionally pass `npm ci && npm run build` from
`docs/`.

## Layout

```text
src/
  geometry/       aircraft, sections, airfoils, bodies, direct meshes
  solver/         lattice generation, Biot-Savart, AIC, cached solutions
  aero/           near/far-field loads and stability derivatives
  gradients/      safe analytical kernels, adjoints, geometry sensitivities
  trim/           constraints, presets, Newton solver
  dynamics/       mass properties and eigenmode analysis
  io/             AVL geometry and mass readers
  visualization/  feature-gated vtkio exports
tests/            feature tests and frozen v0.2.0 parity fixtures
examples/         native Rust workflows and CSV scaling runner
benches/          Criterion benchmark
docs/             Docusaurus site
```

The source forbids unsafe code. Preserve double precision, panel ordering,
coordinate/sign conventions, and the deliberate unit split: geometry
incidence and twist are degrees; flight angles, rates, controls, trim, and
derivatives are radians.

The legacy behavioral oracle is tag `v0.2.0`, commit `556b0ea`. C++, CMake,
Philote, and Philote-dependent OpenMDAO code do not belong in the 0.3 tree.

`CHANGELOG.md` follows Keep a Changelog. Every user-visible change must update
its Unreleased section.
