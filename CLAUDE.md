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

The legacy behavioral oracle is tag `v0.2.0`, commit `556b0ea`. C++ and CMake
do not belong in the 0.3 tree.

Philote lives behind the non-default `philote` feature in `src/philote/`, built
on the `philote-mdo` crate. Nothing outside that module and its two example
servers may reference it, and the feature must stay off by default: the gRPC
stack and async runtime it pulls in are not dependencies of the library itself.
`cargo clippy --all-targets` and `cargo test` with the feature off must stay
clean, which is what keeps the core crate's dependency graph unchanged for
consumers who do not want MDO transport.

`CHANGELOG.md` follows Keep a Changelog. Every user-visible change must update
its Unreleased section.
