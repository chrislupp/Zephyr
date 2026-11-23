# AI Agent Guidelines for Project Zephyr

Zephyr is a safe Rust 2024 vortex-lattice-method library. The crate package is
`zephyr-vlm`, its import name is `zephyr_vlm`, and the pinned toolchain is
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

## Versioning and releases

`CHANGELOG.md` follows Keep a Changelog. Every user-visible change must update
its Unreleased section.

The release workflow owns the version. Do not hand-edit `version` in
`Cargo.toml`, do not add a dated heading to `CHANGELOG.md`, and do not create a
`v*` tag by hand — `.github/workflows/release.yml` does all three.

For the same reason, no document outside `CHANGELOG.md` states the crate's
version. `README.md` and `docs/docs/getting-started/installation.md` install
with `cargo add zephyr-vlm`, which resolves the current release on its own, so
a bump leaves them correct and nothing has to be edited alongside it. Do not
introduce a pinned version into prose, install snippets or examples.

A release runs when a pull request into `main` merges carrying a `release` or
`prerelease` label plus exactly one of `major`, `minor`, or `patch`; a
`prerelease` additionally needs exactly one of `alpha`, `beta`, or `rc`. The
workflow derives the current version from the highest stable `v*` git tag, not
from `Cargo.toml` — the manifest is only a fallback used when no tag exists, so
a stale version there is harmless and self-correcting. It then writes the new
version via `scripts/update_copyright.py`, propagates it into `Cargo.lock`,
rewrites `## [Unreleased]` into a dated release heading with updated compare
links, snapshots the Docusaurus docs, commits, signs a tag, and attaches the
packaged `.crate` to a GitHub release.

The project is pre-1.0: a breaking change takes a `minor` bump.

The workflow publishes to crates.io once it has pushed the tag and created the
GitHub release, authenticating with the `CARGO_REGISTRY_TOKEN` repository
secret. That step runs `cargo publish --locked` against the committed release
tree and so needs no `--allow-dirty`, unlike the earlier `cargo package` step,
which carries that flag because the version and changelog edits are still
uncommitted at that point. Prereleases publish too, keeping their
`-alpha.N`/`-beta.N`/`-rc.N` qualifier.

A published version number can never be reused or overwritten, only yanked, so
the publish step first asks crates.io whether that exact version exists and
skips when it does, which keeps a re-run of the workflow safe.
