# Product

<!-- impeccable:product-schema 1 -->

## Platform

web

## Users

Aerospace engineers and researchers doing aircraft analysis, stability and
control, trim, and multidisciplinary design optimization (MDO). They are fluent
in aerodynamics and in AVL-style tools; their Rust experience varies. They come
to the documentation in the middle of engineering work, with an aircraft model
or an optimization loop already in mind. (Confirmed by the maintainer.)

Rust developers without an aerodynamics background are a secondary audience and
are not the reader the documentation is written for.

## Product Purpose

Zephyr (`zephyr-vlm`, imported as `zephyr_vlm`) is a vortex-lattice-method
library for aircraft geometry, aerodynamic loads, stability and control
derivatives, trim, eigenmodes, analytical and adjoint gradients, AVL input, and
VTK output.

The documentation succeeds when a reader can do two things (confirmed by the
maintainer):

1. Run a first analysis: add the crate, describe a wing, get loads, and find
   the API page they need.
2. Wire Zephyr into an MDO loop: use the analytical and adjoint gradients and
   the Philote servers inside an optimizer.

## Positioning

A safe, Rust-only VLM with analytical and adjoint gradients built in. No C,
C++, or Python runtime dependency; Philote transport is optional. Aerodynamic
results are checked against independent OptVL fixtures.

## Operating Context

- Readers work from AVL geometry and mass files, OptVL results, and MDO
  frameworks reached through Philote (gRPC).
- Results are consumed as coefficients, derivatives, eigenvalues, and VTK files
  opened in ParaView or VisIt.
- The documentation is a Docusaurus site published at
  `chrislupp.github.io/Zephyr/`, with versioned snapshots of released docs and
  a link out to rustdoc on docs.rs for item-level API reference.

## Capabilities and Constraints

- Double precision throughout. Right-handed body frame: x aft, y starboard,
  z up.
- Deliberate unit split: section twist and surface incidence in degrees; flight
  angles, rates, control deflections, trim variables, and derivatives in
  radians.
- Features: `vtk` (default), `parallel`, `philote` (off by default).
- The source forbids unsafe code.
- No document outside `CHANGELOG.md` states the crate's version; install
  snippets use `cargo add zephyr-vlm`.
- Released documentation snapshots under `docs/versioned_docs/` are frozen
  history and are not rewritten.
- Documentation pages contain prose, Rust and shell code, KaTeX mathematics,
  and tables. The only figure is Plate 1 on the contents page, drawn from the
  lattice the library builds for the Quick Start wing.

## Brand Commitments

- Name: Zephyr. Formerly Pegasus.
- Distribution Statement A must remain reachable from every page. The rail
  lists the Copyright and License sheet on every page, and that sheet holds the
  copyright, the licence terms, and the clearance line: approved for public
  release, distribution unlimited, case AFRL-2026-2198. The footer is the
  copyright line alone — no link columns, no clearance line. (Confirmed by the
  maintainer.)
- License: GPL-3.0-or-later. Copyright Christopher A. Lupp.
- The documentation is a research library's reference, not a marketing site
  (confirmed by the maintainer). It should not read as the stock Docusaurus
  template, and it should not keep the earlier near-black, glowing-cyan look.

## Evidence on Hand

- Validation: OptVL 2.1.1 fixtures covering full and symmetric section wings, a
  wing-tail configuration, and 680-panel uCRM geometry (`tests/`).
- Runnable examples in `examples/`, including the Philote servers.
- No logo, photography, or screenshots exist; the favicon is a drawn mark
  (red bar and the letter Z). No benchmark or adoption claims are on hand;
  none may be invented.

## Product Principles

1. The reader is an engineer in the middle of work. Answer the question, then
   get out of the way.
2. Conventions are part of the result. Axes, units, and sign conventions are
   stated where a reader will meet them, not buried.
3. Claims are checkable. Anything the documentation asserts about accuracy or
   behavior points at a fixture, an example, or the source.
4. One name per thing. Terms follow AVL and the crate's own API so the prose,
   the code, and the reader's existing vocabulary agree.
