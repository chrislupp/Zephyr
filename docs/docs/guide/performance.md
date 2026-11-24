---
sidebar_position: 13
---

# Performance

A vortex-lattice solve has two costs that grow differently. Assembling the
influence operators is $O(N^2)$ kernel evaluations; factoring the normalwash
matrix is $O(N^3)$. Storing the operators is $O(N^2)$ doubles — and since there
are two of them, one $N \times N$ and one $3N \times N$, memory is usually the
first wall a large model hits, not time.

Everything in this sheet is about not paying those costs more than once.

## First, reuse

The cheapest optimization is not repeating work. Keep one solver alive across a
sweep: changing alpha, beta, the rates, or dynamic pressure reuses the lattice,
both operators, the factorization, and the six unit responses. See the
invalidation table in [Solver Workflow](./solver-workflow.md).

The second cheapest is not computing what you will not read. `OutputSelection`
in `GradientCalculator` skips whole adjoint solves; the request interface in
[Selective Analysis](./selective-analysis.md) skips assembly stages too.

## The three symmetry accelerations

All three require surfaces explicitly duplicated about a plane, and they are
independent of one another.

### Mirror reuse

```rust
solver.set_mirror_reuse(true);
```

Influence assembly derives the mirrored half of **both** operators by
reflection instead of evaluating the kernel for it, roughly halving the kernel
evaluations. Results are **bit-identical** to having it off: the fast path is
taken only when assembly can prove the mesh satisfies the reflection identity
exactly, and the horseshoe kernel is reflection-equivariant in floating point by
construction, not merely in exact arithmetic.

A mesh that does not qualify — image planes in use, partially duplicated
geometry, or a mirror plane that is not a pure sign flip — evaluates every entry.
`InfluenceMatrixBuilder::mirror_reuse_applies` reports which happened.

Off by default. `Analysis::set_mirror_reuse` is the same setting on the request
interface.

### Symmetric factorization

```rust
solver.set_symmetric_factorization(true);
println!("{}", solver.symmetric_factorization_applies());
```

When the system matrix is exactly mirror-symmetric it has the block form
$\begin{pmatrix} P & Q \\ Q & P \end{pmatrix}$, whose eigenvalues are those of
$P + Q$ together with those of $P - Q$. Factoring those two half-size systems
costs a quarter of the full factorization, and neither half is worse conditioned
than the full matrix. The transpose has the same structure, so adjoint solves
benefit too.

The structure holds for an undeflected mirrored lattice at any flight condition
— angles and rates only enter the right-hand side — and for symmetrically
deflected controls. It does **not** hold for antisymmetric deflections such as
an aileron, for partially duplicated geometry, or for a half-wing model. Those
fall back to the full LU silently, with results bit-identical to having the
setting off. `symmetric_factorization_applies()` tells you which was used.

Results under the half-size path agree with the full factorization to the
rounding of the elimination, not bit for bit. That is the one acceleration here
that is not exact.

Off by default.

### Half-wing mode

`SolverMode::HalfWing` in `analysis::Analysis` goes further: it builds only one
semispan's panels and factors $A_{\text{same}} + A_{\text{reflected}}$, so there
is no second set of unknowns at all. Reflected sources still influence the
stored half, as images.

It carries real restrictions — zero sideslip, zero roll and yaw rate, symmetric
control reflection signs, and no antisymmetric derivatives — and it never
silently falls back. [Selective Analysis](./selective-analysis.md) has the
complete list.

Neither mirror reuse nor symmetric factorization has any effect in half-wing
mode; there are no mirror pairs left to exploit.

## The `parallel` feature

```bash
cargo add zephyr-vlm --features parallel
RAYON_NUM_THREADS=10 cargo run --release --example selective_analysis
```

`parallel` multithreads two things through `rayon`:

- influence assembly, over source panels,
- the geometry reverse sweep in the adjoints, over target rows.

Results are **independent of the thread count** and identical to the serial
build. Each source column writes a disjoint slice of the operators and reads
nothing another column writes; the reverse sweep splits rows into a fixed number
of chunks whose partial column sums are added back in chunk order. That is a
design constraint, not an accident — a derivative that changes with the size of
the thread pool is not usable in an optimizer.

Systems with fewer than 128 circulation unknowns use the serial loops: below
that, the coordination costs more than the work. The LU factorization is serial
either way.

The feature deliberately does **not** enable `faer/rayon`. The gemm
cache-blocking parameters differ between its sequential and parallel paths,
which changes the order in which partial panels accumulate in the LU — enough to
perturb the last unit in the last place, and so enough to break the
thread-independence property above.

Off by default, so a standard build is single-threaded.

## Measuring instead of guessing

`AnalysisResult::work` reports what a request actually did. It is the tool for
confirming that a request avoids work, which is a different question from how
long it took:

| Field | Counts |
| --- | --- |
| `unknowns` | circulation unknowns solved |
| `operator_bytes` | stored AIC and velocity-operator bytes, excluding the LU |
| `influence_evaluations` | kernel evaluations |
| `factorizations`, `primal_solves`, `adjoint_solves` | linear-algebra calls |
| `velocity_products` | velocity-operator products |
| `near_field_panels`, `moment_panels`, `force_derivative_panels` | integration work |
| `trefftz_interactions` | far-field pair interactions |
| `mach_interactions`, `geometry_interactions` | derivative-stage work |
| `alpha_derivative_panels`, `beta_derivative_panels` | per-input derivative work |

A cache hit reports zero work. Re-evaluating an identical request is the
simplest way to see the caching behave:

```rust
let first = analysis.evaluate(condition, request)?;
let second = analysis.evaluate(condition, request)?;
println!("{:?}", first.work);
println!("{:?}", second.work); // zeros
```

For wall-clock scaling there is a Criterion benchmark that solves a
six-chordwise mesh at spanwise counts from 8 to 96 strips, measuring solver
construction, meshing, assembly, factorization, and the solve together:

```bash
cargo bench --bench analysis_scaling
```

`cargo run --release --example analysis_scaling` runs the same sweep as a plain
CSV-emitting program, mirroring the legacy harness iteration for iteration: one
untimed warm-up, then fifty timed repetitions, with mesh construction outside
the timed region.

Measure your own case. This project publishes no timing claims.

## Build settings

The workspace sets `lto = "fat"` and `codegen-units = 1` for release builds.
That matters here for a specific reason: the influence kernel is defined in a
different module from the assembly loop that calls it $N^2$ times, and without
fat LTO it is not reliably inlined. The loop then loses the instruction-level
parallelism that hides its divide and square-root latency.

Cargo ignores a library's profile settings when it is consumed as a dependency,
so those settings govern only this workspace's own examples, benchmarks, and
tests. **If you depend on `zephyr-vlm` and care about assembly throughput, set
them in your own release profile.**

```toml
[profile.release]
lto = "fat"
codegen-units = 1
```

Always benchmark in release mode. The debug build of a kernel called $N^2$ times
is not slower by a constant factor you can reason past.
