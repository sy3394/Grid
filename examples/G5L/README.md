# Interior eigenvalues of the Wilson Dirac operator

Computes complex eigenvalues of `D_W` directly (not `H_W = g5 D_W`).

Two solvers are provided (see "Comparison" below):

- **RefinedArnoldi** (`Grid/algorithms/iterative/RefinedArnoldi.h`) — single-vector
  Arnoldi + refined Ritz extraction. **The default**, and the better method for
  interior modes.
- **gamma5-Block Lanczos** (`Grid/algorithms/iterative/Gamma5BlockLanczos.h`) — the
  gamma5-metric block method (`--g5bl`); kept for comparison and for the manuscript.

Reference: S. Yamamoto, *gamma5-Block Krylov (Block Lanczos) Methods for the
Wilson Dirac Operator* (2026).

## Algorithm (gamma5-Block Lanczos)

`D_W` is self-adjoint in the gamma5 inner product, so the block Krylov recurrence
is three-term (block-tridiagonal `T_m`) at block size 2 — targeting conjugate
eigenvalue pairs together at low storage. Ritz **values** come from `T_m`; Ritz
**vectors** are extracted by **refined Ritz** (the Euclidean-residual minimiser
over the Krylov subspace), and convergence is judged by the raw Euclidean
residual `||D_W u - lambda u|| / ||u||`. A serious breakdown (a gamma5-neutral
residual) is handled by **look-ahead** (the block grows). `thickRestart` locks
converged conjugate pairs and deflates them.

## Files

- `Grid/algorithms/iterative/RefinedArnoldi.h` — the default solver (header only).
- `Grid/algorithms/iterative/Gamma5BlockLanczos.h` — the gamma5-block solver (header only).
- `examples/Example_gamma5_block_lanczos.cc` — demonstrator/comparison driver.
- `examples/Wilson_DW_spectrum.cc` — production driver: D_W eigenvalues in Re in [lo,hi].
- `examples/G5L/dw_spectrum_0to2.slurm` — Frontier production: D_W spectrum in [0,2] (RefinedArnoldi).
- `examples/G5L/g5bl_702_sweep.slurm` — Frontier production: shift-invert cluster sweep (RefinedArnoldi, or g5bl via --g5bl).
- `examples/G5L/compare_solvers.slurm` — Frontier test: RefinedArnoldi vs g5bl head-to-head (--history + --compare).
- `examples/G5L/compare_evals.py`     — compare output vs a reference eval file.
- `examples/G5L/free_wilson_spectrum.py` — analytic free-Wilson spectrum.

## Modes

- **DIRECT** (default): Lanczos on `D_W` — peripheral eigenvalues (doublers, edges).
- **SHIFT-INVERT** (`--shift sigma`): Lanczos on `(D_W - sigma)^{-1}` — eigenvalues
  of `D_W` nearest real `sigma`. Since the Wilson mass is additive,
  `D_W(m) - sigma = D_W(m - sigma)`, inverted by normal-equations CG. Map back:
  `lambda = sigma + 1/theta`. Real `sigma` keeps `(D_W - sigma)` gamma5-Hermitian.
- **SWEEP** (`--shift-sweep lo:hi:n`): n shifts in `[lo,hi]` in one job; accumulate
  and de-duplicate the converged modes of each window.

## Options

```
--config PATH     NERSC gauge config   (else --cold | --weak eps | hot random)
--cold            unit gauge (free field; analytic spectrum)
--weak eps        slight perturbation of the free field (interacting, separated)
--grid X.Y.Z.T    lattice dims
--mass m          Wilson bare mass
--shift sigma     single shift-invert about real sigma
--shift-sweep lo:hi:n   sweep n shifts
--steps N         base step count (default Krylov dim = 2*N)
--krylov N        RefinedArnoldi Krylov dimension (default 2*steps)
--g5bl            use gamma5-block Lanczos instead of the default RefinedArnoldi
--wanted N        (--g5bl only) wanted conjugate pairs
--cycles N        (--g5bl only) thick-restart cycles
--tol eps         convergence/locking tolerance (raw Euclidean residual)
--accept eps      keep sweep modes with raw residual < eps
--stol eps        inner CG tolerance
--out PATH        output .dat
--isolation-only  single non-restarted pass
--dense           exact dense diagonalisation (small lattices)
--compare         head-to-head: g5bl vs single-vector Arnoldi vs block Arnoldi
                  seeded with [v, g5 v] (same operator, equal Krylov dim)
--history         residual & #converged vs Krylov dim for all three methods
```

## Cost reporting

- **Krylov dim** = number of `(D_W - sigma)^{-1}` applications = outer steps
  (the comparison axis; equal for both methods at equal cost).
- **total CG iters** = inner-solve work (the absolute D_W-matvec cost).

## Build (already-built Grid tree)

The example `#include`s the header directly, so no `Algorithms.h` edit is needed.
Place the two source files, then either register the example via
`scripts/filelist` + reconfigure + `make`, or compile it against `libGrid.a`.

## Validation

- Free field (cold): D_W eigenvalues match the analytic free-Wilson spectrum.
- Shift-invert (cold): recovers the eigenvalue nearest `sigma` to ~12 digits (e.g.
  `lambda = (1-cos(pi/8)) + i sin(pi/8) = 0.076120 + 0.382683 i`).

## Why RefinedArnoldi is the default

Run `--compare` / `--history` to see it directly. On weak 8^4, shift-invert, interior
modes, at **equal Krylov dim and equal inner-CG cost**, RefinedArnoldi dominates
Gamma5BlockLanczos (e.g. Krylov dim 80, sigma=0.1: 24 distinct modes converged to 1e-6
and best residual 2e-12, vs g5bl's 5 modes and best 2e-7).  The reasons (established by
the decomposition runs, see git history):

- **Oblique projection hurts g5bl.** The standard g5-Galerkin Ritz vectors are
  Euclidean-suboptimal; refined extraction (used by both solvers here) is what gives
  good vectors.
- **The g5 metric then buys nothing** over a Euclidean Krylov space of the same
  dimension, and its indefinite basis is worse-conditioned.
- **Block seeding is a liability** for non-degenerate D_W: it halves the Krylov depth
  at fixed dimension, so single-vector Arnoldi resolves more interior modes.

The decisive ingredient is **refined extraction**; `RefinedArnoldi` packages it on a
well-conditioned Euclidean basis.
(Block seeding may still help a genuinely near-degenerate cluster; confirm on the
real 32^4 config before deciding.)

## Performance analysis

Run `--history` to write `<out>.{g5bl,blockarnoldi,arnoldi}.hist` (the honest raw
Euclidean residual is `min_refined_raw` for g5bl/blockarnoldi, `min_raw_res` for
arnoldi; see each file header); plot with `~/BNL/G5BL/g5bl_performance.ipynb`.
