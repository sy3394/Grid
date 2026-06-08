# gamma5-Block Lanczos for the Wilson Dirac operator

Computes complex eigenvalues of `D_W` directly (not `H_W = g5 D_W`) using the
gamma5-Block Lanczos algorithm in the indefinite gamma5 inner product
`(u,v) = u^dag g5 v`.

Reference: S. Yamamoto, *gamma5-Block Krylov (Block Lanczos) Methods for the
Wilson Dirac Operator* (2026).

## Algorithm

`D_W` is self-adjoint in the gamma5 inner product, so the block Krylov recurrence
is three-term (block-tridiagonal `T_m`) at block size 2 — targeting conjugate
eigenvalue pairs together at low storage. Ritz **values** come from `T_m`; Ritz
**vectors** are extracted by **refined Ritz** (the Euclidean-residual minimiser
over the Krylov subspace), and convergence is judged by the raw Euclidean
residual `||D_W u - lambda u|| / ||u||`. A serious breakdown (a gamma5-neutral
residual) is handled by **look-ahead** (the block grows). `thickRestart` locks
converged conjugate pairs and deflates them.

## Files

- `Grid/algorithms/iterative/Gamma5BlockLanczos.h` — the algorithm (header only).
- `examples/Example_gamma5_block_lanczos.cc` — driver (double precision).
- `examples/G5L/g5bl_702_sweep.slurm` — Frontier shift-invert sweep.
- `examples/G5L/g5bl_history.slurm`   — Frontier residual-vs-Krylov-dim history.
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
--steps N         Lanczos steps per pass
--wanted N        wanted conjugate pairs
--cycles N        thick-restart cycles
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
- Shift-invert (cold): recovers the eigenvalue nearest `sigma` to ~8 digits.

## Comparison with Euclidean Arnoldi

Grid has no non-Hermitian Arnoldi/Krylov-Schur eigensolver, so the driver also
implements a **block Arnoldi seeded with `[v, g5 v]`** (`blockArnoldiG5`) plus
**refined Ritz extraction** (`refinedRitz`). This spans the *same* block Krylov
subspace as g5bl but in the Euclidean metric — a perfectly-conditioned orthonormal
basis, no oblique-projection penalty, standard restart. On weak 8^4 (shift-invert,
interior modes) it converges *identically* to g5bl, confirming the g5 metric buys
nothing for interior modes: the smooth convergence comes from refined extraction
plus block seeding, not the metric. The Euclidean version is the simpler, more
robust production choice.

## Performance analysis

Run `--history` to write `<out>.{g5bl,blockarnoldi,arnoldi}.hist` (the honest raw
Euclidean residual is `min_refined_raw` for g5bl/blockarnoldi, `min_raw_res` for
arnoldi; see each file header); plot with `~/BNL/G5BL/g5bl_performance.ipynb`.
