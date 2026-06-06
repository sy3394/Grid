# gamma5-Block Lanczos for the Wilson Dirac operator

Computes eigenvalues of `D_W` **directly** (not `H_W = g5 D_W`) using the
gamma5-Block Lanczos algorithm in the indefinite gamma5 inner product
`(u,v) = u^dag g5 v`.

Reference: S. Yamamoto, *gamma5-Block Krylov (Block Lanczos) Methods for the
Wilson Dirac Operator* (2026).

## Files

- `Grid/algorithms/iterative/Gamma5BlockLanczos.h` — the algorithm (header only):
  verified three-term block kernel, block-tridiagonal Ritz extraction, and a
  gamma5-metric **thick restart** (paired locking + deflation; avoids the
  explicit-restart failure modes — noise normalisation and doubler
  re-amplification).
- `examples/Example_gamma5_block_lanczos.cc` — driver (single precision).
- `examples/G5L/g5bl_702_sweep.slurm` — Frontier batch script (shift-invert sweep).
- `examples/G5L/compare_evals.py` — compare output vs a reference eval file.
- `examples/G5L/free_wilson_spectrum.py` — analytic free-Wilson spectrum (validation).

## Modes

- **DIRECT** (default): Lanczos on `D_W`. Converges *peripheral* eigenvalues
  (doublers, spectral edges). Interior near-real modes do **not** converge.
- **SHIFT-INVERT** (`--shift sigma`): Lanczos on `(D_W - sigma)^{-1}`, converging
  the eigenvalues of `D_W` closest to real `sigma` (the near-real interior modes
  around `Re ~ sigma`). Since the Wilson mass is additive,
  `D_W(m) - sigma = D_W(m - sigma)`, inverted by normal-equations CG (HPD
  `D^dag D`, robust where BiCGSTAB fails). Map back: `lambda = sigma + 1/theta`.
- **SWEEP** (`--shift-sweep lo:hi:n`): loop `n` shifts in `[lo,hi]` in ONE job,
  accumulating + de-duplicating the converged modes of each window. Resolves a
  whole interior band; no multiple submissions.

A real `sigma` keeps `(D_W - sigma)` gamma5-Hermitian, so the conjugate-pair
structure and the block recurrence are preserved.

## Verification & diagnostics

- `--check`: verify each eigenpair against the RAW operator,
  `||D_W u - lambda u|| / ||u||`, and use that raw residual as the acceptance
  filter. This is essential in shift-invert: near-`sigma` "ghost" Ritz values
  can have a deceptively small Lanczos (theta-space) residual but a large raw
  residual; `--check` rejects them.
- `--dense`: exact dense diagonalisation of `D_W` (small lattices only) as a
  reference; writes all `N = 12*Volume` eigenvalues to `<out>.dense`.
- Performance counters (printed in shift-invert mode): number of inverter
  solves and total inner CG iterations, quoted against a single common
  Wilson-inverter solve (the baseline).

## Test configurations

- `--cold`            : free field (analytic spectrum).
- `--weak eps`        : slight perturbation of the free field,
  `U_mu = exp(i eps * random algebra)` -- a non-trivial INTERACTING background
  whose spectrum sits near (but shifted from) the free one. Good middle ground:
  a fully random (HOT) config has a dense, unseparated spectrum that converges
  poorly; a weak field keeps modes separated while breaking the free degeneracy.
- `--config PATH`     : a NERSC gauge configuration.

## Build (already-built Grid tree)

The example `#include`s the header directly, so no edit to `Algorithms.h` is
needed. Drop the two source files in place, then either register the example via
`scripts/filelist` + reconfigure and `make`, or compile it by hand against
`libGrid.a` with your platform's flags.

## Validation (machine-checkable, config-independent)

- Free field (cold 4^4, m=0): D_W eigenvalues match the analytic free-Wilson
  spectrum to ~8 digits; thick restart locks 8/8 pairs.
- Shift-invert (cold 8^4, sigma=0.3): recovers the analytic eigenvalue nearest
  sigma, `(0.0761205, +/-0.3826834)`, to 6-7 digits (residual ~1e-6).
- Sweep (cold 8^4, sigma in {0,0.3,0.6}): accumulates the `(0.07612, +/-0.38268)`
  and `(0.36928, +/-0.80380)` bands, matching analytic.

## Run on Frontier

Edit `g5bl_702_sweep.slurm` (account, config path, OUT, SWEEP), then
`sbatch g5bl_702_sweep.slurm`. Compare:

```bash
python3 G5L/compare_evals.py evals_DW_702_2.dat g5bl_702_md3.dat --mdtime 3.0
```

For the 702 reference cluster (`Re in [0.749, 0.812]`, `m = 0`), the default
sweep `0.70:0.85:8` tiles the band; widen `SWEEP` to scan more of the interior.
