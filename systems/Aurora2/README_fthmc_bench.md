# FTHMC kernel benchmark: optimised vs non-optimised rect-FTHMC paths

Measures the FTHMC-specific cost of the optimised `SmearedConfigurationRect`
defaults against the retained non-optimised (`int old`) reference
implementations, in ONE binary built from this tree
(`feature/rect_fthmc-optimise`), for the Lattice 2026 "Code Optimization"
backup slide.  Lives only in this project checkout; nothing to build or
patch anywhere else.

- Driver: `HMC/Benchmark_fthmc_kernels.cc`
- Job script: `systems/Aurora2/fthmc_kernel_bench.pbs`

## DEBUG must stay undefined

`GaugeConfigurationRect.h` ships with `#undef DEBUG`.  With `#define DEBUG`
the default routines self-check against the old ones on every call, so the
"opt" timings would include old-path work.  The driver enforces this at
compile time (`#error` if DEBUG leaks out of the header) — if it builds,
the check is taken into account.

## What is timed

Per kernel (`p`, `r`, `pp`, `pr`, `rp`, `rr` at rho = 0.12; a letter is one
masked step = 8 sub-levels; in two-letter names the **leftmost letter is
applied first** to the thin field, i.e. it is `mask_types[0]`), on 24^3x40
quenched, 5 warm-ups then 20 timed calls per component with
`accelerator_barrier()` + communicator barrier around every sample:

| component    | impl        | what                                                       |
|--------------|-------------|------------------------------------------------------------|
| `lndet`      | opt vs old  | `logDetJacobian()` vs `logDetJacobian(0)`                  |
| `jacforce`   | opt vs old  | `logDetJacobianForce(f)` vs `logDetJacobianForce(0,f)`     |
| `forward`    | single      | F: `set_Field` / `fill_smearedSet` (plain smearers; no old path) |
| `inverse`    | single      | F^-1: driver-side fixed-point inversion (see below)        |
| `xformforce` | single      | gauge-force pullback through F (`smeared_force`; no old overload) |
| `rawderiv`   | single      | raw Wilson deriv on the smeared field — context, not FTHMC cost |

The opt/old pairs are interleaved opt,old,opt,old per sample inside the
process, so both paths see identical thermal/clock state; no second build,
no cross-job variance.  Before timing, the driver prints untimed
consistency checks: `lndet` opt-old difference and the relative L2
difference of the two Jacobian forces.

The inverse has no library implementation: per level,
U_{k+1} = exp(-iQ(U_k)) V on the masked links with exp(-iQ(U_k)) =
U_k W_k^dag taken from the class's own smearer W_k = Stout(U_k), stopping
at relative L2 change 1e-13 (`--fptol/--fpmaxit`); iteration counts are
reported separately (`fp_iters_*`).  The level masks are transcribed from
the `SmearedConfigurationRect` constructor (private there); the driver
validates the transcription by printing |F^-1(F(U)) - U|/|U|
(`fp_relerr`, expect ~1e-13) before timing.

## Building

In the Aurora clone of this branch (`~/src/Grid_rect_opt`), after pulling:

```bash
cd ~/src/Grid_rect_opt
./scripts/filelist                     # regenerate HMC/Make.inc (new target)
autoreconf -fvi                        # Make.inc is included at automake time
cd systems/Aurora2
./config.status                        # refresh Makefiles, same configure flags
make -C HMC Benchmark_fthmc_kernels
```

If `--config` will point at an Aurora-GPU-written SciDAC file, the tree
also needs commit `c5fb9a59` ("IldgIO.h: tolerate wrong stored SciDAC
checksum on read"; currently on `g5-block-lanczos`) — Aurora GPU builds
wrote wrong checksum records; payloads are intact and the GRID_FIELD_NORM
record remains a hard check.  Without `--config` the fixed-seed hot start
needs nothing.

## Running

From a work dir on flare:

```bash
qsub ~/src/Grid_rect_opt/systems/Aurora2/fthmc_kernel_bench.pbs
```

One node, 12 tiles, `--mpi 2.2.3.1`, grid 24.24.24.40, one driver
invocation per (kernel, rep), `NREP=2` passes (drift check only — the
comparison is within-run).  Overridable via `qsub -v`: `GRID MPI NP PPN
KERNELS RHO NWARM NMEAS NREP DEVMEM CONFIG ROOT BIN`.  A thermalised
configuration gives production-representative fixed-point iteration
counts; the hot start is fine for the kernel timings themselves.

Walltime is 2 h for the full 6-kernel matrix; the old-path `jacforce` on
`r`-containing kernels is the unknown.  Smoke-run first:

```bash
qsub -q debug -l walltime=0:59:00 -v KERNELS=r,NREP=1 fthmc_kernel_bench.pbs
```

and scale up (or trim `NMEAS`) from the observed run time.

## Output

Each run prints its own per-kernel table (`FTHMCBENCH-TABLE`), and
`fthmc_kernel_bench_logs/SUMMARY.txt` ends with the slide table — per
kernel, component x impl, median wall time per call in ms (median over 20
calls, median across passes) and the old/opt speedup:

```
kernel component    non-opt[ms]      opt[ms]  speedup
p      lndet             ...          ...       ...x
p      jacforce          ...          ...       ...x
p      forward             -          ...        -
...
```

Only `lndet` and `jacforce` have a non-opt column: `forward`
(`fill_smearedSet`) and the driver-side `inverse` use the plain smearers,
which the optimisation does not touch, and `xformforce` has no old
overload — they are reported single-path for completeness of the FTHMC
cost budget.  Raw per-run medians are in `bench_raw.dat`; full logs per
(kernel, rep) sit next to it.
