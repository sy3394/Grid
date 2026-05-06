# Autocovariance Executables — Usage Guide

Sources: `autocovariance_VS.cc`, `autocovariance_NVS.cc`, `ACC.hpp`

---

## Vocabulary — read this first

This codebase distinguishes three independent concepts that are easy to confuse.
Throughout the source, docs, and outputs:

| Term | Axis | Meaning |
|---|---|---|
| **binning** | MD-chain partition | Splits the Markov chain into `n_bin` consecutive segments; each bin gives one estimate of ρ(t), variance comes from inter-bin sample variance |
| **blocking** / **sparsening** | spatial-lattice partition | Coarsens the per-site autocovariance field G(x,t) into spatial cells (block: average of `l_B^d` sites; sparse: corner site only) |
| **ACC** | observable | The autocorrelation coefficient ρ(t) = G(t)/G(0); **never** a label for the centered/connected form of G |
| **MFCOV** | observable | The unnormalized autocovariance G(t) itself |
| **G_cent / G_conn** | construction of G | Centered form (subtract mean, then product) vs connected form (product, then subtract product of means); they differ at O(1/V) |

In particular, **binning ≠ blocking**: binning lives on the MC-time axis,
blocking lives on the lattice. The two routines `binning_avg_cov` and
`binning_avg_rho` (formerly `binning` / `binning2`) both partition the MD chain
— they differ in what is averaged across bins, not in any spatial operation.

---

## Overview

Two Grid-based executables compute the equal-time spatial autocovariance function
of a Wilson-flowed observable field (energy density `E` or topological charge
density `topo`) from a sequence of pre-flowed lattice configurations:

| Executable | Estimator | Mean subtracted |
|---|---|---|
| `autocovariance_VS` | VS — Volume-Subtracted | per-config spatial mean `⟨A⟩_i = Σ_x A(x,i)/V` |
| `autocovariance_NVS` | NVS — Non-Volume-Subtracted | ensemble mean `⟪a⟫` (scalar, uniform across x) |

Both read pre-computed flowed observable fields from disk (SciDAC/ILDG format,
written by `ComputeWilsonFlow.cc`), accumulate the autocovariance field
`G[t](x) = ⟨(A(x,i)−m_i)(A(x,i+t)−m_{i+t})⟩_i`, and then estimate errors via
one of three methods selected by the input parameters.

For each spatial pre-blocking size `l_B` in `space_block_sizes` the code also
computes a **sparse-sampled** variant (sites with all coordinates divisible by
`l_B`) alongside the standard **block-averaged** variant.

Both executables also compute a second variant `G_conn` using the alternative
mean-subtraction convention (see source comments); its output is tagged
`G_conn` vs the primary `G_cent`.

### Three orthogonal axes

Every output row is determined by three independent choices:

| Axis | Choices | Notes |
|---|---|---|
| **Field construction** | `G_cent` (centered) vs `G_conn` (connected) | how `G(x,t)` is built from raw `A(x,i)`. `G_cent` subtracts means before the product; `G_conn` subtracts the product of means after. They differ at `O(1/V)` only |
| **Spatial coarsening** | `blocked` (cell = mean of `l_B^d` sites) vs `sparse` (cell = corner site) | how `G(x,t)` is reduced to a `V/l_B^d`-cell coarse lattice |
| **Error estimation** | `MS_approx`, `MF_approx`, `binning_avg_cov`, `binning_avg_rho` | which variance formula consumes the coarsened field |

Field construction and spatial coarsening are orthogonal preprocessing steps;
the error-estimation routine is selected at runtime by `MDtime_div_fac` and `R`.

### Combinations actually present in the output

|  | MS | MF | binning_avg_cov | binning_avg_rho |
|---|---|---|---|---|
| cent × blocked  (`G_cent_B`) | ✓ VS, NVS | ✓ VS only | ✓ VS, NVS | ✓ VS, NVS |
| cent × sparse   (`G_cent_s`) | ✓ VS, NVS | — | ✓ VS, NVS | ✓ VS, NVS |
| conn × blocked  (`G_conn_B`) | ✓ VS, NVS | ✓ VS only | ✓ VS, NVS | ✓ VS, NVS |
| conn × sparse   (`G_conn_s`) | ✓ VS, NVS | — | ✓ VS, NVS | ✓ VS, NVS |

> **MF + sparse is intentionally omitted.** `MF_approx` integrates the spatial
> covariance density `Cov[G(x,t), G(y,t)]` over all separations `|x−y| ≤ R`.
> Sparse sampling discards every separation that isn't a multiple of `l_B` —
> exactly the sub-`l_B` structure MF is designed to use. Sparse-MF is strictly
> noisier than blocked-MF with no diagnostic upside.

> **NVS has no MF.** Subtracting a single ensemble-mean scalar from every site
> breaks the spatial covariance decomposition that `MF_approx` relies on, so
> NVS only outputs MS, `binning_avg_cov`, and `binning_avg_rho`.

---

## Input file — `input_ACF.xml`

The file must be in the run directory (hard-coded name `input_ACF.xml`).

```xml
<root>
  <WilsonFlow>
    <tau>4</tau>
    <data_name>E</data_name>
    <path>/path/to/flowed/fields/</path>
  </WilsonFlow>

  <Configurations>
    <conf_prefix>ckpoint_lat</conf_prefix>
    <StartConfiguration>100</StartConfiguration>
    <EndConfiguration>500</EndConfiguration>
  </Configurations>

  <Autocorrelations>
    <MDtime_div_fac>1</MDtime_div_fac>
    <MScut>100</MScut>
    <space_block_sizes>2 4</space_block_sizes>
    <R>16</R>
    <isFullTimeAvg>1</isFullTimeAvg>
  </Autocorrelations>
</root>
```

### Parameter reference

#### `[WilsonFlow]`

| Parameter | Type | Description |
|---|---|---|
| `tau` | `int` | Wilson-flow step τ_W (integer index into flowed-field filenames) |
| `data_name` | `string` | Observable name; used in filename construction and output tags. Typical values: `E` (energy density), `topo5li` (5-link topological charge density), `plaq` (per-site plaquette at τ=0) |
| `path` | `string` | Directory containing flowed-field files |

**File naming convention** (constructed inside the executable):

```
{path}{data_name}_{tau}_{conf_prefix}.{conf_number}
```

Example: `/data/flowed/E_4_ckpoint_lat.200`

#### `[Configurations]`

| Parameter | Type | Description |
|---|---|---|
| `conf_prefix` | `string` | Filename prefix for gauge configs (e.g. `ckpoint_lat`) |
| `StartConfiguration` | `int` | First config index (inclusive) |
| `EndConfiguration` | `int` | Last config index (inclusive) |

Total configs loaded: `N = EndConfiguration − StartConfiguration + 1`.

#### `[Autocorrelations]`

| Parameter | Type | Description |
|---|---|---|
| `MDtime_div_fac` | `int` | Number of MD-chain bins. `1` = single chain (no binning). `n>1` = split chain into `n` bins of `T = N/n` configs each |
| `MScut` (`W`) | `int` | Lag cutoff for Madras–Sokal error estimate. Only used when running MS mode; set to `100` (dummy) for MF mode |
| `space_block_sizes` | `vector<int>` | List of spatial pre-blocking sizes `l_B` to loop over (e.g. `2 4 8`). Each produces independent output |
| `R` | `int` | Summation radius for Master-Field error estimate (in lattice units). Signals the error mode: `R ≥ 0` → MF mode for VS; `R < 0` → MS mode for VS; see table below |
| `isFullTimeAvg` | `int` | `1` = use all source times (full-time average over `i`); `0` = fix source time at `i=0` only. Full-time average gives better statistics but is more expensive |

---

## Error estimation modes

The mode is determined by the combination of `MDtime_div_fac` and `R`:

| Mode | Condition | VS output | NVS output | Notes |
|---|---|---|---|---|
| **Master-Field (MF)** | `MDtime_div_fac=1`, `R ≥ 0` | `auto2_VS_MF.dat` | — | VS only; uses spatial covariance sum out to radius `R` (in units of `l_B` blocks on the coarsened lattice) |
| **Madras–Sokal (MS)** | `MDtime_div_fac=1`, `R < 0` | `auto2_VS_MS.dat` | `auto2_NVS_MS.dat` | Uses four-point variance with inner cutoff `W = MScut` (≥ 100). **Cross-check only** — combining per-site MS into a spatial-average error needs an unverifiable site-independence assumption; see [Method limitations](#method-limitations) |
| **MD-chain binning** | `MDtime_div_fac > 1` | `auto2_VS.dat` | `auto2_NVS.dat` | Splits chain into `n_bin` segments; error from sample variance across bins |

> **VS MF mode** (`R ≥ 0`, `MDtime_div_fac=1`): the assertion `R < 0 || MScut ≥ 100`
> in the source means you must set `MScut` to at least 100 when using MF mode
> (it is unused but must pass the guard).

> **NVS MS mode** (`MDtime_div_fac=1`): NVS has no MF implementation; it always
> falls back to the local Madras–Sokal approximation (sparse spatial input only;
> retained as a per-site cross-check rather than as the primary error estimator,
> since combining per-site results into a spatial-average error requires an
> unverifiable site-independence assumption).

### Choosing `R` for MF mode

`R` is the summation radius in lattice units.  Inside `MF_approx`, it is
converted to blocks: `R_b = R / block_size`.  The variance accumulation loop
runs over all lattice shells `r = 1 … R_b` (Chebyshev ball in 4D: sites with
`max(|Δx|) = r`).  Choose `R` large enough that `σ_ρ(t)` has saturated — see
the `σ_ρ vs R/l_B` plot in `fthmc_utils.plot_sigma_vs_R`.

**Typical values**: `R = 16` for a 32⁴ lattice with `l_B = 2` gives `R_b = 8`
blocks, which usually covers the saturation plateau.

> **Block size `l_B` in MF is for data compression / variance reduction only**
> (cf. Bruno 2023). It does not enter the asymptotic σ_ρ formula in any
> essential way — the error estimate comes from the exponential falloff of
> spatial correlations, integrated over a 4D ball of radius R. Different `l_B`
> values should give the same σ_ρ once the saturation plateau is reached.

---

## Running

Both executables are standard Grid programs.  Pass the lattice geometry and
MPI/SIMD layout via standard Grid flags:

```bash
./autocovariance_VS --grid 32.32.32.32 --mpi 1.1.1.1
./autocovariance_NVS --grid 32.32.32.32 --mpi 1.1.1.1
```

Optional flag:
```bash
--debug    # prints per-config spatial sum for each loaded field (sanity check)
```

The executable reads `input_ACF.xml` from the current working directory and
writes all output to **stdout** via `GridLogMessage`.

---

## Output: stdout log format

All results are printed as tagged lines in stdout.  The relevant tags are:

| Tag | Function | Columns (after tag) |
|---|---|---|
| `{ACC_type} MFACC {name} (Master-Field Approx):` | ACF ratio | `tau  block_size  R  t  ρ(t)` |
| `{ACC_type} MFACC {name} (binning_avg_cov):`     | ACF ratio | `tau  block_size  n_bin  t  ρ(t)` |
| `{ACC_type} MFACC {name} (binning_avg_rho):`     | ACF ratio | `tau  block_size  n_bin  t  ρ(t)` |
| `{ACC_type} Variance {name} (Master-Field Approx):` | MF variance | `tau  block_size  r  t  G(t)  var  cov` |
| `{ACC_type} Variance {name} (binning_avg_cov):`     | binning variance | `tau  block_size  n_bin  t  σ²_full  σ²_no_cov  −2cov_term` |
| `{ACC_type} Variance {name} (binning_avg_rho):`     | binning variance | `tau  block_size  n_bin  t  var  …` |
| `{ACC_type} MFCOV {name} (Master-Field Approx):`    | Raw covariance | `tau  block_size  R  t  G(t)` |
| `{ACC_type} MFACC {name} (Madras-Sokal Approx):`    | MS ACF ratio | `tau  block_size  1  t  ρ(t)` |
| `{ACC_type} Variance {name} (Madras-Sokal Approx):` | MS variance | `tau  block_size  1  t  σ²` |

`{ACC_type}` is `LVS` (VS executable) or `NVS`.
`{name}` is composed from `data_name` and spatial sampling, e.g.
`Blocked E G_cent` (data_name=E, blocked spatial coarsening, centered G).

> **Migration note (2026-04-30):** previous versions of this code used
> `(Binning)` / `(Binning2)` log markers and tagged the centered/connected
> distinction with `ACC` / `ACC2`. The current convention is
> `(binning_avg_cov)` / `(binning_avg_rho)` for the routine, and
> `G_cent` / `G_conn` for the form of G. **Update any post-processing
> scripts that grep for the old markers.**

**Extracting `.dat` files from stdout:**

A post-processing script greps for the relevant tags and reformats the columns
into the compact `auto2_*.dat` layout understood by `fthmc_utils.py`.  The
`tau_kind` column (see below) encodes the observable × sampling × estimator
combination; it is assigned during parsing.

---

## Output data files

After parsing stdout, the following files are written to the run directory:

| File | Estimator | Columns | When produced |
|---|---|---|---|
| `auto2_VS_MF.dat`  | VS, Master-Field   | 9 | `MDtime_div_fac=1`, `R ≥ 0` |
| `auto2_VS_MS.dat`  | VS, Madras–Sokal   | 7 | `MDtime_div_fac=1`, `R < 0` |
| `auto2_VS.dat`     | VS, MD-chain bin   | 9 | `MDtime_div_fac > 1` |
| `auto2_NVS_MS.dat` | NVS, Madras–Sokal  | 7 | `MDtime_div_fac=1` |
| `auto2_NVS.dat`    | NVS, MD-chain bin  | 9 | `MDtime_div_fac > 1` |

### Column layout — `auto2_VS_MF.dat` (9 columns)

```
tau_kind  bin_method  tau  block_size  R  t  G(t)  var  cov
```

| Column | Name | Description |
|---|---|---|
| 0 | `tau_kind` | Observable × sampling × G-form index (see encoding below) |
| 1 | `bin_method` | `0` = MFACC rows (ACF ratio at max R, `var=−1`, `cov=−1`); `1` = Variance rows (covariance-of-covariances for each shell radius `r`) |
| 2 | `tau` | Wilson-flow step |
| 3 | `block_size` | Spatial pre-blocking size l_B |
| 4 | `R` | `bin_method=0`: max summation radius; `bin_method=1`: shell radius `r` from 1 to R_max |
| 5 | `t` | MC-time lag |
| 6 | `G(t)` | `bin_method=0`: ACF ratio ρ(t) = G(t)/G(0); `bin_method=1`: raw covariance sum G(t) |
| 7 | `var` | `bin_method=0`: `−1` (not available); `bin_method=1`: variance term for σ_ρ formula |
| 8 | `cov` | `bin_method=0`: `−1` (not available); `bin_method=1`: covariance term for σ_ρ formula |

> **`bin_method=1` rows are the input to `plot_sigma_vs_R`**: `fthmc_utils` reads
> `G(t)`, `var`, `cov` at each shell radius `r` to compute σ_ρ(t) as a function
> of `R/l_B` using the Bruno (2023) covariance-of-covariances formula.

### Column layout — `auto2_VS_MS.dat` / `auto2_NVS_MS.dat` (7 columns)

```
tau_kind  bin_method  tau  block_size  bin_size=1  t  rho(t)
```

Madras–Sokal output: no separate variance column (σ is computed inline and
appended to the same row internally, but the standard `.dat` layout stores only
the ACF ratio here).

### Column layout — `auto2_VS.dat` / `auto2_NVS.dat` (9 columns)

```
tau_kind  bin_method  tau  block_size  n_bin  t  val  −1  flag
```

`flag=0`: `val` is the ACF ratio ρ(t);  `flag=1`: `val` is the variance σ²(t).

> **Watch out — `bin_method` is overloaded across files.** In `auto2_VS_MF.dat`
> it distinguishes MFACC ρ(t) rows (0) from per-shell variance rows (1). In
> `auto2_VS.dat` / `auto2_NVS.dat` it distinguishes which binning routine
> produced the row: `0` = `binning_avg_cov`, `1` = `binning_avg_rho`.
> See [Method limitations](#method-limitations) for why `binning_avg_rho`
> is the preferred reported number.

---

## `tau_kind` encoding

`tau_kind = 4 * i_obs + 2 * i_blk + i_form`

| `i_obs` | `i_blk` | `i_form` | `tau_kind` | Meaning |
|---|---|---|---|---|
| 0 | 0 | 0 | **0** | Energy density, blocked,  G_cent (centered) |
| 0 | 0 | 1 | **1** | Energy density, blocked,  G_conn (connected) |
| 0 | 1 | 0 | **2** | Energy density, sparse,   G_cent |
| 0 | 1 | 1 | **3** | Energy density, sparse,   G_conn |
| 1 | 0 | 0 | **4** | Topo. charge density, blocked,  G_cent |
| 1 | 0 | 1 | **5** | Topo. charge density, blocked,  G_conn |
| 1 | 1 | 0 | **6** | Topo. charge density, sparse,   G_cent |
| 1 | 1 | 1 | **7** | Topo. charge density, sparse,   G_conn |

`G_cent` = centered (subtract per-config means *before* the product).
`G_conn` = connected (take the product, *then* subtract the product of
per-config means). They differ at `O(1/V)`; cross-checking them is a
finite-volume-bias diagnostic.

> **`auto2_VS_MF.dat` only contains `tau_kind ∈ {0, 1}` (or `{4, 5}` for
> `data_name=topo`).** All other `tau_kind` values use sparse sampling, which
> `MF_approx` does not compute. Sparse and `G_conn` sparse entries appear in
> `auto2_VS_MS.dat`, `auto2_VS.dat`, and the NVS variants.

`fthmc_utils.MasterFieldACF.rho()` selects the right `tau_kind` from its
`observable`, `blocked`, and `connected` keyword arguments:
```python
tau_kind = 4 * (0 if observable == "E" else 1) + 2 * (0 if blocked else 1) + (1 if connected else 0)
```

---

## Reading outputs with `fthmc_utils`

`fthmc_utils.MasterFieldACF` (accessed via `HMCRun.master_field`) wraps all
five file types.  Typical usage:

```python
# Select estimator and get (t, ρ(t), σ(t))
t, rho, sigma = run.master_field.rho(
    estimator="MF",     # "MF" | "MS" | "VS" | "NVS" | "NVS_MS"
    tau=4,              # Wilson-flow step τ_W
    block_size=2,       # l_B
    R=16,               # max summation radius (MF mode, bin_method=0 rows)
    bin_method=0,       # 0 = ACF ratio rows; 1 = variance rows
    observable="E",     # "E" or "topo"
    blocked=True,       # True = block-averaged, False = sparse-sampled
    connected=False,    # False = G_cent (centered), True = G_conn (connected)
)

# σ_ρ vs R plot (uses bin_method=1 rows internally)
run.master_field.plot_sigma_vs_R(
    t_fixed=5,
    tau=4, block_size=2,
    estimator="MF",
    lattice_L=32,       # L for 4D volume V = L^4
)

# Full per-tau analysis suite (lim-sup diagnostic + ACF + fit, all τ_W)
run.master_field.mf_analysis(
    [0, 4, 8, 12, 16],
    fit_ranges={4: [20, 40]},
    block_size=2,
    estimator="MF",
)
```

---

## Typical workflow

### MF mode (single long chain, master-field error)

```xml
<Autocorrelations>
  <MDtime_div_fac>1</MDtime_div_fac>
  <MScut>100</MScut>         <!-- dummy; must be ≥ 100 to pass VS assertion -->
  <space_block_sizes>2 4</space_block_sizes>
  <R>16</R>                  <!-- summation radius in lattice units -->
  <isFullTimeAvg>1</isFullTimeAvg>
</Autocorrelations>
```

Produces `auto2_VS_MF.dat`.  Use `plot_sigma_vs_R` to find the saturation
radius R_sat, then use `mf_analysis` to extract τ_exp.

### Binning mode (multiple shorter chains)

```xml
<Autocorrelations>
  <MDtime_div_fac>6</MDtime_div_fac>   <!-- 6 MD-chain bins -->
  <MScut>50</MScut>
  <space_block_sizes>2</space_block_sizes>
  <R>-1</R>                            <!-- negative → no MF; triggers MS assertion bypass -->
  <isFullTimeAvg>1</isFullTimeAvg>
</Autocorrelations>
```

Produces `auto2_VS.dat` and `auto2_NVS.dat`.  Use `estimator="VS"` or `"NVS"`
in `rho()`. **Use `bin_method=1` (= `binning_avg_rho`) as the reported number;
`bin_method=0` (= `binning_avg_cov`) as a cross-check.**

### MS mode (single chain, Madras–Sokal approximation) — diagnostic only

```xml
<Autocorrelations>
  <MDtime_div_fac>1</MDtime_div_fac>
  <MScut>50</MScut>
  <space_block_sizes>2</space_block_sizes>
  <R>-1</R>                            <!-- triggers MS path in VS code -->
  <isFullTimeAvg>1</isFullTimeAvg>
</Autocorrelations>
```

Produces `auto2_VS_MS.dat` and `auto2_NVS_MS.dat`. Currently unreliable —
retained as a placeholder (see below).

---

## Notes

- **Both G_cent (centered) and G_conn (connected) variants** are always
  computed and written; they differ only in the order of mean subtraction
  relative to the product (`G_cent`: subtract-then-multiply;
  `G_conn`: multiply-then-subtract-product-of-means). They agree at `O(1/V)`.
  `connected=False` (default) selects `G_cent`.
- **Sparse vs blocked** (spatial-lattice coarsening): block-averaging
  preserves all long-wavelength information `k < 1/l_B` and is the safer
  default. Sparse sampling at the same `l_B` keeps only `V/l_B^d` sites
  and discards information between sample points, so it can run out of
  statistics if `l_B` exceeds the correlation length. Block-averaged is
  the primary variant; sparse is a cross-check (and the only acceptable
  spatial input for the local-MS path, where block-averaged input would
  require the unknown intra-block cross-site covariance).
- **`isFullTimeAvg=1`** is strongly recommended: it averages the product
  `A(x,i)·A(x,i+t)` over all valid source times `i`, significantly reducing
  statistical noise at large `t`.
- **Field files must exist for every config** in `[StartConfiguration,
  EndConfiguration]`.  Missing files will cause a read error.

### `binning_avg_cov` vs `binning_avg_rho`

Both run when `MDtime_div_fac > 1` and write to the same `auto2_VS.dat` /
`auto2_NVS.dat` file, distinguished by `bin_method`
(0 = `binning_avg_cov`, 1 = `binning_avg_rho`).

| Routine | `bin_method` | Formula | Comment |
|---|---|---|---|
| `binning_avg_cov` | 0 | `ρ(t) = ⟨G(t)⟩_b / ⟨G(0)⟩_b` (average covariance, then ratio) | Variance via delta-method linearization using inter-bin `Cov(G_b(t), G_b(0))` |
| `binning_avg_rho` | 1 | `ρ(t) = ⟨G(t)/G(0)⟩_b` (per-bin ratio, then average) | **Preferred.** Variance computed directly from inter-bin sample variance of ρ_b — no linearization, no covariance estimate needed |

These are two distinct point estimators of the population ACC. They agree
at leading order in the delta-method linearisation — expanding `G_b(0)`
around its bin-mean and dropping `O(1/n_bin)` corrections gives identical
leading-order variances — but are **not** mathematically equivalent at
finite `n_bin`. **`binning_avg_rho` is the recommended reported number**
because:
1. The variance is computed directly from a sample of `ρ_b` values, with no
   reliance on linearisation.
2. Each MD-chain bin is treated as an independent replicate of the ACC
   estimate — the standard interpretation, justified by the LLN when bin
   width > τ_int.
3. Within-bin spatial coarsening (block-avg vs sparse) only changes
   within-bin precision, not the inter-bin variance that supplies the
   reported error. Block-averaging is the safer default — sparse-sampling
   at large `l_B` discards information between sample points and can run
   out of statistics if `l_B` exceeds the correlation length.

`binning_avg_cov` is run alongside as a cross-check on the linearisation;
disagreement between the two at the working `n_bin` is itself an empirical
diagnostic on the validity of the linearisation.

### Method limitations

The four error estimators implemented here split along two independent axes
— spatial input (full lattice / block-averaged / sparse-sampled) and MD-time
partition (single chain / multi-bin) — plus a prior design choice to use the
autocovariance `G_x(t)` rather than the per-site ACC `ρ_x(t) = G_x(t)/G_x(0)`
as the basic statistical variable (the latter has poor signal-to-noise at
single sites because `G_x(0)` is a small, fluctuating denominator). See §4.2
of `Master_Field_Type_Autocorrelation/main.tex` for the full discussion.

- **(1) `MF_approx` (Master-Field, single chain)** — directly estimates the
  spatial covariance density `Cov[G(x,t), G(x+y,t)]` from a single configuration
  (averaging over `x` via translation invariance) and integrates it over
  `|y| ≤ R`. Saturation at `R_sat` is empirically observable provided
  `R_sat < L/2`; if no plateau is reached the error is unbounded. The block
  size `l_B` is a data-compression / variance-reduction parameter only
  (cf. Bruno 2023), not part of the asymptotic formula. Caveats: with
  `l_B < √(8τ_W)` (Wilson-flow smearing radius), adjacent blocks share
  smeared field and `σ_ρ` is biased low. Always compare blocked-MF for
  several `l_B` (`<space_block_sizes>2 4 8</space_block_sizes>`) and use
  `plot_sigma_vs_R` to confirm saturation.

- **(2) `MS_approx` + sparse (Local Madras–Sokal)** — the MS variance
  formula is the standard, theoretically clean estimator both for the
  volume-summed scalar observable and for `G_x(t)` at any single fixed
  site `x`, provided `T ≫ W` (we use `W ≥ 100`, the four-point inner
  cutoff `Λ` of Lüscher 2005 Eq. E.11; this is *not* the τ_int summation
  window of Madras–Sokal/Wolff automatic windowing, which scales with τ_int).
  The complication arises in combining per-site MS estimates into an error
  on the spatial average: this requires a site-independence assumption not
  ensured by the data. Validating it via `σ_ρ ∝ l_B^{d/2}` Bienaymé scaling
  is circular. **Sparse spatial input is the only acceptable choice** here
  — block-averaged MS would require the unknown intra-block cross-site
  covariance and cannot be obtained from a single configuration.
  Sparse-sampled local MS is retained as a cross-check rather than as a
  primary estimator.

- **(3) Block-first per-block binning ("Block-First")** — multi-bin
  estimator that bins each spatial cell's time series and combines per-cell
  errors under an inter-cell independence assumption. Both block-averaged
  and sparse-sampled inputs are technically permissible at the binning step
  (no intra-block covariance is needed), but the inter-cell independence
  assumption suffers from the same Bienaymé-scaling circularity as method (2).
  Validity reduces to verifying `R_sat ≪ L/2` via the same diagnostic that
  underpins MF, making this method redundant with MF in practice.
  *Currently produced by neither `binning_avg_cov` nor `binning_avg_rho` as
  written* — both apply method (4) below.

- **(4) `binning_avg_cov` / `binning_avg_rho` (preferred)** — bin-first
  per-bin spatial average. The chain is partitioned into `n_bin ≥ 2` bins
  of length `L_bin ≫ τ_int`; within each bin the spatially averaged
  autocovariance `G̃(t)|_b` is computed; inter-bin sample variance of
  `G̃(t)|_b` supplies `Var[G̃(t)]`. The only assumption is the LLN on the
  temporal axis (testable by varying `L_bin`); no spatial-independence
  assumption enters at any stage. Within each bin, **block-averaged spatial
  input is the safer default**: sparse-sampling at the same `l_B` discards
  information between sample points and runs out of statistics if `l_B`
  exceeds the correlation length, whereas block-averaging preserves all
  long-wavelength modes `k < 1/l_B`. Requires `n_bin ≳ 20–50` for the
  sample variance to converge; with smaller `n_bin` the variance estimator
  itself carries `O(1/√n_bin)` fractional uncertainty.

For a single long chain: use `MF_approx` with the saturation check.
For multiple short chains: use `binning_avg_rho` (with `binning_avg_cov` as a
linearisation cross-check).
