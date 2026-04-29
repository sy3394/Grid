# Autocovariance Executables — Usage Guide

Sources: `autocovariance_VS.cc`, `autocovariance_NVS.cc`, `ACC.hpp`

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

Both executables also compute a second estimator `G2` using the alternative
mean-subtraction convention (see source comments); its output is tagged `ACC2`
vs the primary `ACC`.

### Three orthogonal axes

Every output row is determined by three independent choices:

| Axis | Choices | Notes |
|---|---|---|
| **Field construction** | `G` (centered, ACC) vs `G2` (connected, ACC2) | how `G(x,t)` is built from raw `A(x,i)`. `G` subtracts means before the product; `G2` subtracts the product of means after. They differ at `O(1/V)` only |
| **Spatial coarsening** | `blocked` (cell = mean of `l_B^d` sites) vs `sparse` (cell = corner site) | how `G(x,t)` is reduced to a `V/l_B^d`-cell coarse lattice |
| **Error estimation** | `MS_approx`, `MF_approx`, `binning`, `binning2` | which variance formula consumes the coarsened field |

Field construction and spatial coarsening are orthogonal preprocessing steps;
the error-estimation routine is selected at runtime by `MDtime_div_fac` and `R`.

### Combinations actually present in the output

Not every combination is computed — the C++ deliberately skips ones that have
no diagnostic value:

|  | MS | MF | binning | binning2 |
|---|---|---|---|---|
| cen × blocked  (`G_B`)  | ✓ VS, NVS | ✓ VS only | ✓ VS, NVS | ✓ VS, NVS |
| cen × sparse   (`G_s`)  | ✓ VS, NVS | — | ✓ VS, NVS | ✓ VS, NVS |
| conn × blocked (`G2_B`) | ✓ VS, NVS | ✓ VS only | ✓ VS, NVS | ✓ VS, NVS |
| conn × sparse  (`G2_s`) | ✓ VS, NVS | — | ✓ VS, NVS | ✓ VS, NVS |

> **MF + sparse is intentionally omitted.** `MF_approx` integrates the spatial
> covariance density `Cov[G(x,t), G(y,t)]` over all separations `|x−y| ≤ R`.
> Sparse sampling discards every separation that isn't a multiple of `l_B` —
> exactly the sub-`l_B` structure MF is designed to use. Sparse-MF is strictly
> noisier than blocked-MF with no diagnostic upside.

> **NVS has no MF.** Subtracting a single ensemble-mean scalar from every site
> breaks the spatial covariance decomposition that `MF_approx` relies on, so
> NVS only outputs MS, binning, and binning2.

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
| `data_name` | `string` | Observable name; used in filename construction and output tags. Typical values: `E` (energy density), `topo5li` (5-link topological charge density) |
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
| `MDtime_div_fac` | `int` | Number of MD-time bins. `1` = single chain (no binning). `n>1` = split chain into `n` bins of `T = N/n` configs each |
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
| **Madras–Sokal (MS)** | `MDtime_div_fac=1`, `R < 0` | `auto2_VS_MS.dat` | `auto2_NVS_MS.dat` | Uses ACF truncation at lag `W = MScut`; valid when `T ≫ W` |
| **MD-time binning** | `MDtime_div_fac > 1` | `auto2_VS.dat` | `auto2_NVS.dat` | Splits chain into `n_bin` segments; error from sample variance across bins |

> **VS MF mode** (`R ≥ 0`, `MDtime_div_fac=1`): the assertion `R < 0 || MScut ≥ 100`
> in the source means you must set `MScut` to at least 100 when using MF mode
> (it is unused but must pass the guard).

> **NVS MS mode** (`MDtime_div_fac=1`): NVS has no MF implementation; it always
> falls back to the Madras–Sokal approximation.

### Choosing `R` for MF mode

`R` is the summation radius in lattice units.  Inside `MF_approx`, it is
converted to blocks: `R_b = R / block_size`.  The variance accumulation loop
runs over all lattice shells `r = 1 … R_b` (Chebyshev ball in 4D: sites with
`max(|Δx|) = r`).  Choose `R` large enough that `σ_ρ(t)` has saturated — see
the `σ_ρ vs R/l_B` plot in `fthmc_utils.plot_sigma_vs_R`.

**Typical values**: `R = 16` for a 32⁴ lattice with `l_B = 2` gives `R_b = 8`
blocks, which usually covers the saturation plateau.

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
| `{ACC_type} MFACC {name} (Binning):` | ACF ratio | `tau  block_size  n_bin  t  ρ(t)` |
| `{ACC_type} MFACC {name} (Binning2):` | ACF ratio | `tau  block_size  n_bin  t  ρ(t)` |
| `{ACC_type} Variance {name} (Master-Field Approx):` | MF variance | `tau  block_size  r  t  G(t)  var  cov` |
| `{ACC_type} Variance {name} (Binning):` | Binning variance | `tau  block_size  n_bin  t  σ²_full  σ²_no_cov  −2cov_term` |
| `{ACC_type} Variance {name} (Binning2):` | Binning2 variance | `tau  block_size  n_bin  t  var  …` |
| `{ACC_type} MFCOV {name} (Master-Field Approx):` | Raw covariance | `tau  block_size  R  t  G(t)` |
| `{ACC_type} MFACC {name} (Madras-Sokal Approx):` | MS ACF ratio | `tau  block_size  1  t  ρ(t)` |
| `{ACC_type} Variance {name} (Madras-Sokal Approx):` | MS variance | `tau  block_size  1  t  σ²` |

`{ACC_type}` is `LVS` (VS executable) or `NVS`.
`{name}` is composed from `data_name` and spatial sampling, e.g. `Blocked E ACC`.

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
| `auto2_VS_MF.dat` | VS, Master-Field | 9 | `MDtime_div_fac=1`, `R ≥ 0` |
| `auto2_VS_MS.dat` | VS, Madras–Sokal | 7 | `MDtime_div_fac=1`, `R < 0` |
| `auto2_VS.dat` | VS, binning | 9 | `MDtime_div_fac > 1` |
| `auto2_NVS_MS.dat` | NVS, Madras–Sokal | 7 | `MDtime_div_fac=1` |
| `auto2_NVS.dat` | NVS, binning | 9 | `MDtime_div_fac > 1` |

### Column layout — `auto2_VS_MF.dat` (9 columns)

```
tau_kind  bin_method  tau  block_size  R  t  G(t)  var  cov
```

| Column | Name | Description |
|---|---|---|
| 0 | `tau_kind` | Observable × sampling × estimator index (see encoding below) |
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
> `auto2_VS.dat` / `auto2_NVS.dat` it distinguishes the binning routine used:
> `0` = `binning` (ratio-of-averages), `1` = `binning2` (average-of-ratios).
> See [Notes](#notes) for the difference.

---

## `tau_kind` encoding

`tau_kind = 4 * i_obs + 2 * i_blk + i_acc`

| `i_obs` | `i_blk` | `i_acc` | `tau_kind` | Meaning |
|---|---|---|---|---|
| 0 | 0 | 0 | **0** | Energy density, Block-averaged, ACC (G, centered) |
| 0 | 0 | 1 | **1** | Energy density, Block-averaged, ACC2 (G2, connected) |
| 0 | 1 | 0 | **2** | Energy density, Sparse-sampled, ACC (G, centered) |
| 0 | 1 | 1 | **3** | Energy density, Sparse-sampled, ACC2 (G2, connected) |
| 1 | 0 | 0 | **4** | Topo. charge density, Block-averaged, ACC |
| 1 | 0 | 1 | **5** | Topo. charge density, Block-averaged, ACC2 |
| 1 | 1 | 0 | **6** | Topo. charge density, Sparse-sampled, ACC |
| 1 | 1 | 1 | **7** | Topo. charge density, Sparse-sampled, ACC2 |

`ACC` uses estimator `G` (centered: subtract per-config means *before* the product).
`ACC2` uses estimator `G2` (connected: take the product, *then* subtract the
product of per-config means). They differ at `O(1/V)`; cross-checking them is a
finite-volume-bias diagnostic.

> **`auto2_VS_MF.dat` only contains `tau_kind ∈ {0, 1}` (or `{4, 5}` for
> `data_name=topo`).** All other `tau_kind` values use sparse sampling, which
> `MF_approx` does not compute. Sparse and ACC2 entries appear in
> `auto2_VS_MS.dat`, `auto2_VS.dat`, and the NVS variants.

`fthmc_utils.MasterFieldACF.rho()` selects the right `tau_kind` from its
`observable`, `blocked`, and `acc2` keyword arguments:
```python
tau_kind = 4 * (0 if observable == "E" else 1) + 2 * (0 if blocked else 1) + (1 if acc2 else 0)
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
    acc2=False,         # False = ACC (G), True = ACC2 (G2)
)

# σ_ρ vs R plot (uses bin_method=1 rows internally)
run.master_field.plot_sigma_vs_R(
    t_fixed=5,
    tau=4, block_size=2,
    estimator="MF",
    lattice_L=32,       # L for 4D volume V = L^4
)

# ACF plot
run.master_field.plot_acf(tau=4, block_size=2, estimator="MF")

# τ_exp lim-sup diagnostic
run.master_field.plot_tau_diagnostic(tau=4, block_size=2, estimator="MF")
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
radius R_sat, then use `plot_acf` + `plot_tau_diagnostic` to extract τ_exp.

### Binning mode (multiple shorter chains)

```xml
<Autocorrelations>
  <MDtime_div_fac>6</MDtime_div_fac>   <!-- 6 MD-time bins -->
  <MScut>50</MScut>
  <space_block_sizes>2</space_block_sizes>
  <R>-1</R>                            <!-- negative → no MF; triggers MS assertion bypass -->
  <isFullTimeAvg>1</isFullTimeAvg>
</Autocorrelations>
```

Produces `auto2_VS.dat` and `auto2_NVS.dat`.  Use `estimator="VS"` or `"NVS"`
in `rho()`.

### MS mode (single chain, Madras–Sokal approximation)

```xml
<Autocorrelations>
  <MDtime_div_fac>1</MDtime_div_fac>
  <MScut>50</MScut>
  <space_block_sizes>2</space_block_sizes>
  <R>-1</R>                            <!-- triggers MS path in VS code -->
  <isFullTimeAvg>1</isFullTimeAvg>
</Autocorrelations>
```

Produces `auto2_VS_MS.dat` and `auto2_NVS_MS.dat`.

---

## Notes

- **Both G (centered) and G2 (connected) estimators** are always computed and
  written; they differ only in the order of mean subtraction relative to the
  product (`G`: subtract-then-multiply; `G2`: multiply-then-subtract-product-of-means).
  They agree at `O(1/V)`. `acc2=False` (default) selects `G`.
- **Sparse vs blocked**: block-averaging (`blocked=True`) reduces per-cell
  variance by `√(l_B^d)` at the cost of mixing spatial scales within a block.
  Sparse sampling preserves resolution but is noisier. Block-averaged is the
  primary variant; sparse is a cross-check (and the natural choice for `MS`,
  which wants per-site time series).
- **`isFullTimeAvg=1`** is strongly recommended: it averages the product
  `A(x,i)·A(x,i+t)` over all valid source times `i`, significantly reducing
  statistical noise at large `t`.
- **Field files must exist for every config** in `[StartConfiguration,
  EndConfiguration]`.  Missing files will cause a read error.

### `binning` vs `binning2`

Both run when `MDtime_div_fac > 1` and write to the same `auto2_VS.dat` /
`auto2_NVS.dat` file, distinguished by `bin_method` (0 = `binning`,
1 = `binning2`).

| Routine | `bin_method` | Formula | Comment |
|---|---|---|---|
| `binning`  | 0 | `ρ(t) = ⟨G(t)⟩_b / ⟨G(0)⟩_b` (ratio of averages) | typically tighter σ — averages more statistics into the denominator |
| `binning2` | 1 | `ρ(t) = ⟨ G(t)/G(0) ⟩_b` (average of ratios) | needs `bin_size > ~30` for CLT-Gaussianity per bin (see source comment) |

Both are unbiased to leading order; they differ at `O(1/n_bin)`.

### Limitations of each error-estimation method

- **`MS_approx`** — assumes `T ≫ W`; the variance formula for ρ(t) at lag t
  uses ρ(t±k) up to k = W, so the noisy tail at large lag dominates the σ
  estimate. The choice of W is data-dependent (you'd want W ~ a few τ_int but
  τ_int is what you're estimating). The classic 1/√N error scaling itself
  implicitly assumes statistical independence — a circularity discussed in
  §4.3 of `Master_Field_Type_Autocorrelation/main.tex`.
- **`MF_approx`** — requires direct observation of σ_ρ vs R/l_B saturation;
  if no plateau is reached within `R < L/2`, the error is unbounded. With
  `l_B < √(8τ_W)` (Wilson-flow smearing radius), adjacent blocks share
  smeared field, biasing σ_ρ low. Always compare blocked-MF for several
  `l_B` values (`space_block_sizes>2 4 8</space_block_sizes>`) and use
  `plot_sigma_vs_R` to confirm saturation.
- **`binning`** — robust if `T = N/n_bin ≫ τ_int`. With `n_bin` small (~6–10)
  and a chain that's barely long enough, the bin-to-bin variance estimator
  itself has large fractional error (~1/√n_bin).
- **`binning2`** — same caveat as `binning` plus the CLT requirement on
  `bin_size`. Tends to give slightly larger σ.

For a single long chain, use `MF_approx` with the saturation check. For
multiple short chains use `binning`. Cross-check with `binning2` and the
`MS` estimator when in doubt.
