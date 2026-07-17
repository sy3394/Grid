# 5. Inverse_RealPartSite and the Inverse_RealPart CPU bug

**Bug found (Lattice_trace.h):** `toReal` duplicates each complex lane
into TWO adjacent real lanes; the CPU (non-SIMT) branch of
`Inverse_RealPart` read real lane `blane` where complex lane `blane`
lives at `2*blane` — every lane pair inverted the wrong site's matrix and
the upper half was never inverted (residual O(1) on real input). The GPU
(GRID_SIMT) path was unaffected — production GPU results were fine; this
is likely why an earlier attempt to use Inverse_RealPart in the rect code
on a CPU box failed and was reverted to the Eigen `Inverse`.

**Fix:** new `accelerator_inline Inverse_RealPartSite<T,N>` — site-local
real-part LU inverse usable inside coalesced kernels. Under GRID_SIMT the
kernel values are per-thread scalars (no lane handling); on CPU it loops
lanes internally reading `real()` of the complex lane directly (which
structurally avoids the toReal trap). `Inverse_RealPart` is reimplemented
on top of it, deleting the hand-rolled lane code.

**Validation:** `tests/forces/Test_inverse_realpart` — vs the complex
Eigen inverse on exactly real input (residual ~1e-27, rel 2e-32). PASSED
on Aurora PVC 2026-07-14: first GRID_SIMT exercise of the site helper.
