# 2. Fused adjoint kernels

**What:** the per-generator lattice loops of `ComputeNxy` and
`Compute_MpInvJx_dNxxdSy` are replaced by single sitewise kernels using
the hand-expanded `SU3::LieAlgebraProject` overloads (GaugeGroup.h):

- `ComputeNxy(PlaqL,PlaqR,NxAd)`: two-argument overload — the whole 8x8
  N matrix per site in one pass (the T'^b = 2i t^b convention is internal).
- `Compute_MpInvJx_dNxxdSy(...)`: `accelerator_for2d` over (sites x
  generators); UtaU = adj(PlaqL) (i t^a) PlaqR, one-argument
  LieAlgebraProject (internal T'^c = 2i t^c), traceProduct with MpInvJx.

**Convention ledger (this class vs Masked):** ta = i t^a here vs 2i t^a;
dJdX -0.5 vs -1.0; final force scale -1.0 vs -0.5. The factors cancel to
the same net force on both the Fdet1 and Fdet2 chains — verified
numerically by `Test_rect_vs_masked_plq` (forces equal to ~1e-27 rel^2).

**Old versions** kept as `ComputeNxy(int old,...)` /
`Compute_MpInvJx_dNxxdSy(int old,...)` for the checks; to be deleted.
