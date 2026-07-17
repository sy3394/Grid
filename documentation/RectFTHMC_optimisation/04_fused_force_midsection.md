# 4. Fused force midsection (one kernel)

**What:** the sitewise midsection of the default `logDetJacobianForceLevel`
is a single `accelerator_for` ("J_Mab_Inv_dJdX_fusedOpt"):

1. ONE Horner recursion yields BOTH dJdX_b and J = 1 + t2/2 (the XB
   scheme sketched in the old routine; algebraically the same truncation
   sum_{k<=11} X^k/(k+1)! as the reference Taylor loop — FP order differs
   at the last bit only).
2. Mab = 1 - J Nxx in registers.
3. Its real-part 8x8 LU inverse IN-KERNEL via `Inverse_RealPartSite`
   (see 05) — the adjoint rep is real (production choice, as Masked).
4. dJdXe traces (factor -0.5, this class's normalisation) and
   MpInvJx = -MpAdInv J.

Jx, Mab, MpAdInv, nMpInv never exist as lattice fields. Preceded by
`make_adjoint_rep` (ZxAd) and `ComputeNxy` (NxxAd), which stay separate
(ZxAd is also a kernel input; NxxAd is reused later in the level).

**Fusion boundary history:** the previous boundary was the 8x8 inversion —
the Eigen-based `Inverse` is host-only (device round trip per level). The
site-local LU dissolved it. On CPU builds only the LU section serialises
SIMD lanes (data-dependent pivoting); the rest stays vectorised.

**Validation:** per-level force vs `int old` (~1e-27 rel^2), FD level
test (`Test_rect_fd_level0`), all platforms/ranks as in 10.
