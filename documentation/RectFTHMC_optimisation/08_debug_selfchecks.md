# 8. DEBUG self-checks (single-HMC-run validation)

**Enable:** flip `#undef DEBUG` -> `#define DEBUG` at the top of
GaugeConfigurationRect.h and rebuild (a header define — everything
including it recompiles; keep a separate build tree for the DEBUG binary).

Every force/action evaluation then prints, in Masked.h style:
- `DEBUG: BaseSmear_ghost`  padded staple vs unpadded reference
- `DEBUG: ZxAd`             make_adjoint_rep vs generator loop
- `DEBUG: NxxAd`            fused ComputeNxy vs ComputeNxy(old,...)
- `DEBUG: fused J/Mab/Inv/dJdX`  the mega-kernel vs the old lattice
  sequence (Taylor J, complex Eigen inverse, lattice Horner) — compared
  post-mask; expect ~1e-25 absolute vs O(1e3) refs (Horner-J + LU-vs-Eigen
  fingerprint)
- `DEBUG: Plaq L/R line NNNN`    per-term norms (line-tagged)
- `DEBUG: forceLevel` / `logDetJacobianLevel`  full level vs (old,...)
- `DEBUG: ... TOTAL`             totals vs the old top-levels

Triage one-liner (largest diffs first):

    grep ' DEBUG: ' traj.log | awk '{for(i=1;i<NF;i++) if($i=="diff") print $(i+1), $0}' | sort -g -r | head -20

**Cost:** ~3x per evaluation, and the old path's Eigen inverse runs on the
HOST — at production volumes keep the DEBUG trajectory short (few MD
steps). The MD evolution itself is bit-identical to the non-DEBUG build.

Aurora GPU status (2026-07-14): lndet self-checks observed at diff
0 .. 1e-13 in a production-binary run.
