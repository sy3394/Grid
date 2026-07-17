# Rect-FTHMC optimisation — overview

Branch `feature/rect_fthmc-optimise`. Port of the plaquette-kernel FTHMC
optimisations (GaugeConfigurationMasked.h) to the rectangle flow kernel
(SmearedConfigurationRect), plus items found on the way.

Files touched: `Grid/qcd/smearing/GaugeConfigurationRect.h`,
`Grid/qcd/smearing/RectSmearing.h`, `Grid/lattice/Lattice_trace.h`
(one sanctioned bug fix), test drivers under `tests/forces/`.

Method structure: the OPTIMISED routines carry the original names and are
the default; the OLD implementations take a leading `int old` argument and
exist only for consistency checks (delete after sign-off). Normalisation
follows this class's Luscher-leaning conventions (ta = i t^a, dJdX -0.5,
final -1.0); they differ from Masked's (2i, -1.0, -0.5) but cancel to the
same net force (verified numerically).

| # | Item | Doc |
|---|------|-----|
| 1 | Padded staples (Rs + plq) on depth-2 ghost grid | 01 |
| 2 | Fused adjoint kernels (ComputeNxy, Compute_MpInvJx_dNxxdSy) | 02 |
| 3 | Force-term stencils (16 PlaqL/R terms, no Cshift) | 03 |
| 4 | Fused force midsection (one kernel incl. 8x8 inverse) | 04 |
| 5 | Inverse_RealPartSite + CPU lane-bug fix | 05 |
| 6 | Fused lndet kernel + absolute validation | 06 |
| 7 | Compression to the active half (mask-driven) | 07 |
| 8 | DEBUG self-checks (single-HMC-run validation) | 08 |
| 9 | Force normalisation factor 2 (PARKED — see FIXME) | 09 |
| 10| Validation drivers, criteria, Aurora results | 10 |

Status 2026-07-14: implementation complete; validated on CPU (1 rank) and
Aurora PVC (GPU, 1/2/4 ranks); the factor-2 force fix is deliberately
PARKED at the production normalisation (see 09) pending consistency
checks against production; `int old` overloads pending deletion.
