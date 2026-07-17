# 9. Force normalisation factor 2 — PARKED

**Finding:** at the production scale `force = -1.0*(Fdet1+Fdet2)` (this
class; net-equivalent to Masked's -0.5 with its factor conventions), the
Jacobian force is exactly 2x dS/dU for BOTH kernels:

- `Test_rect_fd_level0`: eps->0 limit of dS/dSpred = 1/2 at level grain
  and for the totals, plq AND rect (the historical impression that plq was
  consistent came from a curvature coincidence at eps=0.01; Test_rfthmc
  never asserted or swept eps).
- `Test_rect_numjac`: the lndet is correct ABSOLUTELY (vs a numerical
  Jacobian of the actual map), so the 1/2 belongs to the force.
- Direct FD on `SmearedConfigurationMasked` (in Test_rect_fd_level0):
  same 1/2 — the production plq FTHMC force carries the same factor.
  Equivalent fix there: -0.5 -> -0.25. Provenance in the original
  first-pass commit vs later convention reshuffles: undecided (would need
  FD on 3badbfc3).

**Consequences:** the ACTION is correct, so all ensembles are exact
(Metropolis absorbs force errors); the cost is acceptance / step size,
i.e. autocorrelation through rejections only.

**Current state: PARKED.** The scale is deliberately kept at -1.0 (both
overloads; marked `FIXME EXTRA-FACTOR-2` in GaugeConfigurationRect.h) so
the force matches production during consistency checking. While parked:
`Test_rect_fd_gate` FAILS by design (prints an explanatory NOTE);
`Test_rect_vs_masked_plq` PASSES (both classes carry the factor).

**To apply the fix:** change -1.0 -> -0.5 at both FIXME blocks. Then the
gate passes (quadratic ratios ~3.9-4.0) and the cross-check reports
|dF|^2/|F|^2 = 1/4 vs an unfixed Masked (it prints a HINT for exactly
this signature).
