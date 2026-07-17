# 6. Fused lndet (action side)

**What:** the default `logDetJacobianLevel`:
- Ncb via the fused two-argument `LieAlgebraProject` with PlaqL = 1
  (reproduces Nb = 2 Ta(i T^b U C^dag) exactly — resolves the FIXME the
  old routine carried);
- Zac via `SU3Adjoint::make_adjoint_rep`;
- J Taylor + Mab + `Determinant` + `log` in ONE sitewise kernel (the
  pattern of the completed plaquette version at
  GaugeConfigurationMasked.h:967);
- staple from the padded cell for both kernels;
- with compression (07): kernel runs on the active half and the sum needs
  no mask.

**Absolute validation (not just old-vs-new):**
`tests/forces/Test_rect_numjac` builds the true Jacobian of the
production level-0 map by central differences (the masking makes the
active-active block site-diagonal — staples of active links contain only
inactive same-direction links) and compares log dets. Result: ratio
1 + 3e-8 (rect), 1 + 7e-8 (plq) — the lndet is correct absolutely for
both kernels. This is what pinned the factor-2 (09) to the force side.
