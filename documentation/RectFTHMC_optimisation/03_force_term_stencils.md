# 3. Force-term stencils (the nu loop)

**What:** all 16 PlaqL/PlaqR products of `logDetJacobianForceLevel`
(P1-P6 plaquette, R1-R10 rectangle; one product per derivative slot of
the staple w.r.t. a neighbouring link) are built by GeneralLocalStencil
kernels on the depth-2 padded grid — no CovShift/Cshift chains remain in
the default force path.

**Mask handling — key difference from Masked:** the plaquette class
encodes the mask through per-term checkerboard picks; the type-2 (brick)
mask is not a parity class, so that trick is unavailable. Instead every
kernel reads the padded MASKED link (`gUtmp`, exchanged once per level)
at exactly the position where the reference placed `Utmp`; these entries
are marked `msk` in the constructor shift tables.

**Discipline:** the constructor shift lists transcribe the old routine's
CovShift chains term by term; entry order == kernel read order; each term
carries a comment stating its link product. The `(PlaqR,PlaqL)` argument
swaps of the reference (R1,R3,R6,R8,R9,R10) and all rho signs preserved.

**Memory note:** ~76 stencil entries x padded osites x 16 B per (mu,nu)
pair — GB-scale at production local volumes (same trade the plaquette
class makes). Fallbacks (per-level build/free) noted in the project plan.

**Validation:** per-level force checks vs the `int old` routine, both
kernels, isotropic + anisotropic (decisive for offset mix-ups), CPU+GPU,
1/2/4 ranks: `tests/forces/Test_rect_force_level`.
