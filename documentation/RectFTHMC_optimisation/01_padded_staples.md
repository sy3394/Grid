# 1. Padded staples on the depth-2 ghost grid

**What:** the flow staples are computed by single fused kernels on a
padded (halo-materialised) grid instead of Cshift chains.

- `Rect_Stout::RectStapleStencilRs` / `RectStaplePaddedRs`
  (RectSmearing.h): the Rs staple (sum over nu of upper+lower 1x2
  rectangles, 10 link-reads per nu), writing rho*adj(Stap) directly.
- `SmearedConfigurationRect::BaseSmear_ghost_rect` / `BaseSmear_ghost_plq`
  (GaugeConfigurationRect.h): wrappers; plq version transcribes Masked's
  `BaseSmear_ghost` kernel onto the full grid.

**Why depth 2:** rectangle paths reach two hops in nu; the plaquette class
uses depth 1. One `PaddedCell GhostRect(2,...)` member serves staples and
all force-term stencils; `ExchangePeriodic` runs once per level and the
padded field is reused by every kernel of that level.

**Stencils** are built once in the constructor (`gStencils_rectsmear`,
`gStencils_plqsmear`); entry order must match kernel read order — both
sides carry matching comments.

**Validation:** `tests/forces/Test_rect_staple_padded` vs
`RectStapleUnoptimisedRs`, last-bit agreement, isotropic + anisotropic
volumes, CPU and GPU, 1-2 ranks.
