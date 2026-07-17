# 10. Validation drivers, criteria, results

Drivers (tests/forces/), all run both mask types where applicable:

| Driver | Tests | Criterion |
|---|---|---|
| Test_rect_staple_padded | padded Rs staple vs RectStapleUnoptimisedRs | rel 1e-24 (last-bit) |
| Test_rect_force_level | default vs `int old`: force + lndet, per level and totals | force rel^2 1e-24; lndet rel 1e-12 |
| Test_rect_vs_masked_plq | Rect(mask={1}) vs Masked: configs, force, lndet | configs bit-identical; force ~1e-27 |
| Test_rect_fd_gate | dS vs <F,dU>, 3 schedules, eps pair | quadratic ratio in [3,5.5] |
| Test_rect_fd_level0 | FD per level (no chain) + totals + direct Masked | dS/dSpred sweep (diagnostic) |
| Test_rect_numjac | lndet vs numerical Jacobian of the map | \|ratio-1\| < 1e-6 |
| Test_inverse_realpart | Inverse_RealPart(Site) vs Eigen on real input | residual < 1e-20 |

Batch: `systems/Aurora2/rect_validation.pbs` (single node, debug queue,
1/2/4 ranks, split directions chosen so depth-2 halos cross every
dimension; SUMMARY.txt collects verdicts + the compression
ENABLED/DISABLED line).

Results:
- CPU (macOS arm64, GEN SIMD, single rank): all PASS; compression ENABLED
  on 8.8.8.8, fallback on 8.4.12.16 — identical values.
- Aurora PVC (SYCL, GPU), 2026-07-14: all PASS at np=1,2,4 including both
  compression paths; Test_inverse_realpart PASS = first GRID_SIMT
  validation of the in-kernel LU; DEBUG lndet self-checks at 0..1e-13.
- While the factor 2 is parked (09): fd_gate FAILS by design with a NOTE;
  numjac/xcheck/force_level/staple all PASS.
- Multi-rank on the local Mac is impossible (broken spack OpenMPI for
  np>=2 — fails on an MPI hello-world); Aurora covers it.
