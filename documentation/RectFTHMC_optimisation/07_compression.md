# 7. Compression to the active half (checkerboarding replacement)

**Problem:** the type-2 (brick) mask is red-black in units of 1x2x2x2
bricks — NOT a fine-lattice parity class — so `GridRedBlackCartesian` /
`pickCheckerboard` cannot represent it.

**Design:** geometry-free compression. Every compressed operation is
sitewise (staples, force terms and all Cshifts run on the FULL grid), so
the container needs no geometric metadata: a plain `GridCartesian` with
one extent halved (`HalfGrid`) is pure storage. Per-level tables
(`active_tab[smr]`) map half-grid osites to active full-grid osites,
built in the constructor FROM THE MASKS THEMSELVES — one code path for
both mask types. `pickActive`/`setActive` are whole-vector-word indexed
copies (no lane surgery, unlike pickCheckerboard). The zero-filled
scatter REPLACES the explicit "mask it off": dJdXe was masked anyway, and
inactive-site MpInvJx values only ever multiply the vanishing masked link
(each force term carries exactly one masked-link factor).

**Requirement + fallback:** whole-word copies need the mask uniform
within each SIMD vector word (lane strides = 0 mod 4 perpendicular to
mu). The constructor CHECKS this empirically against the actual masks;
non-uniform layouts (e.g. 8x4x12x16 with GEN simd, or any local extent
giving stride 2) automatically fall back to the validated full-grid path.
A log line reports ENABLED/DISABLED — both are correct.

**Gains where enabled:** the fused force kernel (incl. the LU inverse)
and the lndet kernel run on half the sites; two big temporaries halve;
same gather/scatter overhead as pick/setCheckerboard.

**Validation:** identical results on ENABLED (8.8.8.8) and DISABLED
(8.4.12.16) layouts, CPU + Aurora GPU, 1/2/4 ranks.
