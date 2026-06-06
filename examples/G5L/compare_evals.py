#!/usr/bin/env python3
"""
Compare gamma5-Block Lanczos D_W eigenvalues against a reference eval file.

Reference format (evals_DW_*.dat):   <MD_time>  Re(lambda)  Im(lambda)
Computed  format (g5bl_evals*.dat):   <tag>      Re(lambda)  Im(lambda)

The Wilson bare mass enters D_W only on the diagonal: D_W(m) = D_W(0) + m*I.
So computed eigenvalues differ from the reference by a CONSTANT REAL SHIFT
delta = (m_used - m_ref).  The imaginary parts are mass-independent and must
match exactly.  We therefore:
  1. select reference rows at a chosen MD time (--mdtime, default: the first),
  2. match each computed eigenvalue to its nearest reference eigenvalue in Im,
  3. estimate delta = median(Re_ref - Re_comp) over matched pairs,
  4. report RMS of (Re_comp + delta - Re_ref) and of (Im_comp - Im_ref).

Usage:
  compare_evals.py REFERENCE.dat COMPUTED.dat [--mdtime T] [--n N]
"""
import sys, argparse
import numpy as np

def load(path):
    rows = []
    with open(path) as f:
        for ln in f:
            p = ln.split()
            if len(p) < 3:
                continue
            try:
                rows.append((float(p[0]), float(p[1]), float(p[2])))
            except ValueError:
                continue
    return np.array(rows)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("reference")
    ap.add_argument("computed")
    ap.add_argument("--mdtime", type=float, default=None,
                    help="MD-time slice (col 1) of the reference; if omitted, scan ALL slices "
                         "and pick the best-matching one")
    ap.add_argument("--n", type=int, default=0, help="limit to N lowest-|Im| computed evals")
    args = ap.parse_args()

    ref = load(args.reference)
    com = load(args.computed)
    if ref.size == 0 or com.size == 0:
        print("empty input"); sys.exit(1)

    C = com[:, 1:3]
    if args.n > 0:
        order = np.argsort(np.abs(C[:, 1]))
        C = C[order][:args.n]

    times = np.unique(ref[:, 0])

    def score_slice(t):
        """Return (rms_combined, delta) for reference slice at MD time t."""
        R = ref[np.isclose(ref[:, 0], t)][:, 1:3]
        if len(R) == 0:
            return (1e30, 0.0, R)
        # match each computed eval to nearest reference by Im, estimate real shift
        dl = []
        for c in C:
            j = np.argmin(np.abs(R[:, 1] - c[1]))
            dl.append(R[j, 0] - c[0])
        delta = float(np.median(dl))
        errs = []
        for c in C:
            cs = np.array([c[0] + delta, c[1]])
            j = np.argmin(np.abs(R[:, 0] - cs[0]) + np.abs(R[:, 1] - cs[1]))
            errs.append((cs - R[j]))
        errs = np.array(errs)
        rms = np.sqrt(np.mean(errs[:, 0]**2 + errs[:, 1]**2))
        return (rms, delta, R)

    if args.mdtime is not None:
        t = times[np.argmin(np.abs(times - args.mdtime))]
    else:
        # scan all slices, pick min RMS
        best = min(times, key=lambda tt: score_slice(tt)[0])
        t = best
        print("scanned %d MD-time slices; best match at MD time = %g" % (len(times), t))

    rms0, delta, R = score_slice(t)
    print(f"reference MD time = {t}   ({len(R)} ref evals, {len(C)} computed evals)  "
          f"combined RMS = {rms0:.3e}")

    # match each computed eval to nearest reference by Im, then estimate delta
    matched = []
    for c in C:
        j = np.argmin(np.abs(R[:, 1] - c[1]) + 0.0 * np.abs(R[:, 0] - c[0]))
        matched.append((c, R[j]))
    deltas = [r[0] - c[0] for c, r in matched]
    delta = float(np.median(deltas))

    # refine matching using the estimated real shift, then report
    re_err, im_err = [], []
    print(f"\n estimated real shift delta = m_used - m_ref = {delta:+.6f}\n")
    print(f"{'Re_comp+delta':>14} {'Re_ref':>12} {'Im_comp':>12} {'Im_ref':>12} "
          f"{'dRe':>10} {'dIm':>10}")
    for c, _ in matched:
        cs = np.array([c[0] + delta, c[1]])
        j = np.argmin(np.abs(R[:, 0] - cs[0]) + np.abs(R[:, 1] - cs[1]))
        r = R[j]
        dre, dim = cs[0] - r[0], cs[1] - r[1]
        re_err.append(dre); im_err.append(dim)
        print(f"{cs[0]:14.6f} {r[0]:12.6f} {c[1]:12.6f} {r[1]:12.6f} "
              f"{dre:10.2e} {dim:10.2e}")

    re_err = np.array(re_err); im_err = np.array(im_err)
    print(f"\n RMS dRe = {np.sqrt(np.mean(re_err**2)):.3e}   "
          f"RMS dIm = {np.sqrt(np.mean(im_err**2)):.3e}   "
          f"max|dRe| = {np.max(np.abs(re_err)):.3e}   "
          f"max|dIm| = {np.max(np.abs(im_err)):.3e}")

if __name__ == "__main__":
    main()
