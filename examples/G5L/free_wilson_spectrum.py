#!/usr/bin/env python3
"""
Analytic free Wilson-Dirac spectrum for validation of the gamma5-Block Lanczos.

Grid convention (r=1), unit gauge:
    D_W(p) = [ m + sum_mu (1 - cos p_mu) ]  +  i sum_mu gamma_mu sin p_mu
Eigenvalues:
    lambda_pm(p) = m + sum_mu (1 - cos p_mu)  ±  i sqrt( sum_mu sin^2 p_mu )
each with multiplicity 6 (2 spin x 3 colour) for the +/- pair.

Spatial directions periodic:      p_j = 2 pi n_j / L
Time direction antiperiodic:      p_t = 2 pi (n_t + 1/2) / L   (boundary {1,1,1,-1})

Usage:  free_wilson_spectrum.py L m [n_lowest_by_absIm]
"""
import sys, numpy as np

L = int(sys.argv[1]) if len(sys.argv) > 1 else 4
m = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
nshow = int(sys.argv[3]) if len(sys.argv) > 3 else 16

evals = []
for nt in range(L):
    pt = 2*np.pi*(nt + 0.5)/L           # antiperiodic time
    for nx in range(L):
        px = 2*np.pi*nx/L
        for ny in range(L):
            py = 2*np.pi*ny/L
            for nz in range(L):
                pz = 2*np.pi*nz/L
                p = np.array([px, py, pz, pt])
                re = m + np.sum(1 - np.cos(p))
                im = np.sqrt(np.sum(np.sin(p)**2))
                evals.append((re,  im))
                evals.append((re, -im))

ev = np.array(evals)  # each row counted x6 (spin*colour); we keep distinct (re,im)
# unique up to rounding
uniq = {}
for re, im in ev:
    key = (round(re, 8), round(im, 8))
    uniq[key] = uniq.get(key, 0) + 6   # multiplicity
items = sorted(uniq.items(), key=lambda kv: (abs(kv[0][1]), kv[0][0]))

print(f"# free Wilson spectrum  L={L}  m={m}")
print(f"# {'Re':>12} {'Im':>14} {'mult':>6}")
for (re, im), mult in items[:nshow]:
    print(f"{re:14.8f} {im:14.8f} {mult:6d}")
