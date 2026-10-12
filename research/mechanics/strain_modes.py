"""Guest-Hutchinson style analysis: infinitesimal mechanisms with a deformable periodic lattice.
Which macroscopic strains can the contact framework (and variants) produce at first order?"""
import numpy as np
from pframe import *
from contact_nd import contact_framework, TSTAR

np.set_printoptions(precision=3, suppress=True)

for n in (2, 3, 4):
    for reps in [(1,) * n, (2,) * n]:
        if n == 4 and reps[0] == 2:
            continue
        A, X, B = contact_framework(n)
        A2, X2, B2 = supercell(A, X, B, reps)
        nm, sr, mats, ss = flexible_analysis(A2, X2, B2)
        print(f"n={n} supercell {reps}: mechanisms(non-trivial)={nm}, strain-space dim={sr} of {n*(n+1)//2}, self-stress={ss}")
        diag_part = max([np.abs(np.diag(E)).max() for E in mats] + [0])
        print("   max |normal-strain component| in achievable strains:", round(diag_part, 12))
