"""Theorem 3: maximal hard-sphere packing fraction of the CaB6 structure type.

Model (a = 1): Ca (radius rL) at the cube corners, B (radius rS) at (1/2,1/2,1/2) +- h e_i.
Wyckoff: B at (x,1/2,1/2) with x = 1/2 - h.  gamma = rS/rL.
Claims checked:
  * phi_max = (pi/6)(30 sqrt2 - 41) at gamma = sqrt2 - 1, h = 1 - sqrt2/2 (x = (sqrt2-1)/2);
  * Ca is held by Ca-Ca contact; the Ca-B pair is NOT in contact (gap);
  * phi_max(gamma) = (pi/6)(1 + 6 gamma^3)            for gamma <= sqrt2 - 1,
                   = (pi/6)(sqrt2-1)^3 (6 + gamma^-3)  for gamma >= sqrt2 - 1 (checked up to gamma = 1.5);
  * phi > pi/(3 sqrt2) (fcc) iff ((sqrt2-1)/6)^(1/3) < gamma < (4 + 7 sqrt2)^(-1/3);
  * NaCl at the same gamma: pi(5 - 3 sqrt2)/3 = 0.7931 (so CaB6 is not the densest binary packing there).
Independent check: brute-force maximisation over (h, scale) with ALL pair distances in a 5^3 block.
Run: python cab6_packing.py
"""
import itertools
import numpy as np
import sympy as sp
from scipy.optimize import minimize_scalar

s2 = sp.sqrt(2)
phi_exact = sp.pi / 6 * (30 * s2 - 41)
rS = (s2 - 1) / 2
print("closed form: (4pi/3)(1/8 + 6 rS^3) - (pi/6)(30sqrt2-41) =",
      sp.simplify(sp.Rational(4, 3) * sp.pi * (sp.Rational(1, 8) + 6 * rS ** 3) - phi_exact),
      " value", sp.N(phi_exact, 12))
h = sp.symbols("h", positive=True)
hs = sp.solve(sp.Eq(h * s2, 1 - 2 * h), h)[0]
print("h* =", sp.nsimplify(hs), "=", sp.N(hs), "; Wyckoff x = 1/2 - h* =", sp.nsimplify(sp.Rational(1, 2) - hs))
dCaB = sp.sqrt((sp.Rational(1, 2) - hs) ** 2 + sp.Rational(1, 2))
print("Ca-B distance", sp.N(dCaB, 8), " rL + rS =", sp.N(sp.Rational(1, 2) + rS, 8), " gap =", sp.N(dCaB - sp.Rational(1, 2) - rS, 6))
lo = ((s2 - 1) / 6) ** sp.Rational(1, 3); hi = (4 + 7 * s2) ** sp.Rational(-1, 3)
print("fcc window: gamma in (", sp.N(lo, 6), ",", sp.N(hi, 6), ") ; offsets from sqrt2-1:",
      sp.N(lo - (s2 - 1), 4), sp.N(hi - (s2 - 1), 4))
print("check window ends: phi_left(lo) - fcc =", sp.N(sp.pi / 6 * (1 + 6 * lo ** 3) - sp.pi / (3 * s2), 3),
      " phi_right(hi) - fcc =", sp.N(sp.pi / 6 * (s2 - 1) ** 3 * (6 + hi ** -3) - sp.pi / (3 * s2), 3))
print("NaCl at gamma=sqrt2-1:", sp.nsimplify(sp.simplify(sp.pi / (3 * s2) * (1 + (s2 - 1) ** 3))), "=",
      sp.N(sp.pi / (3 * s2) * (1 + (s2 - 1) ** 3), 6), "; fcc", sp.N(sp.pi / (3 * s2), 6))

# ---------------- brute force with every pair in a block ----------------
shifts = np.array(list(itertools.product(range(-2, 3), repeat=3)), float)


def sites(hh):
    B = []
    for m in range(3):
        for sg in (1, -1):
            p = np.full(3, 0.5); p[m] += sg * hh; B.append(p)
    return np.array(B)


def rL_max(hh, g):
    """largest rL (a = 1) with rS = g rL and no overlaps, over all pairs in a 5^3 block."""
    B = sites(hh)
    allB = (B[None] + shifts[:, None]).reshape(-1, 3)
    Ca = shifts
    lim = [0.5]  # Ca-Ca, nearest distance 1
    dCB = np.linalg.norm(B[:, None] - Ca[None], axis=2).min()
    lim.append(dCB / (1 + g))
    D = np.linalg.norm(B[:, None] - allB[None], axis=2); D[D < 1e-12] = np.inf
    lim.append(D.min() / (2 * g))
    return min(lim)


def phi_bf(g):
    f = lambda hh: -(4 * np.pi / 3) * rL_max(hh, g) ** 3 * (1 + 6 * g ** 3)
    best = min((minimize_scalar(f, bounds=(a, a + 0.05), method="bounded") for a in np.arange(0.0, 0.5, 0.05)),
               key=lambda r: r.fun)
    hstar = 1 - 2 ** 0.5 / 2
    if -f(hstar) > -best.fun:
        return -f(hstar), hstar
    return -best.fun, best.x


def phi_closed(g):
    g0 = 2 ** 0.5 - 1
    return np.pi / 6 * (1 + 6 * g ** 3) if g <= g0 else np.pi / 6 * g0 ** 3 * (6 + g ** -3)


if __name__ == "__main__":
    print("\ngamma | brute force phi (h) | closed form | diff")
    worst = 0
    for g in [0.1, 0.2, 0.3, 0.38, 0.40, 0.41, 0.4142, 2 ** 0.5 - 1, 0.415, 0.42, 0.45, 0.5, 0.7, 1.0, 1.5]:
        pb, hb = phi_bf(g)
        pc = phi_closed(g)
        worst = max(worst, abs(pb - pc))
        print(f"{g:.6f} | {pb:.6f} ({hb:.4f}) | {pc:.6f} | {pb-pc:+.1e}")
    print("max |brute - closed| =", f"{worst:.1e}")
    gg = np.linspace(0.05, 1.5, 3000)
    ph = np.array([phi_closed(g) for g in gg])
    print("argmax gamma (grid) =", gg[ph.argmax()], " max =", ph.max())
