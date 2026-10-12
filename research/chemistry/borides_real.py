"""Compare real MB6 (CaB6 type) and MB12 (UB12 type) structures from COD with the 'dot' ideal.
For each CIF: lattice a, B-B inter-cluster (multiplicity 1) and intra-cluster (multiplicity 4) distances,
effective boron parameter x, and M-B distance. Ideal (all balls touching): inter = intra.
"""
import json, sys, warnings
import numpy as np
from collections import defaultdict
from ase.io import read
warnings.filterwarnings("ignore")

idx = json.load(open("cod_index.json"))
rows = []
for e in idx:
    try:
        at = read(f"cod/{e['cod']}.cif")
    except Exception as ex:
        continue
    sym = np.array(at.get_chemical_symbols())
    if not set(sym) <= {"B", e["M"]} or "B" not in sym:
        continue
    a = at.cell.lengths()[0]
    sup = at.repeat((3, 3, 3))
    P = sup.get_positions(); S = np.array(sup.get_chemical_symbols())
    ctr = P.mean(0)
    Bi = np.where(S == "B")[0]
    i0 = Bi[np.argmin(np.linalg.norm(P[Bi] - ctr, axis=1))]
    d = np.linalg.norm(P[Bi] - P[i0], axis=1)
    d = np.sort(d[d > 0.1])[:5]
    # group: one inter + four intra
    vals = np.round(d, 3)
    uniq, cnt = np.unique(vals, return_counts=True)
    if len(uniq) == 1:
        inter = intra = uniq[0]
    else:
        inter = uniq[cnt == 1][0] if (cnt == 1).any() else np.nan
        intra = uniq[cnt == 4][0] if (cnt == 4).any() else np.nan
    Mi = np.where(S == e["M"])[0]
    dMB = np.min(np.linalg.norm(P[Mi][:, None] - P[Bi][None], axis=2))
    occ = None
    rows.append(dict(M=e["M"], n=e["n"], cod=e["cod"], a=a, inter=inter, intra=intra, dMB=dMB, T=e["T"], P=e["P"]))

# x_eff: MB6 x = inter/(2a); MB12 x from intra: intra = x*sqrt2*a
for r in rows:
    if r["n"] == 6:
        r["x"] = r["inter"] / (2 * r["a"]); r["x_ideal"] = (2 ** 0.5 - 1) / 2
    else:
        r["x"] = r["intra"] / (2 ** 0.5 * r["a"]); r["x_ideal"] = 1 / 6
    r["ratio"] = r["inter"] / r["intra"]
    r["dBB_mean"] = (r["inter"] + 4 * r["intra"]) / 5

json.dump(rows, open("borides_real.json", "w"), indent=1, default=float)
print(f"{'M':3s}{'n':>3s} {'COD':>8s} {'a(A)':>7s} {'B-B inter':>9s} {'intra':>7s} {'inter/intra':>11s} {'x':>7s} {'x_ideal':>7s} {'M-B':>6s} {'M-B/dBB':>7s}  T  P")
for r in sorted(rows, key=lambda r: (r["n"], r["M"], r["cod"])):
    print(f"{r['M']:3s}{r['n']:3d} {r['cod']:>8s} {r['a']:7.4f} {r['inter']:9.3f} {r['intra']:7.3f} {r['ratio']:11.3f} {r['x']:7.4f} {r['x_ideal']:7.4f} {r['dMB']:6.3f} {r['dMB']/r['dBB_mean']:7.3f}  {r['T']} {r['P']}")
for n in (6, 12):
    R = [r for r in rows if r["n"] == n and np.isfinite(r["ratio"]) and not r["P"]]
    if R:
        rat = np.array([r["ratio"] for r in R]); dbb = np.array([r["dBB_mean"] for r in R]); a = np.array([r["a"] for r in R])
        mb = np.array([r["dMB"] / r["dBB_mean"] for r in R])
        k = (1 + 2 ** 0.5) if n == 6 else 3 * 2 ** 0.5
        print(f"\nMB{n}: N={len(R)} entries, metals {sorted(set(r['M'] for r in R))}")
        print(f"  inter/intra: mean {rat.mean():.4f} sd {rat.std():.4f} range {rat.min():.3f}-{rat.max():.3f}  (ideal 1)")
        print(f"  mean B-B {dbb.mean():.3f} A ; a/mean dBB = {np.mean(a/dbb):.4f} (ideal {k:.4f}); M-B/dBB = {mb.mean():.4f} (ideal {1.778824 if n==6 else 1.581139})")
