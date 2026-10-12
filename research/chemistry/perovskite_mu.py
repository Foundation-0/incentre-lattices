"""Thread C: does the octahedral radius-ratio limit mu = rB/rX >= sqrt2-1 (the dot construction's cube value)
separate perovskites from non-perovskites in Bartel et al. 2019 (Sci. Adv. 5, eaav0693) Table S1 (576 ABX3)?"""
import numpy as np, pandas as pd
d = pd.read_csv("data/TableS1.csv")
y = d.exp_label.values == 1
mu = d["rB (Ang)"] / d["rX (Ang)"]; t = d.t.values
rArX = d["rA (Ang)"] / d["rX (Ang)"]
print(f"N={len(d)}, perovskites={y.sum()}")
def acc(pred): return (pred == y).mean()
c = 2**0.5 - 1
print(f"mu >= sqrt2-1          : acc {acc(mu>=c):.3f}; perovskites violating {np.sum(y & (mu<c))}/{y.sum()}; non-perov with mu<c {np.sum(~y&(mu<c))}/{(~y).sum()}")
ths = np.linspace(0.3,0.6,301); a = [acc(mu>=x) for x in ths]
print(f"best single mu threshold: {ths[int(np.argmax(a))]:.3f} acc {max(a):.3f}")
print(f"Goldschmidt 0.825<t<1.059: acc {acc((t>0.825)&(t<1.059)):.3f}")
print(f"t window AND mu>=sqrt2-1  : acc {acc((t>0.825)&(t<1.059)&(mu>=c)):.3f}")
print(f"Bartel tau<4.18 (reported): acc {acc(d.tau.values<4.18):.3f}")
# grid best (t_lo, t_hi, mu_c)
best=(0,)
for lo in np.arange(0.75,0.9,0.005):
  for hi in np.arange(0.98,1.1,0.005):
    for mc in np.arange(0.35,0.5,0.005):
      s=acc((t>lo)&(t<hi)&(mu>=mc))
      if s>best[0]: best=(s,lo,hi,mc)
print("best fitted t-window + mu threshold: acc %.3f lo %.3f hi %.3f mu_c %.3f"%best)
# by anion
for X in ["O","F","Cl","Br","I"]:
    m = d.X==X
    if m.sum(): print(f"  X={X:2s} N={m.sum():3d} perov={y[m].sum():3d}  min mu among perovskites {mu[m&y].min() if (m&y).any() else np.nan:.3f}  acc(mu>=c) {((mu[m]>=c)==y[m]).mean():.3f}")
