"""Thread D: hard-sphere binary packings of CaB6 type (L at 1a, 6 S at (x,1/2,1/2)) and UB12 type
(L at 4a, 48 S at (1/2,x,x)): max packing fraction vs size ratio gamma = rS/rL, optimising x and a.
Where is the 'dot' (all-touching Wythoff) point, and is it the packing optimum?"""
import numpy as np
from scipy.optimize import minimize_scalar
from itertools import product

def mindist_struct(kind, x):
    """Return functions giving contact distances in units of a for given x."""
    if kind == "CaB6":
        L = np.array([[0,0,0]]); S = []
        for i in range(3):
            for s in (x, 1-x):
                p = [0.5,0.5,0.5]; p[i] = s; S.append(p)
        S = np.array(S); nL, nS = 1, 6
    else:
        fcc = np.array([[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0]])
        L = fcc.copy()
        offs = []
        for (i,j) in [(0,1),(0,2),(1,2)]:
            for si, sj in product((1,-1),(1,-1)):
                o = [0,0,0]; o[i]=si*x; o[j]=sj*x; offs.append(o)
        S = np.array([c+np.array([.5,0,0])+o for c in fcc for o in offs]) % 1.0
        nL, nS = 4, 48
    T = np.array(list(product((-1,0,1),repeat=3)))
    def md(A,B,same):
        d = np.linalg.norm((A[:,None,None,:]+T[None,None,:,:]) - B[None,:,None,:], axis=3)
        if same: d[d<1e-9]=np.inf
        return d.min()
    return md(L,L,True), md(L,S,False), md(S,S,True), nL, nS

def best_phi(kind, g):
    def neg(x):
        LL, LS, SS, nL, nS = mindist_struct(kind, x)
        # largest rL (a=1) with rS=g rL: 2rL<=LL, rL(1+g)<=LS, 2 g rL<=SS
        rL = min(LL/2, LS/(1+g), SS/(2*g))
        return -(4/3*np.pi*(nL*rL**3 + nS*(g*rL)**3))
    lo, hi = (0.12, 0.33) if kind=="CaB6" else (0.08, 0.24)
    xs = np.linspace(lo, hi, 400); v = [neg(x) for x in xs]
    x0 = xs[int(np.argmin(v))]
    r = minimize_scalar(neg, bounds=(max(lo,x0-0.01), min(hi,x0+0.01)), method="bounded")
    return -r.fun, r.x

for kind, xid, gid in (("CaB6", (2**.5-1)/2, None), ("UB12", 1/6, None)):
    gs = np.linspace(0.2, 0.8, 121)
    res = [best_phi(kind, g) for g in gs]
    phis = np.array([r[0] for r in res]); i = int(np.argmax(phis))
    print(f"{kind}: max phi {phis[i]:.4f} at gamma {gs[i]:.3f}, x_opt {res[i][1]:.4f} (dot-ideal x {xid:.4f})")
    for g in (0.3,0.35,0.4,0.45,0.5,0.55,0.6):
        p,x = best_phi(kind,g); print(f"   gamma {g:.2f}: phi {p:.4f}, x {x:.4f}")
# the dot-ideal configurations themselves: S all touching, L as large as allowed
for kind, x in (("CaB6",(2**.5-1)/2),("UB12",1/6)):
    LL, LS, SS, nL, nS = mindist_struct(kind, x)
    rS = SS/2; rL = min(LL/2, LS - rS)
    print(f"{kind} at dot-ideal x={x:.4f}: rS={rS:.4f}a, rL={rL:.4f}a (limited by {'L-L' if np.isclose(rL,LL/2) else 'L-S'}), gamma={rS/rL:.4f}, phi={4/3*np.pi*(nL*rL**3+nS*rS**3):.4f}; cage-limited rL would be {LS-rS:.4f}a")
