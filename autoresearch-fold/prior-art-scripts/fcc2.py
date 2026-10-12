import numpy as np, sys, itertools
import os; FOLD_DIR = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))  # run from this folder
sys.path.insert(0,FOLD_DIR)
import verify
from collections import Counter, defaultdict
X0=np.load('final_X.npy'); L0=np.load('final_L.npy'); m=2; R=verify.R; D=2*R
I,J,S=verify.bar_graph(m)
a,b,c=L0; n=np.cross(a,b); n/=np.linalg.norm(n); hc=c@n
z=X0@n; lay=np.round((z-z.min())/(hc/4)).astype(int)%4
e1=a/np.linalg.norm(a); e2=np.cross(n,e1)
u=lambda deg: (D/np.sqrt(3))*(np.cos(np.radians(deg))*e1+np.sin(np.radians(deg))*e2)
def contacts(X,L,tol=1e-7):
    out=defaultdict(set)
    for i in range(len(X)):
        for j in range(len(X)):
            for s in itertools.product((-1,0,1),repeat=3):
                if i==j and s==(0,0,0): continue
                r=np.linalg.norm(X[j]+np.array(s)@L-X[i])
                if r<D*(1-1e-9): return None
                if r<D*(1+tol): out[i].add((j,s))
    return out
def cna(X,L,ct):
    sig=Counter()
    def nb(i,off): return {(j,tuple(np.add(s,off))) for j,s in ct[i]}
    for i in range(len(X)):
        N0=nb(i,(0,0,0))
        for (j,s) in N0:
            com=N0&nb(j,s)
            bl=[(p,q) for p in com for q in com if p<q and q in nb(*p)]
            shared=len(bl)==2 and len(set(bl[0])&set(bl[1]))>0
            sig[(len(com),len(bl),2 if shared else 1)]+=1
    return sig
res=[]
for sd in (None,0,120,240,60,180,300):
  for td in (0,120,240,60,180,300,None):
    X=X0.copy(); L=L0.copy()
    mask=(lay==2)|(lay==3)
    X[mask]+=u(sd) if sd is not None else 0
    if td is not None: L[2]=L[2]+u(td)
    bl=np.linalg.norm(X[J]+S@L-X[I],axis=1)
    if abs(bl-D).max()>1e-9: continue
    ct=contacts(X,L)
    if ct is None: continue
    deg=Counter(len(v) for v in ct.values())
    print('shift',sd,'c-shift',td,'deg',dict(deg),'CNA',dict(cna(X,L,ct)))
    res.append((sd,td))
print(res)
