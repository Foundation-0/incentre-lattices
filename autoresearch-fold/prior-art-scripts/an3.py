import numpy as np, sys
import os; FOLD_DIR = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))  # run from this folder
sys.path.insert(0,FOLD_DIR)
import verify
from collections import Counter
X=np.load('final_X.npy'); L=np.load('final_L.npy'); m=2
I,J,S=verify.bar_graph(m)
a,b,c=L; n=np.cross(a,b); n/=np.linalg.norm(n); hc=c@n
z=X@n; lay=np.round((z-z.min())/(hc/4)).astype(int)%4
print('oct layers',[tuple(sorted(set(lay[6*k:6*k+6]))) for k in range(8)])
e1=a/np.linalg.norm(a); e2=np.cross(n,e1)
rows=[];inplane=0
for i,j,s in zip(I,J,S):
    if i//6==j//6 and not np.any(s): continue  # octahedron edge
    v=X[j]+s@L-X[i]; vz=v@n
    if abs(vz)<1e-6: inplane+=1; continue
    li,lj=lay[i],lay[j]
    if vz<0: v=-v; li,lj=lj,li
    rows.append(((int(li),int(lj)),round(float(np.degrees(np.arctan2(v@e2,v@e1))%360),1)))
print('inter-oct bars in-plane',inplane); print(Counter(rows))
