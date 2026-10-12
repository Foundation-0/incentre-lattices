import numpy as np, itertools, sys
import os; FOLD_DIR = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))  # run from this folder
sys.path.insert(0,FOLD_DIR)
import verify
X=np.load('final_X.npy'); L=np.load('final_L.npy'); m=2; R=verify.R; D=2*R
I,J,S=verify.bar_graph(m)
a,b,c=L; n=np.cross(a,b); n/=np.linalg.norm(n); hc=c@n
z=X@n; lay=np.round((z-z.min())/(hc/4)).astype(int)%4
e1=a/np.linalg.norm(a); e2=np.cross(n,e1)
for i,j,s in zip(I,J,S):
    v=X[j]+s@L-X[i]
    li,lj=lay[i],lay[j]
    if li!=lj and i//6!=j//6 or (li!=lj and np.any(s) ):
        pass
rows=[]
for i,j,s in zip(I,J,S):
    v=X[j]+s@L-X[i]; li,lj=lay[i],lay[j]
    vz=v@n
    if abs(vz)>1e-6 and {li,lj} in ({1,0},{2,3},{0,1},{3,2}) and not ({li,lj}=={1,2} or {li,lj}=={3,0}):
        if vz<0: v=-v; li,lj=lj,li
        ang=np.degrees(np.arctan2(v@e2,v@e1))%360
        rows.append(((li,lj),round(ang,2)))
from collections import Counter
print(Counter(rows))
