import numpy as np, itertools, sys
import os; FOLD_DIR = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))  # run from this folder
sys.path.insert(0,FOLD_DIR)
import verify
d=np.load(FOLD_DIR+'/best/path.npz'); Ls=d['Ls']; Xs=d['Xs']; m=int(d['m'])
R=verify.R; D=2*R
I,J,S=verify.bar_graph(m)
T=len(Ls); N=Xs.shape[1]
dens=N*4/3*np.pi*R**3/np.linalg.det(Ls)
print('T',T,'N',N,'bars',len(I),'dens start/end',dens[0],dens[-1])
print('monotone density?', np.all(np.diff(dens)>=-1e-12), 'min diff',np.diff(dens).min())
# stage A length: count states until fold end (cubic lattice)
cub=[np.allclose(L,L[0,0]*np.eye(3),atol=1e-12) for L in Ls]
nA=np.argmin(cub) ; print('stage A states',nA,'density at end of A',dens[nA-1])
# path length in configuration space (mod nothing)
dl=[np.sqrt(np.sum((Xs[t+1]-Xs[t])**2)+np.sum((Ls[t+1]-Ls[t])**2)) for t in range(T-1)]
print('path length total %.3f, stageA %.3f'%(sum(dl),sum(dl[:nA-1])))
# octahedra rigidity: check all 15 pair distances per octahedron constant
X=Xs[-1];L=Ls[-1]
def octdist(X,c):
    P=X[6*c:6*c+6]; return np.array([np.linalg.norm(P[i]-P[j]) for i,j in itertools.combinations(range(6),2)])
print('octahedron distance drift', max(abs(octdist(Xs[t],c)-octdist(Xs[0],c)).max() for t in range(0,T,50) for c in range(8)))
# end state: layers. hex axis = cross(a,b)
a,b,cv=L
n=np.cross(a,b); n/=np.linalg.norm(n)
h=(X@n)%(abs(cv@n)); 
print('c . n',cv@n, 'layer heights sorted unique',np.unique(np.round(h,4)))
lay=np.round(h/(abs(cv@n)/4)).astype(int)%4
for c in range(8):
    print('oct',c,'layers',sorted(lay[6*c:6*c+6]))
# inter-octahedron bars: layer of endpoints
inter=[(i,j,s) for i,j,s in zip(I,J,S) if i//6!=j//6 or np.any(s)]
from collections import Counter
print('inter bars',len(inter), Counter(tuple(sorted((lay[i],lay[j]))) for i,j,s in inter))
# symmetry along path: count distinct ball environments via sorted distance lists (coarse)
def envs(X,L):
    cl=verify._image_cloud(X,L,3*D) if hasattr(verify,'_image_cloud') else None
for t in [0,nA-1,(nA+T)//2,T-1]:
    X=Xs[t];L=Ls[t]
    ds=[]
    for i in range(N):
        dd=[]
        for j in range(N):
            for s in itertools.product((-1,0,1),repeat=3):
                r=np.linalg.norm(X[j]+np.array(s)@L-X[i])
                if 1e-9<r<2.2*D: dd.append(round(r,6))
        ds.append(tuple(sorted(dd)))
    print('t',t,'dens %.4f'%dens[t],'distinct environments',len(set(ds)))
np.save('final_X.npy',Xs[-1]); np.save('final_L.npy',Ls[-1])
