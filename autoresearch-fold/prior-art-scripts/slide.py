import numpy as np, sys
import os; FOLD_DIR = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))  # run from this folder
sys.path.insert(0,FOLD_DIR)
import verify
d=np.load(FOLD_DIR+'/best/path.npz'); Ls=list(d['Ls']); Xs=list(d['Xs']); m=2
R=verify.R; D=2*R
I,J,S=verify.bar_graph(m)
X0=Xs[-1].copy(); L0=Ls[-1].copy()
a,b,c=L0; n=np.cross(a,b); n/=np.linalg.norm(n); hc=c@n
z=X0@n; lay=np.round((z-z.min())/(hc/4)).astype(int)%4
e1=a/np.linalg.norm(a); e2=np.cross(n,e1)
u=lambda deg: (D/np.sqrt(3))*(np.cos(np.radians(deg))*e1+np.sin(np.radians(deg))*e2)
top=(lay==2)|(lay==3)
def interface_bar(X,L,la,lb):
    for i,j,s in zip(I,J,S):
        if {lay[i],lay[j]}=={la,lb} and i//6!=j//6:
            return i,j,s
def lifted(X,L,lat,mode,ib):
    i,j,s=ib
    def apply(dz):
        Y=X.copy(); M=L.copy(); sh=lat+dz*n
        if mode==1: Y[top]+=sh; M[2]+=sh
        else: M[2]+=sh
        return Y,M
    Y,M=apply(0.0); v=Y[j]+s@M-Y[i]
    sgn=np.sign(v@n); vl=v-(v@n)*n
    dz=sgn*np.sqrt(D**2-vl@vl)-(v@n)
    return apply(dz)
path_L=[];path_X=[]
X,L=X0,L0
for mode,vec,(la,lb) in ((1,u(0),(1,2)),(2,u(120),(3,0))):
    ib=interface_bar(X,L,la,lb)
    T=200
    for t in range(1,T+1):
        Y,M=lifted(X,L,vec*t/T,mode,ib)
        path_X.append(Y);path_L.append(M)
    X,L=path_X[-1],path_L[-1]
bl=[abs(np.linalg.norm(Y[J]+S@M-Y[I],axis=1)-D).max() for Y,M in zip(path_X,path_L)]
md=[verify.min_pair_distance(Y,M) for Y,M in zip(path_X,path_L)]
dens=[48*4/3*np.pi*R**3/np.linalg.det(M) for M in path_L]
print('max bar err',max(bl),'min dist/D',min(md)/D,'density min',min(dens),'end',dens[-1])
np.savez('vcheck/path_fcc.npz',m=2,Ls=np.array(Ls+path_L),Xs=np.array(Xs+path_X))
