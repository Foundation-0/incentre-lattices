import numpy as np, itertools, sys
d=np.load('../best/path.npz'); Ls=d['Ls']; Xs=d['Xs']; m=int(d['m'])
S2=np.sqrt(2); R=S2-1; k=2-S2; B=2*R
T=len(Xs); print('states',T,'N',Xs.shape[1])
# canonical
cells=list(itertools.product(range(m),repeat=3))
dots=[];idx={}
for c in cells:
  for ax in range(3):
    for s in (+1,-1):
      idx[(c,ax,s)]=len(dots)
      p=2*np.array(c)+1.0; p[ax]+=s*k; dots.append(p)
P0=np.array(dots); L0=2*m*np.eye(3)
# bars: (i,j,shift n integer, in lattice coords)
bars=[]
for c in cells:
  for a in range(3):
    for b in range(a+1,3):
      for sa in (1,-1):
        for sb in (1,-1):
          bars.append((idx[(c,a,sa)],idx[(c,b,sb)],(0,0,0)))
  for a in range(3):
    c2=list(c); c2[a]+=1; n=[0,0,0]
    if c2[a]==m: c2[a]=0; n[a]=1
    bars.append((idx[(c,a,1)],idx[(tuple(c2),a,-1)],tuple(n)))
print('bars',len(bars))
# canonical check
def blen(X,L,bars,rowconv=True):
  out=[]
  for i,j,n in bars:
    sh=np.array(n)@L if rowconv else L@np.array(n)
    out.append(np.linalg.norm(X[j]+sh-X[i]))
  return np.array(out)
print('canon bar err',abs(blen(P0,L0,bars)-B).max())
# convention + ordering for state 0
X=Xs[0];L=Ls[0]
print('L0\n',L); print('det',np.linalg.det(L))
# gram invariants of state0 vs canonical: sorted pair distances w/ images
def pairset(X,L,rng=1):
  fr=np.linalg.solve(L.T,X.T).T if True else None
  return None
for row in (True,False):
  e=abs(blen(X,L,bars,row)-B).max(); print('direct order, rowconv',row,'err',e)
print('state0 vs canon max dev',abs(X-P0).max())
for row in (True,False):
  e=max(abs(blen(Xs[t],Ls[t],bars,row)-B).max() for t in range(T)); print('rowconv',row,'max abs bar err over path',e,'rel',e/B)
ROW=True
Lrow=lambda L:L
# bar check w/ fixed integer shifts
E=np.array([abs(blen(Xs[t],Ls[t],bars,ROW)-B).max() for t in range(T)])
print('fixed-shift max rel err',E.max()/B)
# min-image alternative
def minimg_err(X,L):
  worst=0
  for i,j,n in bars:
    best=1e9
    for s in itertools.product((-1,0,1),repeat=3):
      best=min(best,abs(np.linalg.norm(X[j]-X[i]+np.array(s)@L)-B))
    worst=max(worst,best)
  return worst
print('minimg err at samples',[minimg_err(Xs[t],Ls[t])/B for t in (0,T//2,T-1)])
# overlap check
def mind(X,L,exclude=True,rngn=2):
  # all pairs incl images; returns min non-bonded dist, min bonded-pair dist deviations
  shifts=np.array(list(itertools.product(range(-rngn,rngn+1),repeat=3)))@L
  N=len(X); bonded={(i,j,n) for i,j,n in bars}
  bset=set()
  for i,j,n in bars: bset.add((i,j,n)); bset.add((j,i,tuple(-x for x in n)))
  D=X[None,:,None,:]-X[:,None,None,:]+shifts[None,None,:,:]  # i,j,s
  r=np.linalg.norm(D,axis=-1)
  ns=np.array(list(itertools.product(range(-rngn,rngn+1),repeat=3)))
  mask=np.ones_like(r,bool)
  for a in range(N): mask[a,a,np.where((ns==0).all(1))[0][0]]=False
  for (i,j,n) in bset:
    s=np.where((ns==np.array(n)).all(1))[0][0]; mask[i,j,s]=False
  return r[mask].min(), r
ovl=[]
for t in range(T):
  mn,_=mind(Xs[t],Ls[t]); ovl.append(mn)
ovl=np.array(ovl); print('min nonbonded dist/2R over path: min',ovl.min()/B,'at',ovl.argmin())
print('violations (<2R(1-1e-9))',(ovl<B*(1-1e-9)).sum())
# step sizes: displacement in cartesian with unwrapped, and fractional continuity
dX=np.linalg.norm(Xs[1:]-Xs[:-1],axis=2).max(1)
dL=np.abs(Ls[1:]-Ls[:-1]).max((1,2))
print('max atom step/R',dX.max()/R,'max L step',dL.max())
# fractional coord jump
Fr=np.array([np.linalg.solve(Ls[t].T,Xs[t].T).T for t in range(T)])
print('max frac step',np.abs(np.diff(Fr,axis=0)).max())
# interpolated substates
w=10; worst=1e9; bw=0
for t in range(T-1):
  for a in range(1,w):
    u=a/w; X=(1-u)*Xs[t]+u*Xs[t+1]; L=(1-u)*Ls[t]+u*Ls[t+1]
    mn,_=mind(X,L,rngn=1) if False else mind(X,L)
    if mn<worst: worst=mn;wt=(t,a)
    bw=max(bw,abs(blen(X,L,bars,ROW)-B).max())
print('interp min nonbonded/2R',worst/B,wt,'interp max bar dev/B',bw/B)
# final state
X=Xs[-1];L=Ls[-1]; N=len(X)
print('final L\n',L,'\ndet',np.linalg.det(L))
rho=N*4/3*np.pi*R**3/abs(np.linalg.det(L)); print('density',rho,'pi/sqrt18',np.pi/np.sqrt(18),'diff',rho-np.pi/np.sqrt(18))
_,r=mind(X,L)
cnt=(r<B*(1+1e-6)).sum(axis=(1,2))-0  # i->all j incl self? self shift0 r=0
cnt=((r<B*(1+1e-6))&(r>1e-9)).sum(axis=(1,2)); print('contacts per ball',np.unique(cnt,return_counts=True))
print('gap just beyond contact: next dist/2R',np.sort(r[(r>B*(1+1e-6))].ravel())[:1]/B)
np.save('final.npy',np.concatenate([L,X]))
