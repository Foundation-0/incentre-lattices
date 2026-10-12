import numpy as np, itertools, collections
exec(open('chk.py').read().split("# overlap check")[0].split("print('state0")[0])
def contactsets(X,L,tol=1e-6):
  ns=list(itertools.product(range(-2,3),repeat=3)); N=len(X); nb=collections.defaultdict(set)
  for i in range(N):
    for j in range(N):
      for n in ns:
        if i==j and n==(0,0,0): continue
        r=np.linalg.norm(X[j]-X[i]+np.array(n)@L)
        if r<B*(1+tol): nb[(i,0,0,0)].add((j,)+n)
  return nb
def contacts(X,L,tol=1e-6):
  return {(i,j,n) for i in range(len(X)) for j in range(len(X)) for n in itertools.product(range(-2,3),repeat=3) if not(i==j and n==(0,0,0)) and np.linalg.norm(X[j]-X[i]+np.array(n)@L)<B*(1+tol)}
c0=contacts(Xs[0],Ls[0]); c1=contacts(Xs[-1],Ls[-1])
print('initial contact pairs(directed)',len(c0),'per ball',len(c0)/48,' final',len(c1))
print('initial contacts retained at end',len(c0&c1),'new at end',len(c1-c0))
# CNA at final
ctab=collections.defaultdict(set)
for i,j,n in c1: ctab[(i,(0,0,0))].add((j,n))
def nbrs(i,off):  # neighbours of atom i at image off (in absolute image coords)
  return {(j,tuple(np.add(n,off))) for (j,n) in ctab[(i,(0,0,0))]}
sig=collections.Counter()
for i in range(48):
  N0=nbrs(i,(0,0,0)); sigs=[]
  for (j,nj) in N0:
    Nj=nbrs(j,nj); common=N0&Nj
    bonds=[(a,b) for a in common for b in common if a<b and b in nbrs(*a)]
    # chain
    adj=collections.defaultdict(list)
    for a,b in bonds: adj[a].append(b); adj[b].append(a)
    # longest chain via simple path search
    def lp(u,seen):
      return max([1+lp(v,seen|{v}) for v in adj[u] if v not in seen]+[0])
    ch=max([lp(u,{u}) for u in common]+[0])
    sigs.append((len(common),len(bonds),ch))
  sig[tuple(sorted(collections.Counter(sigs).items()))]+=1
for k,v in sig.items(): print(v,'atoms with',k)
# layer analysis: fraction of 12 neighbors forming hcp (cuboctahedron vs anticuboctahedron) -- done by CNA (421 vs 422)
