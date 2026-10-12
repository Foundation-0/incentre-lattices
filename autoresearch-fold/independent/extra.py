exec(open('chk.py').read().split("print('state0")[0].replace("print(","(lambda *a:None)("))
import collections
pairs=collections.Counter(tuple(sorted((i,j))) for i,j,n in bars)
print('max multiplicity of unordered bonded pair',max(pairs.values()),'self-bars',sum(i==j for i,j,n in bars))
deg=collections.Counter(); 
for i,j,n in bars: deg[i]+=1;deg[j]+=1
print('bar degree',collections.Counter(deg.values()))
# does a dot appear bonded to two images of same dot? (covered above) ; wrap bars' shifts
print('wrap bars',[b for b in bars if any(b[2])])
X=Xs[-1];L=Ls[-1]
ds=[np.linalg.norm(X[j]-X[i]+np.array(n)@L) for i,j,n in itertools.product(range(48),range(48),itertools.product(range(-2,3),repeat=3)) if not(i==j and n==(0,0,0))]
ds=np.array(ds); c=ds[ds<B*1.01]; print('contacts count',len(c),'max |d/2R-1|',abs(c/B-1).max())
# edge lengths: min nonbonded over all path states cleared 1e-9?
