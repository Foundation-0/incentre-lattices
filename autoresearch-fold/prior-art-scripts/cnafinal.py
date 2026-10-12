import numpy as np, sys
exec(open('fcc2.py').read().split('res=[]')[0])
d=np.load('vcheck/path_fcc.npz'); X=d['Xs'][-1]; L=d['Ls'][-1]
ct=contacts(X,L); print(dict(cna(X,L,ct)))
