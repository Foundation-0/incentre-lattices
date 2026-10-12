import numpy as np, sys, glob
import os; FOLD_DIR = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))  # run from this folder
exec(open('fcc2.py').read().split('res=[]')[0])
ref=np.load(FOLD_DIR+'/best/path.npz')
for f in sorted(glob.glob('run_s3_f[23]*/out/path.npz')):
    d=np.load(f); X=d['Xs'][-1]; L=d['Ls'][-1]
    same = d['Xs'].shape==ref['Xs'].shape and np.allclose(d['Xs'],ref['Xs'])
    ct=contacts(X,L); print(f, len(d['Xs']), 'identical to best path:',same, dict(cna(X,L,ct)))
