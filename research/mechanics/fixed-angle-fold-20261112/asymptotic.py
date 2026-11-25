#!/usr/bin/env python3
"""Exact reduction of the soft rotations. Requires SymPy. No numerical fit."""
import json
from pathlib import Path
import sympy as S

r=S.sqrt(2); L0=2+r
p=[sign*S.eye(3)[:,i] for i in range(3) for sign in (1,-1)]
C=S.Matrix([[0,-1,-1],[1,0,-1],[1,-1,0]])
bonds=[(i,j,S.zeros(3,1)) for i in range(6) for j in range(i+1,6) if i//2!=j//2]
bonds +=[(2*a,2*a+1,S.eye(3)[:,a]) for a in range(3)]
R0=S.zeros(15,18);R1=S.zeros(15,18);B0=S.zeros(15,6)
M=[S.diag(*[int(k==i) for k in range(3)]) for i in range(3)]
for i,j in [(1,2),(0,2),(0,1)]:
    e=S.zeros(3);e[i,j]=e[j,i]=1/r;M.append(e)
for b,(i,j,n) in enumerate(bonds):
    v=p[j]-p[i]+L0*n;d=v/r
    R0[b,3*i:3*i+3]=-d.T;R0[b,3*j:3*j+3]=d.T
    for k,e in enumerate(M):B0[b,k]=(d.T*e*v)[0]
    if b>=12:
        dp=r*C[b-12,:]
        R1[b,3*i:3*i+3]=-dp;R1[b,3*j:3*j+3]=dp
rot=S.zeros(18,3)
for i,pi in enumerate(p):
    for j in range(3):rot[3*i:3*i+3,j]=S.eye(3)[:,j].cross(pi)
W=S.Matrix.hstack(*R0.T.nullspace());P=W*(W.T*W).inv()*W.T
A2=S.simplify(rot.T*R1.T*P*R1*rot)
f=S.simplify(rot.T*R1.T*P*B0[:,0])
# Angular derivatives at s=0 assembled independently from the vector formula.
G=S.zeros(3)
for i in range(6):
    j=i^1; n=(1 if i%2==0 else -1)*S.eye(3)[:,i//2]
    b=p[j]-p[i]+L0*n; nb=b/r
    for k in range(6):
        if k//2==i//2:continue
        d=p[k]-p[i];nd=d/r;c=(nb.T*nd)[0]
        gb=(nd-c*nb)/r;gd=(nb-c*nd)/r
        row=gb.T*(rot[3*j:3*j+3,:]-rot[3*i:3*i+3,:])+gd.T*(rot[3*k:3*k+3,:]-rot[3*i:3*i+3,:])
        G+=row.T*row
G=S.simplify(G)
assert A2==S.Matrix([[S.Rational(16,3),S.Rational(8,3),S.Rational(8,3)],[S.Rational(8,3),S.Rational(80,9),-S.Rational(56,9)],[S.Rational(8,3),-S.Rational(56,9),S.Rational(80,9)]])
assert S.simplify(f-S.Rational(20,9)*(1+r)*S.Matrix([0,-1,1]))==S.zeros(3,1)
assert S.simplify(G-2*L0**2*S.eye(3))==S.zeros(3)
assert rot.T*rot==4*S.eye(3)
result={'central_rotation_matrix':str(A2),'exx_rotation_forcing':str(f),'angular_rotation_matrix':str(G),
        'scaled_variable':'x=s/sqrt(kappa)',
        'scaled_gain':'sqrt(kappa)*||du/de_xx|| -> (40*(2+sqrt(2))/9)*x / ((136/9)*x^2+2*(2+sqrt(2))^2)',
        'peak_x':'3*(2+sqrt(2))/(2*sqrt(17))','peak_scaled_gain':'5/(3*sqrt(17))',
        'peak_x_decimal':float(3*L0/(2*S.sqrt(17))),'peak_scaled_gain_decimal':float(5/(3*S.sqrt(17)))}
path=Path(__file__).with_name('asymptotic.json');path.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
