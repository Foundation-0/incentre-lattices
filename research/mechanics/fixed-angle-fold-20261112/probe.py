#!/usr/bin/env python3
"""Fixed-angle fold model. Objective analytic forces and full tangent Hessian.

Only NumPy and SciPy are required. No modules from the geometry project are used.
All 24 angular interactions retain the same open-state preferred angle.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import minimize

HERE=Path(__file__).resolve().parent
P=np.array([sign*np.eye(3)[axis] for axis in range(3) for sign in (1,-1)])
C=np.array([[0.,-1.,-1.],[1.,0.,-1.],[1.,-1.,0.]])
BONDS=[(i,j,np.zeros(3)) for i,j in itertools.combinations(range(6),2) if i//2!=j//2]
BONDS += [(2*a,2*a+1,np.eye(3)[a]) for a in range(3)]
ANGLES=[(i,i^1, (1 if i%2==0 else -1)*np.eye(3)[i//2], k)
        for i in range(6) for k in range(6) if i//2!=k//2]
MANDEL=[np.diag(np.eye(3)[i]) for i in range(3)]
for i,j in [(1,2),(0,2),(0,1)]:
    E=np.zeros((3,3)); E[i,j]=E[j,i]=1/math.sqrt(2); MANDEL.append(E)
TRANSLATIONS=np.tile(np.eye(3),(6,1))/math.sqrt(6)
Z=null_space(TRANSLATIONS.T) # orthonormal 18 x 15, no translation
C0=-1/math.sqrt(2)


def geometry(s):
    a=math.sqrt(2)*(1-s*s)/(1+s*s)
    t=2*s/(1+s*s)
    return a,t,(2+a)*np.eye(3)+t*C


def cosine(b,d):
    """c and its first/second derivatives wrt the six Cartesian vector entries."""
    r=np.linalg.norm(b); l=np.linalg.norm(d); n=b/r; m=d/l; c=n@m
    g=np.r_[(m-c*n)/r,(n-c*m)/l]
    Hbb=(-np.outer(m,n)-np.outer(n,m)-c*np.eye(3)+3*c*np.outer(n,n))/(r*r)
    Hdd=(-np.outer(n,m)-np.outer(m,n)-c*np.eye(3)+3*c*np.outer(m,m))/(l*l)
    Hbd=(np.eye(3)-np.outer(n,n))@(np.eye(3)-np.outer(m,m))/(r*l)
    return c,g,np.block([[Hbb,Hbd],[Hbd.T,Hdd]])


def edge_map(i,j,n,v,q=None,affine=False):
    A=np.zeros((3,24 if affine else 18),dtype=complex if q is not None else float)
    A[:,3*i:3*i+3]-=np.eye(3)
    A[:,3*j:3*j+3]+=np.eye(3)*(np.exp(1j*np.dot(q,n)) if q is not None else 1)
    if affine:
        A[:,18:]=np.array([E@v for E in MANDEL]).T
    return A


def terms(s,kappa,u=None,e=None,q=None,affine=True,geometric=True):
    """Cell energy, gradient, Hessian in u and symmetric incremental strain.

    q selects the Hermitian Bloch Hessian. Its energy/gradient return slots are
    unused. Affine vectors refer to the selected unperturbed fold, so the map
    x=(I+E)x_ref+u is linear and its dead-load work has exactly zero Hessian.
    """
    _,_,L=geometry(s)
    u=np.zeros((6,3)) if u is None else np.asarray(u).reshape(6,3)
    E=np.zeros((3,3)) if e is None else sum(x*M for x,M in zip(e,MANDEL))
    dim=24 if affine else 18
    U=0.; g=np.zeros(dim,dtype=complex if q is not None else float)
    H=np.zeros((dim,dim),dtype=g.dtype)
    for i,j,n in BONDS:
        v0=P[j]-P[i]+n@L; v=v0+E@v0+u[j]-u[i]
        l=np.linalg.norm(v); unit=v/l; dl=l-math.sqrt(2)
        A=edge_map(i,j,n,v0,q,affine)
        grad=dl*unit
        Hv=np.outer(unit,unit)+(dl/l)*(np.eye(3)-np.outer(unit,unit))
        U+=.5*dl*dl; g+=A.conj().T@grad; H+=A.conj().T@Hv@A
    for i,j,n,k in ANGLES:
        b0=P[j]-P[i]+n@L; d0=P[k]-P[i]
        b=b0+E@b0+u[j]-u[i]; d=d0+E@d0+u[k]-u[i]
        c,dc,Hc=cosine(b,d); delta=c-C0
        A=np.vstack([edge_map(i,j,n,b0,q,affine),edge_map(i,k,np.zeros(3),d0,q,affine)])
        U+=.5*kappa*delta*delta
        g+=A.conj().T@(kappa*delta*dc)
        H+=A.conj().T@(kappa*(np.outer(dc,dc)+(delta*Hc if geometric else 0)))@A
    return U,g,H


def response(s,kappa):
    U,g,H=terms(s,kappa)
    Huu=Z.T@H[:18,:18]@Z; Hue=Z.T@H[:18,18:]
    eig=np.linalg.eigvalsh(Huu)
    du=-Z@np.linalg.solve(Huu,Hue)
    K=(H[18:,18:]+H[18:,:18]@du)/np.linalg.det(geometry(s)[2])
    K=(K+K.T)/2
    return dict(s=float(s),kappa=float(kappa),energy=float(U),
                holding_joint_force_norm=float(np.linalg.norm(g[:18])),
                holding_strain_force_norm=float(np.linalg.norm(g[18:])),
                internal_eigenvalues=eig.tolist(),relaxed_eigenvalues=np.linalg.eigvalsh(K).tolist(),
                exx_response_norm=float(np.linalg.norm(du[:,0])),
                geometric_stiffness_norm=float(np.linalg.norm(H-terms(s,kappa,geometric=False)[2]))),du,K


def checks():
    rng=np.random.default_rng(12617)
    out=[]
    for s in (0.,.001,.07,.4,.9):
        kappa=.003
        U,g,H=terms(s,kappa)
        a,t,_=geometry(s)
        np.testing.assert_allclose(U,3*kappa*((math.sqrt(2)-a)**2+t*t),atol=5e-16)
        assert np.linalg.norm(g[:18].reshape(6,3).sum(axis=0))<1e-14
        assert np.linalg.norm(H[:,:18]@TRANSLATIONS)<1e-14
        d=rng.normal(size=24); d/=np.linalg.norm(d)
        h=2e-5
        values=[]
        for sign in (-1,1):
            x=sign*h*d
            values.append(terms(s,kappa,x[:18],x[18:]))
        dg=(values[1][1]-values[0][1])/(2*h)
        herror=float(np.max(np.abs(dg-H@d)))
        gerror=float(abs((values[1][0]-values[0][0])/(2*h)-g@d))
        assert herror<3e-9 and gerror<1e-9
        # Objectivity is checked directly on every angle using a random rotation.
        axis=rng.normal(size=3); axis/=np.linalg.norm(axis); theta=.43
        cross=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
        Q=np.eye(3)+math.sin(theta)*cross+(1-math.cos(theta))*(cross@cross)
        b=rng.normal(size=3); d2=rng.normal(size=3)
        np.testing.assert_allclose(cosine(b,d2)[0],cosine(Q@b,Q@d2)[0],atol=5e-16)
        Hq=terms(s,kappa,q=np.zeros(3),affine=False)[2]
        np.testing.assert_allclose(Hq,H[:18,:18],atol=1e-15)
        out.append(dict(s=s,force_derivative_max_error=herror,energy_derivative_error=gerror))
    return out


def scan():
    rows=[]
    for kappa in (1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2):
        points=sorted(set([0.,1.,.25,.5,.75]+np.geomspace(1e-6,1.,61).tolist()+
                      [x*math.sqrt(kappa) for x in (.1,.2,.5,1,2,5,10) if x*math.sqrt(kappa)<=1]))
        for s in points:
            r,_,_=response(s,kappa); rows.append(r)
    return rows


def bloch(s,kappa,n=8,random_count=100):
    # The grid contains Gamma. Refine with small vectors since finite grids alone
    # miss long-wavelength instabilities when acoustic eigenvalues vanish as q².
    qs=[np.array(q)*2*math.pi/n for q in itertools.product(range(n),repeat=3)]
    rng=np.random.default_rng(149)
    qs += [rng.uniform(-math.pi,math.pi,size=3) for _ in range(random_count)]
    directions=[np.eye(3)[i] for i in range(3)]+[np.array([1.,1.,1.])/math.sqrt(3)]
    directions += [v/np.linalg.norm(v) for v in rng.normal(size=(20,3))]
    small=[]
    for radius in (.1,.01,.001,.0001):
        for d in directions: small.append(radius*d)
    qs+=small
    minimum=float('inf'); witness=None; negatives=0
    scaled_min=float('inf')
    for q in qs:
        H=terms(s,kappa,q=q,affine=False)[2]
        if np.linalg.norm(q)<1e-14:
            values=np.linalg.eigvalsh(Z.T@H.real@Z)
        else:
            values=np.linalg.eigvalsh(H)
        val=float(values[0])
        if val<minimum: minimum,witness=val,q.tolist()
        if val<-1e-10: negatives+=1
        if 0<np.linalg.norm(q)<=.10000001:
            scaled_min=min(scaled_min,val/np.dot(q,q))
    return dict(s=s,kappa=kappa,wavevectors=len(qs),grid_n=n,random_wavevectors=random_count,
                minimum_eigenvalue=minimum,witness_q=witness,negative_count_below_minus_1e_10=negatives,
                minimum_small_q_eigenvalue_over_q_squared=scaled_min)



def acoustic(s,kappa,direction):
    """Exact q Taylor coefficients, avoiding subtraction of nearly equal H(q)."""
    _,_,L=geometry(s); direction=np.asarray(direction);direction=direction/np.linalg.norm(direction)
    def maps(i,j,n,v):
        A0=edge_map(i,j,n,v,q=np.zeros(3),affine=False)
        A1=np.zeros_like(A0);A2=np.zeros_like(A0);phase=n@direction
        A1[:,3*j:3*j+3]=1j*phase*np.eye(3)
        A2[:,3*j:3*j+3]=-.5*phase**2*np.eye(3)
        return A0,A1,A2
    H0=np.zeros((18,18),complex);H1=H0.copy();H2=H0.copy()
    locals=[]
    for i,j,n in BONDS:
        v=P[j]-P[i]+n@L;unit=v/np.linalg.norm(v)
        locals.append((maps(i,j,n,v),np.outer(unit,unit)))
    for i,j,n,k in ANGLES:
        b=P[j]-P[i]+n@L;d=P[k]-P[i]
        c,g,Hc=cosine(b,d)
        bm=maps(i,j,n,b);dm=maps(i,k,np.zeros(3),d)
        locals.append((tuple(np.vstack([bm[z],dm[z]]) for z in range(3)),
                       kappa*(np.outer(g,g)+(c-C0)*Hc)))
    for (A0,A1,A2),h in locals:
        H0+=A0.conj().T@h@A0
        H1+=A1.conj().T@h@A0+A0.conj().T@h@A1
        H2+=A2.conj().T@h@A0+A1.conj().T@h@A1+A0.conj().T@h@A2
    Hred=Z.T@H0@Z
    coefficient=TRANSLATIONS.T@(H2-H1@Z@np.linalg.solve(Hred,Z.T@H1))@TRANSLATIONS
    coefficient=(coefficient+coefficient.conj().T)/2
    values=np.linalg.eigvalsh(coefficient)
    # Independent finite-q comparison at a radius large enough for double precision.
    rho=.01
    direct=np.linalg.eigvalsh(terms(s,kappa,q=rho*direction,affine=False)[2])[:3]/rho**2
    return dict(direction=direction.tolist(),coefficient_eigenvalues=values.tolist(),
                direct_q_radius=rho,direct_three_eigenvalues_over_q_squared=direct.tolist())


def acoustic_scan(s,kappa):
    rng=np.random.default_rng(718)
    dirs=[np.eye(3)[i] for i in range(3)]
    dirs += [np.array(x) for x in itertools.product((-1,0,1),repeat=3) if x!=(0,0,0)]
    dirs += list(rng.normal(size=(120,3)))
    rows=[acoustic(s,kappa,d) for d in dirs]
    worst=min(rows,key=lambda x:x['coefficient_eigenvalues'][0])
    return dict(s=s,kappa=kappa,sampled_directions=len(dirs),
                minimum_coefficient_eigenvalue=worst['coefficient_eigenvalues'][0],
                witness=worst,negative_directions_below_minus_1e_10=sum(r['coefficient_eigenvalues'][0]<-1e-10 for r in rows))


def nonlinear(s,kappa,h):
    row,du,K=response(s,kappa)
    U0,g0,_=terms(s,kappa)
    e=np.zeros(6);e[0]=h
    expected=du[:,0]*h
    def fun(v):
        u=Z@v; U,g,H=terms(s,kappa,u,e)
        W=U-U0-g0[:18]@u-g0[18:]@e
        return W,Z.T@(g[:18]-g0[:18]),Z.T@H[:18,:18]@Z
    opt=minimize(lambda x:fun(x)[0],Z.T@expected,jac=lambda x:fun(x)[1],hess=lambda x:fun(x)[2],
                 method='trust-exact',options={'gtol':1e-13,'maxiter':100})
    point=opt.x.copy(); polishing_steps=0
    # Energy subtraction can stop trust-exact near floating-point resolution.
    # Polish the force equation directly and retain the original solver status.
    for _ in range(8):
        _,residual,hessian=fun(point)
        step=np.linalg.solve(hessian,residual)
        if np.linalg.norm(step)<1e-13:break
        candidate=point-step
        if np.linalg.norm(fun(candidate)[1])>=np.linalg.norm(residual):break
        point=candidate;polishing_steps+=1
    u=Z@point; W,grad,Huu=fun(point)
    correction=float(np.linalg.norm(np.linalg.solve(Huu,grad)))
    stationary=correction<max(1e-12,1e-7*np.linalg.norm(u)) and np.linalg.eigvalsh(Huu)[0]>0
    assert stationary, 'Finite-strain state did not pass the force-residual/positive-curvature check'
    pred=.5*np.linalg.det(geometry(s)[2])*K[0,0]*h*h
    return dict(s=s,kappa=kappa,strain=h,solver_success=bool(opt.success),message=str(opt.message),
                gradient_norm=float(np.linalg.norm(grad)),iterations=int(opt.nit),
                stationary_point_verified=bool(stationary),newton_polishing_steps=polishing_steps,
                residual_displacement_error_estimate=float(np.linalg.norm(np.linalg.solve(Huu,grad))),
                internal_hessian_min=float(np.linalg.eigvalsh(Huu)[0]),
                displacement_norm=float(np.linalg.norm(u)),linear_prediction_norm=float(np.linalg.norm(expected)),
                relative_displacement_error=float(np.linalg.norm(u-expected)/np.linalg.norm(expected)),
                energy=float(W),quadratic_energy_prediction=float(pred),
                relative_energy_error=float(abs(W-pred)/max(abs(pred),1e-30)))



def endpoint_controls():
    rows=[]
    for kappa in (1e-8,1e-7,1e-6,1e-5,1e-4):
        q=np.full(3,math.pi)
        H=terms(1.,kappa,q=q,affine=False)[2].real
        eig,vec=np.linalg.eigh(H);v=vec[:,0]
        v*=np.sign(v[np.argmax(np.abs(v))])
        no_geometric=terms(1.,kappa,q=q,affine=False,geometric=False)[2]
        rows.append(dict(kappa=kappa,minimum_eigenvalue=float(eig[0]),
                         eigenvalue_over_kappa=float(eig[0]/kappa),
                         residual=float(np.linalg.norm(H@v-eig[0]*v)),
                         omitted_geometric_stiffness_minimum=float(np.linalg.eigvalsh(no_geometric)[0])))
        assert eig[0]<-.28*kappa and np.linalg.norm(H@v-eig[0]*v)<5e-14
        if kappa==1e-6:mode=v.reshape(6,3);prediction=float(eig[0])
    # Finite 2x2x2 supercell, with real displacement (-1)^(nx+ny+nz)*h*mode.
    # Every interaction is enumerated directly. No Bloch matrix is used here.
    _,_,L=geometry(1.);cells=list(itertools.product(range(2),repeat=3))
    def energy_change(h):
        total=0.
        def displacement(cell,i):
            return (-1)**int(sum(cell))*h*mode[i]
        for cell_tuple in cells:
            cell=np.array(cell_tuple)
            for i,j,n in BONDS:
                v0=P[j]-P[i]+n@L
                v=v0+displacement(cell+n,j)-displacement(cell,i)
                total+=.5*(np.linalg.norm(v)-math.sqrt(2))**2
            for i,j,n,k in ANGLES:
                b0=P[j]-P[i]+n@L;d0=P[k]-P[i]
                b=b0+displacement(cell+n,j)-displacement(cell,i)
                d=d0+displacement(cell,k)-displacement(cell,i)
                c=cosine(b,d)[0];c_ref=cosine(b0,d0)[0]
                dc=c-c_ref
                total+=1e-6*((c_ref-C0)*dc+.5*dc*dc)
        # Joint holding-work sums to zero for this alternating displacement.
        return total
    finite=[]
    for h in (1e-4,3e-5,1e-5):
        measured=(energy_change(h)+energy_change(-h))/(8*h*h)
        finite.append(dict(amplitude=h,measured_curvature=measured,
                           predicted_curvature=prediction,relative_error=abs(measured-prediction)/abs(prediction)))
    assert finite[-1]['relative_error']<.001 and finite[-1]['measured_curvature']<0
    return dict(q=[math.pi]*3,kappa_scan=rows,mode_at_kappa_1e_6=mode.tolist(),
                mode_convention='unit Euclidean norm per cell; multiply by (-1)^(nx+ny+nz)',
                supercell_energy_controls=finite)


def crossover_controls():
    L0=2+math.sqrt(2);xstar=3*L0/(2*math.sqrt(17));gstar=5/(3*math.sqrt(17))
    rows=[]
    for kappa in (1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2):
        s=xstar*math.sqrt(kappa);row,_,_=response(s,kappa)
        row['scaled_gain']=math.sqrt(kappa)*row['exx_response_norm']
        row['relative_asymptotic_error']=abs(row['scaled_gain']-gstar)/gstar
        rows.append(row)
    assert rows[0]['relative_asymptotic_error']<1e-5
    return dict(asymptotic_peak_x=xstar,asymptotic_peak_scaled_gain=gstar,rows=rows)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');parser.add_argument('--quick',action='store_true')
    args=parser.parse_args()
    data={'model':'fixed 24 cosine-angle interactions plus 15 central springs; dead nodal and strain holding forces',
          'checks':checks(),'scan':scan(),'crossover':crossover_controls(),'endpoint':endpoint_controls()}
    if not args.quick:
        selected=[(0.,1e-6),(.001,1e-6),(.003,1e-6),(.01,1e-6),(.1,1e-6),(.5,1e-6),(1.,1e-6),
                  (3*(2+math.sqrt(2))/(2*math.sqrt(17))*.01,1e-4)]
        data['bloch']=[bloch(s,k) for s,k in selected]
        data['acoustic']=[acoustic_scan(s,k) for s,k in selected]
        data['nonlinear']=[nonlinear(s,1e-6,h) for s in (.001,.003,.01) for h in (1e-7,-1e-7,1e-8,-1e-8)]
        data['nonlinear'] += [nonlinear(selected[-1][0],1e-4,h) for h in (1e-4,-1e-4,1e-5,-1e-5,1e-6,-1e-6)]
    path=HERE/('quick.json' if args.quick else 'evidence.json')
    if args.check:
        old=json.loads(path.read_text()); assert old==data,'Numerical output changed'
        print('PASS: output reproduced exactly')
    else:
        path.write_text(json.dumps(data,indent=2)+'\n')
        print(json.dumps({'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}))
    for k in sorted(set(x['kappa'] for x in data['scan'])):
        rows=[x for x in data['scan'] if x['kappa']==k]
        positive=[x for x in rows if x['internal_eigenvalues'][0]>0 and x['relaxed_eigenvalues'][0]>0]
        best=max(positive,key=lambda x:x['exx_response_norm']) if positive else None
        print(k,'stable homogeneous',len(positive),'max gain',None if best is None else (best['s'],best['exx_response_norm']))
    if 'bloch' in data:
        for b in data['bloch']:print('bloch',b)
        for b in data['acoustic']:print('acoustic',b)
        for r in data['nonlinear']:print('nonlinear',r)

if __name__=='__main__':main()
