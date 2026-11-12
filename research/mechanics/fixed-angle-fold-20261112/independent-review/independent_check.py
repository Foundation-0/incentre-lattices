#!/usr/bin/env python3
"""Independent geometry/energy checks; does not call the author's model code."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as s

HERE=Path(__file__).resolve().parent
AUTHOR=HERE.parent


def exact_checks():
    rt=s.sqrt(2)
    p=[sg*s.eye(3)[:,i] for i in range(3) for sg in (1,-1)]
    C=s.Matrix([[0,-1,-1],[1,0,-1],[1,-1,0]])
    z=s.symbols('z',real=True)
    L=(2+rt*(1-z*z)/(1+z*z))*s.eye(3)+2*z/(1+z*z)*C
    edges=[(i,j,s.zeros(3,1)) for i in range(6) for j in range(i+1,6) if i//2!=j//2]
    edges += [(2*a,2*a+1,s.eye(3)[:,a]) for a in range(3)]
    R=s.zeros(15,18)
    B=s.zeros(15,1)
    for index,(i,j,n) in enumerate(edges):
        v=p[j]-p[i]+L.T*n
        assert s.simplify(v.dot(v)-2)==0
        # Differentiate the exact unit bond, not a manually supplied R1.
        d=v/rt
        R[index,3*i:3*i+3]=-d.T
        R[index,3*j:3*j+3]=d.T
        B[index]=(d.T*s.diag(1,0,0)*v)[0]
    R0=R.subs(z,0)
    R1=R.diff(z).subs(z,0)
    B0=B.subs(z,0)
    rotations=s.zeros(18,3)
    for i,pi in enumerate(p):
        for j in range(3):rotations[3*i:3*i+3,j]=s.eye(3)[:,j].cross(pi)
    W=s.Matrix.hstack(*R0.T.nullspace())
    projector=W*(W.T*W).inv()*W.T
    M=s.simplify(rotations.T*R1.T*projector*R1*rotations)
    f=s.simplify(rotations.T*R1.T*projector*B0)
    assert rotations.T*rotations==4*s.eye(3)
    # Endpoint alternating rotations leave every bridging vector fixed exactly
    # to first order, since opposite motif vertices have opposite velocities.
    Hrot=s.zeros(3)
    Gopen=s.zeros(3)
    gauss=s.zeros(3)
    for i in range(6):
        j=i^1
        n=(1 if i%2==0 else -1)*s.eye(3)[:,i//2]
        for k in range(6):
            if i//2==k//2:continue
            d=p[k]-p[i]
            D=rotations[3*k:3*k+3,:]-rotations[3*i:3*i+3,:]
            b=p[j]-p[i]+L.subs(z,1).T*n
            bchange=-rotations[3*j:3*j+3,:]-rotations[3*i:3*i+3,:]
            assert bchange==s.zeros(3)
            assert d.T*D==s.zeros(1,3)
            cref=b.dot(d)/2
            gc=D.T*b/2
            # Direct second-order expansion of c= (b.d)/(sqrt(2)*|d|).
            Hrot+=gc*gc.T-(cref+1/rt)*cref*D.T*D/2
            gauss+=gc*gc.T
            # Open cell uniform rotation: use first variation of two vectors.
            b0=p[j]-p[i]+L.subs(z,0).T*n
            db=rotations[3*j:3*j+3,:]-rotations[3*i:3*i+3,:]
            c0=b0.dot(d)/2
            gc0=(db.T*d+D.T*b0)/2-c0*(db.T*b0+D.T*d)/2
            Gopen+=gc0*gc0.T
    Hrot=s.simplify(Hrot)
    Gopen=s.simplify(Gopen)
    trial=s.Matrix([-rt,1,1])
    quotient=s.simplify((trial.T*Hrot*trial)[0]/(4*trial.dot(trial)))
    eigmin=(2+rt-s.sqrt(38-4*rt))/8
    assert quotient==s.Rational(1,4)-3*rt/8
    # Reduction to the forced rotation axis and its Euclidean norm.
    direction=s.Matrix([0,-1,1])/rt
    central=s.simplify((direction.T*M*direction)[0])
    forcing=s.simplify((direction.T*f)[0])
    angular=s.simplify((direction.T*Gopen*direction)[0])
    x=s.symbols('x',positive=True)
    gain=s.simplify(2*forcing*x/(central*x*x+angular))
    xstar=s.sqrt(angular/central)
    gstar=s.simplify(gain.subs(x,xstar))
    assert s.simplify(xstar-3*(2+rt)/(2*s.sqrt(17)))==0
    assert s.simplify(gstar-5/(3*s.sqrt(17)))==0
    return {
        'central_null_dimension_open':len(R0.nullspace()),
        'left_null_dimension_open':W.cols,
        'central_rotation_matrix':str(M),
        'forcing':str(f),
        'open_angular_rotation_matrix':str(Gopen),
        'peak_s_over_sqrt_kappa':str(s.simplify(xstar)),
        'peak_gain_times_sqrt_kappa':str(gstar),
        'endpoint_rotation_matrix':str(Hrot),
        'endpoint_gauss_newton_rotation_matrix':str(s.simplify(gauss)),
        'endpoint_trial_rayleigh_over_kappa':str(quotient),
        'endpoint_trial_rayleigh_decimal':float(quotient),
        'endpoint_projected_minimum_over_kappa':str(eigmin),
        'endpoint_projected_minimum_decimal':float(eigmin),
        'endpoint_disposition':'negative trial curvature proves local instability for every kappa>0 in the stated dead-load point-joint model',
    }


def high_precision_endpoint():
    """Finite energy for an independently chosen exact negative trial mode."""
    mp.mp.dps=70
    rt=mp.sqrt(2)
    points=[mp.matrix([sg if j==i else 0 for j in range(3)]) for i in range(3) for sg in (1,-1)]
    lattice=mp.matrix([[2,-1,-1],[1,2,-1],[1,-1,2]])
    omega=mp.matrix([-rt,1,1])
    def cross(a,b):return mp.matrix([a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]])
    vel=[cross(omega,p) for p in points]
    norm2=sum(mp.fdot(v,v) for v in vel)
    def dot(a,b):return mp.fdot(a,b)
    def norm(a):return mp.sqrt(dot(a,a))
    def U(h,kappa):
        value=mp.mpf(0)
        for i in range(6):
            for j in range(i+1,6):
                if i//2==j//2:continue
                d=points[j]-points[i]+h*(vel[j]-vel[i])
                value+=(norm(d)-rt)**2/2
        for axis in range(3):
            i,j=2*axis,2*axis+1
            v=points[j]-points[i]+lattice[axis,:].T-h*(vel[j]+vel[i])
            value+=(norm(v)-rt)**2/2
        for i in range(6):
            j=i^1;sign=1 if i%2==0 else -1
            b=points[j]-points[i]+sign*lattice[i//2,:].T-h*(vel[j]+vel[i])
            for k in range(6):
                if i//2==k//2:continue
                d=points[k]-points[i]+h*(vel[k]-vel[i])
                cosine=dot(b,d)/(norm(b)*norm(d))
                value+=kappa*(cosine+1/rt)**2/2
        return value
    kappa=mp.mpf('1e-6')
    predicted=kappa*(mp.mpf(1)/4-3*rt/8)
    rows=[]
    for hs in ['1e-5','1e-7','1e-9']:
        h=mp.mpf(hs)
        curvature=(U(h,kappa)+U(-h,kappa)-2*U(0,kappa))/(h*h*norm2)
        rows.append({'amplitude':hs,'curvature':str(curvature),'predicted':str(predicted),
                     'relative_error':float(abs(curvature-predicted)/abs(predicted))})
    assert rows[-1]['relative_error']<1e-10
    return rows


def main():
    checkpoint=json.loads((AUTHOR/'checkpoint.json').read_text())
    hashes={}
    for item in checkpoint['files']:
        actual=hashlib.sha256((AUTHOR/item['path']).read_bytes()).hexdigest()
        assert actual==item['sha256'],item['path']
        hashes[item['path']]=actual
    result={'source_hashes':hashes,'exact':exact_checks(),'independent_finite_energy':high_precision_endpoint()}
    (HERE/'independent-results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
