"""Test rule on face-transitive (isohedral) P whose duals are the 'missing' coordination polyhedra."""
import numpy as np
from polytopes import analyse, antiprism, prism, perms_signs, classic_ratio, SOLIDS

def bipyramid(n, h):
    return [(np.cos(2*np.pi*k/n), np.sin(2*np.pi*k/n), 0) for k in range(n)] + [(0,0,h),(0,0,-h)]

def polar_dual(V):
    from polytopes import facets
    return [n/h for n,h,idx in facets(np.asarray(V,float))]

print("bipyramids (dual = prisms):")
for n in (3,4,5,6):
    for h in (0.6,1.0,1.7):
        a = analyse(f"{n}-bipyramid h={h}", bipyramid(n,h))
        print(f"  n={n} h={h}: balls={a['n_balls']} eqR={a['equal_R']} touch={a['touch']} ratio={a['ratio']:.4f}  uniform {n}-prism={classic_ratio(prism(n)):.4f}")
print("trapezohedra (dual of antiprisms):")
for n in (3,4,5):
    for s in (0.8,1.0,1.3):
        A = np.array(antiprism(n)); A[:,2]*=s
        T = polar_dual(A)
        a = analyse("trap", T)
        print(f"  n={n} stretch={s}: balls={a['n_balls']} eqR={a['equal_R']} ins_all={a['ins_all']} touch={a['touch']} ratio={a['ratio']:.4f}  uniform antiprism={classic_ratio(antiprism(n)):.4f}")
# Catalan solids (duals of Archimedean): face-transitive, mirror walls -> expect Archimedean radius ratios
arch = {
 "cuboctahedron": perms_signs((1,1,0)),
 "truncated octahedron": perms_signs((0,1,2)),
 "truncated cube": perms_signs((2**0.5-1,1,1)),
 "rhombicuboctahedron": perms_signs((1,1,1+2**0.5)),
 "truncated tetrahedron": [v for v in perms_signs((3,1,1)) if np.prod(np.sign(v))>0],
 "icosidodecahedron": None,
}
phi=(1+5**0.5)/2
arch["icosidodecahedron"]=np.vstack([perms_signs((0,0,phi)) , perms_signs((0.5,phi/2,phi**2/2),even_only=True)])
print("Catalan P (dual Archimedean Q):")
for k,Q in arch.items():
    a = analyse("cat", polar_dual(Q))
    print(f"  Q={k:22s} balls={a['n_balls']:3d} eqR={a['equal_R']} touch={a['touch']} ratio={a['ratio']:.4f}  classic(Q)={classic_ratio(Q):.4f}")
