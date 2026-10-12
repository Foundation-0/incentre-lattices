"""Effective (linear, central-force spring) elastic tensor of periodic frameworks,
including non-affine relaxation. Springs k_b = 1 (or 1/L_b if 'EA')."""
import numpy as np
from pframe import *
from contact_nd import contact_framework, TSTAR


def voigt_strains(n):
    """strain basis in Voigt order: normals then shears (engineering shear gamma=2 eps)."""
    out = []
    for a in range(n):
        E = np.zeros((n, n)); E[a, a] = 1; out.append(E)
    for a in range(n):
        for b in range(a + 1, n):
            E = np.zeros((n, n)); E[a, b] = E[b, a] = 0.5; out.append(E)
    return out


def elastic_tensor(A, X, bonds, kmode="unit"):
    N, n = X.shape
    R = rigidity(A, X, bonds).real
    V = abs(np.linalg.det(A))
    L = np.array([np.linalg.norm(bond_vec(A, X, b)) for b in bonds])
    k = np.ones(len(bonds)) if kmode == "unit" else 1.0 / L
    K = np.diag(k)
    VS = voigt_strains(n)
    Cmat = np.zeros((len(bonds), len(VS)))
    for r, b in enumerate(bonds):
        d = bond_vec(A, X, b); e = d / np.linalg.norm(d)
        for c, E in enumerate(VS):
            Cmat[r, c] = e @ E @ d
    # minimise 1/2 |K^1/2 (R u + C eps)|^2 over u: projected stiffness
    H = R.T @ K @ R
    G = R.T @ K @ Cmat
    u = -np.linalg.pinv(H, rcond=1e-10) @ G
    res = R @ u + Cmat
    Cv = res.T @ K @ res / V
    return Cv


if __name__ == "__main__":
    np.set_printoptions(precision=5, suppress=True)
    for n in (2, 3, 4):
        A, X, B = contact_framework(n)
        A2, X2, B2 = supercell(A, X, B, (2,) * n) if n < 4 else (A, X, B)
        C = elastic_tensor(A2, X2, B2)
        print(f"n={n} Voigt elastic matrix (unit springs, Wythoff point):\n", C)
    print("\nt-scan, n=3, unit springs and EA springs: C11, C12, C44, nu_100 = C12/(C11+C12)")
    for kmode in ("unit", "EA"):
        for t in np.linspace(0.1, 0.95, 18).tolist() + [TSTAR]:
            A, X, B = contact_framework(3, t)
            C = elastic_tensor(A, X, B, kmode)
            print(kmode, round(t, 4), round(C[0, 0], 5), round(C[0, 1], 5), round(C[3, 3], 8), round(C[0, 1] / (C[0, 0] + C[0, 1]), 5))
