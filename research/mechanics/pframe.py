"""Periodic bar-framework rigidity tools (n dimensions).

A framework: lattice matrix A (columns = lattice vectors), sites X (N x n),
bonds [(i, j, m)] meaning site i in cell 0 joined to site j in cell m (int vector).
"""
import numpy as np
import itertools

TOL = 1e-8


def bond_vec(A, X, b):
    i, j, m = b
    return X[j] + A @ np.asarray(m, float) - X[i]


def rigidity(A, X, bonds, q=None):
    """Compatibility (rigidity) matrix R(q): rows = bonds, cols = n*N.
    Bloch convention u(site j, cell m) = u_j exp(i q.(A m))."""
    N, n = X.shape
    q = np.zeros(n) if q is None else np.asarray(q, float)
    R = np.zeros((len(bonds), n * N), complex)
    for r, b in enumerate(bonds):
        i, j, m = b
        d = bond_vec(A, X, b)
        e = d / np.linalg.norm(d)
        ph = np.exp(1j * q @ (A @ np.asarray(m, float)))
        R[r, n * i:n * i + n] -= e
        R[r, n * j:n * j + n] += e * ph
    return R


def rank(M, tol=TOL):
    if M.size == 0:
        return 0
    s = np.linalg.svd(M, compute_uv=False)
    return int((s > tol * max(1.0, s[0])).sum())


def counts(A, X, bonds, q=None):
    R = rigidity(A, X, bonds, q)
    r = rank(R)
    N0 = R.shape[1] - r          # zero modes (incl. trivial translations at q=0)
    S = R.shape[0] - r           # states of self-stress
    return N0, S


def strain_basis(n):
    """basis of symmetric n x n strains"""
    B = []
    for a in range(n):
        for b in range(a, n):
            E = np.zeros((n, n))
            E[a, b] = E[b, a] = 1.0 if a != b else 1.0
            if a != b:
                E /= 2.0
            B.append(((a, b), E))
    return B


def flexible_rigidity(A, X, bonds):
    """q=0 rigidity matrix with extra columns for affine strain (symmetric eps).
    bond extension = e.(u_j - u_i) + e.(eps d)."""
    N, n = X.shape
    R0 = rigidity(A, X, bonds).real
    SB = strain_basis(n)
    C = np.zeros((len(bonds), len(SB)))
    for r, b in enumerate(bonds):
        d = bond_vec(A, X, b)
        e = d / np.linalg.norm(d)
        for k, (_, E) in enumerate(SB):
            C[r, k] = e @ E @ d
    return np.hstack([R0, C]), SB


def flexible_analysis(A, X, bonds):
    """Return (#non-trivial mechanisms incl. strain, dim of achievable strain space,
    basis of achievable strains as n x n matrices, #self-stresses)."""
    N, n = X.shape
    M, SB = flexible_rigidity(A, X, bonds)
    U, s, Vt = np.linalg.svd(M)
    r = int((s > TOL * max(1.0, s[0])).sum())
    null = Vt[r:].T                      # columns = zero modes
    nmech = null.shape[1] - n            # remove rigid translations
    strains = null[n * N:, :]            # strain components of the zero modes
    sr = rank(strains)
    # orthonormal basis for achievable strain subspace
    if sr > 0:
        Us, ss, _ = np.linalg.svd(strains)
        Sbasis = Us[:, :sr]
    else:
        Sbasis = np.zeros((len(SB), 0))
    mats = []
    for c in range(Sbasis.shape[1]):
        E = np.zeros((n, n))
        for k, (_, Eb) in enumerate(SB):
            E += Sbasis[k, c] * Eb
        mats.append(E)
    return nmech, sr, mats, M.shape[0] - r


def supercell(A, X, bonds, reps):
    """Build supercell with reps = tuple of ints along each lattice vector."""
    N, n = X.shape
    cells = list(itertools.product(*[range(k) for k in reps]))
    idx = {c: k for k, c in enumerate(cells)}
    reps = np.array(reps)
    X2 = np.vstack([X + A @ np.array(c, float) for c in cells])
    A2 = A * reps[None, :]
    B2 = []
    for c in cells:
        for (i, j, m) in bonds:
            t = np.array(c) + np.array(m)
            cm = np.mod(t, reps)
            M = (t - cm) // reps
            B2.append((idx[c] * N + i, idx[tuple(cm)] * N + j, tuple(int(v) for v in M)))
    return A2, X2, B2


def random_q(n, rng):
    return rng.uniform(-np.pi, np.pi, n)
