"""Numerical check of Theorem 1 (contact framework of the truncated n-cubic honeycomb).

Independent of the hand proof: builds the periodic (Bloch) rigidity matrix and the
spring-network elastic tensor directly from the geometry, for n = 2..8, and checks
  (a) zero modes / self-stresses at generic q, at q = 0, and on planes q_i = 0;
  (b) the cross-polytope is infinitesimally rigid for n >= 3 with n(n-3)/2 stresses;
  (c) the explicit stress basis (equal tension w_ij on the 4 bars of the (i,j) square,
      sum_j w_ij = 0) spans the self-stresses at generic q;
  (d) the explicit q = 0 stresses (w, tau) with tau_i = sqrt2 sum_j w_ij (tensions);
  (e) generic-q zero modes are exactly the cluster rotations u_{i,s} = s c Omega e_i;
  (f) C11, C12, all shear moduli = 0, nu100 = 1/(n^2-n-1), with FULL internal relaxation
      (no symmetry restriction), for the incentre c = 2 - sqrt2 and for other c.
Run:  python rigidity_nd.py   (the certificate driver is verify.py)
"""
import itertools
import numpy as np

rng = np.random.default_rng(1)
TOL = 1e-9


def framework(n, c):
    """dots (i,s) at s*c*e_i in a cell of side 2 centred at 0; lattice 2Z^n.
    bars: list of (a, b, L) meaning dot a (cell 0) -- dot b (cell L)."""
    verts = [(i, s) for i in range(n) for s in (1, -1)]
    idx = {v: k for k, v in enumerate(verts)}
    pos = np.array([s * c * np.eye(n)[i] for i, s in verts])
    bars, kind = [], []
    for (i, s), (j, t) in itertools.combinations(verts, 2):
        if i != j:
            bars.append((idx[(i, s)], idx[(j, t)], np.zeros(n)))
            kind.append(("intra", i, j, s, t))
    for i in range(n):
        bars.append((idx[(i, 1)], idx[(i, -1)], 2 * np.eye(n)[i]))
        kind.append(("cross", i))
    return verts, pos, bars, kind


def rigidity(n, c, q):
    verts, pos, bars, kind = framework(n, c)
    V = len(verts)
    R = np.zeros((len(bars), n * V), complex)
    for r, (a, b, L) in enumerate(bars):
        d = pos[b] + L - pos[a]
        d = d / np.linalg.norm(d)
        R[r, n * a:n * a + n] -= d
        R[r, n * b:n * b + n] += d * np.exp(1j * q @ L)
    return R


def rank(M):
    if M.size == 0:
        return 0
    s = np.linalg.svd(M, compute_uv=False)
    return int((s > TOL * max(1, s[0])).sum())


def nullspace(M):
    u, s, vh = np.linalg.svd(M)
    r = int((s > TOL * max(1, s[0])).sum())
    return vh[r:].conj().T


def check_n(n, c=2 - np.sqrt(2)):
    out = {}
    verts, pos, bars, kind = framework(n, c)
    V, E = len(verts), len(bars)
    out["dofs"], out["bars"] = n * V, E
    # bar lengths
    lens = [np.linalg.norm(pos[b] + L - pos[a]) for a, b, L in bars]
    out["lengths"] = (min(lens), max(lens))
    # (a) generic q
    zs = []
    for _ in range(3):
        q = rng.uniform(-np.pi, np.pi, n)
        R = rigidity(n, c, q)
        rk = rank(R)
        zs.append((n * V - rk, E - rk))
    out["generic"] = zs
    R0 = rigidity(n, c, np.zeros(n))
    rk0 = rank(R0)
    out["q0"] = (n * V - rk0, E - rk0)
    planes = []
    for k in range(n + 1):
        q = rng.uniform(0.3, 2.8, n)
        q[:k] = 0.0
        R = rigidity(n, c, q)
        rk = rank(R)
        planes.append((k, n * V - rk, E - rk))
    out["planes"] = planes
    # (b) cross-polytope alone
    intra = [r for r, kk in enumerate(kind) if kk[0] == "intra"]
    Rc = R0[intra].real
    rkc = rank(Rc)
    out["crosspolytope"] = (len(intra), rkc, n * V - rkc, len(intra) - rkc)
    # (c) explicit generic-q stress basis: w on edges of K_n with zero vertex sums
    pairs = list(itertools.combinations(range(n), 2))
    Inc = np.zeros((n, len(pairs)))
    for m, (i, j) in enumerate(pairs):
        Inc[i, m] = Inc[j, m] = 1
    W = nullspace(Inc).real  # columns: w vectors, dim C(n,2)-n for n>=3
    S = np.zeros((E, W.shape[1]))
    for col in range(W.shape[1]):
        for r, kk in enumerate(kind):
            if kk[0] == "intra":
                S[r, col] = W[pairs.index((kk[1], kk[2])), col]
    q = rng.uniform(-np.pi, np.pi, n)
    Rq = rigidity(n, c, q)
    resid = np.abs(S.T @ Rq).max() if S.size else 0.0
    left = nullspace(Rq.conj().T)  # all self-stresses at q
    span_ok = rank(np.hstack([left, S.astype(complex)])) == left.shape[1] == W.shape[1] if n >= 3 else None
    out["stress_basis"] = (W.shape[1], resid, span_ok)
    # (d) q = 0 stresses: arbitrary w, tau_i = sqrt2 sum_j w_ij (tensions)  (tension convention: rows of R)
    wq = rng.normal(size=len(pairs))
    s0 = np.zeros(E)
    for r, kk in enumerate(kind):
        if kk[0] == "intra":
            s0[r] = wq[pairs.index((kk[1], kk[2]))]
    # force balance with unit-direction rows: intra bar length c*sqrt2, rows use unit vectors.
    # derive tau from equilibrium numerically and compare with formula sqrt2*sum_j w_ij
    cross = [r for r, kk in enumerate(kind) if kk[0] == "cross"]
    A = R0[cross].real.T
    b = -(R0.real.T @ s0)
    tau, *_ = np.linalg.lstsq(A, b, rcond=None)
    s0[cross] = tau
    out["q0_stress"] = (np.abs(R0.real.T @ s0).max(),
                        np.abs(tau - np.sqrt(2) * (Inc @ wq)).max())
    # (e) generic-q zero modes = rotations of each cross-polytope
    N = nullspace(Rq)
    rots = []
    for i, j in pairs:
        Om = np.zeros((n, n)); Om[i, j], Om[j, i] = 1, -1
        u = np.concatenate([Om @ pos[k] for k in range(V)])
        rots.append(u)
    Rot = np.array(rots).T.astype(complex)
    out["modes_are_rotations"] = (np.abs(Rq @ Rot).max(),
                                  rank(np.hstack([N, Rot])) == N.shape[1] == Rot.shape[1])
    return out


def elastic(n, c):
    """Full relaxation: W(eps) = min_w (1/2V) sum_b (d_b.(eps l_b + w_b - w_a))^2, V = 2^n."""
    verts, pos, bars, kind = framework(n, c)
    V = len(verts)
    sym = [(i, i) for i in range(n)] + list(itertools.combinations(range(n), 2))
    Aeps = np.zeros((len(bars), len(sym)))
    Bw = np.zeros((len(bars), n * V))
    for r, (a, b, L) in enumerate(bars):
        l = pos[b] + L - pos[a]
        d = l / np.linalg.norm(l)
        for m, (i, j) in enumerate(sym):
            Eij = np.zeros((n, n)); Eij[i, j] = Eij[j, i] = 1.0  # tensor strain eps_ij=eps_ji=x
            if i == j:
                Eij[i, i] = 1.0
            Aeps[r, m] = d @ Eij @ l
        Bw[r, n * a:n * a + n] -= d
        Bw[r, n * b:n * b + n] += d
    # Schur complement: K = A^T (I - B B^+) A
    P = np.eye(len(bars)) - Bw @ np.linalg.pinv(Bw)
    K = Aeps.T @ P @ Aeps / 2 ** n  # W = (1/2) x^T K x in variables (eps_ii, eps_ij)
    C11 = K[0, 0]
    C12 = K[0, 1]
    shear = K[n:, n:]
    coupling = np.abs(K[:n, n:]).max()
    Cn = K[:n, :n]
    # nu100 from compliance of the normal block (uniaxial stress along e_1)
    # uniaxial stress: eps_1 = 1, sigma_j = 0 for j >= 2 (valid also for n = 2, where Cn is singular)
    eperp = np.linalg.solve(Cn[1:, 1:], -Cn[1:, 0])
    nu = -eperp[0]
    return C11, C12, np.abs(shear).max(), coupling, nu, Cn


if __name__ == "__main__":
    allok = True
    print("n | dofs bars | generic (zero,SS) | q=0 | crosspoly(bars,rank,zero,SS) | predicted")
    for n in range(2, 9):
        o = check_n(n)
        pred = (n * (n - 1) // 2, n * (n - 3) // 2)
        g = o["generic"]
        ok = all(z == pred for z in g) if n >= 3 else all(z == (2, 0) for z in g)
        okq0 = o["q0"] == (n * (n + 1) // 2, n * (n - 1) // 2) if n >= 3 else True
        okcp = o["crosspolytope"][2] == n * (n + 1) // 2 if n >= 3 else True
        okpl = all(z == pred[0] + k for k, z, s in o["planes"]) if n >= 3 else True
        allok &= ok and okq0 and okcp and okpl
        print(f"{n} | {o['dofs']} {o['bars']} | {g[0]} | {o['q0']} | {o['crosspolytope']} | {pred if n>=3 else '(2,0) for n=2'}"
              f" | lengths {o['lengths'][0]:.12f}..{o['lengths'][1]:.12f} | ok={ok and okq0 and okcp and okpl}")
        print(f"    planes q_1..q_k=0 (k, zero, SS): {o['planes']}")
        print(f"    explicit stress basis dim={o['stress_basis'][0]}, |S^T R(q)|max={o['stress_basis'][1]:.1e}, spans all SS: {o['stress_basis'][2]}")
        print(f"    q=0 stress (w,tau): equilibrium residual {o['q0_stress'][0]:.1e}, |tau - sqrt2*Inc w| = {o['q0_stress'][1]:.1e}")
        print(f"    generic-q zero modes = cluster rotations: residual {o['modes_are_rotations'][0]:.1e}, same span: {o['modes_are_rotations'][1]}")
        if n >= 3:
            allok &= bool(o["stress_basis"][2]) and o["modes_are_rotations"][1]
    print("\nElastic constants, unit springs, full internal relaxation, V = 2^n")
    print("n | c | C11 | 4(n-1)^2/(n^2 2^n) | C12 | 4/(n^2 2^n) | max|shear| | max|normal-shear| | nu100 | 1/(n^2-n-1)")
    for n in range(2, 9):
        for c in (2 - np.sqrt(2), 0.3, 0.75):
            C11, C12, sh, cp, nu, Cn = elastic(n, c)
            p11 = 4 * (n - 1) ** 2 / (n ** 2 * 2 ** n)
            p12 = 4 / (n ** 2 * 2 ** n)
            pn = 1 / (n * n - n - 1)
            ok = abs(C11 - p11) < 1e-12 and abs(C12 - p12) < 1e-12 and sh < 1e-12 and cp < 1e-12 and abs(nu - pn) < 1e-10
            cubic = np.allclose(Cn, p12 * np.ones((n, n)) + (p11 - p12) * np.eye(n), atol=1e-12)
            allok &= ok and cubic
            print(f"{n} | {c:.4f} | {C11:.12f} | {p11:.12f} | {C12:.12f} | {p12:.12f} | {sh:.1e} | {cp:.1e} | {nu:.10f} | {pn:.10f} | ok={ok and cubic}")
    print("\nALL CHECKS PASSED" if allok else "\nSOME CHECK FAILED")
