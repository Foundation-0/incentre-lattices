#!/usr/bin/env python3
"""search.py -- the ONLY file the agent edits.  Generates out/path.npz for verify.py.

Pipeline:  [stage A] analytic cubic [111] checkerboard fold to FOLD_END rad (baseline, m even)
           [stage B] optional numerical continuation (free, possibly sheared lattice):
                     LP picks the first-order motion that keeps every bar length, lets contacts
                     only separate, and shrinks det L fastest; Gauss-Newton retracts onto the
                     bar manifold (touching balls that would overlap become temporary equality
                     constraints = sliding contact); random kicks escape stationary points.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")   # ~2x faster on small dense lstsq
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
import time
import numpy as np
from scipy.optimize import linprog
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import verify
from verify import R, K, D, canonical_state, bar_graph, close_pairs, min_pair_distance

DESC = "steer: from canonical, LP-steer toward explicit m=1 fcc realisation T=[[2,-1,-1],[1,2,-1],[1,-1,2]] (metric+fractional objective), then density descent (skipped if reached); m=1"

# ----------------------------------------------------------------- knobs
M = 1                    # supercell (cubic checkerboard fold needs M even)
FOLD_END = 0.0          # stage A end angle (rad); 0 = skip stage A (start from canonical state)
CONTINUE = True          # stage B on/off
SEED = 1
TIME_BUDGET = 240.0      # seconds of wall clock for the whole script (run.sh kills at 300)
HX = 0.006               # max dot displacement per step (verify limit 0.01)
HL = 0.006               # max lattice Frobenius change per step (verify limit 0.01)
MAX_STEPS = 4000
START_KICKS = 3          # seeded random cone kicks at the start of stage B (makes SEED matter)
NOPROG_STEPS = 40        # descent steps without a new best -> treat as stationary (cycle breaker)
STALL_KICKS = 25         # stop after this many consecutive kicks without new best density
CONTACT_MARGIN = 0.03    # relative distance within which non-bar pairs become LP constraints
TRIM_TO_BEST = True
TARGET = "fcc"           # None or "fcc": steer (stage S) toward an explicit fcc realisation before stage B
TARGET_T = [[2, -1, -1], [1, 2, -1], [1, -1, 2]]   # link lattice in fcc units (rows t_a = 2 e_a + n_a), det 12
W_F = 10.0               # weight of fractional-coordinate mismatch vs metric mismatch
STEER_TOL = 1e-10
SKIP_B_IF_REACHED = True  # target is an exact fcc packing (first-order jammed at m=2): skip density descent
STEER_FRAC = 0.6         # fraction of the remaining time budget given to stage S      # cut path at its max-density state
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "path.npz")
T0 = time.time()


def log(*a):
    print("[%6.1fs]" % (time.time() - T0), *a, flush=True)


# ----------------------------------------------------------------- stage A
def rot(n, th):
    n = np.asarray(n, float) / np.linalg.norm(n)
    Kx = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    return np.eye(3) + np.sin(th) * Kx + (1 - np.cos(th)) * Kx @ Kx


def cubic_state(m, th):
    c = np.cos(th)
    q = 2 * K * (1 - c) / 3
    a = (2 * K * c + q + np.sqrt(4 * R ** 2 - 2 * q ** 2)) / 2
    L = 2 * a * m * np.eye(3)
    X = np.zeros((6 * m ** 3, 3))
    for ix in range(m):
        for iy in range(m):
            for iz in range(m):
                cc = (ix * m + iy) * m + iz
                Rm = rot([1, 1, 1], th if (ix + iy + iz) % 2 == 0 else -th)
                ctr = 2 * a * np.array([ix, iy, iz]) + a
                for ax in range(3):
                    for sg in range(2):
                        e = np.zeros(3)
                        e[ax] = 1 if sg == 0 else -1
                        X[6 * cc + 2 * ax + sg] = ctr + K * Rm @ e
    return L, X


def stage_a(m, th_end):
    L, X = cubic_state(m, 0.0)
    Ls, Xs = [L], [X]
    th = 0.0
    while th < th_end - 1e-15:
        th2 = min(th + 0.004 / m, th_end)   # lattice change scales with m
        L2, X2 = cubic_state(m, th2)
        Ls.append(L2)
        Xs.append(X2)
        th = th2
    return Ls, Xs


# ----------------------------------------------------------------- stage B engine
class Engine:
    def __init__(self, m):
        self.m = m
        self.N = 6 * m ** 3
        self.n = 3 * self.N + 9
        self.I, self.J, self.S = bar_graph(m)
        self.barkeys = {verify._canon_key(i, j, s) for i, j, s in zip(self.I, self.J, self.S)}
        self.rng = np.random.default_rng(SEED)

    def rowsparse(self, X, L, I, J, S, normalize=True):
        """rows d(len)/d(u) for pairs (I,J,S); u = [v (3N), Ldot (9)]; returns csr, lengths, vectors"""
        N = self.N
        r = X[J] + S @ L - X[I]
        ln = np.linalg.norm(r, axis=1)
        rh = r / ln[:, None] if normalize else r
        B = len(I)
        rows = np.repeat(np.arange(B), 3)
        cols_i = (3 * I[:, None] + np.arange(3)).ravel()
        cols_j = (3 * J[:, None] + np.arange(3)).ravel()
        vals = rh.ravel()
        # Ldot entries: S[b,a] * rh[b,c] at column 3N + 3a + c
        ra = np.repeat(np.arange(B), 9)
        ca = 3 * N + np.tile(np.arange(9), B)
        va = (S[:, :, None] * rh[:, None, :]).reshape(B, 9).ravel()
        A = sp.coo_matrix((np.concatenate([-vals, vals, va]),
                           (np.concatenate([rows, rows, ra]), np.concatenate([cols_i, cols_j, ca]))),
                          shape=(B, self.n)).tocsr()
        return A, ln, r

    def near_pairs(self, X, L, rel):
        i, j, s, d = close_pairs(X, L, D * (1 + rel))
        keep = np.array([verify._canon_key(a, b, c) not in self.barkeys for a, b, c in zip(i, j, s)], bool)
        if len(i) == 0:
            return i, j, s, d
        return i[keep], j[keep], s[keep], d[keep]

    def gauge_rows(self, X, L):
        N = self.N
        G = np.zeros((6, self.n))
        for k in range(3):
            G[k, k:3 * N:3] = 1.0
        # angular momentum  sum x_i x v_i + sum_k l_k x ldot_k
        for i in range(N):
            x = X[i]
            for a in range(3):          # v_i[a] contributes (x x e_a)
                e = np.zeros(3); e[a] = 1
                G[3:6, 3 * i + a] = np.cross(x, e)
        for k in range(3):
            for a in range(3):
                e = np.zeros(3); e[a] = 1
                G[3:6, 3 * N + 3 * k + a] = np.cross(L[k], e)
        return G

    def direction(self, X, L, cost_mode="descent", rcost=None, h_try=0.005):
        """LP: returns (u, predicted dlndet) ; u = [v, Ldot] scaled so |v|inf<=1, |Ldot|inf<=m"""
        N, n, m = self.N, self.n, self.m
        Ab, _, _ = self.rowsparse(X, L, self.I, self.J, self.S)
        ci, cj, cs, cd = self.near_pairs(X, L, CONTACT_MARGIN)
        if len(ci):
            Ac, _, _ = self.rowsparse(X, L, ci, cj, cs)
            rhs = -0.5 * (cd - D) / h_try        # allowed approach speed (<=0 when touching: 0)
        G = sp.csr_matrix(self.gauge_rows(X, L))
        Aeq = sp.vstack([Ab, G]).tocsr()
        beq = np.zeros(Aeq.shape[0])
        Linv = np.linalg.inv(L)
        c = np.zeros(n)
        c[3 * N:] = Linv.T.ravel()       # d ln det L = c.u  (we minimise: density up)
        cost = c if cost_mode == "descent" else rcost
        eps = 2e-3 * np.abs(c[3 * N:]).max()
        Aeq2 = sp.hstack([Aeq, -Aeq]).tocsr()
        cc = np.concatenate([cost + eps, -cost + eps])
        ub = np.concatenate([np.ones(3 * N), m * np.ones(9)])
        bounds = [(0, b) for b in np.concatenate([ub, ub])]
        if len(ci):
            Aub = sp.hstack([-Ac, Ac]).tocsr()
            bub = np.maximum(-rhs, 0.0) * 0 + np.where(cd - D > 0, 0.5 * (cd - D) / h_try, 0.0)
            lp = linprog(cc, A_ub=Aub, b_ub=bub, A_eq=Aeq2, b_eq=beq, bounds=bounds, method="highs")
        else:
            lp = linprog(cc, A_eq=Aeq2, b_eq=beq, bounds=bounds, method="highs")
        if lp.status != 0:
            return None, 0.0
        u = lp.x[:n] - lp.x[n:]
        return u, float(c @ u)

    def retract(self, X, L, max_iter=40):
        """Gauss-Newton back to |bar| = D (and sticky equality for any pair that would overlap)."""
        N = self.N
        ei = np.zeros(0, int); ej = np.zeros(0, int); es = np.zeros((0, 3))
        for it in range(max_iter):
            # sticky extras
            pi, pj, ps, pd = close_pairs(X, L, D * (1 + 1e-9))
            if len(pi):
                kp = np.array([verify._canon_key(a, b, c) not in self.barkeys for a, b, c in zip(pi, pj, ps)], bool)
                pi, pj, ps, pd = pi[kp], pj[kp], ps[kp], pd[kp]
                viol = pd < D * (1 - 1e-13)
                for a, b, c in zip(pi[viol], pj[viol], ps[viol]):
                    ei = np.append(ei, a); ej = np.append(ej, b); es = np.vstack([es, c])
            I = np.concatenate([self.I, ei]).astype(int)
            J = np.concatenate([self.J, ej]).astype(int)
            S = np.vstack([self.S, es]) if len(ei) else self.S
            A, ln, _ = self.rowsparse(X, L, I, J, S)
            res = ln - D
            if np.abs(res[:len(self.I)]).max() < 1e-13 * D and (len(pi) == 0 or not viol.any()) \
                    and np.abs(res).max() < 1e-12 * D:
                return X, L, True
            dx = np.linalg.lstsq(A.toarray(), -res, rcond=1e-10)[0]
            X = X + dx[:3 * N].reshape(N, 3)
            L = L + dx[3 * N:].reshape(3, 3)
        return X, L, False

    def step(self, X, L, u, h):
        N = self.N
        X1 = X + h * u[:3 * N].reshape(N, 3)
        L1 = L + h * u[3 * N:].reshape(3, 3)
        return self.retract(X1, L1)


def stage_b(m, L, X, deadline, obj=None, start_kicks=None, tag="B"):
    """obj None: shrink det L (density).  obj(X,L)->(f, grad over u): drive f down (steering)."""
    if start_kicks is None:
        start_kicks = START_KICKS
    val = (lambda X, L: np.linalg.det(L)) if obj is None else (lambda X, L: obj(X, L)[0])
    eng = Engine(m)
    N = eng.N
    Ls, Xs = [L.copy()], [X.copy()]
    capf = 1.0
    best = val(X, L)
    best_k = 0
    kicks = 0
    logt = time.time()
    for k in range(MAX_STEPS):
        if time.time() > deadline:
            log("stage %s: deadline" % tag + ""); break
        if obj is None:
            u, rate = eng.direction(X, L)
        else:
            fv, g = obj(X, L)
            if fv < STEER_TOL:
                log("stage %s: target reached f=%.3e" % (tag, fv)); break
            g = g / max(np.abs(g).max(), 1e-300)
            u, _ = eng.direction(X, L, "kick", g)
            rate = float(g @ u) if u is not None else 0.0
        mode = "descent"
        if k < start_kicks:
            u, _ = eng.direction(X, L, "kick", eng.rng.normal(size=eng.n))
            mode = "kick"
            if u is None:
                continue
        elif u is None or rate > -1e-7 * 0.01 or k - best_k > NOPROG_STEPS:
            # stationary at first order: random kick inside the feasible cone
            kicks += 1
            if kicks > STALL_KICKS:
                log("stage %s: stalled (%d kicks without progress) -> jammed or stuck" % (tag, STALL_KICKS)); break
            rc = eng.rng.normal(size=eng.n)
            u, _ = eng.direction(X, L, "kick", rc)
            best_k = k                     # give the kick NOPROG_STEPS steps to pay off
            if u is None:
                log("stage B: kick LP failed"); break
            mode = "kick"
        sx = np.abs(u[:3 * N]).max()
        sl = np.linalg.norm(u[3 * N:])
        if max(sx, sl) < 1e-12:
            kicks += 1
            if kicks > STALL_KICKS:
                break
            continue
        h = capf * min(HX / max(sx, 1e-12), HL / max(sl, 1e-12))
        ok = False
        for attempt in range(12):
            X1, L1, conv = eng.step(X, L, u, h)
            if conv and np.linalg.det(L1) > 0:
                dx = np.linalg.norm(X1 - X, axis=1).max()
                dl = np.linalg.norm(L1 - L)
                if dx <= 0.0095 and dl <= 0.0095:
                    md = min_pair_distance(X1, L1)
                    mm = min_pair_distance(0.5 * (X + X1), 0.5 * (L + L1))
                    if md >= D * (1 - 1e-11) and mm >= D * (1 - 5e-5):
                        ok = True
                        break
            h *= 0.5
            capf *= 0.5
        if not ok:
            kicks += 1
            log("step failed at k=%d (mode %s)" % (k, mode))
            if kicks > STALL_KICKS:
                break
            continue
        capf = min(1.0, capf * 1.5)
        X, L = X1, L1
        Ls.append(L.copy()); Xs.append(X.copy())
        det = np.linalg.det(L)
        v = val(X, L)
        if v < best - 1e-9 * (1 if obj is None else 1e-3):
            best = v; best_k = k; kicks = 0
        if time.time() - logt > 10:
            logt = time.time()
            dens = N * 4 / 3 * np.pi * R ** 3 / det
            log("%s k=%d density=%.5f best=%.6g mode=%s h=%.2e" % (tag, k, dens, best, mode, h))
    return Ls, Xs


# ----------------------------------------------------------------- steering target
def fcc_target(m, T):
    """explicit fcc realisation: octahedron of cell idx centred at (idx@T) s, dot (a,+-) at +-e_a s, s = D/sqrt2."""
    T = np.asarray(T, float)
    s = D / np.sqrt(2)
    X = np.zeros((6 * m ** 3, 3))
    for ix in range(m):
        for iy in range(m):
            for iz in range(m):
                c = (ix * m + iy) * m + iz
                h = np.array([ix, iy, iz]) @ T
                for a in range(3):
                    e = np.zeros(3); e[a] = 1
                    X[6 * c + 2 * a] = (h + e) * s
                    X[6 * c + 2 * a + 1] = (h - e) * s
    return m * T * s, X


def make_objective(m, Lt, Xt):
    Gt = Lt @ Lt.T
    Ft = Xt @ np.linalg.inv(Lt)
    Ft = Ft - Ft.mean(0)
    w = W_F

    def obj(X, L):
        Li = np.linalg.inv(L)
        dG = (L @ L.T - Gt) / m ** 2
        F = X @ Li
        E = F - F.mean(0) - Ft
        f = float((dG ** 2).sum() + w * (E ** 2).sum())
        gL = 4 * dG @ L / m ** 2 - 2 * w * F.T @ E @ Li.T
        gX = 2 * w * E @ Li.T
        return f, np.concatenate([gX.ravel(), gL.ravel()])
    return obj


# ----------------------------------------------------------------- main
def main():
    m = M
    if FOLD_END > 0:
        Ls, Xs = stage_a(m, FOLD_END)
    else:
        L0, X0 = canonical_state(m)
        Ls, Xs = [L0], [X0]
    log("stage A done: %d states, density %.5f" % (len(Ls), 6 * m ** 3 * 4 / 3 * np.pi * R ** 3 / np.linalg.det(Ls[-1])))
    if TARGET == "fcc":
        Lt, Xt = fcc_target(m, TARGET_T)
        log("fcc target: det %.6f density %.6f, bars ok %s" % (np.linalg.det(Lt), 6 * m ** 3 * 4 / 3 * np.pi * R ** 3 / np.linalg.det(Lt),
            np.abs(np.linalg.norm(verify.bar_vectors(Xt, Lt, *bar_graph(m)), axis=1) / D - 1).max() < 1e-12))
        obj = make_objective(m, Lt, Xt)
        log("stage S start f=%.4g" % obj(Xs[-1], Ls[-1])[0])
        dl = time.time() + STEER_FRAC * (T0 + TIME_BUDGET - time.time())
        L2, X2 = stage_b(m, Ls[-1], Xs[-1], dl, obj=obj, start_kicks=0, tag="S")
        Ls += L2[1:]; Xs += X2[1:]
        log("stage S end f=%.4g density %.5f" % (obj(Xs[-1], Ls[-1])[0], 6 * m ** 3 * 4 / 3 * np.pi * R ** 3 / np.linalg.det(Ls[-1])))
    reached = TARGET == "fcc" and obj(Xs[-1], Ls[-1])[0] < STEER_TOL
    if CONTINUE and not (reached and SKIP_B_IF_REACHED):
        L2, X2 = stage_b(m, Ls[-1], Xs[-1], T0 + TIME_BUDGET, start_kicks=(START_KICKS if TARGET is None else 0))
        Ls += L2[1:]
        Xs += X2[1:]
    Ls = np.array(Ls); Xs = np.array(Xs)
    if TRIM_TO_BEST:
        k = int(np.argmin(np.linalg.det(Ls)))
        Ls, Xs = Ls[:k + 1], Xs[:k + 1]
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    np.savez(OUT, m=m, Ls=Ls, Xs=Xs)
    log("wrote %s (%d states), final density %.5f" % (OUT, len(Ls), 6 * m ** 3 * 4 / 3 * np.pi * R ** 3 / np.linalg.det(Ls[-1])))


if __name__ == "__main__":
    main()
