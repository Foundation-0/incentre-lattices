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

DESC = "free m=1: robust pipeline (1-thread, cycle breaker, 3 start kicks SEED=1), no fold, density descent"

# ----------------------------------------------------------------- knobs
M = 1                    # supercell (cubic checkerboard fold needs M even)
FOLD_END = 0.0
CONTINUE = True          # stage B on/off
SEED = 1
WALL_LIMIT = 3000.0      # safety limit only (s). Stage B is bounded by MAX_STEPS and STALL_KICKS, not by time.
                         # If the limit is hit the run is INCOMPLETE: nothing is written and the exit code is 2,
                         # so a slow machine can never silently produce a different end state.
HX = 0.006               # max dot displacement per step (verify limit 0.01)
HL = 0.006               # max lattice Frobenius change per step (verify limit 0.01)
MAX_STEPS = 4000         # stage B step bound (the search ends here or at a stall, whichever is first)
START_KICKS = 3          # seeded random cone kicks at the start of stage B (makes SEED matter)
NOPROG_STEPS = 40        # descent steps without a new best -> treat as stationary (cycle breaker)
STALL_KICKS = 25         # stop after this many consecutive kicks without new best density
CONTACT_MARGIN = 0.03    # relative distance within which non-bar pairs become LP constraints
TRIM_TO_BEST = True      # cut path at its max-density state
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


def stage_b(m, L, X, deadline):
    """Returns (Ls, Xs, complete); complete is False only if the wall-clock safety limit was hit."""
    eng = Engine(m)
    N = eng.N
    Ls, Xs = [L.copy()], [X.copy()]
    capf = 1.0
    best = np.linalg.det(L)
    best_k = 0
    kicks = 0
    logt = time.time()
    for k in range(MAX_STEPS):
        if time.time() > deadline:
            log("stage B: INCOMPLETE (wall-clock safety limit %.0f s hit at step %d)" % (WALL_LIMIT, k))
            return Ls, Xs, False
        u, rate = eng.direction(X, L)
        mode = "descent"
        if k < START_KICKS:
            u, _ = eng.direction(X, L, "kick", eng.rng.normal(size=eng.n))
            mode = "kick"
            if u is None:
                continue
        elif u is None or rate > -1e-7 * 0.01 or k - best_k > NOPROG_STEPS:
            # stationary at first order: random kick inside the feasible cone
            kicks += 1
            if kicks > STALL_KICKS:
                log("stage B: stalled (%d kicks without progress) -> jammed or stuck" % STALL_KICKS); break
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
        if det < best - 1e-9:
            best = det; best_k = k; kicks = 0
        if time.time() - logt > 10:
            logt = time.time()
            dens = N * 4 / 3 * np.pi * R ** 3 / det
            log("k=%d density=%.5f best=%.5f mode=%s h=%.2e" % (k, dens, N * 4 / 3 * np.pi * R ** 3 / best, mode, h))
    return Ls, Xs, True


# ----------------------------------------------------------------- main
def main():
    m = M
    if FOLD_END > 0:
        Ls, Xs = stage_a(m, FOLD_END)
    else:
        L0, X0 = canonical_state(m)
        Ls, Xs = [L0], [X0]
    log("stage A done: %d states, density %.5f" % (len(Ls), 6 * m ** 3 * 4 / 3 * np.pi * R ** 3 / np.linalg.det(Ls[-1])))
    if CONTINUE:
        L2, X2, complete = stage_b(m, Ls[-1], Xs[-1], T0 + WALL_LIMIT)
        if not complete:
            log("search INCOMPLETE: no path written")
            sys.exit(2)
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
