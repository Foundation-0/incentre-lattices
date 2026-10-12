"""Symbolic (sympy) check of the individual steps in the proof of Theorem 1, n = 2..8.

Step 1 (character decomposition): for an intra bar (i,s)-(j,t) the first-order extension under
  affine strain eps and periodic internal displacement u_{k,s} = a_k + s b_k is
  delta = (1/sqrt2)[ c(eps_ii+eps_jj) + C + sA - tB - st(2c eps_ij + D) ],
  C = b_i^i + b_j^j, A = a_i^i - a_j^i, B = a_i^j - a_j^j, D = b_j^i + b_i^j.
  The cross bar along e_i has delta = 2R eps_ii - 2 b_i^i with 2R = 2 - 2c.
Step 2: sum over (s,t) of delta^2 = 2[(c(e_i+e_j)+C)^2 + A^2 + B^2 + (2c eps_ij + D)^2].
Step 3: min over x of  sum_{i<j}(c(e_i+e_j)+x_i+x_j)^2 + 2 sum_i ((1-c) e_i - x_i)^2
        = 2(n-1)^2/n^2 sum e_i^2 + (2/n^2) sum_{i!=j} e_i e_j   (independent of c).
Run: python proof_steps.py   (or import run() from verify.py)
"""
import itertools
import sympy as sp

def run():
    ok_all = True
    for n in range(2, 9):
        c = sp.symbols("c", positive=True)
        eps = sp.Matrix(n, n, lambda i, j: sp.Symbol(f"e{min(i,j)}{max(i,j)}"))
        a = [sp.Matrix(n, 1, lambda k, _: sp.Symbol(f"a{i}_{k}")) for i in range(n)]
        b = [sp.Matrix(n, 1, lambda k, _: sp.Symbol(f"b{i}_{k}")) for i in range(n)]
        E = sp.eye(n)
        u = lambda i, s: a[i] + s * b[i]
        p = lambda i, s: s * c * E[:, i]
        ok = True
        # step 1 and 2 on all pairs
        for i, j in itertools.combinations(range(n), 2):
            Cc = b[i][i] + b[j][j]; A = a[i][i] - a[j][i]; B = a[i][j] - a[j][j]; D = b[j][i] + b[i][j]
            tot = 0
            for s, t in itertools.product((1, -1), repeat=2):
                l = p(i, s) - p(j, t)
                dhat = l / (c * sp.sqrt(2))
                delta = (dhat.T * (eps * l + u(i, s) - u(j, t)))[0]
                claim = (c * (eps[i, i] + eps[j, j]) + Cc + s * A - t * B - s * t * (2 * c * eps[i, j] + D)) / sp.sqrt(2)
                ok &= sp.simplify(sp.expand(delta - claim)) == 0
                tot += delta ** 2
            claim2 = 2 * ((c * (eps[i, i] + eps[j, j]) + Cc) ** 2 + A ** 2 + B ** 2 + (2 * c * eps[i, j] + D) ** 2)
            ok &= sp.expand(tot - claim2) == 0
        # cross bar: from p(i,+) in cell 0 to p(i,-) in cell 2e_i
        for i in range(n):
            l = 2 * E[:, i] + p(i, -1) - p(i, 1)
            delta = (E[:, i].T * (eps * l + u(i, -1) - u(i, 1)))[0]
            ok &= sp.expand(delta - ((2 - 2 * c) * eps[i, i] - 2 * b[i][i])) == 0
        # step 3
        e = sp.symbols(f"x0:{n}"); x = sp.symbols(f"y0:{n}")
        F = sum((c * (e[i] + e[j]) + x[i] + x[j]) ** 2 for i, j in itertools.combinations(range(n), 2)) \
            + 2 * sum(((1 - c) * e[i] - x[i]) ** 2 for i in range(n))
        sol = sp.solve([sp.diff(F, xi) for xi in x], x, dict=True)[0]
        Fmin = sp.expand(sp.simplify(F.subs(sol)))
        claim3 = sp.Rational(2 * (n - 1) ** 2, n ** 2) * sum(ei ** 2 for ei in e) + sp.Rational(2, n ** 2) * sum(
            e[i] * e[j] for i in range(n) for j in range(n) if i != j)
        ok &= sp.simplify(Fmin - sp.expand(claim3)) == 0
        # Hessian positive definite (minimum, not saddle): eigenvalues n (mult n-1) and 2n
        H = sp.hessian(F, x) / 2
        ok &= sorted(H.eigenvals().items()) == sorted({n: n - 1, 2 * n: 1}.items()) if n > 1 else True
        # relaxed internal coordinate x_k = b_k^k = (2e_k - sum e/n)/n - c e_k
        ok &= all(sp.simplify(sol[x[k]] - ((2 * e[k] - sum(e) / n) / n - c * e[k])) == 0 for k in range(n))
        print(f"n={n}: steps 1-3 verified symbolically: {ok}")
        ok_all &= ok
    print("ALL STEPS VERIFIED" if ok_all else "FAILURE")
    return ok_all


if __name__ == "__main__":
    raise SystemExit(0 if run() else 1)
