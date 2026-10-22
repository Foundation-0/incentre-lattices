"""Symbolic checks: (1) elastic constants of the unit-spring n-D contact network (symmetric-ansatz reduction);
(2) 2D fold end state and gap formula; (3) 3D [111] tilt end spacing."""
import sympy as sp

# (1) energy per cell E = 2 sum_i (e_i - s_i)^2 + sum_{i<j} (s_i + s_j)^2   (unit springs, cell volume 2^n)
for n in range(2, 8):
    e = sp.symbols(f'e0:{n}'); s = sp.symbols(f's0:{n}')
    E = 2 * sum((e[i] - s[i]) ** 2 for i in range(n)) + sum((s[i] + s[j]) ** 2 for i in range(n) for j in range(i + 1, n))
    sol = sp.solve([sp.diff(E, v) for v in s], s, dict=True)[0]
    Emin = sp.expand(E.subs(sol))
    C11 = sp.diff(Emin, e[0], 2) / 2 ** n; C12 = sp.diff(Emin, e[0], e[1]) / 2 ** n
    nu = sp.simplify(C12 / (C11 + (n - 2) * C12))
    print(n, "C11 =", sp.nsimplify(C11), " formula", sp.Rational(4 * (n - 1) ** 2, n ** 2 * 2 ** n),
          " C12 =", C12, " formula", sp.Rational(4, n ** 2 * 2 ** n), " nu =", nu, " 1/(n^2-n-1) =", sp.Rational(1, n * n - n - 1))

# (2) 2D: spacing a(theta) = 2 - 2t(1-cos theta); at 45 deg with t = 2 - sqrt2
t = 2 - sp.sqrt(2); R = sp.sqrt(2) - 1
a45 = 2 - 2 * t * (1 - sp.cos(sp.pi / 4))
print("2D a(45) =", sp.nsimplify(sp.simplify(a45)), " 2R =", sp.simplify(2 * R), " strain =", sp.simplify(a45 / 2 - 1), " -R^2 =", sp.expand(-R ** 2))
dens45 = sp.simplify(4 * sp.pi * R ** 2 / a45 ** 2)
print("   area ratio =", sp.simplify((a45 / 2) ** 2), " = 4R^2 =", sp.expand(4 * R ** 2), "; density 4*pi*R^2/a^2 =", dens45)
assert sp.simplify(dens45 - sp.pi / 4) == 0
# gap between the balls of neighbouring clusters along the bar axis: a(theta) - 2t sin(theta) - 2R = 2 sqrt2 t cos(theta + pi/4)
th = sp.symbols('theta', real=True)
a = 2 - 2 * t * (1 - sp.cos(th)); gap = a - 2 * t * sp.sin(th) - 2 * R
gap_ok = sp.simplify(sp.expand_trig(gap - 2 * sp.sqrt(2) * t * sp.cos(th + sp.pi / 4))) == 0
print("   gap(theta) = 2 sqrt2 t cos(theta + pi/4):", gap_ok, "; gap(45) =", sp.simplify(gap.subs(th, sp.pi / 4)))
assert gap_ok and sp.simplify(gap.subs(th, sp.pi / 4)) == 0
# (3) 3D [111] tilt by 60deg: alpha = 4t/3 + sqrt(L^2 - 2t^2/9)
L = 2 - 2 * t
alpha = sp.Rational(4, 3) * t + sp.sqrt(L ** 2 - 2 * t ** 2 / 9)
print("3D [111] end spacing alpha =", sp.simplify(alpha), "; 8t/3 =", sp.simplify(8 * t / 3), "; equal:", sp.simplify(alpha - 8 * t / 3) == 0)
print("   volume ratio (alpha/2)^3 =", sp.N((alpha / 2) ** 3), "; end density =", sp.N(6 * sp.Rational(4, 3) * sp.pi * R ** 3 / alpha ** 3))
