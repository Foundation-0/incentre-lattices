"""Hinged-squares fold read as a ball framework (the 2D part of the incentre-lattice fold remarks).

Squares of 4 dots at +-c e_i in cells of side 2, cross bars of length 2 - 2c, alternate rotation
+-theta.  Exact: L(theta) = 2R + 2c cos(theta); first new contact at theta = pi/4 for the
incentre c = 2 - sqrt2; there L = 4R, strain -R^2, density pi R^2 -> pi/4.
Engine functions only (numpy + sympy); verify.py is the certificate driver.
"""
import numpy as np
import sympy as sp

R_ = sp.sqrt(2) - 1
c_ = 2 - sp.sqrt(2)
th = sp.symbols("theta", real=True)
R = float(R_)
c = float(c_)


def exact_facts():
    """sympy: gap(theta), gap(pi/4), L(pi/4) - 4R, strain + R^2, densities."""
    L = 2 * R_ + 2 * c_ * sp.cos(th)
    gap = sp.simplify(L - 2 * c_ * sp.sin(th) - 2 * R_)      # gap of the pair that closes (+c e_1 dots of vertically adjacent squares)
    return {
        "L": L,
        "gap": gap,
        "gap_at_pi_over_4": sp.simplify(gap.subs(th, sp.pi / 4)),
        "L_minus_4R": sp.simplify(L.subs(th, sp.pi / 4) - 4 * R_),
        "strain_plus_R2": sp.simplify(L.subs(th, sp.pi / 4) / 2 - 1 + R_ ** 2),
        "density_start": sp.pi * R_ ** 2,
        "density_end": 4 * sp.pi * R_ ** 2 / (4 * R_) ** 2,
    }


def dots2(cc, t, Lp, rng=2):
    P = []
    for i in range(-rng, rng + 1):
        for j in range(-rng, rng + 1):
            s = 1 if (i + j) % 2 == 0 else -1
            for k in range(4):
                a = s * t + k * np.pi / 2
                P.append([Lp * i + cc * np.cos(a), Lp * j + cc * np.sin(a)])
    return np.array(P)


def lock2(cc):
    """Lock angle (degrees, 0.01 deg grid) and the nearest distances / (2r) of a central dot."""
    b = 2 - 2 * cc
    s = cc * 2 ** 0.5
    r2 = min(s, b)
    P0 = dots2(cc, 0, 2.0)
    D0 = np.linalg.norm(P0[:, None] - P0[None], axis=2)
    linked = (np.abs(D0 - s) < 1e-9) | (np.abs(D0 - b) < 1e-9)
    np.fill_diagonal(linked, True)
    for t in np.linspace(0, np.pi / 2, 9001)[1:]:
        P = dots2(cc, t, b + 2 * cc * np.cos(t))
        D = np.linalg.norm(P[:, None] - P[None], axis=2)
        if np.where(linked, np.inf, D).min() <= r2 + 1e-9:
            i0 = np.argmin(np.linalg.norm(P - P[:, :].mean(0), axis=1))
            nn = np.sort(D[i0])[1:6] / r2
            return float(np.degrees(t)), nn
    return None, None
