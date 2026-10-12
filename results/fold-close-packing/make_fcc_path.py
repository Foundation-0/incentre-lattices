#!/usr/bin/env python3
"""Deterministically extend the hcp fold path to a pure fcc end state.

Input : autoresearch-fold/best/path.npz   (793 states, ends in hcp, density pi/sqrt18)
Output: a path whose first 793 states are the input and whose last 400 states are two rigid
        stacking slides (200 states each). Every slide moves the upper layers by one
        close-packed stacking vector and lifts them over the bridge site so that the chosen
        interface bar keeps its length 2R.

Usage : python make_fcc_path.py OUT.npz
No random numbers are used. Requires numpy; imports the bar graph helpers from
autoresearch-fold/verify.py (a generator may use them; the CHECKS are separate programs).
"""
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
FOLD = HERE.parents[1] / "autoresearch-fold"
import importlib.util  # noqa: E402

_spec = importlib.util.spec_from_file_location("fold_verify", FOLD / "verify.py")
verify = importlib.util.module_from_spec(_spec)  # the hash-locked checker, loaded by path
_spec.loader.exec_module(verify)


def build(src=FOLD / "best" / "path.npz", steps=200):
    d = np.load(src)
    Ls, Xs = list(d["Ls"]), list(d["Xs"])
    m = int(d["m"])
    R = verify.R
    D = 2 * R
    I, J, S = verify.bar_graph(m)
    X0, L0 = Xs[-1].copy(), Ls[-1].copy()
    a, b, c = L0
    n = np.cross(a, b)
    n /= np.linalg.norm(n)
    hc = c @ n
    z = X0 @ n
    lay = np.round((z - z.min()) / (hc / 4)).astype(int) % 4  # layer index 0..3
    e1 = a / np.linalg.norm(a)
    e2 = np.cross(n, e1)

    def u(deg):
        return (D / np.sqrt(3)) * (np.cos(np.radians(deg)) * e1 + np.sin(np.radians(deg)) * e2)

    top = (lay == 2) | (lay == 3)

    def interface_bar(la, lb):
        for i, j, s in zip(I, J, S):
            if {lay[i], lay[j]} == {la, lb} and i // 6 != j // 6:
                return i, j, s

    def lifted(X, L, lat, mode, ib):
        i, j, s = ib

        def apply(dz):
            Y = X.copy()
            M = L.copy()
            sh = lat + dz * n
            if mode == 1:
                Y[top] += sh
                M[2] += sh
            else:
                M[2] += sh
            return Y, M

        Y, M = apply(0.0)
        v = Y[j] + s @ M - Y[i]
        sgn = np.sign(v @ n)
        vl = v - (v @ n) * n
        if vl @ vl > D ** 2:
            raise ValueError("interface bar planar component exceeds D: sqrt of negative")
        dz = sgn * np.sqrt(D ** 2 - vl @ vl) - (v @ n)
        return apply(dz)

    path_L, path_X = [], []
    X, L = X0, L0
    for mode, vec, (la, lb) in ((1, u(0), (1, 2)), (2, u(120), (3, 0))):
        ib = interface_bar(la, lb)
        for t in range(1, steps + 1):
            Y, M = lifted(X, L, vec * t / steps, mode, ib)
            path_X.append(Y)
            path_L.append(M)
        X, L = path_X[-1], path_L[-1]
    return m, np.array(Ls + path_L), np.array(Xs + path_X)


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "path_fcc.npz"
    m, Ls, Xs = build()
    np.savez_compressed(out, m=m, Ls=Ls, Xs=Xs)
    print(f"wrote {out}: {len(Xs)} states")
