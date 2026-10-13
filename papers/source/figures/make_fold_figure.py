"""Figure 2 of the paper (v2): the exact motion from the cab framework to fcc.

Everything is computed from the closed form of results/fold-continuous-fcc/PROOF.md:
motif P = {+-e_i}, ball radius 1/sqrt(2), lattice rows L(s) = (2 + a) I + t C with
t = 2s/(1+s^2), a = sqrt(2)(1-s^2)/(1+s^2). Nothing is read from sampled paths.
Run from the repository root:  python papers/source/figures/make_fold_figure.py
(needs numpy and matplotlib; see requirements-dev.txt). Writes fig_motion.pdf next to this file.
"""
import itertools
import os

import numpy as np
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

HERE = os.path.dirname(os.path.abspath(__file__))
SQ2 = np.sqrt(2.0)
C = np.array([[0, -1, -1], [1, 0, -1], [1, -1, 0]], float)
P = np.array([s * np.eye(3)[i] for i in range(3) for s in (1, -1)])
RADIUS = 1 / SQ2
BLUE, ORANGE, INK, MUTED = "#2a78d6", "#eb6834", "#222222", "#8a8a8a"
plt.rcParams.update({"font.size": 9, "font.family": "serif", "mathtext.fontset": "cm",
                     "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.4})


def lattice(s):
    t = 2 * s / (1 + s * s)
    a = SQ2 * (1 - s * s) / (1 + s * s)
    return (2 + a) * np.eye(3) + t * C


def density(s):
    return 2 * SQ2 * np.pi / np.linalg.det(lattice(s))


def view_basis():
    """Fixed orthographic view, the same for every snapshot (the octahedra do not rotate)."""
    d = np.array([1.0, 0.42, 0.25])
    d /= np.linalg.norm(d)
    up = np.array([0.0, 0.0, 1.0])
    x = np.cross(up, d); x /= np.linalg.norm(x)
    y = np.cross(d, x)
    return x, y, d


def draw_state(ax, s, cells=2):
    L = lattice(s)
    x, y, d = view_basis()
    shifts = [np.array(n) for n in itertools.product(range(cells), repeat=3)]
    centres = {}
    for n in shifts:
        for k, p in enumerate(P):
            centres[(tuple(n), k)] = p + n @ L
    # bars: 12 octahedron edges per cell, and the bridge from +e_i to -e_i of the cell n + e_i
    segs = []
    for n in shifts:
        for k1, k2 in itertools.combinations(range(6), 2):
            if abs(P[k1] @ P[k2]) < 1e-12:
                segs.append((centres[(tuple(n), k1)], centres[(tuple(n), k2)], INK, 0.55))
        for i in range(3):
            m = n.copy(); m[i] += 1
            if m[i] < cells:
                segs.append((centres[(tuple(n), 2 * i)], centres[(tuple(m), 2 * i + 1)], ORANGE, 1.1))
    items = []
    for a_, b_, col, lw in segs:
        items.append(((a_ + b_) @ d / 2, "seg", (a_, b_, col, lw)))
    for (_, k), c in centres.items():
        items.append((c @ d, "ball", c))
    items.sort(key=lambda z: z[0])
    for _, kind, obj in items:
        if kind == "seg":
            a_, b_, col, lw = obj
            ax.plot([a_ @ x, b_ @ x], [a_ @ y, b_ @ y], color=col, lw=lw, solid_capstyle="round", zorder=2)
        else:
            c = obj
            ax.add_patch(Circle((c @ x, c @ y), 0.42 * RADIUS, fc="#cde2fb", ec=BLUE, lw=0.6, zorder=2))
    ax.set_aspect("equal")
    ax.axis("off")


def main():
    fig = plt.figure(figsize=(7.0, 4.6))
    snaps = [(0.0, "(a) $s=0$: cab, $\\phi=0.2233$"), (0.5, "(b) $s=1/2$"), (1.0, "(c) $s=1$: fcc, $\\phi=\\pi/\\sqrt{18}$")]
    boxes = []
    for j, (s, title) in enumerate(snaps):
        ax = fig.add_axes([0.01 + j * 0.33, 0.40, 0.31, 0.58])
        draw_state(ax, s)
        ax.set_title(title, fontsize=9)
        boxes.append(ax)
    # common limits so the shrinking of the cell is visible
    xs = [a.get_xlim() for a in boxes]; ys = [a.get_ylim() for a in boxes]
    lo_x = min(v[0] for v in xs); hi_x = max(v[1] for v in xs)
    lo_y = min(v[0] for v in ys); hi_y = max(v[1] for v in ys)
    for a in boxes:
        a.set_xlim(lo_x, hi_x); a.set_ylim(lo_y, hi_y)
    ax = fig.add_axes([0.10, 0.08, 0.84, 0.25])
    ss = np.linspace(0, 1, 401)
    ax.plot(ss, [density(s) for s in ss], color=BLUE)
    for val, lab in [(np.pi / np.sqrt(18), "$\\pi/\\sqrt{18}$ (fcc)"), (27 * np.pi * SQ2 / 256, "0.4686 (cubic tilt locks)")]:
        ax.axhline(val, color=MUTED, lw=0.7, ls=":")
        ax.text(0.01, val + 0.012, lab, color=MUTED, fontsize=8)
    for s, _ in snaps:
        ax.plot([s], [density(s)], "o", color=ORANGE, ms=4)
    ax.set_xlabel("$s$")
    ax.set_ylabel("packing fraction")
    ax.set_xlim(0, 1)
    ax.set_ylim(0.18, 0.80)
    ax.set_title("(d) packing fraction $2\\sqrt{2}\\,\\pi/\\det L(s)$", fontsize=9)
    out = os.path.join(HERE, "fig_motion.pdf")
    fig.savefig(out)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
