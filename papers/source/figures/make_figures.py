"""Generate all figures of the paper as PDFs in this folder.
Run: ../.venv/bin/python make_figures.py
"""
import itertools
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "proofs"))

BLUE, ORANGE, INK, MUTED = "#2a78d6", "#eb6834", "#222222", "#8a8a8a"
plt.rcParams.update({"font.size": 9, "font.family": "serif", "mathtext.fontset": "cm",
                     "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.5})
R = 2 ** 0.5 - 1
c = 2 - 2 ** 0.5


def save(fig, name):
    fig.savefig(os.path.join(HERE, name), bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


# ---------------------------------------------------------------- Figure 1
def fig_construction():
    fig = plt.figure(figsize=(7.0, 2.5))
    ax = fig.add_subplot(1, 3, 1)
    sq = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]])
    ax.add_patch(Polygon(sq, closed=True, fill=False, ec=INK, lw=1.0))
    ax.plot([-1, 1], [-1, 1], color=INK, lw=0.7); ax.plot([-1, 1], [1, -1], color=INK, lw=0.7)
    for k in range(4):
        a = k * np.pi / 2
        d = c * np.array([np.cos(a), np.sin(a)])
        ax.add_patch(Circle(d, R, fill=True, fc="#cde2fb", ec=BLUE, lw=0.8))
        ax.plot(*d, "o", color=BLUE, ms=3)
        for v in sq:  # lines from the incentre to the corners (the user's drawing)
            if np.dot(v, d) > 0:
                ax.plot([d[0], v[0]], [d[1], v[1]], color=ORANGE, lw=0.9)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.1, 1.1); ax.set_ylim(-1.1, 1.1)
    ax.set_title("(a)", fontsize=9)

    ax = fig.add_subplot(1, 3, 2, projection="3d")
    cube = np.array(list(itertools.product((-1, 1), repeat=3)))
    for a_, b_ in itertools.combinations(range(8), 2):
        if np.sum(np.abs(cube[a_] - cube[b_])) == 2:
            ax.plot(*zip(cube[a_], cube[b_]), color=INK, lw=0.6)
    dots = np.array([s * c * e for e in np.eye(3) for s in (1, -1)])
    for v in cube:  # pyramid edges, faint
        ax.plot(*zip((0, 0, 0), v), color=MUTED, lw=0.4, ls=":")
    for p_, q_ in itertools.combinations(range(6), 2):
        if abs(np.dot(dots[p_], dots[q_])) < 1e-9:
            ax.plot(*zip(dots[p_], dots[q_]), color=BLUE, lw=0.9)
    ax.scatter(*dots.T, s=60, color=BLUE, depthshade=False, edgecolor="white", linewidth=0.5)
    ax.set_box_aspect((1, 1, 1)); ax.axis("off"); ax.view_init(18, 28)
    ax.set_title("(b)", fontsize=9)

    ax = fig.add_subplot(1, 3, 3, projection="3d")
    for x0 in (0, 2):
        cc = cube + np.array([x0, 0, 0])
        for a_, b_ in itertools.combinations(range(8), 2):
            if np.sum(np.abs(cube[a_] - cube[b_])) == 2:
                ax.plot(*zip(cc[a_], cc[b_]), color=MUTED, lw=0.5)
        d = dots + np.array([x0, 0, 0])
        for p_, q_ in itertools.combinations(range(6), 2):
            if abs(np.dot(dots[p_], dots[q_])) < 1e-9:
                ax.plot(*zip(d[p_], d[q_]), color=BLUE, lw=0.9)
        ax.scatter(*d.T, s=40, color=BLUE, depthshade=False, edgecolor="white", linewidth=0.5)
    ax.plot([c, 2 - c], [0, 0], [0, 0], color=ORANGE, lw=2.2)
    ax.text(1.0, 0, -1.55, "cross bar", color=INK, fontsize=7, ha="center")
    ax.set_box_aspect((2, 1, 1)); ax.axis("off"); ax.view_init(18, -60)
    ax.set_title("(c)", fontsize=9)
    save(fig, "fig1_construction.pdf")


# ---------------------------------------------------------------- Figure 2
def dots2(t, Lp, rng):
    P, links = [], []
    for i in range(-rng, rng + 1):
        for j in range(-rng, rng + 1):
            s = 1 if (i + j) % 2 == 0 else -1
            for k in range(4):
                a = s * t + k * np.pi / 2
                P.append([Lp * i + c * np.cos(a), Lp * j + c * np.sin(a)])
    return np.array(P)


def fig_fold():
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 3.1))
    for ax, t, lab in [(axs[0], 0.0, r"(a) $\theta=0$: 4.8.8, $\phi=\pi(\sqrt{2}-1)^2\approx0.539$"),
                       (axs[1], np.pi / 4, r"(b) $\theta=\pi/4$: square packing, $\phi=\pi/4$")]:
        Lp = 2 * R + 2 * c * np.cos(t)
        P = dots2(t, Lp, 3)
        D = np.linalg.norm(P[:, None] - P[None], axis=2)
        for i, j in zip(*np.where(np.abs(D - 2 * R) < 1e-9)):
            if i < j:
                ax.plot(*zip(P[i], P[j]), color=INK, lw=0.5, zorder=1)
        for p in P:
            ax.add_patch(Circle(p, R, fc="#cde2fb", ec=BLUE, lw=0.7, zorder=2))
        ax.set_aspect("equal"); ax.axis("off")
        w = 2.05 * Lp
        ax.set_xlim(-w, w); ax.set_ylim(-w, w)
        ax.set_title(lab, fontsize=8.5)
        # one cell outline
        ax.add_patch(Polygon(np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]]) * Lp / 2, closed=True, fill=False,
                             ec=ORANGE, lw=0.9, ls="--", zorder=3))
    save(fig, "fig2_fold.pdf")


# ---------------------------------------------------------------- Figure 3
def fig_modes():
    from rigidity_nd import check_n, elastic
    ns = np.arange(3, 9)
    zero, ss, nu = [], [], []
    for n in ns:
        o = check_n(int(n))
        zero.append(o["generic"][0][0]); ss.append(o["generic"][0][1])
        nu.append(elastic(int(n), c)[4])
    nn = np.linspace(3, 8, 200)
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.5))
    ax = axs[0]
    ax.plot(nn, nn * (nn - 1) / 2, color=BLUE, lw=1.2)
    ax.plot(ns, zero, "o", color=BLUE, ms=6, mec="white", mew=0.6)
    ax.plot(nn, nn * (nn - 3) / 2, color=ORANGE, lw=1.2, ls="--")
    ax.plot(ns, ss, "s", color=ORANGE, ms=6, mec="white", mew=0.6)
    ax.text(7.6, 7.6 * 6.6 / 2 + 1.5, r"zero modes $n(n-1)/2$", ha="right", fontsize=8, color=INK)
    ax.text(6.1, 3.0, r"self-stresses $n(n-3)/2$", ha="left", fontsize=8, color=INK)
    ax.set_xlabel(r"dimension $n$"); ax.set_ylabel("per cell, generic $q$")
    ax.set_title("(a)", fontsize=9); ax.grid(alpha=0.25, lw=0.4)
    ax = axs[1]
    ax.plot(nn, 1 / (nn ** 2 - nn - 1), color=BLUE, lw=1.2)
    ax.plot(ns, nu, "o", color=BLUE, ms=6, mec="white", mew=0.6)
    ax.set_yscale("log")
    ax.set_yticks([1 / 5, 1 / 11, 1 / 19, 1 / 29, 1 / 41, 1 / 55])
    ax.set_yticklabels(["1/5", "1/11", "1/19", "1/29", "1/41", "1/55"])
    ax.minorticks_off()
    ax.set_xlabel(r"dimension $n$"); ax.set_ylabel(r"$\nu_{100}$")
    ax.text(5.2, 0.09, r"$\nu_{100}=1/(n^2-n-1)$", fontsize=8.5, color=INK)
    ax.set_title("(b)", fontsize=9); ax.grid(alpha=0.25, lw=0.4, which="both")
    fig.tight_layout()
    save(fig, "fig3_modes.pdf")


# ---------------------------------------------------------------- Figure 4
def fig_cab6():
    g0 = 2 ** 0.5 - 1
    fcc = np.pi / (3 * 2 ** 0.5)
    phi = lambda g: np.where(g <= g0, np.pi / 6 * (1 + 6 * g ** 3), np.pi / 6 * g0 ** 3 * (6 + g ** -3))
    lo, hi = ((g0) / 6) ** (1 / 3), (4 + 7 * 2 ** 0.5) ** (-1 / 3)
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.6))
    for ax, (a, b), title in [(axs[0], (0.05, 1.0), "(a)"), (axs[1], (0.405, 0.42), "(b)")]:
        g = np.linspace(a, b, 2000)
        ax.plot(g, phi(g), color=BLUE)
        ax.axhline(fcc, color=MUTED, lw=0.8, ls="--")
        ax.plot([g0], [phi(np.array(g0))], "o", color=ORANGE, ms=6, mec="white", mew=0.6, zorder=3)
        ax.set_xlim(a, b)
        ax.set_xlabel(r"$\gamma=r_{\mathrm{B}}/r_{\mathrm{Ca}}$"); ax.set_title(title, fontsize=9)
        ax.grid(alpha=0.25, lw=0.4)
    axs[0].set_ylabel(r"$\phi_{\max}(\gamma)$")
    axs[0].text(0.88, fcc + 0.012, "fcc", color=INK, fontsize=8)
    axs[0].annotate(r"$\gamma=\sqrt{2}-1$" "\n" r"$\phi=\frac{\pi}{6}(30\sqrt{2}-41)$", xy=(g0, phi(np.array(g0))),
                    xytext=(0.55, 0.52), fontsize=8, arrowprops=dict(arrowstyle="-", lw=0.6, color=INK))
    axs[1].axvspan(lo, hi, color="#cde2fb", lw=0)
    axs[1].text((lo + hi) / 2, 0.733, r"$\phi>\phi_{\rm fcc}$", ha="center", fontsize=8)
    axs[1].set_ylim(0.725, 0.75)
    axs[1].set_xticks([0.405, 0.41, 0.415, 0.42])
    axs[1].text(0.4185, fcc + 0.0007, "fcc", color=INK, fontsize=8)
    fig.tight_layout()
    save(fig, "fig4_cab6.pdf")


if __name__ == "__main__":
    fig_construction()
    fig_fold()
    fig_modes()
    fig_cab6()
