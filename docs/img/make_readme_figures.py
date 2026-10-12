"""Regenerate the README figures in docs/img/.

Run from anywhere:  python docs/img/make_readme_figures.py
Needs numpy and matplotlib (see requirements-dev.txt).

Every plotted number is read from the repository's verified data:
  results/*/certificate.json, autoresearch-fold/best/path.npz,
  results/fold-close-packing/path_fcc.npz, and the closed form in
  results/cab6-packing/cab6_packing.py (imported, not copied).
Geometry that is drawn rather than measured (the construction figure) uses the
exact incentre formula r = (sqrt2 - 1)/2 for the unit square.
Each figure is written twice, <name>-light.png and <name>-dark.png, with a
transparent background, for use in a <picture> element.
"""
import contextlib
import importlib.util
import io
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Polygon

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
R_BALL = np.sqrt(2.0) - 1.0
KEPLER = np.pi / np.sqrt(18.0)
DPI = 150

THEMES = {
    "light": dict(fg="#1f2a2e", muted="#5f6f75", teal="#0f766e", amber="#b45309",
                  teal2="#5eaaa0", fill="#0f766e", grid="#1f2a2e"),
    "dark": dict(fg="#e6edf3", muted="#9aa7b2", teal="#2dd4bf", amber="#fbbf24",
                 teal2="#2a8f84", fill="#2dd4bf", grid="#e6edf3"),
}


def load(rel):
    return json.loads((ROOT / rel).read_text())


def style(t):
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "text.color": t["fg"], "axes.labelcolor": t["fg"],
        "axes.edgecolor": t["muted"], "xtick.color": t["muted"], "ytick.color": t["muted"],
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.facecolor": "none", "figure.facecolor": "none", "savefig.facecolor": "none",
        "legend.frameon": False, "axes.titleweight": "bold", "axes.titlesize": 11,
        "axes.titlelocation": "left",
    })


def save(fig, name, theme):
    fig.savefig(HERE / f"{name}-{theme}.png", dpi=DPI, transparent=True, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


def density(Ls, n=48):
    return n * 4.0 / 3.0 * np.pi * R_BALL ** 3 / np.abs(np.linalg.det(Ls))


# ---------------------------------------------------------------- construction
def fig_construction(t):
    r = (np.sqrt(2) - 1) / 2
    c = np.array([0.5, 0.5])
    corners = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float)
    # incentre of each triangle (side, centre): on the axis from the side midpoint to the centre
    inc = [np.array([0.5, r]), np.array([1 - r, 0.5]), np.array([0.5, 1 - r]), np.array([r, 0.5])]
    tris = [(corners[0], corners[1]), (corners[1], corners[2]), (corners[2], corners[3]), (corners[3], corners[0])]

    fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 4.9), gridspec_kw=dict(width_ratios=[1, 1.0]))
    # left: one square
    a.add_patch(Polygon(corners, closed=True, fill=False, ec=t["fg"], lw=2))
    for p, q in ((0, 2), (1, 3)):
        a.plot(*zip(corners[p], corners[q]), color=t["muted"], lw=1.3)
    for (p, q), i in zip(tris, inc):
        a.add_patch(Circle(i, r, fill=True, fc=t["fill"], alpha=0.14, ec="none"))
        a.add_patch(Circle(i, r, fill=False, ec=t["teal"], lw=1.4))
        for v in (p, q, c):
            a.plot(*zip(i, v), color=t["teal"], lw=1.3)
        a.plot(*i, "o", color=t["amber"], ms=7, zorder=5)
    a.plot(*c, "o", color=t["fg"], ms=3.5, zorder=5)
    a.set_title("One square, four incentres")
    a.text(0.5, -0.1, "diagonals, incentres (amber), inscribed circles, bisector segments", ha="center",
           va="top", fontsize=8.5, color=t["muted"])
    # right: tiled, truncated square tiling 4.8.8
    n = 3
    for i in range(n):
        for j in range(n):
            o = np.array([i, j], float)
            b.add_patch(Polygon(corners + o, closed=True, fill=False, ec=t["muted"], lw=0.9))
            for p, q in ((0, 2), (1, 3)):
                b.plot(*zip(corners[p] + o, corners[q] + o), color=t["muted"], lw=0.6, alpha=0.7)
    dots = {}
    for i in range(n):
        for j in range(n):
            for k, d in enumerate(inc):
                dots[(i, j, k)] = d + [i, j]
    bonds = []
    for i in range(n):
        for j in range(n):
            for k in range(4):
                bonds.append((dots[(i, j, k)], dots[(i, j, (k + 1) % 4)]))   # ring around the centre
            if i + 1 < n:
                bonds.append((dots[(i, j, 1)], dots[(i + 1, j, 3)]))        # across a vertical side
            if j + 1 < n:
                bonds.append((dots[(i, j, 2)], dots[(i, j + 1, 0)]))        # across a horizontal side
    # octagons around the interior lattice vertices
    for vx in range(1, n):
        for vy in range(1, n):
            ring = [dots[(vx - 1, vy - 1, 1)], dots[(vx - 1, vy - 1, 2)], dots[(vx - 1, vy, 1)], dots[(vx - 1, vy, 0)],
                    dots[(vx, vy, 3)], dots[(vx, vy, 0)], dots[(vx, vy - 1, 3)], dots[(vx, vy - 1, 2)]]
            ctr = np.mean(ring, axis=0)
            ring.sort(key=lambda p: np.arctan2(p[1] - ctr[1], p[0] - ctr[0]))
            b.add_patch(Polygon(ring, closed=True, fc=t["amber"], alpha=0.2, ec="none"))
    for p, q in bonds:
        b.plot(*zip(p, q), color=t["teal"], lw=2.2, solid_capstyle="round")
    for d in dots.values():
        b.plot(*d, "o", color=t["amber"], ms=4.5, zorder=5)
    b.set_title("Tiled: squares and regular octagons (4.8.8)")
    b.text(n / 2, -0.1, "dots joined across shared sides and around each centre", ha="center", va="top",
           fontsize=8.5, color=t["muted"])
    for ax, lim in ((a, (-0.04, 1.04)), (b, (-0.04, n + 0.04))):
        ax.set_aspect("equal")
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.axis("off")
    fig.subplots_adjust(wspace=0.08)
    return fig


# ---------------------------------------------------------------- fold density
def fig_fold(t):
    hcp = np.load(ROOT / "autoresearch-fold/best/path.npz")
    fcc = np.load(ROOT / "results/fold-close-packing/path_fcc.npz")
    d_hcp = density(hcp["Ls"])
    d_fcc = density(fcc["Ls"])
    n_h = len(d_hcp)
    assert abs(d_hcp[0] - 0.2232661157) < 1e-9 and abs(d_hcp[-1] - KEPLER) < 1e-9
    lock = 27 * np.pi * np.sqrt(2) / 256   # fold-3d-cubic certificate: 0.4685853099
    cert = load("results/fold-3d-cubic/certificate.json")
    assert abs(cert["exact_at_60_degrees"]["density_value"] - lock) < 1e-9
    fig, ax = plt.subplots(figsize=(7.2, 3.9))
    ax.plot(np.arange(n_h), d_hcp, color=t["teal"], lw=2.4, label=f"fold to hcp ({n_h} states)")
    ax.plot(np.arange(n_h - 1, len(d_fcc)), d_fcc[n_h - 1:], color=t["amber"], lw=1.8, ls=(0, (4, 2)),
            label=f"two layer slides to fcc ({len(d_fcc)} states in all)")
    for y, lab, dy in ((d_hcp[0], f"start  {d_hcp[0]:.4f}", 0.02), (lock, f"cubic-path lock  {lock:.4f}\n(reference; not on this path)", 0.02),
                       (KEPLER, "π/√18  = 0.7405", 0.02)):
        ax.axhline(y, color=t["muted"], lw=0.9, ls=":")
        ax.text(len(d_fcc) * 0.02, y + dy, lab, color=t["muted"], fontsize=9, va="bottom")
    ax.set_xlim(0, len(d_fcc))
    ax.set_ylim(0.15, 0.82)
    ax.set_xlabel("state index along the verified path")
    ax.set_ylabel("packing density")
    ax.set_title("The fold: 0.2233 to close packing")
    ax.legend(loc="lower right", fontsize=9, bbox_to_anchor=(1.0, 0.1))
    return fig


# ---------------------------------------------------------------- CaB6
def load_cab6():
    spec = importlib.util.spec_from_file_location("cab6_packing", ROOT / "results/cab6-packing/cab6_packing.py")
    mod = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(mod)
    return mod


def fig_cab6(t):
    mod = load_cab6()
    cert = load("results/cab6-packing/certificate.json")["exact"]
    g0 = np.sqrt(2) - 1
    pmax = cert["phi_max_value"]
    assert abs(mod.phi_closed(g0) - pmax) < 1e-9
    fcc = cert["fcc_density"]
    lo = cert["fcc_window_exact"]["lower_value"]
    hi = cert["fcc_window_exact"]["upper_value"]
    g = np.linspace(0.02, 1.2, 800)
    phi = np.array([mod.phi_closed(x) for x in g])
    fig, ax = plt.subplots(figsize=(7.2, 3.9))
    ax.plot(g, phi, color=t["teal"], lw=2.4, label="CaB$_6$ maximum packing fraction")
    ax.axhline(fcc, color=t["muted"], lw=1.1, ls="--")
    ax.text(1.19, fcc + 0.008, f"fcc  {fcc:.4f}", color=t["muted"], ha="right", fontsize=9)
    ax.plot(g0, pmax, "o", color=t["amber"], ms=8, zorder=5)
    ax.annotate("(π/6)(30√2 − 41) = %.4f\nat γ = √2 − 1 = %.4f" % (pmax, g0),
                (g0, pmax), (0.47, 0.805), color=t["amber"], fontsize=9,
                arrowprops=dict(arrowstyle="-", color=t["amber"], lw=1))
    ax.set_xlim(0, 1.2)
    ax.set_ylim(0.2, 0.85)
    ax.set_xlabel("size ratio  γ = r$_B$ / r$_{Ca}$")
    ax.set_ylabel("packing fraction")
    ax.set_title("CaB$_6$ hard-sphere packing")
    # inset on the fcc window
    ins = ax.inset_axes([0.6, 0.42, 0.38, 0.34])
    gz = np.linspace(0.405, 0.421, 400)
    ins.plot(gz, [mod.phi_closed(x) for x in gz], color=t["teal"], lw=2)
    ins.axhline(fcc, color=t["muted"], lw=1, ls="--")
    ins.axvspan(lo, hi, color=t["amber"], alpha=0.25, lw=0)
    ins.plot(g0, pmax, "o", color=t["amber"], ms=5)
    ins.set_xlim(0.405, 0.421)
    ins.set_ylim(0.735, 0.75)
    ins.set_xticks([lo, hi])
    ins.set_xticklabels([f"{lo:.4f}", f"{hi:.4f}"], fontsize=7)
    ins.set_yticks([])
    ins.set_title("beats fcc only here", fontsize=8, color=t["amber"], loc="left")
    for s in ins.spines.values():
        s.set_color(t["muted"])
    ins.spines["top"].set_visible(True)
    ins.spines["right"].set_visible(True)
    return fig


# ---------------------------------------------------------------- rigidity
def fig_rigidity(t):
    cert = load("results/rigidity-nd/certificate.json")
    rows = cert["counts_by_n"]
    el = {e["n"]: e for e in cert["elastic_by_n"]}
    ns = np.array([r["n"] for r in rows])
    zero = np.array([r["generic_q_zero_modes_selfstresses"][0] for r in rows])
    stress = np.array([r["generic_q_zero_modes_selfstresses"][1] for r in rows])
    nu = np.array([el[n]["nu100_value"] for n in ns])
    nn = np.linspace(2, 8, 100)
    panels = [
        (zero, nn * (nn - 1) / 2, "zero modes per cell", "n(n−1)/2"),
        (stress, nn * (nn - 3) / 2, "self-stresses per cell", "n(n−3)/2"),
        (nu, 1 / (nn ** 2 - nn - 1), "Poisson ratio ν$_{100}$", "1/(n$^2$−n−1)"),
    ]
    fig, axs = plt.subplots(1, 3, figsize=(9.6, 3.1))
    for ax, (y, f, title, form) in zip(axs, panels):
        ax.plot(nn, f, color=t["teal"], lw=2, zorder=2)
        ax.plot(ns[1:], y[1:], "o", color=t["amber"], ms=6, zorder=3)
        ax.plot(ns[:1], y[:1], "o", mfc="none", mec=t["amber"], ms=6, zorder=3)
        ax.set_title(title, fontsize=10)
        right = title.startswith("Poisson")
        ax.text(0.96 if right else 0.04, 0.88, form, transform=ax.transAxes, color=t["teal"], fontsize=10,
                ha="right" if right else "left")
        ax.set_xlabel("dimension n")
        ax.set_xticks(ns)
    axs[2].set_yscale("log")
    axs[0].text(0.04, 0.78, "markers: checked, n = 2..8\nopen: n = 2, outside the theorem",
                transform=axs[0].transAxes, color=t["muted"], fontsize=7.5, va="top")
    fig.tight_layout(w_pad=1.6)
    return fig


# ---------------------------------------------------------------- radius ratio
def fig_radius(t):
    cert = load("results/radius-ratio/certificate.json")["cases_satisfying_theorem"]
    want = ["tetrahedron -> Q tetrahedron (CN4)", "cube -> Q octahedron (CN6)", "octahedron -> Q cube (CN8)",
            "dodecahedron -> Q icosahedron (CN12)", "rhombic dodecahedron -> Q cuboctahedron"]
    labels = ["tetrahedron\nCN4", "cube\nCN6", "octahedron\nCN8", "dodecahedron\nCN12 (icosahedral)",
              "rhombic dodecahedron\nCN12 (cuboctahedral)"]
    by = {c["name"]: c for c in cert}
    hole = np.array([by[w]["hole_over_ball"] for w in want])
    # classical closed forms, computed here independently of the certificate
    cl = np.array([np.sqrt(1.5) - 1, np.sqrt(2) - 1, np.sqrt(3) - 1, 2 * np.sin(2 * np.pi / 5) - 1, 1.0])
    err = np.abs(hole - cl).max()
    assert err < 1e-7
    y = np.arange(len(want))[::-1]
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    ax.barh(y, hole, color=t["teal"], alpha=0.85, height=0.55, label="hole/ball from the construction")
    ax.plot(cl, y, "|", color=t["amber"], ms=26, mew=3, label="classical closed form")
    for yy, h in zip(y, hole):
        ax.text(h + 0.02, yy, f"{h:.4f}", va="center", fontsize=9, color=t["fg"])
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xlim(0, 1.2)
    ax.set_xlabel(f"radius ratio (hole / ball); largest difference {err:.0e}")
    ax.set_title("Radius ratio: construction versus the classical limits")
    ax.legend(loc="upper right", fontsize=8.5)
    return fig


# ---------------------------------------------------------------- fold frames
def fig_frames(t):
    hcp = np.load(ROOT / "autoresearch-fold/best/path.npz")
    Ls, Xs = hcp["Ls"], hcp["Xs"]
    d = density(Ls)
    lock = 27 * np.pi * np.sqrt(2) / 256
    picks = [(0, "start"), (int(np.argmax(d >= lock)), "passing density 0.469 (free continuation)"),
             (len(d) - 1, "end: hcp")]
    rot = np.array([[0.866, 0.0, -0.5], [0.25, 0.866, 0.433], [0.433, -0.5, 0.75]])  # fixed view
    fig, axs = plt.subplots(1, 3, figsize=(9.6, 3.4))
    c0 = np.array([[a, b, c] for a in (0, 1) for b in (0, 1) for c in (0, 1)], float) @ Ls[0] @ rot.T
    mid0 = (c0.min(axis=0) + c0.max(axis=0)) / 2
    half = (c0.max(axis=0) - c0.min(axis=0))[:2].max() / 2 + R_BALL
    for ax, (i, lab) in zip(axs, picks):
        L = Ls[i]
        frac = Xs[i] @ np.linalg.inv(L)
        f = frac - np.floor(frac)
        pos = f @ L
        # periodic images in the 3 positive directions to fill the box edges
        pts = [pos]
        for sx in (0, 1):
            for sy in (0, 1):
                for sz in (0, 1):
                    if (sx, sy, sz) != (0, 0, 0):
                        sel = np.ones(len(f), bool)
                        for k, s in enumerate((sx, sy, sz)):
                            if s:
                                sel &= f[:, k] < 1e-6
                        if sel.any():
                            pts.append((f[sel] + [sx, sy, sz]) @ L)
        P = np.vstack(pts)
        Q = P @ rot.T
        order = np.argsort(Q[:, 2])
        zmin, zmax = Q[:, 2].min(), Q[:, 2].max()
        for k in order:
            sh = (Q[k, 2] - zmin) / (zmax - zmin + 1e-9)
            ax.add_patch(Circle(Q[k, :2], R_BALL, fc=t["fill"], alpha=0.25 + 0.5 * sh, ec=t["teal"], lw=0.7))
        corners = np.array([[a, b, c] for a in (0, 1) for b in (0, 1) for c in (0, 1)], float) @ L @ rot.T
        for a_ in range(8):
            for b_ in range(a_ + 1, 8):
                if bin(a_ ^ b_).count("1") == 1:
                    ax.plot(*zip(corners[a_, :2], corners[b_, :2]), color=t["muted"], lw=0.7, alpha=0.7)
        ax.set_aspect("equal")
        ax.autoscale_view()
        ax.set_xlim(mid0[0] - half, mid0[0] + half)
        ax.set_ylim(mid0[1] - half, mid0[1] + half)
        ax.axis("off")
        ax.set_title(f"{lab}\ndensity {d[i]:.4f}", fontsize=9, loc="center")
    fig.tight_layout()
    return fig


FIGS = {
    "construction": fig_construction,
    "fold-density": fig_fold,
    "cab6-packing": fig_cab6,
    "rigidity": fig_rigidity,
    "radius-ratio": fig_radius,
    "fold-frames": fig_frames,
}

if __name__ == "__main__":
    for theme, t in THEMES.items():
        for name, fn in FIGS.items():
            style(t)
            save(fn(t), name, theme)
            print("wrote", f"{name}-{theme}.png")
