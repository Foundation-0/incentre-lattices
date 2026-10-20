"""Contact numbers and densities at the collision (end) states of the tilt mechanisms."""
import numpy as np
from scipy.spatial.transform import Rotation as Rot
from mechanism import *

def contacts(n, R, ncell=2):
    steps, R = build_state(n, R)
    pts, bonded, index, cells = positions(n, steps, R, ncell=ncell)
    cid = [index[((0,) * n, i, s)] for i in range(n) for s in (1, -1)]
    out = []
    for a in cid:
        d = np.linalg.norm(pts - pts[a], axis=1)
        out.append(int(np.sum(np.abs(d - 2 * RB) < 1e-9)))
    dens = len(cid) * (np.pi ** (n / 2) / __import__('math').gamma(n / 2 + 1)) * RB ** n / abs(np.linalg.det(steps))
    return out, dens, steps

print("2D theta=0:", contacts(2, rot2(0)))
print("2D theta=45:", contacts(2, rot2(np.pi / 4)), "pi/4 =", np.pi / 4)
for name, ax, th in [("[001]", [0, 0, 1], np.pi / 4), ("[111]", [1, 1, 1], np.pi / 3), ("[011]", [0, 1, 1], np.arccos(1 / np.sqrt(3)))]:
    ax = np.array(ax, float) / np.linalg.norm(ax)
    c, dens, steps = contacts(3, Rot.from_rotvec(th * ax).as_matrix())
    print("3D", name, "theta=%.4f deg" % np.degrees(th), "contacts per ball:", c, "density %.5f" % dens, "steps diag", np.round(np.diag(steps), 6))
print("3D theta=0:", contacts(3, np.eye(3))[:2])
# exact value checks
R2 = (np.sqrt(2) - 1) ** 2
print("linear scale factor 1 - R^2 = 2R:", 1 - R2, 2 * (np.sqrt(2) - 1))
