# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["numpy", "scipy", "trimesh", "manifold3d"]
# ///
"""Parametric vessel STLs for HemoCell: bifurcation, stenosis, S-bend.

Usage (from the repository root; uv fetches the dependencies):
    uv run scripts/geometry/make_vessels.py bifurcation -o cases/vessel_geometries/bifurcation.stl
    uv run scripts/geometry/make_vessels.py stenosis --severity 0.5 -o ...
    uv run scripts/geometry/make_vessels.py curved --offset 20 -o ...

Units are micrometres. Every vessel runs along +x from x = 0 to x = length,
with flat inlet/outlet faces normal to x (for a pre-inlet and outlet boundary
later), and lies in the xy plane. The bounding box is moved to start at 0, so
with <refDir> 0 and <refDirN> = length / dx the voxeliser maps a point p [um]
to lattice node p / dx + 1, and a .pos position q [um] to q / dx: pack the box
printed below with `packCells ... --offset dx dx dx`.
"""

import argparse

import numpy as np
import trimesh

SEGMENTS = 48  # around the circumference


def tube(center, radius):
    """Closed tube along a planar (xy) centerline.

    center: (N, 3) points [um]; radius: (N,) radii [um]. The cross-sections lie
    in the plane spanned by z and the in-plane normal, so the end caps are
    normal to the end tangents.
    """
    t = np.gradient(center, axis=0)
    t /= np.linalg.norm(t, axis=1)[:, None]
    u = np.cross([0.0, 0.0, 1.0], t)
    u /= np.linalg.norm(u, axis=1)[:, None]
    z = np.array([0.0, 0.0, 1.0])
    th = np.linspace(0, 2 * np.pi, SEGMENTS, endpoint=False)
    rings = center[:, None, :] + radius[:, None, None] * (
        np.cos(th)[None, :, None] * u[:, None, :] + np.sin(th)[None, :, None] * z
    )
    n, m = len(center), SEGMENTS
    verts = np.vstack([rings.reshape(-1, 3), center[0], center[-1]])
    faces = []
    for i in range(n - 1):
        for j in range(m):
            a, b = i * m + j, i * m + (j + 1) % m
            faces += [[a, b, a + m], [b, b + m, a + m]]
    c0, c1 = n * m, n * m + 1
    for j in range(m):
        faces += [[c0, (j + 1) % m, j], [c1, (n - 1) * m + j, (n - 1) * m + (j + 1) % m]]
    mesh = trimesh.Trimesh(verts, np.array(faces), process=True)
    mesh.fix_normals()
    assert mesh.is_watertight and mesh.volume > 0, "tube is not a closed volume"
    return mesh


def samples(x0, x1, step):
    return np.linspace(x0, x1, max(2, int(round((x1 - x0) / step)) + 1))


def s_curve(x, xa, xb, h):
    """y rising smoothly (zero slope at both ends) from 0 at xa to h at xb."""
    s = np.clip((x - xa) / (xb - xa), 0, 1)
    return h * (1 - np.cos(np.pi * s)) / 2


def stenosis(a):
    R = a.diameter / 2
    x = samples(0, a.length, a.step)
    xc, ls = a.length / 2, a.stenosis_length
    bump = np.where(np.abs(x - xc) < ls / 2, (1 + np.cos(2 * np.pi * (x - xc) / ls)) / 2, 0.0)
    r = R * (1 - a.severity * bump)
    return tube(np.c_[x, 0 * x, 0 * x], r)


def curved(a):
    x = samples(0, a.length, a.step)
    xa = (a.length - a.bend_length) / 2
    y = s_curve(x, xa, xa + a.bend_length, a.offset)
    return tube(np.c_[x, y, 0 * x], np.full_like(x, a.diameter / 2))


def bifurcation(a):
    R = a.diameter / 2
    rd = R / 2 ** (1 / 3)  # Murray's law: R^3 = 2 rd^3
    xj = a.length - a.branch_length - a.outlet_length  # where the daughters start to diverge
    # Parent: keeps R until the daughters start to diverge, then tapers to 0.7 rd
    # at xj + 16, where its end cap lies inside the daughters (no step in the
    # wall). Chosen numerically for the defaults: the junction narrows to 0.89 R
    # (in y) instead of rd; the rest is the union of thinner circular daughters.
    t0, t1 = xj + 2, xj + 16
    xp = samples(0, t1, a.step)
    taper = np.clip((xp - t0) / (t1 - t0), 0, 1)
    rp = R + (0.7 * rd - R) * (1 - np.cos(np.pi * taper)) / 2
    parts = [tube(np.c_[xp, 0 * xp, 0 * xp], rp)]
    xb = samples(xj - 2, a.length, a.step)
    for sign in (1, -1):
        y = sign * s_curve(xb, xj, xj + a.branch_length, a.offset)
        parts.append(tube(np.c_[xb, y, 0 * xb], np.full_like(xb, rd)))
    mesh = trimesh.boolean.union(parts, engine="manifold")
    assert mesh.is_watertight and mesh.body_count == 1, "branches do not form one vessel"
    return mesh


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("kind", choices=["bifurcation", "stenosis", "curved"])
    p.add_argument("-o", "--output", required=True, help="STL file to write")
    p.add_argument("--diameter", type=float, default=20.0, help="inlet diameter [um] (default 20)")
    p.add_argument("--length", type=float, default=100.0, help="extent along x [um] (default 100)")
    p.add_argument("--dx", type=float, default=0.5, help="lattice spacing [um], for the printed refDirN (default 0.5)")
    p.add_argument("--step", type=float, default=0.5, help="centerline sampling [um] (default 0.5)")
    p.add_argument("--severity", type=float, default=0.5, help="stenosis: radius reduction, 0..1 (default 0.5)")
    p.add_argument("--stenosis_length", type=float, default=30.0, help="stenosis: length [um] (default 30)")
    p.add_argument("--offset", type=float, default=15.0,
                   help="curved: lateral shift of the S-bend; bifurcation: half the outlet separation [um] (default 15)")
    p.add_argument("--bend_length", type=float, default=50.0, help="curved: length of the S-bend [um] (default 50)")
    p.add_argument("--branch_length", type=float, default=40.0, help="bifurcation: length of the diverging part [um]")
    p.add_argument("--outlet_length", type=float, default=20.0, help="bifurcation: straight outlets [um] (default 20)")
    a = p.parse_args()

    mesh = {"bifurcation": bifurcation, "stenosis": stenosis, "curved": curved}[a.kind](a)
    mesh.apply_translation(-mesh.bounds[0])
    assert mesh.is_watertight and mesh.is_winding_consistent and mesh.volume > 0
    size = mesh.bounds[1]
    assert abs(size[0] - a.length) < 1e-6, "x extent must equal --length (refDirN)"
    mesh.export(a.output)

    refdirn = a.length / a.dx
    print(f"{a.output}: {len(mesh.faces)} triangles, lumen {mesh.volume:.0f} um^3")
    print(f"bounding box [um]: {size[0]:.2f} x {size[1]:.2f} x {size[2]:.2f}")
    print(f"config: <refDir> 0 </refDir> <refDirN> {refdirn:g} </refDirN> <dx> {a.dx * 1e-6:g} </dx>")
    print(f"pack:   packCells {size[0]:.2f} {size[1]:.2f} {size[2]:.2f} --walls xyz --offset {a.dx:g} {a.dx:g} {a.dx:g} ...")


if __name__ == "__main__":
    main()
