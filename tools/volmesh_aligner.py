"""0.70 mm aligner as a true volumetric tetrahedral mesh (no folded elements), for FEBio and the custom solver.

ofj4d/aligner4d.py extrudes the inner surface along per-vertex offset directions. Where the appliance
wraps embrasures and attachment edges with curvature radius below its 0.70 mm thickness, those prisms
must self-intersect: 4,163 folded prisms and 1,083 slivers on patient 8, overlapping volume
(docs/FINDINGS.md section 31). Smoothing the directions cannot remove this.

Here the shell is the region between the same signed-distance level sets aligner4d uses, 0 <= phi <= 0.70
mm, trimmed at the gingival margin by aligner4d's rule. Its boundary is one closed isosurface of
psi = max(-phi, phi - 0.70, gingival cut), which closes the rim; TetGen meshes the enclosed volume with
that surface preserved (switch Y). Inner-surface faces are those on phi = 0.

Output keeps the aligner.pkl contract of ofj4d/model4d.py: inner-surface nodes first (n_inner),
inner_faces, vert_tooth; then the other surface nodes (outer surface and rim) up to n_layer_nodes, so
loads4d's precision-cut search stays on the surface; interior nodes last.

TetGen (https://wias-berlin.de/software/tetgen/) is AGPL-licensed; the `tetgen` Python wrapper is MIT.

usage: python tools/volmesh_aligner.py [--switches pq1.5/15YQ] [--geom build4d/geom.pkl] [--out build4d/aligner_fe.pkl]
"""
import sys, os, time, pickle, argparse
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [os.path.join(HERE, 'ofj4d'), os.path.join(HERE, 'ofj')]
from scipy import ndimage
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from skimage import measure
from build_aligner import rasterize_tris, fill_interior
import aligner4d as A4
import spec
import tetgen

ap = argparse.ArgumentParser()
ap.add_argument('--switches', default='pq1.5/15YQ')
ap.add_argument('--step', type=int, default=A4.STEP, help='marching-cubes step in voxels (3 = 0.45 mm surface)')
ap.add_argument('--geom', default=os.path.join(HERE, 'build4d/geom.pkl'))
ap.add_argument('--out', default=os.path.join(HERE, 'build4d/aligner_fe.pkl'))
a = ap.parse_args()
t0 = time.time()
THK = spec.ALIGNER_THICKNESS


def quality(X, T):
    P = X[T]
    vol = np.einsum('ij,ij->i', np.cross(P[:, 1]-P[:, 0], P[:, 2]-P[:, 0]), P[:, 3]-P[:, 0])/6
    A = np.stack([0.5*np.linalg.norm(np.cross(P[:, b]-P[:, c0], P[:, c]-P[:, c0]), axis=1)
                  for c0, b, c in [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]], 1).max(1)
    L = np.stack([np.linalg.norm(P[:, i]-P[:, j], axis=1) for i in range(4) for j in range(i+1, 4)], 1).max(1)
    return 3*np.abs(vol)/A/L, vol


# ---- the signed distance of aligner4d.build ------------------------------------------------------
G = pickle.load(open(a.geom, 'rb')); nodes, teeth = G['nodes'], G['teeth']
tri = []
for u, t in teeth.items():
    tri += [t['root_faces'], t['crown_faces']] + ([t['att']['faces']] if t['att'] is not None else [])
tri = np.vstack(tri)
pts = nodes[np.unique(tri)]
lo = pts.min(0) - A4.PAD; hi = pts.max(0) + A4.PAD
shape = tuple(np.ceil((hi - lo)/A4.VOX).astype(int) + 1)
Sv = rasterize_tris(nodes[tri], lo, shape, A4.VOX); solid = fill_interior(Sv); del Sv
d_out = ndimage.distance_transform_edt(~solid, sampling=A4.VOX); D = d_out <= A4.R_CLOSE; del d_out, solid
d_in = ndimage.distance_transform_edt(D, sampling=A4.VOX); closed = d_in > A4.R_CLOSE; del d_in, D
phi = (ndimage.distance_transform_edt(~closed, sampling=A4.VOX)
       - ndimage.distance_transform_edt(closed, sampling=A4.VOX)).astype(np.float32)
del closed

# ---- gingival trim of aligner4d: nearest crown/attachment point's tooth, keep e_o . x >= gingival level
pa, ti = [], []
for u, t in teeth.items():
    cn = np.unique(t['crown_faces']); pa.append(nodes[cn]); ti.append(np.full(len(cn), u))
    if t['att'] is not None:
        an = np.unique(t['att']['faces']); pa.append(nodes[an]); ti.append(np.full(len(an), u))
pa = np.vstack(pa); ti = np.concatenate(ti); tree = cKDTree(pa)
EOt = {u: t['e_o'] for u, t in teeth.items()}; GGt = {u: t['gingival_level'] for u, t in teeth.items()}


def gingival(X):
    _, nn = tree.query(X, workers=-1); vt = ti[nn]
    EO = np.array([EOt[u] for u in vt]); GG = np.array([GGt[u] for u in vt])
    return GG - np.einsum('ij,ij->i', X, EO), vt


band = (phi > -0.4) & (phi < THK + 0.4)
g, _ = gingival(lo + np.argwhere(band)*A4.VOX)
psi = np.full(shape, 1.0, np.float32)
pv = phi[band]
psi[band] = np.maximum(np.maximum(-pv, pv - THK), g.astype(np.float32))
del band, g, pv

# ---- closed shell surface --------------------------------------------------------------------------
verts, faces, _, _ = measure.marching_cubes(psi, level=0.0, spacing=(A4.VOX,)*3, step_size=a.step)
del psi
verts = verts + lo
e = np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
nc, lab = connected_components(coo_matrix((np.ones(len(e), np.int8), (e[:, 0], e[:, 1])), shape=(len(verts),)*2), directed=False)
cnt = np.bincount(lab[faces[:, 0]], minlength=nc)
faces = faces[lab[faces[:, 0]] == cnt.argmax()]
used = np.unique(faces); rm = np.full(len(verts), -1); rm[used] = np.arange(len(used)); verts = verts[used]; faces = rm[faces]
ek = np.sort(np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), 1)
_, ec = np.unique(ek, axis=0, return_counts=True)
assert np.all(ec == 2), f'shell surface is not closed ({np.sum(ec != 2)} bad edges)'
print(f"closed shell surface: {len(verts)} verts, {len(faces)} tris  [{time.time()-t0:.0f} s]", flush=True)

# ---- volume mesh ---------------------------------------------------------------------------------
res = tetgen.TetGen(verts, faces).tetrahedralize(switches=a.switches)
X, T = np.asarray(res[0], float), np.asarray(res[1])[:, :4].astype(np.int64)
assert np.allclose(X[:len(verts)], verts), 'TetGen did not keep the surface nodes first'
q, vol = quality(X, T)
T[vol < 0] = T[vol < 0][:, [0, 1, 3, 2]]
fc = np.sort(np.vstack([T[:, [0, 1, 2]], T[:, [0, 1, 3]], T[:, [0, 2, 3]], T[:, [1, 2, 3]]]), 1)
fu, fcount = np.unique(fc, axis=0, return_counts=True)
assert fcount.max() <= 2, 'non-conforming tetrahedralisation'
bnd = fu[fcount == 1]
assert len(bnd) == len(faces) and np.array_equal(np.unique(bnd), np.arange(len(verts))), 'boundary is not the input surface'
print(f"TetGen {a.switches}: {len(X)} nodes, {len(T)} tets, volume {np.abs(vol).sum():.1f} mm3; quality min {q.min():.4f}, "
      f"q<0.01 {np.sum(q<0.01)}, q<0.03 {np.sum(q<0.03)}, median {np.median(q):.3f}  [{time.time()-t0:.0f} s]", flush=True)

# ---- inner surface: surface faces on phi = 0 -----------------------------------------------------------
ijk = ((verts - lo)/A4.VOX).T
phiv = ndimage.map_coordinates(phi, ijk, order=1, mode='nearest')
del phi
# a surface vertex belongs to the bounding level set it is closest to: inner (phi = 0), outer
# (phi = thickness) or the gingival cut.  A tolerance on |phi| alone fails: marching cubes samples
# every STEP voxels, so vertices miss the level sets by a fraction of 0.45 mm.
gv, _ = gingival(verts)
dist = np.stack([np.abs(phiv), np.abs(phiv - THK), np.abs(gv)], 1)
kind = dist.argmin(1)                                   # 0 inner, 1 outer, 2 rim
inner_f = faces[(kind[faces] == 0).all(1)]
inner_nodes = np.unique(inner_f)
# orient inner faces with the outward normal of the TEETH (as aligner4d: along grad phi)
v = verts[inner_f]; fn = np.cross(v[:, 1]-v[:, 0], v[:, 2]-v[:, 0])
_, nn = tree.query(v.mean(1), workers=-1)
if np.einsum('ij,ij->i', fn, v.mean(1) - pa[nn]).sum() < 0:
    inner_f = inner_f[:, [0, 2, 1]]

# ---- node order: inner surface, other surface, interior ----------------------------------------------
surf = np.arange(len(verts))
other = np.setdiff1d(surf, inner_nodes)
interior = np.arange(len(verts), len(X))
order = np.concatenate([inner_nodes, other, interior])
new = np.empty(len(X), np.int64); new[order] = np.arange(len(X))
Xn, Tn, Fn = X[order], new[T], new[inner_f]
n_inner = len(inner_nodes)
_, vt = gingival(Xn[:n_inner])
_, nnv = tree.query(Xn[:n_inner], workers=-1); vt = ti[nnv]
area = 0.5*np.linalg.norm(np.cross(Xn[Fn[:, 1]]-Xn[Fn[:, 0]], Xn[Fn[:, 2]]-Xn[Fn[:, 0]]), axis=1).sum()
outer_area = 0.5*np.linalg.norm(np.cross(verts[faces[:, 1]]-verts[faces[:, 0]], verts[faces[:, 2]]-verts[faces[:, 0]]), axis=1)[(kind[faces] == 1).all(1)].sum()
print(f"surface vertices: inner {np.sum(kind==0)}, outer {np.sum(kind==1)}, rim {np.sum(kind==2)}")
print(f"inner surface: {n_inner} nodes, {len(Fn)} faces, {area:.0f} mm2; outer surface {outer_area:.0f} mm2; other surface nodes {len(other)}; "
      f"interior {len(interior)}; mean thickness 2 x volume / (inner + outer area) {2*np.abs(vol).sum()/(area + outer_area):.3f} mm", flush=True)

A = dict(nodes=Xn, tets=Tn, inner_faces=Fn, n_inner=n_inner, n_layer_nodes=n_inner + len(other),
         vert_tooth=vt, thickness=THK, inner_verts=Xn[:n_inner],
         volmesh=dict(tool='tools/volmesh_aligner.py', switches=a.switches, step=a.step, geom=a.geom,
                      quality_min=float(q.min()), q_below_001=int(np.sum(q < 0.01)), volume=float(np.abs(vol).sum())))
pickle.dump(A, open(a.out, 'wb'))
print(f"-> {a.out}  [{time.time()-t0:.0f} s]")
