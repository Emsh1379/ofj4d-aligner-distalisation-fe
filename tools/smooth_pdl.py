"""PDL mesh without folded prisms, for finite-strain solvers (FEBio) and the custom solver alike.

ofj4d/geom4d.extrude_layers offsets the root surface along its vertex normals. Where the root is
concave on the scale of the 0.30 mm PDL, neighbouring normals diverge, the prisms twist, and some of
their three tets fold over their neighbours (29 folded prisms on patient 8; docs/FINDINGS.md section 31).
This tool rebuilds the middle and outer layers along extrusion directions smoothed over the root
surface (Laplacian, `--iterations`), keeping each node's PDL thickness; the outer layer (the fixed bone
surface) moves. Before that it removes root-surface micro-triangles by collapsing interior edges shorter
than `--collapse` x the median edge (to their midpoint, only when the link condition keeps the surface
manifold, never on the boundary): extruded 0.3 mm over a triangle a few percent the size of its
neighbours they give needle tets and free nodes with ~1% of the usual support, which invert first under
finite strain. Collapsed root vertices move by at most half a short edge; all others are unchanged.

usage: python tools/smooth_pdl.py [--iterations 20] [--collapse 0.35] [--geom build4d/geom.pkl] [--out build4d/geom_fe.pkl]
"""
import argparse, copy, pickle
import numpy as np
from scipy.sparse import coo_matrix

ap = argparse.ArgumentParser()
ap.add_argument('--iterations', type=int, default=20)
ap.add_argument('--collapse', type=float, default=0.35)
ap.add_argument('--layers', type=int, default=0, help='element layers through the PDL (0 = as the source)')
ap.add_argument('--subdivide', type=int, default=0, help='1->4 midpoint subdivisions of the root surface (mesh convergence)')
ap.add_argument('--geom', default='build4d/geom.pkl'); ap.add_argument('--out', default='build4d/geom_fe.pkl')
a = ap.parse_args()


def vol6(X, T):
    P = X[T]
    return np.einsum('ij,ij->i', np.cross(P[:, 1]-P[:, 0], P[:, 2]-P[:, 0]), P[:, 3]-P[:, 0])


def quality(X, T):
    """altitude / longest edge (0 = flat, regular tet 0.667)"""
    P = X[T]; v = np.abs(vol6(X, T))/6
    A = np.stack([0.5*np.linalg.norm(np.cross(P[:, b]-P[:, c0], P[:, c]-P[:, c0]), axis=1)
                  for c0, b, c in [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]], 1).max(1)
    L = np.stack([np.linalg.norm(P[:, i]-P[:, j], axis=1) for i in range(4) for j in range(i+1, 4)], 1).max(1)
    return 3*v/A/L


def extrude(inner, f, d, h, layers):
    """geom4d.extrude_layers with per-node direction and thickness; canonical (unoriented) tets too."""
    nv = len(inner)
    X = np.vstack([inner + d*(h[:, None]*L/layers) for L in range(layers+1)])
    T = []
    for L in range(layers):
        b, t = f + L*nv, f + (L+1)*nv
        i, j, k = b.T; I, J, K = t.T
        T += [np.stack([i, j, k, K], 1), np.stack([i, j, J, K], 1), np.stack([i, I, J, K], 1)]
    T = np.vstack(T)
    # folded prisms: canonical signs must be (+, -, +) times the prism orientation
    sg = np.sign(vol6(X, T)).reshape(layers, 3, len(f))
    nb = np.cross(inner[f[:, 1]]-inner[f[:, 0]], inner[f[:, 2]]-inner[f[:, 0]])
    sp = np.sign(np.einsum('ij,ij->i', nb, (d[f]*h[f][:, :, None]).mean(1)))
    folded = int((sg != np.array([1, -1, 1])[None, :, None]*sp[None, None, :]).any(1).sum())
    Tpos = T.copy(); s = vol6(X, T); Tpos[s < 0] = Tpos[s < 0][:, [0, 1, 3, 2]]
    return X, Tpos, T, folded


def collapse_short_edges(inner, d, h, f, ratio):
    """Collapse interior root-surface edges shorter than ratio x median edge. Returns new inner, d, h, f, count."""
    X = inner.copy(); D = d.copy(); H = h.copy(); F = f.copy()
    alive = np.ones(len(X), bool); ncol = 0
    e0 = np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    lmin = ratio*np.median(np.linalg.norm(X[e0[:, 0]] - X[e0[:, 1]], axis=1))
    for _ in range(50):
        e = np.sort(np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]), 1)
        eu, ec = np.unique(e, axis=0, return_counts=True)
        bverts = set(np.unique(eu[ec == 1]).tolist())
        L = np.linalg.norm(X[eu[:, 0]] - X[eu[:, 1]], axis=1)
        cand = np.argsort(L); cand = cand[L[cand] < lmin]
        if len(cand) == 0:
            break
        nbr = {}
        for a_, b_ in eu:
            nbr.setdefault(a_, set()).add(b_); nbr.setdefault(b_, set()).add(a_)
        touched = set(); done = 0
        for k in cand:
            a_, b_ = eu[k]
            if a_ in touched or b_ in touched or a_ in bverts or b_ in bverts or ec[k] != 2:
                continue
            if len(nbr[a_] & nbr[b_]) != 2:                      # link condition: keeps the surface manifold
                continue
            X[a_] = 0.5*(X[a_] + X[b_]); D[a_] = D[a_] + D[b_]; D[a_] /= np.linalg.norm(D[a_]); H[a_] = 0.5*(H[a_] + H[b_])
            F[F == b_] = a_; alive[b_] = False
            touched |= {a_, b_} | nbr[a_] | nbr[b_]; done += 1
        F = F[(F[:, 0] != F[:, 1]) & (F[:, 1] != F[:, 2]) & (F[:, 0] != F[:, 2])]
        ncol += done
        if done == 0:
            break
    keep = np.where(alive)[0]; rm = np.full(len(X), -1); rm[keep] = np.arange(len(keep))
    return X[keep], D[keep], H[keep], np.sort(rm[F], 1), ncol


def support_stats(X, T, nv, layers):
    vol = np.abs(vol6(X, T))/6
    sup = np.zeros(len(X)); np.add.at(sup, T.ravel(), np.repeat(vol, 4)/4)
    free = np.arange(nv, layers*nv)
    med = np.median(sup[free])
    return int(np.sum(sup[free] < 0.05*med)), float(sup[free].min()/med)


def subdivide(inner, d, h, f):
    """1->4 midpoint subdivision of the root surface; direction and thickness interpolated."""
    e = np.sort(np.vstack([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]]), 1)
    eu, inv = np.unique(e, axis=0, return_inverse=True)
    nv = len(inner); m = nv + np.arange(len(eu))
    Xm = 0.5*(inner[eu[:, 0]] + inner[eu[:, 1]])
    Dm = d[eu[:, 0]] + d[eu[:, 1]]; Dm /= np.linalg.norm(Dm, axis=1, keepdims=True)
    Hm = 0.5*(h[eu[:, 0]] + h[eu[:, 1]])
    nf = len(f); mab, mbc, mca = m[inv[:nf]], m[inv[nf:2*nf]], m[inv[2*nf:]]
    a_, b_, c_ = f.T
    F = np.vstack([np.stack([a_, mab, mca], 1), np.stack([mab, b_, mbc], 1),
                   np.stack([mca, mbc, c_], 1), np.stack([mab, mbc, mca], 1)])
    return np.vstack([inner, Xm]), np.vstack([d, Dm]), np.concatenate([h, Hm]), np.sort(F, 1)


def smooth(d, f, it):
    n = len(d)
    e = np.vstack([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]]); e = np.vstack([e, e[:, ::-1]])
    W = coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])), shape=(n, n)).tocsr(); W.data[:] = 1
    deg = np.asarray(W.sum(1)).ravel()
    for _ in range(it):
        d = 0.5*d + 0.5*(W @ d)/deg[:, None]
        d /= np.linalg.norm(d, axis=1, keepdims=True)
    return d


G = pickle.load(open(a.geom, 'rb')); G2 = copy.deepcopy(G)
tot = dict(folded0=0, folded1=0, q001_0=0, q001_1=0, q003_0=0, q003_1=0); shifts = []
for u, t in G['teeth'].items():
    pn, pt, lids = t['pdl_local']; nv = len(lids[0]); layers = len(lids) - 1
    assert all(np.array_equal(l, np.arange(nv) + L*nv) for L, l in enumerate(lids)), 'unexpected layer layout'
    nf = len(pt)//(3*layers)
    blk = pt[:nf]
    f = np.sort(np.where(blk < nv, blk, np.iinfo(np.int64).max), 1)[:, :3]
    inner = pn[lids[0]]
    dv = pn[lids[-1]] - inner; h = np.linalg.norm(dv, axis=1); d0 = dv/h[:, None]
    X0, T0, C0, fold0 = extrude(inner, f, d0, h, layers)
    assert np.allclose(X0, pn) and np.array_equal(np.sort(T0, 1), np.sort(pt, 1)), 'could not reproduce the shipped PDL'
    w0, m0 = support_stats(X0, T0, nv, layers)
    inner1, d1, h1, f1, ncol = collapse_short_edges(inner, d0, h, f, a.collapse) if a.collapse > 0 else (inner, d0, h, f, 0)
    for _ in range(a.subdivide):
        inner1, d1, h1, f1 = subdivide(inner1, d1, h1, f1)
    nv1 = len(inner1); layers1 = a.layers or layers
    # 4^n iterations keep the same physical smoothing length at mesh size h/2^n
    X1, T1, C1, fold1 = extrude(inner1, f1, smooth(d1, f1, a.iterations*(4**a.subdivide)), h1, layers1)
    lids1 = [np.arange(nv1) + L*nv1 for L in range(layers1+1)]
    w1, m1 = support_stats(X1, T1, nv1, layers1)
    q0, q1 = quality(X0, T0), quality(X1, T1)
    _, nn_ = __import__('scipy.spatial', fromlist=['cKDTree']).cKDTree(X0[lids[-1]]).query(X1[lids1[-1]])
    sh = np.linalg.norm(X1[lids1[-1]] - X0[lids[-1]][nn_], axis=1); shifts.append(sh)
    v0, v1 = np.abs(vol6(X0, T0)).sum()/6, np.abs(vol6(X1, T1)).sum()/6
    print(f"FDI {t['fdi']}: collapsed edges {ncol}, root vertices {nv} -> {nv1}; weak free nodes (<5% support) {w0} -> {w1}, min support {m0:.4f} -> {m1:.4f}; folded prisms {fold0} -> {fold1}; q<0.01 {np.sum(q0<0.01)} -> {np.sum(q1<0.01)}, "
          f"q<0.03 {np.sum(q0<0.03)} -> {np.sum(q1<0.03)}, min q {q0.min():.4f} -> {q1.min():.4f}; "
          f"outer-layer shift median {np.median(sh)*1e3:.1f} um, max {sh.max()*1e3:.1f} um (nearest-node); volume {v0:.3f} -> {v1:.3f} mm3")
    tot['folded0'] += fold0; tot['folded1'] += fold1
    tot['q001_0'] += int(np.sum(q0 < 0.01)); tot['q001_1'] += int(np.sum(q1 < 0.01))
    tot['q003_0'] += int(np.sum(q0 < 0.03)); tot['q003_1'] += int(np.sum(q1 < 0.03))
    G2['teeth'][u]['pdl_local'] = (X1, T1, lids1)
    G2['teeth'][u]['pdl_vol'] = float(v1)
s = np.concatenate(shifts)*1e3
print(f"all teeth: folded prisms {tot['folded0']} -> {tot['folded1']}, q<0.01 {tot['q001_0']} -> {tot['q001_1']}, "
      f"q<0.03 {tot['q003_0']} -> {tot['q003_1']}; outer-layer shift median {np.median(s):.1f} um, p95 {np.percentile(s, 95):.1f}, max {s.max():.1f}")
G2['pdl_smoothing'] = dict(tool='tools/smooth_pdl.py', iterations=a.iterations, collapse=a.collapse, layers=a.layers, subdivide=a.subdivide, source=a.geom)
pickle.dump(G2, open(a.out, 'wb')); print('->', a.out)
