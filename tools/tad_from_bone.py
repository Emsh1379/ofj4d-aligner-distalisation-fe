"""Miniscrew points from the real alveolar bone (Open-Full-Jaw patient 8 maxilla).

Specification (Mao et al. 2024): the TAD "was located in the buccal [palatal]
interradicular space between the first and second molars 4 mm above the alveolar
crest".  The previous rule (inputs.py) took the most occlusal bone node on each side
of the root midline, which lands on the ridge of the septum; both screws then sat
in the middle of the septum, 3.5 mm apart (FINDINGS 19).

Here every step uses the segmented bone and PDL of the dataset:

  1. interradicular septum -- on each side, the mesiodistal gap between the PDL of
     16 and of 17 at the height in question;
  2. alveolar crest -- vertical rays through the bone across that side of the
     septum; the crest height is the highest bone entry;
  3. screw point -- at 4 mm apical to that crest, a ray from the septum centre
     outward (buccally / palatally) marches through the bone until it first leaves
     it: that exit is the cortical surface, where the screw is inserted.

The TAD is a point (the fixed end of the elastic); nothing about a screw body is
assumed.  Output: build4d/tad_bone.pkl.

Position variants (defaults reproduce the paper's site):
    python tools/tad_from_bone.py [--pair 16 17] [--apical 4.0] [--out build4d/tad_bone.pkl]
  --pair    the two adjacent teeth (FDI) whose interradicular space holds the screw
  --apical  mm apical of that side's alveolar crest
Select a variant for a run with OFJ_TAD_FILE=<out> (ofj4d/loads4d.py).
"""
import sys, os, pickle, argparse
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, 'ofj')); sys.path.insert(0, os.path.join(HERE, 'ofj4d'))
os.chdir(HERE)
from scipy.spatial import cKDTree
from mesh_io import load_msh22
from anatomy import load_axes, split_bodies, TAG_TEETH, TAG_PDL, TAG_BONE
from build_frames import orthonormal_frame
import spec

MSH = 'data/P8/maxilla_volumetric_mesh.msh'
AXES = 'data/P8/teeth_axes_maxilla.json'
CACHE = os.environ.get('OFJ_MSH_CACHE', 'data/cache/p8_maxilla.npz')
APICAL = spec.GROUPS['buccal_tad']['tad_above_crest']        # 4.0 mm (default; --apical)
SLICE = 0.75          # mm, half-height of the PDL slice that locates the septum
STEP = 0.05           # mm, ray-march resolution
EXIT_RUN = 0.5        # mm of continuous non-bone needed to call it an exit


def unit(v):
    return v/np.linalg.norm(v)


class BoneQuery:
    """Point-in-bone test on the tetrahedral bone mesh."""

    def __init__(self, nodes, tets):
        X = nodes[tets]
        self.T0 = X[:, 0]
        M = np.stack([X[:, 1]-X[:, 0], X[:, 2]-X[:, 0], X[:, 3]-X[:, 0]], -1)
        det = np.linalg.det(M)
        good = np.abs(det) > 1e-12
        self.Minv = np.full_like(M, np.nan)
        self.Minv[good] = np.linalg.inv(M[good])
        self.tree = cKDTree(X.mean(1))

    def inside(self, P, k=32):
        P = np.atleast_2d(P)
        _, idx = self.tree.query(P, k=k)
        r = P[:, None, :] - self.T0[idx]
        lam = np.einsum('nkij,nkj->nki', self.Minv[idx], r)
        l0 = 1.0 - lam.sum(-1)
        with np.errstate(invalid='ignore'):
            ok = (lam >= -1e-9).all(-1) & (l0 >= -1e-9)
        return ok.any(1)


def label_teeth(nodes, tets, tags, axes):
    """Tooth bodies and PDL components -> UNN, by the same rule as build_frames."""
    tT = tets[tags == TAG_TEETH]
    lab, nb = split_bodies(tT, len(nodes))
    unns = sorted(axes)
    C = np.array([axes[u]['c'] for u in unns])
    cen = np.array([nodes[np.unique(tT[lab == i])].mean(0) for i in range(nb)])
    D = np.linalg.norm(cen[:, None] - C[None], axis=2)
    body_unn = {}
    for i in np.argsort(D.min(1)):
        for j in np.argsort(D[i]):
            if unns[j] not in body_unn.values():
                body_unn[i] = unns[j]; break
    tooth_cen = {body_unn[i]: cen[i] for i in body_unn}
    tP = tets[tags == TAG_PDL]
    plab, npc = split_bodies(tP, len(nodes))
    pdl_nodes = {}
    us = sorted(tooth_cen); TC = np.array([tooth_cen[u] for u in us])
    for i in range(npc):
        pn = nodes[np.unique(tP[plab == i])]
        u = us[int(np.argmin(np.linalg.norm(TC - pn.mean(0), axis=1)))]
        pdl_nodes.setdefault(u, []).append(pn)
    return tooth_cen, {u: np.vstack(v) for u, v in pdl_nodes.items()}


def compute(pair=(16, 17), apical=APICAL):
    nodes, tets, tags = load_msh22(MSH, cache=CACHE if os.path.exists(os.path.dirname(CACHE)) else None)
    axes = load_axes(AXES)
    tooth_cen, pdl = label_teeth(nodes, tets, tags, axes)
    bone = BoneQuery(nodes, tets[tags == TAG_BONE])

    # mesial (m) and distal (d) tooth of the chosen interradicular space; names keep the 16/17 code
    um, ud = sorted((spec.UNN[pair[0]], spec.UNN[pair[1]]), reverse=True)   # higher UNN is mesial on the right side
    ed16, eb16, eo16 = orthonormal_frame(axes[um]['ax_x'], axes[um]['ax_y'], axes[um]['ax_z'])
    ed17, eb17, eo17 = orthonormal_frame(axes[ud]['ax_x'], axes[ud]['ax_y'], axes[ud]['ax_z'])
    eo = unit(eo16)
    es = unit(eb16 + eb17); es = unit(es - (es@eo)*eo)
    ed = unit(ed16 + ed17); ed = unit(ed - (ed@eo)*eo - (ed@es)*es)
    if (tooth_cen[ud] - tooth_cen[um]) @ ed < 0:
        ed = -ed
    mid = 0.5*(tooth_cen[ud] + tooth_cen[um])
    p16, p17 = pdl[um], pdl[ud]

    def septum_md(h, sign):
        """Mesiodistal centre of the 16/17 gap at height h on one side."""
        for w in (SLICE, 1.5, 2.5, 4.0):
            s16 = p16[(np.abs(p16@eo - h) < w) & ((((p16-mid)@es)*sign) > 0)]
            s17 = p17[(np.abs(p17@eo - h) < w) & ((((p17-mid)@es)*sign) > 0)]
            if len(s16) and len(s17):
                lo, hi = (s16@ed).max(), (s17@ed).min()
                return 0.5*(lo + hi), lo, hi, w
        raise RuntimeError('no 16/17 PDL on this side near that height')

    def point(md, h, s):
        """Frame coordinates -> xyz: md along ed, h along eo, s along es from mid."""
        return mid + ed*(md - mid@ed) + eo*(h - mid@eo) + es*s

    out = dict(frame=dict(ed=ed, es=es, eo=eo, mid=mid), sides={})
    top16 = (p16@eo).max()                    # cervical end of the 16 PDL
    for side, sign in (('buccal', 1.0), ('palatal', -1.0)):
        # --- septum centre near the crest, from the cervical 3 mm of the PDL -----
        md_c, _, _, _ = septum_md(top16 - 1.5, sign)
        # --- crest: highest bone entry of vertical rays across this half ---------
        s_grid = np.arange(0.0, 9.0, 0.1)*sign
        hs = np.arange(top16 + 3.0, top16 - 12.0, -STEP)
        prof = []
        for s in s_grid:
            P = np.array([point(md_c, h, s) for h in hs])
            inn = bone.inside(P)
            prof.append(hs[np.argmax(inn)] if inn.any() else np.nan)
        prof = np.array(prof)
        crest_h = np.nanmax(prof)
        crest_s = s_grid[np.nanargmax(prof)]
        target_h = crest_h - apical
        # --- screw point: march outward at the target height until bone ends -----
        md_t, lo, hi, w = septum_md(target_h, sign)
        s_ray = np.arange(0.0, 20.0, STEP)*sign
        P = np.array([point(md_t, target_h, s) for s in s_ray])
        inn = bone.inside(P)
        if not inn.any():
            raise RuntimeError(f'{side}: ray at the target height never enters bone')
        first_in = int(np.argmax(inn))
        run = int(round(EXIT_RUN/STEP))
        exit_i = None
        for i in range(first_in, len(inn) - run):
            if not inn[i:i+run].any():
                exit_i = i; break
        if exit_i is None:
            raise RuntimeError(f'{side}: no cortical exit found within 20 mm')
        # refine the exit between the last bone sample and the first non-bone one
        a, b = s_ray[exit_i-1], s_ray[exit_i]
        for _ in range(20):
            m = 0.5*(a + b)
            if bone.inside(point(md_t, target_h, m))[0]:
                a = m
            else:
                b = m
        s_exit = 0.5*(a + b)
        tad = point(md_t, target_h, s_exit)
        root16 = float(np.linalg.norm(p16 - tad, axis=1).min())
        root17 = float(np.linalg.norm(p17 - tad, axis=1).min())
        slice_all = np.vstack([p16, p17])
        k = (np.abs(slice_all@eo - target_h) < SLICE) & ((((slice_all-mid)@es)*sign) > 0)
        outer_root = float(((((slice_all[k]-mid)@es)*sign)).max()) if k.any() else np.nan
        out['sides'][side] = dict(
            tad=tad, crest_height=float(crest_h), crest_s=float(crest_s),
            target_height=float(target_h), apical_of_crest=float(crest_h - target_h),
            septum_md=float(md_t), gap=(float(lo), float(hi)), slice=w,
            s_exit=float(abs(s_exit)), outer_root=outer_root,
            beyond_root=float(abs(s_exit) - outer_root),
            clear16=float(md_t - lo), clear17=float(hi - md_t),
            dist_pdl16=root16, dist_pdl17=root17, profile=(s_grid, prof),
            first_bone_s=float(abs(s_ray[first_in])))
    return out, nodes, tets, tags, pdl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pair', type=int, nargs=2, default=[16, 17])
    ap.add_argument('--apical', type=float, default=APICAL)
    ap.add_argument('--out', default='build4d/tad_bone.pkl')
    a = ap.parse_args()
    R, nodes, tets, tags, pdl = compute(tuple(a.pair), a.apical)
    fr = R['frame']
    B, P = R['sides']['buccal'], R['sides']['palatal']
    bar = '='*78
    print(f"\n{bar}\nMINISCREW POINTS FROM THE REAL BONE  (FDI {a.pair[0]}/{a.pair[1]}, {a.apical:g} mm apical)\n{bar}")
    print(f"  {'':<48}{'buccal':>12}{'palatal':>12}")
    rows = [('alveolar crest height (occlusal axis)', 'crest_height', '{:.2f}'),
            ('screw height, apical of crest', 'apical_of_crest', '{:.2f} mm'),
            ('interradicular gap at screw height', None, None),
            ('clearance to mesial tooth PDL, mesiodistal', 'clear16', '{:.2f} mm'),
            ('clearance to distal tooth PDL, mesiodistal', 'clear17', '{:.2f} mm'),
            ('cortical surface, from ridge centre', 's_exit', '{:.2f} mm'),
            ('outermost root/PDL on that side', 'outer_root', '{:.2f} mm'),
            ('screw beyond the roots', 'beyond_root', '{:+.2f} mm'),
            ('nearest PDL of mesial tooth (3D)', 'dist_pdl16', '{:.2f} mm'),
            ('nearest PDL of distal tooth (3D)', 'dist_pdl17', '{:.2f} mm')]
    for label, k, fmt in rows:
        if k is None:
            print(f"  {label:<48}{B['gap'][1]-B['gap'][0]:>9.2f} mm{P['gap'][1]-P['gap'][0]:>9.2f} mm")
            continue
        print(f"  {label:<48}{fmt.format(B[k]):>12}{fmt.format(P[k]):>12}")
    sep = abs((B['tad'] - P['tad']) @ fr['es'])
    print(f"\n  buccal-to-palatal separation of the screws: {sep:.2f} mm "
          f"(3D {np.linalg.norm(B['tad']-P['tad']):.2f} mm)")

    # the rule these replace (inputs.py / slim.pkl)
    slim = pickle.load(open('build4d/slim.pkl', 'rb'))
    eo_old = unit(np.asarray(slim['e_o_16']))
    print(f"\n  previous rule (septum-ridge crest, straight apical):")
    for side, S in (('buccal', B), ('palatal', P)):
        old = np.asarray(slim['crest'][side]) - eo_old*APICAL
        print(f"    {side:<8} old {np.round(old, 2)}  new {np.round(S['tad'], 2)}  "
              f"moved {np.linalg.norm(S['tad']-old):.2f} mm")
    os.makedirs('build4d', exist_ok=True)
    pickle.dump(dict(frame=fr, sides={k: {kk: vv for kk, vv in v.items() if kk != 'profile'}
                                      for k, v in R['sides'].items()},
                     profiles={k: v['profile'] for k, v in R['sides'].items()},
                     source=MSH, rule=__doc__.split('\n\n')[2], pair=tuple(a.pair), apical=a.apical),
                open(a.out, 'wb'))
    print(f"\n  wrote {a.out}")


if __name__ == '__main__':
    main()
