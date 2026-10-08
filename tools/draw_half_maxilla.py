"""Render the half-maxilla FE model: bone, teeth, attachments, aligner, TAD points.

Bone is drawn cropped to the alveolar region (within DISPLAY_CROP of the teeth) for
legibility only; the deformable-bone model itself uses the whole right half.
"""
import sys, os, pickle
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, 'tools'))
sys.path.insert(0, os.path.join(HERE, 'ofj4d')); sys.path.insert(0, os.path.join(HERE, 'ofj'))
os.chdir(HERE)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import cKDTree
from mesh_io import boundary_faces
from extrude import tet_vol
import spec

DISPLAY_CROP = 9.0
C_BONE, C_TOOTH, C_ATT, C_ALIGN = '#dccbab', '#f3efe3', '#c8791f', '#57c0cd'
C_TAD = {'buccal': '#2a78d6', 'palatal': '#1baf7a'}
INK, INK_2 = '#0b0b0b', '#52514e'


def cam(elev, azim):
    e, a = np.radians(elev), np.radians(azim)
    return np.array([np.cos(e)*np.cos(a), np.cos(e)*np.sin(a), np.sin(e)])


def faces_rgb(pts, tris, color, light):
    v = pts[tris]
    n = np.cross(v[:, 1]-v[:, 0], v[:, 2]-v[:, 0])
    L = np.linalg.norm(n, axis=1); k = L > 1e-12
    v, n = v[k], n[k]/L[k, None]
    b = 0.40 + 0.60*np.abs(n @ light)
    return v, np.clip(np.array(to_rgb(color))[None]*b[:, None], 0, 1)


def draw(ax, layers, elev, azim, lim):
    light = cam(elev, azim)
    V, C = [], []
    for pts, tris, col in layers:
        if len(tris):
            v, c = faces_rgb(pts, tris, col, light); V.append(v); C.append(c)
    pc = Poly3DCollection(np.vstack(V), facecolors=np.vstack(C), edgecolors='none')
    ax.add_collection3d(pc)
    ax.view_init(elev, azim)
    (x0, x1), (y0, y1), (z0, z1) = lim
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_zlim(z0, z1)
    ax.set_box_aspect((x1-x0, y1-y0, z1-z0)); ax.set_axis_off()
    ax.computed_zorder = False


def points(ax, sides, hooks, which):
    for side in which:
        p, h = sides[side], hooks[side]
        ax.plot([h[0], p[0]], [h[1], p[1]], [h[2], p[2]], color=C_TAD[side], lw=2.2, zorder=20)
        ax.scatter([p[0]], [p[1]], [p[2]], s=110, c='#ffd23f', edgecolors='#111111',
                   linewidths=1.6, depthshade=False, zorder=21)


def main():
    G = pickle.load(open('build4d/geom.pkl', 'rb'))
    A = pickle.load(open('build4d/aligner.pkl', 'rb'))
    Bd = pickle.load(open('build4d/bone.pkl', 'rb'))
    T = pickle.load(open('build4d/tad_bone.pkl', 'rb'))
    gn, teeth = G['nodes'], G['teeth']
    tooth_tris = np.vstack([np.vstack([t['root_faces'], t['crown_faces']]) for t in teeth.values()])
    att_tris = np.vstack([t['att']['faces'] for t in teeth.values() if t['att'] is not None])
    al_tris = boundary_faces(A['tets'])

    bf = boundary_faces(Bd['tets'])
    keep = ~Bd['socket'][bf].all(1)
    tpts = gn[np.unique(tooth_tris)]
    d, _ = cKDTree(tpts).query(Bd['nodes'][bf].mean(1))
    bone_tris = bf[keep & (d < DISPLAY_CROP)]

    from tad_check import compute
    R = compute()
    sides = {s: T['sides'][s]['tad'] for s in ('buccal', 'palatal')}
    hooks = {s: R['sides'][s]['hook'] for s in ('buccal', 'palatal')}

    allp = np.vstack([tpts, sides['buccal'][None], sides['palatal'][None]])
    lo, hi = allp.min(0) - 3.5, allp.max(0) + 3.5
    lim = list(zip(lo, hi))

    bone = (Bd['nodes'], bone_tris, C_BONE)
    teeth_l = (gn, tooth_tris, C_TOOTH); att_l = (gn, att_tris, C_ATT)
    align_l = (A['nodes'], al_tris, C_ALIGN)

    fig = plt.figure(figsize=(17, 10.5), facecolor='white')
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.78], left=0.01, right=0.99,
                          top=0.90, bottom=0.02, wspace=0.02, hspace=0.06)
    views = [(gs[0, 0], 'Buccal view — bone, teeth, attachments', [bone, teeth_l, att_l], 12, 158, ['buccal']),
             (gs[0, 1], 'Palatal view', [bone, teeth_l, att_l], 8, -22, ['palatal']),
             (gs[1, 0], 'With the clear aligner', [bone, teeth_l, att_l, align_l], 12, 158, ['buccal']),
             (gs[1, 1], 'Occlusal view, from below — both miniscrews', [bone, teeth_l, att_l], -88, -96, ['buccal', 'palatal'])]
    for cell, title, layers, el, az, which in views:
        ax = fig.add_subplot(cell, projection='3d')
        draw(ax, layers, el, az, lim)
        points(ax, sides, hooks, which)
        ax.set_title(title, fontsize=12.5, color=INK, pad=-2)

    ax = fig.add_subplot(gs[:, 2]); ax.set_axis_off()
    Vb = np.abs(tet_vol(Bd['nodes'], Bd['tets']))
    cort_vol = 100*Vb[Bd['cortical']].sum()/Vb.sum()
    sep = np.linalg.norm(sides['buccal'] - sides['palatal'])
    txt = [
        ('Model', None),
        ('dentition', 'FDI 11–17, right half arch, teeth rigid'),
        ('PDL', f'{spec.PDL_THICKNESS:.2f} mm, 3 layers, bilinear\nE₁ {spec.PDL_E1} / E₂ {spec.PDL_E2} MPa, ν {spec.PDL_NU}'),
        ('aligner', f'{spec.ALIGNER_THICKNESS:.2f} mm, E {spec.E_ALIGNER:.0f} MPa, ν {spec.NU_ALIGNER}'),
        ('attachments', f'{len(spec.ATTACHMENTS)}, E {spec.E_ATTACH/1000:.0f} GPa, rigid with tooth'),
        ('Bone (deformable variant)', None),
        ('source', 'Open-Full-Jaw patient 8 maxilla'),
        ('domain', f'right half, {len(Bd["tets"]):,} tets'),
        ('cortical', f'≤ {Bd["materials"]["cortical_t"]} mm from outer surface\nE {Bd["materials"]["E_cortical"]/1000:.1f} GPa, {cort_vol:.0f}% of bone volume'),
        ('cancellous', f'E {Bd["materials"]["E_cancellous"]/1000:.2f} GPa, ν {Bd["materials"]["nu"]}'),
        ('support', f'superior slab z ≥ {Bd["z_fixed"]:.1f} mm fixed\n(dataset authors’ region);\nmidsagittal symmetry'),
        ('PDL–bone', 'outer PDL layer tied into bone,\nre-tied after each remodelling step'),
        ('Miniscrews (points)', None),
        ('placement', '16/17 interradicular, cortical surface,\n4 mm apical to that side’s crest'),
        ('separation', f'buccal–palatal {sep:.1f} mm'),
        ('elastic', f'{spec.ELASTIC_FORCE_GF:.0f} gf to canine precision cut'),
        ('Size', None),
        ('rigid bone', '150,644 dof'),
        ('deformable bone', '304,499 dof'),
    ]
    y = 0.97
    for k, v in txt:
        if v is None:
            y -= 0.012
            ax.text(0.02, y, k, fontsize=12.5, fontweight='semibold', color=INK, va='top', transform=ax.transAxes)
            y -= 0.038
        else:
            ax.text(0.02, y, k, fontsize=10, color=INK_2, va='top', transform=ax.transAxes)
            ax.text(0.34, y, v, fontsize=10, color=INK, va='top', transform=ax.transAxes)
            y -= 0.030*(1 + v.count('\n'))
    ax.text(0.02, 0.005, f'bone drawn within {DISPLAY_CROP:.0f} mm of the teeth for legibility;\n'
            'the model uses the whole right half', fontsize=8.5, color=INK_2, transform=ax.transAxes)

    fig.suptitle('Half-maxilla finite element model — Open-Full-Jaw patient 8, 4D aligner staging',
                 fontsize=15, fontweight='semibold', color=INK, y=0.975)
    os.makedirs('out4d', exist_ok=True)
    fig.savefig('out4d/half_maxilla_model.png', dpi=115, facecolor='white')
    print('wrote out4d/half_maxilla_model.png')


if __name__ == '__main__':
    main()
