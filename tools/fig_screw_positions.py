"""Figure: miniscrew positions of the position study on the patient-8 anatomy.

usage: python tools/fig_screw_positions.py [--out out/positions/screw_positions.png]

Three views in the frame of the 16/17 interradicular space (tools/tad_from_bone.py): occlusal (mesiodistal x
buccopalatal), buccal side (mesiodistal x vertical) and palatal side. Teeth are the crown and root surface points
of the model (grey), bone the outer bone surface of the dataset where the mesh cache exists. Circles = height
series at 16/17 (2 / 4 / 6 / 8 mm), triangles = site series at 4 mm (16/17 / 15/16 / 14/15), each an ordered blue
ramp and labelled directly. The miniscrew is a point in the model; marker size carries no meaning.
"""
import sys, os, pickle, argparse
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, 'ofj4d')); sys.path.insert(0, os.path.join(HERE, 'ofj')); os.chdir(HERE)
os.environ.setdefault('OFJ_GPU', 'off')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import spec

ap = argparse.ArgumentParser()
ap.add_argument('--out', default='out/positions/screw_positions.png')
a = ap.parse_args()
os.makedirs(os.path.dirname(a.out), exist_ok=True)
SURFACE, INK, INK_2, MUTED, GRID, TOOTH, BONE = '#fcfcfb', '#0b0b0b', '#52514e', '#8a8983', '#e4e3de', '#c9c8c2', '#ecebe6'
RAMP4 = ['#86b6ef', '#3987e5', '#1c5cab', '#0d366b']
RAMP3 = ['#86b6ef', '#2a78d6', '#104281']
HEIGHT = [('h1617_a2', '2 mm'), ('baseline', '4 mm'), ('h1617_a6', '6 mm'), ('h1617_a8', '8 mm')]
SITE = [('baseline', '16/17'), ('md1516_a4', '15/16'), ('md1415_a4', '14/15')]


def tad(v):
    return pickle.load(open('build4d/tad_bone.pkl' if v == 'baseline' else f'build4d/tad_positions/tad_{v}.pkl', 'rb'))


fr = tad('baseline')['frame']
ed, es, eo, mid = (np.asarray(fr[k]) for k in ('ed', 'es', 'eo', 'mid'))
to2 = lambda P: np.stack([(P - mid) @ ed, (P - mid) @ es, (P - mid) @ eo], -1)     # mesiodistal, buccal +, occlusal +

G = pickle.load(open('build4d/geom.pkl', 'rb'))
teeth = []
for u, t in G['teeth'].items():
    idx = np.unique(np.vstack([t['crown_faces'], t['root_faces']]))
    teeth.append(to2(G['nodes'][idx]))
teeth = np.vstack(teeth)
bone = None
try:
    from mesh_io import boundary_faces
    from anatomy import TAG_BONE
    c = np.load('data/cache/p8_maxilla.npz'); f, own = boundary_faces(c['tets'], return_owner=True)
    bone = to2(c['nodes'][np.unique(f[c['tags'][own] == TAG_BONE])])
except Exception:
    pass

fig, axes = plt.subplots(1, 3, figsize=(15, 5.6), facecolor=SURFACE, gridspec_kw=dict(width_ratios=[1.15, 1, 1]))
fig.subplots_adjust(left=0.05, right=0.99, top=0.83, bottom=0.12, wspace=0.22)
views = [('Occlusal view', 0, 1, 'Mesiodistal (mm, distal +)', 'Buccopalatal (mm, buccal +)', None),
         ('Buccal side', 0, 2, 'Mesiodistal (mm, distal +)', 'Vertical (mm, occlusal +)', 'buccal'),
         ('Palatal side', 0, 2, 'Mesiodistal (mm, distal +)', 'Vertical (mm, occlusal +)', 'palatal')]
XLIM, YLIM_OCC, YLIM_SIDE = (-18, 14), (-12, 12), (-12, 12)
for ax, (title, i, j, xl, yl, side) in zip(axes, views):
    ax.set_facecolor(SURFACE)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'):
        ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK_2, labelsize=8.5, length=3)
    sel = slice(None)
    if side is not None:
        sgn = 1 if side == 'buccal' else -1
        sel_t = teeth[:, 1]*sgn > -2.0
    else:
        sel_t = np.ones(len(teeth), bool)
    ylim = YLIM_OCC if side is None else YLIM_SIDE
    inwin = lambda P: (P[:, i] > XLIM[0]) & (P[:, i] < XLIM[1]) & (P[:, j] > ylim[0]) & (P[:, j] < ylim[1])
    if bone is not None:
        kb = inwin(bone) & ((bone[:, 1]*(1 if side == 'buccal' else -1) > 2.0) if side else (bone[:, 2] > -8))
        ax.scatter(bone[kb, i], bone[kb, j], s=0.3, color=BONE, zorder=1, rasterized=True)
    kt = sel_t & inwin(teeth)
    ax.scatter(teeth[kt, i], teeth[kt, j], s=0.4, color=TOOTH, zorder=2, rasterized=True)
    for series, marker, ramp, kind in ((HEIGHT, 'o', RAMP4, 'height'), (SITE, '^', RAMP3, 'site')):
        for (v, name), col in zip(series, ramp):
            T = tad(v)
            for sd in (['buccal', 'palatal'] if side is None else [side]):
                p = to2(np.asarray(T['sides'][sd]['tad'])[None])[0]
                if kind == 'site' and v == 'baseline':
                    continue                                   # the paper's site is drawn once, as the 4 mm circle
                ax.scatter(p[i], p[j], s=70, marker=marker, color=col, edgecolors=SURFACE, linewidths=1.5, zorder=4)
                if side is None:                                   # occlusal: heights overlap in plan
                    if kind == 'site' and v != 'baseline':
                        ax.annotate(name, (p[i], p[j]), xytext=(p[i] - 1.0, p[j] + (1.1 if sd == 'buccal' else -1.8)),
                                    fontsize=7.5, color=INK_2, zorder=5)
                    elif kind == 'height' and v == 'h1617_a8':
                        ax.annotate('16/17: 2-8 mm (paper 4 mm)', (p[i], p[j]), xytext=(p[i] + 1.2, p[j] + (1.1 if sd == 'buccal' else -1.8)),
                                    fontsize=7.5, color=INK_2, zorder=5)
                else:
                    if kind == 'height':
                        lab = name if v != 'baseline' else '4 mm (paper)'
                        ax.annotate(lab, (p[i], p[j]), xytext=(p[i] + 0.9, p[j] - 0.25), fontsize=7.5, color=INK_2, zorder=5)
                    elif v != 'baseline':
                        ax.annotate(name, (p[i], p[j]), xytext=(p[i] - 1.1, p[j] - 1.3), fontsize=7.5, color=INK_2, zorder=5)
    ax.set_title(title, color=INK, fontsize=11, loc='left')
    ax.set_xlabel(xl, color=INK_2, fontsize=9); ax.set_ylabel(yl, color=INK_2, fontsize=9)
    ax.set_xlim(*XLIM); ax.set_ylim(*ylim); ax.set_aspect('equal', adjustable='box')
from matplotlib.lines import Line2D
handles = [Line2D([], [], ls='', marker='o', color=RAMP4[1], markersize=8, label='height series at 16/17 (2, 4, 6, 8 mm)'),
           Line2D([], [], ls='', marker='^', color=RAMP3[1], markersize=8, label='site series at 4 mm (15/16, 14/15; 16/17 = the 4 mm circle)'),
           Line2D([], [], ls='', marker='s', color=TOOTH, markersize=7, label='teeth (crown and root surface)')]
if bone is not None:
    handles.append(Line2D([], [], ls='', marker='s', color=BONE, markersize=7, label='outer bone surface'))
fig.legend(handles=handles, loc='upper right', bbox_to_anchor=(0.99, 0.975), ncol=2, frameon=False, fontsize=8.5, labelcolor=INK_2)
fig.suptitle('Miniscrew positions tested, patient 8 right maxilla', x=0.05, ha='left', y=0.96, color=INK, fontsize=12.5)
fig.text(0.05, 0.02, 'Frame of the 16/17 interradicular space. The miniscrew is modelled as a point (the fixed end of the elastic); '
         'darker = further from the crest (circles) or from the paper\'s site (triangles); the paper\'s site (16/17, 4 mm) is shared by both series.', color=MUTED, fontsize=8)
fig.savefig(a.out, dpi=160, facecolor=SURFACE); plt.close(fig)
print('->', a.out)
