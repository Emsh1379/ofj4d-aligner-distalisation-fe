"""Picture of the TAD placement: occlusal view, coronal section, line of action."""
import sys, os
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, 'tools'))
sys.path.insert(0, os.path.join(HERE, 'ofj4d')); sys.path.insert(0, os.path.join(HERE, 'ofj'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tad_check import compute
import spec

SURFACE, INK, INK_2, MUTED, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#8a8983', '#e4e3de'
COL = {'buccal': '#2a78d6', 'palatal': '#1baf7a'}   # validated pair; aqua is
                                                    # also labelled directly


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK_2, labelsize=9, length=3, width=0.8)
    ax.grid(True, color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)


BONE_CACHE = os.environ.get('OFJ_MSH_CACHE',
                            os.path.join(HERE, 'data', 'cache', 'p8_maxilla.npz'))


def bone_outer_points():
    """Outer-surface nodes of the segmented bone (Open-Full-Jaw), if available."""
    try:
        from mesh_io import boundary_faces
        from anatomy import TAG_BONE
        d = np.load(BONE_CACHE)
        nodes, tets, tags = d['nodes'], d['tets'], d['tags']
    except Exception:
        return None
    f, own = boundary_faces(tets, return_owner=True)
    return nodes[np.unique(f[tags[own] == TAG_BONE])]


def main():
    R = compute()
    os.chdir(HERE)
    S, st, fr = R['S'], R['st'], R['frame']
    mid, es, ed, a, p, l, eo = (fr[k] for k in ('mid', 'es', 'ed', 'a', 'p', 'l', 'eo'))
    B, P = R['sides']['buccal'], R['sides']['palatal']
    bone_pts = bone_outer_points()

    fig = plt.figure(figsize=(16, 7.4), facecolor=SURFACE)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.0, 0.85], left=0.045,
                          right=0.985, top=0.80, bottom=0.10, wspace=0.28)

    # ---- occlusal view ---------------------------------------------------------
    ax = fig.add_subplot(gs[0]); style(ax)
    def occ(X):
        X = np.atleast_2d(X) - mid
        return X @ l, -(X @ p)             # lateral right, anterior up
    for u in S.keep:
        T = S.tooth[u]
        cr = T['surf_rest'][:T['surf_n']][::6]
        x, y = occ(cr)
        ax.scatter(x, y, s=0.6, c='#cfcdc4', lw=0, zorder=1)
    h = B['tad'] @ eo
    win = np.vstack([S.tooth[u]['surf_rest'][:S.tooth[u]['surf_n']] for u in S.keep]
                    + [B['tad'][None], P['tad'][None], B['hook'][None], P['hook'][None]])
    wx, wy = occ(win)
    xl = (wx.min()-4, wx.max()+4); yl = (wy.min()-4, wy.max()+4)
    if bone_pts is not None:
        k = np.abs(bone_pts @ eo - h) < 0.35
        x, y = occ(bone_pts[k])
        m = (x > xl[0]) & (x < xl[1]) & (y > yl[0]) & (y < yl[1])
        ax.scatter(x[m], y[m], s=3.0, c='#b8955e', lw=0, zorder=3)
    for u in (2, 3):
        pd = S.tooth[u]['pdl_rest']
        k = np.abs(pd @ eo - h) < 0.6
        x, y = occ(pd[k])
        ax.scatter(x, y, s=1.4, c='#8f8a80', lw=0, zorder=2)
    for side, s in (('buccal', B), ('palatal', P)):
        c = COL[side]
        tx, ty = occ(s['tad']); hx, hy = occ(s['hook']); cx, cy = occ(s['crest'])
        ax.plot([hx[0], tx[0]], [hy[0], ty[0]], color=c, lw=2.0, zorder=3)
        ax.scatter(tx, ty, s=90, c=c, edgecolors=SURFACE, linewidths=1.5, zorder=5)
        ax.scatter(hx, hy, s=60, marker='D', c=c, edgecolors=SURFACE,
                   linewidths=1.2, zorder=5)
        ax.scatter(cx, cy, s=40, facecolors='none', edgecolors=c, linewidths=1.4,
                   zorder=4)
        ax.annotate(f"{side} TAD", (tx[0], ty[0]), xytext=(8 if side == 'buccal' else -8, 6),
                    textcoords='offset points', color=INK, fontsize=9.5,
                    ha='left' if side == 'buccal' else 'right')
        ax.annotate(f"{side} hook", (hx[0], hy[0]), xytext=(8 if side == 'buccal' else -8, -2),
                    textcoords='offset points', color=INK_2, fontsize=8.5,
                    ha='left' if side == 'buccal' else 'right')
    for u, lab in ((2, '17'), (3, '16'), (6, '13'), (8, '11')):
        x, y = occ(S.tooth[u]['crown_c'])
        ax.text(x[0], y[0], lab, color=MUTED, fontsize=8.5, ha='center', va='center')
    ax.set_xlim(*xl); ax.set_ylim(*yl)
    ax.set_aspect('equal'); ax.set_anchor('N')
    ax.set_xlabel('buccal  →  (mm)', color=INK_2, fontsize=9.5)
    ax.set_ylabel('anterior  →  (mm)', color=INK_2, fontsize=9.5)
    ax.set_title('Occlusal view', color=INK, fontsize=12, fontweight='semibold',
                 loc='left', pad=18)
    ax.text(0, 1.012, 'crowns light  ·  bone surface tan  ·  roots dark, at screw height',
            transform=ax.transAxes, color=INK_2, fontsize=8.5, va='bottom')
    from matplotlib.lines import Line2D
    ax.legend(handles=[
        Line2D([], [], marker='o', ls='none', mfc=INK_2, mec=SURFACE, ms=8, label='TAD (screw)'),
        Line2D([], [], marker='o', ls='none', mfc='none', mec=INK_2, ms=7, label='crest height (4 mm occlusal)'),
        Line2D([], [], marker='D', ls='none', mfc=INK_2, mec=SURFACE, ms=7, label='precision cut'),
        Line2D([], [], marker='s', ls='none', mfc='#b8955e', mec='#b8955e', ms=6, label='bone surface')],
        frameon=False, fontsize=8.5, labelcolor=INK_2, loc='lower left', bbox_to_anchor=(1.02, 0.0))

    # ---- coronal section through the septum ------------------------------------
    ax = fig.add_subplot(gs[1]); style(ax)
    def cor(X):
        X = np.atleast_2d(X) - mid
        return X @ es, X @ a                # buccal right, apical up
    dsep = 0.5*((B['gap_lo'] + B['gap_hi'])) - mid @ ed
    if bone_pts is not None:
        k = np.abs((bone_pts - mid) @ ed - dsep) < 0.35
        x, y = cor(bone_pts[k])
        m = (np.abs(x) < 16) & (y > -4) & (y < 16)
        ax.scatter(x[m], y[m], s=3.0, c='#b8955e', lw=0, zorder=3)
    for u in (2, 3):
        pd = S.tooth[u]['pdl_rest']
        k = np.abs((pd - mid) @ ed - dsep) < 1.0
        x, y = cor(pd[k])
        ax.scatter(x, y, s=1.2, c='#8f8a80', lw=0, zorder=2)
    for side, s in (('buccal', B), ('palatal', P)):
        c = COL[side]
        tx, ty = cor(s['tad']); cx, cy = cor(s['crest'])
        ax.annotate('', xy=(tx[0], ty[0]), xytext=(cx[0], cy[0]),
                    arrowprops=dict(arrowstyle='->', color=c, lw=1.6), zorder=4)
        ax.scatter(cx, cy, s=45, facecolors='none', edgecolors=c, linewidths=1.5,
                   zorder=5)
        ax.scatter(tx, ty, s=90, c=c, edgecolors=SURFACE, linewidths=1.5, zorder=5)
        ax.annotate(f"{side}\n{s['apical_of_crest']:.1f} mm above crest",
                    (tx[0], ty[0]), xytext=(10 if side == 'buccal' else -10, 0),
                    textcoords='offset points', color=INK, fontsize=9,
                    ha='left' if side == 'buccal' else 'right', va='center')
        ax.axvline(s['outer']*(1 if side == 'buccal' else -1), color=c, lw=1.0,
                   ls=(0, (4, 3)), zorder=3)
    ax.set_xlim(-16, 16); ax.set_ylim(-4, 16)
    ax.set_aspect('equal'); ax.set_anchor('N')
    ax.set_xlabel('palatal  ←   →  buccal  (mm from 16/17 midline)',
                  color=INK_2, fontsize=9.5)
    ax.set_ylabel('apical  →  (mm)', color=INK_2, fontsize=9.5)
    ax.set_title('Coronal section at the 16/17 septum', color=INK, fontsize=12,
                 fontweight='semibold', loc='left', pad=18)
    ax.text(0, 1.012, 'tan = bone surface  ·  dashed = outermost root/PDL on that side',
            transform=ax.transAxes, color=INK_2, fontsize=8.5, va='bottom')

    # ---- line of action --------------------------------------------------------
    ax = fig.add_subplot(gs[2]); style(ax)
    labels = ['sagittal\n(distal)', 'lateral\n(+ buccal)', 'apical\n(intrusive)']
    x = np.arange(3); w = 0.36
    for i, (side, s) in enumerate((('buccal', B), ('palatal', P))):
        v = np.array(s['comp'])*spec.ELASTIC_FORCE_N
        ax.bar(x + (i-0.5)*w, v, w*0.92, color=COL[side], label=f"{side} TAD", zorder=3)
        for xi, vi in zip(x + (i-0.5)*w, v):
            ax.text(xi, vi + (0.03 if vi >= 0 else -0.09), f"{vi:.2f}", ha='center',
                    color=INK_2, fontsize=8.5)
    ax.axhline(0, color=INK, lw=1.0, zorder=4)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel('force on the precision cut (N)', color=INK_2, fontsize=9.5)
    ax.legend(frameon=False, fontsize=9.5, labelcolor=INK_2, loc='upper right')
    ax.set_title('Elastic line of action', color=INK, fontsize=12,
                 fontweight='semibold', loc='left', pad=18)
    ax.text(0, 1.012, f'{spec.ELASTIC_FORCE_GF:.0f} gf = {spec.ELASTIC_FORCE_N:.2f} N',
            transform=ax.transAxes, color=INK_2, fontsize=8.5, va='bottom')

    fig.suptitle('Miniscrew (TAD) placement — buccal and palatal groups',
                 color=INK, fontsize=15, fontweight='semibold', x=0.045, ha='left', y=0.965)
    fig.text(0.045, 0.905, 'Specification (Mao et al. 2024): interradicular space between '
             'the first and second molars, 4 mm above the alveolar crest; elastic to a '
             'precision cut at the canine mesial cervical region.',
             color=INK_2, fontsize=10, ha='left')
    os.makedirs('out', exist_ok=True)
    fig.savefig('out/tad_placement.png', dpi=150, facecolor=SURFACE)
    print('wrote', os.path.abspath('out/tad_placement.png'))


if __name__ == '__main__':
    main()
