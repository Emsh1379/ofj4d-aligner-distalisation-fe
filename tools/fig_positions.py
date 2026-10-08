"""Figures for the miniscrew-position study: staging curves and anchorage loss per position.

usage: python tools/fig_positions.py [--set old|v2] [--out out/positions]

Small multiples, one PNG per series (height series at 16/17: 2/4/6/8 mm; site series at 4 mm:
16/17, 15/16, 14/15): rows = quantity, columns = buccal / palatal TAD group. Within a panel the positions are
an ORDERED series and wear one validated blue ordinal ramp, light -> dark with distance from the crest or
from the paper's site; the paper's site is the reference line in every panel. Each line is labelled at its end,
and a legend names every series. A thin grey line is the prescribed movement.
  staging_<series>.png   FDI 17 and FDI 16 distal movement, incisor (FDI 11) retraction, vs step
  anchorage_<series>.png FDI 17 and FDI 16 anchorage loss vs step: loss of distal movement already
                         achieved (running maximum minus current), the paper's mesial relapse
"""
import sys, os, pickle, argparse
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, 'ofj4d')); os.chdir(HERE)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import spec

ap = argparse.ArgumentParser()
ap.add_argument('--set', default='v2', choices=['old', 'v2'])
ap.add_argument('--out', default='out/positions')
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)

SURFACE, INK, INK_2, MUTED, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#8a8983', '#e4e3de'
RAMP4 = ['#86b6ef', '#3987e5', '#1c5cab', '#0d366b']        # validated ordinal (light, 4 steps)
RAMP3 = ['#86b6ef', '#2a78d6', '#104281']                   # validated ordinal (light, 3 steps)
SERIES = {
    'height': dict(title='Miniscrew height at the 16/17 site', variants=[('h1617_a2', '2 mm'), ('baseline', '4 mm (paper)'),
                   ('h1617_a6', '6 mm'), ('h1617_a8', '8 mm')], colors=RAMP4),
    'site': dict(title='Miniscrew site at 4 mm above the crest', variants=[('baseline', '16/17 (paper)'), ('md1516_a4', '15/16'),
                 ('md1415_a4', '14/15')], colors=RAMP3),
}
GROUPS = [('buccal_tad', 'Buccal TAD'), ('palatal_tad', 'Palatal TAD')]


def path(g, v):
    if a.set == 'v2':
        return f'build4d/v2/runs/v2_{g}_{v}.pkl'
    return (f'build4d/final/run_{g}.pkl' if v == 'baseline' else f'build4d/tad_positions/run_{g}_{v}_rigid.pkl')


def load(g, v):
    p = path(g, v)
    if not os.path.exists(p):
        return None
    D = pickle.load(open(p, 'rb'))
    steps = np.array([h['step'] for h in D['hist']])
    m = lambda u, k: np.array([h['meas'][u][k] for h in D['hist']])
    return dict(steps=steps, d17=m(2, 'crown_distal'), d16=m(3, 'crown_distal'), inc=-m(8, 'crown_buccal'))


def prescribed(u, steps):
    a_, b_ = spec.V_PATTERN[u]
    inc = spec.STEP_PALATAL if u in (7, 8) else spec.STEP_DISTAL
    return np.clip(steps - a_ + 1, 0, b_ - a_ + 1)*inc


def style(ax):
    ax.set_facecolor(SURFACE)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'):
        ax.spines[sp].set_color(GRID); ax.spines[sp].set_linewidth(0.8)
    ax.tick_params(colors=INK_2, labelsize=8.5, length=3, width=0.8)
    ax.grid(True, color=GRID, lw=0.6, zorder=0); ax.set_axisbelow(True)
    ax.set_xlim(0, 72)


def end_labels(ax, items):
    """Direct labels at line ends, nudged apart so they never overlap."""
    if not items:
        return
    lo, hi = ax.get_ylim(); gap = 0.055*(hi - lo)
    items = sorted(items, key=lambda t: t[0])
    ys = [t[0] for t in items]
    for i in range(1, len(ys)):
        ys[i] = max(ys[i], ys[i-1] + gap)
    for (y0, text), y in zip(items, ys):
        ax.annotate(text, xy=(70, y0), xytext=(71.5, y), fontsize=7.5, color=INK_2, va='center', ha='left',
                    annotation_clip=False)


def figure(kind, series):
    cfg = SERIES[series]
    rows = ([('d17', 2, 'FDI 17 distal (mm)'), ('d16', 3, 'FDI 16 distal (mm)'), ('inc', 8, 'Incisor retraction (mm)')]
            if kind == 'staging' else [('d17', 2, 'FDI 17 anchorage loss (mm)'), ('d16', 3, 'FDI 16 anchorage loss (mm)')])
    H = 3.0*len(rows) + 1.6
    fig, axes = plt.subplots(len(rows), 2, figsize=(11.5, H), facecolor=SURFACE, sharex=True)
    # fixed inch margins: title 0.45 in, legend row 0.35 in, panel titles 0.35 in; x-axis band + footnote 0.8 in
    fig.subplots_adjust(left=0.07, right=0.88, top=1 - 1.15/H, bottom=0.8/H, wspace=0.32, hspace=0.30)
    missing = []
    for c, (g, gname) in enumerate(GROUPS):
        data = {v: load(g, v) for v, _ in cfg['variants']}
        missing += [f'{g}/{v}' for v, d in data.items() if d is None]
        for r, (key, u, ylab) in enumerate(rows):
            ax = axes[r][c]; style(ax)
            labels = []
            if kind == 'staging':
                st = np.arange(0, 71)
                ax.plot(st, prescribed(u, st), color=MUTED, lw=1.0, zorder=2, label='prescribed')
            for (v, vname), col in zip(cfg['variants'], cfg['colors']):
                d = data[v]
                if d is None:
                    continue
                y = d[key] if kind == 'staging' else np.maximum.accumulate(d[key]) - d[key]
                ax.plot(d['steps'], y, color=col, lw=2.0, zorder=3, label=vname, solid_capstyle='round')
                if d['steps'][-1] >= 70:
                    labels.append((float(y[-1]), vname))
            ax.set_ylabel(ylab, color=INK_2, fontsize=9)
            if r == 0:
                ax.set_title(gname, color=INK, fontsize=11, loc='left', pad=8)
            if r == len(rows) - 1:
                ax.set_xlabel('Aligner step', color=INK_2, fontsize=9)
            end_labels(ax, labels)
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc='upper left', bbox_to_anchor=(0.062, 1 - 0.42/H), ncol=len(l), frameon=False, fontsize=8.5, labelcolor=INK_2)
    what = 'Staging curves' if kind == 'staging' else 'Anchorage loss (distal movement lost since its maximum)'
    fig.suptitle(f"{what}: {cfg['title'].lower()}", x=0.07, ha='left', y=1 - 0.12/H, va='top', color=INK, fontsize=12.5)
    note = f"Final model, rigid bone, {'repaired meshes, K_N 1000' if a.set == 'v2' else 'shipped meshes, K_N 2000'}; miniscrew modelled as a point."
    if missing:
        note += f"  Missing runs: {', '.join(missing)}"
    fig.text(0.07, 0.12/H, note, color=MUTED, fontsize=8)
    out = f"{a.out}/{kind}_{series}_{a.set}.png"
    fig.savefig(out, dpi=160, facecolor=SURFACE); plt.close(fig)
    print('->', out)


for series in ('height', 'site'):
    for kind in ('staging', 'anchorage'):
        figure(kind, series)
