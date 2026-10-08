"""Three-group comparison: expression, final position, anchorage, and the paper."""
import sys, os, pickle
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE,'ofj4d')); sys.path.insert(0, os.path.join(HERE,'ofj')); os.chdir(HERE)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import spec

SURFACE, INK, INK_2, MUTED, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#8a8983', '#e4e3de'
BUCCAL, CONTROL, PALATAL = '#2a78d6', '#eb6834', '#1baf7a'
GROUPS = [('buccal TAD',  'buccal_tad',  BUCCAL),
          ('palatal TAD', 'palatal_tad', PALATAL),
          ('control',     'control',     CONTROL)]
TEETH = spec.HALF_ARCH_UNN
OUT = 'out'


DIR = sys.argv[1] if len(sys.argv) > 1 else 'build4d'       # holds run_<group>.pkl
TAG = sys.argv[2] if len(sys.argv) > 2 else ''              # label shown top right
SUFFIX = sys.argv[3] if len(sys.argv) > 3 else ''           # output file suffix


def load():
    return {k: pickle.load(open(f'{DIR}/run_{k}.pkl', 'rb')) for _, k, _ in GROUPS}


def val(D, step, u):
    m = D['hist'][step-1]['meas'][u]
    return -m['crown_buccal'] if u in (7, 8) else m['crown_distal']


def style(ax, ylab=None, xlab=None):
    ax.set_facecolor(SURFACE)
    ax.grid(True, axis='y', color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK_2, labelsize=9, length=3, width=0.8)
    if ylab: ax.set_ylabel(ylab, color=INK_2, fontsize=9.5)
    if xlab: ax.set_xlabel(xlab, color=INK_2, fontsize=9.5)


def title(ax, t, sub=None):
    ax.set_title(t, color=INK, fontsize=12, fontweight='semibold', loc='left', pad=14)
    if sub:
        ax.text(0.0, 1.015, sub, color=INK_2, fontsize=9,
                transform=ax.transAxes, va='bottom')


def main():
    D = load()
    fig = plt.figure(figsize=(15.5, 9.6), facecolor=SURFACE)
    gs = fig.add_gridspec(2, 2, left=0.055, right=0.985, top=0.845, bottom=0.07,
                          hspace=0.42, wspace=0.20)
    x = np.arange(len(TEETH))
    w = 0.26

    # -- A: expression at the end of each tooth's own active phase ---------------
    ax = fig.add_subplot(gs[0, 0]); style(ax, 'expression (%)')
    title(ax, 'Expression at the end of each tooth’s phase',
          'achieved / prescribed, before later relapse')
    vmax = 100.0
    for i, (name, key, col) in enumerate(GROUPS):
        v = []
        for u in TEETH:
            b = spec.V_PATTERN[u][1]
            v.append(100*val(D[key], b, u)/spec.total_prescribed(u))
        vmax = max(vmax, max(v))
        bars = ax.bar(x + (i-1)*w, v, w*0.92, color=col, label=name, zorder=3)
        for xi, vi in zip(x + (i-1)*w, v):
            ax.text(xi, vi + (2 if vi >= 0 else -6), f"{vi:.0f}", ha='center',
                    color=INK_2, fontsize=7.5, zorder=4)
    ax.axhline(0, color=MUTED, lw=1.0, zorder=2)
    ax.set_xticks(x); ax.set_xticklabels([f"FDI {spec.FDI[u]}" for u in TEETH])
    # headroom for the value labels and the legend when a tooth overshoots 100%
    ax.set_ylim(-20, 100 if vmax <= 100 else vmax + 25)
    ax.legend(frameon=False, fontsize=9.5, labelcolor=INK_2, loc='upper right', ncol=3)

    # -- B: where every tooth ends up -------------------------------------------
    ax = fig.add_subplot(gs[0, 1]); style(ax, 'crown movement (mm)')
    title(ax, 'Net movement retained at step 70',
          'positive = distal (incisors: palatal); negative = lost to anchorage')
    for i, (name, key, col) in enumerate(GROUPS):
        v = [val(D[key], 70, u) for u in TEETH]
        ax.bar(x + (i-1)*w, v, w*0.92, color=col, label=name, zorder=3)
    ax.axhline(0, color=INK, lw=1.2, zorder=4)
    ax.set_xticks(x); ax.set_xticklabels([f"FDI {spec.FDI[u]}" for u in TEETH])
    post = [u for u in TEETH if u not in (7, 8)]
    lost = [spec.FDI[u] for u in post if val(D['control'], 70, u) < 0]
    note = ('control is net mesial on every posterior tooth' if len(lost) == len(post)
            else f"control net mesial on FDI {', '.join(map(str, lost))}" if lost
            else 'every posterior tooth retains net distal movement')
    ax.text(0.015, 0.97, note, transform=ax.transAxes, ha='left', va='top',
            color=INK_2, fontsize=9)

    # -- C: anchorage over the protocol -----------------------------------------
    ax = fig.add_subplot(gs[1, 0]); style(ax, 'FDI 11 movement (mm)', 'staging step')
    title(ax, 'Anterior anchorage through the protocol',
          'central incisor; up = palatal (retraction), down = labial (anchorage loss)')
    steps = np.array([h['step'] for h in D['control']['hist']])
    for name, key, col in GROUPS:
        ax.plot(steps, [val(D[key], s, 8) for s in steps], color=col, lw=2.2,
                label=name, zorder=3, solid_capstyle='round')
        ax.plot([70], [val(D[key], 70, 8)], 'o', ms=6, mfc=col, mec=SURFACE,
                mew=1.6, zorder=4)
        ax.annotate(f" {val(D[key],70,8):.2f}", (70, val(D[key], 70, 8)),
                    color=INK_2, fontsize=9, va='center')
    ax.axhline(0, color=MUTED, lw=1.0)
    ax.axvspan(50.5, 70.5, color='#000000', alpha=0.04, lw=0, zorder=1)
    ax.text(60.5, ax.get_ylim()[0], ' incisors\n prescribed', color=INK_2,
            fontsize=8.5, va='bottom')
    ax.set_xlim(0, 76)
    ax.legend(frameon=False, fontsize=9.5, labelcolor=INK_2, loc='upper left')

    # -- D: against the paper ----------------------------------------------------
    ax = fig.add_subplot(gs[1, 1]); style(ax, 'mm')
    title(ax, 'Against Mao et al. 2024',
          'filled = this model, hatched = the paper, where it reports a value')
    metrics = [
        ('FDI 17\nstep 10', {'buccal_tad': 1.11, 'palatal_tad': 1.11, 'control': 0.97},
         lambda k: val(D[k], 10, 2)),
        ('FDI 17 drift\nsteps 20–40', {'control': 0.89},
         lambda k: val(D[k], 20, 2)-val(D[k], 40, 2)),
        ('FDI 17\nstep 70', {'palatal_tad': 1.31, 'control': 0.73},
         lambda k: val(D[k], 70, 2)),
        ('FDI 11 retraction\nstep 70', {'buccal_tad': 1.71, 'palatal_tad': 1.88,
                                        'control': 0.87},
         lambda k: val(D[k], 70, 8)),
    ]
    xm = np.arange(len(metrics)); wm = 0.13
    for i, (name, key, col) in enumerate(GROUPS):
        ours = [f(key) for _, _, f in metrics]
        paper = [p.get(key) for _, p, _ in metrics]
        off = (i-1)*2*wm
        ax.bar(xm + off - wm*0.5, ours, wm*0.95, color=col, zorder=3,
               label=name if i == 0 else None)
        pv = [p if p is not None else 0.0 for p in paper]
        ax.bar(xm + off + wm*0.5, pv, wm*0.95, facecolor='none', edgecolor=col,
               lw=1.3, hatch='////', zorder=3)
    ax.axhline(0, color=INK, lw=1.2, zorder=4)
    ax.set_xticks(xm); ax.set_xticklabels([m[0] for m in metrics], fontsize=9)
    order = sorted((k for _, k, _ in GROUPS), key=lambda k: -val(D[k], 70, 8))
    same = order == ['palatal_tad', 'buccal_tad', 'control']
    ax.text(0.015, 0.97, 'incisor retraction at step 70: '
            + ('palatal > buccal > control, as in the paper' if same
               else ' > '.join(order) + ' (paper: palatal > buccal > control)'),
            transform=ax.transAxes, ha='left', va='top', color=INK_2, fontsize=9)
    ax.set_ylim(top=2.25)

    fig.suptitle('Maxillary whole-arch distalisation — three anchorage groups',
                 color=INK, fontsize=16, fontweight='semibold', x=0.055, ha='left',
                 y=0.965)
    fig.text(0.055, 0.905,
             '70-step staging simulation, Open-Full-Jaw patient 8 maxilla, right half arch. '
             'Reproduction of Mao et al., Bioengineering 2024;11(1):3.',
             color=INK_2, fontsize=10.5, ha='left')
    os.makedirs(OUT, exist_ok=True)
    if TAG:
        fig.text(0.985, 0.905, TAG, color='#b4472f', fontsize=10.5, ha='right',
                 fontweight='semibold')
    p = os.path.join(OUT, f'group_comparison{SUFFIX}.png')
    fig.savefig(p, dpi=150, facecolor=SURFACE)
    print('wrote', os.path.abspath(p))


if __name__ == '__main__':
    main()
