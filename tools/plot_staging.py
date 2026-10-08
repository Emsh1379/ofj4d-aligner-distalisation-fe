"""Staging curves: achieved vs prescribed tooth movement, both anchorage groups."""
import sys, os, pickle
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE,'ofj4d')); sys.path.insert(0, os.path.join(HERE,'ofj')); os.chdir(HERE)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import spec

SURFACE   = '#fcfcfb'
INK       = '#0b0b0b'
INK_2     = '#52514e'
MUTED     = '#8a8983'
GRID      = '#e4e3de'
BUCCAL    = '#2a78d6'      # categorical slot 1
CONTROL   = '#eb6834'      # slot 2
PALATAL   = '#1baf7a'      # slot 3 -- validated all-pairs; carries a contrast WARN,
                           # so every series is also labelled directly (relief rule)

DIR = sys.argv[1] if len(sys.argv) > 1 else 'build4d'       # holds run_<group>.pkl
TAG = sys.argv[2] if len(sys.argv) > 2 else ''              # output file suffix
LABEL = sys.argv[3] if len(sys.argv) > 3 else ''            # model label for the title
RUNS = [('buccal TAD',  f'{DIR}/run_buccal_tad.pkl',  BUCCAL),
        ('palatal TAD', f'{DIR}/run_palatal_tad.pkl', PALATAL),
        ('control',     f'{DIR}/run_control.pkl',     CONTROL)]
OUT = 'out'


def series(D, u):
    """Cumulative crown movement of tooth u, in its prescribed direction."""
    key, sign = ('crown_buccal', -1.0) if u in (7, 8) else ('crown_distal', 1.0)
    return np.array([sign*h['meas'][u][key] for h in D['hist']])


def prescribed(u, steps):
    a, b = spec.V_PATTERN[u]
    inc = spec.STEP_PALATAL if u in (7, 8) else spec.STEP_DISTAL
    return np.clip(steps - a + 1, 0, b - a + 1)*inc


def main():
    data = {}
    for name, path, _ in RUNS:
        data[name] = pickle.load(open(path, 'rb'))
    steps = np.array([h['step'] for h in data['buccal TAD']['hist']])

    fig, axes = plt.subplots(2, 4, figsize=(15, 7.2), facecolor=SURFACE)
    fig.subplots_adjust(left=0.05, right=0.985, top=0.85, bottom=0.09,
                        wspace=0.27, hspace=0.42)

    for k, u in enumerate(spec.HALF_ARCH_UNN):
        ax = axes[k//4][k % 4]
        ax.set_facecolor(SURFACE)
        pres = prescribed(u, steps)
        ax.plot(steps, pres, color=MUTED, lw=1.6, ls=(0, (5, 3)), zorder=2,
                label='prescribed')
        for name, _, colour in RUNS:
            ax.plot(steps, series(data[name], u), color=colour, lw=2.0,
                    zorder=3, label=name, solid_capstyle='round')
        a, b = spec.V_PATTERN[u]
        ax.axvspan(a-0.5, b+0.5, color='#000000', alpha=0.04, lw=0, zorder=1)
        direction = 'palatal' if u in (7, 8) else 'distal'
        ax.set_title(f"FDI {spec.FDI[u]}   ·   steps {a}–{b} {direction}",
                     color=INK, fontsize=11.5, fontweight='semibold', loc='left',
                     pad=8)
        # mark and label the end of this tooth's own active phase -- the point an
        # expression rate refers to, before the later relapse
        want = spec.total_prescribed(u)
        note = []
        for name, _, colour in RUNS:
            got = series(data[name], u)[b-1]
            ax.plot([b], [got], marker='o', ms=6, mfc=colour, mec=SURFACE,
                    mew=1.6, zorder=4)
            note.append(f"{100*got/want:.0f}%")
        ax.text(0.975, 0.055, "at step %d:  %s" % (b, ' / '.join(note)),
                color=INK_2, fontsize=8.5, transform=ax.transAxes, ha='right')
        ax.grid(True, color=GRID, lw=0.8, zorder=0)
        ax.set_axisbelow(True)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        for s in ('left', 'bottom'):
            ax.spines[s].set_color(GRID)
        ax.tick_params(colors=INK_2, labelsize=9, length=3, width=0.8)
        ax.set_xlim(0, 71)
        if k >= 3:
            ax.set_xlabel('staging step', color=INK_2, fontsize=9.5)
        if k % 4 == 0:
            ax.set_ylabel('crown movement (mm)', color=INK_2, fontsize=9.5)

    # last cell: the anchorage story, as the difference between the groups
    ax = axes[1][3]
    ax.set_facecolor(SURFACE)
    for name, _, colour in RUNS:
        if name == 'control':
            continue
        d = series(data[name], 8) - series(data['control'], 8)
        ax.plot(steps, d, color=colour, lw=2.0, label=name)
    ax.axhline(0, color=MUTED, lw=1.0)
    ax.set_title('FDI 11  ·  anchorage saved', color=INK,
                 fontsize=11.5, fontweight='semibold', loc='left', pad=8)
    ax.grid(True, color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK_2, labelsize=9, length=3, width=0.8)
    ax.set_xlim(0, 71)
    ax.set_xlabel('staging step', color=INK_2, fontsize=9.5)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK_2, loc='lower right')
    ax.set_ylabel('mm more palatal with the TAD', color=INK_2, fontsize=9)

    # with a model label the full title would run into the legend on the right
    fig.suptitle(f'70-step staging simulation  ·  {LABEL}' if LABEL else
                 'Maxillary whole-arch distalisation, 70-step staging simulation',
                 color=INK, fontsize=15, fontweight='semibold', x=0.05, ha='left',
                 y=0.955)
    fig.text(0.05, 0.905,
             'Crown movement against the V-pattern prescription. Shaded band = the steps in which '
             'that tooth is prescribed to move; dot = end of its phase, labelled buccal / palatal / control expression.',
             color=INK_2, fontsize=10, ha='left')
    handles = [plt.Line2D([], [], color=c, lw=2.0) for _, _, c in RUNS] + \
              [plt.Line2D([], [], color=MUTED, lw=1.6, ls=(0, (5, 3)))]
    labels = [n for n, _, _ in RUNS] + ['prescribed']
    fig.legend(handles, labels, frameon=False, fontsize=10, labelcolor=INK_2,
               loc='upper right', ncol=4, bbox_to_anchor=(0.985, 0.982))

    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, f'staging_curves{TAG}.png')
    fig.savefig(p, dpi=150, facecolor=SURFACE)
    print('wrote', os.path.abspath(p))


if __name__ == '__main__':
    main()
