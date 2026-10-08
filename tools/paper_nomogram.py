"""Design chart: elastic force needed for a target second-molar distalisation, by miniscrew height.

    python tools/paper_nomogram.py      # -> out/paper/fig9_nomogram.png/.pdf

Uses the linear force law fitted in build4d/v2/force_study.json (anchorage = intercept + slope x force,
R2 >= 0.995 over 100-250 gf). Slope and intercept are interpolated linearly between the three simulated
heights (4, 6, 8 mm). Extrapolation beyond 250 gf is drawn dashed: it is not backed by a run.
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(HERE)
FS = json.load(open('build4d/v2/force_study.json'))
V2 = json.load(open('build4d/v2/compare_positions_v2.json'))
OUT = 'out/paper'
os.makedirs(OUT, exist_ok=True)

INK, MUTED, GRID, ACCENT, CLAY = '#1a1a1a', '#5a5a5a', '#d8d8d8', '#0072B2', '#D55E00'
RAMP = ['#08306b', '#2171b5', '#4292c6', '#9ecae1']
plt.rcParams.update({
    'font.size': 8, 'axes.titlesize': 8.5, 'axes.labelsize': 8, 'legend.fontsize': 7.5,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
    'text.color': INK, 'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.spines.top': False,
    'axes.spines.right': False, 'savefig.dpi': 300, 'savefig.bbox': 'tight', 'legend.frameon': False,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.5, 'axes.axisbelow': True})

H = np.array([4.0, 6.0, 8.0])
SITES = ['baseline', 'h1617_a6', 'h1617_a8']
TARGETS = [0.6, 0.8, 1.0, 1.2]
hh = np.linspace(4, 8, 81)

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=True)
for ax, (g, gl) in zip(axes, [('buccal_tad', 'Buccal miniscrew'), ('palatal_tad', 'Palatal miniscrew')]):
    d = FS['groups'][g]['d17_70']
    slope = np.array([d[s]['slope_per_gf'] for s in SITES])       # mm per gf
    icept = np.array([d[s]['intercept'] for s in SITES])          # mm at 0 gf
    sl = np.interp(hh, H, slope)
    ic = np.interp(hh, H, icept)
    for i, t in enumerate(TARGETS):
        f = (t - ic)/sl
        inside = f <= 250
        ax.plot(hh[inside], f[inside], lw=1.8, color=RAMP[i])
        if (~inside).any():
            ax.plot(hh[~inside], f[~inside], lw=1.4, ls=':', color=RAMP[i])
        j = np.argmax(f > 0)
        ax.annotate('%.1f mm' % t, (hh[j], f[j]), xytext=(3, 3), textcoords='offset points',
                    size=6.8, color=RAMP[i])
    ax.axhspan(250, 400, color=CLAY, alpha=0.06)
    ax.axhline(250, color=CLAY, lw=1, ls='--')
    ax.annotate('above 250 gf: not simulated', (4.05, 256), size=6.8, color=CLAY)
    ax.axhline(150, color=MUTED, lw=1, ls='--')
    paper = V2[g]['metrics']['baseline']['d17_70']
    ax.plot([4], [150], marker='*', ms=11, color=ACCENT, zorder=5)
    ax.annotate('source study\n4 mm, 150 gf\n(%.2f mm)' % paper, (4.05, 150), xytext=(6, -26),
                textcoords='offset points', size=6.8, color=ACCENT)
    ax.set_xlabel('miniscrew height above the alveolar crest (mm)')
    ax.set_title(gl, loc='left')
    ax.set_xticks([4, 5, 6, 7, 8]); ax.set_ylim(60, 330)
    ax.grid(axis='both')
axes[0].set_ylabel('elastic force required (gf)')
fig.text(0.5, -0.04, 'contours: second-molar distalisation achieved at step 70, from the fitted linear force law; '
         'dotted where force exceeds the simulated range', ha='center', color=MUTED, size=7)
for ext in ('png', 'pdf'):
    fig.savefig('%s/fig9_nomogram.%s' % (OUT, ext))
print('-> %s/fig9_nomogram.png' % OUT)

rows = ['| target FDI 17 @70 | buccal 4 mm | buccal 6 mm | buccal 8 mm | palatal 4 mm | palatal 6 mm | palatal 8 mm |',
        '|---|---|---|---|---|---|---|']
for t in TARGETS:
    cells = []
    for g in ('buccal_tad', 'palatal_tad'):
        d = FS['groups'][g]['d17_70']
        for s in SITES:
            f = (t - d[s]['intercept'])/d[s]['slope_per_gf']
            cells.append('%.0f gf' % f if f <= 250 else '>250 gf')
    rows.append('| %.1f mm | ' % t + ' | '.join(cells) + ' |')
open('%s/table5c_nomogram.md' % OUT, 'w').write(
    '**Table 5c. Elastic force required for a target second-molar distalisation (linear force law; '
    '">250 gf" is outside the simulated range)**\n\n' + '\n'.join(rows) + '\n')
print('-> %s/table5c_nomogram.md' % OUT)
print('\n'.join(rows))
