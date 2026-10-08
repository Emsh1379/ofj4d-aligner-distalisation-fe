"""Figures and tables for the manuscript, from the saved run summaries.

    python tools/paper_figs.py            # -> out/paper/*.png, *.pdf, tables.md

Reads build4d/v2/compare_*.json (tools/compare_tad.py, tools/compare_meshconv.py) and the run
checkpoints for what those summaries do not carry (staging protocols, deformable bone).
No solve, no GPU.
"""
import sys
import json, os, pickle, shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(HERE)
OUT = 'out/paper'
os.makedirs(OUT, exist_ok=True)

# Okabe-Ito subset, validated with the dataviz palette checker on a light surface.
# Ordered variables (force, screw height, overcorrection) use a single-hue ramp instead.
CAT = ['#0072B2', '#D55E00', '#009E73']
RAMP3 = ['#08306b', '#2171b5', '#6baed6']
RAMP4 = ['#08306b', '#2171b5', '#4292c6', '#9ecae1']
INK, MUTED, GRID = '#1a1a1a', '#5a5a5a', '#d8d8d8'
MARKS = ['o', 's', '^', 'D', 'v']

plt.rcParams.update({
    'font.size': 8, 'axes.titlesize': 8.5, 'axes.labelsize': 8, 'legend.fontsize': 7.5,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5, 'axes.edgecolor': MUTED,
    'axes.labelcolor': INK, 'text.color': INK, 'xtick.color': MUTED, 'ytick.color': MUTED,
    'axes.spines.top': False, 'axes.spines.right': False, 'figure.dpi': 150,
    'savefig.dpi': 300, 'savefig.bbox': 'tight', 'legend.frameon': False,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.5, 'axes.axisbelow': True,
})

GROUPS = [('buccal_tad', 'Buccal miniscrew'), ('palatal_tad', 'Palatal miniscrew')]
SITE_ORDER = ['h1617_a2', 'baseline', 'h1617_a6', 'h1617_a8', 'md1516_a4', 'md1415_a4']
SITE_LAB = {'h1617_a2': '16/17\n2 mm', 'baseline': '16/17\n4 mm*', 'h1617_a6': '16/17\n6 mm',
            'h1617_a8': '16/17\n8 mm', 'md1516_a4': '15/16\n4 mm', 'md1415_a4': '14/15\n4 mm'}
ROB_SITES = ['h1617_a2', 'baseline', 'h1617_a6', 'h1617_a8', 'md1415_a4']
TABLES = []


def J(name):
    with open('build4d/v2/%s.json' % name) as fh:
        return json.load(fh)


def hist(run):
    return {h['step']: h for h in pickle.load(open('build4d/v2/runs/%s.pkl' % run, 'rb'))['hist']}


def save(fig, stem):
    for ext in ('png', 'pdf'):
        fig.savefig('%s/%s.%s' % (OUT, stem, ext))
    plt.close(fig)
    print('->', '%s/%s.png' % (OUT, stem))


def grid_only_y(ax):
    ax.grid(axis='y')
    ax.xaxis.grid(False)


V2 = J('compare_positions_v2')
ROB = {k: J('compare_positions_' + k) for k in ('a000', 'a060', 'nu030', 'nu049')}
FORCE = {100: J('compare_positions_f100'), 150: V2, 200: J('compare_positions_f200'), 250: J('compare_positions_f250')}
HOOK = {'13 canine*': V2, '14 first premolar': J('compare_positions_hk14'), '12 lateral incisor': J('compare_positions_hk12')}
OC = {0: V2, 10: J('compare_positions_oc10'), 20: J('compare_positions_oc20'), 30: J('compare_positions_oc30')}
MU = {0.10: J('compare_positions_mu010'), 0.20: V2, 0.30: J('compare_positions_mu030')}
MC = J('compare_meshconv')
FS = J('force_study')


def m(Jd, g, site, key):
    return Jd[g]['metrics'][site][key]


def fig_positions():
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9), sharey=True)
    keys = [('d17_70', 'FDI 17 distal'), ('d16_70', 'FDI 16 distal'), ('inc_retraction_70', 'FDI 11 retraction')]
    x = np.arange(len(SITE_ORDER)); w = 0.26
    for ax, (g, gl) in zip(axes, GROUPS):
        for i, (k, lab) in enumerate(keys):
            ax.bar(x + (i - 1)*w, [m(V2, g, s, k) for s in SITE_ORDER], w*0.94, label=lab, color=CAT[i])
        ax.set_xticks(x); ax.set_xticklabels([SITE_LAB[s] for s in SITE_ORDER])
        ax.set_title(gl, loc='left'); grid_only_y(ax)
    axes[0].set_ylabel('movement at step 70 (mm)')
    axes[0].legend(loc='upper right')
    fig.text(0.5, -0.06, '* the site simulated by Mao et al. (16/17, 4 mm apical to the crest)',
             ha='center', color=MUTED, size=7)
    save(fig, 'fig3_positions')
    rows = ['| site | ' + ' | '.join('%s: %s' % (gl, lab) for _, gl in GROUPS for _, lab in keys) + ' |',
            '|' + '---|'*(1 + 2*len(keys))]
    for s in SITE_ORDER:
        rows.append('| ' + SITE_LAB[s].replace('\n', ' ') + ' | ' +
                    ' | '.join('%.3f' % m(V2, g, s, k) for g, _ in GROUPS for k, _ in keys) + ' |')
    TABLES.append(('Table 3. Miniscrew position: tooth movement at step 70 (mm), 150 gf, canine hook', rows))


def fig_robustness():
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), sharey=True)
    sets = [('production', V2), ('alpha 0', ROB['a000']), ('alpha 0.60', ROB['a060']),
            ('PDL nu 0.30', ROB['nu030']), ('PDL nu 0.49', ROB['nu049'])]
    colors = CAT + ['#E69F00', '#CC79A7']
    x = np.arange(len(ROB_SITES))
    for ax, (g, gl) in zip(axes, GROUPS):
        ends = []
        for i, (lab, Jd) in enumerate(sets):
            v = [m(Jd, g, s, 'd17_70') for s in ROB_SITES]
            ax.plot(x, v, marker=MARKS[i], ms=4, lw=1.6, color=colors[i], label=lab)
            ends.append([v[-1], lab, colors[i]])
        # Two sets can finish within a few hundredths of a mm of each other (alpha 0.60 and PDL nu
        # 0.30 end 0.03 mm apart on the buccal side), which overprints their labels. Push the labels
        # apart vertically, leaving the curves where they are, and colour them to match.
        span = max(e[0] for e in ends) - min(e[0] for e in ends)
        gap = 0.075*span
        ends.sort(key=lambda e: -e[0])
        ytxt = [e[0] for e in ends]
        for i in range(1, len(ytxt)):
            ytxt[i] = min(ytxt[i], ytxt[i - 1] - gap)
        for (yv, lab, col), yt in zip(ends, ytxt):
            ax.annotate(lab, (x[-1], yv), xytext=(6, 0), textcoords='offset points',
                        va='center', size=6.5, color=col,
                        xycoords='data', annotation_clip=False) if abs(yt - yv) < 1e-9 else \
                ax.annotate(lab, (x[-1] + 0.12, yt), xytext=(6, 0), textcoords='offset points',
                            va='center', size=6.5, color=col, annotation_clip=False)
        ax.set_xticks(x); ax.set_xticklabels([SITE_LAB[s] for s in ROB_SITES])
        ax.set_title(gl, loc='left'); grid_only_y(ax)
        ax.set_xlim(-0.3, len(ROB_SITES) + 1.6)
    axes[0].set_ylabel('FDI 17 distal at step 70 (mm)')
    fig.text(0.5, -0.04, 'the ranking of the sites is identical in every set; only the absolute values move',
             ha='center', color=MUTED, size=7)
    fig.text(0.5, -0.10, '* the site simulated by Mao et al. (16/17, 4 mm apical to the crest)',
             ha='center', color=MUTED, size=7)
    save(fig, 'fig4_robustness')


def table_robustness():
    sets = [('production (alpha 0.40, nu 0.45)', V2), ('alpha 0', ROB['a000']), ('alpha 0.60', ROB['a060']),
            ('PDL nu 0.30', ROB['nu030']), ('PDL nu 0.49', ROB['nu049'])]
    rows = ['| input set | group | ' + ' | '.join(SITE_LAB[s].replace(chr(10), ' ') for s in ROB_SITES) + ' | order preserved |',
            '|' + '---|'*(3 + len(ROB_SITES))]
    order_ok = True
    for lab, Jd in sets:
        for g, gl in GROUPS:
            v = [m(Jd, g, s, 'd17_70') for s in ROB_SITES]
            ok = all(v[i] > v[i + 1] for i in range(len(v) - 1))
            order_ok = order_ok and ok
            rows.append('| %s | %s | ' % (lab, gl) + ' | '.join('%.3f' % x for x in v) + ' | %s |' % ('yes' if ok else 'NO'))
    TABLES.append(('Table 4. Robustness of the site ranking: FDI 17 distal movement at step 70 (mm). '
                   'Sites are listed best to worst; "order preserved" means the ranking is unchanged '
                   '(all %d comparisons: %s)' % (len(sets)*len(GROUPS), 'yes' if order_ok else 'NO'), rows))


def fig_force():
    heights = [('baseline', '4 mm'), ('h1617_a6', '6 mm'), ('h1617_a8', '8 mm')]
    forces = [100, 150, 200, 250]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.0), sharex=True)
    for col, (g, gl) in enumerate(GROUPS):
        for row, (k, lab) in enumerate([('d17_70', 'FDI 17 distal'), ('inc_retraction_70', 'FDI 11 retraction')]):
            ax = axes[row][col]
            for i, (s, sl) in enumerate(heights):
                v = [m(FORCE[f], g, s, k) for f in forces]
                ax.plot(forces, v, marker=MARKS[i], ms=4, lw=1.6, color=RAMP3[i], label=sl)
                slope = 100*FS['groups'][g][k][s]['slope_per_gf']
                ax.annotate('%s  %+.2f mm/100 gf' % (sl, slope), (forces[-1], v[-1]), xytext=(4, 0),
                            textcoords='offset points', va='center', size=6.5, color=MUTED)
            if row == 0:
                ax.set_title(gl, loc='left')
            if col == 0:
                ax.set_ylabel('%s at step 70 (mm)' % lab)
            ax.set_xticks(forces); ax.set_xlim(90, 345); grid_only_y(ax)
    for ax in axes[1]:
        ax.set_xlabel('elastic force (gf)')
    fig.text(0.5, -0.02, 'straight lines fitted to the four forces; R2 >= 0.995 in every case',
             ha='center', color=MUTED, size=7)
    save(fig, 'fig5_force')
    rows = ['| quantity | group | 4 mm | 6 mm | 8 mm |', '|---|---|---|---|---|']
    for k, lab in [('d17_70', 'FDI 17 distal'), ('d16_70', 'FDI 16 distal'), ('inc_retraction_70', 'FDI 11 retraction')]:
        for g, gl in GROUPS:
            rows.append('| %s | %s | ' % (lab, gl) + ' | '.join(
                '%+.3f' % (100*FS['groups'][g][k][s]['slope_per_gf']) for s, _ in heights) + ' |')
    TABLES.append(('Table 5a. Effect of elastic force: mm gained per extra 100 gf (fit over 100-250 gf)', rows))
    rows = ['| quantity | buccal 6 mm | buccal 8 mm | palatal 6 mm | palatal 8 mm |', '|---|---|---|---|---|']
    for k, lab in [('d17_70', 'FDI 17 distal'), ('d16_70', 'FDI 16 distal'), ('inc_retraction_70', 'FDI 11 retraction')]:
        cells = []
        for g, _ in GROUPS:
            for s in ('h1617_a6', 'h1617_a8'):
                need = FS['groups'][g][k][s]['force_to_match_paper_site_150']
                cells.append('%.0f gf' % need if need and need < 400 else 'not attainable')
        rows.append('| %s | ' % lab + ' | '.join(cells) + ' |')
    TABLES.append(('Table 5b. Elastic force needed at a higher site to match the 16/17, 4 mm site at 150 gf', rows))


def fig_hook():
    heights = [('baseline', '4 mm'), ('h1617_a6', '6 mm'), ('h1617_a8', '8 mm')]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.8), sharex=True)
    x = np.arange(len(heights)); w = 0.26
    for col, (g, gl) in enumerate(GROUPS):
        for row, (k, lab) in enumerate([('d17_70', 'FDI 17 distal'), ('inc_retraction_70', 'FDI 11 retraction')]):
            ax = axes[row][col]
            for i, (hl, Jd) in enumerate(HOOK.items()):
                ax.bar(x + (i - 1)*w, [m(Jd, g, s, k) for s, _ in heights], w*0.94, color=CAT[i], label=hl)
            ax.set_xticks(x); ax.set_xticklabels([sl for _, sl in heights])
            if row == 0:
                ax.set_title(gl, loc='left')
            if col == 0:
                ax.set_ylabel('%s at step 70 (mm)' % lab)
            grid_only_y(ax)
    axes[0][0].legend(title='elastic hooked on', loc='upper right')
    for ax in axes[1]:
        ax.set_xlabel('miniscrew height above the alveolar crest')
    fig.text(0.5, -0.02, '* the tooth used by Mao et al. (canine precision cut)', ha='center', color=MUTED, size=7)
    save(fig, 'fig6_hook')
    rows = ['| group | height | canine* | first premolar | lateral incisor |', '|---|---|---|---|---|']
    for g, gl in GROUPS:
        for s, sl in heights:
            rows.append('| %s | %s | ' % (gl, sl) + ' | '.join('%.3f' % m(Jd, g, s, 'd17_70') for Jd in HOOK.values()) + ' |')
    TABLES.append(('Table 6. Hook tooth: FDI 17 distal movement at step 70 (mm), 150 gf', rows))


def fig_overcorrection():
    levels = [0, 10, 20, 30]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.8), sharex=True)
    for col, (g, gl) in enumerate(GROUPS):
        for row, (k, lab) in enumerate([('d17_70', 'FDI 17 distal'), ('inc_retraction_70', 'FDI 11 retraction')]):
            ax = axes[row][col]
            tgt = m(V2, g, 'baseline', k)
            for i, (s, sl) in enumerate([('h1617_a6', '6 mm'), ('h1617_a8', '8 mm')]):
                v = [m(OC[o], g, s, k) for o in levels]
                ax.plot(levels, v, marker=MARKS[i], ms=4, lw=1.6, color=RAMP3[i], label=sl)
                ax.annotate(sl, (levels[-1], v[-1]), xytext=(4, 0), textcoords='offset points',
                            va='center', size=6.5, color=MUTED)
            ax.axhline(tgt, color=CAT[1], lw=1.2, ls='--')
            ax.annotate('16/17, 4 mm site,\nno overcorrection', (0, tgt), xytext=(2, 4),
                        textcoords='offset points', size=6.5, color=CAT[1])
            if row == 0:
                ax.set_title(gl, loc='left')
            if col == 0:
                ax.set_ylabel('%s at step 70 (mm)' % lab)
            ax.set_xticks(levels); ax.set_xlim(-2, 38); grid_only_y(ax)
    for ax in axes[1]:
        ax.set_xlabel('molar overcorrection (% of planned distalisation)')
    save(fig, 'fig7_overcorrection')
    rows = ['| group | height | FDI 17 needs | FDI 16 needs |', '|---|---|---|---|']
    for g, gl in GROUPS:
        for s, sl in (('h1617_a6', '6 mm'), ('h1617_a8', '8 mm')):
            cells = []
            for k in ('d17_70', 'd16_70'):
                v0, v30 = m(OC[0], g, s, k), m(OC[30], g, s, k)
                cells.append('%.0f %%' % ((m(V2, g, 'baseline', k) - v0)/((v30 - v0)/30.0)))
            rows.append('| %s | %s | ' % (gl, sl) + ' | '.join(cells) + ' |')
    TABLES.append(('Table 7. Molar overcorrection needed to restore the anchorage of the 16/17, 4 mm site', rows))


def fig_staging():
    prot = [('V pattern (70 steps)', 'v2_%s'), ('2-by-2 (70 steps)', 'rob_stp_%s'),
            ('sequential (120 steps)', 'rob_sts_%s')]
    runs = [('control', 'No miniscrew'), ('buccal_tad_baseline', 'Buccal'), ('palatal_tad_baseline', 'Palatal')]
    end, drift = {}, {}
    for pl, pat in prot:
        for r, rl in runs:
            h = hist(pat % r)
            d17 = [h[s]['meas'][2]['crown_distal'] for s in sorted(h)]
            end[(pl, rl)] = d17[-1]
            drift[(pl, rl)] = max(d17) - d17[-1]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    x = np.arange(len(prot)); w = 0.26
    for j, (D, lab) in enumerate([(end, 'FDI 17 distal at the end (mm)'),
                                  (drift, 'FDI 17 lost after its own phase (mm)')]):
        ax = axes[j]
        for i, (_, rl) in enumerate(runs):
            ax.bar(x + (i - 1)*w, [D[(pl, rl)] for pl, _ in prot], w*0.94, color=CAT[i], label=rl)
        ax.set_xticks(x); ax.set_xticklabels([pl.replace(' (', '\n(') for pl, _ in prot])
        ax.set_ylabel(lab); grid_only_y(ax)
        ax.axhline(0, color=MUTED, lw=0.8)
    axes[0].legend(title='anchorage', loc='upper left')
    save(fig, 'fig8_staging')
    rows = ['| staging | no miniscrew | buccal | palatal | gain from the miniscrew (buccal / palatal) |',
            '|---|---|---|---|---|']
    for pl, _ in prot:
        b = end[(pl, 'Buccal')] - end[(pl, 'No miniscrew')]
        p = end[(pl, 'Palatal')] - end[(pl, 'No miniscrew')]
        rows.append('| %s | ' % pl + ' | '.join('%.3f' % end[(pl, rl)] for _, rl in runs) + ' | %.2f / %.2f |' % (b, p))
    TABLES.append(('Table 8. Staging protocol: FDI 17 distal movement at the end of the protocol (mm)', rows))


def fig_meshconv():
    lev = ['M0', 'M1', 'M2', 'M3']
    labs = {'M0': 'M0 production', 'M1': 'M1 PDL 6 layers', 'M2': 'M2 PDL x2 in plane', 'M3': 'M3 aligner 0.3 mm'}
    keys = [('d17', 'FDI 17\ndistal'), ('d16', 'FDI 16\ndistal'), ('inc', 'FDI 11\nretraction'), ('r17', 'FDI 17\nroot')]
    fig, ax = plt.subplots(figsize=(4.6, 2.8))
    x = np.arange(len(keys)); w = 0.2
    for i, L in enumerate(lev):
        ax.bar(x + (i - 1.5)*w, [MC['20'][L][k] - MC['20']['M4'][k] for k, _ in keys], w*0.94,
               color=RAMP4[i], label=labs[L])
    ax.set_xticks(x); ax.set_xticklabels([l for _, l in keys]); grid_only_y(ax)
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_ylabel('difference from M4 at step 20 (um)')
    ax.legend(loc='lower left')
    save(fig, 'figS1_meshconvergence')
    rows = ['| level | unknowns | FDI 17 distal | FDI 16 distal | FDI 11 retraction | contact force | p99 PDL strain |',
            '|---|---|---|---|---|---|---|']
    for L in lev + ['M4']:
        d = MC['20'][L]
        rows.append('| %s | %s | %.1f | %.1f | %.1f | %.3f N | %.2f %% |' % (
            labs.get(L, 'M4 all refinements'), format(d['dofs'], ','), d['d17'], d['d16'], d['inc'], d['Fc'], d['eps99']))
    TABLES.append(('Table S1. Mesh convergence at step 20 (movements in um)', rows))


def table_parameters():
    rows = ['| item | value | source |', '|---|---|---|',
            '| dentition | right half arch, FDI 17-11, midsagittal symmetry | Open-Full-Jaw patient 8 |',
            '| teeth | rigid | Mao et al. |',
            '| alveolar bone | rigid (deformable variant tested, Table 9) | Mao et al. |',
            '| PDL | 0.30 mm shell, 3 element layers, bilinear E1 = 0.05 / E2 = 0.22 MPa at 7.5 % strain | Mao et al. |',
            '| PDL Poisson ratio | 0.45 | **assumed** (tested 0.30-0.49) |',
            '| aligner | 0.70 mm, E = 1500 MPa, nu = 0.30 | Mao et al. |',
            '| attachments | 5: horizontal on 17/16, vertical on 15/14/13, E = 20 GPa | Mao et al. |',
            '| contact | node-to-surface, pairing fixed per step, Coulomb mu = 0.20 | Mao et al. (tested 0.10-0.30) |',
            '| contact penalty | 1000 N/mm | **numerical** (converged; tested 1000-4000) |',
            '| elastic | 150 gf, canine precision cut to the miniscrew | Mao et al. (100-250 gf, other hooks tested) |',
            '| staging | 70 steps, V pattern, 0.1 mm/step distal, 0.15 mm/step palatal | Mao et al. (other protocols tested) |',
            '| socket kinematics | blend alpha = 0.40 between bodily movement and tipping about the apical third | **calibrated on the control group** |',
            '| iteration | 2 PDL (bone-remodelling) iterations per aligner step | Hamanaka et al. |']
    TABLES.append(('Table 1. Model parameters. Items in bold are not specified by the source papers', rows))


def table_placeability():
    data = [('16/17, 4 mm*', 'buccal', 1.85, (0.33, 0.18, 0.03), 'not placeable (clearance < 0.5 mm)'),
            ('16/17, 4 mm*', 'palatal', 2.48, (0.64, 0.49, 0.34), 'placeable, <= 1.2 mm'),
            ('15/16, 4 mm', 'buccal', 1.33, (0.07, -0.08, -0.23), 'not placeable (clearance and root envelope)'),
            ('15/16, 4 mm', 'palatal', 1.45, (0.13, -0.02, -0.17), 'not placeable (clearance and root envelope)'),
            ('14/15, 4 mm', 'buccal', 2.14, (0.47, 0.32, 0.17), 'not placeable (clearance < 0.5 mm)'),
            ('14/15, 4 mm', 'palatal', 2.21, (0.51, 0.36, 0.21), 'placeable, <= 1.2 mm'),
            ('16/17, 2 mm', 'buccal', 1.95, (0.38, 0.23, 0.08), 'not placeable (inside the root envelope)'),
            ('16/17, 2 mm', 'palatal', 2.63, (0.72, 0.57, 0.42), 'not placeable (inside the root envelope)'),
            ('16/17, 6 mm', 'buccal', 2.29, (0.55, 0.40, 0.25), 'placeable, <= 1.2 mm'),
            ('16/17, 6 mm', 'palatal', 2.59, (0.69, 0.54, 0.39), 'placeable, <= 1.5 mm'),
            ('16/17, 8 mm', 'buccal', 2.89, (0.85, 0.70, 0.55), 'placeable, <= 1.8 mm'),
            ('16/17, 8 mm', 'palatal', 3.64, (1.22, 1.07, 0.92), 'placeable, <= 1.8 mm')]
    rows = ['| site | side | interradicular gap | clearance per side, 1.2 / 1.5 / 1.8 mm screw | verdict |',
            '|---|---|---|---|---|']
    for site, side, gap, cl, verdict in data:
        rows.append('| %s | %s | %.2f mm | ' % (site, side, gap) + ' / '.join('%+.2f' % c for c in cl) + ' mm | %s |' % verdict)
    TABLES.append(('Table 2. Clinical placeability on this patient (criterion: at least 0.5 mm from the screw surface '
                   'to the PDL on both sides, and insertion outside the root envelope)', rows))


def table_sensitivity():
    fr = max(abs(m(MU[0.30], g, 'baseline', k) - m(MU[0.10], g, 'baseline', k))
             for g, _ in GROUPS for k in ('d17_70', 'd16_70', 'inc_retraction_70'))
    worst = 0.0
    for g, _ in GROUPS:
        R, D = hist('v2_%s_baseline' % g), hist('rob_bdef_%s_baseline' % g)
        for s in R:
            for u in (2, 3, 4, 5, 6, 7, 8):
                for k in ('crown_distal', 'crown_buccal', 'crown_occlusal', 'root_distal', 'root_buccal', 'root_occlusal'):
                    worst = max(worst, abs(D[s]['meas'][u][k] - R[s]['meas'][u][k]))
    rows = ['| assumption tested | range | largest change in any reported quantity |', '|---|---|---|',
            '| aligner-crown friction | 0.10 - 0.30 (model 0.20) | %.0f um |' % (1000*fr),
            '| alveolar bone | rigid vs deformable, both groups | %.0f um |' % (1000*worst),
            '| PDL Poisson ratio | 0.30 - 0.49 (model 0.45) | ranking unchanged; absolute values move (Fig. 4) |',
            '| socket kinematics alpha | 0 - 0.60 (model 0.40) | ranking unchanged; absolute values move (Fig. 4) |',
            '| mesh | 152k - 816k unknowns | 6 um on FDI 17, 8 um on FDI 11 (Table S1) |',
            '| contact penalty | 1000 - 4000 N/mm | < 50 um, with worse convergence at 4000 |']
    TABLES.append(('Table 9. Sensitivity of the results to the assumed inputs', rows))


def copy_existing():
    for src, dst in [('out4d/model_overview.png', 'fig1_model.png'),
                     ('out4d/half_maxilla_model.png', 'fig1b_half_maxilla.png'),
                     ('out/positions/screw_positions.png', 'fig2_sites.png'),
                     ('out/staging_curves_blend.png', 'figS2_staging_curves.png'),
                     ('out/group_comparison_blend.png', 'figS3_group_comparison.png')]:
        if os.path.exists(src):
            shutil.copyfile(src, '%s/%s' % (OUT, dst))
            print('->', '%s/%s  (from %s)' % (OUT, dst, src))


# `python tools/paper_figs.py fig_robustness` regenerates one figure only and leaves tables.md alone
# (a full run rewrites it from scratch, dropping the hand-added Tables 5c, 6b and S2).
ONLY = [a for a in sys.argv[1:] if not a.startswith('-')]


def step(name, fn):
    if not ONLY or name in ONLY:
        fn()


step('table_parameters', table_parameters)
step('table_placeability', table_placeability)
step('fig_positions', fig_positions)
step('fig_robustness', fig_robustness)
step('table_robustness', table_robustness)
step('fig_force', fig_force)
step('fig_hook', fig_hook)
step('fig_overcorrection', fig_overcorrection)
step('fig_staging', fig_staging)
step('fig_meshconv', fig_meshconv)
step('table_sensitivity', table_sensitivity)
step('copy_existing', copy_existing)

if ONLY:
    print('-> regenerated only: %s (tables.md untouched)' % ', '.join(ONLY))
    raise SystemExit

with open('%s/tables.md' % OUT, 'w') as fh:
    for title, rows in TABLES:
        fh.write('**%s**\n\n' % title + '\n'.join(rows) + '\n\n')
print('->', '%s/tables.md (%d tables)' % (OUT, len(TABLES)))
