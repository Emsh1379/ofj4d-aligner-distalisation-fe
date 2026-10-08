"""Manuscript tables and a data appendix for the other-patient replication (FINDINGS 47).

usage: python tools/replication_tables.py [--out out/paper/replication]

Reads each patient's `build4d/v2/compare_positions_{v2,hk14}.json` (written by tools/compare_tad.py)
and `build4d/v2/runs/v2_control.pkl`, and writes markdown tables plus a long-format CSV appendix.
Patient 8 is this package; the others are ../ofj4d_patients/PN.
"""
import argparse
import csv
import json
import os
import pickle
import subprocess
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAT = os.path.normpath(os.path.join(HERE, '..', 'ofj4d_patients'))
PATIENTS = ['P8', 'P3', 'P4', 'P5', 'P6']
SITES = ['h1617_a2', 'baseline', 'h1617_a6', 'h1617_a8', 'md1516_a4', 'md1415_a4']
LABEL = {'h1617_a2': '16/17, 2 mm', 'baseline': '16/17, 4 mm (paper)', 'h1617_a6': '16/17, 6 mm',
         'h1617_a8': '16/17, 8 mm', 'md1516_a4': '15/16, 4 mm', 'md1415_a4': '14/15, 4 mm'}
RANK = ['h1617_a2', 'baseline', 'h1617_a6', 'h1617_a8', 'md1415_a4']   # 15/16 reported, not ranked
SHORT = {'h1617_a2': '2 mm', 'baseline': 'paper', 'h1617_a6': '6 mm', 'h1617_a8': '8 mm',
         'md1516_a4': '15/16', 'md1415_a4': '14/15'}
GROUPS = ['buccal_tad', 'palatal_tad']
Q = [('d17_70', 'FDI 17 distal @70'), ('inc_retraction_70', 'incisor retraction @70')]


def ws(p):
    return HERE if p == 'P8' else os.path.join(PAT, p)


def summaries(p, tag):
    """compare_tad output for one patient and one job set, regenerating it when missing."""
    f = os.path.join(ws(p), 'build4d/v2', f'compare_positions_{tag}.json')
    if not os.path.exists(f):
        subprocess.run([sys.executable, 'tools/compare_tad.py', '--set', tag, '--json', f],
                       cwd=ws(p), check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if not os.path.exists(f):
        return {}
    d = json.load(open(f))
    return {g: {v: m for v, m in gd['metrics'].items() if m.get('steps', 0) >= 70}
            for g, gd in d.items()}


def control(p):
    f = os.path.join(ws(p), 'build4d/v2/runs/v2_control.pkl')
    if not os.path.exists(f):
        return None
    h = {x['step']: x for x in pickle.load(open(f, 'rb'))['hist']}
    if 70 not in h:
        return None
    return dict(d17_70=h[70]['meas'][2]['crown_distal'],
                inc_retraction_70=-h[70]['meas'][8]['crown_buccal'])


def fmt(x, w=7, sign=False):
    return ('-'.rjust(w) if x is None else f"{x:{'+' if sign else ''}{w}.3f}")


ap = argparse.ArgumentParser()
ap.add_argument('--out', default='out/paper/replication')
a = ap.parse_args()
out = os.path.join(HERE, a.out)
os.makedirs(out, exist_ok=True)

anat = json.load(open(os.path.join(PAT, 'anatomy.json')))
V = {p: summaries(p, 'v2') for p in PATIENTS}
H = {p: summaries(p, 'hk14') for p in PATIENTS}
C = {p: control(p) for p in PATIENTS}
rows, md = [], []

md.append('**Table R1. The five maxillae.** Patient 8 is the source study\'s patient; the other four are '
          'the replication set. Root area is the PDL-covered root surface.\n')
md.append('| patient | 16-26 width (mm) | half-arch length (mm) | root area 17 / 16 / 13 / 11 (mm2) |')
md.append('|---|---|---|---|')
for p in PATIENTS:
    ra = anat[p]['root_area']
    md.append(f"| {p[1:]} | {anat[p]['width_mb16_26']:.1f} | {anat[p]['half_arch_length']:.1f} | "
              f"{ra['17']:.0f} / {ra['16']:.0f} / {ra['13']:.0f} / {ra['11']:.0f} |")

for key, title in Q:
    md.append(f'\n**Table R2{"ab"[[k for k, _ in Q].index(key)]}. {title} (mm) by patient, group and '
              f'miniscrew site.** Control is the no-miniscrew run.\n')
    md.append('| patient | group | control | ' + ' | '.join(LABEL[s] for s in SITES) + ' |')
    md.append('|---' * (len(SITES) + 3) + '|')
    for p in PATIENTS:
        for g in GROUPS:
            c = (C[p] or {}).get(key)
            cells = []
            for s in SITES:
                m = V[p].get(g, {}).get(s)
                cells.append('-' if m is None else f"{m[key]:.3f}")
            md.append(f"| {p[1:]} | {g.replace('_tad', '')} | "
                      f"{'-' if c is None else f'{c:.3f}'} | " + ' | '.join(cells) + ' |')

md.append('\n**Table R3. Site ranking per patient and group.** Ranked best to worst on FDI 17 distal '
          'movement at step 70; the 15/16 site is reported but not ranked, as in the single-patient study.\n')
md.append('| patient | group | ranking | same as patient 8 |')
md.append('|---|---|---|---|')
ref = None
for p in PATIENTS:
    for g in GROUPS:
        got = {s: V[p].get(g, {}).get(s) for s in RANK}
        if any(v is None for v in got.values()):
            md.append(f"| {p[1:]} | {g.replace('_tad', '')} | (incomplete) | - |")
            continue
        order = sorted(RANK, key=lambda s: -got[s]['d17_70'])
        txt = ' > '.join(SHORT[s] for s in order)
        if p == 'P8' and g == 'buccal_tad':
            ref = order
        md.append(f"| {p[1:]} | {g.replace('_tad', '')} | {txt} | "
                  f"{'yes' if order == ref else 'NO'} |")

md.append('\n**Table R4. Palatal minus buccal at each miniscrew site (mm).** Positive favours the palatal '
          'site. The difference rises with screw height in every patient.\n')
md.append('| patient | ' + ' | '.join(LABEL[s] for s in SITES) + ' | mean |')
md.append('|---' * (len(SITES) + 2) + '|')
for p in PATIENTS:
    cells, vals = [], []
    for s in SITES:
        b = V[p].get('buccal_tad', {}).get(s)
        q = V[p].get('palatal_tad', {}).get(s)
        if b is None or q is None:
            cells.append('-')
            continue
        d = q['d17_70'] - b['d17_70']
        vals.append(d)
        cells.append(f"{d:+.3f}")
    m = f"{sum(vals)/len(vals):+.3f}" if vals else '-'
    md.append(f"| {p[1:]} | " + ' | '.join(cells) + f" | {m} |")

md.append('\n**Table R5. First-premolar minus canine precision cut (mm), step 70.** Positive favours the '
          'first-premolar cut.\n')
md.append('| patient | group | FDI 17, paper site | incisor, paper site | FDI 17, 6 mm | incisor, 6 mm |')
md.append('|---|---|---|---|---|---|')
for p in PATIENTS:
    for g in GROUPS:
        cs = []
        for s in ('baseline', 'h1617_a6'):
            for k, _ in Q:
                h = H[p].get(g, {}).get(s)
                v = V[p].get(g, {}).get(s)
                cs.append('-' if (h is None or v is None) else f"{h[k] - v[k]:+.3f}")
        md.append(f"| {p[1:]} | {g.replace('_tad', '')} | {cs[0]} | {cs[1]} | {cs[2]} | {cs[3]} |")

md.append('\n**Table R6. Miniscrew benefit over the control at the paper site (mm), step 70.**\n')
md.append('| patient | buccal, FDI 17 | palatal, FDI 17 | buccal, incisor | palatal, incisor |')
md.append('|---|---|---|---|---|')
for p in PATIENTS:
    c = C[p] or {}
    cs = []
    for k, _ in Q:
        for g in GROUPS:
            m = V[p].get(g, {}).get('baseline')
            cs.append('-' if (m is None or k not in c) else f"{m[k] - c[k]:+.3f}")
    md.append(f"| {p[1:]} | {cs[0]} | {cs[1]} | {cs[2]} | {cs[3]} |")

md.append('\n**Table R7. Molar anchorage does not explain the palatal-buccal split.** The PDL secant '
          'modulus of FDI 17 and 16 was scaled at the paper screw site, holding the geometry fixed '
          '(`OFJ_PDL_E_SCALE`); the 1.0 row is the production run. Palatal minus buccal, step 70 (mm).\n')
md.append('| patient | molar PDL scale | buccal FDI 17 | palatal FDI 17 | palatal - buccal, FDI 17 | palatal - buccal, incisor |')
md.append('|---|---|---|---|---|---|')


def run70(p, name):
    f = os.path.join(ws(p), 'build4d/v2/runs', name + '.pkl')
    if not os.path.exists(f):
        return None
    h = {x['step']: x for x in pickle.load(open(f, 'rb'))['hist']}
    if 70 not in h:
        return None
    return h[70]['meas'][2]['crown_distal'], -h[70]['meas'][8]['crown_buccal']


for p in ('P3', 'P5'):
    for sc, pat in (('0.7', 'rob_pdl070_%s_tad_baseline'), ('1.0', 'v2_%s_tad_baseline'),
                    ('1.4', 'rob_pdl140_%s_tad_baseline')):
        b, q = run70(p, pat % 'buccal'), run70(p, pat % 'palatal')
        if b is None or q is None:
            continue
        md.append(f"| {p[1:]} | {sc} | {b[0]:.3f} | {q[0]:.3f} | {q[0]-b[0]:+.3f} | {q[1]-b[1]:+.3f} |")

open(os.path.join(out, 'tables_replication.md'), 'w').write('\n'.join(md) + '\n')

# long-format appendix: every completed run behind the tables above
with open(os.path.join(out, 'replication_runs.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['patient', 'group', 'site', 'set', 'steps', 'd17_10', 'd17_20', 'd17_70',
                'relapse17_20_40', 'd16_70', 'inc_retraction_70', 'inc_expression_70', 'tip17_10'])
    for p in PATIENTS:
        if C[p]:
            w.writerow([p[1:], 'control', '-', 'v2', 70, '', '', f"{C[p]['d17_70']:.4f}", '', '',
                        f"{C[p]['inc_retraction_70']:.4f}", '', ''])
        for tag, S in (('v2', V[p]), ('hk14', H[p])):
            for g in GROUPS:
                for s, m in sorted(S.get(g, {}).items()):
                    w.writerow([p[1:], g.replace('_tad', ''), LABEL.get(s, s), tag, m.get('steps'),
                                *[f"{m[k]:.4f}" if k in m else '' for k in
                                  ('d17_10', 'd17_20', 'd17_70', 'relapse17_20_40', 'd16_70',
                                   'inc_retraction_70', 'inc_expression_70', 'tip17_10')]])

print(f"wrote {out}/tables_replication.md and replication_runs.csv")
