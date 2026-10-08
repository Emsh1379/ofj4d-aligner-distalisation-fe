"""Miniscrew-position study: summarise the 70-step runs of each position variant per group.

usage: python tools/compare_tad.py [--set old|v2|a000|a060|nu030|nu049] [--bone rigid|deformable] [--json out.json]

Reads build4d/tad_positions/run_<group>_<variant>_<bone>.pkl; the baseline (16/17, 4 mm apical) is
build4d/final/run_<group>.pkl for rigid bone and build4d/final/run_<group>_bone.pkl for deformable
(or build4d/tad_positions/run_<group>_baseline_<bone>.pkl when present).

Per variant (movement of each tooth's crown point along its own axes, run4d.measure):
  FDI 17 distal at steps 10 / 20 / 70, its mesial relapse over steps 20-40 (the paper's anchorage loss),
  FDI 16 distal at 30 / 70, incisor (FDI 11) retraction at 70 and its expression (achieved / prescribed),
  FDI 17 tipping at step 10, capped solves; and each quantity's change from the baseline.
"""
import sys, os, glob, pickle, json, argparse
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, 'ofj4d')); os.chdir(HERE)
import spec

ap = argparse.ArgumentParser()
ap.add_argument('--bone', default='rigid', choices=['rigid', 'deformable'])
ap.add_argument('--set', default='old', help="'old' (build4d/tad_positions, shipped meshes), 'v2' (build4d/v2/runs, repaired meshes, K_N 1000), "
               "or a robustness tag such as a000, a060, nu030, nu049 (build4d/v2/runs/rob_<tag>_...)")
ap.add_argument('--json')
a = ap.parse_args()

LABEL = {'baseline': '16/17, 4 mm (paper)', 'md1516_a4': '15/16, 4 mm', 'md1415_a4': '14/15, 4 mm',
         'h1617_a2': '16/17, 2 mm', 'h1617_a6': '16/17, 6 mm', 'h1617_a8': '16/17, 8 mm'}
ORDER = ['baseline', 'md1516_a4', 'md1415_a4', 'h1617_a2', 'h1617_a6', 'h1617_a8']


def metrics(D):
    h = {x['step']: x for x in D['hist']}
    last = max(h)
    m = lambda s, u, k: h[s]['meas'][u][k] if s in h else float('nan')
    d17 = [m(s, 2, 'crown_distal') for s in range(20, 41)] if all(s in h for s in range(20, 41)) else []   # needs the whole window
    presc11 = sum(spec.prescribed_increment(8, s) for s in range(1, last + 1))
    capped = sum(1 for x in D['hist'] for i in x['info'] if not i['converged'])
    return dict(steps=last,
                d17_10=m(10, 2, 'crown_distal'), d17_20=m(20, 2, 'crown_distal'), d17_70=m(70, 2, 'crown_distal'),
                relapse17_20_40=(d17[0] - min(d17)) if d17 else float('nan'),
                d16_30=m(30, 3, 'crown_distal'), d16_70=m(70, 3, 'crown_distal'),
                inc_retraction_70=-m(70, 8, 'crown_buccal'),
                inc_expression_70=(-m(70, 8, 'crown_buccal')/presc11*100 if presc11 > 0 and 70 in h else float('nan')),
                tip17_10=m(10, 2, 'rot_coronal'), capped=capped)


COLS = [('d17_10', 'FDI17 @10'), ('d17_20', 'FDI17 @20'), ('relapse17_20_40', 'relapse 20-40'), ('d17_70', 'FDI17 @70'),
        ('d16_30', 'FDI16 @30'), ('d16_70', 'FDI16 @70'), ('inc_retraction_70', 'incisor @70'),
        ('inc_expression_70', 'expr. %'), ('tip17_10', 'tip17 @10'), ('capped', 'capped')]
out = {}
for g in ('buccal_tad', 'palatal_tad'):
    paths = {}
    for v in ORDER:
        if a.set == 'old':
            p = f'build4d/tad_positions/run_{g}_{v}_{a.bone}.pkl'
            if v == 'baseline' and not os.path.exists(p):
                p = f'build4d/final/run_{g}.pkl' if a.bone == 'rigid' else f'build4d/final/run_{g}_bone.pkl'
        elif a.set == 'v2':
            p = f'build4d/v2/runs/v2_{g}_{v}.pkl'
        else:
            p = f'build4d/v2/runs/rob_{a.set}_{g}_{v}.pkl'
        if os.path.exists(p):
            paths[v] = p
    if not paths:
        continue
    res = {v: metrics(pickle.load(open(p, 'rb'))) for v, p in paths.items()}
    base = res.get('baseline')
    print(f"\n=== {g}, {a.bone} bone, set {a.set}   (mm, deg; second line: change from the paper's site)")
    print(f"  {'miniscrew':<22}{'steps':>6}" + "".join(f"{lab:>15}" for _, lab in COLS))
    for v in ORDER:
        if v not in res:
            continue
        r = res[v]
        print(f"  {(LABEL[v] + ('' if r['steps'] >= 70 else ' *')):<22}{r['steps']:>6}" + "".join(f"{r[k]:>15.3f}" if k != 'capped' else f"{r[k]:>15d}" for k, _ in COLS))
        if base is not None and v != 'baseline':
            print(f"  {'':<22}{'':>6}" + "".join(f"{r[k]-base[k]:>+15.3f}" if k != 'capped' else f"{'':>15}" for k, _ in COLS))
    if any(r['steps'] < 70 for r in res.values()):
        print("  * run still in progress: quantities at steps it has not reached are nan")
    out[g] = dict(paths=paths, metrics=res)
if a.json:
    json.dump(out, open(a.json, 'w'), indent=1); print('->', a.json)
