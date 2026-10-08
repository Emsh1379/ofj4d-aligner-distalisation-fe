"""Mesh convergence of the staging model: the same 20 steps of the buccal TAD group on refined meshes.

usage: python tools/compare_meshconv.py [--steps 10 20] [--json out.json]

Levels (build4d/v2/jobs_meshconv.txt): M0 repaired meshes (steps 1-20 of v2_buccal_tad_baseline), M1 PDL 6 element
layers, M2 PDL in-plane x2, M3 aligner surface 0.30 mm (TetGen), M4 all refinements. For each quantity of interest
the change of every level from M0, and -- when M4 exists -- from M4 (the finest).
"""
import sys, os, re, pickle, json, argparse
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, 'ofj4d')); os.chdir(HERE)

ap = argparse.ArgumentParser()
ap.add_argument('--steps', type=int, nargs='+', default=[10, 20]); ap.add_argument('--json')
a = ap.parse_args()
LEVELS = [('M0', 'repaired (production)', 'build4d/v2/runs/v2_buccal_tad_baseline.pkl', 'build4d/v2/logs/v2_buccal_tad_baseline.log'),
          ('M1', 'PDL 6 layers', 'build4d/v2/runs/mc_M1_t6.pkl', 'build4d/v2/logs/mc_M1_t6.log'),
          ('M2', 'PDL in-plane x2', 'build4d/v2/runs/mc_M2_h2.pkl', 'build4d/v2/logs/mc_M2_h2.log'),
          ('M3', 'aligner 0.30 mm', 'build4d/v2/runs/mc_M3_al.pkl', 'build4d/v2/logs/mc_M3_al.log'),
          ('M4', 'all refined', 'build4d/v2/runs/mc_M4_all.pkl', 'build4d/v2/logs/mc_M4_all.log')]
QOI = [('d17', 'FDI17 distal'), ('r17', 'FDI17 root distal'), ('tip17', 'FDI17 tip deg'), ('d16', 'FDI16 distal'),
       ('inc', 'incisor retr.'), ('Fc', 'contact N'), ('eps99', 'p99 strain %')]


def qoi(path, log, s):
    D = pickle.load(open(path, 'rb')); h = {x['step']: x for x in D['hist']}
    if s not in h:
        return None
    m = h[s]['meas']
    q = dict(d17=m[2]['crown_distal']*1e3, r17=m[2]['root_distal']*1e3, tip17=m[2]['rot_coronal'], d16=m[3]['crown_distal']*1e3,
             inc=-m[8]['crown_buccal']*1e3, Fc=h[s]['info'][-1]['Fc'], eps99=h[s]['info'][-1]['emax']*100)
    t = None
    if os.path.exists(log):
        for l in open(log):
            mm = re.search(rf'step\s+{s}/\d+ .*\((\d+)s, total', l)
            if mm: t = int(mm.group(1))
    q['sec_per_step'] = t
    dofs = None
    if os.path.exists(log):
        for l in open(log):
            mm = re.search(r'reduced dofs (\d+)', l)
            if mm: dofs = int(mm.group(1)); break
    q['dofs'] = dofs
    return q


out = {}
for s in a.steps:
    res = {k: qoi(p, lg, s) for k, _, p, lg in LEVELS if os.path.exists(p)}
    res = {k: v for k, v in res.items() if v is not None}
    if 'M0' not in res:
        print(f"step {s}: M0 not available yet"); continue
    ref_fine = res.get('M4')
    print(f"\n=== step {s}   (movements in um; change from M0" + (" / from M4" if ref_fine else "") + ")")
    print(f"  {'level':<24}{'dofs':>9}{'s/step':>8}" + "".join(f"{lab:>20}" for _, lab in QOI))
    for k, name, _, _ in LEVELS:
        if k not in res:
            continue
        r = res[k]
        print(f"  {k + ' ' + name:<24}{(r['dofs'] or 0):>9}{(r['sec_per_step'] or 0):>8}" + "".join(f"{r[q]:>20.3f}" for q, _ in QOI))
        if k != 'M0':
            cells = []
            for q, _ in QOI:
                d0 = r[q] - res['M0'][q]
                cells.append(f"{d0:+.3f}" + (f" / {r[q]-ref_fine[q]:+.3f}" if ref_fine and k != 'M4' else ""))
            print(f"  {'':<24}{'':>9}{'':>8}" + "".join(f"{c:>20}" for c in cells))
    out[s] = res
if a.json:
    json.dump(out, open(a.json, 'w'), indent=1); print('->', a.json)
