"""Does the patient 8 result replicate in other patients?  Rankings and effects per patient.

usage: python tools/compare_patients.py [N ...] [--json out.json]
       (default: patient 8 from this package + every ../ofj4d_patients/P* with runs)

For each patient runs tools/compare_tad.py --set v2 and --set hk14 in that patient's workspace (patient 8:
this package), then reports per group:
  * the site ranking on FDI 17 @70 and incisor retraction @70, and whether it equals patient 8's
    order 2 mm > paper > 6 mm > 8 mm > 14/15 (15/16 reported, not ranked, as in docs/FINDINGS.md 34-35);
  * effect of each site vs the paper site, and the miniscrew benefit over control;
  * first premolar minus canine precision cut at the paper site and at 6 mm;
  * palatal minus buccal at the paper site.
Only complete (70-step) runs are used; missing ones are shown as '-'.
"""
import sys, os, glob, json, subprocess, argparse, tempfile, pickle
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAT = os.path.normpath(os.path.join(HERE, '..', 'ofj4d_patients'))
RANK = ['h1617_a2', 'baseline', 'h1617_a6', 'h1617_a8', 'md1415_a4']
LAB = {'h1617_a2': '2 mm', 'baseline': 'paper', 'h1617_a6': '6 mm', 'h1617_a8': '8 mm', 'md1415_a4': '14/15', 'md1516_a4': '15/16'}
Q = [('d17_70', 'FDI 17 @70'), ('inc_retraction_70', 'incisor @70')]


def summary(ws, tag):
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        out = f.name
    subprocess.run([sys.executable, os.path.join(ws, 'tools', 'compare_tad.py'), '--set', tag, '--json', out],
                   cwd=ws, check=True, stdout=subprocess.DEVNULL)
    R = json.load(open(out)); os.unlink(out)
    return {g: {v: m for v, m in d['metrics'].items() if m['steps'] >= 70} for g, d in R.items()}


def control(ws):
    p = os.path.join(ws, 'build4d/v2/runs/v2_control.pkl')
    if not os.path.exists(p):
        return None
    D = pickle.load(open(p, 'rb')); h = {x['step']: x for x in D['hist']}
    return dict(d17_70=h[70]['meas'][2]['crown_distal'], inc_retraction_70=-h[70]['meas'][8]['crown_buccal']) if 70 in h else None


ap = argparse.ArgumentParser()
ap.add_argument('patients', nargs='*')
ap.add_argument('--json')
a = ap.parse_args()
ws = {'P8': HERE}
for w in sorted(glob.glob(os.path.join(PAT, 'P*')), key=lambda s: int(s.rsplit('P', 1)[1])):
    if os.path.basename(w) != 'P8' and glob.glob(os.path.join(w, 'build4d/v2/runs/v2_*.pkl')):
        ws[os.path.basename(w)] = w
if a.patients:
    ws = {k: v for k, v in ws.items() if k[1:] in a.patients}

OUT = {}
for P, w in ws.items():
    OUT[P] = dict(v2=summary(w, 'v2'), hk14=summary(w, 'hk14'), control=control(w))

f = lambda x: f'{x:+.3f}' if x is not None else '     -'
for g in ('buccal_tad', 'palatal_tad'):
    print(f"\n=== {g} ===")
    for k, lab in Q:
        print(f"\n  {lab}: value at the paper site, change from the paper site per site (mm), ranking")
        print(f"  {'patient':<8}{'control':>8}{'paper':>8}" + ''.join(f"{LAB[v]:>8}" for v in ['h1617_a2', 'h1617_a6', 'h1617_a8', 'md1516_a4', 'md1415_a4']) + "   ranking (best first)          same as 2>paper>6>8>14/15")
        for P, R in OUT.items():
            M = R['v2'].get(g, {}); b = M.get('baseline'); c = R['control']
            if b is None:
                print(f"  {P:<8}  no complete paper-site run"); continue
            row = f"  {P:<8}{(c[k] if c else float('nan')):8.3f}{b[k]:8.3f}" + ''.join(
                f"{f(M[v][k]-b[k]) if v in M else '       -':>8}" for v in ['h1617_a2', 'h1617_a6', 'h1617_a8', 'md1516_a4', 'md1415_a4'])
            have = [v for v in RANK if v in M]
            order = sorted(have, key=lambda v: -M[v][k])
            same = order == have if len(have) == len(RANK) else None
            print(row + f"   {' > '.join(LAB[v] for v in order):<36}{'yes' if same else ('NO' if same is False else 'incomplete')}")
    print(f"\n  first premolar minus canine precision cut (mm): FDI 17 / incisor @70")
    for P, R in OUT.items():
        M, H = R['v2'].get(g, {}), R['hk14'].get(g, {})
        cells = []
        for v in ('baseline', 'h1617_a6'):
            if v in M and v in H:
                cells.append(f"{LAB[v]} {H[v]['d17_70']-M[v]['d17_70']:+.3f} / {H[v]['inc_retraction_70']-M[v]['inc_retraction_70']:+.3f}")
        print(f"  {P:<8}" + ('   '.join(cells) if cells else '-'))
print("\n=== palatal minus buccal at the paper site, and miniscrew benefit over control (FDI 17 @70) ===")
for P, R in OUT.items():
    B, Pl, c = R['v2'].get('buccal_tad', {}).get('baseline'), R['v2'].get('palatal_tad', {}).get('baseline'), R['control']
    if B and Pl:
        print(f"  {P:<8}palatal-buccal FDI 17 {Pl['d17_70']-B['d17_70']:+.3f}, incisor {Pl['inc_retraction_70']-B['inc_retraction_70']:+.3f}"
              + (f";  benefit buccal {B['d17_70']-c['d17_70']:+.3f}, palatal {Pl['d17_70']-c['d17_70']:+.3f}" if c else ''))
if a.json:
    json.dump(OUT, open(a.json, 'w'), indent=1); print('->', a.json)
