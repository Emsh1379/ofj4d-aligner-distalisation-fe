"""Rigid vs deformable alveolar bone: the same staging group run twice, differing only in OFJ_BONE.

usage: python tools/compare_bone.py RIGID.pkl DEFORMABLE.pkl [--steps 10 20 40 70] [--json out.json]

Per tooth and step: crown / root movement along the tooth axes (distal, buccal, occlusal), tipping
(rotation about e_b) and total rotation, as stored by ofj4d/run4d.measure; plus the quantities the
anchorage comparison of Mao et al. uses -- second and first molar distalisation, incisor retraction
(palatal = -buccal crown movement of FDI 11) and expression (achieved / prescribed) -- and the
convergence record of both runs.
"""
import sys, os, pickle, json, argparse
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, 'ofj4d')); os.chdir(HERE)
import spec

ap = argparse.ArgumentParser()
ap.add_argument('rigid'); ap.add_argument('deformable')
ap.add_argument('--steps', type=int, nargs='+', default=[10, 20, 40, 70])
ap.add_argument('--json')
a = ap.parse_args()
R, D = pickle.load(open(a.rigid, 'rb')), pickle.load(open(a.deformable, 'rb'))
assert R.get('group') == D.get('group'), f"groups differ: {R.get('group')} vs {D.get('group')}"
hr = {h['step']: h for h in R['hist']}; hd = {h['step']: h for h in D['hist']}
common = sorted(set(hr) & set(hd))
print(f"group {R['group']}: rigid {len(hr)} steps, deformable {len(hd)} steps, compared on {len(common)} steps")
steps = [s for s in a.steps if s in hr and s in hd]
if common and common[-1] not in steps:
    steps.append(common[-1])

KEYS = [('crown_distal', 'crown distal'), ('crown_buccal', 'crown buccal'), ('crown_occlusal', 'crown occlusal'),
        ('root_distal', 'root distal'), ('rot_coronal', 'tip (deg)'), ('rot_total', 'rotation (deg)')]


def prescribed(u, step):
    return sum(spec.prescribed_increment(u, s) for s in range(1, step + 1))


out = dict(group=R['group'], steps={})
for s in steps:
    print(f"\n=== step {s}: rigid / deformable (difference)   [mm unless noted]")
    print(f"  {'FDI':>4} " + " ".join(f"{lab:>27}" for _, lab in KEYS))
    rec = {}
    for u in spec.HALF_ARCH_UNN:
        mr, md = hr[s]['meas'][u], hd[s]['meas'][u]
        cells = []
        for k, _ in KEYS:
            cells.append(f"{mr[k]:8.3f}/{md[k]:8.3f} ({md[k]-mr[k]:+7.3f})")
        print(f"  {spec.FDI[u]:>4} " + " ".join(cells))
        rec[spec.FDI[u]] = {k: (float(mr[k]), float(md[k])) for k, _ in KEYS}
    out['steps'][s] = rec

print("\n=== anchorage summary: rigid / deformable")
summary = {}
for s in steps:
    r17, d17 = hr[s]['meas'][2]['crown_distal'], hd[s]['meas'][2]['crown_distal']
    r16, d16 = hr[s]['meas'][3]['crown_distal'], hd[s]['meas'][3]['crown_distal']
    ri, di = -hr[s]['meas'][8]['crown_buccal'], -hd[s]['meas'][8]['crown_buccal']
    pi = prescribed(8, s)
    line = (f"  step {s:2d}: FDI 17 distal {r17:.3f} / {d17:.3f} mm ({100*(d17-r17)/max(abs(r17),1e-9):+.1f} %), "
            f"FDI 16 distal {r16:.3f} / {d16:.3f} mm, incisor retraction {ri:.3f} / {di:.3f} mm")
    if pi > 0:
        line += f" (expression {100*ri/pi:.1f} % / {100*di/pi:.1f} %)"
    print(line)
    summary[s] = dict(fdi17=(r17, d17), fdi16=(r16, d16), incisor_retraction=(ri, di), incisor_prescribed=pi)
out['summary'] = summary

def conv(hist):
    it = [i['nit'] for h in hist for i in h['info']]
    cap = sum(1 for h in hist for i in h['info'] if not i['converged'])
    fc = [h['info'][-1]['Fc'] for h in hist]; e99 = [h['info'][-1]['emax'] for h in hist]
    return dict(solves=len(it), median_it=float(np.median(it)), capped=cap, Fc_median=float(np.median(fc)),
                eps99_max=float(np.max(e99)))
cr, cd = conv([hr[s] for s in common]), conv([hd[s] for s in common])
print(f"\n=== convergence over steps {common[0]}-{common[-1]}")
for name, c in (('rigid', cr), ('deformable', cd)):
    print(f"  {name:<10}: {c['solves']} solves, median {c['median_it']:.0f} iterations, {c['capped']} at the cap, "
          f"median contact force {c['Fc_median']:.2f} N, max p99 PDL strain {c['eps99_max']*100:.1f} %")
out['convergence'] = dict(rigid=cr, deformable=cd)

# largest per-tooth difference over every common step
worst = (0.0, None)
for s in common:
    for u in spec.HALF_ARCH_UNN:
        for k in ('crown_distal', 'crown_buccal', 'crown_occlusal', 'root_distal', 'root_buccal'):
            d = abs(hd[s]['meas'][u][k] - hr[s]['meas'][u][k])
            if d > worst[0]:
                worst = (d, (s, spec.FDI[u], k))
print(f"\nlargest difference over all compared steps: {worst[0]*1e3:.1f} um (step {worst[1][0]}, FDI {worst[1][1]}, {worst[1][2]})")
out['largest_difference'] = dict(mm=worst[0], at=worst[1])
if a.json:
    json.dump(out, open(a.json, 'w'), indent=1); print('->', a.json)
