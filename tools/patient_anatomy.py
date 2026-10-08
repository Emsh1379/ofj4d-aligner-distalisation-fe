"""Anatomy and miniscrew placeability of the patients built by pipeline/patient_build.sh, side by side.

usage: python tools/patient_anatomy.py [N ...] [--json out.json]     (default: every ../ofj4d_patients/P*)

Per patient: half-arch length and width at the first molars (crown centres, occlusal-plane frame of
ofj4d/arch4d.py), mean PDL root area, aligner size; and for each miniscrew site (build4d/tad_bone.pkl and
build4d/tad_positions/tad_*.pkl) the placeability criterion of docs/FINDINGS.md section 33: largest screw
diameter d in {1.2, 1.5, 1.8} mm leaving >= 0.5 mm to the PDL on both sides of the gap, and the cortical
insertion point outside the root envelope.
"""
import sys, os, glob, pickle, json, argparse
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAT = os.path.normpath(os.path.join(HERE, '..', 'ofj4d_patients'))

SITES = [('baseline', 'tad_bone.pkl', '16/17, 4 mm (paper)'), ('h1617_a2', 'tad_positions/tad_h1617_a2.pkl', '16/17, 2 mm'),
         ('h1617_a6', 'tad_positions/tad_h1617_a6.pkl', '16/17, 6 mm'), ('h1617_a8', 'tad_positions/tad_h1617_a8.pkl', '16/17, 8 mm'),
         ('md1516_a4', 'tad_positions/tad_md1516_a4.pkl', '15/16, 4 mm'), ('md1415_a4', 'tad_positions/tad_md1415_a4.pkl', '14/15, 4 mm')]
CLEAR, DIAMS = 0.5, (1.2, 1.5, 1.8)


def verdict(s):
    gap = s['gap'][1] - s['gap'][0]
    if s['beyond_root'] < 0:
        return gap, None, 'inside root envelope'
    ok = [d for d in DIAMS if (gap - d)/2 >= CLEAR]
    return gap, (max(ok) if ok else None), (f'placeable, d <= {max(ok)} mm' if ok else 'gap too narrow')


def patient(w):
    b = os.path.join(w, 'build4d')
    G = pickle.load(open(os.path.join(b, 'geom_fe.pkl'), 'rb'))
    IN = pickle.load(open(os.path.join(b, 'slim.pkl'), 'rb'))
    T = G['teeth']
    inc = 0.5*(IN['crown_c'][8] + IN['crown_c'][9])
    width16 = float(np.linalg.norm(IN['mb_cusp'][16] - IN['mb_cusp'][26]))
    c = [T[u]['crown_c'] for u in (8, 7, 6, 5, 4, 3, 2)]
    arch_len = float(sum(np.linalg.norm(c[i+1] - c[i]) for i in range(6)) + np.linalg.norm(c[0] - inc))
    out = dict(width_mb16_26=width16, half_arch_length=arch_len,
               root_area={int(T[u]['fdi']): float(T[u]['pdl_root_area']) for u in T},
               sites={})
    for key, f, label in SITES:
        p = os.path.join(b, f)
        if not os.path.exists(p):
            continue
        S = pickle.load(open(p, 'rb'))['sides']
        out['sites'][key] = {}
        for side in ('buccal', 'palatal'):
            gap, d, v = verdict(S[side])
            out['sites'][key][side] = dict(label=label, gap=gap, beyond_root=S[side]['beyond_root'],
                                           nearest_pdl=min(S[side]['dist_pdl16'], S[side]['dist_pdl17']),
                                           max_diameter=d, verdict=v)
    return out


ap = argparse.ArgumentParser()
ap.add_argument('patients', nargs='*')
ap.add_argument('--json')
a = ap.parse_args()
ws = [os.path.join(PAT, f'P{n}') for n in a.patients] or sorted(glob.glob(os.path.join(PAT, 'P*')), key=lambda s: int(s.rsplit('P', 1)[1]))
R = {}
for w in ws:
    if not os.path.exists(os.path.join(w, 'build4d', 'tad_positions', 'tad_md1516_a4.pkl')):
        print(f'{os.path.basename(w)}: not built'); continue
    R[os.path.basename(w)] = patient(w)
print(f"{'patient':<8}{'16-26 width':>12}{'half arch':>11}{'root area 17/16/13/11 (mm2)':>32}")
for k, r in R.items():
    ra = r['root_area']
    print(f"{k:<8}{r['width_mb16_26']:>10.1f}mm{r['half_arch_length']:>9.1f}mm{'':>6}{ra[17]:6.0f}{ra[16]:6.0f}{ra[13]:6.0f}{ra[11]:6.0f}")
print(f"\nplaceability (>= {CLEAR} mm screw-to-PDL each side, insertion outside the roots): gap B/P mm, verdict B | P")
for key, _, label in SITES:
    print(f"  {label}")
    for k, r in R.items():
        s = r['sites'].get(key)
        if s:
            print(f"    {k:<6}{s['buccal']['gap']:5.2f}/{s['palatal']['gap']:4.2f}   {s['buccal']['verdict']:<24} | {s['palatal']['verdict']}")
if a.json:
    json.dump(R, open(a.json, 'w'), indent=1); print('->', a.json)
