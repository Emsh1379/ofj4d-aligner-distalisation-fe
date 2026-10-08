"""Sensitivity of the 70-step control run to the ASSUMED inputs: PDL Poisson's
ratio, normal contact penalty stiffness, attachment size.

Reads the checkpoints written by .runs/sens_queue.sh; incomplete runs are reported
up to the last step they reached.

usage: python tools/sensitivity.py
"""
import sys, os, pickle
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(HERE)

V = 'build4d/variants'
# defaults: the final (blended, alpha 0.40) model; pass BASELINE and PREFIX for others,
# e.g. build4d/run_control.pkl build4d/variants/sens_  for the bodily model
BASE = sys.argv[1] if len(sys.argv) > 1 else 'build4d/final/run_control.pkl'
PREFIX = sys.argv[2] if len(sys.argv) > 2 else f'{V}/sensblend_'
RUNS = [('baseline (nu 0.45, K_N 2000, att 1x)', BASE),
        ('PDL nu 0.30', f'{PREFIX}nu030.pkl'),
        ('PDL nu 0.49', f'{PREFIX}nu049.pkl'),
        ('K_N 1000', f'{PREFIX}kn1000.pkl'),
        ('K_N 4000', f'{PREFIX}kn4000.pkl'),
        ('attachments 0.75x', f'{PREFIX}A075.pkl'),
        ('attachments 1.25x', f'{PREFIX}A125.pkl')]
PAPER = {'17@10': 0.97, 'tip@10': 5.92, 'drift20-40': 0.89, '17@70': 0.73, '16@70': 0.62,
         '11@70': 0.87}


def metrics(h):
    n = len(h)
    m = lambda s: h[s-1]['meas'] if s <= n else None
    out = {}
    out['17@10'] = m(10)[2]['crown_distal'] if m(10) else np.nan
    out['tip@10'] = abs(m(10)[2]['rot_coronal']) if m(10) else np.nan
    out['17@20'] = m(20)[2]['crown_distal'] if m(20) else np.nan
    out['drift20-40'] = (m(20)[2]['crown_distal'] - m(40)[2]['crown_distal']) if m(40) else np.nan
    out['17@70'] = m(70)[2]['crown_distal'] if m(70) else np.nan
    out['16@70'] = m(70)[3]['crown_distal'] if m(70) else np.nan
    out['11@70'] = -m(70)[8]['crown_buccal'] if m(70) else np.nan
    info = [i for e in h for i in e['info']]
    out['capped'] = sum(not i['converged'] for i in info)
    out['final unconv'] = sum(not e['info'][-1]['converged'] for e in h)
    out['steps'] = n
    return out


def main():
    cols = ['17@10', 'tip@10', '17@20', 'drift20-40', '17@70', '16@70', '11@70']
    rows = []
    for name, p in RUNS:
        if not os.path.exists(p):
            rows.append((name, None)); continue
        rows.append((name, metrics(pickle.load(open(p, 'rb'))['hist'])))
    base = rows[0][1]
    if base is None:
        print(f"baseline {RUNS[0][1]} not found: pass it as the first argument, e.g.\n"
              "  python tools/sensitivity.py build4d/variants/cal_a0.40_control.pkl")
        return
    w = 38
    print(f"{'control group, crown movement (mm; tipping deg)':<{w}}"
          + "".join(f"{c:>12}" for c in cols) + f"{'steps':>7}{'capped':>8}{'final unc':>10}")
    print(f"{'paper':<{w}}" + "".join(f"{PAPER.get(c, np.nan):>12.2f}" if c in PAPER
                                     else f"{'-':>12}" for c in cols))
    for name, r in rows:
        if r is None:
            print(f"{name:<{w}}  (not run yet)"); continue
        print(f"{name:<{w}}" + "".join(f"{r[c]:>12.3f}" if np.isfinite(r[c]) else f"{'-':>12}"
                                       for c in cols)
              + f"{r['steps']:>7}{r['capped']:>8}{r['final unconv']:>10}")
    print("\nchange from baseline (mm; deg)")
    for name, r in rows[1:]:
        if r is None:
            continue
        print(f"{name:<{w}}" + "".join(f"{r[c]-base[c]:>+12.3f}" if np.isfinite(r[c]) else f"{'-':>12}"
                                       for c in cols))


if __name__ == '__main__':
    main()
