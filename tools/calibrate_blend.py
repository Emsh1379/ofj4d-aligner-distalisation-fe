"""Calibrate the socket-kinematics blend (OFJ_TIP_BLEND) on the CONTROL group only.

alpha = 0 is bodily translation of each socket, alpha = 1 tipping about the apical
third; the crown centre lands on its prescription for every alpha (arch4d.py).
The two TAD groups are NOT used here, so they remain an out-of-sample test of the
calibrated model.

Metrics are the paper's control-group values, measured exactly as tools/analyse.py
measures them (crown point: cusp tip / incisal edge). Tipping at step 10 is
reported but kept out of the primary objective: no alpha in [0, 1] reaches the
paper's 5.92 deg (alpha = 1 gives 3.8), so including it only drags alpha to the
upper end.

usage: python tools/calibrate_blend.py [ALPHA:RUN.pkl ...]
"""
import sys, os, pickle
import numpy as np
from scipy.interpolate import PchipInterpolator
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(HERE)

V = 'build4d/variants'
DEFAULT = [(0.0, 'build4d/run_control.pkl'), (0.4, f'{V}/cal_a0.40_control.pkl'),
           (0.6, f'{V}/cal_a0.60_control.pkl'), (1.0, f'{V}/tip_run_control.pkl')]
PAPER = {'17@10': 0.97, 'drift20-40': 0.89, '17@70': 0.73, '16@70': 0.62, '11@70': 0.87,
         'tip@10': 5.92}
PRIMARY = ['17@10', 'drift20-40', '17@70', '16@70', '11@70']


def metrics(h):
    m = lambda s: h[s-1]['meas']
    return {'17@10': m(10)[2]['crown_distal'],
            'tip@10': abs(m(10)[2]['rot_coronal']),
            'drift20-40': m(20)[2]['crown_distal'] - m(40)[2]['crown_distal'],
            '17@70': m(70)[2]['crown_distal'],
            '16@70': m(70)[3]['crown_distal'],
            '11@70': -m(70)[8]['crown_buccal']}


def objective(vals, keys):
    return np.sqrt(np.mean([((vals[k] - PAPER[k])/PAPER[k])**2 for k in keys]))


def main():
    pts = ([(float(a), p) for a, _, p in (x.partition(':') for x in sys.argv[1:])]
           if len(sys.argv) > 1 else DEFAULT)
    rows = []
    for a, p in pts:
        if not os.path.exists(p):
            print(f"alpha {a:.2f}: {p} missing -- skipped"); continue
        h = pickle.load(open(p, 'rb'))['hist']
        if len(h) < 70:
            print(f"alpha {a:.2f}: {p} has {len(h)} steps -- skipped"); continue
        rows.append((a, metrics(h)))
    rows.sort()
    keys = list(PAPER)
    print(f"\n{'alpha':<8}" + "".join(f"{k:>12}" for k in keys) + f"{'rms rel err':>13}")
    print(f"{'paper':<8}" + "".join(f"{PAPER[k]:>12.2f}" for k in keys))
    for a, r in rows:
        print(f"{a:<8.2f}" + "".join(f"{r[k]:>12.3f}" for k in keys)
              + f"{objective(r, PRIMARY):>13.3f}")
    if len(rows) < 3:
        print("\nneed at least 3 complete runs to interpolate"); return

    A = np.array([a for a, _ in rows])
    grid = np.linspace(A.min(), A.max(), 1001)
    curves = {k: PchipInterpolator(A, [r[k] for _, r in rows])(grid) for k in keys}
    obj = np.array([objective({k: curves[k][i] for k in keys}, PRIMARY)
                    for i in range(len(grid))])
    obj_t = np.array([objective({k: curves[k][i] for k in keys}, PRIMARY + ['tip@10'])
                      for i in range(len(grid))])
    i = int(obj.argmin()); j = int(obj_t.argmin())
    print(f"\nbest alpha (displacements): {grid[i]:.3f}   rms relative error {obj[i]:.3f}")
    print("  predicted: " + "  ".join(f"{k} {curves[k][i]:.3f}" for k in keys))
    print(f"best alpha incl. tipping:   {grid[j]:.3f}   rms relative error {obj_t[j]:.3f}")
    print("\nper-metric alpha that would match the paper on its own:")
    for k in keys:
        c = curves[k] - PAPER[k]
        s = np.where(np.sign(c[:-1]) != np.sign(c[1:]))[0]
        print(f"  {k:<11} " + (", ".join(f"{grid[x]:.2f}" for x in s) if len(s)
                               else f"none in [{A.min():.1f}, {A.max():.1f}] "
                                    f"(range {curves[k].min():.2f} to {curves[k].max():.2f})"))


if __name__ == '__main__':
    main()
