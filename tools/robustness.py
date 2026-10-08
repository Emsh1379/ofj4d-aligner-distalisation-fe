"""Is the group comparison robust to the most influential assumed input?

Compares the three anchorage groups of the final model (alpha 0.40) at the baseline
PDL Poisson's ratio (0.45) and at 0.49: each group's published metrics, the TAD
gain over the control, and the incisor-retraction ordering.

usage: python tools/robustness.py [BASE_PREFIX] [ALT_CONTROL] [ALT_PREFIX]
  defaults: build4d/final/run_  build4d/variants/sensblend_nu049.pkl
            build4d/variants/sensblend_nu049_
"""
import sys, os, pickle
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(HERE)

BASE = sys.argv[1] if len(sys.argv) > 1 else 'build4d/final/run_'
ALT_CONTROL = sys.argv[2] if len(sys.argv) > 2 else 'build4d/variants/sensblend_nu049.pkl'
ALT = sys.argv[3] if len(sys.argv) > 3 else 'build4d/variants/sensblend_nu049_'
GROUPS = ('control', 'buccal_tad', 'palatal_tad')
PAPER = {'control': {'17@10': 0.97, '17@70': 0.73, '16@70': 0.62, '11@70': 0.87},
         'buccal_tad': {'17@10': 1.11, '11@70': 1.71},
         'palatal_tad': {'17@10': 1.11, '17@70': 1.31, '16@70': 1.00, '11@70': 1.88}}
KEYS = ('17@10', '17@70', '16@70', '11@70')


def metrics(path):
    if not os.path.exists(path):
        return None, 0
    h = pickle.load(open(path, 'rb'))['hist']
    if len(h) < 70:
        return None, len(h)
    m = lambda s: h[s-1]['meas']
    return {'17@10': m(10)[2]['crown_distal'], '17@70': m(70)[2]['crown_distal'],
            '16@70': m(70)[3]['crown_distal'], '11@70': -m(70)[8]['crown_buccal']}, 70


def main():
    sets = {'nu 0.45 (baseline)': {g: f'{BASE}{g}.pkl' for g in GROUPS},
            'nu 0.49': {'control': ALT_CONTROL, 'buccal_tad': f'{ALT}buccal_tad.pkl',
                        'palatal_tad': f'{ALT}palatal_tad.pkl'}}
    res = {}
    for lab, paths in sets.items():
        print(f"\n{lab}")
        print(f"  {'group':<13}" + "".join(f"{k:>9}" for k in KEYS) + "   rms rel err vs paper")
        for g in GROUPS:
            r, n = metrics(paths[g])
            if r is None:
                print(f"  {g:<13}  incomplete ({n} of 70 steps): {paths[g]}"); continue
            res[(lab, g)] = r
            e = np.sqrt(np.mean([((r[k]-v)/v)**2 for k, v in PAPER[g].items()]))
            print(f"  {g:<13}" + "".join(f"{r[k]:>9.3f}" for k in KEYS) + f"   {e:.3f}")
        if all((lab, g) in res for g in GROUPS):
            c = res[(lab, 'control')]
            for g in ('buccal_tad', 'palatal_tad'):
                t = res[(lab, g)]
                print(f"  TAD gain {g:<12} incisor {t['11@70']-c['11@70']:+.3f} mm, "
                      f"FDI 17 @70 {t['17@70']-c['17@70']:+.3f} mm, FDI 16 @70 {t['16@70']-c['16@70']:+.3f} mm")
            order = sorted(GROUPS, key=lambda g: -res[(lab, g)]['11@70'])
            print(f"  incisor retraction ordering: {' > '.join(order)}"
                  f"{'  (as the paper)' if order == ['palatal_tad', 'buccal_tad', 'control'] else '  (paper: palatal_tad > buccal_tad > control)'}")


if __name__ == '__main__':
    main()
