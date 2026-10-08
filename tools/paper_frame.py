"""Tooth movement in the paper's coordinate system (Mao et al. 2024).

Origin at the central incisors' midpoint; X parallel to the line through the upper
first molars' mesiobuccal cusps; Y perpendicular to it in the occlusal plane
(posterior +); Z the occlusal-plane normal (occlusal +, so + = extrusion).
run4d.py reports movement along each tooth's OWN axes; for a proclined incisor a
horizontal retraction then shows up partly as "intrusion". This tool uses the
global frame, as the paper does.

Angles are changes of the crown-point-to-apex axis within the tooth's vertical
plane: along the arch tangent for canine..molars (+ = distal crown tip), along the
labial direction for incisors (+ = labial crown inclination).

usage: python tools/paper_frame.py PREFIX [PREFIX ...]     e.g. build4d/final/run_
"""
import sys, os, pickle
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, 'ofj4d')); sys.path.insert(0, os.path.join(HERE, 'ofj'))
os.chdir(HERE)
from model4d import Static
from arch4d import ArchForm
from run4d import crown_point, root_point

S = Static(); A = ArchForm(S)
EX, EY, EN = A.ex, A.ey, A.en
CP0 = {u: crown_point(S, u) for u in S.keep}
RP0 = {u: root_point(S, u) for u in S.keep}


def plane_dir(u):
    if u in (7, 8):
        return A._palatal(u)                                  # labial, in the occlusal plane
    t = A.frame(A.s0[u])[1]
    d = EX*t[0] + EY*t[1]; d = d - (d @ EN)*EN
    return d/np.linalg.norm(d)                                # distal along the arch


def movement(T, u):
    R, t = T[u][:3, :3], T[u][:3, 3]
    cp = R @ CP0[u] + t; rp = R @ RP0[u] + t
    d = plane_dir(u)
    ang = lambda la: np.degrees(np.arctan2(la @ d, la @ EN))
    dc = cp - CP0[u]; dr = rp - RP0[u]
    return dict(sag=dc @ EY, ver=dc @ EN, crown_pl=dc @ d, root_pl=dr @ d,
                tip=ang(cp - rp) - ang(CP0[u] - RP0[u]))


PAPER = {'control':     {'17 sag@10': 0.97, '17 tip@10': 5.92, '17 sag@70': 0.73, '16 sag@70': 0.62,
                         '11 retr@70': 0.87, '11 vert@70': -1.08, '11 incl@70': +0.65},
         'buccal_tad':  {'17 sag@10': 1.11, '17 tip@10': 6.75, '11 retr@70': 1.71,
                         '11 vert@70': -1.80, '11 incl@70': +1.9},
         'palatal_tad': {'17 sag@10': 1.11, '17 tip@10': 6.75, '17 sag@70': 1.31, '16 sag@70': 1.00,
                         '11 retr@70': 1.88, '11 vert@70': +0.3, '11 incl@70': -6.0}}
COLS = ['17 sag@10', '17 tip@10', '17 sag@70', '16 sag@70', '11 retr@70', '11 vert@70', '11 incl@70']


def main():
    for prefix in sys.argv[1:]:
        print(f"\n=== {prefix}<group>.pkl   (mm, deg; tip + distal crown, incl + labial crown, vert + extrusion)")
        print(f"  {'':<22}" + "".join(f"{c:>11}" for c in COLS) + f"{'13 tip@20':>11}{'13 crown/root@20':>18}")
        for g in ('control', 'buccal_tad', 'palatal_tad'):
            p = f'{prefix}{g}.pkl'
            if not os.path.exists(p):
                continue
            h = pickle.load(open(p, 'rb'))['hist']
            if len(h) < 70:
                print(f"  {g:<22}  ({len(h)} steps, skipped)"); continue
            m10 = lambda u: movement(h[9]['T'], u); m70 = lambda u: movement(h[69]['T'], u)
            c20 = movement(h[19]['T'], 6)
            v = {'17 sag@10': m10(2)['sag'], '17 tip@10': m10(2)['tip'], '17 sag@70': m70(2)['sag'],
                 '16 sag@70': m70(3)['sag'], '11 retr@70': m70(8)['sag'], '11 vert@70': m70(8)['ver'],
                 '11 incl@70': m70(8)['tip']}
            pv = PAPER[g]
            print(f"  {g+' paper':<22}" + "".join(f"{pv[c]:>11.2f}" if c in pv else f"{'-':>11}" for c in COLS)
                  + f"{'mesial':>11}{'-':>18}")
            print(f"  {g+' model':<22}" + "".join(f"{v[c]:>11.3f}" for c in COLS)
                  + f"{c20['tip']:>11.2f}{c20['crown_pl']:>9.2f}/{c20['root_pl']:<8.2f}")


if __name__ == '__main__':
    main()
