"""Staging diagnostic: per-tooth lag behind the prescription, and the rigid wear-in
shift of the appliance, for one or more runs.

A uniform lag on every tooth that equals the wear-in shift means the appliance and
dentition are moving as one block: the common activation is removed by the seat,
and only the differential part is delivered as force.

usage: python tools/lag_seat.py RUN.pkl[:ALIGNER.pkl] [...] [--steps 10,20,70]
"""
import sys, os, pickle
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, 'ofj4d')); sys.path.insert(0, os.path.join(HERE, 'ofj'))
os.chdir(HERE)
import spec

PAPER_STEP10_17 = {'control': (0.97, 5.92), 'buccal_tad': (1.11, 6.75),
                   'palatal_tad': (1.11, 6.75)}


def ap(T, p):
    return p @ T[:3, :3].T + T[:3, 3]


def report(run, aligner, steps):
    os.environ['OFJ_ALIGNER'] = aligner
    from model4d import Static
    from arch4d import ArchForm
    from morph4d import Morpher, kabsch
    S = Static(aligner=aligner); arch = ArchForm(S); mo = Morpher(S, arch)
    D = pickle.load(open(run, 'rb')); h = D['hist']; g = D['group']
    cc = {u: S.tooth[u]['crown_c'] for u in S.keep}
    print(f"\n== {run}  [{g}, aligner {os.path.basename(aligner)}, "
          f"{len(S.aligner_tets)} aligner tets, {len(h)} steps]")
    m10 = h[min(10, len(h))-1]['meas'][2]
    pv = PAPER_STEP10_17.get(g, (np.nan, np.nan))
    ratio = m10['root_distal']/m10['crown_distal'] if abs(m10['crown_distal']) > 1e-9 else np.nan
    print(f"   FDI 17 at step {min(10, len(h))}: crown {m10['crown_distal']:+.3f} mm "
          f"(paper {pv[0]:.2f}), root {m10['root_distal']:+.3f}, root/crown {ratio:.2f}, "
          f"tipping {abs(m10['rot_coronal']):.2f} deg (paper {pv[1]:.2f})")
    print("   step" + "".join(f"{'L'+str(spec.FDI[u]):>7}" for u in S.keep)
          + "   lag sd   seat mm  seat deg")
    for s in steps:
        if s > len(h):
            continue
        T = h[s-1]['T']; lag = []
        for u in S.keep:
            Tp = arch.prescribed_transform(u, s)
            dp = ap(Tp, cc[u]) - cc[u]; da = ap(T[u], cc[u]) - cc[u]
            if u in (7, 8):
                d = -arch._palatal(u)
            elif np.linalg.norm(dp) > 1e-9:
                d = dp/np.linalg.norm(dp)
            else:
                t = arch.frame(arch.s0[u])[1]; d = arch.ex*t[0] + arch.ey*t[1]
            lag.append((dp - da) @ d)
        pre = np.zeros_like(mo.inner_target_rest); act = np.zeros_like(pre)
        for u in S.keep:
            k = mo.inner_tooth == u
            pre[k] = ap(arch.prescribed_transform(u, s), mo.inner_target_rest[k])
            act[k] = ap(T[u], mo.inner_target_rest[k])
        R, t = kabsch(pre, act); c = pre.mean(0); sh = R @ c + t - c
        ang = np.degrees(np.arccos(np.clip((np.trace(R)-1)/2, -1, 1)))
        print(f"   {s:4d}" + "".join(f"{x:7.3f}" for x in lag)
              + f"   {np.std(lag):6.3f}   {np.linalg.norm(sh):6.3f}  {ang:7.2f}")


def main():
    args = sys.argv[1:]
    steps = [10, 20, 30, 40, 50, 60, 70]
    if '--steps' in args:
        i = args.index('--steps'); steps = [int(x) for x in args[i+1].split(',')]
        del args[i:i+2]
    for a in args:
        run, _, al = a.partition(':')
        report(run, al or 'build4d/aligner.pkl', steps)


if __name__ == '__main__':
    main()
