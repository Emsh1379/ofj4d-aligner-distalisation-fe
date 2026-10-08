"""One aligner solve, both staging backends, identical configuration.

Builds step STEP of the final model exactly as ofj4d/run4d.py does (prescribed shape,
pairing, standoff reference, wear-in seat), then solves it with solve4d (custom) and
febio4d (FEBio), and compares tooth motions and wall time.

usage: OFJ_PRESCRIBE=tip OFJ_TIP_BLEND=0.40 python tools/febio_step_test.py [GROUP] [STEP] [--only febio|custom]
"""
import os, sys, time, pickle, argparse
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, 'ofj4d')); sys.path.insert(0, os.path.join(HERE, 'ofj')); os.chdir(HERE)
from model4d import Static, State
from arch4d import ArchForm
from morph4d import Morpher
from loads4d import Elastic
from contact import build_pairs
from run4d import crown_point, root_point
import solve4d, febio4d, spec

ap = argparse.ArgumentParser()
ap.add_argument('group', nargs='?', default='buccal_tad'); ap.add_argument('step', nargs='?', type=int, default=1)
ap.add_argument('--only', choices=['febio', 'custom']); ap.add_argument('--out')
a = ap.parse_args()

S = Static(); state = State(S); arch = ArchForm(S); morph = Morpher(S, arch)
for u in S.keep:
    S.tooth[u]['centroid_cur'] = state.centroid(u)
load = Elastic(S, state, group=a.group)
shape = morph.shape(a.step)
pre = S.prescribed_nodes(state, arch, a.step, shape)
slave = np.unique(S.aligner_inner_faces)
Cp = build_pairs(pre, slave, S.master_tris, k=16, max_gap=solve4d.MAX_GAP)
gap_pre, _ = S.gap_for_pairing(pre, slave, Cp['tri'], Cp['bary'])
cref = dict(slave=slave, tri=Cp['tri'], bary=Cp['bary'], gap_pre=gap_pre)
state.aligner = morph.seat(shape, state, a.step)

res = {}
for name, fn in (('custom', lambda: solve4d.solve_config(S, state, load, cache={}, contact_ref=cref)),
                 ('febio', lambda: febio4d.solve_config(S, state, load, contact_ref=cref, verbose=True))):
    if a.only and a.only != name:
        continue
    t0 = time.perf_counter(); r = fn(); wall = time.perf_counter() - t0
    res[name] = dict(q=r['q'], wall=wall, hist=r['hist'][-1], eps99=float(np.percentile(vm_eq := solve4d.vm_equiv_strain(r['eps']), 99)))
    h = r['hist'][-1]
    print(f"{name}: {wall:.1f} s  iterations {h['it'] + (0 if name == 'febio' else 1)}  |Fc| {h['Fc']:.3f} N  "
          f"p99 PDL strain {res[name]['eps99']*100:.2f}%", flush=True)

def motion(q, u):
    T = S.tooth[u]; c = state.centroid(u)
    cp, rp = crown_point(S, u), root_point(S, u)
    dc = q[:3] + np.cross(q[3:], cp - c); dr = q[:3] + np.cross(q[3:], rp - c)
    return dict(crown_d=dc @ T['e_d']*1e3, crown_b=dc @ T['e_b']*1e3, crown_o=dc @ T['e_o']*1e3,
                root_d=dr @ T['e_d']*1e3, tip=np.degrees(q[3:] @ T['e_b']))

if len(res) == 2:
    print(f"\nstep {a.step}, {a.group}: tooth motion in um (custom / FEBio)")
    print(f"  {'FDI':>4} {'crown distal':>20} {'crown buccal':>20} {'crown occlusal':>20} {'root distal':>20} {'tip deg':>18}")
    for u in S.keep:
        mc, mf = motion(res['custom']['q'][u], u), motion(res['febio']['q'][u], u)
        print(f"  {spec.FDI[u]:>4} " + " ".join(f"{mc[k]:9.2f}/{mf[k]:9.2f}" if k != 'tip' else f"{mc[k]:8.4f}/{mf[k]:8.4f}"
                                              for k in ('crown_d', 'crown_b', 'crown_o', 'root_d', 'tip')))
if a.out:
    pickle.dump(res, open(a.out, 'wb'))
