"""Is the 35-iteration cap costing accuracy at a busy step?

Replays one staging step from the recorded state at the end of the previous step,
first at the cap the run used and then at a much higher one, and compares the
tooth movements.  `hist[k]['T']` holds every tooth's rigid transform at the end of
step k+1, so the state going into step N is hist[N-2]['T'].
"""
import sys, os, pickle, time
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE,'ofj4d')); sys.path.insert(0, os.path.join(HERE,'ofj')); os.chdir(HERE)
from model4d import Static, State
from arch4d import ArchForm
from morph4d import Morpher
from loads4d import Elastic
from solve4d import solve_config, MAX_GAP
from contact import build_pairs
import run4d
import spec

STEP = int(sys.argv[1]) if len(sys.argv) > 1 else 30
GROUP = sys.argv[2] if len(sys.argv) > 2 else 'buccal_tad'
CAPS = tuple(int(x) for x in (sys.argv[3].split(',') if len(sys.argv)>3 else ('35','250')))

D = pickle.load(open(f'build4d/run_{GROUP}.pkl', 'rb'))
T_before = D['hist'][STEP-2]['T']
recorded = D['hist'][STEP-1]

S = Static()
arch = ArchForm(S)
morph = Morpher(S, arch)
slave = np.unique(S.aligner_inner_faces)
ref = dict(crown={u: run4d.crown_point(S, u) for u in S.keep},
           root={u: run4d.root_point(S, u) for u in S.keep})
ref['LA'] = {}
for u in S.keep:
    v = ref['crown'][u] - ref['root'][u]
    ref['LA'][u] = v/np.linalg.norm(v)

print(f"replaying step {STEP} of {GROUP}; the run recorded "
      f"{[i['nit'] for i in recorded['info']]} iterations "
      f"({'capped' if any(not i['converged'] for i in recorded['info']) else 'converged'})")

out = {}
for cap in CAPS:
    state = State(S)
    for u in S.keep:
        state.T[u] = T_before[u].copy()
        S.tooth[u]['centroid_cur'] = state.centroid(u)
    load = Elastic(S, state, group=GROUP)
    shape = morph.shape(STEP)
    pre = S.prescribed_nodes(state, arch, STEP, shape)
    Cp = build_pairs(pre, slave, S.master_tris, k=16, max_gap=MAX_GAP)
    gap_pre, _ = S.gap_for_pairing(pre, slave, Cp['tri'], Cp['bary'])
    cref = dict(slave=slave, tri=Cp['tri'], bary=Cp['bary'], gap_pre=gap_pre)
    cache, E_prev, nits, fcs, dus = {}, None, [], [], []
    t0 = time.time()
    for k in range(spec.PDL_ITERATIONS_PER_STEP):
        state.aligner = morph.seat(shape, state, STEP)
        r = solve_config(S, state, load, E_pdl0=E_prev, x0=cache.get('x0'),
                         cache=cache, contact_ref=cref, maxit=cap)
        E_prev = r['E_pdl']
        for u in S.keep:
            q = r['q'][u]
            state.apply_rigid(u, q[:3], q[3:])
        nits.append(len(r['hist']))
        fcs.append(r['hist'][-1]['Fc']); dus.append(r['hist'][-1]['du'])
    out[cap] = run4d.measure(S, state, ref)
    print(f"  maxit={cap:3d}: {nits} iterations, {time.time()-t0:.0f} s, "
          f"Fc={[f'{f:.1f}' for f in fcs]} N, du={[f'{d:.1e}' for d in dus]}")

keys = ('crown_distal', 'crown_buccal', 'crown_occlusal',
        'root_distal', 'root_buccal', 'root_occlusal')
a, b = out[CAPS[0]], out[CAPS[1]]
worst = max(abs(a[u][k]-b[u][k]) for u in a for k in keys)
scale = max(abs(a[u][k]) for u in a for k in keys)
print(f"\n  cap {CAPS[0]} vs {CAPS[1]}: largest movement difference "
      f"{worst*1e6:.4f} nm, on movements up to {scale*1e3:.1f} um "
      f"({worst/scale:.2e} relative)")

# and does the replay reproduce what the run recorded?
rec = recorded['meas']
drift = max(abs(a[u][k]-rec[u][k]) for u in a for k in keys)
print(f"  replay vs the recorded run at maxit=35: {drift*1e6:.4f} nm")
print("\n  VERDICT: the cap is not costing accuracy"
      if worst < 1e-5 else "\n  VERDICT: the cap IS affecting the result -- rerun with a higher cap")
