"""One solve of staging step 1 -- to verify the deformable-bone variant.

Environment decides the model:  OFJ_BONE=rigid|deformable,  OFJ_BONE_E_SCALE.
    python tools/bone_check.py <group> <out.pkl>
"""
import sys, os, time, pickle
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, 'ofj4d')); sys.path.insert(0, os.path.join(HERE, 'ofj'))
os.chdir(HERE)
from model4d import Static, State
from arch4d import ArchForm
from morph4d import Morpher
from loads4d import Elastic
from solve4d import solve_config, MAX_GAP
from contact import build_pairs

GROUP, OUT = sys.argv[1], sys.argv[2]
S = Static(); st = State(S); arch = ArchForm(S); morph = Morpher(S, arch)
for u in S.keep:
    S.tooth[u]['centroid_cur'] = st.centroid(u)
load = Elastic(S, st, group=GROUP)
slave = np.unique(S.aligner_inner_faces)
shape = morph.shape(1)
pre = S.prescribed_nodes(st, arch, 1, shape)
Cp = build_pairs(pre, slave, S.master_tris, k=16, max_gap=MAX_GAP)
gp, _ = S.gap_for_pairing(pre, slave, Cp['tri'], Cp['bary'])
st.aligner = morph.seat(shape, st, 1)
t0 = time.time()
r = solve_config(S, st, load, verbose=True, cache={},
                 contact_ref=dict(slave=slave, tri=Cp['tri'], bary=Cp['bary'], gap_pre=gp))
wall = time.time() - t0
out = dict(q=r['q'], hist=r['hist'], wall=wall, n_red=S.n_red,
           deformable=getattr(S, 'deformable_bone', False),
           tie_inside=getattr(S, 'tie_inside', None), n_bone_react=len(load.bone_nodes)
           if hasattr(load, 'bone_nodes') else 0)
if out['deformable']:
    ub = r['u'][3*S.bone_off:3*(S.bone_off + S.nB)].reshape(-1, 3)
    out['bone_umax'] = float(np.linalg.norm(ub, axis=1).max())
pickle.dump(out, open(OUT, 'wb'))
print(f"RESULT n_red {S.n_red}  iterations {len(r['hist'])}  wall {wall:.0f} s  "
      f"tie_inside {out['tie_inside']}  bone_umax {out.get('bone_umax')}  "
      f"bone reaction nodes {out['n_bone_react']}")
