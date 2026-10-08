"""Save the scaled linear system of one nonlinear iteration (staging step 1), exactly
as the solver receives it, for benchmarking linear solvers offline.

    OFJ_BONE=rigid|deformable python tools/dump_system.py <group> <out-prefix>
writes <out-prefix>_U.npz (upper triangle, CSR) and <out-prefix>_b.npy.
"""
import sys, os
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, 'ofj4d')); sys.path.insert(0, os.path.join(HERE, 'ofj'))
os.chdir(HERE)
from scipy.sparse import triu, save_npz
from model4d import Static, State
from arch4d import ArchForm
from morph4d import Morpher
from loads4d import Elastic
from contact import build_pairs
import solve4d

GROUP, OUT = sys.argv[1], sys.argv[2]


def capture(self, As, bs, upper=False):
    U = As if upper else triu(As, format='csr')
    save_npz(OUT + '_U.npz', U.tocsr())
    np.save(OUT + '_b.npy', np.asarray(bs))
    print(f"saved system: n {U.shape[0]}, upper nnz {U.nnz}")
    raise SystemExit(0)          # BaseException: not swallowed by the solver fallback


solve4d.Pardiso.solve = capture
S = Static(); st = State(S); arch = ArchForm(S); morph = Morpher(S, arch)
for u in S.keep:
    S.tooth[u]['centroid_cur'] = st.centroid(u)
load = Elastic(S, st, group=GROUP)
slave = np.unique(S.aligner_inner_faces)
shape = morph.shape(1)
pre = S.prescribed_nodes(st, arch, 1, shape)
Cp = build_pairs(pre, slave, S.master_tris, k=16, max_gap=solve4d.MAX_GAP)
gp, _ = S.gap_for_pairing(pre, slave, Cp['tri'], Cp['bary'])
st.aligner = morph.seat(shape, st, 1)
solve4d.solve_config(S, st, load, cache={}, maxit=1,
                     contact_ref=dict(slave=slave, tri=Cp['tri'], bary=Cp['bary'], gap_pre=gp))
