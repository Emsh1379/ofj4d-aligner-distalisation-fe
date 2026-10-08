"""Which tooth does each group's precision cut (elastic hook) actually sit on?"""
import sys, os
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, 'ofj4d')); sys.path.insert(0, os.path.join(HERE, 'ofj'))
os.chdir(HERE)
from scipy.spatial import cKDTree
from model4d import Static, State
from loads4d import Elastic
import spec

S = Static(); st = State(S)
for u in S.keep:
    S.tooth[u]['centroid_cur'] = st.centroid(u)
ids = list(S.keep)
trees = {u: cKDTree(S.tooth[u]['surf_rest']) for u in ids}
can = S.tooth[6]
for side in ('buccal', 'palatal'):
    L = Elastic(S, st, group=f'{side}_tad')
    X = S.aligner_rest[L.hook_nodes]
    D = np.stack([trees[u].query(X)[0] for u in ids], axis=1)
    own = np.array(ids)[D.argmin(1)]
    counts = {spec.FDI[u]: int((own == u).sum()) for u in ids if (own == u).any()}
    c = X.mean(0) - can['crown_c']
    print(f"{side:<8} cut on (FDI: nodes) {counts}  | vs canine crown centre: "
          f"distal {c @ can['e_d']:+.2f}, buccal {c @ can['e_b']:+.2f}, occlusal {c @ can['e_o']:+.2f} mm"
          f" | mean distance to canine surface {D[:, ids.index(6)].mean():.2f} mm, "
          f"to lateral incisor {D[:, ids.index(7)].mean():.2f} mm")
