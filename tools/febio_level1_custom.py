"""Level-1 cross-validation, custom side: PDL-only system of the half arch, solved
with this project's own operators (ofj4d/solve4d.Operators, model4d projection,
MKL PARDISO), and exported for the FEBio side (tools/febio_level1_export.py).

Problem (docs/FINDINGS.md section 28 / pipeline/febio_cudss/README.md, Level 1):
  every tooth rigid (6 dofs about its centroid), PDL linear elastic E1 = 0.05 MPa,
  nu = spec.PDL_NU, outer PDL layer fixed, one force of F N applied at the crown
  point (the most occlusal crown vertex) of every tooth at once.  Load cases: distal
  (e_d), buccal (e_b), intrusive (-e_o).  The teeth are uncoupled (separate PDLs, fixed
  outer layer), so one solve per load case covers all seven.

usage: python tools/febio_level1_custom.py [--repeat N] [--out build4d/febio_val/custom.pkl]
       MKL_NUM_THREADS sets the PARDISO thread count.
"""
import sys, os, time, pickle, argparse
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, 'ofj4d')); sys.path.insert(0, os.path.join(HERE, 'ofj')); os.chdir(HERE)
from model4d import Static, State
from solve4d import Operators, Pardiso
import spec

ap = argparse.ArgumentParser()
ap.add_argument('--repeat', type=int, default=3)
ap.add_argument('--force', type=float, default=1.0)
ap.add_argument('--out', default='build4d/febio_val/custom.pkl')
a = ap.parse_args()

S = Static(); st = State(S)
for u in S.keep:
    S.tooth[u]['centroid_cur'] = st.centroid(u)
nodes = st.nodes(); ndof = 3*S.n_global

CASES = {'distal': lambda T: T['e_d'], 'buccal': lambda T: T['e_b'], 'intrusive': lambda T: -T['e_o']}

teeth = {}
for u in S.keep:
    T = S.tooth[u]
    crown = T['surf_rest'][:T['surf_n']]
    L = T['pdl_layers']
    teeth[u] = dict(fdi=T['fdi'], centroid=T['centroid_cur'].copy(),
                    e_d=T['e_d'], e_b=T['e_b'], e_o=T['e_o'],
                    crown_pt=crown[np.argmax(crown @ T['e_o'])],
                    apex=T['pdl_rest'][np.argmin(T['pdl_rest'] @ T['e_o'])],
                    pdl_nodes=T['pdl_rest'].copy(), pdl_layers=[np.asarray(l) for l in L],
                    pdl_off=T['pdl_off'])
pdl_tets_local = {u: S.pdl_tets[S.pdl_owner == u] - S.tooth[u]['pdl_off'] for u in S.keep}

timings = []
for rep in range(a.repeat):
    t0 = time.perf_counter()
    P = S.build_P(nodes).tocsr()
    tP = time.perf_counter()
    opP = Operators(nodes, S.pdl_tets, spec.PDL_NU)
    tO = time.perf_counter()
    Kg = opP.K(np.full(len(S.pdl_tets), spec.PDL_E1), ndof)
    tK = time.perf_counter()
    K = (P.T @ Kg @ P).tocsr(); del Kg
    tR = time.perf_counter()
    keep = np.where(np.abs(K.diagonal()) > 0)[0]          # dofs the PDL touches
    Kk = K[keep][:, keep].tocsr()
    t1 = time.perf_counter()
    B = np.zeros((Kk.shape[0], len(CASES)))
    col = {c: j for j, c in enumerate(keep)}
    for j, (name, dirf) in enumerate(CASES.items()):
        for u in S.keep:
            Fv = a.force*dirf(S.tooth[u]); tt = teeth[u]
            M = np.cross(tt['crown_pt'] - tt['centroid'], Fv)
            q0 = S.q_col[u]
            for k in range(3):
                B[col[q0+k], j] += Fv[k]; B[col[q0+3+k], j] += M[k]
    solver = Pardiso()                                      # the pipeline's own MKL wrapper
    ts = time.perf_counter()
    x1 = solver.solve(Kk, B[:, 0])[0]                       # analysis + factorisation + 1 solve
    ts1 = time.perf_counter()
    X = np.column_stack([x1] + [solver.solve(Kk, B[:, j])[0] for j in range(1, B.shape[1])])
    t2 = time.perf_counter()
    solver.free()
    tm = dict(build_P=tP-t0, element_ops=tO-tP, assemble_K=tK-tO, project_PtKP=tR-tK, extract=t1-tR,
              solve_first=ts1-ts, solve_extra_rhs=(t2-ts1)/max(B.shape[1]-1, 1),
              assemble=t1-t0, solve=t2-t1, total=t2-t0, to_first_solution=ts1-t0)
    timings.append(tm)
    print(f"rep {rep}: " + "  ".join(f"{k} {v:.3f}" for k, v in tm.items()), flush=True)
    print("TIMING " + __import__('json').dumps(tm), flush=True)

res = {}
for j, name in enumerate(CASES):
    x = np.zeros(K.shape[0]); x[keep] = X[:, j]
    ug = (P @ x).reshape(-1, 3)
    r = np.linalg.norm(Kk @ X[:, j] - B[:, j])/np.linalg.norm(B[:, j])
    res[name] = dict(residual=r, teeth={u: dict(q=x[S.q_col[u]:S.q_col[u]+6].copy(),
                     u_pdl=ug[S.tooth[u]['pdl_off']:S.tooth[u]['pdl_off']+S.tooth[u]['pdl_n']].copy())
                     for u in S.keep})
    print(f"{name}: relative residual {r:.1e}")
info = dict(n_red=int(Kk.shape[0]), nnz=int(Kk.nnz), n_pdl_tets=len(S.pdl_tets),
            threads=os.environ.get('MKL_NUM_THREADS', 'default'))
print(info)
pickle.dump(dict(teeth=teeth, tets=pdl_tets_local, cases=list(CASES), force=a.force, E=spec.PDL_E1,
                 nu=spec.PDL_NU, res=res, timings=timings, info=info), open(a.out, 'wb'))
print('->', a.out)
