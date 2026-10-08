"""Where does a nonlinear iteration actually go?

cProfile cannot answer this: CUDA launches are asynchronous, so device work is
billed to whichever host call happens to synchronise next.  Every timer here
synchronises first.
"""
import sys, os, time
import numpy as np, torch
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE,'ofj4d')); sys.path.insert(0, os.path.join(HERE,'ofj')); os.chdir(HERE)
from model4d import Static, State
from arch4d import ArchForm
from morph4d import Morpher
from loads4d import Elastic
import solve4d
from solve4d import (Operators, contact_ops, K_N, K_T, MAX_GAP, direct_solver,
                     pdl_secant_modulus, vm_equiv_strain)
from contact import build_pairs
import gpu4d, spec

S = Static(); st = State(S); AF = ArchForm(S); M = Morpher(S, AF)
for u in S.keep: S.tooth[u]['centroid_cur'] = st.centroid(u)
load = Elastic(S, st, group='control')
slave = np.unique(S.aligner_inner_faces); shape = M.shape(1)
pre = S.prescribed_nodes(st, AF, 1, shape)
Cp = build_pairs(pre, slave, S.master_tris, k=16, max_gap=MAX_GAP)
gap_pre, _ = S.gap_for_pairing(pre, slave, Cp['tri'], Cp['bary'])
st.aligner = M.seat(shape, st, 1)
nodes = st.nodes(); ndof = 3*S.n_global
P = S.build_P(nodes).tocsr()
opA = Operators(nodes, S.aligner_tets, spec.NU_ALIGNER)
opP = Operators(nodes, S.pdl_tets, spec.PDL_NU)
gap_now, nrm = S.gap_for_pairing(nodes, slave, Cp['tri'], Cp['bary'])
dofs, R, a = contact_ops(dict(slave=slave, tri=Cp['tri'], bary=Cp['bary'],
                              normal=nrm, gap0=gap_now))
gref = gap_now - gap_pre
Pt_tan = np.eye(3)[None] - np.einsum('ni,nj->nij', nrm, nrm)
fr_ext = P.T @ load.vector(S, st, ndof)

g = gpu4d.for_step(Cp['tri'], ndof, dofs, opA.dofs, opP.dofs)
g.refresh(P, opA, opP, spec.E_ALIGNER, R, a, Pt_tan)
direct = direct_solver()
E_pdl = np.full(len(S.pdl_tets), spec.PDL_E1)
active = gref <= 0.0
ksec = np.full(len(gref), K_T)
lam = np.zeros(len(gref)); tmag = np.zeros(len(gref))

T = {}
def tic():
    torch.cuda.synchronize(); return time.time()
def toc(k, t0):
    torch.cuda.synchronize(); T[k] = T.get(k, 0.0) + time.time()-t0

N = 8
for it in range(N):
    t0 = tic(); As, bs, dg = g.system(E_pdl, active, ksec, K_N, gref, fr_ext); toc('GPU assembly -> host upper triangle', t0)
    t0 = tic(); xs, reused = direct.solve(As, bs, upper=True); toc('PARDISO solve (host)', t0)
    t0 = tic()
    r = bs - (As@xs + As.T@xs - As.diagonal()*xs)
    res = np.linalg.norm(r)/np.linalg.norm(bs)
    toc('residual check (host, 2 SpMV)', t0)
    t0 = tic(); x = xs/dg; u, eps_dev, E_dev = g.pdl_update(x, spec.PDL_E1, spec.PDL_E2, spec.PDL_EPS_KNEE); toc('PDL strain + modulus (GPU)', t0)
    t0 = tic()
    rel = np.einsum('nij,nj->ni', R, u[dofs])
    cn = gref + np.einsum('ni,ni->n', nrm, rel)
    lam_new = np.where(cn < 0, -K_N*cn, 0.0)
    lam = lam_new if it == 0 else 0.5*lam + 0.5*lam_new
    tang = rel - np.einsum('ni,ni->n', rel, nrm)[:, None]*nrm
    tmag = 0.5*tmag + 0.5*np.linalg.norm(tang, axis=1)
    active = cn < 0
    toc('contact kinematics (host)', t0)
    t0 = tic()
    E_new = E_dev.cpu().numpy()
    dE = np.abs(E_new-E_pdl).max()/spec.PDL_E2
    E_pdl = 0.5*E_pdl + 0.5*E_new
    ksec = np.clip(np.where(K_T*tmag <= spec.FRICTION*lam, K_T,
                            spec.FRICTION*lam/np.maximum(tmag, 1e-12)), 0, K_T)
    toc('modulus download + secant (host)', t0)

tot = sum(T.values())
print(f"\n  per nonlinear iteration, {N} iterations averaged\n")
for k, v in sorted(T.items(), key=lambda kv: -kv[1]):
    print(f"    {k:<38}{v/N*1e3:7.1f} ms   {100*v/tot:5.1f}%")
print(f"    {'TOTAL':<38}{tot/N*1e3:7.1f} ms")
