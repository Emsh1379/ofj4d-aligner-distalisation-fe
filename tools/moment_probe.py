"""What the appliance does to each tooth: net contact force and moment, split into
normal pressure and friction, crown and attachment, for the first PDL iteration of
one staging step (restarted from a checkpoint).

The readout is the height of the effective line of action below the crown point,
z = M_tip / F.  A force at the crown (z ~ 0) tips the tooth; a line of action
through the centre of resistance (z ~ 70% of crown-apex, FINDINGS 15) translates it.

usage: python tools/moment_probe.py RUN.pkl STEP
"""
import sys, os, pickle
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, 'ofj4d')); sys.path.insert(0, os.path.join(HERE, 'ofj'))
os.chdir(HERE)
from model4d import Static, State
from arch4d import ArchForm
from morph4d import Morpher
from loads4d import Elastic
from solve4d import solve_config, MAX_GAP, K_T
from contact import build_pairs
from run4d import crown_point, root_point
import spec


def main(run, step):
    S = Static(); st = State(S); arch = ArchForm(S); mo = Morpher(S, arch)
    D = pickle.load(open(run, 'rb')); g = D['group']
    if step > 1:
        st.T = {u: D['hist'][step-2]['T'][u].copy() for u in S.keep}
    for u in S.keep:
        S.tooth[u]['centroid_cur'] = st.centroid(u)
    load = Elastic(S, st, group=g)
    shape = mo.shape(step); slave = np.unique(S.aligner_inner_faces)
    pre = S.prescribed_nodes(st, arch, step, shape)
    Cp = build_pairs(pre, slave, S.master_tris, k=16, max_gap=MAX_GAP)
    gap_pre, _ = S.gap_for_pairing(pre, slave, Cp['tri'], Cp['bary'])
    cref = dict(slave=slave, tri=Cp['tri'], bary=Cp['bary'], gap_pre=gap_pre)
    st.aligner = mo.seat(shape, st, step)
    r = solve_config(S, st, load, contact_ref=cref, cache={})

    X = r['nodes']; C = r['C']; lam = r['lam']; n = r['nrm']
    us = r['u'].reshape(-1, 3); tri = C['tri']; w = C['bary']; s = C['slave']
    rel = us[s] - np.einsum('nj,njk->nk', w, us[tri])
    tang = rel - (rel*n).sum(1)[:, None]*n
    tm = np.linalg.norm(tang, axis=1)
    fmag = np.where(lam > 0, np.minimum(K_T*tm, spec.FRICTION*lam), 0.0)
    Ft = fmag[:, None]*tang/np.maximum(tm, 1e-15)[:, None]      # on the tooth
    Fn = -lam[:, None]*n                                         # on the tooth
    p = np.einsum('nj,njk->nk', w, X[tri])

    # which tooth, and crown or attachment
    G = pickle.load(open(os.environ.get('OFJ_GEOM', 'build4d/geom.pkl'), 'rb'))
    owner = np.full(len(X), -1); is_att = np.zeros(len(X), bool)
    for u in S.keep:
        T = S.tooth[u]; base = T['surf_off']
        owner[base:base+T['surf_n']] = u
        t = G['teeth'][u]
        if t['att'] is not None:
            sf = np.vstack([t['crown_faces'], t['att']['faces']])
            su = np.unique(sf)
            is_att[base + np.searchsorted(su, np.unique(t['att']['faces']))] = True
    own = owner[tri[:, 0]]; att = is_att[tri].all(1)

    print(f"\n{run}  step {step}  [{g}]  mu={spec.FRICTION}  K_T={K_T}  "
          f"iterations {len(r['hist'])}  active pairs {int((lam > 0).sum())}")
    print(f"{'tooth':<7}{'part':<11}{'F_dir N':>9}{'|F| N':>8}{'M_tip Nmm':>11}{'z mm':>8}"
          f"{'  (crown-apex mm, CRes ~70%)'}")
    for u in S.keep:
        T = S.tooth[u]; Tm = st.T[u]
        cp = crown_point(S, u) @ Tm[:3, :3].T + Tm[:3, 3]
        rp = root_point(S, u) @ Tm[:3, :3].T + Tm[:3, 3]
        L = np.linalg.norm(cp - rp); LA = (cp - rp)/L
        if u in (7, 8):
            d = -T['e_b']; lab = 'palatal'
        else:
            d = T['e_d']; lab = 'distal'
        d = d - (d @ LA)*LA; d /= np.linalg.norm(d)
        ax = np.cross(LA, d)
        k0 = own == u
        rows = [('normal', k0, Fn), ('friction', k0, Ft), ('total', k0, Fn + Ft),
                ('attachmt', k0 & att, Fn + Ft), ('crown', k0 & ~att, Fn + Ft)]
        for name, k, F in rows:
            if not k.any():
                continue
            Fs = F[k].sum(0); M = np.cross(p[k] - cp, F[k]).sum(0)
            Fd = Fs @ d; Mt = M @ ax
            z = -Mt/Fd if abs(Fd) > 1e-6 else np.nan
            tag = f"FDI {spec.FDI[u]}" if name == 'normal' else ''
            extra = f"   L={L:.1f}, {lab}" if name == 'normal' else ''
            print(f"{tag:<7}{name:<11}{Fd:9.3f}{np.linalg.norm(Fs):8.3f}{Mt:11.2f}{z:8.1f}{extra}")


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]))
