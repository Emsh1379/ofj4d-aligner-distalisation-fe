"""PDL stress at chosen staging steps of a finished run (replayed from its stored tooth positions).

usage: [same OFJ_* environment as the run]  python tools/pdl_stress.py RUN.pkl [--steps 10 20 40 70] [--out PREFIX]

run4d stores each step's tooth transforms, not the fields. For each requested step s this replays the step
exactly as run4d does -- teeth at the end of step s-1, aligner made for step s, wear-in seat, contact reference,
the two PDL (bone-remodelling) iterations -- and evaluates the PDL stress of the FIRST iteration: the
appliance of step s activated against the teeth, i.e. the loaded state that drives that step's movement.
The replay is checked against the stored run: tooth movement after both iterations vs the recorded step.

Stress = D(E_s, nu) eps per element, with E_s the bilinear PDL secant modulus at convergence (solve4d), in kPa.
Hydrostatic = mean normal stress (negative = compression). Per tooth: 5th / 95th percentile hydrostatic,
95th percentile von Mises, over the whole PDL and in its cervical and apical thirds (by depth along the
tooth axis between crown point and root point). Writes PREFIX.json (summary) and PREFIX_step<s>.npz
(element centroids, owner tooth, hydrostatic, von Mises) for figures.
"""
import sys, os, json, pickle, argparse, time
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, 'ofj4d')); sys.path.insert(0, os.path.join(HERE, 'ofj')); os.chdir(HERE)
from model4d import Static, State
from arch4d import ArchForm
from morph4d import Morpher
from loads4d import Elastic
from contact import build_pairs
from fem import isotropic_D
import solve4d, run4d, spec

ap = argparse.ArgumentParser()
ap.add_argument('run'); ap.add_argument('--steps', type=int, nargs='+', default=[10, 20, 40, 70])
ap.add_argument('--out')
a = ap.parse_args()
D = pickle.load(open(a.run, 'rb'))
out_prefix = a.out or a.run[:-4] + '_pdlstress'
S = Static(); arch = ArchForm(S); morph = Morpher(S, arch)
slave = np.unique(S.aligner_inner_faces)
ref = dict(crown={u: run4d.crown_point(S, u) for u in S.keep}, root={u: run4d.root_point(S, u) for u in S.keep})
ref['LA'] = {}
for u in S.keep:
    v = ref['crown'][u] - ref['root'][u]; ref['LA'][u] = v/np.linalg.norm(v)
hist = {h['step']: h for h in D['hist']}
Dhat = isotropic_D(1.0, spec.PDL_NU)
summary = dict(run=a.run, group=D['group'], nu=spec.PDL_NU, steps={})

for s in a.steps:
    if s not in hist:
        print(f"step {s}: not in the run"); continue
    t0 = time.time()
    state = State(S)
    if s >= 2:
        for u in S.keep:
            state.T[u] = hist[s - 1]['T'][u].copy()
    for u in S.keep:
        S.tooth[u]['centroid_cur'] = state.centroid(u)
    load = Elastic(S, state, group=D['group'])
    shape = morph.shape(s)
    pre = S.prescribed_nodes(state, arch, s, shape)
    Cp = build_pairs(pre, slave, S.master_tris, k=16, max_gap=solve4d.MAX_GAP)
    gap_pre, _ = S.gap_for_pairing(pre, slave, Cp['tri'], Cp['bary'])
    cref = dict(slave=slave, tri=Cp['tri'], bary=Cp['bary'], gap_pre=gap_pre)
    cache, E_prev, first = {}, None, None
    for k in range(spec.PDL_ITERATIONS_PER_STEP):
        state.aligner = morph.seat(shape, state, s)
        r = solve4d.solve_config(S, state, load, E_pdl0=E_prev, x0=cache.get('x0'), cache=cache, contact_ref=cref)
        E_prev = r['E_pdl']
        if k == 0:
            nodes = r['nodes']
            sig = np.einsum('kl,nl->nk', Dhat, r['eps'])*r['E_pdl'][:, None]*1e3          # kPa
            hyd = sig[:, :3].mean(1)
            sx, sy, sz, txy, tyz, tzx = sig.T
            vm = np.sqrt(0.5*((sx-sy)**2 + (sy-sz)**2 + (sz-sx)**2) + 3*(txy**2 + tyz**2 + tzx**2))
            cen = nodes[S.pdl_tets].mean(1)
            first = dict(hyd=hyd, vm=vm, cen=cen, Fc=r['hist'][-1]['Fc'], it=len(r['hist']))
        for u in S.keep:
            q = r['q'][u]; state.apply_rigid(u, q[:3], q[3:])
    meas = run4d.measure(S, state, ref)
    rec = hist[s]['meas']
    dev = max(abs(meas[u][k] - rec[u][k]) for u in S.keep
              for k in ('crown_distal', 'crown_buccal', 'crown_occlusal', 'root_distal', 'root_buccal', 'root_occlusal'))
    teeth = {}
    for u in S.keep:
        m = S.pdl_owner == u; T = S.tooth[u]
        cp = ref['crown'][u] @ state.T[u][:3, :3].T + state.T[u][:3, 3]
        rp = ref['root'][u] @ state.T[u][:3, :3].T + state.T[u][:3, 3]
        ax = (cp - rp)/np.linalg.norm(cp - rp)
        depth = ((cp - first['cen'][m]) @ ax)/np.linalg.norm(cp - rp)          # 0 crown point .. 1 root point
        dmin, dmax = depth.min(), depth.max(); rel = (depth - dmin)/max(dmax - dmin, 1e-9)
        def stats(sel):
            h_, v_ = first['hyd'][m][sel], first['vm'][m][sel]
            return dict(hyd_p5=float(np.percentile(h_, 5)), hyd_p95=float(np.percentile(h_, 95)), vm_p95=float(np.percentile(v_, 95)))
        teeth[int(spec.FDI[u])] = dict(all=stats(slice(None)), cervical=stats(rel <= 1/3), apical=stats(rel >= 2/3))
    summary['steps'][s] = dict(first_iteration_Fc=first['Fc'], first_iteration_it=first['it'],
                               replay_vs_recorded_max_mm=float(dev), teeth=teeth)
    np.savez_compressed(f"{out_prefix}_step{s}.npz", centroid=first['cen'], owner=np.array([spec.FDI[u] for u in S.pdl_owner]),
                        hydrostatic_kPa=first['hyd'], von_mises_kPa=first['vm'])
    print(f"step {s}: replay vs recorded run max {dev*1e3:.3f} um; first iteration |Fc| {first['Fc']:.2f} N; {time.time()-t0:.0f} s")
    for fdi, t in teeth.items():
        print(f"   FDI {fdi}: hydrostatic p5/p95 {t['all']['hyd_p5']:7.2f}/{t['all']['hyd_p95']:6.2f} kPa "
              f"(cervical {t['cervical']['hyd_p5']:7.2f}, apical {t['apical']['hyd_p5']:7.2f}), von Mises p95 {t['all']['vm_p95']:6.2f} kPa")
json.dump(summary, open(out_prefix + '.json', 'w'), indent=1)
print('->', out_prefix + '.json')
