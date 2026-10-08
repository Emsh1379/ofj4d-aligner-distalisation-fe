"""Level-1 cross-validation: compare FEBio runs (tools/febio_level1_export.py) with the
custom solution (tools/febio_level1_custom.py).

For each tooth: the rigid motion (translation t at the centroid, rotation vector w)
is recovered from FEBio's inner-PDL-layer node displacements by an exact rigid
(Kabsch) fit, which is valid at any rotation size; crown-point and apex displacement,
rotation and centre of rotation follow from it.  The PDL displacement field is
compared node by node.  FEBio results are scaled by F_custom / F_febio, so a small-load
FEBio run compares against the linear custom solution directly.

usage: python tools/febio_level1_compare.py RUN.feb [RUN.feb ...]  [--json out.json]
"""
import sys, os, re, pickle, json, argparse
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ap = argparse.ArgumentParser()
ap.add_argument('runs', nargs='+')
ap.add_argument('--custom', default=os.path.join(HERE, 'build4d/febio_val/custom.pkl'))
ap.add_argument('--json')
a = ap.parse_args()
D = pickle.load(open(a.custom, 'rb'))


def last_block(path):
    blocks, cur = [], None
    for line in open(path):
        if line.startswith('*'):
            if line.startswith('*Time'):
                cur = {}; blocks.append((float(line.split('=')[1]), cur))
            continue
        v = line.strip().split(',')
        if cur is not None and len(v) >= 4 and v[0].strip().isdigit():
            cur[int(v[0])] = [float(x) for x in v[1:4]]
    return blocks[-1]


def kabsch(X, Y):
    """Rigid map Y ~ R (X - cx) + cy."""
    cx, cy = X.mean(0), Y.mean(0)
    H = (X - cx).T @ (Y - cy)
    U, _, Vt = np.linalg.svd(H)
    s = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, s]) @ U.T
    return R, cx, cy


def rotvec(R):
    th = np.arccos(np.clip((np.trace(R) - 1)/2, -1, 1))
    if th < 1e-14:
        return np.zeros(3)
    return th/(2*np.sin(th))*np.array([R[2, 1]-R[1, 2], R[0, 2]-R[2, 0], R[1, 0]-R[0, 1]])


def febio_log(path):
    out = {}
    if not os.path.exists(path):
        return out
    txt = open(path, errors='ignore').read()
    for key, pat in [('equilibrium_iters', r'Total number of equilibrium iterations\s*\.*\s*:\s*(\d+)'),
                     ('reformations', r'Total number of stiffness reformations\s*\.*\s*:\s*(\d+)'),
                     ('rhs_evals', r'Total number of right hand evaluations\s*\.*\s*:\s*(\d+)'),
                     ('equations', r'Nr of equations\s*\.*\s*:\s*(\d+)'),
                     ('nnz', r'Nr of nonzeroes in stiffness matrix\s*\.*\s*:\s*(\d+)')]:
        m = re.search(pat, txt)
        if m: out[key] = int(m.group(1))
    sec = lambda label: r'\s' + label + r'\s*\.*\s*:\s*[\d:]+\s*\(([\d.eE+-]+) sec\)'
    for key, label in [('total_s', 'Total elapsed time'), ('input_s', 'Input time'), ('solve_s', 'Solve time'),
                       ('linear_solver_s', 'time in linear solver'), ('stiffness_s', 'evaluating stiffness'),
                       ('reform_s', 'reforming stiffness'), ('residual_s', 'evaluating residual'),
                       ('update_s', 'model update'), ('io_s', r'IO-time \(plot, dmp, data\)')]:
        m = re.search(sec(label), txt)
        if m: out[key] = float(m.group(1))
    out['normal_termination'] = 'N O R M A L   T E R M I N A T I O N' in txt
    return out


summary = []
for run in a.runs:
    M = pickle.load(open(run + '.map.pkl', 'rb'))
    case, Ff = M['case'], M['force']
    scale = D['force']/Ff
    t_end, U = last_block(run[:-4] + '_nodes.txt')
    ref = D['res'][case]['teeth']
    rows = []; num = den = 0.0
    for u, T in D['teeth'].items():
        first = M['idmap'][u]; X = T['pdl_nodes']; n = len(X)
        Uf = np.array([U[first+i] for i in range(n)])
        Uc = ref[u]['u_pdl']
        num += np.sum((Uf*scale - Uc)**2); den += np.sum(Uc**2)
        inner = T['pdl_layers'][0]; c = T['centroid']
        R, cx, cy = kabsch(X[inner], X[inner] + Uf[inner])
        w_f = rotvec(R); t_f = (R @ (c - cx) + cy) - c
        q = ref[u]['q']; t_c, w_c = q[:3], q[3:]
        def disp_f(p): return (R @ (p - cx) + cy) - p
        def disp_c(p): return t_c + np.cross(w_c, p - c)
        ed = {'distal': T['e_d'], 'buccal': T['e_b'], 'intrusive': -T['e_o']}[case]
        ax = np.cross(T['e_o'], ed); nax = np.linalg.norm(ax)
        rot_axis = ax/nax if nax > 1e-9 else T['e_o']
        def cor(t, w):          # centre of rotation on the long axis, distance below the crown point
            aa = np.cross(w, T['e_o']); b = t + np.cross(w, T['crown_pt'] - c)
            return float((aa @ b)/(aa @ aa)) if aa @ aa > 1e-30 else float('nan')
        L_ca = float((T['crown_pt'] - T['apex']) @ T['e_o'])
        rows.append(dict(fdi=T['fdi'],
                   crown_c=float(disp_c(T['crown_pt']) @ ed*1e3), crown_f=float(disp_f(T['crown_pt']) @ ed*scale*1e3),
                   apex_c=float(disp_c(T['apex']) @ ed*1e3), apex_f=float(disp_f(T['apex']) @ ed*scale*1e3),
                   rot_c=float(np.degrees(w_c @ rot_axis)), rot_f=float(np.degrees(w_f @ rot_axis)*scale),
                   cor_c=cor(t_c, w_c)/L_ca*100, cor_f=cor(t_f, w_f)/L_ca*100,
                   field_rel_err=float(np.linalg.norm(Uf*scale - Uc)/np.linalg.norm(Uc)),
                   max_node_diff_um=float(np.abs(Uf*scale - Uc).max()*1e3)))
    log = febio_log(run[:-4] + '.log')
    rec = dict(run=os.path.basename(run), case=case, force=Ff, steps=M['steps'], pdl_mat=M['pdl_mat'],
               t_end=t_end, field_rel_err=float(np.sqrt(num/den)), teeth=rows, febio=log)
    summary.append(rec)
    print(f"\n=== {rec['run']}: {case}, F = {Ff:g} N ({M['pdl_mat']}, {M['steps']} steps), "
          f"results x{scale:g} -> per 1 N;  FEBio {log}")
    print(f"  whole PDL displacement field, relative L2 difference: {rec['field_rel_err']:.3e}")
    print(f"  {'FDI':>4} {'crown um c/F':>20} {'d%':>7} {'apex um c/F':>20} {'d%':>7} {'rot deg c/F':>18} {'d%':>7}"
          f" {'CoR % c/F':>14} {'field err':>10}")
    pc = lambda x, y: 100*(y - x)/x if abs(x) > 1e-12 else float('nan')
    for r in rows:
        print(f"  {r['fdi']:>4} {r['crown_c']:9.2f}/{r['crown_f']:9.2f} {pc(r['crown_c'], r['crown_f']):+7.2f}"
              f" {r['apex_c']:9.2f}/{r['apex_f']:9.2f} {pc(r['apex_c'], r['apex_f']):+7.2f}"
              f" {r['rot_c']:8.4f}/{r['rot_f']:8.4f} {pc(r['rot_c'], r['rot_f']):+7.2f}"
              f" {r['cor_c']:6.1f}/{r['cor_f']:6.1f} {r['field_rel_err']:10.2e}")
if a.json:
    json.dump(summary, open(a.json, 'w'), indent=1)
    print('->', a.json)
