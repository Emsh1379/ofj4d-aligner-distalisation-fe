"""Level-1 cross-validation, FEBio side: write the PDL-only half-arch problem solved by
tools/febio_level1_custom.py as a FEBio 4 model.

Same mesh (this project's 3-layer PDL, node for node), same material constants, same
loads: each tooth is a FEBio rigid body with its centre of mass at the tooth centroid,
the inner PDL layer is tied to it (`rigid` node-set BC), the outer PDL layer is fixed,
and the crown-point force is applied as rigid_force + rigid_moment = (p - c) x F about
the centroid -- exactly the load vector of the custom solve.

FEBio has no small-strain solver: `isotropic elastic` is St Venant-Kirchhoff in the
nonlinear `solid` solver.  Use a small --force (e.g. 1e-3 N) to compare against the
linear custom solution, the full 1 N to measure the effect of finite deformation.

usage: python tools/febio_level1_export.py CASE FORCE_N TIME_STEPS OUT.feb [--pdl-mat isotropic|neo]
       CASE = distal | buccal | intrusive
"""
import sys, os, pickle, argparse
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ap = argparse.ArgumentParser()
ap.add_argument('case'); ap.add_argument('force', type=float)
ap.add_argument('steps', type=int); ap.add_argument('out')
ap.add_argument('--custom', default=os.path.join(HERE, 'build4d/febio_val/custom.pkl'))
ap.add_argument('--pdl-mat', default='isotropic', choices=['isotropic', 'neo'])
ap.add_argument('--auto', action='store_true', help='FEBio adaptive time stepper (large loads)')
a = ap.parse_args()

D = pickle.load(open(a.custom, 'rb'))
DIR = {'distal': lambda T: T['e_d'], 'buccal': lambda T: T['e_b'], 'intrusive': lambda T: -T['e_o']}[a.case]
teeth = D['teeth']; keep = list(teeth)
f = lambda v: ','.join(f'{x:.12g}' for x in v)

nodes, elems, sets, idmap = [], [], {}, {}
nid = 1; eid = 1
outer = []
for u in keep:
    T = teeth[u]
    X = T['pdl_nodes']; n = len(X)
    first = nid
    nodes.append((first, X)); idmap[u] = first
    tets = D['tets'][u].copy()
    P = X[tets]
    vol = np.einsum('ij,ij->i', np.cross(P[:, 1]-P[:, 0], P[:, 2]-P[:, 0]), P[:, 3]-P[:, 0])
    tets[vol < 0] = tets[vol < 0][:, [0, 1, 3, 2]]
    elems.append((f"PDL{T['fdi']}", eid, tets + first)); eid += len(tets)
    sets[f"inner{T['fdi']}"] = T['pdl_layers'][0] + first
    outer.append(T['pdl_layers'][-1] + first)
    nid += n
    # a small free-standing tet that carries the rigid material (0.5 mm, at the centroid)
    c = T['centroid']
    R = c + 0.5*np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    nodes.append((nid, R))
    elems.append((f"RB{T['fdi']}", eid, np.array([[nid, nid+1, nid+2, nid+3]]))); eid += 1
    nid += 4
sets['outer'] = np.concatenate(outer)

L = ['<?xml version="1.0" encoding="ISO-8859-1"?>', '<febio_spec version="4.0">',
     '\t<Module type="solid"/>',
     '\t<Control>', '\t\t<analysis>STATIC</analysis>',
     f'\t\t<time_steps>{a.steps}</time_steps>', f'\t\t<step_size>{1.0/a.steps:.12g}</step_size>',
     '\t\t<plot_level>PLOT_NEVER</plot_level>', '\t\t<output_level>OUTPUT_FINAL</output_level>',
     *(['\t\t<time_stepper type="default">', '\t\t\t<max_retries>25</max_retries>',
        '\t\t\t<opt_iter>8</opt_iter>', f'\t\t\t<dtmin>{1e-4/a.steps:.6g}</dtmin>',
        f'\t\t\t<dtmax>{1.0/a.steps:.6g}</dtmax>', '\t\t</time_stepper>'] if a.auto else []),
     '\t\t<solver type="solid">',
     '\t\t\t<dtol>1e-9</dtol>', '\t\t\t<etol>1e-14</etol>', '\t\t\t<rtol>0</rtol>',
     '\t\t\t<lstol>0.9</lstol>', '\t\t\t<max_refs>50</max_refs>',
     '\t\t\t<qn_method type="BFGS">', '\t\t\t\t<max_ups>0</max_ups>', '\t\t\t</qn_method>',
     '\t\t</solver>', '\t</Control>', '\t<Material>']
if a.pdl_mat == 'isotropic':
    L += ['\t\t<material id="1" name="PDL" type="isotropic elastic">']
else:
    L += ['\t\t<material id="1" name="PDL" type="neo-Hookean">']
L += ['\t\t\t<density>1</density>', f"\t\t\t<E>{D['E']}</E>", f"\t\t\t<v>{D['nu']}</v>", '\t\t</material>']
for k, u in enumerate(keep):
    T = teeth[u]
    L += [f'\t\t<material id="{k+2}" name="rb{T["fdi"]}" type="rigid body">', '\t\t\t<density>1</density>',
          f'\t\t\t<center_of_mass>{f(T["centroid"])}</center_of_mass>', '\t\t</material>']
L += ['\t</Material>', '\t<Mesh>', '\t\t<Nodes name="all">']
for first, X in nodes:
    L += [f'\t\t\t<node id="{first+i}">{f(x)}</node>' for i, x in enumerate(X)]
L += ['\t\t</Nodes>']
for name, e0, E in elems:
    L.append(f'\t\t<Elements type="tet4" name="{name}">')
    L += [f'\t\t\t<elem id="{e0+i}">{",".join(map(str, e))}</elem>' for i, e in enumerate(E)]
    L.append('\t\t</Elements>')
for name, ids in sets.items():
    L.append(f'\t\t<NodeSet name="{name}">{",".join(map(str, ids))}</NodeSet>')
L += ['\t</Mesh>', '\t<MeshDomains>']
for name, _, _ in elems:
    mat = 'PDL' if name.startswith('PDL') else 'rb' + name[2:]
    L.append(f'\t\t<SolidDomain name="{name}" mat="{mat}"/>')
L += ['\t</MeshDomains>', '\t<Boundary>',
      '\t\t<bc name="fixed_outer" node_set="outer" type="zero displacement">',
      '\t\t\t<x_dof>1</x_dof>', '\t\t\t<y_dof>1</y_dof>', '\t\t\t<z_dof>1</z_dof>', '\t\t</bc>']
for u in keep:
    fdi = teeth[u]['fdi']
    L += [f'\t\t<bc name="tie{fdi}" node_set="inner{fdi}" type="rigid">', f'\t\t\t<rb>rb{fdi}</rb>', '\t\t</bc>']
L += ['\t</Boundary>', '\t<Rigid>']
for u in keep:
    T = teeth[u]; fdi = T['fdi']
    F = a.force*DIR(T); M = np.cross(T['crown_pt'] - T['centroid'], F)
    for k, ax in enumerate('xyz'):
        L += [f'\t\t<rigid_load name="F{fdi}{ax}" type="rigid_force">', f'\t\t\t<rb>rb{fdi}</rb>',
              f'\t\t\t<dof>R{ax}</dof>', f'\t\t\t<value lc="1">{F[k]:.12g}</value>',
              '\t\t\t<load_type>LOAD</load_type>', '\t\t\t<relative>0</relative>', '\t\t</rigid_load>']
    for k, ax in enumerate('uvw'):
        L += [f'\t\t<rigid_load name="M{fdi}{ax}" type="rigid_moment">', f'\t\t\t<rb>rb{fdi}</rb>',
              f'\t\t\t<dof>R{ax}</dof>', f'\t\t\t<value lc="1">{M[k]:.12g}</value>',
              '\t\t\t<relative>0</relative>', '\t\t</rigid_load>']
L += ['\t</Rigid>', '\t<LoadData>', '\t\t<load_controller id="1" type="loadcurve">',
      '\t\t\t<interpolate>LINEAR</interpolate>', '\t\t\t<points>', '\t\t\t\t<pt>0,0</pt>', '\t\t\t\t<pt>1,1</pt>',
      '\t\t\t</points>', '\t\t</load_controller>', '\t</LoadData>',
      '\t<Output>', '\t\t<logfile>',
      f'\t\t\t<node_data data="ux;uy;uz" delim="," file="{os.path.basename(a.out)[:-4]}_nodes.txt"/>',
      f'\t\t\t<rigid_body_data data="x;y;z;thx;thy;thz" delim="," file="{os.path.basename(a.out)[:-4]}_rb.txt"/>',
      '\t\t</logfile>', '\t</Output>', '</febio_spec>']
open(a.out, 'w').write('\n'.join(L) + '\n')
pickle.dump(dict(idmap=idmap, case=a.case, force=a.force, steps=a.steps, pdl_mat=a.pdl_mat, auto=a.auto),
            open(a.out + '.map.pkl', 'wb'))
print(f"{a.out}: {nid-1} nodes, {eid-1} elements, case {a.case}, F = {a.force:g} N, {a.steps} steps, PDL {a.pdl_mat}")
