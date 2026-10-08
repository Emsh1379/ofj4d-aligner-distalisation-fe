"""Reduced FEBio model of Open-Full-Jaw patient 8 for the cross-validation.

From the dataset's full-maxilla P8_maxilla_tipping.feb (FEBio 2.5 format) keep the
right half arch UNN 2-8 (FDI 17-11, the teeth of this project's model):
  * every tooth becomes its own rigid body (no tooth deformation, 6 dof per tooth);
  * the bone is removed and the PDL nodes it shared are fixed (nodes that also belong
    to a tooth are not fixed, or the rigid tooth would be pinned);
  * the same pressure loads, load curves and solver controls;
  * output: displacement only (plot file) plus a node displacement log.

usage: python tools/febio_reduce_p8.py FULL.feb OUT.feb TIME_STEPS NODE_LOG_FILE
writes OUT.feb and OUT.feb.map (new node id, original node id)
"""
import re, sys
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

full, out, nsteps, nodelog = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
KEEP_UNN = [2, 3, 4, 5, 6, 7, 8]
L = open(full).read().split('\n')


def find(pat, start=0):
    for i in range(start, len(L)):
        if pat in L[i]:
            return i
    raise KeyError(pat)


def section(i0, closing):
    return i0, find(closing, i0)


# ---- nodes, element parts, bone node set, surfaces, loads --------------------------
n0, n1 = section(find('<Nodes'), '</Nodes>')
coords = {}
for l in L[n0+1:n1]:
    m = re.match(r'\s*<node id="(\d+)">([^<]+)</node>', l)
    if m:
        coords[int(m.group(1))] = m.group(2).strip()


def elements(name):
    i0, i1 = section(find(f'name="{name}"'), '</Elements>')
    return np.array([[int(v) for v in re.search(r'>([^<]+)<', l).group(1).split(',')]
                     for l in L[i0+1:i1] if '<elem' in l], dtype=np.int64)


teeth, pdl, bone = elements('Part1'), elements('Part2'), elements('Part3')
surf = {}
for i, l in enumerate(L):
    m = re.match(r'\s*<Surface name="([^"]+)">', l)
    if m:
        j, faces = i + 1, []
        while '</Surface>' not in L[j]:
            faces.append([int(v) for v in re.search(r'>([^<]+)<', L[j]).group(1).split(',')]); j += 1
        surf[m.group(1)] = faces


def components(E):
    un, inv = np.unique(E.ravel(), return_inverse=True)
    ne, nn = len(E), len(un)
    r = np.repeat(np.arange(ne), E.shape[1]); c = ne + inv
    g = coo_matrix((np.ones(len(r)), (r, c)), shape=(ne + nn, ne + nn))
    _, lab = connected_components(g, directed=False)
    return lab[:ne], dict(zip(un.tolist(), lab[ne:].tolist()))


t_elab, t_nlab = components(teeth)
p_elab, p_nlab = components(pdl)

# tooth component of each kept UNN, from its load triangle
tooth_comp = {}
for u in KEEP_UNN:
    tri = surf[f'Pressure_load_tooth{u}'][0]
    labs = {t_nlab.get(n) for n in tri}
    assert len(labs) == 1 and None not in labs, f'load triangle of UNN {u} not on one tooth: {labs}'
    tooth_comp[u] = labs.pop()

# PDL component of each kept tooth: the one sharing most nodes with it
tooth_nodes = {u: set(np.unique(teeth[t_elab == c]).tolist()) for u, c in tooth_comp.items()}
pdl_comp = {}
for u, tn in tooth_nodes.items():
    counts = {}
    for n in tn:
        if n in p_nlab:
            counts[p_nlab[n]] = counts.get(p_nlab[n], 0) + 1
    pdl_comp[u] = max(counts, key=counts.get)
assert len(set(pdl_comp.values())) == len(pdl_comp), 'two teeth share one PDL component'

keep_teeth = {u: teeth[t_elab == c] for u, c in tooth_comp.items()}
keep_pdl = pdl[np.isin(p_elab, list(pdl_comp.values()))]
all_tooth_nodes = set().union(*tooth_nodes.values())
bone_nodes = set(np.unique(bone).tolist())
pdl_nodes = set(np.unique(keep_pdl).tolist())
fixed = sorted((pdl_nodes & bone_nodes) - all_tooth_nodes)

# ---- renumber ----------------------------------------------------------------------
used = sorted(pdl_nodes | all_tooth_nodes)
new = {o: i + 1 for i, o in enumerate(used)}

o = []
o += L[:find('<Control>')]
c0, c1 = section(find('<Control>'), '</Control>')
o += [re.sub(r'<time_steps>\d+</time_steps>', f'<time_steps>{nsteps}</time_steps>', l) for l in L[c0:c1+1]]
g0, g1 = section(find('<Globals>'), '</Globals>')
o += L[g0:g1+1]
o += ['\t\t<Material>',
      '\t\t\t<material id="1" name="PDL_mat" type="neo-Hookean">',
      '\t\t\t\t<density>1</density>', '\t\t\t\t<E>0.0689</E>', '\t\t\t\t<v>0.45</v>',
      '\t\t\t</material>']
for k, u in enumerate(KEEP_UNN):
    o += [f'\t\t\t<material id="{k+2}" name="rigid_tooth{u}" type="rigid body">',
          '\t\t\t\t<density>1</density>', '\t\t\t</material>']
o += ['\t\t</Material>', '\t\t<Geometry>', '\t\t\t<Nodes name="reduced">']
o += [f'\t\t\t\t<node id="{new[n]}">{coords[n]}</node>' for n in used]
o += ['\t\t\t</Nodes>']
eid = 1


def write_elems(E, mat, name):
    global eid
    o.append(f'\t\t\t<Elements mat="{mat}" name="{name}" type="tet4">')
    for e in E:
        o.append(f'\t\t\t\t<elem id="{eid}">' + ','.join(str(new[n]) for n in e) + '</elem>'); eid += 1
    o.append('\t\t\t</Elements>')


write_elems(keep_pdl, 1, 'PDL')
for k, u in enumerate(KEEP_UNN):
    write_elems(keep_teeth[u], k + 2, f'tooth{u}')
o += ['\t\t\t<NodeSet name="PDL_outer">'] + [f'\t\t\t\t<node id="{new[n]}" />' for n in fixed] + ['\t\t\t</NodeSet>']
for u in KEEP_UNN:
    o.append(f'\t\t\t<Surface name="Pressure_load_tooth{u}">')
    for j, f in enumerate(surf[f'Pressure_load_tooth{u}']):
        o.append(f'\t\t\t\t<tri3 id="{j+1}">' + ','.join(str(new[n]) for n in f) + '</tri3>')
    o.append('\t\t\t</Surface>')
o += ['\t\t</Geometry>', '\t\t<Boundary>', '\t\t\t<fix bc="x,y,z" node_set="PDL_outer" />', '\t\t</Boundary>', '\t\t<Loads>']
for u in KEEP_UNN:
    i0, i1 = section(find(f'<surface_load surface="Pressure_load_tooth{u}"'), '</surface_load>')
    o += L[i0:i1+1]
o += ['\t\t</Loads>']
d0, d1 = section(find('<LoadData>'), '</LoadData>')
o += L[d0:d1+1]
o += ['\t\t<Output>', '\t\t\t<plotfile type="febio">', '\t\t\t\t<var type="displacement" />', '\t\t\t</plotfile>',
      '\t\t\t<logfile>', f'\t\t\t\t<node_data data="ux;uy;uz" delim="," file="{nodelog}"></node_data>',
      '\t\t\t</logfile>', '\t\t</Output>', '</febio_spec>']
open(out, 'w').write('\n'.join(o) + '\n')
with open(out + '.map', 'w') as f:
    for n in used:
        f.write(f'{new[n]},{n}\n')
ntooth = sum(len(e) for e in keep_teeth.values())
print(f"reduced model: {len(used)} nodes ({len(used)}/{len(coords)}), PDL {len(keep_pdl)} tets "
      f"(of {len(pdl)}), rigid tooth tets {ntooth} (of {len(teeth)}), bone removed ({len(bone)} tets), "
      f"fixed outer PDL nodes {len(fixed)}, time steps {nsteps}")
print("estimated equations: about", 3*(len(pdl_nodes) - len(fixed) - len(pdl_nodes & all_tooth_nodes)) + 6*len(KEEP_UNN))
