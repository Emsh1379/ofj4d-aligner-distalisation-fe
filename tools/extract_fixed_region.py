"""Coordinates of the dataset authors' fixed bone region, from their FEBio file.

Open-Full-Jaw ships a tipping-scenario FEBio model per patient whose node set
'Bone_fixed_region' is the support the authors used.  That .feb was meshed
separately from the volumetric .msh (its nodes sit a median 0.14 mm away), so the
set is exported as coordinates and ofj4d/bone4d.py reproduces the region on the
.msh -- it is a superior slab, z >= 78.38 mm for patient 8.

    python tools/extract_fixed_region.py <P8_maxilla_tipping.feb> <out.npy>
"""
import re, sys
import numpy as np


def extract(feb, out):
    txt = open(feb, encoding='utf-8', errors='ignore').read()
    m = re.search(r'<NodeSet[^>]*name="Bone_fixed_region"[^>]*>(.*?)</NodeSet>', txt, re.S)
    if m is None:
        raise SystemExit('Bone_fixed_region not found')
    ids = [int(x) for x in re.findall(r'id="(\d+)"', m.group(1))]
    nd = re.search(r'<Nodes[^>]*>(.*?)</Nodes>', txt, re.S).group(1)
    pairs = re.findall(r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*</node>', nd)
    lut = {int(p[0]): (float(p[1]), float(p[2]), float(p[3])) for p in pairs}
    P = np.array([lut[i] for i in ids if i in lut])
    np.save(out, P)
    print(f'{len(P)} fixed nodes, z {P[:, 2].min():.2f} .. {P[:, 2].max():.2f} -> {out}')


if __name__ == '__main__':
    extract(sys.argv[1], sys.argv[2])
