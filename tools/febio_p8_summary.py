"""Summarise the FEBio Open-Full-Jaw patient 8 tipping step: per-tooth load and crown
surface displacement, and peak displacement of teeth, PDL and bone.

usage: python analyse_p8.py MODEL.feb NODE_DISPLACEMENTS.txt
"""
import re, sys
import numpy as np

feb, disp = sys.argv[1], sys.argv[2]
xyz = {}; parts = {}; part_mat = {}; surfaces = {}; pressures = {}; lcurve = {}
mode = cur = load = None
with open(feb) as f:
    for line in f:
        s = line.strip()
        if s.startswith('<Nodes'):
            mode = 'nodes'; continue
        if s.startswith('<Elements'):
            mode = 'elems'; cur = re.search(r'name="([^"]+)"', s).group(1)
            part_mat[cur] = int(re.search(r'mat="(\d+)"', s).group(1)); parts[cur] = []; continue
        if s.startswith('<Surface '):
            mode = 'surf'; cur = re.search(r'name="([^"]+)"', s).group(1); surfaces[cur] = []; continue
        if s.startswith(('</Nodes', '</Elements', '</Surface')):
            mode = None; continue
        if mode == 'nodes' and s.startswith('<node'):
            m = re.match(r'<node id="(\d+)">([^<]+)</node>', s)
            xyz[int(m.group(1))] = [float(v) for v in m.group(2).split(',')]
        elif mode == 'elems' and s.startswith('<elem'):
            parts[cur].append([int(v) for v in re.search(r'>([^<]+)<', s).group(1).split(',')])
        elif mode == 'surf' and s.startswith(('<tri3', '<quad4')):
            surfaces[cur].append([int(v) for v in re.search(r'>([^<]+)<', s).group(1).split(',')])
        m = re.match(r'<surface_load surface="([^"]+)"', s)
        if m:
            load = m.group(1)
        m = re.match(r'<pressure lc="(\d+)">([^<]+)</pressure>', s)
        if m and load:
            pressures[load] = float(m.group(2)); lcurve[load] = int(m.group(1))

# displacements: last data block of the log file
U = {}
block = {}
with open(disp) as f:
    for line in f:
        if line.startswith('*'):
            if line.startswith('*Step') and block:
                U = block; block = {}
            continue
        v = line.strip().split(',')
        if len(v) >= 4:
            block[int(v[0])] = [float(x) for x in v[1:4]]
if block:
    U = block

ids = np.array(sorted(xyz)); X = np.array([xyz[i] for i in ids])
Uarr = np.array([U.get(i, [np.nan]*3) for i in ids]); idx = {n: k for k, n in enumerate(ids)}
mag = np.linalg.norm(Uarr, axis=1)
print(f"nodes {len(ids)}, displacement records {len(U)}; units as in the model (mm, MPa, N)")
names = {1: 'teeth', 2: 'PDL', 3: 'bone'}
for p, el in parts.items():
    nodes = np.unique(np.array(el).ravel()); k = [idx[n] for n in nodes]
    print(f"  {names.get(part_mat[p], p):<6} ({p}, {len(el)} elements): max |u| {np.nanmax(mag[k])*1e3:8.2f} um, "
          f"mean |u| {np.nanmean(mag[k])*1e3:7.2f} um")

LOAD_FRACTION = 0.1          # linear load curve 0 -> 1 over t = 0..1, one step of 0.1
print(f"\nper tooth (load at t = 0.1, i.e. {LOAD_FRACTION:.0%} of the full pressure):")
print(f"  {'surface':<22}{'pressure':>9}{'area mm2':>10}{'force N':>9}{'mean |u| um':>13}{'max |u| um':>12}   mean u (x, y, z) um")
def key(name):
    m = re.search(r'(\d+)$', name); return int(m.group(1)) if m else 0
for sname in sorted(surfaces, key=key):
    faces = np.array(surfaces[sname]); tri = faces[:, :3]
    P = np.array([[xyz[n] for n in f] for f in tri])
    area = 0.5*np.linalg.norm(np.cross(P[:, 1]-P[:, 0], P[:, 2]-P[:, 0]), axis=1).sum()
    nodes = np.unique(faces.ravel()); k = [idx[n] for n in nodes]
    p = pressures.get(sname, np.nan)
    um = np.nanmean(Uarr[k], axis=0)*1e3
    print(f"  {sname:<22}{p:9.3f}{area:10.2f}{LOAD_FRACTION*p*area:9.3f}{np.nanmean(mag[k])*1e3:13.2f}"
          f"{np.nanmax(mag[k])*1e3:12.2f}   ({um[0]:+.2f}, {um[1]:+.2f}, {um[2]:+.2f})")
