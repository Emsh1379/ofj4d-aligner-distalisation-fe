"""Full vs reduced FEBio patient 8 model: displacement of each kept tooth's load
triangle at a given time (default t = 0.1).

usage: python tools/febio_compare_reduced.py FULL.feb FULL_NODES.txt REDUCED.feb.map REDUCED_NODES.txt [TIME]
"""
import re, sys
import numpy as np

full_feb, full_nodes, redmap, red_nodes = sys.argv[1:5]
T = float(sys.argv[5]) if len(sys.argv) > 5 else 0.1
KEEP = {2: 17, 3: 16, 4: 15, 5: 14, 6: 13, 7: 12, 8: 11}      # UNN: FDI


def blocks(path):
    out, cur, t = {}, None, None
    for line in open(path):
        if line.startswith('*Time'):
            t = float(line.split('=')[1]); cur = out.setdefault(t, {})
        elif line[:1].isdigit() and cur is not None:
            v = line.strip().split(',')
            cur[int(v[0])] = np.array([float(x) for x in v[1:4]])
    return out


def at(b, t):
    k = min(b, key=lambda x: abs(x - t))
    return b[k], k


surf = {}
lines = open(full_feb).read().split('\n')
for i, l in enumerate(lines):
    m = re.match(r'\s*<Surface name="Pressure_load_tooth(\d+)">', l)
    if m and int(m.group(1)) in KEEP:
        surf[int(m.group(1))] = [int(v) for v in re.search(r'>([^<]+)<', lines[i+1]).group(1).split(',')]
o2n = {}
for line in open(redmap):
    a, b = line.strip().split(','); o2n[int(b)] = int(a)

uf, tf = at(blocks(full_nodes), T)
ur, tr = at(blocks(red_nodes), T)
print(f"load-triangle displacement at t = {tf:g} (full) / {tr:g} (reduced), um")
print(f"  {'tooth':<14}{'full |u|':>10}{'reduced |u|':>13}{'reduced/full':>14}   direction angle")
for u, fdi in KEEP.items():
    a = np.mean([uf[n] for n in surf[u]], axis=0)
    b = np.mean([ur[o2n[n]] for n in surf[u]], axis=0)
    ang = np.degrees(np.arccos(np.clip(a @ b/(np.linalg.norm(a)*np.linalg.norm(b)), -1, 1)))
    print(f"  UNN {u} / FDI {fdi:<3}{np.linalg.norm(a)*1e3:10.2f}{np.linalg.norm(b)*1e3:13.2f}"
          f"{np.linalg.norm(b)/np.linalg.norm(a):14.2f}   {ang:6.1f} deg")
