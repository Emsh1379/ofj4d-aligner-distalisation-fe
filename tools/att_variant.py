"""Attachment-size variant of the half-arch model.

Rebuilds every attachment on the shipped geometry (build4d/geom.pkl) at a scaled
size and re-wraps the aligner around the new attachments.  build/stage2.pkl is not
needed: geom.pkl carries the crown faces, tooth axes and gingival levels the
attachment builder uses.

The attachment thickness is set by geom4d.PROTRUSION (the outer face sits that far
beyond the most prominent crown point under it), not by the depth in
spec.ATTACHMENTS, so the thickness scale is applied there.

usage: python tools/att_variant.py FOOTPRINT_SCALE THICKNESS_SCALE TAG
  -> build4d/variants/geom_TAG.pkl, build4d/variants/aligner_TAG.pkl
run with OFJ_GEOM=... OFJ_ALIGNER=... to use it.
"""
import sys, os, pickle
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, 'ofj4d')); sys.path.insert(0, os.path.join(HERE, 'ofj'))
os.chdir(HERE)
import numpy as np
import spec
import geom4d
import aligner4d


def main(fp, th, tag):
    G = pickle.load(open('build4d/geom.pkl', 'rb'))
    nodes = G['nodes']
    geom4d.PROTRUSION = geom4d.PROTRUSION*th
    for u in G['keep']:
        t = G['teeth'][u]
        if t['att'] is None:
            continue
        kind, (md, og, dp) = spec.ATTACHMENTS[u]
        size = (md*fp, og*fp, dp*th)
        old = t['att']
        nodes, a_tets, a_faces, a_vol, a_span = geom4d.build_attachment(nodes, t, size, u)
        t['att'] = dict(tets=a_tets, faces=a_faces, volume=a_vol, span=a_span,
                        label=kind, nominal=size)
        print(f"FDI {spec.FDI[u]}: {kind} nominal {size[0]:.2f} x {size[1]:.2f} mm, "
              f"protrusion {geom4d.PROTRUSION:.2f} mm -> span {a_span[0]:.2f} x {a_span[1]:.2f} mm, "
              f"vol {a_vol:.2f} mm3  (shipped: span {old['span'][0]:.2f} x {old['span'][1]:.2f}, "
              f"vol {old['volume']:.2f})")
    G['nodes'] = nodes
    os.makedirs('build4d/variants', exist_ok=True)
    gp = f'build4d/variants/geom_{tag}.pkl'
    pickle.dump(G, open(gp, 'wb'))
    aligner4d.build(geom=gp, out=f'build4d/variants/aligner_{tag}.pkl')


if __name__ == '__main__':
    main(float(sys.argv[1]), float(sys.argv[2]), sys.argv[3])
