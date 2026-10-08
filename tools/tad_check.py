"""Where are the miniscrews?  Check the buccal and palatal TAD placement against the
paper's specification and against the actual root geometry.

Mao et al. 2024, Methods:
  buccal TAD  -- 'located in the buccal interradicular space between the first and
                 second molars 4 mm above the alveolar crest'
  palatal TAD -- 'located in the palatal interradicular space between the first and
                 second molars 4 mm above the alveolar crest'
  elastics    -- 'at the buccal or lingual mesial cervical region of the canine to
                 simulate the precision cut', 150 g

The bone is not in the shipped model (build/model.pkl is not included), so the
screw is measured against the PDL, whose outer surface is the socket wall.
"""
import sys, os
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, 'ofj4d')); sys.path.insert(0, os.path.join(HERE, 'ofj'))
os.chdir(HERE)
from model4d import Static, State
from loads4d import Elastic
import spec


def unit(v):
    return v/np.linalg.norm(v)


def compute():
    S = Static(); st = State(S)
    for u in S.keep:
        S.tooth[u]['centroid_cur'] = st.centroid(u)
    t16, t17 = S.tooth[3], S.tooth[2]

    # local frame at the 16/17 septum
    eo = unit(t16['e_o'])                                   # occlusal
    es = unit(0.5*(t16['e_b'] + t17['e_b'])); es = unit(es - (es@eo)*eo)   # buccal
    ed = unit(0.5*(t16['e_d'] + t17['e_d']))
    ed = unit(ed - (ed@eo)*eo - (ed@es)*es)                 # distal
    mid = 0.5*(t16['centroid_cur'] + t17['centroid_cur'])
    if (t17['centroid_cur'] - t16['centroid_cur']) @ ed < 0:
        ed = -ed

    # anatomical frame for the elastic's line of action
    a = -eo                                                  # apical
    inc = S.tooth[8]['crown_c']
    p = mid - inc; p = unit(p - (p@eo)*eo)                   # posterior (sagittal)
    l = unit(np.cross(a, p))
    if l @ es < 0:
        l = -l                                               # lateral = buccal

    pdl16, pdl17 = t16['pdl_rest'], t17['pdl_rest']
    out = dict(S=S, st=st, frame=dict(eo=eo, es=es, ed=ed, a=a, p=p, l=l, mid=mid),
               sides={})
    for side, sign in (('buccal', 1.0), ('palatal', -1.0)):
        L = Elastic(S, st, group=f'{side}_tad')
        tad, crest = L.tad, L.crest
        hook = st.aligner[L.hook_nodes].mean(0)
        h = tad @ eo

        def side_slice(P, width=1.0):
            k = (np.abs(P@eo - h) < width) & ((((P - mid)@es)*sign) > 0)
            return P[k]

        w = 1.0
        s16, s17 = side_slice(pdl16, w), side_slice(pdl17, w)
        while (len(s16) == 0 or len(s17) == 0) and w < 4.0:
            w += 0.5
            s16, s17 = side_slice(pdl16, w), side_slice(pdl17, w)

        td = tad @ ed
        gap_lo = (s16 @ ed).max() if len(s16) else np.nan      # distal wall of 16
        gap_hi = (s17 @ ed).min() if len(s17) else np.nan      # mesial wall of 17
        both = np.vstack([s16, s17]) if len(s16) and len(s17) else np.empty((0, 3))
        coord = ((both - mid) @ es)*sign
        outer = coord.max() if len(coord) else np.nan
        inner = coord.min() if len(coord) else np.nan
        t_s = ((tad - mid) @ es)*sign

        v = tad - hook
        u = unit(v)

        def frac(t):
            uu = unit(t - hook)
            return uu @ p, uu @ l, uu @ a

        tad_surface = tad + es*sign*(outer - t_s)
        tad_head = tad_surface + es*sign*2.0

        out['sides'][side] = dict(
            tad=tad, crest=crest, hook=hook, n_hook=len(L.hook_nodes),
            apical_of_crest=float((crest - tad) @ eo),
            crest_below_gingival=float(unit(t16['e_o']) @ np.zeros(3) +
                                       t16['gingival_level'] - crest @ unit(t16['e_o'])),
            slice_width=w,
            gap_lo=gap_lo, gap_hi=gap_hi, td=td,
            clear16=td - gap_lo, clear17=gap_hi - td, gap=gap_hi - gap_lo,
            d16=float(np.linalg.norm(pdl16 - tad, axis=1).min()),
            d17=float(np.linalg.norm(pdl17 - tad, axis=1).min()),
            outer=outer, inner=inner, t_s=t_s, beyond_outer=t_s - outer,
            crest_lateral=float(((crest - mid) @ es)*sign),
            length=float(np.linalg.norm(v)),
            comp=(u@p, u@l, u@a),
            comp_surface=frac(tad_surface), comp_head=frac(tad_head),
            tad_surface=tad_surface, tad_head=tad_head,
            hook_rel=((hook - S.tooth[6]['crown_c']) @ S.tooth[6]['e_d'],
                      (hook - S.tooth[6]['crown_c']) @ S.tooth[6]['e_b'],
                      (hook - S.tooth[6]['crown_c']) @ S.tooth[6]['e_o']))
    return out


def main():
    R = compute()
    B, P = R['sides']['buccal'], R['sides']['palatal']
    es = R['frame']['es']
    F = spec.ELASTIC_FORCE_N
    bar = '='*78

    print(f"\n{bar}\nSPECIFICATION CHECK\n{bar}")
    print(f"  {'':<44}{'buccal':>14}{'palatal':>14}")
    rows = [
        ('apical of the alveolar crest (spec 4.0 mm)', 'apical_of_crest', '{:.2f} mm'),
        ('crest below the gingival level of 16', 'crest_below_gingival', '{:.2f} mm'),
        ('between 16 and 17: clearance to 16 root', 'clear16', '{:.2f} mm'),
        ('between 16 and 17: clearance to 17 root', 'clear17', '{:.2f} mm'),
        ('interradicular gap at screw height', 'gap', '{:.2f} mm'),
        ('nearest PDL of 16 (3D)', 'd16', '{:.2f} mm'),
        ('nearest PDL of 17 (3D)', 'd17', '{:.2f} mm'),
    ]
    for label, k, fmt in rows:
        print(f"  {label:<44}{fmt.format(B[k]):>14}{fmt.format(P[k]):>14}")

    print(f"\n{bar}\nBUCCO-PALATAL POSITION  (mm from the 16/17 midline, positive = toward that side)\n{bar}")
    print(f"  {'':<44}{'buccal':>14}{'palatal':>14}")
    for label, k in (('alveolar crest point', 'crest_lateral'),
                     ('screw', 't_s'),
                     ('outermost root/PDL surface on that side', 'outer'),
                     ('screw relative to that surface', 'beyond_outer')):
        print(f"  {label:<44}{B[k]:>11.2f} mm{P[k]:>11.2f} mm")
    sep = abs((B['tad'] - P['tad']) @ es)
    csep = abs((B['crest'] - P['crest']) @ es)
    print(f"\n  buccal-to-palatal separation:  crest points {csep:.2f} mm,  screws {sep:.2f} mm"
          f"  (3D {np.linalg.norm(B['tad']-P['tad']):.2f} mm)")

    print(f"\n{bar}\nPRECISION CUT (hook) vs canine crown centre\n{bar}")
    for name, s in (('buccal', B), ('palatal', P)):
        d, b, o = s['hook_rel']
        print(f"  {name:<9} {d:+.2f} mm distal  {b:+.2f} mm buccal  {o:+.2f} mm occlusal"
              f"   ({s['n_hook']} aligner nodes)")

    print(f"\n{bar}\nELASTIC LINE OF ACTION  (unit components; force = {F:.3f} N)\n{bar}")
    print(f"  {'':<34}{'sagittal':>10}{'lateral':>10}{'apical':>10}{'length':>10}")
    for name, s in (('buccal TAD, as modelled', B), ('palatal TAD, as modelled', P)):
        c = s['comp']
        print(f"  {name:<34}{c[0]:>10.3f}{c[1]:>10.3f}{c[2]:>10.3f}{s['length']:>8.1f}mm")
    print("\n  what-if: same height and mesiodistal position, moved out along the"
          "\n  bucco-palatal axis (geometric estimates -- the bone is not in the model)")
    for name, s in (('buccal', B), ('palatal', P)):
        for tag, key in (('at the root/PDL surface', 'comp_surface'),
                         ('2 mm beyond it (head)', 'comp_head')):
            c = s[key]
            print(f"  {name+' '+tag:<34}{c[0]:>10.3f}{c[1]:>10.3f}{c[2]:>10.3f}")
    print("\n  lateral > 0 pulls buccally, < 0 palatally; apical > 0 is intrusive on the hook")


if __name__ == '__main__':
    main()
