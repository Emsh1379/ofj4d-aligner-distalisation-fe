#!/usr/bin/env bash
# rigid vs deformable-bone verification on one solve of step 1 (buccal TAD, final model)
cd "$(dirname "$0")/../.."
export OFJ_PRESCRIBE=tip OFJ_TIP_BLEND=0.40
for cfg in "rigid:1:rigid" "deformable:1000:stiff1000" "deformable:1:real"; do
  IFS=: read bone scale tag <<< "$cfg"
  s=$(date +%s)
  OFJ_BONE=$bone OFJ_BONE_E_SCALE=$scale python -u tools/bone_check.py buccal_tad build4d/bone_check/$tag.pkl 2>&1 \
    | grep -E "bone model|global nodes|reduced dofs|Traceback|Error"
  echo "  [$tag] $(( $(date +%s)-s )) s"
done
python - <<'PY'
import pickle, numpy as np, sys
sys.path.insert(0, 'ofj4d'); import spec
R = pickle.load(open('build4d/bone_check/rigid.pkl', 'rb'))
for tag in ('stiff1000', 'real'):
    D = pickle.load(open(f'build4d/bone_check/{tag}.pkl', 'rb'))
    print(f"\n{tag}: deformable={D['deformable']}, reduced dofs {D['n_red']}, iterations {len(D['hist'])} (rigid {len(R['hist'])}), "
          f"|Fc| {D['hist'][-1]['Fc']:.4f} N (rigid {R['hist'][-1]['Fc']:.4f}), peak bone displacement {D.get('bone_umax', float('nan'))*1e3:.3f} um, "
          f"ties inside bone {D['tie_inside']:.3f}, TAD reaction nodes {D['n_bone_react']}, wall {D['wall']:.0f} s (rigid {R['wall']:.0f} s)")
    for u in R['q']:
        qr, qd = R['q'][u], D['q'][u]
        tr, td = np.linalg.norm(qr[:3]), np.linalg.norm(qd[:3])
        print(f"   FDI {spec.FDI[u]}: translation {tr*1e3:8.3f} -> {td*1e3:8.3f} um ({100*(td-tr)/tr:+.3f} %), "
              f"|dq_t| {np.linalg.norm(qd[:3]-qr[:3])*1e6:8.2f} nm, rotation diff {np.degrees(np.linalg.norm(qd[3:]-qr[3:])):.5f} deg")
PY
echo VERIFY_DONE
