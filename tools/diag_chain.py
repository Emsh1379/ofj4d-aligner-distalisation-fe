"""Reproduce the real failure: run the staging chain (carrying cache and E_pdl
across steps, as run4d does) and dump the full trace of the bad step."""
import sys, os
import numpy as np
HERE=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE,'ofj4d')); sys.path.insert(0, os.path.join(HERE,'ofj')); os.chdir(HERE)
from model4d import Static, State
from arch4d import ArchForm
from morph4d import Morpher
from loads4d import Elastic
from solve4d import solve_config, MAX_GAP
from contact import build_pairs
import spec
GROUP=sys.argv[1] if len(sys.argv)>1 else 'buccal_tad'
LAST=int(sys.argv[2]) if len(sys.argv)>2 else 16
S=Static(); state=State(S); arch=ArchForm(S); morph=Morpher(S,arch)
for u in S.keep: S.tooth[u]['centroid_cur']=state.centroid(u)
load=Elastic(S,state,group=GROUP)
slave=np.unique(S.aligner_inner_faces)
cache={}; E_prev=None; traces={}
for step in range(1, LAST+1):
    shape=morph.shape(step)
    pre=S.prescribed_nodes(state,arch,step,shape)
    Cp=build_pairs(pre,slave,S.master_tris,k=16,max_gap=MAX_GAP)
    gp,_=S.gap_for_pairing(pre,slave,Cp['tri'],Cp['bary'])
    cref=dict(slave=slave,tri=Cp['tri'],bary=Cp['bary'],gap_pre=gp)
    for k in range(spec.PDL_ITERATIONS_PER_STEP):
        state.aligner=morph.seat(shape,state,step)
        r=solve_config(S,state,load,E_pdl0=E_prev,x0=cache.get('x0'),
                       cache=cache,contact_ref=cref)
        E_prev=r['E_pdl']
        for u in S.keep:
            q=r['q'][u]; state.apply_rigid(u,q[:3],q[3:])
        traces[(step,k)]=r['hist']
    h0,h1=traces[(step,0)],traces[(step,1)]
    print(f"  step {step:2d}: solves {len(h0)}+{len(h1)} it, "
          f"final du {h0[-1]['du']:.1e}/{h1[-1]['du']:.1e}, "
          f"dA {h0[-1]['dA']}/{h1[-1]['dA']}, Fc {h0[-1]['Fc']:.2f}/{h1[-1]['Fc']:.2f} N",
          flush=True)
h=traces[(LAST,0)]
print(f"\n  full trace of step {LAST}, solve 0:")
print(f"  {'it':>3}{'active':>8}{'dA':>7}{'|Fc| N':>11}{'du':>10}{'dF':>10}{'dE':>10}")
for e in h:
    print(f"  {e['it']:>3}{e['active']:>8}{e['dA']:>7}{e['Fc']:>11.3f}"
          f"{e['du']:>10.2e}{e['dF']:>10.2e}{e['dE']:>10.2e}")
