"""Is the PDL model itself responsible for the movement-type discrepancy?

Removes the aligner and the contact entirely and asks the classical question of
orthodontic biomechanics directly of the PDL: where is this tooth's centre of
resistance, and where does it rotate about when a single force is applied at the
crown?  The paper reports "the centre of rotation ... around the apical third".

The teeth are rigid bodies with 6 dofs in the reduced system, so the tooth's
effective 6x6 stiffness is the Schur complement of the PDL stiffness onto those
dofs -- exact, and six linear solves.
"""
import sys, os
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE,'ofj4d')); sys.path.insert(0, os.path.join(HERE,'ofj')); os.chdir(HERE)
from model4d import Static, State
from solve4d import Operators
import spec

TOOTH = int(sys.argv[1]) if len(sys.argv) > 1 else 2          # FDI 17
NU = float(os.environ.get('OFJ_PDL_NU', spec.PDL_NU))

S = Static(); st = State(S)
for u in S.keep:
    S.tooth[u]['centroid_cur'] = st.centroid(u)
nodes = st.nodes(); ndof = 3*S.n_global
P = S.build_P(nodes).tocsr()
opP = Operators(nodes, S.pdl_tets, NU)
K = (P.T @ opP.K(np.full(len(S.pdl_tets), spec.PDL_E1), ndof) @ P).tocsr()

q0 = S.q_col[TOOTH]
qd = np.arange(q0, q0+6)
keep = np.where(np.abs(K.diagonal()) > 0)[0]          # dofs the PDL actually touches
other = np.setdiff1d(keep, qd)

Kqq = K[qd][:, qd].toarray()
Kqf = K[qd][:, other]
Kff = K[other][:, other].tocsc()
print(f"FDI {spec.FDI[TOOTH]}: PDL-only system, {len(other)} free dofs "
      f"+ 6 rigid dofs, nu = {NU}")
lu = splu(Kff)
X = np.column_stack([lu.solve(np.asarray(Kqf.T.todense())[:, i]) for i in range(6)])
Keff = Kqq - Kqf @ X                                   # 6x6, exact Schur complement

T = S.tooth[TOOTH]
c = T['centroid_cur']
e_d, e_b, e_o = T['e_d'], T['e_b'], T['e_o']
crown = T['surf_rest'][:T['surf_n']]
crown_pt = crown[np.argmax(crown @ e_o)]
apex = T['pdl_rest'][np.argmin(T['pdl_rest'] @ e_o)]
L = (crown_pt - apex) @ e_o                            # crown-to-apex along the long axis
print(f"  crown point to apex along the long axis: {L:.2f} mm")

C = np.linalg.inv(Keff)


def motion(F, p):
    """Rigid motion (t, w) from force F applied at point p."""
    M = np.cross(p - c, F)
    q = C @ np.concatenate([F, M])
    return q[:3], q[3:]


def centre_of_rotation(t, w, axis=e_o, origin=None):
    """Point on the long axis where the displacement vanishes, as a distance
    below the crown point measured along the long axis."""
    origin = crown_pt if origin is None else origin
    # u(s) = t + w x (origin + s*axis - c); minimise |u| over s
    a = np.cross(w, axis)
    b = t + np.cross(w, origin - c)
    if a @ a < 1e-30:
        return np.inf
    return -(a @ b)/(a @ a)


# --- centre of resistance: the point where a force gives no rotation ------------
# search along the long axis for the application point with zero rotation
def rot_for(s):
    t, w = motion(e_d, crown_pt + s*e_o*-1.0)          # e_o points occlusally
    return w @ e_b


lo, hi = -5.0, 30.0
for _ in range(60):
    mid = 0.5*(lo+hi)
    if rot_for(lo)*rot_for(mid) <= 0:
        hi = mid
    else:
        lo = mid
s_cres = 0.5*(lo+hi)

t, w = motion(e_d, crown_pt)
s_crot = -centre_of_rotation(t, w)

print(f"\n  distal force of 1 N at the crown point:")
print(f"    crown displacement   {(t + np.cross(w, crown_pt-c)) @ e_d*1e3:8.2f} um")
print(f"    apex  displacement   {(t + np.cross(w, apex-c)) @ e_d*1e3:8.2f} um")
print(f"    rotation             {np.degrees(w @ e_b):8.4f} deg")
print(f"\n  centre of RESISTANCE  {s_cres:6.2f} mm below the crown point "
      f"({100*s_cres/L:5.1f}% of the crown-apex distance)")
print(f"  centre of ROTATION    {s_crot:6.2f} mm below the crown point "
      f"({100*s_crot/L:5.1f}% of the crown-apex distance)")
print(f"\n  the apex is at 100%.  The paper puts the centre of rotation in the")
print(f"  apical third of the ROOT, i.e. well inside this span, not beyond it.")

# translational and rotational compliance, for reference against the literature
print(f"\n  single-tooth stiffness, distal translation at the centroid: "
      f"{1.0/ (C[0:3,0:3] @ e_d @ e_d):6.2f} N/mm")
