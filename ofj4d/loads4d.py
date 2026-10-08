"""Anchorage groups: the elastic from a TAD to a precision cut on the aligner.

Buccal TAD group (Mao et al. 2024):
    'the TAD was located in the buccal interradicular space between the first and
     second molars 4 mm above the alveolar crest ... The application of elastics on
     the CAs were all set at the buccal or lingual mesial cervical region of the
     canine to simulate the precision cut ... the elastic forces for all the test
     groups were set as 150 g'

The TAD is a fixed point in the bone.  The precision cut is a material patch on
the aligner, so it travels with the appliance; the line of action is recomputed
from the CURRENT hook position at every iteration, as a real elastic would.
"""
import sys, os, pickle
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scipy.spatial import cKDTree
import spec
from inputs import load as load_inputs

HOOK_RADIUS = 1.75        # mm, precision-cut patch radius
TAD_PATCH = 0.6           # mm, bone nodes around the TAD point that take the reaction


class Elastic:
    """150 gf elastic, TAD -> precision cut at the canine."""

    def __init__(self, S, state, group='buccal_tad', model='build/model.pkl'):
        self.group = group
        g = spec.GROUPS[group]
        self.enabled = g['elastic']
        if not self.enabled:
            self.hook_nodes = np.array([], np.int64)
            self.tad = None
            return

        side = g['side']
        ua, ub = g['tad_between']                 # (16, 17)
        ta, tb = S.tooth[ua], S.tooth[ub]
        e_o = ta['e_o']/np.linalg.norm(ta['e_o'])
        e_side = 0.5*(ta['e_b'] + tb['e_b']); e_side /= np.linalg.norm(e_side)
        if side == 'palatal':
            e_side = -e_side

        # ---- the miniscrew point -------------------------------------------------
        # Preferred: build4d/tad_bone.pkl from tools/tad_from_bone.py -- a point on
        # the cortical surface of the segmented bone, 4 mm apical to that side's
        # alveolar crest, centred in the 16/17 interradicular gap.  The legacy rule
        # (most occlusal bone node on each side of the root midline, offset straight
        # apically) lands on the septum ridge and put both screws in the middle of
        # the septum, 3.5 mm apart (FINDINGS 19-20).  OFJ_TAD_RULE=legacy restores it.
        tb = os.environ.get('OFJ_TAD_FILE', os.path.join('build4d', 'tad_bone.pkl'))   # position variants: tools/tad_from_bone.py
        if os.path.exists(tb) and os.environ.get('OFJ_TAD_RULE', 'bone') == 'bone':
            T = pickle.load(open(tb, 'rb'))['sides'][side]
            self.tad = np.asarray(T['tad'], float)
            self.crest = self.tad + e_o*g['tad_above_crest']     # reference only
            self.tad_rule = f'bone surface ({os.path.basename(tb)})'
        else:
            IN, _ = load_inputs(model)
            crest = IN['crest'][side]
            # 'above the alveolar crest' = apical, i.e. away from the occlusal plane
            self.tad = crest - e_o*g['tad_above_crest']
            self.crest = crest
            self.tad_rule = 'legacy septum-ridge crest'

        # ---- precision cut: mesial cervical region of the canine, on the aligner --
        can = S.tooth[g['hook_tooth']]
        e_b = can['e_b']; e_d = can['e_d']; e_oc = can['e_o']
        if side == 'palatal':
            e_b = -e_b
        cerv = can['gingival_level'] + 0.35*(can['occl_tip'] - can['gingival_level'])
        seed = (can['crown_c'] - e_d*2.2          # mesial
                + e_b*6.0                          # push outside, snapped back below
                + e_oc*(cerv - can['crown_c'] @ e_oc))
        al = state.aligner
        # restrict to the outer surface layer (ids >= n_inner) ...
        cand = np.arange(S.n_inner, getattr(S, 'nA_layers', S.nA))   # layer nodes only, not prism centroids
        # ... lying over the hook tooth itself, on the requested side of it.  The
        # nearest outer node to the seed is not enough: on the palatal side the arch
        # curves, and the closest appliance surface to a point palatal of the
        # canine belongs to the lateral incisor, so the cut landed there.
        hook_u = g['hook_tooth']
        ids = list(S.keep)
        D = np.stack([cKDTree(S.tooth[u]['surf_rest']).query(S.aligner_rest[cand])[0]
                      for u in ids], axis=1)
        on_tooth = np.array(ids)[D.argmin(1)] == hook_u
        sgn = 1.0 if side == 'buccal' else -1.0
        side_ok = ((S.aligner_rest[cand] - can['crown_c']) @ can['e_b'])*sgn > 0
        cand = cand[on_tooth & side_ok]
        if len(cand) == 0:
            raise RuntimeError(f'no aligner surface over the {side} side of FDI {spec.FDI[hook_u]}')
        j = cand[np.argmin(np.linalg.norm(al[cand] - seed, axis=1))]
        centre = al[j]
        self.hook_nodes = cand[np.linalg.norm(al[cand] - centre, axis=1) <= HOOK_RADIUS]
        self.hook_rest = al[self.hook_nodes].copy()
        self.magnitude = spec.ELASTIC_FORCE_N
        # With deformable bone the elastic pulls on the bone as well: the equal and
        # opposite force goes into the bone at the TAD point.  With rigid bone the
        # reaction is absorbed by the fixed support, as before.
        self.bone_nodes = np.array([], np.int64)
        if getattr(S, 'deformable_bone', False):
            Bx = S.nodes_rest[S.bone_off:S.bone_off + S.nB]
            dist = np.linalg.norm(Bx - self.tad, axis=1)
            near = np.where(dist <= TAD_PATCH)[0]
            if len(near) < 3:
                near = np.argsort(dist)[:3]
            near = near + S.bone_off
            self.bone_nodes = near[~np.isin(near, S.bone_fixed)]

    def direction(self, state):
        c = state.aligner[self.hook_nodes].mean(0)
        v = self.tad - c
        return v/np.linalg.norm(v), c

    def vector(self, S, state, ndof):
        f = np.zeros(ndof)
        if not self.enabled:
            return f
        dirn, c = self.direction(state)
        F = self.magnitude*dirn/len(self.hook_nodes)
        for d in range(3):
            np.add.at(f, 3*self.hook_nodes + d, F[d])
        if len(self.bone_nodes):
            Fb = -self.magnitude*dirn/len(self.bone_nodes)
            for d in range(3):
                np.add.at(f, 3*self.bone_nodes + d, Fb[d])
        return f

    def describe(self, S, state):
        if not self.enabled:
            return f"group '{self.group}': no elastic (control)"
        dirn, c = self.direction(state)
        can = S.tooth[spec.GROUPS[self.group]['hook_tooth']]
        rel = c - can['crown_c']
        return (f"group '{self.group}'   [TAD rule: {self.tad_rule}]\n"
                f"  TAD            {np.round(self.tad, 2)}  "
                f"({spec.GROUPS[self.group]['tad_above_crest']:.0f} mm apical to the "
                f"crest at {np.round(self.crest, 2)})\n"
                f"  precision cut  {np.round(c, 2)}  on {len(self.hook_nodes)} aligner nodes\n"
                f"                 vs canine crown centre: "
                f"{rel @ can['e_d']:+.2f} mm distal, {rel @ can['e_b']:+.2f} mm buccal, "
                f"{rel @ can['e_o']:+.2f} mm occlusal\n"
                f"  force          {self.magnitude:.4f} N ({spec.ELASTIC_FORCE_GF:.0f} gf) "
                f"along {np.round(dirn, 3)}\n"
                f"  free length    {np.linalg.norm(self.tad - c):.2f} mm")
