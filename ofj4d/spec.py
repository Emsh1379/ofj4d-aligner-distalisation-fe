"""Specification of Mao et al., Bioengineering 2024, 11, 3 -- buccal TAD group.

'Effect of Different Anchorage Reinforcement Methods on Long-Term Maxillary Whole
Arch Distalization with Clear Aligner: A 4D Finite Element Study with Staging
Simulation.'  doi:10.3390/bioengineering11010003

Every value below is quoted from that paper (or, where it says 'following the
established methodology described in previous publications', from the companion
method paper Mao et al., Prog Orthod 2023;24:16, doi:10.1186/s40510-023-00468-1).
Anything NOT specified by either paper is marked ASSUMED.
"""
import os
import numpy as np

CITATION = ("Mao B, Tian Y, Xiao Y, Liu J, Liu D, Zhou Y, Li J. Bioengineering 2024;11(1):3. "
            "doi:10.3390/bioengineering11010003")

# ---------------------------------------------------------------- dentition --
# "To simplify the model, only the right half of the maxillary dentition was
#  involved. A contact surface was formed on the midsagittal plane for symmetric
#  analysis, and the front end of the CA, on the midsagittal plane, was
#  constrained in the x-direction."                       [Prog Orthod 2023]
HALF_ARCH_UNN = [2, 3, 4, 5, 6, 7, 8]          # FDI 17,16,15,14,13,12,11
FDI = {2: 17, 3: 16, 4: 15, 5: 14, 6: 13, 7: 12, 8: 11}
UNN = {v: k for k, v in FDI.items()}

# ---------------------------------------------------------------- materials --
# "the elastic modulus values were set at 1500 MPa and 20,000 MPa, respectively,
#  with Poisson's ratios of 0.30 for both"                          [Bioeng 2024]
E_ALIGNER, NU_ALIGNER = 1500.0, 0.30
E_ALIGNER *= float(os.environ.get('OFJ_E_ALIGNER_SCALE', 1.0))   # diagnostic only
E_ATTACH,  NU_ATTACH  = 20000.0, 0.30

# "For the PDL, the nonlinear hyperelastic model was used based on the double
#  linear stress-strain curve of Vollmer's research: when the dependent variable
#  of PDL e < 7.5%, E1 = 0.05 MPa, and when e > 7.5%, E2 = 0.22 MPa"
#                                                          [Prog Orthod 2023]
PDL_E1, PDL_E2, PDL_EPS_KNEE = 0.05, 0.22, 0.075
PDL_NU = float(os.environ.get('OFJ_PDL_NU', 0.45))
# Experiment knob (default 1.0 = the specified model): scales the PDL secant modulus of the
# listed teeth, a stand-in for the anchorage a larger root area provides -- the ligament of a
# bigger root is stiffer in the same proportion at equal thickness. Used to test whether the
# palatal-vs-buccal ordering follows molar anchorage (FINDINGS: other-patient replication).
PDL_E_SCALE = float(os.environ.get('OFJ_PDL_E_SCALE', 1.0))
PDL_E_SCALE_FDI = os.environ.get('OFJ_PDL_E_SCALE_FDI', '17,16')
# ASSUMED - neither paper states Poisson's ratio for PDL.  Overridable so the
# assumption can be tested against the paper's reported step-10 tipping.

# "The PDL was simulated with a 0.30-mm shell element"      [Prog Orthod 2023]
PDL_THICKNESS = 0.30

# "The CA was developed by making an external offset of dental crowns with a
#  thickness of 0.7 mm"                                     [Prog Orthod 2023]
ALIGNER_THICKNESS = 0.70

# Teeth and alveolar bone are RIGID: "the tooth and the alveolar bone were
# assumed to be rigid bodies during the initial displacement since their
# deformation is negligible compared with that of the PDL"   [Prog Orthod 2023]
# The bone is therefore not meshed at all; the outer PDL surface is constrained.

# -------------------------------------------------------------- attachments --
# "Conventional vertical and horizontal rectangular attachments were designed
#  according to the manufacturer's suggestion."  Fig.1a / Fig.4 of both papers
# show five attachments, on 17,16 (horizontal) and 15,14,13 (vertical).
ATTACHMENTS = {                      # unn: (label, length x height x depth mm)
    2: ('horizontal rectangular', (3.0, 2.0, 1.0)),   # 17
    3: ('horizontal rectangular', (3.0, 2.0, 1.0)),   # 16
    4: ('vertical rectangular',   (2.0, 3.0, 1.0)),   # 15
    5: ('vertical rectangular',   (2.0, 3.0, 1.0)),   # 14
    6: ('vertical rectangular',   (2.0, 3.0, 1.0)),   # 13
}

# ------------------------------------------------------------------ contact --
# "The interaction between the aligners and crowns was defined as a small-sliding
#  surface-to-surface contact, with a friction coefficient set at 0.2"
FRICTION = float(os.environ.get('OFJ_MU', 0.20))   # override for sensitivity only

# "stress and displacement converged when the element size was smaller than 0.2 mm"
TARGET_ELEM_SIZE = 0.20

# ------------------------------------------------------------ anchorage arm --
# "(1) buccal TAD group: the TAD was located in the buccal interradicular space
#  between the first and second molars 4 mm above the alveolar crest ...
#  The application of elastics on the CAs were all set at the buccal or lingual
#  mesial cervical region of the canine to simulate the precision cut.
#  ... the elastic forces for all the test groups were set as 150 g"
GROUPS = {
    'control':    dict(elastic=False),
    'buccal_tad': dict(elastic=True, side='buccal',
                       tad_between=(3, 2),       # between 16 and 17
                       tad_above_crest=4.0,
                       hook_tooth=6,             # canine 13
                       hook_site='mesial cervical, buccal'),
    'palatal_tad': dict(elastic=True, side='palatal',
                        tad_between=(3, 2), tad_above_crest=4.0,
                        hook_tooth=6, hook_site='mesial cervical, palatal'),
}
# Study variants, off by default: the precision cut on another tooth (FDI number) and the elastic force.
if os.environ.get('OFJ_HOOK_FDI'):
    for _g in ('buccal_tad', 'palatal_tad'):
        GROUPS[_g]['hook_tooth'] = UNN[int(os.environ['OFJ_HOOK_FDI'])]
ELASTIC_FORCE_GF = float(os.environ.get('OFJ_ELASTIC_GF', 150.0))
GF_TO_N = 9.80665e-3
ELASTIC_FORCE_N = ELASTIC_FORCE_GF * GF_TO_N      # 1.4710 N

# ------------------------------------------------------------------ staging --
# "The common 'V pattern' staging strategy for maxillary whole arch distalization
#  was designed. During the total 70 steps, for the movement of the canine to the
#  second molar, 0.1 mm distal movement along the arch form was prescribed for the
#  target teeth; for the movement of the incisors, 0.15 mm palatal movement was
#  designed in each relevant step."       Schedule read off Fig.4 of Bioeng 2024.
N_STEPS = 70
STEP_DISTAL = 0.10        # mm per step, canine .. second molar, along the arch form
STEP_PALATAL = 0.15       # mm per step, incisors, palatally

V_PATTERN = {             # unn: (first step, last step)  inclusive, 1-based
    2: (1, 20),           # 17
    3: (11, 30),          # 16
    4: (21, 40),          # 15
    5: (31, 50),          # 14
    6: (41, 60),          # 13
    7: (51, 70),          # 12
    8: (51, 70),          # 11
}

# Study variants, off by default. OFJ_STAGING: 'pairs' moves 17+16, then 15+14, then 13, the incisors
# keeping their V-pattern timing (70 steps); 'seq' moves one tooth at a time with no overlap (120 steps).
# OFJ_OVERCORR: planned extra distal movement of the molars (17, 16) as a fraction, 0.2 = 2.4 mm.
STAGINGS = {
    'v': V_PATTERN,
    'pairs': {2: (1, 20), 3: (1, 20), 4: (21, 40), 5: (21, 40), 6: (41, 60), 7: (51, 70), 8: (51, 70)},
    'seq': {2: (1, 20), 3: (21, 40), 4: (41, 60), 5: (61, 80), 6: (81, 100), 7: (101, 120), 8: (101, 120)},
}
STAGING = os.environ.get('OFJ_STAGING', 'v')
V_PATTERN = STAGINGS[STAGING]
N_STEPS = max(b for _, b in V_PATTERN.values())
MOLAR_OVERCORR = float(os.environ.get('OFJ_OVERCORR', 0.0))


def step_mm(unn):
    """Prescribed movement per active step for this tooth."""
    if unn in (7, 8):
        return STEP_PALATAL
    return STEP_DISTAL*(1 + MOLAR_OVERCORR) if unn in (2, 3) else STEP_DISTAL


def active_teeth(step):
    """Teeth prescribed to move at a given 1-based step."""
    return [u for u, (a, b) in V_PATTERN.items() if a <= step <= b]


def prescribed_increment(unn, step):
    """mm of prescribed movement for this tooth at this step (0 if not active)."""
    a, b = V_PATTERN[unn]
    if not (a <= step <= b):
        return 0.0
    return step_mm(unn)


def total_prescribed(unn):
    a, b = V_PATTERN[unn]
    return (b - a + 1) * step_mm(unn)


# ------------------------------------------------------- 4D iteration scheme --
# "each step of the CA was accompanied by two iterations of PDL to simulate the
#  clinical situation that the CA was well worn with enough time."
# Bone-remodelling loop of Hamanaka et al. (AJODO 2017;152:601):
#   phase 1 - outer PDL surface constrained, tooth moves elastically under load
#   phase 2 - tooth position retained, PDL restored to its original thickness
PDL_ITERATIONS_PER_STEP = 2

# "Prior to each iteration, the CA from the previous step was removed and
#  regenerated, followed by the application of a best-fit algorithm to match the
#  inner surface of the CA with the dental crowns, simulating a wear-in process."
WEAR_IN_BEST_FIT = True


def summary():
    lines = [f"Paper: {CITATION}", ""]
    lines.append(f"  half arch          FDI {[FDI[u] for u in HALF_ARCH_UNN]}")
    lines.append(f"  aligner            {ALIGNER_THICKNESS} mm, E={E_ALIGNER} MPa, nu={NU_ALIGNER}")
    lines.append(f"  attachments        {len(ATTACHMENTS)}, E={E_ATTACH} MPa")
    lines.append(f"  PDL                {PDL_THICKNESS} mm shell, bilinear "
                 f"E1={PDL_E1}/E2={PDL_E2} MPa at eps={PDL_EPS_KNEE}")
    lines.append(f"  teeth + bone       rigid (bone not meshed; outer PDL constrained)")
    lines.append(f"  contact            small-sliding surface-to-surface, mu={FRICTION}")
    lines.append(f"  elastic            {ELASTIC_FORCE_GF:.0f} gf = {ELASTIC_FORCE_N:.4f} N")
    lines.append(f"  staging            {N_STEPS} steps, V pattern, "
                 f"{STEP_DISTAL} mm distal / {STEP_PALATAL} mm palatal per step")
    lines.append(f"  iterations         {PDL_ITERATIONS_PER_STEP} PDL per CA step "
                 f"= {N_STEPS*PDL_ITERATIONS_PER_STEP} solves")
    lines.append("")
    lines.append("  V pattern (1-based step ranges) and total prescribed movement:")
    for u in HALF_ARCH_UNN:
        a, b = V_PATTERN[u]
        lines.append(f"    FDI {FDI[u]}   steps {a:3d}-{b:3d}   "
                     f"{total_prescribed(u):.2f} mm "
                     f"{'palatal' if u in (7, 8) else 'distal'}")
    return "\n".join(lines)


if __name__ == '__main__':
    print(summary())
