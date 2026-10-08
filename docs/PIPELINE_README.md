# 4D clear-aligner staging FEM — portable benchmark package

Reproduction of **Mao et al., *Bioengineering* 2024;11(1):3** (`doi:10.3390/bioengineering11010003`),
*"Effect of Different Anchorage Reinforcement Methods on Long-Term Maxillary Whole Arch
Distalization with Clear Aligner: A 4D Finite Element Study with Staging Simulation"* —
**buccal TAD group**, built on the maxilla of
[Open-Full-Jaw](https://github.com/diku-dk/Open-Full-Jaw) patient 8.

This package exists to answer two questions on a different machine:

1. **How fast is it here?**
2. **Is a better linear solver available here?**

```bash
pip install -r requirements.txt
python benchmark.py                 # full benchmark, ~3 min on an idle machine
python benchmark.py --probe-only    # just list CPU/GPU/solvers, instant
python benchmark.py --quick         # skip alternative-solver trials
```

> **Run it on an idle machine.** The numbers in `benchmark.py` are wall-clock and
> contend badly. On the reference laptop a contended run reported 24.5 s per
> nonlinear iteration where an idle run gives ~11 s.

---

## Answers from the AWS machine (4-core Xeon / 31 GiB / Tesla T4, 2026-09-12)

**Yes — a much better linear solver was available, and the GPU turned out to be worth
using too, though not for the reason this README expected.** Both are now wired in.

| per nonlinear iteration | AMG+CG (original) | + PARDISO | + GPU assembly |
|---|---|---|---|
| assembly (host -> device) | 1.41 s | 0.68 s | **0.06 s** |
| linear solve | 28.8 s | 0.86 s | **0.36 s** |
| measured, whole loop | 29.6 s | 1.27 s | **0.45 s** |
| measured staging step | ~25 min | 85 s | **~25 s** |
| 70 steps / group | ~30 h | ~1.6 h | **~30 min** |

Both anchorage groups, full 70-step protocol: **~1 h**, down from ~60 h.

The linear solve is now 80% of what is left and is at its floor: MKL's
CG-on-the-previous-factorisation mode gains nothing here, and sparse Cholesky is
supernodal dense GEMM — the one thing a T4 is *worse* at. See `docs/FINDINGS.md`
§16, which also records two optimisations that were aimed at a phantom because
cProfile misattributes asynchronous GPU work.

Every figure above is measured on the real nonlinear loop, not modelled;
`python benchmark.py` reproduces the table.

### 1. The direct solver (`Pardiso` in `ofj4d/solve4d.py`)

The system is SPD, 150,644 dof, 5.74M nnz, and the Cholesky factor is only 25.5M nnz
(~230 MB) — so the direct solve this README predicted would win does win, by 21x on
the linear solve. It is exact, which is why it succeeds where every AMG reuse scheme
failed: there is no hierarchy left to go stale when the active set moves.

### 2. The GPU (`ofj4d/gpu4d.py`) — bandwidth, not flops

The README's reasoning ("CuPy can only accelerate CG — the 14%, not the 75%", and
pre-Pascal FP64 being hopeless) was right about *arithmetic*, and right about the
AMG-era cost split. Measured here, a Tesla T4 against this host:

| | T4 | host | ratio |
|---|---|---|---|
| FP64 arithmetic | 0.25 TFLOP/s | 0.14 TFLOP/s | 1.8x — **useless** |
| FP64 memory bandwidth | 236 GB/s | 6.8 GB/s | **35x** |

Every per-iteration assembly kernel is bandwidth-bound; the factorisation is
flop-bound. So the work is split: **assemble on the device, factorise on the host**,
17 MB over PCIe per iteration. The assembly kernels go 234 ms -> 5 ms (PDL stiffness),
165 ms -> 17 ms (triple product), 50 ms -> 1 ms (diagonal scaling).

The second gain is structural. Assembling on the device means choosing the sparsity
pattern rather than letting scipy derive it from the values: every contact pair is
assembled every iteration, inactive ones contributing a genuine zero. The pattern then
belongs to the *step* instead of to the active set, so PARDISO's symbolic analysis —
0.59 s of its 0.86 s — runs once per step instead of once per iteration: one analysis
and 32 reuses in a typical solve. (On the host this was a dead end, because scipy
prunes explicit zeros; cuSPARSE's SpGEMM does not.)

Peak device memory is 2.2 GiB of the T4's 14.6 GiB.

Reductions use atomics by default. `OFJ_GPU_DETERMINISTIC=1` switches to sorted
segment sums, which make reruns bit-identical at a cost of ~30% of the runtime, to
remove a variation of ~1e-16 — worth having, not worth defaulting to.

### Running it

Both paths are automatic, and both degrade gracefully:

```bash
pip install pypardiso                        # the solver; torch is already present
OFJ_GPU=off       python ofj4d/run4d.py ...  # assemble on the host
OFJ_LINSOLVE=amg  python ofj4d/run4d.py ...  # the original AMG+CG path
```

Verified against the original: tooth movements agree with the host path to
**0.0001 nm** on movements of 105 um, with identical active sets, and with the AMG
path to 1.9 nm — well inside the nonlinear loop's own 1e-3 tolerance. See
`docs/FINDINGS.md` §8 and §10.

Not a win here: SciPy SuperLU (25 s, worse than AMG+CG). Untested, not installed:
CHOLMOD, UMFPACK, PETSc, pyamgx.

---

## Why it was slow, and where to look (the original diagnosis)

The model is 150,644 DOF / 291k tetrahedra — small by FEA standards. The cost is
**not** the problem size; it is that the nonlinear loop rebuilds an AMG hierarchy
**every iteration**, ~26 iterations per solve, 2 solves per staging step, 70 steps
= **~3,600 AMG setups**.

Measured split per nonlinear iteration on the reference machine (idle):

| stage | time | share |
|---|---|---|
| PDL stiffness rebuild + triple product | 1.1 s | 9% |
| contact operator + assembly + scaling | 0.2 s | 2% |
| **AMG setup** | **~10 s** | **~75%** |
| CG (11–36 iterations) | 1.5–4 s | 14% |

**Everything hinges on that AMG setup.** Things already tried on the reference
machine, all measured, all failures:

| attempt | result |
|---|---|
| Reuse the hierarchy every 6 iterations | **worse** — CG grinds against a stale hierarchy |
| Jacobi smoothing instead of `smooth='energy'` | **2.3× worse** (setup 83 s vs 30 s) |
| One fixed hierarchy built with all contact pairs active | **catastrophic** — 400 iterations / 6418 s once the active set shrinks |
| Faster active-set relaxation (0.5 → 0.85) | 3% |
| GPU via CuPy | can only accelerate CG — the 14%, not the 75% |

The contact springs change the spectrum enough that the preconditioner genuinely
must be rebuilt. So the promising directions were:

1. **A sparse direct solver.** At 150k DOF with 5.7M nnz this is very plausible.
   `pypardiso` (Intel MKL PARDISO) or `sksparse` (CHOLMOD) would replace
   setup+CG entirely, and the matrix is SPD. **This is the first thing to try** —
   `benchmark.py` tests both automatically if installed.
2. **`pyamgx`** (NVIDIA AMGX) — the only GPU route that moves the AMG *setup* onto
   the device. Needs a Pascal-or-newer GPU; pre-Pascal consumer cards run FP64 at
   1/32 rate and will not help.
3. **Refactorisation reuse.** The sparsity pattern is constant within a solve; only
   values change. A solver supporting symbolic-factorisation reuse (CHOLMOD's
   `analyze` / `cholesky_inplace`, PARDISO phase 11 once + phase 22 per iteration)
   would amortise most of the cost.

Direction 1 was right and is now implemented. Direction 3 turned out to rest on a false
premise — the pattern is *not* constant within a solve, because the contact active set
changes every iteration and scipy prunes the pruned-out entries. The machinery to reuse
the symbolic factorisation is in place (`Pardiso` in `ofj4d/solve4d.py`) and fires
whenever the pattern does repeat, but in practice it never does. The same idea *did* pay
off one level up, in the assembly: see `docs/FINDINGS.md` §8.

## Contents

```
benchmark.py             speed + solver-availability benchmark  <- start here
ofj4d/                   the 4D staging model
  spec.py                every paper parameter, quoted, with ASSUMED items flagged
  geom4d.py              half-arch geometry, PDL shells, attachments
  aligner4d.py           0.70 mm aligner (voxel closing + offset)
  model4d.py             reduced system: rigid teeth, projection, near-nullspace
  solve4d.py             nonlinear solve: contact active set + bilinear PDL
  arch4d.py              arch form + 'V pattern' staging schedule
  morph4d.py             aligner staging morph + wear-in seating
  loads4d.py             anchorage groups (control / buccal TAD / palatal TAD)
  run4d.py               the staging loop
  inputs.py              slim runtime inputs
  draw_model.py          render the model
  gpu4d.py               device assembly of the reduced system (optional)
tools/                   analysis and diagnostics added for this study
  analyse.py             validate the runs: completeness, convergence, expression,
                         anchorage, and a side-by-side with the paper
  plot_staging.py        per-tooth staging curves, all groups   -> out/
  plot_compare.py        four-panel group comparison            -> out/
  pdl_diag.py            centre of resistance / centre of rotation of one tooth
  verify_cap.py          replay a step at a higher iteration cap
  diag_chain.py          replay the staging chain and dump a solve's full trace
  breakdown.py           synchronised per-iteration timing breakdown
  lag_seat.py            per-tooth lag behind the prescription + wear-in shift
  moment_probe.py        net contact force/moment per tooth, line of action
  att_variant.py         rebuild attachments at a scaled size + re-wrap the aligner
  sensitivity.py         tabulate the 70-step sensitivity runs against the paper
ofj/                     FE core (tet4 elasticity, contact pairing, meshing)
out/                     result figures (staging curves, group comparison)
out4d/                   rendered model overview (draw_model.py)
build4d/run_<group>.pkl  the 70-step results, one per anchorage group
build4d/geom.pkl         half-arch geometry, prebuilt
build4d/aligner.pkl      aligner mesh, prebuilt
build4d/slim.pkl         0.8 kB extract standing in for the 99.6 MB full model
docs/spec.txt            printed parameter summary
docs/FINDINGS.md         verification status and the defects found so far
```

`build/model.pkl` (99.6 MB, the full Open-Full-Jaw derived model) is **not**
included — `build4d/slim.pkl` supplies the handful of values the runtime needs
(two mesiobuccal cusps, four crown centres, one axis pair, two alveolar crest
points), verified to reproduce the full-model values exactly. To rebuild geometry
from scratch you need the Open-Full-Jaw dataset and the `ofj/build_*.py` stages.

## Source data and model variants

The geometry is Open-Full-Jaw patient 8 (CC BY-NC-SA 4.0). The package ships the
derived half-arch files in `build4d/`; the raw data are only needed for the
bone-based steps below.

```bash
# the dataset zip is stored with Git LFS -- fetch the real file, not the pointer
curl -L -o Patient_8.zip https://media.githubusercontent.com/media/diku-dk/Open-Full-Jaw/main/dataset/Patient_8.zip
unzip Patient_8.zip
mkdir -p data/P8
ln -s $PWD/Patient_8/output/maxilla/volumetric_meshes/maxilla_volumetric_mesh.msh data/P8/
ln -s $PWD/Patient_8/input/maxilla/teeth_axes_maxilla.json data/P8/

python tools/tad_from_bone.py                 # miniscrew points on the real bone -> build4d/tad_bone.pkl
python tools/extract_fixed_region.py Patient_8/output/maxilla/simulation/P8_maxilla_tipping.feb fixed.npy
OFJ_BONE_FIXED=fixed.npy python ofj4d/bone4d.py   # right-half bone -> build4d/bone.pkl
```

Two bone models, differing only in whether the bone deforms (FINDINGS 21):

```bash
python ofj4d/run4d.py --group buccal_tad --steps 70                     # rigid bone, as Mao et al.
OFJ_BONE=deformable python ofj4d/run4d.py --group buccal_tad --steps 70 \
    --out build4d/run_buccal_tad_bone.pkl                                # deformable bone
```

Linear solver: `OFJ_LINSOLVE=auto` uses NVIDIA cuDSS on the GPU when installed, else
PARDISO. On a 2-core host it is 6x faster on the deformable system (FINDINGS 22):

```bash
pip install --target /path/to/cudss_libs "nvmath-python[cu12]" nvidia-cudss-cu12
export OFJ_CUDSS_DIR=/path/to/cudss_libs
```

Miniscrews are points: on the cortical surface in the 16/17 interradicular space,
4 mm apical to the crest on that side (FINDINGS 20). `OFJ_TAD_RULE=legacy` restores
the earlier septum-ridge placement.

Sensitivity and diagnostic knobs (FINDINGS 24). Defaults reproduce the shipped model;
every variant should write to its own `--out`:

```bash
OFJ_PDL_NU=0.30 | OFJ_K_N=4000 | OFJ_MU=0          # ASSUMED PDL ratio, contact penalty, friction
OFJ_E_ALIGNER_SCALE=0.5                              # aligner modulus (diagnostic)
OFJ_NO_ATTACH=1                                      # attachments out of the contact surface
OFJ_SEAT=trans|none                                  # wear-in fit: translation only / none
OFJ_PRESCRIBE=tip  OFJ_TIP_PIVOT=0.83                # sockets tipped about the apical third
OFJ_PRESCRIBE=tip  OFJ_TIP_BLEND=0.5                 # between: 0 bodily .. 1 tipping (tools/calibrate_blend.py)
OFJ_TIP_BLEND_INCISOR=0                              # separate blend for the incisors (default: OFJ_TIP_BLEND)
python tools/paper_frame.py build4d/final/run_       # movement in the paper's global frame
bash pipeline/build_febio.sh $HOME/febio             # FEBio from source with MKL, for cross-validation
python tools/att_variant.py 1.25 1.25 A125           # attachment size variant, then run with
OFJ_GEOM=build4d/variants/geom_A125.pkl OFJ_ALIGNER=build4d/variants/aligner_A125.pkl ...
OFJ_SHELL_LAYERS=4 python ofj4d/aligner4d.py out.pkl # aligner through-thickness refinement
```

## Running the actual simulation

```bash
python ofj4d/run4d.py --group buccal_tad --steps 30
python ofj4d/run4d.py --group control    --steps 30     # for the comparison
```

Checkpoints every step to `build4d/run_<group>.pkl` and resumes automatically.
Run groups **sequentially** unless you have >8 GB free — concurrent AMG setup
spikes OOM'd on a 16 GB machine twice. (On the PARDISO path the peak is far lower —
~230 MB for the factor plus ~0.5 GB transient while the assembly map is built — but
the two groups were still measured one at a time.)

Measured here: **~30 s per step** with the GPU assembly (85 s without), so
`--steps 30` is ~16 min per group and the full 70-step protocol ~37 min.

`--steps 30` is a deliberate truncation: both molars (FDI 17, 16) complete their
prescribed 2 mm by step 30, which covers the paper's step-10 and step-20
anchorage comparisons. 70 steps is the full protocol including premolar, canine
and incisor phases.

## Model at a glance

| | |
|---|---|
| Dentition | right half arch, FDI 17→11, midsagittal symmetry |
| Teeth + bone | **rigid** (bone not meshed; outer PDL surface constrained) |
| PDL | 0.30 mm shell, 3 element layers, bilinear E₁=0.05 / E₂=0.22 MPa at ε=7.5% |
| Aligner | 0.70 mm, E=1500 MPa, generated by rolling-ball closing + offset |
| Attachments | 5 — horizontal on 17/16, vertical on 15/14/13, E=20 GPa, rigid with tooth |
| Contact | node-to-surface, pairing fixed per step, Coulomb μ=0.2 |
| Buccal TAD | 16/17 buccal interradicular space, 4 mm apical to the alveolar crest |
| Elastic | 150 gf to a precision cut at the canine mesial cervical region |
| Staging | 70 steps, 'V pattern', 0.1 mm/step distal, 0.15 mm/step palatal incisors |
| Iteration | 2 PDL (bone-remodelling) iterations per aligner step |
| Linear algebra | GPU assembly + MKL PARDISO; falls back to host, then to AMG+CG |

Units throughout: **N, mm, MPa**.

## Verification status

Passing:

- unstaged aligner is **exactly passive** (0.0000 N)
- a tooth at its prescribed position feels **0.00000000 mm interference on zero
  nodes**, at every step — the invariant that bounds the force
- staged interference equals the prescribed movement (0.1006 mm at step 1)
- single-tooth PDL stiffness 13–16 N/mm pre-knee, 58–72 N/mm post-knee — no
  volumetric locking
- per-tooth force balance (PDL vs contact) within ~1 N
- PDL mesh: median shape quality 0.495; 7 of 167,751 elements near-degenerate
  (report p99 strain, not max — those slivers dominate the maximum)

See `docs/FINDINGS.md` for the defects found and fixed, which is worth reading
before trusting any output.

## Citation

Geometry: Gholamalizadeh et al., *Open-Full-Jaw: An open-access dataset and
pipeline for finite element models of human jaw*, Comput Methods Programs Biomed
224 (2022) 107009 — patient 8, maxilla.

Protocol: Mao B, Tian Y, Xiao Y, Liu J, Liu D, Zhou Y, Li J. *Bioengineering*
2024;11(1):3. Method: Mao et al., *Prog Orthod* 2023;24:16.
Bone remodelling: Hamanaka et al., *AJODO* 2017;152:601.
