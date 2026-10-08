# Final pipeline: setting up a new machine

4D finite-element staging simulation of clear-aligner maxillary whole-arch
distalisation, reproducing Mao et al., Bioengineering 2024;11(1):3, on the right
half maxilla of Open-Full-Jaw patient 8. The final model uses blended socket
kinematics (alpha 0.40), calibrated on the control group only (`docs/FINDINGS.md`
sections 24-25).

## 1. What is in this package

| path | contents |
|---|---|
| `ofj4d/`, `ofj/` | the model and the FE core (tet4 elasticity, contact, meshing helpers) |
| `tools/` | analysis, calibration, sensitivity, figures, diagnostics |
| `pipeline/` | this guide, `requirements.txt`, `smoke_test.py`, `run_final.sh`, `run_sensitivity.sh` |
| `build4d/geom.pkl`, `aligner.pkl`, `slim.pkl`, `tad_bone.pkl` | prebuilt model inputs: everything the rigid-bone model needs |
| `build4d/final/run_<group>.pkl` | **final results** (alpha 0.40), 70 steps, three groups |
| `build4d/run_<group>.pkl` | bodily-socket model (alpha 0), 70 steps, three groups |
| `build4d/variants/` | calibration runs (alpha 0.6, alpha 1), sensitivity checkpoints (resumable) |
| `out/`, `out4d/` | figures |
| `docs/FINDINGS.md` | every verification finding, the diagnosis, the calibration and its limits |

**Not included** (download or rebuild only if needed):

* Open-Full-Jaw raw data. Needed only to rebuild the bone, the bone-based miniscrew
  points or the TAD figure: see section 5.
* `build4d/bone.pkl` (deformable-bone variant): rebuilt from the raw data by `ofj4d/bone4d.py`.
* Attachment-size variants and refined aligners: rebuilt automatically by
  `pipeline/run_sensitivity.sh` (`tools/att_variant.py`, about a minute each).

## 2. Install

```bash
unzip ofj4d_final_pipeline.zip && cd ofj4d_final_pipeline
python3 -m venv .venv && source .venv/bin/activate      # or a conda env, Python 3.10-3.12
pip install -r pipeline/requirements.txt
```

GPU (optional, 3-6x faster): follow the two commands at the end of
`pipeline/requirements.txt` (a CUDA build of PyTorch, then cuDSS in `$HOME/.pylibs_cudss`).
Without a GPU everything runs on the CPU; `OFJ_GPU=off` forces that.

## 3. Check the installation (1-5 minutes)

```bash
python pipeline/smoke_test.py
```

Runs two staging steps of the final model and compares all seven teeth with the
shipped reference. It must print `SMOKE TEST PASSED`.

## 4. Run

```bash
bash pipeline/run_final.sh 2          # three groups, 2 at a time, then analysis + figures
bash pipeline/run_sensitivity.sh 2    # PDL Poisson ratio, contact stiffness, attachment size
```

Both are resumable: after an interruption run the same command again. The final
results are already shipped, so `run_final.sh` only recomputes them if you delete
`build4d/final/`. Measured cost: about 25 s per step with one run on a T4 (40-45 s
each with two), about 85 s per step on a 4-core CPU without a GPU; 70 steps per run.

Individual runs and the models:

```bash
python ofj4d/run4d.py --group control --steps 70 --out build4d/x.pkl             # bodily sockets (alpha 0)
OFJ_PRESCRIBE=tip OFJ_TIP_BLEND=0.40 python ofj4d/run4d.py --group palatal_tad \
    --steps 70 --out build4d/y.pkl                                                # final model
python tools/analyse.py build4d/final/run_        # completeness, convergence, against the paper
python tools/calibrate_blend.py                   # alpha calibration on the control group
python tools/sensitivity.py                       # sensitivity table
```

External cross-validation (optional): `bash pipeline/build_febio.sh $HOME/febio` builds
FEBio from source with Intel MKL PARDISO (about 10-20 min on 4 cores; MKL is fetched
with pip and unpacked by the script) and writes `$HOME/febio/febio.sh -i model.feb`.
Verified 2026-09-13: FEBio 4 (commit dbbcf7b) with MKL 2024.2.2 built in 12 min on
4 cores; one time step of Open-Full-Jaw's `P8_maxilla_tipping.feb` (791,127 tet4)
ran to normal termination in 68 s, 36 s of it in the PARDISO linear solver. MKL 2025+
fails to compile FEBio's ILU0 preconditioner (`mkl_dcsrtrsv`), hence the pin.

Groups: `control`, `buccal_tad`, `palatal_tad`. Knobs (all optional; defaults reproduce
the shipped model): `OFJ_PRESCRIBE`, `OFJ_TIP_BLEND`, `OFJ_PDL_NU`, `OFJ_K_N`, `OFJ_MU`,
`OFJ_E_ALIGNER_SCALE`, `OFJ_NO_ATTACH`, `OFJ_SEAT`, `OFJ_GEOM`, `OFJ_ALIGNER`,
`OFJ_BONE=deformable`, `OFJ_LINSOLVE=auto|cudss|pardiso|amg`, `OFJ_GPU=off`,
`OFJ_CUDSS_DIR` (see `README.md`).

## 5. Raw data (optional)

```bash
curl -L -o Patient_8.zip https://media.githubusercontent.com/media/diku-dk/Open-Full-Jaw/main/dataset/Patient_8.zip
unzip Patient_8.zip
mkdir -p data/P8 data/cache
ln -s $PWD/Patient_8/output/maxilla/volumetric_meshes/maxilla_volumetric_mesh.msh data/P8/
ln -s $PWD/Patient_8/input/maxilla/teeth_axes_maxilla.json data/P8/
```

Then `README.md`, section "Source data and model variants". Mesh caches default to
`data/cache/` (`OFJ_MSH_CACHE`, `OFJ_BONE_FIXED` override). Open-Full-Jaw is
CC BY-NC-SA 4.0; the derived geometry in `build4d/` inherits that licence.

## 6. Final results shipped (alpha 0.40)

| | paper | model |
|---|---|---|
| FDI 17 at step 10: control / buccal / palatal | 0.97 / 1.09-1.13 / 1.09-1.13 mm | 0.93 / 1.04 / 1.06 mm |
| incisor retraction at step 70: control / buccal / palatal | 0.87 / 1.71 / 1.88 mm | 1.24 / 1.88 / 2.06 mm |
| group ordering | palatal > buccal > control | palatal > buccal > control |
| FDI 17 / 16 at step 70, palatal | 1.31 / 1.00 mm | 1.02 / 1.50 mm |
| FDI 17 tipping at step 10 | 5.92-6.87° | 2.2-2.4° |

alpha is calibrated to the paper, not measured. `docs/FINDINGS.md` section 25 gives
the out-of-sample test, the literature on the movement type, and what the calibration
does not fix.
