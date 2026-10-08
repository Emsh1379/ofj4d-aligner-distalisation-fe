# FEBioCuDSS: GPU linear solver for FEBio (draft) and the validation plan

**Status: compiled and verified on 2026-09-13** (FEBio 4.13.0, cuDSS 0.8.0, CUDA 12,
Tesla T4). On one step of Open-Full-Jaw patient 8 (472,719 equations) the solver time
fell from 35.9 s (PARDISO) to 7.05 s, the whole run from ~68 s to ~41 s, with the same
12 iterations / 2 reformations and node displacements equal to PARDISO within
1.0e-12 mm. Run it with a configuration file:

```xml
<febio_config version="3.0">
	<default_linear_solver type="cudss"/>
	<import>/path/to/libFEBioCuDSS.so</import>
</febio_config>
```
```bash
LD_LIBRARY_PATH=$CUDSS/lib:/usr/local/cuda/lib64:$HOME/febio/mkl/lib:$HOME/febio/FEBio/build/lib \
  febio4 -cnf febio_cudss.xml -i model.feb
```

Known gap: the plug-in does not update FEBio's linear-solver call counters. The
validation of this project's solver against FEBio (section 2 below) is still to do.

## Why this and what to expect

FEBio has no GPU solver. On Open-Full-Jaw patient 8 (472,719 equations, one load step)
the run took 68 s, 36 s of it in PARDISO (docs/FINDINGS.md section 28). Only the linear
solver moves to the GPU, so the whole-run gain is bounded below 2x on that model;
expect roughly 1.5-1.8x with a T4, more on models where factorisation dominates.

## 1. Build and verify the plug-in

```bash
bash pipeline/build_febio.sh $HOME/febio                       # FEBio with MKL PARDISO
FEBIO=$HOME/febio/FEBio bash pipeline/febio_cudss/build.sh      # the plug-in
```

Select it in the model's control section (FEBio 4 format):
`<solver><linear_solver type="cudss"/></solver>`, or `<linear_solver type="cudss"/>`
inside `<Control>` for the older 2.5 format. Load it with `-import libFEBioCuDSS.so`.

Verification, before trusting it:
1. run one step of `P8_maxilla_tipping.feb` with `pardiso` and with `cudss`;
2. compare node displacements (add a `<logfile><node_data data="ux;uy;uz"/>` output and use
   `tools/febio_p8_summary.py`): the difference must be at round-off level, as cuDSS and
   PARDISO agreed to 1e-11 in this project's own solver;
3. compare iteration counts and time in linear solver.

Things to check in the draft: the header name for the CUDA runtime API, `FE_SDK_VERSION`
and the plug-in registration macro against `FECore/sdk.h`, and whether
`CompactSymmMatrix(1)` stores the upper triangle row-wise (cuDSS `CUDSS_MVIEW_UPPER`) or
the lower triangle column-wise (identical as CSR for a symmetric matrix).

## 2. Validation of this project's solver against FEBio

**Level 1 -- single tooth, identical load (hours).** For each of the seven half-arch teeth:
the same rigid tooth, the same PDL mesh and material (bilinear 0.05 / 0.22 MPa or a
linear 0.05 MPa in both), the outer PDL surface fixed, and the same 1 N force at the
same crown point and direction. Export `build4d/geom.pkl` nodes/tets to a `.feb`
(tooth as a rigid body, PDL solid, fixed outer surface, nodal load on the crown point).
Compare crown and apex displacement, rotation and centre of rotation per tooth with
the Schur-complement analysis of section 15. Pass: within a few percent, same ranking.

**Level 2 -- one aligner staging step (days).** Export the half arch, PDL, attachments
and the 0.70 mm aligner of step 1 with frictional sliding contact (mu 0.2) and the
same prescribed aligner shape; compare per-tooth movement, contact force and centre of
rotation after the first PDL iteration.
