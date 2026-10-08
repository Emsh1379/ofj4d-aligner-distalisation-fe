# Hand-off, 2026-09-17: the replication is finished

**All four patients are complete** (3, 5, 6, 4 at 17/17; 68 runs, no failures, no solver fallbacks),
finished 2026-09-17 17:13 UTC on a fresh studio restored from the 2026-09-16 zips. Results and the
five-patient comparison: `compare_all5.json` next to this package, `python tools/compare_patients.py 8 3 4 5 6`.
Replicates: the site ranking (9/10 patient-group rankings, the tenth a 6 um tie), first-premolar cut +
palatal screw (5/5), miniscrew over control (+0.62..+0.92 mm). Does NOT replicate: palatal > buccal
(2 palatal, 2 buccal, 1 tie) and absolute movement (control molar +0.40 P3 to -0.66 P6).

**Experiment A finished 2026-09-17 19:46 UTC** (`build4d/v2/jobs_pdlscale.txt` in P3 and P5, 8 runs):
scaling the molar PDL modulus 0.7/1.4 at the paper site with the geometry fixed moves the palatal-buccal
difference by only ~0.08 mm, in **opposite** directions in the two patients, against the 0.377 mm gap it
would have to explain - so molar root area is not the cause, despite r = +0.93. What does track it is the
elastic's line of action: the palatal advantage rises monotonically with screw height in all five
patients (2 < 4 < 6 < 8 mm). Full write-up: FINDINGS 48 and `docs/REPLICATION.md`, Table R7.
The knob (`OFJ_PDL_E_SCALE`, `spec.py`) stays at its default 1.0 = the specified model.

**Restoring on a new studio, in this order:**
1. Unzip both zips into `pipeline_data/` as siblings, and re-link the code into each workspace
   (`for P in 3 4 5 6; do for d in ofj ofj4d tools pipeline; do ln -sfn $PWD/$d ../ofj4d_patients/P$P/$d; done; done`).
2. **Reinstall cuDSS** - it is not in the zips and is worth ~2x: `pip install --target $HOME/.pylibs_cudss
   "nvmath-python[cu12]==1.0.0" "nvidia-cudss-cu12==0.8.0.10"`. Without it the solver silently falls back to
   PARDISO on the CPU: ~60 s per step instead of ~35 s at PAR=2. `python pipeline/smoke_test.py` must still
   print `SMOKE TEST PASSED`.
3. `setsid nohup bash build4d/v2/guard.sh >/dev/null 2>&1 &` (`~/.lightning_studio/on_start.sh` and a
   fallback in `~/.zshrc` both do this; `.lightning_studio/.studiorc` is a platform symlink - do not edit it).

**Throughput** (T4, 8 cores, cuDSS, measured 2026-09-17 on the same jobs): PAR=2 35.6 s/step (~202 steps/h),
PAR=3 45 s (~240/h), PAR=4 57 s (~253/h, GPU peak 12.2 GB of 14.5) - `supervisor.sh` now runs PAR=4 and
`queue.sh` caps each job at `nproc/PAR` MKL threads.

**Restart hardening added 2026-09-17** (the studio rebooted ~every 30 min, and once slept 6 h):
`queue.sh` skips a finished job without starting python, so re-walking all phases takes ~1.4 s instead of
~7 min; `guard.sh` kills the supervisor's whole process group (killing only `supervisor.sh` leaves an
orphaned `xargs` pool that respawns workers), has a 30-min stall watchdog and a flock on its single-instance
check; `supervisor.sh` verifies every patient job exited 0 before writing `SUPERVISOR_DONE`. Zips do not
carry empty directories, so `patient_queue.sh` now creates `build4d/v2/{runs,logs}` itself - without that,
patients 5/6/4 failed instantly and the chain falsely reported itself finished.

---

# Hand-off, 2026-09-16 (final)

State of the package after all queued work finished (`SUPERVISOR_DONE v5`, 2026-09-16 09:38 UTC). Details and every
number: `docs/FINDINGS.md` sections 29-46. **Paused 2026-09-16 12:03 UTC (studio released):** the other-patient replication (FINDINGS 47). Patients 3, 5, 6, 4
are built; patient 3 has 2 of 17 runs complete (control, buccal paper site) and 2 at step 22 of 70 (buccal 2 mm, 6 mm),
patients 5, 6, 4 not started. ~6 h of GPU per patient remains (~37 s per step at 2 in parallel on a T4).
**To resume**: restore `pipeline_data/ofj4d_patients/` from `ofj4d_patients_2026-09-16.zip` next to this package, re-link the code
into each workspace (`for P in 3 4 5 6; do for d in ofj ofj4d tools pipeline; do ln -sfn $PWD/$d ../ofj4d_patients/P$P/$d; done; done`,
from this package; the zip holds each workspace's `build4d/` only - enough to run), then
`setsid nohup bash build4d/v2/guard.sh >/dev/null 2>&1 &` (`~/.lightning_studio/on_start.sh` does this on studio start).
Runs resume from their last step. Summarise with `python tools/compare_patients.py`.

## Finished

| work | result | where |
|---|---|---|
| Solver core verified against FEBio 4 (PDL + rigid teeth, Level 1) | agrees to 0.14 % | FINDINGS 29, `tools/febio_level1_*.py` |
| FEBio (GPU linear solver) runs the aligner step with this model's contact law | step 1 FDI 17 within 1.5 %; ~20x slower, not adopted | FINDINGS 30-31, `ofj4d/febio4d.py`, `pipeline/febio_*` |
| Meshes repaired (PDL folded prisms, aligner volume mesh) | `geom_fe.pkl`, `aligner_fe.pkl` | FINDINGS 31, `build4d/` |
| Rigid vs deformable bone | <= 17 um in both groups: rigid bone adequate | FINDINGS 32, 40 |
| Miniscrew placeability on this patient | only 6 and 8 mm placeable on both sides | FINDINGS 33 |
| Production position study (13 runs, repaired meshes, K_N 1000) + PDL stress | complete | FINDINGS 34, `build4d/v2/runs/`, `build4d/v2/stress/` |
| Robustness of the position ranking (40 runs: alpha 0/0.60, PDL nu 0.30/0.49) | ranking holds in all 10 comparisons | FINDINGS 35 |
| Mesh convergence (M1-M4) | production mesh converged (<= 6 um molar) | FINDINGS 36 |
| Elastic force 100-250 gf (18 runs) | linear in force, R2 >= 0.995 | FINDINGS 37 |
| Precision-cut location (12 runs) + its robustness (8 runs) | first premolar best; palatal gain robust, buccal a tie at nu 0.49 | FINDINGS 38, 44 |
| Friction mu 0.10/0.30 (4 runs) | not sensitive (<= 64 um) | FINDINGS 39 |
| Molar overcorrection (12 runs) | molar-only remedy, more relapse | FINDINGS 41 |
| Staging protocol V / pairs / sequential (6 runs) | sequential best with a miniscrew | FINDINGS 42 |
| PDL stress at 250 gf (replays) | +4-6 % peak compression for +67 % force | FINDINGS 45 |
| Force design chart | `out/paper/fig9_nomogram.*`, Table 5c | FINDINGS 46 |
| Draft manuscript, figures, tables | `docs/manuscript.md` -> `out/paper/manuscript.docx` / `.tex` (pandoc, regenerated 2026-09-16); `out/paper/` | FINDINGS 43 |

Summaries: `python tools/compare_tad.py --set <tag>` with tag `v2`, `a000`, `a060`, `nu030`, `nu049`, `f100`/`f200`/`f250`,
`hk14`/`hk12`, `hka000`/`hka060`/`hknu030`/`hknu049`, `mu010`/`mu030`, `oc10`/`oc20`/`oc30`; saved as
`build4d/v2/compare_positions_<tag>.json`. `python tools/paper_figs.py` and `tools/paper_nomogram.py` regenerate the figures.

## Open (not computational runs)

* `tools/paper_figs.py` rewrites `out/paper/tables.md` from scratch; Tables 5c, 6b and S2 were added by hand
  (and mirrored in `out/paper/page.html`). Re-add them after regenerating, or teach the script to write them.
* A second patient; PDL stress contour maps from `build4d/v2/stress/*.npz`; clinical citations marked [verify] in
  `docs/manuscript.md`. See FINDINGS "Open items".

## Re-running

Every job file is resumable: `bash build4d/v2/supervisor.sh` returns immediately for finished jobs. The guard
(`build4d/v2/guard.sh`) and `~/.lightning_studio/on_start.sh` relaunch the chain after preemption; add phases to
`supervisor.sh` and bump its end marker.

## Rebuilding what is not shipped

* FEBio 4 with MKL: `bash pipeline/build_febio.sh $HOME/febio`; plug-ins: `pipeline/febio_cudss/build.sh` (see its header for CUDA 12
  runtime variables), `pipeline/febio_pdl/build.sh`, `pipeline/febio_pairs/build.sh`.
* cuDSS for the custom solver: see `pipeline/requirements.txt`; `pip install tetgen` for `tools/volmesh_aligner.py` (TetGen is AGPL).
* Refined mesh-convergence meshes (not shipped): `bash build4d/meshconv/rebuild.sh`.
* Raw Open-Full-Jaw data (for the bone mesh, miniscrew points): `pipeline/SETUP.md` section 5.

## Environment of the rerun

`OFJ_GEOM=build4d/geom_fe.pkl OFJ_ALIGNER=build4d/aligner_fe.pkl OFJ_K_N=1000 OFJ_PRESCRIBE=tip OFJ_TIP_BLEND=0.40 OFJ_BONE=rigid`,
miniscrew file via `OFJ_TAD_FILE` (defaults to `build4d/tad_bone.pkl`, the paper's site). `build4d/v2/queue.sh` sets these.
