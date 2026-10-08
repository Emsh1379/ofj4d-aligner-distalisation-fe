# Verification findings

Read this before trusting any output. Four defects were found and fixed; each one
silently produced plausible-looking but wrong numbers.

## 1. Staged interference was clamped to zero (fatal)

`gref = max(0, gap - ADJUST)` was carried over from an earlier force-driven
full-arch model, where an external 250 gf load did the work and the clamp
correctly suppressed spurious pull-in. In a **staging** model the aligner-tooth
interference *is* the load, so the clamp deleted the entire driving mechanism.

Symptom: the control group produced **exactly 0.000 N** and no tooth movement;
the buccal TAD group was driven only by its 150 gf elastic.

## 2. The gap reference was pairing-inconsistent (fatal)

The reference standoff was measured at the step-0 configuration, but the current
gap uses the current pairing. Any slave node that re-paired to a different master
triangle had its reference subtracted from a gap measured against a *different
surface*.

Symptom: one contact pair carried **93 N** against a 0.5 N mean; per-tooth force
balance was out by **113 N**; net force on the unloaded control dentition was
117 N instead of ~0.

## 3. The reference drifted with sliding (fatal)

Even with consistent pairing, anchoring to the original configuration made the
reference decay as the appliance slid 0.2, 0.4, 0.6 mm over the crowns, so force
grew with *distance travelled* rather than with interference.

Symptom: contact force climbed monotonically 4 -> 19 N over six steps while the
tooth tracked its prescribed position essentially perfectly — physically
impossible, since a tooth at its prescribed position must feel nothing.

Fix: the pairing is fixed per step and referenced to the configuration the
appliance was **manufactured** for (teeth at prescribed positions).

## 4. The staged motion was applied twice

`morph.shape(step)` already places the aligner at the prescribed positions; the
wear-in `seat()` then moved it *again* by a rest->actual rigid fit.

Symptom: teeth placed exactly at their prescribed positions still showed
-0.034 mm interference on 5,441 nodes, growing to -1.07 mm by step 30.

Fix: `seat()` fits prescribed->actual. The invariant now holds exactly.

## 5. Non-fatal: the strain "explosion" was a reporting artifact

Peak PDL equivalent strain of 424% came from **7 elements out of 167,751**
(0.004% of the mesh, 0.000% of the volume) at FDI 17's root apex, where the
constant-thickness offset shell thins. p99 strain is 6.8%. The runner now
reports p99.

## 6. Non-fatal: AMG near-nullspace was wrong

`P.T @ rigid_body_modes` is not the near-nullspace of a reduced system —
transpose is not pseudo-inverse. Nullspace quality was 1.5e8; the exact reduced
rigid modes give 1.1e4. CG went from **600+ iterations unconverged to 71
converged**.

## 7. Non-fatal: reported contact force included stale values

`lam` was returned after relaxation, so pairs that had separated retained decayed
values (0.5^k, never zero): 6,878 pairs appeared loaded when 1,560 were active.
`solve_config` now returns the true penalty force, zero off the active set.

## 8. Non-fatal: the linear solver was 95% of the wall time

Measured on a 4-core AWS Xeon (31 GiB, Tesla T4 idle), 150,644 dof / 5.74M nnz:

| per nonlinear iteration | before | after |
|---|---|---|
| AMG setup | 10.4 s | — |
| CG (36 iterations, residual 1e-6) | 18.4 s | — |
| PARDISO (exact, residual 1e-15) | — | 0.86 s |
| PDL stiffness rebuild | 0.96 s | 0.24 s |
| rest of assembly | 0.45 s | 0.44 s |
| **total** | **30.2 s** | **2.06 s** |

Two changes, both in `ofj4d/solve4d.py`:

**(a) MKL PARDISO replaces AMG+CG** (`Pardiso`, selected by `OFJ_LINSOLVE=auto`,
falling back to AMG+CG when `pypardiso` is absent). The reduced system is SPD, so
mtype 2 applies and only the upper triangle is passed; the Cholesky factor is
25.5M nnz / ~230 MB, i.e. the fill-in that would have made this hopeless at scale
simply is not there at 150k dof. Phase split: 0.59 s analysis, 0.20 s numeric
factorisation, 0.07 s back-substitution.

The `AMG_EVERY` reuse scheme and every other preconditioner-reuse attempt failed
because a stale hierarchy makes CG grind. A direct factorisation has no analogue of
that failure: it is exact every time, whatever the active set does.

**(b) The assembly writes into a fixed CSR pattern.** `Operators.K` re-sorted 24.2M
COO entries into CSR on every iteration, although only the per-element modulus
changes. The COO -> CSR map is now built once per configuration (2.0 s) and each
rebuild is a scatter-add: 0.96 s -> 0.24 s, identical pattern, values differing only
by summation order (3.6e-16 relative).

**Verification.** One full solve of step 1 (buccal TAD), both backends, same
configuration:

| | AMG+CG | PARDISO |
|---|---|---|
| tooth rigid motions | — | **1.9 nm** different, on movements of 105 um |
| contact force | 2.986 N | 2.992 N (0.2%) |
| active set | 2127 pairs | 2126 pairs, 7 differ |
| linear residual | 1e-6 | **1e-15** |

Every difference is far inside the nonlinear loop's own convergence tolerance
(`U_TOL = 1e-3`), and it is the PARDISO solution that is the exact one — the
differences are CG's truncation, not PARDISO's error. End-to-end, two staging steps
reproduce the documented invariants unchanged (0.105 mm interference at step 1,
p99 strain 6.8%).

**What did not work.** Reusing the symbolic factorisation across iterations (PARDISO
phase 11 once, then phase 22+23 per iteration, 0.86 s -> 0.29 s) is implemented and
is checked for on every solve, but it never fires: the contact active set changes on
every single iteration, so the matrix pattern changes too. Forcing a constant pattern
by assembling all 10,470 pairs with zeros off the active set does not work either --
scipy prunes explicit zeros in both the triple product and the sparse addition, so
the pattern still tracks the active set. It could only be forced with an epsilon
instead of a zero, which is not worth a numerical fiction for 0.5 s.

## 9. Non-fatal: the 35-iteration cap is a limit cycle, not a truncation

`it=35` (reported as `CAP`) appears from step 2 onward. Raising `MAXIT` to 90
(`OFJ_MAXIT=90`, now honoured) does **not** clear it -- the second step still runs to
the cap. The active-set loop is not converging slowly, it is oscillating between two
nearly-equal active sets.

At step 2 it costs nothing: those 55 extra iterations change the computed tooth
movement by **0.02 nm** out of 204.5 um (1e-7 relative). The physical answer is
settled long before the cap; what oscillates is a handful of marginal contact pairs.
`CAP` in the log should be read as "the active set is hunting", not "this result is
truncated".

**That conclusion does not generalise to every capped solve** -- see §11, which
measures it across a full 70-step run. In the molar phase a few solves stop while
still far from stationary, and there the cap does move the answer, by up to 0.3%.

## 10. Non-fatal: the GPU is worth using, for bandwidth rather than arithmetic

The README ruled the GPU out on the grounds that CuPy could only accelerate CG, and
that pre-Pascal FP64 runs at 1/32 rate. Both statements are true and neither settles
the question, because the per-iteration cost is not arithmetic. Measured, Tesla T4
against this host:

| | T4 | host | ratio |
|---|---|---|---|
| FP64 arithmetic | 0.25 TFLOP/s | 0.14 TFLOP/s | 1.8x |
| FP64 memory bandwidth | 236 GB/s | 6.8 GB/s | **35x** |

Sparse assembly is bandwidth-bound and sparse Cholesky is flop-bound, so the split is
assembly on the device, factorisation on the host (`ofj4d/gpu4d.py`), 17 MB over PCIe
per iteration:

| kernel | host | T4 |
|---|---|---|
| PDL stiffness scatter-add (24.2M entries) | 234 ms | 5 ms |
| triple product P.T K P (SpGEMM) | 165 ms | 17 ms |
| diagonal scaling | 50 ms | 1 ms |
| **whole nonlinear iteration, measured** | **1.27 s** | **0.44 s** |

The larger gain is structural, and it retracts the dead end recorded in §8. Assembling
on the device means the pattern is chosen rather than inferred from the values: every
contact pair is assembled every iteration, inactive ones contributing a genuine zero,
so the pattern belongs to the staging step and not to the active set. scipy prunes
explicit zeros — which is what killed this on the host — but cuSPARSE's SpGEMM derives
its output pattern from the input patterns alone, verified stable across active sets
differing by 1000 pairs. PARDISO's symbolic analysis therefore runs once per step
rather than once per iteration: **one analysis and 32 reuses** in a typical solve,
0.86 s -> 0.29 s per linear solve.

No epsilon is substituted for any zero. The explicit zeros are structural zeros, the
ordinary FEM kind, and the assembled matrix agrees with the host's to 1.1e-15.

**Verification.** Same configuration, host vs device assembly, both with PARDISO:
identical iteration counts, **identical active sets** (0 of 2126 pairs differ), tooth
rigid motions agreeing to 7.4e-11 relative — 0.0000 nm on movements of 105 um — and
contact force to 8e-10. Over three staging steps the two paths agree on contact force
to six decimals (3.516056 N vs 3.516056 N at step 1) and on crown movement to less
than 0.0001 nm. The device path is the same arithmetic in a different summation order.

**Reproducibility.** Reductions use sorted segment sums rather than atomics, so runs
are bit-reproducible; `OFJ_GPU_DETERMINISTIC=0` selects the faster atomic path at the
cost of ~1e-16 run-to-run variation. Peak device memory is 2.2 GiB.

**What is left.** Of the remaining 0.44 s per iteration: PARDISO 0.26 s (flop-bound,
and a GPU sparse Cholesky on a T4 would not beat it), the two SpGEMMs of the triple
product 0.16 s, and host-side einsums for strain and contact kinematics ~0.15 s. The
triple product is the next real target: P is a projection whose free-node block is an
identity, so P.T K P is mostly a submatrix gather rather than a product, and the
reduced matrix's pattern is already fixed for the step.

## 11. The 70-step production run: what converged and what did not

Both anchorage groups, full 70-step protocol, GPU assembly + PARDISO, 38 min and
37 min. 140 steps, 280 solves, **zero errors, zero non-finite values, all 70 steps
and all 7 teeth present in both checkpoints**.

Convergence is not uniform, and the summary statistics mislead unless split:

| | buccal TAD | control |
|---|---|---|
| solves at the 35-iteration cap | 52 of 140 | 38 of 140 |
| final du > 1e-3 (the tolerance) | 50 | 34 |
| final du > 1e-2 (materially unconverged) | **4** | **4** |
| steps affected | 12, 15, 16, 19 | 13, 14, 17, 18 |

All eight sit in the second-molar phase and all are the *first* of the step's two
PDL iterations; the second iteration of those steps converged normally. Two were
replayed from their checkpoints at a 250-iteration cap:

| | buccal TAD step 16 | control step 17 |
|---|---|---|
| converges by 250 iterations? | no (still capped) | nearly (du 0.15 -> 0.005) |
| change in tooth movement vs cap 35 | **0.09 nm** | **3.3 um** |
| relative to the movement | 7e-8 | 2.7e-3 |

So the cap is harmless in one case and worth 0.27% in the other. The results below
are unaffected at that magnitude, but the honest statement is "tooth positions carry
an uncertainty of a few um at the eight affected steps", not "the run converged".

Two diagnostics in the checkpoints are **not** physical and should not be quoted:
peak contact force of 320 N (TAD) and 722 N (control). Both belong to capped first
iterations, and replaying control step 17 from its checkpoint gives 13 N where the
run recorded 722 N. They are transients of a chattering active set at the instant
the cap fell, not forces the dentition ever carried. Median contact force is
6.0-7.0 N.

### Results

Expression at the end of each tooth's own active phase, before later relapse:

| tooth | phase | prescribed | buccal TAD | control |
|---|---|---|---|---|
| FDI 17 | 1-20 | 2.00 mm | 1.606 mm (80%) | 1.386 mm (69%) |
| FDI 16 | 11-30 | 2.00 mm | 1.723 mm (86%) | 1.448 mm (72%) |
| FDI 15 | 21-40 | 2.00 mm | 1.266 mm (63%) | 0.997 mm (50%) |
| FDI 14 | 31-50 | 2.00 mm | 0.993 mm (50%) | 0.614 mm (31%) |
| FDI 13 | 41-60 | 2.00 mm | 0.471 mm (24%) | -0.161 mm (-8%) |
| FDI 12 | 51-70 | 3.00 mm | 0.841 mm (28%) | 0.187 mm (6%) |
| FDI 11 | 51-70 | 3.00 mm | 0.829 mm (28%) | 0.169 mm (6%) |

This **supersedes the open item** recording a 90-100% expression rate: that figure
came from 30-step truncations, in which only the molars had ever been loaded. Over
the full protocol expression falls monotonically down the arch, 80% at the second
molar to 28% at the incisors, which brackets the paper's reported 37-68% rather
than contradicting it.

The buccal TAD is better at every tooth, and the margin widens as the active tooth
moves forward: +11 points at FDI 17, +19 at FDI 14, +32 at the canine, +22 at the
incisors. At the end of the run the TAD group holds net distalisation on every
tooth (+0.09 to +0.84 mm) while the control group is net **mesial** on all five
posterior teeth (-0.01 to -0.72 mm): it has lost the whole distalisation and
drifted forward past its starting position.

`out/staging_curves.png` plots all seven teeth against their prescription.

## 12. Comparison with the paper's published results

Numbers extracted from the open-access full text (PMC10813679); verify against the
PDF before citing. Note the paper's step-70 efficiencies are quoted for the
**palatal** TAD group, while this package models the **buccal** TAD group, so only
the control column is directly comparable at step 70.

**Agrees:**

| quantity | paper | this run |
|---|---|---|
| FDI 17 distal, step 10, control | 0.97 mm | 0.840 mm (-13%) |
| FDI 17 distal, step 10, buccal TAD | 1.09-1.13 mm | 0.959 mm (-13%) |
| buccal TAD : control ratio at step 10 | 1.12-1.16 | **1.14** |
| FDI 17 mesial drift, steps 20-40, control | 0.89 mm | 0.752 mm |
| FDI 17 mesial drift, steps 20-40, TAD groups | 0.61-0.76 mm | 0.549 mm |
| ordering (TAD better than control) | everywhere | everywhere |

The mid-protocol relapse is the notable one: the paper reports "undesired mesial
tipping ... in the already distally shifted second molars" over steps 20-40, worst
in the control group, and this run reproduces both the effect and its magnitude
without having been tuned for it.

**Differs:**

| quantity | paper | this run |
|---|---|---|
| FDI 17 tipping at step 10, control | 5.92° | 1.10° |
| FDI 17 tipping at step 10, buccal TAD | 6.63-6.87° | 1.36° |
| FDI 16 distal at step 70, control | +0.62 mm (31.96%) | **-0.370 mm** |
| FDI 17 efficiency at step 70, control | 37.63% | **-35.9%** |
| FDI 11 retraction at step 70, buccal TAD | 1.71 mm | 0.829 mm |
| FDI 11 retraction at step 70, control | 0.87 mm | 0.169 mm |

Two distinct discrepancies:

1. **Movement type.** The paper's molars tip 5-6x more than this model's for the
   same crown translation. This model does tip (root moves half as far as the
   crown, ratio 0.50 at step 10) but far more bodily. Something makes the teeth
   here much stiffer in rotation than in the paper's model -- a candidate is the
   attachment/PDL interaction, or the rigid-body rotation being reported about a
   different reference than the paper's "tilting".

2. **Cumulative retention.** Through step 20 the agreement is within 13%. Past
   step 40 this model keeps relapsing where the paper's arrests: the paper's
   control keeps +0.62 to +0.75 mm at step 70 where this run ends net *mesial*.
   The relapse rate over steps 20-40 matches; what differs is that it does not
   stop.

The eight materially unconverged solves (§11) are worth 0.3% and cannot explain a
sign change. The documented candidates remain: a real unaligned Open-Full-Jaw
patient-8 arch against the paper's idealised dentition, and the ASSUMED PDL
Poisson ratio and von Mises strain driver. The movement-type discrepancy (1) is
the more promising thread, since a model that tips less transfers reciprocal load
to the anterior segment differently, which is exactly what the late phase depends
on.

**Bottom line:** the early mechanics and every group-to-group comparison reproduce;
the long-horizon absolute movements do not. Relative comparisons between anchorage
groups are the trustworthy output, as the open items already warned.

## 13. Attempts to close the gap with the paper, and what they ruled out

Three hypotheses for the movement-type discrepancy of §12, tested against the
paper's step-10 benchmark (FDI 17: 0.97 mm crown, 5.92°, centre of rotation at the
apical third). All are now exposed as environment knobs.

| variant | crown | root | root/crown | tipping | capped solves | Fc |
|---|---|---|---|---|---|---|
| baseline (nu 0.45, K_T 2000, no clearance) | 0.840 mm | 0.450 mm | 0.54 | 1.10° | 5/20 | 7.0 N |
| `OFJ_PDL_NU=0.30` | 0.745 mm | 0.535 mm | 0.72 | 0.60° | - | - |
| `OFJ_K_T=200` | 0.840 mm | 0.450 mm | 0.54 | 1.10° | - | - |
| `OFJ_CLEARANCE=0.03` | 0.807 mm | 0.115 mm | 0.14 | 1.94° | 19/20 | 2.2 N |
| `OFJ_CLEARANCE=0.045` | 0.791 mm | 0.020 mm | **0.03** | 2.16° | 20/20 | **104.6 N** |
| `OFJ_CLEARANCE=0.06` | — singular — | | | | | |

**PDL Poisson's ratio is not the cause.** 0.30 moves the model *away* from the
paper (more bodily, not less), so of the two plausible values the ASSUMED 0.45 is
the better one. That is worth knowing: the assumption flagged in the open items is
not where the error lives.

**Tangential contact stiffness is not a free parameter at all.** `K_T` changes the
answer by 1.6e-4 mm, because `ksec = min(K_T, mu*lam/|t|)` is Coulomb-capped: once
pairs slide, the couple the appliance can carry is set by `mu = 0.2`, which the
paper states. The knob remains for documentation; it is not a lever.

**Appliance fit controls the centre of rotation but is not a usable fix.** The
reference configuration makes the aligner a perfect glove -- every slave node at
exactly zero gap -- so the shell resists the crown rotating as well as translating.
Introducing thermoforming slack does move the centre of rotation to where the paper
puts it (root/crown 0.54 -> 0.03, i.e. the apex stops moving), which confirms the
mechanism. But it is not usable:

* the tipping *angle* still only reaches 2.16° against the paper's 5.92°, so most
  of the discrepancy survives;
* the contact solution degrades badly -- nearly every solve hits the cap and the
  penalty force chatters to 104 N against 7 N at baseline;
* past ~0.05 mm the appliance is no longer held at all. The aligner is supported
  only by contact with the crowns, so once enough pairs go inactive it retains a
  rigid-body freedom and the reduced system is singular.

So the fit assumption explains the *character* of the discrepancy without fixing
it, and the remaining candidate is structural: `spec.py` quotes the method paper as
simulating the PDL "with a 0.30-mm shell element", whereas this package uses three
layers of solid tetrahedra. If that is literal, the two models differ in exactly
the ratio -- rotational against translational PDL stiffness -- that this
disagreement is about. Testing it means re-meshing, not re-running.

## 14. The palatal TAD group: a comparison that does NOT reproduce

Run to close §12's caveat that the paper's fullest numbers are for the palatal TAD
group while only buccal TAD and control had been simulated. 70 steps, 47 min.

| | paper | this model |
|---|---|---|
| palatal TAD, FDI 17 at step 70 | 1.31 mm (67.53%) | **-0.089 mm** |
| palatal TAD, FDI 16 at step 70 | 1.00 mm (51.55%) | 0.213 mm |
| control, FDI 17 at step 70 | 0.73 mm (37.63%) | -0.717 mm |
| control, FDI 16 at step 70 | 0.62 mm (31.96%) | -0.370 mm |

Incisor retraction at step 70, where the paper's ordering is explicit:

| group | paper | this model |
|---|---|---|
| palatal TAD | 1.88 mm | 0.622 mm |
| buccal TAD | 1.71 mm | 0.829 mm |
| control | 0.87 mm | 0.169 mm |

**The paper's ordering is palatal > buccal > control; this model gives buccal >
palatal > control.** This retracts the claim in §12 that every group-to-group
comparison reproduces. That claim was made on two groups, where only the
TAD-versus-control contrast was in evidence; with the third group the
buccal-versus-palatal contrast is inverted. What survives is narrower: **either TAD
beats control, in the expected direction and by a plausible margin, but the model
does not resolve the difference between the two TAD positions.**

## 15. The PDL is not the problem: where the tooth actually rotates

§13 left the PDL element type (solid tetrahedra here, "0.30-mm shell element" in
the method paper) as the leading suspect for the movement-type discrepancy. It is
not the cause, and this can be settled without re-meshing.

Remove the aligner and the contact entirely and ask the PDL the classical question
directly. The teeth are rigid bodies with 6 dofs, so a tooth's effective stiffness
is the exact Schur complement of the PDL stiffness onto those dofs -- six linear
solves. For FDI 17 under a 1 N distal force at the crown point:

| | |
|---|---|
| crown displacement | +280 um |
| apex displacement | **-104 um** (opposite direction) |
| rotation | 1.07° |
| centre of resistance | 14.4 mm below the crown point, **70.9%** of the crown-apex distance |
| centre of rotation | 15.0 mm below the crown point, **73.7%** |

Both lie inside the apical portion of the root, which is exactly what the paper
describes ("the centre of the rotation ... around the apical third"). **The PDL
model reproduces the classical biomechanics on its own.** Re-meshing it as shell
elements would have been wasted effort.

So the discrepancy is entirely in how the *appliance* loads the crown: driven by
the aligner, the same tooth translates (root/crown 0.54) where a pure crown force
tips it about the apical third. The appliance is delivering a large counter-moment.
That is consistent with §13's clearance result, which reduced the engagement and
moved the centre of rotation back towards the apex. Anyone continuing this should
look at the contact model -- the passive-fit reference that puts every slave node
at zero gap over the whole crown -- and not at the PDL.

(Note for whoever checks: the single-tooth figure quoted in the README, 13-16 N/mm
pre-knee, is not the same measurement as the 26.1 N/mm here, which is pure
translation at the centroid; at the crown point, where rotation contributes, this
model gives 3.6 N/mm.)

## 16. Performance, second round: what cProfile cannot see

**cProfile lies about GPU code.** CUDA launches are asynchronous, so device work is
billed to whichever host call synchronises next. It attributed 160 ms per iteration
to the two SpGEMMs of the triple product; timed with explicit synchronisation they
take **12.9 ms**. Two optimisations were aimed at that phantom before the
measurement was corrected -- including a complete fixed-linear-map replacement for
the triple product (A = P^T Kg P has fixed patterns within a step, so A.data is a
fixed linear combination of Kg.data). It was built and verified correct to 2.2e-16,
and it is **10x slower** than cuSPARSE's SpGEMM: 49.6M triplets of gather-multiply
cannot beat a tuned sparse kernel. Recorded so nobody tries it again.

The real per-iteration split, every timer synchronised:

| | before | after |
|---|---|---|
| PARDISO solve (host) | 374 ms | 355 ms |
| GPU assembly -> host | **207 ms** | **63 ms** |
| residual check (host, 2 SpMV) | 20 ms | 17 ms |
| PDL strain + modulus (GPU) | 6 ms | 6 ms |
| contact kinematics, modulus (host) | 5 ms | 5 ms |
| **total** | **591 ms** | **451 ms** |

**The assembly cost was the reproducibility guarantee, not the arithmetic.**
`torch.segment_reduce` costs O(number of segments), and this matrix has 8.84M of
them averaging 2.7 contributions each. The contact scatter spent 48 ms writing
1.5M values because it scanned all 8.84M segments, 45 ms of it over empty ones.
Two changes:

* the deterministic path now reduces over only the segments a scatter actually
  touches (`unique_consecutive`), then scatters that compacted result;
* the default is now atomics (`OFJ_GPU_DETERMINISTIC=1` restores bit-reproducible
  reruns). Determinism was costing 143 ms per iteration -- 30% of the runtime -- to
  remove a ~1e-16 variation in a model whose nonlinear tolerance is 1e-3 and which
  has eight materially unconverged solves (§11). That is the wrong trade by
  default, though it is the right one to keep available.

Also moved to the device: the post-solve chain (reduced solution -> global
displacement -> element strains -> von Mises equivalent strain -> secant modulus),
worth ~20 ms.

**PARDISO is at its floor.** It is now 80% of the iteration. MKL's
CG-on-the-previous-factorisation mode (`iparm[3]` = 62 / 82 / 102) was tried across
successive systems of the real loop and gives 280-350 ms against 280-370 ms for
full refactorisation -- no gain, because the matrix changes enough between
nonlinear iterations that CG does not converge cheaply and it refactorises anyway.
Sparse Cholesky is supernodal dense GEMM, the one thing a T4 is bad at (1.8x the
host's FP64), so moving it to the device would not help either.

**Result:** 591 -> 451 ms per iteration; measured staging step ~30 s -> **~25 s**;
a 70-step group ~38 min -> **~30 min**. Verified identical: same iteration counts,
same active sets, tooth motions agreeing to 6e-11 (0.0000 nm on 105 um).

**A defect this round exposed.** With a large clearance PARDISO hit a zero pivot,
and the AMG fallback then died: the device assembler emits only the upper triangle,
which AMG+CG cannot use. The fallback now hands assembly back to the host and
redoes the iteration. It never fired in the production runs, but the README's claim
that both paths "degrade gracefully" was false until this was fixed.

## 17. The three-group experiment

All three anchorage groups, full 70-step protocol, one consistent code version:
29 min each, 87 min total (the same three runs would be ~90 h on the code as
shipped). Complete, zero errors, zero non-finite values, 70 steps and 7 teeth in
every checkpoint.

Reruns reproduce the earlier ones with **identical total iteration counts**
(3262 / 3375 / 3413) and movement differences of 0.000-0.007 nm on ~2 mm, so
none of the performance work in §16 changed the physics.

Expression at the end of each tooth's own active phase:

| tooth | phase | prescribed | control | buccal TAD | palatal TAD |
|---|---|---|---|---|---|
| FDI 17 | 1-20 | 2.00 mm | 1.386 (69%) | **1.606 (80%)** | 1.569 (78%) |
| FDI 16 | 11-30 | 2.00 mm | 1.448 (72%) | **1.723 (86%)** | 1.700 (85%) |
| FDI 15 | 21-40 | 2.00 mm | 0.997 (50%) | 1.266 (63%) | **1.300 (65%)** |
| FDI 14 | 31-50 | 2.00 mm | 0.614 (31%) | **0.993 (50%)** | 0.965 (48%) |
| FDI 13 | 41-60 | 2.00 mm | -0.161 (-8%) | **0.471 (24%)** | 0.231 (12%) |
| FDI 12 | 51-70 | 3.00 mm | 0.187 (6%) | **0.841 (28%)** | 0.439 (15%) |
| FDI 11 | 51-70 | 3.00 mm | 0.169 (6%) | **0.829 (28%)** | 0.622 (21%) |

Either TAD beats the control at every tooth, and the margin widens down the arch
(+11 points at the second molar, +19 at the first premolar, +32 at the canine).
The two TAD positions are within a few points of each other through the posterior
phases and separate only in the anterior phases, where the buccal TAD leads.

At step 70 the qualitative split is the same as in §14: both TAD groups hold net
distalisation on most teeth while the control is net **mesial** on every posterior
tooth, having lost the whole distalisation and drifted past its start.

### Against the paper

| | paper | control | buccal TAD | palatal TAD |
|---|---|---|---|---|
| FDI 17 at step 10 | 0.97 / 1.09-1.13 mm | 0.840 | 0.959 | 0.938 |
| FDI 17 tipping at step 10 | 5.92 / 6.63-6.87° | **1.10** | **1.36** | **1.28** |
| FDI 17 mesial drift, steps 20-40 | 0.89 / 0.61-0.76 mm | 0.752 | 0.549 | 0.586 |
| FDI 17 at step 70 | 0.73 (ctrl) / 1.31 (pal) | **-0.717** | 0.088 | **-0.089** |
| FDI 16 at step 70 | 0.62 (ctrl) / 1.00 (pal) | **-0.370** | 0.307 | **0.213** |
| FDI 11 retraction at step 70 | 0.87 / 1.71 / 1.88 mm | 0.169 | 0.829 | 0.622 |

Reproduces: step-10 molar translation within 13% in all three groups, the
group ratio (1.14 against the paper's 1.12-1.16), and the mid-protocol relapse
including its ordering (worst in control).

Does not reproduce: the tipping angle (5x too small -- §15 shows the PDL is not at
fault, the appliance is), the long-horizon retention (this model keeps relapsing
where the paper's arrests), and the **ordering of the two TAD positions** -- the
paper has palatal ahead of buccal on incisor retraction (1.88 vs 1.71 mm), this
model has buccal ahead of palatal (0.829 vs 0.622 mm). With all three groups run,
the honest summary is that the model resolves **TAD against no TAD** but not
**buccal against palatal**.

`out/staging_curves.png` plots all seven teeth and all three groups.

## 18. Convergence: two defects in the contact loop, fixed

§11 recorded eight solves per group ending with du of 0.1-0.3 and contact forces of
320-722 N. Replaying the staging chain with its warm starts intact reproduced them
exactly, and they turned out to be two separate problems.

**A limit cycle seeded by a pairing-inconsistent cache (fixed).** `solve_config`
reused the cached active set, `lam` and `tmag` whenever `npair` matched -- which it
always does -- but `build_pairs` runs every step, and ~0.8% of slave nodes (70-83
of 10,470) re-pair to a different master triangle each step. Their contact state
was carried onto a different surface: the same inconsistency §2 fixed for the gap
reference. The replayed step 16 cycled with period ~12, active set churning by
100-500 pairs, force spiking to 3000 N; the same configuration solved cold
converged to 22.18 N. The cache is now kept only while the pairing is the one it
was computed against (between a step's two PDL iterations):

| buccal TAD | before | after |
|---|---|---|
| step 12, first solve | 35 it, du 1.1e-1, 129.6 N | 27 it, du 9.6e-4, 19.1 N |
| step 15, first solve | 35 it, du 1.9e-1, 42.2 N | 21 it, du 7.6e-4, 21.8 N |
| step 16, first solve | 35 it, du 3.0e-1, 320.7 N | 21 it, du 7.4e-4, 22.2 N |

**Max-norm criteria pinned by marginal entities (fixed).** With the active set
frozen for 30 iterations and the force constant to 1e-8, du sat at 1.4e-3 against
a 1e-3 tolerance and dE at 1.29e-2 against 1e-2: a handful of free aligner nodes and
elements straddling the 7.5% bilinear knee flip every iteration and hold the
maximum. du and dE now use p99.9, as §5 already does for strain; the maxima are
still recorded as `du_max` / `dE_max`. Second solves fell from 35 iterations to
11-18.

The reported 320 N and 722 N peaks in the results from before this fix were
artefacts of the first defect. One solve in the replayed chain (step 13) still
cycles.

## 19. The miniscrews sit in the middle of the septum, not on its buccal and palatal sides

Checked with `tools/tad_check.py` against the paper's specification ("the buccal /
palatal interradicular space between the first and second molars 4 mm above the
alveolar crest") and against the root geometry. `tools/plot_tad.py` draws it.

What matches: both screws are exactly 4.00 mm apical of their crest point, lie
between 16 and 17 mesiodistally, and the precision cuts sit at the canine's mesial
cervical region on the correct side (buccal +4.6 mm, palatal -5.9 mm).

What does not:

| at screw height, mm from the 16/17 midline | buccal | palatal |
|---|---|---|
| screw | +2.96 | +0.55 |
| outermost root/PDL surface on that side | +5.79 | +5.92 |
| screw relative to that surface | **-2.82** (inside) | **-5.36** (inside) |
| nearest PDL of 16 | 0.75 | **0.29** |

The two screws are 3.5 mm apart where the roots span ~11.7 mm. The "palatal" screw
is 0.55 mm from the midline -- effectively a mid-septum point. Cause:
`inputs.py` takes the crest as the most occlusal bone node on each side of the
midline within +/-2.5 mm mesiodistally, which lands on the ridge of the septum,
not on the buccal or palatal alveolar margin; the screw is then offset straight
apically along the occlusal axis and so stays in the septum.

Root proximity does not enter the solution. The TAD is modelled as a point -- the
fixed end of the elastic, with no screw body meshed -- and the schematic
(`draw_model.py`) now draws it as a point as well. The 0.29 mm and 0.75 mm to 16's
PDL are geometry only; what the position does affect is the line of action.

The line of action as modelled (unit components):

| | sagittal | lateral (+buccal) | apical |
|---|---|---|---|
| buccal TAD | 0.906 | -0.272 | 0.323 |
| palatal TAD | 0.952 | -0.037 | 0.304 |

Moving each screw out to its own side (same height and mesiodistal position; the
bone is not in the model, so these are geometric estimates):

| | sagittal | lateral | apical |
|---|---|---|---|
| buccal, at root surface / +2 mm | 0.931 / 0.942 | -0.170 / -0.097 | 0.32 |
| palatal, at root surface / +2 mm | 0.924 / 0.906 | -0.221 / -0.289 | 0.31 |

As modelled, palatal has the larger sagittal component, which is the mechanism the
paper gives for palatal TADs performing best. With the screws on their proper
sides, this patient's arch reverses it. So the placement is a live candidate for
the inverted buccal/palatal ordering (§17), but correcting it would not by itself
move the model toward the paper -- the real patient-8 arch form, rather than the
paper's idealised dentition, also enters through exactly this geometry.

## 20. Miniscrew points placed on the real bone (Open-Full-Jaw source data)

With the dataset downloaded (`dataset/Patient_8.zip`, Git LFS, 236 MB; CC BY-NC-SA
4.0) the geometry was confirmed identical to what the package was built from:
crown centres match `slim.pkl` to 0.00e+00 mm. The maxilla volumetric mesh has
791,321 tets -- teeth 221,403, PDL 101,283, bone 468,635 (one bone tag, no
cortical/cancellous split).

`tools/tad_from_bone.py` replaces the septum-ridge rule of §19 with three steps on
the segmented anatomy: the 16/17 interradicular gap is found on each side from the
real PDL; the crest height from vertical rays through the bone across that side of
the septum; and the screw point by marching a ray outward at 4 mm apical until it
first leaves the bone -- the cortical surface. The TAD stays a point.

| | buccal | palatal |
|---|---|---|
| apical of that side's crest | 4.00 mm | 4.00 mm |
| mesiodistal clearance to 16 / 17 PDL | 0.93 / 0.93 mm | 1.24 / 1.24 mm |
| on the cortical surface, beyond the roots | +0.10 mm | +0.54 mm |
| nearest PDL of 16 (3D) | 1.56 mm (was 0.75) | 2.60 mm (was 0.29) |
| moved from the legacy point | 2.85 mm | 5.88 mm |

The two screws are now **12.1 mm apart** bucco-palatally, against 3.5 mm before.
`loads4d.py` uses these points (`build4d/tad_bone.pkl`); `OFJ_TAD_RULE=legacy`
restores the old rule for comparison. The rigid-bone buccal and palatal groups are
being re-run with them.

## 21. A deformable-bone variant of the half-maxilla model, verified

Mao et al. treat teeth and alveolar bone as rigid (the outer PDL layer is fixed and
bone is not meshed). With the Open-Full-Jaw source data available, the same model
can be solved with the bone deforming, so the effect of that assumption can be
measured rather than assumed. Switched by `OFJ_BONE=deformable`; the rigid model is
unchanged and remains the default.

**Construction** (`ofj4d/bone4d.py`, `model4d.py`, `solve4d.py`, `gpu4d.py`):

| | |
|---|---|
| domain | right half of the patient-8 maxilla bone, 245,653 tets, cut at the midsagittal plane of the aligner's symmetry constraint |
| materials | cortical 13.7 GPa within 1.5 mm of the outer bone surface (63% of the volume), cancellous 1.37 GPa, nu 0.30 -- the values of the package's full-arch model; the segmentation has one bone label |
| support | the dataset authors' own fixed region -- a superior slab, z >= 78.38 mm (11,971 nodes on the right, 9.2 mm from the nearest root); symmetry on the midsagittal cut (1,674 nodes) |
| PDL-bone | the same 0.30 mm three-layer PDL as the rigid model, its outer layer tied into the bone by barycentric interpolation in the containing bone tet; 96.3% of those nodes lie inside bone at rest, the rest (cervical rim, where the segmented PDL is thicker) are pinned to the nearest bone |
| remodelling | bone stress-free at the start of every solve and the ties re-formed as teeth move -- the bone accommodates each new tooth position, the counterpart of restoring the PDL |
| TAD | the elastic's equal and opposite force goes into the bone nodes at the screw point |

Only one thing separates the two variants: whether the bone deforms.

The authors' `.feb` was meshed separately from the `.msh` -- node-by-node matching of
their fixed set recovered 58 of 22,142 nodes and would have left the bone on 30
point supports. The region is a slab, reproduced by height on this mesh
(`tools/extract_fixed_region.py`).

**Verification** (one solve of staging step 1, buccal TAD, `tools/bone_check.py`):

| check | result |
|---|---|
| V0: rigid path after the patches vs the frozen code | identical to 1.0e-11 |
| V1: deformable code, bone E x1000 (rigid limit) vs rigid model | 1.2e-6 relative, <= 0.08 nm; same 26 iterations and force |
| V2: deformable bone, real moduli, vs rigid model | tooth translations 0.10-3.09% larger (FDI 12 most); rotations < 0.005 deg; peak bone displacement 0.50 um |

V1 passing is what establishes the tie: with near-rigid bone it reproduces the fixed
boundary to a millionth. V2 solved at the iteration cap (35 against 26), so part of
its 0.1-3% may be convergence path rather than bone compliance; the 70-step runs
will show whether the difference accumulates.

**Cost.** 304,499 reduced dof against 150,644; PARDISO's numeric factorisation is
~5.4 s per iteration against ~0.8 s on this 2-core machine, so a 70-step group is
roughly 7 h.

## 22. The factorisation moves to the GPU: cuDSS

After the studio moved to a 2-core host, and with the bone meshed, nearly all the
time is the sparse Cholesky factorisation. That is dense supernodal BLAS3 -- the
float64 work §10 found a T4 poor at against 4 cores, but not against 2. Both real
systems were dumped from the loop (`tools/dump_system.py`) and factorised offline:

| per nonlinear iteration | PARDISO (2 cores) | cuDSS float64 | cuDSS float32 + refinement |
|---|---|---|---|
| rigid, 150,644 dof, 25.9M-nnz factor | 0.69 s | **0.12 s** | 0.04 s + 8 refinement steps, 0.28 s |
| deformable, 304,499 dof, 112.3M-nnz factor | 5.66 s | **0.89 s** | 0.37 s + 10 refinement steps, 0.94 s |

float64 wins and is exact (residuals ~1e-15). float32 converges to 1e-10 in 8-10
refinement steps, but those residual products run on the host and cost more than
the cheaper factorisation saves; moving them to the GPU is the next step if more
speed is wanted.

**Integration.** A Studio allows one environment, so `nvmath-python` and
`nvidia-cudss-cu12` are installed with `pip --target` into a separate directory
(`OFJ_CUDSS_DIR`). That directory carries its own, newer cuBLAS; loading it next to
torch's would put two `libcublas.so.12` in one process. The `CuDSS` backend
therefore imports torch first and preloads `libcudss` by absolute path, so it binds
to torch's already-loaded cuBLAS 12.8.4.1 -- verified: one cuBLAS mapped, torch GPU
ops unaffected. `OFJ_LINSOLVE=auto` now prefers cuDSS when present and falls back to
PARDISO; a cuDSS failure mid-run also switches to PARDISO, never to AMG+CG.

**Verification** (one solve, buccal TAD, step 1, against PARDISO):

| | iterations | rigid motions | wall (contended) |
|---|---|---|---|
| rigid bone | 26 vs 26 | 1.6e-11 relative | per-iteration solve 0.1 s vs 0.8 s |
| deformable bone | 35 vs 35 | 9.5e-11 relative, bone displacement identical | 98 s vs 227 s |

The rigid control group ran with PARDISO; the rigid buccal and palatal groups and
all deformable runs use cuDSS. The two agree to 1e-10, so the experiment set stays
consistent.

## 23. The palatal elastic was attached to the lateral incisor, not the canine

The paper puts the elastic "at the buccal or lingual mesial cervical region of the
canine". `loads4d.py` placed the precision cut at the aligner's outer surface
nearest a seed 2.2 mm mesial and 6 mm palatal of the canine crown centre. On the
buccal side that lands on the canine. On the palatal side the arch curves inward,
and the nearest appliance surface to that seed belongs to the **lateral incisor**:
`tools/hook_check.py` found all 112 nodes of the palatal cut over FDI 12 (0.63 mm
from its surface, 3.26 mm from the canine), against all 92 buccal nodes correctly
over FDI 13.

Fix: the cut is now chosen only among outer aligner nodes that lie over the hook
tooth and on the requested side of it. The buccal cut is unchanged node for node;
the palatal cut moves onto the canine's lingual surface.

**Consequence.** Every palatal-TAD result before this fix (§14, §17, and the rigid
palatal runs made before it) loaded the lateral incisor instead of the canine. The
inverted buccal-versus-palatal ordering reported in §14 and §17 was therefore
measured with the palatal elastic on the wrong tooth, and has to be re-evaluated
with the corrected runs before any conclusion is drawn from it.

## 24. Why incisor retraction is half the paper's and the control molar ends mesial

With the corrected runs (§23) the ordering reproduces, palatal > buccal > control,
but the magnitudes do not: incisor retraction 0.92 / 0.76 / 0.17 mm against
1.88 / 1.71 / 0.87, and the control's FDI 17 ends 0.72 mm **mesial** against the
paper's 0.73 mm distal. Both have one cause.

**The dentition moves as one block.** `tools/lag_seat.py` measures each tooth's lag
behind its prescribed position and the rigid wear-in shift of the appliance. In the
control at step 70 every tooth lags by nearly the same amount, 1.89-2.94 mm
(spread 0.33 mm), and the wear-in fit has shifted the whole appliance 2.65 mm
anteriorly and rotated it 3.0°. The TAD groups lag 1.3-2.3 mm. The common part of
the activation is removed by the seat, only the differential part becomes force,
and the anchorage loss accumulates without ever arresting. The fit itself is not
at fault: of the placements tried, the full rigid fit leaves by far the least
initial interference (summed penetration at step 10: 91.5 mm, against 355 mm for
translation only and 499 mm for no fit), and it weights the teeth evenly (9-17%
of the fitted nodes each).

**Every tooth is loaded through its centre of resistance.** `tools/moment_probe.py`
restarts a step from a checkpoint and resolves the net contact force and moment on
each tooth. At control step 10 the line of action sits 64-77% of the crown-apex
distance below the crown point on **all seven teeth** (FDI 17: 5.45 N distal at
13.3 mm of 20.8 mm). §15 put the centre of resistance at ~71%. So the appliance
translates every tooth bodily, the moving molar and the anchorage teeth alike.
That is why the step-10 tipping is 1.1° against the paper's 5.9°. Bodily movement
costs several times the force per mm of crown movement (§15: 3.6 N/mm at the crown
point against 26 N/mm in translation), so each step expresses less of its
prescription and pushes harder on anchorage teeth that also cannot tip. That is the
long-horizon relapse (§12) and the short incisor retraction.

**It is not how the load is transmitted.** Ten-step control restarts, FDI 17 at step 10:

| variant | crown | root | root/crown | tipping | line of action |
|---|---|---|---|---|---|
| baseline | 0.840 mm | 0.450 mm | 0.54 | 1.10° | 13.3 mm |
| aligner E x0.5 | 0.859 | 0.384 | 0.45 | 1.34° | - |
| aligner 4 layers through the 0.7 mm sheet (was 2) | 0.842 | 0.444 | 0.53 | 1.13° | - |
| friction mu = 0 (paper: 0.2) | 0.842 | 0.454 | 0.54 | 1.09° | 13.3 mm |
| attachments removed from contact | 0.849 | 0.430 | 0.51 | 1.18° | 13.2 mm |

Aligner stiffness, its through-thickness resolution, friction and attachments
change the result by at most 0.02 mm and 0.24°. Removing friction shifts load from
tangential traction to normal pressure but leaves the line of action exactly where
it was. (An aligner at E x0.2 did not converge: every solve reached the iteration
cap with the penalty force chattering to 120 N, so it is not evidence either way.
A finer in-plane aligner mesh could not be built without an inverted element.)

What these share is the socket kinematics. The appliance is manufactured with each
tooth's socket **translated** to its prescribed position (`arch4d.py`), and a shell
four orders of magnitude stiffer than the PDL holds the crown at the socket's
orientation, whatever carries the load.

**Test: tip the sockets instead.** `OFJ_PRESCRIBE=tip` gives every socket the same
crown-centre displacement (identical to 0.0000 um) but as a rotation about a pivot
at 83% of the crown-apex distance, the middle of the apical third where the paper
reports the centre of rotation. Ten-step control, FDI 17 at step 10:

| | translate (baseline) | tip | paper |
|---|---|---|---|
| crown | 0.840 mm | 1.059 mm | 0.97 mm |
| root | +0.450 mm | **-0.302 mm** | apical-third centre of rotation |
| tipping | 1.10° | 3.81° | 5.92° |
| mesial drift of FDI 16 (anchorage) | 0.114 mm | 0.051 mm | - |
| wear-in shift of the appliance | 0.121 mm, 1.24° | 0.051 mm, 0.24° | - |

The movement type changes as predicted: the root now moves against the crown, the
step-10 crown movement goes from 13% under the paper to 9% over it, and the
anchorage teeth drift half as far. Two caveats. The pivot is taken from the paper's
reported centre of rotation, so this is partly calibration rather than independent
reproduction. And the first PDL iteration reached the iteration cap in 8 of 10
steps (du up to 4.9e-2, the penalty force chattering to 92 N once), while the second
converged in every step (du <= 9e-4, 0.2-0.8 N). Whether it closes the long-horizon
gap is being tested with 70-step runs of all three groups
(`build4d/variants/tip_run_<group>.pkl`), followed by the sensitivity runs of the
assumed inputs around that model (`senstip_*`).

**Over 70 steps, full tipping overshoots.** Control and palatal TAD finished (the
buccal run was stopped at step 17 for the calibration below). Crown-centre movement
at step 70 as a share of the prescription, against the paper's efficiencies:

| | translate | paper | tip |
|---|---|---|---|
| control FDI 17 / 16 / 11 retraction | -20% / -8% / 2% | 38% / 32% / 29% | 61% / 70% / 74% |
| palatal FDI 17 / 16 / 11 retraction | 7% / 19% / 26% | 68% / 52% / 63% | 93% / 101% / 102% |

The two socket kinematics bracket every published step-70 value, and both keep
palatal ahead of control. Every final solve ended below tolerance (max du 9.7e-4),
but 25-31% of solves reached the cap, against 4-6% with translation.

## 25. A blended socket kinematics, calibrated on the control group

**This section is calibration, not independent reproduction.** At the user's
request the model was tuned to the paper. `OFJ_TIP_BLEND=alpha` scales the tipping
rotation while the crown centre still lands exactly on its prescription: alpha = 0
is bodily translation (verified identical to the translate model, max |dT| 0), and
alpha = 1 is tipping about the apical third. The blend is fitted on the **control
group only** (`tools/calibrate_blend.py`, PCHIP interpolation over alpha = 0, 0.4,
0.6, 1), so the two TAD groups remain an out-of-sample test.

| control | 17 @10 | 17 drift 20-40 | 17 @70 | 16 @70 | 11 retraction @70 | tipping @10 | rms rel. error |
|---|---|---|---|---|---|---|---|
| paper | 0.97 | 0.89 | 0.73 | 0.62 | 0.87 | 5.92° | - |
| alpha 0 (translate) | 0.840 | 0.752 | -0.717 | -0.370 | 0.169 | 1.10° | 1.197 |
| alpha 0.4 | 0.927 | 0.540 | 0.103 | 0.653 | 1.237 | 2.15° | **0.464** |
| alpha 0.6 | 0.970 | 0.422 | 0.526 | 1.170 | 1.786 | 2.70° | 0.671 |
| alpha 1 (tip) | 1.059 | 0.226 | 1.361 | 2.183 | 2.837 | 3.81° | 1.598 |

The rms relative error is over the five displacements; tipping is kept out because no
alpha reaches 5.92°. The optimum is alpha = 0.42 (error 0.461); alpha = 0.40 is used
because its control run already exists and scores 0.464.

**No single alpha fits everything.** Each metric on its own would want a different
blend: incisor retraction 0.26, FDI 16 at step 70 0.39, FDI 17 at step 10 0.60, FDI
17 at step 70 0.70. The mid-protocol relapse of FDI 17 never reaches the paper's
0.89 mm (range 0.23-0.75 mm), nor does the tipping (range 1.1-3.8°). At alpha 0.4 the
step-10 molar and FDI 16 at step 70 fit within 13%, the incisor retraction
overshoots by 42%, and FDI 17 at step 70 stays 86% short. The blend removes most of
the systematic error, but some of the paper's behaviour is not a matter of
socket kinematics.

The palatal and buccal TAD groups at alpha 0.40 (`build4d/variants/blend_run_<group>.pkl`)
test whether a blend fitted without them reproduces their published values and the
group ordering. The control at alpha 0.40 is `cal_a0.40_control.pkl`, linked as
`blend_run_control.pkl`.

**Out-of-sample test: the TAD groups at alpha 0.40.** Neither group was used in the
fit. RMS relative error over each group's published values (crown point, as
`tools/analyse.py`):

| group | published values used | alpha 0 | **alpha 0.4** | alpha 1 |
|---|---|---|---|---|
| control (fitted) | 17 @10, 17 / 16 / 11 @70 | 1.336 | **0.480** | 1.748 |
| buccal TAD | 17 @10, 11 @70 | 0.406 | **0.081** | not run to 70 |
| palatal TAD | 17 @10, 17 / 16 / 11 @70 | 0.606 | **0.278** | 1.209 |

| alpha 0.4 | paper | model |
|---|---|---|
| buccal TAD: FDI 17 @10 / 11 retraction @70 | 1.09-1.13 / 1.71 mm | 1.044 / 1.879 mm |
| palatal TAD: FDI 17 @10 | 1.09-1.13 mm | 1.061 mm |
| palatal TAD: FDI 17 / 16 @70 | 1.31 / 1.00 mm | 1.015 / 1.496 mm |
| palatal TAD: FDI 11 retraction @70 | 1.88 mm | 2.063 mm |
| incisor retraction ordering @70 | palatal 1.88 > buccal 1.71 > control 0.87 | palatal 2.06 > buccal 1.88 > control 1.24 |

At the crown centre, where the prescription is defined, incisor expression matches
the paper within 3-4 points in all three groups (control 31.7% vs 29.0%, buccal
54.1% vs 57.0%, palatal 58.4% vs 62.7%); the crown-point (incisal edge) values run
higher because the tooth tips. The group ordering reproduces. What does not: FDI 16
in the palatal group overshoots by 50% while FDI 17 falls 23% short, the step-10
tipping stays at 2.2-2.4° (paper 5.9-6.8°), and the mid-protocol molar relapse is
0.29-0.54 mm (paper 0.61-0.89).

Convergence: the final PDL iteration of every step is within tolerance in all three
runs (max final du 9.7e-4 to 1.0e-3) and the second iteration never reached the
cap. The first iteration reached it in 8 / 13 / 18 of 70 steps (control / buccal /
palatal), with first-iteration du up to 9.2e-2.

**Literature on the movement type** (the quantity alpha controls). No Mao et al.
paper has such a parameter: the temperature control method expands the aligner in a
1 mm band either side of each contact point and the movement type emerges.
Their models tip even when bodily movement is designed: Mao et al. 2023 (Prog
Orthod, molar distalization) report FDI 17 crown +1.45 mm, root -0.60 mm and 4.63°
at step 20 with the centre of rotation "around the apical third"; the 2024 premolar
study designs "distal body movement" and still reports 4.2-5.1° of tipping. The
2023 tipping, not used in the fit, falls between alpha 0.4 (4.11° at step 20) and
0.6 (5.14°). Clinically the picture is mixed: a CBCT study of >= 2 mm Invisalign
distalization found no significant change in molar angulation (Elfouly et al. 2024),
another reports distal tipping with large scatter (Qiang et al. 2025), and a 2026
systematic review finds translation-dominant movement in 56-60% of molars (Serafin
et al., Angle Orthod). In Burstone's moment-to-force terms the planned FDI 17
kinematics are 100% (alpha 0), 94% (0.4), 89% (0.6) and 71% (1, about controlled
tipping) of the translation value. There is no standard alpha; 0.4 is defensible for
reproducing the Mao et al. models, 0 for the clinically planned movement.

## 26. Sensitivity of the final model to the assumed inputs

Control group, alpha 0.40, 70 steps each (`tools/sensitivity.py`, runs
`build4d/variants/sensblend_*`). Baseline: PDL Poisson's ratio 0.45, K_N 2000,
attachments as specified.

| control | 17 @10 | tipping @10 | 17 drift 20-40 | 17 @70 | 16 @70 | 11 retraction @70 | capped solves |
|---|---|---|---|---|---|---|---|
| paper | 0.97 | 5.92° | 0.89 | 0.73 | 0.62 | 0.87 | - |
| baseline, nu 0.45 | 0.927 | 2.15° | 0.540 | 0.103 | 0.653 | 1.237 | 8 / 140 |
| PDL nu 0.30 | 0.875 (-0.052) | 1.85° (-0.30) | 0.541 (+0.001) | 0.146 (+0.043) | 0.677 (+0.025) | 1.325 (+0.088) | 38 / 140 |
| PDL nu 0.49 | 0.977 (+0.051) | 2.71° (+0.56) | 0.501 (-0.039) | 0.277 (+0.175) | 0.837 (+0.184) | 1.583 (+0.345) | 3 / 140 |

The final PDL iteration converged in every step of both runs. **Lowering the ratio to
0.30 changes nothing by more than 0.1 mm**, so the two values the Mao et al. group
itself used (0.30 in the 2023 molar study, 0.45 in the 2024 premolar study) give
the same answer. **Raising it to 0.49 matters**: step-70 incisor retraction +0.35 mm
(+28%), both molars +0.18 mm, tipping +0.56°. Linear tetrahedra approach volumetric
locking as the ratio nears 0.5, so part of this response may belong to the element
rather than the material; it was not separated here. It moves FDI 17 at step 70 and
the tipping towards the paper but the incisor further away. The assumed 0.45 sits
between the two.

**Contact penalty stiffness.** Same control group and baseline:

| control | 17 @10 | tipping @10 | 17 drift 20-40 | 17 @70 | 16 @70 | 11 retraction @70 | capped / final unconverged |
|---|---|---|---|---|---|---|---|
| baseline, K_N 2000 | 0.927 | 2.15° | 0.540 | 0.103 | 0.653 | 1.237 | 8 / 0 |
| K_N 1000 | 0.927 (+0.000) | 2.16° (+0.01) | 0.540 (-0.001) | 0.104 (+0.001) | 0.653 (+0.000) | 1.242 (+0.004) | **0 / 0** |
| K_N 4000 | 0.932 (+0.005) | 2.16° (+0.01) | 0.523 (-0.017) | 0.142 (+0.039) | 0.692 (+0.039) | 1.289 (+0.052) | **89 / 29** |

**The answer is converged in the penalty.** Halving it changes nothing by more than
0.006 mm, and the solver then never reaches its iteration cap. Doubling it changes
results by at most 0.05 mm but degrades the contact loop badly: 89 of 140 solves hit
the cap, 29 of them the final solve of the step, so its small differences are
partly unconverged iteration rather than physics. K_N = 2000 is adequate; 1000 would
be the numerically cleaner choice with the same result.

**Attachment size.** Both footprint and protrusion scaled (`tools/att_variant.py`):

| control | 17 @10 | tipping @10 | 17 drift 20-40 | 17 @70 | 16 @70 | 11 retraction @70 | capped / final unconverged |
|---|---|---|---|---|---|---|---|
| baseline, 1x | 0.927 | 2.15° | 0.540 | 0.103 | 0.653 | 1.237 | 8 / 0 |
| 0.75x | 0.930 (+0.003) | 2.17° (+0.02) | 0.535 (-0.005) | 0.113 (+0.010) | 0.657 (+0.004) | 1.260 (+0.023) | 29 / 3 |
| 1.25x | 0.926 (-0.001) | 2.15° (-0.00) | 0.525 (-0.016) | 0.138 (+0.035) | 0.684 (+0.031) | 1.279 (+0.042) | 18 / 0 |

Attachment size changes nothing by more than 0.042 mm or 0.02°. The ASSUMED
dimensions do not drive the result.

**Summary of the sensitivity analysis.** Of the assumed inputs, only a PDL Poisson's
ratio near incompressibility matters (nu 0.49: up to +0.35 mm); nu 0.30, the contact
penalty (1000-4000) and the attachment size (0.75-1.25x) each change every reported
quantity by at most 0.1 mm, most by less than 0.05 mm. PDL nu 0.49 is not tested on
the TAD groups: 0.45 is the model value.

## 27. Measuring in the paper's coordinate system, and a split socket kinematics

**Part of the trend mismatch was the measurement.** `run4d.py` reports movement
along each tooth's own axes. The paper uses one global frame: origin at the incisor
midpoint, X parallel to the first molars' mesiobuccal cusps, Y in the occlusal plane,
Z its normal. For a proclined incisor a horizontal retraction then reads partly as
"intrusion", and the molars' distal axes differ from the sagittal Y.
`tools/paper_frame.py` measures in the paper's frame, with explicit angle
definitions (distal crown tip +, labial crown inclination +, extrusion +). The final
model (alpha 0.40), in that frame:

| trend | paper | own tooth axes | paper frame | verdict |
|---|---|---|---|---|
| palatal-TAD incisor, vertical @70 | +0.30 (extrusion) | -0.53 | **+0.22** | reproduces: was a measurement artefact |
| FDI 17 vs 16 retained @70 | 17 > 16 | 16 > 17 in all groups | control 0.25 > 0.22; buccal 1.133 vs 1.129; palatal 1.00 vs 1.05 | mostly an artefact; palatal gap 0.48 -> 0.05 mm |
| palatal-TAD incisor inclination @70 | -6.0° | - | **-6.00°** | reproduces |
| FDI 17 tipping @10 | 5.92-6.75° | 2.15-2.42° | 2.53-2.86° | real: ~40% of the paper |
| control / buccal incisor inclination @70 | +0.65° / +1.9° (labial) | - | -4.0° / -6.3° (lingual) | real |
| control incisor vertical @70 | -1.08 (intrusion) | -0.12 | +0.31 (extrusion) | real, hidden by the tooth-axis measurement |
| canine steps 10-20 | mesial crown tipping | - | distal crown tip 0.35-0.70° (root drifts mesially 0.33-0.42 mm) | real; only alpha 1 gives mesial tipping (control -0.28°) |

**Recalibrated in the paper frame**, over the same four control runs, the optimum
moves from 0.42 to **0.48** (rms relative error 0.436; 0.489 at 0.40). More useful is
where each quantity would match on its own: FDI 17 at step 10 at 0.24, incisor
retraction at 0.27, FDI 16 at step 70 at 0.55, FDI 17 at step 70 at 0.61. No uniform
alpha gives the incisor intrusion (the incisor extrudes 0.29-0.33 mm at every alpha)
or a labial incisor inclination (only alpha 0: +0.32°). The molars want more tipping
freedom than the incisors, and the incisors' lingual inclination is imposed by
applying the blend to them: at alpha 0.40 their sockets are planned to tip 5.2°.

**Split socket kinematics.** `OFJ_TIP_BLEND_INCISOR` gives the incisors their own
blend. The test model plans canine..molars at alpha 0.55 and a bodily palatal
movement of the incisors (alpha 0); verified: incisor transforms identical to
translation, crown centres exactly on the prescription, posterior sockets tipped
4.0-5.6° at step 70. Two things change at once (posterior 0.40 -> 0.55, incisors
0.40 -> 0), so the run tests a candidate model, not either factor alone. Control and
buccal TAD were run (`build4d/variants/split_run_<group>.pkl`).

**Control group, paper frame** (`tools/paper_frame.py`):

| control | paper | bodily (alpha 0) | blended alpha 0.40 | **split 0.55 / 0** |
|---|---|---|---|---|
| FDI 17 @10 | 0.97 mm | 0.898 | 1.020 | 1.066 |
| FDI 17 tipping @10 | 5.92° | 1.36° | 2.53° | 2.98° |
| FDI 17 drift steps 20-40 | 0.89 mm | 0.777 | 0.553 | 0.453 |
| FDI 17 / 16 @70 | 0.73 / 0.62 mm | -0.626 / -0.803 | 0.247 / 0.216 | **0.494 / 0.430** (17 > 16) |
| incisor retraction @70 | 0.87 mm | 0.127 | 1.250 | **0.769** |
| incisor vertical @70 | -1.08 (intrusion) | +0.300 | +0.307 | **-0.371 (intrusion)** |
| incisor inclination @70 | +0.65° (labial) | +0.32° | -4.04° | **+0.65°** |
| canine tip @20 | mesial | +0.81° | +0.35° | +0.19° (still distal) |
| rms relative error, 5 displacements | - | 1.376 | 0.489 | **0.304** |

With bodily incisors the control incisor now **flares labially by +0.65°, the paper's
value, and intrudes** instead of extruding; the retraction drops to within 12% of the
paper, and both molars move towards the paper with the second molar ahead of the first,
as reported. The overall error against the paper falls from 0.489 to 0.304. Still not
reproduced: the molar tipping (2.98° against 5.92°), the size of the intrusion
(0.37 against 1.08 mm), the mid-protocol relapse (0.45 against 0.89 mm, now further
off), and the canine's mesial tipping. Final PDL iteration within tolerance in every
step (max du 1.0e-3); 15 first iterations reached the cap, no second iteration did.

**Buccal TAD group, paper frame:**

| buccal TAD | paper | bodily (alpha 0) | blended alpha 0.40 | **split 0.55 / 0** |
|---|---|---|---|---|
| FDI 17 @10 | 1.09-1.13 mm | 1.009 | 1.137 | 1.185 |
| FDI 17 tipping @10 | 6.63-6.87° | 1.70° | 2.86° | 3.32° |
| FDI 17 drift steps 20-40 | 0.61-0.76 mm | 0.586 | 0.330 | 0.234 |
| FDI 17 / 16 @70 | - | 0.166 / 0.012 | 1.133 / 1.129 | 1.381 / 1.341 (17 > 16) |
| incisor retraction @70 | 1.71 mm | 0.848 | 2.048 | **1.553** |
| incisor vertical @70 | -1.80 (intrusion) | +0.009 | -0.065 | **-0.713 (intrusion)** |
| incisor inclination @70 | +1.9° (labial) | -1.79° | -6.30° | -1.64° (lingual, closer) |
| canine tip @20 | mesial | +0.94° | +0.50° | +0.33° (still distal) |

Final PDL iteration within tolerance in every step (max du 9.8e-4); 21 first
iterations and one second iteration reached the cap.

**Scorecard of the five mismatched trends, split model (control and buccal run; the
palatal group was not run with it):**

| trend | status |
|---|---|
| FDI 17 vs 16 retained @70 | **fixed**: 17 > 16 in both groups (paper frame + split) |
| palatal incisor vertical (extrusion) | **fixed by the measurement** (blended model, paper frame: +0.22 mm against +0.30) |
| incisor inclination | **control fixed** (+0.65°, the paper's value); buccal improved from -6.3° to -1.6° but still lingual |
| incisor intrusion (control, buccal) | direction now right in both; size 0.37 / 0.71 mm against 1.08 / 1.80 |
| molar tipping @10 | improved (2.5 -> 3.0°, 2.9 -> 3.3°) but still about half the paper |
| canine mesial tipping steps 10-20 | **not fixed** (distal 0.19-0.33°) |

The price is the mid-protocol relapse of FDI 17, which moves further from the paper
(control 0.45, buccal 0.23 mm against 0.89 and 0.61-0.76). Incisor retraction is within
12% of the paper in both groups (0.77 / 1.55 against 0.87 / 1.71). Molar tipping and
the canine are the two behaviours no socket-kinematics setting reproduced; both point to
the aligner-crown coupling (the conforming, rotation-locking pocket of section 24)
rather than to the prescription. `pipeline/run_split.sh` runs all three groups,
including the palatal group, on a new machine.

## 28. FEBio built from source and run on the Open-Full-Jaw patient 8 model

Preparation for an external cross-validation against an established solver.
`pipeline/build_febio.sh` built FEBio 4.13.0 (commit dbbcf7b) with Intel MKL 2024.2.2
PARDISO in 12 min on 4 cores. oneMKL 2025+ no longer declares `mkl_dcsrtrsv`, which
FEBio's ILU0 preconditioner calls, so the MKL version is pinned below 2025.

Smoke test: the first of the ten load steps of the dataset's own
`P8_maxilla_tipping.feb`. Full maxilla, 179,715 nodes, 791,127 linear tetrahedra
(221,121 teeth / 98,769 PDL / 471,237 bone), 472,719 equations. Neo-Hookean
materials: teeth E 2000, PDL E 0.0689 (nu 0.45), bone E 1500 (MPa). Superior bone
region fixed, and a pressure on one crown triangle per tooth, scaled so that every one
of the 14 teeth receives 1 N at full load. The step **converged** (energy norm 1.2e-13
against 4.6e-4 required, after one extra stiffness reformation when the 10
quasi-Newton updates ran out) in 68 s, 36 s of it in PARDISO (two factorisations).
FEBio has no GPU solver; the linear solver is 53% of the run, which bounds a GPU
solver plug-in to under 2x on this model.

Result at t = 0.1 (0.1 N per tooth): peak displacement teeth 212 um, PDL 45 um,
bone 0.86 um. Load-point displacement per tooth (UNN / FDI): 2/17 14.4, 3/16 7.7,
4/15 37.4, 5/14 38.1, 6/13 11.8, 7/12 28.5, 8/11 61.2, 9/21 159.1, 10/22 61.7,
11/23 13.1, 12/24 19.6, 13/25 27.7, 14/26 11.4, 15/27 16.1 um; right-side teeth move
+x, left-side -x, i.e. the loads point towards the arch centre.

As a plausibility check only (different load point, a linear PDL of E 0.0689 against
this model's bilinear 0.05 / 0.22, deformable teeth and bone): FDI 17 moves about
0.14 mm per N at FEBio's load point, against about 0.28 mm per N at the cusp tip in
this model's rigid-bone single-tooth analysis (section 15). Same order of magnitude; a
real cross-validation needs the same geometry, load and aligner in both codes.

**GPU speed-up for FEBio: a cuDSS linear-solver plug-in** (`pipeline/febio_cudss/`).
FEBio has no GPU solver, but its `LinearSolver` interface (PreProcess / Factor /
BackSolve) maps one-to-one onto cuDSS analysis / factorisation / solve, and FEBio
loads solver plug-ins. The plug-in takes FEBio's 1-based upper-triangle CSR matrix
directly (cuDSS SYMMETRIC, UPPER view; GENERAL for unsymmetric systems), keeps the
analysis while the sparsity is unchanged, and is selected through a FEBio
configuration file (`<default_linear_solver type="cudss"/>`, `<import>`). Same step,
same model, same machine (Tesla T4):

| | PARDISO (MKL, 4 threads) | cuDSS plug-in (T4) |
|---|---|---|
| equilibrium iterations / reformations / RHS evaluations | 12 / 2 / 15 | 12 / 2 / 15 |
| time in linear solver | 35.9 s | **7.05 s (5.1x)** |
| whole run (wall clock) | ~68 s | ~41 s (1.65x) |
| max node displacement difference | - | **1.0e-12 mm** (4.7e-12 of the peak) |

The solutions are identical to round-off and the nonlinear iteration is unchanged, so
the plug-in is a drop-in replacement. As predicted from the 53% linear-solver share,
the whole run gains 1.65x; the rest (stiffness, residual, model update, about 20 s) is
FEBio's CPU assembly. The plug-in does not yet update FEBio's linear-solver call
counters (the log reports 0 calls). The cross-validation of this project's solver
against FEBio (single-tooth identical-load test, then one aligner step) is not done:
the plan is in `pipeline/febio_cudss/README.md`.

**Reduced FEBio model** (`tools/febio_reduce_p8.py`): the right half arch UNN 2-8
(FDI 17-11), every tooth a rigid body, bone removed and the PDL nodes it shared fixed
(nodes that also belong to a tooth are left free, or the tooth would be pinned), same
pressures, load curves and controls, displacement-only output.

| | full model | reduced model |
|---|---|---|
| nodes / elements | 179,715 / 791,127 (deformable teeth and bone) | 38,534 / 52,236 PDL + 111,535 rigid-tooth tets |
| equations | 472,719 | **423** |
| 1 step, PARDISO | 68 s | **3.1 s** |
| 1 step, cuDSS plug-in | 41 s | 3.5 s |
| 10 steps (full load), PARDISO | ~10 min (estimated) | **15.8 s** |
| 10 steps, cuDSS plug-in | - | 15.9 s |

Only 423 equations remain because the dataset's PDL is a single element thick: nearly
every PDL node lies either on a tooth (rigid) or on the removed bone (fixed), so the
unknowns are the seven rigid bodies and a few free PDL nodes. The 10-step time is now
stiffness evaluation 8.0 s, residual 2.2 s, model update 1.5 s, linear solver 0.05 s,
and ~4 s of input, set-up and plot output. The GPU plug-in no longer helps: 0.40 s
against 0.05 s in the solver, device overhead on a 423-equation system; use PARDISO.
PARDISO and cuDSS agree to 1.0e-12 mm on the reduced model too (peak 790 um at full load).

**What the simplifications cost**, load-triangle displacement at t = 0.1
(`tools/febio_compare_reduced.py`): the reduced model moves each tooth 90-97% as far as
the full model (FDI 17 0.97, 16 0.95, 15 0.96, 14 0.97, 13 0.90, 12 0.92, 11 0.97) in
the same direction within 0.1-2.6°. Rigid teeth plus rigid bone stiffen the response
by 3-10%, which supports the rigid-body assumption of the paper and of this model.
Note for the validation: a one-element-thick PDL is much coarser across the thickness
than this project's three layers, so the single-tooth comparison should use matched
PDL meshes.

## 29. External cross-validation against FEBio, Level 1: the solver core is verified

The Level-1 test planned in §28 / `pipeline/febio_cudss/README.md`, done on 2026-09-14
(8-core host, Tesla T4 unused; FEBio 4.13.0 built by `pipeline/build_febio.sh`, MKL
2024.2 PARDISO). Same problem, node for node, in both codes:

* the half-arch PDL of `build4d/geom.pkl` (7 teeth, 38,276 nodes, 167,751 tet4, 3 layers);
* every tooth rigid: in this code the model4d projection onto 6 dofs about the
  centroid, in FEBio a `rigid body` material with `center_of_mass` at the same
  centroid and the inner PDL layer tied to it by the `rigid` node-set BC;
* outer PDL layer fixed; PDL E = 0.05 MPa, nu = 0.45;
* a 1 N force at the crown point (most occlusal crown vertex) of every tooth at once,
  in FEBio as `rigid_force` + `rigid_moment` = (p - c) x F. Three load cases:
  distal (e_d), buccal (e_b), intrusive (-e_o).

Both systems have **57,456 equations**. Scripts: `tools/febio_level1_custom.py`
(this code's `Operators`, `build_P`, `Pardiso`), `tools/febio_level1_export.py`
(the .feb), `tools/febio_level1_compare.py` (rigid motion from FEBio's inner layer by an
exact Kabsch fit, node-by-node field difference), `tools/febio_level1_bench.py`.
Outputs in `build4d/febio_val/` (`level1_linear.json`, `bench_level1.json`).

FEBio has no small-strain solver (`isotropic elastic` is St Venant-Kirchhoff in the
nonlinear `solid` solver), so FEBio was loaded with 1 mN and its results scaled x1000.

### Accuracy

Largest per-tooth difference over the seven teeth, FEBio relative to this code:

| load case | crown displ. | apex displ. | rotation | centre of rotation | whole PDL field (rel. L2) |
|---|---|---|---|---|---|
| distal, 1 mN | 0.05% | 0.06% | 0.05% | identical to 0.1% of root length | 1.22e-3 |
| buccal, 1 mN | 0.08% | 0.06% | 0.08% | identical | 1.45e-3 |
| intrusive, 1 mN | 0.07% | 0.07% | 0.14% (of rotations of 0.001-0.27°) | identical | 6.2e-4 |
| distal, **0.01 mN** | **0.00%** | **0.00%** | **0.00%** | identical | **1.22e-5** |

The residual difference falls exactly 100x when the load falls 100x: it is FEBio's
finite-deformation kinematics, not a discretisation or implementation difference.
In the linear limit the two codes agree to five significant figures on every tooth,
and FDI 17 reproduces §15 exactly (crown +280.5 um, apex -103.7 um, 1.069°, centre of
rotation 73.7% of the crown-apex distance). **The element operators, the rigid-tooth
projection, the load vector and the PARDISO solve are correct.**

### Speed (idle machine, median of 3 fresh processes)

| | this code | FEBio 4 |
|---|---|---|
| 1 thread, wall to solution | **6.12 s** (element operators 1.34, K 0.91, PtKP 0.19, PARDISO 0.94, pypardiso MKL lookup ~2.7) | **7.02 s** (input 0.51, stiffness 1.98, PARDISO 2.80 for 3 factorisations, residual 0.53, update 0.37) |
| 8 threads, wall to solution | 6.16 s (assembly 3.06, PARDISO 0.43) | 7.73 s (stiffness 3.97, PARDISO 2.22) |
| process wall, 1 / 8 threads | 7.13 / 7.15 s | 7.15 / 7.88 s |
| each further load case | 0.22 s (1 thread) / 0.12 s (8), back-substitution only | a full re-run, ~7 s |

Per factorisation the two are the same library at the same cost (0.94 s vs 2.80/3 s).
The difference is in the work, not in the code: this code solves the linear problem
once, FEBio's Newton loop needs 3 iterations and 3 reformations to *confirm*
convergence of a problem that is linear. Of this code's 6.1 s, ~2.7 s is pypardiso
failing `find_library('mkl_rt')` and globbing the conda environment for it; with
`PYPARDISO_MKL_RT=$CONDA_PREFIX/lib/libmkl_rt.so.3` that overhead measured 0.02 s, so the
real in-solve cost is ~3.4 s against FEBio's 6.4 s (about 1.9x). At 8 threads both codes'
assembly got slower (Python/NumPy and FEBio's OpenMP overhead on 5k-element-per-tooth
loops), and only the factorisation benefits.

### At a full 1 N: what the small-strain assumption costs

Same problem at the full 1 N, where FEBio's finite-deformation answer is the reference
(`build4d/febio_val/level1_neo_1N.json`). St Venant-Kirchhoff (`isotropic elastic`)
fails with negative Jacobians at 0.05-0.25 N, as that material does under large
compression; neo-Hookean with the adaptive time stepper converged in all three cases
(20, 111 and 20 steps). Difference of this code's linear solution from FEBio, per tooth:

| load case | FDI 17 / 16 / 13 (large roots) | FDI 15 / 14 / 12 / 11 (small roots) | centre of rotation | whole field |
|---|---|---|---|---|
| distal | 1-3% | this code **11-26% larger** | within 2.6 points of root length | 32% |
| buccal | 0-11% | this code **30-36% larger** | within 2.1 points | 34% |
| intrusive | 4-6% | this code 9-11% larger | within 0.7 points | 10% |

The linear solution over-predicts movement where the PDL strain is large: the small
premolar and incisor roots move 0.7-2.2 mm across a 0.3 mm PDL at 1 N with E = 0.05 MPa,
far outside small strain. The movement *type* is robust (centres of rotation agree
within 3 points of root length); the magnitude is not. For intrusion the rotations are
0.001-0.3° and their percentages are noise.

Two things make this an upper bound for the production model, not an estimate of its
error: the production PDL is bilinear (E rises to 0.22 MPa above 7.5% strain, so strains
stay far lower than with a linear 0.05 MPa), and the staging loop re-linearises on the
updated geometry every 0.1 mm step instead of taking a 1 N load in one linear jump.
How large the effect is at production strain can only be estimated from this test. The
linear 1 N solution has p99 von Mises equivalent strain of 110% (distal), 143% (buccal) and
38% (intrusive); the 1 mN runs are therefore at ~0.1% strain, where the field difference
was 1.2e-3 -- i.e. the difference grows about one-for-one with p99 strain while small
(it saturates later: 32% at 110%). The production p99 PDL strain is 3.7-9.8% (run logs),
so the first-order ESTIMATE is a 4-10% field difference and roughly half that on the rigid
tooth motions (their ratio to the field difference was ~0.4 at 1 mN). This is an
extrapolation, not a measurement, and the bilinear stiffening lowers it further; Level 2
is what measures it.

Cost at 1 N: FEBio 254 s / 286 s / 1831 s (buccal needed 111 load steps), against ~3 s for
the linear solve here -- not a like-for-like speed figure, since the physics differ.

### What Level 1 does not cover

It validates the PDL/rigid-tooth core that every staging solve is built on. It does not
validate the aligner, the contact model (passive-fit gap reference, penalty normal,
secant Coulomb friction), the bilinear PDL secant iteration or the staging morph -- the
parts §13-§15 and §24 identify as governing movement type. That is Level 2 (one aligner
step in FEBio with sliding contact), still to do. Nor does it compare the production
loop's speed: one staging step here is ~44 nonlinear iterations of the 150,644-dof
contact system (median 39 s per step on the shipped runs), which has no FEBio
counterpart until Level 2 exists.

## 30. FEBio on the GPU as the staging solver: built, partly verified, not adopted

Attempt (2026-09-14) to replace the custom staging solver by FEBio 4 with a GPU linear
solver, for finite-deformation kinematics. Outcome: **the custom solver stays the
production solver**; the FEBio pieces are kept as optional tooling.

**Built and verified**

| piece | verification |
|---|---|
| `pipeline/febio_cudss` cuDSS linear-solver plug-in | Level-1 model (57,456 eq.): node displacements equal PARDISO to **1.5e-14 mm**; linear solver 2.32 s -> 0.95 s, whole run 7.3 s -> 6.5 s (**1.12x**; FEBio's CPU stiffness assembly, 3.5 s, is untouched). Median of 3, idle. |
| `pipeline/febio_pdl` material `pdl bilinear` | the bilinear secant law of `solve4d` on neo-Hookean finite-strain kinematics; in the linear limit it reproduces the Level-1 solution to a field difference of **3.6e-6** |
| `ofj4d/febio4d.py`, `OFJ_SOLVER=febio` in `run4d.py` | writes a staging configuration as FEBio 4 (rigid crown/attachment shells, tied PDL, sliding-elastic contact with mu 0.2, penalty from K_N, nodal elastic, symmetry frame) |

Build notes for this host (CUDA 13 toolkit, cuDSS wheel for CUDA 12): link the plug-in
against a CUDA 12 runtime (`CUDART_LIB`, `CUDART_INCLUDE` in `build.sh`), and at run time
put `~/.pylibs_cudss/nvidia/{cu12,cublas,cuda_runtime}/lib` first on the library path.
libcudss finds libcublas through its own RUNPATH; taking libcublasLt from torch instead
mixes two cuBLAS releases (`undefined symbol cublasLtGetEnvironmentMode`) and FEBio drops
the plug-in with only "Failed loading plugin". `febio4d.py` sets this and raises if a
plug-in fails to load.

**Why FEBio cannot run the aligner model**

The custom contact references every gap to the configuration the aligner was made for
(`gref = gap_now - gap_pre`, §2-§3), because the generated aligner stands off the crowns
by a meshing artefact: median 0.11 mm, p95 0.29-0.35 mm, max 0.8-1.1 mm, i.e. as large as
a staging increment. FEBio's sliding-elastic contact has no per-node reference (only one
global `offset`, not load-curve driven). Moving each aligner through-thickness column by
-gap_pre instead inverts **3,014 of 122,964** aligner tets at step 1 (the standoff jumps by
up to 1.5 mm across 0.45 mm edges, near attachments and interproximally); FEBio stops with
a negative Jacobian at the first iteration at any step size. `febio4d.solve_config` now
detects this and refuses with an explanation.

**Why it would not be the fastest anyway.** Same step-1 configuration, buccal TAD group:
the custom solver on the GPU (torch assembly + cuDSS) solved it in **12.0 s** (35
iterations). FEBio's failed first attempt had already spent 36.5 s evaluating stiffness
on the CPU for 11 linear-solver calls (contended run; the idle Level-1 figure is
3.5 s per 3 reformations on half the elements). A GPU only accelerates FEBio's linear
solver, which was 22% of that time.

**What would make it possible** (not done): (1) an FEBio contact plug-in, a copy of
sliding-elastic with a per-primary-node gap reference, which keeps the model identical;
or (2) regenerating the aligner with its inner surface on the crowns, which removes the
reference for both solvers but changes the geometry and invalidates the shipped runs and
the alpha 0.40 calibration. Neither changes the speed picture above.

## 31. FEBio running the aligner experiment: contact, mesh quality, and what failed on the way

Follows §30 (2026-09-14, same host). Goal: FEBio 4 on the GPU solving the identical staging
configuration as solve4d, for a finite-deformation cross-check of the aligner step (Level 2).
Driver `ofj4d/febio4d.py` (`OFJ_SOLVER=febio`), test `tools/febio_step_test.py GROUP STEP`
(same configuration through both backends). Everything below is step 1, buccal TAD group,
final model (`OFJ_PRESCRIBE=tip OFJ_TIP_BLEND=0.40`).

### What did not work, in order, and why

| attempt | outcome | cause |
|---|---|---|
| move each aligner column by -gap_pre (§30) | 3,014 aligner tets inverted | standoff jumps 1.5 mm over 0.45 mm edges |
| FEBio sliding-elastic + per-point offset map (`pipeline/febio_contact`, superseded) | negative Jacobians at iteration 1 | a closest-point reference is not passive in FEBio's ray projection; with `auto_passive` the start is exact, but partner facets then **jump** during iterations (penetration 0.008 -> 0.81 mm, 277 N at one point) |
| same, frozen segments (`seg_up`), augmented Lagrangian, bilateral / frictionless variants | all fail | the projection jumps remain |
| fixed-pairing contact (below) | converges in energy; 22 PDL tets invert | displacement criterion held up by active-set hunting (§9); PDL slivers |
| + energy/residual criteria, extruded meshes | 2 load steps, then PDL inversion | folded prisms and micro-triangle needles in the PDL (below) |
| PARDISO instead of cuDSS | **identical** failure (same 14 iterations, same time) | the linear solver is not the cause |
| frictionless | fails earlier (t = 0.1 vs 0.2) | friction stabilises; it is not the trigger |
| repaired meshes (below) | 2 load steps, then inversion (aligner first, then apical PDL) | mesh quality was not the trigger |
| linear PDL (E2 = E1) | identical failure | the bilinear law is not the trigger |
| near-contact stiffness band (0.02 mm, no force) | identical failure | loss of seating alone is not the trigger |
| **full Newton** (`max_ups` 0 instead of 10 BFGS updates) | **all 10 load steps converge**, no inverted element | BFGS updates on the non-smooth contact problem build a wrong stiffness; the next increment inverts elements |

A per-iteration rigid-body log showed every tooth within 14 um and 0.05 deg in all converged
iterations, consistent with solve4d's step-1 answer scaled to 20-30 % of the load: the failures were
local PDL collapse, not a runaway tooth. The largest motions in the displacement file written after
a failed attempt (FDI 15 52 um) are the diverged iterate, not a solution.

### Contact: solve4d's own law as an FEBio constraint

`pipeline/febio_pairs` (`aligner-pairs`, FENLConstraint): every aligner inner node keeps its master
triangle, barycentric weights and master normal for the whole solve, as solve4d does (small
sliding). cn = ramp(t) gref + n.(u_s - sum w_i u_mi), lam = K_N max(-cn, 0), secant Coulomb friction
(K_T, mu), symmetric tangent for cn <= 0 (so a passive appliance is seated, as solve4d's `gref <= 0`
active set). At ramp = 1 the gap is exactly solve4d's. FEBio therefore differs from solve4d only in
the continuum (finite deformation, neo-Hookean aligner, finite-strain PDL); **this is not an
independent check of the contact algorithm**. Rigid crown nodes are handled by FEBio's residual and
linear-system assembly (node lists).

Other FEBio-side changes: `pipeline/febio_pdl` tangent now uses the incremental modulus (E1 / E2)
instead of the secant modulus (stress unchanged; Level-1 still 3.6e-6); `pipeline/febio_cudss`
logs zero / perturbed pivots and the true residual of every solve and runs 2 steps of iterative
refinement (`FEBIO_CUDSS_IR_STEPS`, `FEBIO_CUDSS_PIVOT`, `FEBIO_CUDSS_RESCHECK`); convergence on
energy and residual (1e-3), no displacement criterion.

### Mesh quality: defects that also affect the custom solver

Both layered meshes contained elements that overlap their neighbours. Small-strain solve4d never
checks Jacobians, so this was invisible; FEBio exposed it.

**PDL** (`tools/smooth_pdl.py`, writes `build4d/geom_fe.pkl`). The prism extrusion along vertex
normals twists where the root is concave on the scale of 0.30 mm, and a few root-surface
micro-triangles give needle tets whose free nodes carry ~1 % of the usual support -- exactly where
FEBio's first inversions were (FDI 15 nodes 1843/2803, FDI 14 nodes 1811/1812/1859). Fix: collapse
interior root edges shorter than 0.35 x the median edge (73 collapses, link condition, boundary
untouched), then extrude along directions smoothed over the root surface (20 Laplacian iterations),
keeping each node's thickness.

| PDL, 7 teeth | shipped | `geom_fe` |
|---|---|---|
| folded prisms | 29 | **0** |
| tets with quality < 0.01 / < 0.03 | 14 / 46 | **0 / 0** |
| worst quality (altitude / longest edge) | 0.0008 | 0.059 |
| elements / nodes | 167,751 / 38,276 | 166,437 / 37,984 |
| outer (bone) surface shift | - | median 28 um, p95 80 um, max 155 um |
| PDL volume | - | -1.2 % (thickness kept along the smoothed direction) |

**Aligner** (`tools/volmesh_aligner.py`, writes `build4d/aligner_fe.pkl`). A 0.70 mm shell over
embrasures and attachment edges with curvature radius below its thickness cannot be a vertex-offset
extrusion without self-intersection: 4,163 folded prisms, 1,083 slivers, and 1,060.7 mm3 of volume
that counts overlaps (smoothing the directions leaves 100 folded prisms and moves the outer surface
by 1 mm). Rebuilt as a true volume mesh: the region 0 <= phi <= 0.70 mm between aligner4d's own
signed-distance level sets, trimmed by aligner4d's gingival rule, bounded by one closed isosurface
(watertight, verified) and meshed by TetGen with the surface preserved (`pq1.5/15Y`). Inner-surface
nodes are classified by the nearest bounding level set and ordered first, so the aligner.pkl
contract of model4d is kept (`n_layer_nodes` marks the end of the surface nodes for loads4d's
precision-cut search). TetGen is AGPL-licensed (the `tetgen` Python wrapper is MIT).

| aligner | extruded (shipped) | `aligner_fe` |
|---|---|---|
| folded prisms | 4,163 | none (volume mesh) |
| tets with quality < 0.01 | 1,083 | 357 |
| tets / nodes | 122,964 / 31,410 | 134,561 / 31,958 |
| volume | 1,060.7 mm3 (with overlaps) | 1,021.9 mm3 |
| inner (contact) surface | 1,436 mm2, 10,470 nodes | 1,420 mm2, 10,286 nodes |
| mean thickness | 0.70 nominal | 0.677 mm |

No aligner element inverted in any fixed-pairing FEBio run, on either mesh.

**Checks on `geom_fe` + `aligner_fe` with the custom solver** (`OFJ_GEOM`, `OFJ_ALIGNER`):

* unstaged aligner on unmoved teeth: contact force **0.000 N**, tooth motion **0 um / 0 deg** (the §1-3
  passive invariant holds on the new meshes);
* first two staging steps against the shipped reference (`pipeline/smoke_test.py`): largest difference
  **3.4 um** (FDI 11 root, step 2); FDI 17 still +0.099 / +0.195 mm; contact force 1.45 / 2.65 N
  against 1.47 / 2.53 N. The meshing defects above therefore did not materially change solve4d's
  early results; the 70-step runs have **not** been repeated on the new meshes.

### FEBio result on the repaired meshes

With full Newton (now the driver default) FEBio completes step 1: 10 of 10 load steps, 50 Newton
iterations (one stiffness reformation each), normal termination. The FEBio model was checked to be the
same configuration as the custom solve (node coordinates equal to 5e-9 mm). Both on `geom_fe` +
`aligner_fe`, buccal TAD group, step 1 (`build4d/febio_val/step1_final_comparison.pkl`):

| FDI 17, the moving tooth | custom (small strain per solve) | FEBio (finite deformation) |
|---|---|---|
| crown distal | 105.8 um | 104.2 um (-1.5 %) |
| crown buccal / occlusal | -39.3 / 15.9 um | -33.1 / 12.6 um |
| root distal | -14.0 um | -8.0 um |
| tipping about e_b / total rotation | 0.340 / 0.365 deg | 0.318 / 0.344 deg |
| centre of rotation (crown point to root point) | 87.8 % | 90.0 % |
| p99 PDL equivalent strain | 12.1 % | 12.4 % |

The six anchorage teeth move less than 15 um in both codes; the largest per-tooth differences are 8.4 um
(centroid translation) and 0.07 deg (rotation), with sign changes on components of a few um.

**Reading.** At the first staging step the finite-deformation answer and solve4d's small-strain answer
agree to 1.5 % on the prescribed movement and to 2 points on the centre of rotation; the differences on
the anchorage teeth are a few micrometres. This measures what the small-strain continuum costs at
production strain (p99 12 %), which §29 could only estimate (4-10 % on the field). It is **not** an
independent verification of the contact algorithm: both codes use solve4d's contact law. Only step 1
of one group has been compared; later steps (larger accumulated rotations) are not.

**Speed, idle machine** (same configuration, same Tesla T4; FEBio run alone, then the custom solve alone):

| | custom solve4d (GPU assembly + cuDSS) | FEBio 4 (full Newton, cuDSS plug-in) |
|---|---|---|
| one staging solve | **10.0 s** (median of 2 after warm-up; 11.7 s cold) | **197 s** |
| nonlinear iterations | 29 (active-set) | 50 Newton, 50 stiffness reformations |
| time per iteration | 0.35 s | 3.9 s |
| breakdown | assembly on the GPU, factorisation on the GPU | stiffness evaluation 88 s (CPU), linear solver 61 s (GPU, ~1.2 s per factorisation), reformation 26 s, residual 9 s, update 6 s |

FEBio is **~20x slower per solve**. At 2 solves per step and 70 steps, one anchorage group would take
about 7.7 h in FEBio against about 23 min of solve time here. FEBio's cost is dominated by CPU stiffness
evaluation, which a GPU linear-solver plug-in cannot reach; full Newton (needed for convergence) also
reforms the stiffness every iteration.

**Conclusion.** FEBio with the GPU linear solver now runs the aligner experiment (fixed-pairing contact,
repaired meshes, full Newton) and agrees with the custom solver at step 1 to 1.5 % on the prescribed
movement. It is not faster -- about 20x slower -- so the custom solver remains the production solver;
`OFJ_SOLVER=febio` is a finite-deformation cross-check.


## 32. Rigid vs deformable alveolar bone in the buccal TAD group

The question: does the paper's rigid-bone assumption change the buccal TAD result? Same model twice,
differing only in `OFJ_BONE` (§21): final model (`OFJ_PRESCRIBE=tip OFJ_TIP_BLEND=0.40`), shipped
meshes (`geom.pkl`, `aligner.pkl`, so the rigid reference is the shipped
`build4d/final/run_buccal_tad.pkl`), cuDSS on the GPU. Deformable run:
`build4d/final/run_buccal_tad_bone.sh` -> `run_buccal_tad_bone.pkl` / `.log`; comparison
`tools/compare_bone.py`.

**Bone model** rebuilt from the Open-Full-Jaw data (`ofj4d/bone4d.py`, fixed region from
`tools/extract_fixed_region.py`), identical to §21: right half-maxilla 245,653 tets / 63,814 nodes,
cortical 13.7 GPa within 1.5 mm of the outer surface (63.1 % of the volume), cancellous 1.37 GPa,
nu 0.30; 11,971 fixed nodes (superior slab z >= 78.38 mm, 9.2 mm from the nearest root), 1,674
symmetry nodes; no tet below quality 0.05. 304,499 reduced dof against 150,644 rigid. The elastic's
reaction goes into 4 bone nodes at the miniscrew.

**Verification, one solve of step 1** (`tools/bone_check.py`, `build4d/bone_check/`):

| | K_N 2000 (production) | K_N 1000 (§26: same answer, no capped solves) |
|---|---|---|
| iterations rigid / E x1000 / real | 35 / 35 / 32 -- all at or near the cap | 21 / 21 / 21, converged |
| rigid limit (bone E x1000) vs rigid | <= 0.29 um (FDI 12), contact force +5 % | **<= 0.04 nm, relative 8.4e-6**, force identical |
| real bone vs rigid | 0.10-0.74 um, rotations < 0.006 deg | **0.005-0.41 %, <= 41 nm**, rotations < 0.0004 deg |
| peak bone displacement | 0.33 um | 0.33 um |

At the production penalty every solve hunts to the iteration cap (§9), and where the cycle stops
differs between runs by up to ~0.3 um -- larger than the bone effect itself. With the cycle removed the
tie is exact and **deformable bone changes step-1 tooth movement by less than 0.5 %**.

### 70 steps (`tools/compare_bone.py`, `build4d/final/compare_bone_buccal_tad.json`)

| buccal TAD | rigid | deformable | difference |
|---|---|---|---|
| FDI 17 distal, step 10 / 20 / 70 | 1.044 / 1.832 / 0.935 mm | 1.043 / 1.827 / 0.938 mm | -0.1 / -0.3 / +0.4 % |
| FDI 17 mesial relapse, steps 20-40 | 0.322 mm | 0.320 mm | -2 um |
| FDI 16 distal, step 70 | 1.350 mm | 1.354 mm | +4 um |
| incisor retraction, step 70 | 1.879 mm (62.6 %) | 1.882 mm (62.7 %) | +3 um |
| FDI 17 tipping, step 10 / 70 | -2.370 / -2.910 deg | -2.367 / -2.917 deg | < 0.01 deg |
| capped solves / median iterations | 13 of 140 / 23 | 10 of 140 / 22 | |
| median contact force, max p99 PDL strain | 3.63 N, 10.3 % | 3.61 N, 10.3 % | |

Over all 70 steps and seven teeth the largest difference in any crown or root movement is **14 um**
(FDI 12 crown buccal, step 30); every rotation agrees within 0.03 deg. Peak bone displacement is a few
tenths of a micrometre per solve.

**Limitation.** The bone geometry is not remodelled: it is stress-free at the start of every solve
(§21) but keeps its original sockets, so as teeth move the outer PDL layer leaves the socket wall. The
fraction of those nodes inside a bone tet fell from 0.96 (step 1) to 0.86 (10), 0.62 (20), a minimum of
0.60 (30), and recovered to 0.68 (70); the remainder are tied to the nearest bone tet with clamped
weights. The result is therefore "a compliant bone with this tie", not a model of socket remodelling.

**Conclusion.** For the buccal TAD group the rigid-bone assumption of Mao et al. is adequate: deformable
alveolar bone changes no reported quantity by more than 14 um or 0.4 %. The miniscrew position study
(§33) therefore uses rigid bone. The palatal group was not run with deformable bone.

## 33. Miniscrew position study: set-up

Both TAD groups, final model, 70 steps per run, five positions besides the paper's (16/17, 4 mm
apical of the crest). The placement rule is §20's, generalised: `tools/tad_from_bone.py --pair A B
--apical H --out FILE` finds the interradicular gap between the two teeth from the real PDL, the crest
height from vertical rays through the bone, and marches outward at H mm apical to the cortical surface,
separately on the buccal and palatal sides. Its defaults reproduce the shipped `tad_bone.pkl` to
1.4e-14 mm. A run selects a file with `OFJ_TAD_FILE` (`ofj4d/loads4d.py`). Files:
`build4d/tad_positions/tad_<variant>.pkl`; runner `build4d/tad_positions/run_grid.sh BONE`;
summary `tools/compare_tad.py --bone BONE`. The bone model (rigid or deformable) follows the outcome
of §32.

| variant | site | apical of crest | gap at screw height B / P | clearance to each root B / P | beyond the roots B / P |
|---|---|---|---|---|---|
| baseline (paper) | 16/17 | 4 mm | 1.85 / 2.48 mm | 0.93 / 1.24 mm | +0.10 / +0.54 mm |
| `md1516_a4` | 15/16 | 4 mm | 1.33 / 1.45 mm | 0.67 / 0.73 mm | **-0.21 / -1.65 mm** |
| `md1415_a4` | 14/15 | 4 mm | 2.14 / 2.21 mm | 1.07 / 1.11 mm | +0.78 / +1.12 mm |
| `h1617_a2` | 16/17 | 2 mm | 1.95 / 2.63 mm | 0.98 / 1.32 mm | **-1.23 / -0.59 mm** |
| `h1617_a6` | 16/17 | 6 mm | 2.29 / 2.59 mm | 1.15 / 1.29 mm | +0.55 / +1.01 mm |
| `h1617_a8` | 16/17 | 8 mm | 2.89 / 3.64 mm | 1.45 / 1.82 mm | +1.32 / +0.96 mm |

**Clinical caveats.** The TAD is a point: root proximity and screw diameter do not enter the model.
At 15/16 the interradicular space at 4 mm is 1.3-1.5 mm, too narrow for a 1.4-1.6 mm miniscrew without
root contact; at 15/16 and at 2 mm the cortical exit lies inside the outermost root surface (concave or
thin bone). These two variants are mechanical what-ifs, not clinically placeable sites on this patient.

**Line of action at step 0** (elastic from the precision cut to the screw; unit components in the paper
frame: Y sagittal, posterior +; X transverse; Z occlusal-plane normal, occlusal +):

| variant | sagittal B / P | transverse B / P | vertical B / P | elastic length B / P |
|---|---|---|---|---|
| 16/17, 2 mm | 0.903 / 0.908 | 0.282 / 0.224 | -0.325 / -0.353 | 26.6 / 27.6 mm |
| **16/17, 4 mm (paper)** | **0.864 / 0.893** | 0.315 / 0.179 | **-0.393 / -0.412** | 27.2 / 28.3 mm |
| 16/17, 6 mm | 0.832 / 0.868 | 0.316 / 0.165 | -0.456 / -0.469 | 27.9 / 29.2 mm |
| 16/17, 8 mm | 0.795 / 0.836 | 0.327 / 0.169 | -0.512 / -0.522 | 28.9 / 30.1 mm |
| 15/16, 4 mm | 0.815 / 0.859 | 0.291 / 0.206 | -0.501 / -0.469 | 17.5 / 20.1 mm |
| 14/15, 4 mm | 0.751 / 0.782 | 0.180 / 0.195 | -0.635 / -0.592 | 12.9 / 15.9 mm |

A higher or more mesial screw trades sagittal (distalising) for apical (intrusive) pull, and the mesial
sites halve the elastic length; the palatal screw keeps the larger sagittal component at every site.
The staging runs show what that does to anchorage; results follow.

### Clinical placeability on this patient

The model treats the miniscrew as a **point** -- the fixed end of the elastic. Screw diameter, length,
insertion angle and root proximity do not enter the mechanics, so a position that would contact a root
is still solved. Placeability is therefore judged separately, from the anatomy
(`build4d/tad_positions/tad_*.pkl`): the PDL-to-PDL gap at screw height, the clearance left on each side
by a screw centred in it, and whether the cortical insertion point lies beyond the outermost root.
Criterion (an assumption of this study, not taken from the paper): **at least 0.5 mm from the screw
surface to the PDL on both sides, and the insertion point outside the root envelope.**

| site | side | gap | clearance per side, screw 1.2 / 1.5 / 1.8 mm | beyond roots | nearest PDL (3D) | verdict |
|---|---|---|---|---|---|---|
| 16/17, 4 mm (paper) | buccal | 1.85 mm | 0.33 / 0.18 / 0.03 mm | +0.10 mm | 1.56 mm | what-if |
| 16/17, 4 mm (paper) | palatal | 2.48 mm | 0.64 / 0.49 / 0.34 mm | +0.54 mm | 2.60 mm | placeable, d <= 1.2 mm |
| 15/16, 4 mm | buccal | 1.33 mm | 0.07 / -0.08 / -0.23 mm | -0.21 mm | 1.07 mm | what-if (inside root envelope) |
| 15/16, 4 mm | palatal | 1.45 mm | 0.13 / -0.02 / -0.17 mm | -1.65 mm | 2.32 mm | what-if (inside root envelope) |
| 14/15, 4 mm | buccal | 2.14 mm | 0.47 / 0.32 / 0.17 mm | +0.78 mm | 2.27 mm | what-if |
| 14/15, 4 mm | palatal | 2.21 mm | 0.51 / 0.36 / 0.21 mm | +1.12 mm | 2.58 mm | placeable, d <= 1.2 mm |
| 16/17, 2 mm | buccal | 1.95 mm | 0.38 / 0.23 / 0.08 mm | -1.23 mm | 1.16 mm | what-if (inside root envelope) |
| 16/17, 2 mm | palatal | 2.63 mm | 0.72 / 0.57 / 0.42 mm | -0.59 mm | 1.87 mm | what-if (inside root envelope) |
| 16/17, 6 mm | buccal | 2.29 mm | 0.55 / 0.40 / 0.25 mm | +0.55 mm | 1.91 mm | placeable, d <= 1.2 mm |
| 16/17, 6 mm | palatal | 2.59 mm | 0.69 / 0.54 / 0.39 mm | +1.01 mm | 3.00 mm | placeable, d <= 1.5 mm |
| 16/17, 8 mm | buccal | 2.89 mm | 0.85 / 0.70 / 0.55 mm | +1.32 mm | 2.72 mm | placeable, d <= 1.8 mm |
| 16/17, 8 mm | palatal | 3.64 mm | 1.22 / 1.07 / 0.92 mm | +0.96 mm | 3.08 mm | placeable, d <= 1.8 mm |

On this patient even the paper's buccal site is too narrow for a common miniscrew under this criterion;
only the 6 and 8 mm heights are placeable on both sides. The 15/16 site and the 2 mm height are
**mechanical what-ifs**: they show the direction and size of the effect of moving the line of action,
not recommended placements. Any article reporting this study should present them that way.

## 34. Miniscrew position study: results (repaired meshes, K_N 1000)

Production rerun (`build4d/v2/`): control, both paper-site baselines and the 10 position runs on `geom_fe` + `aligner_fe`
with K_N 1000 (§26: same answer, no hunting), final model alpha 0.40, rigid bone (§32), 70 steps. All 13 runs complete.
Summary tool `tools/compare_tad.py --set v2` (`build4d/v2/compare_positions_v2.json`); figures `out/positions/*_v2.png`.
The miniscrew is a point; which positions are clinically placeable on this patient is §33.

**Buccal TAD** (mm; change from the paper's site in brackets)

| miniscrew | FDI 17 @20 | FDI 17 relapse 20-40 | FDI 17 @70 | FDI 16 @70 | incisor retraction @70 | expression | capped solves |
|---|---|---|---|---|---|---|---|
| 16/17, 2 mm | 1.881 (+0.050) | 0.284 (-0.047) | 1.052 (+0.171) | 1.463 (+0.167) | 2.125 (+0.314) | 70.8 % | 0 of 140 |
| 16/17, 4 mm (paper) | 1.831 | 0.332 | 0.881 | 1.295 | 1.811 | 60.4 % | 0 of 140 |
| 16/17, 6 mm | 1.787 (-0.044) | 0.375 (+0.043) | 0.726 (-0.155) | 1.144 (-0.152) | 1.529 (-0.282) | 51.0 % | 0 of 140 |
| 16/17, 8 mm | 1.744 (-0.088) | 0.417 (+0.085) | 0.573 (-0.308) | 0.995 (-0.300) | 1.260 (-0.551) | 42.0 % | 1 of 140 |
| 15/16, 4 mm | 1.759 (-0.072) | 0.401 (+0.070) | 0.637 (-0.243) | 1.056 (-0.239) | 1.364 (-0.447) | 45.5 % | 0 of 140 |
| 14/15, 4 mm | 1.670 (-0.161) | 0.482 (+0.151) | 0.348 (-0.533) | 0.773 (-0.522) | 0.849 (-0.962) | 28.3 % | 0 of 140 |

**Palatal TAD** (mm; change from the paper's site in brackets)

| miniscrew | FDI 17 @20 | FDI 17 relapse 20-40 | FDI 17 @70 | FDI 16 @70 | incisor retraction @70 | expression | capped solves |
|---|---|---|---|---|---|---|---|
| 16/17, 2 mm | 1.884 (+0.022) | 0.287 (-0.021) | 1.025 (+0.073) | 1.527 (+0.088) | 2.193 (+0.189) | 73.1 % | 2 of 140 |
| 16/17, 4 mm (paper) | 1.861 | 0.308 | 0.952 | 1.439 | 2.003 | 66.8 % | 2 of 140 |
| 16/17, 6 mm | 1.834 (-0.027) | 0.334 (+0.026) | 0.862 (-0.090) | 1.336 (-0.103) | 1.794 (-0.209) | 59.8 % | 1 of 140 |
| 16/17, 8 mm | 1.805 (-0.056) | 0.362 (+0.054) | 0.765 (-0.187) | 1.229 (-0.210) | 1.581 (-0.422) | 52.7 % | 1 of 140 |
| 15/16, 4 mm | 1.827 (-0.034) | 0.341 (+0.033) | 0.840 (-0.112) | 1.317 (-0.122) | 1.762 (-0.241) | 58.7 % | 1 of 140 |
| 14/15, 4 mm | 1.757 (-0.105) | 0.407 (+0.098) | 0.612 (-0.339) | 1.065 (-0.374) | 1.271 (-0.733) | 42.4 % | 0 of 140 |

Control (no miniscrew): FDI 17 @70 0.052 mm, incisor retraction @70 1.183 mm.

**Effect of the rerun** against the same runs on the shipped meshes with K_N 2000 (12 runs): FDI 17 @70 changed by -0.066 to -0.028 mm, incisor retraction @70 by -0.068 to -0.029 mm; capped solves 13-19 -> 0-2 per run. Every position effect (difference from the paper's site) and every ranking is unchanged.

**PDL stress, paper-site baselines** (`tools/pdl_stress.py`, first PDL iteration of the step: the activated appliance; 5th-percentile hydrostatic stress = compression, kPa; replays reproduce the stored runs to <= 0.1 um):

| group | step | FDI 17 | FDI 16 | FDI 13 | FDI 12 | FDI 11 |
|---|---|---|---|---|---|---|
| buccal | 10 | -43.4 | -7.8 | -3.6 | -2.0 | -2.2 |
| buccal | 20 | -20.1 | -25.5 | -8.9 | -7.5 | -6.1 |
| buccal | 40 | -5.0 | -5.5 | -10.4 | -8.6 | -7.3 |
| buccal | 70 | -6.6 | -3.8 | -16.2 | -51.5 | -49.6 |
| palatal | 10 | -43.5 | -8.3 | -3.1 | -1.7 | -2.2 |
| palatal | 20 | -19.5 | -26.1 | -8.7 | -7.2 | -5.5 |
| palatal | 40 | -4.8 | -5.9 | -10.3 | -8.1 | -6.2 |
| palatal | 70 | -5.6 | -3.6 | -16.4 | -52.5 | -50.0 |

Compression concentrates in the PDL of the tooth being moved at that phase (FDI 17 at step 10-20, incisors at step 70), most in the cervical third.

**Reading.** Anchorage improves as the elastic's line of action becomes more horizontal and posterior: a lower or more distal miniscrew. Each 2 mm of height and each interradicular space mesially cost molar distalisation and incisor retraction in both groups; the palatal miniscrew is less sensitive to position and outperforms the buccal one at every position. On this patient the positions that give the best anchorage (low, distal) are the least placeable (§33): the clinical choice trades root clearance against line of action.

Robustness of these rankings is 35; mesh convergence is 36; elastic force 37, hook location 38, friction 39, deformable bone 40, overcorrection 41 and staging 42 follow.

---

## 35. Robustness of the position ranking (40 runs)

`build4d/v2/jobs_robustness.txt` + `jobs_rob6mm.txt`, summaries `build4d/v2/compare_positions_{a000,a060,nu030,nu049}.json`.
Five sites (16/17 at 2, 4, 6, 8 mm and 14/15 at 4 mm) x both groups x four changed assumptions:
alpha 0, alpha 0.60, PDL nu 0.30, PDL nu 0.49. All 40 runs complete, 70 steps, exit 0.

**The ranking 2 mm > paper site > 6 mm > 8 mm > 14/15 holds in all 10 group-by-set comparisons**, on
FDI 17 @70, FDI 16 @70, incisor @70 and the 20-40 relapse. Absolute values move by much more than the
differences between sites: FDI 17 @70 at the paper site is -0.034 mm (alpha 0), 0.881 (production),
1.328 (alpha 0.60), 1.639 (nu 0.30), 0.492 (nu 0.49). Ranking is the output to trust; absolute
movement is not.

## 36. Mesh convergence (4 runs, M1-M4)

`build4d/v2/jobs_meshconv.txt`, summary `build4d/v2/compare_meshconv.json`, figure `out/paper/figS1_meshconvergence.png`.
M0 production (151,774 dof) against M1 PDL 6 layers, M2 PDL x2 in plane, M3 aligner 0.30 mm,
M4 all three (815,891 dof), 20 steps each.

M0 vs M4 at step 20: FDI 17 distal 6.1 um (0.33 %), FDI 16 5.4 um (0.46 %), incisor 8.0 um (2.4 %),
contact force 0.15 %, p99 PDL strain +0.30 points (4 %). M1+M2+M3 sum to +5.4 um on FDI 17 against
M4's +6.1, so the refinements are nearly independent. PDL strain is the most mesh-sensitive quantity.
**The production mesh is converged for every quantity this study reports**; the smallest position
effect anywhere (19 um) is 3x the largest discretisation error.

M4 needed `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True` and `OFJ_LINSOLVE=pardiso` on a 16 GB T4
(see 43); 6 of its 20 steps hit the iteration cap on the first PDL iteration, none on the second.

## 37. Elastic force (18 runs)

`build4d/v2/jobs_force.txt`, summaries `compare_positions_f{100,200,250}.json` and `force_study.json`,
figure `out/paper/fig5_force.png`. 100/200/250 gf at 4, 6, 8 mm, both groups (150 gf is the v2 set).

**Every outcome is linear in force over 100-250 gf**: R2 >= 0.995 (worst 0.9948, buccal incisor at
8 mm), largest deviation from a straight line 3.3 um. Extrapolating to 0 gf reproduces the control
run (FDI 17 0.039-0.058 vs 0.052; relapse 0.549-0.551 vs 0.549; incisor 1.15-1.16 vs 1.183). So
anchorage = control + slope(site) x force.

Slope (mm per 100 gf), buccal / palatal: FDI 17 0.554/0.596 (4 mm), 0.456/0.539 (6), 0.357/0.476 (8);
incisor 0.433/0.573, 0.250/0.431, **0.076**/0.288. The force lever collapses for buccal incisor
retraction as the screw rises. Force needed to match the paper site at 150 gf: FDI 17 184 gf (buccal
6 mm), 236 (buccal 8), 167/189 (palatal); incisor 262 (buccal 6), **not attainable** (buccal 8),
198/295 (palatal).

## 38. Precision-cut (hook) location (12 runs)

`build4d/v2/jobs_hook.txt`, summaries `compare_positions_hk14.json`, `hk12.json`, figure `out/paper/fig6_hook.png`.
`OFJ_HOOK_FDI` moves the cut to the mesial cervical region of FDI 14 or 12, same side as the screw.

**First premolar > canine > lateral incisor on every measure, in both groups, at all three heights.**
Palatal: premolar hook +0.237/+0.264/+0.282 mm FDI 17 at 4/6/8 mm and +0.31 to +0.37 mm incisor - the
gain grows with height, so a palatal screw at 8 mm with a premolar cut (1.047 mm) beats the paper
configuration (0.952 mm) at the same force, and is placeable (33). Buccal: premolar +0.045 to +0.064 mm
only; lateral incisor -0.17 to -0.25 mm, consistent with its more transverse line of pull (dir x
-0.518 vs -0.318 at the canine). Cost: 0.12-0.15 deg more FDI 17 tipping at step 10 with the palatal
premolar cut. Convergence improves (0 capped solves vs 1-2).

## 39. Aligner-crown friction (4 runs)

`build4d/v2/jobs_friction.txt`, summaries `compare_positions_mu{010,030}.json`. mu 0.10 and 0.30
against the model's 0.20, paper site, both groups.

Tripling friction changes molars by 0.02-0.03 mm and incisors by 0.05-0.06 mm (largest change
anywhere 64 um), monotone, more friction giving slightly more anchorage. **mu is not a sensitive
assumption.** At mu 0.30, 30-35 of 140 solves per run hit the cap - but on the *first* PDL iteration;
the final iteration converged in 139 of 140 steps (one at 1.03e-3 vs the 1e-3 tolerance), so the
reported positions are converged. Higher friction means more stick-slip switching in the active set.

## 40. Deformable bone on the repaired meshes, both groups (2 runs)

`build4d/v2/jobs_bone.txt`. 32 tested only the buccal group on the shipped meshes; this repeats it on
`geom_fe`/`aligner_fe` with K_N 1000 and adds the palatal group. 305,629 reduced dof, ~60 s/step,
0 capped solves.

Largest difference in any crown or root movement over all 70 steps and 7 teeth: **15.6 um (buccal),
17.1 um (palatal)**; rotations <= 0.03 deg; FDI 17 @70 +1.2 %, incisor @70 +0.9 % in both groups.
Tie fraction 0.96 (step 1) -> 0.59-0.60 (step ~32) -> 0.68-0.69 (step 70), as in 32. **The rigid-bone
assumption is adequate for both groups**, and 32's limitation (sockets are not remodelled) still applies.

## 41. Molar overcorrection (12 runs)

`build4d/v2/jobs_overcorr.txt`, summaries `compare_positions_oc{10,20,30}.json`, figure `out/paper/fig7_overcorrection.png`.
`OFJ_OVERCORR` scales the planned per-step distal movement of FDI 17 and 16 (2.2 / 2.4 / 2.6 mm
instead of 2.00), at 6 and 8 mm, both groups.

Molar response is nearly proportional (~0.0085 mm of FDI 17 per 1 %). To restore the paper site's
anchorage: FDI 17 needs 18 % (buccal 6 mm), 36 % (buccal 8), 11 % / 22 % (palatal 6 / 8); FDI 16
needs 9 / 17 / 6 / 12 %. **Two costs**: incisor retraction *falls* slightly as overcorrection rises
(1.529 -> 1.495 mm at buccal 6 mm), and the 20-40 relapse grows (0.375 -> 0.459 mm). Overcorrection
is a molar-only remedy that buys position with more post-phase drift.

## 42. Staging protocol (6 runs)

`build4d/v2/jobs_staging.txt`, figure `out/paper/fig8_staging.png`. `OFJ_STAGING=pairs` moves 17+16,
then 15+14, then 13, incisors on V timing (70 steps); `seq` moves one tooth at a time, no overlap
(120 steps). Both run for control and both groups; totals prescribed are unchanged.

FDI 17 at the end of the protocol (control / buccal / palatal): V pattern 0.052 / 0.881 / 0.952;
pairs 0.033 / 0.864 / 0.932; **sequential -0.298 / 1.187 / 1.300**. The miniscrew is worth
0.83-0.90 mm under the V pattern and **1.49-1.60 mm under sequential staging**. Sequential is the
best protocol with a miniscrew (incisor expression 63 / 74 % vs 60 / 67 %) and the worst without one
(net mesial second molar, expression 24 % vs 39 %). Pairs is marginally worse than V overall but
drifts less after each phase (0.77-0.81 vs 0.91-0.95 mm). Sequential costs 120 steps instead of 70.

## 43. Code and infrastructure added for these studies

* **`ofj4d/spec.py`**: `OFJ_ELASTIC_GF` (elastic force, default 150), `OFJ_HOOK_FDI` (precision cut on
  another tooth), `OFJ_STAGING` (`v` | `pairs` | `seq`), `OFJ_OVERCORR` (molar overcorrection);
  `step_mm(unn)` replaces the two step constants. `ofj4d/arch4d.py` uses `step_mm`. All default to the
  production model: with no variables set the prescription is identical to the shipped logs to the
  micrometre, and the smoke test reproduces its reference to 0.8100 um as before.
* **`ofj4d/gpu4d.py` (`for_step`)**: the previous system is now fully released before the next step's
  is built. `close()` freed only the per-solve tensors, so at every step change the old projections and
  scatter maps were still resident - at 816k dof that doubled peak device memory and forced M4 onto the
  host assembly path. Memory lifetime only; the smoke test is bit-identical.
* **`build4d/v2/guard.sh`, `supervisor.sh`, `status.sh`**: a resumable chain for an interruptible
  machine. The guard is single-instance (pid file, not flock - a lock fd is inherited by the supervisor
  and every worker, so a dead guard could not be replaced), restarts the supervisor after a preemption,
  clears orphaned workers first, and exits on the supervisor's versioned end marker, which it reads from
  `supervisor.sh` so that adding phases needs no guard change. `~/.lightning_studio/on_start.sh` relaunches
  the guard after every preemption. Four preemptions during these studies; all resumed automatically,
  each costing at most the step in flight.
* **`tools/paper_figs.py`**: regenerates every manuscript figure and table from the saved summaries.
  Draft manuscript in `docs/manuscript.md`; figures and tables in `out/paper/`.

## 44. Robustness of the precision-cut finding (8 runs)

`build4d/v2/jobs_hookrob.txt`. First-premolar cut (`OFJ_HOOK_FDI=14`) at the 6 mm site under alpha 0,
alpha 0.60, PDL nu 0.30 and 0.49, both groups; the canine comparators are the 35 runs
`rob_<set>_*_h1617_a6`. All 8 complete, 70 steps, exit 0 (the queue resumed once after a preemption at
~07:55 and once after a supervisor restart at ~09:13, each from the step checkpoint).

Premolar minus canine at step 70 (mm):

| input set | buccal FDI 17 | buccal FDI 11 | palatal FDI 17 | palatal FDI 11 |
|---|---|---|---|---|
| production (alpha 0.40, nu 0.45) | +0.055 | +0.025 | +0.264 | +0.344 |
| alpha 0 | +0.050 | +0.028 | +0.239 | +0.325 |
| alpha 0.60 | +0.056 | +0.024 | +0.272 | +0.350 |
| PDL nu 0.30 | +0.083 | +0.073 | +0.359 | +0.488 |
| PDL nu 0.49 | +0.012 | -0.014 | +0.116 | +0.165 |

**Palatal: the premolar cut is better on both measures in every set** (+0.12 to +0.36 mm FDI 17, +0.17
to +0.49 mm incisor) - the headline claim is robust. **Buccal: always small, and a tie at nu 0.49**
(+12 um molar, -14 um incisor); state it as "no worse than the canine". As in 35, a near-incompressible
PDL compresses every difference. Values from `tools/compare_tad.py --set hka000` (hka060, hknu030, hknu049;
summaries `build4d/v2/compare_positions_hk*.json`) against the canine runs of 35; 0-4 capped solves per run.

## 45. PDL stress at 250 gf (replays, no new runs)

`tools/pdl_stress.py` on `build4d/v2/runs/rob_f250_{buccal,palatal}_tad_baseline.pkl` at steps
10/20/40/70 with `OFJ_ELASTIC_GF=250` and the v2 environment, against the 150 gf baselines of 34.
Outputs `build4d/v2/stress/f250_*_tad_baseline.json` and `_step<s>.npz` (per-element fields for contour
maps), summary `build4d/v2/stress/force_stress_summary.json`. The replays reproduce the stored runs to
<= 2e-5 mm.

Peak compression (hydrostatic 5th percentile, any tooth, any step): **buccal -51.5 -> -53.7 kPa (x1.04),
palatal -52.5 -> -55.5 kPa (x1.06)**; cervical third -62.5 -> -65.2 and -63.5 -> -67.1 kPa. For the tooth
being moved, every step rises by 4-10 %. **67 % more elastic force adds 4-6 % to peak ligament
compression**, because the ligament is loaded mainly by the aligner's activation (contact 3-16 N against
an elastic of 1.5-2.5 N). Absolute compression exceeds capillary pressure (~4.7 kPa) at both forces, as is
normal for the loaded state of an activated appliance in FE; only the relative result is defensible, and
no root-resorption threshold is implied.

## 46. Force design chart

`tools/paper_nomogram.py` -> `out/paper/fig9_nomogram.png` / `.pdf` and `out/paper/table5c_nomogram.md`.
The linear force law of 37 (slope and intercept per height from `force_study.json`), interpolated between
4, 6 and 8 mm, solved for the force that gives a target FDI 17 @70; drawn dotted and tabled as ">250 gf"
beyond the simulated range. Examples: 1.0 mm needs 172 / 210 / >250 gf buccal and 158 / 176 / 199 gf
palatal at 4 / 6 / 8 mm.

Note: `tools/paper_figs.py` rewrites `out/paper/tables.md` from scratch; Tables 5c, 6b and S2 were added to
that file by hand and are not yet generated by the script.


## 47. Other patients: set-up, build check and placeability (runs in progress)

Replication of the study beyond Open-Full-Jaw patient 8 (manuscript limitation 1).
`bash pipeline/patient_build.sh N` builds patient N in its own workspace `../ofj4d_patients/PN/` (code linked,
own `data/`, `build/`, `build4d/`): stage 1-2 (`ofj/build_frames.py`, `build_attachments.py`), `slim.pkl` from
stage 2 (`pipeline/patient_slim.py`), `geom4d`, `smooth_pdl` -> `geom_fe`, `volmesh_aligner` -> `aligner_fe`,
and the six miniscrew sites of 33. Nothing is re-fitted: alpha 0.40, K_N 1000, 150 gf, rigid bone, as in 34.

**Build check.** Patient 8 rebuilt from the raw dataset this way reproduces `geom.pkl`, `geom_fe.pkl`,
`aligner_fe.pkl`, `slim.pkl`, `tad_bone.pkl` and all five `tad_positions` files exactly (every array
identical; only metadata keys added by later tool versions differ), and a 2-step buccal run in that
workspace matches `v2_buccal_tad_baseline` to 1e-11 mm. Build time ~1-2 min per patient.

**Eligibility.** A patient needs a maxilla with FDI 17-11, 21 and 26 (arch frame, midsagittal plane).
Patients 1, 2, 7, 9, 10: mandible only. 12: no 26. 13, 14: no 16. 15: no 11. 17: has 18 (not in the
paper's design). 11: has 28, which `ofj/build_frames.py` has no crown height for (not attempted further).
16: built, but 8 folded PDL prisms remain after `smooth_pdl` (0 on every other patient) - excluded.
**Built and queued: patients 3, 5, 6, 4** (`../ofj4d_patients/queue_order.txt`).

| patient | 16-26 width (mb cusps) | half-arch length | PDL root area 17 / 16 / 13 / 11 (mm2) | folded PDL prisms before -> after repair | aligner tets |
|---|---|---|---|---|---|
| 8 | 48.3 mm | 47.8 mm | 250 / 313 / 166 / 92 | 29 -> 0 | 134,561 |
| 3 | 56.3 mm | 55.2 mm | 208 / 288 / 178 / 144 | 16 -> 0 | 157,276 |
| 4 | 47.9 mm | 50.1 mm | 229 / 279 / 179 / 155 | 4 -> 0 | 140,155 |
| 5 | 49.9 mm | 51.4 mm | 336 / 370 / 181 / 138 | 3 -> 0 | 145,432 |
| 6 | 53.6 mm | 51.3 mm | 223 / 333 / 196 / 147 | 14 -> 0 | 147,202 |

**Placeability** (`tools/patient_anatomy.py`, criterion of 33: >= 0.5 mm screw-to-PDL on both sides and
insertion outside the root envelope; largest diameter of 1.2 / 1.5 / 1.8 mm; B = buccal, P = palatal;
`../ofj4d_patients/anatomy.json`):

| site | P8 B / P | P3 B / P | P4 B / P | P5 B / P | P6 B / P |
|---|---|---|---|---|---|
| 16/17, 4 mm (paper) | no / 1.2 | 1.8 / 1.8 | no / no (roots) | 1.2 / no | no / no |
| 16/17, 2 mm | no (roots) / no (roots) | 1.8 / 1.5 | no (roots) / no | 1.2 / no (roots) | no / no |
| 16/17, 6 mm | 1.2 / 1.5 | 1.8 / 1.8 | no / no (roots) | no / no | no / 1.8 |
| 16/17, 8 mm | 1.8 / 1.8 | 1.8 / 1.8 | no / no (roots) | no / no | no / no (roots) |
| 15/16, 4 mm | no (roots) / no (roots) | 1.8 / 1.8 | 1.5 / no (roots) | no / no (roots) | no / no (roots) |
| 14/15, 4 mm | no / 1.2 | 1.8 / 1.8 | 1.8 / 1.8 | 1.8 / 1.8 | 1.8 / 1.8 |

Placeability is strongly patient-specific: the 16/17 space - the source study's site - admits a
miniscrew on both sides at some height only in patients 3 and 8, while **14/15 at 4 mm is placeable on both
sides in all four new patients** (gap 3.0-4.2 mm). Section 33's "only 6 and 8 mm are placeable" is a
patient 8 statement, not a general one.

**Runs** (`bash pipeline/patient_queue.sh N`, chained by `build4d/v2/supervisor.sh` v6 under the guard):
17 per patient - control, the six sites x both groups (the v2 set), and the first-premolar cut at the paper
site and 6 mm x both groups. ~32 s per step at 2 in parallel: ~5.5 h per patient. Summary across patients:
`python tools/compare_patients.py` (ranking per patient against patient 8's order, site effects, premolar
minus canine, palatal minus buccal).

**Status 2026-09-16 12:03 UTC (paused, machine released):** patient 3, 2 of 17 complete, 0 capped solves:
control FDI 17 @70 0.396 mm, incisor 0.880 mm (patient 8: 0.052 / 1.183); buccal paper site 1.222 / 1.747 mm
(patient 8: 0.881 / 1.811). Miniscrew benefit on FDI 17 +0.83 mm in both patients. Buccal 2 mm and 6 mm checkpointed
at step 22; everything else not started. The rankings cannot be judged yet. Resume: docs/HANDOFF.md.

## 48. The replication finished, and what splits palatal from buccal is the elastic, not the roots

**Status 2026-09-17 17:13 UTC: complete.** Patients 3, 5, 6, 4 at 17/17 (68 runs, no failed or
unconverged finals, no solver fallbacks, 0 capped solves anywhere in the four workspaces). Tables and the
per-run appendix: `out/paper/replication/` (`python tools/replication_tables.py`); manuscript-facing
summary and drafted edits: `docs/REPLICATION.md`; raw comparison `compare_all5.json`.

**Replicates** (Tables R3, R5, R6): the site ranking `2 mm > paper > 6 mm > 8 mm > 14/15` in 9 of 10
patient-group rankings; the first-premolar precision cut with a palatal screw in 5 of 5; the miniscrew
benefit over control, +0.62 to +0.92 mm on FDI 17, in 5 of 5. The one exception is patient 5's buccal
group, where 8 mm and 14/15 differ by 6 um - at the discretisation error of §36, i.e. a tie.

**Does not replicate:** palatal > buccal (patients 8 and 5 palatal, 3 and 4 buccal, 6 a tie) and the
absolute movements (control FDI 17 from +0.40 mm in patient 3 to -0.66 mm in patient 6; paper-site
distalisation 1.22 to 0.22 mm). Consistent with §14, where that ordering inverted before calibration.

**New: the palatal-buccal difference tracks screw height, monotonically, in all five patients** (Table
R4): 2 mm < 4 mm < 6 mm < 8 mm, about +0.06 mm per 2 mm of height, so three of the four patients that
favour buccal at the paper site cross over by 8 mm. The mesial 14/15 site favours palatal in all five.
Over the 30 paired runs the strongest within-patient predictor is the vertical component of the *buccal*
elastic (r = -0.83): raising the screw steepens that pull and costs it distalisation, while the palatal
line of action changes less. Within-patient spread across the six sites (SD 0.077 mm) is nearly the
between-patient spread of the patient means (SD 0.112 mm).

**Molar root area is not the cause, though it correlates.** Root area of FDI 17 against the patient mean
difference gives r = +0.86 (p = 0.063; r = +0.93 at the paper site alone), but it is collinear with the
geometry it would have to beat (root 17 vs the palatal elastic's transverse component r = +0.92; buccal
free length r = -0.79), and twelve predictors were screened at n = 5. Tested directly by scaling the PDL
secant modulus of FDI 17 and 16 by 0.7 and 1.4 at the paper site with the geometry fixed (new knob
`OFJ_PDL_E_SCALE`, default 1.0 = unchanged and verified against the smoke-test reference;
`build4d/v2/jobs_pdlscale.txt`, 8 runs, Table R7):

| | palatal - buccal, FDI 17 @70 | |
|---|---|---|
| molar PDL scale | patient 3 | patient 5 |
| 0.7 | -0.206 | +0.266 |
| 1.0 (production) | -0.161 | +0.216 |
| 1.4 | -0.126 | +0.185 |

The two patients move in **opposite** directions (+0.080 and -0.081 mm) where the correlation requires
both to rise; neither changes sign over a two-fold change in anchorage; and the effect is ~0.08 mm
against the 0.377 mm gap it would have to explain. The knob itself behaves (a stiffer molar ligament
reduces distalisation in both groups, patient 3 buccal 1.425 -> 1.066 mm), so this is a real null.
**Conclusion: the buccal-palatal ordering is set by the elastic's line of action, not by molar
anchorage.** A larger series is what would settle which anatomical measure, if any, predicts it.

## Open items

- The active-set loop needs 22-35 iterations and does not strictly converge (§9);
  it no longer dominates the cost, which is now ~0.44 s per iteration, 60% of it the
  host factorisation (§10).
- Relative sliding reaches 2-3 mm over 30 steps, which is **not** small compared
  to the 0.45 mm element size, so the small-sliding assumption is strained late in
  the run. The paper makes the same assumption.
- ~~Expression rate here (~90-100% per step) is much higher than the paper's
  37-68% overall.~~ **Resolved by the 70-step run (§11):** that figure was an
  artefact of 30-step truncation. Over the full protocol expression runs 80% at the
  second molar down to 28% at the incisors, consistent with the paper. Movement
  *pattern* and the relative comparison between groups remain the most robust
  outputs.
- PDL Poisson's ratio (0.45) and the choice of von Mises equivalent strain as the
  bilinear driver are ASSUMED — neither paper states them. Poisson's ratio is now
  bounded: 0.30 and 0.49 leave every site ranking unchanged (35).
- Open for the manuscript: other patients - built and running (47); placeability already differs
  (the half arch still judges the right side only - the model is a half arch with midsagittal symmetry, as in the
  source study); PDL stress contour maps (the relative 150 vs 250 gf result is 45, and the
  per-element fields are saved in `build4d/v2/stress/*.npz`); and the clinical citations
  marked [verify] in `docs/manuscript.md`.
