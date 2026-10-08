# Other-patient replication: what to write, and what not to claim

Companion to `docs/manuscript.md`. Numbers: `out/paper/replication/tables_replication.md` (Tables R1-R6),
raw per-run appendix `out/paper/replication/replication_runs.csv`, regenerate both with
`python tools/replication_tables.py`. Method and build check: `docs/FINDINGS.md` section 47.

Four Open-Full-Jaw maxillae (patients 3, 4, 5, 6) were built with `pipeline/patient_build.sh` and run
through the production set of the single-patient study: 17 runs each, 70 steps, alpha 0.40, K_N 1000,
150 gf, rigid bone. **Nothing was re-calibrated to the new anatomies** - the socket-kinematics blend
stays at the value fitted to patient 8's control group, which is what makes these out-of-sample. All 68
runs completed (2026-09-17), no failed or unconverged finals, no solver fallbacks.

## What replicates

| claim | result across the five patients |
|---|---|
| Site ranking `2 mm > paper > 6 mm > 8 mm > 14/15` | 9 of 10 patient-group rankings identical (Table R3) |
| First-premolar cut beats the canine cut, palatal screw | 5 of 5, both teeth (Table R5) |
| A miniscrew beats no miniscrew | 5 of 5, +0.62 to +0.92 mm on FDI 17 (Table R6) |
| Effect sizes of the sites | 2 mm gains +0.11 to +0.17 mm; 14/15 loses 0.21 to 0.59 mm (Table R2a) |

The single exception is patient 5's buccal group, where 8 mm and 14/15 differ by **6 um** (0.214 vs
0.208 mm) and swap. Mesh convergence (FINDINGS 36) puts discretisation error at up to 6 um on this
quantity, so this is a tie, not a reversal - say so rather than reporting it as a failed replication.

## What does not replicate

**Palatal versus buccal.** Patient 8's palatal advantage is not general: two patients favour palatal
(8, 5), two favour buccal (3, 4), one is a tie (6) (Table R4). This was already the study's most fragile
comparison - FINDINGS 14 records that the ordering inverted before calibration - so the honest reading is
that the model does not resolve which of the two sites is better, while it does resolve the height and
mesiodistal ranking within each side.

**Absolute movement.** The control second molar ends between +0.40 mm (patient 3) and -0.66 mm (patient
6); distalisation at the paper site ranges 1.22 to 0.22 mm, a five-fold spread (Table R2a). Incisor
retraction is steadier (1.6-2.2 mm). Relative comparisons are the transferable output, as the
single-patient discussion already argued; the replication now demonstrates it rather than asserting it.

## The one genuinely new result

**The palatal-buccal difference tracks screw height, in the same direction in every patient.** Table R4
read across: 2 mm < 4 mm < 6 mm < 8 mm, monotonically, in all five (slope about +0.06 mm per 2 mm of
height). Three of the patients that favour buccal at the paper site would cross to palatal at 8 mm. The
mesial site (14/15) likewise favours palatal in all five.

Mechanically this is the elastic's line of action, not the anatomy: raising the screw steepens the buccal
elastic (its vertical component is the strongest within-patient predictor, r = -0.83 over the 30 paired
runs) and costs it distalisation, while the palatal pull changes less. Within-patient spread across the
six sites (SD 0.077 mm) is nearly as large as the spread between patient means (SD 0.112 mm).

## What not to claim: molar root size (tested and rejected)

Second-molar root area tracks the palatal-buccal difference across the five patients (r = +0.86 on the
six-site means, p = 0.063; r = +0.93 at the paper site alone). **It is not the cause.** Scaling the molar
PDL modulus 0.7 to 1.4 at the paper site, geometry held fixed (Table R7, `OFJ_PDL_E_SCALE`), fails the
prediction three ways:

* the two patients move in *opposite* directions (P3 +0.080 mm, P5 -0.081 mm), where the correlation
  requires both to rise with anchorage;
* neither changes sign - P3 stays buccal-favouring and P5 palatal-favouring across a two-fold change;
* the effect is ~0.08 mm against the 0.377 mm gap between these patients it would have to explain.

The correlation is an artefact of n = 5 with twelve candidate predictors tested, and root area is
collinear with the elastic geometry it would have to beat (root area vs the palatal elastic's transverse
component r = +0.92, vs buccal free length r = -0.79). Report it, if at all, as a rejected explanation.

The scaling did behave sensibly in the direction that checks the knob: a stiffer molar ligament reduces
distalisation in both groups (P3 buccal 1.425 -> 1.066 mm), so the null result is not a dead control.

What remains is the elastic's line of action, which the height trend above already supports directly.

## Suggested edits to the manuscript

*Limitations item 1* currently reads "One patient. All anatomy, and therefore all placeability, comes
from a single maxilla. The rankings are robust to the model's assumptions but not demonstrated across
anatomies." Replace with:

> **Five patients, one calibration.** The site ranking, the precision-cut result and the miniscrew
> benefit were reproduced on four further Open-Full-Jaw maxillae without re-calibration (Tables R1-R6).
> Absolute movements vary more than five-fold between them, and which of the buccal and palatal sites is
> better is patient-specific, so only the relative results should be read as general. Placeability was
> judged on each patient's own anatomy; the calibration itself still rests on patient 8's control group.

*Conclusion 1* ("the ranking of miniscrew sites is stable under the model's assumed inputs, while
absolute movements are not") can be strengthened to "stable under the model's assumed inputs and across
five maxillae".

*Conclusion 4* (first-premolar cut with a palatal screw) holds in 5 of 5 and can be stated as replicated.

*Discussion* needs a new paragraph for the palatal-buccal split and the height trend, and the abstract's
group ordering should be softened from a general claim to a patient-8 result.

## Reproducing the numbers

```
python tools/compare_patients.py 8 3 4 5 6      # console summary, --json for the raw dict
python tools/replication_tables.py              # Tables R1-R6 and the CSV appendix
```
