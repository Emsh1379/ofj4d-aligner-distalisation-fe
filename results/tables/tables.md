**Table 1. Model parameters. Items in bold are not specified by the source papers**

| item | value | source |
|---|---|---|
| dentition | right half arch, FDI 17-11, midsagittal symmetry | Open-Full-Jaw patient 8 |
| teeth | rigid | Mao et al. |
| alveolar bone | rigid (deformable variant tested, Table 9) | Mao et al. |
| PDL | 0.30 mm shell, 3 element layers, bilinear E1 = 0.05 / E2 = 0.22 MPa at 7.5 % strain | Mao et al. |
| PDL Poisson ratio | 0.45 | **assumed** (tested 0.30-0.49) |
| aligner | 0.70 mm, E = 1500 MPa, nu = 0.30 | Mao et al. |
| attachments | 5: horizontal on 17/16, vertical on 15/14/13, E = 20 GPa | Mao et al. |
| contact | node-to-surface, pairing fixed per step, Coulomb mu = 0.20 | Mao et al. (tested 0.10-0.30) |
| contact penalty | 1000 N/mm | **numerical** (converged; tested 1000-4000) |
| elastic | 150 gf, canine precision cut to the miniscrew | Mao et al. (100-250 gf, other hooks tested) |
| staging | 70 steps, V pattern, 0.1 mm/step distal, 0.15 mm/step palatal | Mao et al. (other protocols tested) |
| socket kinematics | blend alpha = 0.40 between bodily movement and tipping about the apical third | **calibrated on the control group** |
| iteration | 2 PDL (bone-remodelling) iterations per aligner step | Hamanaka et al. |

**Table 2. Clinical placeability on this patient (criterion: at least 0.5 mm from the screw surface to the PDL on both sides, and insertion outside the root envelope)**

| site | side | interradicular gap | clearance per side, 1.2 / 1.5 / 1.8 mm screw | verdict |
|---|---|---|---|---|
| 16/17, 4 mm* | buccal | 1.85 mm | +0.33 / +0.18 / +0.03 mm | not placeable (clearance < 0.5 mm) |
| 16/17, 4 mm* | palatal | 2.48 mm | +0.64 / +0.49 / +0.34 mm | placeable, <= 1.2 mm |
| 15/16, 4 mm | buccal | 1.33 mm | +0.07 / -0.08 / -0.23 mm | not placeable (clearance and root envelope) |
| 15/16, 4 mm | palatal | 1.45 mm | +0.13 / -0.02 / -0.17 mm | not placeable (clearance and root envelope) |
| 14/15, 4 mm | buccal | 2.14 mm | +0.47 / +0.32 / +0.17 mm | not placeable (clearance < 0.5 mm) |
| 14/15, 4 mm | palatal | 2.21 mm | +0.51 / +0.36 / +0.21 mm | placeable, <= 1.2 mm |
| 16/17, 2 mm | buccal | 1.95 mm | +0.38 / +0.23 / +0.08 mm | not placeable (inside the root envelope) |
| 16/17, 2 mm | palatal | 2.63 mm | +0.72 / +0.57 / +0.42 mm | not placeable (inside the root envelope) |
| 16/17, 6 mm | buccal | 2.29 mm | +0.55 / +0.40 / +0.25 mm | placeable, <= 1.2 mm |
| 16/17, 6 mm | palatal | 2.59 mm | +0.69 / +0.54 / +0.39 mm | placeable, <= 1.5 mm |
| 16/17, 8 mm | buccal | 2.89 mm | +0.85 / +0.70 / +0.55 mm | placeable, <= 1.8 mm |
| 16/17, 8 mm | palatal | 3.64 mm | +1.22 / +1.07 / +0.92 mm | placeable, <= 1.8 mm |

**Table 3. Miniscrew position: tooth movement at step 70 (mm), 150 gf, canine hook**

| site | Buccal miniscrew: FDI 17 distal | Buccal miniscrew: FDI 16 distal | Buccal miniscrew: FDI 11 retraction | Palatal miniscrew: FDI 17 distal | Palatal miniscrew: FDI 16 distal | Palatal miniscrew: FDI 11 retraction |
|---|---|---|---|---|---|---|
| 16/17 2 mm | 1.052 | 1.463 | 2.125 | 1.025 | 1.527 | 2.193 |
| 16/17 4 mm* | 0.881 | 1.295 | 1.811 | 0.952 | 1.439 | 2.003 |
| 16/17 6 mm | 0.726 | 1.144 | 1.529 | 0.862 | 1.336 | 1.794 |
| 16/17 8 mm | 0.573 | 0.995 | 1.260 | 0.765 | 1.229 | 1.581 |
| 15/16 4 mm | 0.637 | 1.056 | 1.364 | 0.840 | 1.317 | 1.762 |
| 14/15 4 mm | 0.348 | 0.773 | 0.849 | 0.612 | 1.065 | 1.271 |

**Table 4. Robustness of the site ranking: FDI 17 distal movement at step 70 (mm). Sites are listed best to worst; "order preserved" means the ranking is unchanged (all 10 comparisons: yes)**

| input set | group | 16/17 2 mm | 16/17 4 mm* | 16/17 6 mm | 16/17 8 mm | 14/15 4 mm | order preserved |
|---|---|---|---|---|---|---|---|
| production (alpha 0.40, nu 0.45) | Buccal miniscrew | 1.052 | 0.881 | 0.726 | 0.573 | 0.348 | yes |
| production (alpha 0.40, nu 0.45) | Palatal miniscrew | 1.025 | 0.952 | 0.862 | 0.765 | 0.612 | yes |
| alpha 0 | Buccal miniscrew | 0.118 | -0.034 | -0.171 | -0.306 | -0.509 | yes |
| alpha 0 | Palatal miniscrew | 0.083 | 0.021 | -0.057 | -0.142 | -0.276 | yes |
| alpha 0.60 | Buccal miniscrew | 1.506 | 1.328 | 1.166 | 1.007 | 0.774 | yes |
| alpha 0.60 | Palatal miniscrew | 1.484 | 1.406 | 1.311 | 1.208 | 1.046 | yes |
| PDL nu 0.30 | Buccal miniscrew | 1.917 | 1.639 | 1.387 | 1.134 | 0.802 | yes |
| PDL nu 0.30 | Palatal miniscrew | 1.949 | 1.821 | 1.661 | 1.488 | 1.225 | yes |
| PDL nu 0.49 | Buccal miniscrew | 0.555 | 0.492 | 0.436 | 0.380 | 0.281 | yes |
| PDL nu 0.49 | Palatal miniscrew | 0.514 | 0.492 | 0.464 | 0.432 | 0.377 | yes |

**Table 5a. Effect of elastic force: mm gained per extra 100 gf (fit over 100-250 gf)**

| quantity | group | 4 mm | 6 mm | 8 mm |
|---|---|---|---|---|
| FDI 17 distal | Buccal miniscrew | +0.554 | +0.456 | +0.357 |
| FDI 17 distal | Palatal miniscrew | +0.596 | +0.539 | +0.476 |
| FDI 16 distal | Buccal miniscrew | +0.463 | +0.367 | +0.272 |
| FDI 16 distal | Palatal miniscrew | +0.559 | +0.492 | +0.421 |
| FDI 11 retraction | Buccal miniscrew | +0.433 | +0.250 | +0.076 |
| FDI 11 retraction | Palatal miniscrew | +0.573 | +0.431 | +0.288 |

**Table 5b. Elastic force needed at a higher site to match the 16/17, 4 mm site at 150 gf**

| quantity | buccal 6 mm | buccal 8 mm | palatal 6 mm | palatal 8 mm |
|---|---|---|---|---|
| FDI 17 distal | 184 gf | 236 gf | 167 gf | 189 gf |
| FDI 16 distal | 191 gf | 260 gf | 171 gf | 200 gf |
| FDI 11 retraction | 262 gf | not attainable | 198 gf | 295 gf |

**Table 6. Hook tooth: FDI 17 distal movement at step 70 (mm), 150 gf**

| group | height | canine* | first premolar | lateral incisor |
|---|---|---|---|---|
| Buccal miniscrew | 4 mm | 0.881 | 0.945 | 0.631 |
| Buccal miniscrew | 6 mm | 0.726 | 0.780 | 0.514 |
| Buccal miniscrew | 8 mm | 0.573 | 0.618 | 0.399 |
| Palatal miniscrew | 4 mm | 0.952 | 1.189 | 0.926 |
| Palatal miniscrew | 6 mm | 0.862 | 1.126 | 0.836 |
| Palatal miniscrew | 8 mm | 0.765 | 1.047 | 0.739 |

**Table 7. Molar overcorrection needed to restore the anchorage of the 16/17, 4 mm site**

| group | height | FDI 17 needs | FDI 16 needs |
|---|---|---|---|
| Buccal miniscrew | 6 mm | 18 % | 9 % |
| Buccal miniscrew | 8 mm | 36 % | 17 % |
| Palatal miniscrew | 6 mm | 11 % | 6 % |
| Palatal miniscrew | 8 mm | 22 % | 12 % |

**Table 8. Staging protocol: FDI 17 distal movement at the end of the protocol (mm)**

| staging | no miniscrew | buccal | palatal | gain from the miniscrew (buccal / palatal) |
|---|---|---|---|---|
| V pattern (70 steps) | 0.052 | 0.881 | 0.952 | 0.83 / 0.90 |
| 2-by-2 (70 steps) | 0.033 | 0.864 | 0.932 | 0.83 / 0.90 |
| sequential (120 steps) | -0.298 | 1.187 | 1.300 | 1.49 / 1.60 |

**Table S1. Mesh convergence at step 20 (movements in um)**

| level | unknowns | FDI 17 distal | FDI 16 distal | FDI 11 retraction | contact force | p99 PDL strain |
|---|---|---|---|---|---|---|
| M0 production | 151,774 | 1831.3 | 1165.8 | 321.2 | 9.054 N | 7.35 % |
| M1 PDL 6 layers | 237,238 | 1831.7 | 1166.1 | 321.7 | 9.038 N | 7.57 % |
| M2 PDL x2 in plane | 319,666 | 1835.0 | 1168.8 | 325.8 | 9.036 N | 7.45 % |
| M3 aligner 0.3 mm | 310,697 | 1832.6 | 1166.8 | 321.6 | 9.202 N | 7.36 % |
| M4 all refinements | 815,891 | 1837.4 | 1171.1 | 329.2 | 9.068 N | 7.65 % |

**Table 9. Sensitivity of the results to the assumed inputs**

| assumption tested | range | largest change in any reported quantity |
|---|---|---|
| aligner-crown friction | 0.10 - 0.30 (model 0.20) | 64 um |
| alveolar bone | rigid vs deformable, both groups | 17 um |
| PDL Poisson ratio | 0.30 - 0.49 (model 0.45) | ranking unchanged; absolute values move (Fig. 4) |
| socket kinematics alpha | 0 - 0.60 (model 0.40) | ranking unchanged; absolute values move (Fig. 4) |
| mesh | 152k - 816k unknowns | 6 um on FDI 17, 8 um on FDI 11 (Table S1) |
| contact penalty | 1000 - 4000 N/mm | < 50 um, with worse convergence at 4000 |

**Table 5c. Elastic force required for a target second-molar distalisation at step 70 (linear force law of Table 5a, interpolated between heights; ">250 gf" is outside the simulated range)**

| target FDI 17 @70 | buccal 4 mm | buccal 6 mm | buccal 8 mm | palatal 4 mm | palatal 6 mm | palatal 8 mm |
|---|---|---|---|---|---|---|
| 0.6 mm | 99 gf | 122 gf | 157 gf | 91 gf | 101 gf | 115 gf |
| 0.8 mm | 135 gf | 166 gf | 213 gf | 125 gf | 139 gf | 157 gf |
| 1.0 mm | 172 gf | 210 gf | >250 gf | 158 gf | 176 gf | 199 gf |
| 1.2 mm | 208 gf | >250 gf | >250 gf | 192 gf | 213 gf | 241 gf |

**Table 6b. Robustness of the precision-cut finding at the 6 mm site: first premolar minus canine, step 70 (mm)**

| input set | buccal FDI 17 | buccal FDI 11 | palatal FDI 17 | palatal FDI 11 |
|---|---|---|---|---|
| production (alpha 0.40, nu 0.45) | +0.055 | +0.025 | +0.264 | +0.344 |
| alpha 0 | +0.050 | +0.028 | +0.239 | +0.325 |
| alpha 0.60 | +0.056 | +0.024 | +0.272 | +0.350 |
| PDL nu 0.30 | +0.083 | +0.073 | +0.359 | +0.488 |
| PDL nu 0.49 | +0.012 | -0.014 | +0.116 | +0.165 |

**Table S2. Periodontal ligament stress at 150 and 250 gf, 16/17 4 mm site (kPa; hydrostatic 5th percentile = compression, for the tooth being moved at that step; replays reproduce the stored runs to 2e-5 mm)**

| group | step (tooth) | hydrostatic p5, 150 -> 250 gf | cervical third | von Mises p95 | contact force (N) |
|---|---|---|---|---|---|
| buccal | 10 (FDI 17) | -43.4 -> -47.7 | -58.8 -> -63.8 | 15.9 -> 17.4 | 8.50 -> 10.02 |
| buccal | 20 (FDI 16) | -25.5 -> -27.4 | -37.8 -> -41.0 | 13.4 -> 14.1 | 14.95 -> 16.25 |
| buccal | 40 (FDI 15) | -38.4 -> -42.0 | -43.6 -> -47.0 | 18.2 -> 19.5 | 3.03 -> 3.46 |
| buccal | 70 (FDI 11) | -49.6 -> -53.1 | -53.3 -> -55.1 | 20.6 -> 21.6 | 5.10 -> 6.31 |
| palatal | 10 (FDI 17) | -43.5 -> -47.5 | -58.4 -> -63.1 | 15.8 -> 17.2 | 7.84 -> 8.98 |
| palatal | 20 (FDI 16) | -26.1 -> -28.3 | -38.5 -> -42.0 | 13.8 -> 14.6 | 13.87 -> 14.68 |
| palatal | 40 (FDI 15) | -39.4 -> -43.6 | -44.5 -> -48.5 | 18.3 -> 19.9 | 3.72 -> 4.42 |
| palatal | 70 (FDI 11) | -50.0 -> -53.1 | -55.3 -> -58.4 | 21.1 -> 22.3 | 5.06 -> 6.38 |
| **peak anywhere** | buccal / palatal | -51.5 -> -53.7 (x1.04) / -52.5 -> -55.5 (x1.06) | -62.5 -> -65.2 / -63.5 -> -67.1 | | |



