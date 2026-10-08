#!/usr/bin/env bash
# Second-patient replication: run the core staging set on a patient built by pipeline/patient_build.sh.
#
#   bash pipeline/patient_queue.sh N [PARALLEL]      (default PARALLEL 2; resumable, rerun after a restart)
#
# Same model and settings as the production set of patient 8 (build4d/v2/queue.sh: geom_fe, aligner_fe,
# K_N 1000, alpha 0.40, rigid bone, 150 gf, 70 steps); nothing is re-calibrated to the new anatomy.
# Jobs (17): control; paper site and the 2/6/8 mm, 15/16 and 14/15 variants, both groups (the v2 set);
# first-premolar precision cut at the paper site and at 6 mm, both groups (the hk14 set).
# Summaries, from the workspace: python tools/compare_tad.py --set v2   (and --set hk14)
set -uo pipefail
N=$1; PAR=${2:-2}
PKG=$(cd "$(dirname "$0")/.." && pwd)
W=$PKG/../ofj4d_patients/P$N
[ -f "$W/build4d/geom_fe.pkl" ] || { echo "build patient $N first: bash pipeline/patient_build.sh $N"; exit 1; }
mkdir -p "$W/build4d/v2/runs" "$W/build4d/v2/logs"   # a workspace restored from a zip has no empty dirs
ln -sfn "$PKG/build4d/v2/queue.sh" "$W/build4d/v2/queue.sh"
J=$W/build4d/v2/jobs_patient.txt
{
  echo "# patient $N: v2 production set + first-premolar cut (hk14) at the paper site and 6 mm"
  echo "v2_control control 70"
  for g in buccal_tad palatal_tad; do
    echo "v2_${g}_baseline $g 70 OFJ_TAD_FILE=build4d/tad_bone.pkl"
    for v in h1617_a2 h1617_a6 h1617_a8 md1516_a4 md1415_a4; do
      echo "v2_${g}_$v $g 70 OFJ_TAD_FILE=build4d/tad_positions/tad_$v.pkl"
    done
  done
  for g in buccal_tad palatal_tad; do
    echo "rob_hk14_${g}_baseline $g 70 OFJ_TAD_FILE=build4d/tad_bone.pkl OFJ_HOOK_FDI=14"
    echo "rob_hk14_${g}_h1617_a6 $g 70 OFJ_TAD_FILE=build4d/tad_positions/tad_h1617_a6.pkl OFJ_HOOK_FDI=14"
  done
} > "$J"
export OFJ_CUDSS_DIR=${OFJ_CUDSS_DIR:-$HOME/.pylibs_cudss}
echo "patient $N start $(date -u +%FT%TZ)"
bash "$W/build4d/v2/queue.sh" "$J" "$PAR"
echo "patient $N done $(date -u +%FT%TZ)"
