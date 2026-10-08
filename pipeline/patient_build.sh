#!/usr/bin/env bash
# Build the half-arch model of another Open-Full-Jaw patient in its own workspace.
#
#   bash pipeline/patient_build.sh N [RAW_DIR]
#
# RAW_DIR holds the unzipped Patient_N/ folder (default: ../../ofj_raw). The workspace
# ../ofj4d_patients/PN/ links the code (ofj, ofj4d, tools, pipeline) and keeps its own
# data/, build/ and build4d/, so every relative path in the code resolves to patient N
# and the patient 8 results are never touched. The patient's files sit under data/P8/
# because that directory name is hard-coded in the builders; it is only a path.
#
# Chain: stage 1 frames -> stage 2 attachments (ofj/) -> slim.pkl (inputs) -> geom4d ->
# smooth_pdl (geom_fe) -> volmesh_aligner (aligner_fe) -> miniscrew sites.
set -euo pipefail
N=$1
PKG=$(cd "$(dirname "$0")/.." && pwd)
RAW=${2:-$PKG/../../ofj_raw}
RAW=$(cd "$RAW" && pwd)
W=$PKG/../ofj4d_patients/P$N
mkdir -p "$W"/data/P8 "$W"/data/cache "$W"/build "$W"/build4d/tad_positions "$W"/build4d/v2/runs "$W"/build4d/v2/logs
W=$(cd "$W" && pwd)
for d in ofj ofj4d tools pipeline; do ln -sfn "$PKG/$d" "$W/$d"; done
ln -sfn "$RAW/Patient_$N/output/maxilla/volumetric_meshes/maxilla_volumetric_mesh.msh" "$W/data/P8/"
ln -sfn "$RAW/Patient_$N/input/maxilla/teeth_axes_maxilla.json" "$W/data/P8/"
cd "$W"
export PYTHONUNBUFFERED=1 OFJ_MSH_CACHE=data/cache/maxilla.npz
log=build4d/build.log; : > $log
run() { echo "+ $*" | tee -a $log; "$@" 2>&1 | tee -a $log; }

run python -c "import sys; sys.path.insert(0,'ofj'); import build_frames as B; B.build('data/P8/maxilla_volumetric_mesh.msh','data/P8/teeth_axes_maxilla.json', cache='data/cache/maxilla.npz')"
run python ofj/build_attachments.py
run python pipeline/patient_slim.py
run python ofj4d/geom4d.py
run python tools/smooth_pdl.py
run python tools/volmesh_aligner.py
run python tools/tad_from_bone.py --out build4d/tad_bone.pkl
for v in "16 17 2 h1617_a2" "16 17 6 h1617_a6" "16 17 8 h1617_a8" "14 15 4 md1415_a4" "15 16 4 md1516_a4"; do
  set -- $v
  run python tools/tad_from_bone.py --pair $1 $2 --apical $3 --out build4d/tad_positions/tad_$4.pkl
done
echo "built patient $N in $W" | tee -a $log
