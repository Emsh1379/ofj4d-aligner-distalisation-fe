#!/usr/bin/env bash
# Miniscrew-position study: buccal and palatal TAD groups, final model, 70 steps per run.
#   usage: bash build4d/tad_positions/run_grid.sh BONE [PARALLEL] [VARIANT ...]
#     BONE      rigid | deformable
#     PARALLEL  runs at once (default 2 for rigid, 1 for deformable)
#     VARIANT   default: md1516_a4 md1415_a4 h1617_a2 h1617_a6 h1617_a8
#               (baseline 16/17 4 mm = build4d/tad_bone.pkl; add "baseline" to run it too)
# Resumable: rerun the same command. Results: build4d/tad_positions/run_<group>_<variant>_<bone>.pkl/.log
set -uo pipefail
cd "$(dirname "$0")/../.."
BONE=${1:?rigid or deformable}; shift
PAR=${1:-}; [ -n "$PAR" ] && [[ "$PAR" =~ ^[0-9]+$ ]] && shift || PAR=$([ "$BONE" = rigid ] && echo 2 || echo 1)
VARIANTS=${*:-md1516_a4 md1415_a4 h1617_a2 h1617_a6 h1617_a8}
export BONE OFJ_PRESCRIBE=tip OFJ_TIP_BLEND=0.40
run() {
  local g=$1 v=$2
  local tad=build4d/tad_positions/tad_$v.pkl; [ "$v" = baseline ] && tad=build4d/tad_bone.pkl
  local out=build4d/tad_positions/run_${g}_${v}_${BONE}
  OFJ_BONE=$BONE OFJ_TAD_FILE=$tad python -u ofj4d/run4d.py --group "$g" --steps 70 --out "$out.pkl" 2>&1 \
    | grep --line-buffered -v -E "Resetting the LHS|UserWarning|sparse_csr_tensor|return t.sparse" >> "$out.log"
  echo "RUN_EXIT ${PIPESTATUS[0]}" >> "$out.log"
}
export -f run
for v in $VARIANTS; do for g in buccal_tad palatal_tad; do echo "$g $v"; done; done \
  | xargs -P "$PAR" -L 1 bash -c 'run $0 $1'
echo GRID_DONE >> build4d/tad_positions/grid_${BONE}.done
