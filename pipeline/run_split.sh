#!/usr/bin/env bash
# Split socket kinematics (docs/FINDINGS.md section 27): canine..molars blended at
# alpha 0.55, incisors bodily (alpha 0).  All three groups, then the comparison in the
# paper's coordinate frame.  Resumable: shipped checkpoints continue where they stopped.
#
# usage: bash pipeline/run_split.sh [PARALLEL]
set -euo pipefail
cd "$(dirname "$0")/.."
PAR=${1:-1}
export V=build4d/variants
mkdir -p "$V"
export OFJ_PRESCRIBE=tip OFJ_TIP_BLEND=0.55 OFJ_TIP_BLEND_INCISOR=0

run() {
  local g=$1
  python -u ofj4d/run4d.py --group "$g" --steps 70 --out "$V/split_run_$g.pkl" 2>&1 \
    | grep --line-buffered -v "Resetting the LHS" | tee -a "$V/split_run_$g.log"
}
export -f run
printf '%s\n' control buccal_tad palatal_tad | xargs -P "$PAR" -I{} bash -c 'run {}'

python tools/paper_frame.py build4d/final/run_ "$V/split_run_"
