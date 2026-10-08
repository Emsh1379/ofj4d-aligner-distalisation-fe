#!/usr/bin/env bash
# Final model: Mao et al. 2024 reproduction with blended socket kinematics
# (OFJ_PRESCRIBE=tip, OFJ_TIP_BLEND=0.40; docs/FINDINGS.md sections 24-25).
# Runs the three anchorage groups for 70 steps, then the analysis and figures.
# Resumable: after an interruption, run the same command again.
#
# usage: bash pipeline/run_final.sh [PARALLEL]
#   PARALLEL = groups run at once (default 1; 2 fits a 16 GB T4 with cuDSS)
set -euo pipefail
cd "$(dirname "$0")/.."
PAR=${1:-1}
export OUT=build4d/final
mkdir -p "$OUT"
export OFJ_PRESCRIBE=tip OFJ_TIP_BLEND=0.40

run() {
  local g=$1
  python -u ofj4d/run4d.py --group "$g" --steps 70 --out "$OUT/run_$g.pkl" 2>&1 \
    | grep --line-buffered -v "Resetting the LHS" | tee -a "$OUT/run_$g.log"
}
export -f run
printf '%s\n' control buccal_tad palatal_tad | xargs -P "$PAR" -I{} bash -c 'run {}'

python tools/analyse.py "$OUT/run_"
python tools/plot_compare.py "$OUT" "blended sockets, alpha 0.40" _final
python tools/plot_staging.py "$OUT" _final "blended sockets, alpha 0.40"
