#!/usr/bin/env bash
# Sensitivity of the final model (control group) to the ASSUMED inputs: PDL Poisson's
# ratio, contact penalty stiffness, attachment size.  Resumable: the checkpoints
# shipped in build4d/variants/sensblend_*.pkl continue from the step the original
# machine reached.
#
# usage: bash pipeline/run_sensitivity.sh [PARALLEL]
set -euo pipefail
cd "$(dirname "$0")/.."
PAR=${1:-1}
V=build4d/variants
mkdir -p "$V"
export OFJ_PRESCRIBE=tip OFJ_TIP_BLEND=0.40

# attachment-size variants are rebuilt from build4d/geom.pkl (~1 min each)
for spec in 0.75:A075 1.25:A125; do
  tag=${spec#*:}; scale=${spec%:*}
  [ -f "$V/aligner_$tag.pkl" ] || python tools/att_variant.py "$scale" "$scale" "$tag"
done

# tag  group  settings
cat > "$V/jobs_sens.txt" <<EOF
nu030 control OFJ_PDL_NU=0.30
nu049 control OFJ_PDL_NU=0.49
kn1000 control OFJ_K_N=1000
kn4000 control OFJ_K_N=4000
A075 control OFJ_GEOM=$V/geom_A075.pkl OFJ_ALIGNER=$V/aligner_A075.pkl
A125 control OFJ_GEOM=$V/geom_A125.pkl OFJ_ALIGNER=$V/aligner_A125.pkl
EOF

job() {
  local tag=$1 group=$2; shift 2
  env "$@" python -u ofj4d/run4d.py --group "$group" --steps 70 \
      --out "build4d/variants/sensblend_$tag.pkl" 2>&1 \
    | grep --line-buffered -v "Resetting the LHS" | tee -a "build4d/variants/sensblend_$tag.log"
}
export -f job
xargs -P "$PAR" -L 1 bash -c 'job "$@"' _ < "$V/jobs_sens.txt"

python tools/sensitivity.py      # control group: every assumed input
