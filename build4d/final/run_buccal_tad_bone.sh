#!/usr/bin/env bash
# deformable-bone buccal TAD group, final model, 70 steps (resumable: rerun this script)
cd "$(dirname "$0")/../.."
export OFJ_BONE=deformable OFJ_PRESCRIBE=tip OFJ_TIP_BLEND=0.40
python -u ofj4d/run4d.py --group buccal_tad --steps 70 --out build4d/final/run_buccal_tad_bone.pkl 2>&1 \
  | grep --line-buffered -v -E "Resetting the LHS|UserWarning|sparse_csr_tensor|return t.sparse" >> build4d/final/run_buccal_tad_bone.log
echo "RUN_EXIT ${PIPESTATUS[0]}" >> build4d/final/run_buccal_tad_bone.log
