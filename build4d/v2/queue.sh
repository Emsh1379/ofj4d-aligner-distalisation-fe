#!/usr/bin/env bash
# Resumable job queue for staging runs (production rerun, robustness, mesh convergence).
#   usage: bash build4d/v2/queue.sh JOBFILE [PARALLEL]      (default PARALLEL 2)
# JOBFILE lines:  NAME GROUP STEPS [KEY=VALUE ...]      ('#' comments allowed)
# Every job gets the v2 defaults below; KEY=VALUE overrides them.
# Results: build4d/v2/runs/NAME.pkl, log build4d/v2/logs/NAME.log (ends with RUN_EXIT <code>).
# Resumable: rerun the same command; finished jobs return immediately (run4d resumes from its checkpoint).
set -uo pipefail
cd "$(dirname "$0")/../.."
JOBS=${1:?jobfile}; PAR=${2:-2}
# Cap each job's MKL/PARDISO thread count so PAR jobs share the core budget instead
# of each defaulting to every core (which just adds contention, not throughput,
# once PAR > 1). Override with OFJ_THREADS_PER_JOB=N.
THREADS=${OFJ_THREADS_PER_JOB:-$(( $(nproc) / PAR ))}
[ "$THREADS" -lt 1 ] && THREADS=1
export THREADS
run_job() {
  local name=$1 group=$2 steps=$3; shift 3
  local log=build4d/v2/logs/$name.log out=build4d/v2/runs/$name.pkl
  # Finished already (checkpoint + final step logged + last exit 0): skip without starting
  # python, so a restart reaches unfinished work in seconds instead of minutes.
  if [ -s "$out" ] && grep -qE "\] step +$steps/$steps " "$log" 2>/dev/null \
     && [ "$(grep -o 'RUN_EXIT [0-9]*' "$log" | tail -n 1)" = "RUN_EXIT 0" ]; then
    return 0
  fi
  env OFJ_GEOM=build4d/geom_fe.pkl OFJ_ALIGNER=build4d/aligner_fe.pkl OFJ_K_N=1000 \
      OFJ_PRESCRIBE=tip OFJ_TIP_BLEND=0.40 OFJ_BONE=rigid \
      OMP_NUM_THREADS=$THREADS MKL_NUM_THREADS=$THREADS PYTORCH_CUDA_ALLOC_CONF=${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True} "$@" \
      python -u ofj4d/run4d.py --group "$group" --steps "$steps" --out "$out" 2>&1 \
    | grep --line-buffered -v -E "Resetting the LHS|UserWarning|sparse_csr_tensor|return t.sparse" >> "$log"
  echo "RUN_EXIT ${PIPESTATUS[0]} $(date -u +%FT%TZ)" >> "$log"
}
export -f run_job
grep -v -E '^\s*(#|$)' "$JOBS" | xargs -P "$PAR" -L 1 bash -c 'run_job "$@"' _
echo "QUEUE_DONE $(date -u +%FT%TZ)" >> "${JOBS%.txt}.done"
