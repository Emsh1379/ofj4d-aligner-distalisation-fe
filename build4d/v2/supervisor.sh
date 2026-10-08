#!/usr/bin/env bash
# Chains the GPU work, each phase only if its job file exists:
#   production -> robustness -> mesh convergence -> elastic force -> robustness at 6 mm -> hook location
#   -> friction -> deformable bone -> molar overcorrection -> staging -> hook robustness
#   -> other patients (pipeline/patient_queue.sh, each only if that patient is built)
# Resumable: rerun this script after a restart; finished jobs return immediately.
# guard.sh reads the version from the final marker line below, so bump it whenever phases are added:
# the end of an earlier, shorter chain already in the log is then not taken for the end of this one.
cd "$(dirname "$0")/../.."
while pgrep -f '[r]un_grid.sh rigid' > /dev/null; do sleep 30; done
phase() {
  [ -f "$2" ] || return 0
  echo "$1 start $(date -u +%FT%TZ)"
  bash build4d/v2/queue.sh "$2" "$3"
}
phase production build4d/v2/jobs_production.txt 4
phase robustness build4d/v2/jobs_robustness.txt 4
phase "mesh convergence" build4d/v2/jobs_meshconv.txt 1
phase "elastic force" build4d/v2/jobs_force.txt 4
phase "robustness 6 mm" build4d/v2/jobs_rob6mm.txt 4
phase "hook location" build4d/v2/jobs_hook.txt 4
phase friction build4d/v2/jobs_friction.txt 4
phase "deformable bone" build4d/v2/jobs_bone.txt 1
phase "molar overcorrection" build4d/v2/jobs_overcorr.txt 4
phase "staging pattern" build4d/v2/jobs_staging.txt 4
phase "hook robustness" build4d/v2/jobs_hookrob.txt 4
for n in $(cat ../ofj4d_patients/queue_order.txt 2>/dev/null); do
  [ -f ../ofj4d_patients/P$n/build4d/geom_fe.pkl ] || continue
  echo "patient $n start $(date -u +%FT%TZ)"
  bash pipeline/patient_queue.sh "$n" 4
done
# Experiment A (2026-09-17): molar PDL stiffness sweep at the paper screw site, patients 3 and 5.
for n in 3 5; do
  W=$(cd ../ofj4d_patients/P$n 2>/dev/null && pwd) || continue
  [ -f "$W/build4d/v2/jobs_pdlscale.txt" ] || continue
  echo "pdl-scale patient $n start $(date -u +%FT%TZ)"
  bash "$W/build4d/v2/queue.sh" "$W/build4d/v2/jobs_pdlscale.txt" 4
done

# Only declare the chain finished when every queued patient's jobs really exited 0; otherwise exit
# non-zero without the marker, so guard.sh restarts the chain instead of shutting down for good.
incomplete=""
for n in $(cat ../ofj4d_patients/queue_order.txt 2>/dev/null); do
  W=../ofj4d_patients/P$n
  [ -f $W/build4d/geom_fe.pkl ] || continue
  total=0; ok=0
  for J in "$W/build4d/v2/jobs_patient.txt" "$W/build4d/v2/jobs_pdlscale.txt"; do
    [ -f "$J" ] || continue
    n_j=$(grep -cvE '^\s*(#|$)' "$J" 2>/dev/null); total=$((total + ${n_j:-0}))
    while read -r name _; do
      case "$name" in ''|\#*) continue;; esac
      [ "$(grep -o 'RUN_EXIT [0-9]*' "$W/build4d/v2/logs/$name.log" 2>/dev/null | tail -n 1)" = "RUN_EXIT 0" ] && ok=$((ok+1))
    done < "$J"
  done
  { [ "$total" -gt 0 ] && [ "$ok" -eq "$total" ]; } || incomplete="$incomplete P$n($ok/$total)"
done
if [ -n "$incomplete" ]; then
  echo "SUPERVISOR_INCOMPLETE$incomplete $(date -u +%FT%TZ)"
  exit 1
fi
echo "SUPERVISOR_DONE v14 $(date -u +%FT%TZ)"
