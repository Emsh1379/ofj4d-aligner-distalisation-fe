#!/usr/bin/env bash
# Where the remaining pipeline stands. Safe to run any time:
#   bash build4d/v2/status.sh
cd /teamspace/studios/this_studio/pipeline_data/ofj4d_final_pipeline || exit 1

echo "time        $(date -u +%FT%TZ)   (machine up $(uptime -p | sed 's/^up //'))"

if pgrep -f '[g]uard\.sh' >/dev/null 2>&1; then
  echo "guard       running (pid $(pgrep -f '[g]uard\.sh' | head -n 1))"
else
  echo "guard       NOT RUNNING  ->  setsid nohup bash build4d/v2/guard.sh >/dev/null 2>&1 &"
fi

if pgrep -f '^bash build4d/v2/supervisor\.sh' >/dev/null 2>&1; then
  echo "supervisor  running (pid $(pgrep -f '^bash build4d/v2/supervisor\.sh' | head -n 1))"
else
  echo "supervisor  down (the guard restarts it within ~60 s)"
fi

echo "phase       $(tail -n 1 build4d/v2/supervisor_run.log 2>/dev/null)"

L=build4d/v2/logs
n() { grep -l "RUN_EXIT 0" "$@" 2>/dev/null | wc -l; }
echo "complete    robustness $(n $L/rob_a000_*.log $L/rob_a060_*.log $L/rob_nu030_*.log $L/rob_nu049_*.log)/40   mesh convergence $(n $L/mc_*.log)/4"
echo "            elastic force $(n $L/rob_f*.log)/18   hook location $(n $L/rob_hk12_*.log $L/rob_hk14_*.log)/12"
echo "            hook robustness $(n $L/rob_hka*.log $L/rob_hknu*.log)/8"
echo "            friction $(n $L/rob_mu*.log)/4   deformable bone $(n $L/rob_bdef*.log)/2   overcorrection $(n $L/rob_oc*.log)/12   staging $(n $L/rob_st*.log)/6"

echo "in flight:"
pgrep -af '[r]un4d\.py' 2>/dev/null | sed 's/.*--out build4d\/v2\/runs\///; s/\.pkl//' | while read -r name; do
  step=$(grep -E "^\[.*\] step" "build4d/v2/logs/$name.log" 2>/dev/null | tail -n 1 | grep -oE 'step +[0-9]+/[0-9]+')
  printf '  %-38s %s\n' "$name" "${step:-starting}"
done
[ -z "$(pgrep -f '[r]un4d\.py' 2>/dev/null)" ] && echo "  (none)"

echo "recent guard events:"
tail -n 5 build4d/v2/guard.log 2>/dev/null | sed 's/^/  /'
