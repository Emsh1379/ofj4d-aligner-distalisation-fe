#!/usr/bin/env bash
# Keeps the queued GPU work (build4d/v2/supervisor.sh) running on an interruptible machine:
# starts the supervisor if it is not up, restarts it after a preemption or a stall, and exits once
# the work is done.
#
#   setsid nohup bash build4d/v2/guard.sh >/dev/null 2>&1 &
#
# No tmux: the supervisor is a direct child, so nothing depends on a tmux server
# staying alive. (An earlier version ran both under tmux and killed itself --
# `kill-session -t ofj4d` prefix-matches `ofj4dguard` when no exact match exists.)
#
# Single instance only: run4d.py checkpoints per step, so two supervisors sharing a
# checkpoint would corrupt it. A pid file enforces this; a second run exits immediately.
#
# The supervisor is started with setsid, so it and every queue/xargs/worker under it share one
# process group; restarts kill that whole group, which leaves no orphaned worker pools behind.
set -uo pipefail

ROOT=/teamspace/studios/this_studio/pipeline_data/ofj4d_final_pipeline
export PATH=/home/zeus/miniconda3/envs/cloudspace/bin:$PATH
export HOME=${HOME:-/teamspace/studios/this_studio}
export OFJ_CUDSS_DIR=$HOME/.pylibs_cudss

cd "$ROOT" || exit 1

GUARD_LOG=build4d/v2/guard.log
SUP_LOG=build4d/v2/supervisor_run.log
PIDFILE=build4d/v2/guard.pid
SUP_PAT='^bash build4d/v2/supervisor\.sh'
STALL_SEC=${OFJ_STALL_SEC:-1800}   # no log written for this long while the supervisor is up -> restart

log() { echo "[$(date -u +%FT%TZ)] $*" >> "$GUARD_LOG"; }

# The end marker is read from supervisor.sh on every check, so swapping in a longer chain
# (with a bumped version) needs no guard restart; an older chain's marker in the log never matches.
done_mark() { grep -oE 'SUPERVISOR_DONE v[0-9]+' build4d/v2/supervisor.sh 2>/dev/null | tail -n 1; }

# A pid file rather than flock: a lock fd is inherited by every child (the supervisor,
# its workers, even a plain `sleep`), and each one holds the lock open, so a guard that
# died could not be replaced while the run it started was still alive. Liveness is
# checked against the command line, so a recycled pid does not read as a live guard.
guard_running() {
  local p
  [ -f "$PIDFILE" ] || return 1
  p=$(cat "$PIDFILE" 2>/dev/null)
  [ -n "$p" ] && [ "$p" != "$$" ] || return 1
  kill -0 "$p" 2>/dev/null || return 1
  tr '\0' ' ' < "/proc/$p/cmdline" 2>/dev/null | grep -q 'guard\.sh'
}

# The check-and-claim is serialised by a short flock (held only inside this subshell, so no child
# inherits it): two guards launched at the same instant cannot both pass the check.
if ! ( flock -x 9; guard_running && exit 1; echo $$ > "$PIDFILE"; exit 0 ) 9> build4d/v2/guard.lock; then
  log "guard: another instance (pid $(cat "$PIDFILE")) is running; exiting"
  exit 0
fi

mark0=$(done_mark); [ -n "$mark0" ] && grep -q "$mark0" "$SUP_LOG" 2>/dev/null && exit 0

log "guard: up (boot $(cat /proc/sys/kernel/random/boot_id 2>/dev/null)), waiting for '$(done_mark)'"

# The GPU is not always ready the instant the machine boots.
for _ in $(seq 1 30); do
  nvidia-smi >/dev/null 2>&1 && break
  sleep 10
done

done_count() { grep -l "RUN_EXIT 0" "$@" 2>/dev/null | wc -l; }

patient_progress() {
  local out="" n d
  for n in $(cat ../ofj4d_patients/queue_order.txt 2>/dev/null); do
    d=../ofj4d_patients/P$n/build4d/v2/logs
    out="$out P$n $(ls "$d"/*.log 2>/dev/null | xargs -r grep -l 'RUN_EXIT 0' 2>/dev/null | wc -l)/17"
  done
  echo "${out# }"
}

newest_log() {
  find build4d/v2/logs ../ofj4d_patients/P*/build4d/v2/logs -name '*.log' -printf '%T@\n' 2>/dev/null \
    | sort -n | tail -n 1 | cut -d. -f1
}

stop_supervisor() {
  local p pg
  for p in $(pgrep -f "$SUP_PAT"); do
    pg=$(ps -o pgid= -p "$p" 2>/dev/null | tr -d ' ')
    [ -n "$pg" ] && [ "$pg" != "$(ps -o pgid= -p $$ | tr -d ' ')" ] && kill -TERM -- -"$pg" 2>/dev/null
  done
  sleep 10
  pkill -KILL -f '[r]un4d\.py' 2>/dev/null
  pkill -KILL -f 'xargs -P [0-9]+ -L 1 bash -c run_job' 2>/dev/null
}

sup_since=$(date +%s)   # when the current supervisor was (last known to be) started
beat=0
while true; do
  mark=$(done_mark)
  if [ -n "$mark" ] && grep -q "$mark" "$SUP_LOG" 2>/dev/null; then
    log "guard: $mark - all queued work finished; exiting"
    break
  fi

  # Heartbeat every ~10 min: a throughput record that outlives each preemption.
  if [ $((beat % 10)) -eq 0 ]; then
    log "guard: patients [$(patient_progress)], workers $(pgrep -fc '[r]un4d\.py'), gpu $(nvidia-smi --query-gpu=memory.used,utilization.gpu --format=csv,noheader 2>/dev/null | tr -d ' ')"
  fi
  beat=$((beat + 1))

  if ! pgrep -f "$SUP_PAT" >/dev/null 2>&1; then
    log "guard: supervisor down -> starting"
    # Nothing is driving the queue now, so clear any orphaned workers before
    # restarting: two queues over one checkpoint would corrupt it. run4d.py
    # checkpoints per step, so this costs at most the step in flight.
    if pgrep -f '[r]un4d\.py' >/dev/null 2>&1 || pgrep -f 'xargs -P [0-9]+ -L 1 bash -c run_job' >/dev/null 2>&1; then
      log "guard: orphaned workers found -> clearing"
      stop_supervisor
    fi
    setsid nohup bash build4d/v2/supervisor.sh >> "$SUP_LOG" 2>&1 < /dev/null &
    sup_since=$(date +%s)
    sleep 20
    if pgrep -f "$SUP_PAT" >/dev/null 2>&1; then
      log "guard: supervisor started (pid $(pgrep -f "$SUP_PAT" | head -n 1))"
    else
      log "guard: WARNING supervisor did not come up; retrying next cycle"
    fi
  else
    # Stall watchdog: every running job writes a log line per step (~1-2 min), so a long silence
    # with the supervisor up means something hung. Restart it; checkpoints make this cheap.
    now=$(date +%s); last=$(newest_log); last=${last:-0}
    [ "$last" -lt "$sup_since" ] && last=$sup_since
    if [ $((now - last)) -gt "$STALL_SEC" ]; then
      log "guard: STALL - no log written for $((now - last)) s -> restarting supervisor"
      stop_supervisor
      continue
    fi
  fi

  sleep 60
done
