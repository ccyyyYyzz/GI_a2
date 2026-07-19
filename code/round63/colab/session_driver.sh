#!/usr/bin/env bash
# ROUND63 session driver -- runs one Colab session's worth of shards, N_LANES at
# a time, with a keep-alive touch loop and a session-level status JSON.
#
# Launched (detached, via setsid) inside a single Colab session:
#   setsid bash session_driver.sh <plan_file> --bundle-root <BDIR> \
#          --session pro2a --status-json <BDIR>/results/round63/status/pro2a.json &
#
# It reads a lane plan (one shard_id per line; '#'/blank ignored), then keeps up
# to --n-lanes remote_lane.py workers busy through a simple queue: when a lane
# frees, the next queued shard starts. remote_lane.py pins BLAS and writes each
# shard's heartbeat; this driver aggregates those into a session status JSON
# every --status-interval seconds (per-lane: shard / heartbeat age / cells done).
#
# Every shard is resumable (shard_runner is META-as-truth), so on a VM recycle
# just relaunch this driver with the SAME plan: completed shards' CSV+meta (if
# restored into the bundle) are skipped and in-flight shards resume.
set -uo pipefail

usage() {
  cat <<'EOF'
usage: session_driver.sh <plan_file> --bundle-root DIR [options]

  <plan_file>              file listing shard_ids (one per line; # comments ok)

required:
  --bundle-root DIR        unpacked bundle root (== repo-root layout)

options:
  --session NAME           session label for status/logs (default: basename of plan)
  --n-lanes N              concurrent remote_lane workers (default: 6)
  --status-json PATH       session status JSON path
                           (default: <bundle-root>/results/round63/status/<session>.json)
  --python PATH            python interpreter (default: python3)
  --wall-budget-s S        per-shard wall budget passed through (default: 21600)
  --status-interval S      seconds between status writes (default: 120)
  --keepalive-interval S   seconds between keep-alive touches (default: 60)
  --heartbeat-interval S   seconds between per-lane heartbeats (default: 60)
  --poll-interval S        scheduler poll granularity (default: 5)
  --keepalive-file PATH    file touched by the keep-alive loop
                           (default: <bundle-root>/results/round63/status/<session>.keepalive)
  -h, --help               show this help
EOF
}

# ------------------------------------------------------------------ args ---- #
PLAN=""
BUNDLE_ROOT=""
SESSION=""
N_LANES=6
STATUS_JSON=""
PYTHON="python3"
WALL_BUDGET_S=21600
STATUS_INTERVAL=120
KEEPALIVE_INTERVAL=60
HEARTBEAT_INTERVAL=60
POLL_INTERVAL=5
KEEPALIVE_FILE=""

while [ $# -gt 0 ]; do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    --bundle-root) BUNDLE_ROOT="$2"; shift 2 ;;
    --session) SESSION="$2"; shift 2 ;;
    --n-lanes) N_LANES="$2"; shift 2 ;;
    --status-json) STATUS_JSON="$2"; shift 2 ;;
    --python) PYTHON="$2"; shift 2 ;;
    --wall-budget-s) WALL_BUDGET_S="$2"; shift 2 ;;
    --status-interval) STATUS_INTERVAL="$2"; shift 2 ;;
    --keepalive-interval) KEEPALIVE_INTERVAL="$2"; shift 2 ;;
    --heartbeat-interval) HEARTBEAT_INTERVAL="$2"; shift 2 ;;
    --poll-interval) POLL_INTERVAL="$2"; shift 2 ;;
    --keepalive-file) KEEPALIVE_FILE="$2"; shift 2 ;;
    --*) echo "unknown option: $1" >&2; usage; exit 2 ;;
    *) if [ -z "$PLAN" ]; then PLAN="$1"; shift
       else echo "unexpected arg: $1" >&2; usage; exit 2; fi ;;
  esac
done

[ -n "$PLAN" ] || { echo "error: <plan_file> required" >&2; usage; exit 2; }
[ -f "$PLAN" ] || { echo "error: plan file not found: $PLAN" >&2; exit 2; }
[ -n "$BUNDLE_ROOT" ] || { echo "error: --bundle-root required" >&2; usage; exit 2; }
[ -d "$BUNDLE_ROOT" ] || { echo "error: bundle root not found: $BUNDLE_ROOT" >&2; exit 2; }
[ -z "$SESSION" ] && SESSION="$(basename "$PLAN" .txt | sed 's/^session_//')"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REMOTE_LANE="$HERE/remote_lane.py"
# prefer the bundle-local copy of the worker if present
[ -f "$BUNDLE_ROOT/code/round63/colab/remote_lane.py" ] && \
  REMOTE_LANE="$BUNDLE_ROOT/code/round63/colab/remote_lane.py"
[ -f "$REMOTE_LANE" ] || { echo "error: remote_lane.py not found ($REMOTE_LANE)" >&2; exit 2; }

STATUS_DIR="$BUNDLE_ROOT/results/round63/status"
HB_DIR="$BUNDLE_ROOT/results/round63/heartbeats"
LANE_LOG_DIR="$BUNDLE_ROOT/results/round63/logs"
STATE_DIR="$STATUS_DIR/.state_$SESSION"
mkdir -p "$STATUS_DIR" "$HB_DIR" "$LANE_LOG_DIR" "$STATE_DIR"
[ -z "$STATUS_JSON" ] && STATUS_JSON="$STATUS_DIR/session_${SESSION}.json"
[ -z "$KEEPALIVE_FILE" ] && KEEPALIVE_FILE="$STATUS_DIR/${SESSION}.keepalive"

# ------------------------------------------------------------- read plan ---- #
QUEUE=()
while IFS= read -r line || [ -n "$line" ]; do
  line="${line%%#*}"                       # strip trailing comment
  line="$(echo "$line" | tr -d '[:space:]')"
  [ -n "$line" ] && QUEUE+=("$line")
done < "$PLAN"
N_TOTAL="${#QUEUE[@]}"
[ "$N_TOTAL" -gt 0 ] || { echo "error: no shards in plan $PLAN" >&2; exit 2; }
echo "[driver $SESSION] $N_TOTAL shards, $N_LANES lanes, bundle=$BUNDLE_ROOT" >&2

# ------------------------------------------------------- keep-alive loop ---- #
keepalive_loop() {
  while true; do
    date -u +"[keepalive $SESSION %Y-%m-%dT%H:%M:%SZ]"
    : > "$KEEPALIVE_FILE" 2>/dev/null || true
    sleep "$KEEPALIVE_INTERVAL"
  done
}
keepalive_loop & KA_PID=$!

cleanup() {
  kill "$KA_PID" 2>/dev/null || true
  for pid in "${SLOT_PID[@]:-}"; do
    [ -n "${pid:-}" ] && kill "$pid" 2>/dev/null || true
  done
}
trap cleanup EXIT INT TERM

# --------------------------------------------------- heartbeat field read --- #
hb_field() {   # $1 = heartbeat file, $2 = key -> prints number or 'null'/''
  [ -f "$1" ] || { echo ""; return; }
  grep -oE "\"$2\":[[:space:]]*(-?[0-9]+|null)" "$1" 2>/dev/null \
    | head -1 | grep -oE '(-?[0-9]+|null)' | head -1
}

# ------------------------------------------------------- status writer ------ #
COMPLETED_SHARDS=()
COMPLETED_RCS=()
write_status() {
  local now; now="$(date -u +%s)"
  local tmp="$STATUS_JSON.tmp"
  {
    printf '{\n'
    printf '  "session": "%s",\n' "$SESSION"
    printf '  "ts": %s,\n' "$now"
    printf '  "ts_utc": "%s",\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf '  "n_lanes": %s,\n' "$N_LANES"
    printf '  "plan_file": "%s",\n' "$PLAN"
    printf '  "total_shards": %s,\n' "$N_TOTAL"
    printf '  "queued_remaining": %s,\n' "$((N_TOTAL - QHEAD))"
    printf '  "completed_count": %s,\n' "${#COMPLETED_SHARDS[@]}"
    # completed [{shard,rc}]
    printf '  "completed": ['
    local i first=1
    for i in "${!COMPLETED_SHARDS[@]}"; do
      [ "$first" -eq 1 ] || printf ', '
      printf '{"shard": "%s", "rc": %s}' "${COMPLETED_SHARDS[$i]}" "${COMPLETED_RCS[$i]}"
      first=0
    done
    printf '],\n'
    # lanes
    printf '  "lanes": [\n'
    local s shard hb epoch cells rows age state
    for (( s=0; s<N_LANES; s++ )); do
      shard="${SLOT_SHARD[$s]:-}"
      if [ -z "$shard" ]; then
        printf '    {"lane": %s, "shard": "-", "state": "idle"}' "$s"
      else
        # normalize: plan lines may carry manifest PATHS (e.g. results/.../X.json)
        # while remote_lane names heartbeats by bare shard_id -> join on basename.
        hb="$HB_DIR/$(basename "${shard%.json}").hb.json"
        epoch="$(hb_field "$hb" epoch)"
        cells="$(hb_field "$hb" cells_done)"
        rows="$(hb_field "$hb" rows)"
        if [ -n "$epoch" ] && [ "$epoch" != "null" ]; then
          age=$((now - epoch)); state="running"
        else
          age=-1; state="starting"
        fi
        [ -n "$cells" ] || cells=0
        [ -n "$rows" ] || rows=0
        printf '    {"lane": %s, "shard": "%s", "heartbeat_age_s": %s, "cells_done": %s, "rows": %s, "state": "%s"}' \
          "$s" "$shard" "$age" "$cells" "$rows" "$state"
      fi
      [ "$s" -lt "$((N_LANES - 1))" ] && printf ',\n' || printf '\n'
    done
    printf '  ]\n'
    printf '}\n'
  } > "$tmp" && mv -f "$tmp" "$STATUS_JSON"
}

# ------------------------------------------------------- lane launcher ------ #
declare -a SLOT_PID SLOT_SHARD
for (( s=0; s<N_LANES; s++ )); do SLOT_PID[$s]=""; SLOT_SHARD[$s]=""; done

launch_lane() {   # $1 slot, $2 shard
  local slot="$1" shard="$2"
  local rcfile="$STATE_DIR/slot_${slot}.rc"
  rm -f "$rcfile"
  (
    "$PYTHON" "$REMOTE_LANE" \
      --manifest "$shard" \
      --bundle-root "$BUNDLE_ROOT" \
      --wall-budget-s "$WALL_BUDGET_S" \
      --heartbeat-dir "$HB_DIR" \
      --heartbeat-interval "$HEARTBEAT_INTERVAL" \
      --log-dir "$LANE_LOG_DIR" \
      --python "$PYTHON"
    echo "$?" > "$rcfile"
  ) >> "$LANE_LOG_DIR/lane_slot${slot}.log" 2>&1 &
  SLOT_PID[$slot]=$!
  SLOT_SHARD[$slot]="$shard"
  echo "[driver $SESSION] lane $slot <- $shard (pid ${SLOT_PID[$slot]})" >&2
}

reap_slot() {     # $1 slot ; returns 0 if a shard completed on this slot
  local slot="$1"
  local rcfile="$STATE_DIR/slot_${slot}.rc"
  [ -n "${SLOT_SHARD[$slot]:-}" ] || return 1
  [ -f "$rcfile" ] || return 1
  local rc; rc="$(cat "$rcfile" 2>/dev/null)"; [ -n "$rc" ] || rc=1
  wait "${SLOT_PID[$slot]}" 2>/dev/null || true
  COMPLETED_SHARDS+=("${SLOT_SHARD[$slot]}")
  COMPLETED_RCS+=("$rc")
  echo "[driver $SESSION] lane $slot done: ${SLOT_SHARD[$slot]} rc=$rc" >&2
  SLOT_PID[$slot]=""
  SLOT_SHARD[$slot]=""
  rm -f "$rcfile"
  return 0
}

# ------------------------------------------------------------ scheduler ----- #
QHEAD=0
LAST_STATUS=0
write_status
while :; do
  # reap finished lanes
  for (( s=0; s<N_LANES; s++ )); do reap_slot "$s" || true; done
  # fill free lanes from the queue
  for (( s=0; s<N_LANES; s++ )); do
    if [ -z "${SLOT_SHARD[$s]:-}" ] && [ "$QHEAD" -lt "$N_TOTAL" ]; then
      launch_lane "$s" "${QUEUE[$QHEAD]}"
      QHEAD=$((QHEAD + 1))
    fi
  done
  # done?  queue drained and all lanes idle
  busy=0
  for (( s=0; s<N_LANES; s++ )); do [ -n "${SLOT_SHARD[$s]:-}" ] && busy=1; done
  if [ "$QHEAD" -ge "$N_TOTAL" ] && [ "$busy" -eq 0 ]; then
    break
  fi
  # throttled status write
  now="$(date -u +%s)"
  if [ "$((now - LAST_STATUS))" -ge "$STATUS_INTERVAL" ]; then
    write_status; LAST_STATUS="$now"
  fi
  sleep "$POLL_INTERVAL"
done

write_status
# summarise; non-zero session exit if any shard did not complete (rc!=0)
fail=0
for rc in "${COMPLETED_RCS[@]:-}"; do [ "${rc:-1}" -ne 0 ] && fail=1; done
echo "[driver $SESSION] all lanes drained: ${#COMPLETED_SHARDS[@]}/$N_TOTAL shards, session_fail=$fail" >&2
exit "$fail"
