#!/bin/bash
# ROUND63 M1 self-healing fleet watchdog v2 (detached daemon; launch via setsid).
#
# v2: progress is read from VM-side GROUND TRUTH (heartbeat files + driver
# status) via a per-session exec, NOT the driver status JSON's lane table --
# the M1 plans carry manifest PATHS, and session_driver v1 joins heartbeat
# filenames from the raw plan line, so its per-lane progress is blind (the
# underlying computation and completed-shard reaping are unaffected; verified
# on r63m1_a 2026-07-19 15:13Z: lanes 5/16 cells, fresh heartbeats, CSVs+metas).
#
#   cycle (300 s): rebind-first (guide 3.2) -> per-session exec summarizer ->
#   taxonomy -> verbose log; every 6th cycle (30 min) + on any non-OK state,
#   append a digest block to the D:-side digest log.
#
# Taxonomy: COMPLETE done==total | OK | STALL (no heartbeat progress >900 s)
#           TOKEN/CLI (endpoint assigned but exec fails) | DEAD (endpoint gone).
# Outcome-blind: cell/row counts, ages, rc codes only. Exits on FLEET_COMPLETE.
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
VENVPY=/var/tmp/codex-colab-tools/colab-cli-venv/bin/python
STAGE=/var/tmp/r63_m1_stage/watch
DIGEST=/mnt/d/GI_another/results/round63_m1/watch/M1_WATCH_DIGEST.log
VERBOSE="$STAGE/watch_verbose.log"
CYCLE_S=300
STALL_S=900
mkdir -p "$STAGE" "$(dirname "$DIGEST")"

TAGS=(a b c p1a p1b)
declare -A ACCT=( [a]=pro2 [b]=pro2 [c]=pro2 [p1a]=pro1 [p1b]=pro1 )
declare -A SESS=( [a]=r63m1_a [b]=r63m1_b [c]=r63m1_c [p1a]=r63m1_p1a [p1b]=r63m1_p1b )
declare -A TOTAL=( [a]=11 [b]=8 [c]=8 [p1a]=8 [p1b]=8 )
declare -A EP=(
  [a]=gpu-l4-s-kkb-ass1a0-31p0z6j43zm2y
  [b]=gpu-l4-s-kkb-ass1a0-1qoqoql79dk3e
  [c]=gpu-l4-s-kkb-ass1a1-2qniwt6j3y0pq
  [p1a]=gpu-l4-s-kkb-ass1a2-2ink89ozfh59f
  [p1b]=gpu-l4-s-kkb-usw4a0-64dgi267p0qa
)

cat > "$STAGE/map.json" <<'EOF'
{"gpu-l4-s-kkb-ass1a0-31p0z6j43zm2y": "r63m1_a",
 "gpu-l4-s-kkb-ass1a0-1qoqoql79dk3e": "r63m1_b",
 "gpu-l4-s-kkb-ass1a1-2qniwt6j3y0pq": "r63m1_c",
 "gpu-l4-s-kkb-ass1a2-2ink89ozfh59f": "r63m1_p1a",
 "gpu-l4-s-kkb-usw4a0-64dgi267p0qa": "r63m1_p1b"}
EOF
cat > "$STAGE/rebind.py" <<'PYEOF'
import json, sys
mapping = json.load(open(sys.argv[1]))
from colab_cli.common import state
from colab_cli.state import SessionState
for a in state.client.list_assignments():
    print("FOUND", a.endpoint)
    name = mapping.get(a.endpoint)
    if not name:
        continue
    state.store.add(SessionState(name=name, token=a.runtime_proxy_info.token,
                                 url=a.runtime_proxy_info.url, endpoint=a.endpoint,
                                 variant=a.variant.name, accelerator=a.accelerator.name))
PYEOF
# VM-side ground-truth summarizer: ONE line, counts and rc codes only.
cat > "$STAGE/summ_remote.py" <<'PYEOF'
import glob, json, os, time
B = '/content/round63_m1_bundle'
now = int(time.time())
done = 0; fails = []; inflight = 0
cells = rows = 0; ctot = 0; ages = []
for h in sorted(glob.glob(B + '/results/round63/heartbeats/*.hb.json')):
    try:
        d = json.load(open(h))
    except Exception:
        continue
    rc = d.get("returncode")
    if rc == 0:
        done += 1
    elif rc is None:
        inflight += 1
        cells += int(d.get("cells_done") or 0)
        ctot += int(d.get("cells_total") or 0)
        rows += int(d.get("rows") or 0)
        ages.append(now - int(d.get("epoch", now)))
    else:
        fails.append("%s:rc=%s" % (d.get("shard"), rc))
ncsv = len(glob.glob(B + '/results/round63_m1/shards/*.csv'))
print("M1SUMM done=%d fail=%s inflight=%d cells=%d/%d rows=%d max_hb=%s csv=%d"
      % (done, ";".join(fails) if fails else "-", inflight, cells, ctot, rows,
         max(ages) if ages else -1, ncsv))
PYEOF

log()    { echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) $*" >> "$VERBOSE"; }
digest() { echo "$*" >> "$DIGEST"; }

digest "=== M1 watchdog v2 armed $(date -u +%Y-%m-%dT%H:%M:%SZ) (pid $$, cycle ${CYCLE_S}s, stall ${STALL_S}s; VM ground-truth exec; counts+rc only) ==="
CYCLE=0
while true; do
  CYCLE=$((CYCLE + 1))
  T0=$SECONDS
  # -- rebind + live-endpoint census ------------------------------------------ #
  : > "$STAGE/found.txt"
  for acct in pro1 pro2; do
    HOME="/var/tmp/codex-colab-accounts/$acct" timeout 90 "$VENVPY" \
      "$STAGE/rebind.py" "$STAGE/map.json" >> "$STAGE/found.txt" 2>/dev/null
  done
  # -- per-session ground-truth exec ------------------------------------------- #
  ALERT=0
  LINES=()
  N_COMPLETE=0
  for tag in "${TAGS[@]}"; do
    if ! grep -q "FOUND ${EP[$tag]}" "$STAGE/found.txt"; then
      ALERT=1
      LINES+=("m1_$tag  state=DEAD  endpoint no longer assigned")
      continue
    fi
    OUT=$(HOME="/var/tmp/codex-colab-accounts/${ACCT[$tag]}" timeout 90 "$C" --auth oauth2 \
          exec --session "${SESS[$tag]}" --file "$STAGE/summ_remote.py" --timeout 60 2>&1 \
          | grep "^M1SUMM" | tail -1)
    if [ -z "$OUT" ]; then
      ALERT=1
      LINES+=("m1_$tag  state=TOKEN/CLI  summarizer exec failed after rebind")
      continue
    fi
    done_n=$(echo "$OUT" | sed -nE 's/.*done=([0-9]+).*/\1/p')
    fail_s=$(echo "$OUT" | sed -nE 's/.*fail=([^ ]+).*/\1/p')
    maxhb=$(echo "$OUT"  | sed -nE 's/.*max_hb=(-?[0-9]+).*/\1/p')
    state="OK"
    if [ "${done_n:-0}" -ge "${TOTAL[$tag]}" ]; then
      state="COMPLETE"; N_COMPLETE=$((N_COMPLETE + 1))
    elif [ "${maxhb:--1}" -gt "$STALL_S" ]; then
      state="STALL"; ALERT=1
    fi
    [ "$fail_s" != "-" ] && ALERT=1
    LINES+=("m1_$tag  state=$state  ${OUT#M1SUMM } total=${TOTAL[$tag]}")
  done
  # -- logging ----------------------------------------------------------------- #
  for l in "${LINES[@]}"; do log "cyc$CYCLE $l"; done
  if [ "$ALERT" = 1 ] || [ $(( (CYCLE - 1) % 6 )) = 0 ]; then
    digest "--- $(date -u +%Y-%m-%dT%H:%M:%SZ) cycle $CYCLE $( [ "$ALERT" = 1 ] && echo ALERT || echo digest ) ---"
    for l in "${LINES[@]}"; do digest "  $l"; done
  fi
  if [ "$N_COMPLETE" = "${#TAGS[@]}" ]; then
    digest "=== FLEET_COMPLETE $(date -u +%Y-%m-%dT%H:%M:%SZ) — all ${#TAGS[@]} sessions done==total; watchdog exiting ==="
    log "FLEET_COMPLETE; exiting"
    exit 0
  fi
  EL=$((SECONDS - T0))
  [ "$EL" -lt "$CYCLE_S" ] && sleep $((CYCLE_S - EL))
done
