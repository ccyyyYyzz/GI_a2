#!/bin/bash
# ROUND63 S2 standing watchdog: poll the five session status JSONs every 15 min,
# print compact progress lines, flag stalls (heartbeat age > 600 s) and dead
# sessions (download failure twice in a row).
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
STAGE=/var/tmp/r63_s2_stage/status
mkdir -p "$STAGE"
declare -A ACCT=( [pro1a]=pro1 [pro1b]=pro1 [pro2a]=pro2 [pro2b]=pro2 [pro2c]=pro2 )
declare -A SESS=( [pro1a]=r63_pro1a [pro1b]=r63_pro1b [pro2a]=r63_pro2a [pro2b]=r63_pro2b [pro2c]=r63_pro2c )
declare -A MISS=( [pro1a]=0 [pro1b]=0 [pro2a]=0 [pro2b]=0 [pro2c]=0 )

LAST_DIGEST=0
REBIND=/mnt/d/GI_another/code/round63/colab/live_rebind.sh
while true; do
  # self-heal: refresh runtime-proxy tokens every cycle (they expire ~30-45
  # min and the CLI's 401 path prunes rows instead of refreshing)
  bash "$REBIND" >/dev/null 2>&1
  all_done=1
  CYCLE_OUT=""
  ANOMALY=0
  for tag in pro1a pro1b pro2a pro2b pro2c; do
    a="${ACCT[$tag]}"; s="${SESS[$tag]}"
    f="$STAGE/$tag.json"
    rm -f "$f"
    HOME="/var/tmp/codex-colab-accounts/$a" timeout 240 "$C" --auth oauth2 \
      download --session "$s" \
      "/content/round63_s2_bundle/results/round63/status/$tag.json" "$f" \
      >/dev/null 2>&1
    if [ ! -s "$f" ]; then
      MISS[$tag]=$(( ${MISS[$tag]} + 1 ))
      if [ "${MISS[$tag]}" -ge 2 ]; then
        CYCLE_OUT="$CYCLE_OUT | $tag UNREACHABLEx${MISS[$tag]}"
        ANOMALY=1
      fi
      all_done=0
      continue
    fi
    MISS[$tag]=0
    line=$(python3 - "$f" "$tag" <<'PYEOF'
import json, sys, time
d = json.load(open(sys.argv[1]))
tag = sys.argv[2]
done = d.get("completed_count", "?"); tot = d.get("total_shards", "?")
lanes = d.get("lanes", [])
ages = [l.get("heartbeat_age_s") for l in lanes if l.get("heartbeat_age_s") is not None]
busy = sum(1 for l in lanes if l.get("state") == "running")
cells = sum(l.get("cells_done") or 0 for l in lanes)
mx = max(ages) if ages else -1
stall = " STALL" if (ages and mx > 600) else ""
fin = " ALL_DONE" if (isinstance(done, int) and isinstance(tot, int) and done >= tot and busy == 0) else ""
print(f"S2WATCH {tag}: {done}/{tot} shards, {busy} busy (+{cells} cells in-flight), max_hb {mx:.0f}s{stall}{fin}")
PYEOF
)
    short=$(echo "$line" | sed 's/^S2WATCH //')
    CYCLE_OUT="$CYCLE_OUT | $short"
    echo "$line" | grep -qE "STALL" && ANOMALY=1
    echo "$line" | grep -q "ALL_DONE" || all_done=0
  done
  now=$(date +%s)
  if [ "$all_done" -eq 1 ]; then
    echo "S2WATCH: ALL SESSIONS REPORT COMPLETE$CYCLE_OUT"
    break
  fi
  if [ "$ANOMALY" -eq 1 ] || [ $(( now - LAST_DIGEST )) -ge 3600 ]; then
    echo "S2WATCH_DIGEST $(date -u +%H:%MZ):$CYCLE_OUT"
    LAST_DIGEST=$now
  fi
  sleep 900
done
