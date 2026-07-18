#!/bin/bash
# One-shot S2 status pull (rebind first, then download + summarize all five).
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
STAGE=/var/tmp/r63_s2_stage/status_once
mkdir -p "$STAGE"
bash /mnt/d/GI_another/code/round63/colab/live_rebind.sh >/dev/null 2>&1
declare -A ACCT=( [pro1a]=pro1 [pro1b]=pro1 [pro2a]=pro2 [pro2b]=pro2 [pro2c]=pro2 )
declare -A SESS=( [pro1a]=r63_pro1a [pro1b]=r63_pro1b [pro2a]=r63_pro2a [pro2b]=r63_pro2b [pro2c]=r63_pro2c )
TOT_DONE=0; TOT_ALL=0
for tag in pro1a pro1b pro2a pro2b pro2c; do
  f="$STAGE/$tag.json"; rm -f "$f"
  HOME="/var/tmp/codex-colab-accounts/${ACCT[$tag]}" timeout 240 "$C" --auth oauth2 \
    download --session "${SESS[$tag]}" \
    "/content/round63_s2_bundle/results/round63/status/$tag.json" "$f" >/dev/null 2>&1
  if [ -s "$f" ]; then
    python3 - "$f" "$tag" <<'PYEOF'
import json, sys
d = json.load(open(sys.argv[1]))
lanes = d.get("lanes", [])
busy = sum(1 for l in lanes if l.get("state") == "running")
cells = sum(l.get("cells_done") or 0 for l in lanes)
rows = sum(l.get("rows") or 0 for l in lanes)
ages = [l.get("heartbeat_age_s", 0) for l in lanes]
cur = [l.get("shard") for l in lanes if l.get("state") == "running"]
print(f"{sys.argv[2]}: {d.get('completed_count')}/{d.get('total_shards')} shards done, "
      f"{busy} busy, in-flight {cells} cells/{rows} rows, max_hb {max(ages) if ages else -1}s")
print(f"   running: {', '.join(cur[:6])}")
PYEOF
  else
    echo "$tag: STATUS DOWNLOAD FAILED"
  fi
done
