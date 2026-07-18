#!/bin/bash
# Fetch ALL completed shard CSVs+metas from every session (idempotent).
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
DEST=/mnt/d/GI_another/results/round63/s2_downloads
mkdir -p "$DEST"
bash /mnt/d/GI_another/code/round63/colab/live_rebind.sh >/dev/null 2>&1
declare -A ACCT=( [pro1a]=pro1 [pro1b]=pro1 [pro2a]=pro2 [pro2b]=pro2 [pro2c]=pro2 )
declare -A SESS=( [pro1a]=r63_pro1a [pro1b]=r63_pro1b [pro2a]=r63_pro2a [pro2b]=r63_pro2b [pro2c]=r63_pro2c )
STAGE=/var/tmp/r63_s2_stage/status_once
ok=0; fail=0; skip=0
for tag in pro1a pro1b pro2a pro2b pro2c; do
  f="$STAGE/$tag.json"; rm -f "$f"
  HOME="/var/tmp/codex-colab-accounts/${ACCT[$tag]}" timeout 240 "$C" --auth oauth2 \
    download --session "${SESS[$tag]}" \
    "/content/round63_s2_bundle/results/round63/status/$tag.json" "$f" >/dev/null 2>&1
  [ -s "$f" ] || { echo "STATUS_MISS $tag"; continue; }
  for sh in $(python3 -c "
import json
d=json.load(open('$STAGE/$tag.json'))
print(' '.join(c['shard'] for c in d.get('completed',[]) if c.get('rc')==0))"); do
    for fn in "${sh}.csv" "${sh}_meta.json"; do
      out="$DEST/$fn"
      [ -s "$out" ] && { skip=$((skip+1)); continue; }
      HOME="/var/tmp/codex-colab-accounts/${ACCT[$tag]}" timeout 240 "$C" --auth oauth2 \
        download --session "${SESS[$tag]}" \
        "/content/round63_s2_bundle/results/round63/shards/$fn" "$out" >/dev/null 2>&1
      if [ -s "$out" ]; then ok=$((ok+1)); else echo "DL_FAIL $fn ($tag)"; fail=$((fail+1)); fi
    done
  done
done
echo "== fetched ok=$ok skip=$skip fail=$fail; local shard csvs: $(ls "$DEST" | grep -c 'csv$')"
