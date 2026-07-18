#!/bin/bash
# Early confirmatory pipeline step 1: fetch all 32 S2A_DETAIL shard CSVs (+
# meta JSONs) from whichever session ran them, into results/round63/ locally.
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
DEST=/mnt/d/GI_another/results/round63/s2_downloads
mkdir -p "$DEST"
bash /mnt/d/GI_another/code/round63/colab/live_rebind.sh >/dev/null 2>&1
declare -A ACCT=( [pro1a]=pro1 [pro1b]=pro1 [pro2a]=pro2 [pro2b]=pro2 [pro2c]=pro2 )
declare -A SESS=( [pro1a]=r63_pro1a [pro1b]=r63_pro1b [pro2a]=r63_pro2a [pro2b]=r63_pro2b [pro2c]=r63_pro2c )
STAGE=/var/tmp/r63_s2_stage/status_once

# which session completed which shard (from the freshly pulled status JSONs)
ok=0; fail=0
for tag in pro1a pro1b pro2a pro2b pro2c; do
  f="$STAGE/$tag.json"; rm -f "$f"
  HOME="/var/tmp/codex-colab-accounts/${ACCT[$tag]}" timeout 240 "$C" --auth oauth2 \
    download --session "${SESS[$tag]}" \
    "/content/round63_s2_bundle/results/round63/status/$tag.json" "$f" >/dev/null 2>&1
  [ -s "$f" ] || { echo "STATUS_MISS $tag"; continue; }
  for sh in $(python3 -c "
import json
d=json.load(open('$STAGE/$tag.json'))
print(' '.join(c['shard'] for c in d.get('completed',[])
               if c.get('rc')==0 and str(c.get('shard','')).startswith('S2A_DETAIL')))"); do
    for fn in "${sh}.csv" "${sh}_meta.json"; do
      out="$DEST/$fn"
      [ -s "$out" ] && continue
      HOME="/var/tmp/codex-colab-accounts/${ACCT[$tag]}" timeout 240 "$C" --auth oauth2 \
        download --session "${SESS[$tag]}" \
        "/content/round63_s2_bundle/results/round63/shards/$fn" "$out" >/dev/null 2>&1
      if [ -s "$out" ]; then ok=$((ok+1)); else echo "DL_FAIL $fn ($tag)"; fail=$((fail+1)); fi
    done
  done
done
echo "== downloaded files ok=$ok fail=$fail"
ls "$DEST" | grep -c '^S2A_DETAIL.*csv$' | xargs echo "detail csvs present:"
