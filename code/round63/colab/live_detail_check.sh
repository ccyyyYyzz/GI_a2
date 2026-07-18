#!/bin/bash
# Count completed S2A_DETAIL shards across the five sessions (from the status
# JSONs already pulled by the watchdog cycle; re-pull fresh here).
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
STAGE=/var/tmp/r63_s2_stage/status_once
mkdir -p "$STAGE"
bash /mnt/d/GI_another/code/round63/colab/live_rebind.sh >/dev/null 2>&1
declare -A ACCT=( [pro1a]=pro1 [pro1b]=pro1 [pro2a]=pro2 [pro2b]=pro2 [pro2c]=pro2 )
declare -A SESS=( [pro1a]=r63_pro1a [pro1b]=r63_pro1b [pro2a]=r63_pro2a [pro2b]=r63_pro2b [pro2c]=r63_pro2c )
: > /tmp/r63_completed.txt
for tag in pro1a pro1b pro2a pro2b pro2c; do
  f="$STAGE/$tag.json"; rm -f "$f"
  HOME="/var/tmp/codex-colab-accounts/${ACCT[$tag]}" timeout 240 "$C" --auth oauth2 \
    download --session "${SESS[$tag]}" \
    "/content/round63_s2_bundle/results/round63/status/$tag.json" "$f" >/dev/null 2>&1
  [ -s "$f" ] && python3 -c "
import json,sys
d=json.load(open('$f'))
for s in d.get('completed',[]): print(s)
for l in d.get('lanes',[]):
    if l.get('state')=='running': print('RUNNING::'+str(l.get('shard')))
" >> /tmp/r63_completed.txt
done
echo "== completed S2A_DETAIL: $(grep -c '^S2A_DETAIL' /tmp/r63_completed.txt) / 32"
echo "== running  S2A_DETAIL: $(grep -c '^RUNNING::S2A_DETAIL' /tmp/r63_completed.txt)"
echo "== completed total: $(grep -vc '^RUNNING::' /tmp/r63_completed.txt) / 82"
grep '^RUNNING::S2A_DETAIL' /tmp/r63_completed.txt | sed 's/RUNNING:://' | tr '\n' ' '; echo
