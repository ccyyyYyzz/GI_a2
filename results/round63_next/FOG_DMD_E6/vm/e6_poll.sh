#!/bin/bash
# Rebind then fetch logs + partial JSONs from both sessions. Idempotent.
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
H=/var/tmp/codex-colab-accounts/pro2
BASE=/mnt/d/GI_another/results/round63_next/FOG_DMD_E6
bash "$BASE/vm/e6_rebind.sh" >/dev/null 2>&1
for pair in "e6s1:E6_s1.json:e6_s1.log" "e6s2:E6_s2.json:e6_s2.log"; do
  s="${pair%%:*}"; rest="${pair#*:}"; j="${rest%%:*}"; lg="${rest#*:}"
  HOME=$H timeout 120 $C --auth oauth2 download --session "$s" "/content/$lg" "$BASE/$lg" 2>&1 | tail -1
  HOME=$H timeout 120 $C --auth oauth2 download --session "$s" "/content/$j" "$BASE/$j" 2>&1 | tail -1
done
echo "=== s1 log tail ==="; tail -5 "$BASE/e6_s1.log" 2>/dev/null
echo "=== s2 log tail ==="; tail -5 "$BASE/e6_s2.log" 2>/dev/null
echo "=== DONE markers ==="; grep -l "DONE ->" "$BASE"/e6_s*.log 2>/dev/null
