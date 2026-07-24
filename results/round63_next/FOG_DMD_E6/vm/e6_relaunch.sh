#!/bin/bash
# Kill stale drivers, re-upload FIXED code, relaunch fresh on both sessions.
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
H=/var/tmp/codex-colab-accounts/pro2
BASE=/mnt/d/GI_another/results/round63_next/FOG_DMD_E6
bash "$BASE/vm/e6_rebind.sh" >/dev/null 2>&1
for s in e6s1 e6s2; do
  echo "== redeploy $s =="
  HOME=$H timeout 240 $C --auth oauth2 upload --session "$s" "$BASE/fog_e6.py" /content/fog_e6.py 2>&1 | tail -1
  HOME=$H timeout 240 $C --auth oauth2 upload --session "$s" "$BASE/run_e6.py" /content/run_e6.py 2>&1 | tail -1
  HOME=$H timeout 240 $C --auth oauth2 upload --session "$s" "$BASE/vm/kill.py" /content/kill.py 2>&1 | tail -1
  echo "-- kill stale on $s --"
  HOME=$H timeout 200 $C --auth oauth2 exec --session "$s" --file "$BASE/vm/kill.py" --timeout 60 2>&1 | tail -2
done
echo "== relaunch s1 (A0,A1,A4) =="
HOME=$H timeout 200 $C --auth oauth2 exec --session e6s1 --file "$BASE/vm/launch_s1.py" --timeout 120 2>&1 | tail -3
echo "== relaunch s2 (A2,A3) =="
HOME=$H timeout 200 $C --auth oauth2 exec --session e6s2 --file "$BASE/vm/launch_s2.py" --timeout 120 2>&1 | tail -3
