#!/bin/bash
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
H=/var/tmp/codex-colab-accounts/pro2
BASE=/mnt/d/GI_another/results/round63_next/FOG_DMD_E6
HOME=$H timeout 240 $C --auth oauth2 upload --session e6s1 "$BASE/vm/launch_s1.py" /content/launch_s1.py 2>&1 | tail -1
HOME=$H timeout 240 $C --auth oauth2 upload --session e6s2 "$BASE/vm/launch_s2.py" /content/launch_s2.py 2>&1 | tail -1
echo "== launch s1 (A0,A1,A4) =="
HOME=$H timeout 200 $C --auth oauth2 exec --session e6s1 --file "$BASE/vm/launch_s1.py" --timeout 120 2>&1 | tail -4
echo "== launch s2 (A2,A3) =="
HOME=$H timeout 200 $C --auth oauth2 exec --session e6s2 --file "$BASE/vm/launch_s2.py" --timeout 120 2>&1 | tail -4
