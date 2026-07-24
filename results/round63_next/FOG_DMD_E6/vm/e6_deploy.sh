#!/bin/bash
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
H=/var/tmp/codex-colab-accounts/pro2
BASE=/mnt/d/GI_another/results/round63_next/FOG_DMD_E6
up() {  # session
  local s="$1"
  echo "== upload to $s =="
  for f in fog_e6.py run_e6.py; do
    HOME=$H timeout 240 $C --auth oauth2 upload --session "$s" "$BASE/$f" "/content/$f" 2>&1 | tail -1
  done
  HOME=$H timeout 240 $C --auth oauth2 upload --session "$s" "$BASE/vm/sanity.py" "/content/sanity.py" 2>&1 | tail -1
  echo "== sanity on $s =="
  HOME=$H timeout 240 $C --auth oauth2 exec --session "$s" --file "$BASE/vm/sanity.py" --timeout 300 2>&1 | tail -6
}
up e6s1
up e6s2
