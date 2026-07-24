#!/bin/bash
# Upload the 5 byte-identical frozen/driver .py files + the two VM launchers to both A100
# sessions, then exec the launchers (which place files, verify md5, background-run the tasks).
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
VENVPY=/var/tmp/codex-colab-tools/colab-cli-venv/bin/python
H=/var/tmp/codex-colab-accounts/pro2
RN=/mnt/d/GI_another/results/round63_next

cat > /tmp/cbt_camp_rebind.py <<'PYEOF'
from colab_cli.common import state
from colab_cli.state import SessionState
MAP = {"gpu-a100-s-kkb-ass1b2-xtrulwvrwi44": "cbt_camp1",
       "gpu-a100-s-kkb-ass1c0-3p8qvfl62mp9y": "cbt_camp2"}
try:
    for a in state.client.list_assignments():
        nm = MAP.get(a.endpoint)
        if nm:
            state.store.add(SessionState(name=nm, token=a.runtime_proxy_info.token,
                                         url=a.runtime_proxy_info.url, endpoint=a.endpoint,
                                         variant=a.variant.name, accelerator=a.accelerator.name))
            print("REBOUND", nm)
except Exception as e:
    print("REBIND_ERR", type(e).__name__, str(e)[:120])
PYEOF
HOME=$H timeout 120 "$VENVPY" /tmp/cbt_camp_rebind.py 2>&1

FILES=(
  "$RN/CURVED_BLINDNESS_TEST/curved_blindness_test.py:curved_blindness_test.py"
  "$RN/SCRAMBLE_EXT/scramble_toy.py:scramble_toy.py"
  "$RN/JET_TEST/jet_test.py:jet_test.py"
  "$RN/CBT_LAW_PLANE/cbt_law_plane.py:cbt_law_plane.py"
  "$RN/CBT_PRIVACY_TEST/cbt_privacy_test.py:cbt_privacy_test.py"
)
upload_common () {
  local S=$1
  for pair in "${FILES[@]}"; do
    local src="${pair%%:*}"; local base="${pair##*:}"
    HOME=$H timeout 240 "$C" --auth oauth2 upload --session "$S" "$src" "/content/$base" 2>&1 | tail -1
  done
}

echo "=================== upload cbt_camp1 ==================="
upload_common cbt_camp1
HOME=$H timeout 240 "$C" --auth oauth2 upload --session cbt_camp1 "$RN/CBT_CAMPAIGN/launcher_camp1.py" /content/launcher_camp1.py 2>&1 | tail -1
echo "=================== upload cbt_camp2 ==================="
upload_common cbt_camp2
HOME=$H timeout 240 "$C" --auth oauth2 upload --session cbt_camp2 "$RN/CBT_CAMPAIGN/launcher_camp2.py" /content/launcher_camp2.py 2>&1 | tail -1

echo "=================== exec launcher camp1 ==================="
HOME=$H timeout 300 "$C" --auth oauth2 exec --session cbt_camp1 --file "$RN/CBT_CAMPAIGN/launcher_camp1.py" --timeout 600 2>&1 | tail -20
echo "=================== exec launcher camp2 ==================="
HOME=$H timeout 300 "$C" --auth oauth2 exec --session cbt_camp2 --file "$RN/CBT_CAMPAIGN/launcher_camp2.py" --timeout 600 2>&1 | tail -20
