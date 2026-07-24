#!/bin/bash
# One-shot: rebind both sessions, download both progress logs + list result JSONs present.
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
VENVPY=/var/tmp/codex-colab-tools/colab-cli-venv/bin/python
H=/var/tmp/codex-colab-accounts/pro2
RN=/mnt/d/GI_another/results/round63_next
CAMP=$RN/CBT_CAMPAIGN

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
except Exception as e:
    print("REBIND_ERR", type(e).__name__, str(e)[:120])
PYEOF
HOME=$H timeout 120 "$VENVPY" /tmp/cbt_camp_rebind.py 2>&1

HOME=$H timeout 120 "$C" --auth oauth2 download --session cbt_camp1 \
  /content/results/round63_next/CBT_LAW_PLANE/camp1_progress.log "$CAMP/camp1_progress.log" >/dev/null 2>&1 || true
HOME=$H timeout 120 "$C" --auth oauth2 download --session cbt_camp2 \
  /content/results/round63_next/CBT_PRIVACY_TEST/camp2_progress.log "$CAMP/camp2_progress.log" >/dev/null 2>&1 || true
echo "===== camp1_progress.log ====="; tail -8 "$CAMP/camp1_progress.log" 2>/dev/null || echo "(none)"
echo "===== camp2_progress.log ====="; tail -8 "$CAMP/camp2_progress.log" 2>/dev/null || echo "(none)"
