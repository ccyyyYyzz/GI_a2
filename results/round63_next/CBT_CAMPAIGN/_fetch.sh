#!/bin/bash
# Idempotent fetch: rebind both, download all result JSONs + progress logs into the job dirs.
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
VENVPY=/var/tmp/codex-colab-tools/colab-cli-venv/bin/python
H=/var/tmp/codex-colab-accounts/pro2
RN=/mnt/d/GI_another/results/round63_next
LP=$RN/CBT_LAW_PLANE
PV=$RN/CBT_PRIVACY_TEST
RLP=/content/results/round63_next/CBT_LAW_PLANE
RPV=/content/results/round63_next/CBT_PRIVACY_TEST

cat > /tmp/cbt_camp_rebind.py <<'PYEOF'
from colab_cli.common import state
from colab_cli.state import SessionState
MAP = {"gpu-a100-s-kkb-ass1b2-xtrulwvrwi44": "cbt_camp1",
       "gpu-a100-s-kkb-ass1c0-3p8qvfl62mp9y": "cbt_camp2"}
for a in state.client.list_assignments():
    nm = MAP.get(a.endpoint)
    if nm:
        state.store.add(SessionState(name=nm, token=a.runtime_proxy_info.token,
                                     url=a.runtime_proxy_info.url, endpoint=a.endpoint,
                                     variant=a.variant.name, accelerator=a.accelerator.name))
PYEOF
HOME=$H timeout 120 "$VENVPY" /tmp/cbt_camp_rebind.py 2>&1 || true

dl () { HOME=$H timeout 150 "$C" --auth oauth2 download --session "$1" "$2" "$3" 2>&1 | tail -1; }
echo "=== camp1 (CBT_LAW_PLANE) ==="
dl cbt_camp1 "$RLP/CBT_LAW_PLANE_analytic.json" "$LP/CBT_LAW_PLANE_analytic.json"
for k in 0 1 2 3; do dl cbt_camp1 "$RLP/CBT_LAW_PLANE_mc_corner$k.json" "$LP/CBT_LAW_PLANE_mc_corner$k.json"; done
dl cbt_camp1 "$RLP/camp1_progress.log" "$LP/camp1_progress.log"
echo "=== camp2 (CBT_PRIVACY_TEST) ==="
dl cbt_camp2 "$RPV/CBT_PRIVACY_TEST.json" "$PV/CBT_PRIVACY_TEST.json"
dl cbt_camp2 "$RPV/camp2_progress.log" "$PV/camp2_progress.log"
echo "=== local listing ==="
ls -la "$LP"/CBT_LAW_PLANE_*.json "$PV"/CBT_PRIVACY_TEST.json 2>/dev/null
