#!/bin/bash
# Fetch the A100 full-run artifacts (JSON + 2 figures) that were never downloaded.
# Rebind first (fresh token, guide 3.2), then idempotent downloads into colab_a100/.
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
VENVPY=/var/tmp/codex-colab-tools/colab-cli-venv/bin/python
H=/var/tmp/codex-colab-accounts/pro2
S=cbt_replica
L=/mnt/d/GI_another/results/round63_next/CURVED_BLINDNESS_TEST/colab_a100
R=/content/results/round63_next/CURVED_BLINDNESS_TEST

cat > /tmp/cbt_rebind.py <<'PYEOF'
from colab_cli.common import state
from colab_cli.state import SessionState
MAP = {"gpu-a100-s-kkb-usc1b1-y3aoc1muozy1": "cbt_replica"}
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
HOME=$H timeout 120 "$VENVPY" /tmp/cbt_rebind.py 2>&1

for f in CURVED_BLINDNESS_TEST.json cbt_slopes_coef_kink.png cbt_cusum_model2.png; do
  echo "=== download $f ==="
  HOME=$H timeout 180 "$C" --auth oauth2 download --session "$S" "$R/$f" "$L/$f" 2>&1 | tail -2
done
echo "===== local colab_a100 listing ====="
ls -la "$L"
