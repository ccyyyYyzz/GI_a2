#!/bin/bash
# Rebind the pro2 CLI session (fresh runtime-proxy token per guide 3.2), then report
# whether the A100 endpoint is still assigned and the cbt_replica session is alive.
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
VENVPY=/var/tmp/codex-colab-tools/colab-cli-venv/bin/python
H=/var/tmp/codex-colab-accounts/pro2
S=cbt_replica

cat > /tmp/cbt_rebind.py <<'PYEOF'
from colab_cli.common import state
from colab_cli.state import SessionState
MAP = {"gpu-a100-s-kkb-usc1b1-y3aoc1muozy1": "cbt_replica"}
print("=== list_assignments ===")
try:
    asn = list(state.client.list_assignments())
    if not asn:
        print("NO_ASSIGNMENTS (VM likely reclaimed)")
    for a in asn:
        print("ASSIGN endpoint=%s variant=%s accel=%s" % (a.endpoint, a.variant.name, a.accelerator.name))
        nm = MAP.get(a.endpoint)
        if nm:
            state.store.add(SessionState(name=nm, token=a.runtime_proxy_info.token,
                                         url=a.runtime_proxy_info.url, endpoint=a.endpoint,
                                         variant=a.variant.name, accelerator=a.accelerator.name))
            print("REBOUND", nm)
except Exception as e:
    print("REBIND_ERR", type(e).__name__, str(e)[:200])
PYEOF

HOME=$H timeout 120 "$VENVPY" /tmp/cbt_rebind.py 2>&1
echo "=== sessions ==="
HOME=$H timeout 120 "$C" --auth oauth2 sessions 2>&1 | head -20
