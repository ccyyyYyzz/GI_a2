#!/bin/bash
# Release the two finished pro1 assignments (results fetched, 82/82 local)
# and kill their keep-alive daemons.
set -u
VENVPY=/var/tmp/codex-colab-tools/colab-cli-venv/bin/python
cat > /tmp/r63_release_p1.py <<'PYEOF'
from colab_cli.common import state
targets = {"gpu-l4-s-kkb-ass1a1-37hszw4iqnyz6",   # r63_pro1a
           "gpu-l4-s-kkb-ass1b0-2vtuz2jdjpg73"}   # r63_pro1b
for a in state.client.list_assignments():
    if a.endpoint in targets:
        try:
            state.client.unassign(a.endpoint)
            print("RELEASED", a.endpoint)
        except Exception as e:
            print("RELEASE_FAIL", a.endpoint, type(e).__name__, str(e)[:120])
print("remaining:", [x.endpoint for x in state.client.list_assignments()])
PYEOF
HOME=/var/tmp/codex-colab-accounts/pro1 timeout 240 "$VENVPY" /tmp/r63_release_p1.py
for ep in gpu-l4-s-kkb-ass1a1-37hszw4iqnyz6 gpu-l4-s-kkb-ass1b0-2vtuz2jdjpg73; do
  pkill -f "keep-alive $ep" 2>/dev/null && echo "keepalive killed: $ep" || echo "no keepalive: $ep"
done
