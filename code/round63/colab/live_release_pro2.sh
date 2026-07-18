#!/bin/bash
# Release the three finished pro2 assignments (results already fetched) and
# kill their keep-alive daemons. pro1 sessions untouched.
set -u
VENVPY=/var/tmp/codex-colab-tools/colab-cli-venv/bin/python
cat > /tmp/r63_release.py <<'PYEOF'
from colab_cli.common import state
targets = {"gpu-l4-s-kkb-ass1b1-2jga36pyizwwq",   # r63_pro2a
           "gpu-l4-s-kkb-ass1a2-2liz0l1bktulc",   # r63_pro2b
           "gpu-l4-s-kkb-usw4c2-2qxckgl655b7r"}   # r63_pro2c
for a in state.client.list_assignments():
    if a.endpoint in targets:
        try:
            state.client.unassign(a.endpoint)
            print("RELEASED", a.endpoint)
        except Exception as e:
            print("RELEASE_FAIL", a.endpoint, type(e).__name__, str(e)[:120])
print("remaining:", [x.endpoint for x in state.client.list_assignments()])
PYEOF
HOME=/var/tmp/codex-colab-accounts/pro2 timeout 240 "$VENVPY" /tmp/r63_release.py
for ep in gpu-l4-s-kkb-ass1b1-2jga36pyizwwq gpu-l4-s-kkb-ass1a2-2liz0l1bktulc gpu-l4-s-kkb-usw4c2-2qxckgl655b7r; do
  pkill -f "keep-alive $ep" 2>/dev/null && echo "keepalive killed: $ep" || echo "no keepalive: $ep"
done
