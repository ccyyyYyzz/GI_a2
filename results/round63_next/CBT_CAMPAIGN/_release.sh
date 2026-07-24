#!/bin/bash
# Release both campaign A100 assignments (results fetched) and kill any keep-alive daemons.
set -u
VENVPY=/var/tmp/codex-colab-tools/colab-cli-venv/bin/python
EPS="gpu-a100-s-kkb-ass1b2-xtrulwvrwi44 gpu-a100-s-kkb-ass1c0-3p8qvfl62mp9y"
cat > /tmp/cbt_camp_release.py <<'PYEOF'
from colab_cli.common import state
targets = {"gpu-a100-s-kkb-ass1b2-xtrulwvrwi44",
           "gpu-a100-s-kkb-ass1c0-3p8qvfl62mp9y"}
seen = set()
for a in state.client.list_assignments():
    if a.endpoint in targets:
        seen.add(a.endpoint)
        try:
            state.client.unassign(a.endpoint); print("RELEASED", a.endpoint)
        except Exception as e:
            print("RELEASE_FAIL", a.endpoint, type(e).__name__, str(e)[:120])
for t in targets - seen:
    print("ALREADY_GONE", t)
print("remaining:", [x.endpoint for x in state.client.list_assignments()])
PYEOF
HOME=/var/tmp/codex-colab-accounts/pro2 timeout 240 "$VENVPY" /tmp/cbt_camp_release.py
for ep in $EPS; do pkill -f "keep-alive $ep" 2>/dev/null && echo "keepalive killed: $ep" || echo "no keepalive: $ep"; done
