#!/bin/bash
# Release both e6 VMs (unassign) + kill keep-alives. Run when E6 sweep is complete.
set -u
VPY=/var/tmp/codex-colab-tools/colab-cli-venv/bin/python
H=/var/tmp/codex-colab-accounts/pro2
cat > /tmp/e6_release.py <<'PYEOF'
from colab_cli.common import state
targets = {"gpu-l4-s-kkb-ass1b2-876uvi64s25f",   # e6s1
           "gpu-l4-s-kkb-ass1b0-32rhhwq9z4gca"}  # e6s2
for a in state.client.list_assignments():
    if a.endpoint in targets:
        try:
            state.client.unassign(a.endpoint); print("RELEASED", a.endpoint)
        except Exception as e:
            print("RELEASE_FAIL", a.endpoint, type(e).__name__, str(e)[:120])
print("remaining:", [x.endpoint for x in state.client.list_assignments()])
PYEOF
HOME=$H timeout 240 "$VPY" /tmp/e6_release.py
for ep in gpu-l4-s-kkb-ass1b2-876uvi64s25f gpu-l4-s-kkb-ass1b0-32rhhwq9z4gca; do
  pkill -f "keep-alive $ep" 2>/dev/null && echo "keepalive killed: $ep" || echo "no keepalive: $ep"
done
