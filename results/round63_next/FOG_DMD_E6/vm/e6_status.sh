#!/bin/bash
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
VPY=/var/tmp/codex-colab-tools/colab-cli-venv/bin/python
H=/var/tmp/codex-colab-accounts/pro2
echo "=== pro2 sessions ==="
HOME=$H timeout 120 $C --auth oauth2 sessions 2>&1 | tail -20
echo "=== pro2 live assignments ==="
cat > /tmp/e6_assign.py <<'PYEOF'
from colab_cli.common import state
try:
    a = state.client.list_assignments()
    print("N_ASSIGN", len(a))
    for x in a:
        print("EP", x.endpoint, x.variant.name, x.accelerator.name)
except Exception as e:
    print("ERR", type(e).__name__, str(e)[:200])
PYEOF
HOME=$H timeout 120 $VPY /tmp/e6_assign.py 2>&1 | tail -20
