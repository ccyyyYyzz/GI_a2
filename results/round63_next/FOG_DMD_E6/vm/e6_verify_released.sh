#!/bin/bash
set -u
VPY=/var/tmp/codex-colab-tools/colab-cli-venv/bin/python
H=/var/tmp/codex-colab-accounts/pro2
cat > /tmp/e6_chk.py <<'PYEOF'
from colab_cli.common import state
a = state.client.list_assignments()
print("N_ASSIGN", len(a))
for x in a: print("EP", x.endpoint)
PYEOF
HOME=$H timeout 120 "$VPY" /tmp/e6_chk.py 2>&1 | tail
echo "=== stray keep-alive procs ==="
pgrep -af keep-alive 2>/dev/null || echo "none"
