#!/bin/bash
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
H=/var/tmp/codex-colab-accounts/pro2
for s in e6s1 e6s2; do
  echo "== creating $s on pro2 =="
  HOME=$H timeout 240 $C --auth oauth2 new --session "$s" --gpu L4 2>&1 | tail -4
  echo "== rc=$? =="
done
echo "== sessions =="
HOME=$H timeout 120 $C --auth oauth2 sessions 2>&1 | tail -20
echo "== endpoints =="
VPY=/var/tmp/codex-colab-tools/colab-cli-venv/bin/python
cat > /tmp/e6_eps.py <<'PYEOF'
from colab_cli.common import state
for x in state.client.list_assignments():
    print("EP", x.endpoint, x.variant.name, x.accelerator.name)
PYEOF
HOME=$H timeout 120 $VPY /tmp/e6_eps.py 2>&1 | tail
