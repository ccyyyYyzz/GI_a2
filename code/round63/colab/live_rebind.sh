#!/bin/bash
# Re-bind CLI session names to still-assigned runtimes with FRESH tokens from
# list_assignments (the runtime-proxy token expires ~30-45 min; the CLI's
# 401 path prunes rows instead of refreshing — this repairs them).
set -u
VENVPY=/var/tmp/codex-colab-tools/colab-cli-venv/bin/python

cat > /tmp/r63_map.json <<'EOF'
{"gpu-l4-s-kkb-ass1b0-2vtuz2jdjpg73": "r63_pro1a",
 "gpu-l4-s-kkb-ass1a1-37hszw4iqnyz6": "r63_pro1b",
 "gpu-l4-s-kkb-ass1b1-2jga36pyizwwq": "r63_pro2a",
 "gpu-l4-s-kkb-ass1a2-2liz0l1bktulc": "r63_pro2b",
 "gpu-l4-s-kkb-usw4c2-2qxckgl655b7r": "r63_pro2c"}
EOF

cat > /tmp/r63_rebind.py <<'PYEOF'
import json, sys
mapping = json.load(open('/tmp/r63_map.json'))
from colab_cli.common import state
from colab_cli.state import SessionState
assignments = state.client.list_assignments()
n = 0
for a in assignments:
    name = mapping.get(a.endpoint)
    if not name:
        print('SKIP unknown endpoint', a.endpoint)
        continue
    s = SessionState(
        name=name,
        token=a.runtime_proxy_info.token,
        url=a.runtime_proxy_info.url,
        endpoint=a.endpoint,
        variant=a.variant.name,          # int-enum: .value is 1, .name 'GPU'
        accelerator=a.accelerator.name,
    )
    state.store.add(s)
    print('REBOUND', name, a.endpoint)
    n += 1
print('TOTAL_REBOUND', n)
PYEOF

for acct in pro1 pro2; do
  echo "== rebind $acct"
  HOME="/var/tmp/codex-colab-accounts/$acct" timeout 240 "$VENVPY" /tmp/r63_rebind.py 2>&1 | grep -E "REBOUND|SKIP|TOTAL|Error|error" | head -8
done
