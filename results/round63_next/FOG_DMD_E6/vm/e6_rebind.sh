#!/bin/bash
# Rebind e6s1/e6s2 CLI rows to live runtimes with fresh tokens (token expires ~30-45min;
# CLI 401 path prunes rows instead of refreshing). pro2 only.
set -u
VPY=/var/tmp/codex-colab-tools/colab-cli-venv/bin/python
cat > /tmp/e6_map.json <<'EOF'
{"gpu-l4-s-kkb-ass1b2-876uvi64s25f": "e6s1",
 "gpu-l4-s-kkb-ass1b0-32rhhwq9z4gca": "e6s2"}
EOF
cat > /tmp/e6_rebind.py <<'PYEOF'
import json
mapping = json.load(open('/tmp/e6_map.json'))
from colab_cli.common import state
from colab_cli.state import SessionState
n = 0
for a in state.client.list_assignments():
    name = mapping.get(a.endpoint)
    if not name:
        continue
    s = SessionState(name=name, token=a.runtime_proxy_info.token,
                     url=a.runtime_proxy_info.url, endpoint=a.endpoint,
                     variant=a.variant.name, accelerator=a.accelerator.name)
    state.store.add(s); print('REBOUND', name, a.endpoint); n += 1
print('TOTAL_REBOUND', n)
PYEOF
HOME=/var/tmp/codex-colab-accounts/pro2 timeout 240 "$VPY" /tmp/e6_rebind.py 2>&1 | grep -E "REBOUND|TOTAL|Error|error" | head
