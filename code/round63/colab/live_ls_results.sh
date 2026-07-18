#!/bin/bash
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
cat > /tmp/r63_ls.py <<'PYEOF'
import os
root = '/content/round63_s2_bundle/results/round63'
for dirpath, dirnames, filenames in os.walk(root):
    rel = os.path.relpath(dirpath, root)
    if rel.count(os.sep) > 1:
        dirnames[:] = []
        continue
    print('DIR', rel, '->', sorted(filenames)[:6], '...' if len(filenames) > 6 else '')
PYEOF
HOME=/var/tmp/codex-colab-accounts/pro1 timeout 240 "$C" --auth oauth2 \
  exec --session r63_pro1a --file /tmp/r63_ls.py --timeout 120 2>&1 | grep '^DIR' | head -15
