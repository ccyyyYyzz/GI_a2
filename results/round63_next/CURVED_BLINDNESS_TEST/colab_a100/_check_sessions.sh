#!/bin/bash
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
for H in pro2 pro1 proplus_web; do
  echo "=== account $H sessions ==="
  HOME=/var/tmp/codex-colab-accounts/$H timeout 120 "$C" --auth oauth2 sessions 2>&1 | head -20
done
echo "=== gpu options (new --help) ==="
HOME=/var/tmp/codex-colab-accounts/pro2 "$C" --auth oauth2 new --help 2>&1 | head -50
