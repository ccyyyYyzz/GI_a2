#!/bin/bash
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
H=/var/tmp/codex-colab-accounts/pro2
SESS=cbt_replica
for GPU in A100 H100 L4; do
  echo "=== trying --gpu $GPU ==="
  OUT=$(HOME=$H timeout 240 "$C" --auth oauth2 new --session "$SESS" --gpu "$GPU" 2>&1)
  echo "$OUT"
  if echo "$OUT" | grep -qiE "created|ready|session|assigned|endpoint"; then
    if ! echo "$OUT" | grep -qiE "error|unavailable|not available|failed|no runtime|quota"; then
      echo ">>> SUCCESS_GPU=$GPU"
      break
    fi
  fi
  echo ">>> $GPU not obtained, trying next"
done
echo "=== sessions now ==="
HOME=$H timeout 120 "$C" --auth oauth2 sessions 2>&1 | head -20
