#!/bin/bash
# Create 2 pro2 sessions for the CBT calibration-plane (JOB 2) + privacy (JOB 3) campaigns.
# Try A100 first (preferred), fall back to L4.  Print endpoints for rebind/keepalive/release.
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
H=/var/tmp/codex-colab-accounts/pro2
for S in cbt_camp1 cbt_camp2; do
  echo "=================== create $S ==================="
  for GPU in A100 L4; do
    echo "--- try --gpu $GPU ---"
    OUT=$(HOME=$H timeout 240 "$C" --auth oauth2 new --session "$S" --gpu "$GPU" 2>&1)
    echo "$OUT" | tail -6
    if echo "$OUT" | grep -qiE "created|ready|assigned|endpoint|Hardware"; then
      if ! echo "$OUT" | grep -qiE "error|unavailable|not available|failed|no runtime|quota"; then
        echo ">>> GOT $GPU for $S"; break
      fi
    fi
    echo ">>> $GPU not obtained for $S"
  done
done
echo "=================== sessions ==================="
HOME=$H timeout 120 "$C" --auth oauth2 sessions 2>&1 | head -20
