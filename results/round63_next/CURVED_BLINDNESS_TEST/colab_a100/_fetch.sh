#!/bin/bash
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
H=/var/tmp/codex-colab-accounts/pro2
S=cbt_replica
L=/mnt/d/GI_another/results/round63_next/CURVED_BLINDNESS_TEST/colab_a100
R=/content/results/round63_next/CURVED_BLINDNESS_TEST
HOME=$H timeout 120 "$C" --auth oauth2 download --session "$S" "$R/run.log" "$L/run.log" 2>&1 | tail -1 || true
echo "===== VM run.log ====="
cat "$L/run.log" 2>/dev/null || echo "(no run.log yet)"
