#!/bin/bash
# Self-healing A100 watcher: each cycle re-binds the CLI session (fresh runtime-proxy token,
# guide 3.2) then downloads run.log + colab_run.out and reports.  Exits on VERDICT or crash.
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
VENVPY=/var/tmp/codex-colab-tools/colab-cli-venv/bin/python
H=/var/tmp/codex-colab-accounts/pro2
S=cbt_replica
EP=gpu-a100-s-kkb-usc1b1-y3aoc1muozy1
L=/mnt/d/GI_another/results/round63_next/CURVED_BLINDNESS_TEST/colab_a100
R=/content/results/round63_next/CURVED_BLINDNESS_TEST

cat > /tmp/cbt_rebind.py <<'PYEOF'
import sys
from colab_cli.common import state
from colab_cli.state import SessionState
MAP = {"gpu-a100-s-kkb-usc1b1-y3aoc1muozy1": "cbt_replica"}
try:
    for a in state.client.list_assignments():
        nm = MAP.get(a.endpoint)
        if not nm:
            continue
        state.store.add(SessionState(name=nm, token=a.runtime_proxy_info.token,
                                     url=a.runtime_proxy_info.url, endpoint=a.endpoint,
                                     variant=a.variant.name, accelerator=a.accelerator.name))
        print("REBOUND", nm)
except Exception as e:
    print("REBIND_ERR", type(e).__name__, str(e)[:100])
PYEOF

for i in $(seq 1 70); do
  HOME=$H timeout 120 "$VENVPY" /tmp/cbt_rebind.py >/dev/null 2>&1 || true
  HOME=$H timeout 120 "$C" --auth oauth2 download --session "$S" "$R/run.log" "$L/run.log" >/dev/null 2>&1 || true
  HOME=$H timeout 120 "$C" --auth oauth2 download --session "$S" "$R/colab_run.out" "$L/colab_run.out" >/dev/null 2>&1 || true
  echo "===[cycle $i $(date +%H:%M:%S)]==="
  tail -3 "$L/run.log" 2>/dev/null
  if grep -q "VERDICT:" "$L/run.log" 2>/dev/null; then echo ">>> A100_DONE_VERDICT"; break; fi
  if grep -qE "Traceback|Error:|Killed|MemoryError|CUDA out of memory" "$L/colab_run.out" 2>/dev/null; then echo ">>> A100_CRASH"; break; fi
  sleep 55
done
echo "===== A100 FINAL run.log ====="
cat "$L/run.log" 2>/dev/null
