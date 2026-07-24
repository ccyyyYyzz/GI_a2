#!/bin/bash
# Self-healing watchdog (runs INSIDE wsl): each cycle rebinds both sessions (fresh token,
# guide 3.2), downloads both progress logs, and emits an event only on state change /
# completion / crash / periodic heartbeat.  Exits when BOTH campaigns finish or on a crash.
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
VENVPY=/var/tmp/codex-colab-tools/colab-cli-venv/bin/python
H=/var/tmp/codex-colab-accounts/pro2
CAMP=/mnt/d/GI_another/results/round63_next/CBT_CAMPAIGN
L1=$CAMP/camp1_progress.log
L2=$CAMP/camp2_progress.log

cat > /tmp/cbt_camp_rebind.py <<'PYEOF'
from colab_cli.common import state
from colab_cli.state import SessionState
MAP = {"gpu-a100-s-kkb-ass1b2-xtrulwvrwi44": "cbt_camp1",
       "gpu-a100-s-kkb-ass1c0-3p8qvfl62mp9y": "cbt_camp2"}
try:
    for a in state.client.list_assignments():
        nm = MAP.get(a.endpoint)
        if nm:
            state.store.add(SessionState(name=nm, token=a.runtime_proxy_info.token,
                                         url=a.runtime_proxy_info.url, endpoint=a.endpoint,
                                         variant=a.variant.name, accelerator=a.accelerator.name))
except Exception:
    pass
PYEOF

count () { local n; n=$(grep -c "$1" "$2" 2>/dev/null); echo "${n:-0}"; }
prev_n1=-1; prev_c2=-1
for i in $(seq 1 55); do
  HOME=$H timeout 90 "$VENVPY" /tmp/cbt_camp_rebind.py >/dev/null 2>&1 || true
  HOME=$H timeout 90 "$C" --auth oauth2 download --session cbt_camp1 \
    /content/results/round63_next/CBT_LAW_PLANE/camp1_progress.log "$L1" >/dev/null 2>&1 || true
  HOME=$H timeout 90 "$C" --auth oauth2 download --session cbt_camp2 \
    /content/results/round63_next/CBT_PRIVACY_TEST/camp2_progress.log "$L2" >/dev/null 2>&1 || true

  n1=$(count "wrote .*mc_corner" "$L1")
  c1=$(count "CAMP1_ALL_DONE" "$L1")
  c2=$(count "CAMP2_ALL_DONE" "$L2")
  crash=$(grep -lE "Traceback|CUDA out of memory|MemoryError|Killed|Segmentation fault" "$L1" "$L2" 2>/dev/null | wc -l)

  if [ "$crash" -gt 0 ]; then
    echo "CBT_CAMPAIGN_CRASH :: camp1='$(tail -2 "$L1"|tr '\n' ' ')' camp2='$(tail -2 "$L2"|tr '\n' ' ')'"
    break
  fi
  if [ "$n1" != "$prev_n1" ] || [ "$c2" != "$prev_c2" ]; then
    echo "progress[cyc $i]: camp1 corners_done=$n1/4 all=$c1 | camp2 done=$c2 last='$(tail -1 "$L2" 2>/dev/null)'"
    prev_n1=$n1; prev_c2=$c2
  elif [ $((i % 6)) -eq 0 ]; then
    echo "heartbeat[cyc $i]: camp1 corners=$n1/4 done=$c1 | camp2 done=$c2"
  fi
  if [ "$c1" -gt 0 ] && [ "$c2" -gt 0 ] 2>/dev/null; then
    echo "CBT_CAMPAIGN_BOTH_DONE :: camp1 corners=$n1/4 | $(grep VERDICT "$L2" | tail -1)"
    break
  fi
  sleep 70
done
echo "===== camp1 tail ====="; tail -4 "$L1" 2>/dev/null
echo "===== camp2 tail ====="; tail -4 "$L2" 2>/dev/null
