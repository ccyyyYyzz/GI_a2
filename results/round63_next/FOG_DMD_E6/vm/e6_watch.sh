#!/bin/bash
# Watchdog cycle (per COLAB_USAGE_GUIDE): rebind -> ensure keep-alives -> heartbeat check
# -> idempotent fetch. Exits when both sessions print "DONE ->" or after MAX cycles.
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
H=/var/tmp/codex-colab-accounts/pro2
BASE=/mnt/d/GI_another/results/round63_next/FOG_DMD_E6
EP1=gpu-l4-s-kkb-ass1b2-876uvi64s25f    # e6s1
EP2=gpu-l4-s-kkb-ass1b0-32rhhwq9z4gca   # e6s2
MAX=60
ensure_ka() {  # endpoint session
  pgrep -f "keep-alive $1" >/dev/null 2>&1 || \
    setsid bash -c "HOME=$H $C --auth oauth2 keep-alive $1 $2 >/tmp/e6_ka_$2.log 2>&1" &
}
for i in $(seq 1 $MAX); do
  bash "$BASE/vm/e6_rebind.sh" >/dev/null 2>&1        # 1. rebind (self-heal token)
  ensure_ka "$EP1" e6s1; ensure_ka "$EP2" e6s2         # 2. ensure keep-alives (idempotent)
  bash "$BASE/vm/e6_poll.sh" >/tmp/e6_poll_out.txt 2>&1 # 3. idempotent fetch (rebinds again, harmless)
  s1done=$(grep -c "DONE ->" "$BASE/e6_s1.log" 2>/dev/null); s1done=${s1done:-0}
  s2done=$(grep -c "DONE ->" "$BASE/e6_s2.log" 2>/dev/null); s2done=${s2done:-0}
  s1n=$(grep -c "^    \[" "$BASE/e6_s1.log" 2>/dev/null); s1n=${s1n:-0}
  s2n=$(grep -c "^    \[" "$BASE/e6_s2.log" 2>/dev/null); s2n=${s2n:-0}
  echo "[watch $i] s1_cells=$s1n done=$s1done | s2_cells=$s2n done=$s2done | $(date +%H:%M:%S)"
  if [ "$s1done" -ge 1 ] && [ "$s2done" -ge 1 ]; then echo "BOTH_DONE"; break; fi
  sleep 110
done
echo "=== final s1 tail ==="; tail -3 "$BASE/e6_s1.log" 2>/dev/null
echo "=== final s2 tail ==="; tail -3 "$BASE/e6_s2.log" 2>/dev/null
