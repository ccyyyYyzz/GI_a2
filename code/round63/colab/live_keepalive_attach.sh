#!/bin/bash
# (Re)attach keep-alive daemons for all five sessions; idempotent (skips if a
# matching keep-alive process already runs).
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
STAGE=/var/tmp/r63_s2_stage
declare -A EP=(
  [r63_pro1a]=gpu-l4-s-kkb-ass1b0-2vtuz2jdjpg73
  [r63_pro1b]=gpu-l4-s-kkb-ass1a1-37hszw4iqnyz6
  [r63_pro2a]=gpu-l4-s-kkb-ass1b1-2jga36pyizwwq
  [r63_pro2b]=gpu-l4-s-kkb-ass1a2-2liz0l1bktulc
  [r63_pro2c]=gpu-l4-s-kkb-usw4c2-2qxckgl655b7r
)
declare -A ACCT=(
  [r63_pro1a]=pro1 [r63_pro1b]=pro1
  [r63_pro2a]=pro2 [r63_pro2b]=pro2 [r63_pro2c]=pro2
)
for s in r63_pro1a r63_pro1b r63_pro2a r63_pro2b r63_pro2c; do
  ep="${EP[$s]}"
  if pgrep -f "keep-alive $ep $s" >/dev/null; then
    echo "$s: keep-alive already running"
    continue
  fi
  HOME="/var/tmp/codex-colab-accounts/${ACCT[$s]}" setsid nohup timeout 86400 \
    "$C" --auth oauth2 keep-alive "$ep" "$s" \
    >"$STAGE/ka_$s.log" 2>&1 < /dev/null &
  sleep 2
  pgrep -f "keep-alive $ep $s" >/dev/null \
    && echo "$s: keep-alive attached" || echo "$s: KEEPALIVE_FAIL"
done
